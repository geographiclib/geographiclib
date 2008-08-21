/**
 * \file TransverseMercator.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/TransverseMercator.hpp"
#include "GeographicLib/Constants.hpp"
#include <limits>

/*
 * Implementation taken from the report:
 *
 * JHS 154, ETRS89 - järjestelmään liittyvät karttaprojektiot,
 * tasokoordinaatistot ja karttalehtijako (Map projections, plane coordinates,
 * and map sheet index for ETRS89), Published by JUHTA, Finnish Geodetic
 * Institute, and the National Land Survey of Finland, 34 p (2006).
 *
 * http://www.jhs-suositukset.fi/suomi/jhs154
 * http://docs.jhs-suositukset.fi/jhs-suositukset/JHS154/JHS154.pdf
 *
 * This achieves approximately 1micron accuracy.
 *
 * Other accurate implementations are given in
 * http://www.ign.fr/telechargement/MPro/geodesie/CIRCE/NTG_76.pdf
 * http://www.lantmateriet.se/upload/filer/kartor/geodesi_gps_och_detaljmatning/geodesi/Formelsamling/Gauss_Conformal_Projection.pdf
 */

/*
 *
 * Transformation errors are primarly a function of x and rapidly increase for
 * x > 70e5
 *
 * For x < 70e5
 *
 * Error < 1mm for forward tx: F(ll) - F_e(ll) - this method - exact
 * Error < 1mm for round trip: ll-F^-1(F(ll)) converted to a distance
 *
 * For reverse tx restrict inputs to |x| < 70e5
 * For forward tx restrict inputs to (90-lam)^2 + phi^2 > 36^2 and
 * require output |x| < 7e6
 *
 * Max relative error in scale is 16%
 * Max relative in convergence is 4deg
 *
 * UTM subsets:
 *
 * dx = max |F(ll)-F_e(ll)| + |ll-F^-1(F(ll))| expressed as a distance
 * dang = max discrepancy in meridan convergence
 * dscale = max relative error in scale
 * set                       dx     dang  dscale
 * x<5e5, y<96e5             0.38um 1.07' 0.08%
 * x<4e5, y<95e5             0.36um 0.27' 0.03%
 * x<4e5, phi<84.5, lam<12   0.36um 0.35" 0.022%%
 * x<4e5, phi<84.5, lam<9    0.36um 0.10" 0.016%%
 * 
 */


namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = TRANSVERSEMERCATOR_HPP;

}

namespace GeographicLib {

  TransverseMercator::TransverseMercator(double a, double invf, double k0)
    : _a(a)
    , _f(1 / invf)
    , _k0(k0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(_e2))
    , _e2m(1 - _e2)
    , _e1(sqrt(_e2 / _e2m))
    , _n(_f / (2 - _f))
      // _a1 is the equivalent radius for computing the circumference of ellipse.
      // Relative error is f^6/16384 = 8.8e-20 for WGS84.
    , _a1(_a / (1 + _n) * (_n * _n * (_n * _n + 16) + 64) / 64)
    , _h1(_n * (_n * ((555 - 4 * _n) * _n - 960) + 720) / 1440)
    , _h2(_n * _n * ((96 - 437 * _n) * _n + 30) / 1440)
    , _h3((119 - 148 * _n) * _n * _n * _n / 3360)
    , _h4(4397 * _n * _n * _n * _n / 161280)
    , _h1p(_n * (_n * (_n * (164 * _n + 225) - 480) + 360) / 720)
    , _h2p(_n * _n * (_n * (557 * _n - 864) + 390) / 1440)
    , _h3p((427 - 1236 * _n) * _n * _n * _n / 1680)
    , _h4p(49561 * _n * _n * _n * _n / 161280)
    , _tol(0.1*sqrt(std::numeric_limits<double>::epsilon()))
    , _numit(5)
  {}

  const TransverseMercator
  TransverseMercator::UTM(Constants::WGS84_a, Constants::WGS84_invf,
			  Constants::UTM_k0);

  double TransverseMercator::Convergence(double phi, double l) const {
    double
      cosphi = cos(phi),
      v2 = _e1 * cosphi;
    v2 = 1 + v2 * v2;
    return l * sin(phi) * (1 + v2 * (2 * v2 - 1) * cosphi * cosphi * l * l / 3);
  }

  double TransverseMercator::Scale(double phi, double l) const {
    double cosphi = cos(phi);
    return _k0 * ( 1 + cosphi * cosphi * l * l / 2);
  }
  void TransverseMercator::Forward(double lon0, double lat, double lon,
				   double& x, double& y,
				   double& gamma, double& k) const {
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon - lon0 > 180)
      lon -= lon0 - 360;
    else if (lon - lon0 <= -180)
      lon -= lon0 + 360;
    else
      lon -= lon0;
    // Now lon in (-180, 180]
    // Explicitly enforce the parity
    int
      latsign = lat < 0 ? -1 : 1,
      lonsign = lon < 0 ? -1 : 1;
    lon *= lonsign;
    lat *= latsign;
    bool backside = lon > 90;
    if (backside) {
      if (lat == 0)
	latsign = -1;
      lon = 180 - lon;
    }
    double
      phi = lat * Constants::degree,
      l = lon * Constants::degree;
    double etap, xip;
    if (lat < 90) {
      double
	qp = asinh(tan(phi)),
	qpp = atanh(_e * sin(phi)),
	q = qp - _e * qpp,
	bet = atan(sinh(q));
      etap = cos(bet) * sin(l);
      etap = atanh(etap);
      if (lon < 90)
	xip = asin(sin(bet) * cosh(etap));
      else
	xip = Constants::pi/2;
    } else {
      xip = Constants::pi/2;
      etap = 0;
    }

    double
      xi1  = _h1p * sin(2 * xip) * cosh(2 * etap),
      xi2	 = _h2p * sin(4 * xip) * cosh(4 * etap),
      xi3	 = _h3p * sin(6 * xip) * cosh(6 * etap),
      xi4	 = _h4p * sin(8 * xip) * cosh(8 * etap),
      eta1 = _h1p * cos(2 * xip) * sinh(2 * etap),
      eta2 = _h2p * cos(4 * xip) * sinh(4 * etap),
      eta3 = _h3p * cos(6 * xip) * sinh(6 * etap),
      eta4 = _h4p * cos(8 * xip) * sinh(8 * etap),
      xi = xip + xi1 + xi2 + xi3 + xi4,
      eta = etap + eta1 + eta2 + eta3 + eta4;

    y = _a1 * _k0 * (backside ? Constants::pi - xi : xi) * latsign;
    x = _a1 * _k0 * eta * lonsign;
    gamma = Convergence(phi, l) / Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= latsign * lonsign;
    k = Scale(phi, l);
  }

  void TransverseMercator::Reverse(double lon0, double x, double y,
				   double& lat, double& lon,
				   double& gamma, double& k) const {
    double
      xi = y / (_a1 * _k0),
      eta = x / (_a1 * _k0);
    // Explicitly enforce the parity
    int
      xisign = xi < 0 ? -1 : 1,
      etasign = eta < 0 ? -1 : 1;
    xi *= xisign;
    eta *= etasign;
    bool backside = xi > Constants::pi/2;
    if (backside) {
      xi = Constants::pi - xi;
    }
    
    double
      xi1p  = _h1 * sin(2 * xi) * cosh(2 * eta),
      xi2p  = _h2 * sin(4 * xi) * cosh(4 * eta),
      xi3p  = _h3 * sin(6 * xi) * cosh(6 * eta),
      xi4p  = _h4 * sin(8 * xi) * cosh(8 * eta),
      eta1p = _h1 * cos(2 * xi) * sinh(2 * eta),
      eta2p = _h2 * cos(4 * xi) * sinh(4 * eta),
      eta3p = _h3 * cos(6 * xi) * sinh(6 * eta),
      eta4p = _h4 * cos(8 * xi) * sinh(8 * eta),
      xip = xi - (xi1p + xi2p + xi3p + xi4p),
      etap = eta - (eta1p + eta2p + eta3p + eta4p),
      bet = asin(sin(xip) / cosh(etap));
    double l, phi;
    if (bet < Constants::pi/2) {
      l = asin(tanh(etap) / cos(bet));
      double
	q = asinh(tan(bet)),
	qp = q;
      for (int i = 0; i < _numit; ++i) {
	double
	  t = tanh(qp),
	  dqp = -(qp - _e * atanh(_e * t) - q) *
	  (1 - _e2 * t * t) / _e2m;
	qp += dqp;
	if (std::abs(dqp) < _tol)
	  break;
      }
      phi = atan(sinh(qp));
    } else {
      phi = Constants::pi/2;
      l = 0;
    }
    lat = phi / Constants::degree * xisign;
    lon = l / Constants::degree;
    if (backside)
      lon = 180 - lon;
    lon *= etasign;
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon + lon0 >= 180)
      lon += lon0 - 360;
    else if (lon + lon0 < -180)
      lon += lon0 + 360;
    else
      lon += lon0;
    gamma = Convergence(phi, l) / Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= xisign * etasign;
    k = Scale(phi, l);
  }

} // namespace GeographicLib

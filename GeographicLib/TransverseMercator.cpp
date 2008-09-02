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
 * This achieves approximately 1um accuracy.  This is a straight transcription
 * of the formulas in this paper with the following exceptions:
 *
 *  * use Newton's method instead of plain iteration to solve for latitude in
 *    terms of isometric latitude in the Reverse method
 *
 *  * use more accurate expressions for the convergence and scale.
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
 * Error < 1.4mm for forward tx: F(ll) - F_e(ll) - this method - exact
 * Error < 0.9mm for round trip: ll-F^-1(F(ll)) converted to a distance
 *
 * For reverse tx restrict inputs to |x| < 70e5
 * For forward tx restrict inputs to (90-lam)^2 + phi^2 > 36^2 and
 * require output |x| < 70e5
 *
 * Max relative error in scale is 0.08%%
 * Max relative in convergence is 0.23"
 *
 * UTM subsets:
 *
 * dx = max |F(ll)-F_e(ll)| + |ll-F^-1(F(ll))| expressed as a distance
 * dang = max discrepancy in meridan convergence
 * dscale = max relative error in scale
 * set                       dx     dang  dscale
 * x<5e5, y<96e5             0.22um 9e-8" 1e-6%%
 * x<4e5, y<95e5             0.21um 3e-8" 4e-7%%
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
    , _ep2(_e2 / _e2m)
    , _n(_f / (2 - _f))
      // _a1 is the equivalent radius for computing the circumference of
      // ellipse.  Relative error is f^6/16384 = 8.8e-20 for WGS84.
    , _a1(_a / (1 + _n) * (_n * _n * (_n * _n + 16) + 64) / 64)
    , _h1(_n * (_n * ((555 - 4 * _n) * _n - 960) + 720) / 1440)
    , _h2(_n * _n * ((96 - 437 * _n) * _n + 30) / 1440)
    , _h3((119 - 148 * _n) * _n * _n * _n / 3360)
    , _h4(4397 * _n * _n * _n * _n / 161280)
    , _h1p(_n * (_n * (_n * (164 * _n + 225) - 480) + 360) / 720)
    , _h2p(_n * _n * (_n * (557 * _n - 864) + 390) / 1440)
    , _h3p((427 - 1236 * _n) * _n * _n * _n / 1680)
    , _h4p(49561 * _n * _n * _n * _n / 161280)
#if 0
      // Here are order _n^8 terms with ^ indicating exponentiation.  These
      // need to be converted to Horner form, but are here left in expanded
      // form so that they can be easily truncated to lower order in _n.
    , _a1( _a*(1+_n^2/4+_n^4/64+_n^6/256+25*_n^8/16384)/(1+_n) )
    , _h1( _n/2-2*_n^2/3+37*_n^3/96-_n^4/360-81*_n^5/512+96199*_n^6/604800-
	  5406467*_n^7/38707200+7944359*_n^8/67737600 )
    , _h2( _n^2/48+_n^3/15-437*_n^4/1440+46*_n^5/105-1118711*_n^6/3870720+
	  51841*_n^7/1209600+24749483*_n^8/348364800 )
    , _h3( 17*_n^3/480-37*_n^4/840-209*_n^5/4480+5569*_n^6/90720+9261899*
	  _n^7/58060800-6457463*_n^8/17740800 )
    , _h4( 4397*_n^4/161280-11*_n^5/504-830251*_n^6/7257600+
	   466511*_n^7/2494800+324154477*_n^8/7664025600 )
    , _h5( 4583*_n^5/161280-108847*_n^6/3991680-8005831*_n^7/63866880+
	  22894433*_n^8/124540416 )
    , _h6( 20648693*_n^6/638668800-16363163*_n^7/518918400-
	  2204645983*_n^8/12915302400 )
    , _h7( 219941297*_n^7/5535129600-497323811*_n^8/12454041600 )
    , _h8( 191773887257*_n^8/3719607091200 )
    , _h1p( _n/2-2*_n^2/3+5*_n^3/16+41*_n^4/180-127*_n^5/288+7891*_n^6/37800+
	   72161*_n^7/387072-18975107*_n^8/50803200 )
    , _h2p( 13*_n^2/48-3*_n^3/5+557*_n^4/1440+281*_n^5/630-
	    1983433*_n^6/1935360+13769*_n^7/28800+148003883*_n^8/174182400 )
    , _h3p( 61*_n^3/240-103*_n^4/140+15061*_n^5/26880+167603*_n^6/181440-
	   67102379*_n^7/29030400+79682431*_n^8/79833600 )
    , _h4p( 49561*_n^4/161280-179*_n^5/168+6601661*_n^6/7257600+
	    97445*_n^7/49896-40176129013*_n^8/7664025600 )
    , _h5p( 34729*_n^5/80640-3418889*_n^6/1995840+14644087*_n^7/9123840+
	   2605413599*_n^8/622702080 )
    , _h6p( 212378941*_n^6/319334400-30705481*_n^7/10378368+
	   175214326799*_n^8/58118860800 )
    , _h7p( 1522256789*_n^7/1383782400-16759934899*_n^8/3113510400 )
    , _h8p( 1424729850961*_n^8/743921418240 )
#endif
    , _tol(0.1*sqrt(std::numeric_limits<double>::epsilon()))
    , _numit(5)
  {}

  const TransverseMercator
  TransverseMercator::UTM(Constants::WGS84_a, Constants::WGS84_invf,
			  Constants::UTM_k0);

  void TransverseMercator::Scale(double phi, double l,
				 double& gamma, double& k) const {
    double
      c = cos(phi),
      s = sin(l),
      c2 = c * c,
      s2 = s * s,
      d = 1 - s2 * c2,
      // Accurate to order _ep2^2
      carg = 1 + c2 * c2 * s2 / d * _ep2 *
      (1 + c2 / (3 * d * d) *
       (2 - s2 * (c2 * ((4 * c2 - 1) * s2 - 9) + 8)) * _ep2),
      // Accurate to order _ep2
      cabs = 1 + c2 * c2 * s2 * ((c2 - 2) * s2 + 1) / (2 * d * d) * _ep2;
#if 0
    // To order _ep2^4
    carg = 1+c2^2*s2*_ep2/d-
      ((4*c2^5-c2^4)*s2^3+(-9*c2^4+8*c2^3)*s2^2-2*c2^3*s2)*_ep2^2/(3*d^3)+
      ((28*c2^8-11*c2^7-2*c2^6)*s2^5+(-120*c2^7+75*c2^6+48*c2^5)*s2^4+
       (165*c2^6-336*c2^5+120*c2^4)*s2^3+
       (93*c2^5-60*c2^4)*s2^2)*_ep2^3/(15*d^5)-
      ((836*c2^11-429*c2^10-75*c2^9-17*c2^8)*s2^7+
       (-5264*c2^10+3472*c2^9+1106*c2^8+1264*c2^7)*s2^6+
       (13097*c2^9-13566*c2^8-14460*c2^7+11536*c2^6)*s2^5+
       (-13727*c2^8+52288*c2^7-43064*c2^6+8064*c2^5)*s2^4+
       (-15054*c2^7+20272*c2^6-6048*c2^5)*s2^3+
       (-735*c2^6+504*c2^5)*s2^2)*_ep2^4/(315*d^7);
    // To order _ep2^4
    cabs = 1+((c2^3-2*c2^2)*s2^2+c2^2*s2)*_ep2/(2*d^2)-
      ((3*c2^6+12*c2^5-28*c2^4)*s2^4+(-34*c2^5+124*c2^4-64*c2^3)*s2^3+
       (-61*c2^4+48*c2^3)*s2^2)*_ep2^2/(24*d^4)+
      ((15*c2^9+30*c2^8+140*c2^7-328*c2^6)*s2^6+
       (-147*c2^8-1044*c2^7+4268*c2^6-2688*c2^5)*s2^5+
       (1405*c2^7-8274*c2^6+8480*c2^5-1920*c2^4)*s2^4+
       (3383*c2^6-5280*c2^5+1920*c2^4)*s2^3+
       (280*c2^5-240*c2^4)*s2^2)*_ep2^3/(240*d^6)-
      ((1575*c2^12+2520*c2^11+5880*c2^10+27552*c2^9-64016*c2^8)*s2^8+
       (-17396*c2^11-54680*c2^10-457680*c2^9+1860384*c2^8-1247744*c2^7)*s2^7+
       (114506*c2^10+1432424*c2^9-8037000*c2^8+9297792*c2^7-2874368*c2^6)*s2^6+
       (-1212372*c2^9+10085240*c2^8-16270464*c2^7+
	8397312*c2^6-1032192*c2^5)*s2^5+
       (-3418969*c2^8+8422624*c2^7-6228096*c2^6+1290240*c2^5)*s2^4+
       (-685664*c2^7+985152*c2^6-322560*c2^5)*s2^3)*_ep2^4/(40320*d^8);
#endif
    gamma = atan2(sin(phi) * s * carg, cos(l));
    k = cabs/sqrt(d);
    // Old expressions from JHS 154.  These are Taylor expansions in l as well
    // as ep2.
    // double v2 = 1 + _ep2 * c2;
    // Accurate to _ep2^2 and l^3
    // gamma = l * sin(phi) * (1 + v2 * (2 * v2 - 1) * c2 * l * l / 3)
    // Accurate to _ep2^0 and l^2
    // k = ( 1 + c2 * l * l / 2);
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
      // {q,bet} is {isometric,conformal} latitude
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
    // {xi',eta'} is {northing,easting} for spheroidal transverse mercator (for
    // eta' = 0, xi' = bet). {xi,eta} is {northing,easting} for transverse
    // mercator with constant scale on the central meridian (for eta = 0, xip =
    // rectifying latitude).  Define
    //
    //   zeta = xi + i*eta
    //   zeta' = xi' + i*eta'
    //
    // The conversion from conformal to rectifying latitude can be expresses as
    // a series in _n:
    //
    //   zeta = zeta' + sum(h[j]' * sin(2 * j * zeta'), j = 1..4)
    //
    // where h[j]' = O(_n^j).  The reversion of this series gives
    //
    //   zeta' = zeta - sum(h[j] * sin(2 * j * zeta), j = 1..4)
    //
    // which is used in Reverse.  Extensions of these series to order _n^8 are
    // given above.
    double
      xi1  = _h1p * sin(2 * xip) * cosh(2 * etap),
      xi2  = _h2p * sin(4 * xip) * cosh(4 * etap),
      xi3  = _h3p * sin(6 * xip) * cosh(6 * etap),
      xi4  = _h4p * sin(8 * xip) * cosh(8 * etap),
      eta1 = _h1p * cos(2 * xip) * sinh(2 * etap),
      eta2 = _h2p * cos(4 * xip) * sinh(4 * etap),
      eta3 = _h3p * cos(6 * xip) * sinh(6 * etap),
      eta4 = _h4p * cos(8 * xip) * sinh(8 * etap),
      xi = xip + xi1 + xi2 + xi3 + xi4,
      eta = etap + eta1 + eta2 + eta3 + eta4;

    y = _a1 * _k0 * (backside ? Constants::pi - xi : xi) * latsign;
    x = _a1 * _k0 * eta * lonsign;
    Scale(phi, l, gamma, k);
    gamma /= Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= latsign * lonsign;
    k *= _k0;
  }

  void TransverseMercator::Reverse(double lon0, double x, double y,
				   double& lat, double& lon,
				   double& gamma, double& k) const {
    // This undoes the steps in Forward.  The wrinkles are: (1) Use of the
    // reverted series to express zeta' in terms of zeta. (2) Newton's method
    // to solve for phi in terms of q.
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
    Scale(phi, l, gamma, k);
    gamma /= Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= xisign * etasign;
    k *= _k0;
  }

} // namespace GeographicLib

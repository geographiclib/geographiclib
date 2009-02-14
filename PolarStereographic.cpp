/**
 * \file PolarStereographic.cpp
 * \brief Implementation for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/PolarStereographic.hpp"
#include "GeographicLib/Constants.hpp"
#include <limits>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = POLARSTEREOGRAPHIC_HPP;
}

namespace GeographicLib {

  using namespace std;

  PolarStereographic::PolarStereographic(double a, double invf, double k0)
    throw()
    : _a(a)
    , _f(invf > 0 ? 1 / invf : 0)
    , _k0(k0)
    , _e(sqrt(_f * (2 - _f)))
    , _e2m(1 - sq(_e))
    , _c(sqrt( pow(1 + _e, 1 + _e) * pow(1 - _e, 1 - _e) ))
    , _tol(0.1*sqrt(numeric_limits<double>::epsilon()))
    , _numit(5)
  {}

  const PolarStereographic
  PolarStereographic::UPS(Constants::WGS84_a, Constants::WGS84_invf,
			  Constants::UPS_k0);

  void PolarStereographic::Forward(bool northp, double lat, double lon,
				   double& x, double& y,
				   double& gamma, double& k) const throw() {
    double theta = 90 - (northp ? lat : -lat); //  the colatitude
    double rho;
    theta *= Constants::degree;
    double
      ecos = _e * cos(theta),
      f = pow((1 + ecos)/(1 - ecos), _e/2),
      t2 = 2 * tan(theta/2) * f,	   // Snyder (15-9) (t2 = 2 * t)
      m = sin(theta) / sqrt(1 - sq(ecos)); // Snyder (14-15)
    rho = _a * _k0 * t2 / _c;		   // Snyder (21-33)
    k = m < numeric_limits<double>::epsilon() ? _k0 : rho / (_a * m);
    double
      lam = lon * Constants::degree;
    x = rho * sin(lam);
    y = (northp ? -rho : rho) * cos(lam);
    gamma = northp ? lon : -lon;
  }

  void PolarStereographic::Reverse(bool northp, double x, double y,
				   double& lat, double& lon,
				   double& gamma, double& k) const throw() {
    double
      rho = hypot(x, y),
      t2 = rho * _c / (_a * _k0),
      theta = Constants::pi/2,	// initial estimate of colatitude
      ecos = _e * cos(theta);
    // Solve from theta using Newton's method on Snyder (15-9) which converges
    // more rapidly than the iteration procedure given by Snyder (7-9).  First
    // rewrite as
    // v(theta) = 2 * tan(theta/2) * f - t2 = 0
    for (int i = 0; i < _numit; ++i) {
      double
	f = pow((1 + ecos)/(1 - ecos), _e/2),
	c2 = cos(theta/2),
	v = 2 * tan(theta/2) * f - t2,
	dv = _e2m * f / ((1 - sq(ecos)) * sq(c2)), // dv/dtheta
	dtheta = -v/dv;
      theta += dtheta;
      ecos = _e * cos(theta);
      if (abs(dtheta) < _tol)
	break;
    }
    lat = (northp ? 1 : -1) * (90 - theta / Constants::degree);
    // Result is in [-180, 180).  Assume atan2(0,0) = 0.
    lon = -atan2( -x, northp ? -y : y ) / Constants::degree;
    double m = sin(theta) / sqrt(1 - sq(ecos));
    k = m == 0 ? _k0 : rho / (_a * m);
    gamma = northp ? lon : -lon;
  }

} // namespace GeographicLib

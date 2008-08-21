/**
 * \file PolarStereographic.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/PolarStereographic.hpp"
#include "GeographicLib/Constants.hpp"
#include <cmath>
#include <limits>

/*
 * Implementation taken from the report:
 *
 * J. P.Snyder, Map Projections: A Working Manual,
 * USGS Professional Paper 1395 (1987), pp. 160 - 163
 *
 * http://pubs.er.usgs.gov/usgspubs/pp/pp1395
 */

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = POLARSTEREOGRAPHIC_HPP;
}

namespace GeographicLib {

  PolarStereographic::PolarStereographic(double a, double invf, double k0)
    : _a(a)
    , _f(1 / invf)
    , _k0(k0)
    , _e(sqrt(_f * (2 - _f)))
    , _e2m(1 - _e * _e)
    , _c(sqrt( std::pow(1 + _e, 1 + _e) * std::pow(1 - _e, 1 - _e) ))
    , _tol(0.1*sqrt(std::numeric_limits<double>::epsilon()))
    , _numit(5)
  {}

  const PolarStereographic
  PolarStereographic::UPS(Constants::WGS84_a, Constants::WGS84_invf,
			  Constants::UPS_k0);

  void PolarStereographic::Forward(bool northp, double lat, double lon,
				   double& x, double& y,
				   double& gamma, double& k) const {
    double
      // theta is the colatitude
      theta = 90 - (northp ? lat : -lat);
    double rho;
    theta *= Constants::degree;
    double
      ecos = _e * cos(theta),
      f = std::pow((1 + ecos)/(1 - ecos), _e/2),
      t = 2 * tan(theta/2) * f,
      m = sin(theta) / sqrt(1 - ecos * ecos);
    rho = _a * _k0 * t / _c;
    k = m < std::numeric_limits<double>::epsilon() ? _k0 : rho / (_a * m);
    double
      lam = lon * Constants::degree;
    x = rho * sin(lam);
    y = (northp ? -rho : rho) * cos(lam);
    gamma = lon;
  }

  void PolarStereographic::Reverse(bool northp, double x, double y,
				   double& lat, double& lon,
				   double& gamma, double& k) const {
    double
      rho = hypot(x, y),
      t = rho * _c / (_a * _k0),
      theta = Constants::pi/2,	// initial estimate of colatitude
      ecos = _e * cos(theta);
    for (int i = 0; i < _numit; ++i) {
      double
	f = std::pow((1 + ecos)/(1 - ecos), _e/2),
	c2 = cos(theta/2),
	v = 2 * tan(theta/2) * f - t,
	dv = _e2m * f / ((1 - ecos * ecos) * c2 * c2),
	dtheta = -v/dv;
      theta += dtheta;
      ecos = _e * cos(theta);
      if (std::abs(dtheta) < _tol)
	break;
    }
    lat = (northp ? 1 : -1) * (90 - theta / Constants::degree);
    // Result is in [-180, 180)
    lon = rho == 0 ? 0 : -atan2( -x, northp ? -y : y ) / Constants::degree;
    double
      m = sin(theta) / sqrt(1 - ecos * ecos);
    k = m == 0 ? _k0 : rho / (_a * m);
    gamma = lon;
  }

} // namespace GeographicLib

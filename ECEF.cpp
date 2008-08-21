/**
 * \file ECEF.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/ECEF.hpp"
#include "GeographicLib/Constants.hpp"
#include <cmath>
#include <limits>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = ECEF_HPP;
}

namespace GeographicLib {

  ECEF::ECEF(double a, double invf)
    : _a(a)
    , _f(1 / invf)
    , _e2(_f * (2 - _f))
    , _e12(_e2 / (1 - _e2))
    , _b(_a * (1 - _f))
    , _tol(_a * std::numeric_limits<double>::epsilon())
    , _numit(5)
  {}

  const ECEF ECEF::WGS84(Constants::WGS84_a, Constants::WGS84_invf);

  void ECEF::Forward(double lat, double lon, double h,
		     double& x, double& y, double& z) const {
    double
      phi = lat * Constants::degree,
      lam = lon * Constants::degree,
      sphi = sin(phi),
      n = _a/sqrt(1 - _e2 * sphi * sphi);
    z = ( (1 - _f) * (1 - _f) * n + h) * sphi;
    x = (n + h) * cos(phi);
    y = x * sin(lam);
    x *= cos(lam);
  }

  void ECEF::Reverse(double x, double y, double z,
		     double& lat, double& lon, double& h) const {
    double p = hypot(x, y);
    // Use Pollard's latitude first method, J. Pollard, Iterative vector methods
    // for computing geodetic latitude and height from rectangular coordinates,
    // J. Geodesy 76, 36-40 (2002).  The advantage of this method is that it's
    // simple and convergence (via Newton's method) is fast.
    if (p != 0) {
      double z0 = _b * z / hypot(p, z);
      for (int i = 0; i < _numit; ++i) {
	double
	  t = _b / _a * (z + _e12 * z0) / p,
	  tsq = 1/hypot(1.0, t),
	  dz0 = - (_b * t / tsq - z0) / 
	  (_b * _b * _e12/(_a * p  * tsq * tsq * tsq) - 1);
	z0 = z0 + dz0;
	if (std::abs(dz0) < _tol)
	  break;
      }
      double
	phi = atan2(z + _e12, p),
	sphi = sin(phi);
      h = p * cos(phi) + z * sphi - _a * sqrt(1 - _e2 * sphi * sphi);
      lat = phi / Constants::degree;
      lon =  -atan2(-y, x) / Constants::degree;
    } else {
      lon = 0;
      lat = z >= 0 ? 90 : -90;
      h = std::abs(z) - _a * (1 - _f);
    }
  }

} // namespace GeographicLib

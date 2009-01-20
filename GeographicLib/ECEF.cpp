/**
 * \file ECEF.cpp
 * \brief Implementation for GeographicLib::ECEF class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/ECEF.hpp"
#include "GeographicLib/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = ECEF_HPP;
}

#if defined(_MSC_VER)
#define hypot _hypot
#endif

namespace GeographicLib {

  ECEF::ECEF(double a, double invf)
    : _a(a)
    , _f(1 / invf)
    , _e2(_f * (2 - _f))
    , _e12(_e2 / (1 - _e2))
    , _b(_a * (1 - _f))
    , _tol(_a * std::numeric_limits<double>::epsilon())
    , _maxrad(2 * _a / std::numeric_limits<double>::epsilon())
    , _numit(5)
  {}

  const ECEF ECEF::WGS84(Constants::WGS84_a, Constants::WGS84_invf);

  void ECEF::Forward(double lat, double lon, double h,
		     double& x, double& y, double& z) const {
    double
      phi = lat * Constants::degree,
      lam = lon * Constants::degree,
      sphi = std::sin(phi),
      n = _a/std::sqrt(1 - _e2 * sq(sphi));
    z = ( sq(1 - _f) * n + h) * sphi;
    x = (n + h) * std::cos(phi);
    y = x * std::sin(lam);
    x *= std::cos(lam);
  }

  void ECEF::Reverse(double x, double y, double z,
		     double& lat, double& lon, double& h) const {
    // H. Vermeille, Direct transformation from geocentric coordinates to
    // geodetic coordinates, J. Geodesy 76, 451-454 (2002).
    double
      rad = hypot(x, y),
      p = sq(rad / _a),
      q = (1 - _e2) * sq(z / _a),
      e4 = sq(_e2),
      r = (p + q - e4) / 6;
    double phi;
    h = hypot(rad, z);		// Distance to center of earth
    if (h > _maxrad) {
      // We really far away (> 12 million light years); treat the earth as a
      // point and h, above, is an acceptable approximation to the height.
      // This avoids overflow, e.g., in the computation of disc below.  It's
      // possible that h has overflowed to inf; but that's OK.
      //
      // Treat the case x, y finite, but rad overflows to +inf.  Fix by scaling
      // by _maxrad.
      double rad1 = hypot(x/_maxrad, y/_maxrad);
      phi = atan2(z/_maxrad, rad1);
    }
    else if ( !(e4 * q == 0 && r <= 0) ) {
      double 
	// Avoid possible division by zero when r = 0 by multiplying equations
	// for s and t by r^3 and r, resp.
	S = e4 * p * q / 4,	// S = r^3 * s
	r2 = sq(r),
	r3 = r * r2,
	disc =  S * (2 * r3 + S);
      double u = r;
      if (disc >= 0) {
	double T3 = r3 + S;
	// Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
	// of precision due to cancellation.  The result is unchanged because
	// of the way the T is used in definition of u.
	T3 += T3 < 0e0 ? -std::sqrt(disc) : std::sqrt(disc); // T3 = (r * t)^3
	// N.B. cbrt on negative args returns negative real cube root
	double T = cbrt(T3);	// T = r * t
	// T can be zero; but then r2 / T -> 0.
	u += T + (T != 0 ? r2 / T : 0);
      } else {
	// T is complex, but the way u is defined the result is real.
	double ang = std::atan2(std::sqrt(-disc), r3 + S);
	// There are three possible real solutions for u depending on the
	// multiple of 2*pi here.  We choose multiplier = 1 which leads to a
	// jump in the solution across the line 2 + s = 0; but this
	// nevertheless leads to a continuous (and accurate) solution for k.
	// Other choices of the multiplier lead to poorly conditioned solutions
	// near s = 0 (i.e., near p = 0 or q = 0).
	u += 2 * std::abs(r) * std::cos((2 * Constants::pi + ang) / 3.0);
      }
      double
	v = std::sqrt(sq(u) + e4 * q), // guaranteed positive
	// Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
	// e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
	uv = u < 0 ? e4 * q / (v - u) : u + v, //  u + v, guaranteed positive
	// Need to guard against w going negative due to roundoff in uv - q.
	w = std::max(0.0, _e2 * (uv - q) / (2 * v)),
	// Rearrange expression for k to avoid loss of accuracy due to
	// subtraction.  Division by 0 not possible because uv > 0, w >= 0.
	k = uv / (std::sqrt(uv + sq(w)) + w), // guaranteed positive
	d = k * rad / (k + _e2);
      // Prefer atan2 to Vermeille's half-angle formula.
      phi = std::atan2(z, d);
      h = (k + _e2 - 1) * hypot(d, z) / k;
    } else {			// e4 * q == 0 && r <= 0
      // Very near equatorial plane with rad <= a * e^2.  This leads to k = 0
      // using the general formula and division by 0 in formula for h.  So
      // handle this case directly.  The condition e4 * q == 0 implies abs(z) <
      // 1.e-145 for WGS84 so it's OK to treat these points as though z = 0.
      // (But we do take care that the sign of phi matches the sign of z.)
      phi = std::atan2(std::sqrt( -6 * r), std::sqrt(p * (1 - _e2)));
      if (z < 0) phi = -phi;	// for tiny negative z
      h = _a * (_e2 - 1) / std::sqrt(1 - _e2 * sq(std::sin(phi)));
    }
    lat = phi / Constants::degree;
    // Negative signs return lon in [-180, 180)
    lon = (rad != 0 ? -std::atan2(-y, x) : 0) / Constants::degree;
  }

} // namespace GeographicLib


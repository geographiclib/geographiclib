/**
 * \file Geocentric.cpp
 * \brief Implementation for GeographicLib::Geocentric class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/Geocentric.hpp"
#include "GeographicLib/Constants.hpp"
#include <algorithm>
#include <limits>

#define GEOCENTRIC_CPP "$Id$"

RCSID_DECL(GEOCENTRIC_CPP)
RCSID_DECL(GEOCENTRIC_HPP)

namespace GeographicLib {

  using namespace std;

  Geocentric::Geocentric(double a, double r)
    throw()
    : _a(a)
    , _f(r > 0 ? 1 / r : 0)
    , _e2(_f * (2 - _f))
    , _e4(sq(_e2))
    , _e2m(1 - _e2)
    , _maxrad(2 * _a / numeric_limits<double>::epsilon())
  {}

  const Geocentric Geocentric::WGS84(Constants::WGS84_a(),
				     Constants::WGS84_r());

  void Geocentric::Forward(double lat, double lon, double h,
			   double& x, double& y, double& z) const throw() {
    double
      phi = lat * Constants::degree(),
      lam = lon * Constants::degree(),
      sphi = sin(phi),
      n = _a/sqrt(1 - _e2 * sq(sphi));
    z = ( sq(1 - _f) * n + h) * sphi;
    x = (n + h) * cos(phi);
    y = x * sin(lam);
    x *= cos(lam);
  }

  void Geocentric::Reverse(double x, double y, double z,
			   double& lat, double& lon, double& h) const throw() {
    double rad = hypot(x, y);
    h = hypot(rad, z);		// Distance to center of earth
    double phi;
    if (h > _maxrad)
      // We really far away (> 12 million light years); treat the earth as a
      // point and h, above, is an acceptable approximation to the height.
      // This avoids overflow, e.g., in the computation of disc below.  It's
      // possible that h has overflowed to inf; but that's OK.
      //
      // Treat the case x, y finite, but rad overflows to +inf by scaling by 2.
      phi = atan2(z/2, hypot(x/2, y/2));
    else if (_e4 == 0) {
      // Treat the spherical case.  Dealing with underflow in the general case
      // with _e2 = 0 is difficult.  Origin maps to N pole same as an
      // ellipsoid.
      phi = atan2(h != 0 ? z : 1.0, rad);
      h -= _a;
    } else {
      double
	p = sq(rad / _a),
	q = _e2m * sq(z / _a),
	r = (p + q - _e4) / 6;
      if ( !(_e4 * q == 0 && r <= 0) ) {
	double
	  // Avoid possible division by zero when r = 0 by multiplying
	  // equations for s and t by r^3 and r, resp.
	  S = _e4 * p * q / 4,	// S = r^3 * s
	  r2 = sq(r),
	  r3 = r * r2,
	  disc =  S * (2 * r3 + S);
	double u = r;
	if (disc >= 0) {
	  double T3 = r3 + S;
	  // Pick the sign on the sqrt to maximize abs(T3).  This minimizes
	  // loss of precision due to cancellation.  The result is unchanged
	  // because of the way the T is used in definition of u.
	  T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
	  // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
	  double T = cbrt(T3);	// T = r * t
	  // T can be zero; but then r2 / T -> 0.
	  u += T + (T != 0 ? r2 / T : 0);
	} else {
	  // T is complex, but the way u is defined the result is real.
	  double ang = atan2(sqrt(-disc), r3 + S);
	  // There are three possible real solutions for u depending on the
	  // multiple of 2*pi here.  We choose multiplier = 1 which leads to a
	  // jump in the solution across the line 2 + s = 0; but this
	  // nevertheless leads to a continuous (and accurate) solution for k.
	  // Other choices of the multiplier lead to poorly conditioned
	  // solutions near s = 0 (i.e., near p = 0 or q = 0).
	  u += 2 * abs(r) * cos((2 * Constants::pi() + ang) / 3.0);
	}
	double
	  v = sqrt(sq(u) + _e4 * q), // guaranteed positive
	  // Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
	  // e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
	  uv = u < 0 ? _e4 * q / (v - u) : u + v, //  u+v, guaranteed positive
	  // Need to guard against w going negative due to roundoff in uv - q.
	  w = max(0.0, _e2 * (uv - q) / (2 * v)),
	  // Rearrange expression for k to avoid loss of accuracy due to
	  // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
	  k = uv / (sqrt(uv + sq(w)) + w), // guaranteed positive
	  d = k * rad / (k + _e2);
	// Probably atan2 returns the result for phi more accurately than the
	// half-angle formula that Vermeille uses.  It's certainly simpler.
	phi = atan2(z, d);
	h = (k + _e2 - 1) * hypot(d, z) / k;
      } else {			// e4 * q == 0 && r <= 0
	// Very near equatorial plane with rad <= a * e^2.  This leads to k = 0
	// using the general formula and division by 0 in formula for h.  So
	// handle this case directly.  The condition e4 * q == 0 implies abs(z)
	// < 1.e-145 for WGS84 so it's OK to treat these points as though z =
	// 0.  (But we do take care that the sign of phi matches the sign of
	// z.)
	phi = atan2(sqrt( -6 * r), sqrt(p * _e2m));
	if (z < 0) phi = -phi;	// for tiny negative z
	h = - _a * _e2m / sqrt(1 - _e2 * sq(sin(phi)));
      }
    }
    lat = phi / Constants::degree();
    // Negative signs return lon in [-180, 180).  Assume atan2(0,0) = 0.
    lon = -atan2(-y, x) / Constants::degree();
  }

} // namespace GeographicLib


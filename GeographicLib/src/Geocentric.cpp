/**
 * \file Geocentric.cpp
 * \brief Implementation for GeographicLib::Geocentric class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/Geocentric.hpp"
#include <algorithm>
#include <limits>

#define GEOGRAPHICLIB_GEOCENTRIC_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GEOCENTRIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_GEOCENTRIC_HPP)

namespace GeographicLib {

  using namespace std;

  Geocentric::Geocentric(real a, real r) throw()
    : _a(a)
    , _f(r != 0 ? 1 / r : 0)
    , _e2(_f * (2 - _f))
    , _e2m(sq(1 - _f))          // 1 - _e2
      // Constants with the x suffix are for use by Reverse and support prolate
      // spheroids by interchanging the roles of a and b.
    , _ax(_f >= 0 ? _a : _a * (1 - _f))
    , _e2x(_f >= 0 ? _e2 : - _e2/(1 - _e2))
    , _e4x(sq(_e2x))
    , _e2mx(_f >= 0 ? _e2m : 1/_e2m)
    , _maxrad(2 * _ax / numeric_limits<real>::epsilon())
  {}

  const Geocentric Geocentric::WGS84(Constants::WGS84_a(),
                                     Constants::WGS84_r());

  void Geocentric::Forward(real lat, real lon, real h,
                           real& x, real& y, real& z) const throw() {
    lon = lon >= 180 ? lon - 360 : lon < -180 ? lon + 360 : lon;
    real
      phi = lat * Constants::degree(),
      lam = lon * Constants::degree(),
      sphi = sin(phi),
      cphi = abs(lat) == 90 ? 0 : cos(phi),
      n = _a/sqrt(1 - _e2 * sq(sphi));
    z = ( _e2m * n + h) * sphi;
    x = (n + h) * cphi;
    y = x * (lon == -180 ? 0 : sin(lam));
    x *= (abs(lon) == 90 ? 0 : cos(lam));
  }

  void Geocentric::Reverse(real x, real y, real z,
                           real& lat, real& lon, real& h) const throw() {
    real R = Math::hypot(x, y);
    h = Math::hypot(R, z);      // Distance to center of earth
    real phi;
    if (h > _maxrad)
      // We really far away (> 12 million light years); treat the earth as a
      // point and h, above, is an acceptable approximation to the height.
      // This avoids overflow, e.g., in the computation of disc below.  It's
      // possible that h has overflowed to inf; but that's OK.
      //
      // Treat the case x, y finite, but R overflows to +inf by scaling by 2.
      phi = atan2(z/2, Math::hypot(x/2, y/2));
    else if (_e4x == 0) {
      // Treat the spherical case.  Dealing with underflow in the general case
      // with _e2 = 0 is difficult.  Origin maps to N pole same as an
      // ellipsoid.
      phi = atan2(h != 0 ? z : real(1), R);
      h -= _ax;
    } else {
      // Treat prolate spheroids by swapping R and z here and by switching
      // the arguments to phi = atan2(...) at the end.
      if (_f < 0) swap(R, z);
      real
        p = sq(R / _ax),
        q = _e2mx * sq(z / _ax),
        r = (p + q - _e4x) / 6;
      if ( !(_e4x * q == 0 && r <= 0) ) {
        real
          // Avoid possible division by zero when r = 0 by multiplying
          // equations for s and t by r^3 and r, resp.
          S = _e4x * p * q / 4, // S = r^3 * s
          r2 = sq(r),
          r3 = r * r2,
          disc =  S * (2 * r3 + S);
        real u = r;
        if (disc >= 0) {
          real T3 = r3 + S;
          // Pick the sign on the sqrt to maximize abs(T3).  This minimizes
          // loss of precision due to cancellation.  The result is unchanged
          // because of the way the T is used in definition of u.
          T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
          // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
          real T = Math::cbrt(T3); // T = r * t
          // T can be zero; but then r2 / T -> 0.
          u += T + (T != 0 ? r2 / T : 0);
        } else {
          // T is complex, but the way u is defined the result is real.
          real ang = atan2(sqrt(-disc), r3 + S);
          // There are three possible real solutions for u depending on the
          // multiple of 2*pi here.  We choose multiplier = 1 which leads to a
          // jump in the solution across the line 2 + s = 0; but this
          // nevertheless leads to a continuous (and accurate) solution for k.
          // Other choices of the multiplier lead to poorly conditioned
          // solutions near s = 0 (i.e., near p = 0 or q = 0).
          u += 2 * abs(r) * cos((2 * Constants::pi() + ang) / real(3));
        }
        real
          v = sqrt(sq(u) + _e4x * q), // guaranteed positive
          // Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
          // e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
          uv = u < 0 ? _e4x * q / (v - u) : u + v, // u+v, guaranteed positive
          // Need to guard against w going negative due to roundoff in uv - q.
          w = max(real(0), _e2x * (uv - q) / (2 * v)),
          // Rearrange expression for k to avoid loss of accuracy due to
          // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
          k = uv / (sqrt(uv + sq(w)) + w), // guaranteed positive
          d = k * R / (k + _e2x);
        // Probably atan2 returns the result for phi more accurately than the
        // half-angle formula that Vermeille uses.  It's certainly simpler.
        phi = _f >= 0 ? atan2(z, d) : atan2(d, z);
        h = (k + _e2x - 1) * Math::hypot(d, z) / k;
      } else {                  // e4 * q == 0 && r <= 0
        // Very near equatorial plane with R <= a * e^2.  This leads to k = 0
        // using the general formula and division by 0 in formula for h.  So
        // handle this case directly.  The condition e4 * q == 0 implies abs(z)
        // < 1.e-145 for WGS84 so it's OK to treat these points as though z =
        // 0.  (But we do take care that the sign of phi matches the sign of
        // z.)
        phi = _f >= 0 ?
          atan2(sqrt( -6 * r), sqrt(p * _e2mx)) :
          atan2(sqrt(p * _e2mx), sqrt( -6 * r));
        if (z < 0) phi = -phi;  // for tiny negative z (not for prolate)
        h = - _a * (_f >= 0 ? _e2m : real(1)) / sqrt(1 - _e2 * sq(sin(phi)));
      }
    }
    lat = phi / Constants::degree();
    // Negative signs return lon in [-180, 180).  Assume atan2(0,0) = 0.
    lon = -atan2(-y, x) / Constants::degree();
  }

} // namespace GeographicLib


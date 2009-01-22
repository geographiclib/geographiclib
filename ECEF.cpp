/**
 * \file ECEF.cpp
 * \brief Implementation for GeographicLib::ECEF class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 *
 * For the reverse transformation we use H. Vermeille,
 * <a href="http://dx.doi.org/10.1007/s00190-002-0273-6">
 * Direct transformation from geocentric coordinates to geodetic
 * coordinates</a>, J. Geodesy 76, 451&ndash;454 (2002).

 * This provides a closed-form solution but can't directly be applied close to
 * the center of the earth.  Several changes have been made to remove this
 * restriction and to improve the numerical accuracy.  Now the method is
 * accurate for all inputs (provided that \e h is finite).

 * The problems encountered near the center of the ellipsoid are:
 * - There's a potential division by zero in the definition of \e s. The
 *   equations are easily reformulated to avoid this problem.
 * - \e t<sup>3</sup> may be negative.  This is OK; we just take the real root.
 * - The solution for \e t may be complex.  However this leads to 3 real roots
 *   for \e u/\e r.  It's then just a matter of picking the one that computes
 *   the geodetic result which minimizes |\e h| and which avoids large
 *   round-off errors.
 * - Some of the equations result in a large loss of accuracy due to
 *   subtracting nearly equal quantities.  E.g., \e k= sqrt(\e u + \e v + \e
 *   w<sup>2</sup>) - \e w is inaccurate if \e u + \e v is small; we can fix
 *   this by writing \e k= (\e u + \e v)/(sqrt(\e u + \e v + \e w<sup>2</sup>)
 *   + \e w).
 *
 * The error is computed as follows.  Write a version of ECEF::WGS84.Forward
 * which uses long doubles (including using long doubles for the WGS84
 * parameters).  Generate random (long double) geodetic coordinates (\e lat0,
 * \e lon0, \e h0) and use the "long double" WGS84.Forward to obtain the
 * corresponding (long double) ECEF coordinates (\e x0, \e y0, \e z0).  [We
 * restrict \e h0 so that \e h0 >= - \e a (1 - \e e<sup>2</sup>) / sqrt(1 - \e
 * e<sup>2</sup> sin<sup>2</sup>\e lat0), which ensures that (\e lat0, \e lon0,
 * \e h0) is the principal geodetic inverse of (\e x0, \e y0, \e z0).]  Because
 * the forward calculation is numerically stable and because long doubles (on
 * Linux systems using G++) provide 11 bits additional accuracy (about 3.3
 * decimal digits), we regard this set of test data as exact.
 *
 * Apply the double version of WGS84.Reverse to (\e x0, \e y0, \e z0) to
 * compute the approximate geodetic coordinates (\e lat1, \e lon1, \e h1).
 * Convert (\e lat1 - \e lat0, \e lon1 - \e lon0) to a distance, \e ds, on the
 * surface of the ellipsoid and define \e err = hypot(\e ds, \e h1 - \e h0).
 * For |\e h0| < 5000 km, we have \e err < 7 nm.
 *
 * This methodology is not very useful very far from the globe, because the
 * absolute errors in the approximate geodetic height become large, or within
 * 50 km of the center of the earth, because of errors in computing the
 * approximate geodetic latitude.  To illustrate the second issue, consider
 * points with ECEF coordinates satisfying hypot(\e x, \e y) = \e a \e
 * e<sup>2</sup> (= 42.7 km), \e z = 0.  (This is the center of meridional
 * curvature for \e lat = 0.)  The geodetic latitude for these points is \e lat
 * = 0.  However, if we move 1 nm towards the center of the earth, the geodetic
 * latitude becomes 0.04", a distance of 1.4 m from the equator.  If, instead,
 * we move 1 nm up, the geodetic latitude becomes 7.45", a distance of 229 m
 * from the equator.
 *
 * To obtain a practical measure of the error for the general case we define
 * - <i>err</i><sub>h</sub> = |\e h1 - \e h0| / max(1, \e h0 / \e a)
 * - for \e h0 > 0, <i>err</i><sub>out</sub> = \e ds
 * - for \e h0 < 0, apply the long double version of WGS84.Forward to (\e lat1,
 *   \e lon1, \e h1) to give (\e x1, \e y1, \e z1) and compute
 *   <i>err</i><sub>in</sub> = hypot(\e x1 - \e x0, \e y1 - \e y0, \e z1 - \e
 *   z0).
 * .
 * We then find <i>err</i><sub>h</sub> < 8 nm, <i>err</i><sub>out</sub> < 4 nm,
 * and <i>err</i><sub>in</sub> < 7 nm.
 **********************************************************************/

#include "GeographicLib/ECEF.hpp"
#include "GeographicLib/Constants.hpp"
#include <algorithm>
#include <limits>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = ECEF_HPP;
}

namespace GeographicLib {

  using namespace std;

  ECEF::ECEF(double a, double invf)
    : _a(a)
    , _f(1 / invf)
    , _e2(_f * (2 - _f))
    , _e12(_e2 / (1 - _e2))
    , _b(_a * (1 - _f))
    , _tol(_a * numeric_limits<double>::epsilon())
    , _maxrad(2 * _a / numeric_limits<double>::epsilon())
  {}

  const ECEF ECEF::WGS84(Constants::WGS84_a, Constants::WGS84_invf);

  void ECEF::Forward(double lat, double lon, double h,
		     double& x, double& y, double& z) const {
    double
      phi = lat * Constants::degree,
      lam = lon * Constants::degree,
      sphi = sin(phi),
      n = _a/sqrt(1 - _e2 * sq(sphi));
    z = ( sq(1 - _f) * n + h) * sphi;
    x = (n + h) * cos(phi);
    y = x * sin(lam);
    x *= cos(lam);
  }

  void ECEF::Reverse(double x, double y, double z,
		     double& lat, double& lon, double& h) const {
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
	// Other choices of the multiplier lead to poorly conditioned solutions
	// near s = 0 (i.e., near p = 0 or q = 0).
	u += 2 * abs(r) * cos((2 * Constants::pi + ang) / 3.0);
      }
      double
	v = sqrt(sq(u) + e4 * q), // guaranteed positive
	// Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
	// e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
	uv = u < 0 ? e4 * q / (v - u) : u + v, //  u + v, guaranteed positive
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
      // handle this case directly.  The condition e4 * q == 0 implies abs(z) <
      // 1.e-145 for WGS84 so it's OK to treat these points as though z = 0.
      // (But we do take care that the sign of phi matches the sign of z.)
      phi = atan2(sqrt( -6 * r), sqrt(p * (1 - _e2)));
      if (z < 0) phi = -phi;	// for tiny negative z
      h = _a * (_e2 - 1) / sqrt(1 - _e2 * sq(sin(phi)));
    }
    lat = phi / Constants::degree;
    // Negative signs return lon in [-180, 180)
    lon = (rad != 0 ? -atan2(-y, x) : 0) / Constants::degree;
  }

} // namespace GeographicLib


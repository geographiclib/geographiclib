/**
 * \file Geodesic.cpp
 * \brief Implementation for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * This is a reformulation of the geodesic problem.  The notation is as
 * follows:
 * - at a general point (no suffix or 1 or 2 as suffix)
 *   - phi = latitude
 *   - beta = latitude on auxiliary sphere
 *   - omega = longitude on auxiliary sphere
 *   - lambda = longitude
 *   - alpha = azimuth of great circle
 *   - sigma = arc length along greate circle
 *   - s = distance
 *   - tau = scaled distance (= sigma at multiples of pi/2)
 * - at northwards equator crossing
 *   - beta = phi = 0
 *   - omega = lambda = 0
 *   - alpha = alpha0
 *   - sigma = s = 0
 * - a 12 suffix means a difference, e.g., s12 = s2 - s1.
 * - s and c prefixes mean sin and cos
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"

#define GEOGRAPHICLIB_GEODESIC_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GEODESIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_GEODESIC_HPP)

namespace GeographicLib {

  using namespace std;

  // Underflow guard.  We require
  //   eps2 * epsilon() > 0
  //   eps2 + epsilon() == epsilon()
  const Math::real Geodesic::eps2 = sqrt(numeric_limits<real>::min());
  const Math::real Geodesic::tol0 = numeric_limits<real>::epsilon();
  const Math::real Geodesic::tol1 = 100 * tol0;
  const Math::real Geodesic::tol2 = sqrt(numeric_limits<real>::epsilon());
  const Math::real Geodesic::xthresh = 1000 * tol2;

  Geodesic::Geodesic(real a, real r)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _f1(1 - _f)
    , _e2(_f * (2 - _f))
    , _ep2(_e2 / sq(_f1))       // e2 / (1 - e2)
    , _n(_f / ( 2 - _f))
    , _b(_a * _f1)
    , _etol2(tol2 / max(real(0.1), sqrt(abs(_e2))))
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    A3coeff();
    C3coeff();
  }

  const Geodesic Geodesic::WGS84(Constants::WGS84_a(), Constants::WGS84_r());

  Math::real Geodesic::SinSeries(real sinx, real cosx,
                                 const real c[], int n) throw() {
    // Evaluate y = sum(c[i] * sin(2 * i * x), i, 1, n) using Clenshaw
    // summation.  (c[0] is unused.)
    // Approx operation count = (n + 5) mult and (2 * n + 2) add
    real
      ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
      y0 = n & 1 ? c[n--] : 0, y1 = 0;        // Accumulators for sum
    // Now n is even
    while (n) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + c[n--];
      y0 = ar * y1 - y0 + c[n--];
    }
    return 2 * sinx * cosx * y0; // sin(2 * x) * y0
  }

  GeodesicLine Geodesic::Line(real lat1, real lon1, real azi1) const throw() {
    return GeodesicLine(*this, lat1, lon1, azi1);
  }

  Math::real Geodesic::Direct(real lat1, real lon1, real azi1, real s12,
                              real& lat2, real& lon2, real& azi2, real& m12,
                              bool arcmode) const throw() {
    GeodesicLine l(*this, lat1, lon1, azi1);
    return l.Position(s12, lat2, lon2, azi2, m12, arcmode);
  }

  Math::real Geodesic::Inverse(real lat1, real lon1, real lat2, real lon2,
                               real& s12, real& azi1, real& azi2, real& m12)
    const throw() {
    lon1 = AngNormalize(lon1);
    real lon12 = AngNormalize(AngNormalize(lon2) - lon1);
    // If very close to being on the same meridian, then make it so.
    // Not sure this is necessary...
    lon12 = AngRound(lon12);
    // Make longitude difference positive.
    int lonsign = lon12 >= 0 ? 1 : -1;
    lon12 *= lonsign;
    if (lon12 == 180)
      lonsign = 1;
    // If really close to the equator, treat as on equator.
    lat1 = AngRound(lat1);
    lat2 = AngRound(lat2);
    // Swap points so that point with higher (abs) latitude is point 1
    int swapp = abs(lat1) >= abs(lat2) ? 1 : -1;
    if (swapp < 0) {
      lonsign *= -1;
      swap(lat1, lat2);
    }
    // Make lat1 <= 0
    int latsign = lat1 < 0 ? 1 : -1;
    lat1 *= latsign;
    lat2 *= latsign;
    // Now we have
    //
    //     0 <= lon12 <= 180
    //     -90 <= lat1 <= 0
    //     lat1 <= lat2 <= -lat1
    //
    // longsign, swapp, latsign register the transformation to bring the
    // coordinates to this canonical form.  In all cases, 1 means no change was
    // made.  We make these transformations so that there are few cases to
    // check, e.g., on verifying quadrants in atan2.  In addition, this
    // enforces some symmetries in the results returned.

    real phi, sbet1, cbet1, sbet2, cbet2;

    phi = lat1 * Constants::degree();
    // Ensure cbet1 = +epsilon at poles
    sbet1 = _f1 * sin(phi);
    cbet1 = lat1 == -90 ? eps2 : cos(phi);
    SinCosNorm(sbet1, cbet1);

    phi = lat2 * Constants::degree();
    // Ensure cbet2 = +epsilon at poles
    sbet2 = _f1 * sin(phi);
    cbet2 = abs(lat2) == 90 ? eps2 : cos(phi);
    SinCosNorm(sbet2, cbet2);

    real
      lam12 = lon12 * Constants::degree(),
      slam12 = lon12 == 180 ? 0 : sin(lam12),
      clam12 = cos(lam12);      // lon12 == 90 isn't interesting

    real sig12, calp1, salp1, calp2, salp2;
    // index zero elements of these arrays are unused
    real C1a[nC1 + 1], C2a[nC2 + 1], C3a[nC3];

    bool meridian = lat1 == -90 || slam12 == 0;

    if (meridian) {

      // Endpoints are on a single full meridian, so the geodesic might lie on
      // a meridian.

      calp1 = clam12; salp1 = slam12; // Head to the target longitude
      calp2 = 1; salp2 = 0;           // At the target we're heading north

      real
        // tan(bet) = tan(sig) * cos(alp),
        ssig1 = sbet1, csig1 = calp1 * cbet1,
        ssig2 = sbet2, csig2 = calp2 * cbet2;

      // sig12 = sig2 - sig1
      sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, real(0)),
                    csig1 * csig2 + ssig1 * ssig2);
      {
        real dummy;
        Lengths(_n, sig12, ssig1, csig1, ssig2, csig2,
                cbet1, cbet2, s12, m12, dummy, C1a, C2a);
      }
      // Add the check for sig12 since zero length geodesics might yield m12 <
      // 0.  Test case was
      //
      //    echo 20.001 0 20.001 0 | Geod -i
      //
      // In fact, we will have sig12 > pi/2 for meridional geodesic which is
      // not a shortest path.
      if (sig12 < 1 || m12 >= 0) {
        m12 *= _a;
        s12 *= _b;
        sig12 /= Constants::degree();
      } else
        // m12 < 0, i.e., prolate and too close to anti-podal
        meridian = false;
    }

    if (!meridian &&
        sbet1 == 0 &&   // and sbet2 == 0
         // Mimic the way Lambda12 works with calp1 = 0
        (_f <= 0 || lam12 <= Constants::pi() - _f * Constants::pi())) {

      // Geodesic runs along equator
      calp1 = calp2 = 0; salp1 = salp2 = 1;
      s12 = _a * lam12;
      m12 = _b * sin(lam12 / _f1);
      sig12 = lon12 / _f1;

    } else if (!meridian) {

      // Now point1 and point2 belong within a hemisphere bounded by a
      // meridian and geodesic is neither meridional or equatorial.

      // Figure a starting point for Newton's method
      sig12 = InverseStart(sbet1, cbet1, sbet2, cbet2,
                           lam12,
                           salp1, calp1, salp2, calp2,
                           C1a, C2a);

      if (sig12 >= 0) {
        // Short lines (InverseStart sets salp2, calp2)
        s12 = m12 = sig12 * _a * sqrt(1 - _e2 * sq(cbet1));
        sig12 /= Constants::degree();
      } else {

        // Newton's method
        real ssig1, csig1, ssig2, csig2, eps;
        real ov = 0;
        unsigned numit = 0;
        for (unsigned trip = 0; numit < maxit; ++numit) {
          real dv;
          real v = Lambda12(sbet1, cbet1, sbet2, cbet2, salp1, calp1,
                            salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                            eps, trip < 1, dv, C1a, C2a, C3a) - lam12;

          if (abs(v) <= eps2 || !(trip < 1)) {
            if (abs(v) > max(tol1, ov))
              numit = maxit;
            break;
          }
          real
            dalp1 = -v/dv;
          real
            sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
            nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
          calp1 = calp1 * cdalp1 - salp1 * sdalp1;
          salp1 = max(real(0), nsalp1);
          SinCosNorm(salp1, calp1);
          // In some regimes we don't get quadratic convergence because slope
          // -> 0.  So use convergence conditions based on epsilon instead of
          // sqrt(epsilon).  The first criterion is a test on abs(v) against
          // 100 * epsilon.  The second takes credit for an anticipated
          // reduction in abs(v) by v/ov (due to the latest update in alp1) and
          // checks this against epsilon.
          if (abs(v) < tol1 || sq(v) < ov * tol0) ++trip;
          ov = abs(v);
        }

        {
          real dummy;
          Lengths(eps, sig12, ssig1, csig1, ssig2, csig2,
                  cbet1, cbet2, s12, m12, dummy, C1a, C2a);
        }
        m12 *= _a;
        s12 *= _b;
        sig12 /= Constants::degree();

        if (numit >= maxit) {
          // Signal failure to converge by negating the distance and azimuths.
          s12 *= -1; sig12 *= -1; m12 *= -1;
          salp1 *= -1; calp1 *= -1;
          salp2 *= -1; calp2 *= -1;
        }
      }
    }

    // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
    if (swapp < 0) {
      swap(salp1, salp2);
      swap(calp1, calp2);
    }

    // minus signs give range [-180, 180). 0- converts -0 to +0.
    azi1 = 0 - atan2(- swapp * lonsign * salp1,
                     + swapp * latsign * calp1) / Constants::degree();
    azi2 = 0 - atan2(- swapp * lonsign * salp2,
                     + swapp * latsign * calp2) / Constants::degree();
    // Returned value in [0, 180], unless it's negated to signal convergence
    // failure
    return sig12;
  }

  void Geodesic::Lengths(real eps, real sig12,
                         real ssig1, real csig1, real ssig2, real csig2,
                         real cbet1, real cbet2,
                         real& s12b, real& m12a, real& m0,
                         // Scratch areas of the right size
                         real C1a[], real C2a[]) const throw() {
    // Return m12a = (reduced length)/_a; also calculate s12b = distance/_b,
    // and m0 = coefficient of secular term in expression for reduced length.
    C1f(eps, C1a);
    C2f(eps, C2a);
    real
      A1m1 = A1m1f(eps),
      AB1 = (1 + A1m1) * (SinSeries(ssig2, csig2, C1a, nC1) -
                          SinSeries(ssig1, csig1, C1a, nC1)),
      A2m1 = A2m1f(eps),
      AB2 = (1 + A2m1) * (SinSeries(ssig2, csig2, C2a, nC2) -
                          SinSeries(ssig1, csig1, C2a, nC2));
    m0 = A1m1 - A2m1;
    // Missing a factor of _a.
    // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure accurate
    // cancellation in the case of coincident points.
    m12a = (sqrt(1 - _e2 * sq(cbet2)) * (csig1 * ssig2) -
            sqrt(1 - _e2 * sq(cbet1)) * (ssig1 * csig2))
      - _f1 * csig1 * csig2 * ( m0 * sig12 + (AB1 - AB2) );
    // Missing a factor of _b
    s12b =  (1 + A1m1) * sig12 + AB1;
  }

  Math::real Geodesic::Astroid(real x, real y) throw() {
    // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
    // This solution is adapted from Geocentric::Reverse.
    real k;
    real
      p = sq(x),
      q = sq(y),
      r = (p + q - 1) / 6;
    if ( !(q == 0 && r <= 0) ) {
      real
        // Avoid possible division by zero when r = 0 by multiplying equations
        // for s and t by r^3 and r, resp.
        S = p * q / 4,            // S = r^3 * s
        r2 = sq(r),
        r3 = r * r2,
        // The discrimant of the quadratic equation for T3.  This is zero on
        // the evolute curve p^(1/3)+q^(1/3) = 1
        disc =  S * (S + 2 * r3);
      real u = r;
      if (disc >= 0) {
        real T3 = S + r3;
        // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
        // of precision due to cancellation.  The result is unchanged because
        // of the way the T is used in definition of u.
        T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
        // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
        real T = Math::cbrt(T3); // T = r * t
        // T can be zero; but then r2 / T -> 0.
        u += T + (T != 0 ? r2 / T : 0);
      } else {
        // T is complex, but the way u is defined the result is real.
        real ang = atan2(sqrt(-disc), -(S + r3));
        // There are three possible cube roots.  We choose the root which
        // avoids cancellation.  Note that disc < 0 implies that r < 0.
        u += 2 * r * cos(ang / 3);
      }
      real
        v = sqrt(sq(u) + q),    // guaranteed positive
        // Avoid loss of accuracy when u < 0.
        uv = u < 0 ? q / (v - u) : u + v, // u+v, guaranteed positive
        w = (uv - q) / (2 * v);           // positive?
      // Rearrange expression for k to avoid loss of accuracy due to
      // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
      k = uv / (sqrt(uv + sq(w)) + w);   // guaranteed positive
    } else {               // q == 0 && r <= 0
      // y = 0 with |x| <= 1.  Handle this case directly.
      // for y small, positive root is k = abs(y)/sqrt(1-x^2)
      k = 0;
    }
    return k;
  }

  Math::real Geodesic::InverseStart(real sbet1, real cbet1,
                                    real sbet2, real cbet2,
                                    real lam12,
                                    real& salp1, real& calp1,
                                    // Only updated if return val >= 0
                                    real& salp2, real& calp2,
                                    // Scratch areas of the right size
                                    real C1a[], real C2a[]) const throw() {
    // Return a starting point for Newton's method in salp1 and calp1 (function
    // value is -1).  If Newton's method doesn't need to be used, return also
    // salp2 and calp2 and function value is sig12.
    real
      sig12 = -1,               // Return value
      // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
      sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
      cbet12 = cbet2 * cbet1 + sbet2 * sbet1,
      sbet12a = sbet2 * cbet1 + cbet2 * sbet1;

    bool shortline = cbet12 >= 0 && sbet12 < real(0.5) &&
      lam12 <= Constants::pi() / 6;
    real
      omg12 = shortline ? lam12 / sqrt(1 - _e2 * sq(cbet1)) : lam12,
      somg12 = sin(omg12), comg12 = cos(omg12);

    salp1 = cbet2 * somg12;
    calp1 = comg12 >= 0 ?
      sbet12 + cbet2 * sbet1 * sq(somg12) / (1 + comg12) :
      sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);

    real
      ssig12 = Math::hypot(salp1, calp1),
      csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

    if (shortline && ssig12 < _etol2) {
      // really short lines
      salp2 = cbet1 * somg12;
      calp2 = sbet12 - cbet1 * sbet2 * sq(somg12) / (1 + comg12);
      // Set return value
      sig12 = atan2(ssig12, csig12);
    } else if (csig12 >= 0 ||
               ssig12 >= 3 * abs(_f) * Constants::pi() * sq(cbet1)) {
      // Nothing to do, zeroth order spherical approximation is OK

    } else {
      // Scale lam12 and bet2 to x, y coordinate system where antipodal point
      // is at origin and singular point is at y = 0, x = -1.
      real x, y, lamscale, betscale;
      if (_f >= 0) {            // In fact f == 0 does not get here
        // x = dlong, y = dlat
        {
          real
            mu = sq(sbet1),
            k2 = mu * _ep2,
            eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
          lamscale = _f * cbet1 * A3f(eps) * Constants::pi();
        }
        betscale = lamscale * cbet1;

        x = (lam12 - Constants::pi()) / lamscale;
        y = sbet12a / betscale;
      } else {                  // _f < 0
        // x = dlat, y = dlong
        real
          cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
          bet12a = atan2(sbet12a, cbet12a);
        real m0, dummy;
        // In the case of lon12 = 180, this repeats a calculation made in
        // Inverse.
        Lengths(_n, Constants::pi() + bet12a, sbet1, -cbet1, sbet2, cbet2,
                cbet1, cbet2, dummy, x, m0, C1a, C2a);
        x = -1 + x/(_f1 * cbet1 * cbet2 * m0 * Constants::pi());
        betscale = x < -real(0.01) ? sbet12a / x :
          -_f * sq(cbet1) * Constants::pi();
        lamscale = betscale / cbet1;
        y = (lam12 - Constants::pi()) / lamscale;
      }

      if (y > -tol1 && x >  -1 - xthresh) {
        // strip near cut
        if (_f >= 0) {
          salp1 = min(real( 1), -x); calp1 = - sqrt(1 - sq(salp1));
        } else {
          calp1 = max(real(-1),  x); salp1 =   sqrt(1 - sq(calp1));
        }
      } else {
        // Estimate omega12, by solving the astroid problem.
        real k = Astroid(x, y);
        // estimate omg12a = pi - omg12
        real
          omg12a = lamscale * ( _f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k ),
          somg12 = sin(omg12a), comg12 = -cos(omg12a);
        // Update spherical estimate of alp1 using omg12 instead of lam12
        salp1 = cbet2 * somg12;
        calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
      }
    }
    SinCosNorm(salp1, calp1);
    return sig12;
  }

  Math::real Geodesic::Lambda12(real sbet1, real cbet1, real sbet2, real cbet2,
                                real salp1, real calp1,
                                real& salp2, real& calp2,
                                real& sig12,
                                real& ssig1, real& csig1,
                                real& ssig2, real& csig2,
                                real& eps, bool diffp, real& dlam12,
                                // Scratch areas of the right size
                                real C1a[], real C2a[], real C3a[]) const
    throw() {

    if (sbet1 == 0 && calp1 == 0)
      // Break degeneracy of equatorial line.  This case has already been
      // handled.
      calp1 = -eps2;

    real
      // sin(alp1) * cos(bet1) = sin(alp0),
      salp0 = salp1 * cbet1,
      calp0 = Math::hypot(calp1, salp1 * sbet1); // calp0 > 0

    real somg1, comg1, somg2, comg2, omg12, lam12, mu, k2;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
    ssig1 = sbet1; somg1 = salp0 * sbet1;
    csig1 = comg1 = calp1 * cbet1;
    SinCosNorm(ssig1, csig1);
    SinCosNorm(somg1, comg1);

    // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
    // about this case, since this can yield singularities in the Newton
    // iteration.
    // sin(alp2) * cos(bet2) = sin(alp0),
    salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
    // calp2 = sqrt(1 - sq(salp2))
    //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
    // and subst for calp0 and rearrange to give (choose positive sqrt
    // to give alp2 in [0, pi/2]).
    calp2 = cbet2 != cbet1 || abs(sbet2) != -sbet1 ?
      sqrt(sq(calp1 * cbet1) + (cbet1 < -sbet1 ?
                                (cbet2 - cbet1) * (cbet1 + cbet2) :
                                (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
      abs(calp1);
    // tan(bet2) = tan(sig2) * cos(alp2)
    // tan(omg2) = sin(alp0) * tan(sig2).
    ssig2 = sbet2; somg2 = salp0 * sbet2;
    csig2 = comg2 = calp2 * cbet2;
    SinCosNorm(ssig2, csig2);
    SinCosNorm(somg2, comg2);

    // sig12 = sig2 - sig1, limit to [0, pi]
    sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, real(0)),
                  csig1 * csig2 + ssig1 * ssig2);

    // omg12 = omg2 - omg1, limit to [0, pi]
    omg12 = atan2(max(comg1 * somg2 - somg1 * comg2, real(0)),
                  comg1 * comg2 + somg1 * somg2);
    real B312, h0;
    mu = sq(calp0);
    k2 = mu * _ep2;
    eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
    C3f(eps, C3a);
    B312 = (SinSeries(ssig2, csig2, C3a, nC3-1) -
            SinSeries(ssig1, csig1, C3a, nC3-1));
    h0 = -_f * A3f(eps),
    lam12 = omg12 + salp0 * h0 * (sig12 + B312);

    if (diffp) {
      if (calp2 == 0)
        dlam12 = - 2 * sqrt(1 - _e2 * sq(cbet1)) / sbet1;
      else {
        real dummy1, dummy2;
        Lengths(eps, sig12, ssig1, csig1, ssig2, csig2,
                cbet1, cbet2, dummy1, dlam12, dummy2, C1a, C2a);
        dlam12 /= calp2 * cbet2;
      }
    }

    return lam12;
  }

  GeodesicLine::GeodesicLine(const Geodesic& g,
                             real lat1, real lon1, real azi1) throw()
    : _a(g._a)
    , _r(g._r)
    , _b(g._b)
    , _f1(g._f1)
  {
    azi1 = Geodesic::AngNormalize(azi1);
    // Guard against underflow in salp0
    azi1 = Geodesic::AngRound(azi1);
    lon1 = Geodesic::AngNormalize(lon1);
    _lat1 = lat1;
    _lon1 = lon1;
    _azi1 = azi1;
    // alp1 is in [0, pi]
    real
      alp1 = azi1 * Constants::degree(),
      // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
      // problems directly than to skirt them.
      salp1 = azi1 == -180 ? 0 : sin(alp1),
      calp1 = azi1 ==   90 ? 0 : cos(alp1);
    real cbet1, sbet1, phi;
    phi = lat1 * Constants::degree();
    // Ensure cbet1 = +epsilon at poles
    sbet1 = _f1 * sin(phi);
    cbet1 = abs(lat1) == 90 ? Geodesic::eps2 : cos(phi);
    Geodesic::SinCosNorm(sbet1, cbet1);

    // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    _salp0 = salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    _calp0 = Math::hypot(calp1, salp1 * sbet1);
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    _ssig1 = sbet1; _somg1 = _salp0 * sbet1;
    _csig1 = _comg1 = sbet1 != 0 || calp1 != 0 ? cbet1 * calp1 : 1;
    Geodesic::SinCosNorm(_ssig1, _csig1); // sig1 in (-pi, pi]
    Geodesic::SinCosNorm(_somg1, _comg1);

    real mu = sq(_calp0);
    _k2 = mu * g._ep2;
    real eps = _k2 / (2 * (1 + sqrt(1 + _k2)) + _k2);
    _A1m1 =  Geodesic::A1m1f(eps);
    _A2m1 =  Geodesic::A2m1f(eps);

    Geodesic::C1f(eps, _C1a);
    Geodesic::C1pf(eps, _C1pa);
    Geodesic::C2f(eps, _C2a);
    g.C3f(eps, _C3a);

    _B11 = Geodesic::SinSeries(_ssig1, _csig1, _C1a, nC1);
    _B21 = Geodesic::SinSeries(_ssig1, _csig1, _C2a, nC2);
    {
      real s = sin(_B11), c = cos(_B11);
      // tau1 = sig1 + B11
      _stau1 = _ssig1 * c + _csig1 * s;
      _ctau1 = _csig1 * c - _ssig1 * s;
    }
    // Not necessary because C1pa reverts C1a
    //    _B11 = -SinSeries(_stau1, _ctau1, _C1pa, nC1p);

    _A3c = -g._f * _salp0 * g.A3f(eps);
    _B31 = Geodesic::SinSeries(_ssig1, _csig1, _C3a, nC3-1);
  }

  Math::real GeodesicLine::Position(real s12, real& lat2, real& lon2,
                                    real& azi2, real& m12, bool arcmode)
  const throw() {
    if (!Init())
      // Uninitialized
      return 0;
    // Avoid warning about uninitialized B12.
    real sig12, ssig12, csig12, B12 = 0;
    if (arcmode) {
      // Interpret s12 as spherical arc length
      sig12 = s12 * Constants::degree();
      real s12a = abs(s12);
      s12a -= 180 * floor(s12a / 180);
      ssig12 = s12a ==  0 ? 0 : sin(sig12);
      csig12 = s12a == 90 ? 0 : cos(sig12);
    } else {
      // Interpret s12 as distance
      real
        tau12 = s12 / (_b * (1 + _A1m1)),
        s = sin(tau12),
        c = cos(tau12);
      // tau2 = tau1 + tau12
      B12 = - Geodesic::SinSeries(_stau1 * c + _ctau1 * s,
                                  _ctau1 * c - _stau1 * s,
                                  _C1pa, nC1p);
      sig12 = tau12 - (B12 - _B11);
      ssig12 = sin(sig12);
      csig12 = cos(sig12);
    }

    real omg12, lam12, lon12;
    real ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
    csig2 = _csig1 * csig12 - _ssig1 * ssig12;
    if (arcmode)
      B12 = Geodesic::SinSeries(ssig2, csig2, _C1a, nC1);
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = _calp0 * ssig2;
    // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
    cbet2 = Math::hypot(_salp0, _calp0 * csig2);
    if (cbet2 == 0)
      // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
      cbet2 = csig2 = Geodesic::eps2;
    // tan(omg2) = sin(alp0) * tan(sig2)
    somg2 = _salp0 * ssig2; comg2 = csig2;  // No need to normalize
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = _salp0; calp2 = _calp0 * csig2; // No need to normalize
    // omg12 = omg2 - omg1
    omg12 = atan2(somg2 * _comg1 - comg2 * _somg1,
                  comg2 * _comg1 + somg2 * _somg1);
    lam12 = omg12 + _A3c *
      ( sig12 + (Geodesic::SinSeries(ssig2, csig2, _C3a, nC3-1)  - _B31));
    lon12 = lam12 / Constants::degree();
    // Can't use AngNormalize because longitude might have wrapped multiple
    // times.
    lon12 = lon12 - 360 * floor(lon12/360 + real(0.5));
    lat2 = atan2(sbet2, _f1 * cbet2) / Constants::degree();
    lon2 = Geodesic::AngNormalize(_lon1 + lon12);
    // minus signs give range [-180, 180). 0- converts -0 to +0.
    azi2 = 0 - atan2(-salp2, calp2) / Constants::degree();

    real
      B22 = Geodesic::SinSeries(ssig2, csig2, _C2a, nC2),
      AB1 = (1 + _A1m1) * (B12 - _B11),
      AB2 = (1 + _A2m1) * (B22 - _B21);
    // Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
    // accurate cancellation in the case of coincident points.
    m12 = _b * ((sqrt(1 + _k2 * sq( ssig2)) * (_csig1 * ssig2) -
                 sqrt(1 + _k2 * sq(_ssig1)) * (_ssig1 * csig2))
                - _csig1 * csig2 * ( (_A1m1 - _A2m1) * sig12 + (AB1 - AB2) ));
    if (arcmode)
      s12 = _b * ((1 + _A1m1) * sig12 + AB1);

    return arcmode ? s12 : sig12 /  Constants::degree();
  }

  void GeodesicLine::Scale(real a12, real& M12, real& M21) const throw() {
    if (!Init())
      // Uninitialized
      return;
    real sig12 = a12 * Constants::degree(), ssig12, csig12;
    {
      real a12a = abs(a12);
      a12a -= 180 * floor(a12a / 180);
      ssig12 = a12a ==  0 ? 0 : sin(sig12);
      csig12 = a12a == 90 ? 0 : cos(sig12);
    }
    // sig2 = sig1 + sig12
    real
      ssig2 = _ssig1 * csig12 + _csig1 * ssig12,
      csig2 = _csig1 * csig12 - _ssig1 * ssig12,
      ssig1sq = sq(_ssig1),
      ssig2sq = sq( ssig2),
      w1 = sqrt(1 + _k2 * ssig1sq),
      w2 = sqrt(1 + _k2 * ssig2sq),
      B12 = Geodesic::SinSeries(ssig2, csig2, _C1a, nC1),
      B22 = Geodesic::SinSeries(ssig2, csig2, _C2a, nC2),
      AB1 = (1 + _A1m1) * (B12 - _B11),
      AB2 = (1 + _A2m1) * (B22 - _B21),
      J12 = (_A1m1 - _A2m1) * sig12 + (AB1 - AB2);
    M12 = csig12 + (_k2 * (ssig2sq - ssig1sq) * ssig2/ (w1 + w2)
                    - csig2 * J12) * _ssig1 / w1;
    M21 = csig12 - (_k2 * (ssig2sq - ssig1sq) * _ssig1/ (w1 + w2)
                    - _csig1 * J12) * ssig2 / w2;
  }

  Math::real Geodesic::A3f(real eps) const throw() {
    // Evaluation sum(_A3c[k] * eps^k, k, 0, nA3x-1) by Horner's method
    real v = 0;
    for (int i = nA3x; i; )
      v = eps * v + _A3x[--i];
    return v;
  }

  void Geodesic::C3f(real eps, real c[]) const throw() {
    // Evaluation C3 coeffs by Horner's method
    for (int j = nC3x, k = nC3; k; ) {
      real t = 0;
      for (int i = nC3 - k; i; --i)
        t = eps * t + _C3x[--j];
      c[k--] = t;
    }

    real mult = 1;
    for (int k = 0; k < nC3; ) {
      mult *= eps;
      c[++k] *= mult;
    }
  }

  // Generated by Maxima on 2010-05-03 08:35:50-04:00

  // The scale factor A1-1 = mean value of I1-1
  Math::real Geodesic::A1m1f(real eps) throw() {
    real
      eps2 = sq(eps),
      t;
    switch (nA1/2) {
    case 0:
      t = 0;
      break;
    case 1:
      t = eps2/4;
      break;
    case 2:
      t = eps2*(eps2+16)/64;
      break;
    case 3:
      t = eps2*(eps2*(eps2+4)+64)/256;
      break;
    case 4:
      t = eps2*(eps2*(eps2*(25*eps2+64)+256)+4096)/16384;
      break;
    default:
      STATIC_ASSERT(nA1 >= 0 && nA1 <= 8, "Bad value of nA1");
      t = 0;
    }
    return (t + eps) / (1 - eps);
  }

  // The coefficients C1[l] in the Fourier expansion of B1
  void Geodesic::C1f(real eps, real c[]) throw() {
    real
      eps2 = sq(eps),
      d = eps;
    switch (nC1) {
    case 0:
      break;
    case 1:
      c[1] = -d/2;
      break;
    case 2:
      c[1] = -d/2;
      d *= eps;
      c[2] = -d/16;
      break;
    case 3:
      c[1] = d*(3*eps2-8)/16;
      d *= eps;
      c[2] = -d/16;
      d *= eps;
      c[3] = -d/48;
      break;
    case 4:
      c[1] = d*(3*eps2-8)/16;
      d *= eps;
      c[2] = d*(eps2-2)/32;
      d *= eps;
      c[3] = -d/48;
      d *= eps;
      c[4] = -5*d/512;
      break;
    case 5:
      c[1] = d*((6-eps2)*eps2-16)/32;
      d *= eps;
      c[2] = d*(eps2-2)/32;
      d *= eps;
      c[3] = d*(9*eps2-16)/768;
      d *= eps;
      c[4] = -5*d/512;
      d *= eps;
      c[5] = -7*d/1280;
      break;
    case 6:
      c[1] = d*((6-eps2)*eps2-16)/32;
      d *= eps;
      c[2] = d*((64-9*eps2)*eps2-128)/2048;
      d *= eps;
      c[3] = d*(9*eps2-16)/768;
      d *= eps;
      c[4] = d*(3*eps2-5)/512;
      d *= eps;
      c[5] = -7*d/1280;
      d *= eps;
      c[6] = -7*d/2048;
      break;
    case 7:
      c[1] = d*(eps2*(eps2*(19*eps2-64)+384)-1024)/2048;
      d *= eps;
      c[2] = d*((64-9*eps2)*eps2-128)/2048;
      d *= eps;
      c[3] = d*((72-9*eps2)*eps2-128)/6144;
      d *= eps;
      c[4] = d*(3*eps2-5)/512;
      d *= eps;
      c[5] = d*(35*eps2-56)/10240;
      d *= eps;
      c[6] = -7*d/2048;
      d *= eps;
      c[7] = -33*d/14336;
      break;
    case 8:
      c[1] = d*(eps2*(eps2*(19*eps2-64)+384)-1024)/2048;
      d *= eps;
      c[2] = d*(eps2*(eps2*(7*eps2-18)+128)-256)/4096;
      d *= eps;
      c[3] = d*((72-9*eps2)*eps2-128)/6144;
      d *= eps;
      c[4] = d*((96-11*eps2)*eps2-160)/16384;
      d *= eps;
      c[5] = d*(35*eps2-56)/10240;
      d *= eps;
      c[6] = d*(9*eps2-14)/4096;
      d *= eps;
      c[7] = -33*d/14336;
      d *= eps;
      c[8] = -429*d/262144;
      break;
    default:
      STATIC_ASSERT(nC1 >= 0 && nC1 <= 8, "Bad value of nC1");
    }
  }

  // The coefficients C1p[l] in the Fourier expansion of B1p
  void Geodesic::C1pf(real eps, real c[]) throw() {
    real
      eps2 = sq(eps),
      d = eps;
    switch (nC1p) {
    case 0:
      break;
    case 1:
      c[1] = d/2;
      break;
    case 2:
      c[1] = d/2;
      d *= eps;
      c[2] = 5*d/16;
      break;
    case 3:
      c[1] = d*(16-9*eps2)/32;
      d *= eps;
      c[2] = 5*d/16;
      d *= eps;
      c[3] = 29*d/96;
      break;
    case 4:
      c[1] = d*(16-9*eps2)/32;
      d *= eps;
      c[2] = d*(30-37*eps2)/96;
      d *= eps;
      c[3] = 29*d/96;
      d *= eps;
      c[4] = 539*d/1536;
      break;
    case 5:
      c[1] = d*(eps2*(205*eps2-432)+768)/1536;
      d *= eps;
      c[2] = d*(30-37*eps2)/96;
      d *= eps;
      c[3] = d*(116-225*eps2)/384;
      d *= eps;
      c[4] = 539*d/1536;
      d *= eps;
      c[5] = 3467*d/7680;
      break;
    case 6:
      c[1] = d*(eps2*(205*eps2-432)+768)/1536;
      d *= eps;
      c[2] = d*(eps2*(4005*eps2-4736)+3840)/12288;
      d *= eps;
      c[3] = d*(116-225*eps2)/384;
      d *= eps;
      c[4] = d*(2695-7173*eps2)/7680;
      d *= eps;
      c[5] = 3467*d/7680;
      d *= eps;
      c[6] = 38081*d/61440;
      break;
    case 7:
      c[1] = d*(eps2*((9840-4879*eps2)*eps2-20736)+36864)/73728;
      d *= eps;
      c[2] = d*(eps2*(4005*eps2-4736)+3840)/12288;
      d *= eps;
      c[3] = d*(eps2*(8703*eps2-7200)+3712)/12288;
      d *= eps;
      c[4] = d*(2695-7173*eps2)/7680;
      d *= eps;
      c[5] = d*(41604-141115*eps2)/92160;
      d *= eps;
      c[6] = 38081*d/61440;
      d *= eps;
      c[7] = 459485*d/516096;
      break;
    case 8:
      c[1] = d*(eps2*((9840-4879*eps2)*eps2-20736)+36864)/73728;
      d *= eps;
      c[2] = d*(eps2*((120150-86171*eps2)*eps2-142080)+115200)/368640;
      d *= eps;
      c[3] = d*(eps2*(8703*eps2-7200)+3712)/12288;
      d *= eps;
      c[4] = d*(eps2*(1082857*eps2-688608)+258720)/737280;
      d *= eps;
      c[5] = d*(41604-141115*eps2)/92160;
      d *= eps;
      c[6] = d*(533134-2200311*eps2)/860160;
      d *= eps;
      c[7] = 459485*d/516096;
      d *= eps;
      c[8] = 109167851*d/82575360;
      break;
    default:
      STATIC_ASSERT(nC1p >= 0 && nC1p <= 8, "Bad value of nC1p");
    }
  }

  // The scale factor A2-1 = mean value of I2-1
  Math::real Geodesic::A2m1f(real eps) throw() {
    real
      eps2 = sq(eps),
      t;
    switch (nA2/2) {
    case 0:
      t = 0;
      break;
    case 1:
      t = eps2/4;
      break;
    case 2:
      t = eps2*(9*eps2+16)/64;
      break;
    case 3:
      t = eps2*(eps2*(25*eps2+36)+64)/256;
      break;
    case 4:
      t = eps2*(eps2*(eps2*(1225*eps2+1600)+2304)+4096)/16384;
      break;
    default:
      STATIC_ASSERT(nA2 >= 0 && nA2 <= 8, "Bad value of nA2");
      t = 0;
    }
    return t * (1 - eps) - eps;
  }

  // The coefficients C2[l] in the Fourier expansion of B2
  void Geodesic::C2f(real eps, real c[]) throw() {
    real
      eps2 = sq(eps),
      d = eps;
    switch (nC2) {
    case 0:
      break;
    case 1:
      c[1] = d/2;
      break;
    case 2:
      c[1] = d/2;
      d *= eps;
      c[2] = 3*d/16;
      break;
    case 3:
      c[1] = d*(eps2+8)/16;
      d *= eps;
      c[2] = 3*d/16;
      d *= eps;
      c[3] = 5*d/48;
      break;
    case 4:
      c[1] = d*(eps2+8)/16;
      d *= eps;
      c[2] = d*(eps2+6)/32;
      d *= eps;
      c[3] = 5*d/48;
      d *= eps;
      c[4] = 35*d/512;
      break;
    case 5:
      c[1] = d*(eps2*(eps2+2)+16)/32;
      d *= eps;
      c[2] = d*(eps2+6)/32;
      d *= eps;
      c[3] = d*(15*eps2+80)/768;
      d *= eps;
      c[4] = 35*d/512;
      d *= eps;
      c[5] = 63*d/1280;
      break;
    case 6:
      c[1] = d*(eps2*(eps2+2)+16)/32;
      d *= eps;
      c[2] = d*(eps2*(35*eps2+64)+384)/2048;
      d *= eps;
      c[3] = d*(15*eps2+80)/768;
      d *= eps;
      c[4] = d*(7*eps2+35)/512;
      d *= eps;
      c[5] = 63*d/1280;
      d *= eps;
      c[6] = 77*d/2048;
      break;
    case 7:
      c[1] = d*(eps2*(eps2*(41*eps2+64)+128)+1024)/2048;
      d *= eps;
      c[2] = d*(eps2*(35*eps2+64)+384)/2048;
      d *= eps;
      c[3] = d*(eps2*(69*eps2+120)+640)/6144;
      d *= eps;
      c[4] = d*(7*eps2+35)/512;
      d *= eps;
      c[5] = d*(105*eps2+504)/10240;
      d *= eps;
      c[6] = 77*d/2048;
      d *= eps;
      c[7] = 429*d/14336;
      break;
    case 8:
      c[1] = d*(eps2*(eps2*(41*eps2+64)+128)+1024)/2048;
      d *= eps;
      c[2] = d*(eps2*(eps2*(47*eps2+70)+128)+768)/4096;
      d *= eps;
      c[3] = d*(eps2*(69*eps2+120)+640)/6144;
      d *= eps;
      c[4] = d*(eps2*(133*eps2+224)+1120)/16384;
      d *= eps;
      c[5] = d*(105*eps2+504)/10240;
      d *= eps;
      c[6] = d*(33*eps2+154)/4096;
      d *= eps;
      c[7] = 429*d/14336;
      d *= eps;
      c[8] = 6435*d/262144;
      break;
    default:
      STATIC_ASSERT(nC2 >= 0 && nC2 <= 8, "Bad value of nC2");
    }
  }

  // The scale factor A3 = mean value of I3
  void Geodesic::A3coeff() throw() {
    switch (nA3) {
    case 0:
      break;
    case 1:
      _A3x[0] = 1;
      break;
    case 2:
      _A3x[0] = 1;
      _A3x[1] = -1/real(2);
      break;
    case 3:
      _A3x[0] = 1;
      _A3x[1] = (_n-1)/2;
      _A3x[2] = -1/real(4);
      break;
    case 4:
      _A3x[0] = 1;
      _A3x[1] = (_n-1)/2;
      _A3x[2] = (-_n-2)/8;
      _A3x[3] = -1/real(16);
      break;
    case 5:
      _A3x[0] = 1;
      _A3x[1] = (_n-1)/2;
      _A3x[2] = (_n*(3*_n-1)-2)/8;
      _A3x[3] = (-3*_n-1)/16;
      _A3x[4] = -3/real(64);
      break;
    case 6:
      _A3x[0] = 1;
      _A3x[1] = (_n-1)/2;
      _A3x[2] = (_n*(3*_n-1)-2)/8;
      _A3x[3] = ((-_n-3)*_n-1)/16;
      _A3x[4] = (-2*_n-3)/64;
      _A3x[5] = -3/real(128);
      break;
    case 7:
      _A3x[0] = 1;
      _A3x[1] = (_n-1)/2;
      _A3x[2] = (_n*(3*_n-1)-2)/8;
      _A3x[3] = (_n*(_n*(5*_n-1)-3)-1)/16;
      _A3x[4] = ((-10*_n-2)*_n-3)/64;
      _A3x[5] = (-5*_n-3)/128;
      _A3x[6] = -5/real(256);
      break;
    case 8:
      _A3x[0] = 1;
      _A3x[1] = (_n-1)/2;
      _A3x[2] = (_n*(3*_n-1)-2)/8;
      _A3x[3] = (_n*(_n*(5*_n-1)-3)-1)/16;
      _A3x[4] = (_n*((-5*_n-20)*_n-4)-6)/128;
      _A3x[5] = ((-5*_n-10)*_n-6)/256;
      _A3x[6] = (-15*_n-20)/1024;
      _A3x[7] = -25/real(2048);
      break;
    default:
      STATIC_ASSERT(nA3 >= 0 && nA3 <= 8, "Bad value of nA3");
    }
  }

  // The coefficients C3[l] in the Fourier expansion of B3
  void Geodesic::C3coeff() throw() {
    switch (nC3) {
    case 0:
      break;
    case 1:
      break;
    case 2:
      _C3x[0] = 1/real(4);
      break;
    case 3:
      _C3x[0] = (1-_n)/4;
      _C3x[1] = 1/real(8);
      _C3x[2] = 1/real(16);
      break;
    case 4:
      _C3x[0] = (1-_n)/4;
      _C3x[1] = 1/real(8);
      _C3x[2] = 3/real(64);
      _C3x[3] = (2-3*_n)/32;
      _C3x[4] = 3/real(64);
      _C3x[5] = 5/real(192);
      break;
    case 5:
      _C3x[0] = (1-_n)/4;
      _C3x[1] = (1-_n*_n)/8;
      _C3x[2] = (3*_n+3)/64;
      _C3x[3] = 5/real(128);
      _C3x[4] = ((_n-3)*_n+2)/32;
      _C3x[5] = (3-2*_n)/64;
      _C3x[6] = 3/real(128);
      _C3x[7] = (5-9*_n)/192;
      _C3x[8] = 3/real(128);
      _C3x[9] = 7/real(512);
      break;
    case 6:
      _C3x[0] = (1-_n)/4;
      _C3x[1] = (1-_n*_n)/8;
      _C3x[2] = ((3-_n)*_n+3)/64;
      _C3x[3] = (2*_n+5)/128;
      _C3x[4] = 3/real(128);
      _C3x[5] = ((_n-3)*_n+2)/32;
      _C3x[6] = ((-3*_n-2)*_n+3)/64;
      _C3x[7] = (_n+3)/128;
      _C3x[8] = 5/real(256);
      _C3x[9] = (_n*(5*_n-9)+5)/192;
      _C3x[10] = (9-10*_n)/384;
      _C3x[11] = 7/real(512);
      _C3x[12] = (7-14*_n)/512;
      _C3x[13] = 7/real(512);
      _C3x[14] = 21/real(2560);
      break;
    case 7:
      _C3x[0] = (1-_n)/4;
      _C3x[1] = (1-_n*_n)/8;
      _C3x[2] = (_n*((-5*_n-1)*_n+3)+3)/64;
      _C3x[3] = (_n*(2*_n+2)+5)/128;
      _C3x[4] = (11*_n+12)/512;
      _C3x[5] = 21/real(1024);
      _C3x[6] = ((_n-3)*_n+2)/32;
      _C3x[7] = (_n*(_n*(2*_n-3)-2)+3)/64;
      _C3x[8] = ((2-9*_n)*_n+6)/256;
      _C3x[9] = (_n+5)/256;
      _C3x[10] = 27/real(2048);
      _C3x[11] = (_n*((5-_n)*_n-9)+5)/192;
      _C3x[12] = ((-6*_n-10)*_n+9)/384;
      _C3x[13] = (21-4*_n)/1536;
      _C3x[14] = 3/real(256);
      _C3x[15] = (_n*(10*_n-14)+7)/512;
      _C3x[16] = (7-10*_n)/512;
      _C3x[17] = 9/real(1024);
      _C3x[18] = (21-45*_n)/2560;
      _C3x[19] = 9/real(1024);
      _C3x[20] = 11/real(2048);
      break;
    case 8:
      _C3x[0] = (1-_n)/4;
      _C3x[1] = (1-_n*_n)/8;
      _C3x[2] = (_n*((-5*_n-1)*_n+3)+3)/64;
      _C3x[3] = (_n*((2-2*_n)*_n+2)+5)/128;
      _C3x[4] = (_n*(3*_n+11)+12)/512;
      _C3x[5] = (10*_n+21)/1024;
      _C3x[6] = 243/real(16384);
      _C3x[7] = ((_n-3)*_n+2)/32;
      _C3x[8] = (_n*(_n*(2*_n-3)-2)+3)/64;
      _C3x[9] = (_n*((-6*_n-9)*_n+2)+6)/256;
      _C3x[10] = ((1-2*_n)*_n+5)/256;
      _C3x[11] = (69*_n+108)/8192;
      _C3x[12] = 187/real(16384);
      _C3x[13] = (_n*((5-_n)*_n-9)+5)/192;
      _C3x[14] = (_n*(_n*(10*_n-6)-10)+9)/384;
      _C3x[15] = ((-77*_n-8)*_n+42)/3072;
      _C3x[16] = (12-_n)/1024;
      _C3x[17] = 139/real(16384);
      _C3x[18] = (_n*((20-7*_n)*_n-28)+14)/1024;
      _C3x[19] = ((-7*_n-40)*_n+28)/2048;
      _C3x[20] = (72-43*_n)/8192;
      _C3x[21] = 127/real(16384);
      _C3x[22] = (_n*(75*_n-90)+42)/5120;
      _C3x[23] = (9-15*_n)/1024;
      _C3x[24] = 99/real(16384);
      _C3x[25] = (44-99*_n)/8192;
      _C3x[26] = 99/real(16384);
      _C3x[27] = 429/real(114688);
      break;
    default:
      STATIC_ASSERT(nC3 >= 0 && nC3 <= 8, "Bad value of nC3");
    }
  }
} // namespace GeographicLib

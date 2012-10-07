/**
 * \file Geodesic.cpp
 * \brief Implementation for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
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
 *   - sigma = arc length along great circle
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

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <iostream>

#if defined(_MSC_VER)
// Squelch warnings about potentially uninitialized local variables
#  pragma warning (disable: 4701)
#endif

namespace GeographicLib {

  using namespace std;

  // Underflow guard.  We require
  //   tiny_ * epsilon() > 0
  //   tiny_ + epsilon() == epsilon()
  const Math::real Geodesic::tiny_ = sqrt(numeric_limits<real>::min());
  const Math::real Geodesic::tol0_ = numeric_limits<real>::epsilon();
  // Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse case
  // 52.784459512564 0 -52.784459512563990912 179.634407464943777557
  // which otherwise failed for Visual Studio 10 (Release and Debug)
  const Math::real Geodesic::tol1_ = 200 * tol0_;
  const Math::real Geodesic::tol2_ = sqrt(numeric_limits<real>::epsilon());
  const Math::real Geodesic::xthresh_ = 1000 * tol2_;

  Geodesic::Geodesic(real a, real f)
    : _a(a)
    , _f(f <= 1 ? f : 1/f)
    , _f1(1 - _f)
    , _e2(_f * (2 - _f))
    , _ep2(_e2 / Math::sq(_f1))       // e2 / (1 - e2)
    , _n(_f / ( 2 - _f))
    , _b(_a * _f1)
    , _c2((Math::sq(_a) + Math::sq(_b) *
           (_e2 == 0 ? 1 :
            (_e2 > 0 ? Math::atanh(sqrt(_e2)) : atan(sqrt(-_e2))) /
            sqrt(abs(_e2))))/2) // authalic radius squared
      // The sig12 threshold for "really short"
    , _etol2(10 * tol2_ / max(real(0.1), sqrt(abs(_e2))))
  {
    if (!(Math::isfinite(_a) && _a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(Math::isfinite(_b) && _b > 0))
      throw GeographicErr("Minor radius is not positive");
    A3coeff();
    C3coeff();
    C4coeff();
  }

  const Geodesic Geodesic::WGS84(Constants::WGS84_a<real>(),
                                 Constants::WGS84_f<real>());

  Math::real Geodesic::SinCosSeries(bool sinp,
                                    real sinx, real cosx,
                                    const real c[], int n) throw() {
    // Evaluate
    // y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
    //            sum(c[i] * cos((2*i+1) * x), i, 0, n-1) :
    // using Clenshaw summation.  N.B. c[0] is unused for sin series
    // Approx operation count = (n + 5) mult and (2 * n + 2) add
    c += (n + sinp);            // Point to one beyond last element
    real
      ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
      y0 = n & 1 ? *--c : 0, y1 = 0;          // accumulators for sum
    // Now n is even
    n /= 2;
    while (n--) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + *--c;
      y0 = ar * y1 - y0 + *--c;
    }
    return sinp
      ? 2 * sinx * cosx * y0    // sin(2 * x) * y0
      : cosx * (y0 - y1);       // cos(x) * (y0 - y1)
  }

  GeodesicLine Geodesic::Line(real lat1, real lon1, real azi1, unsigned caps)
    const throw() {
    return GeodesicLine(*this, lat1, lon1, azi1, caps);
  }

  Math::real Geodesic::GenDirect(real lat1, real lon1, real azi1,
                                 bool arcmode, real s12_a12, unsigned outmask,
                                 real& lat2, real& lon2, real& azi2,
                                 real& s12, real& m12, real& M12, real& M21,
                                 real& S12) const throw() {
    return GeodesicLine(*this, lat1, lon1, azi1,
                        // Automatically supply DISTANCE_IN if necessary
                        outmask | (arcmode ? NONE : DISTANCE_IN))
      .                         // Note the dot!
      GenPosition(arcmode, s12_a12, outmask,
                  lat2, lon2, azi2, s12, m12, M12, M21, S12);
  }

  Math::real Geodesic::GenInverse(real lat1, real lon1, real lat2, real lon2,
                                  unsigned outmask,
                                  real& s12, real& azi1, real& azi2,
                                  real& m12, real& M12, real& M21, real& S12)
    const throw() {
    outmask &= OUT_ALL;
    lon1 = Math::AngNormalize(lon1);
    real lon12 = Math::AngNormalize(Math::AngNormalize(lon2) - lon1);
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

    real phi, sbet1, cbet1, sbet2, cbet2, s12x, m12x;

    phi = lat1 * Math::degree<real>();
    // Ensure cbet1 = +epsilon at poles
    sbet1 = _f1 * sin(phi);
    cbet1 = lat1 == -90 ? tiny_ : cos(phi);
    SinCosNorm(sbet1, cbet1);

    phi = lat2 * Math::degree<real>();
    // Ensure cbet2 = +epsilon at poles
    sbet2 = _f1 * sin(phi);
    cbet2 = abs(lat2) == 90 ? tiny_ : cos(phi);
    SinCosNorm(sbet2, cbet2);

    // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
    // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
    // a better measure.  This logic is used in assigning calp2 in Lambda12.
    // Sometimes these quantities vanish and in that case we force bet2 = +/-
    // bet1 exactly.  An example where is is necessary is the inverse problem
    // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
    // which failed with Visual Studio 10 (Release and Debug)

    if (cbet1 < -sbet1) {
      if (cbet2 == cbet1)
        sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
    } else {
      if (abs(sbet2) == -sbet1)
        cbet2 = cbet1;
    }

    real
      lam12 = lon12 * Math::degree<real>(),
      slam12 = lon12 == 180 ? 0 : sin(lam12),
      clam12 = cos(lam12);      // lon12 == 90 isn't interesting

    real a12, sig12, calp1, salp1, calp2, salp2;
    // index zero elements of these arrays are unused
    real C1a[nC1_ + 1], C2a[nC2_ + 1], C3a[nC3_];

    bool meridian = lat1 == -90 || slam12 == 0;

    if (meridian) {

      // Endpoints are on a single full meridian, so the geodesic might lie on
      // a meridian.

      calp1 = clam12; salp1 = slam12; // Head to the target longitude
      calp2 = 1; salp2 = 0;           // At the target we're heading north

      real
        // tan(bet) = tan(sig) * cos(alp)
        ssig1 = sbet1, csig1 = calp1 * cbet1,
        ssig2 = sbet2, csig2 = calp2 * cbet2;

      // sig12 = sig2 - sig1
      sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, real(0)),
                    csig1 * csig2 + ssig1 * ssig2);
      {
        real dummy;
        Lengths(_n, sig12, ssig1, csig1, ssig2, csig2,
                cbet1, cbet2, s12x, m12x, dummy,
                (outmask & GEODESICSCALE) != 0U, M12, M21, C1a, C2a);
      }
      // Add the check for sig12 since zero length geodesics might yield m12 <
      // 0.  Test case was
      //
      //    echo 20.001 0 20.001 0 | Geod -i
      //
      // In fact, we will have sig12 > pi/2 for meridional geodesic which is
      // not a shortest path.
      if (sig12 < 1 || m12x >= 0) {
        m12x *= _a;
        s12x *= _b;
        a12 = sig12 / Math::degree<real>();
      } else
        // m12 < 0, i.e., prolate and too close to anti-podal
        meridian = false;
    }

    real omg12;
    if (!meridian &&
        sbet1 == 0 &&   // and sbet2 == 0
        // Mimic the way Lambda12 works with calp1 = 0
        (_f <= 0 || lam12 <= Math::pi<real>() - _f * Math::pi<real>())) {

      // Geodesic runs along equator
      calp1 = calp2 = 0; salp1 = salp2 = 1;
      s12x = _a * lam12;
      m12x = _b * sin(lam12 / _f1);
      if (outmask & GEODESICSCALE)
        M12 = M21 = cos(lam12 / _f1);
      a12 = lon12 / _f1;
      sig12 = omg12 = lam12 / _f1;

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
        real wm = sqrt(1 - _e2 * Math::sq((cbet1 + cbet2) / 2));
        s12x = sig12 * _a * wm;
        m12x = Math::sq(wm) * _a / _f1 * sin(sig12 * _f1 / wm);
        if (outmask & GEODESICSCALE)
          M12 = M21 = cos(sig12 * _f1 / wm);
        a12 = sig12 / Math::degree<real>();
        omg12 = lam12 / wm;
      } else {

        // Newton's method
        real ssig1, csig1, ssig2, csig2, eps;
        real ov = 0;
        unsigned numit = 0;
        // Bracketing values for bisection method
        real salp1a = 0, calp1a = 1, salp1b = 0, calp1b = -1;
        for (unsigned trip = 0; numit < maxit_; ++numit) {
          real dv;
          real v = Lambda12(sbet1, cbet1, sbet2, cbet2, salp1, calp1,
                            salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                            eps, omg12, trip < 1, dv, C1a, C2a, C3a) - lam12;
          // Update bracketing values
          if (v >= 0 && calp1 > calp1b) {
            salp1b = salp1; calp1b = calp1;
          } else if (v <= 0 && calp1 < calp1a) {
            salp1a = salp1; calp1a = calp1;
          }
          if (!(abs(v) > tiny_) || !(trip < 1)) {
            if (!(abs(v) <= max(tol1_, ov)))
              numit = maxit_;
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
          // 200 * epsilon.  The second takes credit for an anticipated
          // reduction in abs(v) by v/ov (due to the latest update in alp1) and
          // checks this against epsilon.
          if (!(abs(v) >= tol1_ && Math::sq(v) >= ov * tol0_)) ++trip;
          ov = abs(v);
        }
        if (numit >= maxit_) {
          // Resort to the safer bisection method
          for (unsigned i = 0; i < bisection_; ++i) {
            ++numit;
            salp1 = (salp1a + salp1b)/2;
            calp1 = (calp1a + calp1b)/2;
            SinCosNorm(salp1, calp1);
            if (abs(salp1a - salp1b) < tol0_ && calp1a - calp1b < tol0_)
              break;
            real
              dummy,
              v = Lambda12(sbet1, cbet1, sbet2, cbet2, salp1, calp1,
                           salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                           eps, omg12, false, dummy, C1a, C2a, C3a) - lam12;
            if (false)
              cerr << i << " "
                   << calp1a << " " << calp1 << " " << calp1b << " "
                   << v << "\n";
            // Be more tolerant on error.  It is approximately 1 ulp for a
            // number in [0, pi].
            if (abs(v) <= 2 * tol0_) break;
            if (v > 0) {
              salp1b = salp1; calp1b = calp1;
            } else {
              salp1a = salp1; calp1a = calp1;
            }
          }
          cerr << numit - maxit_ << "\n";
        }

        if (numit >= maxit_ + bisection_) {
          // Signal failure.
          if (outmask & DISTANCE)
            s12 = Math::NaN<real>();
          if (outmask & AZIMUTH)
            azi1 = azi2 = Math::NaN<real>();
          if (outmask & REDUCEDLENGTH)
            m12 = Math::NaN<real>();
          if (outmask & GEODESICSCALE)
            M12 = M21 = Math::NaN<real>();
          if (outmask & AREA)
            S12 = Math::NaN<real>();
          return Math::NaN<real>();
        }

        {
          real dummy;
          Lengths(eps, sig12, ssig1, csig1, ssig2, csig2,
                  cbet1, cbet2, s12x, m12x, dummy,
                  (outmask & GEODESICSCALE) != 0U, M12, M21, C1a, C2a);
        }
        m12x *= _a;
        s12x *= _b;
        a12 = sig12 / Math::degree<real>();
        omg12 = lam12 - omg12;
      }
    }

    if (outmask & DISTANCE)
      s12 = 0 + s12x;           // Convert -0 to 0

    if (outmask & REDUCEDLENGTH)
      m12 = 0 + m12x;           // Convert -0 to 0

    if (outmask & AREA) {
      real
        // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
        salp0 = salp1 * cbet1,
        calp0 = Math::hypot(calp1, salp1 * sbet1); // calp0 > 0
      real alp12;
      if (calp0 != 0 && salp0 != 0) {
        real
          // From Lambda12: tan(bet) = tan(sig) * cos(alp)
          ssig1 = sbet1, csig1 = calp1 * cbet1,
          ssig2 = sbet2, csig2 = calp2 * cbet2,
          k2 = Math::sq(calp0) * _ep2,
          // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
          A4 = Math::sq(_a) * calp0 * salp0 * _e2;
        SinCosNorm(ssig1, csig1);
        SinCosNorm(ssig2, csig2);
        real C4a[nC4_];
        C4f(k2, C4a);
        real
          B41 = SinCosSeries(false, ssig1, csig1, C4a, nC4_),
          B42 = SinCosSeries(false, ssig2, csig2, C4a, nC4_);
        S12 = A4 * (B42 - B41);
      } else
        // Avoid problems with indeterminate sig1, sig2 on equator
        S12 = 0;

      if (!meridian &&
          omg12 < real(0.75) * Math::pi<real>() && // Long difference too big
          sbet2 - sbet1 < real(1.75)) {            // Lat difference too big
        // Use tan(Gamma/2) = tan(omg12/2)
        // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
        // with tan(x/2) = sin(x)/(1+cos(x))
        real
          somg12 = sin(omg12), domg12 = 1 + cos(omg12),
          dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
        alp12 = 2 * atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                           domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
      } else {
        // alp12 = alp2 - alp1, used in atan2 so no need to normalize
        real
          salp12 = salp2 * calp1 - calp2 * salp1,
          calp12 = calp2 * calp1 + salp2 * salp1;
        // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
        // salp12 = -0 and alp12 = -180.  However this depends on the sign
        // being attached to 0 correctly.  The following ensures the correct
        // behavior.
        if (salp12 == 0 && calp12 < 0) {
          salp12 = tiny_ * calp1;
          calp12 = -1;
        }
        alp12 = atan2(salp12, calp12);
      }
      S12 += _c2 * alp12;
      S12 *= swapp * lonsign * latsign;
      // Convert -0 to 0
      S12 += 0;
    }

    // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
    if (swapp < 0) {
      swap(salp1, salp2);
      swap(calp1, calp2);
      if (outmask & GEODESICSCALE)
        swap(M12, M21);
    }

    salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
    salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

    if (outmask & AZIMUTH) {
      // minus signs give range [-180, 180). 0- converts -0 to +0.
      azi1 = 0 - atan2(-salp1, calp1) / Math::degree<real>();
      azi2 = 0 - atan2(-salp2, calp2) / Math::degree<real>();
    }

    // Returned value in [0, 180]
    return a12;
  }

  void Geodesic::Lengths(real eps, real sig12,
                         real ssig1, real csig1, real ssig2, real csig2,
                         real cbet1, real cbet2,
                         real& s12b, real& m12a, real& m0,
                         bool scalep, real& M12, real& M21,
                         // Scratch areas of the right size
                         real C1a[], real C2a[]) const throw() {
    // Return m12a = (reduced length)/_a; also calculate s12b = distance/_b,
    // and m0 = coefficient of secular term in expression for reduced length.
    C1f(eps, C1a);
    C2f(eps, C2a);
    real
      A1m1 = A1m1f(eps),
      AB1 = (1 + A1m1) * (SinCosSeries(true, ssig2, csig2, C1a, nC1_) -
                          SinCosSeries(true, ssig1, csig1, C1a, nC1_)),
      A2m1 = A2m1f(eps),
      AB2 = (1 + A2m1) * (SinCosSeries(true, ssig2, csig2, C2a, nC2_) -
                          SinCosSeries(true, ssig1, csig1, C2a, nC2_)),
      cbet1sq = Math::sq(cbet1),
      cbet2sq = Math::sq(cbet2),
      w1 = sqrt(1 - _e2 * cbet1sq),
      w2 = sqrt(1 - _e2 * cbet2sq),
      // Make sure it's OK to have repeated dummy arguments
      m0x = A1m1 - A2m1,
      J12 = m0x * sig12 + (AB1 - AB2);
    m0 = m0x;
    // Missing a factor of _a.
    // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure accurate
    // cancellation in the case of coincident points.
    m12a = (w2 * (csig1 * ssig2) - w1 * (ssig1 * csig2))
      - _f1 * csig1 * csig2 * J12;
    // Missing a factor of _b
    s12b = (1 + A1m1) * sig12 + AB1;
    if (scalep) {
      real csig12 = csig1 * csig2 + ssig1 * ssig2;
      J12 *= _f1;
      M12 = csig12 + (_e2 * (cbet1sq - cbet2sq) * ssig2 / (w1 + w2)
                      - csig2 * J12) * ssig1 / w1;
      M21 = csig12 - (_e2 * (cbet1sq - cbet2sq) * ssig1 / (w1 + w2)
                      - csig1 * J12) * ssig2 / w2;
    }
  }

  Math::real Geodesic::Astroid(real x, real y) throw() {
    // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
    // This solution is adapted from Geocentric::Reverse.
    real k;
    real
      p = Math::sq(x),
      q = Math::sq(y),
      r = (p + q - 1) / 6;
    if ( !(q == 0 && r <= 0) ) {
      real
        // Avoid possible division by zero when r = 0 by multiplying equations
        // for s and t by r^3 and r, resp.
        S = p * q / 4,            // S = r^3 * s
        r2 = Math::sq(r),
        r3 = r * r2,
        // The discrimant of the quadratic equation for T3.  This is zero on
        // the evolute curve p^(1/3)+q^(1/3) = 1
        disc = S * (S + 2 * r3);
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
        v = sqrt(Math::sq(u) + q),    // guaranteed positive
        // Avoid loss of accuracy when u < 0.
        uv = u < 0 ? q / (v - u) : u + v, // u+v, guaranteed positive
        w = (uv - q) / (2 * v);           // positive?
      // Rearrange expression for k to avoid loss of accuracy due to
      // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
      k = uv / (sqrt(uv + Math::sq(w)) + w);   // guaranteed positive
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
      cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
#if defined(__GNUC__) && __GNUC__ == 4 && \
  (__GNUC_MINOR__ < 6 || defined(__MINGW32__))
    // Volatile declaration needed to fix inverse cases
    // 88.202499451857 0 -88.202499451857 179.981022032992859592
    // 89.262080389218 0 -89.262080389218 179.992207982775375662
    // 89.333123580033 0 -89.333123580032997687 179.99295812360148422
    // which otherwise fail with g++ 4.4.4 x86 -O3 (Linux)
    // and g++ 4.4.0 (mingw) and g++ 4.6.1 (tdm mingw).
    real sbet12a;
    {
      volatile real xx1 = sbet2 * cbet1;
      volatile real xx2 = cbet2 * sbet1;
      sbet12a = xx1 + xx2;
    }
#else
    real sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
#endif
    bool shortline = cbet12 >= 0 && sbet12 < real(0.5) &&
      lam12 <= Math::pi<real>() / 6;
    real
      omg12 = (!shortline ? lam12 :
               lam12 / sqrt(1 - _e2 * Math::sq((cbet1 + cbet2) / 2))),
      somg12 = sin(omg12), comg12 = cos(omg12);

    salp1 = cbet2 * somg12;
    calp1 = comg12 >= 0 ?
      sbet12 + cbet2 * sbet1 * Math::sq(somg12) / (1 + comg12) :
      sbet12a - cbet2 * sbet1 * Math::sq(somg12) / (1 - comg12);

    real
      ssig12 = Math::hypot(salp1, calp1),
      csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

    if (shortline && ssig12 < _etol2) {
      // really short lines
      salp2 = cbet1 * somg12;
      calp2 = sbet12 - cbet1 * sbet2 * Math::sq(somg12) / (1 + comg12);
      SinCosNorm(salp2, calp2);
      // Set return value
      sig12 = atan2(ssig12, csig12);
    } else if (csig12 >= 0 ||
               ssig12 >= 3 * abs(_f) * Math::pi<real>() * Math::sq(cbet1)) {
      // Nothing to do, zeroth order spherical approximation is OK
    } else {
      // Scale lam12 and bet2 to x, y coordinate system where antipodal point
      // is at origin and singular point is at y = 0, x = -1.
      real y, lamscale, betscale;
      // Volatile declaration needed to fix inverse case
      // 56.320923501171 0 -56.320923501171 179.664747671772880215
      // which otherwise fails with g++ 4.4.4 x86 -O3
      volatile real x;
      if (_f >= 0) {            // In fact f == 0 does not get here
        // x = dlong, y = dlat
        {
          real
            k2 = Math::sq(sbet1) * _ep2,
            eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
          lamscale = _f * cbet1 * A3f(eps) * Math::pi<real>();
        }
        betscale = lamscale * cbet1;

        x = (lam12 - Math::pi<real>()) / lamscale;
        y = sbet12a / betscale;
      } else {                  // _f < 0
        // x = dlat, y = dlong
        real
          cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
          bet12a = atan2(sbet12a, cbet12a);
        real m12a, m0, dummy;
        // In the case of lon12 = 180, this repeats a calculation made in
        // Inverse.
        Lengths(_n, Math::pi<real>() + bet12a, sbet1, -cbet1, sbet2, cbet2,
                cbet1, cbet2, dummy, m12a, m0, false,
                dummy, dummy, C1a, C2a);
        x = -1 + m12a/(_f1 * cbet1 * cbet2 * m0 * Math::pi<real>());
        betscale = x < -real(0.01) ? sbet12a / x :
          -_f * Math::sq(cbet1) * Math::pi<real>();
        lamscale = betscale / cbet1;
        y = (lam12 - Math::pi<real>()) / lamscale;
      }

      if (y > -tol1_ && x > -1 - xthresh_) {
        // strip near cut
        if (_f >= 0) {
          salp1 = min(real(1), -real(x)); calp1 = - sqrt(1 - Math::sq(salp1));
        } else {
          calp1 = max(real(x > -tol1_ ? 0 : -1), real(x));
          salp1 = sqrt(1 - Math::sq(calp1));
        }
      } else {
        // Estimate alp1, by solving the astroid problem.
        //
        // Could estimate alpha1 = theta + pi/2, directly, i.e.,
        //   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
        //   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
        //
        // However, it's better to estimate omg12 from astroid and use
        // spherical formula to compute alp1.  This reduces the mean number of
        // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
        // (min 0 max 5).  The changes in the number of iterations are as
        // follows:
        //
        // change percent
        //    1       5
        //    0      78
        //   -1      16
        //   -2       0.6
        //   -3       0.04
        //   -4       0.002
        //
        // The histogram of iterations is (m = number of iterations estimating
        // alp1 directly, n = number of iterations estimating via omg12, total
        // number of trials = 148605):
        //
        //  iter    m      n
        //    0   148    186
        //    1 13046  13845
        //    2 93315 102225
        //    3 36189  32341
        //    4  5396      7
        //    5   455      1
        //    6    56      0
        //
        // Because omg12 is near pi, estimate work with omg12a = pi - omg12
        real k = Astroid(x, y);
        real
          omg12a = lamscale * ( _f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k ),
          somg12 = sin(omg12a), comg12 = -cos(omg12a);
        // Update spherical estimate of alp1 using omg12 instead of lam12
        salp1 = cbet2 * somg12;
        calp1 = sbet12a - cbet2 * sbet1 * Math::sq(somg12) / (1 - comg12);
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
                                real& eps, real& domg12,
                                bool diffp, real& dlam12,
                                // Scratch areas of the right size
                                real C1a[], real C2a[], real C3a[]) const
    throw() {

    if (sbet1 == 0 && calp1 == 0)
      // Break degeneracy of equatorial line.  This case has already been
      // handled.
      calp1 = -tiny_;

    real
      // sin(alp1) * cos(bet1) = sin(alp0)
      salp0 = salp1 * cbet1,
      calp0 = Math::hypot(calp1, salp1 * sbet1); // calp0 > 0

    real somg1, comg1, somg2, comg2, omg12, lam12;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
    ssig1 = sbet1; somg1 = salp0 * sbet1;
    csig1 = comg1 = calp1 * cbet1;
    SinCosNorm(ssig1, csig1);
    // SinCosNorm(somg1, comg1); -- don't need to normalize!

    // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
    // about this case, since this can yield singularities in the Newton
    // iteration.
    // sin(alp2) * cos(bet2) = sin(alp0)
    salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
    // calp2 = sqrt(1 - sq(salp2))
    //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
    // and subst for calp0 and rearrange to give (choose positive sqrt
    // to give alp2 in [0, pi/2]).
    calp2 = cbet2 != cbet1 || abs(sbet2) != -sbet1 ?
      sqrt(Math::sq(calp1 * cbet1) +
           (cbet1 < -sbet1 ?
            (cbet2 - cbet1) * (cbet1 + cbet2) :
            (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
      abs(calp1);
    // tan(bet2) = tan(sig2) * cos(alp2)
    // tan(omg2) = sin(alp0) * tan(sig2).
    ssig2 = sbet2; somg2 = salp0 * sbet2;
    csig2 = comg2 = calp2 * cbet2;
    SinCosNorm(ssig2, csig2);
    // SinCosNorm(somg2, comg2); -- don't need to normalize!

    // sig12 = sig2 - sig1, limit to [0, pi]
    sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, real(0)),
                  csig1 * csig2 + ssig1 * ssig2);

    // omg12 = omg2 - omg1, limit to [0, pi]
    omg12 = atan2(max(comg1 * somg2 - somg1 * comg2, real(0)),
                  comg1 * comg2 + somg1 * somg2);
    real B312, h0;
    real k2 = Math::sq(calp0) * _ep2;
    eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
    C3f(eps, C3a);
    B312 = (SinCosSeries(true, ssig2, csig2, C3a, nC3_-1) -
            SinCosSeries(true, ssig1, csig1, C3a, nC3_-1));
    h0 = -_f * A3f(eps);
    domg12 = salp0 * h0 * (sig12 + B312);
    lam12 = omg12 + domg12;

    if (diffp) {
      if (calp2 == 0)
        dlam12 = - 2 * sqrt(1 - _e2 * Math::sq(cbet1)) / sbet1;
      else {
        real dummy;
        Lengths(eps, sig12, ssig1, csig1, ssig2, csig2,
                cbet1, cbet2, dummy, dlam12, dummy,
                false, dummy, dummy, C1a, C2a);
        dlam12 /= calp2 * cbet2;
      }
    }

    return lam12;
  }

  Math::real Geodesic::A3f(real eps) const throw() {
    // Evaluation sum(_A3c[k] * eps^k, k, 0, nA3x_-1) by Horner's method
    real v = 0;
    for (int i = nA3x_; i; )
      v = eps * v + _A3x[--i];
    return v;
  }

  void Geodesic::C3f(real eps, real c[]) const throw() {
    // Evaluation C3 coeffs by Horner's method
    // Elements c[1] thru c[nC3_ - 1] are set
    for (int j = nC3x_, k = nC3_ - 1; k; ) {
      real t = 0;
      for (int i = nC3_ - k; i; --i)
        t = eps * t + _C3x[--j];
      c[k--] = t;
    }

    real mult = 1;
    for (int k = 1; k < nC3_; ) {
      mult *= eps;
      c[k++] *= mult;
    }
  }

  void Geodesic::C4f(real k2, real c[]) const throw() {
    // Evaluation C4 coeffs by Horner's method
    // Elements c[0] thru c[nC4_ - 1] are set
    for (int j = nC4x_, k = nC4_; k; ) {
      real t = 0;
      for (int i = nC4_ - k + 1; i; --i)
        t = k2 * t + _C4x[--j];
      c[--k] = t;
    }

    real mult = 1;
    for (int k = 1; k < nC4_; ) {
      mult *= k2;
      c[k++] *= mult;
    }
  }

  // Generated by Maxima on 2010-09-04 10:26:17-04:00

  // The scale factor A1-1 = mean value of (d/dsigma)I1 - 1
  Math::real Geodesic::A1m1f(real eps) throw() {
    real
      eps2 = Math::sq(eps),
      t;
    switch (nA1_/2) {
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
      STATIC_ASSERT(nA1_ >= 0 && nA1_ <= 8, "Bad value of nA1_");
      t = 0;
    }
    return (t + eps) / (1 - eps);
  }

  // The coefficients C1[l] in the Fourier expansion of B1
  void Geodesic::C1f(real eps, real c[]) throw() {
    real
      eps2 = Math::sq(eps),
      d = eps;
    switch (nC1_) {
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
      STATIC_ASSERT(nC1_ >= 0 && nC1_ <= 8, "Bad value of nC1_");
    }
  }

  // The coefficients C1p[l] in the Fourier expansion of B1p
  void Geodesic::C1pf(real eps, real c[]) throw() {
    real
      eps2 = Math::sq(eps),
      d = eps;
    switch (nC1p_) {
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
      STATIC_ASSERT(nC1p_ >= 0 && nC1p_ <= 8, "Bad value of nC1p_");
    }
  }

  // The scale factor A2-1 = mean value of (d/dsigma)I2 - 1
  Math::real Geodesic::A2m1f(real eps) throw() {
    real
      eps2 = Math::sq(eps),
      t;
    switch (nA2_/2) {
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
      STATIC_ASSERT(nA2_ >= 0 && nA2_ <= 8, "Bad value of nA2_");
      t = 0;
    }
    return t * (1 - eps) - eps;
  }

  // The coefficients C2[l] in the Fourier expansion of B2
  void Geodesic::C2f(real eps, real c[]) throw() {
    real
      eps2 = Math::sq(eps),
      d = eps;
    switch (nC2_) {
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
      STATIC_ASSERT(nC2_ >= 0 && nC2_ <= 8, "Bad value of nC2_");
    }
  }

  // The scale factor A3 = mean value of (d/dsigma)I3
  void Geodesic::A3coeff() throw() {
    switch (nA3_) {
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
      STATIC_ASSERT(nA3_ >= 0 && nA3_ <= 8, "Bad value of nA3_");
    }
  }

  // The coefficients C3[l] in the Fourier expansion of B3
  void Geodesic::C3coeff() throw() {
    switch (nC3_) {
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
      STATIC_ASSERT(nC3_ >= 0 && nC3_ <= 8, "Bad value of nC3_");
    }
  }

  // The coefficients C4[l] in the Fourier expansion of I4
  void Geodesic::C4coeff() throw() {
    switch (nC4_) {
    case 0:
      break;
    case 1:
      _C4x[0] = 2/real(3);
      break;
    case 2:
      _C4x[0] = (10-_ep2)/15;
      _C4x[1] = -1/real(20);
      _C4x[2] = 1/real(180);
      break;
    case 3:
      _C4x[0] = (_ep2*(4*_ep2-7)+70)/105;
      _C4x[1] = (4*_ep2-7)/140;
      _C4x[2] = 1/real(42);
      _C4x[3] = (7-4*_ep2)/1260;
      _C4x[4] = -1/real(252);
      _C4x[5] = 1/real(2100);
      break;
    case 4:
      _C4x[0] = (_ep2*((12-8*_ep2)*_ep2-21)+210)/315;
      _C4x[1] = ((12-8*_ep2)*_ep2-21)/420;
      _C4x[2] = (3-2*_ep2)/126;
      _C4x[3] = -1/real(72);
      _C4x[4] = (_ep2*(8*_ep2-12)+21)/3780;
      _C4x[5] = (2*_ep2-3)/756;
      _C4x[6] = 1/real(360);
      _C4x[7] = (3-2*_ep2)/6300;
      _C4x[8] = -1/real(1800);
      _C4x[9] = 1/real(17640);
      break;
    case 5:
      _C4x[0] = (_ep2*(_ep2*(_ep2*(64*_ep2-88)+132)-231)+2310)/3465;
      _C4x[1] = (_ep2*(_ep2*(64*_ep2-88)+132)-231)/4620;
      _C4x[2] = (_ep2*(16*_ep2-22)+33)/1386;
      _C4x[3] = (8*_ep2-11)/792;
      _C4x[4] = 1/real(110);
      _C4x[5] = (_ep2*((88-64*_ep2)*_ep2-132)+231)/41580;
      _C4x[6] = ((22-16*_ep2)*_ep2-33)/8316;
      _C4x[7] = (11-8*_ep2)/3960;
      _C4x[8] = -1/real(495);
      _C4x[9] = (_ep2*(16*_ep2-22)+33)/69300;
      _C4x[10] = (8*_ep2-11)/19800;
      _C4x[11] = 1/real(1925);
      _C4x[12] = (11-8*_ep2)/194040;
      _C4x[13] = -1/real(10780);
      _C4x[14] = 1/real(124740);
      break;
    case 6:
      _C4x[0] = (_ep2*(_ep2*(_ep2*((832-640*_ep2)*_ep2-1144)+1716)-3003)+
                30030)/45045;
      _C4x[1] = (_ep2*(_ep2*((832-640*_ep2)*_ep2-1144)+1716)-3003)/60060;
      _C4x[2] = (_ep2*((208-160*_ep2)*_ep2-286)+429)/18018;
      _C4x[3] = ((104-80*_ep2)*_ep2-143)/10296;
      _C4x[4] = (13-10*_ep2)/1430;
      _C4x[5] = -1/real(156);
      _C4x[6] = (_ep2*(_ep2*(_ep2*(640*_ep2-832)+1144)-1716)+3003)/540540;
      _C4x[7] = (_ep2*(_ep2*(160*_ep2-208)+286)-429)/108108;
      _C4x[8] = (_ep2*(80*_ep2-104)+143)/51480;
      _C4x[9] = (10*_ep2-13)/6435;
      _C4x[10] = 5/real(3276);
      _C4x[11] = (_ep2*((208-160*_ep2)*_ep2-286)+429)/900900;
      _C4x[12] = ((104-80*_ep2)*_ep2-143)/257400;
      _C4x[13] = (13-10*_ep2)/25025;
      _C4x[14] = -1/real(2184);
      _C4x[15] = (_ep2*(80*_ep2-104)+143)/2522520;
      _C4x[16] = (10*_ep2-13)/140140;
      _C4x[17] = 5/real(45864);
      _C4x[18] = (13-10*_ep2)/1621620;
      _C4x[19] = -1/real(58968);
      _C4x[20] = 1/real(792792);
      break;
    case 7:
      _C4x[0] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(512*_ep2-640)+832)-1144)+1716)-
                3003)+30030)/45045;
      _C4x[1] = (_ep2*(_ep2*(_ep2*(_ep2*(512*_ep2-640)+832)-1144)+1716)-
                3003)/60060;
      _C4x[2] = (_ep2*(_ep2*(_ep2*(128*_ep2-160)+208)-286)+429)/18018;
      _C4x[3] = (_ep2*(_ep2*(64*_ep2-80)+104)-143)/10296;
      _C4x[4] = (_ep2*(8*_ep2-10)+13)/1430;
      _C4x[5] = (4*_ep2-5)/780;
      _C4x[6] = 1/real(210);
      _C4x[7] = (_ep2*(_ep2*(_ep2*((640-512*_ep2)*_ep2-832)+1144)-1716)+
                3003)/540540;
      _C4x[8] = (_ep2*(_ep2*((160-128*_ep2)*_ep2-208)+286)-429)/108108;
      _C4x[9] = (_ep2*((80-64*_ep2)*_ep2-104)+143)/51480;
      _C4x[10] = ((10-8*_ep2)*_ep2-13)/6435;
      _C4x[11] = (5-4*_ep2)/3276;
      _C4x[12] = -1/real(840);
      _C4x[13] = (_ep2*(_ep2*(_ep2*(128*_ep2-160)+208)-286)+429)/900900;
      _C4x[14] = (_ep2*(_ep2*(64*_ep2-80)+104)-143)/257400;
      _C4x[15] = (_ep2*(8*_ep2-10)+13)/25025;
      _C4x[16] = (4*_ep2-5)/10920;
      _C4x[17] = 1/real(2520);
      _C4x[18] = (_ep2*((80-64*_ep2)*_ep2-104)+143)/2522520;
      _C4x[19] = ((10-8*_ep2)*_ep2-13)/140140;
      _C4x[20] = (5-4*_ep2)/45864;
      _C4x[21] = -1/real(8820);
      _C4x[22] = (_ep2*(8*_ep2-10)+13)/1621620;
      _C4x[23] = (4*_ep2-5)/294840;
      _C4x[24] = 1/real(41580);
      _C4x[25] = (5-4*_ep2)/3963960;
      _C4x[26] = -1/real(304920);
      _C4x[27] = 1/real(4684680);
      break;
    case 8:
      _C4x[0] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((8704-7168*_ep2)*_ep2-10880)+
                14144)-19448)+29172)-51051)+510510)/765765;
      _C4x[1] = (_ep2*(_ep2*(_ep2*(_ep2*((8704-7168*_ep2)*_ep2-10880)+14144)-
                19448)+29172)-51051)/1021020;
      _C4x[2] = (_ep2*(_ep2*(_ep2*((2176-1792*_ep2)*_ep2-2720)+3536)-4862)+
                7293)/306306;
      _C4x[3] = (_ep2*(_ep2*((1088-896*_ep2)*_ep2-1360)+1768)-2431)/175032;
      _C4x[4] = (_ep2*((136-112*_ep2)*_ep2-170)+221)/24310;
      _C4x[5] = ((68-56*_ep2)*_ep2-85)/13260;
      _C4x[6] = (17-14*_ep2)/3570;
      _C4x[7] = -1/real(272);
      _C4x[8] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(7168*_ep2-8704)+10880)-14144)+
                19448)-29172)+51051)/9189180;
      _C4x[9] = (_ep2*(_ep2*(_ep2*(_ep2*(1792*_ep2-2176)+2720)-3536)+4862)-
                7293)/1837836;
      _C4x[10] = (_ep2*(_ep2*(_ep2*(896*_ep2-1088)+1360)-1768)+2431)/875160;
      _C4x[11] = (_ep2*(_ep2*(112*_ep2-136)+170)-221)/109395;
      _C4x[12] = (_ep2*(56*_ep2-68)+85)/55692;
      _C4x[13] = (14*_ep2-17)/14280;
      _C4x[14] = 7/real(7344);
      _C4x[15] = (_ep2*(_ep2*(_ep2*((2176-1792*_ep2)*_ep2-2720)+3536)-4862)+
                7293)/15315300;
      _C4x[16] = (_ep2*(_ep2*((1088-896*_ep2)*_ep2-1360)+1768)-2431)/4375800;
      _C4x[17] = (_ep2*((136-112*_ep2)*_ep2-170)+221)/425425;
      _C4x[18] = ((68-56*_ep2)*_ep2-85)/185640;
      _C4x[19] = (17-14*_ep2)/42840;
      _C4x[20] = -7/real(20400);
      _C4x[21] = (_ep2*(_ep2*(_ep2*(896*_ep2-1088)+1360)-1768)+2431)/42882840;
      _C4x[22] = (_ep2*(_ep2*(112*_ep2-136)+170)-221)/2382380;
      _C4x[23] = (_ep2*(56*_ep2-68)+85)/779688;
      _C4x[24] = (14*_ep2-17)/149940;
      _C4x[25] = 1/real(8976);
      _C4x[26] = (_ep2*((136-112*_ep2)*_ep2-170)+221)/27567540;
      _C4x[27] = ((68-56*_ep2)*_ep2-85)/5012280;
      _C4x[28] = (17-14*_ep2)/706860;
      _C4x[29] = -7/real(242352);
      _C4x[30] = (_ep2*(56*_ep2-68)+85)/67387320;
      _C4x[31] = (14*_ep2-17)/5183640;
      _C4x[32] = 7/real(1283568);
      _C4x[33] = (17-14*_ep2)/79639560;
      _C4x[34] = -1/real(1516944);
      _C4x[35] = 1/real(26254800);
      break;
    default:
      STATIC_ASSERT(nC4_ >= 0 && nC4_ <= 8, "Bad value of nC4_");
    }
  }

} // namespace GeographicLib

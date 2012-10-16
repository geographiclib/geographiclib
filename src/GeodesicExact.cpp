/**
 * \file GeodesicExact.cpp
 * \brief Implementation for GeographicLib::GeodesicExact class
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
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

#include <GeographicLib/GeodesicExact.hpp>
#include <GeographicLib/GeodesicLineExact.hpp>

#if defined(_MSC_VER)
// Squelch warnings about potentially uninitialized local variables
#  pragma warning (disable: 4701)
#endif

namespace GeographicLib {

  using namespace std;

  // Underflow guard.  We require
  //   tiny_ * epsilon() > 0
  //   tiny_ + epsilon() == epsilon()
  const Math::real GeodesicExact::tiny_ = sqrt(numeric_limits<real>::min());
  const Math::real GeodesicExact::tol0_ = numeric_limits<real>::epsilon();
  // Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse case
  // 52.784459512564 0 -52.784459512563990912 179.634407464943777557
  // which otherwise failed for Visual Studio 10 (Release and Debug)
  const Math::real GeodesicExact::tol1_ = 200 * tol0_;
  const Math::real GeodesicExact::tol2_ = sqrt(numeric_limits<real>::epsilon());
  const Math::real GeodesicExact::xthresh_ = 1000 * tol2_;

  GeodesicExact::GeodesicExact(real a, real f)
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
    C4coeff();
  }

  const GeodesicExact GeodesicExact::WGS84(Constants::WGS84_a<real>(),
                                           Constants::WGS84_f<real>());

  Math::real GeodesicExact::SinCosSeries(real sinx, real cosx,
                                         const real c[], int n) throw() {
    // Evaluate
    // y = sum(c[i] * cos((2*i+1) * x), i, 0, n-1) :
    // using Clenshaw summation.
    // Approx operation count = (n + 5) mult and (2 * n + 2) add
    c += n ;                    // Point to one beyond last element
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
    return cosx * (y0 - y1);    // cos(x) * (y0 - y1)
  }

  GeodesicLineExact GeodesicExact::Line(real lat1, real lon1, real azi1,
                                        unsigned caps) const throw() {
    return GeodesicLineExact(*this, lat1, lon1, azi1, caps);
  }

  Math::real GeodesicExact::GenDirect(real lat1, real lon1, real azi1,
                                      bool arcmode, real s12_a12,
                                      unsigned outmask,
                                      real& lat2, real& lon2, real& azi2,
                                      real& s12, real& m12,
                                      real& M12, real& M21,
                                      real& S12) const throw() {
    return GeodesicLineExact(*this, lat1, lon1, azi1,
                        // Automatically supply DISTANCE_IN if necessary
                        outmask | (arcmode ? NONE : DISTANCE_IN))
      .                         // Note the dot!
      GenPosition(arcmode, s12_a12, outmask,
                  lat2, lon2, azi2, s12, m12, M12, M21, S12);
  }

  Math::real GeodesicExact::GenInverse(real lat1, real lon1,
                                       real lat2, real lon2,
                                       unsigned outmask,
                                       real& s12, real& azi1, real& azi2,
                                       real& m12, real& M12, real& M21,
                                       real& S12) const throw() {
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
    // Initialize for the meridian.  No longitude calculation is done in this
    // case to let the parameter default to 0.
    EllipticFunction E(-_ep2);

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
      dn1 = (_f >= 0 ? sqrt(1 + _ep2 * Math::sq(sbet1)) :
             sqrt(1 - _e2 * Math::sq(cbet1)) / _f1),
      dn2 = (_f >= 0 ? sqrt(1 + _ep2 * Math::sq(sbet2)) :
             sqrt(1 - _e2 * Math::sq(cbet2)) / _f1);

    real
      lam12 = lon12 * Math::degree<real>(),
      slam12 = lon12 == 180 ? 0 : sin(lam12),
      clam12 = cos(lam12);      // lon12 == 90 isn't interesting

    real a12, sig12, calp1, salp1, calp2, salp2;

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
        Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                cbet1, cbet2, s12x, m12x, dummy,
                (outmask & GEODESICSCALE) != 0U, M12, M21);
      }
      // Add the check for sig12 since zero length geodesics might yield m12 <
      // 0.  Test case was
      //
      //    echo 20.001 0 20.001 0 | Geod -i
      //
      // In fact, we will have sig12 > pi/2 for meridional geodesic which is
      // not a shortest path.
      if (sig12 < 1 || m12x >= 0) {
        m12x *= _b;
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
      sig12 = InverseStart(E, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                           lam12,
                           salp1, calp1, salp2, calp2);

      if (sig12 >= 0) {
        // Short lines (InverseStart sets salp2, calp2)
        real dnm = (dn1 + dn2) / 2;
        s12x = sig12 * _b * dnm;
        m12x = Math::sq(dnm) * _b * sin(sig12 / dnm);
        if (outmask & GEODESICSCALE)
          M12 = M21 = cos(sig12 / dnm);
        a12 = sig12 / Math::degree<real>();
        omg12 = lam12 / (_f1 * dnm);
      } else {

        // Newton's method.  This is a straightforward solution of f(alp1) =
        // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
        // root in the interval (0, pi) and its derivative is positive at the
        // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
        // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
        // maintained which brackets the root and with each evaluation of
        // f(alp) the range is shrunk if possible.  Newton's method is
        // restarted whenever the derivative of f is negative (because the new
        // value of alp1 is guaranteed to be further from the solution) or if
        // the new estimate of alp1 lies outside (0,pi); in this case, the new
        // starting guess is taken to be (alp1a + alp1b) / 2.
        real ssig1, csig1, ssig2, csig2;
        real ov = 0;
        unsigned numit = 0;
        // Bracketing range
        real salp1a = tiny_, calp1a = 1, salp1b = tiny_, calp1b = -1;
        for (unsigned trip = 0; numit < maxit_; ++numit) {
          // For the WGS84 test set: mean = 1.62, sd = 1.13, max = 16
          real dv;
          real v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                            salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                            E, omg12, trip < 1, dv) - lam12;
          // Update bracketing values
          if (v >= 0 && calp1/salp1 > calp1b/salp1b) {
            salp1b = salp1; calp1b = calp1;
          } else if (v <= 0 && calp1/salp1 < calp1a/salp1a) {
            salp1a = salp1; calp1a = calp1;
          }
          if (!(abs(v) > tiny_) || !(trip < 1)) {
            if (!(abs(v) <= max(tol1_, ov)))
              numit = maxit_;
            break;
          }
          if (dv >= 0) {
            real
              dalp1 = -v/dv;
            real
              sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
              nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
            if (nsalp1 > 0 && abs(dalp1) < Math::pi<real>()) {
              calp1 = calp1 * cdalp1 - salp1 * sdalp1;
              salp1 = nsalp1;
              SinCosNorm(salp1, calp1);
              // In some regimes we don't get quadratic convergence because
              // slope -> 0.  So use convergence conditions based on epsilon
              // instead of sqrt(epsilon).  The first criterion is a test on
              // abs(v) against 200 * epsilon.  The second takes credit for an
              // anticipated reduction in abs(v) by v/ov (due to the latest
              // update in alp1) and checks this against epsilon.
              if (!(abs(v) >= tol1_ && Math::sq(v) >= ov * tol0_)) ++trip;
              ov = abs(v);
              continue;
            }
          }
          // Either dv was not postive or updated value was outside legal
          // range.  Use the midpoint of the bracket as the next estimate.
          // This mechanism is not needed for the WGS84 ellipsoid, but it does
          // catch problems with more eccentric ellipsoids.  Its efficacy is
          // such for the WGS84 test set with the starting guess set to alp1 =
          // 90deg: mean = 4.86, sd = 3.42, max = 22
          salp1 = (salp1a + salp1b)/2;
          calp1 = (calp1a + calp1b)/2;
          SinCosNorm(salp1, calp1);
          trip = 0;
          ov = 0;
        }
        if (numit >= maxit_) {
          // Resort to the safer bisection method
          for (unsigned i = 0; i < bisection_; ++i) {
            ++numit;
            salp1 = (salp1a + salp1b)/2;
            calp1 = (calp1a + calp1b)/2;
            SinCosNorm(salp1, calp1);
            if ( (abs(salp1 - salp1b) < tol0_ && calp1 - calp1b < tol0_) ||
                 (abs(salp1a - salp1) < tol0_ && calp1a - calp1 < tol0_) )
              break;
            real
              dummy,
              v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                           salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                           E, omg12, false, dummy) - lam12;
            // Be more tolerant on error.  It is approximately 1 ulp for a
            // number in [0, pi].
            if (abs(v) <= 2 * tol0_) break;
            if (v > 0) {
              salp1b = salp1; calp1b = calp1;
            } else {
              salp1a = salp1; calp1a = calp1;
            }
          }
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
          Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                  cbet1, cbet2, s12x, m12x, dummy,
                  (outmask & GEODESICSCALE) != 0U, M12, M21);
        }
        m12x *= _b;
        s12x *= _b;
        a12 = sig12 / Math::degree<real>();
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
          B41 = SinCosSeries(ssig1, csig1, C4a, nC4_),
          B42 = SinCosSeries(ssig2, csig2, C4a, nC4_);
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

  void GeodesicExact::Lengths(const EllipticFunction& E,
                              real sig12,
                              real ssig1, real csig1, real dn1,
                              real ssig2, real csig2, real dn2,
                              real cbet1, real cbet2,
                              real& s12b, real& m12b, real& m0,
                              bool scalep, real& M12, real& M21) const throw() {
    // Return m12b = (reduced length)/_b; also calculate s12b = distance/_b,
    // and m0 = coefficient of secular term in expression for reduced length.

    // It's OK to have repeated dummy arguments,
    // e.g., s12b = m0 = M12 = M21 = dummy
    m0 = - E.k2() * E.D() / (Math::pi<real>() / 2);
    real J12 = m0 *
      (sig12 + E.deltaD(ssig2, csig2, dn2) - E.deltaD(ssig1, csig1, dn1));
    // Missing a factor of _a.
    // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure accurate
    // cancellation in the case of coincident points.
    m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;
    // Missing a factor of _b
    s12b = E.E() / (Math::pi<real>() / 2) *
      (sig12 + E.deltaE(ssig2, csig2, dn2) - E.deltaE(ssig1, csig1, dn1));
    if (scalep) {
      real csig12 = csig1 * csig2 + ssig1 * ssig2;
      real t = _ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
      M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
      M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
    }
  }

  Math::real GeodesicExact::Astroid(real x, real y) throw() {
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

  Math::real GeodesicExact::InverseStart(EllipticFunction& E,
                                         real sbet1, real cbet1, real dn1,
                                         real sbet2, real cbet2, real dn2,
                                         real lam12,
                                         real& salp1, real& calp1,
                                         // Only updated if return val >= 0
                                         real& salp2, real& calp2)
    const throw() {
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
      omg12 = (!shortline ? lam12 : lam12 / (_f1 * (dn1 + dn2) / 2)),
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
    } else if (abs(_n) > real(0.1) || // Skip astroid calc if too eccentric
               csig12 >= 0 ||
               ssig12 >= 6 * abs(_n) * Math::pi<real>() * Math::sq(cbet1)) {
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
          real k2 = Math::sq(sbet1) * _ep2;
          E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
          lamscale = _e2/_f1 * cbet1 * 2 * E.H();
        }
        betscale = lamscale * cbet1;

        x = (lam12 - Math::pi<real>()) / lamscale;
        y = sbet12a / betscale;
      } else {                  // _f < 0
        // x = dlat, y = dlong
        real
          cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
          bet12a = atan2(sbet12a, cbet12a);
        real m12b, m0, dummy;
        // In the case of lon12 = 180, this repeats a calculation made in
        // Inverse.
        Lengths(E, Math::pi<real>() + bet12a,
                sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                cbet1, cbet2, dummy, m12b, m0, false,
                dummy, dummy);
        x = -1 + m12b/(cbet1 * cbet2 * m0 * Math::pi<real>());
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
    if (salp1 > 0)              // Sanity check on starting guess
      SinCosNorm(salp1, calp1);
    else {
      salp1 = 1; calp1 = 0;
    }
    return sig12;
  }

  Math::real GeodesicExact::Lambda12(real sbet1, real cbet1, real dn1,
                                     real sbet2, real cbet2, real dn2,
                                     real salp1, real calp1,
                                     real& salp2, real& calp2,
                                     real& sig12,
                                     real& ssig1, real& csig1,
                                     real& ssig2, real& csig2,
                                     EllipticFunction& E,
                                     real& omg12,
                                     bool diffp, real& dlam12) const
    throw() {

    if (sbet1 == 0 && calp1 == 0)
      // Break degeneracy of equatorial line.  This case has already been
      // handled.
      calp1 = -tiny_;

    real
      // sin(alp1) * cos(bet1) = sin(alp0)
      salp0 = salp1 * cbet1,
      calp0 = Math::hypot(calp1, salp1 * sbet1); // calp0 > 0

    real somg1, comg1, somg2, comg2, cchi1, cchi2, lam12;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
    ssig1 = sbet1; somg1 = salp0 * sbet1;
    csig1 = comg1 = calp1 * cbet1;
    // Without normalization we have schi1 = somg1.
    cchi1 = _f1 * dn1 * comg1;
    SinCosNorm(ssig1, csig1);
    // SinCosNorm(somg1, comg1); -- don't need to normalize!
    // SinCosNorm(schi1, cchi1); -- don't need to normalize!

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
    // Without normalization we have schi2 = somg2.
    cchi2 = _f1 * dn2 * comg2;
    SinCosNorm(ssig2, csig2);
    // SinCosNorm(somg2, comg2); -- don't need to normalize!
    // SinCosNorm(schi2, cchi2); -- don't need to normalize!

    // sig12 = sig2 - sig1, limit to [0, pi]
    sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, real(0)),
                  csig1 * csig2 + ssig1 * ssig2);

    // omg12 = omg2 - omg1, limit to [0, pi]
    omg12 = atan2(max(comg1 * somg2 - somg1 * comg2, real(0)),
                  comg1 * comg2 + somg1 * somg2);
    real k2 = Math::sq(calp0) * _ep2;
    E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
    real chi12 = atan2(max(cchi1 * somg2 - somg1 * cchi2, real(0)),
                       cchi1 * cchi2 + somg1 * somg2);
    lam12 = chi12 -
      _e2/_f1 * salp0 * E.H() / (Math::pi<real>() / 2) *
      (sig12 + E.deltaH(ssig2, csig2, dn2) - E.deltaH(ssig1, csig1, dn1) );

    if (diffp) {
      if (calp2 == 0)
        dlam12 = - 2 * _f1 * dn1 / sbet1;
      else {
        real dummy;
        Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                cbet1, cbet2, dummy, dlam12, dummy,
                false, dummy, dummy);
        dlam12 *= _f1 / (calp2 * cbet2);
      }
    }

    return lam12;
  }

  void GeodesicExact::C4f(real k2, real c[]) const throw() {
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

  // The coefficients C4[l] in the Fourier expansion of I4
  void GeodesicExact::C4coeff() throw() {
    switch (nC4_) {
    case 30:
      _C4x[0] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(25874463385137300111360.L)-
        real(24601948792425629614080.L)*_ep2)*_ep2-
        real(27260595352198226903040.L))+real(28775072871764795064320.L))-
        real(30435173229751225548800.L))+real(32261283623536299081728.L))-
        real(34277613850007317774336.L))+real(36513110405442577629184.L))-
        real(39002640660359117012992.L))+real(41788543564670482513920.L))-
        real(44922684332020768702464.L))+real(48469212042443460968448.L))-
        real(52508313045980416049152.L))+real(57141399491213982171136.L))-
        real(62498405693515292999680.L))+real(68748246262866822299648.L))-
        real(76114129791031124688896.L))+real(84896529382303946768384.L))-
        real(95508595555091940114432.L))+real(108532494948968113766400.L))-
        real(124812369191313330831360.L))+real(145614430723198885969920.L))-
        real(172917136483798677089280.L))+real(209970808587469822179840.L))-
        real(262463510734337277724800.L))+real(341202563954638461042240.L))-
        real(469153525437627883933080.L))+real(703730288156441825899620.L))-
        real(1231528004273773195324335.L))+real(12315280042737731953243350.L))/
        real(18472920064106597929865025.L);
      _C4x[1] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(25874463385137300111360.L)-
        real(24601948792425629614080.L)*_ep2)*_ep2-
        real(27260595352198226903040.L))+real(28775072871764795064320.L))-
        real(30435173229751225548800.L))+real(32261283623536299081728.L))-
        real(34277613850007317774336.L))+real(36513110405442577629184.L))-
        real(39002640660359117012992.L))+real(41788543564670482513920.L))-
        real(44922684332020768702464.L))+real(48469212042443460968448.L))-
        real(52508313045980416049152.L))+real(57141399491213982171136.L))-
        real(62498405693515292999680.L))+real(68748246262866822299648.L))-
        real(76114129791031124688896.L))+real(84896529382303946768384.L))-
        real(95508595555091940114432.L))+real(108532494948968113766400.L))-
        real(124812369191313330831360.L))+real(145614430723198885969920.L))-
        real(172917136483798677089280.L))+real(209970808587469822179840.L))-
        real(262463510734337277724800.L))+real(341202563954638461042240.L))-
        real(469153525437627883933080.L))+real(703730288156441825899620.L))-
        real(1231528004273773195324335.L))/real(24630560085475463906486700.L);
      _C4x[2] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(6468615846284325027840.L)-
        real(6150487198106407403520.L)*_ep2)*_ep2-
        real(6815148838049556725760.L))+real(7193768217941198766080.L))-
        real(7608793307437806387200.L))+real(8065320905884074770432.L))-
        real(8569403462501829443584.L))+real(9128277601360644407296.L))-
        real(9750660165089779253248.L))+real(10447135891167620628480.L))-
        real(11230671083005192175616.L))+real(12117303010610865242112.L))-
        real(13127078261495104012288.L))+real(14285349872803495542784.L))-
        real(15624601423378823249920.L))+real(17187061565716705574912.L))-
        real(19028532447757781172224.L))+real(21224132345575986692096.L))-
        real(23877148888772985028608.L))+real(27133123737242028441600.L))-
        real(31203092297828332707840.L))+real(36403607680799721492480.L))-
        real(43229284120949669272320.L))+real(52492702146867455544960.L))-
        real(65615877683584319431200.L))+real(85300640988659615260560.L))-
        real(117288381359406970983270.L))+real(175932572039110456474905.L))/
        real(7389168025642639171946010.L);
      _C4x[3] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(3234307923142162513920.L)-
        real(3075243599053203701760.L)*_ep2)*_ep2-
        real(3407574419024778362880.L))+real(3596884108970599383040.L))-
        real(3804396653718903193600.L))+real(4032660452942037385216.L))-
        real(4284701731250914721792.L))+real(4564138800680322203648.L))-
        real(4875330082544889626624.L))+real(5223567945583810314240.L))-
        real(5615335541502596087808.L))+real(6058651505305432621056.L))-
        real(6563539130747552006144.L))+real(7142674936401747771392.L))-
        real(7812300711689411624960.L))+real(8593530782858352787456.L))-
        real(9514266223878890586112.L))+real(10612066172787993346048.L))-
        real(11938574444386492514304.L))+real(13566561868621014220800.L))-
        real(15601546148914166353920.L))+real(18201803840399860746240.L))-
        real(21614642060474834636160.L))+real(26246351073433727772480.L))-
        real(32807938841792159715600.L))+real(42650320494329807630280.L))-
        real(58644190679703485491635.L))/real(4222381728938650955397720.L);
      _C4x[4] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(404288490392770314240.L)-real(384405449881650462720.L)*
        _ep2)*_ep2-real(425946802378097295360.L))+
        real(449610513621324922880.L))-real(475549581714862899200.L))+
        real(504082556617754673152.L))-real(535587716406364340224.L))+
        real(570517350085040275456.L))-real(609416260318111203328.L))+
        real(652945993197976289280.L))-real(701916942687824510976.L))+
        real(757331438163179077632.L))-real(820442391343444000768.L))+
        real(892834367050218471424.L))-real(976537588961176453120.L))+
        real(1074191347857294098432.L))-real(1189283277984861323264.L))+
        real(1326508271598499168256.L))-real(1492321805548311564288.L))+
        real(1695820233577626777600.L))-real(1950193268614270794240.L))+
        real(2275225480049982593280.L))-real(2701830257559354329520.L))+
        real(3280793884179215971560.L))-real(4100992355224019964450.L))+
        real(5331290061791225953785.L))/real(586441906797034854916350.L);
      _C4x[5] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(202144245196385157120.L)-real(192202724940825231360.L)*_ep2)*
        _ep2-real(212973401189048647680.L))+real(224805256810662461440.L))-
        real(237774790857431449600.L))+real(252041278308877336576.L))-
        real(267793858203182170112.L))+real(285258675042520137728.L))-
        real(304708130159055601664.L))+real(326472996598988144640.L))-
        real(350958471343912255488.L))+real(378665719081589538816.L))-
        real(410221195671722000384.L))+real(446417183525109235712.L))-
        real(488268794480588226560.L))+real(537095673928647049216.L))-
        real(594641638992430661632.L))+real(663254135799249584128.L))-
        real(746160902774155782144.L))+real(847910116788813388800.L))-
        real(975096634307135397120.L))+real(1137612740024991296640.L))-
        real(1350915128779677164760.L))+real(1640396942089607985780.L))-
        real(2050496177612009982225.L))/real(319877403707473557227100.L);
      _C4x[6] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(50536061299096289280.L)-real(48050681235206307840.L)*_ep2)*_ep2-
        real(53243350297262161920.L))+real(56201314202665615360.L))-
        real(59443697714357862400.L))+real(63010319577219334144.L))-
        real(66948464550795542528.L))+real(71314668760630034432.L))-
        real(76177032539763900416.L))+real(81618249149747036160.L))-
        real(87739617835978063872.L))+real(94666429770397384704.L))-
        real(102555298917930500096.L))+real(111604295881277308928.L))-
        real(122067198620147056640.L))+real(134273918482161762304.L))-
        real(148660409748107665408.L))+real(165813533949812396032.L))-
        real(186540225693538945536.L))+real(211977529197203347200.L))-
        real(243774158576783849280.L))+real(284403185006247824160.L))-
        real(337728782194919291190.L))+real(410099235522401996445.L))/
        real(86120839459704419253450.L);
      _C4x[7] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(25268030649548144640.L)-real(24025340617603153920.L)*_ep2)*_ep2-
        real(26621675148631080960.L))+real(28100657101332807680.L))-
        real(29721848857178931200.L))+real(31505159788609667072.L))-
        real(33474232275397771264.L))+real(35657334380315017216.L))-
        real(38088516269881950208.L))+real(40809124574873518080.L))-
        real(43869808917989031936.L))+real(47333214885198692352.L))-
        real(51277649458965250048.L))+real(55802147940638654464.L))-
        real(61033599310073528320.L))+real(67136959241080881152.L))-
        real(74330204874053832704.L))+real(82906766974906198016.L))-
        real(93270112846769472768.L))+real(105988764598601673600.L))-
        real(121887079288391924640.L))+real(142201592503123912080.L))-
        real(168864391097459645595.L))/real(45931114378509023601840.L);
      _C4x[8] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1579251915596759040.L)-real(1501583788600197120.L)*_ep2)*_ep2-
        real(1663854696789442560.L))+real(1756291068833300480.L))-
        real(1857615553573683200.L))+real(1969072486788104192.L))-
        real(2092139517212360704.L))+real(2228583398769688576.L))-
        real(2380532266867621888.L))+real(2550570285929594880.L))-
        real(2741863057374314496.L))+real(2958325930324918272.L))-
        real(3204853091185328128.L))+real(3487634246289915904.L))-
        real(3814599956879595520.L))+real(4196059952567555072.L))-
        real(4645637804628364544.L))+real(5181672935931637376.L))-
        real(5829382052923092048.L))+real(6624297787412604600.L))-
        real(7617942455524495290.L))+real(8887599531445244505.L))/
        real(3039559039754273620710.L);
      _C4x[9] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(789625957798379520.L)-real(750791894300098560.L)*_ep2)*_ep2-
        real(831927348394721280.L))+real(878145534416650240.L))-
        real(928807776786841600.L))+real(984536243394052096.L))-
        real(1046069758606180352.L))+real(1114291699384844288.L))-
        real(1190266133433810944.L))+real(1275285142964797440.L))-
        real(1370931528687157248.L))+real(1479162965162459136.L))-
        real(1602426545592664064.L))+real(1743817123144957952.L))-
        real(1907299978439797760.L))+real(2098029976283777536.L))-
        real(2322818902314182272.L))+real(2590836467965818688.L))-
        real(2914691026461546024.L))+real(3312148893706302300.L))-
        real(3808971227762247645.L))/real(1599767915660144010900.L);
      _C4x[10] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(197406489449594880.L)-
        real(187697973575024640.L)*_ep2)*_ep2-real(207981837098680320.L))+
        real(219536383604162560.L))-real(232201944196710400.L))+
        real(246134060848513024.L))-real(261517439651545088.L))+
        real(278572924846211072.L))-real(297566533358452736.L))+
        real(318821285741199360.L))-real(342732882171789312.L))+
        real(369790741290614784.L))-real(400606636398166016.L))+
        real(435954280786239488.L))-real(476824994609949440.L))+
        real(524507494070944384.L))-real(580704725578545568.L))+
        real(647709116991454672.L))-real(728672756615386506.L))+
        real(828037223426575575.L))/real(418986835053847240950.L);
      _C4x[11] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(98703244724797440.L)-
        real(93848986787512320.L)*_ep2)*_ep2-real(103990918549340160.L))+
        real(109768191802081280.L))-real(116100972098355200.L))+
        real(123067030424256512.L))-real(130758719825772544.L))+
        real(139286462423105536.L))-real(148783266679226368.L))+
        real(159410642870599680.L))-real(171366441085894656.L))+
        real(184895370645307392.L))-real(200303318199083008.L))+
        real(217977140393119744.L))-real(238412497304974720.L))+
        real(262253747035472192.L))-real(290352362789272784.L))+
        real(323854558495727336.L))-real(364336378307693253.L))/
        real(218601826984615951800.L);
      _C4x[12] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(12337905590599680.L)-
        real(11731123348439040.L)*_ep2)*_ep2-real(12998864818667520.L))+
        real(13721023975260160.L))-real(14512621512294400.L))+
        real(15383378803032064.L))-real(16344839978221568.L))+
        real(17410807802888192.L))-real(18597908334903296.L))+
        real(19926330358824960.L))-real(21420805135736832.L))+
        real(23111921330663424.L))-real(25037914774885376.L))+
        real(27247142549139968.L))-real(29801562163121840.L))+
        real(32781718379434024.L))-real(36294045348659098.L))+
        real(40481819811965917.L))/real(28418237508000073734.L);
      _C4x[13] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(6168952795299840.L)-
        real(5865561674219520.L)*_ep2)*_ep2-real(6499432409333760.L))+
        real(6860511987630080.L))-real(7256310756147200.L))+
        real(7691689401516032.L))-real(8172419989110784.L))+
        real(8705403901444096.L))-real(9298954167451648.L))+
        real(9963165179412480.L))-real(10710402567868416.L))+
        real(11555960665331712.L))-real(12518957387442688.L))+
        real(13623571274569984.L))-real(14900781081560920.L))+
        real(16390859189717012.L))-real(18147022674329549.L))/
        real(14735382411555593788.L);
      _C4x[14] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(3810874761609205290.L);
      _C4x[15] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(1966903102766041440.L);
      _C4x[16] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(24097471856640.L)-real(22912350289920.L)*_ep2)*_ep2-
        real(25388407848960.L))+real(26798874951680.L))-real(28344963891200.L))+
        real(30045661724672.L))-real(31923515582464.L))+real(34005483990016.L))-
        real(36324039716608.L))+real(38918613982080.L))-real(41837510030736.L))+
        real(45140471348952.L))-real(48902177294698.L))+real(53217075291289.L))/
        real(63328319596633910.L);
      _C4x[17] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(32568850078268868.L);
      _C4x[18] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(8362272317393358.L);
      _C4x[19] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(4288344778150440.L);
      _C4x[20] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/real(549117319153410.L);
      _C4x[21] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(280943744683140.L);
      _C4x[22] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(71796734752358.L);
      _C4x[23] = (_ep2*(_ep2*(_ep2*(_ep2*((real(511580160.L)-real(486420480.L)*
        _ep2)*_ep2-real(538986240.L))+real(568929920.L))-real(601752800.L))+
        real(637857968.L))-real(677724091.L))/real(1594007062032.L);
      _C4x[24] = (_ep2*(_ep2*(_ep2*((real(31973760.L)-real(30401280.L)*_ep2)*
        _ep2-real(33686640.L))+real(35558120.L))-real(37609550.L))+
        real(39866123.L))/real(101658613650.L);
      _C4x[25] = (_ep2*(_ep2*((real(3197376.L)-real(3040128.L)*_ep2)*_ep2-
        real(3368664.L))+real(3555812.L))-real(3760955.L))/real(10365191980.L);
      _C4x[26] = (_ep2*((real(61488.L)-real(58464.L)*_ep2)*_ep2-real(64782.L))+
        real(68381.L))/real(203091570.L);
      _C4x[27] = ((real(3416.L)-real(3248.L)*_ep2)*_ep2-real(3599.L))/
        real(11488008.L);
      _C4x[28] = (real(61.L)-real(58.L)*_ep2)/real(208742.L);
      _C4x[29] = -1/real(3660.L);
      _C4x[30] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(24601948792425629614080.L)*_ep2-
        real(25874463385137300111360.L))+real(27260595352198226903040.L))-
        real(28775072871764795064320.L))+real(30435173229751225548800.L))-
        real(32261283623536299081728.L))+real(34277613850007317774336.L))-
        real(36513110405442577629184.L))+real(39002640660359117012992.L))-
        real(41788543564670482513920.L))+real(44922684332020768702464.L))-
        real(48469212042443460968448.L))+real(52508313045980416049152.L))-
        real(57141399491213982171136.L))+real(62498405693515292999680.L))-
        real(68748246262866822299648.L))+real(76114129791031124688896.L))-
        real(84896529382303946768384.L))+real(95508595555091940114432.L))-
        real(108532494948968113766400.L))+real(124812369191313330831360.L))-
        real(145614430723198885969920.L))+real(172917136483798677089280.L))-
        real(209970808587469822179840.L))+real(262463510734337277724800.L))-
        real(341202563954638461042240.L))+real(469153525437627883933080.L))-
        real(703730288156441825899620.L))+real(1231528004273773195324335.L))/
        real(221675040769279175158380300.L);
      _C4x[31] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(6150487198106407403520.L)*_ep2-
        real(6468615846284325027840.L))+real(6815148838049556725760.L))-
        real(7193768217941198766080.L))+real(7608793307437806387200.L))-
        real(8065320905884074770432.L))+real(8569403462501829443584.L))-
        real(9128277601360644407296.L))+real(9750660165089779253248.L))-
        real(10447135891167620628480.L))+real(11230671083005192175616.L))-
        real(12117303010610865242112.L))+real(13127078261495104012288.L))-
        real(14285349872803495542784.L))+real(15624601423378823249920.L))-
        real(17187061565716705574912.L))+real(19028532447757781172224.L))-
        real(21224132345575986692096.L))+real(23877148888772985028608.L))-
        real(27133123737242028441600.L))+real(31203092297828332707840.L))-
        real(36403607680799721492480.L))+real(43229284120949669272320.L))-
        real(52492702146867455544960.L))+real(65615877683584319431200.L))-
        real(85300640988659615260560.L))+real(117288381359406970983270.L))-
        real(175932572039110456474905.L))/real(44335008153855835031676060.L);
      _C4x[32] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(3075243599053203701760.L)*_ep2-
        real(3234307923142162513920.L))+real(3407574419024778362880.L))-
        real(3596884108970599383040.L))+real(3804396653718903193600.L))-
        real(4032660452942037385216.L))+real(4284701731250914721792.L))-
        real(4564138800680322203648.L))+real(4875330082544889626624.L))-
        real(5223567945583810314240.L))+real(5615335541502596087808.L))-
        real(6058651505305432621056.L))+real(6563539130747552006144.L))-
        real(7142674936401747771392.L))+real(7812300711689411624960.L))-
        real(8593530782858352787456.L))+real(9514266223878890586112.L))-
        real(10612066172787993346048.L))+real(11938574444386492514304.L))-
        real(13566561868621014220800.L))+real(15601546148914166353920.L))-
        real(18201803840399860746240.L))+real(21614642060474834636160.L))-
        real(26246351073433727772480.L))+real(32807938841792159715600.L))-
        real(42650320494329807630280.L))+real(58644190679703485491635.L))/
        real(21111908644693254776988600.L);
      _C4x[33] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(384405449881650462720.L)*_ep2-
        real(404288490392770314240.L))+real(425946802378097295360.L))-
        real(449610513621324922880.L))+real(475549581714862899200.L))-
        real(504082556617754673152.L))+real(535587716406364340224.L))-
        real(570517350085040275456.L))+real(609416260318111203328.L))-
        real(652945993197976289280.L))+real(701916942687824510976.L))-
        real(757331438163179077632.L))+real(820442391343444000768.L))-
        real(892834367050218471424.L))+real(976537588961176453120.L))-
        real(1074191347857294098432.L))+real(1189283277984861323264.L))-
        real(1326508271598499168256.L))+real(1492321805548311564288.L))-
        real(1695820233577626777600.L))+real(1950193268614270794240.L))-
        real(2275225480049982593280.L))+real(2701830257559354329520.L))-
        real(3280793884179215971560.L))+real(4100992355224019964450.L))-
        real(5331290061791225953785.L))/real(2638988580586656847123575.L);
      _C4x[34] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(192202724940825231360.L)*_ep2-
        real(202144245196385157120.L))+real(212973401189048647680.L))-
        real(224805256810662461440.L))+real(237774790857431449600.L))-
        real(252041278308877336576.L))+real(267793858203182170112.L))-
        real(285258675042520137728.L))+real(304708130159055601664.L))-
        real(326472996598988144640.L))+real(350958471343912255488.L))-
        real(378665719081589538816.L))+real(410221195671722000384.L))-
        real(446417183525109235712.L))+real(488268794480588226560.L))-
        real(537095673928647049216.L))+real(594641638992430661632.L))-
        real(663254135799249584128.L))+real(746160902774155782144.L))-
        real(847910116788813388800.L))+real(975096634307135397120.L))-
        real(1137612740024991296640.L))+real(1350915128779677164760.L))-
        real(1640396942089607985780.L))+real(2050496177612009982225.L))/
        real(1343485095571388940353820.L);
      _C4x[35] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(48050681235206307840.L)*_ep2-real(50536061299096289280.L))+
        real(53243350297262161920.L))-real(56201314202665615360.L))+
        real(59443697714357862400.L))-real(63010319577219334144.L))+
        real(66948464550795542528.L))-real(71314668760630034432.L))+
        real(76177032539763900416.L))-real(81618249149747036160.L))+
        real(87739617835978063872.L))-real(94666429770397384704.L))+
        real(102555298917930500096.L))-real(111604295881277308928.L))+
        real(122067198620147056640.L))-real(134273918482161762304.L))+
        real(148660409748107665408.L))-real(165813533949812396032.L))+
        real(186540225693538945536.L))-real(211977529197203347200.L))+
        real(243774158576783849280.L))-real(284403185006247824160.L))+
        real(337728782194919291190.L))-real(410099235522401996445.L))/
        real(344483357838817677013800.L);
      _C4x[36] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(24025340617603153920.L)*_ep2-real(25268030649548144640.L))+
        real(26621675148631080960.L))-real(28100657101332807680.L))+
        real(29721848857178931200.L))-real(31505159788609667072.L))+
        real(33474232275397771264.L))-real(35657334380315017216.L))+
        real(38088516269881950208.L))-real(40809124574873518080.L))+
        real(43869808917989031936.L))-real(47333214885198692352.L))+
        real(51277649458965250048.L))-real(55802147940638654464.L))+
        real(61033599310073528320.L))-real(67136959241080881152.L))+
        real(74330204874053832704.L))-real(82906766974906198016.L))+
        real(93270112846769472768.L))-real(105988764598601673600.L))+
        real(121887079288391924640.L))-real(142201592503123912080.L))+
        real(168864391097459645595.L))/real(177162869745677662464240.L);
      _C4x[37] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(3003167577200394240.L)*_ep2-real(3158503831193518080.L))+
        real(3327709393578885120.L))-real(3512582137666600960.L))+
        real(3715231107147366400.L))-real(3938144973576208384.L))+
        real(4184279034424721408.L))-real(4457166797539377152.L))+
        real(4761064533735243776.L))-real(5101140571859189760.L))+
        real(5483726114748628992.L))-real(5916651860649836544.L))+
        real(6409706182370656256.L))-real(6975268492579831808.L))+
        real(7629199913759191040.L))-real(8392119905135110144.L))+
        real(9291275609256729088.L))-real(10363345871863274752.L))+
        real(11658764105846184096.L))-real(13248595574825209200.L))+
        real(15235884911048990580.L))-real(17775199062890489010.L))/
        real(22796692798157052155325.L);
      _C4x[38] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(750791894300098560.L)*_ep2-real(789625957798379520.L))+
        real(831927348394721280.L))-real(878145534416650240.L))+
        real(928807776786841600.L))-real(984536243394052096.L))+
        real(1046069758606180352.L))-real(1114291699384844288.L))+
        real(1190266133433810944.L))-real(1275285142964797440.L))+
        real(1370931528687157248.L))-real(1479162965162459136.L))+
        real(1602426545592664064.L))-real(1743817123144957952.L))+
        real(1907299978439797760.L))-real(2098029976283777536.L))+
        real(2322818902314182272.L))-real(2590836467965818688.L))+
        real(2914691026461546024.L))-real(3312148893706302300.L))+
        real(3808971227762247645.L))/real(5865815690753861373300.L);
      _C4x[39] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(187697973575024640.L)*_ep2-real(197406489449594880.L))+
        real(207981837098680320.L))-real(219536383604162560.L))+
        real(232201944196710400.L))-real(246134060848513024.L))+
        real(261517439651545088.L))-real(278572924846211072.L))+
        real(297566533358452736.L))-real(318821285741199360.L))+
        real(342732882171789312.L))-real(369790741290614784.L))+
        real(400606636398166016.L))-real(435954280786239488.L))+
        real(476824994609949440.L))-real(524507494070944384.L))+
        real(580704725578545568.L))-real(647709116991454672.L))+
        real(728672756615386506.L))-real(828037223426575575.L))/
        real(1508352606193850067420.L);
      _C4x[40] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(93848986787512320.L)*
        _ep2-real(98703244724797440.L))+real(103990918549340160.L))-
        real(109768191802081280.L))+real(116100972098355200.L))-
        real(123067030424256512.L))+real(130758719825772544.L))-
        real(139286462423105536.L))+real(148783266679226368.L))-
        real(159410642870599680.L))+real(171366441085894656.L))-
        real(184895370645307392.L))+real(200303318199083008.L))-
        real(217977140393119744.L))+real(238412497304974720.L))-
        real(262253747035472192.L))+real(290352362789272784.L))-
        real(323854558495727336.L))+real(364336378307693253.L))/
        real(775042841127274738200.L);
      _C4x[41] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11731123348439040.L)*_ep2-
        real(12337905590599680.L))+real(12998864818667520.L))-
        real(13721023975260160.L))+real(14512621512294400.L))-
        real(15383378803032064.L))+real(16344839978221568.L))-
        real(17410807802888192.L))+real(18597908334903296.L))-
        real(19926330358824960.L))+real(21420805135736832.L))-
        real(23111921330663424.L))+real(25037914774885376.L))-
        real(27247142549139968.L))+real(29801562163121840.L))-
        real(32781718379434024.L))+real(36294045348659098.L))-
        real(40481819811965917.L))/real(99463831278000258069.L);
      _C4x[42] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(5865561674219520.L)*_ep2-
        real(6168952795299840.L))+real(6499432409333760.L))-
        real(6860511987630080.L))+real(7256310756147200.L))-
        real(7691689401516032.L))+real(8172419989110784.L))-
        real(8705403901444096.L))+real(9298954167451648.L))-
        real(9963165179412480.L))+real(10710402567868416.L))-
        real(11555960665331712.L))+real(12518957387442688.L))-
        real(13623571274569984.L))+real(14900781081560920.L))-
        real(16390859189717012.L))+real(18147022674329549.L))/
        real(51007092963077055420.L);
      _C4x[43] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(1466390418554880.L)*_ep2-
        real(1542238198824960.L))+real(1624858102333440.L))-
        real(1715127996907520.L))+real(1814077689036800.L))-
        real(1922922350379008.L))+real(2043104997277696.L))-
        real(2176350975361024.L))+real(2324738541862912.L))-
        real(2490791294853120.L))+real(2677600641967104.L))-
        real(2888990166332928.L))+real(3129739346860672.L))-
        real(3405892818642496.L))+real(3725195270390230.L))-
        real(4097714797429253.L))/real(13065856325517275280.L);
      _C4x[44] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(6687470549404540896.L);
      _C4x[45] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(91649401159680.L)*_ep2-real(96389887426560.L))+
        real(101553631395840.L))-real(107195499806720.L))+
        real(113379855564800.L))-real(120182646898688.L))+
        real(127694062329856.L))-real(136021935960064.L))+
        real(145296158866432.L))-real(155674455928320.L))+
        real(167350040122944.L))-real(180561885395808.L))+
        real(195608709178792.L))-real(212868301165156.L))/
        real(854932314554557785.L);
      _C4x[46] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(109201438497725028.L);
      _C4x[47] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(2864043786240.L)*_ep2-real(3012183982080.L))+
        real(3173550981120.L))-real(3349859368960.L))+real(3543120486400.L))-
        real(3755707715584.L))+real(3990439447808.L))-real(4250685498752.L))+
        real(4540504964576.L))-real(4864826747760.L))+real(5229688753842.L))-
        real(5642558918619.L))/real(27874241057977860.L);
      _C4x[48] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(14219248474919880.L);
      _C4x[49] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(1812087153206253.L);
      _C4x[50] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(89501368320.L)*
        _ep2-real(94130749440.L))+real(99173468160.L))-real(104683105280.L))+
        real(110722515200.L))-real(117365866112.L))+real(124701232744.L))-
        real(132833921836.L))+real(141890780143.L))/real(923100875387460.L);
      _C4x[51] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(234971131916808.L);
      _C4x[52] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(119550529652400.L);
      _C4x[53] = (_ep2*(_ep2*(_ep2*(_ep2*(real(60802560.L)*_ep2-
        real(63947520.L))+real(67373280.L))-real(71116240.L))+real(75219100.L))-
        real(79732246.L))/real(660780988725.L);
      _C4x[54] = (_ep2*(_ep2*(_ep2*(real(15200640.L)*_ep2-real(15986880.L))+
        real(16843320.L))-real(17779060.L))+real(18804775.L))/
        real(167916110076.L);
      _C4x[55] = (_ep2*(_ep2*(real(760032.L)*_ep2-real(799344.L))+
        real(842166.L))-real(888953.L))/real(8529845940.L);
      _C4x[56] = (_ep2*(real(9744.L)*_ep2-real(10248.L))+real(10797.L))/
        real(111050744.L);
      _C4x[57] = (real(406.L)*_ep2-real(427.L))/real(4696695.L);
      _C4x[58] = real(29.L)/real(340380.L);
      _C4x[59] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(6468615846284325027840.L)-
        real(6150487198106407403520.L)*_ep2)*_ep2-
        real(6815148838049556725760.L))+real(7193768217941198766080.L))-
        real(7608793307437806387200.L))+real(8065320905884074770432.L))-
        real(8569403462501829443584.L))+real(9128277601360644407296.L))-
        real(9750660165089779253248.L))+real(10447135891167620628480.L))-
        real(11230671083005192175616.L))+real(12117303010610865242112.L))-
        real(13127078261495104012288.L))+real(14285349872803495542784.L))-
        real(15624601423378823249920.L))+real(17187061565716705574912.L))-
        real(19028532447757781172224.L))+real(21224132345575986692096.L))-
        real(23877148888772985028608.L))+real(27133123737242028441600.L))-
        real(31203092297828332707840.L))+real(36403607680799721492480.L))-
        real(43229284120949669272320.L))+real(52492702146867455544960.L))-
        real(65615877683584319431200.L))+real(85300640988659615260560.L))-
        real(117288381359406970983270.L))+real(175932572039110456474905.L))/
        real(369458401282131958597300500.L);
      _C4x[60] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(3234307923142162513920.L)-
        real(3075243599053203701760.L)*_ep2)*_ep2-
        real(3407574419024778362880.L))+real(3596884108970599383040.L))-
        real(3804396653718903193600.L))+real(4032660452942037385216.L))-
        real(4284701731250914721792.L))+real(4564138800680322203648.L))-
        real(4875330082544889626624.L))+real(5223567945583810314240.L))-
        real(5615335541502596087808.L))+real(6058651505305432621056.L))-
        real(6563539130747552006144.L))+real(7142674936401747771392.L))-
        real(7812300711689411624960.L))+real(8593530782858352787456.L))-
        real(9514266223878890586112.L))+real(10612066172787993346048.L))-
        real(11938574444386492514304.L))+real(13566561868621014220800.L))-
        real(15601546148914166353920.L))+real(18201803840399860746240.L))-
        real(21614642060474834636160.L))+real(26246351073433727772480.L))-
        real(32807938841792159715600.L))+real(42650320494329807630280.L))-
        real(58644190679703485491635.L))/real(105559543223466273884943000.L);
      _C4x[61] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(404288490392770314240.L)-real(384405449881650462720.L)*
        _ep2)*_ep2-real(425946802378097295360.L))+
        real(449610513621324922880.L))-real(475549581714862899200.L))+
        real(504082556617754673152.L))-real(535587716406364340224.L))+
        real(570517350085040275456.L))-real(609416260318111203328.L))+
        real(652945993197976289280.L))-real(701916942687824510976.L))+
        real(757331438163179077632.L))-real(820442391343444000768.L))+
        real(892834367050218471424.L))-real(976537588961176453120.L))+
        real(1074191347857294098432.L))-real(1189283277984861323264.L))+
        real(1326508271598499168256.L))-real(1492321805548311564288.L))+
        real(1695820233577626777600.L))-real(1950193268614270794240.L))+
        real(2275225480049982593280.L))-real(2701830257559354329520.L))+
        real(3280793884179215971560.L))-real(4100992355224019964450.L))+
        real(5331290061791225953785.L))/real(10262733368948109961036125.L);
      _C4x[62] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(202144245196385157120.L)-real(192202724940825231360.L)*_ep2)*
        _ep2-real(212973401189048647680.L))+real(224805256810662461440.L))-
        real(237774790857431449600.L))+real(252041278308877336576.L))-
        real(267793858203182170112.L))+real(285258675042520137728.L))-
        real(304708130159055601664.L))+real(326472996598988144640.L))-
        real(350958471343912255488.L))+real(378665719081589538816.L))-
        real(410221195671722000384.L))+real(446417183525109235712.L))-
        real(488268794480588226560.L))+real(537095673928647049216.L))-
        real(594641638992430661632.L))+real(663254135799249584128.L))-
        real(746160902774155782144.L))+real(847910116788813388800.L))-
        real(975096634307135397120.L))+real(1137612740024991296640.L))-
        real(1350915128779677164760.L))+real(1640396942089607985780.L))-
        real(2050496177612009982225.L))/real(4478283651904629801179400.L);
      _C4x[63] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(50536061299096289280.L)-real(48050681235206307840.L)*_ep2)*_ep2-
        real(53243350297262161920.L))+real(56201314202665615360.L))-
        real(59443697714357862400.L))+real(63010319577219334144.L))-
        real(66948464550795542528.L))+real(71314668760630034432.L))-
        real(76177032539763900416.L))+real(81618249149747036160.L))-
        real(87739617835978063872.L))+real(94666429770397384704.L))-
        real(102555298917930500096.L))+real(111604295881277308928.L))-
        real(122067198620147056640.L))+real(134273918482161762304.L))-
        real(148660409748107665408.L))+real(165813533949812396032.L))-
        real(186540225693538945536.L))+real(211977529197203347200.L))-
        real(243774158576783849280.L))+real(284403185006247824160.L))-
        real(337728782194919291190.L))+real(410099235522401996445.L))/
        real(1033450073516453031041400.L);
      _C4x[64] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(25268030649548144640.L)-real(24025340617603153920.L)*_ep2)*_ep2-
        real(26621675148631080960.L))+real(28100657101332807680.L))-
        real(29721848857178931200.L))+real(31505159788609667072.L))-
        real(33474232275397771264.L))+real(35657334380315017216.L))-
        real(38088516269881950208.L))+real(40809124574873518080.L))-
        real(43869808917989031936.L))+real(47333214885198692352.L))-
        real(51277649458965250048.L))+real(55802147940638654464.L))-
        real(61033599310073528320.L))+real(67136959241080881152.L))-
        real(74330204874053832704.L))+real(82906766974906198016.L))-
        real(93270112846769472768.L))+real(105988764598601673600.L))-
        real(121887079288391924640.L))+real(142201592503123912080.L))-
        real(168864391097459645595.L))/real(492119082626882395734000.L);
      _C4x[65] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3158503831193518080.L)-real(3003167577200394240.L)*_ep2)*_ep2-
        real(3327709393578885120.L))+real(3512582137666600960.L))-
        real(3715231107147366400.L))+real(3938144973576208384.L))-
        real(4184279034424721408.L))+real(4457166797539377152.L))-
        real(4761064533735243776.L))+real(5101140571859189760.L))-
        real(5483726114748628992.L))+real(5916651860649836544.L))-
        real(6409706182370656256.L))+real(6975268492579831808.L))-
        real(7629199913759191040.L))+real(8392119905135110144.L))-
        real(9291275609256729088.L))+real(10363345871863274752.L))-
        real(11658764105846184096.L))+real(13248595574825209200.L))-
        real(15235884911048990580.L))+real(17775199062890489010.L))/
        real(59705623995173231835375.L);
      _C4x[66] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(789625957798379520.L)-real(750791894300098560.L)*_ep2)*_ep2-
        real(831927348394721280.L))+real(878145534416650240.L))-
        real(928807776786841600.L))+real(984536243394052096.L))-
        real(1046069758606180352.L))+real(1114291699384844288.L))-
        real(1190266133433810944.L))+real(1275285142964797440.L))-
        real(1370931528687157248.L))+real(1479162965162459136.L))-
        real(1602426545592664064.L))+real(1743817123144957952.L))-
        real(1907299978439797760.L))+real(2098029976283777536.L))-
        real(2322818902314182272.L))+real(2590836467965818688.L))-
        real(2914691026461546024.L))+real(3312148893706302300.L))-
        real(3808971227762247645.L))/real(14664539226884653433250.L);
      _C4x[67] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(197406489449594880.L)-
        real(187697973575024640.L)*_ep2)*_ep2-real(207981837098680320.L))+
        real(219536383604162560.L))-real(232201944196710400.L))+
        real(246134060848513024.L))-real(261517439651545088.L))+
        real(278572924846211072.L))-real(297566533358452736.L))+
        real(318821285741199360.L))-real(342732882171789312.L))+
        real(369790741290614784.L))-real(400606636398166016.L))+
        real(435954280786239488.L))-real(476824994609949440.L))+
        real(524507494070944384.L))-real(580704725578545568.L))+
        real(647709116991454672.L))-real(728672756615386506.L))+
        real(828037223426575575.L))/real(3631219237133342754900.L);
      _C4x[68] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(98703244724797440.L)-
        real(93848986787512320.L)*_ep2)*_ep2-real(103990918549340160.L))+
        real(109768191802081280.L))-real(116100972098355200.L))+
        real(123067030424256512.L))-real(130758719825772544.L))+
        real(139286462423105536.L))-real(148783266679226368.L))+
        real(159410642870599680.L))-real(171366441085894656.L))+
        real(184895370645307392.L))-real(200303318199083008.L))+
        real(217977140393119744.L))-real(238412497304974720.L))+
        real(262253747035472192.L))-real(290352362789272784.L))+
        real(323854558495727336.L))-real(364336378307693253.L))/
        real(1808433295963641055800.L);
      _C4x[69] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(12337905590599680.L)-
        real(11731123348439040.L)*_ep2)*_ep2-real(12998864818667520.L))+
        real(13721023975260160.L))-real(14512621512294400.L))+
        real(15383378803032064.L))-real(16344839978221568.L))+
        real(17410807802888192.L))-real(18597908334903296.L))+
        real(19926330358824960.L))-real(21420805135736832.L))+
        real(23111921330663424.L))-real(25037914774885376.L))+
        real(27247142549139968.L))-real(29801562163121840.L))+
        real(32781718379434024.L))-real(36294045348659098.L))+
        real(40481819811965917.L))/real(226054161995455131975.L);
      _C4x[70] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(6168952795299840.L)-
        real(5865561674219520.L)*_ep2)*_ep2-real(6499432409333760.L))+
        real(6860511987630080.L))-real(7256310756147200.L))+
        real(7691689401516032.L))-real(8172419989110784.L))+
        real(8705403901444096.L))-real(9298954167451648.L))+
        real(9963165179412480.L))-real(10710402567868416.L))+
        real(11555960665331712.L))-real(12518957387442688.L))+
        real(13623571274569984.L))-real(14900781081560920.L))+
        real(16390859189717012.L))-real(18147022674329549.L))/
        real(113349095473504567600.L);
      _C4x[71] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(28476866350486369200.L);
      _C4x[72] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(14330294034438301920.L);
      _C4x[73] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(96389887426560.L)-real(91649401159680.L)*_ep2)*_ep2-
        real(101553631395840.L))+real(107195499806720.L))-
        real(113379855564800.L))+real(120182646898688.L))-
        real(127694062329856.L))+real(136021935960064.L))-
        real(145296158866432.L))+real(155674455928320.L))-
        real(167350040122944.L))+real(180561885395808.L))-
        real(195608709178792.L))+real(212868301165156.L))/
        real(1804857108504066435.L);
      _C4x[74] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(227502996870260475.L);
      _C4x[75] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(57388143354660300.L);
      _C4x[76] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(28965135782244200.L);
      _C4x[77] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/real(3655965309100335.L);
      _C4x[78] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(1846201750774920.L);
      _C4x[79] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(466212563327000.L);
      _C4x[80] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(235478315982000.L);
      _C4x[81] = (_ep2*(_ep2*(_ep2*((real(1470792960.L)-real(1398458880.L)*
        _ep2)*_ep2-real(1549585440.L))+real(1635673520.L))-real(1730039300.L))+
        real(1833841658.L))/real(29735144492625.L);
      _C4x[82] = (_ep2*(_ep2*((real(3197376.L)-real(3040128.L)*_ep2)*_ep2-
        real(3368664.L))+real(3555812.L))-real(3760955.L))/real(65300709474.L);
      _C4x[83] = (_ep2*((real(799344.L)-real(760032.L)*_ep2)*_ep2-
        real(842166.L))+real(888953.L))/real(16491035484.L);
      _C4x[84] = ((real(133224.L)-real(126672.L)*_ep2)*_ep2-real(140361.L))/
        real(2776268600.L);
      _C4x[85] = (real(3843.L)-real(3654.L)*_ep2)/real(80887525.L);
      _C4x[86] = -real(203.L)/real(4538400.L);
      _C4x[87] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(3075243599053203701760.L)*_ep2-
        real(3234307923142162513920.L))+real(3407574419024778362880.L))-
        real(3596884108970599383040.L))+real(3804396653718903193600.L))-
        real(4032660452942037385216.L))+real(4284701731250914721792.L))-
        real(4564138800680322203648.L))+real(4875330082544889626624.L))-
        real(5223567945583810314240.L))+real(5615335541502596087808.L))-
        real(6058651505305432621056.L))+real(6563539130747552006144.L))-
        real(7142674936401747771392.L))+real(7812300711689411624960.L))-
        real(8593530782858352787456.L))+real(9514266223878890586112.L))-
        real(10612066172787993346048.L))+real(11938574444386492514304.L))-
        real(13566561868621014220800.L))+real(15601546148914166353920.L))-
        real(18201803840399860746240.L))+real(21614642060474834636160.L))-
        real(26246351073433727772480.L))+real(32807938841792159715600.L))-
        real(42650320494329807630280.L))+real(58644190679703485491635.L))/
        real(1034483523589969484072441400.L);
      _C4x[88] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(384405449881650462720.L)*_ep2-
        real(404288490392770314240.L))+real(425946802378097295360.L))-
        real(449610513621324922880.L))+real(475549581714862899200.L))-
        real(504082556617754673152.L))+real(535587716406364340224.L))-
        real(570517350085040275456.L))+real(609416260318111203328.L))-
        real(652945993197976289280.L))+real(701916942687824510976.L))-
        real(757331438163179077632.L))+real(820442391343444000768.L))-
        real(892834367050218471424.L))+real(976537588961176453120.L))-
        real(1074191347857294098432.L))+real(1189283277984861323264.L))-
        real(1326508271598499168256.L))+real(1492321805548311564288.L))-
        real(1695820233577626777600.L))+real(1950193268614270794240.L))-
        real(2275225480049982593280.L))+real(2701830257559354329520.L))-
        real(3280793884179215971560.L))+real(4100992355224019964450.L))-
        real(5331290061791225953785.L))/real(57471306866109415781802300.L);
      _C4x[89] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(192202724940825231360.L)*_ep2-
        real(202144245196385157120.L))+real(212973401189048647680.L))-
        real(224805256810662461440.L))+real(237774790857431449600.L))-
        real(252041278308877336576.L))+real(267793858203182170112.L))-
        real(285258675042520137728.L))+real(304708130159055601664.L))-
        real(326472996598988144640.L))+real(350958471343912255488.L))-
        real(378665719081589538816.L))+real(410221195671722000384.L))-
        real(446417183525109235712.L))+real(488268794480588226560.L))-
        real(537095673928647049216.L))+real(594641638992430661632.L))-
        real(663254135799249584128.L))+real(746160902774155782144.L))-
        real(847910116788813388800.L))+real(975096634307135397120.L))-
        real(1137612740024991296640.L))+real(1350915128779677164760.L))-
        real(1640396942089607985780.L))+real(2050496177612009982225.L))/
        real(18808791337999445164953480.L);
      _C4x[90] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(48050681235206307840.L)*_ep2-real(50536061299096289280.L))+
        real(53243350297262161920.L))-real(56201314202665615360.L))+
        real(59443697714357862400.L))-real(63010319577219334144.L))+
        real(66948464550795542528.L))-real(71314668760630034432.L))+
        real(76177032539763900416.L))-real(81618249149747036160.L))+
        real(87739617835978063872.L))-real(94666429770397384704.L))+
        real(102555298917930500096.L))-real(111604295881277308928.L))+
        real(122067198620147056640.L))-real(134273918482161762304.L))+
        real(148660409748107665408.L))-real(165813533949812396032.L))+
        real(186540225693538945536.L))-real(211977529197203347200.L))+
        real(243774158576783849280.L))-real(284403185006247824160.L))+
        real(337728782194919291190.L))-real(410099235522401996445.L))/
        real(3617075257307585608644900.L);
      _C4x[91] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(24025340617603153920.L)*_ep2-real(25268030649548144640.L))+
        real(26621675148631080960.L))-real(28100657101332807680.L))+
        real(29721848857178931200.L))-real(31505159788609667072.L))+
        real(33474232275397771264.L))-real(35657334380315017216.L))+
        real(38088516269881950208.L))-real(40809124574873518080.L))+
        real(43869808917989031936.L))-real(47333214885198692352.L))+
        real(51277649458965250048.L))-real(55802147940638654464.L))+
        real(61033599310073528320.L))-real(67136959241080881152.L))+
        real(74330204874053832704.L))-real(82906766974906198016.L))+
        real(93270112846769472768.L))-real(105988764598601673600.L))+
        real(121887079288391924640.L))-real(142201592503123912080.L))+
        real(168864391097459645595.L))/real(1515726774490797778860720.L);
      _C4x[92] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1501583788600197120.L)*_ep2-real(1579251915596759040.L))+
        real(1663854696789442560.L))-real(1756291068833300480.L))+
        real(1857615553573683200.L))-real(1969072486788104192.L))+
        real(2092139517212360704.L))-real(2228583398769688576.L))+
        real(2380532266867621888.L))-real(2550570285929594880.L))+
        real(2741863057374314496.L))-real(2958325930324918272.L))+
        real(3204853091185328128.L))-real(3487634246289915904.L))+
        real(3814599956879595520.L))-real(4196059952567555072.L))+
        real(4645637804628364544.L))-real(5181672935931637376.L))+
        real(5829382052923092048.L))-real(6624297787412604600.L))+
        real(7617942455524495290.L))-real(8887599531445244505.L))/
        real(83587873593242524569525.L);
      _C4x[93] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(750791894300098560.L)*_ep2-real(789625957798379520.L))+
        real(831927348394721280.L))-real(878145534416650240.L))+
        real(928807776786841600.L))-real(984536243394052096.L))+
        real(1046069758606180352.L))-real(1114291699384844288.L))+
        real(1190266133433810944.L))-real(1275285142964797440.L))+
        real(1370931528687157248.L))-real(1479162965162459136.L))+
        real(1602426545592664064.L))-real(1743817123144957952.L))+
        real(1907299978439797760.L))-real(2098029976283777536.L))+
        real(2322818902314182272.L))-real(2590836467965818688.L))+
        real(2914691026461546024.L))-real(3312148893706302300.L))+
        real(3808971227762247645.L))/real(38127801989900098926450.L);
      _C4x[94] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(187697973575024640.L)*_ep2-real(197406489449594880.L))+
        real(207981837098680320.L))-real(219536383604162560.L))+
        real(232201944196710400.L))-real(246134060848513024.L))+
        real(261517439651545088.L))-real(278572924846211072.L))+
        real(297566533358452736.L))-real(318821285741199360.L))+
        real(342732882171789312.L))-real(369790741290614784.L))+
        real(400606636398166016.L))-real(435954280786239488.L))+
        real(476824994609949440.L))-real(524507494070944384.L))+
        real(580704725578545568.L))-real(647709116991454672.L))+
        real(728672756615386506.L))-real(828037223426575575.L))/
        real(8896487130976689749505.L);
      _C4x[95] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(93848986787512320.L)*
        _ep2-real(98703244724797440.L))+real(103990918549340160.L))-
        real(109768191802081280.L))+real(116100972098355200.L))-
        real(123067030424256512.L))+real(130758719825772544.L))-
        real(139286462423105536.L))+real(148783266679226368.L))-
        real(159410642870599680.L))+real(171366441085894656.L))-
        real(184895370645307392.L))+real(200303318199083008.L))-
        real(217977140393119744.L))+real(238412497304974720.L))-
        real(262253747035472192.L))+real(290352362789272784.L))-
        real(323854558495727336.L))+real(364336378307693253.L))/
        real(4219677690581829130200.L);
      _C4x[96] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11731123348439040.L)*_ep2-
        real(12337905590599680.L))+real(12998864818667520.L))-
        real(13721023975260160.L))+real(14512621512294400.L))-
        real(15383378803032064.L))+real(16344839978221568.L))-
        real(17410807802888192.L))+real(18597908334903296.L))-
        real(19926330358824960.L))+real(21420805135736832.L))-
        real(23111921330663424.L))+real(25037914774885376.L))-
        real(27247142549139968.L))+real(29801562163121840.L))-
        real(32781718379434024.L))+real(36294045348659098.L))-
        real(40481819811965917.L))/real(506361322869819495624.L);
      _C4x[97] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(5865561674219520.L)*_ep2-
        real(6168952795299840.L))+real(6499432409333760.L))-
        real(6860511987630080.L))+real(7256310756147200.L))-
        real(7691689401516032.L))+real(8172419989110784.L))-
        real(8705403901444096.L))+real(9298954167451648.L))-
        real(9963165179412480.L))+real(10710402567868416.L))-
        real(11555960665331712.L))+real(12518957387442688.L))-
        real(13623571274569984.L))+real(14900781081560920.L))-
        real(16390859189717012.L))+real(18147022674329549.L))/
        real(245246224751764428080.L);
      _C4x[98] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(1466390418554880.L)*_ep2-
        real(1542238198824960.L))+real(1624858102333440.L))-
        real(1715127996907520.L))+real(1814077689036800.L))-
        real(1922922350379008.L))+real(2043104997277696.L))-
        real(2176350975361024.L))+real(2324738541862912.L))-
        real(2490791294853120.L))+real(2677600641967104.L))-
        real(2888990166332928.L))+real(3129739346860672.L))-
        real(3405892818642496.L))+real(3725195270390230.L))-
        real(4097714797429253.L))/real(59801419336021375320.L);
      _C4x[99] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(29321986255081448544.L);
      _C4x[100] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(45824700579840.L)*_ep2-real(48194943713280.L))+
        real(50776815697920.L))-real(53597749903360.L))+real(56689927782400.L))-
        real(60091323449344.L))+real(63847031164928.L))-real(68010967980032.L))+
        real(72648079433216.L))-real(77837227964160.L))+real(83675020061472.L))-
        real(90280942697904.L))+real(97804354589396.L))-
        real(106434150582578.L))/real(1804857108504066435.L);
      _C4x[101] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(445905873865710531.L);
      _C4x[102] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(5728087572480.L)*_ep2-real(6024367964160.L))+
        real(6347101962240.L))-real(6699718737920.L))+real(7086240972800.L))-
        real(7511415431168.L))+real(7980878895616.L))-real(8501370997504.L))+
        real(9081009929152.L))-real(9729653495520.L))+real(10459377507684.L))-
        real(11285117837238.L))/real(220944351915442155.L);
      _C4x[103] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(54863374834603720.L);
      _C4x[104] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(6824468576987292.L);
      _C4x[105] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(3400897961953800.L);
      _C4x[106] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(848506865255140.L);
      _C4x[107] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(423860968767600.L);
      _C4x[108] = (_ep2*(_ep2*(_ep2*(_ep2*(real(699229440.L)*_ep2-
        real(735396480.L))+real(774792720.L))-real(817836760.L))+
        real(865019650.L))-real(916920829.L))/real(26491310547975.L);
      _C4x[109] = (_ep2*(_ep2*(_ep2*(real(349614720.L)*_ep2-real(367698240.L))+
        real(387396360.L))-real(408918380.L))+real(432509825.L))/
        real(13256044023222.L);
      _C4x[110] = (_ep2*(_ep2*(real(760032.L)*_ep2-real(799344.L))+
        real(842166.L))-real(888953.L))/real(28859312097.L);
      _C4x[111] = (_ep2*(real(633360.L)*_ep2-real(666120.L))+real(701805.L))/
        real(24098011448.L);
      _C4x[112] = (real(6786.L)*_ep2-real(7137.L))/real(258840080.L);
      _C4x[113] = real(87.L)/real(3328160.L);
      _C4x[114] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(404288490392770314240.L)-real(384405449881650462720.L)*
        _ep2)*_ep2-real(425946802378097295360.L))+
        real(449610513621324922880.L))-real(475549581714862899200.L))+
        real(504082556617754673152.L))-real(535587716406364340224.L))+
        real(570517350085040275456.L))-real(609416260318111203328.L))+
        real(652945993197976289280.L))-real(701916942687824510976.L))+
        real(757331438163179077632.L))-real(820442391343444000768.L))+
        real(892834367050218471424.L))-real(976537588961176453120.L))+
        real(1074191347857294098432.L))-real(1189283277984861323264.L))+
        real(1326508271598499168256.L))-real(1492321805548311564288.L))+
        real(1695820233577626777600.L))-real(1950193268614270794240.L))+
        real(2275225480049982593280.L))-real(2701830257559354329520.L))+
        real(3280793884179215971560.L))-real(4100992355224019964450.L))+
        real(5331290061791225953785.L))/real(665025122307837525475140900.L);
      _C4x[115] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(202144245196385157120.L)-real(192202724940825231360.L)*_ep2)*
        _ep2-real(212973401189048647680.L))+real(224805256810662461440.L))-
        real(237774790857431449600.L))+real(252041278308877336576.L))-
        real(267793858203182170112.L))+real(285258675042520137728.L))-
        real(304708130159055601664.L))+real(326472996598988144640.L))-
        real(350958471343912255488.L))+real(378665719081589538816.L))-
        real(410221195671722000384.L))+real(446417183525109235712.L))-
        real(488268794480588226560.L))+real(537095673928647049216.L))-
        real(594641638992430661632.L))+real(663254135799249584128.L))-
        real(746160902774155782144.L))+real(847910116788813388800.L))-
        real(975096634307135397120.L))+real(1137612740024991296640.L))-
        real(1350915128779677164760.L))+real(1640396942089607985780.L))-
        real(2050496177612009982225.L))/real(120913658601425004631843800.L);
      _C4x[116] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(50536061299096289280.L)-real(48050681235206307840.L)*_ep2)*_ep2-
        real(53243350297262161920.L))+real(56201314202665615360.L))-
        real(59443697714357862400.L))+real(63010319577219334144.L))-
        real(66948464550795542528.L))+real(71314668760630034432.L))-
        real(76177032539763900416.L))+real(81618249149747036160.L))-
        real(87739617835978063872.L))+real(94666429770397384704.L))-
        real(102555298917930500096.L))+real(111604295881277308928.L))-
        real(122067198620147056640.L))+real(134273918482161762304.L))-
        real(148660409748107665408.L))+real(165813533949812396032.L))-
        real(186540225693538945536.L))+real(211977529197203347200.L))-
        real(243774158576783849280.L))+real(284403185006247824160.L))-
        real(337728782194919291190.L))+real(410099235522401996445.L))/
        real(17051926213021475012183100.L);
      _C4x[117] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(25268030649548144640.L)-real(24025340617603153920.L)*_ep2)*_ep2-
        real(26621675148631080960.L))+real(28100657101332807680.L))-
        real(29721848857178931200.L))+real(31505159788609667072.L))-
        real(33474232275397771264.L))+real(35657334380315017216.L))-
        real(38088516269881950208.L))+real(40809124574873518080.L))-
        real(43869808917989031936.L))+real(47333214885198692352.L))-
        real(51277649458965250048.L))+real(55802147940638654464.L))-
        real(61033599310073528320.L))+real(67136959241080881152.L))-
        real(74330204874053832704.L))+real(82906766974906198016.L))-
        real(93270112846769472768.L))+real(105988764598601673600.L))-
        real(121887079288391924640.L))+real(142201592503123912080.L))-
        real(168864391097459645595.L))/real(5846374701607362861319920.L);
      _C4x[118] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1579251915596759040.L)-real(1501583788600197120.L)*_ep2)*_ep2-
        real(1663854696789442560.L))+real(1756291068833300480.L))-
        real(1857615553573683200.L))+real(1969072486788104192.L))-
        real(2092139517212360704.L))+real(2228583398769688576.L))-
        real(2380532266867621888.L))+real(2550570285929594880.L))-
        real(2741863057374314496.L))+real(2958325930324918272.L))-
        real(3204853091185328128.L))+real(3487634246289915904.L))-
        real(3814599956879595520.L))+real(4196059952567555072.L))-
        real(4645637804628364544.L))+real(5181672935931637376.L))-
        real(5829382052923092048.L))+real(6624297787412604600.L))-
        real(7617942455524495290.L))+real(8887599531445244505.L))/
        real(279422320297410724989555.L);
      _C4x[119] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(789625957798379520.L)-real(750791894300098560.L)*_ep2)*_ep2-
        real(831927348394721280.L))+real(878145534416650240.L))-
        real(928807776786841600.L))+real(984536243394052096.L))-
        real(1046069758606180352.L))+real(1114291699384844288.L))-
        real(1190266133433810944.L))+real(1275285142964797440.L))-
        real(1370931528687157248.L))+real(1479162965162459136.L))-
        real(1602426545592664064.L))+real(1743817123144957952.L))-
        real(1907299978439797760.L))+real(2098029976283777536.L))-
        real(2322818902314182272.L))+real(2590836467965818688.L))-
        real(2914691026461546024.L))+real(3312148893706302300.L))-
        real(3808971227762247645.L))/real(114383405969700296779350.L);
      _C4x[120] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(197406489449594880.L)-
        real(187697973575024640.L)*_ep2)*_ep2-real(207981837098680320.L))+
        real(219536383604162560.L))-real(232201944196710400.L))+
        real(246134060848513024.L))-real(261517439651545088.L))+
        real(278572924846211072.L))-real(297566533358452736.L))+
        real(318821285741199360.L))-real(342732882171789312.L))+
        real(369790741290614784.L))-real(400606636398166016.L))+
        real(435954280786239488.L))-real(476824994609949440.L))+
        real(524507494070944384.L))-real(580704725578545568.L))+
        real(647709116991454672.L))-real(728672756615386506.L))+
        real(828037223426575575.L))/real(24510729850650063595575.L);
      _C4x[121] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(98703244724797440.L)-
        real(93848986787512320.L)*_ep2)*_ep2-real(103990918549340160.L))+
        real(109768191802081280.L))-real(116100972098355200.L))+
        real(123067030424256512.L))-real(130758719825772544.L))+
        real(139286462423105536.L))-real(148783266679226368.L))+
        real(159410642870599680.L))-real(171366441085894656.L))+
        real(184895370645307392.L))-real(200303318199083008.L))+
        real(217977140393119744.L))-real(238412497304974720.L))+
        real(262253747035472192.L))-real(290352362789272784.L))+
        real(323854558495727336.L))-real(364336378307693253.L))/
        real(10850599775781846334800.L);
      _C4x[122] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(12337905590599680.L)-
        real(11731123348439040.L)*_ep2)*_ep2-real(12998864818667520.L))+
        real(13721023975260160.L))-real(14512621512294400.L))+
        real(15383378803032064.L))-real(16344839978221568.L))+
        real(17410807802888192.L))-real(18597908334903296.L))+
        real(19926330358824960.L))-real(21420805135736832.L))+
        real(23111921330663424.L))-real(25037914774885376.L))+
        real(27247142549139968.L))-real(29801562163121840.L))+
        real(32781718379434024.L))-real(36294045348659098.L))+
        real(40481819811965917.L))/real(1229734641255275917944.L);
      _C4x[123] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(6168952795299840.L)-
        real(5865561674219520.L)*_ep2)*_ep2-real(6499432409333760.L))+
        real(6860511987630080.L))-real(7256310756147200.L))+
        real(7691689401516032.L))-real(8172419989110784.L))+
        real(8705403901444096.L))-real(9298954167451648.L))+
        real(9963165179412480.L))-real(10710402567868416.L))+
        real(11555960665331712.L))-real(12518957387442688.L))+
        real(13623571274569984.L))-real(14900781081560920.L))+
        real(16390859189717012.L))-real(18147022674329549.L))/
        real(567569834425511962128.L);
      _C4x[124] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(132805749434540976360.L);
      _C4x[125] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(62832827689460246880.L);
      _C4x[126] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(48194943713280.L)-real(45824700579840.L)*_ep2)*_ep2-
        real(50776815697920.L))+real(53597749903360.L))-real(56689927782400.L))+
        real(60091323449344.L))-real(63847031164928.L))+real(68010967980032.L))-
        real(72648079433216.L))+real(77837227964160.L))-real(83675020061472.L))+
        real(90280942697904.L))-real(97804354589396.L))+
        real(106434150582578.L))/real(3748549379200753365.L);
      _C4x[127] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(900911867606231481.L);
      _C4x[128] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(6024367964160.L)-real(5728087572480.L)*_ep2)*_ep2-
        real(6347101962240.L))+real(6699718737920.L))-real(7086240972800.L))+
        real(7511415431168.L))-real(7980878895616.L))+real(8501370997504.L))-
        real(9081009929152.L))+real(9729653495520.L))-real(10459377507684.L))+
        real(11285117837238.L))/real(435576008061871677.L);
      _C4x[129] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(105807937181021460.L);
      _C4x[130] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/real(12903406973295300.L);
      _C4x[131] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(6315953357914200.L);
      _C4x[132] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(1550279460578940.L);
      _C4x[133] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(762949743781680.L);
      _C4x[134] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(47035592197425.L);
      _C4x[135] = (_ep2*(_ep2*((real(73539648.L)-real(69922944.L)*_ep2)*_ep2-
        real(77479272.L))+real(81783676.L))-real(86501965.L))/
        real(4648223228922.L);
      _C4x[136] = (_ep2*((real(18384912.L)-real(17480736.L)*_ep2)*_ep2-
        real(19369818.L))+real(20445919.L))/real(1150249725009.L);
      _C4x[137] = ((real(222040.L)-real(211120.L)*_ep2)*_ep2-real(233935.L))/
        real(13770292256.L);
      _C4x[138] = (real(27755.L)-real(26390.L)*_ep2)/real(1708344528.L);
      _C4x[139] = -real(2639.L)/real(169736160.L);
      _C4x[140] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(192202724940825231360.L)*_ep2-
        real(202144245196385157120.L))+real(212973401189048647680.L))-
        real(224805256810662461440.L))+real(237774790857431449600.L))-
        real(252041278308877336576.L))+real(267793858203182170112.L))-
        real(285258675042520137728.L))+real(304708130159055601664.L))-
        real(326472996598988144640.L))+real(350958471343912255488.L))-
        real(378665719081589538816.L))+real(410221195671722000384.L))-
        real(446417183525109235712.L))+real(488268794480588226560.L))-
        real(537095673928647049216.L))+real(594641638992430661632.L))-
        real(663254135799249584128.L))+real(746160902774155782144.L))-
        real(847910116788813388800.L))+real(975096634307135397120.L))-
        real(1137612740024991296640.L))+real(1350915128779677164760.L))-
        real(1640396942089607985780.L))+real(2050496177612009982225.L))/
        real(1625616965641380617828122200.L);
      _C4x[141] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(48050681235206307840.L)*_ep2-real(50536061299096289280.L))+
        real(53243350297262161920.L))-real(56201314202665615360.L))+
        real(59443697714357862400.L))-real(63010319577219334144.L))+
        real(66948464550795542528.L))-real(71314668760630034432.L))+
        real(76177032539763900416.L))-real(81618249149747036160.L))+
        real(87739617835978063872.L))-real(94666429770397384704.L))+
        real(102555298917930500096.L))-real(111604295881277308928.L))+
        real(122067198620147056640.L))-real(134273918482161762304.L))+
        real(148660409748107665408.L))-real(165813533949812396032.L))+
        real(186540225693538945536.L))-real(211977529197203347200.L))+
        real(243774158576783849280.L))-real(284403185006247824160.L))+
        real(337728782194919291190.L))-real(410099235522401996445.L))/
        real(125047458895490816756009400.L);
      _C4x[142] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(24025340617603153920.L)*_ep2-real(25268030649548144640.L))+
        real(26621675148631080960.L))-real(28100657101332807680.L))+
        real(29721848857178931200.L))-real(31505159788609667072.L))+
        real(33474232275397771264.L))-real(35657334380315017216.L))+
        real(38088516269881950208.L))-real(40809124574873518080.L))+
        real(43869808917989031936.L))-real(47333214885198692352.L))+
        real(51277649458965250048.L))-real(55802147940638654464.L))+
        real(61033599310073528320.L))-real(67136959241080881152.L))+
        real(74330204874053832704.L))-real(82906766974906198016.L))+
        real(93270112846769472768.L))-real(105988764598601673600.L))+
        real(121887079288391924640.L))-real(142201592503123912080.L))+
        real(168864391097459645595.L))/real(30964132678883440339583280.L);
      _C4x[143] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(3003167577200394240.L)*_ep2-real(3158503831193518080.L))+
        real(3327709393578885120.L))-real(3512582137666600960.L))+
        real(3715231107147366400.L))-real(3938144973576208384.L))+
        real(4184279034424721408.L))-real(4457166797539377152.L))+
        real(4761064533735243776.L))-real(5101140571859189760.L))+
        real(5483726114748628992.L))-real(5916651860649836544.L))+
        real(6409706182370656256.L))-real(6975268492579831808.L))+
        real(7629199913759191040.L))-real(8392119905135110144.L))+
        real(9291275609256729088.L))-real(10363345871863274752.L))+
        real(11658764105846184096.L))-real(13248595574825209200.L))+
        real(15235884911048990580.L))-real(17775199062890489010.L))/
        real(2390613184766736202688415.L);
      _C4x[144] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(750791894300098560.L)*_ep2-real(789625957798379520.L))+
        real(831927348394721280.L))-real(878145534416650240.L))+
        real(928807776786841600.L))-real(984536243394052096.L))+
        real(1046069758606180352.L))-real(1114291699384844288.L))+
        real(1190266133433810944.L))-real(1275285142964797440.L))+
        real(1370931528687157248.L))-real(1479162965162459136.L))+
        real(1602426545592664064.L))-real(1743817123144957952.L))+
        real(1907299978439797760.L))-real(2098029976283777536.L))+
        real(2322818902314182272.L))-real(2590836467965818688.L))+
        real(2914691026461546024.L))-real(3312148893706302300.L))+
        real(3808971227762247645.L))/real(419405821888901088190950.L);
      _C4x[145] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(187697973575024640.L)*_ep2-real(197406489449594880.L))+
        real(207981837098680320.L))-real(219536383604162560.L))+
        real(232201944196710400.L))-real(246134060848513024.L))+
        real(261517439651545088.L))-real(278572924846211072.L))+
        real(297566533358452736.L))-real(318821285741199360.L))+
        real(342732882171789312.L))-real(369790741290614784.L))+
        real(400606636398166016.L))-real(435954280786239488.L))+
        real(476824994609949440.L))-real(524507494070944384.L))+
        real(580704725578545568.L))-real(647709116991454672.L))+
        real(728672756615386506.L))-real(828037223426575575.L))/
        real(79886823216933540607800.L);
      _C4x[146] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(93848986787512320.L)*
        _ep2-real(98703244724797440.L))+real(103990918549340160.L))-
        real(109768191802081280.L))+real(116100972098355200.L))-
        real(123067030424256512.L))+real(130758719825772544.L))-
        real(139286462423105536.L))+real(148783266679226368.L))-
        real(159410642870599680.L))+real(171366441085894656.L))-
        real(184895370645307392.L))+real(200303318199083008.L))-
        real(217977140393119744.L))+real(238412497304974720.L))-
        real(262253747035472192.L))+real(290352362789272784.L))-
        real(323854558495727336.L))+real(364336378307693253.L))/
        real(32207335842400083565200.L);
      _C4x[147] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11731123348439040.L)*_ep2-
        real(12337905590599680.L))+real(12998864818667520.L))-
        real(13721023975260160.L))+real(14512621512294400.L))-
        real(15383378803032064.L))+real(16344839978221568.L))-
        real(17410807802888192.L))+real(18597908334903296.L))-
        real(19926330358824960.L))+real(21420805135736832.L))-
        real(23111921330663424.L))+real(25037914774885376.L))-
        real(27247142549139968.L))+real(29801562163121840.L))-
        real(32781718379434024.L))+real(36294045348659098.L))-
        real(40481819811965917.L))/real(3381770263452008774346.L);
      _C4x[148] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(5865561674219520.L)*_ep2-
        real(6168952795299840.L))+real(6499432409333760.L))-
        real(6860511987630080.L))+real(7256310756147200.L))-
        real(7691689401516032.L))+real(8172419989110784.L))-
        real(8705403901444096.L))+real(9298954167451648.L))-
        real(9963165179412480.L))+real(10710402567868416.L))-
        real(11555960665331712.L))+real(12518957387442688.L))-
        real(13623571274569984.L))+real(14900781081560920.L))-
        real(16390859189717012.L))+real(18147022674329549.L))/
        real(1464470313517679013392.L);
      _C4x[149] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(1466390418554880.L)*_ep2-
        real(1542238198824960.L))+real(1624858102333440.L))-
        real(1715127996907520.L))+real(1814077689036800.L))-
        real(1922922350379008.L))+real(2043104997277696.L))-
        real(2176350975361024.L))+real(2324738541862912.L))-
        real(2490791294853120.L))+real(2677600641967104.L))-
        real(2888990166332928.L))+real(3129739346860672.L))-
        real(3405892818642496.L))+real(3725195270390230.L))-
        real(4097714797429253.L))/real(324636276395544608880.L);
      _C4x[150] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(146609931275407242720.L);
      _C4x[151] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(91649401159680.L)*_ep2-real(96389887426560.L))+
        real(101553631395840.L))-real(107195499806720.L))+
        real(113379855564800.L))-real(120182646898688.L))+
        real(127694062329856.L))-real(136021935960064.L))+
        real(145296158866432.L))-real(155674455928320.L))+
        real(167350040122944.L))-real(180561885395808.L))+
        real(195608709178792.L))-real(212868301165156.L))/
        real(16799054625307079895.L);
      _C4x[152] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(1948125662430568929.L);
      _C4x[153] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(2864043786240.L)*_ep2-real(3012183982080.L))+
        real(3173550981120.L))-real(3349859368960.L))+real(3543120486400.L))-
        real(3755707715584.L))+real(3990439447808.L))-real(4250685498752.L))+
        real(4540504964576.L))-real(4864826747760.L))+real(5229688753842.L))-
        real(5642558918619.L))/real(456317722731484614.L);
      _C4x[154] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(215534686850228900.L);
      _C4x[155] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(358005473280.L)*_ep2-real(376522997760.L))+real(396693872640.L))-
        real(418732421120.L))+real(442890060800.L))-real(469463464448.L))+
        real(498804930976.L))-real(531335687344.L))+real(567563120572.L))-
        real(608103343470.L))/real(51255199921700775.L);
      _C4x[156] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(12260380047715800.L);
      _C4x[157] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(2947444900359960.L);
      _C4x[158] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(1423280516177520.L);
      _C4x[159] = (_ep2*(_ep2*(_ep2*(_ep2*(real(1398458880.L)*_ep2-
        real(1470792960.L))+real(1549585440.L))-real(1635673520.L))+
        real(1730039300.L))-real(1833841658.L))/real(172463838057225.L);
      _C4x[160] = (_ep2*(_ep2*(_ep2*(real(69922944.L)*_ep2-real(73539648.L))+
        real(77479272.L))-real(81783676.L))+real(86501965.L))/
        real(8386476831018.L);
      _C4x[161] = (_ep2*(_ep2*(real(17480736.L)*_ep2-real(18384912.L))+
        real(19369818.L))-real(20445919.L))/real(2044888400016.L);
      _C4x[162] = (_ep2*(real(14567280.L)*_ep2-real(15320760.L))+
        real(16141515.L))/real(1666205362976.L);
      _C4x[163] = (real(237510.L)*_ep2-real(249795.L))/real(26621702228.L);
      _C4x[164] = real(1131.L)/real(124473184.L);
      _C4x[165] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(50536061299096289280.L)-real(48050681235206307840.L)*_ep2)*_ep2-
        real(53243350297262161920.L))+real(56201314202665615360.L))-
        real(59443697714357862400.L))+real(63010319577219334144.L))-
        real(66948464550795542528.L))+real(71314668760630034432.L))-
        real(76177032539763900416.L))+real(81618249149747036160.L))-
        real(87739617835978063872.L))+real(94666429770397384704.L))-
        real(102555298917930500096.L))+real(111604295881277308928.L))-
        real(122067198620147056640.L))+real(134273918482161762304.L))-
        real(148660409748107665408.L))+real(165813533949812396032.L))-
        real(186540225693538945536.L))+real(211977529197203347200.L))-
        real(243774158576783849280.L))+real(284403185006247824160.L))-
        real(337728782194919291190.L))+real(410099235522401996445.L))/
        real(1921183686667086184705962600.L);
      _C4x[166] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(25268030649548144640.L)-real(24025340617603153920.L)*_ep2)*_ep2-
        real(26621675148631080960.L))+real(28100657101332807680.L))-
        real(29721848857178931200.L))+real(31505159788609667072.L))-
        real(33474232275397771264.L))+real(35657334380315017216.L))-
        real(38088516269881950208.L))+real(40809124574873518080.L))-
        real(43869808917989031936.L))+real(47333214885198692352.L))-
        real(51277649458965250048.L))+real(55802147940638654464.L))-
        real(61033599310073528320.L))+real(67136959241080881152.L))-
        real(74330204874053832704.L))+real(82906766974906198016.L))-
        real(93270112846769472768.L))+real(105988764598601673600.L))-
        real(121887079288391924640.L))+real(142201592503123912080.L))-
        real(168864391097459645595.L))/real(256157824888944824627461680.L);
      _C4x[167] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3158503831193518080.L)-real(3003167577200394240.L)*_ep2)*_ep2-
        real(3327709393578885120.L))+real(3512582137666600960.L))-
        real(3715231107147366400.L))+real(3938144973576208384.L))-
        real(4184279034424721408.L))+real(4457166797539377152.L))-
        real(4761064533735243776.L))+real(5101140571859189760.L))-
        real(5483726114748628992.L))+real(5916651860649836544.L))-
        real(6409706182370656256.L))+real(6975268492579831808.L))-
        real(7629199913759191040.L))+real(8392119905135110144.L))-
        real(9291275609256729088.L))+real(10363345871863274752.L))-
        real(11658764105846184096.L))+real(13248595574825209200.L))-
        real(15235884911048990580.L))+real(17775199062890489010.L))/
        real(14126350637257986652249725.L);
      _C4x[168] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(789625957798379520.L)-real(750791894300098560.L)*_ep2)*_ep2-
        real(831927348394721280.L))+real(878145534416650240.L))-
        real(928807776786841600.L))+real(984536243394052096.L))-
        real(1046069758606180352.L))+real(1114291699384844288.L))-
        real(1190266133433810944.L))+real(1275285142964797440.L))-
        real(1370931528687157248.L))+real(1479162965162459136.L))-
        real(1602426545592664064.L))+real(1743817123144957952.L))-
        real(1907299978439797760.L))+real(2098029976283777536.L))-
        real(2322818902314182272.L))+real(2590836467965818688.L))-
        real(2914691026461546024.L))+real(3312148893706302300.L))-
        real(3808971227762247645.L))/real(1982645703474805144175400.L);
      _C4x[169] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(197406489449594880.L)-
        real(187697973575024640.L)*_ep2)*_ep2-real(207981837098680320.L))+
        real(219536383604162560.L))-real(232201944196710400.L))+
        real(246134060848513024.L))-real(261517439651545088.L))+
        real(278572924846211072.L))-real(297566533358452736.L))+
        real(318821285741199360.L))-real(342732882171789312.L))+
        real(369790741290614784.L))-real(400606636398166016.L))+
        real(435954280786239488.L))-real(476824994609949440.L))+
        real(524507494070944384.L))-real(580704725578545568.L))+
        real(647709116991454672.L))-real(728672756615386506.L))+
        real(828037223426575575.L))/real(320999780562587499533160.L);
      _C4x[170] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(98703244724797440.L)-
        real(93848986787512320.L)*_ep2)*_ep2-real(103990918549340160.L))+
        real(109768191802081280.L))-real(116100972098355200.L))+
        real(123067030424256512.L))-real(130758719825772544.L))+
        real(139286462423105536.L))-real(148783266679226368.L))+
        real(159410642870599680.L))-real(171366441085894656.L))+
        real(184895370645307392.L))-real(200303318199083008.L))+
        real(217977140393119744.L))-real(238412497304974720.L))+
        real(262253747035472192.L))-real(290352362789272784.L))+
        real(323854558495727336.L))-real(364336378307693253.L))/
        real(114189645259418478094800.L);
      _C4x[171] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(12337905590599680.L)-
        real(11731123348439040.L)*_ep2)*_ep2-real(12998864818667520.L))+
        real(13721023975260160.L))-real(14512621512294400.L))+
        real(15383378803032064.L))-real(16344839978221568.L))+
        real(17410807802888192.L))-real(18597908334903296.L))+
        real(19926330358824960.L))-real(21420805135736832.L))+
        real(23111921330663424.L))-real(25037914774885376.L))+
        real(27247142549139968.L))-real(29801562163121840.L))+
        real(32781718379434024.L))-real(36294045348659098.L))+
        real(40481819811965917.L))/real(10848016299644755419006.L);
      _C4x[172] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(6168952795299840.L)-
        real(5865561674219520.L)*_ep2)*_ep2-real(6499432409333760.L))+
        real(6860511987630080.L))-real(7256310756147200.L))+
        real(7691689401516032.L))-real(8172419989110784.L))+
        real(8705403901444096.L))-real(9298954167451648.L))+
        real(9963165179412480.L))-real(10710402567868416.L))+
        real(11555960665331712.L))-real(12518957387442688.L))+
        real(13623571274569984.L))-real(14900781081560920.L))+
        real(16390859189717012.L))-real(18147022674329549.L))/
        real(4326844108120415266840.L);
      _C4x[173] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(895209125818016951760.L);
      _C4x[174] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(381185821316058831072.L);
      _C4x[175] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(96389887426560.L)-real(91649401159680.L)*_ep2)*_ep2-
        real(101553631395840.L))+real(107195499806720.L))-
        real(113379855564800.L))+real(120182646898688.L))-
        real(127694062329856.L))+real(136021935960064.L))-
        real(145296158866432.L))+real(155674455928320.L))-
        real(167350040122944.L))+real(180561885395808.L))-
        real(195608709178792.L))+real(212868301165156.L))/
        real(41511713495593528005.L);
      _C4x[176] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(4604660656654072014.L);
      _C4x[177] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(1037085733480646850.L);
      _C4x[178] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(473056650359593300.L);
      _C4x[179] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(376522997760.L)-real(358005473280.L)*_ep2)*_ep2-
        real(396693872640.L))+real(418732421120.L))-real(442890060800.L))+
        real(469463464448.L))-real(498804930976.L))+real(531335687344.L))-
        real(567563120572.L))+real(608103343470.L))/real(109033788924345285.L);
      _C4x[180] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(25356695098684950.L);
      _C4x[181] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(5942175013025160.L);
      _C4x[182] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(2803431319743600.L);
      _C4x[183] = (_ep2*(_ep2*(_ep2*((real(1470792960.L)-real(1398458880.L)*
        _ep2)*_ep2-real(1549585440.L))+real(1635673520.L))-real(1730039300.L))+
        real(1833841658.L))/real(332549888694075.L);
      _C4x[184] = (_ep2*(_ep2*((real(367698240.L)-real(349614720.L)*_ep2)*_ep2-
        real(387396360.L))+real(408918380.L))-real(432509825.L))/
        real(79290326402352.L);
      _C4x[185] = (_ep2*((real(1414224.L)-real(1344672.L)*_ep2)*_ep2-
        real(1489986.L))+real(1572763.L))/real(292126914288.L);
      _C4x[186] = ((real(1178520.L)-real(1120560.L)*_ep2)*_ep2-real(1241655.L))/
        real(234094968352.L);
      _C4x[187] = (real(12627.L)-real(12006.L)*_ep2)/real(2420154748.L);
      _C4x[188] = -real(29.L)/real(5657872.L);
      _C4x[189] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(24025340617603153920.L)*_ep2-real(25268030649548144640.L))+
        real(26621675148631080960.L))-real(28100657101332807680.L))+
        real(29721848857178931200.L))-real(31505159788609667072.L))+
        real(33474232275397771264.L))-real(35657334380315017216.L))+
        real(38088516269881950208.L))-real(40809124574873518080.L))+
        real(43869808917989031936.L))-real(47333214885198692352.L))+
        real(51277649458965250048.L))-real(55802147940638654464.L))+
        real(61033599310073528320.L))-real(67136959241080881152.L))+
        real(74330204874053832704.L))-real(82906766974906198016.L))+
        real(93270112846769472768.L))-real(105988764598601673600.L))+
        real(121887079288391924640.L))-real(142201592503123912080.L))+
        real(168864391097459645595.L))/real(4433500815385583503167606000.L);
      _C4x[190] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1501583788600197120.L)*_ep2-real(1579251915596759040.L))+
        real(1663854696789442560.L))-real(1756291068833300480.L))+
        real(1857615553573683200.L))-real(1969072486788104192.L))+
        real(2092139517212360704.L))-real(2228583398769688576.L))+
        real(2380532266867621888.L))-real(2550570285929594880.L))+
        real(2741863057374314496.L))-real(2958325930324918272.L))+
        real(3204853091185328128.L))-real(3487634246289915904.L))+
        real(3814599956879595520.L))-real(4196059952567555072.L))+
        real(4645637804628364544.L))-real(5181672935931637376.L))+
        real(5829382052923092048.L))-real(6624297787412604600.L))+
        real(7617942455524495290.L))-real(8887599531445244505.L))/
        real(65198541402729169164229500.L);
      _C4x[191] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(750791894300098560.L)*_ep2-real(789625957798379520.L))+
        real(831927348394721280.L))-real(878145534416650240.L))+
        real(928807776786841600.L))-real(984536243394052096.L))+
        real(1046069758606180352.L))-real(1114291699384844288.L))+
        real(1190266133433810944.L))-real(1275285142964797440.L))+
        real(1370931528687157248.L))-real(1479162965162459136.L))+
        real(1602426545592664064.L))-real(1743817123144957952.L))+
        real(1907299978439797760.L))-real(2098029976283777536.L))+
        real(2322818902314182272.L))-real(2590836467965818688.L))+
        real(2914691026461546024.L))-real(3312148893706302300.L))+
        real(3808971227762247645.L))/real(12963452676566033634993000.L);
      _C4x[192] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(187697973575024640.L)*_ep2-real(197406489449594880.L))+
        real(207981837098680320.L))-real(219536383604162560.L))+
        real(232201944196710400.L))-real(246134060848513024.L))+
        real(261517439651545088.L))-real(278572924846211072.L))+
        real(297566533358452736.L))-real(318821285741199360.L))+
        real(342732882171789312.L))-real(369790741290614784.L))+
        real(400606636398166016.L))-real(435954280786239488.L))+
        real(476824994609949440.L))-real(524507494070944384.L))+
        real(580704725578545568.L))-real(647709116991454672.L))+
        real(728672756615386506.L))-real(828037223426575575.L))/
        real(1666729629844204324499100.L);
      _C4x[193] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(93848986787512320.L)*
        _ep2-real(98703244724797440.L))+real(103990918549340160.L))-
        real(109768191802081280.L))+real(116100972098355200.L))-
        real(123067030424256512.L))+real(130758719825772544.L))-
        real(139286462423105536.L))+real(148783266679226368.L))-
        real(159410642870599680.L))+real(171366441085894656.L))-
        real(184895370645307392.L))+real(200303318199083008.L))-
        real(217977140393119744.L))+real(238412497304974720.L))-
        real(262253747035472192.L))+real(290352362789272784.L))-
        real(323854558495727336.L))+real(364336378307693253.L))/
        real(500677675368219480877200.L);
      _C4x[194] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11731123348439040.L)*_ep2-
        real(12337905590599680.L))+real(12998864818667520.L))-
        real(13721023975260160.L))+real(14512621512294400.L))-
        real(15383378803032064.L))+real(16344839978221568.L))-
        real(17410807802888192.L))+real(18597908334903296.L))-
        real(19926330358824960.L))+real(21420805135736832.L))-
        real(23111921330663424.L))+real(25037914774885376.L))-
        real(27247142549139968.L))+real(29801562163121840.L))-
        real(32781718379434024.L))+real(36294045348659098.L))-
        real(40481819811965917.L))/real(41723139614018290073100.L);
      _C4x[195] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(5865561674219520.L)*_ep2-
        real(6168952795299840.L))+real(6499432409333760.L))-
        real(6860511987630080.L))+real(7256310756147200.L))-
        real(7691689401516032.L))+real(8172419989110784.L))-
        real(8705403901444096.L))+real(9298954167451648.L))-
        real(9963165179412480.L))+real(10710402567868416.L))-
        real(11555960665331712.L))+real(12518957387442688.L))-
        real(13623571274569984.L))+real(14900781081560920.L))-
        real(16390859189717012.L))+real(18147022674329549.L))/
        real(14977537297339899000600.L);
      _C4x[196] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(1466390418554880.L)*_ep2-
        real(1542238198824960.L))+real(1624858102333440.L))-
        real(1715127996907520.L))+real(1814077689036800.L))-
        real(1922922350379008.L))+real(2043104997277696.L))-
        real(2176350975361024.L))+real(2324738541862912.L))-
        real(2490791294853120.L))+real(2677600641967104.L))-
        real(2888990166332928.L))+real(3129739346860672.L))-
        real(3405892818642496.L))+real(3725195270390230.L))-
        real(4097714797429253.L))/real(2840567418461015327700.L);
      _C4x[197] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(1124009473111455527520.L);
      _C4x[198] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(22912350289920.L)*_ep2-real(24097471856640.L))+
        real(25388407848960.L))-real(26798874951680.L))+real(28344963891200.L))-
        real(30045661724672.L))+real(31923515582464.L))-real(34005483990016.L))+
        real(36324039716608.L))-real(38918613982080.L))+real(41837510030736.L))-
        real(45140471348952.L))+real(48902177294698.L))-real(53217075291289.L))/
        real(28738878573872442465.L);
      _C4x[199] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(12075159064652286750.L);
      _C4x[200] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(2864043786240.L)*_ep2-real(3012183982080.L))+
        real(3173550981120.L))-real(3349859368960.L))+real(3543120486400.L))-
        real(3755707715584.L))+real(3990439447808.L))-real(4250685498752.L))+
        real(4540504964576.L))-real(4864826747760.L))+real(5229688753842.L))-
        real(5642558918619.L))/real(2592714333701617125.L);
      _C4x[201] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(1133656469796658500.L);
      _C4x[202] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(125808217989629175.L);
      _C4x[203] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(56564935220143350.L);
      _C4x[204] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(12855667095487125.L);
      _C4x[205] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(5898622460094000.L);
      _C4x[206] = (_ep2*(_ep2*(_ep2*(_ep2*(real(699229440.L)*_ep2-
        real(735396480.L))+real(774792720.L))-real(817836760.L))+
        real(865019650.L))-real(916920829.L))/real(341076808917000.L);
      _C4x[207] = (_ep2*(_ep2*(_ep2*(real(69922944.L)*_ep2-real(73539648.L))+
        real(77479272.L))-real(81783676.L))+real(86501965.L))/
        real(31780333254384.L);
      _C4x[208] = (_ep2*(_ep2*(real(17480736.L)*_ep2-real(18384912.L))+
        real(19369818.L))-real(20445919.L))/real(7449236314344.L);
      _C4x[209] = (_ep2*(real(2913456.L)*_ep2-real(3064152.L))+real(3228303.L))/
        real(1170474841760.L);
      _C4x[210] = (real(17342.L)*_ep2-real(18239.L))/real(6600422040.L);
      _C4x[211] = real(8671.L)/real(3140118960.L);
      _C4x[212] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1579251915596759040.L)-real(1501583788600197120.L)*_ep2)*_ep2-
        real(1663854696789442560.L))+real(1756291068833300480.L))-
        real(1857615553573683200.L))+real(1969072486788104192.L))-
        real(2092139517212360704.L))+real(2228583398769688576.L))-
        real(2380532266867621888.L))+real(2550570285929594880.L))-
        real(2741863057374314496.L))+real(2958325930324918272.L))-
        real(3204853091185328128.L))+real(3487634246289915904.L))-
        real(3814599956879595520.L))+real(4196059952567555072.L))-
        real(4645637804628364544.L))+real(5181672935931637376.L))-
        real(5829382052923092048.L))+real(6624297787412604600.L))-
        real(7617942455524495290.L))+real(8887599531445244505.L))/
        real(1256158564359248659230821700.L);
      _C4x[213] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(789625957798379520.L)-real(750791894300098560.L)*_ep2)*_ep2-
        real(831927348394721280.L))+real(878145534416650240.L))-
        real(928807776786841600.L))+real(984536243394052096.L))-
        real(1046069758606180352.L))+real(1114291699384844288.L))-
        real(1190266133433810944.L))+real(1275285142964797440.L))-
        real(1370931528687157248.L))+real(1479162965162459136.L))-
        real(1602426545592664064.L))+real(1743817123144957952.L))-
        real(1907299978439797760.L))+real(2098029976283777536.L))-
        real(2322818902314182272.L))+real(2590836467965818688.L))-
        real(2914691026461546024.L))+real(3312148893706302300.L))-
        real(3808971227762247645.L))/real(132227217300973543076928600.L);
      _C4x[214] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(197406489449594880.L)-
        real(187697973575024640.L)*_ep2)*_ep2-real(207981837098680320.L))+
        real(219536383604162560.L))-real(232201944196710400.L))+
        real(246134060848513024.L))-real(261517439651545088.L))+
        real(278572924846211072.L))-real(297566533358452736.L))+
        real(318821285741199360.L))-real(342732882171789312.L))+
        real(369790741290614784.L))-real(400606636398166016.L))+
        real(435954280786239488.L))-real(476824994609949440.L))+
        real(524507494070944384.L))-real(580704725578545568.L))+
        real(647709116991454672.L))-real(728672756615386506.L))+
        real(828037223426575575.L))/real(11963414898659511040293540.L);
      _C4x[215] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(98703244724797440.L)-
        real(93848986787512320.L)*_ep2)*_ep2-real(103990918549340160.L))+
        real(109768191802081280.L))-real(116100972098355200.L))+
        real(123067030424256512.L))-real(130758719825772544.L))+
        real(139286462423105536.L))-real(148783266679226368.L))+
        real(159410642870599680.L))-real(171366441085894656.L))+
        real(184895370645307392.L))-real(200303318199083008.L))+
        real(217977140393119744.L))-real(238412497304974720.L))+
        real(262253747035472192.L))-real(290352362789272784.L))+
        real(323854558495727336.L))-real(364336378307693253.L))/
        real(2837173493753243724970800.L);
      _C4x[216] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(12337905590599680.L)-
        real(11731123348439040.L)*_ep2)*_ep2-real(12998864818667520.L))+
        real(13721023975260160.L))-real(14512621512294400.L))+
        real(15383378803032064.L))-real(16344839978221568.L))+
        real(17410807802888192.L))-real(18597908334903296.L))+
        real(19926330358824960.L))-real(21420805135736832.L))+
        real(23111921330663424.L))-real(25037914774885376.L))+
        real(27247142549139968.L))-real(29801562163121840.L))+
        real(32781718379434024.L))-real(36294045348659098.L))+
        real(40481819811965917.L))/real(198602144562727060747956.L);
      _C4x[217] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(6168952795299840.L)-
        real(5865561674219520.L)*_ep2)*_ep2-real(6499432409333760.L))+
        real(6860511987630080.L))-real(7256310756147200.L))+
        real(7691689401516032.L))-real(8172419989110784.L))+
        real(8705403901444096.L))-real(9298954167451648.L))+
        real(9963165179412480.L))-real(10710402567868416.L))+
        real(11555960665331712.L))-real(12518957387442688.L))+
        real(13623571274569984.L))-real(14900781081560920.L))+
        real(16390859189717012.L))-real(18147022674329549.L))/
        real(62239988324501358069160.L);
      _C4x[218] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(10577732005888161839340.L);
      _C4x[219] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(3821632208578948793568.L);
      _C4x[220] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(24097471856640.L)-real(22912350289920.L)*_ep2)*_ep2-
        real(25388407848960.L))+real(26798874951680.L))-real(28344963891200.L))+
        real(30045661724672.L))-real(31923515582464.L))+real(34005483990016.L))-
        real(36324039716608.L))+real(38918613982080.L))-real(41837510030736.L))+
        real(45140471348952.L))-real(48902177294698.L))+real(53217075291289.L))/
        real(90474247362191022575.L);
      _C4x[221] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(35581468710508738290.L);
      _C4x[222] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(7212459873751771275.L);
      _C4x[223] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(2997891553462274700.L);
      _C4x[224] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/real(318068981891729145.L);
      _C4x[225] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(137371985534633850.L);
      _C4x[226] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(30110829152540955.L);
      _C4x[227] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(13370210909546400.L);
      _C4x[228] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(750368979617400.L);
      _C4x[229] = (_ep2*(_ep2*((real(367698240.L)-real(349614720.L)*_ep2)*_ep2-
        real(387396360.L))+real(408918380.L))-real(432509825.L))/
        real(340167270759888.L);
      _C4x[230] = (_ep2*((real(18384912.L)-real(17480736.L)*_ep2)*_ep2-
        real(19369818.L))+real(20445919.L))/real(15551914410648.L);
      _C4x[231] = ((real(5106920.L)-real(4855760.L)*_ep2)*_ep2-real(5380505.L))/
        real(3979614461984.L);
      _C4x[232] = (real(383019.L)-real(364182.L)*_ep2)/real(276777697544.L);
      _C4x[233] = -real(8671.L)/real(6147020752.L);
      _C4x[234] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(750791894300098560.L)*_ep2-real(789625957798379520.L))+
        real(831927348394721280.L))-real(878145534416650240.L))+
        real(928807776786841600.L))-real(984536243394052096.L))+
        real(1046069758606180352.L))-real(1114291699384844288.L))+
        real(1190266133433810944.L))-real(1275285142964797440.L))+
        real(1370931528687157248.L))-real(1479162965162459136.L))+
        real(1602426545592664064.L))-real(1743817123144957952.L))+
        real(1907299978439797760.L))-real(2098029976283777536.L))+
        real(2322818902314182272.L))-real(2590836467965818688.L))+
        real(2914691026461546024.L))-real(3312148893706302300.L))+
        real(3808971227762247645.L))/real(2807883849744202885339483800.L);
      _C4x[235] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(187697973575024640.L)*_ep2-real(197406489449594880.L))+
        real(207981837098680320.L))-real(219536383604162560.L))+
        real(232201944196710400.L))-real(246134060848513024.L))+
        real(261517439651545088.L))-real(278572924846211072.L))+
        real(297566533358452736.L))-real(318821285741199360.L))+
        real(342732882171789312.L))-real(369790741290614784.L))+
        real(400606636398166016.L))-real(435954280786239488.L))+
        real(476824994609949440.L))-real(524507494070944384.L))+
        real(580704725578545568.L))-real(647709116991454672.L))+
        real(728672756615386506.L))-real(828037223426575575.L))/
        real(133708754749723946920927800.L);
      _C4x[236] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(93848986787512320.L)*
        _ep2-real(98703244724797440.L))+real(103990918549340160.L))-
        real(109768191802081280.L))+real(116100972098355200.L))-
        real(123067030424256512.L))+real(130758719825772544.L))-
        real(139286462423105536.L))+real(148783266679226368.L))-
        real(159410642870599680.L))+real(171366441085894656.L))-
        real(184895370645307392.L))+real(200303318199083008.L))-
        real(217977140393119744.L))+real(238412497304974720.L))-
        real(262253747035472192.L))+real(290352362789272784.L))-
        real(323854558495727336.L))+real(364336378307693253.L))/
        real(22196710274657730318889200.L);
      _C4x[237] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11731123348439040.L)*_ep2-
        real(12337905590599680.L))+real(12998864818667520.L))-
        real(13721023975260160.L))+real(14512621512294400.L))-
        real(15383378803032064.L))+real(16344839978221568.L))-
        real(17410807802888192.L))+real(18597908334903296.L))-
        real(19926330358824960.L))+real(21420805135736832.L))-
        real(23111921330663424.L))+real(25037914774885376.L))-
        real(27247142549139968.L))+real(29801562163121840.L))-
        real(32781718379434024.L))+real(36294045348659098.L))-
        real(40481819811965917.L))/real(1220819065106175167538906.L);
      _C4x[238] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(5865561674219520.L)*_ep2-
        real(6168952795299840.L))+real(6499432409333760.L))-
        real(6860511987630080.L))+real(7256310756147200.L))-
        real(7691689401516032.L))+real(8172419989110784.L))-
        real(8705403901444096.L))+real(9298954167451648.L))-
        real(9963165179412480.L))+real(10710402567868416.L))-
        real(11555960665331712.L))+real(12518957387442688.L))-
        real(13623571274569984.L))+real(14900781081560920.L))-
        real(16390859189717012.L))+real(18147022674329549.L))/
        real(319986763503612864426152.L);
      _C4x[239] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(1466390418554880.L)*_ep2-
        real(1542238198824960.L))+real(1624858102333440.L))-
        real(1715127996907520.L))+real(1814077689036800.L))-
        real(1922922350379008.L))+real(2043104997277696.L))-
        real(2176350975361024.L))+real(2324738541862912.L))-
        real(2490791294853120.L))+real(2677600641967104.L))-
        real(2888990166332928.L))+real(3129739346860672.L))-
        real(3405892818642496.L))+real(3725195270390230.L))-
        real(4097714797429253.L))/real(47288684261617664693520.L);
      _C4x[240] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(15254414277941182159200.L);
      _C4x[241] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(91649401159680.L)*_ep2-real(96389887426560.L))+
        real(101553631395840.L))-real(107195499806720.L))+
        real(113379855564800.L))-real(120182646898688.L))+
        real(127694062329856.L))-real(136021935960064.L))+
        real(145296158866432.L))-real(155674455928320.L))+
        real(167350040122944.L))-real(180561885395808.L))+
        real(195608709178792.L))-real(212868301165156.L))/
        real(1314537594027128386825.L);
      _C4x[242] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(119302571558764593090.L);
      _C4x[243] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(2864043786240.L)*_ep2-real(3012183982080.L))+
        real(3173550981120.L))-real(3349859368960.L))+real(3543120486400.L))-
        real(3755707715584.L))+real(3990439447808.L))-real(4250685498752.L))+
        real(4540504964576.L))-real(4864826747760.L))+real(5229688753842.L))-
        real(5642558918619.L))/real(22570756781387895990.L);
      _C4x[244] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(8833359604051943100.L);
      _C4x[245] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(358005473280.L)*_ep2-real(376522997760.L))+real(396693872640.L))-
        real(418732421120.L))+real(442890060800.L))-real(469463464448.L))+
        real(498804930976.L))-real(531335687344.L))+real(567563120572.L))-
        real(608103343470.L))/real(1777444310571427575.L);
      _C4x[246] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(366118097194114650.L);
      _C4x[247] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(76921782036743280.L);
      _C4x[248] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(32874989177590560.L);
      _C4x[249] = (_ep2*(_ep2*(_ep2*(_ep2*(real(699229440.L)*_ep2-
        real(735396480.L))+real(774792720.L))-real(817836760.L))+
        real(865019650.L))-real(916920829.L))/real(1782126326591325.L);
      _C4x[250] = (_ep2*(_ep2*(_ep2*(real(69922944.L)*_ep2-real(73539648.L))+
        real(77479272.L))-real(81783676.L))+real(86501965.L))/
        real(156547567512336.L);
      _C4x[251] = (_ep2*(_ep2*(real(17480736.L)*_ep2-real(18384912.L))+
        real(19369818.L))-real(20445919.L))/real(34763102800272.L);
      _C4x[252] = (_ep2*(real(4855760.L)*_ep2-real(5106920.L))+real(5380505.L))/
        real(8661513829024.L);
      _C4x[253] = (real(1820910.L)*_ep2-real(1915095.L))/real(2938727906276.L);
      _C4x[254] = real(4669.L)/real(6870199664.L);
      _C4x[255] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(197406489449594880.L)-
        real(187697973575024640.L)*_ep2)*_ep2-real(207981837098680320.L))+
        real(219536383604162560.L))-real(232201944196710400.L))+
        real(246134060848513024.L))-real(261517439651545088.L))+
        real(278572924846211072.L))-real(297566533358452736.L))+
        real(318821285741199360.L))-real(342732882171789312.L))+
        real(369790741290614784.L))-real(400606636398166016.L))+
        real(435954280786239488.L))-real(476824994609949440.L))+
        real(524507494070944384.L))-real(580704725578545568.L))+
        real(647709116991454672.L))-real(728672756615386506.L))+
        real(828037223426575575.L))/real(3103450570769908452217324200.L);
      _C4x[256] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(98703244724797440.L)-
        real(93848986787512320.L)*_ep2)*_ep2-real(103990918549340160.L))+
        real(109768191802081280.L))-real(116100972098355200.L))+
        real(123067030424256512.L))-real(130758719825772544.L))+
        real(139286462423105536.L))-real(148783266679226368.L))+
        real(159410642870599680.L))-real(171366441085894656.L))+
        real(184895370645307392.L))-real(200303318199083008.L))+
        real(217977140393119744.L))-real(238412497304974720.L))+
        real(262253747035472192.L))-real(290352362789272784.L))+
        real(323854558495727336.L))-real(364336378307693253.L))/
        real(269865267023470300192810800.L);
      _C4x[257] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(12337905590599680.L)-
        real(11731123348439040.L)*_ep2)*_ep2-real(12998864818667520.L))+
        real(13721023975260160.L))-real(14512621512294400.L))+
        real(15383378803032064.L))-real(16344839978221568.L))+
        real(17410807802888192.L))-real(18597908334903296.L))+
        real(19926330358824960.L))-real(21420805135736832.L))+
        real(23111921330663424.L))-real(25037914774885376.L))+
        real(27247142549139968.L))-real(29801562163121840.L))+
        real(32781718379434024.L))-real(36294045348659098.L))+
        real(40481819811965917.L))/real(10344835235899694840724414.L);
      _C4x[258] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(6168952795299840.L)-
        real(5865561674219520.L)*_ep2)*_ep2-real(6499432409333760.L))+
        real(6860511987630080.L))-real(7256310756147200.L))+
        real(7691689401516032.L))-real(8172419989110784.L))+
        real(8705403901444096.L))-real(9298954167451648.L))+
        real(9963165179412480.L))-real(10710402567868416.L))+
        real(11555960665331712.L))-real(12518957387442688.L))+
        real(13623571274569984.L))-real(14900781081560920.L))+
        real(16390859189717012.L))-real(18147022674329549.L))/
        real(2122017484287116890405008.L);
      _C4x[259] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(261332202498413410148400.L);
      _C4x[260] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(73060615752244609288800.L);
      _C4x[261] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(96389887426560.L)-real(91649401159680.L)*_ep2)*_ep2-
        real(101553631395840.L))+real(107195499806720.L))-
        real(113379855564800.L))+real(120182646898688.L))-
        real(127694062329856.L))+real(136021935960064.L))-
        real(145296158866432.L))+real(155674455928320.L))-
        real(167350040122944.L))+real(180561885395808.L))-
        real(195608709178792.L))+real(212868301165156.L))/
        real(5604081321905126280675.L);
      _C4x[262] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(461512579451010399585.L);
      _C4x[263] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(80383572396872682210.L);
      _C4x[264] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(29289560792382758700.L);
      _C4x[265] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(376522997760.L)-real(358005473280.L)*_ep2)*_ep2-
        real(396693872640.L))+real(418732421120.L))-real(442890060800.L))+
        real(469463464448.L))-real(498804930976.L))+real(531335687344.L))-
        real(567563120572.L))+real(608103343470.L))/
        real(5536441369291862925.L);
      _C4x[266] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(1079084918045811600.L);
      _C4x[267] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(215816983609162320.L);
      _C4x[268] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(88243392003006240.L);
      _C4x[269] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(4596010000156575.L);
      _C4x[270] = (_ep2*(_ep2*((real(73539648.L)-real(69922944.L)*_ep2)*_ep2-
        real(77479272.L))+real(81783676.L))-real(86501965.L))/
        real(389309082366204.L);
      _C4x[271] = (_ep2*((real(18384912.L)-real(17480736.L)*_ep2)*_ep2-
        real(19369818.L))+real(20445919.L))/real(83625172990128.L);
      _C4x[272] = ((real(15320760.L)-real(14567280.L)*_ep2)*_ep2-
        real(16141515.L))/real(60630596803168.L);
      _C4x[273] = (real(7015.L)-real(6670.L)*_ep2)/real(24421561548.L);
      _C4x[274] = -real(667.L)/real(2169536736.L);
      _C4x[275] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(93848986787512320.L)*
        _ep2-real(98703244724797440.L))+real(103990918549340160.L))-
        real(109768191802081280.L))+real(116100972098355200.L))-
        real(123067030424256512.L))+real(130758719825772544.L))-
        real(139286462423105536.L))+real(148783266679226368.L))-
        real(159410642870599680.L))+real(171366441085894656.L))-
        real(184895370645307392.L))+real(200303318199083008.L))-
        real(217977140393119744.L))+real(238412497304974720.L))-
        real(262253747035472192.L))+real(290352362789272784.L))-
        real(323854558495727336.L))+real(364336378307693253.L))/
        real(6798034583591228038190329200.L);
      _C4x[276] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11731123348439040.L)*_ep2-
        real(12337905590599680.L))+real(12998864818667520.L))-
        real(13721023975260160.L))+real(14512621512294400.L))-
        real(15383378803032064.L))+real(16344839978221568.L))-
        real(17410807802888192.L))+real(18597908334903296.L))-
        real(19926330358824960.L))+real(21420805135736832.L))-
        real(23111921330663424.L))+real(25037914774885376.L))-
        real(27247142549139968.L))+real(29801562163121840.L))-
        real(32781718379434024.L))+real(36294045348659098.L))-
        real(40481819811965917.L))/real(135960691671824560763806584.L);
      _C4x[277] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(5865561674219520.L)*_ep2-
        real(6168952795299840.L))+real(6499432409333760.L))-
        real(6860511987630080.L))+real(7256310756147200.L))-
        real(7691689401516032.L))+real(8172419989110784.L))-
        real(8705403901444096.L))+real(9298954167451648.L))-
        real(9963165179412480.L))+real(10710402567868416.L))-
        real(11555960665331712.L))+real(12518957387442688.L))-
        real(13623571274569984.L))+real(14900781081560920.L))-
        real(16390859189717012.L))+real(18147022674329549.L))/
        real(19367619896271304952109200.L);
      _C4x[278] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(1466390418554880.L)*_ep2-
        real(1542238198824960.L))+real(1624858102333440.L))-
        real(1715127996907520.L))+real(1814077689036800.L))-
        real(1922922350379008.L))+real(2043104997277696.L))-
        real(2176350975361024.L))+real(2324738541862912.L))-
        real(2490791294853120.L))+real(2677600641967104.L))-
        real(2888990166332928.L))+real(3129739346860672.L))-
        real(3405892818642496.L))+real(3725195270390230.L))-
        real(4097714797429253.L))/real(1860436393976800229389800.L);
      _C4x[279] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(432101356020418117793760.L);
      _C4x[280] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(45824700579840.L)*_ep2-real(48194943713280.L))+
        real(50776815697920.L))-real(53597749903360.L))+real(56689927782400.L))-
        real(60091323449344.L))+real(63847031164928.L))-real(68010967980032.L))+
        real(72648079433216.L))-real(77837227964160.L))+real(83675020061472.L))-
        real(90280942697904.L))+real(97804354589396.L))-
        real(106434150582578.L))/real(14321541155979767161725.L);
      _C4x[281] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(2094074084991999568185.L);
      _C4x[282] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(5728087572480.L)*_ep2-real(6024367964160.L))+
        real(6347101962240.L))-real(6699718737920.L))+real(7086240972800.L))-
        real(7511415431168.L))+real(7980878895616.L))-real(8501370997504.L))+
        real(9081009929152.L))-real(9729653495520.L))+real(10459377507684.L))-
        real(11285117837238.L))/real(660293630402882746725.L);
      _C4x[283] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(110494480661211147900.L);
      _C4x[284] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(9701954399520978840.L);
      _C4x[285] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(3545564730721952400.L);
      _C4x[286] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(669717782469702120.L);
      _C4x[287] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(486420480.L)*_ep2-
        real(511580160.L))+real(538986240.L))-real(568929920.L))+
        real(601752800.L))-real(637857968.L))+real(677724091.L))/
        real(11313255385000800.L);
      _C4x[288] = (_ep2*(_ep2*(_ep2*(_ep2*(real(30401280.L)*_ep2-
        real(31973760.L))+real(33686640.L))-real(35558120.L))+real(37609550.L))-
        real(39866123.L))/real(562776734713050.L);
      _C4x[289] = (_ep2*(_ep2*(_ep2*(real(15200640.L)*_ep2-real(15986880.L))+
        real(16843320.L))-real(17779060.L))+real(18804775.L))/
        real(228641842024596.L);
      _C4x[290] = (_ep2*(_ep2*(real(760032.L)*_ep2-real(799344.L))+
        real(842166.L))-real(888953.L))/real(9457608850074.L);
      _C4x[291] = (_ep2*(real(48720.L)*_ep2-real(51240.L))+real(53985.L))/
        real(509500813472.L);
      _C4x[292] = (real(3654.L)*_ep2-real(3843.L))/real(32562082064.L);
      _C4x[293] = real(203.L)/real(1560543968.L);
      _C4x[294] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(12337905590599680.L)-
        real(11731123348439040.L)*_ep2)*_ep2-real(12998864818667520.L))+
        real(13721023975260160.L))-real(14512621512294400.L))+
        real(15383378803032064.L))-real(16344839978221568.L))+
        real(17410807802888192.L))-real(18597908334903296.L))+
        real(19926330358824960.L))-real(21420805135736832.L))+
        real(23111921330663424.L))-real(25037914774885376.L))+
        real(27247142549139968.L))-real(29801562163121840.L))+
        real(32781718379434024.L))-real(36294045348659098.L))+
        real(40481819811965917.L))/real(3694584012821319585973005000.L);
      _C4x[295] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*((real(6168952795299840.L)-
        real(5865561674219520.L)*_ep2)*_ep2-real(6499432409333760.L))+
        real(6860511987630080.L))-real(7256310756147200.L))+
        real(7691689401516032.L))-real(8172419989110784.L))+
        real(8705403901444096.L))-real(9298954167451648.L))+
        real(9963165179412480.L))-real(10710402567868416.L))+
        real(11555960665331712.L))-real(12518957387442688.L))+
        real(13623571274569984.L))-real(14900781081560920.L))+
        real(16390859189717012.L))-real(18147022674329549.L))/
        real(273672889838616265627630000.L);
      _C4x[296] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(18199921245425219635335000.L);
      _C4x[297] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(3287727708851007417996000.L);
      _C4x[298] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(48194943713280.L)-real(45824700579840.L)*_ep2)*_ep2-
        real(50776815697920.L))+real(53597749903360.L))-real(56689927782400.L))+
        real(60091323449344.L))-real(63847031164928.L))+real(68010967980032.L))-
        real(72648079433216.L))+real(77837227964160.L))-real(83675020061472.L))+
        real(90280942697904.L))-real(97804354589396.L))+
        real(106434150582578.L))/real(90287976852915923410875.L);
      _C4x[299] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(11380837418434780261875.L);
      _C4x[300] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(6024367964160.L)-real(5728087572480.L)*_ep2)*_ep2-
        real(6347101962240.L))+real(6699718737920.L))-real(7086240972800.L))+
        real(7511415431168.L))-real(7980878895616.L))+real(8501370997504.L))-
        real(9081009929152.L))+real(9729653495520.L))-real(10459377507684.L))+
        real(11285117837238.L))/real(3178432071815118811875.L);
      _C4x[301] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(480410785483526730000.L);
      _C4x[302] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/
        real(38667209563308249000.L);
      _C4x[303] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(13103174004841998000.L);
      _C4x[304] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(2316217627118535000.L);
      _C4x[305] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(848494153875060000.L);
      _C4x[306] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(40043729200736250.L);
      _C4x[307] = (_ep2*(_ep2*((real(73539648.L)-real(69922944.L)*_ep2)*_ep2-
        real(77479272.L))+real(81783676.L))-real(86501965.L))/
        real(3102996427476660.L);
      _C4x[308] = (_ep2*((real(1414224.L)-real(1344672.L)*_ep2)*_ep2-
        real(1489986.L))+real(1572763.L))/real(47288044250370.L);
      _C4x[309] = ((real(235704.L)-real(224112.L)*_ep2)*_ep2-real(248331.L))/
        real(6368760168400.L);
      _C4x[310] = (real(88389.L)-real(84042.L)*_ep2)/real(1963302006800.L);
      _C4x[311] = -real(2001.L)/real(39013599200.L);
      _C4x[312] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(5865561674219520.L)*_ep2-
        real(6168952795299840.L))+real(6499432409333760.L))-
        real(6860511987630080.L))+real(7256310756147200.L))-
        real(7691689401516032.L))+real(8172419989110784.L))-
        real(8705403901444096.L))+real(9298954167451648.L))-
        real(9963165179412480.L))+real(10710402567868416.L))-
        real(11555960665331712.L))+real(12518957387442688.L))-
        real(13623571274569984.L))+real(14900781081560920.L))-
        real(16390859189717012.L))+real(18147022674329549.L))/
        real(7980301467694050305701690800.L);
      _C4x[313] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(_ep2*(real(1466390418554880.L)*_ep2-
        real(1542238198824960.L))+real(1624858102333440.L))-
        real(1715127996907520.L))+real(1814077689036800.L))-
        real(1922922350379008.L))+real(2043104997277696.L))-
        real(2176350975361024.L))+real(2324738541862912.L))-
        real(2490791294853120.L))+real(2677600641967104.L))-
        real(2888990166332928.L))+real(3129739346860672.L))-
        real(3405892818642496.L))+real(3725195270390230.L))-
        real(4097714797429253.L))/real(275182809230829320886265200.L);
      _C4x[314] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(34323877280404517443878240.L);
      _C4x[315] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(91649401159680.L)*_ep2-real(96389887426560.L))+
        real(101553631395840.L))-real(107195499806720.L))+
        real(113379855564800.L))-real(120182646898688.L))+
        real(127694062329856.L))-real(136021935960064.L))+
        real(145296158866432.L))-real(155674455928320.L))+
        real(167350040122944.L))-real(180561885395808.L))+
        real(195608709178792.L))-real(212868301165156.L))/
        real(1462665225017237959256175.L);
      _C4x[316] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(76206087353839288633515.L);
      _C4x[317] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(2864043786240.L)*_ep2-real(3012183982080.L))+
        real(3173550981120.L))-real(3349859368960.L))+real(3543120486400.L))-
        real(3755707715584.L))+real(3990439447808.L))-real(4250685498752.L))+
        real(4540504964576.L))-real(4864826747760.L))+real(5229688753842.L))-
        real(5642558918619.L))/real(9153884366827542178200.L);
      _C4x[318] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(2445977199233270379600.L);
      _C4x[319] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(177482491895584862910.L);
      _C4x[320] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(55033330820336391600.L);
      _C4x[321] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(9005454134236864080.L);
      _C4x[322] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(3082347853531581600.L);
      _C4x[323] = (_ep2*(_ep2*(_ep2*(_ep2*(real(699229440.L)*_ep2-
        real(735396480.L))+real(774792720.L))-real(817836760.L))+
        real(865019650.L))-real(916920829.L))/real(136949553866517975.L);
      _C4x[324] = (_ep2*(_ep2*(_ep2*(real(349614720.L)*_ep2-real(367698240.L))+
        real(387396360.L))-real(408918380.L))+real(432509825.L))/
        real(50268542125121892.L);
      _C4x[325] = (_ep2*(_ep2*(real(1344672.L)*_ep2-real(1414224.L))+
        real(1489986.L))-real(1572763.L))/real(145917393686856.L);
      _C4x[326] = (_ep2*(real(373520.L)*_ep2-real(392840.L))+real(413885.L))/
        real(31334300028528.L);
      _C4x[327] = (real(1334.L)*_ep2-real(1403.L))/real(88348590306.L);
      _C4x[328] = real(667.L)/real(35525324448.L);
      _C4x[329] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*((real(1542238198824960.L)-real(1466390418554880.L)*
        _ep2)*_ep2-real(1624858102333440.L))+real(1715127996907520.L))-
        real(1814077689036800.L))+real(1922922350379008.L))-
        real(2043104997277696.L))+real(2176350975361024.L))-
        real(2324738541862912.L))+real(2490791294853120.L))-
        real(2677600641967104.L))+real(2888990166332928.L))-
        real(3129739346860672.L))+real(3405892818642496.L))-
        real(3725195270390230.L))+real(4097714797429253.L))/
        real(8571434909745461439457371600.L);
      _C4x[330] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*((real(771119099412480.L)-real(733195209277440.L)*_ep2)*
        _ep2-real(812429051166720.L))+real(857563998453760.L))-
        real(907038844518400.L))+real(961461175189504.L))-
        real(1021552498638848.L))+real(1088175487680512.L))-
        real(1162369270931456.L))+real(1245395647426560.L))-
        real(1338800320983552.L))+real(1444495083166464.L))-
        real(1564869673430336.L))+real(1702946409321248.L))-
        real(1862597635195115.L))/real(552995800628739447706927200.L);
      _C4x[331] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(96389887426560.L)-real(91649401159680.L)*_ep2)*_ep2-
        real(101553631395840.L))+real(107195499806720.L))-
        real(113379855564800.L))+real(120182646898688.L))-
        real(127694062329856.L))+real(136021935960064.L))-
        real(145296158866432.L))+real(155674455928320.L))-
        real(167350040122944.L))+real(180561885395808.L))-
        real(195608709178792.L))+real(212868301165156.L))/
        real(16233778238154283029275325.L);
      _C4x[332] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(654807861707063517147240.L);
      _C4x[333] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(64890869178177465663240.L);
      _C4x[334] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(14887243941012374038800.L);
      _C4x[335] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/
        real(953146715735548337850.L);
      _C4x[336] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(265994432298292559400.L);
      _C4x[337] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/real(39764824222371008880.L);
      _C4x[338] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(12580545683673344160.L);
      _C4x[339] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(521514967754315925.L);
      _C4x[340] = (_ep2*(_ep2*((real(73539648.L)-real(69922944.L)*_ep2)*_ep2-
        real(77479272.L))+real(81783676.L))-real(86501965.L))/
        real(35994758558729256.L);
      _C4x[341] = (_ep2*((real(18384912.L)-real(17480736.L)*_ep2)*_ep2-
        real(19369818.L))+real(20445919.L))/real(6425769670135992.L);
      _C4x[342] = ((real(1178520.L)-real(1120560.L)*_ep2)*_ep2-real(1241655.L))/
        real(302898233609104.L);
      _C4x[343] = (real(63135.L)-real(60030.L)*_ep2)/real(12241188012398.L);
      _C4x[344] = -real(69.L)/real(10854960248.L);
      _C4x[345] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(_ep2*(real(733195209277440.L)*_ep2-
        real(771119099412480.L))+real(812429051166720.L))-
        real(857563998453760.L))+real(907038844518400.L))-
        real(961461175189504.L))+real(1021552498638848.L))-
        real(1088175487680512.L))+real(1162369270931456.L))-
        real(1245395647426560.L))+real(1338800320983552.L))-
        real(1444495083166464.L))+real(1564869673430336.L))-
        real(1702946409321248.L))+real(1862597635195115.L))/
        real(18325136703593745146426104800.L);
      _C4x[346] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(_ep2*(real(22912350289920.L)*_ep2-real(24097471856640.L))+
        real(25388407848960.L))-real(26798874951680.L))+real(28344963891200.L))-
        real(30045661724672.L))+real(31923515582464.L))-real(34005483990016.L))+
        real(36324039716608.L))-real(38918613982080.L))+real(41837510030736.L))-
        real(45140471348952.L))+real(48902177294698.L))-real(53217075291289.L))/
        real(69413396604521761918280700.L);
      _C4x[347] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(7699637270417539977489960.L);
      _C4x[348] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(2864043786240.L)*_ep2-real(3012183982080.L))+
        real(3173550981120.L))-real(3349859368960.L))+real(3543120486400.L))-
        real(3755707715584.L))+real(3990439447808.L))-real(4250685498752.L))+
        real(4540504964576.L))-real(4864826747760.L))+real(5229688753842.L))-
        real(5642558918619.L))/real(589611863049991800078060.L);
      _C4x[349] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(111397652937920178152400.L);
      _C4x[350] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(6113285831959034166900.L);
      _C4x[351] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(1502934058355377170600.L);
      _C4x[352] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(201909322991176933020.L);
      _C4x[353] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(58275401270348939040.L);
      _C4x[354] = (_ep2*(_ep2*(_ep2*(_ep2*(real(699229440.L)*_ep2-
        real(735396480.L))+real(774792720.L))-real(817836760.L))+
        real(865019650.L))-real(916920829.L))/real(2229926069018454300.L);
      _C4x[355] = (_ep2*(_ep2*(_ep2*(real(69922944.L)*_ep2-real(73539648.L))+
        real(77479272.L))-real(81783676.L))+real(86501965.L))/
        real(143414853066284904.L);
      _C4x[356] = (_ep2*(_ep2*(real(17480736.L)*_ep2-real(18384912.L))+
        real(19369818.L))-real(20445919.L))/real(24041241696888108.L);
      _C4x[357] = (_ep2*(real(14567280.L)*_ep2-real(15320760.L))+
        real(16141515.L))/real(13922873979342608.L);
      _C4x[358] = (real(420210.L)*_ep2-real(441945.L))/real(287878973257084.L);
      _C4x[359] = real(667.L)/real(336503767688.L);
      _C4x[360] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*((real(24097471856640.L)-real(22912350289920.L)*_ep2)*_ep2-
        real(25388407848960.L))+real(26798874951680.L))-real(28344963891200.L))+
        real(30045661724672.L))-real(31923515582464.L))+real(34005483990016.L))-
        real(36324039716608.L))+real(38918613982080.L))-real(41837510030736.L))+
        real(45140471348952.L))-real(48902177294698.L))+real(53217075291289.L))/
        real(2438425448462070926742183300.L);
      _C4x[361] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(12048735928320.L)-real(11456175144960.L)*_ep2)*_ep2-
        real(12694203924480.L))+real(13399437475840.L))-real(14172481945600.L))+
        real(15022830862336.L))-real(15961757791232.L))+real(17002741995008.L))-
        real(18162019858304.L))+real(19459306991040.L))-real(20918755015368.L))+
        real(22570235674476.L))-real(24451088647349.L))/
        real(139338597054975481528124760.L);
      _C4x[362] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(7322598944330543323550100.L);
      _C4x[363] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(1067261384598783642298800.L);
      _C4x[364] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/
        real(48156916134335359469580.L);
      _C4x[365] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(10132684457944317053400.L);
      _C4x[366] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/
        real(1197499072302510197220.L);
      _C4x[367] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(310175522890566933600.L);
      _C4x[368] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(10813942549756160100.L);
      _C4x[369] = (_ep2*(_ep2*((real(367698240.L)-real(349614720.L)*_ep2)*_ep2-
        real(387396360.L))+real(408918380.L))-real(432509825.L))/
        real(3206015908868885112.L);
      _C4x[370] = (_ep2*((real(18384912.L)-real(17480736.L)*_ep2)*_ep2-
        real(19369818.L))+real(20445919.L))/real(100042586416082772.L);
      _C4x[371] = ((real(15320760.L)-real(14567280.L)*_ep2)*_ep2-
        real(16141515.L))/real(54344121016143728.L);
      _C4x[372] = (real(127673.L)-real(121394.L)*_ep2)/real(306451810241412.L);
      _C4x[373] = -real(203.L)/real(358213688184.L);
      _C4x[374] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (_ep2*(real(11456175144960.L)*_ep2-real(12048735928320.L))+
        real(12694203924480.L))-real(13399437475840.L))+real(14172481945600.L))-
        real(15022830862336.L))+real(15961757791232.L))-real(17002741995008.L))+
        real(18162019858304.L))-real(19459306991040.L))+real(20918755015368.L))-
        real(22570235674476.L))+real(24451088647349.L))/
        real(5172417617949847420362207000.L);
      _C4x[375] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(2864043786240.L)*_ep2-real(3012183982080.L))+
        real(3173550981120.L))-real(3349859368960.L))+real(3543120486400.L))-
        real(3755707715584.L))+real(3990439447808.L))-real(4250685498752.L))+
        real(4540504964576.L))-real(4864826747760.L))+real(5229688753842.L))-
        real(5642558918619.L))/real(139795070755401281631411000.L);
      _C4x[376] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(13960641343994190068454000.L);
      _C4x[377] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(485217412565651727988950.L);
      _C4x[378] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(83824935061175713805400.L);
      _C4x[379] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(8467165157694516546000.L);
      _C4x[380] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(1926847945229279436000.L);
      _C4x[381] = (_ep2*(_ep2*(_ep2*(_ep2*(real(699229440.L)*_ep2-
        real(735396480.L))+real(774792720.L))-real(817836760.L))+
        real(865019650.L))-real(916920829.L))/real(60213998288414982375.L);
      _C4x[382] = (_ep2*(_ep2*(_ep2*(real(69922944.L)*_ep2-real(73539648.L))+
        real(77479272.L))-real(81783676.L))+real(86501965.L))/
        real(3249194574308196696.L);
      _C4x[383] = (_ep2*(_ep2*(real(17480736.L)*_ep2-real(18384912.L))+
        real(19369818.L))-real(20445919.L))/real(466865403275052936.L);
      _C4x[384] = (_ep2*(real(971152.L)*_ep2-real(1021384.L))+real(1076101.L))/
        real(15719373847644880.L);
      _C4x[385] = (real(2262.L)*_ep2-real(2379.L))/real(23216046230410.L);
      _C4x[386] = real(377.L)/real(2550915658280.L);
      _C4x[387] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(3012183982080.L)-real(2864043786240.L)*_ep2)*_ep2-
        real(3173550981120.L))+real(3349859368960.L))-real(3543120486400.L))+
        real(3755707715584.L))-real(3990439447808.L))+real(4250685498752.L))-
        real(4540504964576.L))+real(4864826747760.L))-real(5229688753842.L))+
        real(5642558918619.L))/real(5467984338975552987240047400.L);
      _C4x[388] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(1506091991040.L)-real(1432021893120.L)*_ep2)*_ep2-
        real(1586775490560.L))+real(1674929684480.L))-real(1771560243200.L))+
        real(1877853857792.L))-real(1995219723904.L))+real(2125342749376.L))-
        real(2270252482288.L))+real(2432413373880.L))-real(2614844376921.L))/
        real(280409453280797589089233200.L);
      _C4x[389] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/
        real(6668273584116528033219570.L);
      _C4x[390] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(886149313503857545942800.L);
      _C4x[391] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/
        real(73398225966986180573040.L);
      _C4x[392] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(14258674794696667826400.L);
      _C4x[393] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(391022331742319334525.L);
      _C4x[394] = (_ep2*(_ep2*((real(367698240.L)-real(349614720.L)*_ep2)*_ep2-
        real(387396360.L))+real(408918380.L))-real(432509825.L))/
        real(94458727981674003948.L);
      _C4x[395] = (_ep2*((real(18384912.L)-real(17480736.L)*_ep2)*_ep2-
        real(19369818.L))+real(20445919.L))/real(2467717131596708376.L);
      _C4x[396] = ((real(222040.L)-real(211120.L)*_ep2)*_ep2-real(233935.L))/
        real(16617623781796016.L);
      _C4x[397] = (real(16653.L)-real(15834.L)*_ep2)/real(734049170812418.L);
      _C4x[398] = -real(2639.L)/real(75507103485088.L);
      _C4x[399] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(1432021893120.L)*_ep2-real(1506091991040.L))+
        real(1586775490560.L))-real(1674929684480.L))+real(1771560243200.L))-
        real(1877853857792.L))+real(1995219723904.L))-real(2125342749376.L))+
        real(2270252482288.L))-real(2432413373880.L))+real(2614844376921.L))/
        real(11527102120002517108235775600.L);
      _C4x[400] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(179002736640.L)*_ep2-real(188261498880.L))+real(198346936320.L))-
        real(209366210560.L))+real(221445030400.L))-real(234731732224.L))+
        real(249402465488.L))-real(265667843672.L))+real(283781560286.L))-
        real(304051671735.L))/real(140574416097591672051655800.L);
      _C4x[401] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(12765340110744758702365200.L);
      _C4x[402] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(812339825229211917423240.L);
      _C4x[403] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(129252960165980064566880.L);
      _C4x[404] = (_ep2*(_ep2*(_ep2*(_ep2*(real(699229440.L)*_ep2-
        real(735396480.L))+real(774792720.L))-real(817836760.L))+
        real(865019650.L))-real(916920829.L))/real(3022496942656846747950.L);
      _C4x[405] = (_ep2*(_ep2*(_ep2*(real(69922944.L)*_ep2-real(73539648.L))+
        real(77479272.L))-real(81783676.L))+real(86501965.L))/
        real(128011635218407626972.L);
      _C4x[406] = (_ep2*(_ep2*(real(58464.L)*_ep2-real(61488.L))+real(64782.L))-
        real(68381.L))/real(50021293208041386.L);
      _C4x[407] = (_ep2*(real(48720.L)*_ep2-real(51240.L))+real(53985.L))/
        real(21108873452551696.L);
      _C4x[408] = (real(2030.L)*_ep2-real(2135.L))/real(476140002689136.L);
      _C4x[409] = real(29.L)/real(3895943914464.L);
      _C4x[410] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        ((real(188261498880.L)-real(179002736640.L)*_ep2)*_ep2-
        real(198346936320.L))+real(209366210560.L))-real(221445030400.L))+
        real(234731732224.L))-real(249402465488.L))+real(265667843672.L))-
        real(283781560286.L))+real(304051671735.L))/
        real(6059117781026964120995728200.L);
      _C4x[411] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(94130749440.L)-
        real(89501368320.L)*_ep2)*_ep2-real(99173468160.L))+
        real(104683105280.L))-real(110722515200.L))+real(117365866112.L))-
        real(124701232744.L))+real(132833921836.L))-real(141890780143.L))/
        real(281819431675672749813754800.L);
      _C4x[412] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/
        real(12240641981872654789890360.L);
      _C4x[413] = (_ep2*(_ep2*(_ep2*(_ep2*((real(11766343680.L)-
        real(11187671040.L)*_ep2)*_ep2-real(12396683520.L))+
        real(13085388160.L))-real(13840314400.L))+real(14670733264.L))-
        real(15587654093.L))/real(1494694488073256644093920.L);
      _C4x[414] = (_ep2*(_ep2*(_ep2*((real(735396480.L)-real(699229440.L)*_ep2)*
        _ep2-real(774792720.L))+real(817836760.L))-real(865019650.L))+
        real(916920829.L))/real(28597471072830165384450.L);
      _C4x[415] = (_ep2*(_ep2*((real(3197376.L)-real(3040128.L)*_ep2)*_ep2-
        real(3368664.L))+real(3555812.L))-real(3760955.L))/
        real(44858778153459082956.L);
      _C4x[416] = (_ep2*((real(799344.L)-real(760032.L)*_ep2)*_ep2-
        real(842166.L))+real(888953.L))/real(4590049143423607182.L);
      _C4x[417] = ((real(666120.L)-real(633360.L)*_ep2)*_ep2-real(701805.L))/
        real(1730927623109239072.L);
      _C4x[418] = (real(35685.L)-real(33930.L)*_ep2)/real(45550726923927344.L);
      _C4x[419] = -real(377.L)/real(266222834155040.L);
      _C4x[420] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*
        (real(89501368320.L)*_ep2-real(94130749440.L))+real(99173468160.L))-
        real(104683105280.L))+real(110722515200.L))-real(117365866112.L))+
        real(124701232744.L))-real(132833921836.L))+real(141890780143.L))/
        real(12709369004105339375747137200.L);
      _C4x[421] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(22375342080.L)*_ep2-
        real(23532687360.L))+real(24793367040.L))-real(26170776320.L))+
        real(27680628800.L))-real(29341466528.L))+real(31175308186.L))-
        real(33208480459.L))/real(282430422313451986127714160.L);
      _C4x[422] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(11187671040.L)*_ep2-
        real(11766343680.L))+real(12396683520.L))-real(13085388160.L))+
        real(13840314400.L))-real(14670733264.L))+real(15587654093.L))/
        real(23514096214810988669282400.L);
      _C4x[423] = (_ep2*(_ep2*(_ep2*(_ep2*(real(30401280.L)*_ep2-
        real(31973760.L))+real(33686640.L))-real(35558120.L))+real(37609550.L))-
        real(39866123.L))/real(14996234830874355018675.L);
      _C4x[424] = (_ep2*(_ep2*(_ep2*(real(15200640.L)*_ep2-real(15986880.L))+
        real(16843320.L))-real(17779060.L))+real(18804775.L))/
        real(2211209528003434308636.L);
      _C4x[425] = (_ep2*(_ep2*(real(760032.L)*_ep2-real(799344.L))+
        real(842166.L))-real(888953.L))/real(38511631837505387088.L);
      _C4x[426] = (_ep2*(real(633360.L)*_ep2-real(666120.L))+real(701805.L))/
        real(12707541818436120992.L);
      _C4x[427] = (real(6786.L)*_ep2-real(7137.L))/real(59715892003929140.L);
      _C4x[428] = real(1131.L)/real(4746558335788640.L);
      _C4x[429] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*((real(23532687360.L)-
        real(22375342080.L)*_ep2)*_ep2-real(24793367040.L))+
        real(26170776320.L))-real(27680628800.L))+real(29341466528.L))-
        real(31175308186.L))+real(33208480459.L))/
        real(13300502446156750509502818000.L);
      _C4x[430] = (_ep2*(_ep2*(_ep2*(_ep2*((real(511580160.L)-real(486420480.L)*
        _ep2)*_ep2-real(538986240.L))+real(568929920.L))-real(601752800.L))+
        real(637857968.L))-real(677724091.L))/
        real(24607775108523127677156000.L);
      _C4x[431] = (_ep2*(_ep2*(_ep2*((real(31973760.L)-real(30401280.L)*_ep2)*
        _ep2-real(33686640.L))+real(35558120.L))-real(37609550.L))+
        real(39866123.L))/real(245868501296893495073625.L);
      _C4x[432] = (_ep2*(_ep2*((real(3197376.L)-real(3040128.L)*_ep2)*_ep2-
        real(3368664.L))+real(3555812.L))-real(3760955.L))/
        real(5553735558706300124016.L);
      _C4x[433] = (_ep2*((real(799344.L)-real(760032.L)*_ep2)*_ep2-
        real(842166.L))+real(888953.L))/real(394968131170694783856.L);
      _C4x[434] = ((real(44408.L)-real(42224.L)*_ep2)*_ep2-real(46787.L))/
        real(7388105708393093600.L);
      _C4x[435] = (real(793.L)-real(754.L)*_ep2)/real(50589875285720700.L);
      _C4x[436] = -real(29.L)/real(827888081823600.L);
      _C4x[437] = (_ep2*(_ep2*(_ep2*(_ep2*(_ep2*(real(486420480.L)*_ep2-
        real(511580160.L))+real(538986240.L))-real(568929920.L))+
        real(601752800.L))-real(637857968.L))+real(677724091.L))/
        real(1207968338105057534196391200.L);
      _C4x[438] = (_ep2*(_ep2*(_ep2*(_ep2*(real(30401280.L)*_ep2-
        real(31973760.L))+real(33686640.L))-real(35558120.L))+real(37609550.L))-
        real(39866123.L))/real(6163103765842130276512200.L);
      _C4x[439] = (_ep2*(_ep2*(_ep2*(real(15200640.L)*_ep2-real(15986880.L))+
        real(16843320.L))-real(17779060.L))+real(18804775.L))/
        real(473713073766689229096624.L);
      _C4x[440] = (_ep2*(_ep2*(real(760032.L)*_ep2-real(799344.L))+
        real(842166.L))-real(888953.L))/real(5156528379172959678120.L);
      _C4x[441] = (_ep2*(real(126672.L)*_ep2-real(133224.L))+real(140361.L))/
        real(236123858440243271456.L);
      _C4x[442] = (real(522.L)*_ep2-real(549.L))/real(317029885123849720.L);
      _C4x[443] = real(87.L)/real(19640659160215120.L);
      _C4x[444] = (_ep2*(_ep2*(_ep2*((real(31973760.L)-real(30401280.L)*_ep2)*
        _ep2-real(33686640.L))+real(35558120.L))-real(37609550.L))+
        real(39866123.L))/real(314842811527382016891612600.L);
      _C4x[445] = (_ep2*(_ep2*((real(3197376.L)-real(3040128.L)*_ep2)*_ep2-
        real(3368664.L))+real(3555812.L))-real(3760955.L))/
        real(2469355384528486406993040.L);
      _C4x[446] = (_ep2*((real(799344.L)-real(760032.L)*_ep2)*_ep2-
        real(842166.L))+real(888953.L))/real(91391237018107987486680.L);
      _C4x[447] = ((real(10248.L)-real(9744.L)*_ep2)*_ep2-real(10797.L))/
        real(246171682203657878752.L);
      _C4x[448] = (real(549.L)-real(522.L)*_ep2)/real(3503517496453777544.L);
      _C4x[449] = -real(29.L)/real(61429295671311120.L);
      _C4x[450] = (_ep2*(_ep2*(_ep2*(real(3040128.L)*_ep2-real(3197376.L))+
        real(3368664.L))-real(3555812.L))+real(3760955.L))/
        real(131077415411399860093650960.L);
      _C4x[451] = (_ep2*(_ep2*(real(58464.L)*_ep2-real(61488.L))+real(64782.L))-
        real(68381.L))/real(190242983180551320890640.L);
      _C4x[452] = (_ep2*(real(9744.L)*_ep2-real(10248.L))+real(10797.L))/
        real(4526544605418280586848.L);
      _C4x[453] = (real(58.L)*_ep2-real(61.L))/real(5469777315892122084.L);
      _C4x[454] = real(29.L)/real(703302752481745680.L);
      _C4x[455] = (_ep2*((real(61488.L)-real(58464.L)*_ep2)*_ep2-real(64782.L))+
        real(68381.L))/real(10478285093218993340819760.L);
      _C4x[456] = ((real(3416.L)-real(3248.L)*_ep2)*_ep2-real(3599.L))/
        real(42336505427147447841696.L);
      _C4x[457] = (real(61.L)-real(58.L)*_ep2)/real(104211770430232260620.L);
      _C4x[458] = -real(29.L)/real(10232365536106966560.L);
      _C4x[459] = (_ep2*(real(3248.L)*_ep2-real(3416.L))+real(3599.L))/
        real(2416376017304170372096800.L);
      _C4x[460] = (real(58.L)*_ep2-real(61.L))/real(3028040121935050591600.L);
      _C4x[461] = real(29.L)/real(201751358211920378400.L);
      _C4x[462] = (real(61.L)-real(58.L)*_ep2)/
        real(178874588293945079492880.L);
      _C4x[463] = -1/real(209087771237808392160.L);
      _C4x[464] = 1/real(12769026871558087949280.L);
      break;
    default:
      STATIC_ASSERT(nC4_ == 30, "Bad value of nC4_");
    }
  }

} // namespace GeographicLib

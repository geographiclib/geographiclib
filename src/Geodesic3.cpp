/**
 * \file Geodesic3.cpp
 * \brief Implementation for GeographicLib::Triaxial::Geodesic3 class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <iostream>
#include <iomanip>
#include <GeographicLib/Triaxial/Geodesic3.hpp>

namespace GeographicLib {
  namespace Triaxial {

  using namespace std;

  Geodesic3::Geodesic3(const Ellipsoid3& t)
    : _t(t)
    , _umbalt(false)
    , _ellipthresh(1/real(8))
  {}

  Geodesic3::Geodesic3(real a, real b, real c)
    : Geodesic3(Ellipsoid3(a, b, c))
  {}

  Geodesic3::Geodesic3(real b, real e2, real k2, real kp2)
    : Geodesic3(Ellipsoid3(b, e2, k2, kp2))
  {}

  GeodesicLine3
  Geodesic3::Inverse(Angle bet1, Angle omg1, Angle bet2, Angle omg2,
                     real& s12, Angle& alp1, Angle& alp2) const {
    using TL = GeodesicLine3;
    string msg;
    bet1.round();
    omg1.round();
    bet2.round();
    omg2.round();

    if (!_umbline)
      // Initialize _umbline
      _umbline = make_shared<GeodesicLine3>(GeodesicLine3(*this));

    // In triaxial + oblate cases, [bet, omg] are initially put into [-90,90] x
    // [-180,180].  For prolate case, maybe we instead put [bet, omg] into
    // [-180,180] x [0,180].
    // What about rotating coordinates to make
    //   oblate: omg1 = 0
    //   prolate: bet1 = -90
    // this eliminates many of the special cases

    ang bet10(bet1), omg10(omg1);
    // If prolate put omg in [ 0, 180], bet in [-180, 180]
    //       else put bet in [-90, 90], omg in [-180, 180]
    bool flip1 = t().AngNorm(bet1, omg1, prolate()),
      flip2 = t().AngNorm(bet2, omg2, prolate());
    bool swap12;
    {
      // For prolate, swap based on omg, switch 1 & 2 because poles are at
      // 0/180, instead of +/-90.
      ang tmp1(prolate() ? omg2 : bet1), tmp2(prolate() ? omg1 : bet2);
      tmp1.setquadrant(0U); tmp2.setquadrant(0U);
      ang tmp12 = tmp2 - tmp1; // |bet2| - |bet1|
      swap12 = tmp12.s() > 0; // is |bet2| > |bet1|
      if (!biaxial() && tmp12.s() == 0) {
        // don't need to do this if biaxial
        tmp1 = omg1; tmp2 = omg2;
        tmp1.setquadrant(0U); tmp2.setquadrant(0U);
        tmp12 = tmp2 - tmp1;
        swap12 = tmp12.s() < 0; // is |omg2| < |omg1|
      }
      // N.B. No swapping if bet1 = +0 and bet2 = -0.
    }
    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
    }
    if (oblate()) {
      // Rotate, subtracting omg1 from omg[12], so omg1 = 0
      omg2 -= omg1;
      omg2 = omg2.base();
      omg1 = ang::cardinal(0);
    } else if (prolate()) {
      // Rotate, subtracting bet1 + 90 from bet[12], so bet1 = -90
      bet2 -= bet1 + ang::cardinal(1);
      bet2 = bet2.base();
      bet1 = ang::cardinal(-1);
    }

    bool swapomg = swapomg_ &&
      // Don't want the swap to make a new umbilical point
      !(bet1.c() == 0 && omg2.s() == 0) &&
      fabs(omg2.s()) < fabs(omg1.s());
    if (swapomg) swap(omg1, omg2);
    // Now |bet1| >= |bet2|
    bool flipz = bet1.s() > 0;
    if (flipz) {                // Not needed for prolate, bet1 already -90
      bet1.reflect(true);
      bet2.reflect(true);
    }
    // Now bet1 <= 0
    bool flipy = prolate() ? signbit(bet2.c()) :
      signbit(omg1.s()) || (omg1.s() == 0 && signbit(omg2.s()));
    if (flipy) {
      if (prolate())
        bet2.reflect(false, true);
      else {
        omg1.reflect(true);
        omg2.reflect(true);
      }
    }
    // For oblate omg2 >= 0, for prolate bet2 in [-90, 90]
    bool flipx = signbit(omg1.c());
    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
    }
    bool flipomg = bet2.c() == 0 && signbit(omg2.s());
    // Eliminate coordinate ambiguity bet2 = +/90 and omg2 < 0 for point 2
    // Point 1 is already done with flipy.  (Maybe skip for oblate?)
    if (flipomg) omg2.reflect(true);

    bool umb1 = bet1.c() == 0 && omg1.s() == 0,
      umb2 = bet2.c() == 0 && omg2.s() == 0;

    // Now bet1 <= 0, bet1 <= bet2 <= -bet1, 0 <= omg1 <= 90
    //
    // Oblate: omg1 = 0, omg2 >= 0; rotation results in meridional geodesics
    //   following triaxial middle ellipse (omg2 = 0/180); minor ellipse folded
    //   into middle ellipse.
    //
    // Prolate: bet1 = -90, bet2 in [-90,90], omg2 >= 0; rotation results in
    //   meridional geodesics following triaxial middle ellipse (bet2 = +/-90);
    //   major ellipse folded into middle ellipse.
    //
    // Distinguish two general categories of solution
    //
    // A The geodesic follows one of the principal ellipses (this means, of
    //   course, that both points need to be on the same principal ellipse) and
    //   the initial and final azimuth are easily determined.  The subcases
    //   are:
    //
    //   a Minor ellipse, omg = +/-90: the ellipse is followed regardless of
    //     the points; however we treat this using B.d (and the solution along
    //     the ellipse is found on the first iteration).
    //     Oblate: this becomes the middle ellipse case, A.c.
    //     Prolate: same as triaxial case.
    //
    //   b Major ellipse, bet = 0: the ellipse is followed provided that both
    //     points are close enough "equatorial geodesics".  From point 1 we
    //     follow an equatorial geodesic +/- 180deg in arc distance (= psi
    //     variable).
    //     Oblate: same as triaxial case, except conjugate point is trivially
    //       found (but don't use this)
    //     Prolate: treated by middle ellipse cases A.c.2 and A.c.3.
    //
    //     1 If point 2 is within the resulting span of longitudes, the
    //       geodesic is equatorial (compute with ArcPos0 with the thr
    //       variable).
    //
    //     2 Otherwise, do search with strategy B.
    //
    //   c Middle ellipse, bet = +/-90 or omg = 0,180: the ellipse is followed
    //     provided that both points are close enough "meridional geodesics".
    //     The four umbilical points divide the middle ellipse into 4 segments.
    //     Oblate: same as triaxial case, but geodesic always follows ellipse.
    //     Prolate: this becomes the major ellipse case, A.b
    //
    //     There are several subcases:
    //
    //     1 opposite umbilical points: multiple shortest geodesic, two of
    //       which follow the middle ellipse.  But it's more useful to return
    //       the one going through bet = 0, omg = 90.  Then all the others can
    //       be generated by multipling tan(alp1) and tan(alp2) by a constant.
    //       Oblate/Prolate: opposite poles with tht12 = 180
    //
    //     2 bet1 = -90, bet2 = -90: geodesic follows ellipse (compute with
    //       ArcPos0).  Points may be adjacent umbilical points.
    //       Oblate: same pole
    //       Prolate: meridional geodesic not crossing a pole
    //
    //     3 bet1 = -90, bet2 = 90: geodesic may follow the ellipse if if they
    //       are close enough.  See case A.b for strategy.  Points may be
    //       adjacent umbilical points.
    //       Oblate: opposite poles with tht12 < 180
    //       Prolate: meridional geodesic crossing a pole + opposite meridians,
    //       non-meridional
    //
    //     4 |bet2| != 90: geodesic follows the ellipse.  Compute with Hybrid.
    //       Oblate: meridional geodesic
    //       Prolate: same (omg2 = 0) or opposite (omg2 = 180) poles
    //
    // B The geodesic does not follow a principal ellipse
    //
    //   a point 1 is an umbilical point, gam = 0, so we know alp2
    //     Oblate: remaining meridional cases (one end a pole)
    //     Prolate: remaining meridional cases (one end a pole)
    //
    //   b bet1 = -90: do a search with alp1 bracket = [-90+eps, 90-eps].
    //     Oblate: already treated by B.a
    //     Prolate: this treate the general case
    //
    //   c omg1 = 0: do a search with alp1 bracket = [eps, 180-eps], or
    //     [-180+eps, -eps], depending on the sign of omg2.
    //     Oblate: this treats the general case
    //     Prolate: already handled by B.a
    //
    //   d general case, bracket search by two neighboring umbilical
    //     directions.  This treats case A.a.
    //     Oblate/Prolate: already handled by B.b or B.c

    // Set up variables for search, "h" variables are for hybridalt
    real fa = Math::NaN(), fb = Math::NaN();
    ang alpa, alpb;

    // and for the final result
    ang bet2a, omg2a;

    // Much of the initial logic for the inverse solution uses umbilical
    // geodesics; use the cached value for this.
    TL::fline lf = _umbline->_f;
    TL::fline::fics fic;
    TL::fline::disttx d{Math::NaN(), Math::NaN(), 0};

    real aa = k2() * Math::sq(bet2.c()), bb = kp2() * Math::sq(omg2.s());

    if constexpr (debug_)
      cout << "COORDS " << real(bet1) << " " << real(omg1) << " "
           << real(bet2) << " " << real(omg2) << "\n";
    //    bet1.setn(0); omg1.setn(0); bet2.setn(0); omg2.setn(0);
    // flag for progress
    bool done = false, backside = false;
    if (bet1.c() * omg1.s() == 0 && bet2.c() * omg2.s() == 0) {
      // Case A.c, both points on middle ellipse
      if (umb1 && umb2 && bet2.s() > 0 && omg2.c() < 0) {
        // Case A.c.1, process opposite umbilical points
        // For oblate/prolate this gives 0/90
        alp1 = biaxial() ? ang(k(), kp(), 0, true) :
          ang(exp(lf.deltashift()/2), 1);
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        bool betp = k2() < kp2();
        //        if (biaxial()) betp = !betp;
        //        betp = !betp;
        d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2, betp);
        if constexpr (debug_) msg = "A.c opposite umbilics";
        backside = signbit(bet2a.c());
        done = true;
      } else if (bet1.c() == 0 && bet2.c() == 0) {
        // Case A.c.{2,3}, bet1 = -90, bet2 = +/-90
        if (bet2.s() < 0) {
          // Case A.c.2, bet1 = bet2 = -90
          // If oblate, bet1 = -90, omg1 = 0, need alp1 = omg2 to follow omg2
          // meridian.
          alp1 = ang::cardinal(oblate() ? 2 : 1);
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          ang omg12 = omg2 - omg1;
          if (omg12.s() == 0 && omg12.c() < 0) {
            // adjacent E/W umbilical points
            // Should be able to get ArcPos0 to return this?
            d = oblate() ?
              TL::fline::disttx{ (biaxial() ? 1 : -1 ) * Math::pi()/2,
                                 -Math::pi()/2, 0 } :
              prolate() ?
              TL::fline::disttx{ Math::pi()/2, -Math::pi()/2, 0 } :
              TL::fline::disttx{ -BigValue(), BigValue(), 0 };
            if constexpr (debug_) msg = "A.c.2 adjacent EW umbilics";
            alp2 = ang::cardinal(prolate() ? 1 : 0);
          } else {
            d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
            if constexpr (debug_) msg = "A.c.2 bet1/2 = -90";
          }
          // not needed
          // if (omg2a.s() < 0) alp2.reflect(true);
          done = true;
        } else {
          // Case A.c.3, bet1 = -90, bet2 = 90
          // need to see how far apart the points are
          // If point 1 is at [-90, 0], direction is 0 else -90.
          // For prolate use -90.
          alp1 = ang::cardinal(omg1.s() == 0 && !prolate() ? 0 : -1);
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          // If point 1 is [-90, 0] and point 2 is [90, 0]
          if (omg1.s() == 0 && omg2.s() == 0) {
            // adjacent N/S umbilical points
            // Should be able to get ArcPos0 to return this?
            d = biaxial() ?
              TL::fline::disttx{ Math::pi()/2, -Math::pi()/2, 0 } :
              TL::fline::disttx{ BigValue(), -BigValue(), 0 };
            alp2 = ang::cardinal(oblate() ? 0 :  1);
            if constexpr (debug_) msg = "A.c.3 adjacent NS umbilics";
            done = true;
          } else {
            // FIX ME for oblate
            if (omg1.s() == 0)
              omg2a = ang::cardinal(2);
            else
              // Compute conjugate point along the middle ellipse
              d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
            // XXX FIX HERE for prolate case -90 -1 90 177
            omg2a -= omg2;
            if (omg2a.s() >= -numeric_limits<real>::epsilon()/2) {
              // Includes all cases where point 1 = [-90, 0], omg1.s() == 0.
              ang omg12 = omg2 + omg1;
              // FIX ME for oblate
              d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
              if (!biaxial() && signbit(omg2a.s()))
                alp2.reflect(true);
              if constexpr (debug_) msg = "A.c.3 bet1/2 = -/+90 meridional";
              done = true;
            } else {
              alpa = ang::cardinal(-1) + ang::eps();
              fa = omg2a.radians0();
              (void) lf.ArcPos0(fic, ang::cardinal(-2), bet2a, omg2a, alp2);
              omg2a -= omg2;
              alpb = -alpa;
              fb = omg2a.radians0();
              if constexpr (debug_)
                msg = "A.c.3 general bet1/2 = -/+90, non-meridional";
              done = false;     // A marker
            }
            if constexpr (debug_)
              cout << "ALP/F "
                   << real(alpa) << " " << fa << " "
                   << real(alpb) << " " << fb << "\n";
          }
        }
      } else {
        // Case A.c.4, other meridional cases, invoke Hybrid with the following
        // value of alp1
        // If oblate, bet1 = -90, omg1 = 0, need alp1 = omg2 to follow omg2
        // meridian.
        // If prolate, bet1 = -90, omg1 = 0, need alp1 = -bet2 to follow bet2
        // meridian.
        alp1 = oblate() ? omg2 :
          (prolate() ? -bet2 :
           ang::cardinal(bet1.c() == 0 ?
                         // TODO: CHECK omg2.c() < 1 test; CHANGE TO < 0
                         (omg2.c() < 0 ? 1 :
                          (omg1.s() == 0 && !prolate() ? 0 : -1)) :
                         (omg2.c() > 0 ? 0 : 2)));
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        d = prolate() ?
          lf.ArcPos0(fic, (omg2-omg1).base(), bet2a, omg2a, alp2, false) :
          lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
        if (prolate() && signbit(bet2a.c()))
          alp2.reflect(true,true);
        if constexpr (debug_) msg = "A.c.4 other meridional";
        backside = signbit(bet2a.c());
        done = true;
      }
    } else if (bet1.s() == 0 && bet2.s() == 0) {
      // Case A.b, both points on equator
      ang omg12 = (omg2 - omg1).base();
      int E = signbit(omg12.s()) ? -1 : 1; // Fix #1 for triaxial sphere
      // set direction for probe as +/-90 based on sign of omg12
      alp1 = ang::cardinal(E);
      bet1.reflect(true);
      lf = TL::fline(this->t(), gamma(bet1, omg1, alp1));
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      (void) lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
      omg2a -= omg2;
      if (E * omg2a.s() >= -numeric_limits<real>::epsilon()/2) {
        // geodesic follows the equator
        d = lf.ArcPos0(fic, omg12.flipsign(E), bet2a, omg2a, alp2, false);
        if constexpr (debug_) msg = "A.b.1 bet1/2 = 0 equatorial";
        done = true;
      } else {
        // geodesic does not follow the equator
        alpb = ang::cardinal(-1) - ang::eps();
        alpa = -alpb;
        (E > 0 ? fa : fb) = omg2a.radians0();
        (void) lf.ArcPos0(fic, ang::cardinal(-2), bet2a, omg2a, alp2);
        omg2a -= omg2;
        (E > 0 ? fb : fa) = omg2a.radians0();
        // Fix for
        //
        //   echo 0 20 0 -50 | Geod3Solve -i -t 1.5 1 0.5
        //
        // With an ellipsoid this eccentric (a > 2*c?), the E and W conjugate
        // points from a given point on the equator can be less than 180deg
        // apart.  So the assumption that relevant angle differences always lie
        // in [-180,180], and so can be computed with radians0(), fails.
        //
        // In this case the E/W conjugate points are at omega = 104.3 and
        // -40.6.  So that fa, the differenc in omega, should be -205.7d =
        // 104.3-360-(-50).  This is reduced to [-180,180] by radians0() and
        // becomes +154.3d.  Fix by checking signs on fa and fb.
        //
        // findroot which calls HybridA also reduces the differences in omega
        // to [-180,180].  But the first iteration in HybridA uses the
        // bisecting alpha which will result (??? TO CHECK) is a omega
        // difference "in range".
        //
        // If this turns out to be false, we'll need to rethink.  Perhaps use
        // omg1 + 180 instead of omg2 as the base value.  This would
        // necessitate a change in the signature for findroot.  So hold off on
        // this for now.
        //
        // The Octave routine triaxial.distance also failed for equatorial
        // geodesics with this ellipsoid.  But here the problem was accepting
        // the equatorial geodesic becase m12 >= 0.  In fact there may be two
        // intervening conjugate points with an ellipsoid this eccentric.
        // Failure case is t.distance([0,20],[0,-175]).
        if (fa > 0) fa -= 2*Math::pi();
        if (fb < 0) fb += 2*Math::pi();
        if constexpr (debug_) msg = "A.b.2 general bet1/2 = 0 non-equatorial";
        done = false;           // A marker
      }
    } else if (umb1) {
      // Case B.a, umbilical point to general point
      alp2 = ang(kp() * omg2.s(), k() * bet2.c());
      // RETHINK THIS.  If we know alp2, we can compute delta.  This should be
      // enough to find alp1.
      fic = TL::fline::fics(lf, bet2, omg2, alp2);
      //      bool betp = k2() > kp2();   // This could be betb = !prolate();
      bool betp = aa > bb;
      real delta = (lf.transpolar() ? -1 : 1) * fic.delta;
      alp1 = oblate() ?
        // For oblate
        // delta =
        //   atan2(bet1.s() * fabs(alp1.s()), bet0.c() * alp1.c())) - omg1
        // Here omg1 = -90 (because of pi/2 shift in ext. vs int. omg)
        // bet1.s() = -1, bet0.c() = 1, alp1.s() > 0, so
        // alp1 = 90 - delta
        ang::cardinal(1) - ang::radians(delta) :
        (prolate() ?
         // For prolate
         // delta =
         //   bet1 - atan2(omg1.s() * fabs(alp1.c()), omg0.c() * alp1.s())
         // Here bet1 = -90, omg1.s() = -1, alp1.c() > 0, omg0.c() = 1
         // delta = bet1 + atan2(alp1.c(), alp1.s()) = -90 + 90 - alp1
         // alp1 = -delta
         -ang::radians(delta) :
         // For triaxial case at an umbilic point
         // delta = f.deltashift()/2 - log(fabs(alp1.t()));
         // alp1.t() = exp(f.deltashift()/2 - delta)
         ang(exp(lf.deltashift()/2 - delta), 1));
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      d = lf.ArcPos0(fic, (betp ? bet2 - bet1 : omg2 - omg1).base(),
                     bet2a, omg2a, alp2, betp);
      if constexpr (debug_) msg = "B.a umbilic to general";
      done = true;
    } else if (bet1.c() == 0) {
      // Case B.b, bet1 = -90 to general point
      // subsumed by B.a for oblate
      // general handling for prolate
      if (!signbit(omg2.s())) {
        alpa = ang::cardinal(-1) + ang::eps();
        alpb = -alpa;
        fa = -omg2.radians();
        fb = (ang::cardinal(2) - omg2).radians0();
      } else {
        alpa = ang::cardinal(1) + ang::eps();
        alpb = -alpa;
        fa = (ang::cardinal(2) - omg2).radians0();
        fb = -omg2.radians();
      }
      if constexpr (debug_) msg = "B.b general bet1 = -90";
      done = false;             // A marker
    } else if (omg1.s() == 0) {
      // Case B.c, omg1 = 0 to general point
      // subsumed by B.a for prolate
      // general handling of oblate
      if (omg2.s() > 0) {
        alpa = ang::eps();
        alpb = ang::cardinal(2) - alpa;
        fa = -omg2.radians0();
        fb = (ang::cardinal(2)-omg2).radians0();
      } else {
        // alpb = -ang::eps();
        // alpa = ang::cardinal(-2) - alpb;
        alpa = ang(-numeric_limits<real>::epsilon()/(1<<20), -1, 0, true);
        alpb = ang(-numeric_limits<real>::epsilon()/(1<<20),  1, 0, true);
        fa = (ang::cardinal(2)-omg2).radians0();
        fb = -omg2.radians0();
      }
      if constexpr (debug_) msg = "B.c general omg1 = 0";
      done = false;             // A marker
    } else {
      // Case B.d, general case
      real f[4];
      alpa = ang( kp() * fabs(omg1.s()), k() * fabs(bet1.c()));
      alpb = alpa;

      fic = TL::fline::fics(lf, bet1, omg1, alpb);
      unsigned qb = 0U, qa = 3U; // qa = qb - 1 (mod 4)
      if constexpr (debug_) msg = "B.d general";
      for (; !done && qb <= 4U; ++qb, ++qa) {
        if (qb) {
          alpb.setquadrant(qb);
          fic.setquadrant(lf, qb);
        }
        if (qb < 4U) {
          f[qb] = lf.Hybrid0(fic, bet2, omg2);
          if constexpr (debug_)
            cout << "f[qb] " << qb << " " << f[qb] << "\n";
          if (fabs(f[qb]) < 2*numeric_limits<real>::epsilon()) {
            alp1 = alpb;
            d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
            if constexpr (debug_) msg = "B.d accidental umbilic";
            backside = signbit(bet2a.c()); // qb == 1U || qb == 2U;
            done = true;
            break;
          }
        }
        if (qb && (f[qa & 3U] < 0 && f[qb & 3U] > 0) &&
            // Fix #2 for triaxial sphere
            // Expect f[qb] - f[qa] <= pi.  The following condition catches
            // cases where f[qb] = pi and f[qa] = -pi.  This can happen with e2
            // == 0 and bet2 == - bet1.  Here "4" is a standin for pi+eps.
            f[qb & 3U] - f[qa & 3U] < 4) {
          break;
        }
      }
      if constexpr (debug_)
        cout << "fDD " << done << " " << qa << " " << qb << " "
             << f[qa & 3U] << " " << f[qb & 3U] << "\n";
      if (!done) {
        fa = f[qa & 3U]; fb = f[qb & 3U];
        alpa.setquadrant(qa);
        done = false;           // A marker
      }
    }

    int countn = 0, countb = 0;
    if (!done) {
      // Iterative search for the solution
      if constexpr (debug_)
        cout << "X " << done << " " << msg << "\n";
      alp1 = findroot(
                      [this, &bet1, &omg1, &bet2, &omg2]
                      (const ang& alp) -> real
                      {
                        return HybridA(bet1, omg1, alp, bet2, omg2, true);
                      },
                      alpa, alpb,
                      fa, fb,
                      &countn, &countb);
      if constexpr (debug_)
        cout << "ALP1 " << real(alp1) << "\n";
      lf = TL::fline(this->t(), gamma(bet1, omg1, alp1));
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      // Let aa = k2() * Math::sq(bet2.c()), bb = kp2() * Math::sq(omg2.s());
      // Figure sin/cos alp2 and call lf.Hybrid with betp = |salp2| <
      // |calp2|.  For !transpolar
      //   gammax = aa * salp2^2 - bb * calp2^2
      //   salp2^2 = (bb + gammax) / (aa + bb)
      // and use betp = salp2^2 < 1/2
      // For transpolar
      //   gammax = - aa * salp2^2 + bb * calp2^2
      //   salp2^2 = (bb + gammax) / (aa + bb)
      // and use betp = salp2^2 < 1/2
      // For transpolar
      //   gammax = bb * calp2^2 - aa * salp2^2
      //   calp2^2 = (aa + gammax) / (aa + bb)
      // and use betp = calp2^2 > 1/2
      bool betp = true;
      if (hybridalt_ && (swapomg_ || omg1.s() <= fabs(omg2.s())))
        betp = (lf.transpolar() ?
                2 * (aa + lf.gammax()) > (aa + bb) :
                2 * (bb + lf.gammax()) < (aa + bb));
      d = lf.Hybrid(fic, betp ? bet2 : omg2, bet2a, omg2a, alp2, betp);
      // Don't set backside for !betp and bet2.c() == 0.  Not sure why.
      backside = (betp || bet2.c() != 0) && signbit(bet2a.c());
    }

    if (!biaxial() && backside) alp2.reflect(true, true);
    alp2.round();

    TL::gline lg(this->t(), lf.gm());
    TL::gline::gics gic(lg, fic);
    s12 = lg.dist(gic, d);

    if constexpr (debug_)
      cout << "FLIPS "
           << flip1 << flip2 << swap12 << swapomg
           << flipz << flipy << flipx << flipomg
           << lf.transpolar() << "\n";
    if constexpr (debug_)
      cout << "A "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    // Undo switches in reverse order flipomg flipx flipy flipz swapomg swap12
    // flip2 flip1
    if (flipomg) {
      omg2.reflect(true);
      alp2.reflect(true, true);
    }

    if constexpr (debug_)
      cout << "B "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
      alp1.reflect(true);
      alp2.reflect(true);
    }

    if constexpr (debug_)
      cout << "C "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipy) {
      if (prolate()) {
        bet2.reflect(false, true);
        alp1.reflect(false, true);
        alp2.reflect(false, true);
      } else {
        omg1.reflect(true);
        // This was omg2.reflect(true, true); (by mistake?)
        omg2.reflect(true);
        alp1.reflect(true);
        alp2.reflect(true);
      }
    }

    if constexpr (debug_)
      cout << "D "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipz) {
      bet1.reflect(true);
      bet2.reflect(true);
      alp1.reflect(false, true);
      alp2.reflect(false, true);
    }

    if constexpr (debug_)
      cout << "E "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (swapomg) {
      swap(omg1, omg2);
      if (lf.transpolar()) {
        real
          calp1 = copysign(hypot(k()/kp()*bet1.c(), lf.nu()), alp1.c()),
          calp2 = copysign(hypot(k()/kp()*bet2.c(), lf.nu()), alp2.c()),
          salp1 = (lf.nu() < lf.nup() ?
                   (omg1.s() - lf.nu()) * (omg1.s() + lf.nu()) :
                   (lf.nup() - omg1.c()) * (lf.nup() + omg1.c())),
          salp2 = (lf.nu() < lf.nup() ?
                   (omg2.s() - lf.nu()) * (omg2.s() + lf.nu()) :
                   (lf.nup() - omg2.c()) * (lf.nup() + omg2.c()));
        salp1 = -copysign(signbit(salp1) ? 0 : sqrt(salp1), alp2.s());
        salp2 = -copysign(signbit(salp2) ? 0 : sqrt(salp2), alp1.s());
        alp1 = ang(salp1, calp1);
        alp2 = ang(salp2, calp2);
      } else {
        real
          salp1 = -copysign(hypot(kp()/k()*omg1.s(), lf.nu()), alp1.s()),
          salp2 = -copysign(hypot(kp()/k()*omg2.s(), lf.nu()), alp2.s()),
          calp1 = (lf.nu() < lf.nup() ?
                   (bet1.c() - lf.nu()) * (bet1.c() + lf.nu()) :
                   (lf.nup() - bet1.s()) * (lf.nup() + bet1.s())),
          calp2 = (lf.nu() < lf.nup() ?
                   (bet2.c() - lf.nu()) * (bet2.c() + lf.nu()) :
                   (lf.nup() - bet2.s()) * (lf.nup() + bet2.s()));
        calp1 = copysign(signbit(calp1) ? 0 : sqrt(calp1), alp1.c());
        calp2 = copysign(signbit(calp2) ? 0 : sqrt(calp2), alp2.c());
        alp1 = ang(salp1, calp1);
        alp2 = ang(salp2, calp2);
      }
    }
    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
      swap(alp1, alp2);
      swap(umb1, umb2);
      // points not swapped if umb1 == true
      alp1 += ang::cardinal(2);
      if (umb2 && !biaxial())
        alp2 += ang::cardinal((signbit(alp2.s()) ? -1 : 1) * bet2.s());
      else
        alp2 += ang::cardinal(2);
    }

    if constexpr (debug_)
      cout << "F "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flip1) t().Flip(bet1, omg1, alp1);
    if (flip2) t().Flip(bet2, omg2, alp2);
    alp1.setn(); alp2.setn();
    if constexpr (debug_)
      cout << "G "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";

    fic = TL::fline::fics(lf, bet10, omg10, alp1);
    gic = TL::gline::gics(lg, fic);
    gic.s13 = signbit(s12) ? 0 : s12;

    // clang needs std::move instead of move.
    return TL(std::move(lf), std::move(fic), std::move(lg), std::move(gic));
  }

  Math::real Geodesic3::HybridA(ang bet1, ang omg1, ang alp1,
                                ang bet2a, ang omg2b,
                                bool betp) const {
    ang b1{bet1}, o1{omg1}, a1{alp1};
    // a1 -= ang(1e-8);
    gamblk gam = gamma(b1, o1, a1);
    GeodesicLine3::fline l(this->t(), gam);
    GeodesicLine3::fline::fics ic(l, b1, o1, a1);
    real dang = l.Hybrid0(ic, bet2a, omg2b, betp);
    return dang;
  }

  // Solve f(alp1) = 0 where alp1 is an azimuth and f(alp1) is the difference
  // in lontitude on bet2 and the target longitude.
  Angle Geodesic3::findroot(const function<real(const ang&)>& f,
                            ang xa,  ang xb,
                            real fa, real fb,
                            int* countn, int* countb) {
    // Implement root finding method of Chandrupatla (1997)
    // https://doi.org/10.1016/s0965-9978(96)00051-8
    // Here we follow Scherer (2013), Section 6.1.7.3
    // https://doi.org/10.1007/978-3-319-00401-3

    // Here the independent variable is an ang, but the computations on this
    // variable essentially involve its conversion to radians.  There's no need
    // to worry about the angle wrapping around because (xb-xa).radians() is in
    // (0,pi).

    // require xa and xb to be normalized (the result is normalized)
    // require fa and fb to have opposite signs

    ang xm;                  // The return value
    int cntn = 0, cntb = 0;
    bool trip = false;
    const bool debug = debug_;
    // 25 iterations with line 498534 of testset.txt
    //   43 -2 -43 141
    // This is a near conjugate case m12 = 0.0003857
    // 27 iterations with line 360115 of testpro.txt
    //   -75 29 75 -169
    // Converge failures with line 40045 of testsph[bc]
    // echo 80 -90 -80 90 | ./Geod3Solve -e 1 0 2 1 -i
    // echo 80 -90 -80 90 | ./Geod3Solve -e 1 0 1 2 -i
    //
    //  Offset for debugging output
    real x0 = real(0);
    // If fa and fb have the same signs, assume that the root is at one of the
    // endpoints if corresponding f is small.  Otherwise, it's an error.
    if (fa * fb >= 0) {
      if constexpr (debug)
        cout << "FA FB " << fa/numeric_limits<real>::epsilon() << " "
             << fb/numeric_limits<real>::epsilon() << " " << (fa == fb) << "\n";
      if (fa == fb && fabs(fa) <= 512*numeric_limits<real>::epsilon())
        // If both fa and fb have the same sign and are small (but not too
        // small!), return mean of the endpoints.  This is the case of
        // antipodal points on a triaxial sphere, case A.c.3 general bet1/2 =
        // -/+90, non-meridional.  The mean angle corresponds to the "minor"
        // ellipse where the call to Hybrid (to compute alp2) gives a well
        // defined result.
        return ang(xa.s() + xb.s(), xa.c() + xb.c());
      else if (fmin(fabs(fa), fabs(fb)) > 2*numeric_limits<real>::epsilon())
        // neither endpoint small enough
        throw GeographicLib::GeographicErr
          ("Bad inputs Geodesic3::findroot");
      else
        // return best endpoint
        return fabs(fa) < fabs(fb) ? xa : xb;
    }
    // tp = 1 - t
    for (real t = 1/real(2), tp = t, ab = 0, ft = 0, fc = 0;
         cntn < maxit_ ||
           (throw_ && (throw GeographicLib::GeographicErr
                       ("Convergence failure Geodesic3::findroot"), false));) {
      ang xt = 2*t == 1 ?
        ang(xa.s() + xb.s(), xa.c() + xb.c()) :
        (t < tp ? xa - ang::radians(t * ab) :
         xb + ang::radians(tp * ab)),
        xc;
      if (trip) {
        xm = xt;
        if constexpr (debug)
          cout << "BREAKA\n";
        break;
      }
      ++cntn;
      ft = f(xt);
      if constexpr (debug)
        cout << "H " << cntn << " " << real(xt)-x0 << " " << ft << "\n";
      if (!(fabs(ft) >= numeric_limits<real>::epsilon())) {
        xm = xt;
        if constexpr (debug)
          cout << "BREAKB\n";
        break;
      }
      if (signbit(ft) == signbit(fa)) {
        xc = xa; xa = xt;
        fc = fa; fa = ft;
      } else {
        xc = xb; xb = xa; xa = xt;
        fc = fb; fb = fa; fa = ft;
      }
      xm = fabs(fb) < fabs(fa) ? xb : xa;
      // ordering is b - a - c
      ab = (xa-xb).radians0();
      real
        ca = (xc-xa).radians0(),
        cb = ca+ab,
        // Scherer has a fabs(cb).  This should be fabs(ab).
        tl = numeric_limits<real>::epsilon() / fabs(ab);
      // Backward tests to deal with NaNs
      if constexpr (debug)
        cout << "R " << cntn << " " << ab << " " << cb << "\n";
      trip = !(2 * tl < 1);
      if (trip && debug)
        cout << "TRIP " << ab << "\n";
      // Increase the amount away from the boundary to make the next iteration.
      // Otherwise we get two equal values of f near the boundary and a
      // bisection is triggered.
      tl = fmin(1/real(32), 16*tl);
      real
        xi = ab / cb,
        xip = ca / cb,          // 1 - xi
        phi = (fa-fb) / (fc-fb),
        phip = (fc-fa) / (fc-fb); // 1 - phi
      if (!trip && Math::sq(phip) < xip && Math::sq(phi) < xi) {
        t = fa/(fb-fa) * fc/(fb-fc) - ca/ab * fa/(fc-fa) * fb/(fc-fb);
        if constexpr (debug)
          cout << "J1 " << cntn << " "  << t << " " << 1 - t << "\n";
        // This equation matches the pseudocode in Scherer.  His Eq (6.40)
        // reads t = fa/(fb-fa) * fc/(fb-fc) + ca/cb * fc/(fc-fa) * fb/(fb-fa);
        // this is wrong.
        tp = fb/(fb-fa) * fc/(fc-fa) + cb/ab * fa/(fc-fa) * fb/(fc-fb);
        t = fmax(tl, t);
        tp = fmax(tl, tp);
        // t = fmin(1 - tl, fmax(tl, t));
        // tp = fmin(1 - tl, fmax(tl, tp));
      } else {
        t = tp = 1/real(2);
        ++cntb;
        if constexpr (debug)
          cout << "J2 " << cntn << " " << t << " " << 1 - t << "\n";
      }
    }
    if (countn) *countn += cntn;
    if (countb) *countb += cntb;
    return xm;
  }

  Geodesic3::gamblk::gamblk(const Geodesic3& tg,
                            ang bet, ang omg, ang alp) {
    real a = tg.k() * bet.c() * alp.s(), b = tg.kp() * omg.s() * alp.c();
    gamma = (a - b) * (a + b);
    // This direct test case
    // -30 -86 58.455576621187896848 -1.577754271270003
    // fails badly with reverse direct if gamma is not set to zero here.
    // Neighboring values of alp as double are
    // 58.455576621187890, 58.455576621187895, 58.455576621187900
    // 30 86 90 180
    // dgam/dalp = 2*alp.c()*alp.s() * hypot(tg.k * bet.c(), tg.kp * omg.s())
    real maxdiff = 0;
    if (!(tg.k2() == 0  || tg.kp2() == 0)) {
      // Force small gamma to zero for triaxial case
      real
        alpdiff = 2 * alp.c() * alp.s()
        * (tg.k2() * Math::sq(bet.c())+tg.kp2() * Math::sq(omg.s())),
        betdiff = -2 * bet.c() * bet.s() * tg.k2() * Math::sq(alp.s()),
        omgdiff = -2 * omg.c() * omg.s() * tg.kp2() * Math::sq(alp.c());
      maxdiff = fmax( fabs(alpdiff), fmax( fabs(betdiff), fabs(omgdiff) ) );
    }
    if (fabs(gamma) <= 3 * maxdiff * numeric_limits<real>::epsilon()) {
      // Set gamma = 0 if a change of alp, bet, or omg by epsilon would include
      // gamma = 0.
      gamma = 0;
      // If (_umbalt and not oblate) or prolate, set gamma = -0
      if ((tg.umbalt() && tg.kp2() > 0) || tg.k2() == 0) gamma = -gamma;
    }
    transpolar = signbit(gamma);
    gammax = fabs(gamma);
    kx2 = !transpolar ? tg.k2() : tg.kp2();
    kxp2 = transpolar ? tg.k2() : tg.kp2();
    kx = !transpolar ? tg.k() : tg.kp();
    kxp = transpolar ? tg.k() : tg.kp();
    // gammap = sqrt(kx2 - gammax)
    real gammap =
      (!transpolar ?
       hypot(kx * hypot(bet.s(), alp.c()*bet.c()),
             kxp * omg.s()*alp.c()) :
       hypot(kxp *  bet.c()*alp.s(),
             kx * hypot(omg.c(), alp.s()*omg.s())));
    // for gam == 0, we have nu = 0, nup = 1
    nu = sqrt(gammax) / kx;
    nup = gammap / kx;
  }

  Geodesic3::gamblk::gamblk(const Geodesic3& tg, bool neg)
    : transpolar(neg)
    , gamma(transpolar ? -real(0) : real(0))
    , nu(0)
    , nup(1)
    , gammax(0)
    , kx2(!transpolar ? tg.k2() : tg.kp2())
    , kxp2(transpolar ? tg.k2() : tg.kp2())
    , kx(!transpolar ? tg.k() : tg.kp())
    , kxp(transpolar ? tg.k() : tg.kp())
  {}

  Geodesic3::gamblk Geodesic3::gamma(ang bet, ang omg, ang alp) const {
    return gamblk(*this, bet, omg, alp);
  }

  GeodesicLine3 Geodesic3::Line(Angle bet1, Angle omg1, Angle alp1) const {
    return GeodesicLine3(*this, bet1, omg1, alp1);
  }

  GeodesicLine3 Geodesic3::Direct(Angle bet1, Angle omg1, Angle alp1, real s12,
                                  Angle& bet2, Angle& omg2, Angle& alp2)
    const {
    GeodesicLine3 l(*this, bet1, omg1, alp1);
    l.Position(s12, bet2, omg2, alp2);
    return l;
  }

  GeodesicLine3 Geodesic3::Inverse(real bet1, real omg1, real bet2, real omg2,
                                   real& s12, real& alp1, real& alp2) const {
    ang alp1a, alp2a;
    GeodesicLine3 l = Inverse(ang(bet1), ang(omg1), ang(bet2), ang(omg2),
                              s12, alp1a, alp2a);
    alp1 = real(alp1a); alp2 = real(alp2a);
    return l;
  }

  GeodesicLine3 Geodesic3::Line(real bet1, real omg1, real alp1) const {
    return Line(ang(bet1), ang(omg1), ang(alp1));
  }

  GeodesicLine3 Geodesic3::Direct(real bet1, real omg1, real alp1, real s12,
                                  real& bet2, real& omg2, real& alp2)
    const {
    ang bet2a, omg2a, alp2a;
    GeodesicLine3 l = Direct(ang(bet1), ang(omg1), ang(alp1), s12,
                            bet2a, omg2a, alp2a);
    bet2 = real(bet2a); omg2 = real(omg2a); alp2 = real(alp2a);
    return l;
  }

  } // namespace Triaxial
} // namespace GeographicLib

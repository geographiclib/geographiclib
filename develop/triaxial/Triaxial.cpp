/**
 * \file Triaxial.cpp
 * \brief Implementation for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "Triaxial.hpp"
#include "TriaxialLine.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <string>

namespace GeographicLib {

  using namespace std;

  Triaxial::Triaxial()
    : Triaxial(1, 0, 1, 0)
  {}

  Triaxial::Triaxial(Math::real a, Math::real b, Math::real c)
    : _a(a)
    , _b(b)
    , _c(c)
    , _axes({_a, _b, _c})
    , _umbalt(false)
    , _oblpro(false)
    , _merid(false)
    , _debug(false)
    , _ellipthresh(1/real(8))
  {
    real s = (_a - _c) * (_a + _c);
    _e2 = s / Math::sq(_b);
    if (s == 0) {
      // The sphere is a nonuniform limit, we can pick any values in [0,1]
      // s.t. k2 + kp2 = 1.  Here we choose to treat the sphere as an
      // oblate ellipsoid.
      _kp2 = 0; _k2 = 1 - _kp2;
    } else {
      _kp2 = (_a - _b) * (_a + _b) / s;
      _k2  = (_b - _c) * (_b + _c) / s;
    }
    _k = sqrt(_k2); _kp = sqrt(_kp2);
    if (! (isfinite(_a) && isfinite(_b) && isfinite(_c) &&
           _a >= _b && _b >= _c && _c >= 0 && _b > 0) )
      throw GeographicErr("Bad semiaxes for triaxial ellipsoid");
  }

  Triaxial::Triaxial(Math::real b, Math::real e2,
                     Math::real k2, Math::real kp2)
    : _b(b)
    , _e2(e2)
    , _k2(k2)
    , _kp2(kp2)
    , _umbalt(false)
    , _oblpro(false)
    , _merid(false)
    , _debug(false)
    , _ellipthresh(1/real(8))
  {
    real ksum = _k2 + _kp2;
    _k2 /= ksum;
    _kp2 /= ksum;
    _k = sqrt(_k2);
    _kp = sqrt(_kp2);
    _a = _b * sqrt(1 + _e2 * _kp2);
    _c = _b * sqrt(1 - _e2 * _k2);
    _axes = vec3({_a, _b, _c});
    if (! (isfinite(_a) && isfinite(_b) && isfinite(_c) &&
           _a >= _b && _b >= _c && _c >= 0 && _b > 0) )
      throw GeographicErr("Bad semiaxes for triaxial ellipsoid");
  }

  void Triaxial::Norm(vec3& r) const {
    vec3 rn{r[0] / _axes[0], r[1] / _axes[1], r[2] / _axes[2]};
    real ra = Math::hypot3(rn[0], rn[1], rn[2]);
    r = {r[0] / ra, r[1] / ra, r[2] / ra};
  }

  void Triaxial::Norm(vec3& r, vec3& v) const {
    Norm(r);
    vec3 axes2 = {Math::sq(_axes[0]), Math::sq(_axes[1]), Math::sq(_axes[2])},
      up = {r[0] / axes2[0], r[1] / axes2[1], r[2] / axes2[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      uv = up[0] * v[0] + up[1] * v[1] + up[2] * v[2],
      f = uv/u2;
    v = {v[0] - f * up[0], v[1] - f * up[1], v[2] - f * up[2]};
    normvec(v);
  }

  void Triaxial::cart2toellip(vec3 r, Angle& bet, Angle& omg)
    const {
    real xi = r[0]/_a, eta = r[1]/_b, zeta = r[2]/_c,
      g = _k2 * Math::sq(xi)
      + (_k2 - _kp2) * Math::sq(eta)
      - _kp2 * Math::sq(zeta);
    if (fabs(r[0]) == _a * _kp2 && r[1] == 0 && fabs(r[2]) == _c * _k2)
      g = 0;
    real h = hypot(g, 2 * _k * _kp * eta),
      so, co, sb, cb;
    if (h == 0) {
      so = 0;
      cb = 0;
    } else if (g < 0) {
      so = copysign(sqrt( (h - g)/2 ) / _kp, eta);
      cb = fabs(eta / so);
    } else {
      cb = sqrt( (h + g)/2 ) / _k;
      so = eta / cb;
    }
    real tz = hypot(_k, _kp * so),
      tx = hypot(_k * cb, _kp);
    sb = tz == 0 ? -1 : zeta / tz;
    co = tx == 0 ?  1 : xi / tx;
    bet = ang(sb, cb, 0, true); omg = ang(so, co, 0, true);
  }

  void Triaxial::cart2toellip(Angle bet, Angle omg,
                              vec3 v, Angle& alp) const {
    real tz = hypot(_k, _kp * omg.s()), tx = hypot(_k * bet.c(), _kp);
    // At oblate pole tx = 0; at prolate pole, tz = 0
    if (tx == 0 || tz == 0 || !(bet.c() == 0 && omg.s() == 0)) {
      // Not a triaxial umbilical point
      vec3
        N = tx == 0 ?
        vec3{-omg.c() * bet.s(), -omg.s() * bet.s(), tx * bet.s()} :
        (tz == 0 ?
         vec3{tz, -bet.s(), bet.c()} :
         vec3{-_a * _k2 * bet.c() * bet.s() * omg.c() / tx,
              -_b * bet.s() * omg.s(),
               _c * bet.c() * tz}),
        E = tx == 0 ?
        vec3{-omg.s(),  omg.c(), tx} :
        (tz == 0 ?
         vec3{tz * omg.c(),  bet.c() * omg.c(), bet.s() * omg.c()} :
         vec3{-_a * tx * omg.s(),
               _b * bet.c() * omg.c(),
               _c * _kp2 * bet.s() * omg.c() * omg.s() / tz});
      normvec(N); normvec(E);
      alp = ang(v[0] * E[0] + v[1] * E[1] + v[2] * E[2],
                v[0] * N[0] + v[1] * N[1] + v[2] * N[2], 0, true);
    } else {                    // bet.c() == 0 && omg.s() == 0
      // Special treatment at umbilical points
      real w = bet.s() * omg.c(),
        upx = omg.c() * tx/_a, upz = bet.s() * tz/_c;
      Math::norm(upx, upz);
      // compute direction cosines of v wrt the plane y = 0; angle = 2*alp
      real s2a = -v[1] * w, c2a = (upz*v[0] - upx*v[2]) * w;
      // Unnecessary: Math::norm(s2a, c2a)
      // We have
      //   2*alp = atan2(s2a, c2a), h2 = hypot(s2a, c2a)
      //   alp = atan2(sa, ca)
      //   tan(2*alp) = 2*tan(alp)/(1-tan(alp)^2)
      // for alp in [-pi/2, pi/2]
      // c2a>0
      //   [sa, ca] = [s2a / sqrt(2*(1+c2a)), sqrt((1+c2a)/2)]
      //           -> [s2a, h2+c2a]
      // c2a<0
      //   [sa, ca] = sign(s2a)*[sqrt((1-c2a)/2), s2a / sqrt(2*(1-c2a))]
      //           -> [sign(s2a) * (h2-c2a), abs(s2a)]
      // for northern umbilical points, we want to flip alp to alp + pi; so
      // multiply [sa, ca] by -bet.s().
      real flip = -bet.s();
      if (c2a >= 0)
        alp = ang(flip * s2a, flip * (1 + c2a));
      else
        alp = ang(flip * copysign(1 - c2a, s2a), flip * fabs(s2a));
    }
  }

  void Triaxial::cart2toellip(vec3 r, vec3 v,
                              Angle& bet, Angle& omg, Angle& alp)
    const {
    cart2toellip(r, bet, omg);
    cart2toellip(bet, omg, v, alp);
  }

  void Triaxial::elliptocart2(Angle bet, Angle omg,
                              vec3& r) const {
    real tx = hypot(_k * bet.c(), _kp), tz = hypot(_k, _kp * omg.s());
    r = vec3{ _a * omg.c() * tx,
              _b * bet.c() * omg.s(),
              _c * bet.s() * tz };
    // Norm(r); r is already normalized
  }

  void Triaxial::elliptocart2(Angle bet, Angle omg,
                              Angle alp,
                              vec3& r, vec3& v) const {
    elliptocart2(bet, omg, r);
    real tx = hypot(_k * bet.c(), _kp), tz = hypot(_k, _kp * omg.s());
    if (bet.c() == 0 && omg.s() == 0 && !(_k == 0 || _kp == 0)) {
      // umbilical point (not oblate or prolate)
      real sa2 = 2 * alp.s() * alp.c(),
        ca2 = (alp.c() - alp.s()) * (alp.c() + alp.s());
      // sign on 2nd component is -sign(cos(bet)*sin(omg)).  negative sign
      // gives normal convention of alpha measured clockwise.
      v = vec3{_a*_k/_b * omg.c() * ca2,
               -omg.c() * bet.s() * sa2,
               -_c*_kp/_b * bet.s() * ca2};
    } else {
      vec3 N, E;
      if (tx == 0) {
        // At an oblate pole tx -> cos(bet)
        N = vec3{-omg.c() * bet.s(), -omg.s() * bet.s(), 0};
        E = vec3{-omg.s()          ,  omg.c()          , 0};
      } else if (tz == 0) {
        // At a prolate pole tz -> sin(omg)
        N = vec3{0, -bet.s()          , bet.c()          };
        E = vec3{0,  bet.c() * omg.c(), bet.s() * omg.c()};
      } else {
        // The general case
        N = vec3{ -_a * _k2 * bet.c() * bet.s() * omg.c() / tx,
                  -_b * bet.s() * omg.s(),
                   _c * bet.c() * tz};
        E = vec3{ -_a * tx * omg.s(),
                   _b * bet.c() * omg.c(),
                   _c * _kp2 * bet.s() * omg.c() * omg.s() / tz};
      }
      normvec(N);
      normvec(E);
      v = vec3{alp.c() * N[0] + alp.s() * E[0],
               alp.c() * N[1] + alp.s() * E[1],
               alp.c() * N[2] + alp.s() * E[2]};
    }
    // normvec(v); v is already normalized
  }

  TriaxialLine Triaxial::Inverse(Angle bet1, Angle omg1,
                                 Angle bet2, Angle omg2,
                                 Angle& alp1, Angle& alp2, real& s12) const {
    typedef TriaxialLine TL;
    string msg;
    bet1.round();
    omg1.round();
    bet2.round();
    omg2.round();
    bool oblate = _kp2 == 0, prolate = _k2 == 0;

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
    bool flip1 = AngNorm(bet1, omg1, prolate),
      flip2 = AngNorm(bet2, omg2, prolate);
    // oblate, set these to bet[12].c() == 0
    // prolate, set these to omg[12].s() == 0
    bool umb1 = (prolate || bet1.c() == 0) && (oblate || omg1.s() == 0),
      umb2 = (prolate || bet2.c() == 0) && (oblate || omg2.s() == 0);
    bool swap12;
    {
      // For prolate, swap based on omg, switch 1 & 2 because poles are at
      // 0/180, instead of +/-90.
      ang tmp1(prolate ? omg2 : bet1), tmp2(prolate ? omg1 : bet2);
      tmp1.setquadrant(0U); tmp2.setquadrant(0U);
      ang tmp12 = tmp2 - tmp1; // |bet2| - |bet1|
      swap12 = tmp12.s() > 0; // is |bet2| > |bet1|
      if (!oblate && !prolate && tmp12.s() == 0) {
        // don't need to do this if oblate or prolate
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
      swap(umb1, umb2);
    }
    if (oblate) {
      // Rotate, subtracting omg1 from omg[12], so omg1 = 0
      omg2 -= omg1;
      omg1 = ang::cardinal(0);
    } else if (prolate) {
      // Rotate, subtracting bet1 + 90 from bet[12], so bet1 = -90
      bet2 -= bet1 + ang::cardinal(1);
      bet1 = ang::cardinal(-1);
    }
    // Now |bet1| >= |bet2|
    bool flipz = bet1.s() > 0;
    if (flipz) {                // Not needed for prolate, bet1 already -90
      bet1.reflect(true);
      bet2.reflect(true);
    }
    // Now bet1 <= 0
    bool flipy = prolate ? signbit(bet2.c()) :
      signbit(omg1.s()) || (omg1.s() == 0 && signbit(omg2.s()));
    if (flipy) {
      if (prolate)
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

    // Reset umb[12] w/o treating prolate/oblate specially.
    umb1 = bet1.c() == 0 && omg1.s() == 0;
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
    //       found
    //     Prolate: treated by middle ellipse cases A.c.2 and A.c.3.
    //
    //     1 If point 2 is within the resulting span of longitudes, the
    //       geodesic is equatorial (compute with ArcPos0).
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
    //       Oblate/prolate: opposite poles
    //
    //     2 bet1 = -90, bet2 = -90: geodesic follows ellipse (compute with
    //       ArcPos0).  Points may be adjacent umbilical points.
    //       Oblate: same pole
    //       Prolate: meridional geodesic, same side of poles

    //     3 bet1 = -90, bet2 = 90: geodesic may follow the ellipse if if they
    //       are close enough.  See case A.b for strategy.  Points may be
    //       adjacent umbilical points.
    //       Oblate: treated by A.c.1
    //       Prolate: meridional geodesic, opposite sides of poles
    //
    //     4 |bet2| != 90: geodesic follows the ellipse.  Compute with Hybrid.
    //       Oblate: meridinal geodesic
    //       Prolate: same pole
    //
    // B The geodesic does not follow a principal ellipse
    //
    //   a point 1 is an umbilical point, gam = 0, so we know alp2
    //     Oblate: remaining meridional cases
    //     Prolate: remaining meridional cases
    //
    //   b bet1 = -90: do a search with alp1 bracket = [-90+eps, 90-eps].
    //     Oblate: already treated by B.a
    //     Prolate: same as triaxial case
    //
    //   c omg1 = 0: do a search with alp1 bracket = [eps, 180-eps], or
    //     [-180+eps, -eps], depending on the sign of omg2.
    //
    //   d general case, bracket search by two neighboring umbilical
    //     directions.  This treats case A.a.
    //     Oblate: there are only two distinct meridional directions, we can
    //       automatically bracket the search in alp1 by [0,180] (because omg1
    //       = 0, omg2 > 0).
    //     Prolate: ditto but bracket alp1 by [-90,90] (because bet1 = -90,
    //       bet2 in (-90,90)).

    // Set up variables for search
    real fa, fb;
    ang alpa, alpb;

    // and for the final result
    ang bet2a, omg2a;

    TL::fline lf(*this);
    TL::fline::fics fic;
    TL::fline::disttx d;

    if (_debug)
      cout << "COORDS " << real(bet1) << " " << real(omg1) << " "
           << real(bet2) << " " << real(omg2) << "\n";
    //    bet1.setn(0); omg1.setn(0); bet2.setn(0); omg2.setn(0);
    // flag for progress
    bool done = false, backside = false;
    if (bet1.c() * omg1.s() == 0 && bet2.c() * omg2.s() == 0) {
      // Case A.c, both points on middle ellipse
      lf = TL::fline(*this, gamblk(*this, (_umbalt && _kp2 > 0) || _k2 == 0));
      if (umb1 && umb2 && bet2.s() > 0 && omg2.c() < 0) {
        // Case A.c.1, process opposite umbilical points
        // For oblate/prolate this gives 0/90
        alp1 = oblate || prolate ? ang(_kp, _k, 0, true) :
          ang(exp(lf.deltashift()/2), 1);
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        bool betp = _k2 > _kp2;
        d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2, betp);
        if (oblate || prolate) alp2 += ang::cardinal(oblate ? 2 : -1);
        if (_debug) msg = "A.c opposite umbilics";
        backside = signbit(bet2a.c());
        done = true;
      } else if (bet1.c() == 0 && bet2.c() == 0) {
        // Case A.c.{2,3}, bet1 = -90, bet2 = +/-90
        if (bet2.s() < 0) {
          // Case A.c.2, bet1 = bet2 = -90
          // If oblate, bet1 = -90, omg1 = 0, need alp1 = omg2 to follow omg2
          // meridian.
          alp1 = ang::cardinal(oblate ? 2 : 1);
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          ang omg12 = omg2 - omg1;
          if (omg12.s() == 0 && omg12.c() < 0) {
            // adjacent E/W umbilical points
            // Should be able to get ArcPos0 to return this?
            d = oblate ?
              TL::fline::disttx{ -lf.fpsi().Max(), lf.ftht().Max(), 0 } :
              prolate ?
              TL::fline::disttx{ lf.fpsi().Max(), -lf.ftht().Max(), 0 } :
              TL::fline::disttx{ -BigValue(), BigValue(), 0 };
            if (_debug) msg = "A.c.2 adjacent EW umbilics";
            alp2 = ang::cardinal(prolate ? 1 : 0);
          } else {
            d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
            if (_debug) msg = "A.c.2 bet1/2 = -90";
          }
          if (omg2a.s() < 0) alp2.reflect(true); // Is this needed?
          done = true;
        } else {
          // Case A.c.3, bet1 = -90, bet2 = 90
          // need to see how far apart the points are
          // If point 1 is at [-90, 0], direction is 0 else -90.
          // XXX Maybe alp1 needs fixing

          alp1 = ang::cardinal(omg1.s() == 0 && (!prolate || omg2.s() == 0) ?
                               0 : -1);
          if (0)
            cout << "FIC " << real(bet1) << " " << real(omg1) << " "
                 << real(alp1) << "\n";
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          // If point 1 is [-90, 0] and point 2 is [90, 0]
          if (omg1.s() == 0 && omg2.s() == 0) {
            // adjacent N/S umbilical points
            // Should be able to get ArcPos0 to return this?
            d = oblate ?
              TL::fline::disttx{ lf.fpsi().Max(), -lf.ftht().Max(), 0 } :
              prolate ?
              TL::fline::disttx{ -lf.fpsi().Max(), lf.ftht().Max(), 0 } :
              TL::fline::disttx{ BigValue(), -BigValue(), 0 };
            alp2 = ang::cardinal(oblate ? 0 : (prolate ? 2 : 1));
            if (_debug) msg = "A.c.3 adjacent NS umbilics";
            done = true;
          } else {
            // FIX ME for oblate
            if (omg1.s() == 0)
              omg2a = ang::cardinal(2);
            else {
              d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
            }
            // XXX FIX HERE for prolate case -90 -1 90 177
            omg2a -= omg2;
            if (omg2a.s() > 0) {
              ang omg12 = omg2 + omg1;
              // FIX ME for oblate
              d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
              if (_debug)
                cout << "APX "
                     << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
                     << real(bet2a) << " " << real(omg2a) << " "
                     << real(alp2) << "\n";
              if (omg2a.s() < 0) alp2.reflect(true); // Is this needed?
              if (_debug) msg = "A.c.3 bet1/2 = -/+90 meridional";
              done = true;
            } else {
              alpa = ang::cardinal(-1) + ang::eps();
              if (false && prolate) {
                fa = -omg2.radians0();
                alpb = -alpa;
                fb = (ang::cardinal(2) - omg2).radians0();
              } else {
                fa = omg2a.radians0();
                fic.setquadrant(lf, 0U);
                (void) lf.ArcPos0(fic, ang::cardinal(2),
                                  bet2a, omg2a, alp2);
                omg2a -= omg2;
                alpb = -alpa;
                fb = omg2a.radians0();
              }
              if (0)
                cout << "ALP/F "
                     << real(alpa) << " " << fa << " "
                     << real(alpb) << " " << fb << "\n";
              if (_debug) msg = "A.c.3 general bet1/2 = -/+90, non-meridional";
            }
          }
        }
      } else {
        // Case A.c.4, other meridional cases, invoke Hybrid with the following
        // value of alp1
        // If oblate, bet1 = -90, omg1 = 0, need alp1 = omg2 to follow omg2
        // meridian.
        // If prolate, bet1 = -90, omg1 = 0, need alp1 = -bet2
        alp1 = oblate ? omg2 :
          (prolate ? -bet2 :
           ang::cardinal(bet1.c() == 0 ?
                         // TODO: CHECK omg2.c() < 1 test; CHANGE TO < 0
                         (omg2.c() < 0 ? 1 :
                          (omg1.s() == 0 && !prolate ? 0 : -1)) :
                         (omg2.c() > 0 ? 0 : 2)));
        if (_debug)
          cout << "ALP1 " << real(alp1) << " "
               << bet1.c() << " " << omg2.c() << "\n";
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        d = prolate ?
          lf.ArcPos0(fic, (omg2-omg1).base(), bet2a, omg2a, alp2, false) :
          lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
        if (prolate)
          alp2 += ang::cardinal(1);
        if (_debug) msg = "A.c.4 other meridional";
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
      lf = TL::fline(*this, gamma(bet1, omg1, alp1));
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      (void) lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
      omg2a -= omg2;
      if (E * omg2a.s() >= 0) {
        // geodesic follows the equator
        d = lf.ArcPos0(fic, omg12.flipsign(E), bet2a, omg2a, alp2, false);
        if (_debug) msg = "A.b.1 bet1/2 = 0 equatorial";
        done = true;
      } else {
        // geodesic does not follow the equator
        alpb = ang::cardinal(-1) - ang::eps();
        alpa = -alpb;
        (E > 0 ? fa : fb) = omg2a.radians0();
        alp1.setquadrant(E > 0 ? 3U : 0U);
        fic.setquadrant(lf, E > 0 ? 3U : 0U);
        (void) lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
        omg2a -= omg2;
        (E > 0 ? fb : fa) = omg2a.radians0();
        if (_debug) msg = "A.b.2 general bet1/2 = 0 non-equatorial";
      }
    } else if (umb1) {
      // Case B.a, umbilical point to general point
      lf = TL::fline(*this, gamblk(*this, (_umbalt && _kp2 > 0) || _k2 == 0));
      alp2 = ang(_kp * omg2.s(), _k * bet2.c());
      // RETHINK THIS.  If we know alp2, we can compute delta.  This should be
      // enough to find alp1.
      if (0)
        cout << "ALP2 " << real(alp2) << " " << real(bet1 - bet2) << "\n";
      fic = TL::fline::fics(lf, bet2, omg2, alp2);
      bool betp = _k2 > _kp2;   // This could be betb = !prolate;
      if (0) {
        (void) lf.ArcPos0(fic, (betp ? bet1 - bet2 :  omg1 - omg2).base(),
                          bet2a, omg2a, alp1, betp);
        //      (void) lf.ArcPos0(fic, bet1 - bet2, bet2a, omg2a, alp1);
        if (0)
          cout << "APOUT "
               << real(bet2a) << " " << real(omg2a) << " " << real(alp1) << "\n";
        if (alp1.s() < 0) alp1 += ang::cardinal(1);
        if (prolate) alp1 += ang::cardinal(1);
        if (oblate) alp1 += omg2;
      } else {
        // cout << "DELTA " << fic.delta/Math::degree() << "\n";
        real delta = (lf.transpolar() ? -1 : 1) * fic.delta;
        alp1 = oblate ?
          // For oblate
          // delta =
          //   atan2(bet1.s() * fabs(alp1.s()), bet0.c() * alp1.c())) - omg1
          // Here omg1 = -90 (because of pi/2 shift in ext. vs int. omg)
          // bet1.s() = -1, bet0.c() = 1, alp1.s() > 0, so
          // alp1 = 90 - delta
          ang::cardinal(1) - ang::radians(delta) :
          (prolate ?
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
      }
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      d = lf.ArcPos0(fic, (betp ? bet2 - bet1 : omg2 - omg1).base(),
                     bet2a, omg2a, alp2, betp);
      if (_debug) msg = "B.a umbilic to general";
      done = true;
    } else if (bet1.c() == 0) {
      // Case B.b, bet1 = -90 to general point
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
      if (_debug) msg = "B.b general bet1 = -90";
    } else if (omg1.s() == 0) {
      // Case B.c, omg1 = 0 to general point
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
      if (_debug) msg = "B.c general omg1 = 0";
    } else {
      // Case B.d, general case
      real f[4];
      alpa = ang( _kp * fabs(omg1.s()), _k * fabs(bet1.c()));
      alpb = alpa;

      lf = TL::fline(*this, gamblk(*this, (_umbalt && _kp2 > 0) || _k2 == 0));
      fic = TL::fline::fics(lf, bet1, omg1, alpb);
      unsigned qb = 0U, qa = 3U; // qa = qb - 1 (mod 4)
      if (_debug) msg = "B.d general";
      for (; !done && qb <= 4U; ++qb, ++qa) {
        if (qb) {
          alpb.setquadrant(qb);
          fic.setquadrant(lf, qb);
        }
        if (qb < 4U) {
          f[qb] = lf.Hybrid0(fic, bet2, omg2);
          if (_debug)
            cout << "f[qb] " << qb << " " << f[qb] << "\n";
          if (fabs(f[qb]) < numeric_limits<real>::epsilon()) {
            alp1 = alpb;
            d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
            if (_debug) msg = "B.d accidental umbilic";
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
      if (_debug)
        cout << "fDD " << done << " " << qa << " " << qb << " "
             << f[qa & 3U] << " " << f[qb & 3U] << "\n";
      if (!done) {
        fa = f[qa & 3U]; fb = f[qb & 3U];
        alpa.setquadrant(qa);
      }
    }

    int countn = 0, countb = 0;
    if (!done) {
      // Iterative search for the solution
      if (_debug)
        cout << "X " << done << " " << msg << "\n";
      alp1 = findroot(
                      [this, &bet1, &omg1, &bet2, &omg2]
                      (const ang& alp) -> Math::real
                      {
                        return HybridA(bet1, omg1, alp, bet2, omg2);
                      },
                      alpa,  alpb,
                      fa, fb,
                      &countn, &countb);
      if (_debug)
        cout << "ALP1 " << real(alp1) << "\n";
      alp1.round();
      lf = TL::fline(*this, gamma(bet1, omg1, alp1));
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
      backside = signbit(bet2a.c());
    }

    if (backside) alp2.reflect(true, true);
    alp2.round();

    TL::gline lg(*this, lf.gm());
    TL::gline::gics gic(lg, fic);
    s12 = lg.dist(gic, d);

    if (_debug)
      cout << "FLIPS " << flip1 << flip2 << flipz << flipy << flipx << flipomg
           << "\n";
    if (_debug)
      cout << "A "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    // Undo switches in reverse order flipz, swap12, flip1
    if (flipomg) {
      omg2.reflect(true);
      alp2.reflect(true, true);
    }

    if (_debug)
      cout << "B "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
      alp1.reflect(true);
      alp2.reflect(true);
    }

    if (_debug)
      cout << "C "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipy) {
      if (prolate) {
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

    if (_debug)
      cout << "D "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipz) {
      bet1.reflect(true);
      bet2.reflect(true);
      alp1.reflect(false, true);
      alp2.reflect(false, true);
    }

    if (_debug)
      cout << "E "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
      swap(alp1, alp2);
      swap(umb1, umb2);
      // points not swapped if umb1 == true
      alp1 += ang::cardinal(2);
      if (umb2 && !(prolate || oblate))
        alp2 += ang::cardinal((signbit(alp2.s()) ? -1 : 1) * bet2.s());
      else
        alp2 += ang::cardinal(2);
    }

    if (_debug)
      cout << "F "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flip1) Flip(bet1, omg1, alp1);
    if (flip2) Flip(bet2, omg2, alp2);
    alp1.setn(); alp2.setn();
    if (_debug)
      cout << "G "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";

    fic = TL::fline::fics(lf, bet10, omg10, alp1);
    gic = TL::gline::gics(lg, fic);
    gic.s13 = 0 + fmax(real(0), s12);

    if (_debug)
      cout << countn << " " << countb << " "
           << lf.gamma() << " "
           << lf.fbet().NCoeffs() << " " << lf.fomg().NCoeffs() << " "
           << lg.gbet().NCoeffs() << " " << lg.gomg().NCoeffs() << " MSG "
           << msg << "\n";
    // clang needs std::move instead of move.
    return TL(std::move(lf), std::move(fic), std::move(lg), std::move(gic));
  }

  Math::real Triaxial::HybridA(Angle bet1, Angle omg1,
                               Angle alp1,
                               Angle bet2a, Angle omg2b) const {
    ang b1{bet1}, o1{omg1}, a1{alp1};
    // a1 -= ang(1e-8);
    gamblk gam = gamma(b1, o1, a1);
    TriaxialLine::fline l(*this, gam);
    real domg;
    TriaxialLine::fline::fics ic(l, b1, o1, a1);
    domg = l.Hybrid0(ic, bet2a, omg2b);
    if (_debug)
      cout << "HA " << signbit(gam.gamma) << " " << gam.gamma << " "
           << real(bet1) << " " << real(omg1) << " "
           << real(alp1) << " " << real(bet2a) << " " << real(omg2b) << " "
           << domg << "\n";
    return domg;
  }

  // Solve f(alp1) = 0 where alp1 is an azimuth and f(alp1) is the difference
  // in lontitude on bet2 and the target longitude.
  Angle Triaxial::findroot(const function<Math::real(const Angle&)>& f,
                           Angle xa,  Angle xb,
                           Math::real fa, Math::real fb,
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
    int cntn = 0, cntb = 0, maxcnt = 100;
    bool trip = false;
    bool debug = false;
    // 25 iterations with line 498534 of testset.txt
    //   43 -2 -43 141
    // This is a near conjugate case m12 = 0.0003857
    // 27 iterations with line 360115 of testpro.txt
    //   -75 29 75 -169
    // Converge failures with line 40045 of testsph[bc]
    // echo 80 -90 -80 90 | ./Geod3Solve -e 1 0 2 1 -i
    // echo 80 -90 -80 90 | ./Geod3Solve -e 1 0 1 2 -i
    real x0 = -88.974388706914264728; // Offset for debugging output
    if (debug)
      cout << "H " << real(xa) << " " << fa << " "
           << real(xb) << " " << fb << " "
           << f(ang(x0)) << "\n";
    // tp = 1 - t
    for (Math::real t = 1/Math::real(2), tp = t, ab = 0, ft = 0, fc = 0;
         cntn < maxcnt ||
           (throw GeographicLib::GeographicErr
            ("Convergence failure Triaxial::findroot"), false)
           || GEOGRAPHICLIB_PANIC("Convergence failure Triaxial::findroot");) {
      ang xt = 2*t == 1 ?
        ang(xa.s() + xb.s(), xa.c() + xb.c()) :
        (t < tp ? xa - ang::radians(t * ab) :
         xb + ang::radians(tp * ab)),
        xc;
      if (trip) {
        xm = xt;
        if (debug)
          cout << "BREAKA\n";
        break;
      }
      ++cntn;
      ft = f(xt);
      if (debug)
        cout << "H " << cntn << " " << real(xt)-x0 << " " << ft << "\n";
      if (!(fabs(ft) >= numeric_limits<Math::real>::epsilon())) {
        xm = xt;
        if (debug)
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
      Math::real
        ca = (xc-xa).radians0(),
        cb = ca+ab,
        // Scherer has a fabs(cb).  This should be fabs(ab).
        tl = numeric_limits<Math::real>::epsilon() / fabs(ab);
      // Backward tests to deal with NaNs
      if (debug)
        cout << "R " << cntn << " " << ab << " " << cb << "\n";
      trip = !(2 * tl < 1);
      if (trip && debug)
        cout << "TRIP " << ab << "\n";
      // Increase the amount away from the boundary to make the next iteration.
      // Otherwise we get two equal values of f near the boundary and a
      // bisection is triggered.
      tl = fmin(1/real(32), 16*tl);
      Math::real
        xi = ab / cb,
        xip = ca / cb,          // 1 - xi
        phi = (fa-fb) / (fc-fb),
        phip = (fc-fa) / (fc-fb); // 1 - phi
      if (!trip && Math::sq(phip) < xip && Math::sq(phi) < xi) {
        t = fa/(fb-fa) * fc/(fb-fc) - ca/ab * fa/(fc-fa) * fb/(fc-fb);
        if (debug)
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
        t = tp = 1/Math::real(2);
        ++cntb;
        if (debug)
          cout << "J2 " << cntn << " " << t << " " << 1 - t << "\n";
      }
    }
    if (countn) *countn += cntn;
    if (countb) *countb += cntb;
    return xm;
  }

  Triaxial::gamblk::gamblk(const Triaxial& t,
                           Angle bet, Angle omg, Angle alp) {
    real a = t.k() * bet.c() * alp.s(), b = t.kp() * omg.s() * alp.c();
    gamma = (a - b) * (a + b);
    // This direct test case
    // -30 -86 58.455576621187896848 -1.577754271270003
    // fails badly with reverse direct if gamma is not set to zero here.
    // Neighboring values of alp as double are
    // 58.455576621187890, 58.455576621187895, 58.455576621187900
    // 30 86 90 180
    // dgam/dalp = 2*alp.c()*alp.s() * hypot(t.k * bet.c(), t.kp * omg.s())
    real
      alpdiff = 2 * alp.c() * alp.s()
      * (t.k2() * Math::sq(bet.c())+t.kp2() * Math::sq(omg.s())),
      betdiff = -2 * bet.c() * bet.s() * t.k2() * Math::sq(alp.s()),
      omgdiff = -2 * omg.c() * omg.s() * t.kp2() * Math::sq(alp.c()),
      maxdiff = fmax( fabs(alpdiff), fmax( fabs(betdiff), fabs(omgdiff) ) );
    // cout << "GAMDIFF " << gamma/ numeric_limits<real>::epsilon() << " "
    //      << maxdiff << "\n";
    if (fabs(gamma) <= 2 * maxdiff * numeric_limits<real>::epsilon()) {
      // Set gamma = 0 if a change of alp, bet, or omg by epsilon would include
      // gamma = 0.
      gamma = 0;
      // If (_umbalt and not oblate) or prolate, set gamma = -0
      if ((t.umbalt() && t.kp2() > 0) || t.k2() == 0) gamma = -gamma;
    }
    transpolar = signbit(gamma);
    gammax = fabs(gamma);
    kx2 = !transpolar ? t.k2() : t.kp2();
    kxp2 = transpolar ? t.k2() : t.kp2();
    kx = !transpolar ? t.k() : t.kp();
    kxp = transpolar ? t.k() : t.kp();
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

  Triaxial::gamblk::gamblk(const Triaxial& t, bool neg)
    : transpolar(neg)
    , gamma(transpolar ? -real(0) : real(0))
    , nu(0)
    , nup(1)
    , gammax(0)
    , kx2(!transpolar ? t.k2() : t.kp2())
    , kxp2(transpolar ? t.k2() : t.kp2())
    , kx(!transpolar ? t.k() : t.kp())
    , kxp(transpolar ? t.k() : t.kp())
  {}

  Triaxial::gamblk Triaxial::gamma(Angle bet, Angle omg, Angle alp) const {
    return gamblk(*this, bet, omg, alp);
#if 0
    real a = _k * bet.c() * alp.s(), b = _kp * omg.s() * alp.c(),
      gam = (a - b) * (a + b);
    // This direct test case
    // -30 -86 58.455576621187896848 -1.577754271270003
    // fails badly with reverse direct if gam is not set to zero here.
    // Neighboring values of alp as double are
    // 58.455576621187890, 58.455576621187895, 58.455576621187900
    // 30 86 90 180
    // dgam/dalp = 2*alp.c()*alp.s() * hypot(_k * bet.c(), _kp * omg.s())
    real
      alpdiff = 2 * alp.c() * alp.s()
      * (_k2 * Math::sq(bet.c())+_kp2 * Math::sq(omg.s())),
      betdiff = -2 * bet.c() * bet.s() * _k2 * Math::sq(alp.s()),
      omgdiff = -2 * omg.c() * omg.s() * _kp2 * Math::sq(alp.c()),
      maxdiff = fmax( fabs(alpdiff), fmax( fabs(betdiff), fabs(omgdiff) ) );
    // cout << "GAMDIFF " << gam/ numeric_limits<real>::epsilon() << " "
    //      << maxdiff << "\n";
    if (fabs(gam) <= 2 * maxdiff * numeric_limits<real>::epsilon()) {
      // Set gam = 0 if a change of alp, bet, or omg by epsilon would include
      // gam = 0.
      gam = 0;
      // If (_umbalt and not oblate) or prolate, set gam = -0
      if ((_umbalt && _kp2 > 0) || _k2 == 0) gam = -gam;
    }
    real gamp = false ? 0 :
      (!signbit(gam) ? // sqrt(k2 - gamma)
       hypot(_k * hypot(bet.s(), alp.c()*bet.c()),
             _kp * omg.s()*alp.c()) :
       // sqrt(kp2 + gamma)
       hypot(_k *  bet.c()*alp.s(),
             _kp * hypot(omg.c(), alp.s()*omg.s()))),
      // for gam == 0, we have nu = 0, nup = 1
      nu = sqrt(fabs(gam)) / (!signbit(gam) ? _k : _kp),
      nup = gamp / (!signbit(gam) ? _k : _kp);
    return gamblk(gam, nu, nup);
#endif
  }

  TriaxialLine Triaxial::Line(Angle bet1, Angle omg1, Angle alp1) const {
    return TriaxialLine(*this, bet1, omg1, alp1);
  }

  TriaxialLine Triaxial::Direct(Angle bet1, Angle omg1, Angle alp1, real s12,
                                Angle& bet2, Angle& omg2, Angle& alp2) const {
    TriaxialLine l(*this, bet1, omg1, alp1);
    l.Position(s12, bet2, omg2, alp2);
    return l;
  }

  TriaxialLine Triaxial::Inverse(real bet1, real omg1,
                                 real bet2, real omg2,
                                 real& alp1, real& alp2, real& s12) const {
    ang alp1a, alp2a;
    TriaxialLine l = Inverse(ang(bet1), ang(omg1), ang(bet2), ang(omg2),
                             alp1a, alp2a, s12);
    alp1 = real(alp1a); alp2 = real(alp2a);
    return l;
  }

  TriaxialLine Triaxial::Line(real bet1, real omg1, real alp1) const {
    return Line(ang(bet1), ang(omg1), ang(alp1));
  }

  TriaxialLine Triaxial::Direct(real bet1, real omg1, real alp1, real s12,
                                real& bet2, real& omg2, real& alp2) const {
    ang bet2a, omg2a, alp2a;
    TriaxialLine l = Direct(ang(bet1), ang(omg1), ang(alp1), s12,
                            bet2a, omg2a, alp2a);
    bet2 = real(bet2a); omg2 = real(omg2a); alp2 = real(alp2a);
    return l;
  }

  Math::real Triaxial::EuclideanInverse(Angle bet1, Angle omg1,
                                        Angle bet2, Angle omg2,
                                        Angle& alp1, Angle& alp2) const {
    vec3 r1, r2, v1, v2;
    elliptocart2(bet1, omg1, r1);
    elliptocart2(bet2, omg2, r2);
    real s12 = EuclideanInverse(r1, r2, v1, v2);
    cart2toellip(bet1, omg1, v1, alp1);
    cart2toellip(bet2, omg2, v2, alp2);
    return s12;
  }

  Math::real Triaxial::EuclideanInverse(vec3 r1, vec3 r2,
                                        vec3& v1, vec3& v2) const {
    vec3 dr{r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]};
    real s12 = Math::hypot3(dr[0], dr[1], dr[2]);
    if (s12 == 0) {
      ang bet, omg;
      cart2toellip(r1, bet, omg);
      // head north/south in southern/northern hemisphere
      // east on the equator
      ang alp = ang::cardinal(2 * (bet.s() > 0 ? 1 : 0) +
                              1 * (bet.s() == 0 ? 1 : 0));
      elliptocart2(bet, omg, alp, dr, v1);
      v2 = v1;
    } else {
      v1 = v2 = dr;
      Norm(r1, v1);
      Norm(r2, v2);
    }
    return s12;
  }

  pair<Math::real, Math::real>
  Triaxial::EuclideanDiff(Angle bet1, Angle omg1, Angle alp1,
                          Angle bet2, Angle omg2, Angle alp2) const {
    vec3 r1, v1, r2, v2;
    elliptocart2(bet1, omg1, alp1, r1, v1);
    elliptocart2(bet2, omg2, alp2, r2, v2);
    return EuclideanDiff(r1, v1, r2, v2);
  }

  pair<Math::real, Math::real>
  Triaxial::EuclideanDiff(vec3 r1, vec3 v1, vec3 r2, vec3 v2) {
    return pair<real, real>
      (Math::hypot3(r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]),
       Math::hypot3(v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]));
  }

} // namespace GeographicLib

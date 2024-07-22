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
    , _newumb(false)
    , _debug(false)
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
    real ksum = _k2 + _kp2;
    if (! (isfinite(_a) && isfinite(_b) && isfinite(_c) &&
           _a >= _b && _b >= _c && _c >= 0 && _b > 0 &&
           fabs(ksum - 1) <= numeric_limits<real>::epsilon()) )
      throw GeographicErr("Bad semiaxes for triaxial ellipsoid");
  }

  Triaxial::Triaxial(Math::real b, Math::real e2,
                     Math::real k2, Math::real kp2)
    : _b(b)
    , _e2(e2)
    , _k2(k2)
    , _kp2(kp2)
    , _k(sqrt(_k2))
    , _kp(sqrt(_kp2))
    , _umbalt(false)
    , _newumb(false)
    , _debug(false)
  {
    _a = _b * sqrt(1 + _e2 * _kp2);
    _c = _b * sqrt(1 - _e2 * _k2);
    _axes = vec3({_a, _b, _c});
    real ksum = _k2 + _kp2;
    if (! (isfinite(_a) && isfinite(_b) && isfinite(_c) &&
           _a >= _b && _b >= _c && _c >= 0 && _b > 0 &&
           fabs(ksum - 1) <= numeric_limits<real>::epsilon()) )
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

  void Triaxial:: cart2toellip(Angle bet, Angle omg,
                               vec3 v, Angle& alp) const {
    real tz = hypot(_k, _kp * omg.s()),
      tx = hypot(_k * bet.c(), _kp);
    if (!(bet.c() == 0 && omg.s() == 0)) {
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
    ang bet10(bet1), omg10(omg1);
    bool flip1 = AngNorm(bet1, omg1), flip2 = AngNorm(bet2, omg2);
    bool umb1 = bet1.c() == 0 && omg1.s() == 0,
      umb2 = bet2.c() == 0 && omg2.s() == 0;
    bool swap12;
    {
      ang tmp1(bet1), tmp2(bet2);
      tmp1.setquadrant(0U); tmp2.setquadrant(0U);
      ang tmp12 = tmp2 - tmp1; // |bet2| - |bet1|
      swap12 = tmp12.s() > 0; // is |bet2| > |bet1|
      if (tmp12.s() == 0) {
        tmp1 = omg1; tmp2 = omg2;
        tmp1.setquadrant(0U); tmp2.setquadrant(0U);
        tmp12 = tmp2 - tmp1;
        swap12 = tmp12.s() < 0; // is |omg2| < |omg1|
      }
      // N.B. No swapping if bet1 = +0 and bet1 = -0.
    }
    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
      swap(umb1, umb2);
    }
    // Now |bet1| > |bet2|
    bool flipz = bet1.s() > 0;
    if (flipz) {
      bet1.reflect(true);
      bet2.reflect(true);
    }
    bool flipy = omg1.s() < 0 || (omg1.s() == 0 && omg2.s() < 0);
    if (flipy) {
      omg1.reflect(true);
      omg2.reflect(true);
    }
    bool flipx = signbit(omg1.c());
    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
    }
    bool flipomg = bet2.c() == 0 && signbit(omg2.s());
    // Eliminate coordinate ambiguity bet2 = +/90 and omg2 < 0 for point 2
    // Point 1 is already done with flipy
    if (flipomg) omg2.reflect(true);

    // Now bet1 <= 0, bet1 <= bet2 <= -bet1, 0 <= omg1 <= 90
    //
    // Distinguish two general categories of solution
    //
    // A The geodesic follows one of the principal ellipse (this means, of
    //   course, that both points need to be on the same principal ellipse) and
    //   the initial and final azimuth are easily determined.  The subcases
    //   are:
    //
    //   a Minor ellipse, omg = +/-90: the ellipse is followed regardless of
    //     the points; however we treat this using B (and the solution along
    //     the ellipse is found on the first iteration).
    //
    //   b Major ellipse, bet = 0: the ellipse is followed provided that both
    //     points are close enough "equatorial geodesics".  From point 1 we
    //     follow an equatorial geodesic +/- 180deg in arc distance (= psi
    //     variable).  If point 2 is within the resulting span of longitudes,
    //     the geodesic is equatorial (compute with ArcPos0).  Other do search
    //     iwith strategy B.
    //
    //   c Middle ellipse, bet = +/-90 or omg = 0,180: the ellipse is followed
    //     provided that both points are close enough "meridional geodesics".
    //     The four umbilical points divide the middle ellipse into 4 segments.
    //     There are several subcases:
    //
    //     1 opposite umbilical points: multiple shortest geodesic, two of
    //       which follow the middle ellipse.  But it's more useful to return
    //       the one going through bet = 0, omg = 90.  Then all the others can
    //       be generated by multipling tan(alp1) and tan(alp2) by a constant.
    //
    //     2 bet1 = -90, bet2 = -90: geodesic follows ellipse (compute with
    //       ArcPos0).  Points may be adjacent umbilical points.
    //
    //     3 bet1 = -90, bet2 = 90: geodesic may follow the ellipse if if they
    //       are close enough.  See case A.b for strategy.  Points may be
    //       adjacent umbilical points.
    //
    //     4 |bet2| != 90: geodesic follows the ellipse.  Compute with Hybrid.
    //
    // B The geodesic does not follow a principal ellipse
    //
    //   a point 1 is an umbilical point, gam = 0, so we know alp2
    //
    //   b bet1 = -90: do a search with alp1 bracket = [-90+eps, 90-eps].

    // Set up variables for search
    real fa, fb;
    ang alpa, alpb;

    // and for the final result
    ang bet2a, omg2a;

    TL::fline lf;
    TL::fline::fics fic;
    TL::fline::disttx d;

    // flag for progress
    bool done = false, backside = false;
    if (bet1.c() * omg1.s() == 0 && bet2.c() * omg2.s() == 0) {
      // both points on middle ellipse
      lf = TL::fline(*this, gamblk{}, 0.5, 1.5);
      if (umb1 && umb2 && bet2.s() > 0 && omg2.c() < 0) {
        // process opposite umbilical points
        fic = TL::fline::fics(lf, bet1, omg1, ang{_kp, _k});
        alp1 = ang(_kp * exp(lf.df), _k);
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2, true);
        if (_debug) msg = "opposite umbilics";
        backside = signbit(bet2a.c());
        done = true;
      } else if (bet1.c() == 0 && bet2.c() == 0) {
        // bet1 = -90, bet2 = +/-90
        if (bet2.s() < 1) {
          // bet1 = bet2 = -90
          alp1 = ang::cardinal(1);
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          ang omg12 = omg2 - omg1;
          if (omg12.s() == 0 && omg12.c() < 0) {
            // adjacent E/W umbilical points
            // Should be able to get ArcPos0 to return this?
            d = TL::fline::disttx{-BigValue(), BigValue(), 0 };
            if (_debug) msg = "adjacent EW umbilics";
            alp2 = ang::cardinal(0);
          } else {
            d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
            if (_debug) msg = "bet1/2 = -90";
          }
          if (omg2a.s() < 0) alp2.reflect(true); // Is this needed?
          done = true;
        } else {
          // bet1 = -90, bet2 = 90
          // need to see how far apart the points are
          // If point 1 is at [-90, 0], direction is 0 else -90.
          alp1 = ang::cardinal(omg1.s() == 0 ? 0 : -1);
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          // If point 1 is [-90, 0] and point 2 is [90, 0]
          if (omg1.s() == 0 && omg2.s() == 0) {
            // adjacent N/S umbilical points
            // Should be able to get ArcPos0 to return this?
            d = TL::fline::disttx{BigValue(), -BigValue(), 0};
            alp2 = ang::cardinal(1);
            if (_debug) msg = "adjacent NS umbilics";
            done = true;
          } else {
            d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
            omg2a -= omg2;
            if (omg2a.s() > 0) {
              ang omg12 = omg2 + omg1;
              d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
              if (omg2a.s() < 0) alp2.reflect(true); // Is this needed?
              if (_debug) msg = "bet1/2 = -/+90 meridional";
              done = true;
            } else {
              alpa = ang::cardinal(-1) + ang::eps();
              fa = omg2a.radians0();
              alpb.setquadrant(0U);
              fic.setquadrant(lf, 0U);
              (void) lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
              omg2a -= omg2;
              alpb = -alpa;
              fb = omg2a.radians0();
              if (_debug) msg = "general bet1/2 = -/+90, non-meridional";
            }
          }
        }
      } else {
        // other meridional cases, invoke Hybrid with the following value of
        // alp1
        alp1 = ang::cardinal(bet1.c() == 0 ?
                             (omg2.c() < 1 ? 1 : (omg1.s() == 0 ? 0 : -1)) :
                             (omg2.c() > 0 ? 0 : 2));
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
        if (_debug) msg = "other meridional";
        done = true;
      }
    } else if (bet1.s() == 0 && bet2.s() == 0) {
      // both points on equator
      ang omg12 = (omg2 - omg1).base();
      int eE = omg12.s() > 0 ? 1 : -1;
      // set direction for probe as +/-90 based on sign of omg12
      alp1 = ang::cardinal(eE);
      bet1.reflect(true);
      lf = TL::fline(*this, gamma(bet1, omg1, alp1), 0.5, 1.5);
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      (void) lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
      omg2a -= omg2;
      if (eE * omg2a.s() >= 0) {
        // geodesic follows the equator
        d = lf.ArcPos0(fic, omg12.flipsign(eE), bet2a, omg2a, alp2, false);
        if (_debug) msg = "bet1/2 = 0 equatorial";
        done = true;
      } else {
        // geodesic does not follow the equator
        alpb = ang::cardinal(-1) - ang::eps();
        alpa = -alpb;
        (eE > 0 ? fa : fb) = omg2a.radians0();
        alp1.setquadrant(eE > 0 ? 3U : 0U);
        fic.setquadrant(lf, eE > 0 ? 3U : 0U);
        (void) lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
        omg2a -= omg2;
        (eE > 0 ? fb : fa) = omg2a.radians0();
        if (_debug) msg = "general bet1/2 = 0 non-equatorial";
      }
    } else if (umb1) {
      // umbilical point to general point
      lf = TL::fline(*this, gamblk{}, 0.5, 1.5);
      alp2 = ang(_kp * omg2.s(), _k * bet2.c());
      fic = TL::fline::fics(lf, bet2, omg2, alp2);
      (void) lf.ArcPos0(fic, bet1 - bet2, bet2a, omg2a, alp1);
      if (alp1.s() < 0) alp1 += ang::cardinal(1);
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      d = lf.ArcPos0(fic, bet2 - bet1, bet2a, omg2a, alp2);
      if (_debug) msg = "umbilic to general";
      done = true;
    } else if (bet1.c() == 0) {
      // bet1 = -90 to general point
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
      if (_debug) msg = "general bet = -90";
    } else if (omg1.s() == 0) {
      // omg1 = 0 to general point
      if (omg2.s() > 0) {
        alpa = ang::eps();
        alpb = ang::cardinal(2) - alpa;
        fa = -omg2.radians();
        fb = (ang::cardinal(2)-omg2).radians0();
      } else {
        alpa = ang(-numeric_limits<real>::epsilon()/(1<<20), -1);
        alpb = ang(-numeric_limits<real>::epsilon()/(1<<20),  1);
        fa = (ang::cardinal(2)-omg2).radians0();
        fb = -omg2.radians();
      }
      if (_debug) msg = "general omg1 = 0";
    } else {
      // general case
      real f[4];
      alpa = ang( _kp * fabs(omg1.s()), _k * fabs(bet1.c()));
      alpb = alpa;

      lf = TL::fline(*this, gamblk{}, 0.5, 1.5);
      fic = TL::fline::fics(lf, bet1, omg1, alpb);
      unsigned qb = 0U, qa = 3U; // qa = qb - 1 (mod 4)
      for (; !done && qb <= 4U; ++qb, ++qa) {
        if (qb) {
          alpb.setquadrant(qb);
          fic.setquadrant(lf, qb);
        }
        if (qb < 4U) {
          f[qb] = lf.Hybrid0(fic, bet2, omg2);
          if (fabs(f[qb]) < numeric_limits<real>::epsilon()) {
            alp1 = alpb;
            d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
            if (_debug) msg = "accidental umbilic";
            backside = signbit(bet2a.c()); // qb == 1U || qb == 2U;
            done = true;
            break;
          }
        }
        if (qb && (f[qa & 3U] < 0 && f[qb & 3U] > 0)) {
          break;
        }
      }
      if (!done) {
        fa = f[qa & 3U]; fb = f[qb & 3U];
        alpa.setquadrant(qa);
      }
      if (_debug) msg = "general";
    }

    int countn = 0, countb = 0;
    if (!done) {
      alp1 = findroot(
                      [this, &bet1, &omg1, &bet2, &omg2]
                      (const ang& alp) -> Math::real
                      {
                        return HybridA(*this,
                                       bet1, omg1, alp,
                                       bet2, omg2);
                      },
                      alpa,  alpb,
                      fa, fb,
                      &countn, &countb);
      alp1.round();
      lf = TL::fline(*this, gamma(bet1, omg1, alp1), 0.5, 1.5);
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
    cerr << "A " << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
               << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    // Undo switches in reverse order flipz, swap12, flip1
    if (flipomg) {
      omg2.reflect(true);
      alp2.reflect(true, true);
    }

    if (_debug)
    cerr << "B " << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
               << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
      alp1.reflect(true);
      alp2.reflect(true);
    }

    if (_debug)
    cerr << "C " << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
               << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipy) {
      omg1.reflect(true);
      omg2.reflect(true, true);
      alp1.reflect(true);
      alp2.reflect(true);
    }

    if (_debug)
    cerr << "D " << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
               << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flipz) {
      bet1.reflect(true);
      bet2.reflect(true);
      alp1.reflect(false, true);
      alp2.reflect(false, true);
    }

    if (_debug)
    cerr << "E " << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
         << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
      swap(alp1, alp2);
      swap(umb1, umb2);
      // points not swapped if umb1 == true
      alp1 += ang::cardinal(2);
      if (umb2)
        alp2 += ang::cardinal((signbit(alp2.s()) ? -1 : 1) * bet2.s());
      else
        alp2 += ang::cardinal(2);
    }

    if (_debug)
    cerr << "F " << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
         << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    if (flip1) Flip(bet1, omg1, alp1);
    if (flip2) Flip(bet2, omg2, alp2);
    alp1.setn(); alp2.setn();
    if (_debug)
    cerr << "G " << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
         << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";

    fic = TL::fline::fics(lf, bet10, omg10, alp1);
    gic = TL::gline::gics(lg, fic);
    gic.s13 = fmax(real(0), s12);

    if (_debug)
      cerr << countn << " " << countb << " "
           << lf.gamma() << " "
           << lf.fbet().NCoeffs() << " " << lf.fomg().NCoeffs() << " "
           << lg.gbet().NCoeffs() << " " << lg.gomg().NCoeffs() << " "
           << msg << "\n";
    return TL(move(lf), move(fic), move(lg), move(gic));
  }

  Math::real Triaxial::HybridA(const Triaxial& t,
                               Angle bet1, Angle omg1,
                               Angle alp1,
                               Angle bet2, Angle omg2) {
    ang b1{bet1}, o1{omg1}, a1{alp1};
    gamblk gam = t.gamma(b1, o1, a1);
    TriaxialLine::fline l(t, gam, 0.5, 1.5);
    TriaxialLine::fline::fics ic(l, b1, o1, a1);
    return l.Hybrid0(ic, bet2, omg2);
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
    bool trip = false, correct = false;
    for (Math::real t = 1/Math::real(2), ab = 0, ft = 0, fm = 0, fc = 0;
         cntn < maxcnt ||
           (throw GeographicLib::GeographicErr
            ("Convergence failure Triaxial::findroot"), false)
           || GEOGRAPHICLIB_PANIC("Convergence failure Triaxial::findroot");) {
      // These inverse problems use lots of iterations
      //  22  48  90   1 -48.5628 -5.7915 0.7706
      //  56 115 -89 179 113.5952 179.8512 1.6130
      // -51  89  90   1 -65.9888 10.0598 2.0530
      // Need to figure out why.
      ang xt = 2*t == 1 ?
        ang(xa.s() + xb.s(), xa.c() + xb.c()) :
        (2*t < 1 ? xa - ang::radians(t * ab) :
         xb + ang::radians((1 - t) * ab)),
        xc;
      if (trip) {
        if (correct) xm = xt;   // Update result if not a bisection step
        break;
      }
      ++cntn;
      ft = f(xt);
      if (signbit(ft) == signbit(fa)) {
        xc = xa; xa = xt;
        fc = fa; fa = ft;
      } else {
        xc = xb; xb = xa; xa = xt;
        fc = fb; fb = fa; fa = ft;
      }
      if (fabs(fb) < fabs(fa)) {
        xm = xb; fm = fb;
      } else {
        xm = xa; fm = fa;
      }
      // ordering is b - a - c
      ab = (xa-xb).radians0();
      Math::real
        ca = (xc-xa).radians0(),
        cb = ca+ab,
        // Scherer has a fabs(cb).  This should be fabs(ab).
        tl = numeric_limits<Math::real>::epsilon() / fabs(ab);
      // Backward tests to deal with NaNs
      trip =  !(2 * tl < 1 &&
                fabs(fm) > numeric_limits<Math::real>::epsilon());
      // If trip update xm one more time, then Hybrid solution is called once
      // more outside this route to update bet2, omg2, alp2, etc.
      if (trip) tl = 0;
      Math::real
        xi = ab / cb,
        phi = (fa-fb) / (fc-fb);
      if ( 2 * tl < 1 && xi / (1 + sqrt(1 - xi)) < phi && phi < sqrt(xi) ) {
        t = fa/(fb-fa) * fc/(fb-fc) - ca/ab * fa/(fc-fa) * fb/(fc-fb);
        // This equation matches the pseudocode in Scherer.  His Eq (6.40)
        // reads t = fa/(fb-fa) * fc/(fb-fc) + ca/cb * fc/(fc-fa) * fb/(fb-fa);
        // this is wrong.
        t = fmin(1 - tl, fmax(tl, t));
        correct = true;
      } else {
        t = 1/Math::real(2);
        ++cntb;
        correct = false;
      }
    }
    if (countn) *countn += cntn;
    if (countb) *countb += cntb;
    return xm;
  }

  Triaxial::gamblk Triaxial::gamma(Angle bet, Angle omg, Angle alp) const {
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
    // cout << "GAMDIFF " << gam/ numeric_limits<real>::epsilon() << " " << maxdiff << "\n";
    if (fabs(gam) <= 2 * maxdiff * numeric_limits<real>::epsilon())
      // Set gam = 0 if a change of alp, bet, or omg by epsilon would include
      // gam = 0.
      gam = 0;
    real gamp = gam == 0 ? 0 :
      (gam > 0 ? // sqrt(k2 - gamma)
       hypot(_k * hypot(bet.s(), alp.c()*bet.c()),
             _kp * omg.s()*alp.c()) :
       // sqrt(kp2 + gamma)
       hypot(_k *  bet.c()*alp.s(),
             _kp * hypot(omg.c(), alp.s()*omg.s()))),
      // for gam == 0, we have nu = nup = 0
      nu = sqrt(fabs(gam)) / (gam > 0 ? _k : _kp),
      nup = gamp / (gam > 0 ? _k : _kp);
    return gamblk{gam, nu, nup};
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
  Triaxial::EuclideanDiff(vec3 r1, vec3 v1, vec3 r2, vec3 v2) const {
    return pair<real, real>
      (Math::hypot3(r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]),
       Math::hypot3(v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]));
  }

} // namespace GeographicLib

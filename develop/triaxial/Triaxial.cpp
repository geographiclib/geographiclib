/**
 * \file Triaxial.cpp
 * \brief Implementation for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "Triaxial.hpp"

namespace GeographicLib {

  using namespace std;

  Triaxial::Triaxial()
    : Triaxial(1, 0, 1, 0)
  {}

  Triaxial::Triaxial(real a, real b, real c)
    : _a(a)
    , _b(b)
    , _c(c)
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
    _oblate = _kp2 == 0;
    _prolate = _k2 == 0;
    _biaxial = _oblate || _prolate;
  }

  Triaxial::Triaxial(real b, real e2, real k2, real kp2)
    : _b(b)
    , _e2(e2)
    , _k2(k2)
    , _kp2(kp2)
  {
    real ksum = _k2 + _kp2;
    _k2 /= ksum;
    _kp2 /= ksum;
    _k = sqrt(_k2);
    _kp = sqrt(_kp2);
    _a = _b * sqrt(1 + _e2 * _kp2);
    _c = _b * sqrt(1 - _e2 * _k2);
    if (! (isfinite(_a) && isfinite(_b) && isfinite(_c) &&
           _a >= _b && _b >= _c && _c >= 0 && _b > 0) )
      throw GeographicErr("Bad semiaxes for triaxial ellipsoid");
    _oblate = _kp2 == 0;
    _prolate = _k2 == 0;
    _biaxial = _oblate || _prolate;
  }

  void Triaxial::Norm(vec3& r) const {
    vec3 rn{r[0] / _a, r[1] / _b, r[2] / _c};
    real ra = Math::hypot3(rn[0], rn[1], rn[2]);
    r = {r[0] / ra, r[1] / ra, r[2] / ra};
  }

  void Triaxial::Norm(vec3& r, vec3& v) const {
    Norm(r);
    vec3 up = {r[0] / Math::sq(_a), r[1] / Math::sq(_b), r[2] / Math::sq(_c)};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      uv = up[0] * v[0] + up[1] * v[1] + up[2] * v[2],
      f = uv/u2;
    v = {v[0] - f * up[0], v[1] - f * up[1], v[2] - f * up[2]};
    normvec(v);
  }

  void Triaxial::cart2toellipint(vec3 r, Angle& bet, Angle& omg, vec3 axes)
    const {
    real a = axes[0], b = axes[1], c = axes[2];
    real xi = r[0]/a, eta = r[1]/b, zeta = r[2]/c,
      g = _k2 * Math::sq(xi)
      + (_k2 - _kp2) * Math::sq(eta)
      - _kp2 * Math::sq(zeta);
    if (fabs(r[0]) == a * _kp2 && r[1] == 0 && fabs(r[2]) == c * _k2)
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
    bet = ang(sb, cb, 0, true);
    omg = ang(so, co, 0, true);
  }

  void Triaxial::cart2toellip(vec3 r, Angle& bet, Angle& omg) const {
    cart2toellipint(r, bet, omg, {_a, _b, _c});
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
    r = { _a * omg.c() * tx,
          _b * bet.c() * omg.s(),
          _c * bet.s() * tz };
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
      v = {_a*_k/_b * omg.c() * ca2,
           -omg.c() * bet.s() * sa2,
           -_c*_kp/_b * bet.s() * ca2};
    } else {
      vec3 N, E;
      if (tx == 0) {
        // At an oblate pole tx -> |cos(bet)|
        real scb = signbit(bet.c()) ? -1 : 1;
        N = {-omg.c() * bet.s() * scb, -omg.s() * bet.s(), 0};
        E = {-omg.s()                ,  omg.c() * scb    , 0};
      } else if (tz == 0) {
        // At a prolate pole tz -> |sin(omg)|
        real sso = signbit(omg.s()) ? -1 : 1;
        N = {0, -bet.s() * sso    , bet.c()                };
        E = {0,  bet.c() * omg.c(), bet.s() * omg.c() * sso};
      } else {
        // The general case
        N = { -_a * _k2 * bet.c() * bet.s() * omg.c() / tx,
              -_b * bet.s() * omg.s(),
               _c * bet.c() * tz};
        E = { -_a * tx * omg.s(),
               _b * bet.c() * omg.c(),
               _c * _kp2 * bet.s() * omg.c() * omg.s() / tz};
      }
      normvec(N);
      normvec(E);
      v = {alp.c() * N[0] + alp.s() * E[0],
           alp.c() * N[1] + alp.s() * E[1],
           alp.c() * N[2] + alp.s() * E[2]};
    }
    // normvec(v); v is already normalized
  }

} // namespace GeographicLib

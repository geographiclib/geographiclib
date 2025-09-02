/**
 * \file Triaxial.cpp
 * \brief Implementation for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
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

  Triaxial::Triaxial(real a, real b, real c)
    : _a(a)
    , _b(b)
    , _c(c)
    , _axes{_a, _b, _c}
    , _axes2{Math::sq(_a), Math::sq(_b), Math::sq(_c)}
    , _linecc2{(_a - _c) * (_a + _c), (_b - _c) * (_b + _c), 0}
    , _umbalt(false)
    , _biaxp(true)
    , _debug(false)
    , _hybridalt(true)
    , _swapomg(false)
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
    if (_k2 == 0) _umbalt = true;
    _oblate = _kp2 == 0;
    _prolate = _k2 == 0;
    _biaxial = _oblate || _prolate;
  }

  Triaxial::Triaxial(real b, real e2, real k2, real kp2)
    : _b(b)
    , _e2(e2)
    , _k2(k2)
    , _kp2(kp2)
    , _umbalt(false)
    , _biaxp(true)
    , _debug(false)
    , _hybridalt(true)
    , _swapomg(false)
    , _ellipthresh(1/real(8))
  {
    real ksum = _k2 + _kp2;
    _k2 /= ksum;
    _kp2 /= ksum;
    _k = sqrt(_k2);
    _kp = sqrt(_kp2);
    _a = _b * sqrt(1 + _e2 * _kp2);
    _c = _b * sqrt(1 - _e2 * _k2);
    _axes = {_a, _b, _c};
    _axes2 = {Math::sq(_a), Math::sq(_b), Math::sq(_c)};
    _linecc2 = {(_a - _c) * (_a + _c), (_b - _c) * (_b + _c), 0};
    if (! (isfinite(_a) && isfinite(_b) && isfinite(_c) &&
           _a >= _b && _b >= _c && _c >= 0 && _b > 0) )
      throw GeographicErr("Bad semiaxes for triaxial ellipsoid");
    if (_k2 == 0) _umbalt = true;
    _oblate = _kp2 == 0;
    _prolate = _k2 == 0;
    _biaxial = _oblate || _prolate;
  }

  void Triaxial::Norm(vec3& r) const {
    vec3 rn{r[0] / _axes[0], r[1] / _axes[1], r[2] / _axes[2]};
    real ra = Math::hypot3(rn[0], rn[1], rn[2]);
    r = {r[0] / ra, r[1] / ra, r[2] / ra};
  }

  void Triaxial::Norm(vec3& r, vec3& v) const {
    Norm(r);
    vec3 up = {r[0] / _axes2[0], r[1] / _axes2[1], r[2] / _axes2[2]};
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
    cart2toellipint(r, bet, omg, _axes);
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
    r = { _axes[0] * omg.c() * tx,
          _axes[1] * bet.c() * omg.s(),
          _axes[2] * bet.s() * tz };
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

  template<int n>
  pair<Math::real, Math::real> Triaxial::funp<n>::operator()(real p) const {
    real fp = 0, fv = -1, fcorr = 0;
    for (int k = 0; k < 3; ++k) {
      if (_r[k] == 0) continue;
      real g = _r[k] / (p + _l[k]);
      if constexpr (n == 2) g *= g;
      real ga = round(g/_d) * _d, gb = g - ga;
      fv = fv + ga; fcorr = fcorr + gb;
      fp = fp - n * g / (p + _l[k]);
    }
    return pair<real, real>(fv + fcorr, fp);
  }

  Math::real Triaxial::cartsolve(const function<pair<real, real>(real)>& f,
                                 real p0, real pscale) {
    // Solve
    //   f(p) = 0
    // Initial guess is p0; pscale is a scale factor for p; scale factor for f
    // = 1.  This assumes that there's a single solution with p >= 0 and that
    // for p > 0, f' < 0 and f'' > 0
    const real eps = numeric_limits<real>::epsilon(),
      tol = Math::sq(cbrt(eps)), tol2 = Math::sq(tol),
      ptol = pscale * sqrt(eps);
    real p = p0;
    int i = 50;
    real od = -1;
    while (i > 0) {
      pair<real, real> fx = f(p);
      real fv = fx.first, fp = fx.second;
      // We're done if f(p) <= 0 on initial guess; this can happens when z = 0.
      // However, since Newton converges from below, any negative f(p)
      // indicates convergence.
      if (!(fv > tol2)) break;
      real d = -fv/fp;          // d is positive
      p = p + d;
      // converged if fv <= 8*eps (after first iteration) or
      // d <= max(eps, |p|) * tol and d <= od.
      // N.B. d is always positive.
      if ( (fv <= 8 * eps || d <= fmax(ptol, p) * tol) && d <= od )
        // The condition d <= od means that this won't trip on the first
        // iteration
        break;
      od = d;
    }
    return p;
  }

  Math::real Triaxial::cubic(vec3 r2) const {
    // Solve sum(r2[i]/(z + lineq2[i]), i,0,2) - 1 = 0 with lineq2[2] = 0.
    // This has three real roots with just one satisifying q >= 0.
    // Express as a cubic equation z^3 + a*z^2 + b*z + c = 0.
    real c = - _linecc2[0]*_linecc2[1] * r2[2],
      b = _linecc2[0]*_linecc2[1]
      - (_linecc2[1] * r2[0] + _linecc2[0] * r2[1] +
         (_linecc2[0] + _linecc2[1]) * r2[2]),
      a = _linecc2[0] + _linecc2[1] - (r2[0] + r2[1] + r2[2]);
    bool recip = b > 0;
    if (recip) {
      // If b positive there a cancellation in p = (3*b - a^2) / 3, so
      // transform to a polynomial in 1/t.  The resulting coefficients are
      real ax = b/c, bx = a/c, cx = 1/c;
      a = ax; b = bx; c = cx;
    }
    // Reduce cubic to w^3 + p*w + q = 0, where z = w - a/3.
    // See https://dlmf.nist.gov/1.11#iii
    real p = (3*b - Math::sq(a)) / 3,
      q = (2*a*Math::sq(a) - 9*a*b + 27*c) / 27;

    // Now switch to https://dlmf.nist.gov/4.43
    // We have 3 real roots, so 4*p^3 + 27*q^2 <= 0
    real A = sqrt(fmax(real(0), -4*p/3)),
      alp = atan2(q, sqrt(fmax(real(0),
                               -(4*p*Math::sq(p)/27 + Math::sq(q)))))/3;
    // alp is in [-pi/3, pi/3]
    // z = A*sin(alp + 2*pi/3 * k) - a/3 for k = -1, 0, 1
    // for the single positive solution we pick k = 1 which gives the
    // algebraically largest result
    real t = A/2 * (cos(alp) * sqrt(real(3)) - sin(alp)) - a/3;

    return recip ? 1/t : t;
  }

  void Triaxial::carttoellip(vec3 r, Angle& bet, Angle& omg, real& H) const {
    // tol2 = eps^(4/3)
    real tol2 = Math::sq(Math::sq(cbrt(numeric_limits<real>::epsilon())));
    vec3 r2 = {Math::sq(r[0]), Math::sq(r[1]), Math::sq(r[2])};
    real qmax = r2[0] + r2[1] + r2[2],
      qmin = fmax(fmax(r2[2], r2[1] + r2[2] - _linecc2[1]),
                  r2[0] + r2[1] + r2[2] - _linecc2[0]),
      q = qmin;
    do {                       // Executed once (provides the ability to break)
      const funp<1> f(r2, _linecc2);
      pair<real, real> fx = f(q);
      if (!( fx.first > tol2 ))
        break;                  // negative means converged
      q = fmax(qmin, fmin(qmax, cubic(r2)));
      fx = f(q);
      if (!( fabs(fx.first) > tol2 ))
        break;                  // test abs(fv) here
      q = fmax(qmin, q - fx.first/fx.second);
      q = cartsolve(f, q, Math::sq(_b));
    } while (false);
    vec3 axes = {sqrt(_linecc2[0] + q), sqrt(_linecc2[1] + q), sqrt(q)};
    cart2toellipint(r, bet, omg, axes);
    H = axes[2] - _c;
  }

  void Triaxial::elliptocart(Angle bet, Angle omg, real H, vec3& r) const {
    vec3 ax;
    real shift = H * (2*_c + H);
    for (int k = 0; k < 2; ++k)
      ax[k] = sqrt(_axes2[k] + shift);
    ax[2] = _c + H;
    real tx = hypot(_k * bet.c(), _kp), tz = hypot(_k, _kp * omg.s());
    r = { ax[0] * omg.c() * tx,
          ax[1] * bet.c() * omg.s(),
          ax[2] * bet.s() * tz };
  }

  template<int n>
  void Triaxial::cart2togeneric(vec3 r, Angle& phi, Angle& lam) const {
    static_assert(n >= 0 && n <= 2, "Bad coordinate conversion");
    if constexpr (n == 2) {
      r[0] /= _axes2[0];
      r[1] /= _axes2[1];
      r[2] /= _axes2[2];
    } else if constexpr (n == 1) {
      r[0] /= _axes[0];
      r[1] /= _axes[1];
      r[2] /= _axes[2];
    } // else n == 0, r is unchanged
    phi = ang(r[2], hypot(r[0], r[1]));
    lam = ang(r[1], r[0]);
  }

  template<int n>
  void Triaxial::generictocart2(Angle phi, Angle lam, vec3& r) const {
    static_assert(n >= 0 && n <= 2, "Bad coordinate conversion");
    r = {phi.c() * lam.c(),
      phi.c() * lam.s(),
      phi.s()};
    if constexpr (n == 2) {
      r[0] *= _axes2[0];
      r[1] *= _axes2[1];
      r[2] *= _axes2[2];
    } else if constexpr (n == 1) {
      r[0] *= _axes[0];
      r[1] *= _axes[1];
      r[2] *= _axes[2];
    } // else n == 0, r is unchanged
    if constexpr (n != 1) {
      real d = Math::hypot3(r[0] / _axes[0], r[1] / _axes[1], r[2] / _axes[2]);
      r[0] /= d; r[1] /= d; r[2] /= d;
    } // else n == 1, d = 1 and r is already on the surface of the ellipsoid
  }

  void Triaxial::cart2togeod(vec3 r, Angle& phi, Angle& lam) const {
    cart2togeneric<2>(r, phi, lam);
  }
  void Triaxial::geodtocart2(Angle phi, Angle lam, vec3& r) const {
    generictocart2<2>(phi, lam, r);
  }

  void Triaxial::cart2toparam(vec3 r, Angle& phip, Angle& lamp) const {
    cart2togeneric<1>(r, phip, lamp);
  }
  void Triaxial::paramtocart2(Angle phip, Angle lamp, vec3& r) const {
    generictocart2<1>(phip, lamp, r);
  }

  void Triaxial::cart2togeocen(vec3 r, Angle& phipp, Angle& lampp) const {
    cart2togeneric<0>(r, phipp, lampp);
  }
  void Triaxial::geocentocart2(Angle phipp, Angle lampp, vec3& r) const {
    generictocart2<0>(phipp, lampp, r);
  }

  void Triaxial::cart2tocart(vec3 r2, real h, vec3& r) const {
    vec3 rn = {r2[0] / _axes2[0], r2[1] / _axes2[1], r2[2] / _axes2[2]};
    real d = h / Math::hypot3(rn[0], rn[1], rn[2]);
    r = r2;
    for (int k = 0; k < 3; ++k)
      r[k] += rn[k] * d;
  }

  void Triaxial::carttocart2(vec3 r, vec3& r2, real& h) const {
    const real eps = numeric_limits<real>::epsilon(), ztol = _b * eps/8;
    for (int k = 0; k < 3; ++k)
      if (fabs(r[k]) <= ztol) r[k] = copysign(real(0), r[k]);
    vec3 s = {r[0] * _axes[0], r[1] * _axes[1], r[2] * _axes[2]};
    real p = fmax(fmax(fabs(s[2]), hypot(s[1], s[2]) - _linecc2[1]),
                  Math::hypot3(s[0], s[1], s[2]) - _linecc2[0]);
    const funp<2> f(s, _linecc2);
    p = cartsolve(f, p, Math::sq(_b));
    r2 = r;
    for (int k = 0; k < 3; ++k)
      r2[k] *= _axes2[k] / (p + _linecc2[k]);
    // Deal with case p == 0 (when r2[2] is indeterminate).
    if (p == 0) {
      if (_linecc2[0] == 0)     // sphere
        r2[0] = r[0];
      if (_linecc2[1] == 0)     // sphere or prolate
        r2[1] = r[1];
      r2[2] = _axes[2] * r[2] *
        sqrt(1 - Math::sq(r2[0]/_axes[0]) + Math::sq(r2[1]/_axes[1]));
    }
    h = (p - _axes2[2]) * Math::hypot3(r2[0] / _axes2[0],
                                       r2[1] / _axes2[1],
                                       r2[2] / _axes2[2]);
  }

  void Triaxial::carttogeod(vec3 r, Angle& phi, Angle& lam, real& h) const {
    vec3 r2;
    carttocart2(r, r2, h);
    cart2togeod(r2, phi, lam);
  }

  void Triaxial::geodtocart(Angle phi, Angle lam, real h, vec3& r) const {
    vec3 r2;
    geodtocart2(phi, lam, r2);
    cart2tocart(r2, h, r);
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

    if (!_umbline)
      // Initialize _umbline
      _umbline = make_shared<TriaxialLine>(TriaxialLine(*this));

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
    bool flip1 = AngNorm(bet1, omg1, _prolate),
      flip2 = AngNorm(bet2, omg2, _prolate);
    bool swap12;
    {
      // For prolate, swap based on omg, switch 1 & 2 because poles are at
      // 0/180, instead of +/-90.
      ang tmp1(_prolate ? omg2 : bet1), tmp2(_prolate ? omg1 : bet2);
      tmp1.setquadrant(0U); tmp2.setquadrant(0U);
      ang tmp12 = tmp2 - tmp1; // |bet2| - |bet1|
      swap12 = tmp12.s() > 0; // is |bet2| > |bet1|
      if (!_biaxial && tmp12.s() == 0) {
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
    if (_oblate) {
      // Rotate, subtracting omg1 from omg[12], so omg1 = 0
      omg2 -= omg1;
      omg2 = omg2.base();
      omg1 = ang::cardinal(0);
    } else if (_prolate) {
      // Rotate, subtracting bet1 + 90 from bet[12], so bet1 = -90
      bet2 -= bet1 + ang::cardinal(1);
      bet2 = bet2.base();
      bet1 = ang::cardinal(-1);
    }

    bool swapomg = _swapomg &&
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
    bool flipy = _prolate ? signbit(bet2.c()) :
      signbit(omg1.s()) || (omg1.s() == 0 && signbit(omg2.s()));
    if (flipy) {
      if (_prolate)
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

    real aa = _k2 * Math::sq(bet2.c()), bb = _kp2 * Math::sq(omg2.s());

    if (_debug)
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
        alp1 = _biaxial ? ang(_k, _kp, 0, true) :
          ang(exp(lf.deltashift()/2), 1);
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        bool betp = _k2 < _kp2;
        //        if (_biaxial) betp = !betp;
        //        betp = !betp;
        d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2, betp);
        if (_debug) msg = "A.c opposite umbilics";
        backside = signbit(bet2a.c());
        done = true;
      } else if (bet1.c() == 0 && bet2.c() == 0) {
        // Case A.c.{2,3}, bet1 = -90, bet2 = +/-90
        if (bet2.s() < 0) {
          // Case A.c.2, bet1 = bet2 = -90
          // If oblate, bet1 = -90, omg1 = 0, need alp1 = omg2 to follow omg2
          // meridian.
          alp1 = ang::cardinal(_oblate ? 2 : 1);
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          ang omg12 = omg2 - omg1;
          if (omg12.s() == 0 && omg12.c() < 0) {
            // adjacent E/W umbilical points
            // Should be able to get ArcPos0 to return this?
            d = _oblate ?
              TL::fline::disttx{ (_biaxial ? 1 : -1 ) * Math::pi()/2,
                                 -Math::pi()/2, 0 } :
              _prolate ?
              TL::fline::disttx{ Math::pi()/2, -Math::pi()/2, 0 } :
              TL::fline::disttx{ -BigValue(), BigValue(), 0 };
            if (_debug) msg = "A.c.2 adjacent EW umbilics";
            alp2 = ang::cardinal(_prolate ? 1 : 0);
          } else {
            d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
            if (_debug) msg = "A.c.2 bet1/2 = -90";
          }
          // not needed
          // if (omg2a.s() < 0) alp2.reflect(true);
          done = true;
        } else {
          // Case A.c.3, bet1 = -90, bet2 = 90
          // need to see how far apart the points are
          // If point 1 is at [-90, 0], direction is 0 else -90.
          // For prolate use -90.
          alp1 = ang::cardinal(omg1.s() == 0 && !_prolate ? 0 : -1);
          if (0)
            cout << "FIC " << real(bet1) << " " << real(omg1) << " "
                 << real(alp1) << "\n";
          fic = TL::fline::fics(lf, bet1, omg1, alp1);
          // If point 1 is [-90, 0] and point 2 is [90, 0]
          if (omg1.s() == 0 && omg2.s() == 0) {
            // adjacent N/S umbilical points
            // Should be able to get ArcPos0 to return this?
            d = _biaxial ?
              TL::fline::disttx{ Math::pi()/2, -Math::pi()/2, 0 } :
              TL::fline::disttx{ BigValue(), -BigValue(), 0 };
            alp2 = ang::cardinal(_oblate ? 0 :  1);
            if (_debug) msg = "A.c.3 adjacent NS umbilics";
            done = true;
          } else {
            // FIX ME for oblate
            if (omg1.s() == 0)
              omg2a = ang::cardinal(2);
            else
              // Compute conjugate point along the middle ellipse
              d = lf.ArcPos0(fic, ang::cardinal(2), bet2a, omg2a, alp2);
            if (_debug)
              cout << "QQ "
                   << real(fic.tht1) << " " << real(fic.phi1) << " "
                   << real(fic.alp1) << " "
                   << real(bet2a) << " " << real(omg2a) << " "
                   << real(alp2) << "\n";
            // XXX FIX HERE for prolate case -90 -1 90 177
            omg2a -= omg2;
            if (omg2a.s() >= -numeric_limits<real>::epsilon()/2) {
              // Includes all cases where point 1 = [-90, 0], omg1.s() == 0.
              ang omg12 = omg2 + omg1;
              // FIX ME for oblate
              d = lf.ArcPos0(fic, omg12.base(), bet2a, omg2a, alp2, false);
              if (_debug)
                cout << "APX "
                     << real(bet1) << " " << real(omg1) << " "
                     << real(alp1) << " "
                     << real(bet2a) << " " << real(omg2a) << " "
                     << real(alp2) << " " << real(omg12.base()) << "\n";
              if (!_biaxial && signbit(omg2a.s()))
                alp2.reflect(true);
              if (_debug) msg = "A.c.3 bet1/2 = -/+90 meridional";
              done = true;
            } else {
              alpa = ang::cardinal(-1) + ang::eps();
              fa = omg2a.radians0();
              (void) lf.ArcPos0(fic, ang::cardinal(-2), bet2a, omg2a, alp2);
              omg2a -= omg2;
              alpb = -alpa;
              fb = omg2a.radians0();
              if (_debug) msg = "A.c.3 general bet1/2 = -/+90, non-meridional";
              done = false;     // A marker
            }
            if (_debug)
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
        alp1 = _oblate ? omg2 :
          (_prolate ? -bet2 :
           ang::cardinal(bet1.c() == 0 ?
                         // TODO: CHECK omg2.c() < 1 test; CHANGE TO < 0
                         (omg2.c() < 0 ? 1 :
                          (omg1.s() == 0 && !_prolate ? 0 : -1)) :
                         (omg2.c() > 0 ? 0 : 2)));
        if (_debug)
          cout << "ALP1 " << real(alp1) << " "
               << bet1.c() << " " << omg2.c() << "\n";
        fic = TL::fline::fics(lf, bet1, omg1, alp1);
        d = _prolate ?
          lf.ArcPos0(fic, (omg2-omg1).base(), bet2a, omg2a, alp2, false) :
          lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
        if (_debug) cout << "AAA " << real(alp2) << " " << real(bet2a) << "\n";
        if (_prolate && signbit(bet2a.c()))
          alp2.reflect(true,true);
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
      if (E * omg2a.s() >= -numeric_limits<real>::epsilon()/2) {
        // geodesic follows the equator
        d = lf.ArcPos0(fic, omg12.flipsign(E), bet2a, omg2a, alp2, false);
        if (_debug) msg = "A.b.1 bet1/2 = 0 equatorial";
        done = true;
      } else {
        // geodesic does not follow the equator
        alpb = ang::cardinal(-1) - ang::eps();
        alpa = -alpb;
        (E > 0 ? fa : fb) = omg2a.radians0();
        (void) lf.ArcPos0(fic, ang::cardinal(-2), bet2a, omg2a, alp2);
        omg2a -= omg2;
        (E > 0 ? fb : fa) = omg2a.radians0();
        if (_debug) msg = "A.b.2 general bet1/2 = 0 non-equatorial";
        done = false;           // A marker
      }
    } else if (umb1) {
      // Case B.a, umbilical point to general point
      alp2 = ang(_kp * omg2.s(), _k * bet2.c());
      // RETHINK THIS.  If we know alp2, we can compute delta.  This should be
      // enough to find alp1.
      if (0)
        cout << "ALP2 " << real(alp2) << " " << real(bet1 - bet2) << "\n";
      fic = TL::fline::fics(lf, bet2, omg2, alp2);
      //      bool betp = _k2 > _kp2;   // This could be betb = !_prolate;
      bool betp = aa > bb;
      // cout << "DELTA " << fic.delta/Math::degree() << "\n";
      real delta = (lf.transpolar() ? -1 : 1) * fic.delta;
      alp1 = _oblate ?
        // For oblate
        // delta =
        //   atan2(bet1.s() * fabs(alp1.s()), bet0.c() * alp1.c())) - omg1
        // Here omg1 = -90 (because of pi/2 shift in ext. vs int. omg)
        // bet1.s() = -1, bet0.c() = 1, alp1.s() > 0, so
        // alp1 = 90 - delta
        ang::cardinal(1) - ang::radians(delta) :
        (_prolate ?
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
      if (_debug) msg = "B.a umbilic to general";
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
      if (_debug) msg = "B.b general bet1 = -90";
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
      if (_debug) msg = "B.c general omg1 = 0";
      done = false;             // A marker
    } else {
      // Case B.d, general case
      real f[4];
      alpa = ang( _kp * fabs(omg1.s()), _k * fabs(bet1.c()));
      alpb = alpa;

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
          if (fabs(f[qb]) < 2*numeric_limits<real>::epsilon()) {
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
        done = false;           // A marker
      }
    }

    int countn = 0, countb = 0;
    if (!done) {
      // Iterative search for the solution
      if (_debug)
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
      if (_debug)
        cout << "ALP1 " << real(alp1) << "\n";
      lf = TL::fline(*this, gamma(bet1, omg1, alp1));
      fic = TL::fline::fics(lf, bet1, omg1, alp1);
      // Let aa = _k2 * Math::sq(bet2.c()), bb = _kp2 * Math::sq(omg2.s());
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
      if (_hybridalt && (swapomg || omg1.s() <= fabs(omg2.s())))
        betp = (lf.transpolar() ?
                2 * (aa + lf.gammax()) > (aa + bb) :
                2 * (bb + lf.gammax()) < (aa + bb));
      d = lf.Hybrid(fic, betp ? bet2 : omg2, bet2a, omg2a, alp2, betp);
      // Don't set backside for !betp and bet2.c() == 0.  Not sure why.
      backside = (betp || bet2.c() != 0) && signbit(bet2a.c());
    }

    if (!_biaxial && backside) alp2.reflect(true, true);
    alp2.round();

    TL::gline lg(*this, lf.gm());
    TL::gline::gics gic(lg, fic);
    s12 = lg.dist(gic, d);

    if (_debug)
      cout << "FLIPS " << flip1 << flip2 << swap12 << swapomg << flipz << flipy << flipx << flipomg
           << lf.transpolar() << "\n";
    if (_debug)
      cout << "A "
           << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
           << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    // Undo switches in reverse order flipomg flipx flipy flipz swapomg swap12
    // flip2 flip1
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
      if (_prolate) {
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
    if (swapomg) {
      swap(omg1, omg2);
      if (lf.transpolar()) {
        real
          calp1 = copysign(hypot(_k/_kp*bet1.c(), lf.nu()), alp1.c()),
          calp2 = copysign(hypot(_k/_kp*bet2.c(), lf.nu()), alp2.c()),
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
          salp1 = -copysign(hypot(_kp/_k*omg1.s(), lf.nu()), alp1.s()),
          salp2 = -copysign(hypot(_kp/_k*omg2.s(), lf.nu()), alp2.s()),
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
      if (_debug)
        cout << "E " << lf.transpolar() << " " << lf.nu() << " " << lf.nup() << " "
             <<  (omg2.s() - lf.nu()) * (omg2.s() + lf.nu()) << " "
             << (lf.nup() - omg2.c()) * (lf.nup() + omg2.c()) << " "
             << real(bet1) << " " << real(omg1) << " " << real(alp1) << " "
             << real(bet2) << " " << real(omg2) << " " << real(alp2) << "\n";
    }
    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
      swap(alp1, alp2);
      swap(umb1, umb2);
      // points not swapped if umb1 == true
      alp1 += ang::cardinal(2);
      if (umb2 && !_biaxial)
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
    gic.s13 = signbit(s12) ? 0 : s12;

    if (_debug)
      cout << countn << " " << countb << " "
           << lf.gamma() << " "
           << lf.fbet().NCoeffs() << " " << lf.fomg().NCoeffs() << " "
           << lg.gbet().NCoeffs() << " " << lg.gomg().NCoeffs() << " MSG "
           << msg << "\n";
    // clang needs std::move instead of move.
    return TL(std::move(lf), std::move(fic), std::move(lg), std::move(gic));
  }

  Math::real Triaxial::HybridA(Angle bet1, Angle omg1, Angle alp1,
                               Angle bet2a, Angle omg2b,
                               bool betp) const {
    ang b1{bet1}, o1{omg1}, a1{alp1};
    // a1 -= ang(1e-8);
    gamblk gam = gamma(b1, o1, a1);
    TriaxialLine::fline l(*this, gam);
    TriaxialLine::fline::fics ic(l, b1, o1, a1);
    real dang = l.Hybrid0(ic, bet2a, omg2b, betp);
    if (_debug)
      cout << "HA " << signbit(gam.gamma) << " " << gam.gamma << " "
           << real(bet1) << " " << real(omg1) << " "
           << real(alp1) << " " << real(bet2a) << " " << real(omg2b) << " "
           << dang << "\n";
    // make -j10 > /dev/null && echo 0 30 0 -135 | ./Geod3Solve -i $SET --debug | head
    return dang;
  }

  // Solve f(alp1) = 0 where alp1 is an azimuth and f(alp1) is the difference
  // in lontitude on bet2 and the target longitude.
  Angle Triaxial::findroot(const function<real(const Angle&)>& f,
                           Angle xa,  Angle xb,
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
    int cntn = 0, cntb = 0, maxcnt = 200;
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
    //
    //  Offset for debugging output
    real x0 = -88.965874767468835986, x0r = -88.966;
    // return ang(x0);
    if (debug) {
      for (int i = -1000; i <= 1000; ++i) {
        real x = x0r+i/real(1000), fx = f(ang(x));
        cout << "PTS " << setprecision(17) << x << " " << fx << "\n";
      }
      real fx0 = f(ang(x0));
      cout << "H0 " << real(xa) << " " << fa << " "
           << real(xb) << " " << fb << " "
           << fx0 << "\n";
    }
    // If fa and fb have the same signs, assume that the root is at one of the
    // endpoints if corresponding f is small.  Otherwise, it's an error.
    if (fa * fb >= 0) {
      if (debug)
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
        throw GeographicLib::GeographicErr("Bad inputs Triaxial::findroot");
      else
        // return best endpoint
        return fabs(fa) < fabs(fb) ? xa : xb;
    }
    // tp = 1 - t
    for (real t = 1/real(2), tp = t, ab = 0, ft = 0, fc = 0;
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
      if (!(fabs(ft) >= numeric_limits<real>::epsilon())) {
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
      real
        ca = (xc-xa).radians0(),
        cb = ca+ab,
        // Scherer has a fabs(cb).  This should be fabs(ab).
        tl = numeric_limits<real>::epsilon() / fabs(ab);
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
      real
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
        t = tp = 1/real(2);
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
    real maxdiff = 0;
    if (!(t.k2() == 0  || t.kp2() == 0)) {
      // Force small gamma to zero for triaxial case
      real
        alpdiff = 2 * alp.c() * alp.s()
        * (t.k2() * Math::sq(bet.c())+t.kp2() * Math::sq(omg.s())),
        betdiff = -2 * bet.c() * bet.s() * t.k2() * Math::sq(alp.s()),
        omgdiff = -2 * omg.c() * omg.s() * t.kp2() * Math::sq(alp.c());
      maxdiff = fmax( fabs(alpdiff), fmax( fabs(betdiff), fabs(omgdiff) ) );
    }
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

} // namespace GeographicLib

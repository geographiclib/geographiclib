/**
 * \file TriaxialCartesian.cpp
 * \brief Implementation for GeographicLib::TriaxialCartesian class
 *
 * Copyright (c) Charles Karney (2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "TriaxialCartesian.hpp"

namespace GeographicLib {

  using namespace std;

  TriaxialCartesian::TriaxialCartesian(const Triaxial& t)
    : _t(t)
    , _b(_t.b())
    , _axes{_t.a(), _t.b(), _t.c()}
    , _axes2{Math::sq(_t.a()), Math::sq(_t.b()), Math::sq(_t.c())}
    , _linecc2{(_t.a() - _t.c()) * (_t.a() + _t.c()),
               (_t.b() - _t.c()) * (_t.b() + _t.c()), 0}
  {}

  template<int n>
  pair<Math::real, Math::real> TriaxialCartesian::funp<n>::operator()(real p)
    const {
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

  Math::real
  TriaxialCartesian::cartsolve(const function<pair<real, real>(real)>& f,
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

  Math::real TriaxialCartesian::cubic(vec3 r2) const {
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

  void TriaxialCartesian::carttoellip(vec3 r, Angle& bet, Angle& omg, real& H)
    const {
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
    _t.cart2toellipint(r, bet, omg, axes);
    H = axes[2] - _t.c();
  }

  void TriaxialCartesian::elliptocart(Angle bet, Angle omg, real H, vec3& r)
    const {
    vec3 ax;
    real shift = H * (2*_t.c() + H);
    for (int k = 0; k < 2; ++k)
      ax[k] = sqrt(_axes2[k] + shift);
    ax[2] = _t.c() + H;
    real tx = hypot(_t.k() * bet.c(), _t.kp()),
      tz = hypot(_t.k(), _t.kp() * omg.s());
    r = { ax[0] * omg.c() * tx,
          ax[1] * bet.c() * omg.s(),
          ax[2] * bet.s() * tz };
  }

  template<int n>
  void TriaxialCartesian::cart2togeneric(vec3 r, Angle& phi, Angle& lam)
    const {
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
  void TriaxialCartesian::generictocart2(Angle phi, Angle lam, vec3& r)
    const {
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

  void TriaxialCartesian::cart2togeod(vec3 r, Angle& phi, Angle& lam) const {
    cart2togeneric<2>(r, phi, lam);
  }
  void TriaxialCartesian::geodtocart2(Angle phi, Angle lam, vec3& r) const {
    generictocart2<2>(phi, lam, r);
  }

  void TriaxialCartesian::cart2toparam(vec3 r, Angle& phip, Angle& lamp)
    const {
    cart2togeneric<1>(r, phip, lamp);
  }
  void TriaxialCartesian::paramtocart2(Angle phip, Angle lamp, vec3& r)
    const {
    generictocart2<1>(phip, lamp, r);
  }

  void TriaxialCartesian::cart2togeocen(vec3 r, Angle& phipp, Angle& lampp)
    const {
    cart2togeneric<0>(r, phipp, lampp);
  }
  void TriaxialCartesian::geocentocart2(Angle phipp, Angle lampp, vec3& r)
    const {
    generictocart2<0>(phipp, lampp, r);
  }

  void TriaxialCartesian::cart2tocart(vec3 r2, real h, vec3& r) const {
    vec3 rn = {r2[0] / _axes2[0], r2[1] / _axes2[1], r2[2] / _axes2[2]};
    real d = h / Math::hypot3(rn[0], rn[1], rn[2]);
    r = r2;
    for (int k = 0; k < 3; ++k)
      r[k] += rn[k] * d;
  }

  void TriaxialCartesian::carttocart2(vec3 r, vec3& r2, real& h) const {
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

  void TriaxialCartesian::carttogeod(vec3 r, Angle& phi, Angle& lam, real& h)
    const {
    vec3 r2;
    carttocart2(r, r2, h);
    cart2togeod(r2, phi, lam);
  }

  void TriaxialCartesian::geodtocart(Angle phi, Angle lam, real h, vec3& r)
    const {
    vec3 r2;
    geodtocart2(phi, lam, r2);
    cart2tocart(r2, h, r);
  }

} // namespace GeographicLib

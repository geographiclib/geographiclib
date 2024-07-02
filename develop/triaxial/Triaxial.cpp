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
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

namespace GeographicLib {

  using namespace std;

  Triaxial::Triaxial()
    : Triaxial(1, 0, 1, 0)
  {}

  Triaxial::Triaxial(Math::real a, Math::real b, Math::real c)
    : a(a)
    , b(b)
    , c(c)
    , axes({a, b, c})
    , umbalt(false)
  {
    real s = (a - c) * (a + c);
    e2 = s / Math::sq(b);
    if (s == 0) {
      // The sphere is a nonuniform limit, we can pick any values in [0,1]
      // s.t. k2 + kp2 = 1.  Here we choose to treat the sphere as an
      // oblate ellipsoid.
      kp2 = 0; k2 = 1 - kp2;
    } else {
      kp2 = (a - b) * (a + b) / s;
      k2  = (b - c) * (b + c) / s;
    }
    k = sqrt(k2); kp = sqrt(kp2);
    real ksum = k2 + kp2;
    if (! (isfinite(a) && isfinite(b) && isfinite(c) &&
           a >= b && b >=c && c > 0 &&
           fabs(ksum - 1) <= numeric_limits<real>::epsilon()) )
      throw GeographicErr("Bad semiaxes for triaxial ellipsoid");
  }

  Triaxial::Triaxial(Math::real b, Math::real e2,
                     Math::real k2, Math::real kp2)
    : b(b)
    , e2(e2)
    , k2(k2)
    , kp2(kp2)
    , k(sqrt(k2))
    , kp(sqrt(kp2))
    , umbalt(false)
  {
    a = b * sqrt(1 + e2 * k2);
    c = b * sqrt(1 - e2 * kp2);
    axes = vec3({a, b, c});
    real ksum = k2 + kp2;
    if (! (isfinite(a) && isfinite(b) && isfinite(c) &&
           a >= b && b >=c && c > 0 &&
           fabs(ksum - 1) <= numeric_limits<real>::epsilon()) )
      throw GeographicErr("Bad semiaxes for triaxial ellipsoid");
  }

  void Triaxial::Norm(vec3& r) const {
    vec3 rn{r[0] / axes[0], r[1] / axes[1], r[2] / axes[2]};
    real ra = hypot3(rn[0], rn[1], rn[2]);
    r = {r[0] / ra, r[1] / ra, r[2] / ra};
  }

  void Triaxial::Norm(vec3& r, vec3& v) const {
    Norm(r);
    vec3 axes2 = {Math::sq(axes[0]), Math::sq(axes[1]), Math::sq(axes[2])},
      up = {r[0] / axes2[0], r[1] / axes2[1], r[2] / axes2[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      uv = up[0] * v[0] + up[1] * v[1] + up[2] * v[2],
      f = uv/u2;
    v = {v[0] - f * up[0], v[1] - f * up[1], v[2] - f * up[2]};
    normvec(v);
  }

  geod_fun::geod_fun(real kap, real kapp, real eps, real mu,
                     real epspow, real nmaxmult)
    : geod_fun(kap, kapp, eps, mu,
               (mu > 0 ? mu / (kap + mu) :
                (mu < 0 ? -mu / kap :
                 kapp)) < Triaxial::EllipticThresh(), epspow, nmaxmult)
  {}

  geod_fun::geod_fun(real kap, real kapp, real eps, real mu, bool tx,
                     real epspow, real nmaxmult)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _tx(tx)
    , _tol(pow(numeric_limits<real>::epsilon(), epspow))
    , _invp(false)
  {
    real k2 = 0, kp2 = 1;
    if (_tx) {
      k2 = _mu > 0 ? _kap / (_kap + _mu) :
        (_mu < 0 ? (_kap + _mu) / _kap :
         _kap);                    // _mu == 0
      kp2 = _mu > 0 ? _mu / (_kap + _mu) :
        (_mu < 0 ? -_mu / _kap :
         _kapp);                // _mu == 0
    }
    _ell = EllipticFunction(k2, 0, kp2, 1);
    if (_mu > 0) {
      _fun = _tx ?
        TrigfunExt(
                   [kap = _kap, kapp = _kapp,
                    eps = _eps, mu = _mu, ell = _ell]
                   (real u) -> real
                   { real sn, cn, dn; (void) ell.am(u, sn, cn, dn);
                     return fup(cn, kap, kapp, eps, mu); },
                   _ell.K(), false, epspow, nmaxmult) :
        TrigfunExt(
                   [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                   (real phi) -> real
                   { return fphip(cos(phi), kap, kapp, eps, mu); },
                   Math::pi()/2, false, epspow, nmaxmult);
    } else if (_mu < 0) {
      _fun = _tx ?
        TrigfunExt(
                   [kap = _kap, kapp = _kapp,
                    eps = _eps, mu = _mu, ell = _ell]
                   (real v) -> real
                   { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                     return fvp(dn, kap, kapp, eps, mu); },
                   _ell.K(), false, epspow, nmaxmult) :
        TrigfunExt(
                   [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                   (real psi) -> real
                   { return fpsip(sin(psi), cos(psi), kap, kapp, eps, mu); },
                   Math::pi()/2);
    } else {                    // _mu == 0
      // N.B. Don't compute the inverse of _fun so not really necessary to
      // supply epspow and nmaxmult args to TrigfunExt.
      _fun = _tx ?
        TrigfunExt(
                   [kap = _kap, kapp = _kapp,
                    eps = _eps, ell = _ell]
                   (real v) -> real
                   { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                    return dfvp(cn, dn, kap, kapp, eps); },
                   2 * _ell.K(), true, epspow, nmaxmult) :
        TrigfunExt(
                   [kap = _kap, kapp = _kapp, eps = _eps]
                   (real phi) -> real
                   { return dfp(cos(phi), kap, kapp, eps); },
                   Math::pi(), true, epspow, nmaxmult);
    }
    _nmax = nmaxmult ? int(ceil(nmaxmult * _fun.NCoeffs())) : 1 << 16;
    _max = _mu == 0 ? _fun(_tx ? _ell.K() : Math::pi()/2) :
        _fun.Max();
  }

  void geod_fun::ComputeInverse() {
    if (!_invp) {
      if (_mu == 0) {
        _chiinv = _tx ?
          // chi maps [-1,1]*K -> [-1,1]*pi/2
          // dchi/dx = sech(lam(am(x)) - dfv(x)) *
          //           (sec(am(x)) * dn(x) - dfvp(x))
          // cosh(lam(am(x))) = 1/cn, sinh(lam(am(x))) = sn/cn
          // cosh(lam(x) - dfv(x))
          //  = cosh(lam(am(x))) * cosh(dfv(x)) -
          //    sinh(lam(am(x))) * sinh(dfv(x))
          //  = 1/cn * cosh(dfv(x)) - sn/cn * sinh(dfv(x))
          //  = 1/cn * (cosh(dfv(x)) - sn * sinh(dfv(x)))
          // (sec(am(x)) * dn - dfvp(x)) = 1/cn * (dn - cn * dfvp(x))
          // dchi/dx = (dn - cn * dfvp(x)) /
          //           (cosh(dfv(x)) - sn * sinh(dfv(x)))
          Trigfun::InverseInit(
                               [fun = _fun, ell = _ell]
                               (real x) -> pair<real, real>
                               { real sn, cn, dn, f = fun(x);
                                 (void) ell.am(x, sn, cn, dn);
                                 return pair<real,real>
                                   ( gd(lam(ell.am(x)) - f),
                                     (dn - cn * fun.deriv(x)) /
                                     (cosh(f) - sn * sinh(f)) ); },
                               _ell.K(), Math::pi()/2,
                               -Math::pi()/2, Math::pi()/2, &_countn, &_countb,
                               _tol, _nmax) :
          // chi maps [-1,1]*pi/2 -> [-1,1]*pi/2
          // dchi/dx = sech(lam(x) - df(x)) * (sec(x) - dfp(x))
          // cosh(lam(x)) = sec(x), sinh(lam(x)) = tan(x)
          // cosh(lam(x) - df(x))
          //  = cosh(lam(x)) * cosh(df(x)) - sinh(lam(x)) * sinh(df(x))
          //  = sec(x) * cosh(df(x)) - tan(x) * sinh(df(x))
          //  = sec(x) * (cosh(df(x)) - sin(x) * sinh(df(x)))
          // (sec(x) - dfp(x)) = sec(x) * (1 - cos(x) * dfp(x))
          // dchi/dx = (1 - cos(x) * dfp(x)) /
          //           (cosh(df(x)) - sin(x) * sinh(df(x)))
          Trigfun::InverseInit(
                               [fun = _fun]
                               (real x) -> pair<real, real>
                               { real f = fun(x);
                                 return pair<real, real>
                                   ( gd(lam(x) - f),
                                     (1 - cos(x) * fun.deriv(x)) /
                                     (cosh(f) - sin(x) * sinh(f)) ); },
                               Math::pi()/2, Math::pi()/2,
                               -Math::pi()/2, Math::pi()/2, &_countn, &_countb,
                               _tol, _nmax);
      }
      else
        _fun.ComputeInverse();
    }
    _invp = true;
  }

  Math::real geod_fun::root(real z, real x0, int* countn, int* countb) const {
    if (_mu != 0) return Math::NaN();
    if (!isfinite(z)) return z; // Deals with +/-inf and nan
    real d = Max()
      + 2 * numeric_limits<real>::epsilon() * fmax(real(1), fabs(z)),
      xa = z - d,
      xb = z + d;
    x0 = fmin(xb, fmax(xa, x0));
    // Solve z = u - _fun(_tx ? _ell.F(gd(u)) : gd(u)) for u
    // N.B. use default tol for root, because we want accurate answers here

    return _tx ?
      Trigfun::root(
                    [fun = _fun, ell = _ell]
                    (real u) -> pair<real, real>
                    { real phi = gd(u), sch = 1/cosh(u);
                      return pair<real, real>
                        (u - fun(ell.F(phi)),
                         1 - fun.deriv(ell.F(phi)) * sch /
                         sqrt(ell.kp2() + ell.k2() * Math::sq(sch))); },
                    z,
                    x0, xa, xb,
                    Math::pi()/2, Math::pi()/2, 1, countn, countb) :
      Trigfun::root(
                    [fun = _fun]
                    (real u) -> pair<real, real>
                    { real phi = gd(u), sch = 1/cosh(u);
                      return pair<real, real>
                        (u - fun(phi),
                         1 - fun.deriv(phi) * sch); },
                    z,
                    x0, xa, xb,
                    Math::pi()/2, Math::pi()/2, 1, countn, countb);
  }
  // Approximate inverse using _chiinv or _fun.inv0
  Math::real geod_fun::inv0(real z) const {
    if (!_invp) return Math::NaN();
    return _mu == 0 ? _fun(_chiinv(gd(z))) : _fun.inv0(z);
  }
  // Accurate inverse by direct Newton (not using _finv)
  Math::real geod_fun::inv1(real z, int* countn, int* countb) const {
    return _mu == 0 ? root(z, z, countn, countb) :
      _fun.inv1(z, countn, countb);
  }
  // Accurate inverse correcting result from _finv
  Math::real geod_fun::inv(real z, int* countn, int* countb) const {
    if (!_invp) return Math::NaN();
    return _mu == 0 ? root(z, inv0(z), countn, countb) :
      _fun.inv(z, countn, countb);
  }

  dist_fun::dist_fun(real kap, real kapp, real eps, real mu)
    : dist_fun(kap, kapp, eps, mu,
               (mu > 0 ? mu / (kap + mu) :
                (mu < 0 ? -mu / kap :
                 kapp)) < Triaxial::EllipticThresh())
  {}

  dist_fun::dist_fun(real kap, real kapp, real eps, real mu, bool tx)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _tx(tx)
  {
    real k2 = 0, kp2 = 1;
    if (_tx) {
      k2 = _mu > 0 ? _kap / (_kap + _mu) :
        (_mu < 0 ? (_kap + _mu) / _kap :
         _kap);                    // _mu == 0
      kp2 = _mu > 0 ? _mu / (_kap + _mu) :
        (_mu < 0 ? -_mu / _kap :
         _kapp);                // _mu == 0
    }
    _ell = EllipticFunction(k2, 0, kp2, 1);
    if (_mu > 0) {
      _fun = _tx ?
        TrigfunExt(
                   [kap = _kap, kapp = _kapp,
                    eps = _eps, mu = _mu, ell = _ell]
                   (real u) -> real
                   { real sn, cn, dn; (void) ell.am(u, sn, cn, dn);
                     return gup(cn, dn, kap, kapp, eps, mu); },
                   _ell.K()) :
        TrigfunExt(
                   [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                   (real phi) -> real
                   { return gphip(cos(phi), kap, kapp, eps, mu); },
                   Math::pi()/2);
    } else if (_mu < 0) {
      _fun = _tx ?
        TrigfunExt(
                   [kap = _kap, kapp = _kapp,
                    eps = _eps, mu = _mu, ell = _ell]
                   (real v) -> real
                   { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                     return gvp(cn, dn, kap, kapp, eps, mu); },
                   _ell.K()) :
        TrigfunExt(
                   [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                   (real psi) -> real
                   { return gpsip(sin(psi), cos(psi), kap, kapp, eps, mu); },
                   Math::pi()/2);
    } else {                    // _mu == 0
      _fun = _tx ?
        TrigfunExt(
                   [kap = _kap, kapp = _kapp, eps = _eps, ell = _ell]
                   (real v) -> real
                   { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                     return g0vp(cn, kap, kapp, eps); },
                   2*_ell.K(), true) :
        TrigfunExt(
                   [kap = _kap, kapp = _kapp, eps = _eps, ell = _ell]
                   (real phi) -> real
                   { return g0p(cos(phi), kap, kapp, eps); },
                   Math::pi(), true);
    }
    _max = _mu == 0 ? _fun(_tx ? _ell.K() : Math::pi()/2) :
      _fun.Max();
  }

  Math::real dist_fun::gfderiv(real u) const {
    real sn = 0, cn = 0, dn = 0;
    if (_mu != 0 && _tx)
      (void) _ell.am(u, sn, cn, dn);
    if (_mu > 0)
      return _tx ? gfup(cn, _kap, _mu) : gfphip(cos(u), _kap, _mu);
    else if (_mu < 0)
      return _tx ? gfvp(cn, _kap, _mu) : gfpsip(cos(u), _kap, _mu);
    else                      // _mu == 0
      return gf0up(u, _kap, _kapp);
  }

  void Triaxial::cart2toellip(const vec3& r, Angle& bet, Angle& omg)
    const {
    real xi = r[0]/a, eta = r[1]/b, zeta = r[2]/c,
      g = k2 * Math::sq(xi)
      + (k2 - kp2) * Math::sq(eta)
      - kp2 * Math::sq(zeta);
    if (fabs(r[0]) == a * kp2 && r[1] == 0 && fabs(r[2]) == c * k2)
      g = 0;
    real h = hypot(g, 2 * k * kp * eta),
      so, co, sb, cb;
    if (h == 0) {
      so = 0;
      cb = 0;
    } else if (g < 0) {
      so = copysign(sqrt( (h - g)/2 ) / kp, eta);
      cb = fabs(eta / so);
    } else {
      cb = sqrt( (h + g)/2 ) / k;
      so = eta / cb;
    }
    real tz = hypot(k, kp * so),
      tx = hypot(k * cb, kp);
    sb = tz == 0 ? -1 : zeta / tz;
    co = tx == 0 ?  1 : xi / tx;
    bet = ang(sb, cb, 0, true); omg = ang(so, co, 0, true);
  }

  void Triaxial:: cart2toellip(const Angle& bet, const Angle& omg,
                               const vec3& v, Angle& alp) const {
    real tz = hypot(k, kp * omg.s()),
      tx = hypot(k * bet.c(), kp);
    if (!(bet.c() == 0 && omg.s() == 0)) {
      vec3
        N = tx == 0 ?
        vec3{-omg.c() * bet.s(), -omg.s() * bet.s(), tx * bet.s()} :
        (tz == 0 ?
         vec3{tz, -bet.s(), bet.c()} :
         vec3{-a * k2 * bet.c() * bet.s() * omg.c() / tx,
            -b * bet.s() * omg.s(),
            c * bet.c() * tz}),
        E = tx == 0 ?
        vec3{-omg.s(),  omg.c(), tx} :
        (tz == 0 ?
         vec3{tz * omg.c(),  bet.c() * omg.c(), bet.s() * omg.c()} :
         vec3{-a * tx * omg.s(),
            b * bet.c() * omg.c(),
            c * kp2 * bet.s() * omg.c() * omg.s() / tz});
      normvec(N); normvec(E);
      alp = ang(v[0] * E[0] + v[1] * E[1] + v[2] * E[2],
                v[0] * N[0] + v[1] * N[1] + v[2] * N[2], 0, true);
    } else {                    // bet.c() == 0 && omg.s() == 0
      // Special treatment at umbilical points
      real w = bet.s() * omg.c(),
        upx = omg.c() * tx/a, upz = bet.s() * tz/c;
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

  void Triaxial::cart2toellip(const vec3& r, const vec3& v,
                              Angle& bet, Angle& omg, Angle& alp)
    const {
    cart2toellip(r, bet, omg);
    cart2toellip(bet, omg, v, alp);
  }

  void Triaxial::elliptocart2(const Angle& bet, const Angle& omg,
                              vec3& r) const {
    real tx = hypot(k * bet.c(), kp), tz = hypot(k, kp * omg.s());
    r = vec3{ a * omg.c() * tx,
              b * bet.c() * omg.s(),
              c * bet.s() * tz };
    // Norm(r); r is already normalized
  }

  void Triaxial::elliptocart2(const Angle& bet, const Angle& omg,
                              const Angle& alp,
                              vec3& r, vec3& v) const {
    elliptocart2(bet, omg, r);
    real tx = hypot(k * bet.c(), kp), tz = hypot(k, kp * omg.s());
    if (bet.c() == 0 && omg.s() == 0 && !(k == 0 || kp == 0)) {
      // umbilical point (not oblate or prolate)
      real sa2 = 2 * alp.s() * alp.c(),
        ca2 = (alp.c() - alp.s()) * (alp.c() + alp.s());
      // sign on 2nd component is -sign(cos(bet)*sin(omg)).  negative sign
      // gives normal convention of alpha measured clockwise.
      v = vec3{a*k/b * omg.c() * ca2,
               -omg.c() * bet.s() * sa2,
               -c*kp/b * bet.s() * ca2};
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
        N = vec3{ -a * k2 * bet.c() * bet.s() * omg.c() / tx,
                  -b * bet.s() * omg.s(),
                  c * bet.c() * tz};
        E = vec3{ -a * tx * omg.s(),
                  b * bet.c() * omg.c(),
                  c * kp2 * bet.s() * omg.c() * omg.s() / tz};
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
                                 Angle bet2, Angle omg2)
    const {
    bool debug = false;
    bet1.rnd();
    omg1.rnd();
    bet2.rnd();
    omg2.rnd();
    bool flip1 = AngNorm(bet1, omg1);
    (void) AngNorm(bet2, omg2);
    bool umb1 = bet1.c() == 0 && omg1.s() == 0,
      umb2 = bet2.c() == 0 && omg2.s() == 0;
    bool swap12;
    /*
    if (umb1 ^ umb2)
      // NOT SURE ABOUT THIS
      swap12 = umb1;          // If one umbilic point, make is point 2
      else */
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
    ang alp1, alp2, bet2a, omg2a;

    TriaxialLineF lf;
    TriaxialLineF::fics fic;
    TriaxialLineF::disttx d;

    // flag for progress
    bool done = false;
    if (bet1.c() * omg1.s() == 0 && bet2.c() * omg2.s() == 0) {
      // both points on middle ellipse
      lf = TriaxialLineF(*this, Triaxial::gamblk{}, 0.5, 1.5);
      if (umb1 && umb2 && bet2.s() > 0 && omg2.c() < 0) {
        // process opposite umbilical points
        fic = TriaxialLineF::fics(lf, bet1, omg1, ang{kp, k});
        alp1 = ang(kp * exp(lf.df), k);
        fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
        d = lf.ArcPos0(fic, Math::pi(), bet2a, omg2a, alp2, true);
        if (debug) cout << "opposite umbilics\n";
        done = true;
      } else if (bet1.c() == 0 && bet2.c() == 0) {
        // bet1 = -90, bet2 = +/-90
        if (bet2.s() < 1) {
          // bet1 = bet2 = -90
          alp1 = ang::cardinal(1);
          fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
          ang omg12 = omg2 - omg1;
          d = omg12.s() == 0 && omg12.c() < 0 ?
            // adjacent E/W umbilical points
            TriaxialLineF::disttx{-Triaxial::BigValue(),
                                  Triaxial::BigValue(), 0} :
            lf.ArcPos0(fic, omg12.radians0(), bet2a, omg2a, alp2, false);
          if (omg2a.s() < 0) alp2.reflect(true); // Is this needed?
          if (debug) {
            if (omg12.s() == 0 && omg12.c() < 0)
              cout << "adj EW umbilics\n";
            else
              cout << "bet1 = bet2 = -90\n";
          }
          done = true;
        } else {
          // bet1 = -90, bet2 = 90
          // need to see how far apart the points are
          // If point 1 is at [-90, 0], direction is 0 else -90.
          alp1 = ang::cardinal(omg1.s() == 0 ? 0 : -1);
          fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
          // If point 1 is [-90, 0] and point 2 is [90, 0]
          if (omg1.s() == 0 && omg2.s() == 0) {
            // adjacent N/S umbilical points
            d = TriaxialLineF::disttx{Triaxial::BigValue(),
                                      -Triaxial::BigValue(), 0};
            if (debug) cout << "bet1 = -90,  bet2 = 90, adj umb\n";
            done = true;
          } else {
            d = lf.ArcPos0(fic, Math::pi(), bet2a, omg2a, alp2);
            omg2a -= omg2;
            if (omg2a.s() > 0) {
              omg2a = omg2 + omg1;
              d = lf.ArcPos0(fic, omg2a.radians0(), bet2a, omg2a, alp2, false);
              if (omg2a.s() < 0) alp2.reflect(true); // Is this needed?
              if (debug) cout << "bet1 = -90,  bet2 = 90, merid\n";
              done = true;
            } else {
              alpa = ang::cardinal(-1) + ang::eps();
              fa = omg2a.radians0();
              alpb.setquadrant(0U);
              fic.setquadrant(lf, 0U);
              (void) lf.ArcPos0(fic, Math::pi(), bet2a, omg2a, alp2);
              omg2a -= omg2;
              alpb = -alpa;
              fb = omg2a.radians0();
              if (debug) cout << "bet1 = -90,  bet2 = 90, non-merid\n";
            }
          }
        }
      } else {
        // other meridional cases, invoke Hybrid with the following value of
        // alp1
        alp1 = ang::cardinal(bet1.c() == 0 ?
                             (omg2.c() < 1 ? 1 : (omg1.s() == 0 ? 0 : -1)) :
                             (omg2.c() > 0 ? 0 : 2));
        fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
        d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
        if (debug) cout << "other merid\n";
        done = true;
      }
    } else if (bet1.s() == 0 && bet2.s() == 0) {
      // both points on equator
      ang omg12 = omg2 - omg1;
      int eE = omg12.s() > 0 ? 1 : -1;
      // set direction for probe as +/-90 based on sign of omg12
      alp1 = ang::cardinal(eE);
      bet1.reflect(true);
      lf = TriaxialLineF(*this, gamblk(*this, bet1, omg1, alp1), 0.5, 1.5);
      fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
      (void) lf.ArcPos0(fic, Math::pi(), bet2a, omg2a, alp2);
      omg2a -= omg2;
      if (eE * omg2a.s() >= 0) {
        // geodesic follows the equator
        d = lf.ArcPos0(fic, eE * omg12.radians0(), bet2a, omg2a, alp2, false);
        if (debug) cout << "bet1 = bet2 = 0 equatorial\n";
        done = true;
      } else {
        // geodesic does not follow the equator
        alpb = ang::cardinal(-1) - ang::eps();
        alpa = -alpb;
        (eE > 0 ? fa : fb) = omg2a.radians0();
        alp1.setquadrant(eE > 0 ? 3U : 0U);
        fic.setquadrant(lf, eE > 0 ? 3U : 0U);
        (void) lf.ArcPos0(fic, Math::pi(), bet2a, omg2a, alp2);
        omg2a -= omg2;
        (eE > 0 ? fb : fa) = omg2a.radians0();
        if (debug) cout << "bet1 = bet2 = 0 non-equatorial\n";
      }
    } else if (umb1) {
      // umbilical point to general point
      lf = TriaxialLineF(*this, Triaxial::gamblk{}, 0.5, 1.5);
      alp2 = ang(kp * omg2.s(), k * bet2.c());
      fic = TriaxialLineF::fics(lf, bet2, omg2, alp2);
      (void) lf.ArcPos0(fic, (bet1 - bet2).radians0(), bet2a, omg2a, alp1);
      if (alp1.s() < 0) alp1 += ang::cardinal(1);
      fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
      d = lf.ArcPos0(fic, (bet2 - bet1).radians0(), bet2a, omg2a, alp2);
      if (debug) cout << "umb to general\n";
      done = true;
    } else if (bet1.c() == 0) {
      // bet1 = -90 to general point
      if (omg2.s() > 0) {
        alpa = ang::cardinal(-1) + ang::eps();
        alpb = -alpa;
        fa = -omg2.radians0();
        fb = (ang::cardinal(2) - omg2).radians0();
      } else {
        alpa = ang::cardinal(1) + ang::eps();
        alpb = -alpa;
        fa = (ang::cardinal(2) - omg2).radians0();
        fb = -omg2.radians0();
      }
      if (debug) cout << "bet = -90 to general\n";
    } else if (omg1.s() == 0) {
      // omg1 = 0 to general point
      if (omg2.s() > 0) {
        alpa = ang::eps();
        alpb = ang::cardinal(2) - alpa;
        fa = -omg2.radians0();
        fb = (ang::cardinal(2)-omg2).radians0();
      } else {
        alpa = ang(-numeric_limits<real>::epsilon()/(1<<20), -1);
        alpb = ang(-numeric_limits<real>::epsilon()/(1<<20),  1);
        fa = (ang::cardinal(2)-omg2).radians0();
        fb = -omg2.radians0();
      }
      if (debug) cout << "omg1 = 0 to general\n";
    } else {
      // general case
      real f[4];
      alpa = ang( kp * fabs(omg1.s()), k * fabs(bet1.c()));
      alpb = alpa;

      lf = TriaxialLineF(*this, Triaxial::gamblk{}, 0.5, 1.5);
      fic = TriaxialLineF::fics(lf, bet1, omg1, alpb);
      {
        unsigned qb = 0U, qa = 3U; // qa = qb - 1 (mod 4)
        for (; qb <= 4U; ++qb, ++qa) {
          if (qb) {
            alpb.setquadrant(qb);
            fic.setquadrant(lf, qb);
          }
          if (qb < 4U) {
            f[qb] = lf.Hybrid0(fic, bet2, omg2);
            if (fabs(f[qb]) < numeric_limits<real>::epsilon()) {
              alp1 = alpb;
              d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
              done = true;
              break;
            }
          }
          if (qb && (f[qa & 3U] < 0 && f[qb & 3U] > 0)) {
            break;
          }
        }
        if (qb > 4U) std::cout << "ERROR\n";
        fa = f[qa & 3U]; fb = f[qb & 3U];
        alpa.setquadrant(qa & 3U);
      }
      if (debug) cout << "general\n";
    }

    if (!done) {
      int countn = 0, countb = 0;
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
      lf = TriaxialLineF(*this, gamblk(*this, bet1, omg1, alp1), 0.5, 1.5);
      fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
      d = lf.Hybrid(fic, bet2, bet2a, omg2a, alp2);
    }

    TriaxialLineG lg(*this, lf.gm());
    TriaxialLineG::gics gic(lg, fic);
    real s13 = lg.dist(gic, d);
    (void) AngNorm(bet2a, omg2a, alp2);

    // Undo switches in reverse order flipz, swap12, flip1
    if (flipomg) {
      omg2.reflect(true);
      alp2.reflect(true, true);
    }

    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
      alp1.reflect(true);
      alp2.reflect(true);
    }

    if (flipy) {
      omg1.reflect(true);
      omg2.reflect(true, true);
      alp1.reflect(true);
      alp2.reflect(true);
    }

    if (flipz) {
      bet1.reflect(true);
      bet2.reflect(true);
      alp1.reflect(false, true);
      alp2.reflect(false, true);
    }

    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
      swap(alp1, alp2);
      swap(umb1, umb2);
      alp1 += ang::cardinal(umb1 ? 1 : 2);
      alp2 += ang::cardinal(umb2 ? 1 : 2);
    }

    if (flip1)
      Flip(bet1, omg1, alp1);

    if (flip1 || swap12 || flipz || flipy || flipx) {
      fic = TriaxialLineF::fics(lf, bet1, omg1, alp1);
      gic = TriaxialLineG::gics(lg, fic);
    }
    gic.s13 = fmax(real(0), s13);

    return TriaxialLine(move(lf), move(fic), move(lg), move(gic));
  }

  Math::real Triaxial::HybridA(const Triaxial& t,
                               const Angle& bet1, const Angle& omg1,
                               const Angle& alp1,
                               const Angle& bet2, const Angle& omg2) {
    ang b1{bet1}, o1{omg1}, a1{alp1};
    gamblk gam(t, b1, o1, a1);
    TriaxialLineF l(t, gam, 0.5, 1.5);
    TriaxialLineF::fics ic(l, b1, o1, a1);
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
           ||GEOGRAPHICLIB_PANIC;) {
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

} // namespace GeographicLib

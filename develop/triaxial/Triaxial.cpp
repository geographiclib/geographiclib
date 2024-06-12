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
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

namespace GeographicLib {

  using namespace std;

  Triaxial::Triaxial()
    : Triaxial(Math::NaN(), Math::NaN(), Math::NaN())
  {}

  Triaxial::Triaxial(Math::real a, Math::real b, Math::real c)
    : a(a)
    , b(b)
    , c(c)
    , axes({a, b, c})
    , axesn({a/b, real(1), c/b})
    , axes2n({Math::sq(a/b), real(1), Math::sq(c/b)})
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

  void Triaxial::Norm(vec6& y) const {
    real ra = hypot3(y[0] / axesn[0], y[1], y[2] / axesn[2]);
    y[0] /= ra; y[1] /= ra; y[2] /= ra;
    vec3 up = {y[0] / axes2n[0], y[1], y[2] / axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      uv = up[0] * y[3+0] + up[1] * y[3+1] + up[2] * y[3+2],
      f = uv/u2;
    y[3+0] -= f * up[0]; y[3+1] -= f * up[1]; y[3+2] -= f * up[2];
    f = hypot3(y[3+0], y[3+1], y[3+2]);
    y[3+0] /= f; y[3+1] /= f; y[3+2] /= f;
  }

  Triaxial::vec6 Triaxial::Accel(const vec6& y) const {
    vec3 up = {y[0] / axes2n[0], y[1], y[2] / axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      f = -(Math::sq(y[3+0]) / axes2n[0] +
            Math::sq(y[3+1]) +
            Math::sq(y[3+2]) / axes2n[2]) / u2;
    return vec6{y[3+0], y[3+1], y[3+2], f * up[0], f * up[1], f * up[2]};
  }

  int Triaxial::Direct(const vec3& r1, const vec3& v1, real s12,
                       vec3& r2, vec3& v2, real eps) const {
    using namespace boost::numeric::odeint;
    if (eps == 0)
      eps = pow(numeric_limits<real>::epsilon(), real(7)/8);
    // Normalize all distances to b.
    vec6 y{r1[0]/b, r1[1]/b, r1[2]/b, v1[0], v1[1], v1[2]};
    Norm(y);
    int n;
    s12 /= b;
    auto fun = [this](const vec6& y, vec6& yp, real /*t*/) -> void {
      yp = y; Norm(yp); yp = Accel(yp);
    };
    if (1) {
      bulirsch_stoer<vec6, real> stepper(eps, real(0));
      n = integrate_adaptive(stepper, fun, y, real(0), s12, 1/real(10));
    } else {
      n = int(round(1/eps));
      runge_kutta4<vec6, real> stepper;
      integrate_n_steps(stepper, fun, y, real(0), s12/n, n);
    }
    Norm(y);
    r2 = {b*y[0], b*y[1], b*y[2]};
    v2 = {y[3], y[4], y[5]};
    return n;
  }

  void Triaxial::Direct(const vec3& r1, const vec3& v1, real ds,
                        long nmin, long nmax,
                        vector<vec3>& r2, vector<vec3>& v2, real eps) const {
    using namespace boost::numeric::odeint;
    if (nmin > 0 || nmax < 0) return;
    if (eps == 0)
      eps = pow(numeric_limits<real>::epsilon(), real(7)/8);
    r2.clear(); v2.clear();
    auto fun = [this](const vec6& y, vec6& yp, real /*t*/) -> void {
      yp = y; Norm(yp); yp = Accel(yp);
    };
    auto observer = [&r2, &v2, b = b](const vec6& y, real /*t*/) -> void {
      r2.push_back( {b*y[0], b*y[1], b*y[2]} );
      v2.push_back( {y[3], y[4], y[5]} );
    };
    bulirsch_stoer<vec6, real> stepper(eps, real(0));
    ds /= b;
    if (nmin < 0) {
      vec6 y{r1[0]/b, r1[1]/b, r1[2]/b, v1[0], v1[1], v1[2]}; Norm(y);
      integrate_n_steps(stepper, fun, y, real(0), -ds, abs(nmin), observer);
      reverse(r2.begin(), r2.end());
      reverse(v2.begin(), v2.end());
    }
    if (nmax > 0) {
      if (nmin < 0) { r2.pop_back(); v2.pop_back(); }
      vec6 y{r1[0]/b, r1[1]/b, r1[2]/b, v1[0], v1[1], v1[2]}; Norm(y);
      integrate_n_steps(stepper, fun, y, real(0), ds, nmax, observer);
    }
    for (int i = 0; i <= nmax - nmin; ++i)
      Norm(r2[i], v2[i]);
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
    real xa = z - Max(), xb = z + Max();
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
    return _mu == 0 ? root(z, z, countn, countb) : _fun.inv1(z, countn, countb);
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

  void Triaxial::cart2toellip(const vec3& r, AuxAngle& bet, AuxAngle& omg)
    const {
    real xi = r[0]/a, eta = r[1]/b, zeta = r[2]/c,
      g = k2 * Math::sq(xi)
      + (k2 - kp2) * Math::sq(eta)
      - kp2 * Math::sq(zeta);
    if (fabs(r[0]) == a * kp2 && r[1] == 0 && fabs(r[2]) == c * k2)
      g = 0;
    real h = hypot(g, 2 * k * kp * eta);
    if (h == 0) {
      omg.y() = 0;
      bet.x() = 0;
    } else if (g < 0) {
      omg.y() = copysign(sqrt( (h - g)/2 ) / kp, eta);
      bet.x() = fabs(eta / omg.y());
    } else {
      bet.x() = sqrt( (h + g)/2 ) / k;
      omg.y() = eta / bet.x();
    }
    real tz = hypot(k, kp * omg.y()),
      tx = hypot(k * bet.x(), kp);
    bet.y() = tz == 0 ? -1 : zeta / tz;
    omg.x() = tx == 0 ?  1 : xi / tx;
    bet.normalize(); omg.normalize();
  }

  void Triaxial::cart2toellip(const vec3& r, const vec3& v,
                              AuxAngle& bet, AuxAngle& omg, AuxAngle& alp)
    const {
    cart2toellip(r, bet, omg);
    real tz = hypot(k, kp * omg.y()),
      tx = hypot(k * bet.x(), kp);
    if (!(bet.x() == 0 && omg.y() == 0)) {
      vec3
        N = tx == 0 ?
        vec3{-omg.x() * bet.y(), -omg.y() * bet.y(), tx * bet.y()} :
        (tz == 0 ?
         vec3{tz, -bet.y(), bet.x()} :
         vec3{-a * k2 * bet.x() * bet.y() * omg.x() / tx,
            -b * bet.y() * omg.y(),
            c * bet.x() * tz}),
        E = tx == 0 ?
        vec3{-omg.y(),  omg.x(), tx} :
        (tz == 0 ?
         vec3{tz * omg.x(),  bet.x() * omg.x(), bet.y() * omg.x()} :
         vec3{-a * tx * omg.y(),
            b * bet.x() * omg.x(),
            c * kp2 * bet.y() * omg.x() * omg.y() / tz});
      normvec(N); normvec(E);
      alp.x() = v[0] * N[0] + v[1] * N[1] + v[2] * N[2];
      alp.y() = v[0] * E[0] + v[1] * E[1] + v[2] * E[2];
    } else {                    // bet.x() == 0 && omg.y() == 0
      // Special treatment at umbilical points
      real w = bet.y() * omg.x(),
        upx = r[0]/Math::sq(a), upz = r[2]/Math::sq(c);
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
      // multiply [sa, ca] by -bet.y().
      real flip = -bet.y();
      if (c2a >= 0) {
        alp.y() = flip * s2a;
        alp.x() = flip * (1 + c2a);
      } else {
        alp.y() = flip * copysign(1 + c2a, s2a);
        alp.x() = flip * fabs(s2a);
      }
    }
    alp.normalize();
  }

  void Triaxial::elliptocart2(const AuxAngle& bet, const AuxAngle& omg,
                              vec3& r) const {
    AuxAngle betn = bet.normalized(), omgn = omg.normalized();
    real tx = hypot(k * betn.x(), kp), tz = hypot(k, kp * omgn.y());
    r = vec3{ a * omgn.x() * tx,
              b * betn.x() * omgn.y(),
              c * betn.y() * tz };
    // Norm(r); r is already normalized
  }

  void Triaxial::elliptocart2(const AuxAngle& bet, const AuxAngle& omg,
                              const AuxAngle& alp,
                              vec3& r, vec3& v) const {
    elliptocart2(bet, omg, r);
    AuxAngle betn = bet.normalized(), omgn = omg.normalized();
    real tx = hypot(k * betn.x(), kp), tz = hypot(k, kp * omgn.y());
    if (betn.x() == 0 && omgn.y() == 0 && !(k == 0 || kp == 0)) {
      // umbilical point (not oblate or prolate)
      real sa2 = 2 * alp.y() * alp.x(),
        ca2 = (alp.x() - alp.y()) * (alp.x() + alp.y());
      // sign on 2nd component is -sign(cos(bet)*sin(omg)).  negative sign
      // gives normal convention of alpha measured clockwise.
      v = vec3{a*k/b * omgn.x() * ca2,
               -omgn.x() * betn.y() * sa2,
               -c*kp/b * betn.y() * ca2};
    } else {
      vec3 N, E;
      if (tx == 0) {
        // At an oblate pole tx -> cos(bet)
        N = vec3{-omgn.x() * betn.y(), -omgn.y() * betn.y(), 0};
        E = vec3{-omgn.y()           ,  omgn.x()           , 0};
      } else if (tz == 0) {
        // At a prolate pole tz -> sin(omg)
        N = vec3{0, -betn.y()           , betn.x()           };
        E = vec3{0,  betn.x() * omgn.x(), betn.y() * omgn.x()};
      } else {
        // The general case
        N = vec3{ -a * k2 * betn.x() * betn.y() * omgn.x() / tx,
                  -b * betn.y() * omgn.y(),
                  c * betn.x() * tz};
        E = vec3{ -a * tx * omgn.y(),
                  b * betn.x() * omgn.x(),
                  c * kp2 * betn.y() * omgn.x() * omgn.y() / tz};
      }
      normvec(N);
      normvec(E);
      v = vec3{alp.x() * N[0] + alp.y() * E[0],
               alp.x() * N[1] + alp.y() * E[1],
               alp.x() * N[2] + alp.y() * E[2]};
    }
    // normvec(v); v is already normalized
  }

  TriaxialLine Triaxial::Inverse(AuxAngle bet1, AuxAngle omg1,
                                 AuxAngle bet2, AuxAngle omg2)
    const {
    bet1.normalize().round();
    omg1.normalize().round();
    bet2.normalize().round();
    omg2.normalize().round();
    bool flip1 = AngNorm(bet1, omg1);
    (void) AngNorm(bet2, omg2);
    bool swap12;
    {
      AuxAngle bet1m(bet1), bet2m(bet2);
      bet1m.y() = fabs(bet1m.y()); bet2m.y() = fabs(bet2m.y());
      AuxAngle bet12 = bet2m - bet1m; // |bet2| - |bet1|
      swap12 = bet12.y() > 0;     // is |bet2| > |bet1|
      // N.B. No swapping if bet1 = +0 and bet1 = -0.
    }
    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
    }
    // Now |bet1| > |bet2|
    bool swapsign = !signbit(bet1.y());
    if (swapsign) {
      bet1.y() *= -1;
      bet2.y() *= -1;
    }
    real f[4];
    AuxAngle alpa( kp * fabs(omg1.y()), k * fabs(bet1.x()), true ),
      alpb(alpa);

    TriaxialLineF l(*this, Triaxial::gamblk{}, 0.5, 1.5);
    TriaxialLineF::ics fic(l, bet1, omg1, alpb);
    TriaxialLineF::disttx d;
    bool done = false, equatorial = false;
    real fa, fb;
    if (alpb.y() == 0 && bet1.y() != 0) {
      AuxAngle omg12 = omg2 - omg1;
      if (omg12.y() == 0) {
        if (omg12.x() < 0) {
          alpb.setquadrant(2U);
          fic.setquadrant(l, 2U);
          done = true;
        }
      } else if (omg12.y() > 0) {
        alpa = AuxAngle(numeric_limits<real>::epsilon()/(1<<20), 1, false);
        alpb = AuxAngle(numeric_limits<real>::epsilon()/(1<<20), -1, false);
        fa = - omg12.radians();
        fb = fa + Math::pi();
      } else {
        alpa = AuxAngle(-numeric_limits<real>::epsilon()/(1<<20), -1, false);
        alpb = AuxAngle(-numeric_limits<real>::epsilon()/(1<<20), 1, false);
        fb = omg12.radians();
        fa = fb - Math::pi();
      }
    } else {
      unsigned q = 0U;
      for (; q <= 4U; ++q) {
        if (q) {
          alpb.setquadrant(q);
          fic.setquadrant(l, q);
        }
        if (q < 4U) {
          f[q] = l.Hybrid0(fic, bet2, omg2);
          if (0)
            cout << "UMB " << q << " "
                 << alpb.quadrant() << " "
                 << alpb.degrees() << " " << scientific
                 << f[q]/Math::degree() << "\n";
          if (fabs(f[q]) < numeric_limits<real>::epsilon()) {
            done = true;
            break;
          }
        }
        if (q && (f[(q+3U) & 3U] < 0 && f[q & 3U] > 0)) {
          break;
        }
      }
      if (q > 4U) std::cout << "ERROR\n";
      fa = f[(q+3U) & 3U]; fb = f[q & 3U];
      alpa.setquadrant((q+3U) & 3U);
      if (!done && bet1.y() == 0 && bet2.y() == 0 && (q == 1U || q == 3U)) {
        // cout << "QQ " << q << "\n";
        // Possibly equatorial geodesics
        // The change in f(beta) between a point and the conjugate point is
        int eE = q&2U ? -1 : 1;
        AuxAngle alpx = AuxAngle(eE, 0);
        l = TriaxialLineF(*this, gamblk(*this, bet1, omg1, alpx), 0.5, 1.5);
        fic = TriaxialLineF::ics(l, bet1, omg1, alpx);
        AuxAngle bet2x, omg2x, alp2x;
        l.ArcPos0(fic, Math::pi(), bet2x, omg2x, alp2x);
        // cout << "OMGa " << omg2.degrees() << " " << omg2x.degrees() << "\n";
        omg2x -= omg2;
        // cout << "OMGb " << eE << " " << omg2x.degrees() << "\n";
        if (eE * omg2x.y() >= 0) {
          equatorial = done = true;
          // We have
          //   delta = l.fbet()(v2) - l.fomg()(u2)
          //   u2 = l.fomg().fwd(eE * omg2)
          //   v2 = l.fbet().fwd(psi2)
          // =>
          //   v2 = l.fbet.inv(fic.delta +  l.fomg()(u2));
          real u2=  l.fomg().fwd(eE * ((omg1 - AuxAngle::cardinal(1U)).radians() +
                                       (omg2 - omg1).radians())),
            v2 = l.fbet().inv(fic.delta + l.fomg()(u2));
          d.omgw2 = u2;
          d.betw2 = v2;
          d.ind2 = 0;
          if (0) {
          cout << "POSQ " << omg2.degrees()-90 << " " << l.fomg().rev(u2) << " " << l.fbet().rev(v2) << "\n";
          cout << "POSP " << u2 << " " << v2 << " " << fic.delta << "\n";
          }
          alpb = alpx;
        } else {
          alpx.x() = -numeric_limits<real>::epsilon()/(1<<20);
          (eE > 0 ? alpa : alpb) = alpx;
          (eE > 0 ? fa : fb) = omg2x.radians();
        }
      }
    }
    AuxAngle alp1;
    if (done)
      alp1 = alpb;
    else {
      int countn = 0, countb = 0;
      if (0)
        cout << "ROOT "
             << alpa.degrees() << " " << fa << " "
             << alpb.degrees() << " " << fb << "\n";
      alp1 = findroot(
                      [this, &bet1, &omg1, &bet2, &omg2]
                      (const AuxAngle& alp) -> Math::real
                      {
                        return HybridA(*this,
                                       bet1, omg1, alp,
                                       bet2, omg2);
                      },
                      alpa,  alpb,
                      fa, fb,
                      &countn, &countb);
      l = TriaxialLineF(*this, gamblk(*this, bet1, omg1, alp1), 0.5, 1.5);
      fic = TriaxialLineF::ics(l, bet1, omg1, alp1);
    }

    AuxAngle bet2a, omg2a, alp2;
    if (equatorial)
      alp2 = alp1;
    else
      d = l.Hybrid(fic, bet2, bet2a, omg2a, alp2);
    TriaxialLineG ld(*this, l.gm());
    TriaxialLineG::ics dic(ld, fic);
    real s13 = ld.dist(dic, d);
    (void) AngNorm(bet2a, omg2a, alp2);

    //    cout << "\nFLIP " << swapsign << swap12 << flip1 << "\n";
    // Undo switches in reverse order swapsign, swap12, flip1
    if (swapsign) {
      bet1.y() *= -1;
      bet2.y() *= -1;
      alp1.x() *= -1;
      alp2.x() *= -1;
    }

    if (swap12) {
      swap(bet1, bet2);
      swap(omg1, omg2);
      swap(alp1, alp2);
      alp1 += AuxAngle::cardinal(2U);
      alp2 += AuxAngle::cardinal(2U);
    }

    if (flip1)
      Flip(bet1, omg1, alp1);

    if (flip1 || swap12 || swapsign) {
      fic = TriaxialLineF::ics(l, bet1, omg1, alp1);
      dic = TriaxialLineG::ics(ld, fic);
    }
    dic.s13 = s13;
    return TriaxialLine(l, fic, ld, dic);
  }

  Math::real Triaxial::HybridA(const Triaxial& t,
                               const AuxAngle& bet1, const AuxAngle& omg1,
                               const AuxAngle& alp1,
                               const AuxAngle& bet2, const AuxAngle& omg2) {
    AuxAngle b1{bet1}, o1{omg1}, a1{alp1};
    gamblk gam(t, b1, o1, a1);
    TriaxialLineF l(t, gam, 0.5, 1.5);
    TriaxialLineF::ics ic(l, b1, o1, a1);
    return l.Hybrid0(ic, bet2, omg2);
  }

  // Solve f(alp1) = 0 where alp1 is an azimuth and f(alp1) is the difference in
  // lontitude on bet2 and the target longitude.
  AuxAngle Triaxial::findroot(const function<Math::real(const AuxAngle&)>& f,
                              AuxAngle xa,  AuxAngle xb,
                              Math::real fa, Math::real fb,
                              int* countn, int* countb) {
    // Implement root finding method of Chandrupatla (1997)
    // https://doi.org/10.1016/s0965-9978(96)00051-8
    // Here we follow Scherer (2013), Section 6.1.7.3
    // https://doi.org/10.1007/978-3-319-00401-3

    // Here the independent variable is an AuxAngle, but the computations on this
    // variable essentially involve its conversion to radians.  There's no need
    // to worry about the angle wrapping around because (xb-xa).radians() is in
    // (0,pi).

    // require xa and xb to be normalized (the result is normalized)
    // require fa and fb to have opposite signs

    AuxAngle xm;                  // The return value
    int cntn = 0, cntb = 0;
    if (0) {
      cout << "START\n";
      int num = 360;
      for (int i = -num/2; i <= num/2; ++i) {
        // { int i = 45;
      Math::real x = 360*i/Math::real(num),
      ff = f(AuxAngle::degrees(x));
      cout << "DAT " << x << " " << ff/Math::degree() << "\n";
      }
      return 0;
    }
    bool trip = false, correct = false;
    for (Math::real t = 1/Math::real(2), ab = 0, ft = 0, fm = 0, fc = 0;
         cntn < 50 || GEOGRAPHICLIB_PANIC;) {
      AuxAngle xt = 2*t == 1 ?
        AuxAngle(xa.y() + xb.y(),
                 xa.x() + xb.x(), true) :
        (2*t < 1 ? xa - AuxAngle::radians(t * ab) :
         xb + AuxAngle::radians((1 - t) * ab)),
        xc;
      if (trip) {
        if (correct) xm = xt;   // Update result if not a bisection step
        break;
      }
      ++cntn;
      ft = f(xt);
      if (0)
        cout << scientific
             << "RT " << cntn << " "
             << xt.degrees() << " " << ft/Math::degree() << " " << ab << " " << t << "\n";
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
      ab = (xa-xb).radians();
      Math::real
        ca = (xc-xa).radians(),
        cb = ca+ab,
        // Scherer has a fabs(cb).  This should be fabs(ab).
        tl = numeric_limits<Math::real>::epsilon() / fabs(ab);
      // Backward tests to deal with NaNs
      trip =  !(2 * tl < 1 && fabs(fm) > numeric_limits<Math::real>::epsilon());
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

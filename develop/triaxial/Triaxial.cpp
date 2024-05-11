/**
 * \file Triaxial.cpp
 * \brief Implementation for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "Triaxial.hpp"
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

#include "kissfft.hh"

namespace GeographicLib {

  using namespace std;

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
    real ra = sqrt(Math::sq(rn[0]) + Math::sq(rn[1]) + Math::sq(rn[2]));
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
    f = sqrt(Math::sq(v[0]) + Math::sq(v[1]) + Math::sq(v[2]));
    v = {v[0]/f, v[1]/f, v[2]/f};
  }

  void Triaxial::Norm(vec6& y) const {
    real ra = sqrt(Math::sq(y[0] / axesn[0]) +
                   Math::sq(y[1]) +
                   Math::sq(y[2] /axesn[2]));
    y[0] /= ra; y[1] /= ra; y[2] /= ra;
    vec3 up = {y[0] / axes2n[0], y[1], y[2] / axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      uv = up[0] * y[3+0] + up[1] * y[3+1] + up[2] * y[3+2],
      f = uv/u2;
    y[3+0] -= f * up[0]; y[3+1] -= f * up[1]; y[3+2] -= f * up[2];
    f = sqrt(Math::sq(y[3+0]) + Math::sq(y[3+1]) + Math::sq(y[3+2]));
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

  geod_fun::geod_fun(real kap, real kapp, real eps, real mu) {
    real kp2 = mu > 0 ? mu / (kap + mu) :
      (mu < 0 ? -mu / kap :
       kapp);                   // mu == 0
    geod_fun f(kap, kapp, eps, mu, kp2 < Triaxial::EllipticThresh());
    *this = f;
  }

  geod_fun::geod_fun(real kap, real kapp, real eps, real mu, bool tx)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _tx(tx)
      /*
    , _ell(_tx && _mu != 0 ?
           (_mu > 0 ? _kap / (_kap + _mu) : (_kap + _mu) / _kap) : 0,
           0,
           _tx && _mu != 0 ?
           (_mu > 0 ? _mu / (_kap + _mu) : -_mu / _kap) : 1,
           1)
    , _fun(_mu > 0 ?
           (_tx ?
            TrigfunExt([kap = _kap, kapp = _kapp,
                        eps = _eps, mu = _mu, ell = _ell]
                       (real u) -> real
            { real sn, cn, dn;
              (void) ell.am(u, sn, cn, dn);
              return fup(cn, kap, kapp, eps, mu); },
                       ell.K()) :
            TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                       (real phi) -> real
            { return fphip(cos(phi), kap, kapp, eps, mu); },
                       Math::pi()/2)) :
           (_tx ?
            TrigfunExt([kap = _kap, kapp = _kapp,
                        eps = _eps, mu = _mu, ell = _ell]
                       (real v) -> real
                       { real sn, cn, dn;
                         (void) ell.am(v, sn, cn, dn);
                         return fvp(dn, kap, kapp, eps, mu); },
                       ell.K()) :
            TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                       (real psi) -> real
            { return fpsip(sin(psi), cos(psi), kap, kapp, eps, mu); },
                       Math::pi()/2)))
      */

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
        TrigfunExt([kap = _kap, kapp = _kapp,
                    eps = _eps, mu = _mu, ell = _ell]
                   (real u) -> real
        { real sn, cn, dn;
          (void) ell.am(u, sn, cn, dn);
          return fup(cn, kap, kapp, eps, mu); },
                   _ell.K()) :
        TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                   (real phi) -> real
        { return fphip(cos(phi), kap, kapp, eps, mu); },
                   Math::pi()/2);
    } else if (_mu < 0) {
      _fun = _tx ?
        TrigfunExt([kap = _kap, kapp = _kapp,
                    eps = _eps, mu = _mu, ell = _ell]
                   (real v) -> real
        { real sn, cn, dn;
          (void) ell.am(v, sn, cn, dn);
          return fvp(dn, kap, kapp, eps, mu); },
                   _ell.K()) :
        TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                   (real psi) -> real
        { return fpsip(sin(psi), cos(psi), kap, kapp, eps, mu); },
                   Math::pi()/2);
    } else {                    // _mu == 0
      _fun = _tx ?
        TrigfunExt([kap = _kap, kapp = _kapp,
                    eps = _eps, ell = _ell]
                   (real v) -> real
        { real sn, cn, dn;
          (void) ell.am(v, sn, cn, dn);
          return dfvp(cn, dn, kap, kapp, eps); },
                   2 * _ell.K(), true) :
        TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps]
                   (real phi) -> real
        { return dfp(cos(phi), kap, kapp, eps); },
                   Math::pi(), true);
    }
    if (_mu == 0) {
      _chiinv = _tx ?
        // chi maps [-1,1]*K -> [-1,1]*pi/2
        // dchi/dx = sech(lam(am(x)) - dfv(x)) * (sec(am(x)) * dn(x) - dfvp(x))
        // cosh(lam(am(x))) = 1/cn, sinh(lam(am(x))) = sn/cn
        // cosh(lam(x) - dfv(x))
        //  = cosh(lam(am(x))) * cosh(dfv(x)) - sinh(lam(am(x))) * sinh(dfv(x))
        //  = 1/cn * cosh(dfv(x)) - sn/cn * sinh(dfv(x))
        //  = 1/cn * (cosh(dfv(x)) - sn * sinh(dfv(x)))
        // (sec(am(x)) * dn - dfvp(x)) = 1/cn * (dn - cn * dfvp(x))
        // dchi/dx = (dn - cn * dfvp(x)) /
        //           (cosh(dfv(x)) - sn * sinh(dfv(x)))
        Trigfun::InverseInit(
                             [fun = _fun, ell = _ell]
                             (real x) -> real
                             { return gd(lam(ell.am(x)) - fun(x)); },
                             [fun = _fun, ell = _ell]
                             (real x) -> real
                             { real sn, cn, dn, f = fun(x);
                               (void) ell.am(x, sn, cn, dn);
                               return (dn - cn * fun.deriv(x)) /
                                 (cosh(f) - sn * sinh(f)); },
                             _ell.K(), Math::pi()/2,
                             -Math::pi()/2, Math::pi()/2, &_countn, &_countb) :
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
                             (real x) -> real
                             { return gd(lam(x) - fun(x)); },
                             [fun = _fun]
                             (real x) -> real
                             { real f = fun(x);
                               return (1 - cos(x) * fun.deriv(x)) /
                                 (cosh(f) - sin(x) * sinh(f)); },
                             Math::pi()/2, Math::pi()/2,
                             -Math::pi()/2, Math::pi()/2, &_countn, &_countb);
    }
  }

  geodu_fun::geodu_fun(real kap, real kapp, real eps) {
    geodu_fun f(kap, kapp, eps, kapp < Triaxial::EllipticThresh());
    *this = f;
  }

  geodu_fun::geodu_fun(real kap, real kapp, real eps, bool tx)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _tx(tx)
    , _countn(0)
    , _countb(0)
    , ell(_tx ? _kap : 0, 0,
          _tx ? _kapp : 1, 1)
    , fun(_tx ?
          TrigfunExt([kap = _kap, kapp = _kapp,
                      eps = _eps, ell = ell]
                     (real v) -> real
          { real sn, cn, dn;
            (void) ell.am(v, sn, cn, dn);
            return dfvp(cn, dn, kap, kapp, eps); },
                     2 * ell.K(), true) :
          TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps]
                     (real phi) -> real
          { return dfp(cos(phi), kap, kapp, eps); },
                     Math::pi(), true)
          )
  {
    _chiinv = _tx ?
      // chi maps [-1,1]*K -> [-1,1]*pi/2
      // dchi/dx = sech(lam(am(x)) - dfv(x)) * (sec(am(x)) * dn(x) - dfvp(x))
      // cosh(lam(am(x))) = 1/cn, sinh(lam(am(x))) = sn/cn
      // cosh(lam(x) - dfv(x))
      //  = cosh(lam(am(x))) * cosh(dfv(x)) - sinh(lam(am(x))) * sinh(dfv(x))
      //  = 1/cn * cosh(dfv(x)) - sn/cn * sinh(dfv(x))
      //  = 1/cn * (cosh(dfv(x)) - sn * sinh(dfv(x)))
      // (sec(am(x)) * dn - dfvp(x)) = 1/cn * (dn - cn * dfvp(x))
      // dchi/dx = (dn - cn * dfvp(x)) /
      //           (cosh(dfv(x)) - sn * sinh(dfv(x)))
      Trigfun::InverseInit(
                           [fun = fun, ell = ell]
                           (real x) -> real
                           { return gd(lam(ell.am(x)) - fun(x)); },
                           [fun = fun, ell = ell]
                           (real x) -> real
                           { real sn, cn, dn, f = fun(x);
                             (void) ell.am(x, sn, cn, dn);
                             return (dn - cn * fun.deriv(x)) /
                               (cosh(f) - sn * sinh(f)); },
                           ell.K(), Math::pi()/2,
                           -Math::pi()/2, 0, &_countn, &_countb) :
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
                           [fun = fun]
                           (real x) -> real
                           { return gd(lam(x) - fun(x)); },
                           [fun = fun]
                           (real x) -> real
                           { real f = fun(x);
                             return (1 - cos(x) * fun.deriv(x)) /
                               (cosh(f) - sin(x) * sinh(f)); },
                           Math::pi()/2, Math::pi()/2,
                           -Math::pi()/2, 0, &_countn, &_countb);
  }

  dist_fun::dist_fun(real kap, real kapp, real eps, real mu) {
    real kp2 = mu > 0 ? mu / (kap + mu) :
      (mu < 0 ? -mu / kap :
       kapp);                   // mu == 0
    dist_fun g(kap, kapp, eps, mu, kp2 < Triaxial::EllipticThresh());
    *this = g;
  }

  dist_fun::dist_fun(real kap, real kapp, real eps, real mu, bool tx)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _tx(tx)
    , _ell(_tx ?
          (_mu == 0 ? _kap : (_mu > 0 ? _kap / (_kap + _mu) :
                               (_kap + _mu) / _kap)) : 0,
          0,
          _tx ?
          (_mu == 0 ? _kapp : (_mu > 0 ? _mu / (_kap + _mu) :
                               -_mu / _kap)) : 1,
          1)
    , _fun(_mu == 0 ?
          (_tx ?
           TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, ell = _ell]
                      (real v) -> real
           { real sn, cn, dn;
             (void) ell.am(v, sn, cn, dn);
             return g0vp(cn, kap, kapp, eps);
           },
                      2*_ell.K(), true) :
           TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, ell = _ell]
                      (real phi) -> real
           {
             return g0p(cos(phi), kap, kapp, eps);
           },
                      Math::pi(), true)) :
          (_mu > 0 ?
           (_tx ?
            TrigfunExt([kap = _kap, kapp = _kapp,
                        eps = _eps, mu = _mu, ell = _ell]
                       (real u) -> real
            { real sn, cn, dn;
              (void) ell.am(u, sn, cn, dn);
              return gup(cn, dn, kap, kapp, eps, mu); },
                       _ell.K()) :
            TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                       (real phi) -> real
            { return gphip(cos(phi), kap, kapp, eps, mu); },
                       Math::pi()/2)) :
           (_tx ?
            TrigfunExt([kap = _kap, kapp = _kapp,
                        eps = _eps, mu = _mu, ell = _ell]
                       (real v) -> real
                       { real sn, cn, dn;
                         (void) ell.am(v, sn, cn, dn);
                         return gvp(cn, dn, kap, kapp, eps, mu); },
                       _ell.K()) :
            TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                       (real psi) -> real
            { return gpsip(sin(psi), cos(psi), kap, kapp, eps, mu); },
                       Math::pi()/2))))

  {}
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

} // namespace GeographicLib

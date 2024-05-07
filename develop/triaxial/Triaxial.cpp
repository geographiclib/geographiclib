/**
 * \file Triaxial.cpp
 * \brief Implementation for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2022) <karney@alum.mit.edu> and licensed under
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

  geod_fun::geod_fun(real kap, real kapp, real eps, real mu)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
  {
    real k2 = _mu > 0 ? _kap / (_kap + _mu) : (_kap + _mu) / _kap,
      kp2 = mu > 0 ? mu / (kap + mu) : -mu / kap;
    _tx = kp2 < 0.125;
    ell = EllipticFunction(_tx ? k2 : 0, 0, _tx ? kp2 : 1, 1);
    fun = _mu > 0 ?
      (_tx ?
       TrigfunExt([kap = _kap, kapp = _kapp,
                   eps = _eps, mu = _mu, ell = ell]
                  (real u) -> real
       { return fupf(u, kap, kapp, eps, mu, ell); },
                  ell.K()) :
       TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                  (real phi) -> real
       { return fphipf(phi, kap, kapp, eps, mu); },
                  Math::pi()/2)) :
      (_tx ?
       TrigfunExt([kap = _kap, kapp = _kapp,
                   eps = _eps, mu = _mu, ell = ell]
                  (real v) -> real
       { return fvpf(v, kap, kapp, eps, mu, ell); },
                  ell.K()) :
       TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                  (real psi) -> real
       { return fpsipf(psi, kap, kapp, eps, mu); },
                  Math::pi()/2));
  }

  geod_fun::geod_fun(real kap, real kapp, real eps, real mu, bool tx)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _tx(tx)
    , ell(_tx && _mu != 0 ?
           (_mu > 0 ? _kap / (_kap + _mu) : (_kap + _mu) / _kap) : 0,
           0,
           _tx && _mu != 0 ?
           (_mu > 0 ? _mu / (_kap + _mu) : -_mu / _kap) : 1,
           1)
    , fun(_mu > 0 ?
           (_tx ?
            TrigfunExt([kap = _kap, kapp = _kapp,
                        eps = _eps, mu = _mu, ell = ell]
                       (real u) -> real
                       { return fupf(u, kap, kapp, eps, mu, ell); },
                       ell.K()) :
            TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                       (real phi) -> real
                       { return fphipf(phi, kap, kapp, eps, mu); },
                       Math::pi()/2)) :
           (_tx ?
            TrigfunExt([kap = _kap, kapp = _kapp,
                        eps = _eps, mu = _mu, ell = ell]
                       (real v) -> real
                       { return fvpf(v, kap, kapp, eps, mu, ell); },
                       ell.K()) :
            TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                       (real psi) -> real
                       { return fpsipf(psi, kap, kapp, eps, mu); },
                       Math::pi()/2)))

  {}

  geodu_fun::geodu_fun(real kap, real kapp, real eps, bool tx)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _tx(tx)
    , ell(_tx ? _kap : 0, 0,
          _tx ? _kapp : 1, 1)
    , fun(_tx ?
          TrigfunExt([kap = _kap, kapp = _kapp,
                        eps = _eps, ell = ell]
                       (real v) -> real
                       { return dfvpf(v, kap, kapp, eps, ell); },
                       ell.K()) :
            TrigfunExt([kap = _kap, kapp = _kapp, eps = _eps]
                       (real phi) -> real
                       { return dfpf(phi, kap, kapp, eps); },
                       Math::pi()/2))
  {}

} // namespace GeographicLib

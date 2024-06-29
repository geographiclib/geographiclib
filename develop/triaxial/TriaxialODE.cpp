/**
 * \file Triaxial.cpp
 * \brief Implementation for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "TriaxialODE.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

namespace GeographicLib {

  using namespace std;

  TriaxialODE::TriaxialODE(const Triaxial& t,
                           const Triaxial::vec3& r1, const Triaxial::vec3& v1)
    : _t(t)
    , _b(t.b)
    , _axesn({t.a/t.b, real(1), t.c/t.b})
    , _axes2n({Math::sq(t.a/t.b), real(1), Math::sq(t.c/t.b)})
    , _r1(r1)
    , _v1(v1)
  {
    _t.Norm(_r1, _v1);
  }

  TriaxialODE::TriaxialODE(const Triaxial& t, const AuxAngle& bet1,
                           const AuxAngle& omg1, const AuxAngle& alp1)
    : _t(t)
    , _b(t.b)
    , _axesn({t.a/t.b, real(1), t.c/t.b})
    , _axes2n({Math::sq(t.a/t.b), real(1), Math::sq(t.c/t.b)})
  {
    _t.elliptocart2(bet1, omg1, alp1, _r1, _v1);
  }

  TriaxialODE::TriaxialODE(const Triaxial& t, real bet1, real omg1, real alp1)
    : TriaxialODE(t, AuxAngle::degrees(bet1), AuxAngle::degrees(omg1),
                  AuxAngle::degrees(alp1))
  {}

  void TriaxialODE::Norm(vec6& y) const {
    real ra = Triaxial::hypot3(y[0] / _axesn[0], y[1], y[2] / _axesn[2]);
    y[0] /= ra; y[1] /= ra; y[2] /= ra;
    vec3 up = {y[0] / _axes2n[0], y[1], y[2] / _axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      uv = up[0] * y[3+0] + up[1] * y[3+1] + up[2] * y[3+2],
      f = uv/u2;
    y[3+0] -= f * up[0]; y[3+1] -= f * up[1]; y[3+2] -= f * up[2];
    f = Triaxial::hypot3(y[3+0], y[3+1], y[3+2]);
    y[3+0] /= f; y[3+1] /= f; y[3+2] /= f;
  }

  void TriaxialODE::Norm(vec10& y) const {
    vec6 y6{y[0], y[1], y[2], y[3+0], y[3+1], y[3+2]};
    Norm(y6);
    for (int i = 0; i < 6; ++i) y[i] = y6[i];
  }

  TriaxialODE::vec6 TriaxialODE::Accel(const vec6& y) const {
    vec3 up = {y[0] / _axes2n[0], y[1], y[2] / _axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      f = -(Math::sq(y[3+0]) / _axes2n[0] +
            Math::sq(y[3+1]) +
            Math::sq(y[3+2]) / _axes2n[2]) / u2;
    return vec6{y[3+0], y[3+1], y[3+2], f * up[0], f * up[1], f * up[2]};
  }

  TriaxialODE::vec10 TriaxialODE::Accel(const vec10& y) const {
    vec3 up = {y[0] / _axes2n[0], y[1], y[2] / _axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      f = -(Math::sq(y[3+0]) / _axes2n[0] +
            Math::sq(y[3+1]) +
            Math::sq(y[3+2]) / _axes2n[2]) / u2,
      K = 1/Math::sq(u2 * _axesn[0] * _axesn[2]);
    return vec10{ y[3+0], y[3+1], y[3+2], f * up[0], f * up[1], f * up[2],
                  y[7], -K * y[6], y[9], -K * y[8] };
  }

  int TriaxialODE::Position(real s12, vec3& r2, vec3& v2, real eps) const {
    using namespace boost::numeric::odeint;
    if (eps == 0)
      eps = pow(numeric_limits<real>::epsilon(), real(7)/8);
    // Normalize all distances to b.
    vec6 y{_r1[0]/_b, _r1[1]/_b, _r1[2]/_b, _v1[0], _v1[1], _v1[2]};
    Norm(y);
    int n;
    s12 /= _b;
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
    r2 = {_b*y[0], _b*y[1], _b*y[2]};
    v2 = {y[3], y[4], y[5]};
    return n;
  }

  int TriaxialODE::Position(real s12, vec3& r2, vec3& v2,
                              real& m12, real& M12, real& M21,
                              real eps) const {
    using namespace boost::numeric::odeint;
    if (eps == 0)
      eps = pow(numeric_limits<real>::epsilon(), real(7)/8);
    // Normalize all distances to b.
    vec10 y{_r1[0]/_b, _r1[1]/_b, _r1[2]/_b, _v1[0], _v1[1], _v1[2], 0, 1, 1, 0};
    Norm(y);
    int n;
    s12 /= _b;
    auto fun = [this](const vec10& y, vec10& yp, real /*t*/) -> void {
      yp = y; Norm(yp); yp = Accel(yp);
    };
    if (1) {
      bulirsch_stoer<vec10, real> stepper(eps, real(0));
      n = integrate_adaptive(stepper, fun, y, real(0), s12, 1/real(10));
    } else {
      n = int(round(1/eps));
      runge_kutta4<vec10, real> stepper;
      integrate_n_steps(stepper, fun, y, real(0), s12/n, n);
    }
    Norm(y);
    r2 = {_b*y[0], _b*y[1], _b*y[2]};
    v2 = {y[3], y[4], y[5]};
    m12 = _b*y[6]; M12 = y[8]; M21 = y[7]; // AG Eq 29: dm12/ds2 = M21
    return n;
  }

  void TriaxialODE::Position(real ds, long nmin, long nmax,
                             vector<vec3>& r2, vector<vec3>& v2,
                             real eps) const {
    using namespace boost::numeric::odeint;
    r2.clear(); v2.clear();
    if (nmin > 0 || nmax < 0) return;
    if (nmin == 0 && nmax == 0) {
      r2.push_back(_r1); v2.push_back(_v1);
      return;
    }
    if (eps == 0)
      eps = pow(numeric_limits<real>::epsilon(), real(7)/8);
    auto fun = [this](const vec6& y, vec6& yp, real /*t*/) -> void {
      yp = y; Norm(yp); yp = Accel(yp);
    };
    auto observer = [&r2, &v2, b = _b](const vec6& y, real /*t*/) -> void {
      r2.push_back( {b*y[0], b*y[1], b*y[2]} );
      v2.push_back( {y[3], y[4], y[5]} );
    };
    bulirsch_stoer<vec6, real> stepper(eps, real(0));
    ds /= _b;
    if (nmin < 0) {
      vec6 y{_r1[0]/_b, _r1[1]/_b, _r1[2]/_b, _v1[0], _v1[1], _v1[2]}; Norm(y);
      integrate_n_steps(stepper, fun, y, real(0), -ds, abs(nmin), observer);
      reverse(r2.begin(), r2.end());
      reverse(v2.begin(), v2.end());
    }
    if (nmax > 0) {
      if (nmin < 0) { r2.pop_back(); v2.pop_back(); }
      vec6 y{_r1[0]/_b, _r1[1]/_b, _r1[2]/_b, _v1[0], _v1[1], _v1[2]}; Norm(y);
      integrate_n_steps(stepper, fun, y, real(0), ds, nmax, observer);
    }
    for (int i = 0; i <= nmax - nmin; ++i)
      _t.Norm(r2[i], v2[i]);
  }

  void TriaxialODE::Position(real ds, long nmin, long nmax,
                             vector<vec3>& r2, vector<vec3>& v2,
                             vector<real>& m12,
                             vector<real>& M12, vector<real>& M21,
                             real eps) const {
    using namespace boost::numeric::odeint;
    r2.clear(); v2.clear();
    m12.clear(); M12.clear(); M21.clear();
    if (nmin > 0 || nmax < 0) return;
    if (nmin == 0 && nmax == 0) {
      r2.push_back(_r1); v2.push_back(_v1);
      m12.push_back(0); M12.push_back(1); M21.push_back(0);
      return;
    }
    if (eps == 0)
      eps = pow(numeric_limits<real>::epsilon(), real(7)/8);
    auto fun = [this](const vec10& y, vec10& yp, real /*t*/) -> void {
      yp = y; Norm(yp); yp = Accel(yp);
    };
    auto observer = [&r2, &v2, &m12, &M12, &M21, b = _b]
      (const vec10& y, real /*t*/) -> void {
      r2.push_back( {b*y[0], b*y[1], b*y[2]} );
      v2.push_back( {y[3], y[4], y[5]} );
      m12.push_back( b*y[6] );
      M12.push_back( y[8] );
      M21.push_back( y[7] );
    };
    bulirsch_stoer<vec10, real> stepper(eps, real(0));
    ds /= _b;
    if (nmin < 0) {
      vec10 y{_r1[0]/_b, _r1[1]/_b, _r1[2]/_b, _v1[0], _v1[1], _v1[2],
        0, 1, 1, 0};
      Norm(y);
      integrate_n_steps(stepper, fun, y, real(0), -ds, abs(nmin), observer);
      reverse(r2.begin(), r2.end());
      reverse(v2.begin(), v2.end());
      reverse(m12.begin(), m12.end());
      reverse(M12.begin(), M12.end());
      reverse(M21.begin(), M21.end());
    }
    if (nmax > 0) {
      if (nmin < 0) {
        r2.pop_back(); v2.pop_back();
        m12.pop_back(); M12.pop_back(); M21.pop_back();
      }
      vec10 y{_r1[0]/_b, _r1[1]/_b, _r1[2]/_b, _v1[0], _v1[1], _v1[2],
        0, 1, 1, 0};
      Norm(y);
      integrate_n_steps(stepper, fun, y, real(0), ds, nmax, observer);
    }
    for (int i = 0; i <= nmax - nmin; ++i)
      _t.Norm(r2[i], v2[i]);
  }

  int TriaxialODE::Position(real s12,
                            AuxAngle& bet2, AuxAngle& omg2, AuxAngle& alp2,
                            real eps) const {
    vec3 r2, v2;
    int n = Position(s12, r2, v2, eps);
    _t.cart2toellip(r2, v2, bet2, omg2, alp2);
    return n;
  }

  int TriaxialODE::Position(real s12,
                            AuxAngle& bet2, AuxAngle& omg2, AuxAngle& alp2,
                            real& m12, real& M12, real& M21,
                            real eps) const {
    vec3 r2, v2;
    int n = Position(s12, r2, v2, m12, M12, M21, eps);
    _t.cart2toellip(r2, v2, bet2, omg2, alp2);
    return n;
  }
    
  void TriaxialODE::Position(real ds, long nmin, long nmax,
                             vector<AuxAngle>& bet2, vector<AuxAngle>& omg2,
                             vector<AuxAngle>& alp2, real eps) const {
    vector<vec3> r2, v2;
    Position(ds, nmin, nmax, r2, v2, eps);
    size_t n = r2.size();
    bet2.resize(n); omg2.resize(n); alp2.resize(n);
    for (size_t i = 0; i < n; ++i)
      _t.cart2toellip(r2[i], v2[i], bet2[i], omg2[i], alp2[i]);
  }
  void TriaxialODE::Position(real ds, long nmin, long nmax,
                             vector<AuxAngle>& bet2, vector<AuxAngle>& omg2,
                             vector<AuxAngle>& alp2,
                             vector<real>& m12,
                             vector<real>& M12, vector<real>& M21,
                             real eps) const {
    vector<vec3> r2, v2;
    Position(ds, nmin, nmax, r2, v2, m12, M12, M21, eps);
    size_t n = r2.size();
    bet2.resize(n); omg2.resize(n); alp2.resize(n);
    for (size_t i = 0; i < n; ++i)
      _t.cart2toellip(r2[i], v2[i], bet2[i], omg2[i], alp2[i]);
  }

} // namespace GeographicLib

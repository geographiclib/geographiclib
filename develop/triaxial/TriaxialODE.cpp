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
#include <limits>
#include <functional>

namespace GeographicLib {

  using namespace std;

  TriaxialODE::TriaxialODE(const Triaxial& t,
                           Triaxial::vec3 r1, Triaxial::vec3 v1,
                           bool extended, real eps)
    : _t(t)
    , _b(t.b())
    , _eps(eps <= 0 ? pow(numeric_limits<real>::epsilon(), real(7)/8) : eps)
    , _axesn({t.a()/t.b(), real(1), t.c()/t.b()})
    , _axes2n({Math::sq(_axesn[0]), real(1), Math::sq(_axesn[2])})
    , _r1(r1)
    , _v1(v1)
    , _extended(extended)
    , _dir(0)
    , _nsteps(0)
    , _step6(_eps, real(0))
    , _step10(_eps, real(0))
  {
    _t.Norm(_r1, _v1);
  }

  TriaxialODE::TriaxialODE(const Triaxial& t, Angle bet1, Angle omg1,
                           Angle alp1,
                           bool extended, real eps)
    : _t(t)
    , _b(t.b())
    , _eps(eps <= 0 ? pow(numeric_limits<real>::epsilon(), real(7)/8) : eps)
    , _axesn({t.a()/t.b(), real(1), t.c()/t.b()})
    , _axes2n({Math::sq(_axesn[0]), real(1), Math::sq(_axesn[2])})
    , _extended(extended)
    , _dir(0)
    , _nsteps(0)
    , _step6(_eps, real(0))
    , _step10(_eps, real(0))
  {
    _t.elliptocart2(bet1, omg1, alp1, _r1, _v1);

  }

  TriaxialODE::TriaxialODE(const Triaxial& t, real bet1, real omg1, real alp1,
                           bool extended, real eps)
    : TriaxialODE(t, Angle(bet1), Angle(omg1), Angle(alp1), extended, eps)
  {}

  void TriaxialODE::Norm(vec6& y) const {
    real ra = Math::hypot3(y[0] / _axesn[0], y[1], y[2] / _axesn[2]);
    y[0] /= ra; y[1] /= ra; y[2] /= ra;
    vec3 up = {y[0] / _axes2n[0], y[1], y[2] / _axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      uv = up[0] * y[3+0] + up[1] * y[3+1] + up[2] * y[3+2],
      f = uv/u2;
    y[3+0] -= f * up[0]; y[3+1] -= f * up[1]; y[3+2] -= f * up[2];
    f = Math::hypot3(y[3+0], y[3+1], y[3+2]);
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

  bool TriaxialODE::Position(real s12, vec3& r2, vec3& v2) {
    if (_extended) {
      real m12, M12, M21;
      return Position(s12, r2, v2, m12, M12, M21);
    }
    s12 /= _b;
    if (s12 == 0) {
      r2 = _r1; v2 = _v1;
    }
    auto fun = [this](const vec6& y, vec6& yp, real /*t*/) -> void {
      yp = y; Norm(yp); yp = Accel(yp);
    };
    if (_dir == 0) {
      _dir = s12 < 0 ? -1 : 1;
      vec6 y{_r1[0] / _b, _r1[1] / _b, _r1[2] / _b, _v1[0], _v1[1], _v1[2]};
      _step6.initialize(y, real(0), _dir / real(4));
      (void) _step6.do_step(fun);
      ++_nsteps;
    }
    if (_dir * s12 < _dir * _step6.previous_time())
      return false;
    while (_dir * _step6.current_time() < _dir * s12) {
      (void) _step6.do_step(fun);
      ++_nsteps;
    }
    vec6 y;
    _step6.calc_state(s12, y);
    Norm(y);
    r2 = {_b * y[0], _b * y[1], _b * y[2]};
    v2 = {y[3], y[4], y[5]};
    return true;
  }

  bool TriaxialODE::Position(real s12, vec3& r2, vec3& v2,
                             real& m12, real& M12, real& M21) {
    if (!_extended) return false;
    s12 /= _b;
    if (s12 == 0) {
      r2 = _r1; v2 = _v1;
      m12 = 0; M12 = M21 = 1;
    }
    auto fun = [this](const vec10& y, vec10& yp, real /*t*/) -> void {
      yp = y; Norm(yp); yp = Accel(yp);
    };
    if (_dir == 0) {
      _dir = s12 < 0 ? -1 : 1;
      vec10 y{_r1[0] / _b, _r1[1] / _b, _r1[2] / _b, _v1[0], _v1[1], _v1[2],
        0, 1, 1, 0};
      _step10.initialize(y, real(0), _dir / real(4));
      (void) _step10.do_step(fun);
      ++_nsteps;
    }
    if (_dir * s12 < _dir * _step10.previous_time())
      return false;
    while (_dir * _step10.current_time() < _dir * s12) {
      (void) _step10.do_step(fun);
      ++_nsteps;
    }
    vec10 y;
    _step10.calc_state(s12, y);
    Norm(y);
    r2 = {_b * y[0], _b * y[1], _b * y[2]};
    v2 = {y[3], y[4], y[5]};
    m12 = _b * y[6]; M12 = y[8]; M21 = y[7]; // AG Eq 29: dm12/ds2 = M21
    return true;
  }

  bool TriaxialODE::Position(real s12,
                             Angle& bet2, Angle& omg2, Angle& alp2) {
    vec3 r2, v2;
    if (!Position(s12, r2, v2)) return false;
    _t.cart2toellip(r2, v2, bet2, omg2, alp2);
    return true;
  }

  bool TriaxialODE::Position(real s12,
                              Angle& bet2, Angle& omg2, Angle& alp2,
                              real& m12, real& M12, real& M21) {
    vec3 r2, v2;
    if (!Position(s12, r2, v2, m12, M12, M21)) return false;
    _t.cart2toellip(r2, v2, bet2, omg2, alp2);
    return true;
  }

  void TriaxialODE::Reset() { _dir = 0; _nsteps = 0; }

} // namespace GeographicLib

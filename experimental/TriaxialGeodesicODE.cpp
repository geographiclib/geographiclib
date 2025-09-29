/**
 * \file TriaxialGeodesicODE.cpp
 * \brief Implementation for GeographicLib::experiemental::TriaxialGeodesicODE
 *   class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <iostream>
#include <iomanip>

#if defined(_MSC_VER)
// Squelch warning triggered by boost:
//   4127: conditional expression is constant
//   4244: conversion from '_Ty' to 'double'
//   4267: conversion from 'size_t' to 'unsigned short'
#  pragma warning (disable: 4127 4244 4267)
#endif
#include "TriaxialGeodesicODE.hpp"

namespace GeographicLib {
namespace experimental {

  using namespace std;
  using namespace Triaxial;

  TriaxialGeodesicODE::TriaxialGeodesicODE(const Ellipsoid3& t,
                                           bool extended, bool dense,
                                           bool normp, real eps)
    : TriaxialGeodesicODE(t, vec3{t.a(), 0, 0}, vec3{0, 0, 1},
                          extended, dense, normp, eps)
  {}

  TriaxialGeodesicODE::TriaxialGeodesicODE(const Ellipsoid3& t,
                                           Ellipsoid3::vec3 R1,
                                           Ellipsoid3::vec3 V1,
                                           bool extended, bool dense,
                                           bool normp, real eps)
    : _t(t)
    , _b(t.b())
    , _eps(eps <= 0 ? pow(numeric_limits<real>::epsilon(),
                          dense ? real(3)/4 : real(9)/10) : eps)
    , _axesn({t.a()/t.b(), real(1), t.c()/t.b()})
    , _axes2n({Math::sq(_axesn[0]), real(1), Math::sq(_axesn[2])})
    , _r1(R1)
    , _v1(V1)
    , _extended(extended)
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
    , _dense(dense)
#else
    , _dense(false && dense)
#endif
    , _normp(normp)
    , _dir(0)
    , _nsteps(0)
    , _intsteps(0)
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
      // eps_abs, eps_rel, factor_x, factor_dxdt, max_dt, control_interpolation
    , _dstep6(_eps, real(0), real(1), real(1), real(0), true)
    , _dstep10(_eps, real(0), real(1), real(1), real(0), true)
#endif
      // eps_abs, eps_rel, factor_x = 1, factor_dxdt = 1, max_dt = 0
    , _step6(_eps, real(0))
    , _step10(_eps, real(0))
  {
    _t.Norm(_r1, _v1);
    _t.cart2toellip(_r1, _v1, _bet1, _omg1, _alp1);
  }

  TriaxialGeodesicODE::TriaxialGeodesicODE(const Ellipsoid3& t,
                                           Angle bet1, Angle omg1,
                                           Angle alp1,
                                           bool extended, bool dense,
                                           bool normp, real eps)
    : _t(t)
    , _b(t.b())
    , _eps(eps <= 0 ? pow(numeric_limits<real>::epsilon(),
                          dense ? real(3)/4 : real(9)/10) : eps)
    , _axesn({t.a()/t.b(), real(1), t.c()/t.b()})
    , _axes2n({Math::sq(_axesn[0]), real(1), Math::sq(_axesn[2])})
    , _bet1(bet1)
    , _omg1(omg1)
    , _alp1(alp1)
    , _extended(extended)
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
    , _dense(dense)
#else
    , _dense(false && dense)
#endif
    , _normp(normp)
    , _dir(0)
    , _nsteps(0)
    , _intsteps(0)
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
      // eps_abs, eps_rel, factor_x, factor_dxdt, max_dt, control_interpolation
    , _dstep6(_eps, real(0), real(1), real(1), real(0), true)
    , _dstep10(_eps, real(0), real(1), real(1), real(0), true)
#endif
      // eps_abs, eps_rel, factor_x = 1, factor_dxdt = 1, max_dt = 0
    , _step6(_eps, real(0))
    , _step10(_eps, real(0))
  {
    _t.elliptocart2(_bet1, _omg1, _alp1, _r1, _v1);
  }

  void TriaxialGeodesicODE::Norm6(vec6& y) const {
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

  void TriaxialGeodesicODE::Accel6(const vec6& y, vec6& yp) const {
    vec3 up = {y[0] / _axes2n[0], y[1], y[2] / _axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      f = -(Math::sq(y[3+0]) / _axes2n[0] +
            Math::sq(y[3+1]) +
            Math::sq(y[3+2]) / _axes2n[2]) / u2;
    yp = vec6{y[3+0], y[3+1], y[3+2], f * up[0], f * up[1], f * up[2]};
    ++_intsteps;
  }

  void TriaxialGeodesicODE::Accel6N(const vec6& y, vec6& yp) const {
    real ra = Math::hypot3(y[0] / _axesn[0], y[1], y[2] / _axesn[2]);
    // merge correction to r with computation of U and put U in yp[3,4,5]
    yp[3+0] = y[0] / (ra * _axes2n[0]);
    yp[3+1] = y[1] / ra;
    yp[3+2] = y[2] / (ra * _axes2n[2]);
    real u2 = Math::sq(yp[3+0]) + Math::sq(yp[3+1]) + Math::sq(yp[3+2]),
      uv = yp[3+0] * y[3+0] + yp[3+1] * y[3+1] + yp[3+2] * y[3+2],
      f = uv/u2;
    yp[0] = y[3+0] - f * yp[3+0];
    yp[1] = y[3+1] - f * yp[3+1];
    yp[2] = y[3+2] - f * yp[3+2];
    f = Math::hypot3(yp[0], yp[1], yp[2]);
    // the corrected v (normal to U and with unit magnitude) in yp[0,1,2]
    yp[0] /= f; yp[1] /= f; yp[2] /= f;
    f = -(Math::sq(yp[0]) / _axes2n[0] +
          Math::sq(yp[1]) +
          Math::sq(yp[2]) / _axes2n[2]) / u2;
    yp[3+0] *= f; yp[3+1] *= f; yp[3+2] *= f;
    ++_intsteps;
  }

  void TriaxialGeodesicODE::Norm10(vec10& y) const {
    vec6 y6{y[0], y[1], y[2], y[3+0], y[3+1], y[3+2]};
    Norm6(y6);
    for (int i = 0; i < 6; ++i) y[i] = y6[i];
  }

  void TriaxialGeodesicODE::Accel10(const vec10& y, vec10& yp) const {
    vec3 up = {y[0] / _axes2n[0], y[1], y[2] / _axes2n[2]};
    real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]),
      f = -(Math::sq(y[3+0]) / _axes2n[0] +
            Math::sq(y[3+1]) +
            Math::sq(y[3+2]) / _axes2n[2]) / u2,
      K = 1/Math::sq(u2 * _axesn[0] * _axesn[2]);
    yp =vec10{
      y[3+0], y[3+1], y[3+2],
      f * up[0], f * up[1], f * up[2],
      y[7], -K * y[6], y[9], -K * y[8] };
    ++_intsteps;
  }

  void TriaxialGeodesicODE::Accel10N(const vec10& y, vec10& yp) const {
    real ra = Math::hypot3(y[0] / _axesn[0], y[1], y[2] / _axesn[2]);
    // merge correction to r with computation of U and put U in yp[3,4,5]
    yp[3+0] = y[0] / (ra * _axes2n[0]);
    yp[3+1] = y[1] / ra;
    yp[3+2] = y[2] / (ra * _axes2n[2]);
    real u2 = Math::sq(yp[3+0]) + Math::sq(yp[3+1]) + Math::sq(yp[3+2]),
      uv = yp[3+0] * y[3+0] + yp[3+1] * y[3+1] + yp[3+2] * y[3+2],
      f = uv/u2,
      K = 1/Math::sq(u2 * _axesn[0] * _axesn[2]);
    yp[0] = y[3+0] - f * yp[3+0];
    yp[1] = y[3+1] - f * yp[3+1];
    yp[2] = y[3+2] - f * yp[3+2];
    f = Math::hypot3(yp[0], yp[1], yp[2]);
    // the corrected v (normal to U and with unit magnitude) in yp[0,1,2]
    yp[0] /= f; yp[1] /= f; yp[2] /= f;
    f = -(Math::sq(yp[0]) / _axes2n[0] +
          Math::sq(yp[1]) +
          Math::sq(yp[2]) / _axes2n[2]) / u2;
    yp[3+0] *= f; yp[3+1] *= f; yp[3+2] *= f;
    yp[6+0] = y[6+1]; yp[6+1] = -K * y[6+0];
    yp[8+0] = y[8+1]; yp[8+1] = -K * y[8+0];
    ++_intsteps;
  }

  pair<Math::real, Math::real>
  TriaxialGeodesicODE::Position(real s12, vec3& R2, vec3& V2) {
    real m12, M12, M21;
    return Position(s12, R2, V2, m12, M12, M21);
  }

  pair<Math::real, Math::real>
  TriaxialGeodesicODE::Position(real s12, vec3& R2, vec3& V2,
                                real& m12, real& M12, real& M21) {
    const auto
      fun6 = [this](const vec6& y, vec6& yp, real /*t*/) -> void {
        return _normp ? Accel6N(y, yp) : Accel6(y, yp);
      };
    const auto
      fun10 = [this](const vec10& y, vec10& yp, real /*t*/) -> void {
        return _normp ? Accel10N(y, yp) : Accel10(y, yp);
      };
    s12 /= _b;
    if (s12 == 0) {
      R2 = _r1; V2 = _v1;
      if (_extended) {
        m12 = copysign(real(0), s12); M12 = M21 = 1;
      }
      return pair<real, real>(0, 0);
    } else if (!isfinite(s12)) {
      R2 = V2 = {Math::NaN(), Math::NaN(), Math::NaN()};
      if (_extended) {
        m12 = M12 = M21 = Math::NaN();
      }
      return pair<real, real>(Math::NaN(), Math::NaN());
    }
    if (_dir == 0) {
      _dir = signbit(s12) ? -1 : 1;
      if (_extended)
        _y10 = {_r1[0] / _b, _r1[1] / _b, _r1[2] / _b,
          _dir * _v1[0], _dir * _v1[1], _dir * _v1[2],
          0, 1, 1, 0};
      else
        _y6 = vec6{_r1[0] / _b, _r1[1] / _b, _r1[2] / _b,
          _dir * _v1[0], _dir * _v1[1], _dir * _v1[2]};
      _s = 0;
      if (_dense) {
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
        if (_extended) {
          // x0, t0, dt0
          _dstep10.initialize(_y10, _s, 1/real(4));
          (void) _dstep10.do_step(fun10);
        } else {
          _dstep6.initialize(_y6, _s, 1/real(4));
          (void) _dstep6.do_step(fun6);
        }
        ++_nsteps;
#endif
      }
    }
    s12 *= _dir;
    if (_dense) {
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
      if (s12 < (_extended ?
                 _dstep10.previous_time() :
                 _dstep6.previous_time())) {
        s12 *= _dir;
        Reset();
        return Position(s12, R2, V2, m12, M12, M21);
      }
      while ((_extended ? _dstep10.current_time() : _dstep6.current_time())
             < s12) {
        if (_extended)
          (void) _dstep10.do_step(fun10);
        else
          (void) _dstep6.do_step(fun6);
        ++_nsteps;
      }
      if (_extended)
        _dstep10.calc_state(s12, _y10);
      else
        _dstep6.calc_state(s12, _y6);
#endif
    } else {
      if (s12 < _s) {
        s12 *= _dir;
        Reset();
        return Position(s12, R2, V2, m12, M12, M21);
      } else if (s12 > _s) {
        _nsteps += long(_extended ?
                        integrate_adaptive(_step10, fun10, _y10, _s, s12,
                                           fmax(s12 - _s, 1/real(4))) :
                        integrate_adaptive(_step6, fun6, _y6, _s, s12,
                                           fmax(s12 - _s, 1/real(4))));
        _s = s12;
      }
    }
    real errr, errv;
    if (_extended) {
      vec10 t = _y10; Norm10(t);
      R2 = {_b * t[0], _b * t[1], _b * t[2]};
      V2 = {_dir * t[3], _dir * t[4], _dir * t[5]};
      m12 = _dir * _b * t[6];
      M12 = t[8]; M21 = t[7];     // AG Eq 29: dm12/ds2 = M21
      for (int i = 0; i < 6; ++i)
        t[i] -= _y10[i];
      errr = _b * Math::hypot3(t[0], t[1], t[2]);
      errv = Math::hypot3(t[3+0], t[3+1], t[3+2]);
    } else {
      vec6 t = _y6; Norm6(t);
      R2 = {_b * t[0], _b * t[1], _b * t[2]};
      V2 = {_dir * t[3], _dir * t[4], _dir * t[5]};
      for (int i = 0; i < 6; ++i)
        t[i] -= _y6[i];
      errr = _b * Math::hypot3(t[0], t[1], t[2]);
      errv = Math::hypot3(t[3+0], t[3+1], t[3+2]);
    }
    return pair<real, real>(errr, errv);
  }

  pair<Math::real, Math::real>
  TriaxialGeodesicODE::Position(real s12,
                                Angle& bet2, Angle& omg2, Angle& alp2) {
    vec3 R2, V2;
    auto err = Position(s12, R2, V2);
    _t.cart2toellip(R2, V2, bet2, omg2, alp2);
    return err;
  }

  pair<Math::real, Math::real>
  TriaxialGeodesicODE::Position(real s12,
                                Angle& bet2, Angle& omg2, Angle& alp2,
                                real& m12, real& M12, real& M21) {
    vec3 R2, V2;
    auto err = Position(s12, R2, V2, m12, M12, M21);
    _t.cart2toellip(R2, V2, bet2, omg2, alp2);
    return err;
  }

  vector<size_t> TriaxialGeodesicODE::sort_indices(const vector<real>& v) {
    // Return a vector of indices into v in "sorted" order:
    //   nans and infs in any order
    //   zeros in any order
    //   negative values in descending order (away from zero)
    //   positive values in ascending order
    //
    // Adapted from Lukasz Wiklendt https://stackoverflow.com/a/12399290/837055
    // Initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    const auto cmp = [](real x, real y) -> bool
    {
      const auto t = [](real x) -> int
      {
        // assign categoried to numbers
        // 0 = nans and infs
        // 1 = zeros
        // 2 = negative
        // 3 = positive
        return !isfinite(x) ? 0 :
          x == 0 ? 1 :
          signbit(x) ? 2 : 3;
      };
      int tx = t(x), ty = t(y);
      // the "less-than" comparison: compare on category first and if
      // categories are equal on absolute values of finite non-zero numbers
      return tx < ty || (tx == ty && tx >= 2 && fabs(x) < fabs(y));
    };
    sort(idx.begin(), idx.end(),
         [&v, &cmp](size_t ix, size_t iy) -> bool
         { return cmp(v[ix], v[iy]); });
    return idx;
  }

  void TriaxialGeodesicODE::Position(const vector<real>& s12,
                                     vector<vec3>& R2, vector<vec3>& V2) {
    vector<real> m12, M12, M21;
    Position(s12, R2, V2, m12, M12, M21);
  }

  void TriaxialGeodesicODE::Position(const vector<real>& s12,
                                     vector<vec3>& R2, vector<vec3>& V2,
                                     vector<real>& m12,
                                     vector<real>& M12, vector<real>& M21) {
    size_t n = s12.size();
    R2.resize(n); V2.resize(n);
    if (_extended) {
      m12.resize(n); M12.resize(n); M21.resize(n);
    }
    vector<size_t> idx = sort_indices(s12);
    for (size_t i = 0; i < n; ++i) {
      size_t k = idx[i];
      if (_extended)
        (void) Position(s12[k], R2[k], V2[k], m12[k], M12[k], M21[k]);
      else
        (void) Position(s12[k], R2[k], V2[k]);
    }
  }

  void TriaxialGeodesicODE::Position(const vector<real>& s12,
                                     vector<Angle>& bet2, vector<Angle>& omg2,
                                     vector<Angle>& alp2) {
    vector<vec3> R2, V2;
    Position(s12, R2, V2);
    size_t n = s12.size();
    bet2.resize(n); omg2.resize(n); alp2.resize(n);
    for (size_t i = 0; i < s12.size(); ++i)
      _t.cart2toellip(R2[i], V2[i], bet2[i], omg2[i], alp2[i]);
  }

  void TriaxialGeodesicODE::Position(const vector<real>& s12,
                                     vector<Angle>& bet2, vector<Angle>& omg2,
                                     vector<Angle>& alp2,
                                     vector<real>& m12,
                                     vector<real>& M12, vector<real>& M21) {
    vector<vec3> R2, V2;
    Position(s12, R2, V2, m12, M12, M21);
    size_t n = s12.size();
    bet2.resize(n); omg2.resize(n); alp2.resize(n);
    for (size_t i = 0; i < s12.size(); ++i)
      _t.cart2toellip(R2[i], V2[i], bet2[i], omg2[i], alp2[i]);
  }

  void TriaxialGeodesicODE::Reset() { _dir = 0; }

  void TriaxialGeodesicODE::Reset(vec3 r1, vec3 v1) {
    _r1 = r1;
    _v1 = v1;
    _t.Norm(_r1, _v1);
    _t.cart2toellip(_r1, _v1, _bet1, _omg1, _alp1);
    _nsteps = 0; _intsteps = 0;
    Reset();
  }

  void TriaxialGeodesicODE::Reset(Angle bet1, Angle omg1, Angle alp1) {
    _bet1 = bet1; _omg1 = omg1; _alp1 = alp1;
    _t.elliptocart2(_bet1, _omg1, _alp1, _r1, _v1);
    Reset();
  }

  pair<Math::real, Math::real> TriaxialGeodesicODE::CurrentDistance() const {
    if (_dir == 0)
      return pair<real, real>(0, 0);
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
    else if (_dense) {
      return pair<real, real>(_dir * _dstep10.previous_time(),
                              _dir * _dstep10.current_time());
    }
#endif
    else
      return pair<real, real>(_dir * _s, _dir * _s);
  }

} // namespace experimental
} // namespace GeographicLib

/**
 * \file TriaxialLine.cpp
 * \brief Implementation for GeographicLib::TriaxialLine class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "TriaxialLine.hpp"
#include <iostream>

namespace GeographicLib {

  using namespace std;

  TriaxialLine::TriaxialLine(const Triaxial& t,
                             const AuxAngle& bet1,
                             const AuxAngle& omg1,
                             const AuxAngle& alp1)
    : _t(t)
    , _bet1(bet1.normalized())
      // omg1 - 90
    , _omg1(AuxAngle(-omg1.x(), omg1.y()).normalized())
    , _alp1(alp1.normalized())
    , _gam(_t.gamma(_bet1.x(), _omg1.x(), _alp1.y(), _alp1.x()))
    , _f(_t.k2, _t.kp2, _t.e2, _gam, 0.5, 1.5) // 0.333, 4)
    , _umbalt(false)
    , _distinit(false)
  {
    _f.fbet().ComputeInverse();
    _f.fomg().ComputeInverse();
    const real eps = numeric_limits<real>::epsilon();
    if (_bet1.y() == 0 && _alp1.x() == 0)
      _alp1.x() = - Math::sq(eps);
    _E = signbit(_alp1.y()) ? -1 : 1;
    _N = signbit(_alp1.x()) ? -1 : 1;
    if (_gam > 0) {
      _flip = signbit(bet1.x()) ? -1 : 1;
      _bet0 = _flip;
      _alp0 = 1;
      _psi1 = AuxAngle( _t.k2 * bet1.y(),
                        _flip * _alp1.x() *
                        hypot(_t.k * _bet1.x(), _t.kp * _omg1.x()) );
      _v0 = fbet().fwd(_psi1.radians());
      _u0 = fomg().fwd(_E * _omg1.radians());
      _delta = fbet()(_v0) - fomg()(_u0);
    } else if (_gam < 0) {
      _flip = _omg1.x() < 0 ? -1 : 1;
      _omg0 = _flip;
      _alp0 = _N;
      _psi1 = AuxAngle( _t.kp2 * _omg1.y(),
                        _flip * alp1.y() *
                        hypot(_t.k * _bet1.x(), _t.kp * _omg1.x()) );
      _v0 = fomg().fwd(_psi1.radians());
      _u0 = fbet().fwd(_N * _bet1.radians());
      _delta = fbet()(_u0) - fomg()(_v0);
    } else {                    // _gam == 0
      _alp0 = _umbalt ? _N : 1;
                                // N.B. factor of k*kp omitted
      _df = fbet().Max() - fomg().Max();
      _deltashift = 2*_df - log(_t.k2/_t.kp2);
      _deltamax = -2*log(eps/2); // = asinh(2/eps^2)

      if (abs(_bet1.x()) < 8*eps && abs(_omg1.x()) < 8*eps) {
        _bet0 = (int(round(_bet1.y())) + _N) % 2 ? -1 : 1;
        _omg0 = (int(round(_omg1.y())) + _E) % 2 ? -1 : 1;
        _delta = _deltashift/2 - log(abs(_alp1.tan()));
      } else {
        _bet0 = signbit(_bet1.x()) ? -1 : 1;
        _omg0 = signbit(_omg1.x()) ? -1 : 1;
        /*
        _delta = _N * fbet()(fbet().fwd(atan2(_bet1.y() * _bet0,
                                              _bet1.x() * _bet0))) -
          _E * fomg()(fomg().fwd(atan2(_omg1.y() * _omg0,
                                       _omg1.x() * _omg0)));
        */
        _delta = _N * fbet()(asinh(_bet1.tan())) -
          _E * fomg()(asinh(_omg1.tan()));

      }
    }

    bool debug = true;

    if (debug) {
      distinit();
      real u = 1;
      if (1) {
        /*
          real du = sqrt(sqrt(numeric_limits<real>::epsilon()));
      cout << "gam " << _gam << "\n";
      cout << "fbet " << fbet().txp() << " " << fbet()(u) << " "
           << fbet().deriv(u) - (fbet()(u+du/2) - fbet()(u-du/2))/du << " "
           << fbet().inv(fbet()(u)) - u << "\n";
      cout << "fomg " << fomg().txp() << " " << fomg()(u) << " "
           << fomg().deriv(u) - (fomg()(u+du/2) - fomg()(u-du/2))/du << " "
           << fomg().inv(fomg()(u)) - u << "\n";
      cout << "gbet " << gbet().txp() << " " << gbet()(u) << " "
           << gbet().deriv(u) - (gbet()(u+du/2) - gbet()(u-du/2))/du << "\n";
      cout << "gomg " << gomg().txp() << " " << gomg()(u) << " "
           << gomg().deriv(u) - (gomg()(u+du/2) - gomg()(u-du/2))/du << "\n";
        */
        cout << "cnts " << _gam << " "
             << fbet().NCoeffs() << " " << fbet().NCoeffsInv() << " "
             << fomg().NCoeffs() << " " << fomg().NCoeffsInv() << " "
             << gbet().NCoeffs() << " " << gomg().NCoeffs() << "\n";
        cout << "invs "
             << fbet().inv(fbet()(u)) - u << " "
             << fomg().inv(fomg()(u)) - u << "\n";
      }
      if (0) {
      cout << "gam " << _gam << "\n";
      cout << "fbet " << " " << fbet().txp() << gbet().txp() << " "
           << gbet().deriv(u) / fbet().deriv(u) - gbet().gfderiv(u)
           << "\n";
      cout << "fomg " << " " << fomg().txp() << gomg().txp() << " "
           << gomg().deriv(u) / fomg().deriv(u) - gomg().gfderiv(u)
           << "\n";
      }
      if (0) {
        cout << "gam " << _gam << "\n";
        geod_fun fbeta(_t.k2, _t.kp2, _t.e2, -_gam, false, 0.333, 4);
        geod_fun fomga(_t.kp2, _t.k2, -_t.e2, _gam, false, 0.333, 4);
        dist_fun gbeta(_t.k2, _t.kp2, _t.e2, -_gam, false);
        dist_fun gomga(_t.kp2, _t.k2, -_t.e2, _gam, false);
        geod_fun fbetb(_t.k2, _t.kp2, _t.e2, -_gam, true, 0.333, 4);
        geod_fun fomgb(_t.kp2, _t.k2, -_t.e2, _gam, true, 0.333, 4);
        dist_fun gbetb(_t.k2, _t.kp2, _t.e2, -_gam, true);
        dist_fun gomgb(_t.kp2, _t.k2, -_t.e2, _gam, true);
        fbeta.ComputeInverse();
        fbetb.ComputeInverse();
        fomga.ComputeInverse();
        fomgb.ComputeInverse();
        real x = fbeta(fbeta.fwd(u));
        cout << "fbet "
             << fbeta(fbeta.fwd(u)) -  fbetb(fbetb.fwd(u)) << " "
             << fbeta.rev(fbeta.inv(x)) -  fbetb.rev(fbetb.inv(x)) << "\n";
        x = fomga(fomga.fwd(u));
        cout << "fomg "
             << fomga(fomga.fwd(u)) -  fomgb(fomgb.fwd(u)) << " "
             << fomga.rev(fomga.inv(x)) -  fomgb.rev(fomgb.inv(x)) << "\n";
        cout << "gbet "
             << gbeta(fbeta.fwd(u)) -  gbetb(fbetb.fwd(u)) << "\n";
        cout << "gomg "
             << gomga(fomga.fwd(u)) -  gomgb(fomgb.fwd(u)) << "\n";
      }
    }
  }
  void TriaxialLine::distinit() {
    if (_distinit) return;
    _g = TriaxialLineG(_t.k2, _t.kp2, _t.e2, _gam);
    if (_gam > 0) {
      _sig1 = gbet()(_v0) + gomg()(_u0);
      _s0 = 0;
    } else if (_gam < 0) {
      _sig1 = gbet()(_u0) + gomg()(_v0);
      _s0 = 0;
    } else {                    // _gam == 0
      /*
      _sig1 = _N * gbet()(gbet().fwd(atan2(_bet1.y() * _bet0,
                                           _bet1.x() * _bet0))) -
        _E * gomg()(gomg().fwd(atan2(_omg1.y() * _omg0,
                                     _omg1.x() * _omg)));
      */
      _sig1 = _N * gbet()(asinh(_bet1.tan())) -
        _E * gomg()(asinh(_omg1.tan()));
      _s0 = gbet().Max() + gomg().Max();
    }

    _distinit = true;
  }

  void TriaxialLine::Position(real s12,
                              AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a)
  const {
    // Compute points at distance s12
    real sig2 = _sig1 + s12/_t.b;
    real bet2, omg2, alp2;
    if (_gam > 0) {
      real u2, v2;
      if (fbet().NCoeffsInv() <= fomg().NCoeffsInv())
          solve2(-_delta, sig2, fomg(), fbet(), gomg(), gbet(), u2, v2);
        else
          solve2( _delta, sig2, fbet(), fomg(), gbet(), gomg(), v2, u2);
      omg2 = _E * fomg().rev(u2);
      omg2a = AuxAngle::radians(omg2);
      AuxAngle psi2 = AuxAngle::radians(fbet().rev(v2));
      real sgam = sqrt(_gam/_t.k2);
      bet2a = AuxAngle(_bet0 * (_flip * sqrt((_t.k2 - _gam)/_t.k2)) * psi2.y(),
                       _bet0 * hypot(psi2.x(), sgam * psi2.y()) );
      bet2 = bet2a.radians();
      alp2a = AuxAngle(_alp0 * _E * sqrt(_gam + _t.kp2 * Math::sq(cos(omg2))),
                       _alp0 * (_flip * sqrt(_t.k2 - _gam)) * psi2.x() );
      alp2 = alp2a.radians();
    } else if (_gam < 0) {
      real u2, v2;
      if (fomg().NCoeffsInv() <= fbet().NCoeffsInv())
        solve2( _delta, sig2, fbet(), fomg(), gbet(), gomg(), u2, v2);
      else
        solve2(-_delta, sig2, fomg(), fbet(), gomg(), gbet(), v2, u2);
      bet2 = _N * fbet().rev(u2);
      bet2a = AuxAngle::radians(bet2);
      AuxAngle psi2 = AuxAngle::radians(fomg().rev(v2));
      real sgam = sqrt(-_gam/_t.kp2);
      omg2a = AuxAngle(_omg0 * (_flip * sqrt((_t.kp2 + _gam)/_t.kp2)) *
                       psi2.y(),
                       _omg0 * hypot(psi2.x(), sgam * psi2.y()));
      omg2 = omg2a.radians();
      alp2a = AuxAngle(_alp0 * (_N * _flip * sqrt(_t.kp2 + _gam)) * psi2.x(),
                       _alp0 * sqrt(-_gam + _t.k2 * Math::sq(cos(bet2))));
      alp2 = alp2a.radians();
    } else {                    // gam == 0
      pair<real, real> sig2n = remx(sig2, 2*_s0);  // reduce to [-s0, s0)
      real u2, v2,
        deltax = max(-_deltamax,
                     min(_deltamax, _delta + sig2n.second * _deltashift));
      solve2u( deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2);
      /*
      if (fbet().NCoeffsInv() <= fomg().NCoeffsInv())
        solve2u( deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2);
      else
        solve2u(-deltax, sig2n.first, fomg(), fbet(), gomg(), gbet(), v2, u2);
      */
      bet2 = fbet().rev(u2); omg2 = fomg().rev(v2);
      bet2a = AuxAngle::lam(u2); omg2a = AuxAngle::lam(v2);
      if (_umbalt) {
        real Ex = _E * (1 - 2*fmod(sig2n.second, real(2)));
        omg2 = (_omg0 < 0 ? Math::pi() : 0) + Ex * omg2;
        bet2 = (_bet0 < 0 ? Math::pi() : 0) +
          _N * (bet2 + sig2n.second * Math::pi());
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2 = atan2( _alp0 * (_N * sqrt(_t.kp2)) * (Ex / cosh(v2)),
                      _alp0 * sqrt(_t.k2) / cosh(u2) );
      } else {
        real Nx = _N * (1 - 2*fmod(sig2n.second, real(2)));
        omg2 = (_omg0 < 0 ? Math::pi() : 0) +
          _E * (omg2 + sig2n.second * Math::pi());
        bet2 = (_bet0 < 0 ? Math::pi() : 0) + Nx * bet2;
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2 = atan2( _alp0 * (_E * sqrt(_t.kp2)) / cosh(v2),
                      _alp0 * sqrt(_t.k2) * (Nx / cosh(u2)) );
      }
    }
    swap(omg2a.x(), omg2a.y());
    omg2a.x() *= -1;
    omg2 = omg2 / Math::degree() + 90;
    bet2 = bet2 / Math::degree();
    alp2 = alp2 / Math::degree();
  }

  void TriaxialLine::solve2(real f0, real g0,
                            const geod_fun& fx, const geod_fun& fy,
                            const dist_fun& gx, const dist_fun& gy,
                            real& x, real& y) {
    x = fx(0) - fy(0) - f0; y = gx(0) + gy(0) - g0;
  }
  void TriaxialLine::solve2u(real f0, real g0,
                             const geod_fun& fx, const geod_fun& fy,
                             const dist_fun& gx, const dist_fun& gy,
                             real& x, real& y) {
    x = fx(0) - fy(0) - f0; y = gx(0) + gy(0) - g0;
  }

  TriaxialLineF::TriaxialLineF(real k2, real kp2, real e2, real gam,
                               real epspow, real nmaxmult)
    : _k2(k2)
    , _kp2(kp2)
    , _e2(e2)
    , _gam(gam)
    , _fbet(_k2 , _kp2,  _e2, -_gam, epspow, nmaxmult)
    , _fomg(_kp2, _k2 , -_e2,  _gam, epspow, nmaxmult)
  {}

  TriaxialLineG::TriaxialLineG(real k2, real kp2, real e2, real gam)
    : _k2(k2)
    , _kp2(kp2)
    , _e2(e2)
    , _gam(gam)
    , _gbet(_k2 , _kp2,  _e2, -_gam)
    , _gomg(_kp2, _k2 , -_e2,  _gam)
  {}

} // namespace GeographicLib

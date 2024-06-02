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
                             AuxAngle bet1, AuxAngle omg1, AuxAngle alp1)
    : _t(t)
    , _gm(_t, bet1, omg1, alp1)
    // This normalizes bet1, omg1, and alp1 and tweaks their values
      /*
    , _gamma(_t.gamma(bet1, omg1, alp1))
    , _rtgam(sqrt(fabs(_gamma)))
    , _gamp(_t.gammap(_gamma, bet1, omg1, alp1))
    , _rtgamp(sqrt(_gamp))
      */
    , _bet1(bet1)
      // omg1 - 90
    , _omg1(AuxAngle(-omg1.x(), omg1.y()))
    , _alp1(alp1)
    , _ibet(0)
      // omg1 - 90 crosses omg1 = -180 if cos(omg1) < 0 and sin(omg1) < 0
    , _iomg(signbit(omg1.x()) && signbit(omg1.y()) ? -1 : 0)
    , _ialp(0)
    , _f(_t, _gm, 0.5, 1.5) // 0.333, 4
    , _umbalt(false)
    , _distinit(false)
  {
    const real eps = numeric_limits<real>::epsilon();
    if (_bet1.y() == 0 && _alp1.x() == 0)
      _alp1.x() = - Math::sq(eps);
    _eE = signbit(_alp1.y()) ? -1 : 1;
    _nN = signbit(_alp1.x()) ? -1 : 1;
    if (_gm.gam > 0) {
      _flip = signbit(bet1.x()) ? -1 : 1;
      _bet0 = _flip;
      _alp0 = 1;
      _psi1 = AuxAngle(_t.k * bet1.y(),
                       _flip * _alp1.x() *
                       hypot(_t.k * _bet1.x(), _t.kp * _omg1.x()),
                       false);
      _v0 = fbet().fwd(_psi1.radians());
      _u0 = fomg().fwd(_eE * _omg1.radians());
      _delta = fbet()(_v0) - fomg()(_u0);
      if (0) {
        cout << "TX " << fbet().txp() << " " << fomg().txp() << "\n";
        cout << "QQ " << _psi1.radians() << " "
             << _eE * _omg1.radians() << " "
             << _v0 << " " << _u0 << "\n";
        cout << "PP " << fbet()(_v0) << " " << fomg()(_u0) << " " << _delta << "\n";
      }
    } else if (_gm.gam < 0) {
      _flip = signbit(_omg1.x()) ? -1 : 1;
      _omg0 = _flip;
      _alp0 = _nN;
      _psi1 = AuxAngle(_t.kp * _omg1.y(),
                       _flip * alp1.y() *
                       hypot(_t.k * _bet1.x(), _t.kp * _omg1.x()),
                       false);
      _v0 = fomg().fwd(_psi1.radians());
      _u0 = fbet().fwd(_nN * _bet1.radians());
      _delta = fbet()(_u0) - fomg()(_v0);
    } else if (_gm.gam == 0) {
      _alp0 = _umbalt && _nN < 0 ? _eE : 0;
                                // N.B. factor of k*kp omitted
      _df = fbet().Max() - fomg().Max();
      _deltashift = 2*_df - log(_t.k2/_t.kp2);
      _deltamax = -2*log(eps/2); // = asinh(2/eps^2)
      //      _deltamax *= 0.75;
      if (abs(_bet1.x()) < 8*eps && abs(_omg1.x()) < 8*eps) {
        //        _bet0 = (int(round(_bet1.y())) + _nN) / 2 ? -1 : 1;
        //        _omg0 = (int(round(_omg1.y())) + _eE) / 2 ? -1 : 1;
        _bet0 = (int(round(_bet1.y())) + _nN) / 2; // -1, 0, or +1
        _omg0 = (int(round(_omg1.y())) + _eE) / 2;
        _delta = _deltashift/2 - log(abs(_alp1.tan()));
      } else {
        _bet0 = signbit(_bet1.x()) ? (signbit(_bet1.y()) ? -1 : 1) : 0;
        _omg0 = signbit(_omg1.x()) ? (signbit(_omg1.y()) ? -1 : 1) : 0;
        /*
        _delta = _nN * fbet()(fbet().fwd(atan2(_bet1.y() * _bet0,
                                              _bet1.x() * _bet0))) -
          _eE * fomg()(fomg().fwd(atan2(_omg1.y() * _omg0,
                                       _omg1.x() * _omg0)));
        */
        _delta = _nN * fbet()(asinh(_bet1.tan())) -
          _eE * fomg()(asinh(_omg1.tan()));

      }
      /*
      cout << "YYY " << _gm.gam << " "
           << _df << " " << _deltashift << " "
           << _deltamax << " " << _delta << " "
           << _bet1.degrees() << " " << _omg1.degrees() << " "
           << _bet0 << " " << _omg0 << "\n";
      */
    } else {
      // gamma = NaN
    }

    distinit();
    /*
    cout << _df  << " " << _deltashift << " "
         << _deltamax << " " << _delta << " "
         << _sig1 << " " << _s0 << "\n";
    */
    //    cout << _eE << " " << _nN << " " << _omg0 << " " << _bet0 << " " << _alp0 << "\n";
    bool debug = false;

    if (debug) {
      distinit();
      real u = 1;
      if (1) {
        /*
          real du = sqrt(sqrt(numeric_limits<real>::epsilon()));
          cout << "gam " << _gm.gam << "\n";
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
        cout << "cnts " << _gm.gam << " "
             << fbet().NCoeffs() << " " << fbet().NCoeffsInv() << " "
             << fomg().NCoeffs() << " " << fomg().NCoeffsInv() << " "
             << gbet().NCoeffs() << " " << gomg().NCoeffs() << "\n";
        cout << "invs "
             << fbet().inv(fbet()(u)) - u << " "
             << fomg().inv(fomg()(u)) - u << "\n";
      }
      if (0) {
        cout << "gam " << _gm.gam << "\n";
        cout << "fbet " << " " << fbet().txp() << gbet().txp() << " "
             << gbet().deriv(u) / fbet().deriv(u) - gbet().gfderiv(u)
             << "\n";
        cout << "fomg " << " " << fomg().txp() << gomg().txp() << " "
             << gomg().deriv(u) / fomg().deriv(u) - gomg().gfderiv(u)
             << "\n";
      }
      if (0) {
        cout << "gam " << _gm.gam << "\n";
        geod_fun fbeta(_t.k2, _t.kp2, _t.e2, -_gm.gam, false, 0.333, 4);
        geod_fun fomga(_t.kp2, _t.k2, -_t.e2, _gm.gam, false, 0.333, 4);
        dist_fun gbeta(_t.k2, _t.kp2, _t.e2, -_gm.gam, false);
        dist_fun gomga(_t.kp2, _t.k2, -_t.e2, _gm.gam, false);
        geod_fun fbetb(_t.k2, _t.kp2, _t.e2, -_gm.gam, true, 0.333, 4);
        geod_fun fomgb(_t.kp2, _t.k2, -_t.e2, _gm.gam, true, 0.333, 4);
        dist_fun gbetb(_t.k2, _t.kp2, _t.e2, -_gm.gam, true);
        dist_fun gomgb(_t.kp2, _t.k2, -_t.e2, _gm.gam, true);
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
  TriaxialLine::TriaxialLine(const Triaxial& t, real bet1, real omg1, real alp1)
      : TriaxialLine(t,
                     AuxAngle::degrees(bet1),
                     AuxAngle::degrees(omg1),
                     AuxAngle::degrees(alp1))
    {
      _ibet = int(round((bet1 - _bet1.degrees()) / Math::td));
      _iomg += int(round((omg1 - AuxAngle::degrees(omg1).degrees()) /
                         Math::td));
      _ialp = int(round((alp1 - _alp1.degrees()) / Math::td));
    }

  void TriaxialLine::distinit() {
    if (_distinit) return;
    _g = TriaxialLineG(_t, _gm);
    if (_gm.gam > 0) {
      _sig1 = gbet()(_v0) + gomg()(_u0);
      _s0 = 0;
    } else if (_gm.gam < 0) {
      _sig1 = gbet()(_u0) + gomg()(_v0);
      _s0 = 0;
    } else if (_gm.gam == 0) {
      /*
      cout << "DEBUG "
           << atan2(_bet1.y() * _bet0, _bet1.x() * _bet0) << " "
           << fbet().fwd(atan2(_bet1.y() * _bet0,
                               _bet1.x() * _bet0)) << " "
           << gbet()(fbet().fwd(atan2(_bet1.y() * _bet0,
                                      _bet1.x() * _bet0))) << "\n"
           << "DEBUG "
           << atan2(_omg1.y() * _omg0, _omg1.x() * _omg0) << " "
           << fomg().fwd(atan2(_omg1.y() * _omg0,
                               _omg1.x() * _omg0)) << " "
           << gomg()(fomg().fwd(atan2(_omg1.y() * _omg0,
                                      _omg1.x() * _omg0))) << "\n";
      */
      _sig1 = _nN * gbet()(fbet().fwd(atan2(_bet1.y() * ( _bet0 ? -1 : 1),
                                            _bet1.x() * ( _bet0 ? -1 : 1)))) +
        _eE * gomg()(fomg().fwd(atan2(_omg1.y() * ( _omg0 ? -1 : 1),
                                      _omg1.x() * ( _omg0 ? -1 : 1))));
      /*
      cout << "QQQ " << _bet1.radians() << " " << _bet0 << " " << _nN << " "
           << gbet()(fbet().fwd(atan2(_bet1.y() * _bet0,
                                      _bet1.x() * _bet0))) << " "
           << _eE  << " "
           << gomg()(fomg().fwd(atan2(_omg1.y() * _omg0,
                                      _omg1.x() * _omg0))) << "\n";
      */
      _s0 = gbet().Max() + gomg().Max();
      // cout << "QQ " << " " << _sig1 << " " << _s0 << "\n";
      /*
        _sig1 = _nN * gbet()(asinh(_bet1.tan())) -
        _eE * gomg()(asinh(_omg1.tan()));
      */
    } else {
      // gamma = NaN
    }


    _distinit = true;
  }

  void TriaxialLine::Position(real s12,
                              AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a,
                              int* ibet2, int* iomg2, int* ialp2,
                              int* countn, int* countb)
  const {
    // Compute points at distance s12
    real sig2 = _sig1 + s12/_t.b;
    // cout << "SIG0 " << _sig1 << " " << s12/_t.b << " " << sig2 << "\n";
    real bet2, omg2, alp2;
    int Ex, Nx;
    if (_gm.gam > 0) {
      if (0) {
        cout << "AA " << _delta << " " << _sig1 << " " << _s0 << "\n";
      }
      real u2, v2;
      if (fbet().NCoeffsInv() <= fomg().NCoeffsInv())
        solve2(-_delta, sig2, fomg(), fbet(), gomg(), gbet(), u2, v2,
               countn, countb);
      else
        solve2( _delta, sig2, fbet(), fomg(), gbet(), gomg(), v2, u2,
                countn, countb);
      omg2 = _eE * fomg().rev(u2);
      omg2a = AuxAngle::radians(omg2);
      AuxAngle psi2 = AuxAngle::radians(fbet().rev(v2));
      bet2a = AuxAngle(_bet0 * _flip * _gm.nu.x() * psi2.y(),
                       _bet0 * hypot(psi2.x(), _gm.nu.y() * psi2.y()),
                       false);
      bet2 = bet2a.radians();
      alp2a = AuxAngle(_alp0 * _eE *
                       hypot(_t.k * _gm.nu.y(), _t.kp * omg2a.x()),
                       _alp0 * _flip * _t.k * _gm.nu.x() * psi2.x(),
                       false);
      alp2 = alp2a.radians();
    } else if (_gm.gam < 0) {
      real u2, v2;
      if (fomg().NCoeffsInv() <= fbet().NCoeffsInv())
        solve2( _delta, sig2, fbet(), fomg(), gbet(), gomg(), u2, v2,
                countn, countb);
      else
        solve2(-_delta, sig2, fomg(), fbet(), gomg(), gbet(), v2, u2,
               countn, countb);
      bet2 = _nN * fbet().rev(u2);
      bet2a = AuxAngle::radians(bet2);
      AuxAngle psi2 = AuxAngle::radians(fomg().rev(v2));
      omg2a = AuxAngle(_omg0 * _flip * _gm.nu.x() * psi2.y(),
                       _omg0 * hypot(psi2.x(), _gm.nu.y() * psi2.y()),
                       false);
      omg2 = omg2a.radians();
      alp2a = AuxAngle(_alp0 * _nN * _flip * _t.kp * _gm.nu.x() * psi2.x(),
                       _alp0 * hypot(_t.kp * _gm.nu.y(), _t.k * bet2a.x()),
                       false);
      alp2 = alp2a.radians();
    } else if (_gm.gam == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_s0);  // reduce to [-s0, s0)
      real u2, v2,
        deltax = max(-_deltamax,
                     min(_deltamax, _delta + sig2n.second * _deltashift));
      solve2u(deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2,
              countn, countb);
      /*
      cout << "AAA " << _sig1 << " " << sig2 << " " << _s0 << " " << deltax << " " << sig2n.first << " "
           << u2 << " " << v2 << "\n";
      */
      /*
      if (fbet().NCoeffsInv() <= fomg().NCoeffsInv())
        solve2u( deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2,
                countn, countb);
      else
        solve2u(-deltax, sig2n.first, fomg(), fbet(), gomg(), gbet(), v2, u2,
                countn, countb);
      */
      bet2 = fbet().rev(u2); omg2 = fomg().rev(v2);
      bet2a = AuxAngle::lam(u2); omg2a = AuxAngle::lam(v2);
      // cout << "UU0 " << bet2 / Math::degree() << " " << bet2a.degrees() << "\n";
      // cout << "UU0 " << omg2 / Math::degree() << " " << omg2a.degrees() << "\n";
      // cout << "SIG2N " << sig2 << " " << _s0 << " " << sig2n.first << " " << sig2n.second << "\n";
      int parity = fmod(sig2n.second, real(2)) ? -1 : 1;
      if (0)
        cout << "DD " << s12 << " " << u2 << " " << v2 << " " << sig2n.first << " " << sig2n.second << " " << bet2a.radians() << " " << omg2a.radians() << "\n";
      if (_umbalt) {
        Ex = _eE * parity;
        omg2a.y() *= Ex * ( _omg0 ? -1 : 1);
        omg2a.x() *= ( _omg0 ? -1 : 1);
        omg2 = _omg0 * Math::pi() + Ex * omg2;
        bet2a.y() *= _nN * parity * ( _bet0 ? -1 : 1);
        bet2a.x() *= parity * ( _bet0 ? -1 : 1);
        bet2 = //(_bet0 < 0 ? Math::pi() : 0) +
          _bet0 * Math::pi() +
          _nN * (bet2 + sig2n.second * Math::pi());
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = AuxAngle((_alp0 ? -1 : 1) * (_nN * _t.kp) * (Ex / cosh(v2)),
                         (_alp0 ? -1 : 1) * _t.k / cosh(u2),
                         false);
        alp2 = alp2a.radians();
        /*
        cout << "DATU " << _alp0 << " " << _nN << " " << _eE << " "
             << ( _omg0 ? -1 : 1) << " "
             << Ex << " " << u2 << " " << v2 << " " << alp2 << "\n";
        */
      } else {
        if (0)
          cout << "BB " << bet2a.radians() << " " << bet2 << " "
         << omg2a.radians() << " " << omg2 << " " << sig2n.second << " " << parity << "\n";
        Nx = _nN * parity;
        omg2a.y() *= _eE * parity * ( _omg0 ? -1 : 1);
        omg2a.x() *= parity * ( _omg0 ? -1 : 1);
        omg2 = _omg0 * Math::pi() +
          _eE * (omg2 + sig2n.second * Math::pi());
        bet2a.y() *= Nx * ( _bet0 ? -1 : 1);
        bet2a.x() *= ( _bet0 ? -1 : 1);
        bet2 = // (_bet0 < 0 ? Math::pi() : 0) +
          _bet0 * Math::pi() +
          Nx * bet2;
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        // _alp0 = 0
        alp2a = AuxAngle((_eE * _t.kp) / cosh(v2), _t.k * (Nx / cosh(u2)),
                         false);
        alp2 = alp2a.radians();
      }
    } else {
      // gamma = NaN
    }
    if (0)
    cout << "AA " << bet2a.radians() << " " << bet2 << " "
         << omg2a.radians() << " " << omg2 << " "
         << alp2a.radians() << " " << alp2 << "\n";

    swap(omg2a.x(), omg2a.y());
    omg2a.x() *= -1;
    omg2 = omg2 / Math::degree() + 90;
    bet2 = bet2 / Math::degree();
    alp2 = alp2 / Math::degree();
    //    cout << "UUU " << _bet0 << " " <<  _ibet << " " << bet2 << " " << bet2a.degrees() << "\n";
    if (ibet2)
      *ibet2 = _ibet + int(round((bet2 - bet2a.degrees()) / Math::td)) -
        (_gm.gam > 0 && _flip < 0 ?
         (signbit(bet2a.y()) ? -1 : 1) - (signbit(_bet1.y()) ? -1 : 1) : 0)/2;
    if (iomg2)
      *iomg2 = _iomg + int(round((omg2 - omg2a.degrees()) / Math::td)) +
        (_gm.gam < 0 && _flip < 0 ?
         (signbit(omg2a.x()) ? -1 : 1) + (signbit(_omg1.y()) ? -1 : 1) : 0)/2;
    if (ialp2)
      *ialp2 = _ialp + int(round((alp2 - alp2a.degrees()) / Math::td)) -
        (_gm.gam < 0 && _nN < 0 ?
         (signbit(alp2a.y()) ? -1 : 1) - _eE : 0)/2 +
        (_gm.gam == 0 && _umbalt ? (_alp0 == _nN * Ex ? _alp0 : 0) : 0);
    //    cout << "AAA " << *ibet2 << " " << bet2a.degrees() << " " << _gm.gam << "\n";
    //    cout << "BBB " << _flip << " " <<       _bet0 << " " <<       _alp0 << " " << _psi1.degrees() << " " <<       _v0<< " " <<       _u0<< " " <<       _delta<< "\n";
  }
  void TriaxialLine::Position(real s12, real& bet2, real& omg2, real& alp2,
                              bool unroll,
                              int* countn, int* countb) const {
    AuxAngle bet2a, omg2a, alp2a;
    int ibet2 = 0, iomg2 = 0, ialp2 = 0;
    Position(s12, bet2a, omg2a, alp2a, &ibet2, &iomg2, &ialp2,
             countn, countb);
    if (unroll) {
      bet2 = bet2a.degrees() + ibet2 * Math::td;
      omg2 = omg2a.degrees() + iomg2 * Math::td;
      alp2 = alp2a.degrees() + ialp2 * Math::td;
    } else {
      Triaxial::AngNorm(bet2a, omg2a, alp2a);
      bet2 = bet2a.degrees();
      omg2 = omg2a.degrees();
      alp2 = alp2a.degrees();
    }
  }

  void TriaxialLine::solve2(real f0, real g0,
                            const geod_fun& fx, const geod_fun& fy,
                            const dist_fun& gx, const dist_fun& gy,
                            real& x, real& y,
                            int* countn, int* countb) {
    // Return x and y, s.t.
    // fx(x) - fy(y) = f0
    // gx(x) + gy(y) = g0
    // We assume that fy.inv() is available
    // fx(x) = fxs*x +/- fxm,
    // fy(y) = fys*y +/- fym,
    // gx(x) = gxs*x +/- gxm,
    // gy(y) = gys*y +/- gym;
    real fxm = fx.Max(), fym = fy.Max(), gxm = gx.Max(), gym = gy.Max(),
      fxs = fx.Slope(), fys = fy.Slope(), gxs = gx.Slope(), gys = gy.Slope(),
      // solve
      //   x = (  fys*g0 + gys*f0 ) / den +/- Dx
      //   y = (- gxs*f0 + fxs*g0 ) / den +/- Dy
      // where
      den = fxs * gys + fys * gxs,
      qf = fxm + fym, qg = gxm + gym,
      Dx = (qf * gys + qg * fys) / fabs(den);
    // Dy = (qf * gxs + qg * fxs) / fabs(den);
    real x0 = (fys * g0 + gys * f0) / den, // Initial guess
      xp = x0 + Dx, xm = x0 - Dx;
    newt2(f0, g0, fx, fy, gx, gy, x0, xm, xp,
          fx.HalfPeriod(), fx.HalfPeriod() * fxs,
          x, y, countn, countb);
  }

  void TriaxialLine::solve2u(real d0, real s0,
                             const geod_fun& fbet, const geod_fun& fomg,
                             const dist_fun& gbet, const dist_fun& gomg,
                             real& u, real& v,
                             int* countn, int* countb) {
    // Return u and v, s.t.
    // fbet(u) - fomg(v) = d0
    // gbet(u) + gomg(v) = s0
    // specialized for umbilics
    //
    // fbet, fomg, gbet, gomg are increasing functions defined in [-1, 1]*pi2
    // Assume fbet(0) = fomg(0) = gbet(0) = gomg(0) = 0
    real pi2 = -2*log(numeric_limits<real>::epsilon()/2),
      sbet = gbet.Max(), somg = gomg.Max(), stot = sbet + somg,
      dbet = fbet.Max(), domg = fomg.Max(), del  = dbet - domg;
    pi2 *= 0.75;
    if (abs(s0) >= stot) {
      // close to umbilic points we have
      // fbet(u) = u -/+ dbet
      // fomg(v) = v -/+ domg
      // fbet(u) - fomg(v) = d0
      // constrain u+v = +/- 2 * pi2
      // u = +/- pi2 + (d0 +/- (dbet - domg))/2
      // v = +/- pi2 - (d0 +/- (dbet - domg))/2
      real t0 = copysign(pi2, s0),
        t1 = (d0 + (1 - 2 * signbit(s0)) * del) / 2;
      u = t0 + t1; v = t0 - t1;
    } else if ((1 - 2 * signbit(d0)) * s0 < sbet - somg) {
      // Use u as independent variable if
      //   d0 < 0 ? (s0 > -sbet + somg) :
      //            (s0 <  sbet - somg)
      // or
      //   sign(d0) * s0 < sbet - somg
      newt2(d0, s0, fbet, fomg, gbet, gomg, 0, -pi2, pi2, 2, 2,
            u, v, countn, countb);
    } else {
      // Otherwise, use v is the independent variable
      newt2(-d0, s0, fomg, fbet, gomg, gbet, 0, -pi2, pi2, 2, 2,
            v, u, countn, countb);
    }
  }

  void TriaxialLine::newt2(real f0, real g0,
                           const geod_fun& fx, const geod_fun& fy,
                           const dist_fun& gx, const dist_fun& gy,
                           real x0, real xa, real xb,
                           real xscale, real zscale,
                           real& x, real& y,
                           int* countn, int* countb) {
    x = Trigfun::root(
                      [&fx, &fy, &gx, &gy, f0, g0]
                      (real x) -> pair<real, real>
                      { real y = fy.inv(fx(x) - f0);
                        return pair<real, real>(gx(x) + gy(y),
                                                fx.deriv(x) *
                                                (gx.gfderiv(x) +
                                                 gy.gfderiv(y))); },
                      g0, x0, xa, xb, xscale, zscale, 1,
                      countn, countb);
    y = fy.inv(fx(x) - f0);
  }

  void TriaxialLine::Hybrid(const AuxAngle& bet2, int dir,
                            AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a,
                            real& s12)
    const {
    //    cout << "GG " << _gm.gam << "\n";
    if (_gm.gam > 0) {
      real spsi = _t.k * bet2.y(),
        // In evaluating equivalent expressions, choose the one
        // with minimum cancelation
        cpsi = 0 + _t.k
        * sqrt(fmax(0,
                    _gm.nu.y() < fabs(bet2.y()) ?
                    (bet2.x() - _gm.nu.y()) * (bet2.x() + _gm.nu.y()) :
                    (_gm.nu.x() - bet2.y()) * (_gm.nu.x() + bet2.y())));
      AuxAngle psi2 = AuxAngle(spsi, dir * cpsi, true),
        psi12 = psi2 - _psi1;
      psi12.y() = fmax(real(0), psi12.y()) + 0; // convert -180deg to 180deg
      real tau12 = psi12.radians(),
        psi2r = _psi1.radians() + tau12,
        v2= fbet().fwd(psi2r),
        u2 = fomg().inv(fbet()(v2) - _delta),
        omg2 = _eE * fomg().rev(u2);
      // cout << "PSI2 " << psi2.degrees() << " " << psi12.degrees() << " "
      //   << 1/psi12.y()  << " " << tau12 / Math::degree() << "\n";
      omg2a = AuxAngle::radians(omg2);
      bet2a = AuxAngle(_bet0 * _flip * _gm.nu.x() * psi2.y(),
                       _bet0 * hypot(psi2.x(), _gm.nu.y() * psi2.y()),
                       true);
      alp2a = AuxAngle(_alp0 * _eE
                       * hypot(_t.k * _gm.nu.y(), _t.kp * omg2a.x()),
                       _alp0 * _flip * _t.k * _gm.nu.x() * psi2.x(),
                       true);
      real sig2 = gbet()(v2) + gomg()(u2);
      s12 = (sig2 - _sig1)  * _t.b;
    } else if (_gm.gam <= 0) {
      bet2a = AuxAngle(bet2.y(), bet2.x() * dir * _nN, false);
      AuxAngle bet12 = bet2a - _bet1;
      bet12.y() = fmax(real(0), bet12.y()) + 0; // convert -180deg to 180deg
      real tau12 =  bet12.radians(),
        bet2r = _bet1.radians() + _nN * tau12;
      //      cout << "TAU12 " << tau12/Math::degree() << " "
      //           << bet2r/Math::degree() << "\n";
      if (_gm.gam < 0) {
        real
          u2 = fbet().fwd(_nN * bet2r),
          v2 = fomg().inv(fbet()(u2) - _delta),
          psi2r = fomg().rev(v2);
        AuxAngle psi2 = AuxAngle::radians(psi2r);
        omg2a = AuxAngle(_omg0 * _flip * _gm.nu.x() * psi2.y(),
                         _omg0 * hypot(psi2.x(), _gm.nu.y() * psi2.y()),
                         true);
        alp2a = AuxAngle(_alp0 * _nN * _flip * _t.kp * _gm.nu.x() * psi2.x(),
                         _alp0 * hypot(_t.kp * _gm.nu.y(), _t.k * bet2.x()),
                         true);
        real sig2 = gbet()(u2) + gomg()(v2);
        s12 = (sig2 - _sig1)  * _t.b;
      } else {                  // _gm.gam == 0
        // Could simplify this.  bet2 is in [-270,90]
        pair<real, real> bet2n =
          remx(_nN * bet2r, Math::pi());  // reduce to [-pi/2, pi/2)
        int parity = fmod(bet2n.second, real(2)) ? -1 :  1;
        int alp0 = _nN < 0 ? _eE : 0;
        real deltax = fmax(-_deltamax,
                           min(_deltamax, _delta + bet2n.second * _deltashift)),
          u2 = fbet().fwd(bet2n.first),
          v2 = fomg().inv(fbet()(u2) - deltax),
          omg2 = _eE * parity * fomg().rev(v2);
        omg2a = AuxAngle::radians(omg2);
        alp2a = AuxAngle((alp0 ? -1 : 1) * _nN * _t.kp * _eE * parity  /
                         cosh(v2),
                         (alp0 ? -1 : 1) * _t.k / cosh(u2),
                         true);
        if (0)
        cout << "AAA "
             << alp2a.degrees() << " " << alp0 << " "
             << (_alp0 ? -1 : 1) * _nN * _t.kp * _eE * parity  / cosh(v2) << " "
             << (_alp0 ? -1 : 1) * _t.k / cosh(u2) << "\n";
        real sig2 = gbet()(u2) + gomg()(v2) + bet2n.second * 2 * _s0;
        s12 = (sig2 - _sig1)  * _t.b;
      }
    } else {
      // gamma = NaN
    }
    omg2a += AuxAngle::cardinal(1U);
    if (_gm.gam <= 0 && dir * _nN < 0) {
      //      omg2a.y() *= -1;
      //      alp2a += AuxAngle::cardinal(2U);
    }
  }

  TriaxialLineF::TriaxialLineF(const Triaxial& t, Triaxial::gamblk gam,
                               real epspow, real nmaxmult)
    : _t(t)
    , _gm(gam)
    , _fbet(_t.k2 , _t.kp2,  _t.e2, -_gm.gam, epspow, nmaxmult)
    , _fomg(_t.kp2, _t.k2 , -_t.e2,  _gm.gam, epspow, nmaxmult)
    {
      _fbet.ComputeInverse();
      _fomg.ComputeInverse();
    }

  TriaxialLineG::TriaxialLineG(const Triaxial& t, Triaxial::gamblk gam)
    : _t(t)
    , _gm(gam)
    , _gbet(_t.k2 , _t.kp2,  _t.e2, -_gm.gam)
    , _gomg(_t.kp2, _t.k2 , -_t.e2,  _gm.gam)
  {}

} // namespace GeographicLib

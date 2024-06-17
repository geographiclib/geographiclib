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

  TriaxialLine::TriaxialLine(TriaxialLineF f, TriaxialLineF::ics fic,
                             TriaxialLineG g, TriaxialLineG::ics gic)
    : _t(f.t())
    , _f(f)
    , _fic(fic)
    , _g(g)
    , _gic(gic)
  {}

  TriaxialLine::TriaxialLine(const Triaxial& t,
                             AuxAngle bet1, AuxAngle omg1, AuxAngle alp1) {
    bet1.normalize().round();
    omg1.normalize().round();
    alp1.normalize().round();
    _t = t;
    Triaxial::gamblk gam(t, bet1, omg1, alp1);
    _f = TriaxialLineF(t, gam, 0.5, 1.5);
    _fic = TriaxialLineF::ics(_f, bet1, omg1, alp1);
    _g = TriaxialLineG(t, gam);
    _gic = TriaxialLineG::ics(_g, _fic);
  }

  TriaxialLine::TriaxialLine(const Triaxial& t, real bet1, real omg1, real alp1)
      : TriaxialLine(t,
                     AuxAngle::degrees(bet1),
                     AuxAngle::degrees(omg1),
                     AuxAngle::degrees(alp1))
    {
      _fic.ibet = int(round((bet1 - _fic.bet1.degrees()) / Math::td));
      _fic.iomg += int(round((omg1 - AuxAngle::degrees(omg1).degrees()) /
                         Math::td));
      _fic.ialp = int(round((alp1 - _fic.alp1.degrees()) / Math::td));
    }

  void TriaxialLine::diag() const {
    cout << _f.fbet().HalfPeriod() << " " << _f.fbet().Slope() << " "
         << _f.fbet().Max() << " "
         << _f.fbet().NCoeffs() << "\n";
  }
  void TriaxialLine::pos1(AuxAngle& bet1, AuxAngle& omg1, AuxAngle& alp1,
                          int* ibet1, int* iomg1, int* ialp1) const {
    bet1 = _fic.bet1; omg1 = _fic.omg1 + AuxAngle::cardinal(1U);
    alp1 = _fic.alp1;
    if (ibet1) *ibet1 = _fic.ibet;
    if (iomg1) *iomg1 = _fic.iomg + (omg1.quadrant() == 2U ? 1 : 0);
    if (ialp1) *ialp1 = _fic.ialp;
  }

  void TriaxialLine::pos1(real& bet1, real& omg1, real& alp1,
                          bool unroll) const {
    AuxAngle bet1a, omg1a, alp1a;
    int ibet1 = 0, iomg1 = 0, ialp1 = 0;
    pos1(bet1a, omg1a, alp1a, &ibet1, &iomg1, &ialp1);
    if (unroll) {
      bet1 = bet1a.degrees() + ibet1 * Math::td;
      omg1 = omg1a.degrees() + iomg1 * Math::td;
      alp1 = alp1a.degrees() + ialp1 * Math::td;
    } else {
      (void) Triaxial::AngNorm(bet1a, omg1a, alp1a);
      bet1 = bet1a.degrees();
      omg1 = omg1a.degrees();
      alp1 = alp1a.degrees();
    }
  }

  void TriaxialLine::Position(real s12,
                              AuxAngle& bet2a, AuxAngle& omg2a,
                              AuxAngle& alp2a,
                              int* ibet2, int* iomg2, int* ialp2,
                              int* countn, int* countb)
  const {
    if (0) {
      cout << "THERE0 " << _f.gamma() << " " << _gic.s0 << "\n";
      cout << "AA " << _fic.delta << " " << _gic.sig1 << " " << _gic.s0 << "\n";
    }
    // Compute points at distance s12
    real sig2 = _gic.sig1 + s12/_t.b;
    // cout << "SIG0 " << _gic.sig1 << " " << s12/_t.b << " " << sig2 << "\n";
    // cout << "SIG1 " << s12 << " " << _t.b << "\n";
    real bet2, omg2, alp2;
    int Ex, Nx;
    if (_f.gamma() > 0) {
      real u2, v2;
      if (fbet().NCoeffsInv() <= fomg().NCoeffsInv())
        solve2(-_fic.delta, sig2, fomg(), fbet(), gomg(), gbet(), u2, v2,
               countn, countb);
      else
        solve2( _fic.delta, sig2, fbet(), fomg(), gbet(), gomg(), v2, u2,
                countn, countb);
      omg2 = _fic.eE * fomg().rev(u2);
      omg2a = AuxAngle::radians(omg2);
      AuxAngle psi2 = AuxAngle::radians(fbet().rev(v2));
      bet2a = AuxAngle(_fic.bet0 * _fic.flip * _f.gm().nup * psi2.y(),
                       _fic.bet0 * hypot(psi2.x(), _f.gm().nu * psi2.y()),
                       false);
      bet2 = bet2a.radians();
      alp2a = AuxAngle(_fic.alp0 * _fic.eE *
                       hypot(_t.k * _f.gm().nu, _t.kp * omg2a.x()),
                       _fic.alp0 * _fic.flip * _t.k * _f.gm().nup * psi2.x(),
                       false);
      alp2 = alp2a.radians();
      if (0) {
      cout << "POSQ " << omg2/Math::degree() << " " << fomg().rev(u2) << " " << fbet().rev(v2) << "\n";
      cout << "POSP " << u2 << " " << v2 << " " << _fic.delta << "\n";
      }
    } else if (_f.gamma() < 0) {
      real u2, v2;
      if (fomg().NCoeffsInv() <= fbet().NCoeffsInv()) {
        solve2( _fic.delta, sig2, fbet(), fomg(), gbet(), gomg(), u2, v2,
                countn, countb);
      }
      else {
        solve2(-_fic.delta, sig2, fomg(), fbet(), gomg(), gbet(), v2, u2,
               countn, countb);
      }
      bet2 = _fic.nN * fbet().rev(u2);
      bet2a = AuxAngle::radians(bet2);
      AuxAngle psi2 = AuxAngle::radians(fomg().rev(v2));
      omg2a = AuxAngle(_fic.omg0 * _fic.flip * _f.gm().nup * psi2.y(),
                       _fic.omg0 * hypot(psi2.x(), _f.gm().nu * psi2.y()),
                       false);
      omg2 = omg2a.radians();
      alp2a = AuxAngle(_fic.alp0 * _fic.nN * _fic.flip * _t.kp * _f.gm().nup * psi2.x(),
                       _fic.alp0 * hypot(_t.kp * _f.gm().nu, _t.k * bet2a.x()),
                       false);
      alp2 = alp2a.radians();
    } else if (_f.gamma() == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_gic.s0);  // reduce to [-s0, s0)
      real u2, v2,
        deltax = fmax(-_fic.deltamax,
                     fmin(_fic.deltamax, _fic.delta + sig2n.second * _fic.deltashift));
      solve2u(deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2,
              countn, countb);
      /*
      cout << "AAA " << _gic.sig1 << " " << sig2 << " " << _gic.s0 << " " << deltax << " " << sig2n.first << " "
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
      // cout << "SIG2N " << sig2 << " " << _gic.s0 << " " << sig2n.first << " " << sig2n.second << "\n";
      int parity = fmod(sig2n.second, real(2)) ? -1 : 1;
      if (0)
        cout << "DD " << s12 << " " << u2 << " " << v2 << " " << sig2n.first << " " << sig2n.second << " " << bet2a.radians() << " " << omg2a.radians() << "\n";
      if (_fic.umbalt) {
        Ex = _fic.eE * parity;
        omg2a.y() *= Ex * ( _fic.omg0 ? -1 : 1);
        omg2a.x() *= ( _fic.omg0 ? -1 : 1);
        omg2 = _fic.omg0 * Math::pi() + Ex * omg2;
        bet2a.y() *= _fic.nN * parity * ( _fic.bet0 ? -1 : 1);
        bet2a.x() *= parity * ( _fic.bet0 ? -1 : 1);
        bet2 = //(_fic.bet0 < 0 ? Math::pi() : 0) +
          _fic.bet0 * Math::pi() +
          _fic.nN * (bet2 + sig2n.second * Math::pi());
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = AuxAngle((_fic.alp0 ? -1 : 1) * (_fic.nN * _t.kp) * (Ex / cosh(v2)),
                         (_fic.alp0 ? -1 : 1) * _t.k / cosh(u2),
                         false);
        alp2 = alp2a.radians();
        /*
        cout << "DATU " << _fic.alp0 << " " << _fic.nN << " " << _fic.eE << " "
             << ( _fic.omg0 ? -1 : 1) << " "
             << Ex << " " << u2 << " " << v2 << " " << alp2 << "\n";
        */
      } else {
        if (0)
          cout << "BB " << bet2a.radians() << " " << bet2 << " "
         << omg2a.radians() << " " << omg2 << " " << sig2n.second << " " << parity << "\n";
        Nx = _fic.nN * parity;
        omg2a.y() *= _fic.eE * parity * ( _fic.omg0 ? -1 : 1);
        omg2a.x() *= parity * ( _fic.omg0 ? -1 : 1);
        omg2 = _fic.omg0 * Math::pi() +
          _fic.eE * (omg2 + sig2n.second * Math::pi());
        bet2a.y() *= Nx * ( _fic.bet0 ? -1 : 1);
        bet2a.x() *= ( _fic.bet0 ? -1 : 1);
        bet2 = // (_fic.bet0 < 0 ? Math::pi() : 0) +
          _fic.bet0 * Math::pi() +
          Nx * bet2;
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        // _fic.alp0 = 0
        alp2a = AuxAngle((_fic.eE * _t.kp) / cosh(v2), _t.k * (Nx / cosh(u2)),
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
    //    cout << "UUU " << _fic.bet0 << " " <<  _fic.ibet << " " << bet2 << " " << bet2a.degrees() << "\n";
    if (ibet2)
      *ibet2 = _fic.ibet + int(round((bet2 - bet2a.degrees()) / Math::td)) -
        (_f.gamma() > 0 && _fic.flip < 0 ?
         (signbit(bet2a.y()) ? -1 : 1) - (signbit(_fic.bet1.y()) ? -1 : 1) : 0)/2;
    if (iomg2)
      *iomg2 = _fic.iomg + int(round((omg2 - omg2a.degrees()) / Math::td)) +
        (_f.gamma() < 0 && _fic.flip < 0 ?
         (signbit(omg2a.x()) ? -1 : 1) + (signbit(_fic.omg1.y()) ? -1 : 1) : 0)/2;
    if (ialp2)
      *ialp2 = _fic.ialp + int(round((alp2 - alp2a.degrees()) / Math::td)) -
        (_f.gamma() < 0 && _fic.nN < 0 ?
         (signbit(alp2a.y()) ? -1 : 1) - _fic.eE : 0)/2 +
        (_f.gamma() == 0 && _fic.umbalt ? (_fic.alp0 == _fic.nN * Ex ? _fic.alp0 : 0) : 0);
    //    cout << "AAA " << *ibet2 << " " << bet2a.degrees() << " " << _f.gamma() << "\n";
    //    cout << "BBB " << _fic.flip << " " <<       _fic.bet0 << " " <<       _fic.alp0 << " " << _fic.psi1.degrees() << " " <<       _fic.v0<< " " <<       _fic.u0<< " " <<       _fic.delta<< "\n";
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
      (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
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
    if (fabs(s0) >= stot) {
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
    // cout << "NEWT " << f0 << " " << g0 << "\n";
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
    //    cout << "GG " << _f.gamma() << "\n";
    if (_f.gamma() > 0) {
      real spsi = _t.k * bet2.y(),
        // In evaluating equivalent expressions, choose the one
        // with minimum cancelation
        cpsi = 0 + _t.k
        * sqrt(fmax(0,
                    _f.gm().nu < fabs(bet2.y()) ?
                    (bet2.x() - _f.gm().nu) * (bet2.x() + _f.gm().nu) :
                    (_f.gm().nup - bet2.y()) * (_f.gm().nup + bet2.y())));
      AuxAngle psi2 = AuxAngle(spsi, dir * cpsi, true),
        psi12 = psi2 - _fic.psi1;
      psi12.y() = fmax(real(0), psi12.y()) + 0; // convert -180deg to 180deg
      real tau12 = psi12.radians(),
        psi2r = _fic.psi1.radians() + tau12,
        v2= fbet().fwd(psi2r),
        u2 = fomg().inv(fbet()(v2) - _fic.delta),
        omg2 = _fic.eE * fomg().rev(u2);
      // cout << "PSI2 " << psi2.degrees() << " " << psi12.degrees() << " "
      //   << 1/psi12.y()  << " " << tau12 / Math::degree() << "\n";
      omg2a = AuxAngle::radians(omg2);
      bet2a = AuxAngle(_fic.bet0 * _fic.flip * _f.gm().nup * psi2.y(),
                       _fic.bet0 * hypot(psi2.x(), _f.gm().nu * psi2.y()),
                       true);
      alp2a = AuxAngle(_fic.alp0 * _fic.eE
                       * hypot(_t.k * _f.gm().nu, _t.kp * omg2a.x()),
                       _fic.alp0 * _fic.flip * _t.k * _f.gm().nup * psi2.x(),
                       true);
      real sig2 = gbet()(v2) + gomg()(u2);
      s12 = (sig2 - _gic.sig1)  * _t.b;
    } else if (_f.gamma() <= 0) {
      bet2a = AuxAngle(bet2.y(), bet2.x() * dir * _fic.nN, false);
      AuxAngle bet12 = bet2a - _fic.bet1;
      bet12.y() *= _fic.nN;
      bet12.y() = fmax(real(0), bet12.y()) + 0; // convert -180deg to 180deg
      real tau12 =  bet12.radians(),
        bet2r = _fic.bet1.radians() + _fic.nN * tau12;
      //      cout << "TAU12 " << tau12/Math::degree() << " "
      //           << bet2r/Math::degree() << "\n";
      if (_f.gamma() < 0) {
        real
          u2 = fbet().fwd(_fic.nN * bet2r),
          v2 = fomg().inv(fbet()(u2) - _fic.delta),
          psi2r = fomg().rev(v2);
        AuxAngle psi2 = AuxAngle::radians(psi2r);
        omg2a = AuxAngle(_fic.omg0 * _fic.flip * _f.gm().nup * psi2.y(),
                         _fic.omg0 * hypot(psi2.x(), _f.gm().nu * psi2.y()),
                         true);
        alp2a = AuxAngle(_fic.alp0 * _fic.nN * _fic.flip * _t.kp * _f.gm().nup * psi2.x(),
                         _fic.alp0 * hypot(_t.kp * _f.gm().nu, _t.k * bet2.x()),
                         true);
        real sig2 = gbet()(u2) + gomg()(v2);
        s12 = (sig2 - _gic.sig1)  * _t.b;
      } else {                  // _f.gamma() == 0
        // Could simplify this.  bet2 is in [-270,90]
        pair<real, real> bet2n =
          remx(_fic.nN * bet2r, Math::pi());  // reduce to [-pi/2, pi/2)
        int parity = fmod(bet2n.second, real(2)) ? -1 :  1;
        int alp0 = _fic.nN < 0 ? _fic.eE : 0;
        real deltax = fmax(-_fic.deltamax,
                           fmin(_fic.deltamax, _fic.delta + bet2n.second * _fic.deltashift)),
          u2 = fbet().fwd(bet2n.first),
          v2 = fomg().inv(fbet()(u2) - deltax),
          omg2 = _fic.eE * parity * fomg().rev(v2);
        omg2a = AuxAngle::radians(omg2);
        alp2a = AuxAngle((alp0 ? -1 : 1) * _fic.nN * _t.kp * _fic.eE * parity  /
                         cosh(v2),
                         (alp0 ? -1 : 1) * _t.k / cosh(u2),
                         true);
        if (0)
        cout << "AAA "
             << alp2a.degrees() << " " << alp0 << " "
             << (_fic.alp0 ? -1 : 1) * _fic.nN * _t.kp * _fic.eE * parity  / cosh(v2) << " "
             << (_fic.alp0 ? -1 : 1) * _t.k / cosh(u2) << "\n";
        real sig2 = gbet()(u2) + gomg()(v2) + bet2n.second * 2 * _gic.s0;
        s12 = (sig2 - _gic.sig1)  * _t.b;
      }
    } else  {
      // gamma = NaN
    }
    omg2a += AuxAngle::cardinal(1U);
    if (_f.gamma() <= 0 && dir * _fic.nN < 0) {
      //      omg2a.y() *= -1;
      //      alp2a += AuxAngle::cardinal(2U);
    }
  }

  TriaxialLineF::disttx
  TriaxialLineF::Hybrid(const ics& fic,
                        const AuxAngle& bet2,
                        AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a)
    const {
    disttx ret{Math::NaN(), Math::NaN(), 0};
    if (gamma() > 0) {
      real spsi = _t.k * bet2.y(),
        // In evaluating equivalent expressions, choose the one
        // with minimum cancelation
        cpsi = 0 + _t.k
        * sqrt(fmax(0,
                    gm().nu < fabs(bet2.y()) ?
                    (bet2.x() - gm().nu) * (bet2.x() + gm().nu) :
                    (gm().nup - bet2.y()) * (gm().nup + bet2.y())));
      if (spsi == 0 && cpsi == 0) cpsi = 1;
      AuxAngle psi2 = AuxAngle(spsi, cpsi, true),
        psi12 = psi2 - fic.psi1;
      psi12.y() = fmax(real(0), psi12.y()) + 0; // convert -180deg to 180deg
      real tau12 = psi12.radians(),
        psi2r = fic.psi1.radians() + tau12,
        v2 = fbet().fwd(psi2r),
        u2 = fomg().inv(fbet()(v2) - fic.delta),
        omg2 = fic.eE * fomg().rev(u2);
      omg2a = AuxAngle::radians(omg2);
      bet2a = AuxAngle(fic.bet0 * fic.flip * gm().nup * psi2.y(),
                       fic.bet0 * hypot(psi2.x(), gm().nu * psi2.y()),
                       true);
      alp2a = AuxAngle(fic.alp0 * fic.eE
                       * hypot(_t.k * gm().nu, _t.kp * omg2a.x()),
                       fic.alp0 * fic.flip * _t.k * gm().nup * psi2.x(),
                       true);
      ret.betw2 = v2;
      ret.omgw2 = u2;
    } else if (gamma() <= 0) {
      bet2a = AuxAngle(bet2.y(), bet2.x() * fic.nN, false);
      AuxAngle bet12 = bet2a - fic.bet1;
      bet12.y() *= fic.nN;
      bet12.y() = fmax(real(0), bet12.y()) + 0; // convert -180deg to 180deg
      real tau12 = bet12.radians(),
        bet2r = fic.bet1.radians() + fic.nN * tau12;
      /*
      cout << "TAU12 " << bet2.degrees() << " " << bet2a.degrees() << " "
           << bet12.degrees()
           << fic.nN << " " << fic.bet1.degrees() << " "
           << tau12/Math::degree() << " "
           << fic.alp1.degrees() << "\n";
      */
      if (gamma() < 0) {
        real
          u2 = fbet().fwd(fic.nN * bet2r),
          v2 = fomg().inv(fbet()(u2) - fic.delta),
          psi2r = fomg().rev(v2);
        AuxAngle psi2 = AuxAngle::radians(psi2r);
        omg2a = AuxAngle(fic.omg0 * fic.flip * gm().nup * psi2.y(),
                         fic.omg0 * hypot(psi2.x(), gm().nu * psi2.y()),
                         true);
        alp2a = AuxAngle(fic.alp0 * fic.nN * fic.flip * _t.kp * gm().nup * psi2.x(),
                         fic.alp0 * hypot(_t.kp * gm().nu, _t.k * bet2.x()),
                         true);
        ret.betw2 = u2;
        ret.omgw2 = v2;
      } else {                  // gamma() == 0
        // Could simplify this.  bet2 is in [-270,90]
        // reduce to [-pi/2, pi/2)
        pair<real, real> bet2n =
          TriaxialLine::remx(fic.nN * bet2r, Math::pi());
        if (0)
          cout << "\nZERO1 " << fic.nN * bet2r << " " << bet2n.first << " "
               << bet2n.second << "\n";
        int parity = fmod(bet2n.second, real(2)) ? -1 :  1;
        int alp0 = fic.nN < 0 ? fic.eE : 0;
        real deltax = fmax(-2*fic.deltamax,
                           fmin(2*fic.deltamax,
                               fic.delta + bet2n.second * fic.deltashift)),
          u2 = fbet().fwd(bet2n.first);
        real
          v2 = fomg().inv(fbet()(u2) - deltax),
          omg2 = fic.eE * parity * fomg().rev(v2);
        omg2a = AuxAngle::radians(omg2);
        omg2a.x() *= fic.omg0 ? -1 : 1;
        omg2a.y() *= fic.omg0 ? -1 : 1;
        alp2a = AuxAngle((alp0 ? -1 : 1) * fic.nN * _t.kp * fic.eE * parity  /
                         cosh(v2),
                         (alp0 ? -1 : 1) * _t.k / cosh(u2),
                         true);
        ret.betw2 = u2;
        ret.omgw2 = v2;
        ret.ind2 = int(bet2n.second);
        if (0)
          cout << "ZERO2 " << fic.delta << " " << deltax << " " << u2 << " "
               << v2 << " " << ret.ind2 << "\n"
               << "ZERO3 " << omg2a.degrees() << " "
               << fic.omg1.degrees() << " "
               << scientific
               << (omg2a-fic.omg1).radians() << "\n";
        // cout << "HHH " << omg2/Math::degree() << "\n";
      }
    } else {
      // gamma = NaN
    }
    omg2a += AuxAngle::cardinal(1U);
    // cout << "HHH " << omg2a.degrees() << "\n";
    return ret;
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

  TriaxialLineG::TriaxialLineG(const Triaxial& t, const Triaxial::gamblk& gam)
    : _t(t)
    , _gm(gam)
    , _gbet(_t.k2 , _t.kp2,  _t.e2, -_gm.gam)
    , _gomg(_t.kp2, _t.k2 , -_t.e2,  _gm.gam)
  {}

  Math::real TriaxialLineF::Hybrid0(const ics& fic,
                                    const AuxAngle& bet2, const AuxAngle& omg2)
  const {
    AuxAngle bet2a, omg2a, alp2a, omg2b(omg2);
    (void) Hybrid(fic, bet2, bet2a, omg2a, alp2a);
    (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
    omg2a -= omg2b;
    return omg2a.radians();
  }

  //      [bet2, omg2, alp2, betw2, omgw2, ind2] = obj.arcdist0(tau12, 1);
  TriaxialLineF::disttx
  TriaxialLineF::ArcPos0(const ics& fic, real tau12,
                         AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a,
                         bool betp)
    const {
    disttx ret{Math::NaN(), Math::NaN(), 0};
    if (gamma() > 0) {
      AuxAngle psi2;
      real psi2r, omg2, u2, v2;
      if (betp) {
        psi2 = AuxAngle::radians(tau12) + fic.psi1;
        psi2r = fic.psi1.radians() + tau12;
        v2 = fbet().fwd(psi2r);
        u2 = fomg().inv(fbet()(v2) - fic.delta);
        omg2 = fic.eE * fomg().rev(u2);
        omg2a = AuxAngle::radians(omg2);
      } else {
        omg2a = fic.eE > 0 ?
          fic.omg1 + AuxAngle::radians(tau12) :
          fic.omg1 - AuxAngle::radians(tau12);
        u2 = fomg().fwd(fic.eE * omg2a.radians());
        v2 = fbet().inv(fomg()(u2) + fic.delta);
        psi2r = fbet().rev(v2);
      }
      bet2a = AuxAngle(fic.bet0 * fic.flip * gm().nup * psi2.y(),
                       fic.bet0 * hypot(psi2.x(), gm().nu * psi2.y()),
                       true);
      alp2a = AuxAngle(fic.alp0 * fic.eE
                       * hypot(_t.k * gm().nu, _t.kp * omg2a.x()),
                       fic.alp0 * fic.flip * _t.k * gm().nup * psi2.x(),
                       true);
      ret.betw2 = v2;
      ret.omgw2 = u2;
    } else if (gamma() < 0) {
      AuxAngle psi2;
      real psi2r, bet2r, u2, v2;
      if (betp) {
        bet2a = fic.nN > 0 ?
          fic.bet1 + AuxAngle::radians(tau12) :
          fic.bet1 - AuxAngle::radians(tau12);
        bet2r = fic.bet1.radians() + fic.nN * tau12;
        u2 = fbet().fwd(fic.nN * bet2r);
        v2 = fomg().inv(fbet()(u2) - fic.delta);
        psi2r = fomg().rev(v2);
      } else {
        psi2 = AuxAngle::radians(tau12) + fic.psi1;
        v2 = fomg().fwd(psi2.radians());
        u2 = fbet().inv(fomg()(v2) + fic.delta);
        bet2r = fic.nN * fbet().rev(u2);
      }
      psi2 = AuxAngle::radians(psi2r);
      omg2a = AuxAngle(fic.omg0 * fic.flip * gm().nup * psi2.y(),
                       fic.omg0 * hypot(psi2.x(), gm().nu * psi2.y()),
                       true);
      alp2a = AuxAngle(fic.alp0 * fic.nN * fic.flip *
                       _t.kp * gm().nup * psi2.x(),
                       fic.alp0 * hypot(_t.kp * gm().nu, _t.k * bet2a.x()),
                       true);
      ret.betw2 = u2;
      ret.omgw2 = v2;
    } else if (gamma() == 0) {
      real u2, v2;
      int ii;
      if (betp) {
        bet2a = fic.nN > 0 ?
          fic.bet1 + AuxAngle::radians(tau12) :
          fic.bet1 - AuxAngle::radians(tau12);
        real bet2r = fic.bet1.radians() + fic.nN * tau12;

        // Could simplify this.  bet2 is in [-270,90]
        // reduce to [-pi/2, pi/2)
        pair<real, real> bet2n =
          TriaxialLine::remx(fic.nN * bet2r, Math::pi());
        if (0)
          cout << "\nZERO1 " << fic.nN * bet2r << " " << bet2n.first << " "
               << bet2n.second << "\n";
        int parity = fmod(bet2n.second, real(2)) ? -1 : 1,
          alp0 = fic.nN < 0 ? fic.eE : 0;
        real deltax = fmax(-2*fic.deltamax,
                           fmin(2*fic.deltamax,
                               fic.delta + bet2n.second * fic.deltashift));
        u2 = fbet().fwd(bet2n.first);
        v2 = fomg().inv(fbet()(u2) - deltax);
        real omg2 = fic.eE * parity * fomg().rev(v2);
        omg2a = AuxAngle::radians(omg2);
        omg2a.x() *= fic.omg0 ? -1 : 1;
        omg2a.y() *= fic.omg0 ? -1 : 1;
        alp2a = AuxAngle((alp0 ? -1 : 1) * fic.nN * _t.kp * fic.eE * parity /
                         cosh(v2),
                         (alp0 ? -1 : 1) * _t.k / cosh(u2),
                         true);
        ii = int(bet2n.second);
        if (0)
            cout << "UU " << fbet()(u2) - deltax << " " << bet2n.first << " " << u2 << " " << v2 << " " << fic.delta << " " << deltax << " " << fbet()(u2) << " " << bet2n.second << " " << fic.deltashift  << "\n";
      } else {
        omg2a = fic.eE > 0 ?
          fic.omg1 + AuxAngle::radians(tau12) :
          fic.omg1 - AuxAngle::radians(tau12);
        real omg2r = fic.omg1.radians() + fic.eE * tau12;
        pair<real, real> omg2n =
          TriaxialLine::remx(fic.eE * omg2r, Math::pi());
        int parity = fmod(omg2n.second, real(2)) ? -1 : 1,
          alp0 =  0;
        real deltax = fmax(-2*fic.deltamax,
                      fmin(2*fic.deltamax,
                           fic.delta + omg2n.second * fic.deltashift));
        v2 = fomg().fwd(omg2n.first);
        u2 = fbet().inv(fomg()(v2) + deltax);
        real bet2 = fic.nN * parity * fbet().rev(u2);
        bet2a = AuxAngle::radians(bet2);
        alp2a = AuxAngle((alp0 ? -1 : 1) * fic.eE * _t.kp / cosh(v2),
                         (alp0 ? -1 : 1) * _t.k * fic.nN * parity / cosh(u2));
        ii = int(omg2n.second);
      }
      ret.betw2 = u2;
      ret.omgw2 = v2;
      ret.ind2 = ii;
    } else {
      // gamma == NaN
    }
    omg2a += AuxAngle::cardinal(1U);
    return ret;
  }

  TriaxialLineF::ics::ics()
    : bet1(AuxAngle::NaN())
    , omg1(AuxAngle::NaN())
    , alp1(AuxAngle::NaN())
    , u0(Math::NaN())
    , v0(Math::NaN())
    , df(Math::NaN())
    , deltashift(Math::NaN())
    , deltamax(Math::NaN())
    , delta(Math::NaN())
    ,ibet(0)
    ,iomg(0)
    ,ialp(0)
    ,nN(0)
    ,eE(0)
    ,flip(0)
    ,bet0(0)
    ,omg0(0)
    ,alp0(0)
    ,umbalt(false)
  {}

  TriaxialLineF::ics::ics(const TriaxialLineF& f,
                         const AuxAngle& bet10, const AuxAngle& omg10,
                         const AuxAngle& alp10)
    : bet1(bet10)
      // omg10 - 90
    , omg1(omg10 - AuxAngle::cardinal(1U))
    , alp1(alp10)
    , ibet(0)
      // omg10 - 90 crossed omg10 = -180 if quadant of omg10 == 2
    , iomg(omg10.quadrant() == 2U ? -1 : 0)
    , ialp(0)
    , umbalt(false)
  {
    const real eps = numeric_limits<real>::epsilon();
    const Triaxial::gamblk& gm = f.gm();
    const Triaxial& t = f.t();
    if (bet1.y() == 0 && fabs(alp1.x()) <= Math::sq(eps))
      alp1.x() = - Math::sq(eps);
    eE = signbit(alp1.y()) ? -1 : 1;
    nN = signbit(alp1.x()) ? -1 : 1;
    if (gm.gam > 0) {
      flip = signbit(bet1.x()) ? -1 : 1;
      bet0 = flip;
      alp0 = 1;
      psi1 = AuxAngle(t.k * bet1.y(),
                       flip * alp1.x() *
                       hypot(t.k * bet1.x(), t.kp * omg1.x()),
                       false);
      v0 = f.fbet().fwd(psi1.radians());
      u0 = f.fomg().fwd(eE * omg1.radians());
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gm.gam < 0) {
      flip = signbit(omg1.x()) ? -1 : 1;
      omg0 = flip;
      alp0 = nN;
      psi1 = AuxAngle(t.kp * omg1.y(),
                       flip * alp1.y() *
                       hypot(t.k * bet1.x(), t.kp * omg1.x()),
                       false);
      v0 = f.fomg().fwd(psi1.radians());
      u0 = f.fbet().fwd(nN * bet1.radians());
      delta = f.fbet()(u0) - f.fomg()(v0);
    } else if (gm.gam == 0) {
      alp0 = umbalt && nN < 0 ? eE : 0;
      // N.B. factor of k*kp omitted
      df = f.fbet().Max() - f.fomg().Max();
      deltashift = 2*df - log(t.k2/t.kp2);
      deltamax = -2*log(eps); // consistent with geod_fun::lam + lamaux
      if (fabs(bet1.x()) < 8*eps && fabs(omg1.x()) < 8*eps) {
        //        bet0 = (int(round(bet1.y())) + nN) / 2 ? -1 : 1;
        //        omg0 = (int(round(omg1.y())) + eE) / 2 ? -1 : 1;
        bet0 = (int(round(bet1.y())) + nN) / 2; // -1, 0, or +1
        omg0 = (int(round(omg1.y())) + eE) / 2;
        delta = deltashift/2 - log(fabs(alp1.tan()));
      } else {
        bet0 = signbit(bet1.x()) ? (signbit(bet1.y()) ? -1 : 1) : 0;
        omg0 = signbit(omg1.x()) ? (signbit(omg1.y()) ? -1 : 1) : 0;
        delta = nN * f.fbet()(geod_fun::lamaux(bet1)) -
          eE * f.fomg()(geod_fun::lamaux(omg1));
        //        cout << "OOOA " << omg1.degrees() << " " << omg0 << "\n";
        if (0) {
          cout << "DEL2 " << omg1.degrees() << " " << geod_fun::lamaux(omg1) << "\n";
          cout << "DEL " <<
            geod_fun::lamaux(bet1) << " " <<
            geod_fun::lamaux(omg1) << " " <<
            f.fbet()(geod_fun::lamaux(bet1)) << " " <<
            f.fomg()(geod_fun::lamaux(omg1)) << " " <<
            delta << " " << deltamax << "\n";
        }
      }
    } else {
      // gamma = NaN
    }
  }

  void TriaxialLineF::ics::setquadrant(const TriaxialLineF& f, unsigned q) {
    const real eps = numeric_limits<real>::epsilon();
    real gam = f.gm().gam;
    alp1.setquadrant(q);
    if (bet1.y() == 0 && fabs(alp1.x()) <= Math::sq(eps))
      alp1.x() = - Math::sq(eps);
    int oE = eE, oN = nN;
    eE = signbit(alp1.y()) ? -1 : 1;
    nN = signbit(alp1.x()) ? -1 : 1;
    if (gam > 0) {
      psi1.x() *= flip * nN;
      v0 = f.fbet().fwd(psi1.radians());
      u0 *= eE/oE;
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gam < 0) {
      alp0 = nN;
      psi1.x() *= flip * eE;
      v0 = f.fomg().fwd(psi1.radians());
      u0 *= nN/oN;
      delta = f.fbet()(u0) - f.fomg()(v0);
    } else if (gam == 0) {
      alp0 = umbalt && nN < 0 ? eE : 0;
      if (fabs(bet1.x()) < 8*eps && fabs(omg1.x()) < 8*eps) {
        //        bet0 = (int(round(bet1.y())) + nN) / 2 ? -1 : 1;
        //        omg0 = (int(round(omg1.y())) + eE) / 2 ? -1 : 1;
        bet0 = (int(round(bet1.y())) + nN) / 2; // -1, 0, or +1
        omg0 = (int(round(omg1.y())) + eE) / 2;
      } else {
        delta = nN * f.fbet()(geod_fun::lamaux(bet1)) -
          eE * f.fomg()(geod_fun::lamaux(omg1));
      }
    } else {
      // gamma = NaN
    }
  }

  TriaxialLineG::ics::ics()
    : s0(Math::NaN())
    , sig1(Math::NaN())
  {}

  TriaxialLineG::ics::ics(const TriaxialLineG& g,
                          const TriaxialLineF::ics& fic)
  {
    if (g.gamma() > 0) {
      sig1 = g.gbet()(fic.v0) + g.gomg()(fic.u0);
      s0 = 0;
    } else if (g.gamma() < 0) {
      sig1 = g.gbet()(fic.u0) + g.gomg()(fic.v0);
      s0 = 0;
    } else if (g.gamma() == 0) {
      sig1 = fic.nN * g.gbet()(geod_fun::lamaux(fic.bet1)) +
        fic.eE * g.gomg()(geod_fun::lamaux(fic.omg1));
      s0 =  g.gbet().Max() + g.gomg().Max();
      if (0)
        cout << "GICS "
             << fic.nN << " " << g.gbet()(geod_fun::lamaux(fic.bet1)) << " "
             << fic.eE << " " << g.gomg()(geod_fun::lamaux(fic.omg1)) << " "
             << sig1 << " " << g.gbet().Max() << " " << g.gomg().Max() << " "
             << s0 << "\n";
    } else {
      // gamma = NaN
    }
  }

  Math::real TriaxialLineG::dist(ics ic, TriaxialLineF::disttx d) const {
    real sig2 = gbet()(d.betw2) + gomg()(d.omgw2) + d.ind2 * 2*ic.s0;
    return (sig2 - ic.sig1) * t().b;
  }
} // namespace GeographicLib

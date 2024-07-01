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
                             Angle bet1, Angle omg1, Angle alp1) {
    bet1.rnd();
    omg1.rnd();
    alp1.rnd();
    _t = t;
    Triaxial::gamblk gam(t, bet1, omg1, alp1);
    _f = TriaxialLineF(t, gam, 0.5, 1.5);
    _fic = TriaxialLineF::ics(_f, bet1, omg1, alp1);
    _g = TriaxialLineG(t, gam);
    _gic = TriaxialLineG::ics(_g, _fic);
  }

  TriaxialLine::TriaxialLine(const Triaxial& t, real bet1, real omg1, real alp1)
      : TriaxialLine(t,
                     ang::degrees(bet1),
                     ang::degrees(omg1),
                     ang::degrees(alp1))
    {
      _fic.ibet = int(round((bet1 - _fic.bet1.degrees()) / Math::td));
      _fic.iomg += int(round((omg1 - ang::degrees(omg1).degrees()) /
                         Math::td));
      _fic.ialp = int(round((alp1 - _fic.alp1.degrees()) / Math::td));
    }

  void TriaxialLine::pos1(Angle& bet1, Angle& omg1, Angle& alp1,
                          int* ibet1, int* iomg1, int* ialp1) const {
    bet1 = _fic.bet1; omg1 = _fic.omg1 + ang::cardinal(1U);
    alp1 = _fic.alp1;
    if (ibet1) *ibet1 = _fic.ibet;
    if (iomg1) *iomg1 = _fic.iomg + (omg1.quadrant() == 2U ? 1 : 0);
    if (ialp1) *ialp1 = _fic.ialp;
  }

  void TriaxialLine::pos1(real& bet1, real& omg1, real& alp1,
                          bool unroll) const {
    ang bet1a, omg1a, alp1a;
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
                              Angle& bet2a, Angle& omg2a,
                              Angle& alp2a,
                              int* ibet2, int* iomg2, int* ialp2,
                              int* countn, int* countb)
  const {
    // Compute points at distance s12
    real sig2 = _gic.sig1 + s12/_t.b;
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
      omg2a = ang::radians(omg2);
      ang psi2 = ang::radians(fbet().rev(v2));
      bet2a = ang::aux(_fic.bet0 * _fic.flip * _f.gm().nup * psi2.y(),
                       _fic.bet0 * hypot(psi2.x(), _f.gm().nu * psi2.y()),
                       false);
      bet2 = bet2a.radians0();
      alp2a = ang::aux(_fic.alp0 * _fic.eE *
                       hypot(_t.k * _f.gm().nu, _t.kp * omg2a.x()),
                       _fic.alp0 * _fic.flip * _t.k * _f.gm().nup * psi2.x(),
                       false);
      alp2 = alp2a.radians0();
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
      bet2a = ang::radians(bet2);
      ang psi2 = ang::radians(fomg().rev(v2));
      omg2a = ang::aux(_fic.omg0 * _fic.flip * _f.gm().nup * psi2.y(),
                       _fic.omg0 * hypot(psi2.x(), _f.gm().nu * psi2.y()),
                       false);
      omg2 = omg2a.radians0();
      alp2a = ang::aux(_fic.alp0 * _fic.nN * _fic.flip * _t.kp * _f.gm().nup * psi2.x(),
                       _fic.alp0 * hypot(_t.kp * _f.gm().nu, _t.k * bet2a.x()),
                       false);
      alp2 = alp2a.radians0();
    } else if (_f.gamma() == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_g.s0);  // reduce to [-s0, s0)
      if (sig2n.first - _g.s0 >= -5 * numeric_limits<real>::epsilon()) {
        sig2n.first = -_g.s0;
        ++sig2n.second;
      }
      real u2, v2,
        deltax = Triaxial::clamp(_fic.delta + sig2n.second * _f.deltashift, 1);
      solve2u(deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2,
              countn, countb);
      bet2 = fbet().rev(u2); omg2 = fomg().rev(v2);
      bet2a = ang::lam(u2); omg2a = ang::lam(v2);
      int parity = fmod(sig2n.second, real(2)) ? -1 : 1;
      if (_fic.umbalt) {
        Ex = _fic.eE * parity;
        omg2a.reflect(Ex * ( _fic.omg0 ? -1 : 1) < 0, _fic.omg0);
        omg2 = _fic.omg0 * Math::pi() + Ex * omg2;
        bet2a.reflect(_fic.nN * parity * ( _fic.bet0 ? -1 : 1) < 0,
                      parity * ( _fic.bet0 ? -1 : 1) < 0);
        bet2 = //(_fic.bet0 < 0 ? Math::pi() : 0) +
          _fic.bet0 * Math::pi() +
          _fic.nN * (bet2 + sig2n.second * Math::pi());
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = ang::aux((_fic.alp0 ? -1 : 1) * (_fic.nN * _t.kp) * (Ex / cosh(v2)),
                         (_fic.alp0 ? -1 : 1) * _t.k / cosh(u2),
                         false);
        alp2 = alp2a.radians0();
      } else {
        Nx = _fic.nN * parity;
        omg2a.reflect(_fic.eE * parity * ( _fic.omg0 ? -1 : 1) < 0,
                      parity * ( _fic.omg0 ? -1 : 1) < 0);
        omg2 = _fic.omg0 * Math::pi() +
          _fic.eE * (omg2 + sig2n.second * Math::pi());
        bet2a.reflect(Nx * ( _fic.bet0 ? -1 : 1) < 0, _fic.bet0);
        bet2 = // (_fic.bet0 < 0 ? Math::pi() : 0) +
          _fic.bet0 * Math::pi() +
          Nx * bet2;
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        // _fic.alp0 = 0
        alp2a = ang::aux((_fic.eE * _t.kp) / cosh(v2), _t.k * (Nx / cosh(u2)),
                         false);
        alp2 = alp2a.radians0();
      }
    } else {
      // gamma = NaN
    }
    bet2a.rnd();
    omg2a.rnd();
    alp2a.rnd();

    swap(omg2a.x(), omg2a.y());
    omg2a.reflect(false, true);
    omg2 = omg2 / Math::degree() + 90;
    bet2 = bet2 / Math::degree();
    alp2 = alp2 / Math::degree();
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
  }

  void TriaxialLine::Position(real s12, real& bet2, real& omg2, real& alp2,
                              bool unroll,
                              int* countn, int* countb) const {
    ang bet2a, omg2a, alp2a;
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
    real pi2 = Triaxial::BigValue(),
      sbet = gbet.Max(), somg = gomg.Max(), stot = sbet + somg,
      dbet = fbet.Max(), domg = fomg.Max(), del  = dbet - domg;
    bool debug = false;
    if (debug)
      cout << "HEREQ " << pi2 << " " << d0 << " " << s0 << " " << stot << " " << sbet-somg << " "
         << fabs(s0) - stot << " "
         << fabs((1 - 2 * signbit(d0)) * s0 - (sbet - somg)) << "\n";
    if (fabs(s0) - stot >= -5 * numeric_limits<real>::epsilon()) {
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
      if (debug) cout << "UV1 " << u << " " << v << "\n";
    } else if (fabs(d0) > 2*pi2/3 &&
               fabs((1 - 2 * signbit(d0)) * s0 - (sbet - somg)) <=
               5 * numeric_limits<real>::epsilon()) {
      if (d0 > 0) {
        u = 2*d0/3; v = -1*d0/3;
      } else {
        u = 1*d0/3; v = -2*d0/3;
      }
      if (debug) cout << "UV2 " << u << " " << v << "\n";
    } else if ((1 - 2 * signbit(d0)) * s0 < sbet - somg) {
      // Use u as independent variable if
      //   d0 < 0 ? (s0 > -sbet + somg) :
      //            (s0 <  sbet - somg)
      // or
      //   sign(d0) * s0 < sbet - somg
      newt2(d0, s0, fbet, fomg, gbet, gomg, 0, -pi2, pi2, pi2, pi2,
            u, v, countn, countb);
      if (debug) cout << "UV3 " << u << " " << v << "\n";
    } else {
      // Otherwise, use v is the independent variable
      newt2(-d0, s0, fomg, fbet, gomg, gbet, 0, -pi2, pi2, pi2, pi2,
            v, u, countn, countb);
      if (debug) cout << "UV4 " << u << " " << v << "\n";
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
    // TO DO: supplement with 1 iteration of 2d Newton
  }

  void TriaxialLine::Hybrid(const Angle& bet2,
                            Angle& bet2a, Angle& omg2a, Angle& alp2a,
                            real& s12)
    const {
    TriaxialLineF::disttx d = _f.Hybrid(_fic, bet2, bet2a, omg2a, alp2a);
    s12 = _g.dist(_gic, d);
  }

  TriaxialLineF::disttx
  TriaxialLineF::Hybrid(const ics& fic,
                        const Angle& bet2,
                        Angle& bet2a, Angle& omg2a, Angle& alp2a)
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
      ang psi2 = ang::aux(spsi, cpsi, true),
        psi12 = psi2 - fic.psi1;
      // convert -180deg to 180deg
      if (signbit(psi12.y())) psi12 = ang(0, psi12.x(), psi12.n());
      real tau12 = psi12.radians0(),
        psi2r = fic.psi1.radians0() + tau12,
        v2 = fbet().fwd(psi2r),
        u2 = fomg().inv(fbet()(v2) - fic.delta),
        omg2 = fic.eE * fomg().rev(u2);
      omg2a = ang::radians(omg2);
      bet2a = ang::aux(fic.bet0 * fic.flip * gm().nup * psi2.y(),
                       fic.bet0 * hypot(psi2.x(), gm().nu * psi2.y()),
                       true);
      alp2a = ang::aux(fic.alp0 * fic.eE
                       * hypot(_t.k * gm().nu, _t.kp * omg2a.x()),
                       fic.alp0 * fic.flip * _t.k * gm().nup * psi2.x(),
                       true);
      ret.betw2 = v2;
      ret.omgw2 = u2;
    } else if (gamma() <= 0) {
      bet2a = ang::aux(bet2.y(), bet2.x() * fic.nN, false);
      ang bet12 = bet2a - fic.bet1;
      bet12.reflect(fic.nN < 0, false);
      // convert -180deg to 180deg
      if (signbit(bet12.y())) bet12 = ang(0, bet12.x(), bet12.n());
      real tau12 = bet12.radians0(),
        bet2r = fic.bet1.radians0() + fic.nN * tau12;
      if (gamma() < 0) {
        real
          u2 = fbet().fwd(fic.nN * bet2r),
          v2 = fomg().inv(fbet()(u2) - fic.delta),
          psi2r = fomg().rev(v2);
        ang psi2 = ang::radians(psi2r);
        omg2a = ang::aux(fic.omg0 * fic.flip * gm().nup * psi2.y(),
                         fic.omg0 * hypot(psi2.x(), gm().nu * psi2.y()),
                         true);
        alp2a = ang::aux(fic.alp0 * fic.nN * fic.flip * _t.kp * gm().nup * psi2.x(),
                         fic.alp0 * hypot(_t.kp * gm().nu, _t.k * bet2.x()),
                         true);
        ret.betw2 = u2;
        ret.omgw2 = v2;
      } else {                  // gamma() == 0
        // Could simplify this.  bet2 is in [-270,90]
        // reduce to [-pi/2, pi/2)
        pair<real, real> bet2n =
          TriaxialLine::remx(fic.nN * bet2r, Math::pi());
        int parity = fmod(bet2n.second, real(2)) ? -1 :  1;
        int alp0 = fic.nN < 0 ? fic.eE : 0;
        real deltax = Triaxial::clamp(fic.delta + bet2n.second * deltashift, 2),
          u2 = fbet().fwd(bet2n.first);
        real
          v2 = fomg().inv(fbet()(u2) - deltax),
          omg2 = fic.eE * parity * fomg().rev(v2);
        omg2a = ang::radians(omg2);
        omg2a.reflect(fic.omg0, fic.omg0);
        alp2a = ang::aux((alp0 ? -1 : 1) * fic.nN * _t.kp * fic.eE * parity  /
                         cosh(v2),
                         (alp0 ? -1 : 1) * _t.k / cosh(u2),
                         true);
        ret.betw2 = u2;
        ret.omgw2 = v2;
        ret.ind2 = int(bet2n.second);
      }
    } else {
      // gamma = NaN
    }
    omg2a += ang::cardinal(1U);
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
      df = _gm.gam == 0 ? _fbet.Max() - _fomg.Max() : 0;
      deltashift = _gm.gam == 0 ? 2*df - log(t.k2/t.kp2) : 0;
    }

  TriaxialLineG::TriaxialLineG(const Triaxial& t, const Triaxial::gamblk& gam)
    : _t(t)
    , _gm(gam)
    , _gbet(_t.k2 , _t.kp2,  _t.e2, -_gm.gam)
    , _gomg(_t.kp2, _t.k2 , -_t.e2,  _gm.gam)
    , s0(_gm.gam == 0 ? _gbet.Max() + _gomg.Max() : 0)
  {}

  Math::real TriaxialLineF::Hybrid0(const ics& fic,
                                    const Angle& bet2, const Angle& omg2)
  const {
    ang bet2a, omg2a, alp2a, omg2b(omg2);
    (void) Hybrid(fic, bet2, bet2a, omg2a, alp2a);
    (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
    omg2a -= omg2b;
    return omg2a.radians0();
  }

  //      [bet2, omg2, alp2, betw2, omgw2, ind2] = obj.arcdist0(tau12, 1);
  TriaxialLineF::disttx
  TriaxialLineF::ArcPos0(const ics& fic, real tau12,
                         Angle& bet2a, Angle& omg2a, Angle& alp2a,
                         bool betp)
    const {
    disttx ret{Math::NaN(), Math::NaN(), 0};
    if (gamma() > 0) {
      ang psi2;
      real psi2r, omg2, u2, v2;
      if (betp) {
        psi2 = ang::radians(tau12) + fic.psi1;
        psi2r = fic.psi1.radians0() + tau12;
        v2 = fbet().fwd(psi2r);
        u2 = fomg().inv(fbet()(v2) - fic.delta);
        omg2 = fic.eE * fomg().rev(u2);
        omg2a = ang::radians(omg2);
      } else {
        omg2a = fic.eE > 0 ?
          fic.omg1 + ang::radians(tau12) :
          fic.omg1 - ang::radians(tau12);
        u2 = fomg().fwd(fic.eE * fic.omg1.radians0() + tau12);
        v2 = fbet().inv(fomg()(u2) + fic.delta);
        psi2r = fbet().rev(v2);
      }
      bet2a = ang::aux(fic.bet0 * fic.flip * gm().nup * psi2.y(),
                       fic.bet0 * hypot(psi2.x(), gm().nu * psi2.y()),
                       true);
      alp2a = ang::aux(fic.alp0 * fic.eE
                       * hypot(_t.k * gm().nu, _t.kp * omg2a.x()),
                       fic.alp0 * fic.flip * _t.k * gm().nup * psi2.x(),
                       true);
      ret.betw2 = v2;
      ret.omgw2 = u2;
    } else if (gamma() < 0) {
      ang psi2;
      real psi2r, bet2r, u2, v2;
      if (betp) {
        bet2a = fic.nN > 0 ?
          fic.bet1 + ang::radians(tau12) :
          fic.bet1 - ang::radians(tau12);
        bet2r = fic.bet1.radians0() + fic.nN * tau12;
        u2 = fbet().fwd(fic.nN * bet2r);
        v2 = fomg().inv(fbet()(u2) - fic.delta);
        psi2r = fomg().rev(v2);
      } else {
        psi2 = ang::radians(tau12) + fic.psi1;
        v2 = fomg().fwd(psi2.radians0());
        u2 = fbet().inv(fomg()(v2) + fic.delta);
        bet2r = fic.nN * fbet().rev(u2);
      }
      psi2 = ang::radians(psi2r);
      omg2a = ang::aux(fic.omg0 * fic.flip * gm().nup * psi2.y(),
                       fic.omg0 * hypot(psi2.x(), gm().nu * psi2.y()),
                       true);
      alp2a = ang::aux(fic.alp0 * fic.nN * fic.flip *
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
          fic.bet1 + ang::radians(tau12) :
          fic.bet1 - ang::radians(tau12);
        real bet2r = fic.bet1.radians0() + fic.nN * tau12;

        // Could simplify this.  bet2 is in [-270,90]
        // reduce to [-pi/2, pi/2)
        pair<real, real> bet2n =
          TriaxialLine::remx(fic.nN * bet2r, Math::pi());
        int parity = fmod(bet2n.second, real(2)) ? -1 : 1,
          alp0 = fic.nN < 0 ? fic.eE : 0;
        real deltax = Triaxial::clamp(fic.delta + bet2n.second * deltashift, 2);
        u2 = fbet().fwd(bet2n.first);
        v2 = fomg().inv(fbet()(u2) - deltax);
        real omg2 = fic.eE * parity * fomg().rev(v2);
        omg2a = ang::radians(omg2);
        omg2a.reflect(fic.omg0, fic.omg0);
        alp2a = ang::aux((alp0 ? -1 : 1) * fic.nN * _t.kp * fic.eE * parity /
                         cosh(v2),
                         (alp0 ? -1 : 1) * _t.k / cosh(u2),
                         true);
        ii = int(bet2n.second);
      } else {
        omg2a = fic.eE > 0 ?
          fic.omg1 + ang::radians(tau12) :
          fic.omg1 - ang::radians(tau12);
        real omg2r = fic.omg1.radians0() + fic.eE * tau12;
        pair<real, real> omg2n =
          TriaxialLine::remx(fic.eE * omg2r, Math::pi());
        int parity = fmod(omg2n.second, real(2)) ? -1 : 1,
          alp0 =  0;
        real deltax = Triaxial::clamp(fic.delta + omg2n.second * deltashift, 2);
        v2 = fomg().fwd(omg2n.first);
        u2 = fbet().inv(fomg()(v2) + deltax);
        real bet2 = fic.nN * parity * fbet().rev(u2);
        bet2a = ang::radians(bet2);
        alp2a = ang((alp0 ? -1 : 1) * fic.eE * _t.kp / cosh(v2),
                         (alp0 ? -1 : 1) * _t.k * fic.nN * parity / cosh(u2));
        ii = int(omg2n.second);
      }
      ret.betw2 = u2;
      ret.omgw2 = v2;
      ret.ind2 = ii;
    } else {
      // gamma == NaN
    }
    omg2a += ang::cardinal(1U);
    return ret;
  }

  TriaxialLineF::ics::ics()
    : bet1(ang::NaN())
    , omg1(ang::NaN())
    , alp1(ang::NaN())
    , u0(Math::NaN())
    , v0(Math::NaN())
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
                         const Angle& bet10, const Angle& omg10,
                         const Angle& alp10)
    : bet1(bet10)
      // omg10 - 90
    , omg1(omg10 - ang::cardinal(1U))
    , alp1(alp10)
    , ibet(0)
      // omg10 - 90 crossed omg10 = -180 if quadant of omg10 == 2
    , iomg(omg10.quadrant() == 2U ? -1 : 0)
    , ialp(0)
    , umbalt(false)
  {
    alp1.rnd();
    const real eps = numeric_limits<real>::epsilon();
    const Triaxial::gamblk& gm = f.gm();
    const Triaxial& t = f.t();
    if (bet1.y() == 0 && fabs(alp1.x()) <= Math::sq(eps))
      alp1 = ang(alp1.y(), - Math::sq(eps), alp1.n());
    eE = signbit(alp1.y()) ? -1 : 1;
    nN = signbit(alp1.x()) ? -1 : 1;
    if (gm.gam > 0) {
      flip = signbit(bet1.x()) ? -1 : 1;
      bet0 = flip;
      alp0 = 1;
      psi1 = ang::aux(t.k * bet1.y(),
                       flip * alp1.x() *
                       hypot(t.k * bet1.x(), t.kp * omg1.x()),
                       false);
      v0 = f.fbet().fwd(psi1.radians0());
      u0 = f.fomg().fwd(eE * omg1.radians0());
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gm.gam < 0) {
      flip = signbit(omg1.x()) ? -1 : 1;
      omg0 = flip;
      alp0 = nN;
      // Need Angle(0, 0) to be treated like Angle(0, 1) here.
      psi1 = ang::aux(t.kp * omg1.y(),
                       flip * alp1.y() *
                       hypot(t.k * bet1.x(), t.kp * omg1.x()),
                       false);
      v0 = f.fomg().fwd(psi1.radians0());
      u0 = f.fbet().fwd(nN * bet1.radians0());
      delta = f.fbet()(u0) - f.fomg()(v0);
    } else if (gm.gam == 0) {
      alp0 = umbalt && nN < 0 ? eE : 0;
      // N.B. factor of k*kp omitted
      if (fabs(bet1.x()) < 8*eps && fabs(omg1.x()) < 8*eps) {
        //        bet0 = (int(round(bet1.y())) + nN) / 2 ? -1 : 1;
        //        omg0 = (int(round(omg1.y())) + eE) / 2 ? -1 : 1;
        bet0 = (int(round(bet1.y())) + nN) / 2; // -1, 0, or +1
        omg0 = (int(round(omg1.y())) + eE) / 2;
        delta = f.deltashift/2 - log(fabs(alp1.tan()));
      } else {
        bet0 = signbit(bet1.x()) ? (signbit(bet1.y()) ? -1 : 1) : 0;
        omg0 = signbit(omg1.x()) ? (signbit(omg1.y()) ? -1 : 1) : 0;

        ang btmp(bet1 - ang::cardinal(bet0 ? 2U : 0U)),
        otmp(omg1 - ang::cardinal(omg0 ? 2U : 0U));

        delta = nN * f.fbet()(geod_fun::lamaux(btmp)) -
          eE * f.fomg()(geod_fun::lamaux(otmp));
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
      alp1 = ang(alp1.y(), -Math::sq(eps), alp1.n());
    int oE = eE, oN = nN;
    eE = signbit(alp1.y()) ? -1 : 1;
    nN = signbit(alp1.x()) ? -1 : 1;
    if (gam > 0) {
      psi1.reflect(false, flip * nN < 0);
      v0 = f.fbet().fwd(psi1.radians0());
      u0 *= eE/oE;
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gam < 0) {
      alp0 = nN;
      psi1.reflect(false, flip * eE < 0);
      v0 = f.fomg().fwd(psi1.radians0());
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
    : sig1(Math::NaN())
  {}

  TriaxialLineG::ics::ics(const TriaxialLineG& g,
                          const TriaxialLineF::ics& fic)
  {
    if (g.gamma() > 0) {
      sig1 = g.gbet()(fic.v0) + g.gomg()(fic.u0);
    } else if (g.gamma() < 0) {
      sig1 = g.gbet()(fic.u0) + g.gomg()(fic.v0);
    } else if (g.gamma() == 0) {
      ang btmp(fic.bet1 - ang::cardinal(fic.bet0 ? 2U : 0U)),
        otmp(fic.omg1 - ang::cardinal(fic.omg0 ? 2U : 0U));
      sig1 = fic.nN * g.gbet()(geod_fun::lamaux(btmp)) +
        fic.eE * g.gomg()(geod_fun::lamaux(otmp));
    } else {
      // gamma = NaN
    }
  }

  Math::real TriaxialLineG::dist(ics ic, TriaxialLineF::disttx d) const {
    real sig2 = gbet()(d.betw2) + gomg()(d.omgw2) + d.ind2 * 2*s0;
    return (sig2 - ic.sig1) * t().b;
  }
} // namespace GeographicLib

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

  TriaxialLine::TriaxialLine(TriaxialLineF f, TriaxialLineF::fics fic,
                             TriaxialLineG g, TriaxialLineG::gics gic)
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
    _fic = TriaxialLineF::fics(_f, bet1, omg1, alp1);
    _g = TriaxialLineG(t, gam);
    _gic = TriaxialLineG::gics(_g, _fic);
  }

  TriaxialLine::TriaxialLine(const Triaxial& t, real bet1, real omg1,
                             real alp1)
      : TriaxialLine(t, ang(bet1), ang(omg1), ang(alp1))
    {}

  void TriaxialLine::pos1(Angle& bet1, Angle& omg1, Angle& alp1) const {
    bet1 = _fic.bet1; omg1 = _fic.omg1 + ang::cardinal(1);
    alp1 = _fic.alp1;
  }

  void TriaxialLine::pos1(real& bet1, real& omg1, real& alp1,
                          bool unroll) const {
    ang bet1a, omg1a, alp1a;
    pos1(bet1a, omg1a, alp1a);
    if (!unroll) {
      (void) Triaxial::AngNorm(bet1a, omg1a, alp1a);
      bet1a.setn(); omg1a.setn(); alp1a.setn();
    }
    bet1 = real(bet1a);
    omg1 = real(omg1a);
    alp1 = real(alp1a);
  }

  void TriaxialLine::Position(real s12,
                              Angle& bet2a, Angle& omg2a, Angle& alp2a,
                              int* countn, int* countb)
  const {
    // Compute points at distance s12
    real sig2 = _gic.sig1 + s12/_t.b;
    real bet2, omg2;
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
      // Already normalized
      bet2a = ang(_f.gm().nup * psi2.s(),
                  _fic.bet0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                  0, true).rebase(_fic.bet0);
      alp2a = ang(_fic.eE * hypot(_t.k * _f.gm().nu, _t.kp * omg2a.c()),
                  _fic.bet0.c() * _t.k * _f.gm().nup * psi2.c())
        .rebase(_fic.alp0);
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
      // Already normalized
      omg2a = ang(_f.gm().nup * psi2.s(),
                  _fic.omg0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                  0, true).rebase(_fic.omg0);
      alp2a = ang(_fic.omg0.c() * _t.kp * _f.gm().nup * psi2.c(),
                  _fic.nN * hypot(_t.kp * _f.gm().nu, _t.k * bet2a.c()))
        .rebase(_fic.alp0);
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
      if (_t.umbalt) {
        Ex = _fic.eE * parity;
        omg2a = omg2a.flipsign(Ex).rebase(_fic.omg0);
        bet2a += ang::cardinal(2 * sig2n.second);
        bet2a = bet2a.flipsign(_fic.nN) + _fic.bet0;
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = ang(_fic.nN * _t.kp * Ex / cosh(v2),
                    _t.k / cosh(u2)).rebase(_fic.alp0);
      } else {
        Nx = _fic.nN * parity;
        omg2a += ang::cardinal(2 * sig2n.second);
        omg2a = omg2a.flipsign(_fic.eE) + _fic.omg0;
        bet2a = bet2a.reflect((_fic.bet0.c() * Nx) < 0, _fic.bet0.c() < 0)
          .rebase(_fic.bet0);
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = ang(_fic.eE * _t.kp / cosh(v2),
                    _t.k * Nx / cosh(u2)).rebase(_fic.alp0);
      }
    } else {
      // gamma = NaN
    }
    bet2a.rnd();
    omg2a.rnd();
    alp2a.rnd();
    omg2a += ang::cardinal(1);
  }

  void TriaxialLine::Position(real s12, real& bet2, real& omg2, real& alp2,
                              bool unroll,
                              int* countn, int* countb) const {
    ang bet2a, omg2a, alp2a;
    Position(s12, bet2a, omg2a, alp2a, countn, countb);
    if (!unroll) {
      (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
      bet2a.setn(); omg2a.setn(); alp2a.setn();
    }
    bet2 = real(bet2a);
    omg2 = real(omg2a);
    alp2 = real(alp2a);
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
      cout << "HEREQ "
           << pi2 << " " << d0 << " " << s0 << " "
           << stot << " " << sbet-somg << " "
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
  TriaxialLineF::Hybrid(const fics& fic,
                        const Angle& bet2,
                        Angle& bet2a, Angle& omg2a, Angle& alp2a)
    const {
    ang tau12;
    if (gamma() > 0) {
      real spsi = _t.k * bet2.s(),
        // In evaluating equivalent expressions, choose the one
        // with minimum cancelation
        cpsi = 0 + _t.k
        * sqrt(fmax(0,
                    gm().nu < fabs(bet2.s()) ?
                    (bet2.c() - gm().nu) * (bet2.c() + gm().nu) :
                    (gm().nup - bet2.s()) * (gm().nup + bet2.s())));
      // Need Angle(0, 0) to be treated like Angle(0, 1) here.
      ang psi2 = ang(spsi, cpsi),
        psi12 = psi2 - fic.psi1;
      // convert -180deg to 180deg
      if (signbit(psi12.s()))
        psi12 = ang(0, copysign(real(1), psi12.c()), 0, true);
      tau12 = psi12;
    } else if (gamma() <= 0) {
      bet2a = bet2; bet2a.reflect(false, fic.nN < 0);
      ang bet12 = bet2a - fic.bet1;
      bet12.reflect(fic.nN < 0);
      // convert -180deg to 180deg
      if (signbit(bet12.s()))
        bet12 = ang(0, copysign(real(1), bet12.c()), 0, true);
      tau12 = bet12;
    } else {
      tau12 = ang::NaN();
    }
    return ArcPos0(fic, tau12.base(), bet2a, omg2a, alp2a, true);
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

  Math::real TriaxialLineF::Hybrid0(const fics& fic,
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
  TriaxialLineF::ArcPos0(const fics& fic, const Angle& tau12,
                         Angle& bet2a, Angle& omg2a, Angle& alp2a,
                         bool betp)
    const {
    disttx ret{Math::NaN(), Math::NaN(), 0};
    if (gamma() > 0) {
      ang psi2;
      real u2, v2;
      if (betp) {
        psi2 = tau12 + fic.psi1;
        v2 = fbet().fwd(psi2.radians());
        u2 = fomg().inv(fbet()(v2) - fic.delta);
        omg2a = ang::radians(fic.eE * fomg().rev(u2));
      } else {
        omg2a = fic.omg1 + tau12.flipsign(fic.eE);
        u2 = fomg().fwd(fic.eE * omg2a.radians());
        v2 = fbet().inv(fomg()(u2) + fic.delta);
        psi2 = ang::radians(fbet().rev(v2));
      }
      // Already normalized
      bet2a = ang(gm().nup * psi2.s(),
                  fic.bet0.c() * hypot(psi2.c(), gm().nu * psi2.s()),
                  0, true).rebase(fic.bet0);
      alp2a = ang(fic.eE * hypot(_t.k * gm().nu, _t.kp * omg2a.c()),
                  fic.bet0.c() * _t.k * gm().nup * psi2.c()).rebase(fic.alp0);
      ret.betw2 = v2;
      ret.omgw2 = u2;
    } else if (gamma() < 0) {
      ang psi2;
      real u2, v2;
      if (betp) {
        bet2a = fic.bet1 + tau12.flipsign(fic.nN);
        u2 = fbet().fwd(fic.nN * bet2a.radians());
        v2 = fomg().inv(fbet()(u2) - fic.delta);
        psi2 = ang::radians(fomg().rev(v2));
      } else {
        psi2 = tau12 + fic.psi1;
        v2 = fomg().fwd(psi2.radians());
        u2 = fbet().inv(fomg()(v2) + fic.delta);
        bet2a = ang::radians(fic.nN * fbet().rev(u2));
      }
      // Already normalized
      omg2a = ang(gm().nup * psi2.s(),
                  fic.omg0.c() * hypot(psi2.c(), gm().nu * psi2.s()),
                  0, true).rebase(fic.omg0);
      alp2a = ang(fic.omg0.c() * _t.kp * gm().nup * psi2.c(),
                  fic.nN * hypot(_t.kp * gm().nu, _t.k * bet2a.c()))
        .rebase(fic.alp0);
      ret.betw2 = u2;
      ret.omgw2 = v2;
    } else if (gamma() == 0) {
      real u2, v2;
      int ii;
      if (betp) {
        bet2a = fic.bet1 + tau12.flipsign(fic.nN);
        pair<real, real> bet2n =
          TriaxialLine::remx(fic.nN * (bet2a - fic.bet0).radians(),
                             Math::pi());
        int parity = fmod(bet2n.second, real(2)) ? -1 : 1;
        real deltax = Triaxial::clamp(fic.delta + bet2n.second * deltashift,
                                      2);
        u2 = fbet().fwd(bet2n.first);
        v2 = fomg().inv(fbet()(u2) - deltax);
        omg2a = ang::radians(fic.eE * parity * fomg().rev(v2))
          .rebase(fic.omg0);
        // umbalt definition of alp0
        ang alp0x(fic.alp1.nearest(2U));
        alp2a = ang(fic.nN * _t.kp * fic.eE * parity / cosh(v2),
                    _t.k / cosh(u2)).rebase(alp0x);
        ii = int(bet2n.second);
      } else {
        omg2a = fic.omg1 + tau12.flipsign(fic.eE);
        pair<real, real> omg2n =
          TriaxialLine::remx(fic.eE * (omg2a - fic.omg0).radians(),
                             Math::pi());
        int parity = fmod(omg2n.second, real(2)) ? -1 : 1;
        real deltax = Triaxial::clamp(fic.delta + omg2n.second * deltashift,
                                      2);
        v2 = fomg().fwd(omg2n.first);
        u2 = fbet().inv(fomg()(v2) + deltax);
        real bet2 = fic.nN * parity * fbet().rev(u2);
        bet2a = ang::radians(bet2);
        // !umbalt definition of alp0
        ang alp0x(fic.alp1.nearest(1U));
        alp2a = ang(fic.eE * _t.kp / cosh(v2),
                    _t.k * fic.nN * parity / cosh(u2)).rebase(alp0x);
        ii = int(omg2n.second);
      }
      ret.betw2 = u2;
      ret.omgw2 = v2;
      ret.ind2 = ii;
    } else {
      // gamma == NaN
    }
    omg2a += ang::cardinal(1);
    return ret;
  }

  TriaxialLineF::fics::fics()
    : bet1(ang::NaN())
    , omg1(ang::NaN())
    , alp1(ang::NaN())
    , bet0(0)
    , omg0(0)
    , alp0(0)
    , u0(Math::NaN())
    , v0(Math::NaN())
    , delta(Math::NaN())
    , nN(0)
    , eE(0)
  {}

  TriaxialLineF::fics::fics(const TriaxialLineF& f,
                            const Angle& bet10, const Angle& omg10,
                            const Angle& alp10)
    : bet1(bet10)
      // omg10 - 90
    , omg1(omg10 - ang::cardinal(1))
    , alp1(alp10)
  {
    const real eps = numeric_limits<real>::epsilon();
    const Triaxial::gamblk& gm = f.gm();
    const Triaxial& t = f.t();
    if (bet1.s() == 0 && fabs(alp1.c()) <= Math::sq(eps))
      alp1 = ang(alp1.s(), - Math::sq(eps), alp1.n(), true);
    eE = signbit(alp1.s()) ? -1 : 1;
    nN = signbit(alp1.c()) ? -1 : 1;
    if (gm.gam > 0) {
      bet0 = bet1.nearest(2U);
      alp0 = alp1.nearest(1U);
      psi1 = ang(t.k * bet1.s(),
                 bet0.c() * alp1.c() *
                 hypot(t.k * bet1.c(), t.kp * omg1.c()));
      v0 = f.fbet().fwd(psi1.radians());
      u0 = f.fomg().fwd(eE * omg1.radians());
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gm.gam < 0) {
      omg0 = omg1.nearest(2U);
      alp0 = alp1.nearest(2U);
      // Need Angle(0, 0) to be treated like Angle(0, 1) here.
      psi1 = ang(t.kp * omg1.s(),
                 omg0.c() * alp1.s() *
                 hypot(t.k * bet1.c(), t.kp * omg1.c()));
      v0 = f.fomg().fwd(psi1.radians());
      u0 = f.fbet().fwd(nN * bet1.radians());
      delta = f.fbet()(u0) - f.fomg()(v0);
    } else if (gm.gam == 0) {
      alp0 = alp1.nearest(t.umbalt ? 2U : 1U);
      // N.B. factor of k*kp omitted
      // bet0, omg0 are the middle of the initial umbilical segment
      if (fabs(bet1.c()) < 8*eps && fabs(omg1.c()) < 8*eps) {
        bet0 = bet1.nearest(1U) + ang::cardinal(nN);
        omg0 = omg1.nearest(1U) + ang::cardinal(eE);
        delta = f.deltashift/2 - log(fabs(alp1.t()));
      } else {
        bet0 = bet1.nearest(2U);
        omg0 = omg1.nearest(2U);
        delta = nN * f.fbet()(geod_fun::lamang(bet1 - bet0)) -
          eE * f.fomg()(geod_fun::lamang(omg1 - omg0));
      }
    } else {
      // gamma = NaN
    }
  }

  void TriaxialLineF::fics::setquadrant(const TriaxialLineF& f, unsigned q) {
    const real eps = numeric_limits<real>::epsilon();
    real gam = f.gm().gam;
    const Triaxial& t = f.t();
    alp1.setquadrant(q);
    if (bet1.s() == 0 && fabs(alp1.c()) <= Math::sq(eps))
      alp1 = ang(alp1.s(), -Math::sq(eps), alp1.n(), true);
    int oE = eE, oN = nN;
    eE = signbit(alp1.s()) ? -1 : 1;
    nN = signbit(alp1.c()) ? -1 : 1;
    if (gam > 0) {
      alp0 = alp1.nearest(1U);
      psi1.reflect(false, nN != oN);
      v0 = f.fbet().fwd(psi1.radians());
      u0 *= eE/oE;
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gam < 0) {
      alp0 = alp1.nearest(2U);
      psi1.reflect(false, eE != oE);
      v0 = f.fomg().fwd(psi1.radians());
      u0 *= nN/oN;
      delta = f.fbet()(u0) - f.fomg()(v0);
    } else if (gam == 0) {
      // Only expect to invoke setquadrant in this case
      alp0 = alp1.nearest(t.umbalt ? 2U : 1U);
      if (fabs(bet1.c()) < 8*eps && fabs(omg1.c()) < 8*eps)
        delta = f.deltashift/2 - log(fabs(alp1.t()));
      else
        delta = nN * f.fbet()(geod_fun::lamang(bet1 - bet0)) -
          eE * f.fomg()(geod_fun::lamang(omg1 - omg0));
    } else {
      // gamma = NaN
    }
  }

  TriaxialLineG::gics::gics()
    : sig1(Math::NaN())
  {}

  TriaxialLineG::gics::gics(const TriaxialLineG& g,
                          const TriaxialLineF::fics& fic)
  {
    if (g.gamma() > 0) {
      sig1 = g.gbet()(fic.v0) + g.gomg()(fic.u0);
    } else if (g.gamma() < 0) {
      sig1 = g.gbet()(fic.u0) + g.gomg()(fic.v0);
    } else if (g.gamma() == 0) {
      sig1 = fic.nN * g.gbet()(geod_fun::lamang(fic.bet1 - fic.bet0)) +
          fic.eE * g.gomg()(geod_fun::lamang(fic.omg1 - fic.omg0));
    } else {
      // gamma = NaN
    }
  }

  Math::real TriaxialLineG::dist(gics ic, TriaxialLineF::disttx d) const {
    real sig2 = gbet()(d.betw2) + gomg()(d.omgw2) + d.ind2 * 2*s0;
    return (sig2 - ic.sig1) * t().b;
  }
} // namespace GeographicLib

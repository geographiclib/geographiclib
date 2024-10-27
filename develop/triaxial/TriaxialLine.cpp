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

  TriaxialLine::TriaxialLine(fline f, fline::fics fic,
                             gline g, gline::gics gic)
    : _t(f.t())
    , _f(f)
    , _fic(fic)
    , _g(g)
    , _gic(gic)
  {}

  TriaxialLine::TriaxialLine(const Triaxial& t,
                             Angle bet1, Angle omg1, Angle alp1) {
    bet1.round();
    omg1.round();
    alp1.round();
    _t = t;
    Triaxial::gamblk gam = t.gamma(bet1, omg1, alp1);
    _f = fline(t, gam);
    _fic = fline::fics(_f, bet1, omg1, alp1);
    _g = gline(t, gam);
    _gic = gline::gics(_g, _fic);
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
    real sig2 = _gic.sig1 + s12/_t._b;
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
      alp2a = ang(_fic.eE * hypot(_t._k * _f.gm().nu, _t._kp * omg2a.c()),
                  _fic.bet0.c() * _t._k * _f.gm().nup * psi2.c())
        .rebase(_fic.alp0);
    } else if (_f.gamma() < 0) {
      real u2, v2;
      if (fomg().NCoeffsInv() <= fbet().NCoeffsInv())
        solve2( _fic.delta, sig2, fbet(), fomg(), gbet(), gomg(), u2, v2,
                countn, countb);
      else
        solve2(-_fic.delta, sig2, fomg(), fbet(), gomg(), gbet(), v2, u2,
               countn, countb);
      bet2 = _fic.nN * fbet().rev(u2);
      bet2a = ang::radians(bet2);
      ang psi2 = ang::radians(fomg().rev(v2));
      // Already normalized
      omg2a = ang(_f.gm().nup * psi2.s(),
                  _fic.omg0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                  0, true).rebase(_fic.omg0);
      alp2a = ang(_fic.omg0.c() * _t._kp * _f.gm().nup * psi2.c(),
                  _fic.nN * hypot(_t._kp * _f.gm().nu, _t._k * bet2a.c()))
        .rebase(_fic.alp0);
    } else if (_f.gamma() == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_g.s0);  // reduce to [-s0, s0)
      if (sig2n.first - _g.s0 >= -5 * numeric_limits<real>::epsilon()) {
        sig2n.first = -_g.s0;
        ++sig2n.second;
      }
      real u2, v2,
        deltax = clamp(_fic.delta + sig2n.second * _f.deltashift, 1);
      solve2u(deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2,
              countn, countb);
      bet2 = fbet().rev(u2); omg2 = fomg().rev(v2);
      bet2a = anglam(u2, _t._kp);
      omg2a = anglam(v2, _t._k);
      int parity = fmod(sig2n.second, real(2)) ? -1 : 1;
      if (_t._umbalt) {
        Ex = _fic.eE * parity;
        omg2a = omg2a.flipsign(Ex).rebase(_fic.omg0);
        bet2a += ang::cardinal(2 * sig2n.second);
        bet2a = bet2a.flipsign(_fic.nN) + _fic.bet0;
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = ang(_fic.nN * _t._kp * Ex / mcosh(v2, _t._k),
                    _t._k / mcosh(u2, _t._kp)).
          rebase(_fic.alp0);
      } else {
        Nx = _fic.nN * parity;
        omg2a += ang::cardinal(2 * sig2n.second);
        omg2a = omg2a.flipsign(_fic.eE) + _fic.omg0;
        bet2a = bet2a.reflect((_fic.bet0.c() * Nx) < 0, _fic.bet0.c() < 0)
          .rebase(_fic.bet0);
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = ang(_fic.eE * _t._kp / mcosh(v2, _t._k),
                    _t._k * Nx / mcosh(u2, _t._kp)).
          rebase(_fic.alp0);
      }
    } else {
      // gamma = NaN
    }
    bet2a.round();
    omg2a.round();
    alp2a.round();
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
                            const ffun& fx, const ffun& fy,
                            const gfun& gx, const gfun& gy,
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
                             const ffun& fbet, const ffun& fomg,
                             const gfun& gbet, const gfun& gomg,
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
               7 * numeric_limits<real>::epsilon()) {
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
                           const ffun& fx, const ffun& fy,
                           const gfun& gx, const gfun& gy,
                           real x0, real xa, real xb,
                           real xscale, real zscale,
                           real& x, real& y,
                           int* countn, int* countb) {
    // Find [x,y] s.t.
    //   fx(x) - fy(y) - f0 = 0
    //   gx(x) + gy(y) - g0 = 0
    // Assume fy.inv is known, then
    //   y = fy.inv(fx(x) - f0)
    // and we consider the 1d problem
    //   gx(x) + gy( fy.inv(fx(x) - f0) ) - g0 = 0
    // d/dx of the LHS is
    //   fx'(x) * ( gx'(x)/fx'(x) + gy'(y)/fy'(y) )
    //   = fx'(x) * (gfx'(x) + gfy'(y))
    // where gfx'(y) = gx'(x)/fx'(x)
    //       gfy'(y) = gy'(y)/fy'(y)
    //    cout << "BBX0\n";
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
    //    cout << "BBX1\n";
    y = fy.inv(fx(x) - f0);
    // Do one round of 2d Newton
    // Trial solution z0 = [x0, y0]'
    // let t0 = [fx(x0) - fy(y0) - f0, gx(x0) + gy(y0) - g0]'
    // Updated solution is
    // z1 = z0 - M . t0
    // M = 1/(gfx'(x) + gfy'(y) *
    //   [ gfy'(y)/fx'(x), 1/fx'(x)]
    //   [-gfx'(x)/fy'(y), 1/fy'(y)]
    int nfix = 1;
    for (int i = 0; i < nfix; ++i) {
      real tf = fx(x) - fy(y) - f0, tg = gx(x) + gy(y) - g0,
        fxp = fx.deriv(x), fyp = fy.deriv(y),
        gfxp = gx.gfderiv(x), gfyp =  gy.gfderiv(y),
        den = gfxp + gfyp;
      x -= ( gfyp * tf + tg) / (fxp * den);
      y -= (-gfxp * tf + tg) / (fyp * den);
      if (countn)
        ++*countn;
    }
  }

  void TriaxialLine::Hybrid(Angle bet2,
                            Angle& bet2a, Angle& omg2a, Angle& alp2a,
                            real& s12)
    const {
    fline::disttx d = _f.Hybrid(_fic, bet2, bet2a, omg2a, alp2a);
    s12 = _g.dist(_gic, d);
  }

  void TriaxialLine::Offset(real s13, bool reverse) {
    ang bet2, omg2, alp2;
    Position(s13, bet2, omg2, alp2);
    if (reverse) {
      alp2.reflect(true, true);
      // TODO: check if point 2 is an umbilical point
    }
    _fic = fline::fics(_f, bet2, omg2, alp2);
    _gic = gline::gics(_g, _fic);
  }

  void TriaxialLine::Optimize() {
    _f.ComputeInverse();
  }

  TriaxialLine::fline::disttx
  TriaxialLine::fline::Hybrid(const fics& fic,
                              Angle bet2,
                              Angle& bet2a, Angle& omg2a, Angle& alp2a)
    const {
    ang tau12;
    if (gamma() > 0) {
      real spsi = _t._k * bet2.s(),
        // In evaluating equivalent expressions, choose the one
        // with minimum cancelation
        cpsi = 0 + _t._k
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
      ang bet2b = bet2; bet2b.reflect(false, fic.nN < 0);
      ang bet12 = bet2b - fic.bet1;
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

  TriaxialLine::fline::fline(const Triaxial& t, Triaxial::gamblk gam)
    : _t(t)
    , _gm(gam)
    , _fbet(_t._k2 , _t._kp2,  _t._e2, -_gm.gamma, t)
    , _fomg(_t._kp2, _t._k2 , -_t._e2,  _gm.gamma, t)
    , _invp(false)
    {
      df = _gm.gamma == 0 ? _fbet.Max() - _fomg.Max() : 0;
      deltashift = _gm.gamma == 0 ? 2*df : 0;
    }

  void TriaxialLine::fline::ComputeInverse() {
    if (!_invp) {
      _fbet.ComputeInverse();
      _fomg.ComputeInverse();
      _invp = true;
    }
  }

  TriaxialLine::gline::gline(const Triaxial& t, const Triaxial::gamblk& gam)
    : _t(t)
    , _gm(gam)
    , _gbet(_t._k2 , _t._kp2,  _t._e2, -_gm.gamma, _t)
    , _gomg(_t._kp2, _t._k2 , -_t._e2,  _gm.gamma, _t)
    , s0(_gm.gamma == 0 ? _gbet.Max() + _gomg.Max() : 0)
  {}

  Math::real TriaxialLine::fline::Hybrid0(const fics& fic,
                                          Angle bet2, Angle omg2)
  const {
    bool debug = false;
    ang bet2a, omg2a, alp2a, omg2b(omg2);
    (void) Hybrid(fic, bet2, bet2a, omg2a, alp2a);
    (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
    if (debug)
      cout << "H0 "
           << real(fic.bet1) << " "
           << real(fic.omg1) << " "
           << real(fic.alp1) << " "
           << real(bet2a) << " " << real(omg2a) << " "
           << real(alp2a) << " " << real(omg2b) << "\n"
           << fic.omg1.s() << " " << fic.omg1.c() << " "
           << omg2a.s() << " " << omg2a.c() << " "
           << omg2b.s() << " " << omg2b.c() << " ";
    omg2a -= omg2b;
    if (debug)
      cout << omg2a.s() << " " << omg2a.c() << " " << omg2a.c() + 1 << "\n";
    // This test has been superceded by Fix #2 for triaxial sphere in
    // Triaxial.hpp
    if (false && _t._e2 < 128 * numeric_limits<real>::epsilon() &&
        ( 1 + omg2a.c() <= 1 - numeric_limits<real>::epsilon() &&
          fabs(omg2a.s()) <= 16 * numeric_limits<real>::epsilon() ) &&
        (signbit(omg2a.s()) ^ signbit(fic.alp1.s()))) {
      if (debug)
        cout << "HERE\n";
      // This is the only place where the limit e2 -> 0 is treated specially.
      // We're following a great circle path the opposite latitude and we need
      // to ensure that the sign of the longitude difference (+/-pi) matches
      // the sign of alp1.
      omg2a = ang(copysign(real(0), fic.alp1.s()), real(-1));
    }
    return omg2a.radians0();
  }

  TriaxialLine::fline::disttx
  TriaxialLine::fline::ArcPos0(const fics& fic, Angle tau12,
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
    if (_t._debug)
      cout << "AP " << real(tau12) << " "
           << real(omg2a) << " "
           << fic.eE * omg2a.radians() << " " << u2 << " "
           << fomg()(u2) + fic.delta << " " << v2 << "\n";
      }
      // Already normalized
      bet2a = ang(gm().nup * psi2.s(),
                  fic.bet0.c() * hypot(psi2.c(), gm().nu * psi2.s()),
                  0, true).rebase(fic.bet0);
      alp2a = ang(fic.eE * hypot(_t._k * gm().nu, _t._kp * omg2a.c()),
                  fic.bet0.c() * _t._k * gm().nup * psi2.c()).rebase(fic.alp0);
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
      alp2a = ang(fic.omg0.c() * _t._kp * gm().nup * psi2.c(),
                  fic.nN * hypot(_t._kp * gm().nu, _t._k * bet2a.c()))
        .rebase(fic.alp0);
      ret.betw2 = u2;
      ret.omgw2 = v2;
    } else if (gamma() == 0) {
      real u2, v2;
      int ii;
      if (betp) {
        bet2a = fic.bet1 + tau12.flipsign(fic.nN);
        pair<real, real> bet2n =
          remx(fic.nN * (bet2a - fic.bet0).radians(), Math::pi());
        int parity = fmod(bet2n.second, real(2)) ? -1 : 1;
        real deltax = clamp(fic.delta + bet2n.second * deltashift, 2);
        u2 = fbet().fwd(bet2n.first);
        v2 = fomg().inv(fbet()(u2) - deltax);
        omg2a = ang::radians(fic.eE * parity * fomg().rev(v2))
          .rebase(fic.omg0);
        // umbalt definition of alp0
        ang alp0x(fic.alp1.nearest(2U));
        alp2a = ang(fic.nN * _t._kp * fic.eE *
                    parity / mcosh(v2, _t._k),
                    _t._k / mcosh(u2, _t._kp)).rebase(alp0x);
        ii = int(bet2n.second);
        // Move forward from umbilical point
        bet2a += ang::eps();
      } else {
        omg2a = fic.omg1 + tau12.flipsign(fic.eE);
        pair<real, real> omg2n =
          remx(fic.eE * (omg2a - fic.omg0).radians(), Math::pi());
        int parity = fmod(omg2n.second, real(2)) ? -1 : 1;
        real deltax = clamp(fic.delta + omg2n.second * deltashift, 2);
        v2 = fomg().fwd(omg2n.first);
        u2 = fbet().inv(fomg()(v2) + deltax);
        real bet2 = fic.nN * parity * fbet().rev(u2);
        bet2a = ang::radians(bet2);
        // !umbalt definition of alp0
        ang alp0x(fic.alp1.nearest(1U));
        alp2a = ang(fic.eE * _t._kp / mcosh(v2, _t._k),
                    _t._k * fic.nN *
                    parity / mcosh(u2, _t._kp)).rebase(alp0x);
        ii = int(omg2n.second);
        // Move forward from umbilical point
        omg2a += ang::eps();
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

  TriaxialLine::fline::fics::fics(const fline& f,
                                  Angle bet10, Angle omg10,
                                  Angle alp10)
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
    if (gm.gamma > 0) {
      bet0 = bet1.nearest(2U);
      alp0 = alp1.nearest(1U);
      psi1 = ang(t._k * bet1.s(),
                 bet0.c() * alp1.c() *
                 hypot(t._k * bet1.c(), t._kp * omg1.c()));
      v0 = f.fbet().fwd(psi1.radians());
      u0 = f.fomg().fwd(eE * omg1.radians());
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gm.gamma < 0) {
      omg0 = omg1.nearest(2U);
      alp0 = alp1.nearest(2U);
      // Need Angle(0, 0) to be treated like Angle(0, 1) here.
      psi1 = ang(t._kp * omg1.s(),
                 omg0.c() * alp1.s() *
                 hypot(t._k * bet1.c(), t._kp * omg1.c()));
      v0 = f.fomg().fwd(psi1.radians());
      u0 = f.fbet().fwd(nN * bet1.radians());
      delta = f.fbet()(u0) - f.fomg()(v0);
    } else if (gm.gamma == 0) {
      alp0 = alp1.nearest(t._umbalt ? 2U : 1U);
      // N.B. factor of k*kp omitted
      // bet0, omg0 are the middle of the initial umbilical segment
      if (fabs(bet1.c()) < 8*eps && fabs(omg1.c()) < 8*eps) {
        bet0 = bet1.nearest(1U) + ang::cardinal(nN);
        omg0 = omg1.nearest(1U) + ang::cardinal(eE);
        delta = f.deltashift/2 - log(fabs(alp1.t()));
      } else {
        bet0 = bet1.nearest(2U);
        omg0 = omg1.nearest(2U);
        delta = nN * f.fbet()(lamang(bet1 - bet0, t._kp)) -
          eE * f.fomg()(lamang(omg1 - omg0, t._k));
      }
    } else {
      // gamma = NaN
    }
  }

  void TriaxialLine::fline::fics::setquadrant(const fline& f, unsigned q) {
    const real eps = numeric_limits<real>::epsilon();
    real gam = f.gm().gamma;
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
      alp0 = alp1.nearest(t._umbalt ? 2U : 1U);
      if (fabs(bet1.c()) < 8*eps && fabs(omg1.c()) < 8*eps)
        delta = f.deltashift/2 - log(fabs(alp1.t()));
      else
        delta = nN * f.fbet()(lamang(bet1 - bet0, t._kp)) -
          eE * f.fomg()(lamang(omg1 - omg0, t._k));
    } else {
      // gamma = NaN
    }
  }

  TriaxialLine::gline::gics::gics(const gline& g, const fline::fics& fic) {
    if (g.gamma() > 0) {
      sig1 = g.gbet()(fic.v0) + g.gomg()(fic.u0);
    } else if (g.gamma() < 0) {
      sig1 = g.gbet()(fic.u0) + g.gomg()(fic.v0);
    } else if (g.gamma() == 0) {
      const Triaxial& t = g.t();
      sig1 = fic.nN *
        g.gbet()(lamang(fic.bet1 - fic.bet0, t._kp)) +
        fic.eE * g.gomg()(lamang(fic.omg1 - fic.omg0, t._k));
    } else {
      // gamma = NaN
    }
  }

  Math::real TriaxialLine::gline::dist(gics ic, fline::disttx d) const {
    real sig2 = gbet()(d.betw2) + gomg()(d.omgw2) + d.ind2 * 2*s0;
    return (sig2 - ic.sig1) * t()._b;
  }

  TriaxialLine::ffun::ffun(real kap, real kapp, real eps, real mu,
                           const Triaxial& t)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _sqrtkapp(sqrt(_kapp))
    // Restrict oblpro treatment to mu == 0
    , _oblpro(t._oblpro && _mu == 0)
    , _invp(false)
  {
    // mu in [-kap, kapp], eps in (-inf, 1/kap)
    if (_oblpro && _kapp == 0 && _mu <= 0) { // mu > 0 not allowed
      // _kap == 1
      // Include scale = 1 in TrigfunExt constructor because this function gets
      // added atan(sqrt(-mu) * tan(phi)).
      _tx = false;
      // f multiplied by sqrt(-mu)
      _fun = TrigfunExt(
                        [kap = _kap, kapp = _kapp,
                         eps = _eps, mu = _mu]
                        (real psi) -> real
                        { return dfpsioblp(sin(psi), cos(psi), eps, mu); },
                        Math::pi()/2, false, 1);
    } else if (_oblpro && _kap == 0 && _mu >= 0) { // mu < 0 not allowed
      // _kapp == 1
      _tx = false;
      // f multiplied by sqrt(mu)
      _fun = TrigfunExt(
                        [eps = _eps, mu = _mu]
                        (real phi) -> real
                        { return fphioblp(phi, eps, mu); },
                        Math::pi()/2, false);
    } else if (_mu > 0) {
      _tx = _mu / (_kap + _mu) < t._ellipthresh;
      if (_tx) {
        _ell = EllipticFunction(_kap / (_kap + _mu), 0, _mu / (_kap + _mu), 1);
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp,
                           eps = _eps, mu = _mu, ell = _ell]
                          (real u) -> real
                          { real sn, cn, dn; (void) ell.am(u, sn, cn, dn);
                            return fup(cn, kap, kapp, eps, mu); },
                          _ell.K(), false);
      } else
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                          (real phi) -> real
                          { return fphip(cos(phi), kap, kapp, eps, mu); },
                          Math::pi()/2, false);
    } else if (_mu < 0) {
      _tx = -_mu / _kap < t._ellipthresh;
      if (_tx) {
        _ell = EllipticFunction((_kap + _mu) / _kap, 0, -_mu / _kap, 1);
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp,
                           eps = _eps, mu = _mu, ell = _ell]
                          (real v) -> real
                          { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                            return fvp(dn, kap, kapp, eps, mu); },
                          _ell.K(), false);
      } else
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return fpsip(sin(psi), cos(psi),
                                         kap, kapp, eps, mu); },
                          Math::pi()/2, false);
    } else if (_mu == 0) {
      _tx = _kapp < t._ellipthresh;
      // f multiplied by sqrt(kap*kapp)
      // Include scale = 1 in TrigfunExt constructor because this function gets
      // added to u.
      if (_tx) {
        _ell = EllipticFunction(_kap, 0, _kapp, 1);
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp,
                           eps = _eps, ell = _ell]
                          (real v) -> real
                          { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                            return dfvp(cn, dn, kap, kapp, eps); },
                          2 * _ell.K(), true, 1);
      } else
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp, eps = _eps]
                          (real phi) -> real
                          { return dfp(cos(phi), kap, kapp, eps); },
                          Math::pi(), true, 1);
    } else {
      // _mu == NaN
      _tx = false;
    }
    _nmax = _fun.NCoeffs();
    // N.B. _max < 0 for _mu == 0 && eps < 0
    _max = _mu == 0 ?
      _fun(_tx ? _ell.K() : Math::pi()/2) : _fun.Max();
  }

  void TriaxialLine::ffun::ComputeInverse() {
    if (!_invp) {
      if (_mu == 0) {
        _countn = _countb = 0;
        // Include scale = 1 in TrigfunExt constructor because _dfinv gets
        // added to u.
        _dfinv = Trigfun(
                         [this]
                         (real phi, real x0) -> real
                         { real u = lam(phi, _sqrtkapp);
                           return
                             root(u, lam(x0, _sqrtkapp), &_countn, &_countb,
                                  sqrt(numeric_limits<real>::epsilon()))
                             - u; },
                         true, true, Math::pi(),
                         int(ceil(real(1.5) * NCoeffs())),
                         sqrt(numeric_limits<real>::epsilon()), 1);
        /*
        cout << "HERE " << _kapp << " " << NCoeffs() << " "
             << _dfinv.NCoeffs() << " "
             << _countn << " " << _countb << "\n";
        _invp = true;
        _countn = 0; _countb = 0;
        for (int i = 0; i <= 20; ++i) {
          real u = i/real(5), z = (*this)(u),
            u1 = inv2(z, &_countn, &_countb);
          cout << u << " " << z << " " << u1-u << "\n";
        }
        cout << _countn << " " << _countb << "\n";
        */
      } else
        _fun.ComputeInverse();
    }
    _invp = true;
  }

  Math::real TriaxialLine::ffun::root(real z, real x0,
                                      int* countn, int* countb,
                                      real tol) const {
    if (_mu != 0) return Math::NaN();
    if (!isfinite(z)) return z; // Deals with +/-inf and nan
    real d = fabs(Max())
      + 2 * numeric_limits<real>::epsilon() * fmax(real(1), fabs(z)),
      xa = z - d,
      xb = z + d;
    x0 = fmin(xb, fmax(xa, x0));
    // Solve z = u - _fun(_tx ? _ell.F(gd(u)) : gd(u)) for u
    // N.B. use default tol for root, because we want accurate answers here
    //    cout << "BBA\n";
    return Trigfun::root(
                         [this]
                         (real u) -> pair<real, real>
                         { return pair<real, real>((*this)(u), deriv(u)); },
                         z,
                         x0, xa, xb,
                         Math::pi()/2, Math::pi()/2, 1, countn, countb, tol);
  }

  // Approximate inverse using _chiinv or _fun.inv0
  Math::real TriaxialLine::ffun::inv0(real z) const {
    if (!_invp) return Math::NaN();
    return _mu == 0 ?
      z + _dfinv(gd(z, _sqrtkapp)) : _fun.inv0(z);
  }

  // Accurate inverse by direct Newton (not using _finv)
  Math::real TriaxialLine::ffun::inv1(real z, int* countn, int* countb) const {
    //    cout << "BBZ\n";
    return _mu == 0 ? root(z, z, countn, countb) :
      _fun.inv1(z, countn, countb);
  }

  // Accurate inverse correcting result from _finv
  Math::real TriaxialLine::ffun::inv2(real z, int* countn, int* countb) const {
    if (!_invp) return Math::NaN();
    return _mu == 0 ? root(z, inv0(z), countn, countb) :
      _fun.inv2(z, countn, countb);
  }

  TriaxialLine::gfun::gfun(real kap, real kapp, real eps, real mu,
                           const Triaxial& t)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _sqrtkapp(sqrt(_kapp))
    // Restrict oblpro treatment to mu == 0
    , _oblpro(t._oblpro && _mu == 0)
    , _invp(false)
  {
    // mu in [-kap, kapp], eps in (-inf, 1/kap)
    // Restrict oblpro treatment to mu == 0
    if (_oblpro && _kapp == 0 && _mu <= 0) { // mu > 0 not allowed
      // _kap == 1
      _tx = false;
      _fun = TrigfunExt(
                        [kap = _kap, kapp = _kapp,
                         eps = _eps, mu = _mu]
                        (real psi) -> real
                        { return gpsioblp(sin(psi), cos(psi), eps, mu); },
                        Math::pi()/2, false);
    } else if (_oblpro && _kap == 0 && _mu >= 0) { // mu < 0 not allowed
      // _kapp == 1
      _tx = false;
      _fun = TrigfunExt(
                        [eps = _eps, mu = _mu]
                        (real phi) -> real
                        { return gphioblp(phi, eps, mu); },
                        Math::pi()/2, false);
    } else if (_mu > 0) {
      _tx = _mu / (_kap + _mu) < t._ellipthresh;
      if (_tx) {
        _ell = EllipticFunction(_kap / (_kap + _mu), 0, _mu / (_kap + _mu), 1);
        _fun =TrigfunExt(
                         [kap = _kap, kapp = _kapp,
                          eps = _eps, mu = _mu, ell = _ell]
                         (real u) -> real
                         { real sn, cn, dn; (void) ell.am(u, sn, cn, dn);
                           return gup(cn, dn, kap, kapp, eps, mu); },
                         _ell.K());
      } else
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                          (real phi) -> real
                          { return gphip(cos(phi), kap, kapp, eps, mu); },
                          Math::pi()/2);
    } else if (_mu < 0) {
      _tx = -_mu / _kap < t._ellipthresh;
      if (_tx) {
        _ell = EllipticFunction((_kap + _mu) / _kap, 0, -_mu / _kap, 1);
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp,
                           eps = _eps, mu = _mu, ell = _ell]
                          (real v) -> real
                          { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                            return gvp(cn, dn, kap, kapp, eps, mu); },
                          _ell.K());
      } else
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return gpsip(sin(psi), cos(psi),
                                            kap, kapp, eps, mu); },
                          Math::pi()/2);
    } else if (_mu == 0) {
      _tx = _kapp < t._ellipthresh;
      if (_tx) {
        _ell = EllipticFunction(_kap, 0, _kapp, 1);
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp, eps = _eps, ell = _ell]
                          (real v) -> real
                          { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                            return g0vp(cn, kap, kapp, eps); },
                          2*_ell.K(), true);
      } else
        _fun = TrigfunExt(
                          [kap = _kap, kapp = _kapp, eps = _eps, ell = _ell]
                          (real phi) -> real
                          { return g0p(cos(phi), kap, kapp, eps); },
                          Math::pi(), true);
    } else {
      // _mu == NaN
      _tx = false;
    }
    _max = _mu == 0 ? _fun(_tx ? _ell.K() : Math::pi()/2) :
      _fun.Max();
  }

  Math::real TriaxialLine::gfun::gfderiv(real u) const {
    // return g'(u)/f'(u)
    real sn = 0, cn = 0, dn = 0;
    if (_oblpro && _kapp == 0 && _mu <= 0) // mu > 0 not allowed
      return gfpsioblp(sin(u), cos(u), _mu);
    else if (_oblpro && _kap == 0 && _mu >= 0) // mu < 0 not allowed
      return gfphioblp(u, _mu);
    else if (_mu != 0 && _tx)
      (void) _ell.am(u, sn, cn, dn);
    if (_mu > 0)
      return _tx ? gfup(cn, _kap, _mu) : gfphip(cos(u), _kap, _mu);
    else if (_mu < 0)
      return _tx ? gfvp(dn, _kap, _mu) :
        gfpsip(sin(u), cos(u), _kap, _mu);
    else                      // _mu == 0
      // This includes factor of sqrt(kap * kapp) because of adjustment of
      // definition of f for mu == 0.
      return gf0up(u, _kap, _kapp);
  }

  Math::real TriaxialLine::gfun::root(real z, real x0,
                                      int* countn, int* countb,
                                      real tol) const {
    if (_mu != 0) return Math::NaN();
    if (!isfinite(z)) return z; // Deals with +/-inf and nan
    real d = fabs(Max())
      + 2 * numeric_limits<real>::epsilon() * fmax(real(1), fabs(z)),
      xa = z - d,
      xb = z + d;
    x0 = fmin(xb, fmax(xa, x0));
    // Solve z = u - _fun(_tx ? _ell.F(gd(u)) : gd(u)) for u
    // N.B. use default tol for root, because we want accurate answers here
    //    cout << "BBD\n";
    return Trigfun::root(
                         [this]
                         (real u) -> pair<real, real>
                         { return pair<real, real>((*this)(u), deriv(u)); },
                         z,
                         x0, xa, xb,
                         Math::pi()/2, Math::pi()/2, 1, countn, countb, tol);
  }

  // Approximate inverse using _chiinv or _fun.inv0
  Math::real TriaxialLine::gfun::inv0(real z) const {
    if (!_invp) return Math::NaN();
    return _mu == 0 ? Math::NaN() : _fun.inv0(z);
  }

  // Accurate inverse by direct Newton (not using _finv)
  Math::real TriaxialLine::gfun::inv1(real z, int* countn, int* countb) const {
    return _mu == 0 ? Math::NaN() : _fun.inv1(z, countn, countb);
  }

  // Accurate inverse correcting result from _finv
  Math::real TriaxialLine::gfun::inv2(real z, int* countn, int* countb) const {
    if (!_invp) return Math::NaN();
    return _mu == 0 ? root(z, inv0(z), countn, countb) :
      _fun.inv2(z, countn, countb);
  }
  void TriaxialLine::gfun::ComputeInverse() {}

  // _mu > 0 && !_tx
  Math::real TriaxialLine::ffun::fphip(real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }
  // This is non-negative
  Math::real TriaxialLine::gfun::gphip(real c, real kap, real kapp,
                                          real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return c2 * sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }
  Math::real TriaxialLine::gfun::gfphip(real c, real kap, real /* mu */) {
    real c2 = kap * Math::sq(c);
    return c2;
  }

  // _mu > 0 && _tx
  Math::real TriaxialLine::ffun::fup(real cn, real kap, real kapp,
                                     real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }
  // This is non-negative
  Math::real TriaxialLine::gfun::gup(real cn, real /* dn */,
                                     real kap, real kapp,
                                     real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return c2 * sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }
  Math::real TriaxialLine::gfun::gfup(real cn, real kap, real /* mu */) {
    real c2 = kap * Math::sq(cn);
    return c2;
  }

  // _mu == 0 && !_tx
  Math::real TriaxialLine::ffun::dfp(real c,
                                     real kap, real kapp, real eps) {
    // function dfp = dfpf(phi, kappa, epsilon)
    // return derivative of sqrt(kap * kapp) * Delta f
    // s = sqrt(1 - kap * sin(phi)^2)
    real c2 = kap * Math::sq(c), s = sqrt(kapp + c2);
    return eps*kap * sqrt(kapp) * c / (s * (1 + sqrt(1 - eps*c2)));
  }
  Math::real TriaxialLine::gfun::g0p(real c, real kap, real kapp, real eps) {
    real c2 = kap * Math::sq(c);
    return sqrt( kap * (1 - eps * c2) / (kapp + c2) ) * c;
  }

  // _mu == 0 && _tx
  Math::real TriaxialLine::ffun::dfvp(real cn, real /* dn */,
                                      real kap, real kapp, real eps) {
    // function dfvp = dfvpf(v, kap, eps)
    // return derivative of sqrt(kap * kapp) * Delta f_v
    return eps*kap * sqrt(kapp) * cn /
      (1  + sqrt(1 - eps*kap * Math::sq(cn)));
  }
  Math::real TriaxialLine::gfun::g0vp(real cn, real kap, real /* kapp */,
                                      real eps) {
    real c2 = kap * Math::sq(cn);
    return sqrt( kap * (1 - eps * c2) ) * cn;
  }

  // _mu == 0 (_tx ignored)
  Math::real TriaxialLine::gfun::gf0up(real u, real kap, real kapp) {
    // Subst tan(phi) = sinh(u) /sqrt(kapp) in
    // kap * cos(phi)^2 gives kap*kapp/(kapp + sinh(u)^2)
    // Divide by sqrt(kap * kappp) to account of factor removed from f
    // functions.
    return sqrt(kap * kapp) / ( kapp + Math::sq(sinh(u)) );
  }

  // _mu < 0 && !_tx
  Math::real TriaxialLine::ffun::fpsip(real s, real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * c2) ) ;
  }
  // This is positive
  Math::real TriaxialLine::gfun::gpsip(real s, real c, real kap, real kapp,
                                          real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt(c2 * (1 - eps * c2) / (kapp + c2)) ;
  }
  Math::real TriaxialLine::gfun::gfpsip(real s, real c, real kap, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return c2;
  }

  // _mu < 0 && _tx
  Math::real TriaxialLine::ffun::fvp(real dn, real kap, real kapp,
                                     real eps, real /* mu */) {
    real c2 = kap * Math::sq(dn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * kap) );
  }
  // This is positive
  Math::real TriaxialLine::gfun::gvp(real /* cn */, real dn,
                                     real kap, real kapp,
                                     real eps, real /* mu */) {
    real dn2 = Math::sq(dn), c2 = kap * dn2;
    return dn2 * sqrt( kap * (1 - eps * c2) / (kapp + c2) );
  }
  Math::real TriaxialLine::gfun::gfvp(real dn, real kap, real /* mu */) {
    real dn2 = Math::sq(dn), c2 = kap * dn2;
    return c2;
  }

  // oblate/prolate variants for kap = 0, kapp = 1, mu > 0
  Math::real TriaxialLine::ffun::fphioblp(real /* phi */, real /* eps */,
                                           real /* mu */) {
    // Multiply by f functions by sqrt(abs(mu))
    // return 1 / sqrt(mu);
    return 1;
  }
  Math::real TriaxialLine::gfun::gphioblp(real /* phi */, real /* eps */,
                                          real /* mu */) {
    return 0;
  }
  Math::real TriaxialLine::gfun::gfphioblp(real /* phi */, real /* mu */) {
    return 0;
  }

  // oblate/prolate variants for kap = 1, kapp = 0, mu <= 0
  Math::real TriaxialLine::ffun::dfpsioblp(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // Multiply by f functions by sqrt(abs(mu))
    return sqrt(-mu) * eps / (1 + sqrt(1 - eps * c2));
  }
  Math::real TriaxialLine::gfun::gpsioblp(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    return sqrt(1 - eps * c2);
  }
  Math::real TriaxialLine::gfun::gfpsioblp(real s, real c, real mu) {
    // cos(phi)^2 = cos(psi)^2 - mu *sin(psi)^2
    // f' = sqrt(1-eps*cos(phi)^2)/cos(phi)^2
    // g' = sqrt(1-eps*cos(phi)^2)
    // g'/f' = cos(phi)^2
    return gfpsip(s, c, 1, mu) / sqrt(-mu);
  }

} // namespace GeographicLib

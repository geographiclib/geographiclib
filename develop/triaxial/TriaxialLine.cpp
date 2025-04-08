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
#include <iomanip>

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
    _ficx = fline::ficsx(_f, bet1, omg1, alp1);
    _g = gline(t, gam);
    _gic = gline::gics(_g, _fic);
    _gicx = gline::gicsx(_g, _fic);
    if (0) {
      cout << "FIC " << real(_fic.psi1) << " "
           << _fic.u0 << " " << _fic.v0 << " " << _fic.delta << " "
           << _fic.N << " " << _fic.E << "\n";
      cout << "GIC " << _gic.sig1 << " " << _gic.s13 << "\n";
    }
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
    if (_t._combine)
      return Positionx(s12, bet2a, omg2a, alp2a, countn, countb);
    // Compute points at distance s12
    real sig2 = _gic.sig1 + s12/_t._b;
    real bet2, omg2;
    int Ex, Nx;
    if (0)
      cout << "POS " << gamma() << " "
           << fbet().NCoeffs() << " " << fomg().NCoeffs() << "\n";
    if (gamma() > 0 || _t._kp2 == 0) {
      real u2, v2;
      if (_t._kp2 == 0 && (_t._oblpro || gamma() == 0)) {
        // If oblate, treat via biaxial machinery for
        // _t._oblpro: all cases
        // !_t._oblpro: meridional only

        // gomg()(x) == 0, so the g equation becomes gbet()(v2) = sig2
        v2 = gbet().inv(sig2);
        // fomg().inv(x) is just x, but keep it general
        u2 = fomg().inv(fbet()(v2) - _fic.delta);
      } else {
        // The general triaxial machinery.  If !_t._oblpro, this is used for
        // non-meridional geodesics on an oblate ellipsoid.
        if (fbet().NCoeffs() <= fomg().NCoeffs())
          solve2(-_fic.delta, sig2, fomg(), fbet(), gomg(), gbet(), u2, v2,
                 countn, countb);
        else
          solve2(_fic.delta, sig2, fbet(), fomg(), gbet(), gomg(), v2, u2,
                 countn, countb);
      }
      omg2 = _fic.E * fomg().rev(u2);
      omg2a = ang::radians(omg2);
      ang psi2 = ang::radians(fbet().rev(v2));
      // Already normalized
      bet2a = ang(_f.gm().nup * psi2.s(),
                  _fic.bet0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                  0, true).rebase(_fic.bet0);
      alp2a = ang(_fic.E * hypot(_t._k * _f.gm().nu, _t._kp * omg2a.c()),
                  _fic.bet0.c() * _t._k * _f.gm().nup * psi2.c())
        .rebase(_fic.alp0);
    } else if (gamma() < 0 || _t._k2 == 0) {
      real u2, v2;
      if (_t._k2 == 0 && (_t._oblpro || gamma() == 0)) {
        // If prolate, treat via biaxial machinery for
        // _t._oblpro: all cases
        // !_t._oblpro: meridional only

        // gbet()(x) == 0, so the g equation becomes gomg()(v2) = sig2
        v2 = gomg().inv(sig2);
        // fbet().inv(x) is just x, but keep it general
        u2 = fbet().inv(fomg()(v2) + _fic.delta);
      } else {
        // The general triaxial machinery.  If !_t._oblpro, this is used for
        // non-meridional geodesics on a prolate ellipsoid.
        if (fomg().NCoeffs() <= fbet().NCoeffs())
          solve2( _fic.delta, sig2, fbet(), fomg(), gbet(), gomg(), u2, v2,
                  countn, countb);
        else
          solve2(-_fic.delta, sig2, fomg(), fbet(), gomg(), gbet(), v2, u2,
                 countn, countb);
      }
      bet2 = _fic.N * fbet().rev(u2);
      bet2a = ang::radians(bet2);
      ang psi2 = ang::radians(fomg().rev(v2));
      // Already normalized
      omg2a = ang(_f.gm().nup * psi2.s(),
                  _fic.omg0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                  0, true).rebase(_fic.omg0);
      alp2a = ang(_fic.omg0.c() * _t._kp * _f.gm().nup * psi2.c(),
                  _fic.N * hypot(_t._kp * _f.gm().nu, _t._k * bet2a.c()))
        .rebase(_fic.alp0);
    } else if (gamma() == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_g.s0);  // reduce to [-s0, s0)
      if (sig2n.first - _g.s0 >= -5 * numeric_limits<real>::epsilon()) {
        sig2n.first = -_g.s0;
        ++sig2n.second;
      }
      real u2, v2,
        deltax = clamp(_fic.delta + sig2n.second * _f.deltashift(), 1);
      solve2u(deltax, sig2n.first, fbet(), fomg(), gbet(), gomg(), u2, v2,
              countn, countb);
      bet2 = fbet().rev(u2); omg2 = fomg().rev(v2);
      bet2a = anglam(u2, _t._kp);
      omg2a = anglam(v2, _t._k);
      int parity = fmod(sig2n.second, real(2)) != 0 ? -1 : 1;
      if (signbit(gamma())) {
        // if t._k2 == 0 then meridional prolate
        Ex = _fic.E * parity;
        omg2a = omg2a.flipsign(Ex).rebase(_fic.omg0);
        bet2a += ang::cardinal(2 * sig2n.second);
        bet2a = bet2a.flipsign(_fic.N) + _fic.bet0;
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = ang(_fic.N * _t._kp * Ex / mcosh(v2, _t._k),
                    _t._k / mcosh(u2, _t._kp)).
          rebase(_fic.alp0);
      } else {
        // if t._kp2 == 0 then meridional oblate
        Nx = _fic.N * parity;
        omg2a += ang::cardinal(2 * sig2n.second);
        omg2a = omg2a.flipsign(_fic.E) + _fic.omg0;
        bet2a = bet2a.reflect((_fic.bet0.c() * Nx) < 0, _fic.bet0.c() < 0)
          .rebase(_fic.bet0);
        // replace cos(bet)/cos(omg) by sech(u)/sech(v)
        alp2a = ang(_fic.E * _t._kp / mcosh(v2, _t._k),
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

  void TriaxialLine::Positionx(real s12,
                               Angle& bet2a, Angle& omg2a, Angle& alp2a,
                               int* countn, int* countb)
  const {
    // Buggy case...
    // echo 90 41 131.136116008822007127 4.28521283894037481354 |
    //   ./Geod3Solve -e 1 3 0 1 --combine --oblpro
    // => -39.36419427 154.00000000 169.90799437
    // latitude should be -90

    // Compute points at distance s12
    real sig2 = _gicx.sig1 + s12/_t._b;
    Angle &phi2a = bet2a, &tht2a = omg2a;
    if (_f.gammax() > 0 || _f.kxp2() == 0) {
      real u2, v2;
      if (_f.kxp2() == 0 && (_t._oblpro || _f.gammax() == 0)) {
        // If oblate, treat via biaxial machinery for
        // _t._oblpro: all cases
        // !_t._oblpro: meridional only

        // gtht()(x) == 0, so the g equation becomes gpsi()(v2) = sig2
        v2 = gpsi().inv(sig2);
        // ftht().inv(x) is just x, but keep it general
        u2 = ftht().inv(fpsi()(v2) - _ficx.delta);
      } else {
        // The general triaxial machinery.  If !_t._oblpro, this is used for
        // non-meridional geodesics on an oblate ellipsoid.
        if (fpsi().NCoeffs() <= ftht().NCoeffs())
          solve2(-_ficx.delta, sig2, ftht(), fpsi(), gtht(), gpsi(), u2, v2,
                 countn, countb);
        else
          solve2(_ficx.delta, sig2, fpsi(), ftht(), gpsi(), gtht(), v2, u2,
                 countn, countb);
      }
      tht2a = ang::radians(_ficx.E * ftht().rev(u2));
      ang psi2 = ang::radians(fpsi().rev(v2));
      // Already normalized
      phi2a = ang(_f.gm().nup * psi2.s(),
                  _ficx.phi0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                  0, true).rebase(_ficx.phi0);
      alp2a = ang(_ficx.E * hypot(_f.kx() * _f.gm().nu, _f.kxp() * tht2a.c()),
                  _ficx.phi0.c() * _f.kx() * _f.gm().nup * psi2.c())
        .rebase(_ficx.alp0);
    } else if (_f.gammax() == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_g.s0);  // reduce to [-s0, s0)
      if (sig2n.first - _g.s0 >= -5 * numeric_limits<real>::epsilon()) {
        sig2n.first = -_g.s0;
        ++sig2n.second;
      }
      real u2, v2,
        deltax = clamp(_ficx.delta + sig2n.second * _f.deltashift(), 1);
      solve2u(deltax, sig2n.first, fpsi(), ftht(), gpsi(), gtht(), u2, v2,
              countn, countb);
      // phi2 = fpsi().rev(u2); tht2 = ftht().rev(v2);
      phi2a = anglam(u2, _f.kxp());
      tht2a = anglam(v2, _f.kx());
      int parity = fmod(sig2n.second, real(2)) != 0 ? -1 : 1;
      // if t._kp2 == 0 then meridional oblate
      int Nx = _ficx.N * parity;
      tht2a += ang::cardinal(2 * sig2n.second);
      tht2a = tht2a.flipsign(_ficx.E) + _ficx.tht0;
      phi2a = phi2a.reflect((_ficx.phi0.c() * Nx) < 0, _ficx.phi0.c() < 0)
        .rebase(_ficx.phi0);
      // replace cos(phi)/cos(tht) by sech(u)/sech(v)
      alp2a = ang(_ficx.E * _f.kxp() / mcosh(v2, _f.kx()),
                  _f.kx() * Nx / mcosh(u2, _f.kxp())).
        rebase(_ficx.alp0);
    } else {
      // gamma = NaN
    }
    phi2a.round();
    tht2a.round();
    alp2a.round();
    if (_f.transpolar()) {
      swap(bet2a, omg2a);
      alp2a.reflect(false, false, true);
    }
    omg2a += ang::cardinal(1);
  }

  void TriaxialLine::Positionx(real s12, real& bet2, real& omg2, real& alp2,
                               bool unroll,
                               int* countn, int* countb) const {
    ang bet2a, omg2a, alp2a;
    Positionx(s12, bet2a, omg2a, alp2a, countn, countb);
    if (!unroll) {
      (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
      bet2a.setn(); omg2a.setn(); alp2a.setn();
    }
    bet2 = real(bet2a);
    omg2 = real(omg2a);
    alp2 = real(alp2a);
  }

  void TriaxialLine::solve2(real f0, real g0,
                            const hfun& fx, const hfun& fy,
                            const hfun& gx, const hfun& gy,
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
    if (0) {
      cout << "FEQ " << fx(x) << " " << fy(y) << " " << f0 << " "
           << fx(x) - fy(y) - f0 << "\n"
           << "GEQ " << gx(x) << " " << gy(y) << " " << g0 << " "
           << gx(x) + gy(y) - g0 << "\n";
      cout << "FF "
           << fx.HalfPeriod() << " " << fx.Slope() << " " << fx.Max() << " "
           << fy.HalfPeriod() << " " << fy.Slope() << " " << fy.Max() << "\n";
      cout << "GG "
           << gx.HalfPeriod() << " " << gx.Slope() << " " << gx.Max() << " "
           << gy.HalfPeriod() << " " << gy.Slope() << " " << gy.Max() << "\n";
    }
  }

  void TriaxialLine::solve2u(real d0, real s0,
                             const hfun& fbet, const hfun& fomg,
                             const hfun& gbet, const hfun& gomg,
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
                           const hfun& fx, const hfun& fy,
                           const hfun& gx, const hfun& gy,
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
    if (0) {
      auto fun = [&fx, &fy, &gx, &gy, f0]
        (real x) -> pair<real, real>
        { real y = fy.inv(fx(x) - f0);
          return pair<real, real>(gx(x) + gy(y),
                                  fx.deriv(x) *
                                  (gx.gfderiv(x) +
                                   gy.gfderiv(y))); };
      int num = 100;
      real dx = (xb - xa) / num;
      for (int i = 0; i <= num; ++i) {
        real x1 = xa + i * dx;
        auto p = fun(x1);
        cout << "NEWT " << x1 << " " << p.first-g0 << " " << p.second << "\n";
      }
    }
    x = Trigfun::root(
                      [&fx, &fy, &gx, &gy, f0]
                      (real x) -> pair<real, real>
                      { real y = fy.inv(fx(x) - f0);
                        return pair<real, real>(gx(x) + gy(y),
                                                fx.deriv(x) *
                                                (gx.gfderiv(x) +
                                                 gy.gfderiv(y))); },
                      g0, x0, xa, xb, xscale, zscale, 1,
                      countn, countb, 0, Trigfun::NEWT2);
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
    _g.ComputeInverse();
  }

  TriaxialLine::fline::disttx
  TriaxialLine::fline::Hybrid(const fics& fic,
                              Angle bet2,
                              Angle& bet2a, Angle& omg2a, Angle& alp2a)
    const {
    // XXX fix for oblate/prolate
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
      // cout << "GRR " << real(fic.bet1) << " " << real(bet2) << "\n";
      ang bet2b = bet2; bet2b.reflect(false, fic.N < 0);
      ang bet12 = bet2b - fic.bet1;
      bet12.reflect(fic.N < 0);
      // convert -180deg to 180deg
      if (signbit(bet12.s()))
        bet12 = ang(0, copysign(real(1), bet12.c()), 0, true);
      tau12 = bet12;
    } else {
      tau12 = ang::NaN();
    }
    return ArcPos0(fic, tau12.base(), bet2a, omg2a, alp2a, true);
  }

  TriaxialLine::fline::fline(const Triaxial& t, bool neg)
    : _t(t)
    , _gm(t, neg)
  {}

  TriaxialLine::fline::fline(const Triaxial& t, Triaxial::gamblk gam)
    : _t(t)
    , _gm(gam)
    // , _transpolar(signbit(_gm.gamma))
    // , _gammax(fabs(_gm.gamma))
    // , _kx2(!_transpolar ? t._k2 : t._kp2)
    // , _kxp2(_transpolar ? t._k2 : t._kp2)
    // , _kx(!_transpolar ? t._k : t._kp)
    // , _kxp(_transpolar ? t._k : t._kp)
    , _fbet(false, _t._k2 , _t._kp2,  _t._e2, -_gm.gamma, t)
    , _fomg(false, _t._kp2, _t._k2 , -_t._e2,  _gm.gamma, t)
    , _ftht( _gm.transpolar ? _fbet : _fomg)
    , _fpsi(!_gm.transpolar ? _fbet : _fomg)
    , _invp(false)
    {
      // Only needed for umbilical lines
      _deltashift = _t._k2 > 0 && _t._kp2 > 0 && _gm.gamma == 0 ?
        2 * (_fbet.Max() - _fomg.Max()) : Math::NaN();
    }

  void TriaxialLine::fline::ComputeInverse() {
    if (!_invp) {
      _fbet.ComputeInverse();
      _fomg.ComputeInverse();
      _ftht.ComputeInverse();
      _fpsi.ComputeInverse();
      _invp = true;
    }
  }

  TriaxialLine::gline::gline(const Triaxial& t, bool neg)
    : _t(t)
    , _gm(t, neg)
  {}

  TriaxialLine::gline::gline(const Triaxial& t, const Triaxial::gamblk& gm)
    : _t(t)
    , _gm(gm)
    // , _transpolar(signbit(_gm.gamma))
    // , _gammax(fabs(_gm.gamma))
    // , _kx2(!_transpolar ? t._k2 : t._kp2)
    // , _kxp2(_transpolar ? t._k2 : t._kp2)
    // , _kx(!_transpolar ? t._k : t._kp)
    // , _kxp(_transpolar ? t._k : t._kp)
    , _gbet(true, _t._k2 , _t._kp2,  _t._e2, -_gm.gamma, _t)
    , _gomg(true, _t._kp2, _t._k2 , -_t._e2,  _gm.gamma, _t)
    , _gtht( _gm.transpolar ? _gbet : _gomg)
    , _gpsi(!_gm.transpolar ? _gbet : _gomg)
    , _invp(false)
    , s0(_gm.gammax == 0 ? _gbet.Max() + _gomg.Max() : 0)
  {}

  void TriaxialLine::gline::ComputeInverse() {
    if (!_invp) {
      _gbet.ComputeInverse();
      _gomg.ComputeInverse();
      _gtht.ComputeInverse();
      _gpsi.ComputeInverse();
      _invp = true;
    }
  }

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
    // XXX fix for oblate/prolate
    disttx ret{Math::NaN(), Math::NaN(), 0};
    bool debug = false;
    real fxx = (_t.oblpro() ? sqrt(fabs(gamma())) : 1);
    if (debug)
      cout << "AP0 " << setprecision(17)
           << _t.oblpro() << " " << real(tau12) << " "
           << gamma() << " " << sqrt(fabs(gamma())) << " "
           << fic.delta / fxx  << "\n"
           << setprecision(6);
    if (gamma() > 0 || _t._kp2 == 0) {
      ang psi2;
      // ang u2a, v2a, psi2a;
      real u2, v2, u2x = 0;
      if (betp) {
        psi2 = tau12 + fic.psi1;
        v2 = fbet().fwd(psi2.radians())+0e-10;
        u2 = fomg().inv(fbet()(v2) - fic.delta);
        omg2a = ang::radians(fic.E * fomg().rev(u2));
      } else {
        omg2a = fic.omg1 + tau12.flipsign(fic.E);
        u2 = fomg().fwd(fic.E * omg2a.radians());
        // u2a = fic.E < 0 ? -omg2a : omg2a;
        u2x = fomg()(u2) + fic.delta,
        v2 = fbet().inv(u2x);
        // v2a = fbet().inv(fomg()(u2a) + fic.deltaa);
        psi2 = ang::radians(fbet().rev(v2));
        // psi2a = v2a;
        // psi2 = psi2a;
        if (0)
          cout << "AP " << real(psi2) << " " << real(tau12) << " "
               << real(omg2a) << " "
               << fic.E * omg2a.radians() << " " << u2 << " "
               << fomg()(u2) + fic.delta << " " << v2/Math::degree() << "\n";
      }
      // Already normalized
      bet2a = ang(gm().nup * psi2.s(),
                  fic.bet0.c() * hypot(psi2.c(), gm().nu * psi2.s()),
                  0, true).rebase(fic.bet0);
      // For oblate, medidional, and !betp: ...
      // psi2 = modang(u2x, 1/sqrt(gam))
      //      = atan2(u2x.s(), u2x.c()*sqrt(gam))
      // psi2.c() = u2x.c()*sqrt(gam)/abs(u2x.s())
      // nu = sqrt(gam)
      // alp2a = atan2(fic.E * hypot(_t._k * gm().nu, _t._kp * omg2a.c()),
      //            fic.bet0.c() * _t._k * gm().nup * psi2.c())
      //       = atan2(fic.E * gm().nu, fic.bet0.c() * psi2.c())
      //       = atan2(fic.E * sqrt(gam),
      //               fic.bet0.c() * u2x.c()*sqrt(gam)/abs(u2x.s())
      //       = atan2(fic.E * abs(u2x.s()),
      //               fic.bet0.c() * u2x.c())
      alp2a = (_t._kp2 == 0 && !betp && gamma() == 0 ?
        ang(fic.E * fabs(sin(u2x)), fic.bet0.c() * cos(u2x)) :
        ang(fic.E * hypot(_t._k * gm().nu, _t._kp * omg2a.c()),
            fic.bet0.c() * _t._k * gm().nup * psi2.c())).rebase(fic.alp0);
      if (0) {
        cout << "AP2 " << real(bet2a) << " "
             << real(fic.omg1) << " " << real(omg2a) << " "
             << real(alp2a) << " "
             << fic.E * hypot(_t._k * gm().nu, _t._kp * omg2a.c()) << " "
             << fic.bet0.c() * _t._k * gm().nup * psi2.c() << "\n";
        cout << "AP3 "
             << fic.u0/Math::degree()+90 << " " << fic.v0/Math::degree()+90 << " "
             << u2/Math::degree() << " " << v2/Math::degree() << " "
             << fbet().inv(fomg()(fic.u0) + fic.delta)/Math::degree() << " "
             << fic.delta/Math::degree() << "\n";
      }
      if (false) {
        for (int i = -360; i <= 360; i += 10) {
          ang a{real(i)},
            fba = fbet()(a),
            foa = fomg()(a),
            ifba = fbet()(fba),
            ifoa = fomg()(foa);
          real p = a.radians(), d = Math::degree(),
            fb = fbet()(fbet().fwd(p)),
            fo = fomg()(fomg().fwd(p)),
            ifb = fbet().inv(fb),
            ifo = fomg().inv(fo);
            cout << "QQ " << i << " "
                 << fb/d << " " << ifb/d << " "
                 << fo/d << " " << ifo/d << "\n";
            cout << "PP " << i << " "
                 << real(fba) << " " << real(ifba) << " "
                 << real(foa) << " " << real(ifoa) << "\n";
        }
      }
      ret.betw2 = v2;
      ret.omgw2 = u2;
    } else if (gamma() < 0 || _t._k2 == 0) {
      ang psi2;
      real u2, v2;
      if (betp) {
        bet2a = fic.bet1 + tau12.flipsign(fic.N);
        u2 = fbet().fwd(fic.N * bet2a.radians());
        if (gamma() == 0 && tau12 == tau12.nearest(2U)) {
          // special case for finding conjugate points on meridional geodesics,
          // gamma = 0, tau12 = multiple of pi.
          real npi = (tau12.ncardinal() + fic.psi1.nearest(2U).ncardinal())
            * Math::pi()/2;
          // Don't worry about the case where we start at a pole -- this is
          // already handled in Triaxial::Inverse.
          //
          // Solve F(tpsi2) = tpsi2 - Deltaf(n*pi + atan(tpsi2))
          //                = (tan(psi1) - Deltaf(psi1)) = c
          //
          // for tpsi2, starting guess tpsi2 = tan(psi1).
          // Test case prolate 1 3 0 1
          // p1 = [-90, -1]
          // [sx,a1x,a2x] = t.distance(p1,[90,178.9293750483]) -> a1x = 90
          // [sx,a1x,a2x] = t.distance(p1,[90,178.9293750484]) -> a1x = 90.0023
          //   omg1 = -91 -> omg2 = 178.9293750483-90 = 88.9293750483
          real c = fic.psi1.t() - fomg().df(fic.psi1.radians()),
            l = exp(Triaxial::BigValue()),
            tpsi2 = Trigfun::root([this, npi] (real tpsi) -> pair<real, real>
                                  {
                                    real psi = atan(tpsi);
                                    return pair<real, real>
                                      (tpsi - fomg().df(npi + psi),
                                       1 - fomg().dfp(psi) /
                                       (1 + Math::sq(tpsi)));
                                  },
                                  c, fic.psi1.t(), -l, l, 1, 1, 1,
                                  nullptr, nullptr, 0, Trigfun::ARCPOS0);
          psi2 = ang(tpsi2, 1) + tau12 + fic.psi1.nearest(2U);
          if (0)
            cout << "RR psi1/2 "
                 << real(fic.psi1) << " "
                 << real(psi2) << "\n";
          v2 = psi2.radians();
        } else {
          v2 = fomg().inv(fbet()(u2) - fic.delta);
          psi2 = ang::radians(fomg().rev(v2));
          if (debug) {
            cout << "ARG2 " << setprecision(17) << _t.oblpro() << " "
                 << fic.delta/fxx << " "
                 << fbet()(u2)/fxx << " "
                 << (fbet()(u2) - fic.delta)/fxx << " "
                 << v2/Math::degree() << " " << real(psi2) << " " << real(fic.psi1) << "\n"
                 << setprecision(6);
            real psi20 = 14 * Math::degree(),
              v20 = fomg().fwd(psi20),
              fomg0 = fomg()(v20),
              v21 = fomg().inv(fomg0),
              psi21 = fomg().rev(v21);
            cout << "ARG3 " << setprecision(17) << _t.oblpro() << " "
                 << psi20/Math::degree() << " " << v20 << " " << fomg0/fxx << "\n" << setprecision(6);
            cout << "ARG3 " << setprecision(17) << _t.oblpro() << " "
                 << psi21/Math::degree() << " " << v21 << " " << fomg0/fxx << "\n" << setprecision(6);
          }
        }
        if (0)
          cout << "QQX " << real(tau12) << " "
               << gamma() << " "
               << real(bet2a) << " "
               << u2 << " " << v2 << " " << real(psi2) << "\n";
      } else {
        psi2 = tau12 + fic.psi1;
        v2 = fomg().fwd(psi2.radians());
        u2 = fbet().inv(fomg()(v2) + fic.delta);
        bet2a = ang::radians(fic.N * fbet().rev(u2));
        if (0)
          cout << "FOOF1 " << real(tau12) << " " << real(fic.psi1) << " "
               << real(psi2) << " " << v2/Math::degree() << " "
               << u2/Math::degree() << " " << real(bet2a) << "\n";
      }
      // Already normalized
      omg2a = ang(gm().nup * psi2.s(),
                  fic.omg0.c() * hypot(psi2.c(), gm().nu * psi2.s()),
                  0, true).rebase(fic.omg0);
      alp2a = ang(fic.omg0.c() * _t._kp * gm().nup * psi2.c(),
                  fic.N * hypot(_t._kp * gm().nu, _t._k * bet2a.c()))
        .rebase(fic.alp0);
      if (0)
        cout << "QQY " << fic.omg0.c() * _t._kp * gm().nup * psi2.c() << " "
             << fic.N * hypot(_t._kp * gm().nu, _t._k * bet2a.c()) << "\n"
             << "FF " << fic.omg0.c() << " " << _t._kp << " "
             << gm().nup << " " << psi2.c() << " "
             << gm().nu << "\n"
             << "ALP2 " << real(alp2a) << "\n";
      ret.betw2 = u2;
      ret.omgw2 = v2;
    } else if (gamma() == 0) {
      real u2, v2;
      int ii;
      if (betp) {
        bet2a = fic.bet1 + tau12.flipsign(fic.N);
        pair<real, real> bet2n =
          remx(fic.N * (bet2a - fic.bet0).radians(), Math::pi());
        int parity = fmod(bet2n.second, real(2)) != 0 ? -1 : 1;
        real deltax = clamp(fic.delta + bet2n.second * _deltashift, 2);
        u2 = fbet().fwd(bet2n.first);
        v2 = fomg().inv(fbet()(u2) - deltax);
        omg2a = ang::radians(fic.E * parity * fomg().rev(v2))
          .rebase(fic.omg0);
        // umbalt definition of alp0
        ang alp0x(fic.alp1.nearest(2U));
        // Conflict XXX
        // testset -50 180 20 0 want fic.N multiplying s()
        // testspha -20 90 20 -90 wants fic.N multiplying c()
        alp2a = ang(_t._kp * fic.E *
                    parity / mcosh(v2, _t._k),
                    fic.N * _t._k / mcosh(u2, _t._kp)).rebase(alp0x);
        ii = int(bet2n.second);
        // Move forward from umbilical point
        bet2a += ang::eps().flipsign(fic.N);
      } else {
        omg2a = fic.omg1 + tau12.flipsign(fic.E);
        pair<real, real> omg2n =
          remx(fic.E * (omg2a - fic.omg0).radians(), Math::pi());
        int parity = fmod(omg2n.second, real(2)) != 0 ? -1 : 1;
        real deltax = clamp(fic.delta + omg2n.second * _deltashift, 2);
        v2 = fomg().fwd(omg2n.first);
        u2 = fbet().inv(fomg()(v2) + deltax);
        real bet2 = fic.N * parity * fbet().rev(u2);
        bet2a = ang::radians(bet2);
        // !umbalt definition of alp0
        ang alp0x(fic.alp1.nearest(1U));
        alp2a = ang(fic.E * _t._kp / mcosh(v2, _t._k),
                    _t._k * fic.N *
                    parity / mcosh(u2, _t._kp)).rebase(alp0x);
        ii = int(omg2n.second);
        // Move forward from umbilical point
        omg2a += ang::eps().flipsign(fic.E);
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
    E = signbit(alp1.s()) ? -1 : 1;
    N = signbit(alp1.c()) ? -1 : 1;
    if (gm.gamma > 0 || t._kp2 == 0) {
      bet0 = bet1.nearest(2U);
      alp0 = alp1.nearest(1U);
      psi1 = ang(t._k * bet1.s(),
                 bet0.c() * alp1.c() *
                 hypot(t._k * bet1.c(), t._kp * omg1.c()));
      // k = 1, kp = 0, sqrt(-mu) = bet1.c() * fabs(alp1.s())
      // modang(psi1, sqrt(-mu)) = atan2(bet1.s() * fabs(alp1.s()),
      //                                 bet0.c() * alp1.c());
      // assume fbet().fwd(x) = x in this case
      v0 = f.fbet().fwd(psi1.radians());
      u0 = f.fomg().fwd(E * omg1.radians());
      // Only used for biaxial cases when fwd rev is the identity
      // v0a = psi1;
      // u0a = E < 0 ? -omg1 : omg1;
      delta = (t._kp2 == 0 && (t._oblpro || gm.gamma == 0) ?
               atan2(bet1.s() * fabs(alp1.s()), bet0.c() * alp1.c())
               - sqrt(gm.gamma) * f.fbet().df(v0)
               : f.fbet()(v0)) - f.fomg()(u0);
      // deltaa = f.fbet()(v0a) - f.fomg()(u0a);
      if (0) {
        real d = Math::degree();
        cout << "V0 " << v0/d << " "
             << atan2(0 * psi1.s(), psi1.c())/d << "\n";
        cout << "INFIC " << real(alp1) << " " << real(psi1) << " "
             << cos(psi1.radians()) << " " << psi1.c() << " "
             << real(bet1) << " " << real(bet0) << " "
             << v0/d << " " << u0/d << " "
             << (t._kp2 == 0 && (t._oblpro || gm.gamma == 0) ?
                 atan2(bet1.s() * fabs(alp1.s()), bet0.c() * alp1.c())
                 - sqrt(gm.gamma) * f.fbet().df(v0)
                 : f.fbet()(v0))/d << " "
             << f.fomg()(u0)/d << " " << delta/d << "\n";
      }
    } else if (gm.gamma < 0 || t._k2 == 0) {
      omg0 = omg1.nearest(2U);
      alp0 = alp1.nearest(2U);
      // Need Angle(0, 0) to be treated like Angle(0, 1) here.
      psi1 = ang(t._kp * omg1.s(),
                 omg0.c() * alp1.s() *
                 hypot(t._k * bet1.c(), t._kp * omg1.c()));
      // For k = 0, kp = 1, sqrt(-mu) = omg1.c() * fabs(alp1.c())
      // modang(psi1, sqrt(-mu)) = atan2(omg1.s() * fabs(alp1.c()),
      //                                 omg0.c() * alp1.s());
      // assume fomg().fwd(x) = x in this case
      if (0)
        cout << "FICS " << real(omg1) << " " << real(omg0) << " "
             << real(psi1) << " "
             << t._kp * omg1.s() << " "
             << omg0.c() * alp1.s() * hypot(t._k * bet1.c(),
                                            t._kp * omg1.c()) << "\n";
      v0 = f.fomg().fwd(psi1.radians());
      u0 = f.fbet().fwd(N * bet1.radians());
      delta = f.fbet()(u0) -
        (t._k2 == 0 && (t._oblpro || gm.gamma == 0) ?
         atan2(omg1.s() * fabs(alp1.c()), omg0.c() * alp1.s())
         - sqrt(-gm.gamma) * f.fomg().df(v0) :
         f.fomg()(v0));
    } else if (gm.gamma == 0) {
      alp0 = alp1.nearest(signbit(f.gamma()) ? 2U : 1U);
      // N.B. factor of k*kp omitted
      // bet0, omg0 are the middle of the initial umbilical segment
      if (fabs(bet1.c()) < 8*eps && fabs(omg1.c()) < 8*eps) {
        bet0 = bet1.nearest(1U) + ang::cardinal(N);
        omg0 = omg1.nearest(1U) + ang::cardinal(E);
        delta = f.deltashift()/2 - log(fabs(alp1.t()));
      } else {
        bet0 = bet1.nearest(2U);
        omg0 = omg1.nearest(2U);
        delta = N * f.fbet()(lamang(bet1 - bet0, t._kp)) -
          E * f.fomg()(lamang(omg1 - omg0, t._k));
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
    int oE = E, oN = N;
    E = signbit(alp1.s()) ? -1 : 1;
    N = signbit(alp1.c()) ? -1 : 1;
    if (gam > 0 || t._kp2 == 0) {
      alp0 = alp1.nearest(1U);
      psi1.reflect(false, N != oN);
      v0 = f.fbet().fwd(psi1.radians());
      u0 *= E/oE;
      // This isn't right if _oblpro
      delta = f.fbet()(v0) - f.fomg()(u0);
    } else if (gam < 0 || t._k2 == 0) {
      alp0 = alp1.nearest(2U);
      psi1.reflect(false, E != oE);
      v0 = f.fomg().fwd(psi1.radians());
      u0 *= N/oN;
      // This isn't right if _oblpro
      delta = f.fbet()(u0) - f.fomg()(v0);
    } else if (gam == 0) {
      // Only expect to invoke setquadrant in this case
      alp0 = alp1.nearest(t._umbalt ? 2U : 1U);
      if (fabs(bet1.c()) < 8*eps && fabs(omg1.c()) < 8*eps)
        delta = f.deltashift()/2 - log(fabs(alp1.t()));
      else
        delta = N * f.fbet()(lamang(bet1 - bet0, t._kp)) -
          E * f.fomg()(lamang(omg1 - omg0, t._k));
    } else {
      // gamma = NaN
    }
  }

  TriaxialLine::fline::ficsx::ficsx(const fline& f,
                                    Angle bet10, Angle omg10,
                                    Angle alp10)
    : tht1(omg10 - ang::cardinal(1))
    , phi1(bet10)
    , alp1(alp10)
  {
    if (f.transpolar()) {
      swap(tht1, phi1);
      alp1.reflect(false, false, true);
    }
    const real eps = numeric_limits<real>::epsilon();
    const Triaxial& t = f.t();
    if (!f.transpolar() && phi1.s() == 0 && fabs(alp1.c()) <= Math::sq(eps))
      alp1 = ang(alp1.s(), - Math::sq(eps), alp1.n(), true);
    E = signbit(alp1.s()) ? -1 : 1;
    N = signbit(alp1.c()) ? -1 : 1;
    if (f.gammax() > 0 || f.kxp2() == 0) {
      phi0 = phi1.nearest(2U);
      alp0 = alp1.nearest(1U);
      psi1 = ang(f.kx() * phi1.s(),
                 phi0.c() * alp1.c() *
                 hypot(f.kx() * phi1.c(), f.kxp() * tht1.c()));
      // k = 1, kp = 0, sqrt(-mu) = bet1.c() * fabs(alp1.s())
      // modang(psi1, sqrt(-mu)) = atan2(bet1.s() * fabs(alp1.s()),
      //                                 bet0.c() * alp1.c());
      // assume fbet().fwd(x) = x in this case
      v0 = f.fpsi().fwd(psi1.radians());
      u0 = f.ftht().fwd(E * tht1.radians());
      // Only used for biaxial cases when fwd rev is the identity
      // v0a = psi1;
      // u0a = E < 0 ? -omg1 : omg1;
      delta = (f.kxp2() == 0 && (t._oblpro || f.gammax() == 0) ?
               atan2(phi1.s() * fabs(alp1.s()), phi0.c() * alp1.c())
               - sqrt(f.gammax()) * f.fpsi().df(v0)
               : f.fpsi()(v0)) - f.ftht()(u0);
      // deltaa = f.fbet()(v0a) - f.fomg()(u0a);
    } else if (f.gammax() == 0) {
      alp0 = alp1.nearest(signbit(f.gamma()) ? 2U : 1U);
      // N.B. factor of k*kp omitted
      // bet0, omg0 are the middle of the initial umbilical segment
      if (fabs(phi1.c()) < 8*eps && fabs(tht1.c()) < 8*eps) {
        phi0 = phi1.nearest(1U) + ang::cardinal(N);
        tht0 = tht1.nearest(1U) + ang::cardinal(E);
        delta = f.deltashift()/2 - log(fabs(alp1.t()));
      } else {
        phi0 = phi1.nearest(2U);
        tht0 = tht1.nearest(2U);
        delta = N * f.fpsi()(lamang(phi1 - phi0, t._kp)) -
          E * f.ftht()(lamang(tht1 - tht0, t._k));
      }
    } else {
      // gamma = NaN
    }
  }

  void TriaxialLine::fline::ficsx::setquadrant(const fline& f, unsigned q) {
    const real eps = numeric_limits<real>::epsilon();
    const Triaxial& t = f.t();
    alp1.setquadrant(q);
    if (phi1.s() == 0 && fabs(alp1.c()) <= Math::sq(eps))
      alp1 = ang(alp1.s(), -Math::sq(eps), alp1.n(), true);
    int oE = E, oN = N;
    E = signbit(alp1.s()) ? -1 : 1;
    N = signbit(alp1.c()) ? -1 : 1;
    if (f.gammax() > 0 || f.kxp2() == 0) {
      alp0 = alp1.nearest(1U);
      psi1.reflect(false, N != oN);
      v0 = f.fpsi().fwd(psi1.radians());
      u0 *= E/oE;
      // This isn't right if _oblpro
      delta = f.fpsi()(v0) - f.ftht()(u0);
    } else if (f.gammax() == 0) {
      // Only expect to invoke setquadrant in this case
      alp0 = alp1.nearest(t._umbalt ? 2U : 1U);
      if (fabs(phi1.c()) < 8*eps && fabs(tht1.c()) < 8*eps)
        delta = f.deltashift()/2 - log(fabs(alp1.t()));
      else
        delta = N * f.fpsi()(lamang(phi1 - phi0, f.kxp())) -
          E * f.ftht()(lamang(tht1 - tht0, f.kx()));
    } else {
      // gamma = NaN
    }
  }

  TriaxialLine::gline::gicsx::gicsx(const gline& g, const fline::fics& fic) {
    const Triaxial& t = g.t();
    if (g.gamma() > 0 || t._kp2 == 0) {
      sig1 = g.gbet()(fic.v0) + g.gomg()(fic.u0);
    } else if (g.gamma() < 0 || t._k2 == 0) {
      sig1 = g.gbet()(fic.u0) + g.gomg()(fic.v0);
    } else if (g.gamma() == 0) {
      sig1 = fic.N *
        g.gbet()(lamang(fic.bet1 - fic.bet0, t._kp)) +
        fic.E * g.gomg()(lamang(fic.omg1 - fic.omg0, t._k));
    } else {
      // gamma = NaN
    }
  }

  TriaxialLine::gline::gics::gics(const gline& g, const fline::fics& fic) {
    const Triaxial& t = g.t();
    if (g.gamma() > 0 || t._kp2 == 0) {
      sig1 = g.gbet()(fic.v0) + g.gomg()(fic.u0);
    } else if (g.gamma() < 0 || t._k2 == 0) {
      sig1 = g.gbet()(fic.u0) + g.gomg()(fic.v0);
    } else if (g.gamma() == 0) {
      sig1 = fic.N *
        g.gbet()(lamang(fic.bet1 - fic.bet0, t._kp)) +
        fic.E * g.gomg()(lamang(fic.omg1 - fic.omg0, t._k));
    } else {
      // gamma = NaN
    }
  }

  Math::real TriaxialLine::gline::dist(gics ic, fline::disttx d) const {
    real sig2 = gbet()(d.betw2) + gomg()(d.omgw2) + d.ind2 * 2*s0;
    return (sig2 - ic.sig1) * t()._b;
  }

  void TriaxialLine::inversedump(ostream& os) const {
    os << "[b, e2, k2, kp2, gam] = deal("
       << _t.b() << ", " << _t.e2() << ", "
       << _t.k2() << ", " << _t.kp2() << ", "
       << gamma() << ");\n";
    os << "tx = ["
       << fbet().txp() << ", " << gbet().txp() << ", "
       << fomg().txp() << ", " << gomg().txp() << "];\n" ;
    _f.inversedump(os);
    _g.inversedump(os);
  }

  void TriaxialLine::fline::inversedump(ostream& os) const {
    _fbet.inversedump(os, "fbet");
    _fomg.inversedump(os, "fomg");
  }

  void TriaxialLine::gline::inversedump(ostream& os) const {
    _gbet.inversedump(os, "gbet");
    _gomg.inversedump(os, "gomg");
  }

  TriaxialLine::hfun::hfun(bool distp, real kap, real kapp, real eps, real mu,
                           const Triaxial& t)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _sqrtkap(sqrt(_kap))
    , _sqrtkapp(sqrt(_kapp))
    , _distp(distp)
      // If oblpro extend special treatment of oblate/prolate cases to mu != 0.
    , _oblpro(t._oblpro)
    , _merid(t._merid)
    , _invp(false)
    , _biaxr(_kap == 0 && (_mu == 0 || (_oblpro && _mu > 0)))
    , _biaxl(_kapp == 0 && (_mu == 0 || (_oblpro && _mu < 0)))
    , _umb(_kap != 0 && _kapp != 0 && _mu == 0)
  {
    (void) _merid;
    // mu in [-kap, kapp], eps in (-inf, 1/kap)
    if (!_distp) {
      if (_biaxr) {
        // oblate/prolate rotating coordinate
        // _kapp == 1, mu < 0 not allowed
        _tx = false;
        // f multiplied by sqrt(mu)
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real tht) -> real
                          // This is a trivial case f' = 1
                          { return fthtoblp(tht, eps, mu); },
                          Math::pi()/2, false);
      } else if (_biaxl) {
        // oblate/prolate librating coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < t._ellipthresh;
        _tx = false;
        // f multiplied by sqrt(-mu)
        /*
          if (_tx) {
          _ell = EllipticFunction(1 + _mu, 0, -_mu, 1);
          _fun = TrigfunExt(
          [kap = _kap, kapp = _kapp,
          eps = _eps, mu = _mu, ell = _ell]
          (real v) -> real
          { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
          return dfvoblp(dn, eps, mu); },
          _ell.K(), false, 1);
          } else
        */
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return dfpsioblp(sin(psi), cos(psi), eps, mu); },
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
                            (real tht) -> real
                            { return fthtp(cos(tht), kap, kapp, eps, mu); },
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
      } else if (_umb) {
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
                            (real tht) -> real
                            { return dfp(cos(tht), kap, kapp, eps); },
                            Math::pi(), true, 1);
      } else {
        // _mu == NaN
        _tx = false;
      }
    } else {
      if (_biaxr) {
        // oblate/prolate symmetry coordinate
        // _kapp == 1, mu < 0 not allowed
        _tx = false;
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real tht) -> real
                          // degenerate f' = 0
                          { return gthtoblp(tht, eps, mu); },
                          Math::pi()/2, false);
      } else if (_biaxl) {
        // oblate/prolate non-symmetry coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < t._ellipthresh;
        _tx = false;
        /*
          if (_tx) {
          _ell = EllipticFunction(1 + _mu, 0, -_mu, 1);
          // Never completed this
          } else
        */
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return gpsioblp(sin(psi), cos(psi), eps, mu); },
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
                            (real tht) -> real
                            { return gthtp(cos(tht), kap, kapp, eps, mu); },
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
      } else if (_umb) {
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
                            [kap = _kap, kapp = _kapp, eps = _eps]
                            (real tht) -> real
                            { return g0p(cos(tht), kap, kapp, eps); },
                            Math::pi(), true);
      } else {
        // _mu == NaN
        _tx = false;
      }
    }
    // N.B. _max < 0 for _umb && eps < 0
    _max = !_distp ?
      ( _umb ? _fun(_tx ? _ell.K() : Math::pi()/2) :
        (_biaxl ? Math::pi()/2 + sqrt(-_mu) * _fun.Max() : _fun.Max()) ) :
      ( _umb ? _fun(_tx ? _ell.K() : Math::pi()/2) : _fun.Max() );
  }

  Math::real TriaxialLine::hfun::operator()(real u) const {
    if (!_distp) {
      if (_biaxl) {
        // This is sqrt(-mu) * f(u)
        return modang(u, sqrt(-_mu)) - sqrt(-_mu) * _fun(u);
      } else if (_umb) {
        // This is sqrt(kap * kapp) * f(u)
        real phi = gd(u, _sqrtkapp);
        return u - _fun(_tx ? _ell.F(phi) : phi);
      } else
        return _fun(u);
    } else {
      if (_umb) {
        real phi = gd(u, _sqrtkapp);
        return _fun(_tx ? _ell.F(phi) : phi);
      } else
        return _fun(u);
    }
  }

  // THIS ISN"T USED ?
  // Should implement an Angle equivalant of _ell.F(phi)
  Angle TriaxialLine::hfun::operator()(const Angle& ang) const {
    if (_distp) return Angle::NaN();
    if (_biaxr)
      return ang;
    else if (_biaxl && _mu == 0)
      return ang.modang(sqrt(-_mu));
    real u = ang.radians();
    if (_biaxl)
      // This is sqrt(-mu) * f(u)
      return ang.modang(sqrt(-_mu)) - ang::radians(sqrt(-_mu) * _fun(u));
    else if (_umb) {
      // This is sqrt(kap * kapp) * f(u)
      real phi = gd(u, _sqrtkapp);
      return ang::radians(u - _fun(_tx ? _ell.F(phi) : phi));
    } else
      return ang::radians(_fun(u));
  }

  Math::real TriaxialLine::hfun::deriv(real u) const {
    if (!_distp) {
      if (_biaxl) {
        // This is sqrt(-mu) * f'(u)
        /*
          if (_tx) {
          // DLMF (22.6.1): sn^2 + cn^2 =  k^2*sn^2 + dn^2 = 1
          // dn^2 = cn^2 + k'^2*sn^2 = cn^2 - mu*sn^2
          // k'^2 = -mu
          real sn, cn, dn;
          (void) _ell.am(u, sn, cn, dn);
          return sqrt(-_mu) / Math::sq(dn) - _fun.deriv(u);
          } else
        */
        // f0(x) = atan(sqrt(-mu) * tanx(u))
        // f0'(x) = sqrt(-mu) / (cos(u)^2 - mu * sin(u)^2)
        // **HERE**
        return sqrt(-_mu) / (Math::sq(cos(u)) - _mu * Math::sq(sin(u)))
          - sqrt(-_mu) * _fun.deriv(u);
      } else if (_umb) {
        // This is sqrt(kap * kapp) * f'(u)
        real phi = gd(u, _sqrtkapp),
          t = _kapp + Math::sq(sinh(u));
        // dphi/du = _sqrtkapp * cosh(u) / t;
        // for tx w = F(phi, sqrtkap)
        //   dw/dphi = 1/sqrt(kapp + kap*cos(phi)^2)
        //           = sqrt(t)/(sqrtkapp *cosh(u))
        //   N.B. cos(phi)^2 = kapp/t
        //   dw/du = dw/dphi*dphi/du = 1/sqrt(t)
        return 1 - _fun.deriv(_tx ? _ell.F(phi) : phi) /
          ( _tx ? sqrt(t) : t / (_sqrtkapp * cosh(u)) );
      } else
        return _fun.deriv(u);
    } else {
      if (_umb) {
        real phi = gd(u, _sqrtkapp),
          t = _kapp + Math::sq(sinh(u));
        // See comments in ffun::deriv
        return _fun.deriv(_tx ? _ell.F(phi) : phi) /
          ( _tx ? sqrt(t) : t / (_sqrtkapp * cosh(u)) );
      } else
        return _fun.deriv(u);
    }
  }

  Math::real TriaxialLine::hfun::gfderiv(real u) const {
    // return g'(u)/f'(u)
    real sn = 0, cn = 0, dn = 0;
    if (_biaxr)
      return gfthtoblp(u, _mu);
    else if (_biaxl) // mu > 0 not allowed
      return gfpsioblp(sin(u), cos(u), _mu);
    else if (_umb)
      // This includes factor of sqrt(kap * kapp) because of adjustment of
      // definition of f for umbilical geodesics.
      return gf0up(u, _kap, _kapp);
    else {
      if (_tx)
        (void) _ell.am(u, sn, cn, dn);
      if (_mu > 0)
        return _tx ? gfup(cn, _kap, _mu) : gfthtp(cos(u), _kap, _mu);
      else if (_mu < 0)
        return _tx ? gfvp(dn, _kap, _mu) :
          gfpsip(sin(u), cos(u), _kap, _mu);
      else
        return Math::NaN();
    }
  }

  void TriaxialLine::hfun::ComputeInverse() {
    if (!_distp) {
      if (!_invp) {
        if (_biaxl) {
          if (_mu == 0) return; // _fun == 0 and there's an analytic inverse
          // now _mu < 0
          _countn = _countb = 0;
          // Include scale = 1 in TrigfunExt constructor because _dfinv gets
          // added to u.
          // Ars are fun, odd, sym, halfp, nmax, tol, scale
          // **HERE**
          _dfinv = Trigfun(
                           [this]
                           (real z, real u1) -> real
                           {
                             real u0 = modang(z/Slope(), 1/sqrt(-_mu));
                             return root(z, u0 + u1, &_countn, &_countb,
                                         sqrt(numeric_limits<real>::epsilon()))
                               - u0;
                           },
                           true, false, Slope() * HalfPeriod(),
                           int(ceil(real(1.5) * NCoeffs())),
                           sqrt(numeric_limits<real>::epsilon()), HalfPeriod());
        } else if (_umb) {
          _countn = _countb = 0;
          // Include scale = 1 in TrigfunExt constructor because _dfinv gets
          // added to z.
          _dfinv = Trigfun(
                           [this]
                           (real phi, real u1) -> real
                           {
                             real z = lam(phi, _sqrtkapp);
                             return root(z, z + u1, &_countn, &_countb,
                                         sqrt(numeric_limits<real>::epsilon()))
                               - z;
                           },
                           true, true, Math::pi(),
                           int(ceil(real(1.5) * NCoeffs())),
                           sqrt(numeric_limits<real>::epsilon()), 1);
        } else
          _fun.ComputeInverse();
      }
      _invp = true;
    } else {
      if (!(_invp || _biaxr || _umb)) {
        // If _umb, the inverse isn't periodic
        _fun.ComputeInverse();
        /*
          real u = 1, z = _fun(u);
          cout << "HERE " << u << " " << z << " " << _fun.inv0(z) << "\n";
          for (int i = -100; i <= 100; ++i) {
          real u = real(i)/10;
          cout << "DD " << u << " " << _fun(u) << "\n";
          }
        */
        _invp = true;
      }
    }
  }

  Math::real TriaxialLine::hfun::root(real z, real u0,
                                      int* countn, int* countb,
                                      real tol) const {
    if (!_distp) {
      if (!isfinite(z)) return z; // Deals with +/-inf and nan
      if (_biaxl) {
        // z = 1/sqrt(-mu)
        // here f(u) = atan(sqrt(-mu)*tan(u))/sqrt(-mu)-Deltaf(u)
        // let fx(u) = sqrt(-mu) * f(u); fun(u) = sqrt(-mu)*Deltaf(u)
        // z = fx(u) = modang(u, sqrt(-mu)) - fun(u)
        // fun(u) is monotonically increasing/decreasing quasilinear function
        // period pi.  Combined function is quasilinear
        // z = s * u +/- m; period = pi
        // Inverting:
        // u = z/s +/- m/s; period s*pi
        // u = fxinv(z) = modang(z/x, 1/sqrt(-mu)) + funinv(z)
        // funinv(z) periodic function of z period (1-s)*pi
        real d = Max(),
          ua = (z - d) / Slope(),
          ub = (z + d) / Slope();
        u0 = fmin(ub, fmax(ua, u0));
        return Trigfun::root(
                             [this]
                             (real u) -> pair<real, real>
                             { return pair<real, real>((*this)(u), deriv(u)); },
                             z,
                             u0, ua, ub,
                             HalfPeriod(), HalfPeriod()/Slope(), 1,
                             countn, countb, tol, Trigfun::FFUNROOT);
      } else if (_umb) {
        real d = fabs(Max())
          + 2 * numeric_limits<real>::epsilon() * fmax(real(1), fabs(z)),
          ua = z - d,
          ub = z + d;
        u0 = fmin(ub, fmax(ua, u0));
        return Trigfun::root(
                             [this]
                             (real u) -> pair<real, real>
                             { return pair<real, real>((*this)(u), deriv(u)); },
                             z,
                             u0, ua, ub,
                             Math::pi()/2, Math::pi()/2, 1, countn, countb, tol,
                             Trigfun::FFUNROOT);
      } else
        return Math::NaN();
    } else {
      // This function isn't neeed.  General inversion mechanisms in Trigfun
      // suffice.  NO, the trigfun for _umb is not invertible.
      if (!(isfinite(z) && _umb))
        return Math::NaN();       // Deals with +/-inf and nan
      // Now we're dealing with _umb.
      if (fabs(z) >= Max())
        return copysign(Triaxial::BigValue(), z);
      real ua = -Triaxial::BigValue(), ub = -ua;
      u0 = fmin(ub, fmax(ua, u0));
      // Solve z = _fun(_tx ? _ell.F(gd(u)) : gd(u)) for u
      return Trigfun::root(
                           [this]
                           (real u) -> pair<real, real>
                           { return pair<real, real>((*this)(u), deriv(u)); },
                           z,
                           u0, ua, ub,
                           Math::pi()/2, Math::pi()/2, 1, countn, countb, tol,
                           Trigfun::GFUNROOT);
    }
  }

  // Approximate inverse using _dfinv _fun.inv0
  Math::real TriaxialLine::hfun::inv0(real z) const {
    if (_distp) {
      if (!_invp) return Math::NaN();
      // For the inverse in the umbilical case, just use gd(z, _sqrtkapp) and
      // not F(gd(z, _sqrt(kapp)))
      return _umb ? z + _dfinv(gd(z, _sqrtkapp)) :
        (_biaxl ? modang(z/Slope(), 1/sqrt(-_mu)) + (_mu == 0 ? 0 : _dfinv(z)) :
         _fun.inv0(z));
    } else {
      return _invp ? _fun.inv0(z) :
        (_umb ?
         // In limit _eps -> 0
         //   g(u) = atan(_sqrtkap/_sqrtkapp * tanh(u))
         // at u = 0, dg/du = _sqrtkap/_sqrtkapp
         //    u = inf, g = atan(_sqrtkap/_sqrtkapp)
         //
         // For _eps finite
         // at u = 0, dg/du = _sqrtkap/_sqrtkapp * sqrt(1 - _eps*_kap)
         //    u = inf, g = _max
         // Note: at u = 0
         //   du/dphi = _sqrtkapp, so
         //   dg/dphi = _sqrtkap * sqrt(1 - _eps*_kap) -- OK
         //
         // Approximate g(u) for _eps finite by
         // g(u) = _max / atan(_sqrtkap/_sqrtkapp) *
         //   atan(_sqrtkap/_sqrtkapp *
         //        tanh(sqrt(1 - _eps*_kap) * atan(_sqrtkap/_sqrtkapp) / _max
         //             * u))
         // Values at +/- inf and slope at origin match.
         //
         // Solve z = g(u) gives u = u0:
         atanh(_sqrtkapp/_sqrtkap * tan(atan(_sqrtkap/_sqrtkapp) / _max * z)) /
         (sqrt(1 - _eps*_kap) * atan(_sqrtkap/_sqrtkapp) / _max)
         : Math::NaN());
    }
  }

  // Accurate inverse by direct Newton (not using _finv)
  Math::real TriaxialLine::hfun::inv1(real z, int* countn, int* countb) const {
    if (!_distp)
      return _umb ? root(z, z, countn, countb) :
        (_biaxl ? (_mu == 0 ?
                   // In this case _fun.Slope() = 0 and Slope() = 1
                   modang(z/Slope(), 1/sqrt(-_mu)) :
                   root(z, modang(z/Slope(), 1/sqrt(-_mu)), countn, countb)) :
         _fun.inv1(z, countn, countb));
    else {
      if (_biaxr) return Math::NaN();
      return _umb ? root(z, inv0(z), countn, countb) :
        _fun.inv1(z, countn, countb);
    }
  }

  // Accurate inverse correcting result from _finv
  Math::real TriaxialLine::hfun::inv2(real z, int* countn, int* countb) const {
    if (!_invp) return Math::NaN();
    if (!_distp)
      return _biaxl && _mu == 0 ? inv1(z) :
        (_umb || _biaxl ? root(z, inv0(z), countn, countb) :
         _fun.inv2(z, countn, countb));
    else
      return _umb ? root(z, inv0(z), countn, countb) :
        _fun.inv1(z, countn, countb);
  }

  Angle TriaxialLine::hfun::inv(const Angle& z, int* countn, int* countb)
    const {
    if (!_distp) return ang::NaN();
    if (_biaxr)
      return z;
    else if (_biaxl && _mu == 0)
      return z.modang(1/sqrt(-_mu));
    else
      return ang::radians(inv(z.radians(), countn, countb));
  }

  // _mu > 0 && !_tx
  Math::real TriaxialLine::hfun::fthtp(real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }
  // This is non-negative
  Math::real TriaxialLine::hfun::gthtp(real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return c2 * sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }
  Math::real TriaxialLine::hfun::gfthtp(real c, real kap, real /* mu */) {
    real c2 = kap * Math::sq(c);
    return c2;
  }

  // _mu > 0 && _tx
  Math::real TriaxialLine::hfun::fup(real cn, real kap, real kapp,
                                     real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }
  // This is non-negative
  Math::real TriaxialLine::hfun::gup(real cn, real /* dn */,
                                     real kap, real kapp,
                                     real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return c2 * sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }
  Math::real TriaxialLine::hfun::gfup(real cn, real kap, real /* mu */) {
    real c2 = kap * Math::sq(cn);
    return c2;
  }

  // _mu == 0 && !_tx
  Math::real TriaxialLine::hfun::dfp(real c,
                                     real kap, real kapp, real eps) {
    // function dfp = dfpf(phi, kappa, epsilon)
    // return derivative of sqrt(kap * kapp) * Delta f
    // s = sqrt(1 - kap * sin(phi)^2)
    real c2 = kap * Math::sq(c), s = sqrt(kapp + c2);
    return eps*kap * sqrt(kapp) * c / (s * (1 + sqrt(1 - eps*c2)));
  }
  Math::real TriaxialLine::hfun::g0p(real c, real kap, real kapp, real eps) {
    real c2 = kap * Math::sq(c);
    return sqrt( kap * (1 - eps * c2) / (kapp + c2) ) * c;
  }

  // _mu == 0 && _tx
  Math::real TriaxialLine::hfun::dfvp(real cn, real /* dn */,
                                      real kap, real kapp, real eps) {
    // function dfvp = dfvpf(v, kap, eps)
    // return derivative of sqrt(kap * kapp) * Delta f_v
    return eps*kap * sqrt(kapp) * cn /
      (1  + sqrt(1 - eps*kap * Math::sq(cn)));
  }
  Math::real TriaxialLine::hfun::g0vp(real cn, real kap, real /* kapp */,
                                      real eps) {
    real c2 = kap * Math::sq(cn);
    return sqrt( kap * (1 - eps * c2) ) * cn;
  }

  // _mu == 0 (_tx ignored)
  Math::real TriaxialLine::hfun::gf0up(real u, real kap, real kapp) {
    // Subst tan(phi) = sinh(u) /sqrt(kapp) in
    // kap * cos(phi)^2 gives kap*kapp/(kapp + sinh(u)^2)
    // Divide by sqrt(kap * kappp) to account of factor removed from f
    // functions.
    return sqrt(kap * kapp) / ( kapp + Math::sq(sinh(u)) );
  }

  // _mu < 0 && !_tx
  Math::real TriaxialLine::hfun::fpsip(real s, real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * c2) ) ;
  }
  // This is positive
  Math::real TriaxialLine::hfun::gpsip(real s, real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt(c2 * (1 - eps * c2) / (kapp + c2)) ;
  }
  Math::real TriaxialLine::hfun::gfpsip(real s, real c, real kap, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return c2;
  }

  // _mu < 0 && _tx
  Math::real TriaxialLine::hfun::fvp(real dn, real kap, real kapp,
                                     real eps, real /* mu */) {
    real c2 = kap * Math::sq(dn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * kap) );
  }
  // This is positive
  Math::real TriaxialLine::hfun::gvp(real /* cn */, real dn,
                                     real kap, real kapp,
                                     real eps, real /* mu */) {
    real dn2 = Math::sq(dn), c2 = kap * dn2;
    return dn2 * sqrt( kap * (1 - eps * c2) / (kapp + c2) );
  }
  Math::real TriaxialLine::hfun::gfvp(real dn, real kap, real /* mu */) {
    real dn2 = Math::sq(dn), c2 = kap * dn2;
    return c2;
  }

  // oblate/prolate variants for kap = 0, kapp = 1, mu > 0
  Math::real TriaxialLine::hfun::fthtoblp(real /* tht */, real /* eps */,
                                          real /* mu */) {
    // Multiply by f functions by sqrt(abs(mu))
    // return 1 / sqrt(mu);
    return 1;
  }
  Math::real TriaxialLine::hfun::gthtoblp(real /* tht */, real /* eps */,
                                          real /* mu */) {
    return 0;
  }
  Math::real TriaxialLine::hfun::gfthtoblp(real /* tht */, real /* mu */) {
    return 0;
  }

  // oblate/prolate variants for kap = 1, kapp = 0, mu <= 0, !_tx
  Math::real TriaxialLine::hfun::dfpsioblp(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // f functions are multiplied by sqrt(abs(mu)) but don't include this
    // factor here; instead include it in operator()(). etc.  This was we can
    // still use this function in the limit mu -> 0 to determine the conjugate
    // point on a meridian.
    return eps / (1 + sqrt(1 - eps * c2));
  }
  Math::real TriaxialLine::hfun::gpsioblp(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // return gpsip(s, c, 1, 0, eps, mu) but with the factor c2 canceled
    return sqrt(1 - eps * c2);
  }
  Math::real TriaxialLine::hfun::gfpsioblp(real s, real c, real mu) {
    // cos(phi)^2 = cos(psi)^2 - mu *sin(psi)^2
    // f' = sqrt(1-eps*cos(phi)^2)/cos(phi)^2
    // g' = sqrt(1-eps*cos(phi)^2)
    // g'/f' = cos(phi)^2
    // Adjust by sqrt(-mu) to accommodate this factor in dfpsioblp
    return gfpsip(s, c, 1, mu) / sqrt(-mu);
  }

#if 0
  // oblate/prolate variants for kap = 1, kapp = 0, mu <= 0, _tx
  Math::real TriaxialLine::hfun::dfvoblp(real dn, real eps, real mu) {
    real c2 = Math::sq(dn);
    // Multiply by f functions by sqrt(abs(mu))
    return sqrt(-mu) * eps * dn / (1 + sqrt(1 - eps * c2));
  }
  // This is positive
  Math::real TriaxialLine::hfun::gvoblp(real /* cn */, real dn,
                                        real eps, real /* mu */) {
    // return gvp(dn, cn, 1, 0, epd, mu) but with cancelation
    real c2 = Math::sq(dn);
    return dn * sqrt(1 - eps * c2);
  }
  Math::real TriaxialLine::hfun::gfvoblp(real dn, real mu) {
    return gfvp(dn, 1, mu) / sqrt(-mu);
  }
#endif

  void TriaxialLine::hfun::inversedump(ostream& os, const string& name) const {
    os << "% " << name << "\n";
    int ndiv = 20;
    real ds = 1/real(100);
    int num = int(round(fmin(3*HalfPeriod(), real(30)) * ndiv));
    if (_distp)
      os << "% u g(u) inv0(g) inv1(g) inv2(g) "
         << "u0-u u1-u u2-u  g'(u) delg delg-g' g'/f'\n";
    else
      os << "% u f(u) inv0(f) inv1(f) inv2(f) "
         << "u0-u u1-u u2-u  f'(u) delf delf-f' phi fwd(phi) fwd(phi)-u\n";
    os << name << " = [\n";
    for (int i = -num; i <= num; ++i) {
      real u = i/real(ndiv),
        f = (*this)(u),
        f1 = deriv(u),
        df = ((*this)(u + ds/2) - (*this)(u - ds/2)) / ds,
        u0 = inv0(f),
        u1 = inv1(f),
        u2 = inv2(f),
        gf1 = gfderiv(u),
        phi = rev(u),
        ux = fwd(phi);
      if (_distp)
        os << u << " " << f << " "
           << u0 << " " << u1 << " " << u2 << " "
           << u0-u << " " << u1-u << " " << u2-u << " "
           << f1 << " " << df << " " << df-f1 << " " << gf1 << ";\n";
      else
        os << u << " " << f << " "
           << u0 << " " << u1 << " " << u2 << " "
           << u0-u << " " << u1-u << " " << u2-u << " "
           << f1 << " " << df << " " << df-f1 << " "
           << phi << " " << ux << " " << ux-u << ";\n";
    }
    os << "];\n";
  }

} // namespace GeographicLib

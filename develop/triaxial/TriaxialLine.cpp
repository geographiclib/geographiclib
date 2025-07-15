/**
 * \file TriaxialLine.cpp
 * \brief Implementation for GeographicLib::TriaxialLine class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

// Set _oblpro = false, _merid = true

#include "TriaxialLine.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>

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

  TriaxialLine::TriaxialLine(const Triaxial& t)
    : _f(t, Triaxial::gamblk(t, (t._umbalt && t._kp2 > 0) || t._k2 == 0))
  {
    // Not worth it...
    // _f.ComputeInverse();
  }

  TriaxialLine::TriaxialLine(const Triaxial& t,
                             Angle bet1, Angle omg1, Angle alp1)
    : _t(t)
  {
    bet1.round();
    omg1.round();
    alp1.round();
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
    _fic.pos1(_f.transpolar(), bet1, omg1, alp1);
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
                              int* countn, int* countb) const {
    // Compute points at distance s12
    // cout << "XX " << _fic.delta << " " << _gic.sig1 << "\n";
    real sig2 = _gic.sig1 + s12/_t._b;
    Angle &phi2a = bet2a, &tht2a = omg2a;
    if (_f.gammax() > 0) {
      real u2, v2;
      // The general triaxial machinery.  This is used for non-meridional
      // geodesics on biaxial ellipsoids.
      solve2(_fic.delta, sig2, fpsi(), ftht(), gpsi(), gtht(), u2, v2,
             countn, countb);
      tht2a = ang::radians(ftht().rev(v2));
      ang psi2 = ang::radians(fpsi().rev(u2));
      // Already normalized
      phi2a = ang(_f.gm().nup * psi2.s(),
                  _fic.phi0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                  0, true).rebase(_fic.phi0);
      alp2a = ang(_fic.Ex * hypot(_f.kx() * _f.gm().nu, _f.kxp() * tht2a.c()),
                  _fic.phi0.c() * _f.kx() * _f.gm().nup * psi2.c());
      if (_t.debug())
        cout << real(tht2a) << " " << real(psi2) << " ";
    } else if (_f.gammax() == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_g.s0);  // reduce to [-s0, s0)
      real u2, v2;
      ang psi2;
      if (_f.kxp2() == 0) {
        // gtht()(x) == 0, so the g equation becomes gpsi()(v2) = sig2
        v2 = gpsi().inv(sig2n.first);
        phi2a = ang::radians(v2);
        psi2 = phi2a + ang::cardinal(2 * sig2n.second);
        int parity = fmod(sig2n.second, real(2)) != 0 ? -1 : 1;
        int Ny = _fic.Nx * parity;
        // cout << "ZZ " << sig2 << " " << sig2n.first << " " << sig2n.second << " " << Ny << "\n";
        phi2a = phi2a.reflect(signbit(_fic.phi0.c() * Ny),
                              signbit(_fic.phi0.c())).rebase(_fic.phi0);
        // ftht().inv(x) is just x, but keep it general
        // u2 = ftht().inv(fpsi()(v2) - _fic.delta);
        tht2a = ang::radians(- _fic.delta) + ang::cardinal(2 * sig2n.second);
        // cout << "YY " << sig2 << " "
        // << sig2n.first << " " << sig2n.second << "\n";
        alp2a = ang(_fic.Ex * real(0), _fic.Nx * parity, 0, true);
        /*
        ang psi2 = ang::radians(fpsi().rev(v2));
        // Already normalized
        phi2a = ang(_f.gm().nup * psi2.s(),
                    _fic.phi0.c() * hypot(psi2.c(), _f.gm().nu * psi2.s()),
                    0, true).rebase(_fic.phi0);
        alp2a = ang(_fic.Ex * hypot(_f.kx() * _f.gm().nu, _f.kxp() * tht2a.c()),
                    _fic.phi0.c() * _f.kx() * _f.gm().nup * psi2.c())
          .rebase(_fic.alp0);
        cerr << "FOO " << sig2 << " " << _g.s0 << "\n";
        cerr << "FOO " << sig2n.first << " " << sig2n.second << "\n";
        */
      } else {
        if (sig2n.first - _g.s0 >= -5 * numeric_limits<real>::epsilon()) {
          sig2n.first = -_g.s0;
          ++sig2n.second;
        }
        real deltax = clamp(_fic.delta + sig2n.second * _f.deltashift(), 1);
        solve2u(deltax, sig2n.first, fpsi(), ftht(), gpsi(), gtht(), u2, v2,
                countn, countb);
        // phi2 = fpsi().rev(u2); tht2 = ftht().rev(v2);
        phi2a = anglam(u2, _f.kxp());
        psi2 = phi2a + ang::cardinal(2 * sig2n.second);
        tht2a = anglam(v2, _f.kx());
        int parity = fmod(sig2n.second, real(2)) != 0 ? -1 : 1;
        // if t._kp2 == 0 then meridional oblate
        int Ny = _fic.Nx * parity;
        tht2a += ang::cardinal(2 * sig2n.second);
        tht2a = tht2a + _fic.tht0;
        phi2a = phi2a.reflect(signbit(_fic.phi0.c() * Ny),
                              signbit(_fic.phi0.c())).rebase(_fic.phi0);
        // replace cos(phi)/cos(tht) by sech(u)/sech(v)
        alp2a = ang(_fic.Ex * _f.kxp() / mcosh(v2, _f.kx()),
                    _f.kx() * Ny / mcosh(u2, _f.kxp()));
      }
      if (_t.debug())
        cout << real(tht2a) << " " << real(psi2) << " ";
    } else {
      // gamma = NaN
    }
    phi2a.round();
    tht2a.round();
    alp2a.round();
    tht2a = tht2a.flipsign(_fic.Ex);
    if (_f.transpolar()) {
      swap(bet2a, omg2a);
      alp2a.reflect(false, false, true);
    }
    alp2a = alp2a.rebase(_fic.alp0);
    omg2a += ang::cardinal(1);
    if (0) {
      // Angle::rebase debug
      ang
        p1 = ang(+0.0, -1, 0, true),
        p2 = ang(-0.0, -1, 1, true),
        m1 = ang(-0.0, -1, 0, true),
        m2 = ang(+0.0, -1, -1, true);
      cout << "DD " << real(p1) << " " << real(p2) << " "
           << real(m1) << " " << real(m2) << "\n";
      cout << "PP "
           << real(p1.rebase(p1)) << " "
           << real(p2.rebase(p1)) << " " // bad
           << real(m1.rebase(p1)) << " " // bad
           << real(m2.rebase(p1)) << " "
           << real(p1.rebase(p2)) << " "
           << real(p2.rebase(p2)) << " "
           << real(m1.rebase(p2)) << " "
           << real(m2.rebase(p2)) << "\n";
      cout << "MM "
           << real(p1.rebase(m1)) << " "
           << real(p2.rebase(m1)) << " "
           << real(m1.rebase(m1)) << " "
           << real(m2.rebase(m1)) << " "
           << real(p1.rebase(m2)) << " "
           << real(p2.rebase(m2)) << " " // bad
           << real(m1.rebase(m2)) << " " // bad
           << real(m2.rebase(m2)) << "\n";
    }
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
                            const hfun& fx, const hfun& fy,
                            const hfun& gx, const hfun& gy,
                            real& x, real& y,
                            int* countn, int* countb) {
    // In the biaxial limit gtht()(x) == 0 and gtht.inv is ill-defined, so x
    // should the tht and y should be psi

    // Return x and y, s.t.
    //   fx(x) - fy(y) = f0
    //   gx(x) + gy(y) = g0
    //
    // We use x as the control variable and assume that gy.inv() is available
    // Given a guess for x, compute y = gy.inv(g0 - gx(x)) and solve the 1d
    // problem:
    //
    //   fx(x) - fy(gy.inv(g0 - gx(x))) - f0 = 0
    //
    // fx(x) = fxs*x +/- fxm,
    // fy(y) = fys*y +/- fym,
    // gx(x) = gxs*x +/- gxm,
    // gy(y) = gys*y +/- gym;
    real fxm = fx.MaxPlus(), fym = fy.MaxPlus(),
      gxm = gx.MaxPlus(), gym = gy.MaxPlus(),
      fxs = fx.Slope(), fys = fy.Slope(), gxs = gx.Slope(), gys = gy.Slope(),
      // solve
      //   x = ( fys*g0 + gys*f0 ) / den +/- Dx
      //   y = ( fxs*g0 - gxs*f0 ) / den +/- Dy
      // where
      den = fxs * gys + fys * gxs, // den > 0
      qf = fxm + fym, qg = gxm + gym,
      Dx = (qf * gys + qg * fys) / den,
      Dy = (qf * gxs + qg * fxs) / den,
      x0 = (fys * g0 + gys * f0) / den, // Initial guess
      y0 = (fxs * g0 - gxs * f0) / den,
      xp = x0 + Dx, xm = x0 - Dx,
      yp = y0 + Dy, ym = y0 - Dy,
      mm = 0;
    if (1) {
      real DxA = (-qf * gys + qg * fys) / den,
        DyA = (-qf * gxs + qg * fxs) / den;
      real
        fA = fx(x0 - Dx) - fy(y0 - DyA) - f0,
        gA = gx(x0 - Dx) + gy(y0 - DyA) - g0,
        fB = fx(x0 - DxA) - fy(y0 - Dy) - f0,
        gB = gx(x0 - DxA) + gy(y0 - Dy) - g0,
        fC = fx(x0 + Dx) - fy(y0 + DyA) - f0,
        gC = gx(x0 + Dx) + gy(y0 + DyA) - g0,
        fD = fx(x0 + DxA) - fy(y0 + Dy) - f0,
        gD = gx(x0 + DxA) + gy(y0 + Dy) - g0;
      if (!( fabs(DxA) <= Dx && fabs(DyA) <= Dy ))
        throw GeographicLib::GeographicErr("Bad Dx/Dy");
      if (!( fA <= 0 && gA <= 0 )) {
        cout << scientific << "midA " << fA << " " << gA << "\n";
        cout << "DA " <<  Dx << " " << DxA << " " << Dy << " " << DyA << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints A TriaxialLine::newt2");
      }
      if (!( fB >= 0 && gB <= 0 )) {
        cout << scientific << "midB " << fB << " " << gB << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints B TriaxialLine::newt2");
      }
      if (!( fC >= 0 && gC >= 0 )) {
        cout << scientific << "midC " << fC << " " << gC << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints C TriaxialLine::newt2");
      }
      //      cout << scientific << "midD " << fD << " " << gD << "\n";
      if (!( fD <= 0 && gD >= 0 )) {
        cout << scientific << "midD " << fD << " " << gD << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints B TriaxialLine::newt2");
      }
      }
    newt2(f0, g0, fx, fy, gx, gy,
          xm-mm*Dx, xp+mm*Dx, fx.HalfPeriod(),
          ym-mm*Dy, yp+mm*Dy, fy.HalfPeriod(),
          (fx.HalfPeriod() * fxs + fy.HalfPeriod() * fys) / 2,
          (gx.HalfPeriod() * gxs + gy.HalfPeriod() * gys) / 2,
          x, y, countn, countb);
    if (0)
      cout << "FEQ " << fx(x) << " " << fy(y) << " " << f0 << " "
           << fx(x) - fy(y) - f0 << "\n"
           << "GEQ " << gx(x) << " " << gy(y) << " " << g0 << " "
           << gx(x) + gy(y) - g0 << "\n";
    if (0) {
      cout << "FF "
           << fx.HalfPeriod() << " " << fx.Slope() << " " << fx.Max() << " "
           << fy.HalfPeriod() << " " << fy.Slope() << " " << fy.Max() << "\n";
      cout << "GG "
           << gx.HalfPeriod() << " " << gx.Slope() << " " << gx.Max() << " "
           << gy.HalfPeriod() << " " << gy.Slope() << " " << gy.Max() << "\n";
    }
  }

  void TriaxialLine::solve2u(real d0, real s0,
                             const hfun& fx, const hfun& fy,
                             const hfun& gx, const hfun& gy,
                             real& u, real& v,
                             int* countn, int* countb) {
    // Return u and v, s.t.
    // fx(u) - fy(v) = d0
    // gx(u) + gy(v) = s0
    // specialized for umbilics
    //
    // fx, fy, gx, gy are increasing functions defined in [-1, 1]*pi2
    // Assume fx(0) = fy(0) = gx(0) = gy(0) = 0
    real pi2 = Triaxial::BigValue(),
      sbet = gx.Max(), somg = gy.Max(), stot = sbet + somg,
      dbet = fx.Max(), domg = fy.Max(), del  = dbet - domg;
    bool debug = false;
    if (debug)
      cout << "HEREQ "
           << pi2 << " " << d0 << " " << s0 << " "
           << stot << " " << sbet-somg << " "
           << fabs(s0) - stot << " "
           << fabs((1 - 2 * signbit(d0)) * s0 - (sbet - somg)) << "\n";
    if (fabs(s0) - stot >= -5 * numeric_limits<real>::epsilon()) {
      // close to umbilic points we have
      // fx(u) = u -/+ dbet
      // fy(v) = v -/+ domg
      // fx(u) - fy(v) = d0
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
    } else {
      real mm = 2;
      newt2(d0, s0, fx, fy, gx, gy,
            -mm * pi2, mm * pi2, pi2,
            -mm * pi2, mm * pi2, pi2,
            pi2, pi2, u, v, countn, countb);
      if (debug) cout << "UV2D " << u << " " << v << "\n";
    }
  }

  void TriaxialLine::newt2(real f0, real g0,
                           const hfun& fx, const hfun& fy,
                           const hfun& gx, const hfun& gy,
                           real xa, real xb, real xscale,
                           real ya, real yb, real yscale,
                           real fscale, real gscale,
                           real& x, real& y,
                           int* countn, int* countb) {
    // solve
    //   f = fx(x) - fy(y) - f0 = 0
    //   g = gx(x) + gy(y) - g0 = 0
    // for x, y.  Assume that
    //   fx and fy are increasing functions
    //   gx and gy are non-decreasing functions
    // The solution is bracketed by x in [xa, xb], y in [ya, yb]
    const bool debug = false, check = true;
    const int maxit = 150;
    const real tol = numeric_limits<real>::epsilon(),
      // Relax [fg]tol to /10 instead of /100.  Otherwise solution resorts to
      // too much bisection.  Example:
      //   echo 63 -173 -61.97921997838416712 -4.64409746197940890408 |
      //     ./Geod3Solve $PROX
      ftol = tol * fscale/10,
      gtol = tol * gscale/10,
      xtol = pow(tol, real(0.75)) * xscale,
      ytol = pow(tol, real(0.75)) * yscale,
      // If tolmult == 0, umbilical case
      //   echo 70 0 180 0.46138515584816834054 | ./Geod3Solve $SET
      // take countn/countb = 55/54 iterations to converge.  With tolmult == 1,
      // it needs 5/3 iterations.
      tolmult = 1;
    int cntn = 0, cntb = 0;
    real oldf = Math::infinity(), oldg = oldf, olddx = oldf, olddy = oldf;
    zset xset(zvals(xa, fx(xa), gx(xa)),
              zvals(xb, fx(xb), gx(xb))),
      yset(zvals(ya, fy(ya), gy(ya)),
           zvals(yb, fy(yb), gy(yb)));
    // x = xset.bisect(), y = yset.bisect();
    auto p = zsetsbisect(xset, yset, f0, g0);
    x = p.first; y = p.second;
    if (check) {
      // A necessary condition for a root is
      //   f01 <= 0 <= f10
      //   g00 <= 0 <= g11
      real
        f01 = xset.min().fz - yset.max().fz - f0,
        f10 = xset.max().fz - yset.min().fz - f0,
        g00 = xset.min().gz + yset.min().gz - g0,
        g11 = xset.max().gz + yset.max().gz - g0;
      // Allow equality on the initial points
      if (!( f01 <= 0 && 0 <= f10 && g00 <= 0 && 0 <= g11 ))
        throw GeographicLib::GeographicErr
          ("Bad initial points TriaxialLine::newt2");
      zvals xv(x, fx(x), gx(x)), yv(y, fy(y), gy(y));
      if (0) {
        real
          fA = xset.min().fz - yv.fz - f0, gA = xset.min().gz + yv.gz - g0,
          fC = xset.max().fz - yv.fz - f0, gC = xset.max().gz + yv.gz - g0,
          fB = xv.fz - yset.min().fz - f0, gB = xv.gz + yset.min().gz - g0,
          fD = xv.fz - yset.max().fz - f0, gD = xv.gz + yset.max().gz - g0;
        // y=-x            y=x
        // g=0             f=0
        //   \    ind=-2   /
        //    \    f<0    /
        //     \   g>0   /
        //      \  yb   /
        //       \  D  /
        // ind=-4 \   / ind=4
        //   f<0   \ /   f>0
        //   g<0    X    g>0
        //   xa    / \   xb
        //    A   /   \   C
        //       /  B  \.
        //      /  f>0  \.
        //     /   g<0   \.
        //    /    ya     \.
        //   /    ind=2    \.
        // Allow equality on the initial midpoints points
        if (!( fA <= 0 && gA <= 0 )) {
          cout << scientific << "midA " << fA << " " << gA << "\n";
          throw GeographicLib::GeographicErr
            ("Bad initial midpoints A TriaxialLine::newt2");
        }
        if (!( fB >= 0 && gB <= 0 )) {
          cout << scientific << "midB " << fB << " " << gB << "\n";
          throw GeographicLib::GeographicErr
            ("Bad initial midpoints B TriaxialLine::newt2");
        }
        if (!( fC >= 0 && gC >= 0 )) {
          cout << scientific << "midC " << fC << " " << gC << "\n";
          throw GeographicLib::GeographicErr
            ("Bad initial midpoints C TriaxialLine::newt2");
        }
        //      cout << scientific << "midD " << fD << " " << gD << "\n";
        if (!( fD <= 0 && gD >= 0 )) {
          cout << "XY " << xv.z << " " << yset.max().z << "\n";
          cout << scientific << "midD " << fD << " " << gD << "\n";
          throw GeographicLib::GeographicErr
            ("Bad initial midpoints B TriaxialLine::newt2");
        }
      }
    }
    bool bis = false;
    int ibis = -1, i = 0;
    for (; i < maxit ||
           (throw GeographicLib::GeographicErr
            ("Convergence failure TriaxialLine::newt2"), false)
           || GEOGRAPHICLIB_PANIC("Convergence failure Trigfun::root"); ++i) {
      ++cntn;
      zvals xv(x, fx(x), gx(x)), yv(y, fy(y), gy(y));
      // zsetsinsert updates xv and yv to enforce monotonicity of f and g
      zsetsinsert(xset, yset, xv, yv, f0, g0);
      real f = xv.fz - yv.fz - f0, g = xv.gz + yv.gz - g0;
      if ((fabs(f) <= ftol && fabs(g) <= gtol) || isnan(f) || isnan(g)) {
        if (debug)
          cout << "break0 " << scientific << f << " " << g << "\n";
        break;
      }
      // Update bounds as follows
      //
      // y=-x            y=x
      // g=0             f=0
      //   \    ind=-2   /
      //    \    f<0    /
      //     \   g>0   /
      //      \  yb   /
      //       \  D  /
      // ind=-4 \   / ind=4
      //   f<0   \ /   f>0
      //   g<0    X    g>0
      //   xa    / \   xb
      //    A   /   \   C
      //       /  B  \.
      //      /  f>0  \.
      //     /   g<0   \.
      //    /    ya     \.
      //   /    ind=2    \.
      //
      // Secant method: data at points A, B, C, D allow piecewise lineear
      // approximations to fx, fy, gx, gy.
      //
      // Sort x(A), x(B), x(C), x(D).  We need to approximate fx, gx in the
      // interval [x(A), x(C)].  x(B) and X(D) can be arbitrarily ordered
      // relative to [x(A), x(C)] to give:
      //
      // x(A), x(C)
      // x(A), x(B), x(C)
      // x(A), x(D), x(C)
      // x(A), x(B), x(D), x(C)
      // x(A), x(D), x(B), x(C)
      //
      // Thus piecewise linear approximations to fx and gx consist of nx = 1,
      // 2, or 3 pieces.
      //
      // Similarly the approximations to fy and gy consist of ny = 1, 2, or 3
      // pieces.
      //
      // Find secant solution to
      //
      //  fx(x) - fy(y) - f0 = 0
      //  gx(x) + gy(y) - g0 = 0
      //
      // by solving the nx x ny systems of linear equations and accepting the
      // solution were [x, y] are in the corresponing intervals.

      // given x in [x(A), x(C)], test x(B), x(D) to find the tightest bracket
      // x in [x(P), x(Q)].  Linearly interpolate fx using fx(x(P)), fx(x(Q))
      // and similarly for gx.
      //
      // Procedure is similar for fy and gy.

      // Rethink using f = fx(x) - fy(y) - f0, g = gx(x) + gy(y) - g0, i.e.,
      // sums of x-dependent and y-dependent terms.  Let

      // fx:fx0 + (fx1-fx0)*(x-x0)/(x1-x0);
      // fy:fy0 + (fy1-fy0)*(y-y0)/(y1-y0);
      // gx:gx0 + (gx1-gx0)*(x-x0)/(x1-x0);
      // gy:gy0 + (gy1-gy0)*(y-y0)/(y1-y0);
      // solve([fx-fy-f0, gx+gy-g0],[x,y])

      // Positions of the roots
      //   [x-x0 = (-(fx0-fy0-f0)*gyp-(gx0+gy0-g0)*fyp)/(fxp*gyp+fyp*gxp),
      //    y-y0 = ( (fx0-fy0-f0)*gxp-(gx0+gy0-g0)*fxp)/(fxp*gyp+fyp*gxp)];
      //   [x-x1 = (-(fx1-fy1-f0)*gyp-(gx1+gy1-g0)*fyp)/(fxp*gyp+fyp*gxp),
      //    y-y1 = ( (fx1-fy1-f0)*gxp-(gx1+gy1-g0)*fxp)/(fxp*gyp+fyp*gxp)];
      //   [x-x1 = (-(fx1-fy0-f0)*gyp-(gx1+gy0-g0)*fyp)/(fxp*gyp+fyp*gxp),
      //    y-y0 = ( (fx1-fy0-f0)*gxp-(gx1+gy0-g0)*fxp)/(fxp*gyp+fyp*gxp)];
      //   [x-x0 = (-(fx0-fy1-f0)*gyp-(gx0+gy1-g0)*fyp)/(fxp*gyp+fyp*gxp),
      //    y-y1 = ( (fx0-fy1-f0)*gxp-(gx0+gy1-g0)*fxp)/(fxp*gyp+fyp*gxp)];

      // Conditions for root to be in the [x0,x1] x [y0,y1] rectangle:
      //   x-x0 > 0
      //   -f00*gyp/fyp-g00 = -f01*gyp/fyp-g01 > 0
      //   y-y0 > 0
      //    f00*gxp/fxp-g00 = -f10*gxy/fxy-g10 > 0
      //   x-x1 < 0
      //   -f10*gyp/fyp-g10 = -f11*gyp/fyp-g11 < 0
      //   y-y1 < 0
      //    f01*gxp/fxp-g01 =  f11*gxp/fxp-g11 < 0
      //
      // where
      //   f00 = fx0-fy0-f0
      //   f01 = fx0-fy1-f0
      //   f10 = fx1-fy0-f0
      //   f11 = fx1+fy1-f0
      //   g00 = gx0+gy0-g0
      //   g01 = gx0+gy1-g0
      //   g10 = gx1+gy0-g0
      //   g11 = gx1+gy1-g0
      //   fxp = (fx1-fx0)/(x1-x0)
      //   fyp = (fy1-fy0)/(y1-y0)
      //   gxp = (gx1-gx0)/(x1-x0)
      //   gyp = (gy1-gy0)/(y1-y0)

      // A necessary condition for a root is
      //   f01 <= 0 <= f10
      //   g00 <= 0 <= g11
      //
      // The condition g00 < g11, f01 < f10 holds for all rectangles.  So can
      // use condition g00 > 0 or f01 > 0 to exit x loop.

      // Find rectangle with x in [x0, x1], y in [y0, y1].  New test point x =
      // (x0 + x1)/2, y = (y0+y1)/2.  But note the actual solution may not be
      // in [x0, x1] x [y0, y1].

      // Maintain a sorted list of x, y values with x in [xa, xb], y in [ya,
      // yb].

      // zvals structure to hold [z, fz(z), gz(z)] (z = x or y).  zset holds a
      // set of zvals with operator< operating on z alone.  Additional data
      // members: za, zb (so we don't bother removing elements from the set)

      // Static function zbisect takes two zsets and returns bisection point of
      // x-y rectangle contains the solution for piecewise linear
      // approximation.

      real
        fxp = fx.deriv(x), fyp = fy.deriv(y),
        gfxp = gx.gfderiv(x), gfyp =  gy.gfderiv(y),
        den = gfxp + gfyp,
        dx = -( gfyp * f + g) / (fxp * den),
        dy = -(-gfxp * f + g) / (fyp * den),
        xn = x + dx, yn = y + dy,
        dxa = 0, dya = 0,
        xa = xset.min().z, xb = xset.max().z,
        ya = yset.min().z, yb = yset.max().z;
      if (check) {
        if (!( fxp > 0 && fyp > 0 && gfxp >= 0 && gfyp >= 0 && den > 0 )) {
          cout << "DERIVS " << x << " " << y << " "
               << fxp << " " << fyp << " "
               << gfxp << " " << gfyp << "\n";
          throw GeographicLib::GeographicErr
            ("Bad derivatives TriaxialLine::newt2");
        }
      }
      bool cond1 = i < ibis + 10 ||
        ((2*fabs(f) < oldf || 2*fabs(g) < oldg) ||
         (2*fabs(dx) < olddx || 2*fabs(dy) < olddy)),
        cond2 = xn >= xa-xtol*tolmult && xn <= xb+xtol*tolmult &&
        yn >= ya-ytol*tolmult && yn <= yb+ytol*tolmult;
      if (cond1 && cond2) {
        oldf = fabs(f); oldg = fabs(g); olddx = fabs(dx); olddy = fabs(dy);
        x = xn; y = yn;
        bis = false;
        if (!(fabs(dx) > xtol || fabs(dy) > ytol) /* && i > ibis + 2 */) {
          if (debug)
            cout << "break1 " << scientific << dx << " " << dy << "\n";
          break;
        }
      } else {
        // xn = xset.bisect(); yn = yset.bisect();
        p = zsetsbisect(xset, yset, f0, g0);
        xn = p.first; yn = p.second;
        ++cntb;
        if (x == xn && y == yn) {
          if (debug)
            cout << "break2\n";
          break;
        }
        dxa = xn - x; dya = yn - y;
        x = xn; y = yn;
        bis = true;
        ibis = i;
      }
      (void) bis;
      if (debug)
        cout << "AA " << scientific << setprecision(4)
             << x-xa << " " << xb-x << " "
             << y-ya << " " << yb-y << "\n";
      if (debug)
        cout << "CC " << i << " "
             << bis << " " << cond1 << " " << cond2 << " "
             << scientific << setprecision(2) << f << " " << g << " "
             << dx << " " << dy << " " << dxa << " " << dya << " "
             << xb-xa << " " << yb-ya << "\n";
      if (debug)
        cout << "BOX " << i << " " << bis << " "
             << xset.num() << " " << yset.num() << " "
             << scientific << setprecision(3)
             << xset.max().z - xset.min().z << " "
             << yset.max().z - yset.min().z << "\n";
      // cout << "BOX " << xset.num() << " " << yset.num() << "\n";
    }
    if (countn)
      *countn += cntn;
    if (countb)
      *countb += cntb;
    // cout << "CNT " << cntn << " " << cntb << "\n";
  }

  int TriaxialLine::zset::insert(zvals& t, int flag) {
    // Inset t into list.  flag = -/+ 1 indicates new min/max.
    // Return -1 if t was already present; othersize return index of newly
    // inserted value.
    // Value of t.fz and t.gz is adjusted to ensuire monotonicity:
    // if a.z == b.z then a.fz == b.gz && a.fz == b.gz
    // if a.z < b.z then a.fz <= b.gz && a.fz <= b.gz
    using std::isnan;
    int ind = -1;
    if (isnan(t.z)) return ind;
    if (t < min()) {
      t.fz = fmin(t.fz, min().fz);
      t.gz = fmin(t.gz, min().gz);
    } else if (t == min()) {
      // Check if t is "other" endpoint and collapse bracket to zero
      if (flag > 0) _s.resize(1);
      t = min();
    } else if (t == max()) {
      // Check if t is "other" endpoint and collapse bracket to zero
      if (flag < 0) { _s[0] = _s.back(); _s.resize(1); }
      t = max();
    } else if (max() < t) {
      t.fz = fmax(t.fz, max().fz);
      t.gz = fmax(t.gz, max().gz);
    }
    if (!(min() < t && t < max())) // Not in range
      return ind;
    // Now min() < t < max()
    auto p = std::lower_bound(_s.begin(), _s.end(), t);
    if (p == _s.end()) return ind; // Can't happen
    // Fix components of t
    if (*p == t)                   // z components match
      t = *p;                      // set fz and gz values
    else {
      t.fz = Math::clamp(t.fz, (p-1)->fz, p->fz);
      t.gz = Math::clamp(t.gz, (p-1)->gz, p->gz);
    }
    if (flag < 0) {
      _s.erase(_s.begin(), p);
      if (!(*p == t)) {
        _s.insert(_s.begin(), t);
        ind = 0;
      }
    } else if (flag > 0) {
      if (!(*p == t)) {
        _s.erase(p, _s.end());
        _s.push_back(t);
        ind = _s.size() - 1;
      } else
        _s.erase(p+1, _s.end());
    } else if (!(*p == t)) {
      ind = p - _s.begin();
      _s.insert(p, t);
    }
    // else it's a duplicate and not a new end value
    return ind;
  }

  void TriaxialLine::zsetsinsert(zset& xset, zset& yset,
                                 zvals& xfg, zvals& yfg,
                                 real f0, real g0) {
    bool debug = false;
    if (debug)
      cout << "BOX0 " << xfg.z << " " << yfg.z << "\n";
    int xind = xset.insert(xfg), yind = yset.insert(yfg);
    if (debug) {
      cout << "BOXA " << xset.num() << " " << yset.num() << " "
           << xind << " " << yind << " "
           << xset.min().z << " " << xset.max().z << " "
           << yset.min().z << " " << yset.max().z << " "
           << xset.max().z - xset.min().z << " "
           << yset.max().z - yset.min().z << "\n";
      zsetsdiag(xset, yset, f0, g0);
    }
    if (xind < 0 && yind < 0) return;
    zvals xa = xset.min(), xb = xset.max(),
      ya = yset.min(), yb = yset.max();
    for (int i = 0; i < xset.num(); ++i) {
      const zvals& x = xset.val(i);
      for (int j = 0; j < yset.num(); ++j) {
        if (i == xind || j == yind) {
          const zvals& y = yset.val(j);
          real f = x.fz - y.fz - f0, g = x.gz + y.gz - g0;
          // Update bounds as follows
          //
          // y=-x            y=x
          // g=0             f=0
          //   \    ind=-2   /
          //    \    f<0    /
          //     \   g>0   /
          //      \  yb   /
          //       \  D  /
          // ind=-4 \   / ind=4
          //   f<0   \ /   f>0
          //   g<0    X    g>0
          //   xa    / \   xb
          //    A   /   \   C
          //       /  B  \.
          //      /  f>0  \.
          //     /   g<0   \.
          //    /    ya     \.
          //   /    ind=2    \.
          if (false) {
            // Problem with
            //   echo -11 165 0.865544228453134485 1.70392636779413409327 |
            //     ./Geod3Solve $HU
            // Pattern is
            // BOXF ---+
            // BOXF ---+
            // BOXF ++++

            // BOXG ..++
            // BOXG --..
            // BOXG --..

            // Allowing equality sets
            // xa/yb xa/yb ../yb ../yb
            // xa/.. xa/.. xa/yb xb/ya
            // ../ya ../ya xb/ya xb/ya

            // Disallowing equality sets
            //             ../yb ../yb
            // xa/.. xa/..
            // ../ya ../ya

            if (f <= 0) {
              if (g <= 0 && xa < x) xa = x;
              if (g >= 0 && y < yb) yb = y;
            }
            if (f >= 0) {
              if (g <= 0 && ya < y) ya = y;
              if (g >= 0 && x < xb) xb = x;
            }
          } else {
            if (f < 0) {
              if (g < 0 && xa < x) xa = x;
              if (g > 0 && y < yb) yb = y;
            }
            if (f > 0) {
              if (g < 0 && ya < y) ya = y;
              if (g > 0 && x < xb) xb = x;
            }
          }
        }
      }
    }
    xset.insert(xa, -1);
    xset.insert(xb, +1);
    yset.insert(ya, -1);
    yset.insert(yb, +1);
    if (debug) {
      cout << "BOXB " << xset.num() << " " << yset.num() << " "
           << xset.min().z << " " << xset.max().z << " "
           << yset.min().z << " " << yset.max().z << " "
           << xa.z << " " << xb.z << " "
           << ya.z << " " << yb.z << "\n";
      zsetsdiag(xset, yset, f0, g0);
    }
    real
      f01 = xset.min().fz - yset.max().fz - f0,
      f10 = xset.max().fz - yset.min().fz - f0,
      g00 = xset.min().gz + yset.min().gz - g0,
      g11 = xset.max().gz + yset.max().gz - g0;
    // Allow equality on the initial points
    if (false && !( f01 <= 0 && 0 <= f10 && g00 <= 0 && 0 <= g11 )) {
      cout << f01 << " " << f10 << " "
           << g00 << " " << g11 << "\n";
      throw GeographicLib::GeographicErr
        ("Bad corner point TriaxialLine::zsetsinsert");
    }
  }

  void TriaxialLine::zsetsdiag(const zset& xset, const zset& yset,
                               real f0, real g0) {
    ostringstream fs, gs;
    for (int j = yset.num() - 1; j >= 0; --j) {
      const zvals& y = yset.val(j);
      fs << "BOXF ";
      gs << "BOXG ";
      for (int i = 0; i < xset.num(); ++i) {
        const zvals& x = xset.val(i);
        real f = x.fz - y.fz - f0, g = x.gz + y.gz - g0;
        fs << (f == 0 ? '.' : f < 0 ? '-' : '+');
        gs << (g == 0 ? '.' : g < 0 ? '-' : '+');
      }
      fs << "\n";
      gs << "\n";
    }
    cout << fs.str() << gs.str();
  }

  pair<Math::real, Math::real>
  TriaxialLine::zsetsbisect(const zset& xset, const zset& yset,
                            real f0, real g0) {
    if (true)
      return pair<real, real>(xset.bisect(), yset.bisect());
    else {
      // A necessary condition for a root is
      //   f01 <= 0 <= f10
      //   g00 <= 0 <= g11
      int cnt = 0;
      real xgap = -1, ygap = -1, x = Math::NaN(), y = Math::NaN();
      for (int i = 0; i < max(1, xset.num() - 1); ++i) {
        int i1 = min(i + 1, xset.num() - 1);
        real xgap1 = xset.val(i1).z - xset.val(i).z,
          xmean = (xset.val(i1).z + xset.val(i).z) / 2;
        for (int j = 0; j < max(1, yset.num() - 1); ++j) {
          int j1 = min(j + 1, yset.num() - 1);
          real ygap1 = yset.val(j1).z - yset.val(j).z,
            ymean = (yset.val(j1).z + yset.val(j).z) / 2;
          real
            f01 = xset.val(i).fz - yset.val(j1).fz - f0,
            f10 = xset.val(i1).fz - yset.val(j).fz - f0,
            g00 = xset.val(i).gz + yset.val(j).gz - g0,
            g11 = xset.val(i1).gz + yset.val(j1).gz - g0;
          if (f01 <= 0 && 0 <= f10 && g00 <= 0 && 0 <= g11) {
            ++cnt;
            if (xgap1 > xgap) { xgap = xgap1; x = xmean; }
            if (ygap1 > ygap) { ygap = ygap1; y = ymean; }
          }
        }
      }
      if (cnt == 0)
        throw GeographicLib::GeographicErr
          ("No legal box TriaxialLine::zsetsbisect");
      // cout << "FOO " << cnt << "\n";
      return pair<real, real>(x, y);
    }
  }

  void TriaxialLine::Hybrid(Angle betomg2,
                            Angle& bet2a, Angle& omg2a, Angle& alp2a,
                            real& s12, bool betp)
    const {
    fline::disttx d = _f.Hybrid(_fic, betomg2, bet2a, omg2a, alp2a, betp);
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
  TriaxialLine::fline::Hybrid(const fics& fic, Angle betomg2,
                              Angle& bet2a, Angle& omg2a, Angle& alp2a,
                              bool betp) const {
    // Is the control variable psi or tht?
    bool psip = !transpolar() ? betp : !betp;
    if (!betp) betomg2 -= ang::cardinal(1);
    ang tau12;
    if (false && _t.debug())
      cout << "HHQ " << transpolar() << " " << psip << " " << gammax() << "\n";
    if (psip) {
      ang phi2 = betomg2;
      if (gammax() > 0) {
        real spsi = phi2.s(),
          // In evaluating equivalent expressions, choose the one with minimum
          // cancelation.  Need the 0 + x to convert -0 to +0.  (Note sqrt(-0) =
          // -0 and fmax(+0, -0) may be -0.)
          cpsi = nu() < nup() ?
          (phi2.c() - nu()) * (phi2.c() + nu()) :
          (nup() - phi2.s()) * (nup() + phi2.s());
        // Return nan if geodesic can't reach phi2 -- given by sqrt(neg) -- but
        // allow a little slop.
        cpsi = !(cpsi > -numeric_limits<real>::epsilon()) ? Math::NaN() :
          (signbit(cpsi) ? 0 : sqrt(cpsi));
        // Need Angle(0, 0) to be treated like Angle(0, 1) here.
        ang psi12 = (ang(spsi, cpsi) - fic.psi1).base();
        // convert -180deg to 180deg
        if (signbit(psi12.s()))
         psi12 = ang(0, copysign(real(1), psi12.c()), 0, true);
        tau12 = psi12;
      } else if (gammax() == 0) {
        tau12 = (fic.Nx > 0 ? phi2 - fic.phi1 :
                 phi2 + fic.phi1 + ang::cardinal(2)).base();
      } else
        tau12 = ang::NaN();
    } else {
      ang tht2 = betomg2.flipsign(fic.Ex);
      if (gammax() >= 0) {
        // cout << "GRR " << real(fic.bet1) << " " << real(betomg2) << "\n";
        // FIX THIS!! betp shouldn't be appearing here.
        // Test case
        // echo -88 21 88 -111 | ./Geod3Solve -i $SET --hybridalt
        ang tht2b = tht2; tht2b.reflect(false, betp && fic.Ex < 0);
        ang tht12 = tht2b - fic.tht1;
        // convert -180deg to 180deg
        if (signbit(tht12.s()))
          tht12 = ang(0, copysign(real(1), tht12.c()), 0, true);
        tau12 = tht12;
        if (false && _t.debug())
          cout << "BBB " << real(tht2) << " " << real(tht2b) << " "
               << real(tht12) << " " << real(fic.tht1) <<  "\n";
      } else
        tau12 = ang::NaN();
    }
    if (false && _t.debug())
      cout << "HERE " << transpolar() << " " << psip << " " << betp << " "
         << real(tau12.base()) << "\n";
    disttx ret = ArcPos0(fic, tau12.base(), bet2a, omg2a, alp2a, betp);
    if (false && _t.debug())
      cout << "HH0X " << real(fic.psi1) << " "
           << real(tau12.base()) << " " << real(omg2a) << "\n";
    return ret;
  }

  TriaxialLine::fline::fline(const Triaxial& t, bool neg)
    : _t(t)
    , _gm(t, neg)
  {}

  TriaxialLine::fline::fline(const Triaxial& t, Triaxial::gamblk gam)
    : _t(t)
    , _gm(gam)
    , _fpsi(false, _gm.kx2 , _gm.kxp2,
            +(_gm.transpolar ? -1 : 1) * _t._e2,
            -(_gm.transpolar ? -1 : 1) * _gm.gamma, t)
    , _ftht(false, _gm.kxp2 , _gm.kx2,
            -(_gm.transpolar ? -1 : 1) * _t._e2,
            +(_gm.transpolar ? -1 : 1) * _gm.gamma, t)
    , _invp(false)
  {
    // Only needed for umbilical lines
    _deltashift = _gm.gamma == 0 ?
      (_t.k2() > 0 && _t.kp2() > 0 ? 2 * (_fpsi.Max() - _ftht.Max()) : 0) :
      Math::NaN();
  }

  void TriaxialLine::fline::ComputeInverse() {
    if (!_invp) {
      _fpsi.ComputeInverse();
      _ftht.ComputeInverse();
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
    , _gpsi(true, _gm.kx2 , _gm.kxp2,
            +(_gm.transpolar ? -1 : 1) * _t._e2,
            -(_gm.transpolar ? -1 : 1) * _gm.gamma, t)
    , _gtht(true, _gm.kxp2 , _gm.kx2,
            -(_gm.transpolar ? -1 : 1) * _t._e2,
            +(_gm.transpolar ? -1 : 1) * _gm.gamma, t)
    , _invp(false)
    , s0(_gm.gammax == 0 ? _gpsi.Max() + _gtht.Max() : 0)
  {}

  void TriaxialLine::gline::ComputeInverse() {
    if (!_invp) {
      _gpsi.ComputeInverse();
      _gtht.ComputeInverse();
      _invp = true;
    }
  }

  Math::real TriaxialLine::fline::Hybrid0(const fics& fic,
                                          Angle bet2, Angle omg2,
                                          bool betp) const {
    ang bet2a, omg2a, alp2a;
    (void) Hybrid(fic, betp ? bet2 : omg2, bet2a, omg2a, alp2a, betp);
    if (false && _t.debug())
      cout << "HH1X " << real(bet2) << " " << real(omg2) << " "
           << real(bet2a) << " " << real(omg2a) << "\n";
    bool angnorm = true || betp;
    if (angnorm)
      (void) Triaxial::AngNorm(bet2a, omg2a, alp2a, !betp);
    if (betp) {
      omg2a -= omg2;
      return omg2a.radians0();
    } else {
      bet2a -= bet2;
      return angnorm ? bet2a.radians0() : bet2.radians();
    }
  }

  TriaxialLine::fline::disttx
  TriaxialLine::fline::ArcPos0(const fics& fic, Angle tau12,
                               Angle& bet2a, Angle& omg2a, Angle& alp2a,
                               bool betp) const {
    // XXX fix for biaxial
    disttx ret{Math::NaN(), Math::NaN(), 0};
    bool psip = transpolar() ? !betp : betp;
    Angle &phi2a = bet2a, &tht2a = omg2a;
    if (gammax() > 0) {
      ang psi2;
      real u2, v2, u2x = 0;
      if (psip) {
        psi2 = tau12 + fic.psi1;
        v2 = fpsi().fwd(psi2.radians());
        u2 = ftht().inv(fpsi()(v2) - fic.delta);
        tht2a = ang::radians(ftht().rev(u2));
      } else {
        tht2a = fic.tht1 + tau12;
        u2 = ftht().fwd(tht2a.radians());
        u2x = ftht()(u2) + fic.delta;
        v2 = fpsi().inv(u2x);
        psi2 = ang::radians(fpsi().rev(v2));
        if (false && _t.debug())
          cout << "FOO "
               << real(fic.tht1) << " " << real(tau12) << " "
               << real(tht2a) << " "
               << u2 << " " << u2x << " " << v2 << " "
               << real(psi2) << "\n";
      }
      // Already normalized
      phi2a = ang(nup() * psi2.s(),
                  fic.phi0.c() * hypot(psi2.c(), nu() * psi2.s()),
                  0, true).rebase(fic.phi0);
      if (!transpolar() && kxp2() == 0 && !psip && gammax() == 0)
        alp2a = ang(fic.Ex * fabs(sin(u2x)), fic.phi0.c() * cos(u2x));
      else {
        real s = fic.Ex * hypot(kx() * nu(), kxp() * tht2a.c()),
          c = fic.phi0.c() * kx() * nup() * psi2.c();
        if (s == 0 && c == 0)
          (transpolar() ? s : c) = 1;
        alp2a = ang(s, c);
      }
      ret.phiw2 = v2;
      ret.thtw2 = u2;
    } else if (gammax() == 0) {
      real u2, v2;
      int ii;
      if (psip) {
        phi2a = fic.phi1 + tau12.flipsign(fic.Nx);
        // remainder with result in [-pi/2, pi/2)
        pair<real, real> phi2n =
          // remx(fic.Nx * (phi2a - fic.phi0).radians(), Math::pi());
          remx((phi2a - fic.phi0).flipsign(fic.Nx));
        if (false && t().debug())
          cout << "PHI2 " << real((phi2a - fic.phi0).flipsign(fic.Nx)) << " "
               << phi2n.first << " " << phi2n.second << "\n";
        u2 = fpsi().fwd(phi2n.first);
        int parity = fmod(phi2n.second, real(2)) != 0 ? -1 : 1;
        if (kxp() == 0) {
          // v2 is independent on u2
          v2 = 0;
          tht2a = ang::radians(-fic.Ex * fic.delta) +
            ang::cardinal(2*fic.Ex*phi2n.second);
          alp2a = fic.alp1.nearest(2U) + ang::cardinal(parity < 0 ? 2 : 0);
          /*
          alp2a = ang::cardinal(2 * (phi2n.second +
                                     (fic.phi1.c() == 0 &&
                                      signbit(fic.alp1.c()) ? 1 : 0)))
            .rebase(fic.alp0);
          */
          if (false && t().debug()) {
            cout << "ALP " << real(fic.alp1.nearest(2U)) << " "
                 << parity << " " << real(alp2a) << "\n";
            cout << "FIC "
                 << real(fic.phi0) << " " << real(fic.phi1) << " "
                 << real(fic.tht0) << " " << real(fic.tht1) << " "
                 << real(fic.alp0) << " " << real(fic.alp1) << " "
                 << fic.Nx << " " << fic.Ex << " "
                 << fic.u0 << " " << fic.v0 << " "
                 << fic.delta/Math::pi() << "\n";
            cout << "WWW " << fic.Ex << " " << parity << " "
                 << real(fic.alp0) << " "
                 << real(alp2a)<< "\n";
          }
        } else {
          real deltax = clamp(fic.delta + phi2n.second * _deltashift, 2);
          v2 = ftht().inv(fpsi()(u2) - deltax);
          tht2a = ang::radians(parity * ftht().rev(v2))
            .rebase(fic.tht0);
          // Conflict XXX
          // testset -50 180 20 0 want fic.Nx multiplying s()
          // testspha -20 90 20 -90 wants fic.Nx multiplying c()
          alp2a = ang(kxp() * fic.Ex * parity / mcosh(v2, kx()),
                      fic.Nx * kx() / mcosh(u2, kxp()));
          // Move forward from umbilical point
          phi2a += ang::eps().flipsign(fic.Nx);
        }
        ii = int(phi2n.second);
      } else {
        tht2a = fic.tht1 + tau12;
        // remainder with result in [-pi/2, pi/2)
        pair<real, real> tht2n =
          // remx(fic.Ex * (tht2a - fic.tht0).radians(), Math::pi());
          remx(tht2a - fic.tht0);
        v2 = ftht().fwd(tht2n.first);
        u2 = 0;
        int parity = fmod(tht2n.second, real(2)) != 0 ? -1 : 1;
        if (kxp() == 0) {
          if (fic.phi1.c() != 0 && tau12 == tau12.nearest(2U)) {
            // STILL TO DO...
            // Special case for finding conjugate points on meridional
            // geodesics, gamma = 0, tau12 = multiple of pi.  This is
            // specialized for prolate ellipsoids for now.
            real npi = (tau12.ncardinal() + fic.phi1.nearest(2U).ncardinal())
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
            // [sx,a1x,a2x] = t.distance(p1,[90,178.9293750483])
            //    -> a1x = 90
            // [sx,a1x,a2x] = t.distance(p1,[90,178.9293750484])
            //    -> a1x = 90.0023
            //   omg1 = -91 -> omg2 = 178.9293750483-90 = 88.9293750483
            real c = fic.Nx * (fic.phi1.t() - fpsi().df(fic.phi1.radians())),
              l = exp(Triaxial::BigValue()),
              tpsi2 = Trigfun::root([this, npi] (real tpsi) -> pair<real, real>
                                    {
                                      real psi = atan(tpsi);
                                      return pair<real, real>
                                        (tpsi - fpsi().df(npi + psi),
                                         1 - fpsi().dfp(psi) /
                                         (1 + Math::sq(tpsi)));
                                    },
                                    c, fic.psi1.t(), -l, l, 1, 1, 1,
                                    nullptr, nullptr, 0, Trigfun::ARCPOS0);
            v2 = atan(tpsi2);
            phi2a = (ang(tpsi2, 1) + fic.phi1.nearest(2U))
              .flipsign(parity*fic.Nx).rebase(fic.phi0);
            if (false && t().debug())
              cout << "PSI2 " << real(fic.phi1) << " " << real(tau12) << " " << npi << " " << c << " " << real(phi2a) << " " << parity << "\n";
            alp2a = ang::cardinal(fic.Nx * parity == 1 ? 0 : 2);
          } else {
            u2 = v2 == 0 ? 0 : copysign(Math::pi()/2, tht2n.first);
            if (false && t().debug()) {
              cout << "FIC "
                   << real(fic.phi0) << " " << real(fic.phi1) << " "
                   << real(fic.tht0) << " " << real(fic.tht1) << " "
                   << real(fic.alp0) << " " << real(fic.alp1) << " "
                   << real(fic.alp1.nearest(2U)) << " "
                   << fic.Nx << " " << fic.Ex << " "
                   << fic.u0 << " " << fic.v0 << " "
                   << fic.delta/Math::pi() << "\n";
              cout << "YY " << real(tau12) << "\n";
              cout << "XX " << real(tht2a) << " "
                   << tht2n.first/Math::degree() << " "
                   << tht2n.second << " "
                   << v2/Math::degree() << " "
                   << parity<< " " << u2/Math::degree() << "\n";
            }
            phi2a = ang::cardinal(fabs(v2) == 0
                                  ? 0 : copysign(real(1), v2 * fic.Nx))
              .rebase(fic.phi0);
            alp2a = fabs(v2) == 0 ? ang::cardinal(2 /* * fic.Ex * parity*/ ) :
              fic.alp1.nearest(2U) +
              ang::cardinal(parity == 1 ? 0 : 2) +
              ang::radians(v2).flipsign(parity * phi2a.s());
          }
          if (false && t().debug())
            cout << "HERE " << real(phi2a) << " "
                 << real((tht2a - fic.tht0).flipsign(fic.Ex)) << " "
                 << real(fic.tht0) << " "
                 << real(fic.tht1) << " " << real(tau12) << " "
                 << fic.Ex << " "
                 << v2/Math::degree() << " " << u2/Math::degree() << "\n";
        } else {
          real deltax = clamp(fic.delta + tht2n.second * _deltashift, 2);
          u2 = fpsi().inv(ftht()(v2) + deltax);
          real phi2 = fic.Nx * parity * fpsi().rev(u2);
          phi2a = ang::radians(phi2);
          alp2a = ang(fic.Ex * kxp() / mcosh(v2, kx()),
                      kx() * fic.Nx * parity / mcosh(u2, kxp()));
          tht2a += ang::eps();
        }
        ii = int(tht2n.second);
        // Move forward from umbilical point
      }
      ret.phiw2 = u2;
      ret.thtw2 = v2;
      ret.ind2 = ii;
    } else {
      // gamma == NaN
    }

    tht2a = tht2a.flipsign(fic.Ex);
    if (transpolar()) {
      swap(bet2a, omg2a);
      alp2a.reflect(false, false, true);
    }
    omg2a += ang::cardinal(1);
    alp2a = alp2a.rebase(fic.alp0);
    if (false && t().debug())
      cout << "ARCPOS0 " << psip << " "
           << real(bet2a) << " " << real(omg2a) << " "
           << real(alp2a) << "\n";
    return ret;
  }

  TriaxialLine::fline::fics::fics(const fline& f,
                                  Angle bet10, Angle omg10,
                                  Angle alp10)
    : tht1(omg10 - ang::cardinal(1))
    , phi1(bet10)
    , alp1(alp10)
  {
    alp0 = alp1.nearest(f.transpolar() ? 2U : 1U);
    if (f.transpolar()) {
      swap(tht1, phi1);
      alp1.reflect(false, false, true);
    }
    const real eps = numeric_limits<real>::epsilon();
    const Triaxial& t = f.t();
    if (!f.transpolar() && phi1.s() == 0 && fabs(alp1.c()) <= Math::sq(eps))
      alp1 = ang(alp1.s(), - Math::sq(eps), alp1.n(), true);
    Ex = signbit(alp1.s()) ? -1 : 1;
    Nx = signbit(alp1.c()) ? -1 : 1;
    tht1 = tht1.flipsign(Ex);
    if (f.gammax() > 0) {
      phi0 = phi1.nearest(2U);
      psi1 = ang(f.kx() * phi1.s(),
                 phi0.c() * alp1.c() *
                 hypot(f.kx() * phi1.c(), f.kxp() * tht1.c()));
      // k = 1, kp = 0, sqrt(-mu) = bet1.c() * fabs(alp1.s())
      // modang(psi1, sqrt(-mu)) = atan2(bet1.s() * fabs(alp1.s()),
      //                                 bet0.c() * alp1.c());
      // assume fbet().fwd(x) = x in this case
      v0 = f.fpsi().fwd(psi1.radians());
      u0 = f.ftht().fwd(tht1.radians());
      // Only used for biaxial cases when fwd rev is the identity
      // v0a = psi1;
      // u0a = Ex < 0 ? -omg1 : omg1;
      delta = f.fpsi()(v0) - f.ftht()(u0);
      // deltaa = f.fbet()(v0a) - f.fomg()(u0a);
    } else if (f.gammax() == 0) {
      if (f.kxp2() == 0) {
        // meridonal geodesic on biaxial ellipsoid
        // N.B. factor of sqrt(f.gammax()) omitted
        // Implicitly assume phi1 in [-90, 90] for now

        // phi1, tht1, alp1 are starting conditions.

        // If alp1 = 0/180, alp1.s() == 0, and phi1 != +/-90, phi1.c() != 0,
        // then, at baseline equator crossing, phi = phi1.nearest(2U), tht =
        // tht0 = tht1, alp = alp1
        phi0 = phi1.nearest(2U); // phi0 = 0
        // psi1 not used, perhaps it should be?
        // psi1 = ang(phi1.s(), phi0.c() * alp1.c() * fabs(phi1.c())));
        tht0 = tht1;

        // Othersise at pole, phi1.c() == 0.  Define baseline equator crossing
        // by alp = alp1.nearest(2U),
        // phi = phi1 + cardinal(signbit(alp1.c()) ? -1 : 1)
        // tht0 = tht1 +
        //   (alp1.nearest(2U)-alp1).flipsign(phi1.s())
        if (phi1.c() == 0) {
          // phi0 = phi1 + ang::cardinal(signbit(alp1.c()) ? -1 : 1);
          tht0 += (alp1.nearest(2U)-alp1).flipsign(phi1.s() * Ex);
        }
        v0 = f.fpsi().fwd((phi1-phi0).radians());
        u0 = 0;
        delta = -f.ftht()(tht0.radians());
       if (0)
          cout << "AA " << real(tht1) << " " << real(tht0) << " "
               << real(alp1) << " " << signbit(phi1.s()) << " "
               << real(alp1.flipsign(phi1.s())) << " "
               << delta/Math::degree() << " " << Ex << " " << Nx << "\n";
      } else {
        // N.B. factor of k*kp omitted
        // phi0, tht0 are the middle of the initial umbilical segment
        if (fabs(phi1.c()) < 8*eps && fabs(tht1.c()) < 8*eps) {
          phi0 = phi1.nearest(1U) + ang::cardinal(Nx);
          tht0 = tht1.nearest(1U) + ang::cardinal(1);
          delta = f.deltashift()/2 - log(fabs(alp1.t()));
        } else {
          phi0 = phi1.nearest(2U);
          tht0 = tht1.nearest(2U);
          delta = Nx * f.fpsi()(lamang(phi1 - phi0, t._kp)) -
            f.ftht()(lamang(tht1 - tht0, t._k));
        }
      }
    } else {
      // gamma = NaN
    }
  }

  void TriaxialLine::fline::fics::pos1(bool transpolar,
                                       Angle& bet10, Angle& omg10,
                                       Angle& alp10) const {
    bet10 = phi1; omg10 = tht1.flipsign(Ex);
    alp10 = alp1;
    if (transpolar) {
      swap(bet10, omg10);
      alp10.reflect(false, false, true);
    }
    omg10 += ang::cardinal(1);
  }

  void TriaxialLine::fline::fics::setquadrant(const fline& f, unsigned q) {
    ang bet1, omg1, alp1;
    pos1(f.transpolar(), bet1, omg1, alp1);
    alp1.setquadrant(q);
    *this = fics(f, bet1, omg1, alp1);
  }

  TriaxialLine::gline::gics::gics(const gline& g, const fline::fics& fic) {
    if (g.gammax() > 0) {
      sig1 = g.gpsi()(fic.v0) + g.gtht()(fic.u0);
    } else if (g.gammax() == 0) {
      sig1 = g.kxp2() == 0 ? fic.Nx * g.gpsi()(fic.v0) :
        fic.Nx * g.gpsi()(lamang(fic.phi1 - fic.phi0, g.kxp())) +
        g.gtht()(lamang(fic.tht1 - fic.tht0, g.kx()));
    } else {
      // gamma = NaN
    }
  }

  Math::real TriaxialLine::gline::dist(gics ic, fline::disttx d) const {
    real sig2 = gpsi()(d.phiw2) + gtht()(d.thtw2) + d.ind2 * 2*s0;
    if (false && t().debug()) {
      cout << "FOO " << gtht()(0) << " " << gtht()(Math::pi()/2) << "\n";
      cout << "DDX " << d.phiw2/Math::degree() << " "
           << d.thtw2/Math::degree() << " "
           << gpsi()(d.phiw2) << " " << gtht()(d.thtw2) << " "
           << d.ind2 * 2*s0 << " " << ic.sig1 << "\n"
           << d.ind2 << " " << s0 << "\n";
    }
    return (sig2 - ic.sig1) * t()._b;
  }

  void TriaxialLine::inversedump(ostream& os) const {
    os << "[b, e2, k2, kp2, gam] = deal("
       << _t.b() << ", " << _t.e2() << ", "
       << _t.k2() << ", " << _t.kp2() << ", "
       << gamma() << ");\n";
    os << "tx = ["
       << fpsi().txp() << ", " << gpsi().txp() << ", "
       << ftht().txp() << ", " << gtht().txp() << "];\n" ;
    _f.inversedump(os);
    _g.inversedump(os);
  }

  void TriaxialLine::fline::inversedump(ostream& os) const {
    _fpsi.inversedump(os, "fpsi");
    _ftht.inversedump(os, "ftht");
  }

  void TriaxialLine::gline::inversedump(ostream& os) const {
    _gpsi.inversedump(os, "gpsi");
    _gtht.inversedump(os, "gtht");
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
    , _umb(!t._biaxial && _mu == 0)
    , _meridr(_kap == 0 && _mu == 0)
    , _meridl(_kapp == 0 && _mu == 0)
    , _invp(false)
  {
    // mu in [-kap, kapp], eps in (-inf, 1/kap)
    if (!_distp) {
      if (_meridr) {
        // biaxial rotating coordinate
        // _kapp == 1, mu < 0 not allowed
        _tx = false;
        // f multiplied by sqrt(mu)
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real tht) -> real
                          // This is a trivial case f' = 1
                          { return fthtbiax(tht, eps, mu); },
                          Math::pi()/2, false);
      } else if (_meridl) {
        // biaxial librating coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < t._ellipthresh;
        _tx = false;
        // f explicitly multiplied by sqrt(-mu) in operator()() for _biaxl
        // For _meridl operator() ignores this
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return dfpsibiax(sin(psi), cos(psi), eps, mu); },
                          Math::pi()/2, false);
      } else if (_mu > 0) {
        _tx = _mu / (_kap + _mu) < t._ellipthresh;
        if (_tx) {
          _ell = EllipticFunction(_kap / (_kap + _mu), 0,
                                  _mu / (_kap + _mu), 1);
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
        // Include scale = 1 in TrigfunExt constructor because this function
        // gets added to u.
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
    } else {                    // _distp
      if (_meridr) {
        // biaxial symmetry coordinate
        // _kapp == 1, mu < 0 not allowed
        _tx = false;
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real tht) -> real
                          // degenerate f' = 0
                          { return gthtbiax(tht, eps, mu); },
                          Math::pi()/2, false);
      } else if (_meridl) {
        // biaxial non-symmetry coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < t._ellipthresh;
        _tx = false;
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return gpsibiax(sin(psi), cos(psi), eps, mu); },
                          Math::pi()/2, false);
      } else if (_mu > 0) {
        _tx = _mu / (_kap + _mu) < t._ellipthresh;
        if (_tx) {
          _ell = EllipticFunction(_kap / (_kap + _mu), 0,
                                  _mu / (_kap + _mu), 1);
          _fun = TrigfunExt(
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
      ( _umb ? _fun(_tx ? _ell.K() : Math::pi()/2) : _fun.Max() ) :
      ( _umb ? _fun(_tx ? _ell.K() : Math::pi()/2) :
        (_meridl ? _fun(Math::pi()/2) : _fun.Max()) );
  }

  Math::real TriaxialLine::hfun::operator()(real u) const {
    if (!_distp) {
      if (_meridl)
        return 0;
      else if (_umb) {
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
    real u = ang.radians();
    if (_umb) {
      // This is sqrt(kap * kapp) * f(u)
      real phi = gd(u, _sqrtkapp);
      return ang::radians(u - _fun(_tx ? _ell.F(phi) : phi));
    } else
      return ang::radians(_fun(u));
  }

  Math::real TriaxialLine::hfun::deriv(real u) const {
    if (!_distp) {
      if (_meridl)
        return 0;
      else if (_umb) {
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
    if (_meridr)
      return gfthtbiax(u, _mu);
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
        // _meridl and _meridr cases fall through to here.  In these cases one
        // of g'/f' is singular, so we can't use our general 2-d solution
        // machinery.  This could easily be fixed (by multiplying by one of the
        // f').  However, we'll just do the simple 1d solution for this case.
        return Math::NaN();
    }
  }

  void TriaxialLine::hfun::ComputeInverse() {
    if (!_distp) {
      if (!_invp) {
        if (_umb) {
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
      if (!(_invp || _umb)) {
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
      if (_umb) {
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
      return _umb ? z + _dfinv(gd(z, _sqrtkapp)) : _fun.inv0(z);
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
      return _umb ? root(z, z, countn, countb) : _fun.inv1(z, countn, countb);
    else
      return _umb ? root(z, inv0(z), countn, countb) :
        _fun.inv1(z, countn, countb);
  }

  // Accurate inverse correcting result from _finv
  Math::real TriaxialLine::hfun::inv2(real z, int* countn, int* countb) const {
    if (!_invp) return Math::NaN();
    if (!_distp)
      return _umb ? root(z, inv0(z), countn, countb) :
        _fun.inv2(z, countn, countb);
    else
      return _umb ? root(z, inv0(z), countn, countb) :
        _fun.inv1(z, countn, countb);
  }

  Angle TriaxialLine::hfun::inv(const Angle& z, int* countn, int* countb)
    const {
    if (!_distp) return ang::NaN();
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

  // biaxial variants for kap = 0, kapp = 1, mu > 0
  Math::real TriaxialLine::hfun::fthtbiax(real /* tht */, real /* eps */,
                                          real /* mu */) {
    // Multiply by f functions by sqrt(abs(mu))
    // return 1 / sqrt(mu);
    return 1;
  }
  Math::real TriaxialLine::hfun::gthtbiax(real /* tht */, real /* eps */,
                                          real /* mu */) {
    return 0;
  }
  Math::real TriaxialLine::hfun::gfthtbiax(real /* tht */, real /* mu */) {
    return 0;
  }

  // biaxial variants for kap = 1, kapp = 0, mu <= 0, !_tx
  Math::real TriaxialLine::hfun::dfpsibiax(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // f functions are multiplied by sqrt(abs(mu)) but don't include this
    // factor here; instead include it in operator()(). etc.  This was we can
    // still use this function in the limit mu -> 0 to determine the conjugate
    // point on a meridian.
    return eps / (1 + sqrt(1 - eps * c2));
  }
  Math::real TriaxialLine::hfun::gpsibiax(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // return gpsip(s, c, 1, 0, eps, mu) but with the factor c2 canceled
    return sqrt(1 - eps * c2);
  }
  Math::real TriaxialLine::hfun::gfpsibiax(real s, real c, real mu) {
    // cos(phi)^2 = cos(psi)^2 - mu *sin(psi)^2
    // f' = sqrt(1-eps*cos(phi)^2)/cos(phi)^2
    // g' = sqrt(1-eps*cos(phi)^2)
    // g'/f' = cos(phi)^2
    // Adjust by sqrt(-mu) to accommodate this factor in dfpsibiax
    return gfpsip(s, c, 1, mu) / sqrt(-mu);
  }

#if 0
  // biaxial variants for kap = 1, kapp = 0, mu <= 0, _tx
  Math::real TriaxialLine::hfun::dfvbiax(real dn, real eps, real mu) {
    real c2 = Math::sq(dn);
    // Multiply by f functions by sqrt(abs(mu))
    return sqrt(-mu) * eps * dn / (1 + sqrt(1 - eps * c2));
  }
  // This is positive
  Math::real TriaxialLine::hfun::gvbiax(real /* cn */, real dn,
                                        real eps, real /* mu */) {
    // return gvp(dn, cn, 1, 0, epd, mu) but with cancelation
    real c2 = Math::sq(dn);
    return dn * sqrt(1 - eps * c2);
  }
  Math::real TriaxialLine::hfun::gfvbiax(real dn, real mu) {
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

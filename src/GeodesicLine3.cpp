/**
 * \file GeodesicLine3.cpp
 * \brief Implementation for GeographicLib::Triaxial::GeodesicLine3 class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <iostream>
#include <iomanip>
// Needed by zsetsdiag
#include <sstream>
#include <GeographicLib/Triaxial/GeodesicLine3.hpp>

namespace GeographicLib {
  namespace Triaxial {

  using namespace std;

  GeodesicLine3::GeodesicLine3(fline f, fline::fics fic,
                               gline g, gline::gics gic)
    : _tg(f.tg())
    , _f(f)
    , _fic(fic)
    , _g(g)
    , _gic(gic)
  {}

  GeodesicLine3::GeodesicLine3(const Geodesic3& tg)
    : _f(tg, Geodesic3::gamblk(tg, (tg._umbalt &&
                                    tg.kp2() > 0) || tg.k2() == 0))
  {}

  GeodesicLine3::GeodesicLine3(const Geodesic3& tg,
                               Angle bet1, Angle omg1, Angle alp1)
    : _tg(tg)
  {
    bet1.round();
    omg1.round();
    alp1.round();
    Geodesic3::gamblk gam = _tg.gamma(bet1, omg1, alp1);
    _f = fline(tg, gam);
    _fic = fline::fics(_f, bet1, omg1, alp1);
    _g = gline(tg, gam);
    _gic = gline::gics(_g, _fic);
  }

  GeodesicLine3::GeodesicLine3(const Geodesic3& tg,
                               real bet1, real omg1, real alp1)
    : GeodesicLine3(tg, ang(bet1), ang(omg1), ang(alp1))
  {}

  void GeodesicLine3::pos1(Angle& bet1, Angle& omg1, Angle& alp1) const {
    _fic.pos1(_f.transpolar(), bet1, omg1, alp1);
  }

  void GeodesicLine3::pos1(real& bet1, real& omg1, real& alp1,
                           bool unroll) const {
    ang bet1a, omg1a, alp1a;
    pos1(bet1a, omg1a, alp1a);
    if (!unroll) {
      (void) Ellipsoid3::AngNorm(bet1a, omg1a, alp1a);
      bet1a.setn(); omg1a.setn(); alp1a.setn();
    }
    bet1 = real(bet1a);
    omg1 = real(omg1a);
    alp1 = real(alp1a);
  }

  void GeodesicLine3::Position(real s12,
                               Angle& bet2a, Angle& omg2a, Angle& alp2a) const {
    // Compute points at distance s12
    real sig2 = _gic.sig1 + s12/_tg.b();
    ang &phi2a = bet2a, &tht2a = omg2a;
    int *countn = nullptr, *countb = nullptr;
    if (_f.gammax() > 0) {
      real u2, v2;
      if constexpr (false && biaxspecial(_tg, _g.gammax())) {
        u2 = gpsi().inv(sig2, countn, countb);
        v2 = ftht().inv(fpsi()(u2) - _fic.delta, countn, countb);
      } else
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
    } else if (_f.gammax() == 0) {
      pair<real, real> sig2n = remx(sig2, 2*_g.s0);  // reduce to [-s0, s0)
      real u2, v2;
      ang psi2;
      if (_f.kxp2() == 0) {
        // meridr
        //  ftht = tht
        //  gtht = 0
        // meridr
        //  fpsi = 0
        //  gpsi = int(sqrt(1-eps*cos(psi)^2), psi)
        solve2(_fic.delta, sig2n.first,
               fpsi(), ftht(), gpsi(), gtht(), u2, v2,
               countn, countb);
        phi2a = ang::radians(u2);
        // u2 in in [-pi/2, pi/2].  With long doubles cos(pi/2) < 0 which leads
        // to a switch in sheets in ellipsoidal coordinates.  Fix by setting
        // cos(phi2a) = +eps
        if (signbit(phi2a.c())) // Not triggered with doubles and quads
          phi2a = ang(copysign(real(1), phi2a.s()),
                      numeric_limits<real>::epsilon()/(1<<11), 0, true);
        psi2 = phi2a + ang::cardinal(2 * sig2n.second);
        int parity = fmod(sig2n.second, real(2)) != 0 ? -1 : 1;
        int Ny = _fic.Nx * parity;
        phi2a = phi2a.reflect(signbit(_fic.phi0.c() * Ny),
                              signbit(_fic.phi0.c())).rebase(_fic.phi0);
        tht2a = ang::radians(- _fic.delta) + ang::cardinal(2 * sig2n.second);
        alp2a = ang(_fic.Ex * real(0), _fic.Nx * parity, 0, true);
      } else {
        if (sig2n.first - _g.s0 >= -5 * numeric_limits<real>::epsilon()) {
          sig2n.first = -_g.s0;
          ++sig2n.second;
        }
        real deltax = bigclamp(_fic.delta + sig2n.second * _f.deltashift(), 1);
        solve2u(deltax, sig2n.first, fpsi(), ftht(), gpsi(), gtht(), u2, v2,
                countn, countb);
        // phi2 = fpsi().rev(u2); tht2 = ftht().rev(v2);
        phi2a = anglam(u2, _f.kxp());
        psi2 = phi2a + ang::cardinal(2 * sig2n.second);
        tht2a = anglam(v2, _f.kx());
        int parity = fmod(sig2n.second, real(2)) != 0 ? -1 : 1;
        // if_tg.t().kp2 == 0 then meridional oblate
        int Ny = _fic.Nx * parity;
        tht2a += ang::cardinal(2 * sig2n.second);
        tht2a = tht2a + _fic.tht0;
        phi2a = phi2a.reflect(signbit(_fic.phi0.c() * Ny),
                              signbit(_fic.phi0.c())).rebase(_fic.phi0);
        // replace cos(phi)/cos(tht) by sech(u)/sech(v)
        alp2a = ang(_fic.Ex * _f.kxp() / mcosh(v2, _f.kx()),
                    _f.kx() * Ny / mcosh(u2, _f.kxp()));
      }
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
  }

  void GeodesicLine3::Position(real s12,
                               real& bet2, real& omg2, real& alp2,
                               bool unroll) const {
    ang bet2a, omg2a, alp2a;
    Position(s12, bet2a, omg2a, alp2a);
    if (!unroll) {
      (void) Ellipsoid3::AngNorm(bet2a, omg2a, alp2a);
      bet2a.setn(); omg2a.setn(); alp2a.setn();
    }
    bet2 = real(bet2a);
    omg2 = real(omg2a);
    alp2 = real(alp2a);
  }

  void GeodesicLine3::solve2(real f0, real g0,
                             const hfun& fx, const hfun& fy,
                             const hfun& gx, const hfun& gy,
                             real& x, real& y,
                             int* countn, int* countb) {
    // Return x and y, s.t.
    //   fx(x) - fy(y) = f0
    //   gx(x) + gy(y) = g0
    //
    // fx(x) = fxs*x +/- fxm,
    // fy(y) = fys*y +/- fym,
    // gx(x) = gxs*x +/- gxm,
    // gy(y) = gys*y +/- gym;
    const bool check = false;
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
    if constexpr (check) {
      real DxA = (-qf * gys + qg * fys) / den,
        DyA = (-qf * gxs + qg * fxs) / den;
      real
        fA = (fx(x0 - Dx) - fy(y0 - DyA)) - f0,
        gA = (gx(x0 - Dx) + gy(y0 - DyA)) - g0,
        fB = (fx(x0 - DxA) - fy(y0 - Dy)) - f0,
        gB = (gx(x0 - DxA) + gy(y0 - Dy)) - g0,
        fC = (fx(x0 + Dx) - fy(y0 + DyA)) - f0,
        gC = (gx(x0 + Dx) + gy(y0 + DyA)) - g0,
        fD = (fx(x0 + DxA) - fy(y0 + Dy)) - f0,
        gD = (gx(x0 + DxA) + gy(y0 + Dy)) - g0;
      if (!( fabs(DxA) <= Dx && fabs(DyA) <= Dy ))
        throw GeographicLib::GeographicErr("Bad Dx/Dy");
      if (!( fA <= 0 && gA <= 0 )) {
        cout << scientific << "midA " << fA << " " << gA << "\n";
        cout << "DA " <<  Dx << " " << DxA << " " << Dy << " " << DyA << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints A GeodesicLine3::newt2");
      }
      if (!( fB >= 0 && gB <= 0 )) {
        cout << scientific << "midB " << fB << " " << gB << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints B GeodesicLine3::newt2");
      }
      if (!( fC >= 0 && gC >= 0 )) {
        cout << scientific << "midC " << fC << " " << gC << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints C GeodesicLine3::newt2");
      }
      if (!( fD <= 0 && gD >= 0 )) {
        cout << scientific << "midD " << fD << " " << gD << "\n";
        throw GeographicLib::GeographicErr
          ("Bad initial midpoints D GeodesicLine3::newt2");
      }
    }
    newt2(f0, g0, fx, fy, gx, gy,
          xm-mm*Dx, xp+mm*Dx, fx.HalfPeriod(),
          ym-mm*Dy, yp+mm*Dy, fy.HalfPeriod(),
          (fx.HalfPeriod() * fxs + fy.HalfPeriod() * fys) / 2,
          (gx.HalfPeriod() * gxs + gy.HalfPeriod() * gys) / 2,
          x, y, countn, countb);
  }

  void GeodesicLine3::solve2u(real d0, real s0,
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
    real pi2 = Geodesic3::BigValue(),
      sbet = gx.Max(), somg = gy.Max(), stot = sbet + somg,
      dbet = fx.Max(), domg = fy.Max(), del  = dbet - domg;
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
    } else if (fabs(d0) > 2*pi2/3 &&
               fabs((1 - 2 * signbit(d0)) * s0 - (sbet - somg)) <=
               7 * numeric_limits<real>::epsilon()) {
      if (d0 > 0) {
        u = 2*d0/3; v = -1*d0/3;
      } else {
        u = 1*d0/3; v = -2*d0/3;
      }
    } else {
      real mm = 2;
      newt2(d0, s0, fx, fy, gx, gy,
            -mm * pi2, mm * pi2, pi2,
            -mm * pi2, mm * pi2, pi2,
            pi2, pi2, u, v, countn, countb);
    }
  }

  void GeodesicLine3::newt2(real f0, real g0,
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
    static const real tol = numeric_limits<real>::epsilon();
    const bool debug = Geodesic3::debug_, check = false;
    // Relax [fg]tol to /10 instead of /100.  Otherwise solution resorts to
    // too much bisection.  Example:
    //   echo 63 -173 -61.97921997838416712 -4.64409746197940890408 |
    //     ./Geod3Solve $PROX
    const real
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
    auto p = zsetsbisect(xset, yset, f0, g0, false);
    x = p.first; y = p.second;
    if constexpr (check) {
      // A necessary condition for a root is
      //   f01 <= 0 <= f10
      //   g00 <= 0 <= g11
      real
        f01 = (xset.min().fz - yset.max().fz) - f0,
        f10 = (xset.max().fz - yset.min().fz) - f0,
        g00 = (xset.min().gz + yset.min().gz) - g0,
        g11 = (xset.max().gz + yset.max().gz) - g0;
      // Allow equality on the initial points
      if (!( f01 <= 0 && 0 <= f10 && g00 <= 0 && 0 <= g11 ))
        throw GeographicLib::GeographicErr
          ("Bad initial points GeodesicLine3::newt2");
      zvals xv(x, fx(x), gx(x)), yv(y, fy(y), gy(y));
    }
    // degen is a flag detected degeneracy to a 1d system.  Only the case gy(y)
    // = const is treated.  Then g = gx(x) + gy(y) - g0 = 0 converges as a 1d
    // Newton system which fixes the value of x.  Once x is found (detected by
    // xset.min().z and xset.max().z being consecutive numbers), g0 is adjusted
    // to that the g equation is satisfied exactly.  This allows f = fx(x) -
    // fy(y) - f0 = 0 to be solved reliably for y.
    bool bis = false, degen = false;
    int ibis = -1, i = 0;
    for (; i < maxit_ ||
           (Geodesic3::throw_ && (throw GeographicLib::GeographicErr
                                  ("Convergence failure GeodesicLine3::newt2"),
                                  false));
         ++i) {
      ++cntn;
      if (!degen && nextafter(xset.min().z, xset.max().z) == xset.max().z &&
          yset.min().gz == yset.max().gz) {
        degen = true;
        real ga = (xset.min().gz + yset.min().gz) - g0,
          gb = (xset.max().gz + yset.max().gz) - g0;
        x = gb < -ga ? xset.max().z : xset.min().z;
        // Adjust g0 so that g = 0 and dx = 0
        if (gb < -ga) {
          g0 = xset.max().gz + yset.max().gz;
          zvals xx = xset.max();
          xset.insert(xx, -1);
        } else {
          g0 = xset.min().gz + yset.min().gz;
          zvals xx = xset.min();
          xset.insert(xx, +1);
        }
      }
      zvals xv(x, fx(x), gx(x)), yv(y, fy(y), gy(y));
      // zsetsinsert updates xv and yv to enforce monotonicity of f and g
      zsetsinsert(xset, yset, xv, yv, f0, g0);
      real f = (xv.fz - yv.fz) - f0, g = (xv.gz + yv.gz) - g0;
      if ((fabs(f) <= ftol && fabs(g) <= gtol) || isnan(f) || isnan(g)) {
        if constexpr (debug)
          cout << "break0 " << scientific << f << " " << g << "\n";
        break;
      }
      real
        fxp = fx.deriv(x), fyp = fy.deriv(y),
        gxp = gx.deriv(x), gyp = gy.deriv(y),
        den = fxp * gyp + fyp * gxp,
        dx = -( gyp * f + fyp * g) / den, // if degen, dx = 0
        dy = -(-gxp * f + fxp * g) / den, // if degen dy = f / fyp
        xn = x + dx, yn = y + dy,
        dxa = 0, dya = 0;
      xa = xset.min().z; xb = xset.max().z;
      ya = yset.min().z; yb = yset.max().z;
      if constexpr (debug) {
        bool bb = gyp == 0 &&
          nextafter(xset.min().z, xset.max().z) == xset.max().z;
        cout << "DERIV " << i << " " << fxp << " " << fyp << " " << gxp << " " << gyp << " " << bb << "\n";
        cout << "DY " << i << " " << -gxp * f << " " << fxp * g << "\n";
        cout << "FG " << i << " " << f << " " << g << "\n";
        cout << "DXY " << i << " " << degen << " " << dx << " " << dy << "\n";
      }
      if constexpr (check) {
        if (!( fxp >= 0 && fyp >= 0 && gxp >= 0 && gyp >= 0 && den > 0 )) {
          cout << "DERIVS " << x << " " << y << " "
               << fxp << " " << fyp << " "
               << gxp << " " << gyp << " " << den << "\n";
          throw GeographicLib::GeographicErr
            ("Bad derivatives GeodesicLine3::newt2");
        }
      }
      bool cond1 = den > 0 &&
        (i < ibis + 2 ||
         ((2*fabs(f) < oldf || 2*fabs(g) < oldg) ||
          (2*fabs(dx) < olddx || 2*fabs(dy) < olddy))),
        cond2 = xn >= xa-xtol*tolmult && xn <= xb+xtol*tolmult &&
        yn >= ya-ytol*tolmult && yn <= yb+ytol*tolmult;
      if (cond1 && cond2) {
        oldf = fabs(f); oldg = fabs(g); olddx = fabs(dx); olddy = fabs(dy);
        x = xn; y = yn;
        bis = false;
        if (!(fabs(dx) > xtol || fabs(dy) > ytol)) {
          if constexpr (debug)
            cout << "break1 " << scientific << dx << " " << dy << " "
                 << f << " " << g << "\n";
          break;
        }
      } else {
        // xn = xset.bisect(); yn = yset.bisect();
        p = zsetsbisect(xset, yset, f0, g0, false);
        xn = p.first; yn = p.second;
        ++cntb;
        if (x == xn && y == yn) {
          if constexpr (debug)
            cout << "break2 " << f << " " << g << "\n";
          break;
        }
        dxa = xn - x; dya = yn - y;
        x = xn; y = yn;
        bis = true;
        ibis = i;
      }
      (void) bis;
      if constexpr (debug)
        cout << "AA " << scientific << setprecision(4)
             << x-xa << " " << xb-x << " "
             << y-ya << " " << yb-y << "\n";
      if constexpr (debug)
        cout << "CC " << i << " "
             << bis << " " << cond1 << " " << cond2 << " "
             << scientific << setprecision(2) << f << " " << g << " "
             << dx << " " << dy << " " << dxa << " " << dya << " "
             << xb-xa << " " << yb-ya << "\n";
      if constexpr (debug)
        cout << "BOX " << i << " " << bis << " "
             << xset.num() << " " << yset.num() << " "
             << scientific << setprecision(3)
             << xset.max().z - xset.min().z << " "
             << yset.max().z - yset.min().z << "\n";
    }
    if (countn)
      *countn += cntn;
    if (countb)
      *countb += cntb;
    if constexpr (debug) {
      cout << "CNT " << cntn << " " << cntb << "\n";
      cout << "XY " << setprecision(18) << x << " " << y << "\n";
    }
  }

  int GeodesicLine3::zset::insert(zvals& t, int flag) {
    // Inset t into list.  flag = -/+ 1 indicates new min/max.
    // Return -1 if t was already present; othersize return index of newly
    // inserted value.
    // If monotonic, value of t.fz and t.gz is adjusted to ensuire monotonicity:
    //   if a.z == b.z then a.fz == b.gz && a.fz == b.gz
    //   if a.z < b.z then a.fz <= b.gz && a.fz <= b.gz
    using std::isnan;
    // Need this flag to be set otherwise conditions for setting the bounds on
    // the allowed results are not met and newt2 fails to converge.
    const bool monotonic = true;
    int ind = -1;
    if (isnan(t.z)) return ind;
    if (t < min()) {
      if constexpr (monotonic) {
        t.fz = fmin(t.fz, min().fz);
        t.gz = fmin(t.gz, min().gz);
      }
    } else if (t == min()) {
      // Check if t is "other" endpoint and collapse bracket to zero
      if (flag > 0) _s.resize(1);
      t = min();
    } else if (t == max()) {
      // Check if t is "other" endpoint and collapse bracket to zero
      if (flag < 0) { _s[0] = _s.back(); _s.resize(1); }
      t = max();
    } else if (max() < t) {
      if constexpr (monotonic) {
        t.fz = fmax(t.fz, max().fz);
        t.gz = fmax(t.gz, max().gz);
      }
    }
    if (!(min() < t && t < max())) // Not in range
      return ind;
    // Now min() < t < max()
    auto p = std::lower_bound(_s.begin(), _s.end(), t);
    if (p == _s.end()) return ind; // Can't happen
    // Fix components of t
    bool ins = !(*p == t);
    if constexpr (monotonic) {
      if (ins) {
        t.fz = Math::clamp(t.fz, (p-1)->fz, p->fz);
        t.gz = Math::clamp(t.gz, (p-1)->gz, p->gz);
      } else                      // z components match
        t = *p;                   // set fz and gz values
    }
    if (flag < 0) {
      _s.erase(_s.begin(), p);
      if (ins) {
        _s.insert(_s.begin(), t);
        ind = 0;
      }
    } else if (flag > 0) {
      if (ins) {
        _s.erase(p, _s.end());
        _s.push_back(t);
        ind = int(_s.size()) - 1;
      } else
        _s.erase(p+1, _s.end());
    } else if (ins) {
      ind = int(p - _s.begin());
      _s.insert(p, t);
    }
    // else it's a duplicate and not a new end value
    return ind;
  }

  void GeodesicLine3::zsetsinsert(zset& xset, zset& yset,
                                  zvals& xfg, zvals& yfg,
                                  real f0, real g0) {
    const bool debug = Geodesic3::debug_;
    real x0 = 0, y0 = 0;
    int xind = xset.insert(xfg), yind = yset.insert(yfg);
    if constexpr (debug) {
      cout << "BOXA " << xset.num() << " " << yset.num() << " "
           << xind << " " << yind << " "
           << xset.min().z - x0 << " " << xset.max().z - x0 << " "
           << yset.min().z - y0 << " " << yset.max().z - y0 << " "
           << xset.max().z - xset.min().z << " "
           << yset.max().z - yset.min().z << "\n";
      zsetsdiag(xset, yset, f0, g0);
    }
    if (xind < 0 && yind < 0) return;
    // We update two sets of bounds xa[0], xb[0], etc. use the values at the
    // newly inserted row and column.  xa[1], xb[1], etc. use just the value at
    // the center point.
    zvals xa[2] = {xset.min(), xset.min()},
      xb[2] = {xset.max(), xset.max()},
      ya[2] = {yset.min(), yset.min()},
      yb[2] = {yset.max(), yset.max()};
    for (int i = 0; i < xset.num(); ++i) {
      const zvals& x = xset.val(i);
      for (int j = 0; j < yset.num(); ++j) {
        if (i == xind || j == yind) {
          const zvals& y = yset.val(j);
          real f = (x.fz - y.fz) - f0, g = (x.gz + y.gz) - g0;
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

          // If the pattern is
          //
          // BOXF --+
          // BOXF -++
          // BOXF -++
          //
          // BOXG ..+
          // BOXG ..+
          // BOXG ---
          //
          // The, allowing equality (tight bounds) collapses the set to 1x1
          // (the middle element).  This is clearly wrong since now we can't
          // get f = 0.  We detect this by checking the values of f and g at
          // the corners:
          //
          // f <= 0 at NW corner f >= 0 at SE corner
          // g <= 0 at SW corner g >= 0 at NE corner
          //
          // If theses requirements fail, use bounds updated with just the
          // center point i == xind && j == yind.
          if (f <= 0) {
            if (g <= 0 && xa[0] < x) xa[0] = x;
            if (g >= 0 && y < yb[0]) yb[0] = y;
            if (i == xind && j == yind) {
              if (g <= 0 && xa[1] < x) xa[1] = x;
              if (g >= 0 && y < yb[1]) yb[1] = y;
            }
          }
          if (f >= 0) {
            if (g <= 0 && ya[0] < y) ya[0] = y;
            if (g >= 0 && x < xb[0]) xb[0] = x;
            if (i == xind && j == yind) {
              if (g <= 0 && ya[1] < y) ya[1] = y;
              if (g >= 0 && x < xb[1]) xb[1] = x;
            }
          }
        }
      }
    }
    int k = 0;
    for (; k < 2; ++k) {
      // Check signs at corners
      if ((xa[k].fz - yb[k].fz) - f0 <= 0 &&
          (xb[k].fz - ya[k].fz) - f0 >= 0 &&
          (xa[k].gz + ya[k].gz) - g0 <= 0 &&
          (xb[k].gz + yb[k].gz) - g0 >= 0)
        // corner constraints met
        break;
    }
    int kk = min(k, 1);
    if constexpr (debug)
      cout << "BOXK " << k << " " << xind << " " << yind << "\n";
    xset.insert(xa[kk], -1);
    xset.insert(xb[kk], +1);
    yset.insert(ya[kk], -1);
    yset.insert(yb[kk], +1);
    if constexpr (debug) {
      cout << "BOXB " << xset.num() << " " << yset.num() << " "
           << xset.min().z - x0 << " " << xset.max().z - x0 << " "
           << yset.min().z - y0 << " " << yset.max().z - y0 << " "
           << xa[kk].z - x0 << " " << xb[kk].z - x0 << " "
           << ya[kk].z - y0 << " " << yb[kk].z - y0 << "\n";
      zsetsdiag(xset, yset, f0, g0);
    }
    if (k == 2)
      throw GeographicLib::GeographicErr
        ("Bad corner points GeodesicLine3::zsetsinsert");
  }

  void GeodesicLine3::zsetsdiag(const zset& xset, const zset& yset,
                                real f0, real g0) {
    ostringstream fs, gs;
    for (int j = yset.num() - 1; j >= 0; --j) {
      const zvals& y = yset.val(j);
      fs << "BOXF ";
      gs << "BOXG ";
      for (int i = 0; i < xset.num(); ++i) {
        const zvals& x = xset.val(i);
        real f = (x.fz - y.fz) - f0, g = (x.gz + y.gz) - g0;
        fs << (f == 0 ? '.' : f < 0 ? '-' : '+');
        gs << (g == 0 ? '.' : g < 0 ? '-' : '+');
      }
      fs << "\n";
      gs << "\n";
    }
    cout << fs.str() << gs.str();
  }

  pair<Math::real, Math::real>
  GeodesicLine3::zsetsbisect(const zset& xset, const zset& yset,
                             real f0, real g0, bool secant) {
    if constexpr (true)
      return pair<real, real>(xset.bisect(), yset.bisect());
    else if (secant && xset.num() <= 2 && yset.num() <= 2) {
      // Use secant solution
      real
        dx = xset.max().z - xset.min().z,
        dy = yset.max().z - yset.min().z,
        fx1 = (xset.max().fz - xset.min().fz) / (dx == 0 ? 1 : dx),
        fy1 = (yset.max().fz - yset.min().fz) / (dy == 0 ? 1 : dy),
        gx1 = (xset.max().gz - xset.min().gz) / (dx == 0 ? 1 : dx),
        gy1 = (yset.max().gz - yset.min().gz) / (dy == 0 ? 1 : dy),
        den = fx1*gy1 + fy1*gx1,
        f00 = (xset.min().fz + yset.min().fz) - f0,
        g00 = (xset.min().gz + yset.min().gz) - g0,
        x = -gy1*f00 - fy1*g00,
        y =  gx1*f00 - fx1*g00;
      if (den <= 0) den = 1;
      x /= den;
      y /= den;
      if (x <= 0 || x >= 1) x = 1/real(2);
      if (y <= 0 || y >= 1) y = 1/real(2);
      x = xset.min().z + x * dx;
      y = yset.min().z + y * dy;
      return pair<real, real>(x, y);
    } else {
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
            f01 = (xset.val(i ).fz - yset.val(j1).fz) - f0,
            f10 = (xset.val(i1).fz - yset.val(j ).fz) - f0,
            g00 = (xset.val(i ).gz + yset.val(j ).gz) - g0,
            g11 = (xset.val(i1).gz + yset.val(j1).gz) - g0;
          if (f01 <= 0 && 0 <= f10 && g00 <= 0 && 0 <= g11) {
            ++cnt;
            if (xgap1 > xgap) { xgap = xgap1; x = xmean; }
            if (ygap1 > ygap) { ygap = ygap1; y = ymean; }
          }
        }
      }
      if (cnt == 0)
        throw GeographicLib::GeographicErr
          ("No legal box GeodesicLine3::zsetsbisect");
      return pair<real, real>(x, y);
    }
  }

  void GeodesicLine3::Hybrid(ang betomg2,
                             ang& bet2a, ang& omg2a, ang& alp2a,
                             real& s12, bool betp) const {
    fline::disttx d = _f.Hybrid(_fic, betomg2, bet2a, omg2a, alp2a, betp);
    s12 = _g.dist(_gic, d);
  }

  GeodesicLine3::fline::disttx
  GeodesicLine3::fline::Hybrid(const fics& fic, ang betomg2,
                               ang& bet2a, ang& omg2a, ang& alp2a,
                               bool betp) const {
    // Is the control variable psi or tht?
    bool psip = !transpolar() ? betp : !betp;
    if (!betp) betomg2 -= ang::cardinal(1);
    ang tau12;
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
        // FIX THIS!! betp shouldn't be appearing here.
        // Test case
        // echo -88 21 88 -111 | ./Geod3Solve -i $SET --hybridalt
        ang tht2b = tht2; tht2b.reflect(false, betp && fic.Ex < 0);
        ang tht12 = tht2b - fic.tht1;
        // convert -180deg to 180deg
        if (signbit(tht12.s()))
          tht12 = ang(0, copysign(real(1), tht12.c()), 0, true);
        tau12 = tht12;
      } else
        tau12 = ang::NaN();
    }
    disttx ret = ArcPos0(fic, tau12.base(), bet2a, omg2a, alp2a, betp);
    return ret;
  }

  GeodesicLine3::fline::fline(const Geodesic3& tg, bool neg)
    : _tg(tg)
    , _gm(tg, neg)
  {}

  GeodesicLine3::fline::fline(const Geodesic3& tg, Geodesic3::gamblk gam)
    : _tg(tg)
    , _gm(gam)
    , _fpsi(false, _gm.kx2 , _gm.kxp2,
            +(_gm.transpolar ? -1 : 1) * _tg.e2(),
            -(_gm.transpolar ? -1 : 1) * _gm.gamma, _tg)
    , _ftht(false, _gm.kxp2 , _gm.kx2,
            -(_gm.transpolar ? -1 : 1) * _tg.e2(),
            +(_gm.transpolar ? -1 : 1) * _gm.gamma, _tg)
  {
    // Only needed for umbilical lines
    _deltashift = _gm.gamma == 0 ?
      (_tg.k2() > 0 && _tg.kp2() > 0 ? 2 * (_fpsi.Max() - _ftht.Max()) : 0) :
      Math::NaN();
  }

  GeodesicLine3::gline::gline(const Geodesic3& tg, bool neg)
    : _tg(tg)
    , _gm(tg, neg)
  {}

  GeodesicLine3::gline::gline(const Geodesic3& tg, const Geodesic3::gamblk& gm)
    : _tg(tg)
    , _gm(gm)
    , _gpsi(true, _gm.kx2 , _gm.kxp2,
            +(_gm.transpolar ? -1 : 1) * _tg.e2(),
            -(_gm.transpolar ? -1 : 1) * _gm.gamma, _tg)
    , _gtht(true, _gm.kxp2 , _gm.kx2,
            -(_gm.transpolar ? -1 : 1) * _tg.e2(),
            +(_gm.transpolar ? -1 : 1) * _gm.gamma, _tg)
    , s0(_gm.gammax == 0 ? _gpsi.Max() + _gtht.Max() : 0)
  {}

  Math::real GeodesicLine3::fline::Hybrid0(const fics& fic, ang bet2, ang omg2,
                                           bool betp) const {
    ang bet2a, omg2a, alp2a;
    (void) Hybrid(fic, betp ? bet2 : omg2, bet2a, omg2a, alp2a, betp);
    bool angnorm = true || betp;
    if (angnorm)
      (void) Ellipsoid3::AngNorm(bet2a, omg2a, alp2a, !betp);
    if (betp) {
      omg2a -= omg2;
      return omg2a.radians0();
    } else {
      bet2a -= bet2;
      return angnorm ? bet2a.radians0() : bet2.radians();
    }
  }

  GeodesicLine3::fline::disttx
  GeodesicLine3::fline::ArcPos0(const fics& fic, ang tau12,
                                ang& bet2a, ang& omg2a, ang& alp2a,
                                bool betp) const {
    // XXX fix for biaxial
    disttx ret{Math::NaN(), Math::NaN(), 0};
    bool psip = transpolar() ? !betp : betp;
    ang &phi2a = bet2a, &tht2a = omg2a;
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
      }
      // Already normalized
      phi2a = ang(nup() * psi2.s(),
                  fic.phi0.c() * hypot(psi2.c(), nu() * psi2.s()),
                  0, true).rebase(fic.phi0);
      real s = fic.Ex * hypot(kx() * nu(), kxp() * tht2a.c()),
        c = fic.phi0.c() * kx() * nup() * psi2.c();
      if (s == 0 && c == 0)
        (transpolar() ? s : c) = 1;
      alp2a = ang(s, c);
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
        u2 = fpsi().fwd(phi2n.first);
        int parity = fmod(phi2n.second, real(2)) != 0 ? -1 : 1;
        if (kxp() == 0) {
          // v2 is independent on u2
          v2 = 0;
          tht2a = ang::radians(-fic.Ex * fic.delta) +
            ang::cardinal(2*fic.Ex*phi2n.second);
          alp2a = fic.alp1.nearest(2U) + ang::cardinal(parity < 0 ? 2 : 0);
        } else {
          real deltax = bigclamp(fic.delta + phi2n.second * _deltashift, 2);
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
            // already handled in Geodesic3::Inverse.
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
              l = exp(Geodesic3::BigValue()),
              tpsi2 = Trigfun::root(Trigfun::ARCPOS0,
                                    [this, npi] (real tpsi) -> pair<real, real>
                                    {
                                      real psi = atan(tpsi);
                                      return pair<real, real>
                                        (tpsi - fpsi().df(npi + psi),
                                         1 - fpsi().dfp(psi) /
                                         (1 + Math::sq(tpsi)));
                                    },
                                    c, fic.psi1.t(), -l, l);
            v2 = atan(tpsi2);
            phi2a = (ang(tpsi2, 1) + fic.phi1.nearest(2U))
              .flipsign(parity*fic.Nx).rebase(fic.phi0);
            alp2a = ang::cardinal(fic.Nx * parity == 1 ? 0 : 2);
          } else {
            u2 = v2 == 0 ? 0 : copysign(Math::pi()/2, tht2n.first);
            phi2a = ang::cardinal(fabs(v2) == 0
                                  ? 0 : copysign(real(1), v2 * fic.Nx))
              .rebase(fic.phi0);
            alp2a = fabs(v2) == 0 ? ang::cardinal(2) :
              fic.alp1.nearest(2U) +
              ang::cardinal(parity == 1 ? 0 : 2) +
              ang::radians(v2).flipsign(parity * phi2a.s());
          }
        } else {
          real deltax = bigclamp(fic.delta + tht2n.second * _deltashift, 2);
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
    return ret;
  }

  GeodesicLine3::fline::fics::fics(const fline& f, ang bet10, ang omg10,
                                   ang alp10)
    : tht1(omg10 - ang::cardinal(1))
    , phi1(bet10)
    , alp1(alp10)
  {
    static const real eps = numeric_limits<real>::epsilon();
    alp0 = alp1.nearest(f.transpolar() ? 2U : 1U);
    if (f.transpolar()) {
      swap(tht1, phi1);
      alp1.reflect(false, false, true);
    }
    const Geodesic3& tg = f.tg();
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
      u0 = f.fpsi().fwd(psi1.radians());
      v0 = f.ftht().fwd(tht1.radians());
      delta = (biaxspecial(tg, f.gammax()) ?
               atan2(phi1.s() * fabs(alp1.s()), phi0.c() * alp1.c())
               - sqrt(f.gammax()) * f.fpsi().df(u0)
               : f.fpsi()(u0)) - f.ftht()(v0);
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
        u0 = f.fpsi().fwd((phi1-phi0).radians());
        v0 = 0;
        delta = -f.ftht()(tht0.radians());
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
          delta = Nx * f.fpsi()(lamang(phi1 - phi0, tg.kp())) -
            f.ftht()(lamang(tht1 - tht0, tg.k()));
        }
      }
    } else {
      // gamma = NaN
    }
  }

  void GeodesicLine3::fline::fics::pos1(bool transpolar,
                                        ang& bet10, ang& omg10, ang& alp10)
    const {
    bet10 = phi1; omg10 = tht1.flipsign(Ex);
    alp10 = alp1;
    if (transpolar) {
      swap(bet10, omg10);
      alp10.reflect(false, false, true);
    }
    omg10 += ang::cardinal(1);
  }

  void GeodesicLine3::fline::fics::setquadrant(const fline& f, unsigned q) {
    ang bet1, omg1, alp1x;  // TODO: Fix so fics uses tau1
    pos1(f.transpolar(), bet1, omg1, alp1x);
    alp1x.setquadrant(q);
    *this = fics(f, bet1, omg1, alp1x);
  }

  GeodesicLine3::gline::gics::gics(const gline& g, const fline::fics& fic) {
    if (g.gammax() > 0) {
      sig1 = g.gpsi()(fic.u0) + g.gtht()(fic.v0);
    } else if (g.gammax() == 0) {
      sig1 = g.kxp2() == 0 ? fic.Nx * g.gpsi()(fic.u0) :
        fic.Nx * g.gpsi()(lamang(fic.phi1 - fic.phi0, g.kxp())) +
        g.gtht()(lamang(fic.tht1 - fic.tht0, g.kx()));
    } else {
      // gamma = NaN
    }
  }

  Math::real GeodesicLine3::gline::dist(gics ic, fline::disttx d) const {
    real sig2 = gpsi()(d.phiw2) + gtht()(d.thtw2) + d.ind2 * 2*s0;
    return (sig2 - ic.sig1) * tg().b();
  }

  GeodesicLine3::hfun::hfun(bool distp, real kap, real kapp, real eps, real mu,
                            const Geodesic3& tg)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _sqrtmu(sqrt(fabs(_mu)))
    , _sqrtkap(sqrt(_kap))
    , _sqrtkapp(sqrt(_kapp))
    , _distp(distp)
    , _umb(!tg.biaxial() && _mu == 0)
    , _meridr(_kap == 0 && _mu == 0)
    , _meridl(_kapp == 0 && _mu == 0)
    , _biaxr(biaxspecial(tg, _mu) && _kap  == 0)
    , _biaxl(biaxspecial(tg, _mu) && _kapp == 0)
  {
    // mu in [-kap, kapp], eps in (-inf, 1/kap)
    if (!_distp) {
      if (_meridr || _biaxr) {
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
      } else if (_meridl || _biaxl) {
        // biaxial librating coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < tg._ellipthresh;
        _tx = false;
        // f explicitly multiplied by sqrt(-mu) in operator()() for _biaxl
        // For _meridl operator() ignores this
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return dfpsibiax(sin(psi), cos(psi), eps, mu); },
                          Math::pi()/2, false);
      } else if (_mu > 0) {
        _tx = _mu / (_kap + _mu) < tg._ellipthresh;
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
        _tx = -_mu / _kap < tg._ellipthresh;
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
        _tx = _kapp < tg._ellipthresh;
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
                            [
#if defined(_MSC_VER)
                             // Visual Studio requires the capture of this.
                             // I can't see why -- probably a bug?
                             this,
#endif
                             kap = _kap, kapp = _kapp, eps = _eps]
                            (real tht) -> real
                            { return dfp(cos(tht), kap, kapp, eps); },
                            Math::pi(), true, 1);
      } else {
        // _mu == NaN
        _tx = false;
      }
    } else {                    // _distp
      if (_meridr || _biaxr) {
        // biaxial symmetry coordinate
        // _kapp == 1, mu < 0 not allowed
        _tx = false;
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real tht) -> real
                          // degenerate f' = 0
                          { return gthtbiax(tht, eps, mu); },
                          Math::pi()/2, false);
      } else if (_meridl || _biaxl) {
        // biaxial non-symmetry coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < tg._ellipthresh;
        _tx = false;
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return gpsibiax(sin(psi), cos(psi), eps, mu); },
                          Math::pi()/2, false);
      } else if (_mu > 0) {
        _tx = _mu / (_kap + _mu) < tg._ellipthresh;
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
        _tx = -_mu / _kap < tg._ellipthresh;
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
        _tx = _kapp < tg._ellipthresh;
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
    _max = _umb ? _fun(_tx ? _ell.K() : Math::pi()/2) :
      !_distp ? ( _biaxl ? Math::pi()/2 + _sqrtmu * _fun.Max() :
                  _fun.Max() ) :
      _meridl ? _fun(Math::pi()/2) : _fun.Max();
  }

  Math::real GeodesicLine3::hfun::operator()(real u) const {
    if (!_distp) {
      if (_biaxl)
        // This is sqrt(-mu) * f(u)
        return modang(u, _sqrtmu) - _sqrtmu * _fun(u);
      else if (_meridl)
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

  Math::real GeodesicLine3::hfun::deriv(real u) const {
    if (!_distp) {
      if (_biaxl)
        // This is sqrt(-mu) * f'(u)
        return _sqrtmu / (Math::sq(cos(u)) - _mu * Math::sq(sin(u)))
          - _sqrtmu * _fun.deriv(u);
      else if (_meridl)
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

  Math::real GeodesicLine3::hfun::root(real z, real u0,
                                       int* countn, int* countb, real tol)
    const {
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
        return Trigfun::root(Trigfun::FFUNROOT,
                             [this]
                             (real u) -> pair<real, real>
                             { return pair<real, real>((*this)(u), deriv(u)); },
                             z,
                             u0, ua, ub,
                             HalfPeriod(), HalfPeriod()/Slope(), 1,
                             countn, countb, tol);
      } else if (_umb) {
        real d = fabs(Max())
          + 2 * numeric_limits<real>::epsilon() * fmax(real(1), fabs(z)),
          ua = z - d,
          ub = z + d;
        u0 = fmin(ub, fmax(ua, u0));
        return Trigfun::root(Trigfun::FFUNROOT,
                             [this]
                             (real u) -> pair<real, real>
                             { return pair<real, real>((*this)(u), deriv(u)); },
                             z,
                             u0, ua, ub,
                             Math::pi()/2, Math::pi()/2, 1,
                             countn, countb, tol);
      } else
        return Math::NaN();
    } else {
      // This function isn't neeed.  General inversion mechanisms in Trigfun
      // suffice.  NO, the trigfun for _umb is not invertible.
      if (!(isfinite(z) && _umb))
        return Math::NaN();       // Deals with +/-inf and nan
      // Now we're dealing with _umb.
      if (fabs(z) >= Max())
        return copysign(Geodesic3::BigValue(), z);
      real ua = -Geodesic3::BigValue(), ub = -ua;
      u0 = fmin(ub, fmax(ua, u0));
      // Solve z = _fun(_tx ? _ell.F(gd(u)) : gd(u)) for u
      return Trigfun::root(Trigfun::GFUNROOT,
                           [this]
                           (real u) -> pair<real, real>
                           { return pair<real, real>((*this)(u), deriv(u)); },
                           z,
                           u0, ua, ub,
                           Math::pi()/2, Math::pi()/2, 1, countn, countb, tol);
    }
  }

  // Accurate inverse by direct Newton (not using _finv)
  Math::real GeodesicLine3::hfun::inv(real z, int* countn, int* countb)
    const {
    if (!_distp)
      return _umb ? root(z, z, countn, countb) :
        _biaxl ? root(z, modang(z/Slope(), 1/_sqrtmu), countn, countb) :
        _fun.inv1(z, countn, countb);
    else {                      // distp
      if (_biaxr) return Math::NaN();
      else if (_umb) {
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
        real u0 = atanh(_sqrtkapp/_sqrtkap *
                          tan(atan(_sqrtkap/_sqrtkapp) / _max * z)) /
          (sqrt(1 - _eps*_kap) * atan(_sqrtkap/_sqrtkapp) / _max);
        return root(z, u0, countn, countb);
      } else
        return _fun.inv1(z, countn, countb);
    }
  }

  // _mu > 0 && !_tx
  Math::real GeodesicLine3::hfun::fthtp(real c, real kap, real kapp,
                                        real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }
  // This is nonnegative
  Math::real GeodesicLine3::hfun::gthtp(real c, real kap, real kapp,
                                        real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return c2 * sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }

  // _mu > 0 && _tx
  Math::real GeodesicLine3::hfun::fup(real cn, real kap, real kapp,
                                      real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }
  // This is nonnegative
  Math::real GeodesicLine3::hfun::gup(real cn, real /* dn */,
                                      real kap, real kapp, real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return c2 * sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }

  // _mu == 0 && !_tx
  Math::real GeodesicLine3::hfun::dfp(real c, real kap, real kapp, real eps) {
    // function dfp = dfpf(phi, kappa, epsilon)
    // return derivative of sqrt(kap * kapp) * Delta f
    // s = sqrt(1 - kap * sin(phi)^2)
    real c2 = kap * Math::sq(c), s = sqrt(kapp + c2);
    return eps*kap * sqrt(kapp) * c / (s * (1 + sqrt(1 - eps*c2)));
  }
  Math::real GeodesicLine3::hfun::g0p(real c, real kap, real kapp, real eps) {
    real c2 = kap * Math::sq(c);
    return sqrt( kap * (1 - eps * c2) / (kapp + c2) ) * c;
  }

  // _mu == 0 && _tx
  Math::real GeodesicLine3::hfun::dfvp(real cn, real /* dn */,
                                       real kap, real kapp, real eps) {
    // function dfvp = dfvpf(v, kap, eps)
    // return derivative of sqrt(kap * kapp) * Delta f_v
    return eps*kap * sqrt(kapp) * cn /
      (1  + sqrt(1 - eps*kap * Math::sq(cn)));
  }
  Math::real GeodesicLine3::hfun::g0vp(real cn, real kap, real /* kapp */,
                                       real eps) {
    real c2 = kap * Math::sq(cn);
    return sqrt( kap * (1 - eps * c2) ) * cn;
  }

  // _mu < 0 && !_tx
  Math::real GeodesicLine3::hfun::fpsip(real s, real c, real kap, real kapp,
                                        real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * c2) ) ;
  }
  // This is positive
  Math::real GeodesicLine3::hfun::gpsip(real s, real c, real kap, real kapp,
                                        real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt(c2 * (1 - eps * c2) / (kapp + c2)) ;
  }

  // _mu < 0 && _tx
  Math::real GeodesicLine3::hfun::fvp(real dn, real kap, real kapp,
                                      real eps, real /* mu */) {
    real c2 = kap * Math::sq(dn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * kap) );
  }
  // This is positive
  Math::real GeodesicLine3::hfun::gvp(real /* cn */, real dn,
                                      real kap, real kapp,
                                      real eps, real /* mu */) {
    real dn2 = Math::sq(dn), c2 = kap * dn2;
    return dn2 * sqrt( kap * (1 - eps * c2) / (kapp + c2) );
  }

  // biaxial variants for kap = 0, kapp = 1, mu > 0
  Math::real
  GeodesicLine3::hfun::fthtbiax(real /* tht */, real /* eps */,
                                real /* mu */) {
    // Multiply by f functions by sqrt(abs(mu))
    // return 1 / sqrt(mu);
    return 1;
  }
  Math::real
  GeodesicLine3::hfun::gthtbiax(real /* tht */, real /* eps */,
                                real /* mu */) {
    return 0;
  }

  // biaxial variants for kap = 1, kapp = 0, mu <= 0, !_tx
  Math::real
  GeodesicLine3::hfun::dfpsibiax(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // f functions are multiplied by sqrt(abs(mu)) but don't include this
    // factor here; instead include it in operator()(). etc.  This was we can
    // still use this function in the limit mu -> 0 to determine the conjugate
    // point on a meridian.
    return eps / (1 + sqrt(1 - eps * c2));
  }
  Math::real
  GeodesicLine3::hfun::gpsibiax(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // return gpsip(s, c, 1, 0, eps, mu) but with the factor c2 canceled
    return sqrt(1 - eps * c2);
  }

  } // namespace Triaxial
} // namespace GeographicLib

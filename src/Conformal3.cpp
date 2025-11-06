/**
 * \file Conformal3.cpp
 * \brief Implementation for GeographicLib::Triaxial::Conformal3 class
 *
 * Copyright (c) Charles Karney (2014-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <iostream>
#include <iomanip>
#include <GeographicLib/Triaxial/Conformal3.hpp>
#include <GeographicLib/Trigfun.hpp>

namespace GeographicLib {
  namespace Triaxial {

  using namespace std;

  Conformal3::Conformal3(const Ellipsoid3& t)
    : _t(t)
    , _ex(kp2() * (1 - e2() * k2 ()),   - e2() * kp2(),
          k2 () * (1 + e2() * kp2()), 1 + e2() * kp2())
    , _ey(k2 () * (1 + e2() * kp2()),     e2() * k2 (),
          kp2() * (1 - e2() * k2 ()), 1 - e2() * k2 ())
    , _x( Math::sq(a()) / b() * _ex.Pi() )
    , _y( Math::sq(c()) / b() * _ey.Pi() )
  {
    _s = EquivSphere(_x, _y, _exs, _eys);
  }
  // Proving the equivalence of this and the Nyrtsov formulations is
  // straightforward but laborious.  Working from this formulation:
  // * switch the angle argument of Pi to its complement
  // * sin(phi)^2 -> cos(phi')^2 = 1 - sin(phi')^2
  // * modulus of Pi then becomes imaginary
  // * Use DLMF 19.7.5 to express as sum of Pi term and F term

  Conformal3::Conformal3(real a, real b, real c)
    : Conformal3(Ellipsoid3(a, b, c))
  {}
  Conformal3::Conformal3(real b, real e2, real k2, real kp2)
    : Conformal3(Ellipsoid3(b, e2, k2, kp2))
  {}

  Math::real Conformal3::Pi(const EllipticFunction& ell, ang phi) {
    if (ell.kp2() == 0 && (signbit(phi.c()) || phi.n() != 0))
      return Math::NaN();
    real p = ell.Pi(phi.s(), phi.c(), ell.Delta(phi.s(), phi.c()));
    if (phi.n() != 0) p += 4 * ell.Pi() * phi.n();
    return p;
  }
  Angle Conformal3::Piinv(const EllipticFunction& ell, real x) {
    real y, n;
    if (ell.kp2() == 0) {
      // ell.Pi() == inf
      y = x; n = 0;
    } else {
      y = remainder(x, 2 * ell.Pi());
      n = 2 * round((x - y) / (2 * ell.Pi()));
    }
    // Now x = n * Pi() + y where y in [-Pi(), Pi()].  Pi() is the quarter
    // period for the elliptic integral which corresponds to pi/2 in angle
    // space.
    if (y == 0)
      return ang::cardinal( n == 0 ? y : n ); // Preserve the sign of +/-0
    else if (fabs(y) == ell.Pi())
      return ang::cardinal(copysign(real(1), y) + n);
    else {
      // solve Pi(phi) = y for phi
      // Pi'(phi) = 1/(sqrt(1 - ell.k2() * Math::sq(sin(phi)))
      //               * (1 - alpha2 * Math::sq(sin(phi))))
      // For k2 in [0,1]
      //    1 - ell.k2() * Math::sq(sin(phi))
      //    = ell.kp2() + ell.k2() * Math::sq(cos(phi))
      // For alpha2 > 0
      //    1 - alpha2 * Math::sq(sin(phi))
      //    alphap2 + alpha2 * Math::sq(sin(phi))
      int countn = 0, countb = 0;
      auto Pif = [&ell]
        (real phi) -> pair<real, real>
        {
          real s = sin(phi), c = cos(phi),
          f = ell.Pi(s, c, ell.Delta(s, c)),
          fp = 1 / (sqrt(ell.kp2() + ell.k2() * c*c) *
                    (ell.alpha2() >= 0 ?
                     ell.alphap2() + ell.alpha2() * s*s :
                     1 - ell.alpha2() * c*c));
          return pair<real, real>(f, fp);
        };
      real z = Trigfun::root(Trigfun::PIINV,
                             Pif, fabs(y), fabs(y) * Math::pi()/(2*ell.Pi()),
                             0, Math::pi()/2,
                             1,1,1,
                             &countn, &countb);
      (void) countn; (void) countb;
      // cout << "CNT " << countn << " " << countb << "\n";
      return ang::radians(copysign(z, y)) + ang::cardinal(n);
    }
  }
  Math::real Conformal3::F(const EllipticFunction& ell, ang phi) {
    if (ell.kp2() == 0 && (signbit(phi.c()) || phi.n() != 0))
      return Math::NaN();
    real p = ell.F(phi.s(), phi.c(), ell.Delta(phi.s(), phi.c()));
    if (phi.n() != 0) p += 4 * ell.K() * phi.n();
    return p;
  }
  Angle Conformal3::Finv(const EllipticFunction& ell, real x) {
    real y, n;
    if (ell.kp2() == 0) {
      // ell.K() == inf
      y = x; n = 0;
    } else {
      y = remainder(x, 2 * ell.K());
      n = 2 * round((x - y) / (2 * ell.K()));
    }
    // Now x = n * K() + y where y in [-K(), K()].  K() is the quarter
    // period for the elliptic integral which corresponds to pi/2 in angle
    // space.
    if (y == 0)
      return ang::cardinal( n == 0 ? y : n ); // Preserve the sign of +/-0
    else if (fabs(y) == ell.K())                // inf == inf is true
      return ang::cardinal(copysign(real(1), y) + n);
    else {
      // solve F(phi) = y for phi
      // F'(phi) = 1/sqrt(1 - ell.k2() * Math::sq(sin(phi)))
      // For k2 in [0,1]
      //    1 - ell.k2() * Math::sq(sin(phi))
      //    = ell.kp2() + ell.k2() * Math::sq(cos(phi))
      int countn = 0, countb = 0;
      auto Ff = [&ell]
        (real phi) -> pair<real, real>
        {
          real s = sin(phi), c = cos(phi),
          f = ell.F(s, c, ell.Delta(s, c)),
          fp = 1 / sqrt(ell.kp2() + ell.k2() * c*c);
          return pair<real, real>(f, fp);
        };
      real z = Trigfun::root(Trigfun::FINV,
                             Ff, fabs(y), fabs(y) * Math::pi()/(2*ell.K()),
                             0, Math::pi()/2,
                             1,1,1,
                             &countn, &countb);
      (void) countn; (void) countb;
      // cout << "CNT " << countn << " " << countb << "\n";
      return ang::radians(copysign(z, y)) + ang::cardinal(n);
    }
  }

  Math::real Conformal3::x(Angle omg) const {
    omg = omegashift(omg, +1);
    if (k2() == 0 && (signbit(omg.c()) || omg.n() != 0))
        return Math::NaN();
    return Math::sq(a()) / b() * Pi(_ex, omg.modang(b() / a()));
  }
  Angle Conformal3::omega(real x) const {
    ang omg = Piinv(_ex, (x * b()) / Math::sq(a())).modang(a() / b());
    return omegashift(omg, -1);
  }

  Math::real Conformal3::y(Angle bet) const {
    if (kp2() == 0 && (signbit(bet.c()) || bet.n() != 0))
        return Math::NaN();
    return  Math::sq(c()) / b() * Pi(_ey, bet.modang(b() / c()));
  }
  Angle Conformal3::beta(real y) const {
    return Piinv(_ey, (y * b()) / Math::sq(c())).modang(c() / b());
  }

  void Conformal3::Forward(Angle bet, Angle omg, real& xa, real& ya) const {
    xa = x(omg); ya = y(bet);
  }
  void Conformal3::Reverse(real x, real y, Angle& bet, Angle& omg) const {
    bet = beta(y); omg = omega(x);
  }
  void Conformal3::Forward(Angle bet, Angle omg, real& x, real& y,
                           real& m) const {
    Forward(bet, omg, x, y); m = 1/invscale(bet, omg);
  }
  void Conformal3::Reverse(real x, real y, Angle& bet, Angle& omg,
                           real& m) const {
    Reverse(x, y, bet, omg); m = 1/invscale(bet, omg);
  }

  Ellipsoid3 Conformal3::EquivSphere(real x, real y,
                                     EllipticFunction& ellx,
                                     EllipticFunction& elly) {
    // find b, k2, kp2 s.t.
    // b*K(kp2) = x
    // b*K(k2)  = y
    // x*K(k2) - y*K(kp2) = 0

    static const real
      N = (log(real(4)) - log(Math::pi())) / (Math::pi()/2 - log(real(4))),
      B = exp(N * Math::pi()/2) - pow(real(4), N);
    real k2 = 1/real(2);
    bool swapxy = x < y;
    int countn = 0, countb = 0;
    if (swapxy) swap(x, y);
    do {
      // Now x >= y, k2 <= 1/2
      real s =  x + y, nx = x/s, ny = y/s;
      if (nx == ny) break;      // k2 = 1/2
      if (ny == 0) { k2 = 0; break; }
      // Find initial guess assume K(k2) = pi/2, so K(kp2) = nx/ny * pi/2.
      // Invert using approximate k(K) given in
      // https://arxiv.org/abs/2505.17159v4
      real KK = nx/ny * Math::pi()/2;
      k2 = 16/pow(exp(N*KK) - B, 2/N);
      // Alternatively using KK = 1/2*log(16/kp) A+S 17.3.26
      k2 = fmin(1/real(2), 16*exp(-2*KK)); // Make sure guess is sane
      static const real logk2min = 2*log(numeric_limits<real>::epsilon());
      // Solve for log(k2) to preserve relative accuracy for tiny k2.
      real logk2 = log(k2);
      if (logk2 > logk2min) {
        auto ksolve = [nx, ny]
          (real logk2) -> pair<real, real>
          {
            real k2 = exp(logk2), kp2 = 1 - k2;
            EllipticFunction elly(k2), ellx(kp2, 0, k2, 1);
            real f = nx * elly.K() - ny * ellx.K(),
            fp = (nx * (elly.E() - kp2 * elly.K()) +
                  ny * (ellx.E() - k2  * ellx.K())) / (2 * k2 * kp2);
            return pair<real, real>(f, k2*fp);
          };
        logk2 = Trigfun::root(Trigfun::KINV, ksolve, 0, logk2,
                              logk2min, -log(real(2)),
                              1, 1, 1,
                              &countn, &countb);
      }
      // otherwise accept the asymptotic result
      k2 = exp(logk2);
    } while (false);
    // b*K(kp2) = x
    // b*K(k2)  = y

    elly = EllipticFunction(k2);
    ellx = EllipticFunction(1-k2, 0, k2, 1);
    real b = y / elly.K();
    real kp2 = 1 - k2;
    if (swapxy) {
      swap(x, y);
      swap(k2, kp2);
      swap(ellx, elly);
    }
    return Ellipsoid3(b, 0, k2, kp2);
  }
  Math::real Conformal3::sphericalscale(real ma, real mb) const {
    real m;
    if (ma == 0 || mb == 0) {
      // Let bet = pi/2 - db, omg = pi - do
      // bet' = pi/2 - c/b*db, omg' = pi/2 - a/b*do
      // Pi(pi/2 - d, a2, k2) = Pi(a2, k2) - d/(ap2*kp)
      // x = x0 - (a^2/b) * (a/b*do) / ((1+e2*kp2)*k*sqrt(1+e2*kp2))
      //   = x0 - (a^2/b) * (a/b*do) / (a^2/b^2*k*a/b)
      //   = x0 - b * do / k
      // y = y0 - (c^2/b) * (c/b*db) / ((1-e2*k2)*kp*sqrt(1-e2*k2))
      //   = y0 - (c^2/b) * (c/b*db) / (c^2/b^2*kp*c/b)
      //   = y0 - b * db / kp
      // for triaxial -> plane
      // mtp = 1/sqrt(k2*db^2 + kp2*do^2)
      // on equivalent sphere
      // x = x0 - b * do / k = x0 - bs * dos / ks
      // y = y0 - b * db / kp = y0 - bs * dbs / kps
      // implies
      // dos = b/bs * ks/k * do
      // dbs = b/bs * kps/kp * db
      // for plane -> sphere
      // mps = sqrt(k2s*dbs^2 + kp2s*dos^2)
      //     = b/bs * sqrt(k2s*kp2s/kp2 * db^2 + kp2s*k2s/k2 * do^2)
      //     = ks*kps/(k*kp) b/bs * sqrt(k2*db^2 + kp2*do^2)
      // Product of scales = ks*kps/(k*kp) * b/bs

      // Oblate case see Aux Lat paper Eqs (50) + (51)
      // sig(1) = sinh(e*atanh(e)) = sg
      // tan(chi) = tan(phi)*sqrt(1+sg^2) - sg*tan(phi)
      //          = tan(phi) * (sqrt(1+sg^2) - sg)
      //          = tan(phi) / (sqrt(1+sg^2) + sg)
      // cos(chi) = cos(phi) * (sqrt(1+sg^2) + sg)
      // Triaxial -> plane scale beta = pi/2 - db
      // mtp = 1/db
      // beta' = pi/2 - c/b*db = phi, chi = betas
      // cos(phi) = c/b * db
      // cos(chi) = c/b * (sqrt(1+sg^2) + sg) * db
      // mps = cos(chi) = c/b * (sqrt(1+sg^2) + sg) * db
      // Product of scales = c/b * (sqrt(1+sg^2) + sg)

      // Prolate case, replace sg = sinh(-e*atan(e)); c/b -> a/b
      // sg = sinh(e*atan(e))
      // Product of scales = a/b / (sqrt(1+sg^2) + sg)
      if (kp2() == 0) {
        // oblate pole
        real e = sqrt(e2()), sg = sinh(e * atanh(e));
        m = (c()/b()) * (hypot(sg, real(1)) + sg);
      } else if (k2() == 0) {
        // prolate pole
        real e = sqrt(e2()), sg = sinh(e * atan(e));
        m = (a()/b()) / (hypot(sg, real(1)) + sg);
      } else {
        // triaxial umbilical
        m =  sqrt(_s.k2() * _s.kp2()/(k2() * kp2())) * (b()/_s.b());
      }
    }
    else
      m = mb/ma;
    return m;
  }

  void Conformal3::ForwardSphere(Angle bet, Angle omg,
                                 vec3& r, vec3& v, real& m) const {
    real x, y;
    Forward(bet, omg, x, y, m);
    real ma = invscale(bet, omg);
    ang omgs = omegashift(Finv(_exs, x/_s.b()), -1),
      bets = Finv(_eys, y/_s.b()),
      alp{};                    // alp = 0, due North
    real mb = sqrt(_s.k2 () * Math::sq(bets.c()) +
                   _s.kp2() * Math::sq(omgs.s()));
    m = sphericalscale(ma, mb);
    _s.elliptocart2(bets, omgs, alp, r, v);
  }

  void Conformal3::ReverseSphere(vec3 r, vec3 v, Angle& bet, Angle& omg,
                                 Angle& gam, real& m) const {
    ang bets, omgs, alp;
    _s.cart2toellip(r, v, bets, omgs, alp);
    Ellipsoid3::AngNorm(bets, omgs, alp, _s.k2() == 0);
    gam = -alp;
    real x = _s.b() * F(_exs, omegashift(omgs, +1)),
      y = _s.b() * F(_eys, bets);
    real mb = sqrt(_s.k2() * Math::sq(bets.c()) +
                   _s.kp2() * Math::sq(omgs.s()));
    omg = omega(x);
    bet = beta(y);
    real ma = invscale(bet, omg);
    m = sphericalscale(ma, mb);
  }

  void Conformal3::ForwardOther(const Conformal3& alt, Angle bet, Angle omg,
                                Angle& betalt, Angle& omgalt,
                                Angle& gam, real& m) const {
    vec3 r, v;
    real ma, mb;
    ForwardSphere(bet, omg, r, v, ma);
    real f = alt._s.b()/_s.b();
    r[0] *= f; r[1] *= f; r[2] *= f;
    alt.ReverseSphere(r, v, betalt, omgalt, gam, mb);
    m = (ma/mb) * f;
  }

  void Conformal3::ReverseOther(const Conformal3& alt,
                                Angle betalt, Angle omgalt,
                                Angle& bet, Angle& omg,
                                Angle& gam, real& m) const {
    alt.ForwardOther(*this, betalt, omgalt, bet, omg, gam, m);
    m = 1/m; gam = -gam;
  }

  } // namespace Triaxial
} // namespace GeographicLib

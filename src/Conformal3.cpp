/**
 * \file Conformal3.cpp
 * \brief Implementation for GeographicLib::Triaxial::Conformal3 class
 *
 * Copyright (c) Charles Karney (2014-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <iostream>
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
    _s0 = Ellipsoid3(_s.b(), 0, 1, 0);
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

  Math::real Conformal3::Pi(const EllipticFunction& ell, Angle phi) {
    real p = ell.Pi(phi.s(), phi.c(), ell.Delta(phi.s(), phi.c()));
    if (phi.n() != 0) p += 4 * ell.Pi() * phi.n();
    return p;
  }
  Angle Conformal3::Piinv(const EllipticFunction& ell, real x) {
    real y = remainder(x, 2 * ell.Pi()),
      n = 2 * round((x - y) / (2 * ell.Pi()));
    // Now x = n * Pi() + y where y in [-Pi(), Pi()].  Pi() is the quarter
    // period for the elliptic integral which corresponds to pi/2 in angle
    // space.
    if (y == 0)
      return Angle::cardinal( n == 0 ? y : n ); // Preserve the sign of +/-0
    else if (fabs(y) == ell.Pi())
      return Angle::cardinal(copysign(real(1), y) + n);
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
      real z = Trigfun::root(Pif, fabs(y), fabs(y) * Math::pi()/(2*ell.Pi()),
                             0, Math::pi()/2,
                             1,1,1,
                             &countn, &countb,
                             0, Trigfun::PIINV);
      (void) countn; (void) countb;
      // cout << "CNT " << countn << " " << countb << "\n";
      return Angle::radians(copysign(z, y)) + Angle::cardinal(n);
    }
  }
  Math::real Conformal3::F(const EllipticFunction& ell, Angle phi) {
    real p = ell.F(phi.s(), phi.c(), ell.Delta(phi.s(), phi.c()));
    if (phi.n() != 0) p += 4 * ell.K() * phi.n();
    return p;
  }
  Angle Conformal3::Finv(const EllipticFunction& ell, real x) {
    real y = remainder(x, 2 * ell.K()),
      n = 2 * round((x - y) / (2 * ell.K()));
    // Now x = n * K() + y where y in [-K(), K()].  K() is the quarter
    // period for the elliptic integral which corresponds to pi/2 in angle
    // space.
    if (y == 0)
      return Angle::cardinal( n == 0 ? y : n ); // Preserve the sign of +/-0
    else if (fabs(y) == ell.K())
      return Angle::cardinal(copysign(real(1), y) + n);
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
      real z = Trigfun::root(Ff, fabs(y), fabs(y) * Math::pi()/(2*ell.K()),
                             0, Math::pi()/2,
                             1,1,1,
                             &countn, &countb,
                             0, Trigfun::FINV);
      (void) countn; (void) countb;
      // cout << "CNT " << countn << " " << countb << "\n";
      return Angle::radians(copysign(z, y)) + Angle::cardinal(n);
    }
  }

  Math::real Conformal3::x(Angle omg) const {
    omg = omegashift(omg, +1);
    return Math::sq(a()) / b() * Pi(_ex, omg.modang(b() / a()));
  }
  Angle Conformal3::omega(real x) const {
    Angle omg = Piinv(_ex, (x * b()) / Math::sq(a())).modang(a() / b());
    return omegashift(omg, -1);
  }

  Math::real Conformal3::y(Angle bet) const {
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
    Forward(bet, omg, x, y); m = scale(bet, omg);
  }
  void Conformal3::Reverse(real x, real y, Angle& bet, Angle& omg,
                           real& m) const {
    Reverse(x, y, bet, omg); m = scale(bet, omg);
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
      if (y == 0) { k2 = 0; break; }
      // Find initial guess assume K(k2) = pi/2, so K(kp2) = nx/ny * pi/2.
      // Invert using approximate k(K) given in
      // https://arxiv.org/abs/2505.17159v4
      real KK = nx/ny * Math::pi()/2;
      k2 = 16/pow(exp(N*KK) - B, 2/N);
      // Alternatively using KK = 1/2*log(16/kp) A+S 17.3.26
      k2 = fmin(1/real(2), 16*exp(-2*KK)); // Make sure guess is sane
      const bool logp = true;
      if constexpr (logp) {
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
          logk2 = Trigfun::root(ksolve, 0, logk2,
                                logk2min, -log(real(2)),
                                1, 1, 1,
                                &countn, &countb, 0, Trigfun::KINV);
        }
        // otherwise accept the asymptotic result
        k2 = exp(logk2);
      } else {
        auto ksolve = [nx, ny]
          (real k2) -> pair<real, real>
          {
            real kp2 = 1 - k2;
            EllipticFunction elly(k2), ellx(kp2, 0, k2, 1);
            real f = nx * elly.K() - ny * ellx.K(),
            fp = (nx * (elly.E() - kp2 * elly.K()) +
                  ny * (ellx.E() - k2  * ellx.K())) / (2 * k2 * kp2);
            return pair<real, real>(f, fp);
          };
        k2 = Trigfun::root(ksolve, 0, k2,
                           Math::sq(numeric_limits<real>::epsilon()),
                           1/real(2),
                           1, 1, 1,
                           &countn, &countb, 0, Trigfun::KINV);
      }
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
    if (0)
      cout << "CNT " << countn << " " << countb << " "
           << x/b - ellx.K() << " " << y/b - elly.K() << " "
           << (x * elly.K() - y * ellx.K())/b << "\n";

    return Ellipsoid3(b, 0, k2, kp2);
  }
  void Conformal3::ForwardSphere(Angle bet, Angle omg, Angle& phi, Angle& lam,
                                 Angle& gamma, real& m) const {
    real x, y;
    Forward(bet, omg, x, y, m);
    Angle omgs = omegashift(Finv(_exs, x/_s.b()), -1),
      bets = Finv(_eys, y/_s.b()),
      alp = Angle{};
    m *= sqrt(_s.k2() * Math::sq(bets.c()) + _s.kp2() * Math::sq(omgs.s()));
    if (0)
    cout << "QQ " << real(bet) << " " << real(omg) << " "
         << x << " " << y << " "
         << real(bets) << " " << real(omgs) << "\n";
    Ellipsoid3::vec3 r, v;
    _s.elliptocart2(bets, omgs, alp, r, v);
    _s0.cart2toellip(r, v, phi, lam, gamma);
  }
  void Conformal3::ReverseSphere(Angle phi, Angle lam, Angle& bet, Angle& omg,
                                 Angle& gamma, real& m) const {
    Ellipsoid3::vec3 r, v;
    Angle alp = Angle{};
    _s0.elliptocart2(phi, lam, alp, r, v);
    Angle bets, omgs;
    _s.cart2toellip(r, v, bets, omgs, alp);
    gamma = -alp;
    real x = _s.b() * F(_exs, omegashift(omgs, +1)),
      y = _s.b() * F(_eys, bets);
    real ms = sqrt(_s.k2() * Math::sq(bets.c()) +
                   _s.kp2() * Math::sq(omgs.s()));
    omg = omega(x);
    bet = beta(y);
    m = ms * scale(bet, omg);
  }
  } // namespace Triaxial
} // namespace GeographicLib

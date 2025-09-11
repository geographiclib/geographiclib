/**
 * \file TriaxialConformal.cpp
 * \brief Implementation for GeographicLib::TriaxialConformal class
 *
 * Copyright (c) Charles Karney (2014-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <iostream>
#include <GeographicLib/TriaxialConformal.hpp>
#include <GeographicLib/Trigfun.hpp>

namespace GeographicLib {

  using namespace std;

  TriaxialConformal::TriaxialConformal(const Triaxial& t)
    : _t(t)
    , _ex(kp2() * (1 - e2() * k2 ()),   - e2() * kp2(),
          k2 () * (1 + e2() * kp2()), 1 + e2() * kp2())
    , _ey(k2 () * (1 + e2() * kp2()),     e2() * k2 (),
          kp2() * (1 - e2() * k2 ()), 1 - e2() * k2 ())
    , _exalt(kp2() * (1 - e2() * k2 ()), kp2(),
             k2 () * (1 + e2() * kp2()), k2 ())
    , _eyalt(k2 () * (1 + e2() * kp2()), k2 (),
             kp2() * (1 - e2() * k2 ()), kp2())
  {}

  // Proving the equivalence of this and the Nyrtsov formulations is
  // straightforward but laborious.  Working from this formulation:
  // * switch the angle argument of Pi to its complement
  // * sin(phi)^2 -> cos(phi')^2 = 1 - sin(phi')^2
  // * modulus of Pi then becomes imaginary
  // * Use DLMF 19.7.5 to express as sum of Pi term and F term

  // Pi = int( 1/(1-a2*sin(t)^2)/sqrt(1-k2*sin(t)^2), t)
  // t = pi/2 - tx
  // Pi = -int( 1/(1-a2*cos(tx)^2)/sqrt(1-k2*cos(tx)^2), tx)
  //    = -int( 1/(1-a2*(1-sin(tx)^2))/sqrt(1-k2*(1-sin(tx)^2)), tx)
  //    = -int( 1/((1-a2)+a2*sin(tx)^2)/sqrt((1-k2)+k2*sin(tx)^2), tx)
  //    = - 1/(1-a2)/sqrt(1-k2) *
  //       int( 1/(1-a2x*sin(tx)^2)/sqrt(1+k2x*sin(tx)^2), tx)
  //    where a2x = -a2/(1-a2) k2x = k2/(1-k2)
  // So Pi(a2, k2) = Pi(t, a2, k) + Pi(tx, a2x, -k2x)/(1-a2)/sqrt(1-k2)
  // DMLMF 19.7.5
  //   a2xx = (a2x + k2x)/(1 + k2x)
  //   k2xx = k2x/(1+k2x);  kp2xx = 1/(1+k2x)
  // tan(txx)^2 = (1+k2x)*tan(tx)^2
  // Pi(tx, a2x, -k2x) = sqrt(kp2xx)/a2xx * (k2xx* F(txx,k2xx) + kp2xx*a2x*Pi(txx,a2xx,k2xx))
  //   [  a2xx = (a2x + k2x)/(1 + k2x),k2xx = k2x/(1+k2x),  kp2xx = 1/(1+k2x),(1+k2x)],a2x = -a2/(1-a2) k2x = k2/(1-k2);
  //   => a2xx =
  //     a2xx = (k2-a2)/(1-a2), k2xx = k2, kp2xx = 1 - k2, (1+k2x) = 1/(1-k2)
  //     [    a2xx = (k2-a2)/(1-a2), k2xx = k2,kp2xx = 1 - k2],
  //     k2 = kp2y * (1 - e2 * k2y),  a2= - e2 * kp2y;
  //   a2xx = kp2y
  //   k2xx = kp2y*(1-e2*k2y)
  //     tan(t)^2 = (1+e2*kp) * tan(ty)^2
  // tan2 multiplier = (1+e2*kp2y)/(1-k2),k2 = kp2y * (1 - e2 * k2y) = 1/k2y
  TriaxialConformal::TriaxialConformal(real a, real b, real c)
    : TriaxialConformal(Triaxial(a, b, c))
  {}
  TriaxialConformal::TriaxialConformal(real b, real e2, real k2, real kp2)
    : TriaxialConformal(Triaxial(b, e2, k2, kp2))
  {}

  Math::real TriaxialConformal::Pi(const EllipticFunction& ell, Angle phi) {
    real p = ell.Pi(phi.s(), phi.c(), ell.Delta(phi.s(), phi.c()));
    if (phi.n() != 0) p += 4 * ell.Pi() * phi.n();
    return p;
  }
  Angle TriaxialConformal::Piinv(const EllipticFunction& ell, real x) {
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
      real z = Trigfun::root([&ell]
                             (real phi) -> pair<real, real>
                             {
                               real s = sin(phi), c = cos(phi),
                                 f = ell.Pi(s, c, ell.Delta(s, c)),
                                 fp = 1 / (sqrt(ell.kp2() + ell.k2() * c*c) *
                                           (ell.alpha2() >= 0 ?
                                            ell.alphap2() + ell.alpha2() * s*s :
                                            1 - ell.alpha2() * c*c));
                               return pair<real, real>(f, fp);
                             },
                             fabs(y),
                             fabs(y) * Math::pi()/(2*ell.Pi()),
                             0, Math::pi()/2,
                             1,1,1,
                             &countn, &countb,
                             0, Trigfun::PIINV);
      (void) countn; (void) countb;
      // cout << "CNT " << countn << " " << countb << "\n";
      return Angle::radians(copysign(z, y)) + Angle::cardinal(n);
    }
  }

  Math::real TriaxialConformal::x() const {
    return Math::sq(a()) / b() * _ex.Pi();
  }
  Math::real TriaxialConformal::x(Angle omg) const {
    omg = omegashift(omg, +1);
    return Math::sq(a()) / b() * Pi(_ex, omg.modang(b() / a()));
  }
  Math::real TriaxialConformal::x2() const {
    return b() * e2() * k2() * _exalt.Pi() +
      Math::sq(c())/b() * _exalt.K();
  }
  Math::real TriaxialConformal::x2(Angle omg) const {
    Angle phi = omg.modang(1/sqrt(k2()));
    return b() * e2() * k2() * Pi(_exalt, phi) +
      Math::sq(c())/b() *
      (_exalt.F(phi.s(), phi.c(), _exalt.Delta(phi.s(), phi.c())) +
       (phi.n() != 0 ? 4 * _exalt.K() * phi.n() : 0))  - x2();
  }

  Angle TriaxialConformal::omega(real x) const {
    Angle omg = Piinv(_ex, (x * b()) / Math::sq(a())).modang(a() / b());
    return omegashift(omg, -1);
  }

  Math::real TriaxialConformal::y() const {
    return Math::sq(c()) / b() * _ey.Pi();
  }
  Math::real TriaxialConformal::y2() const {
    return -b() * e2() * kp2() * _eyalt.Pi() +
      Math::sq(a())/b() * _eyalt.K();
  }
  Math::real TriaxialConformal::y(Angle bet) const {
    return  Math::sq(c()) / b() * Pi(_ey, bet.modang(b() / c()));
  }
  Math::real TriaxialConformal::y2(Angle bet) const {
    bet += Angle::cardinal(1);
    Angle phi = bet.modang(1/sqrt(kp2()));
    return -b() * e2() * kp2() * Pi(_eyalt, phi) +
      Math::sq(a())/b() *
      (_eyalt.F(phi.s(), phi.c(), _eyalt.Delta(phi.s(), phi.c())) +
       4 * _eyalt.K() * phi.n()) - y2();
  }
  Angle TriaxialConformal::beta(real y) const {
    return Piinv(_ey, (y * b()) / Math::sq(c())).modang(c() / b());
  }

  void TriaxialConformal::Forward(Angle bet, Angle omg, real& xa, real& ya)
    const {
    xa = x(omg); ya = y(bet);
  }
  void TriaxialConformal::Reverse(real x, real y, Angle& bet, Angle& omg)
    const {
    bet = beta(y); omg = omega(x);
  }
  void TriaxialConformal::Forward(Angle bet, Angle omg, real& x, real& y,
                                  real& m) const {
    Forward(bet, omg, x, y); m = scale(bet, omg);
  }
  void TriaxialConformal::Reverse(real x, real y, Angle& bet, Angle& omg,
                                  real& m) const {
    Reverse(x, y, bet, omg); m = scale(bet, omg);
  }
} // namespace GeographicLib

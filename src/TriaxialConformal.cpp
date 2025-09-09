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
  {}

  TriaxialConformal::TriaxialConformal(real a, real b, real c)
    : TriaxialConformal(Triaxial(a, b, c))
  {}
  TriaxialConformal::TriaxialConformal(real b, real e2, real k2, real kp2)
    : TriaxialConformal(Triaxial(b, e2, k2, kp2))
  {}

  Math::real TriaxialConformal::Pi(const EllipticFunction& ell, Angle phi) {
    return
      4 * ell.Pi() * phi.n() + ell.Pi(phi.s(), phi.c(),
                                      ell.Delta(phi.s(), phi.c()));
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
    real x = Math::sq(a()) / b() * Pi(_ex, omg.modang(b() / a()));
    return xshift(x, +1);
  }
  Angle TriaxialConformal::omega(real x) const {
    x = xshift(x, -1);
    Angle omg = Piinv(_ex, (x * b()) / Math::sq(a())).modang(a() / b());
    return omegashift(omg, -1);
  }

  Math::real TriaxialConformal::y() const {
    return Math::sq(c()) / b() * _ey.Pi();
  }
  Math::real TriaxialConformal::y(Angle bet) const {
    return Math::sq(c()) / b() * Pi(_ey, bet.modang(b() / c()));
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

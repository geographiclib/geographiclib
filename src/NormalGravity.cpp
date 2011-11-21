/**
 * \file NormalGravity.cpp
 * \brief Implementation for GeographicLib::NormalGravity class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/NormalGravity.hpp>

#define GEOGRAPHICLIB_NORMALGRAVITY_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_NORMALGRAVITY_CPP)
RCSID_DECL(GEOGRAPHICLIB_NORMALGRAVITY_HPP)

#include <iostream>
#include <iomanip>

namespace GeographicLib {

  using namespace std;

  NormalGravity::NormalGravity(real a, real GM, real omega, real J2, bool flatp)
    : _a(a)
    , _GM(GM)
    , _omega(omega)
    , _J2(J2)
    , _omega2(Math::sq(_omega))
    , _aomega2(Math::sq(_omega * _a))
    {
      if (!(Math::isfinite(_a) && _a > 0))
        throw GeographicErr("Major radius is not positive");
      if (!(Math::isfinite(_GM) && _GM > 0))
        throw GeographicErr("Gravitational constants is not positive");
      if (!(Math::isfinite(_J2) && _J2 > 0))
        throw GeographicErr(flatp ? "Flattening is not positive"
                            : "Dynamical form factor is not positive");
      if (!(Math::isfinite(_omega) && _omega != 0))
        throw GeographicErr("Angular velocity is not non-zero");
      real K = 2 * _aomega2 * _a / (15 * _GM);
      if (flatp) {
        _f = _J2;
        _e2 = _f * (2 - _f);
        _ep2 = _e2 / (1 - _e2);
        _q0 = qf(_ep2);
        _J2 = _e2 * ( 1 - K * sqrt(_e2) / _q0) / 3; // H+M, Eq 2-90
      } else {
        _e2 = 3 * _J2;          // See Moritz (1980), p 398.
        for (int j = 0; j < maxit_; ++j) {
          real e2a = _e2;
          real q0 = qf(_e2 / (1 - _e2));
          _e2 = 3 * _J2 + K * _e2 * sqrt(_e2) / q0;
          if (_e2 == e2a)
            break;
        }
        _f = _e2 / (1 + sqrt(1 - _e2));
        _ep2 = _e2 / (1 - _e2);
        _q0 = qf(_ep2);
      }
      _earth = Geocentric(_a, _f);
      _b = _a * (1 - _f);
      _E = a * sqrt(_e2);                               // H+M, Eq 2-54
      _U0 = _GM / _E * atan(sqrt(_ep2)) + _aomega2 / 3; // H+M, Eq 2-61
      _m = _aomega2 * _b / _GM;                         // H+M, Eq 2-70
      real
        Q = _m * sqrt(_ep2) * qpf(_ep2) / (3 * _q0),
        G = (1 - _m - Q / 2);
      _gammae = _GM / (_a * _b) * G;       // H+M, Eq 2-73
      _gammap = _GM / (_a * _a) * (1 + Q); // H+M, Eq 2-74
      // k = b * gammap / (a * gammae) - 1
      _k = (_m + 3 * Q / 2 - _e2 * (1 + Q)) / G;
      // f* = (gammap - gammae) / gammae
      _fstar = (_m + 3 * Q / 2 - _f * (1 + Q)) / G;
    }

  const NormalGravity
  NormalGravity::GRS80(Constants::GRS80_a<real>(), Constants::GRS80_GM<real>(),
                       Constants::GRS80_omega<real>(),
                       Constants::GRS80_J2<real>(), false);

  const NormalGravity
  NormalGravity::WGS84(Constants::WGS84_a<real>(), Constants::WGS84_GM<real>(),
                       Constants::WGS84_omega<real>(),
                       Constants::WGS84_f<real>(), true);

  Math::real NormalGravity::qf(real ep2) throw() {
    // Compute
    //
    //  ((1 + 3/e'^2) * atan(e')  - 3/e')/2
    //
    // See H+M, Eq 2-57, with E/u = e'.  This suffers from two levels of
    // cancelation.  The e'^-1 and e'^1 terms drop out, so that the leading
    // term is O(e'^3).
    real ep = sqrt(ep2);
    if (abs(ep2) >  real(0.5))  // Use the closed expression
      return ((1 + 3 / ep2) * atan(ep)  - 3 / ep)/2;
    else {
      real ep2n = 1, q = 0;     // The series expansion H+M, Eq 2-86
      for (int n = 1; ; ++n) {
        ep2n *= -ep2;
        real
          t = (ep2n * n) / ((2 * n + 1) * (2 * n + 3)),
          qn = q + t;
        if (qn == q)
          break;
        q = qn;
      }
      q *= -2 * ep;
      return q;
    }
  }

  Math::real NormalGravity::qpf(real ep2) throw() {
    // Compute
    //
    //  3*(1 + 1/e'^2) * (1 - atan(e')/e') - 1
    //
    // See H+M, Eq 2-67, with E/u = e'.  This suffers from two levels of
    // cancelation.  The e'^-2 and e'^0 terms drop out, so that the leading
    // term is O(e'^2).
    if (abs(ep2) > real(0.5)) { // Use the closed expression
      real ep = sqrt(ep2);
      return 3 * (1 + 1 / ep2) * (1 - atan(ep) / ep) - 1;
    } else {
      real ep2n = 1, qp = 0;    // The series expansion H+M, Eq 2-101c
      for (int n = 1; ; ++n) {
        ep2n *= -ep2;
        real
          t = ep2n / ((2 * n + 1) * (2 * n + 3)),
          qpn = qp + t;
        if (qpn == qp)
          break;
        qp = qpn;
      }
      qp *= -6;
      return qp;
    }
  }
      
  Math::real NormalGravity::Jn(int n) const throw() {
    // Note Jn(0) = -1; Jn(2) = _J2; Jn(odd) = 0
    if (n & 1 || n < 0)
      return 0;
    n /= 2;
    real e2n = 1;            // Perhaps this should just be e2n = pow(-_e2, n);
    for (int j = n; j--;)
      e2n *= -_e2;
    return                      // H+M, Eq 2-92
      -3 * e2n * (1 - n + 5 * n * _J2 / _e2) / ((2 * n + 1) * (2 * n + 3));
  }

  Math::real NormalGravity::SurfaceGravity(real lat) const throw() {
    real
      phi = lat * Math::degree<real>(),
      sphi2 = abs(lat) == 90 ? 1 : Math::sq(sin(phi));
    // H+M, Eq 2-78
    return _gammae * (1 + _k * sphi2) / sqrt(1 - _e2 * sphi2);
  }

  Math::real NormalGravity::V0(real X, real Y, real Z,
                               real& GammaX, real& GammaY, real& GammaZ)
    const throw() {
    // See H+M, Sec 6-2
    real
      p = Math::hypot(X, Y),
      clam = p ? X/p : 1,
      slam = p ? Y/p : 0,
      r = Math::hypot(p, Z),
      Q = Math::sq(r) - Math::sq(_E),
      t2 = Math::sq(2 * _E * Z),
      disc = sqrt(Math::sq(Q) + t2),
      // This is H+M, Eq 6-8a, but generalized to deal with Q negative
      // accurately.
      u = sqrt((Q >= 0 ? (Q + disc) : t2 / (disc - Q)) / 2),
      uE = Math::hypot(u, _E),
      // H+M, Eq 6-8b
      sbet = Z * uE,
      cbet = p * u,
      s = Math::hypot(cbet, sbet);
    cbet = s ? cbet/s : 0;
    sbet = s ? sbet/s : 1;
    real
      invw = uE / Math::hypot(u, _E * sbet), // H+M, Eq 2-63
      ep = _E/u,
      ep2 = Math::sq(ep),
      q = qf(ep2) / _q0,
      qp = qpf(ep2) / _q0,
      // H+M, Eqs 2-62 + 6-9, but omitting last (rotational) term .
      Vres = (_GM / _E * atan(_E / u)
              + _aomega2 * q * (Math::sq(sbet) - 1/real(3)) / 2),
      // H+M, Eq 6-10
      gamu = - invw * (_GM
                       + (_aomega2 * _E * qp
                          * (Math::sq(sbet) - 1/real(3)) / 2)) / Math::sq(uE),
      gamb = _aomega2 * q  * sbet * cbet * invw / uE,
      t = u * invw / uE;
    // H+M, Eq 6-12
    GammaX = t * cbet * gamu - invw * sbet * gamb;
    GammaY = GammaX * slam;
    GammaX *= clam;
    GammaZ = invw * sbet * gamu + t * cbet * gamb;
    return Vres;
  }

  Math::real NormalGravity::Phi(real X, real Y, real& fX, real& fY)
    const throw() {
    fX = _omega2 * X;
    fY = _omega2 * Y;
    // N.B. fZ = 0;
    return _omega2 * (Math::sq(X) + Math::sq(Y)) / 2;
  }

  Math::real NormalGravity::U(real X, real Y, real Z,
                              real& gammaX, real& gammaY, real& gammaZ)
    const throw() {
    real fX, fY;
    real Ures = V0(X, Y, Z, gammaX, gammaY, gammaZ) + Phi(X, Y, fX, fY);
    gammaX += fX;
    gammaY += fY;
    return Ures;
  }
  
  Math::real NormalGravity::Gravity(real lat, real h, real& gammay, real& gammaz)
    const throw() {
    real X, Y, Z;
    real M[Geocentric::dim2_];
    _earth.IntForward(lat, 0, h, X, Y, Z, M);
    real gammaX, gammaY, gammaZ;
    U(X, Y, Z, gammaX, gammaY, gammaZ);
    // gammax = M[0] * gammaX + M[3] * gammaY + M[6] * gammaZ;
    gammay = M[1] * gammaX + M[4] * gammaY + M[7] * gammaZ;
    gammaz = M[2] * gammaX + M[5] * gammaY + M[8] * gammaZ;
    return Math::hypot(gammay, gammaz);
  }

  void NormalGravity::DumpConstants() const {
    cout << setprecision(16)
         << "a=" << _a << "\n"
         << "GM=" << _GM << "\n"
         << "omega=" << _omega << "\n"
         << "J(2)=" << _J2 << "\n"
         << "J(4)=" << Jn(4) << "\n"
         << "J(6)=" << Jn(6) << "\n"
         << "J(8)=" << Jn(8) << "\n"
         << "J(10)=" << Jn(10) << "\n"
         << "omega2=" << _omega2 << "\n"
         << "aomega2=" << _aomega2 << "\n"
         << "e2=" << _e2 << "\n"
         << "ep2=" << _ep2 << "\n"
         << "f=" << _f << "\n"
         << "1/f=" << 1/_f << "\n"
         << "b=" << _b << "\n"
         << "E=" << _E << "\n"
         << "U0=" << _U0 << "\n"
         << "gammae=" << _gammae << "\n"
         << "gammap=" << _gammap << "\n"
         << "q0=" << _q0 << "\n"
         << "m=" << _m << "\n"
         << "k=" << _k << "\n"
         << "fstar=" << _fstar << "\n"
         << "1/fstar=" << 1/_fstar << "\n";
    for (int n = 2; n <= 20; n += 2)
      cout << "C(" << n << ",0)=" << -Jn(n)/sqrt(2*n+1) << "\n";
  }

} // namespace GeographicLib

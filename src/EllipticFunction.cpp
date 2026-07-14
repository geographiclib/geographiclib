/**
 * \file EllipticFunction.cpp
 * \brief Implementation for GeographicLib::EllipticFunction class
 *
 * Copyright (c) Charles Karney (2008-2024) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/Trigfun.hpp>
#include <type_traits>

namespace GeographicLib {

  using namespace std;

  /*
   * Implementation of methods given in
   *
   *   B. C. Carlson
   *   Computation of elliptic integrals
   *   Numerical Algorithms 10, 13-26 (1995)
   */

  template<typename T>
  T EllipticFunction::RFt(T x, T y, T z) {
    if (x == real(0))
      return RF(y, z);
    else if (y == real(0))
      return RF(z, x);
    else if (z == real(0))
      return RF(x, y);
    else if (y == z)
      return RC(x, y);
    else if (z == x)
      return RC(y, z);
    else if (x == y)
      return RC(z, x);
    // Carlson, eqs 2.2 - 2.7
    static const real tolRF =
      pow(3 * numeric_limits<real>::epsilon() * real(0.01), 1/real(8));
    T
      A0 = (x + y + z) / real(3),
      An = A0,
      x0 = x,
      y0 = y,
      z0 = z;
    real
      Q = fmax(fmax(Abs(A0-x), Abs(A0-y)), Abs(A0-z)) / tolRF,
      mul = 1;
    while (Q >= mul * Abs(An)) {
      // Max 6 trips
      T lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0);
      An = (An + lam) / real(4);
      x0 = (x0 + lam) / real(4);
      y0 = (y0 + lam) / real(4);
      z0 = (z0 + lam) / real(4);
      mul *= 4;
    }
    T
      X = (A0 - x) / (mul * An),
      Y = (A0 - y) / (mul * An),
      Z = - (X + Y),
      E2 = X*Y - Z*Z,
      E3 = X*Y*Z;
    // https://dlmf.nist.gov/19.36.E1
    // Polynomial is
    // (1 - E2/10 + E3/14 + E2^2/24 - 3*E2*E3/44
    //    - 5*E2^3/208 + 3*E3^2/104 + E2^2*E3/16)
    // convert to Horner form...
    return (E3 * (real(6930) * E3 + E2 * (real(15015) * E2 - real(16380)) +
                  real(17160)) +
            E2 * ((real(10010) - real(5775) * E2) * E2 - real(24024)) +
            real(240240)) /
      (real(240240) * sqrt(An));
  }
  Math::real EllipticFunction::RF(real x, real y, real z) {
    return RFt(x, y, z);
  }
  Math::cmplx EllipticFunction::RF(cmplx x, cmplx y, cmplx z) {
    return RFt(x, y, z);
  }

  template<typename T>
  T EllipticFunction::RFt(T x, T y) {
    // Carlson, eqs 2.36 - 2.38
    static const real tolRG0 =
      real(2.7) * sqrt((numeric_limits<real>::epsilon() * real(0.01)));
    T xn = sqrt(x), yn = sqrt(y);
    if (Abs(xn) < Abs(yn)) swap(xn, yn);
    while (Abs(xn-yn) > tolRG0 * Abs(xn)) {
      // Max 4 trips
      T t = (xn + yn) / real(2);
      yn = sqrt(xn * yn);
      xn = t;
    }
    return Math::pi() / (xn + yn);
  }
  Math::real EllipticFunction::RF(real x, real y) {
    return RFt(x, y);
  }
  Math::cmplx EllipticFunction::RF(cmplx x, cmplx y) {
    return RFt(x, y);
  }

  template<typename T>
  T EllipticFunction::RCt(T x, T y) {
    if constexpr (is_same_v<T, Math::real>) {
      // Defined only for y != 0 and x >= 0.
      return ( !(x >= y) ?        // x < y  and catch nans
               // https://dlmf.nist.gov/19.2.E18
               atan(sqrt((y - x) / x)) / sqrt(y - x) :
               ( x == y ? 1 / sqrt(y) :
                 asinh( y > 0 ?
                        // https://dlmf.nist.gov/19.2.E19
                        // atanh(sqrt((x - y) / x))
                        sqrt((x - y) / y) :
                        // https://dlmf.nist.gov/19.2.E20
                        // atanh(sqrt(x / (x - y)))
                        sqrt(-x / y) ) / sqrt(x - y) ) );
    } else {
      static const real tolRC =
        pow(3 * numeric_limits<real>::epsilon() * real(0.01), 1/real(8));
      if (y.imag() == 0 && signbit(y.real()))
        // Carlson, eq 2.14
        return sqrt(x / (x - y)) * RC(x - y, -y);
      // Carlson, eqs 2.9 - 2.13
      cmplx x0 = x, y0 = y, A0 = (x + real(2) * y) / real(3), An = A0;
      real Q = abs(A0 - x) / tolRC,
        mul = 1;
      while (Q >= mul * abs(An)) {
        cmplx lam = real(2) * sqrt(x0) * sqrt(y0) + y0;
        An = (An + lam) / real(4);
        x0 = (x0 + lam) / real(4);
        y0 = (y0 + lam) / real(4);
        mul *= 4;
      }
      cmplx s = (y - A0) / (mul * An);
      // series is
      // 1 + 3/10*s^2 + 1/7*s^3 + 3/8*s^4 + 9/22*s^5 + 159/208*s^6 + 9/8*s^7
      // Write in Horner form
      return (s*s*(s*(s*(s*(s*(real(90090)*s + real(61215)) + real(32760)) +
                         real(30030)) + real(11440)) + real(24024)) +
              real(80080)) / (real(80080) * sqrt(An));
    }
  }

  Math::real EllipticFunction::RC(real x, real y) {
    return RCt(x, y);
  }
  Math::cmplx EllipticFunction::RC(cmplx x, cmplx y) {
    return RCt(x, y);
  }

  template<typename T>
  T EllipticFunction::RGt(T x, T y, T z) {
    return (x == real(0) ? RG(y, z) :
            (y == real(0) ? RG(z, x) :
             (z == real(0) ? RG(x, y) :
              // Carlson, eq 1.7
              (z * RF(x, y, z) - (x-z) * (y-z) * RD(x, y, z) / real(3)
               + sqrt(x * y / z)) / real(2) )));
  }
  Math::real EllipticFunction::RG(real x, real y, real z) {
    return RGt(x, y, z);
  }
  Math::cmplx EllipticFunction::RG(cmplx x, cmplx y, cmplx z) {
    return RGt(x, y, z);
  }

  template<typename T>
  T EllipticFunction::RGt(T x, T y) {
    if (x == real(0))
      // Carlson, top of p. 21
      return sqrt(y) / real(2);
    else if (y == real(0))
      return sqrt(x) / real(2);
    // Carlson, eqs 2.36 - 2.39
    static const real tolRG0 =
      real(2.7) * sqrt((numeric_limits<real>::epsilon() * real(0.01)));
    T
      x0 = sqrt(x),
      y0 = sqrt(y);
    if (Abs(x) < Abs(y)) swap(x0, y0);
    T
      xn = x0,
      yn = y0,
      s = real(0);
    real
      mul = real(0.25);
    while (Abs(xn-yn) > tolRG0 * Abs(xn)) {
      // Max 4 trips
      T t = (xn + yn) / real(2);
      yn = sqrt(xn * yn);
      xn = t;
      mul *= 2;
      t = xn - yn;
      s += mul * t * t;
    }
    return (Math::sq( (x0 + y0) / real(2) ) - s) *
      Math::pi() / (real(2) * (xn + yn));
  }
  Math::real EllipticFunction::RG(real x, real y) {
    return RGt(x, y);
  }
  Math::cmplx EllipticFunction::RG(cmplx x, cmplx y) {
    return RGt(x, y);
  }

  template<typename T>
  T EllipticFunction::RJt(T x, T y, T z, T p) {
    if (x == p)
      return RD(y, z, p);
    else if (y == p)
      return RD(z, x, p);
    else if (z == p)
      return RD(x, y, p);
    // Carlson, eqs 2.17 - 2.25
    static const real
      tolRD = pow(real(0.2) * (numeric_limits<real>::epsilon() * real(0.01)),
                  1/real(8));
    if constexpr (is_same_v<T, Math::cmplx>) {
      if (p.imag() == 0 && signbit(p.real()) &&
          (x.imag() == 0 && y.imag() == 0 && z.imag() == 0)) {
        return RJ(x.real(), y.real(), z.real(), p.real());
      }
    }
    if constexpr (is_same_v<T, Math::real>) {
      // Carlson, eq 2.26
      if (signbit(p)) {
        if (!((z - y) * (y - x) < 0)) { // Backwards test to catch nan
          real q = -p, pmy = (z - y) * (y - x) / (y + q);
          p = pmy + y;
          return (pmy * RJ(x, y, z, p) - 3 * RF(x, y, z)
                  + 3* sqrt( x*y*z / (x*z + p*q) ) * RC(x*z + p*q, p*q)) /
            (y + q);
        } else
          // rotate the positions of the symmetric args so that eventually y
          // lies between x and z.
          return RJ(y, z, x, p);
      }
    }
    T
      A0 = (x + y + z + real(2) * p) / real(5),
      An = A0,
      delta = (p-x) * (p-y) * (p-z),
      x0 = x,
      y0 = y,
      z0 = z,
      p0 = p,
      s = real(0);
    real
      Q = fmax(fmax(Abs(A0-x), Abs(A0-y)),
               fmax(Abs(A0-z), Abs(A0-p))) / tolRD,
      mul = 1,
      mul3 = 1;
    while (Q >= mul * Abs(An)) {
      // Max 7 trips
      T
        lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0),
        d0 = (sqrt(p0)+sqrt(x0)) * (sqrt(p0)+sqrt(y0)) * (sqrt(p0)+sqrt(z0)),
        e0 = delta/(mul3 * Math::sq(d0));
      s += RC(real(1), real(1) + e0) / (mul * d0);
      An = (An + lam) / real(4);
      x0 = (x0 + lam) / real(4);
      y0 = (y0 + lam) / real(4);
      z0 = (z0 + lam) / real(4);
      p0 = (p0 + lam) / real(4);
      mul *= 4;
      mul3 *= 64;
    }
    T
      X = (A0 - x) / (mul * An),
      Y = (A0 - y) / (mul * An),
      Z = (A0 - z) / (mul * An),
      P = -(X + Y + Z) / real(2),
      P2 = P*P,
      XYZ = X*Y*Z,
      E2 = X*Y + X*Z + Y*Z - real(3)*P2,
      E3 = XYZ + real(2)*P * (E2 + real(2)*P2),
      E4 = (real(2)*XYZ + P * (E2 + real(3)*P2)) * P,
      E5 = XYZ*P2;
    // https://dlmf.nist.gov/19.36.E2
    // Polynomial is
    // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
    //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
    //    - 9*(E3*E4+E2*E5)/68)
    return ((real(471240) - real(540540) * E2) * E5 +
            (real(612612) * E2 - real(540540) * E3 - real(556920)) * E4 +
            E3 * (real(306306) * E3 +
                  E2 * (real(675675) * E2 - real(706860)) + real(680680)) +
            E2 * ((real(417690) - real(255255) * E2) * E2 - real(875160)) +
            real(4084080)) /
      ((4084080 * mul) * An * sqrt(An)) + real(6) * s;
  }
  Math::real EllipticFunction::RJ(real x, real y, real z, real p) {
    return RJt(x, y, z, p);
  }
  Math::cmplx EllipticFunction::RJ(cmplx x, cmplx y, cmplx z, cmplx p) {
    return RJt(x, y, z, p);
  }

  template<typename T>
  T EllipticFunction::RDt(T x, T y, T z) {
    if (x == real(0))
      // Carlson eqs. 2.40 - 2.41
      return y == z ? (3 * Math::pi() / 4) / (sqrt(y) * y) :
        real(3) / (z * (y - z)) * (real(2) * RG(y, z) - z * RF(y, z));
    else if (y == real(0))
      // Put 0 in x position
      return RD(y, x, z);
    // Carlson, eqs 2.28 - 2.34
    static const real
      tolRD = pow(real(0.2) * (numeric_limits<real>::epsilon() * real(0.01)),
                  1/real(8));
    T
      A0 = (x + y + real(3) * z) / real(5),
      An = A0,
      x0 = x,
      y0 = y,
      z0 = z,
      s = real(0);
    real
      Q = fmax(fmax(Abs(A0-x), Abs(A0-y)), Abs(A0-z)) / tolRD,
      mul = 1;
    while (Q >= mul * Abs(An)) {
      // Max 7 trips
      T lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0);
      s += real(1) / (mul * sqrt(z0) * (z0 + lam));
      An = (An + lam) / real(4);
      x0 = (x0 + lam) / real(4);
      y0 = (y0 + lam) / real(4);
      z0 = (z0 + lam) / real(4);
      mul *= 4;
    }
    T
      X = (A0 - x) / (mul * An),
      Y = (A0 - y) / (mul * An),
      Z = -(X + Y) / real(3),
      Z2 = Z*Z,
      XY = X*Y,
      E2 = XY - real(6)*Z2,
      E3 = (real(3)*XY - real(8)*Z2)*Z,
      E4 = real(3) * (XY - Z2) * Z2,
      E5 = XY*Z2*Z;
    // https://dlmf.nist.gov/19.36.E2
    // Polynomial is
    // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
    //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
    //    - 9*(E3*E4+E2*E5)/68)
    return ((real(471240) - real(540540) * E2) * E5 +
            (real(612612) * E2 - real(540540) * E3 - real(556920)) * E4 +
            E3 * (real(306306) * E3 + E2 * (real(675675) * E2 - real(706860)) +
                  real(680680)) +
            E2 * ((real(417690) - real(255255) * E2) * E2 - real(875160)) +
            real(4084080)) /
      ((4084080 * mul) * An * sqrt(An)) + real(3) * s;
  }
  Math::real EllipticFunction::RD(real x, real y, real z) {
    return RDt(x, y, z);
  }
  Math::cmplx EllipticFunction::RD(cmplx x, cmplx y, cmplx z) {
    return RDt(x, y, z);
  }

  void EllipticFunction::Reset(real k2, real alpha2,
                               real kp2, real alphap2) {
    // Accept nans here (needed for GeodesicExact)
    if (k2 > 1)
      throw GeographicErr("Parameter k2 is not in (-inf, 1]");
    if (alpha2 > 1)
      throw GeographicErr("Parameter alpha2 is not in (-inf, 1]");
    if (kp2 < 0)
      throw GeographicErr("Parameter kp2 is not in [0, inf)");
    if (alphap2 < 0)
      throw GeographicErr("Parameter alphap2 is not in [0, inf)");
    _k2 = k2;
    _kp2 = kp2;
    _alpha2 = alpha2;
    _alphap2 = alphap2;
    _eps = _k2/Math::sq(sqrt(_kp2) + 1);
    // Values of complete elliptic integrals for k = 0,1 and alpha = 0,1
    //         K     E     D
    // k = 0:  pi/2  pi/2  pi/4
    // k = 1:  inf   1     inf
    //                    Pi    G     H
    // k = 0, alpha = 0:  pi/2  pi/2  pi/4
    // k = 1, alpha = 0:  inf   1     1
    // k = 0, alpha = 1:  inf   inf   pi/2
    // k = 1, alpha = 1:  inf   inf   inf
    //
    // Pi(0, k) = K(k)
    // G(0, k) = E(k)
    // H(0, k) = K(k) - D(k)
    // Pi(0, k) = K(k)
    // G(0, k) = E(k)
    // H(0, k) = K(k) - D(k)
    // Pi(alpha2, 0) = pi/(2*sqrt(1-alpha2))
    // G(alpha2, 0) = pi/(2*sqrt(1-alpha2))
    // H(alpha2, 0) = pi/(2*(1 + sqrt(1-alpha2)))
    // Pi(alpha2, 1) = inf
    // H(1, k) = K(k)
    // G(alpha2, 1) = H(alpha2, 1) = RC(1, alphap2)
    _kKc = K_static(_k2, _kp2);
    if (_k2 != 0) {
      // Complete elliptic integral E(k), Carlson eq. 4.2
      // https://dlmf.nist.gov/19.25.E1
      _eEc = _kp2 != 0 ? 2 * RG(_kp2, 1) : 1;
      // D(k) = (K(k) - E(k))/k^2, Carlson eq.4.3
      // https://dlmf.nist.gov/19.25.E1
      _dDc = _kp2 != 0 ? RD(0, _kp2, 1) / 3 : Math::infinity();
    } else {
      _eEc = _kKc; _dDc = _kKc/2;
    }
    if (_alpha2 != 0) {
      // https://dlmf.nist.gov/19.25.E2
      real rj = (_kp2 != 0 && _alphap2 != 0) ? RJ(0, _kp2, 1, _alphap2) :
        Math::infinity(),
        // Only use rc if _kp2 = 0.
        rc = _kp2 != 0 ? 0 :
        (_alphap2 != 0 ? RC(1, _alphap2) : Math::infinity());
      // Pi(alpha^2, k)
      _pPic = _kp2 != 0 ? _kKc + _alpha2 * rj / 3 : Math::infinity();
      // G(alpha^2, k)
      _gGc = _kp2 != 0 ? _kKc + (_alpha2 - _k2) * rj / 3 :  rc;
      // H(alpha^2, k)
      _hHc = _kp2 != 0 ? _kKc - (_alphap2 != 0 ? _alphap2 * rj : 0) / 3 : rc;
    } else {
      _pPic = _kKc; _gGc = _eEc;
      // Hc = Kc - Dc but this involves large cancellations if k2 is close to
      // 1.  So write (for alpha2 = 0)
      //   Hc = int(cos(phi)^2/sqrt(1-k2*sin(phi)^2),phi,0,pi/2)
      //      = 1/sqrt(1-k2) * int(sin(phi)^2/sqrt(1-k2/kp2*sin(phi)^2,...)
      //      = 1/kp * D(i*k/kp)
      // and use D(k) = RD(0, kp2, 1) / 3
      // so Hc = 1/kp * RD(0, 1/kp2, 1) / 3
      //       = kp2 * RD(0, 1, kp2) / 3
      // using https://dlmf.nist.gov/19.20.E18
      // Equivalently
      //   RF(x, 1) - RD(0, x, 1)/3 = x * RD(0, 1, x)/3 for x > 0
      // For k2 = 0 and alpha2 = 0, we have
      //   Hc = int(cos(phi)^2,...) = pi/4
      // For k2 = 1 and alpha2 = 0, we have
      //   Hc = int(cos(phi),...) = 1
      _hHc = _kp2 == 1 ? Math::pi()/4 :
        (_kp2 == 0 ? 1 : _kp2 * RD(0, 1, _kp2) / 3);
    }
  }
  Math::real EllipticFunction::K_static(real k2, real kp2) {
    // Complete elliptic integral K(k), Carlson eq. 4.1
    // https://dlmf.nist.gov/19.25.E1
    return k2 == 0 ? Math::pi()/2 : kp2 == 0 ? Math::infinity() : RF(kp2, 1);
  }

  /*
   * Implementation of methods given in
   *
   *   R. Bulirsch
   *   Numerical Calculation of Elliptic Integrals and Elliptic Functions
   *   Numericshe Mathematik 7, 78-90 (1965)
   */

  void EllipticFunction::sncndn_static(real x, real& sn, real& cn, real& dn,
                                       real k2, real kp2) {
    // Bulirsch's sncndn routine, p 89.
    static const real tolJAC =
      sqrt(numeric_limits<real>::epsilon() * real(0.01));
    if (kp2 != 0) {
      real mc = kp2, d = 0;
      if (signbit(kp2)) {
        // This implements DLMF Eqs 22.17.2 - 22.17.4.  We *do* need this to
        // treat complex arguments with k2 < 0, kp2 > 1 since k2 and kp2 get
        // interchanged for the imaginary.  Then k2 > 1, kp2 < 0, and, e.g.,
        // sn(z, k) = sn(z*k,1/k)/k, etc.
        d = k2;
        mc /= -d;
        d = sqrt(d);
        x *= d;
      }
      real c = 0;           // To suppress warning about uninitialized variable
      real m[num_], n[num_];
      unsigned l = 0;
      for (real a = 1;
           l < num_ ||
             GEOGRAPHICLIB_PANIC
             ("Convergence failure in EllipticFunction::sncndn");
           ++l) {
        // This converges quadratically.  Max 5 trips
        m[l] = a;
        n[l] = mc = sqrt(mc);
        c = (a + mc) / 2;
        if (!(fabs(a - mc) > tolJAC * a)) {
          ++l;
          break;
        }
        mc *= a;
        a = c;
      }
      x *= c;
      sn = sin(x);
      cn = cos(x);
      dn = 1;
      if (sn != 0) {
        real a = cn / sn;
        c *= a;
        while (l--) {
          real b = m[l];
          a *= c;
          c *= dn;
          dn = (n[l] + a) / (b + a);
          a = c / b;
        }
        a = 1 / sqrt(c*c + 1);
        sn = signbit(sn) ? -a : a;
        cn = c * sn;
        if (signbit(kp2)) {
          // See DLMF Eqs 22.17.2 - 22.17.4
          swap(cn, dn);
          sn /= d;
        }
      }
    } else {
      sn = tanh(x);
      dn = cn = 1 / cosh(x);
    }
  }

  void EllipticFunction::sncndn(real x, real& sn, real& cn, real& dn) const {
    // Bulirsch's sncndn routine, p 89.
    static const real tolJAC =
      sqrt(numeric_limits<real>::epsilon() * real(0.01));
    if (_kp2 != 0) {
      real mc = _kp2, d = 0;
      if (signbit(_kp2)) {
        // This implements DLMF Eqs 22.17.2 - 22.17.4.  But this only
        // accommodates kp2 < 0 or k2 > 1 and these are outside the advertised
        // ranges for the contructor for this class.
        d = _k2;
        mc /= -d;
        d = sqrt(d);
        x *= d;
      }
      real c = 0;           // To suppress warning about uninitialized variable
      real m[num_], n[num_];
      unsigned l = 0;
      for (real a = 1;
           l < num_ ||
             GEOGRAPHICLIB_PANIC
             ("Convergence failure in EllipticFunction::sncndn");
           ++l) {
        // This converges quadratically.  Max 5 trips
        m[l] = a;
        n[l] = mc = sqrt(mc);
        c = (a + mc) / 2;
        if (!(fabs(a - mc) > tolJAC * a)) {
          ++l;
          break;
        }
        mc *= a;
        a = c;
      }
      x *= c;
      sn = sin(x);
      cn = cos(x);
      dn = 1;
      if (sn != 0) {
        real a = cn / sn;
        c *= a;
        while (l--) {
          real b = m[l];
          a *= c;
          c *= dn;
          dn = (n[l] + a) / (b + a);
          a = c / b;
        }
        a = 1 / sqrt(c*c + 1);
        sn = signbit(sn) ? -a : a;
        cn = c * sn;
        if (signbit(_kp2)) {
          // See DLMF Eqs 22.17.2 - 22.17.4
          swap(cn, dn);
          sn /= d;
        }
      }
    } else {
      sn = tanh(x);
      dn = cn = 1 / cosh(x);
    }
  }

  Math::real EllipticFunction::am(real x, real& sn, real& cn, real& dn) const {
    static const real tolJAC =
      pow(numeric_limits<real>::epsilon(), real(0.75));
    static const bool useSala = false;
    // Special cases of k2 = 0 and 1.
    if (_k2 == 0) {
      sn = sin(x); cn = cos(x); dn = 1;
      return x;
    } else if (_kp2 == 0) {
      sn = tanh(x); cn = dn = 1 / cosh(x);
      return atan(sinh(x));     // gd(x)
    }
    // Do argument reduction
    real y = remainder(x, 2 * K()), phi;
    long n = long(rint((x - y) / (2 * K())));
    // Now x = 2*n * K() + y where y in [-K(), K()].  K() is the quarter
    // period for the elliptic integral which corresponds to pi/2 in angle
    // space.  Thus am(x) = am(y) + n*pi with am(y) in [-pi/2, pi/2].
    if (y == 0) {
      sn = y; cn = 1; dn = 1;
      phi = y;
    } else if (fabs(y) == K()) {
      sn = copysign(real(1), y); cn = 0; dn = sqrt(_kp2);
      phi = copysign(Math::pi()/2, y);
    } else {
      if constexpr (useSala) {
        // This implements DLMF Sec 22.20(ii).
        // See also Sala (1989), https://doi.org/10.1137/0520100, Sec 5.
        //
        // But this method does extra calls to sin and asin and I'm doubtful
        // about its speed and accuracy.  So prefer Bulirsch's vgsncndn (the else
        // clause).
        real k2 = _k2, kp2 = _kp2;
        if (_k2 < 0) {
          // Sala Eq. 5.8
          k2 = -_k2 / _kp2; kp2 = 1 / _kp2;
          y *= sqrt(_kp2);
        }
        real a[num_], b, c[num_];
        a[0] = 1; b = sqrt(kp2); c[0] = sqrt(k2);
        int l = 1;
        for (; l < num_ ||
               GEOGRAPHICLIB_PANIC
               ("Convergence failure in EllipticFunction::am");) {
          a[l] = (a[l-1] + b) / 2;
          c[l] = (a[l-1] - b) / 2;
          b = sqrt(a[l-1] * b);
          if (!(c[l] > tolJAC * a[l])) break;
          ++l;
        }
        // Now a[l] = pi/(2*K)
        // Need to initialize phi1 to stop Visual Studio complaining
        phi = a[l] * y * real(1 << l);
        real phi1 = 0;
        for (; l > 0; --l) {
          phi1 = phi;
          phi = (phi + asin(c[l] * sin(phi) / a[l])) / 2;
        }
        if (_k2 < 0)
          // For k2 < 0, see Sala Eq. 5.8
          phi = phi1 - phi;
        sn = sin(phi); cn = cos(phi);
        // Since abs(phi) <= pi/2, cn cannot be negative
        if (signbit(cn)) cn = 0;
        dn = Delta(sn, cn);
      } else {
        sncndn(y, sn, cn, dn);
        if (signbit(cn)) cn = 0;
        phi = atan2(sn, cn);
      }
    }
    if (n % 2 != 0) {
      cn = -cn; sn = -sn;
    }
    return phi + n * Math::pi();
  }
  Math::real EllipticFunction::am(real x) const {
    real sn, cn, dn;
    return am(x, sn, cn, dn);
  }

#if GEOGRAPHICLIB_COMPLEX_JACOBI_AM
  Math::cmplx EllipticFunction::am(cmplx z, cmplx& sn, cmplx& cn, cmplx& dn)
    const {
    // From Lee (1976), modulus for u is k, v is k'
    // Also https://www.peliti.org/Notes/elliptic.pdf
    // sn(u+i*v) = (sn(u)*dn(v) + i*cn(u)*dn(u)*sn(v)*cn(v))/
    //             (1 - dn(u)^2*sn(v)^2)
    // cn(u+i*v) = (cn(u)*cn(v) - i*sn(u)*dn(u)*sn(v)*dn(v)/
    //             (1 - dn(u)^2*sn(v)^2)
    // dn(u+i*v) = (dn(u)*cn(v)*dn(v) - i*k2*sn(u)*cn(u)*sn(v))/
    //             (1 - dn(u)^2*sn(v)^2)
    // Denominator = cn(v)^2 + k2*sn(u)^2*sn(v)^2
    // From Sala (1995)
    // u
    // am(u+i*v) = atan(sn(u)*dn(v)/(cn(u)*cn(v))) +
    //             i*atanh(dn(u)*sn(v))
    // Do argument reduction

    if (!signbit(_k2)) {
      real x = z.real(), y = z.imag(),
        m = _k2, mp = _kp2,
        K = K_static(m, mp), Kp = K_static(mp, m),
        u = remainder(x, 2 * K), v = remainder(y, 2 * Kp);
      long s = long(rint((x - u) / (2 * K))),
        t = long(rint((y - v) / (2 * Kp)));
      // Now x = 2*s * K  + u where u in [-K , K ].
      //     y = 2*t * Kp + v where v in [-Kp, Kp]
      real snu, cnu, dnu, snv, cnv, dnv;
      sncndn_static(u, snu, cnu, dnu, m, mp); if (signbit(cnu)) cnu = 0;
      sncndn_static(v, snv, cnv, dnv, mp, m); if (signbit(cnv)) cnv = 0;
      if (s % 2 != 0) {
        cnu = -cnu; snu = -snu;
      }
      if (t % 2 != 0) {
        cnv = -cnv; snv = -snv;
      }
      // Lee: 1-dnu^2*snv^2 = cnv^2 + k2*snu^2*snv^2 = dnv^2 - k2*cnu^2*snv^2
      real den = Math::sq(cnv) + m * Math::sq(snu*snv); // N.B. m >= 0
      if (s % 2 != 0) {
        cnu = -cnu; snu = -snu;
      }
      sn = cmplx(snu*dnv, cnu*dnu*snv*cnv);
      cn = cmplx(cnu*cnv, -snu*dnu*snv*dnv);
      dn = cmplx(dnu*cnv*dnv, -_k2*snu*snu*snv);
      sn /= den; cn /= den; dn /= den;
      // Sala Eq. 4.13
      return {atan2(snu*dnv, cnu*cnv) + s * Math::pi(), atanh(dnu*snv)};
    } else {
      real m = -_k2/_kp2, mp = 1/_kp2,
        K  = K_static(m, mp), Kp = K_static(mp, m);
      z = z/sqrt(mp);
      real x = z.real(), y = z.imag();
      // Sala Sec 4.1 am(z, -m/m') = pi/2 - am(K - z/k', m)
      // sn(z, -m/m') = cn(K - z/k', m)
      // cn(K + z, m) = -k' * sd(z, m)
      // cn(K - z/k', m) = k' * sd(z/k', m)
      //    sn(z, -m/m') = k' * sd(z/k', m) (DLMF 22.17.6)
      // cn(z, -m/m') = sn(K - z/k, m)
      // sn(K + z, m) = cd(z, m)
      // sn(K - z/k, m) = cd(z/k, m)
      //    cn(z, -m/m') = cd(z/k, m) (DLMF 22.17.7)
      // sn^2 + cn^2 = m' * sd^2 + cd^2 = 1 (DLMF 22.6.4)
      // d(am(z, -m/m'))/dz = dn(z, -m/m') = 1/k'*dn(K - z/k', m)
      // dn(K + z, m) = k'*nd(z)
      // dn(K -z/k', m) = -k'nd(z/k')
      //    dn(z, -m/m') = nd(z/k')  (DLMF 22.17.8)
      // -m/m'*sn^2 + dn^2 = m * sd^2 + nd^2 = 1 (DLMF 22.6.4)
      // DLMF 22.27.5 - 22.17.8
        // -m/(1-m)  = _k2 => m = -k2/kp2; mp=1/kp2
      real
        u = remainder(x, 2 * K), v = remainder(y, 2 * Kp);
      long s = long(rint((x - u) / (2 * K))),
        t = long(rint((y - v) / (2 * Kp)));
      // Now x = 2*s * K  + u where u in [-K , K ].
      //     y = 2*t * Kp + v where v in [-Kp, Kp]
      real snu, cnu, dnu, snv, cnv, dnv;
      sncndn_static(u, snu, cnu, dnu, m, mp); if (signbit(cnu)) cnu = 0;
      sncndn_static(v, snv, cnv, dnv, mp, m); if (signbit(cnv)) cnv = 0;
      if (s % 2 != 0) {
        cnu = -cnu; snu = -snu;
      }
      if (t % 2 != 0) {
        cnv = -cnv; snv = -snv;
      }
      // Need sd, cd, nd from Lee p 113
      // k'*sd = sqrt(mp)*{snu*dnu*cnv, cnu*snv*dnv} / den
      // cd = {cnu*dnu*dnv, -mp*snu*snv*cnv} / den
      // nd = {dnu*cnv*dnv, m*snu*cnu*snv} / den
      // den = dnu^2 + dnv^2 - 1
      // phi = {atan2(sqrt(mp)*snu*snv, snu*dnv), xxx}
      // Lee: dnu^2 + dnv^2 - 1
      //    = m*cnu^2 + mp*cnv^2 = cnu^2*dnv^2+mp*snu^2*cnv^2
      real den = m * Math::sq(cnu) + mp * Math::sq(cnv);
      sn = cmplx(snu*dnu*cnv, cnu*snv*dnv); sn *= sqrt(mp);
      cn = cmplx(cnu*dnu*dnv, -mp*snu*snv*cnv);
      dn = cmplx(dnu*cnv*dnv, m*snu*cnu*snv);
      sn /= den; cn /= den; dn /= den;
      cmplx phi = atan(sn/cn) + s * Math::pi();
      return phi;
    }
  }
  Math::cmplx EllipticFunction::am(cmplx x) const {
    cmplx sn, cn, dn;
    return am(x, sn, cn, dn);
  }
#endif

  template<typename T>
  T EllipticFunction::Ft(T sn, T cn, T dn) const {
    // Carlson, eq. 4.5 and
    // https://dlmf.nist.gov/19.25.E5
    bool negs = signbit(Re(sn)), negc = signbit(Re(cn));
    T cn2 = cn*cn, dn2 = dn*dn,
      sna = negs ? -sn : sn,
      fi = cn2 != real(0) ? sna * RF(cn2, dn2, T(1)) : K();
    // Enforce usual trig-like symmetries
    if (negc) fi = 2 * K() - fi;
    if (negs) fi = -fi;
    return fi;
  }
  Math::real EllipticFunction::F(real sn, real cn, real dn) const {
    return Ft(sn, cn, dn);
  }
  Math::cmplx EllipticFunction::F(cmplx sn, cmplx cn, cmplx dn) const {
    return Ft(sn, cn, dn);
  }

  template<typename T>
  T EllipticFunction::Et(T sn, T cn, T dn) const {
    bool negs = signbit(Re(sn)), negc = signbit(Re(cn));
    T cn2 = cn*cn, dn2 = dn*dn, sn2 = sn*sn,
      cna = negc ? -cn : cn,
      sna = negs ? -sn : sn,
      ei = cn2 != real(0) ?
      sna * ( _k2 <= 0 ?
              // Carlson, eq. 4.6 and
              // https://dlmf.nist.gov/19.25.E9
              RF(cn2, dn2, T(1)) -
              _k2 * sn2 * RD(cn2, dn2, T(1)) / real(3) :
              ( _kp2 >= 0 ?
                // https://dlmf.nist.gov/19.25.E10
                _kp2 * RF(cn2, dn2, T(1)) +
                _k2 * _kp2 * sn2 * RD(cn2, T(1), dn2) / real(3) +
                _k2 * cna / dn :
                // https://dlmf.nist.gov/19.25.E11
                - _kp2 * sn2 * RD(dn2, T(1), cn2) / real(3) +
                dn / cna ) ) :
      E();
    // Enforce usual trig-like symmetries
    if (negc) ei = 2 * E() - ei;
    if (negs) ei = -ei;
    return ei;
  }
  Math::real EllipticFunction::E(real sn, real cn, real dn) const {
    return Et(sn, cn, dn);
  }
  Math::cmplx EllipticFunction::E(cmplx sn, cmplx cn, cmplx dn) const {
    return Et(sn, cn, dn);
  }

  Math::real EllipticFunction::D(real sn, real cn, real dn) const {
    // Carlson, eq. 4.8 and
    // https://dlmf.nist.gov/19.25.E13
    real
      cn2 = cn*cn, dn2 = dn*dn, sn2 = sn*sn,
      di = cn2 != 0 ? fabs(sn) * sn2 * RD(cn2, dn2, 1) / 3 : D();
    // Enforce usual trig-like symmetries
    if (signbit(cn))
      di = 2 * D() - di;
    return copysign(di, sn);
  }

  template<typename T>
  T EllipticFunction::Pit(T sn, T cn, T dn) const {
    // Carlson, eq. 4.7 and
    // https://dlmf.nist.gov/19.25.E14
    bool negs = signbit(Re(sn)), negc = signbit(Re(cn));
    T cn2 = cn*cn, dn2 = dn*dn, sn2 = sn*sn,
      sna = negs ? -sn : sn,
      pii = cn2 != real(0) ?
      sna * (RF(cn2, dn2, T(1)) +
             (_alpha2 * sn2 == real(0) ? T(0) :
              _alpha2 * sn2 *
              RJ(cn2, dn2, real(1), cn2 + _alphap2 * sn2) / real(3))) :
      Pi();
    // Enforce usual trig-like symmetries
    if (negc) pii = 2 * Pi() - pii;
    if (negs) pii = -pii;
    return pii;
  }
  Math::real EllipticFunction::Pi(real sn, real cn, real dn) const {
    return Pit(sn, cn, dn);
  }
  Math::cmplx EllipticFunction::Pi(cmplx sn, cmplx cn, cmplx dn) const {
    return Pit(sn, cn, dn);
  }

  Math::real EllipticFunction::G(real sn, real cn, real dn) const {
    real
      cn2 = cn*cn, dn2 = dn*dn, sn2 = sn*sn,
      gi = cn2 != 0 ? fabs(sn) * (RF(cn2, dn2, 1) +
                                  (_alpha2 - _k2) * sn2 *
                                  RJ(cn2, dn2, 1, cn2 + _alphap2 * sn2) / 3) :
      G();
    // Enforce usual trig-like symmetries
    if (signbit(cn))
      gi = 2 * G() - gi;
    return copysign(gi, sn);
  }

  Math::real EllipticFunction::H(real sn, real cn, real dn) const {
    real
      cn2 = cn*cn, dn2 = dn*dn, sn2 = sn*sn,
      // WARNING: large cancellation if k2 = 1, alpha2 = 0, and phi near pi/2
      hi = cn2 != 0 ? fabs(sn) * (RF(cn2, dn2, 1) -
                                  _alphap2 * sn2 *
                                  RJ(cn2, dn2, 1, cn2 + _alphap2 * sn2) / 3) :
      H();
    // Enforce usual trig-like symmetries
    if (signbit(cn))
      hi = 2 * H() - hi;
    return copysign(hi, sn);
  }

  Math::real EllipticFunction::deltaF(real sn, real cn, real dn) const {
    // Function is periodic with period pi
    if (signbit(cn)) { cn = -cn; sn = -sn; }
    return F(sn, cn, dn) * (Math::pi()/2) / K() - atan2(sn, cn);
  }

  Math::real EllipticFunction::deltaE(real sn, real cn, real dn) const {
    // Function is periodic with period pi
    if (signbit(cn)) { cn = -cn; sn = -sn; }
    return E(sn, cn, dn) * (Math::pi()/2) / E() - atan2(sn, cn);
  }

  Math::real EllipticFunction::deltaPi(real sn, real cn, real dn) const {
    // Function is periodic with period pi
    if (signbit(cn)) { cn = -cn; sn = -sn; }
    return Pi(sn, cn, dn) * (Math::pi()/2) / Pi() - atan2(sn, cn);
  }

  Math::real EllipticFunction::deltaD(real sn, real cn, real dn) const {
    // Function is periodic with period pi
    if (signbit(cn)) { cn = -cn; sn = -sn; }
    return D(sn, cn, dn) * (Math::pi()/2) / D() - atan2(sn, cn);
  }

  Math::real EllipticFunction::deltaG(real sn, real cn, real dn) const {
    // Function is periodic with period pi
    if (signbit(cn)) { cn = -cn; sn = -sn; }
    return G(sn, cn, dn) * (Math::pi()/2) / G() - atan2(sn, cn);
  }

  Math::real EllipticFunction::deltaH(real sn, real cn, real dn) const {
    // Function is periodic with period pi
    if (signbit(cn)) { cn = -cn; sn = -sn; }
    return H(sn, cn, dn) * (Math::pi()/2) / H() - atan2(sn, cn);
  }

  template<typename T>
  T EllipticFunction::Ft(T phi) const {
    if (_k2 == 0)
      return phi;
    else if (_kp2 == 0 && Im(phi) == 0 &&
             fabs(Re(phi)) <= Math::pi()/2 && !signbit(cos(Re(phi))))
      // Use this only for phi real and |phi| < pi/2
      return asinh(tan(phi));
    T sn = sin(phi), cn = cos(phi), Fv = Ft(sn, cn, Delta(sn, cn));
    real n;
    if constexpr (is_same_v<T, Math::real>)
      n = rint( (phi - atan2(sn, cn)) / (2 * Math::pi()) );
    else {
      real phir = phi.real();
      n = rint( (phir - atan2(sin(phir), cos(phir))) / (2 * Math::pi()) );
    }
    return n == 0 ? Fv : Fv + 4 * n * K();
  }
  Math::real EllipticFunction::F(real phi) const {
    return Ft(phi);
  }
  Math::cmplx EllipticFunction::F(cmplx phi) const {
    return Ft(phi);
  }

  template<typename T>
  T EllipticFunction::Et(T phi) const {
    if (_k2 == 0)
      return phi;
    // else if (_kp2 == 0)
    // Despite DLMF Eq 19.6.9 this is probably wrong, since
    // sqrt(1 - k^2*sin(phi)^2) -> abs(cos(phi)) in the limit k -> 1.
    //      return sin(phi);
    T sn = sin(phi), cn = cos(phi), Ev = Et(sn, cn, Delta(sn, cn));
    real n;
    if constexpr (is_same_v<T, Math::real>)
      n = rint( (phi - atan2(sn, cn)) / (2 * Math::pi()) );
    else {
      real phir = phi.real();
      n = rint( (phir - atan2(sin(phir), cos(phir))) / (2 * Math::pi()) );
    }
    return n == 0 ? Ev : Ev + 4 * n * E();
  }
  Math::real EllipticFunction::E(real phi) const {
    return Et(phi);
  }
  Math::cmplx EllipticFunction::E(cmplx phi) const {
    return Et(phi);
  }

  Math::real EllipticFunction::Ed(real ang) const {
    // ang - Math::AngNormalize(ang) is (nearly) an exact multiple of 360
    real n = rint((ang - Math::AngNormalize(ang))/Math::td);
    real sn, cn;
    Math::sincosd(ang, sn, cn);
    return E(sn, cn, Delta(sn, cn)) + 4 * E() * n;
  }

  template<typename T>
  T EllipticFunction::Pit(T phi) const {
    T sn = sin(phi), cn = cos(phi), Piv = Pit(sn, cn, Delta(sn, cn));
    real n;
    if constexpr (is_same_v<T, Math::real>)
      n = rint( (phi - atan2(sn, cn)) / (2 * Math::pi()) );
    else {
      real phir = phi.real();
      n = rint( (phir - atan2(sin(phir), cos(phir))) / (2 * Math::pi()) );
    }
    return n == 0 ? Piv : Piv + 4 * n * Pi();
  }
  Math::real EllipticFunction::Pi(real phi) const {
    return Pit(phi);
  }
  Math::cmplx EllipticFunction::Pi(cmplx phi) const {
    return Pit(phi);
  }

  Math::real EllipticFunction::D(real phi) const {
    real sn = sin(phi), cn = cos(phi), dn = Delta(sn, cn);
    return fabs(phi) < Math::pi() ? D(sn, cn, dn) :
      (deltaD(sn, cn, dn) + phi) * D() / (Math::pi()/2);
  }

  Math::real EllipticFunction::G(real phi) const {
    real sn = sin(phi), cn = cos(phi), dn = Delta(sn, cn);
    return fabs(phi) < Math::pi() ? G(sn, cn, dn) :
      (deltaG(sn, cn, dn) + phi) * G() / (Math::pi()/2);
  }

  Math::real EllipticFunction::H(real phi) const {
    real sn = sin(phi), cn = cos(phi), dn = Delta(sn, cn);
    return fabs(phi) < Math::pi() ? H(sn, cn, dn) :
      (deltaH(sn, cn, dn) + phi) * H() / (Math::pi()/2);
  }

  Math::real EllipticFunction::Einv(real x) const {
    static const real tolJAC =
      sqrt(numeric_limits<real>::epsilon() / real(100));
    if (k2() == 0 || x == 0) return x; // Preserve sign of 0
    real y = remainder(x, 2 * E()), n = rint((x - y) / (2 * E()));
    // Now x = 2 * n * E() + y
    // Linear approximation
    real phi = Math::pi() * y / (2 * E()); // phi in [-pi/2, pi/2)
    // First order correction
    phi -= _eps * sin(2 * phi) / 2;
    // For kp2 close to zero use asin(y/E()) or
    // J. P. Boyd, Applied Math. and Computation 218, 7005-7013 (2012)
    // https://doi.org/10.1016/j.amc.2011.12.021
    for (int i = 0;
         i < num_ ||
           GEOGRAPHICLIB_PANIC("Convergence failure in EllipticFunction::Einv");
         ++i) {
      real
        sn = sin(phi),
        cn = cos(phi),
        dn = Delta(sn, cn),
        err = (E(sn, cn, dn) - y)/dn;
      phi -= err;
      if (!(fabs(err) > tolJAC))
        break;
    }
    return n * Math::pi() + phi;
  }

  Math::real EllipticFunction::Piinv(real x) const {
    // exp(big) is close to max()
    static const real big = log(numeric_limits<real>::max()) - 1;
    // Use Math::tauf for ell.k2() == 1?
    static const bool usetauf = true;
    real y, n;
    if (kp2() == 0) {
      // for k^2 = 1, Pi() == inf
      if constexpr (usetauf) {
        // See tests Conformal3Proj[12] for the improvement this gives.
        //
        // y = Pi(phi; alpha2,1) = asinh(taup(tan(phi), alpha))/alphap2
        // inverse of Pi = Piinv, inverse of taup is tau
        // phi = Piinv(y; alpha2,1)
        // tan(phi) = tau(sinh(alphap2 * y), alpha)
        // This method preserves precision for x large, phi close to pi/2
        real es = copysign(sqrt(fabs(alpha2())), alpha2()),
          t = Math::tauf(sinh(alphap2() * x), es);
        return atan(t);
      } else {
        // if we don't want to rely on Math::tau, we can use the general method
        // for inverting Pi.
        y = x; n = 0;
      }
    } else if ((k2() == 0 && alpha2() == 0) ||
               x == 0)          // Preserve the sign of +/-0
      return x;
    else {
      y = remainder(x, 2 * Pi());
      n = 2 * rint((x - y) / (2 * Pi()));
    }
    // Now x = n * Pi() + y where y in [-Pi(), Pi()].  Pi() is the quarter
    // period for the elliptic integral which corresponds to pi/2 in angle
    // space.
    if (y == 0)
      return  n * Math::pi() / 2;
    else if (fabs(y) == Pi())
      return (copysign(real(1), y) + n) * Math::pi() / 2;
    else {
      // solve Pi(phi) = y for phi
      // Pi'(phi) = 1/(sqrt(1 - k2() * Math::sq(sin(phi)))
      //               * (1 - alpha2 * Math::sq(sin(phi))))
      // For k2 in [0,1]
      //    1 - k2() * Math::sq(sin(phi))
      //    = kp2() + k2() * Math::sq(cos(phi))
      // For alpha2 > 0
      //    1 - alpha2 * Math::sq(sin(phi))
      //    alphap2 + alpha2 * Math::sq(sin(phi))
      //
      // To preserve relative precision in sin(phi) and cos(phi) we let phi =
      // atan(exp(q)) and solved for q in [-inf, inf].
      // d/dq Pi(atan(exp(q))) = tan(phi)*cos(phi)^2 * Pi'
      auto Pif = [this]
        (real q) -> pair<real, real>
        {
          real t = exp(q), sc = hypot(real(1), t), s = t/sc, c = 1/sc,
          d = Delta(s, c),
          f = Pi(s, c, d),
          fp = t*c*c / (d * (alpha2() >= 0 ?
                             alphap2() + alpha2() * s*s :
                             1 - alpha2() * c*c));
          return {f, fp};
        };
      real z = Trigfun::root(Trigfun::PIINV,
                             Pif, fabs(y), 0,
                             -big, big,
                             1,1,1);
      return n * Math::pi() + atan(copysign(exp(z), y));
    }
  }

  Math::real EllipticFunction::deltaEinv(real stau, real ctau) const {
    // Function is periodic with period pi
    if (signbit(ctau)) { ctau = -ctau; stau = -stau; }
    real tau = atan2(stau, ctau);
    return Einv( tau * E() / (Math::pi()/2) ) - tau;
  }

  /// \cond SKIP
  // Instantiate
#define GEOGRAPHICLIB_ELLIPTIC_INSTANTIATE(T)                           \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::RFt(T, T, T);       \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::RFt(T, T);          \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::RCt(T, T);          \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::RGt(T, T, T);       \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::RGt(T, T);          \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::RJt(T, T, T, T);    \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::RDt(T, T, T);       \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::Ft(T) const;        \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::Ft(T, T, T) const;  \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::Et(T) const;        \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::Et(T, T, T) const;  \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::Pit(T) const;       \
  template T GEOGRAPHICLIB_EXPORT EllipticFunction::Pit(T, T, T) const;

  GEOGRAPHICLIB_ELLIPTIC_INSTANTIATE(Math::real)
  GEOGRAPHICLIB_ELLIPTIC_INSTANTIATE(Math::cmplx)

#undef GEOGRAPHICLIB_ELLIPTIC_INSTANTIATE
  /// \endcond

} // namespace GeographicLib

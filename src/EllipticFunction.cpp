/**
 * \file EllipticFunction.cpp
 * \brief Implementation for GeographicLib::EllipticFunction class
 *
 * Copyright (c) Charles Karney (2008-2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/EllipticFunction.hpp>
#include <algorithm>            // For swap

namespace GeographicLib {

  using namespace std;

  const Math::real EllipticFunction::tol_ =
    numeric_limits<real>::epsilon() * real(0.01);
  const Math::real EllipticFunction::tolRF_ = pow(3 * tol_, 1/real(8));
  const Math::real EllipticFunction::tolRD_ = pow(real(0.2) * tol_, 1/real(8));
  const Math::real EllipticFunction::tolRG0_ = real(2.7) * sqrt(tol_);
  const Math::real EllipticFunction::tolJAC_ = sqrt(tol_);

  /*
   * Implementation of methods given in
   *
   *   B. C. Carlson
   *   Computation of elliptic integrals
   *   Numerical Algorithms 10, 13-26 (1995)
   */

  Math::real EllipticFunction::RF(real x, real y, real z) throw() {
    // Carlson, eqs 2.2 - 2.7
    real
      A0 = (x + y + z)/3,
      An = A0,
      Q = max(max(abs(A0-x), abs(A0-y)), abs(A0-z)) / tolRF_,
      x0 = x,
      y0 = y,
      z0 = z,
      mul = 1;
    while (Q >= mul * abs(An)) {
      // Max 6 trips
      real lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0);
      An = (An + lam)/4;
      x0 = (x0 + lam)/4;
      y0 = (y0 + lam)/4;
      z0 = (z0 + lam)/4;
      mul *= 4;
    }
    real
      X = (A0 - x) / (mul * An),
      Y = (A0 - y) / (mul * An),
      Z = - (X + Y),
      E2 = X*Y - Z*Z,
      E3 = X*Y*Z;
    // http://dlmf.nist.gov/19.36.E1
    // Polynomial is
    // (1 - E2/10 + E3/14 + E2^2/24 - 3*E2*E3/44
    //    - 5*E2^3/208 + 3*E3^2/104 + E2^2*E3/16)
    // convert to Horner form...
    return (E3 * (6930 * E3 + E2 * (15015 * E2 - 16380) + 17160) +
            E2 * ((10010 - 5775 * E2) * E2 - 24024) + 240240) /
      (240240 * sqrt(An));
  }


  Math::real EllipticFunction::RF(real x, real y) throw() {
    // Carlson, eqs 2.36 - 2.38
    real xn = sqrt(x), yn = sqrt(y);
    if (xn < yn) swap(xn, yn);
    while (abs(xn-yn) > tolRG0_ * xn) {
      // Max 4 trips
      real t = (xn + yn) /2;
      yn = sqrt(xn * yn);
      xn = t;
    }
    return Math::pi<real>() / (xn + yn);
  }

  Math::real EllipticFunction::RC(real x, real y) throw() {
    return ( x < y ?
             // http://dlmf.nist.gov/19.2.E18
             atan(sqrt((y - x) / x)) / sqrt(y - x) :
             Math::atanh( y > 0 ?
                          // http://dlmf.nist.gov/19.2.E19
                          sqrt((x - y) / x) :
                          // http://dlmf.nist.gov/19.2.E20
                          sqrt(x / (x - y)) ) / sqrt(x - y) );
  }

  Math::real EllipticFunction::RG(real x, real y, real z) throw() {
    if (z == 0)
      swap(y, z);
    // Carlson, eq 1.7
    return (z * RF(x, y, z) - (x-z) * (y-z) * RD(x, y, z) / 3
            + sqrt(x * y / z)) / 2;
  }

  Math::real EllipticFunction::RG(real x, real y) throw() {
    // Carlson, eqs 2.36 - 2.39
    if (x < y) swap(x, y);
    real
      x0 = sqrt(x),
      y0 = sqrt(y),
      xn = x0,
      yn = y0,
      s = 0,
      mul = real(0.25);
    while (abs(xn-yn) > tolRG0_ * xn) {
      // Max 4 trips
      real t = (xn + yn) /2;
      yn = sqrt(xn * yn);
      xn = t;
      mul *= 2;
      t = xn - yn;
      s += mul * t * t;
    }
    return (Math::sq( (x0 + y0)/2 ) - s) * Math::pi<real>() / (2 * (xn + yn));
  }

  Math::real EllipticFunction::RJ(real x, real y, real z, real p) throw() {
    // Carlson, eqs 2.17 - 2.25
    real
      A0 = (x + y + z + 2*p)/5,
      An = A0,
      delta = (p-x) * (p-y) * (p-z),
      Q = max(max(abs(A0-x), abs(A0-y)), max(abs(A0-z), abs(A0-p))) / tolRD_,
      x0 = x,
      y0 = y,
      z0 = z,
      p0 = p,
      mul = 1,
      mul3 = 1,
      s = 0;
    while (Q >= mul * abs(An)) {
      // Max 7 trips
      real
        lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0),
        d0 = (sqrt(p0)+sqrt(x0)) * (sqrt(p0)+sqrt(y0)) * (sqrt(p0)+sqrt(z0)),
        e0 = delta/(mul3 * Math::sq(d0));
      s += RC(1, 1 + e0)/(mul * d0);
      An = (An + lam)/4;
      x0 = (x0 + lam)/4;
      y0 = (y0 + lam)/4;
      z0 = (z0 + lam)/4;
      p0 = (p0 + lam)/4;
      mul *= 4;
      mul3 *= 64;
    }
    real
      X = (A0 - x) / (mul * An),
      Y = (A0 - y) / (mul * An),
      Z = (A0 - z) / (mul * An),
      P = -(X + Y + Z) / 2,
      E2 = X*Y + X*Z + Y*Z - 3*P*P,
      E3 = X*Y*Z + 2*P * (E2 + 2*P*P),
      E4 = (2*X*Y*Z + P * (E2 + 3*P*P)) * P,
      E5 = X*Y*Z*P*P;
    // http://dlmf.nist.gov/19.36.E2
    // Polynomial is
    // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
    //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
    //    - 9*(E3*E4+E2*E5)/68)
    return ((471240 - 540540 * E2) * E5 +
            (612612 * E2 - 540540 * E3 - 556920) * E4 +
            E3 * (306306 * E3 + E2 * (675675 * E2 - 706860) + 680680) +
            E2 * ((417690 - 255255 * E2) * E2 - 875160) + 4084080) /
      (4084080 * mul * An * sqrt(An)) + 6 * s;
  }

  Math::real EllipticFunction::RD(real x, real y, real z) throw() {
    // Carlson, eqs 2.28 - 2.34
    real
      A0 = (x + y + 3*z)/5,
      An = A0,
      Q = max(max(abs(A0-x), abs(A0-y)), abs(A0-z)) / tolRD_,
      x0 = x,
      y0 = y,
      z0 = z,
      mul = 1,
      s = 0;
    while (Q >= mul * abs(An)) {
      // Max 7 trips
      real lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0);
      s += 1/(mul * sqrt(z0) * (z0 + lam));
      An = (An + lam)/4;
      x0 = (x0 + lam)/4;
      y0 = (y0 + lam)/4;
      z0 = (z0 + lam)/4;
      mul *= 4;
    }
    real
      X = (A0 - x) / (mul * An),
      Y = (A0 - y) / (mul * An),
      Z = -(X + Y) / 3,
      E2 = X*Y - 6*Z*Z,
      E3 = (3*X*Y - 8*Z*Z)*Z,
      E4 = 3 * (X*Y - Z*Z) * Z*Z,
      E5 = X*Y*Z*Z*Z;
    // http://dlmf.nist.gov/19.36.E2
    // Polynomial is
    // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
    //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
    //    - 9*(E3*E4+E2*E5)/68)
    return ((471240 - 540540 * E2) * E5 +
            (612612 * E2 - 540540 * E3 - 556920) * E4 +
            E3 * (306306 * E3 + E2 * (675675 * E2 - 706860) + 680680) +
            E2 * ((417690 - 255255 * E2) * E2 - 875160) + 4084080) /
      (4084080 * mul * An * sqrt(An)) + 3 * s;
  }

  EllipticFunction::EllipticFunction(real m) throw()
    : _m(m)
    , _m1(1 - m)
      // Don't initialize _kc, _ec, _kec since this constructor might be called
      // before the static real constants tolRF_, etc., are initialized.
    , _init(false)
  {}

  bool EllipticFunction::Init() const throw() {
    // Complete elliptic integral K(m), Carlson eq. 4.1
    _kc = _m1 ? RF(_m1, 1) : Math::infinity<real>();
    // Complete elliptic integral E(m), Carlson eq. 4.2
    _ec = _m1 ? 2 * RG(_m1, 1) : 1;
    // D(m) = (K(m) - E(m))/m, Carlson eq.4.3
    _dc = _m1 ? RD(real(0), _m1, 1) / 3 : Math::infinity<real>();
    return _init = true;
  }

  /*
   * Implementation of methods given in
   *
   *   R. Bulirsch
   *   Numerical Calculation of Elliptic Integrals and Elliptic Functions
   *   Numericshe Mathematik 7, 78-90 (1965)
   */

  void EllipticFunction::sncndn(real x, real& sn, real& cn, real& dn)
    const throw() {
    // Bulirsch's sncndn routine, p 89.
    if (_m1 != 0) {
      real mc = _m1, d = 0;
      if (_m1 < 0) {
        d = 1 - mc;
        mc /= -d;
        d = sqrt(d);
        x *= d;
      }
      real c = 0;           // To suppress warning about uninitialized variable
      real m[num_], n[num_];
      unsigned l = 0;
      for (real a = 1; l < num_; ++l) {
        // Max 5 trips
        m[l] = a;
        n[l] = mc = sqrt(mc);
        c = (a + mc) / 2;
        if (!(abs(a - mc) > tolJAC_ * a)) {
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
        a = 1 / sqrt(c * c + 1);
        sn = sn < 0 ? -a : a;
        cn = c * sn;
        if (_m1 < 0) {
          swap(cn, dn);
          sn /= d;
        }
      }
    } else {
      sn = tanh(x);
      dn = cn = 1 / cosh(x);
    }
  }

  Math::real EllipticFunction::E(real sn, real cn, real dn) const throw() {
    real
      cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn,
      ei = ( _m <= 0 ?
             // Carlson, eq. 4.6
             // http://dlmf.nist.gov/19.25.E9
             RF(cn2, dn2, 1) - (_m / 3) * sn2 * RD(cn2, dn2, 1) :
             ( _m <= 1 ?
               // http://dlmf.nist.gov/19.25.E10
               _m1 * RF(cn2, dn2, 1) +
               (_m * _m1 / 3) * sn2 * RD(cn2, 1, dn2) +
               _m * abs(cn) / dn :
               // http://dlmf.nist.gov/19.25.E11
               - (_m1 / 3) * sn2 * RD(dn2, 1, cn2) + dn / abs(cn) ) );
    ei *= abs(sn);
    // Enforce usual trig-like symmetries
    if (cn < 0)
      ei = 2 * E() - ei;
    if (sn < 0)
      ei = -ei;
    return ei;
  }

  Math::real EllipticFunction::E(real phi) const throw() {
    real sn, cn, n = Math::sincosp(phi, sn, cn);
    return 4 * n * E() +
      E(sn, cn, sqrt(_m < 0 ? 1 - _m * sn * sn : _m1 + _m * cn * cn));
  }

  Math::real EllipticFunction::Ed(real ang) const throw() {
    real sn, cn, n = Math::sincospd(ang, sn, cn);
    return 4 * n * E() +
      E(sn, cn, sqrt(_m < 0 ? 1 - _m * sn * sn : _m1 + _m * cn * cn));
  }

} // namespace GeographicLib

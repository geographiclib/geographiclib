/**
 * \file EllipticFunction.cpp
 * \brief Implementation for GeographicLib::EllipticFunction class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/EllipticFunction.hpp"
#include <limits>
#include <algorithm>

#define GEOGRAPHICLIB_ELLIPTICFUNCTION_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_ELLIPTICFUNCTION_CPP)
RCSID_DECL(GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real_t EllipticFunction::tol =
    numeric_limits<real_t>::epsilon() * real_t(0.01L);
  const Math::real_t EllipticFunction::tolRF = pow(3 * tol, 1/real_t(6));
  const Math::real_t EllipticFunction::tolRD =
    pow(real_t(0.25L) * tol, 1/real_t(6));
  const Math::real_t EllipticFunction::tolRG0 = real_t(2.7L) * sqrt(tol);
  const Math::real_t EllipticFunction::tolJAC = sqrt(tol);
  const Math::real_t EllipticFunction::tolJAC1 = sqrt(6 * tol);

  /*
   * Implementation of methods given in
   *
   *   B. C. Carlson
   *   Computation of elliptic integrals
   *   Numerical Algorithms 10, 13-26 (1995)
   */

  Math::real_t EllipticFunction::RF(real_t x, real_t y, real_t z) throw() {
    // Carlson, eqs 2.2 - 2.7
    real_t
      a0 = (x + y + z)/3,
      an = a0,
      q = max(max(abs(a0-x), abs(a0-y)), abs(a0-z)) / tolRF,
      x0 = x,
      y0 = y,
      z0 = z,
      mul = 1;
    while (q >= mul * abs(an)) {
      // Max 6 trips
      real_t ln = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0);
      an = (an + ln)/4;
      x0 = (x0 + ln)/4;
      y0 = (y0 + ln)/4;
      z0 = (z0 + ln)/4;
      mul *= 4;
    }
    real_t
      xx = (a0 - x) / (mul * an),
      yy = (a0 - y) / (mul * an),
      zz = - xx - yy,
      e2 = xx * yy - zz * zz,
      e3 = xx * yy * zz;
    return (1 - e2 / 10 + e3 / 14 + e2 * e2 / 24 - 3 * e2 * e3 / 44) / sqrt(an);
  }

  Math::real_t EllipticFunction::RD(real_t x, real_t y, real_t z) throw() {
    // Carlson, eqs 2.28 - 2.34
    real_t
      a0 = (x + y + 3 * z)/5,
      an = a0,
      q = max(max(abs(a0-x), abs(a0-y)), abs(a0-z)) / tolRD,
      x0 = x,
      y0 = y,
      z0 = z,
      mul = 1,
      s = 0;
    while (q >= mul * abs(an)) {
      // Max 7 trips
      real_t ln = sqrt(x0)*sqrt(y0) +
        sqrt(y0)*sqrt(z0) +
        sqrt(z0)*sqrt(x0);
      s += 1/(mul * sqrt(z0) * (z0 + ln ));
      an = (an + ln)/4;
      x0 = (x0 + ln)/4;
      y0 = (y0 + ln)/4;
      z0 = (z0 + ln)/4;
      mul *= 4;
    }
    real_t
      xx = (a0 - x) / (mul * an),
      yy = (a0 - y) / (mul * an),
      zz = -(xx + yy) / 3,
      e2 = xx * yy - 6 * zz * zz,
      e3 = (3 * xx * yy - 8 * zz * zz)*zz,
      e4 = 3 * (xx * yy - zz * zz) * zz * zz,
      e5 = xx * yy * zz * zz * zz;
    return (1 - 3 * e2 / 14 + e3 / 6 + 9 * e2 * e2 / 88 - 3 * e4 / 22
            - 9 * e2 * e3 / 52 + 3 * e5 / 26) / (mul * an * sqrt(an))
      + 3 * s;
  }

  Math::real_t EllipticFunction::RG0(real_t x, real_t y) throw() {
    // Carlson, eqs 2.36 - 2.39
    real_t
      x0 = sqrt(x),
      y0 = sqrt(y),
      xn = x0,
      yn = y0,
      s = 0,
      mul = real_t(0.25L);
    while (abs(xn-yn) >= tolRG0 * abs(xn)) {
      // Max 4 trips
      real_t t = (xn + yn) /2;
      yn = sqrt(xn * yn);
      xn = t;
      mul *= 2;
      t = xn - yn;
      s += mul * t * t;
    }
    x0 = (x0 + y0)/2;
    return  (x0 * x0 - s) * Constants::pi() / (2 * (xn + yn));
  }

  EllipticFunction::EllipticFunction(real_t m) throw()
    : _m(m)
    , _m1(1 - m)
      // Don't initialize _kc, _ec, _kec since this constructor might be called
      // before the static real_t constants tolRF, etc., are initialized.
    , _init(false)
  {}

  bool EllipticFunction::Init() const throw() {
    // Complete elliptic integral K(m), Carlson eq. 4.1
    _kc = RF(real_t(0), _m1, real_t(1));
    // Complete elliptic integral E(m), Carlson eq. 4.2
    _ec = 2 * RG0(_m1, real_t(1));
    // K - E, Carlson eq.4.3
    _kec = _m / 3 * RD(real_t(0), _m1, real_t(1));
    return _init = true;
  }

  /*
   * Implementation of methods given in
   *
   *   R. Bulirsch
   *   Numerical Calculation of Elliptic Integrals and Elliptic Functions
   *   Numericshe Mathematik 7, 78-90 (1965)
   */

  void EllipticFunction::sncndn(real_t x,
                                real_t& sn, real_t& cn, real_t& dn)
    const throw() {
    // Bulirsch's sncndn routine, p 89.
    //
    // Assume _m1 is in [0, 1].  See Bulirsch article for code to treat
    // negative _m1.
    if (_m1 != 0) {
      real_t mc = _m1;
      real_t c;
      real_t m[num], n[num];
      unsigned l = 0;
      for (real_t a = 1; l < num; ++l) {
        // Max 5 trips
        m[l] = a;
        n[l] = mc = sqrt(mc);
        c = (a + mc) / 2;
        if (abs(a - mc) <= tolJAC * a) {
          ++l;
          break;
        }
        mc = a * mc;
        a = c;
      }
      x = c * x;
      sn = sin(x);
      cn = cos(x);
      dn = 1;
      if (sn != 0) {
        real_t a = cn / sn;
        c = a * c;
        while (l--) {
          real_t b = m[l];
          a = c * a;
          c = dn * c;
          dn = (n[l] + a) / (b + a);
          a = c / b;
        }
        a = 1 / sqrt(c * c + 1);
        sn = sn < 0 ? -a : a;
        cn = c * sn;
      }
    } else {
      sn = tanh(x);
      dn = cn = 1 / cosh(x);
    }
  }

  Math::real_t EllipticFunction::E(real_t sn, real_t cn, real_t dn)
    const throw() {
    real_t ei;
    cn *= cn; dn *= dn;
    // Carlson, eq. 4.6
    ei = abs(sn) * (RF(cn, dn, real_t(1)) -
                    (_m/3) * sn*sn * RD(cn, dn, real_t(1)));
    // Enforce usual trig-like symmetries
    if (cn < 0) {
      ei = 2 * E() - ei;
    }
    if (sn < 0)
      ei = -ei;
    return ei;
  }

} // namespace GeographicLib

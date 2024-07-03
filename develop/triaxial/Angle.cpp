/**
 * \file Angle.cpp
 * \brief Implementation for the GeographicLib::Angle class.
 *
 * This file is an implementation of the methods described in
 * - C. F. F. Karney,
 *   <a href="https://doi.org/10.1080/00396265.2023.2217604">
 *   On auxiliary latitudes,</a>
 *   Survey Review 56(395), 165--180 (2024);
 *   preprint
 *   <a href="https://arxiv.org/abs/2212.05818">arXiv:2212.05818</a>.
 * .
 * Copyright (c) Charles Karney (2022-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "Angle.hpp"
#include <iostream>

namespace GeographicLib {

  using namespace std;

  Angle::Angle(real s, real c, real num, bool normp)
    : _s(s)
    , _c(c)
    , _n(num)
  {
    if (!normp) {
      real h = hypot(_s, _c);
      if (h == 0) {
        _s = 0; _c = 1;
      } else if (isfinite(h)) {
        _s /= h; _c /= h;
      } else if (isnan(h) || (isinf(_s) && isinf(_c)))
        _s = _c = Math::NaN();
      else if (isinf(_s)) {
        // infinite, finite
        _s = copysign(real(1), _s);
        _c = copysign(real(0), _c);
      } else {
        // isinf(cos); finite, infinite
        _s = copysign(real(0), _s);
        _c = copysign(real(1), _c);
      }
    }
  }

  Angle::Angle(Math::real deg) {
    Math::sincosd(deg, _s, _c);
    _n = round( (deg - Math::atan2d(_s, _c)) / Math::td );
  }

  Angle::operator Math::real() const {
    return Math::td * _n + Math::atan2d(_s, _c);
  }

  Angle Angle::radians(Math::real rad) {
    real sn = sin(rad), cs = cos(rad);
    return Angle(sn, cs, round( (rad - atan2(sn, cs)) / (2 * Math::pi()) ),
                 true);
  }

  Math::real Angle::radians() const {
    return 2 * Math::pi() * _n + atan2(_s, _c);
  }

  Math::real Angle::radians0() const { return atan2(_s, _c); }

  Angle Angle::lam(Math::real psi) {
    return Angle(sinh(psi), 1, 0);
  }

  Math::real Angle::lam() const {
    return asinh(t());
  }

  Angle Angle::NaN() {
    return Angle(Math::NaN(), Math::NaN(), Math::NaN(), true);
  }

  Angle Angle::cardinal(real q) {
    if (!isfinite(q)) return Angle::NaN();
    q = round(q);
    int iq = int(remainder(q, real(4)));
    // iq is in [-2, 2];
    // We could fold iq = -2 to iq = 2; but this way work too.
    real s, c, z = 0;
    switch (iq) {
    case -2: s = -z; c = -1; break;
    case -1: s = -1; c =  z; break;
    case  1: s =  1; c =  z; break;
    case  2: s =  z; c = -1; break;
    default: s =  z; c =  1; break; // iq = 0
    }
    return Angle(s, c, round((q - iq) / 4));
  }

  Angle Angle::cardinaldir(unsigned ind) const {
    real s, c;
    if (ind == 0U) {
      if (fabs(_c) >= fabs(_s)) {
        s = copysign(real(0), _s); c = copysign(real(1), _c);
      } else {
        s = copysign(real(1), _s); c = copysign(real(0), _c);
      }
    } else if ((ind & 1U) == 0U) { // ind nonzero and even
      s = copysign(real(0), _s); c = copysign(real(1), _c);
    } else {                    // ind odd
      s = copysign(real(1), _s); c = copysign(real(0), _c);
    }
    return Angle(s, c, _n, true);
  }

  Math::real Angle::ncardinal(unsigned ind) const {
    int iq;
    if (ind == 0U)
      iq = (signbit(_s) ? -1 : 1) * (signbit(_c) ?
                                     ( -_c >= fabs(_s) ? 2 : 1 ) :
                                     (  _c >= fabs(_s) ? 0 : 1 ));
    else if ((ind & 1U) == 0U)  // ind nonzero and even
      iq = signbit(_c) ? (signbit(_s) ? -2 : 2) : 0;
    else                        // ind odd
      iq = signbit(_s) ? -1 : 1;
    return 4 * _n + iq;
  }

  Angle Angle::eps() {
    return Angle(numeric_limits<real>::epsilon() / (1 << 20), 1, 0, true);
  }

  // Angle Angle::operator+() const { return *this; }
  Angle Angle::operator-() const {
    return Angle(-_s, _c, -_n, true);
  }

  Angle& Angle::operator+=(const Angle& p) {
    real q = ncardinal() + p.ncardinal();
    real c = _c * p._c - _s * p._s;
    _s = _s * p._c + _c * p._s;
    _c = c;
    _n += p._n;
    q -= ncardinal();
    _n += round(q / 4);
    return *this;
  }

  Angle Angle::operator+(const Angle& p) const {
    Angle t = *this; t += p;
    return t;
  }

  Angle& Angle::operator-=(const Angle& p) {
    *this += -p;
    return *this;
  }

  Angle Angle::operator-(const Angle& p) const {
    Angle t = *this; t -= p;
    return t;
  }

  bool Angle::zerop(real mult) const {
    return _c > 0 && fabs(_s) <= mult * numeric_limits<real>::epsilon();
  }

  bool Angle::operator==(const Angle& p) const {
    Angle t = *this; t -= p;
    return t.zerop();
  }

  Math::real Angle::rnd(real x) {
    // This value of z more-or-less matches the value z = 1/16 in
    // Math::AngRound (where the argument is in degrees).
    static const real z = 1/real(1024);
    GEOGRAPHICLIB_VOLATILE real y = fabs(x);
    GEOGRAPHICLIB_VOLATILE real w = z - y;
    // The compiler mustn't "simplify" z - (z - y) to y
    y = w > 0 ? z - w : y;
    return copysign(y, x);
  }

  Angle& Angle::rnd() {
    _s = rnd(_s); _c = rnd(_c);
    return *this;
  }

  Angle Angle::rnded() const {
    Angle t = *this;
    return t.rnd();
  }

  Angle Angle::base() const {
    return Angle(_s, _c, 0, true);
  }

  Angle Angle::rebase(const Angle& c) const {
    return (*this - c).base() + c;
  }

  Angle& Angle::renormalize() {
    real h = hypot(_s, _c); _s /= h; _c /= h;
    return *this;
  }

  Angle& Angle::setn(Math::real n) {
    _n = round(n);
    return *this;
  }

  Angle& Angle::setquadrant(unsigned q) {
    _s = copysign(_s, real(             q  & 2U ? -1 : 1 ));
    _c = copysign(_c, real( ((q >> 1) ^ q) & 1U ? -1 : 1 ));
    return *this;
  }

  unsigned Angle::quadrant() const {
    return 2U * signbit(_s) + (signbit(_c) ^ signbit(_s));
  }

  Angle& Angle::reflect(bool flips, bool flipc, bool swapp) {
    if (flips) _s *= -1;
    if (flipc) _c *= -1;
    if (swapp) swap(_s, _c);
    return *this;
  }

  Angle Angle::flipsign(int mult) const {
    return mult < 0 ? -*this : *this;
  }

  /*
  Angle Angle::aux(int id, Math::real s, Math::real c, bool donorm) {
    using std::hypot; using std::fabs;
    // donorm = normalization required
    // normq = already normalized
    bool normq = fabs(hypot(s, c) - 1) < 1/real(10000);
    if (donorm && normq)
      cerr << id
           << " angauxOKISH normalization requested but already normalized\n";
    else if (!donorm && !normq)
      cerr << id
           << " angauxBAD normalization not requested but needed\n";
    else if (donorm)
      cerr << id
           << " angauxOKA normalization requested and needed\n";
    else
      cerr << id
           << " angauxOKB normalization not requested nor needed\n";
    return Angle(s, c, 0, !donorm);
  }
  */
} // namespace GeographicLib

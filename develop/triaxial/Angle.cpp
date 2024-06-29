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
      if (isfinite(h)) {
        // h == 0 give _s  = _c = NaN
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
    iq = iq == -2 ? 2 : iq;     // Now iq in (-2, 2];
    real s, c;
    switch (iq) {
    case -1: s = -1; c =  0; break;
    case  0: s =  0; c =  1; break;
    case  1: s =  1; c =  0; break;
    default: s =  0; c = -1; break; // case 2
    }
    return Angle(s, c, round((q - iq) / 4));
  }

  Angle Angle::cardinal() const {
    real s, c;
    if (fabs(_c) >= fabs(_s)) {
      s = copysign(real(0), _s); c = copysign(real(1), _c); 
    } else {
      s = copysign(real(1), _s); c = copysign(real(0), _c);
    }
    return Angle(s, c, _n, true);
  }

  Angle Angle::eps() {
    return Angle(numeric_limits<real>::epsilon() / (1 << 20), 1, 0, true);
  }

  Math::real Angle::quadrant() const {
    int iq = (signbit(_s) ? -1 : 1) * (signbit(_c) ?
                                       ( -_c >= fabs(_s) ? 2 : 1 ) :
                                       (  _c >= fabs(_s) ? 0 : 1 ));
    return 4 * _n + iq;
  }

  // Angle Angle::operator+() const { return *this; }
  Angle Angle::operator-() const {
    return Angle(-_s, _c, -_n, true);
  }

  Angle& Angle::operator+=(const Angle& p) {
    real q = quadrant() + p.quadrant();
    real c = _c * p._c - _s * p._s;
    _s = _s * p._c + _c * p._s;
    _c = c;
    _n += p._n;
    q -= quadrant();
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

  Angle Angle::rounded() const {
    return Angle(rnd(_s), rnd(_c), _n, true);
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

  Angle& Angle::setsigns(unsigned q) {
    _s = copysign(_s, real(             q  & 2U ? -1 : 1 ));
    _c = copysign(_c, real( ((q >> 1) ^ q) & 1U ? -1 : 1 ));
    return *this;
  }
} // namespace GeographicLib

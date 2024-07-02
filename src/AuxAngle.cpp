/**
 * \file AuxAngle.cpp
 * \brief Implementation for the GeographicLib::AuxAngle class.
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

#include <GeographicLib/AuxAngle.hpp>

namespace GeographicLib {

  using namespace std;

  AuxAngle AuxAngle::NaN() {
    return AuxAngle(Math::NaN(), Math::NaN(), false);
  }

  AuxAngle AuxAngle::cardinal(unsigned n) {
    n &= 3U;
    AuxAngle t;                 // Initialized to (0, 1)
    switch (n) {
    case 1U: t._y =  1; t._x =  0; break;
    case 3U: t._y = -1; t._x =  0; break;
    case 2U:            t._x = -1; break;
    default:                       break; // q = 0U, return (0,1)
    }
    return t;
  }

  AuxAngle& AuxAngle::scale() {
    using std::isnan; using std::isinf; // Needed for Centos 7, ubuntu 14
    if (isnan(_x) || isnan(_y) ||       // either component is nan
        (isinf(_x) && isinf(_y)) ||     // both components are infinite
        (_x == 0 && _y == 0))           // both components are zero
      _x = _y = Math::NaN();
    else if (isinf(_x)) {
      _y = 0 * _y; _x = copysign(real(1), _x);
    } else if (isinf(_y)) {
      _x = 0 * _x; _y = copysign(real(1), _y);
    } else {
      int e;
      (void) frexp(fmax(fabs(_x), fabs(_y)), &e);
      if (e > 0) --e;
      if (e) {
        _x = ldexp(_x, -e);
        _y = ldexp(_y, -e);
      }
    }
    return *this;
  }

  AuxAngle& AuxAngle::normalize() {
    real r = hypot(_y, _x); _y /= r; _x /= r;
    return *this;
  }

  AuxAngle AuxAngle::normalized() const {
    AuxAngle t = *this;
    return t.normalize();
  }

  AuxAngle& AuxAngle::complement() {
    swap(_x, _y);
    return *this;
  }

  AuxAngle& AuxAngle::operator+=(const AuxAngle& p) {
    real x = _x * p._x - _y * p._y;
    _y = _y * p._x + _x * p._y;
    _x = x;
    return *this;
  }

  AuxAngle AuxAngle::operator+(const AuxAngle& p) const {
    AuxAngle t = *this; t += p;
    return t;
  }

  AuxAngle& AuxAngle::operator-=(const AuxAngle& p) {
    real x = _x * p._x + _y * p._y;
    _y = _y * p._x - _x * p._y;
    _x = x;
    return *this;
  }

  AuxAngle AuxAngle::operator-(const AuxAngle& p) const {
    AuxAngle t = *this; t -= p;
    return t;
  }

  bool AuxAngle::zerop(real eps) const {
    return _x > 0 && fabs(_y) <= eps * _x;
  }

  bool AuxAngle::operator==(const AuxAngle& p) const {
    AuxAngle t = *this; t -= p;
    return t.zerop();
  }

  unsigned AuxAngle::quadrant() const {
    return 2U * signbit(_y) + (signbit(_x) ^ signbit(_y));
  }

  AuxAngle& AuxAngle::setquadrant(unsigned q) {
    _y = copysign(_y, q & 2U ? -real(1) : real(1));
    _x = copysign(_x, ((q >> 1) ^ q) & 1U ? -real(1) : real(1));
    return *this;
  }

  AuxAngle AuxAngle::copyquadrant(const AuxAngle& p) const {
    return AuxAngle(copysign(y(), p.y()), copysign(x(), p.x()), false);
  }

  Math::real AuxAngle::rnd(real x) {
    // This value of z more-or-less matches the value z = 1/16 in
    // Math::AngRound (where the argument is in degrees).
    static const real z = 1/real(1024);
    GEOGRAPHICLIB_VOLATILE real y = fabs(x);
    GEOGRAPHICLIB_VOLATILE real w = z - y;
    // The compiler mustn't "simplify" z - (z - y) to y
    y = w > 0 ? z - w : y;
    return copysign(y, x);
  }

  AuxAngle& AuxAngle::round() {
    _y = rnd(_y); _x = rnd(_x);
    return *this;
  }

  AuxAngle AuxAngle::rounded() const {
    AuxAngle t = *this;
    return t.round();
  }

} // namespace GeographicLib

/**
 * \file Angle.hpp
 * \brief Header for the GeographicLib::Angle class
 *
 * The class provide an accurate representation of angle via 3 numbers, its
 * sine and cosine, and the number of turns.
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ANGLE_HPP)
#define GEOGRAPHICLIB_ANGLE_HPP 1

#include <GeographicLib/Math.hpp>

namespace GeographicLib {

  /**
   * \brief An accurate representation of angles.
   *
   * The class provide an accurate representation of angle via 3 numbers, its
   * sine and cosine, and the number of turns.  The angle is then
   *   2*pi*n + atan2(sin, cos)
   *
   * N.B. n is stored as a real.  This allows it to be inf or nan.
   *
   * Example of use:
   * xx include example-Angle.cpp
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT Angle {
  private:
    typedef Math::real real;
    real _s, _c, _n;
    static real rnd(real x);
  public:
    real s() const { return _s; }
    real c() const { return _c; }
    real t() const { return _s/_c; }
    real n() const {
      return _n + 0;            // Convert -0 to +0
    }
    /**
     * The default constructor.
     *
     * This sets the angle to 0.
     **********************************************************************/
    Angle() : _s(0), _c(1), _n(0) {}
    /**
     * The general constructor.
     *
     * @param[in] y the \e y coordinate.
     * @param[in] x the \e x coordinate.
     *
     * \note the \e y coordinate is specified \e first.
     * \warning either \e x or \e y can be infinite, but not both.
     *
     * The point (\e x, \e y) is scaled so that it lies reasonably close to the
     * unit circle.
     *
     * normp means "already normalized"
     **********************************************************************/
    Angle(real s, real c, real num = 0, bool normp = false);
    explicit Angle(real deg);
    explicit operator real() const;
    static Angle degrees(real deg);
    real degrees() const;
    real degrees0() const;
    static Angle radians(real rad);
    real radians() const;
    real radians0() const;
    static Angle lam(real q);
    real lam() const;
    static Angle NaN();
    static Angle cardinal(real q);
    // ind == 0 => q has closest cardinal direction
    // ind even => q has closest even (N/S) cardinal direction
    // ind odd  => q has closest odd  (E/W) cardinal direction
    Angle nearest(unsigned ind = 0U) const;
    real ncardinal() const;
    static Angle eps();

    // Angle operator+() const;
    Angle operator-() const;
    Angle& operator+=(const Angle& p);
    Angle& operator-=(const Angle& p);
    Angle operator+(const Angle& p) const;
    Angle operator-(const Angle& p) const;
    bool zerop(real mult = 0) const;
    bool operator==(const Angle& p) const;
    Angle& round();
    Angle base() const;
    Angle rebase(const Angle& c) const;
    Angle& renormalize();
    Angle& setn(real n = 0);
    Angle& setquadrant(unsigned q);
    unsigned quadrant() const;
    Angle& reflect(bool flips, bool flipc = false, bool swapp = false);
    Angle flipsign(int mult) const;
    // Scale the sine component by m
    Angle modang(real m) const;
  };

  inline Angle::Angle(real s, real c, real num, bool normp)
    : _s(s)
    , _c(c)
    , _n(num)
  {
    using std::isfinite; using std::isnan; using std::isinf;
    using std::hypot; using std::copysign; using std::rint;
    _n = rint(_n);
    if (!normp) {
      real h = hypot(_s, _c);
      if (h == 0) {
        // If y is +/-0 and x = -0, +/-π is returned.
        // If y is +/-0 and x = +0, +/-0 is returned.
        // So retain the sign of _s = +/-0
        _c = copysign(real(1), _c);
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

  inline Angle::Angle(Math::real deg) {
    using std::rint;
    Math::sincosd(deg, _s, _c);
    _n = rint( (deg - Math::atan2d(_s, _c)) / Math::td );
  }

  inline Angle::operator Math::real() const {
    real d = Math::atan2d(_s, _c);
    // Preserve sign of +/-0
    return _n == 0 ? d : d + Math::td * _n;
  }

  inline Angle Angle::degrees(Math::real deg) {
    return Angle(deg);
  }

  inline Math::real Angle::degrees() const {
    return real(*this);
  }

  inline Math::real Angle::degrees0() const {
    return Math::atan2d(_s, _c);
  }

  inline Angle Angle::radians(Math::real rad) {
    using std::sin; using std::cos; using std::atan2; using std::rint;
    real sn = sin(rad), cs = cos(rad);
    return Angle(sn, cs, rint( (rad - atan2(sn, cs)) / (2 * Math::pi()) ),
                 true);
  }

  inline Math::real Angle::radians() const {
    real r = radians0();
    // Preserve sign of +/-0
    return _n == 0 ? r : r + 2 * Math::pi() * _n;
  }

  inline Math::real Angle::radians0() const {
    using std::atan2;
    return atan2(_s, _c);
  }

  inline Angle Angle::lam(Math::real psi) {
    using std::sinh;
    return Angle(sinh(psi), 1, 0);
  }

  inline Math::real Angle::lam() const {
    using std::asinh;
    return asinh(t());
  }

  inline Angle Angle::NaN() {
    return Angle(Math::NaN(), Math::NaN(), 0, true);
  }

  inline Math::real Angle::ncardinal() const {
    using std::signbit; using std::fabs;
    int iq = (signbit(_s) ? -1 : 1) * (signbit(_c) ?
                                       ( -_c >= fabs(_s) ? 2 : 1 ) :
                                       (  _c >= fabs(_s) ? 0 : 1 ));
    return 4 * _n + iq;
  }

  inline Angle Angle::eps() {
    return Angle(std::numeric_limits<real>::epsilon() / (1 << 20), 1, 0, true);
  }

  // Angle Angle::operator+() const { return *this; }
  inline Angle Angle::operator-() const {
    return Angle(-_s, _c, -_n, true);
  }

  inline Angle& Angle::operator+=(const Angle& p) {
    using std::rint;
    real q = ncardinal() + p.ncardinal();
    real c = _c * p._c - _s * p._s;
    _s = _s * p._c + _c * p._s;
    _c = c;
    _n += p._n;
    q -= ncardinal();
    _n += rint(q / 4);
    return *this;
  }

  inline Angle Angle::operator+(const Angle& p) const {
    Angle t = *this; t += p;
    return t;
  }

  inline Angle& Angle::operator-=(const Angle& p) {
    *this += -p;
    return *this;
  }

  inline Angle Angle::operator-(const Angle& p) const {
    Angle t = *this; t -= p;
    return t;
  }

  inline bool Angle::zerop(real mult) const {
    using std::fabs;
    return _n == 0 &&_c > 0 &&
      fabs(_s) <= mult * std::numeric_limits<real>::epsilon();
  }

  inline bool Angle::operator==(const Angle& p) const {
    Angle t = *this; t -= p;
    return t.zerop();
  }

  inline Angle& Angle::round() {
    _s = rnd(_s); _c = rnd(_c);
    return *this;
  }

  inline Angle Angle::base() const {
    return Angle(_s, _c, 0, true);
  }

  inline Angle Angle::rebase(const Angle& c) const {
    return (*this - c).base() + c;
  }

  inline Angle& Angle::renormalize() {
    using std::hypot;
    real h = hypot(_s, _c); _s /= h; _c /= h;
    return *this;
  }

  inline Angle& Angle::setn(Math::real n) {
    using std::rint;
    _n = rint(n);
    return *this;
  }

  inline Angle& Angle::setquadrant(unsigned q) {
    using std::copysign;
    _s = copysign(_s, real(             q  & 2U ? -1 : 1 ));
    _c = copysign(_c, real( ((q >> 1) ^ q) & 1U ? -1 : 1 ));
    return *this;
  }

  inline unsigned Angle::quadrant() const {
    using std::signbit;
    return 2U * signbit(_s) + (signbit(_c) ^ signbit(_s));
  }

  inline Angle& Angle::reflect(bool flips, bool flipc, bool swapp) {
    using std::swap;
    if (flips) _s *= -1;
    if (flipc) _c *= -1;
    if (swapp) swap(_s, _c);
    return *this;
  }

  inline Angle Angle::flipsign(int mult) const {
    return mult < 0 ? -*this : *this;
  }

  inline Angle Angle::modang(real m) const {
    using std::signbit;
    return signbit(m) ? Angle::NaN() :
      // Avoid nans if m == inf.
      Angle( _s * (m > 1 ? 1 : m),
             _c / (m > 1 ? m : 1),
             _n );
  }

  inline Angle Angle::cardinal(real q) {
    using std::isfinite; using std::rint; using std::remainder;
    if (!isfinite(q)) return Angle::NaN();
    q = rint(q);
    int iq = int(remainder(q, real(4)));
    // iq is in [-2, 2];
    // We could fold iq = -2 to iq = 2; but this way works too.
    real s, c, z = 0;
    switch (iq) {
    case -2: s = -z; c = -1; break;
    case -1: s = -1; c =  z; break;
    case  1: s =  1; c =  z; break;
    case  2: s =  z; c = -1; break;
    default: s =  z; c =  1; break; // iq = 0
    }
    return Angle(s, c, rint((q - iq) / 4));
  }

  inline Angle Angle::nearest(unsigned ind) const {
    using std::fabs; using std::copysign;
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

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ANGLE_HPP

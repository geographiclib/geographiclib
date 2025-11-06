/**
 * \file Angle.hpp
 * \brief Header for the GeographicLib::AngleT class
 *
 * The class provide an accurate representation of angle via 3 numbers, its
 * sine and cosine, and the number of turns.
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ANGLE_HPP)
#define GEOGRAPHICLIB_ANGLE_HPP 1

#include <string>
#include <GeographicLib/Math.hpp>

namespace GeographicLib {

  /**
   * \brief An accurate representation of angles.
   *
   * @tparam T the working floating point type.
   *
   * This class provides an accurate representation of angle via 3 numbers, its
   * sine = \e s and cosine = \e c, and the number of turns = \e n.  The angle
   * is then 2\e n &pi; + atan2(\e s, \e c).  This representation offers
   * several advantages:
   * - the cardinal directors (multiples of 90&deg;) are exactly represented (a
   *    benefit shared by representing angles as degrees);
   * - angles very close to any cardinal direction are accurately represented;
   * - there's no loss of precision with large angles (outside the "normal"
   *   range [&minus;180&deg;, +180&deg;]);
   * - various operations, such as adding a multiple of 90&deg; to an angle are
   *   performed exactly.
   * .
   * This representation does not favor degrees over radians.  However, the
   * one-argument constructor, AngleT(T), does the conversion from degrees and
   * the cast to T, AngleT::operator T(), returns the angle in degrees.  There
   * are alternatives, radians(T) and radians() const, to allow these
   * conversions with radians.
   *
   * N.B. \e n is stored as a real.  This allows it to be inf or nan.
   *
   * Example of use:
   * \include example-Angle.cpp
   **********************************************************************/
  template<typename T = Math::real>
  class AngleT {
    // No GEOGRAPHICLIB_EXPORT because this is a template class (like
    // PolygonAreaT).  Not sure why Accumulator needs GEOGRAPHICLIB_EXPORT.
  private:
    T _s, _c, _n;
    static T rnd(T x);
  public:
    /** \name Creating AngleT objects.
     **********************************************************************/
    ///@{
    /**
     * The default constructor.
     *
     * This sets the angle to 0.
     **********************************************************************/
    AngleT() : _s(0), _c(1), _n(0) {}
    /**
     * The general constructor.
     *
     * @param[in] s the sine component.
     * @param[in] c the cosine component.
     * @param[in] num the number of turns (default 0).
     * @param[in] normp are \e s and \e c normalized (default false).
     *
     * \warning either \e s or \e c can be infinite, but not both.
     *
     * By default, the point (\e s, \e c) is scaled to lie on the unit circle.
     * Setting \e normp = true skips this step; in this case (\e s, \e c)
     * should already lie on the unit circle.
     **********************************************************************/
    AngleT(T s, T c, T num = 0, bool normp = false);
    /**
     * The 1-argument constructor.
     *
     * @param[in] deg the angle in degrees.
     *
     * \note This is an explicit constructor to avoid accidental conversions.
     **********************************************************************/
    explicit AngleT(T deg);
    /**
     * The convert an angle in degrees to an AngleT.
     *
     * @param[in] deg the angle in degrees.
     * @return the AngleT.
     *
     * \note This mimics the behavior of AngleT(T deg);
     **********************************************************************/
    static AngleT degrees(T deg);
    /**
     * Convert an angle in radians to an AngleT.
     *
     * @param[in] rad the angle in radians.
     * @return the AngleT.
     *
     * \note This is the radians analog of degrees(T).
     **********************************************************************/
    static AngleT radians(T rad);
    /**
     * Convert an lambertian to an AngleT.
     *
     * @param[in] q the lambertian of the angle.
     * @return the AngleT.
     *
     * This sets the angle to atan(sinh(\e q)).
     **********************************************************************/
    static AngleT lam(T q);
    /**
     * Not an angle.
     *
     * @return the AngleT equivalent to not-a-number.
     **********************************************************************/
    static AngleT NaN();
    /**
     * A cardinal direction.
     *
     * @param[in] q the number of quarter turns.
     * @return the AngleT equivalent to \e q quarter turns.
     *
     * \e q is rounded to an integer and \e q = &plusmn;0 are distinguished.
     * NaN() is returned is \e q is not finite.
     **********************************************************************/
    static AngleT cardinal(T q);
    /**
     * Return a tiny angle.
     *
     * @return a tiny angle.
     *
     * This allows angles extremely close to the cardinal directions to be
     * generated.  The round() function will flush this angle to 0.
     **********************************************************************/
    static AngleT eps();
    ///@}

    /** \name Inspector functions.
     **********************************************************************/
    ///@{
    /**
     * @return the sine of the angle.
     **********************************************************************/
    T s() const { return _s; }
    /**
     * @return the cosine of the angle.
     **********************************************************************/
    T c() const { return _c; }
    /**
     * @return the tangent of the angle.
     **********************************************************************/
    T t() const { return _s/_c; }
    /**
     * @return the number of turns.
     **********************************************************************/
    T n() const {
      return _n + 0;            // Convert -0 to +0
    }
    /**
     * @return the number of turns treating &minus;180&deg; as +180&deg; less 1
     *   turn.
     **********************************************************************/
    T n0() const;
    ///@}

    /** \name Converting AngleT into other representations
     **********************************************************************/
    ///@{
    /**
     * Convert an AngleT to degrees via a type conversion.
     *
     * @return the angle in degrees.
     *
     * \note This is an explicit type conversion to avoid accidental
     * conversions.
     **********************************************************************/
    explicit operator T() const;
    /**
     * Convert an AngleT to degrees.
     *
     * @return the angle in degrees.
     *
     * \note This mimics the behavior of AngleT::operator T().
     **********************************************************************/
    T degrees() const;
    /**
     * Convert an AngleT to degrees ignoring the number of turns.
     *
     * @return the angle in degrees assuming n() is zero.
     **********************************************************************/
    T degrees0() const;
    /**
     * Convert an AngleT to radians.
     *
     * @return the angle in radians.
     *
     * \note This is the radians analog of degrees().
     **********************************************************************/
    T radians() const;
    /**
     * Convert an AngleT to radians ignoring the number of turns.
     *
     * @return the angle in radians assuming n() is zero.
     *
     * \note This is the radians analog of degrees0().
     **********************************************************************/
    T radians0() const;
    /**
     * Return the lambertian of the AngleT.
     *
     * @return the lambertian.
     *
     * The lambertian of &phi; is asinh tan &phi;.
     **********************************************************************/
    T lam() const;
    /**
     * Return the nearest cardinal direction as an AngleT.
     *
     * @param[in] ind an indicator.
     * @return the nearest cardinal direction as an AngleT.
     *
     * If \e ind == 0 (the default) the closest cardinal direction is returned.
     * Otherwise, if \e ind is even, the closest even (N/S) cardinal direction
     * is returned; or, if \e ind is odd, the closest odd (E/W) cardinal
     * direction is returned.
     **********************************************************************/
    AngleT nearest(unsigned ind = 0U) const;
    /**
     * Return the nearest cardinal direction as an integer.
     *
     * @return the nearest cardinal direction as an integer.
     *
     * \note This is the reverse of cardinal(T).
     **********************************************************************/
    T ncardinal() const;
    unsigned quadrant() const;
    ///@}

    /** \name Elementary arithmetic operations on AngleT
     **********************************************************************/
    ///@{
    /**
     * Return the negated AngleT.
     *
     * @return minus the AngleT
     **********************************************************************/
    AngleT operator-() const;
    /**
     * Implement the += operator.
     *
     * @param[in] p the AngleT to be added.
     * @return the current AngleT after the addition.
     **********************************************************************/
    AngleT& operator+=(const AngleT& p);
    /**
     * Implement the -= operator.
     *
     * @param[in] p the AngleT to be subtracted.
     * @return the current AngleT after the subtraction.
     **********************************************************************/
    AngleT& operator-=(const AngleT& p);
    /**
     * Implement the + operator.
     *
     * @param[in] p the AngleT to be added.
     * @return the result of the addition; the current AngleT is not modified.
     **********************************************************************/
    AngleT operator+(const AngleT& p) const;
    /**
     * Implement the - operator.
     *
     * @param[in] p the AngleT to be subtracted.
     * @return the result of the subtraction; the current AngleT is not
     *   modified.
     **********************************************************************/
    AngleT operator-(const AngleT& p) const;
    /**
     * Test for a zero angle.
     *
     * @param[in] mult multiplier of machine epsilon used in test (default 0).
     * @return true if this AngleT is withing \e mult &epsilon; of zero.
     **********************************************************************/
    bool zerop(T mult = 0) const;
    /**
     * Implement the == operator.
     *
     * @param[in] p the AngleT to be compared against
     * @return *this == \e p.
     **********************************************************************/
    bool operator==(const AngleT& p) const;
    ///@}

    /** \name Operations which modify a AngleT
     **********************************************************************/
    ///@{
    /**
     * "Round" the AngleT the == operator.
     *
     * @return the AngleT with tiny values of s() and c() set to &plusmn;0.
     *
     * This ensures that the smallest gaps between sine and cosine values is
     * &epsilon;/2048.
     **********************************************************************/
    AngleT& round();
    /**
     * Renormalize the sine and cosine values
     *
     * @return the modified AngleT.
     *
     * During arithmetic operations on AngleT object, no effort is mode to
     * ensure that (s(), c()) remains on the unit circle.  This function
     * corrects this.
     **********************************************************************/
    AngleT& renormalize();
    /**
     * Reduce the angle to [&minus;180&deg;, +180&deg;]
     *
     * @return the modified AngleT, obtained by setting n() to zero.
     **********************************************************************/
    AngleT& setn(T n = 0);
    /**
     * Reduce the angle to (&minus;180&deg;, +180&deg;]
     *
     * @return the modified AngleT.
     *
     * This differs from setn(T) by treating &minus;180&deg; as +180&deg; less
     * 1 turn.
     **********************************************************************/
    AngleT& setn0(T n = 0);
    /**
     * Set the quadrant of an AngleT the angle to (&minus;180&deg;, +180&deg;]
     *
     * @param[in] q the quadrant.
     * @return the modified AngleT.
     *
     * This sets the signs of s() and c() according to e q.
     *
     * \note Only the low two bits of \e q are used.  n() is unchanged.
     **********************************************************************/
    AngleT& setquadrant(unsigned q);
    /**
     * Reflect the angle is various ways
     *
     * @param[in] flips change the sign of s()
     * @param[in] flipc change the sign of c()
     * @param[in] swapp swap s() and c()
     * @return the modified AngleT.
     *
     * \note The operations are carried out in the order of the parameters.
     **********************************************************************/
    AngleT& reflect(bool flips, bool flipc = false, bool swapp = false);
    ///@}

    /** \name Operations which return a new AngleT
     **********************************************************************/
    ///@{
    /**
     * Return an AngleT in [&minus;180&deg;, +180&deg;].
     *
     * @return the new AngleT.
     *
     * This returns the AngleT with n() set to zero.
     **********************************************************************/
    AngleT base() const;
    /**
     * Return an AngleT centered about another AngleT
     *
     * @param[in] c the center AngleT
     * @return the new AngleT.
     *
     * This returns the result of adjusting n() so that the new AngleT is with
     * &plusmn;180&deg; of \e c.
     **********************************************************************/
    AngleT rebase(const AngleT& c) const;
    /**
     * Return an AngleT with the sign optionally flipped
     *
     * @param[in] mult
     * @return the new AngleT.
     *
     * return signbit(\e mult) ? -*this : *this.
     **********************************************************************/
    AngleT flipsign(T mult) const;
    /**
     * The "reduced latitude" operation.
     *
     * @param[in] m
     * @return the atan(m * tan(*this))
     *
     * However the quadrant of the result tracking that of *this through
     * multiples turns.
     **********************************************************************/
    AngleT modang(T m) const;
    ///@}

    /** \name Converting AngleT to and from a string representation
     **********************************************************************/
    ///@{
    /**
     * Interpret two strings as latitude and longitude.
     *
     * @param[in] stra the first string
     * @param[in] strb the second string
     * @param[out] lat the latitude
     * @param[out] lon the longitude
     * @param[in] longfirst (default false) whether the longitude is given
     *   first.
     *
     * In the absence of hemisphere indicators (N/S for latitude and E/W for
     * longitude), it is assumed that the first string is the latitude.
     * Setting \e longfirst = true uses the opposite convention.  The
     * hemisphere indicators can also be used to set the signs of the angles.
     **********************************************************************/
    static void DecodeLatLon(const std::string& stra, const std::string& strb,
                             AngleT& lat, AngleT& lon,
                             bool longfirst = false);
    /**
     * Interpret a string as azimuth
     *
     * @param[in] azistr the string representing the azimuth
     * @return the azimuth
     *
     * The hemisphere indicators E/W can be used to set the sign of the
     * azimuth.
     **********************************************************************/
    static AngleT DecodeAzimuth(const std::string& azistr);
    /**
     * Create a string for a latitude-longitude pair.
     *
     * @param[in] lat the latitude.
     * @param[in] lon the longitude.
     * @param[in] prec the precision relative to 1&deg;.
     * @param[in] dms (default false) whether to use degrees/minutes/seconds as
     *   opposed to decimal degrees
     * @param[in] dmssep (default NULL) the separator to use with the DMS
     *   representation instead of d ' ".
     * @param[in] longfirst (default false) whether to list the longitude
     * first.
     * @return string representation
     *
     * With dms = true the hemisphere indicators N/S and E/W are used to
     * indicator the signs of the latitude and longitude.
     **********************************************************************/
    static std::string LatLonString(AngleT lat, AngleT lon, int prec,
                                    bool dms = false, char dmssep = '\0',
                                    bool longfirst = false);
    /**
     * Create a string for an azimuth.
     *
     * @param[in] azi the azimuth.
     * @param[in] prec the precision relative to 1&deg;.
     * @param[in] dms (default false) whether to use degrees/minutes/seconds as
     *   opposed to decimal degrees
     * @param[in] dmssep (default NULL) the separator to use with the DMS
     *   representation instead of d ' ".
     * @return string representation
     *
     * With dms = true the hemisphere indicators and E/W is used to
     * indicator the sign of the azimuth.
     **********************************************************************/
    static std::string AzimuthString(AngleT azi, int prec,
                                     bool dms = false, char dmssep = '\0');
    ///@}
  };

  template<typename T>
  inline AngleT<T>::AngleT(T s, T c, T num, bool normp)
    : _s(s)
    , _c(c)
    , _n(num)
  {
    using std::isfinite, std::isnan, std::isinf, std::hypot,
      std::copysign, std::rint;
    _n = rint(_n);
    if (!normp) {
      // Cannot just use Math::norm because of all the special cases
      T h = hypot(_s, _c);
      if (h == 0) {
        // If y is +/-0 and x = -0, +/-pi is returned.
        // If y is +/-0 and x = +0, +/-0 is returned.
        // So retain the sign of _s = +/-0
        _c = copysign(T(1), _c);
      } else if (isfinite(h)) {
        _s /= h; _c /= h;
      } else if (isnan(h) || (isinf(_s) && isinf(_c)))
        _s = _c = Math::NaN();
      else if (isinf(_s)) {
        // infinite, finite
        _s = copysign(T(1), _s);
        _c = copysign(T(0), _c);
      } else {
        // isinf(cos); finite, infinite
        _s = copysign(T(0), _s);
        _c = copysign(T(1), _c);
      }
    }
  }

  template<typename T>
  inline AngleT<T>::AngleT(T deg) {
    using std::rint;
    Math::sincosd(deg, _s, _c);
    _n = rint( (deg - Math::atan2d(_s, _c)) / Math::td );
  }

  template<typename T>
  inline AngleT<T>::operator T() const {
    T d = degrees0();
    // Preserve sign of +/-0
    return _n == 0 ? d : d + Math::td * _n;
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::degrees(T deg) {
    return AngleT<T>(deg);
  }

  template<typename T>
  inline T AngleT<T>::degrees() const {
    return T(*this);
  }

  template<typename T>
  inline T AngleT<T>::degrees0() const {
    return Math::atan2d(_s, _c);
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::radians(T rad) {
    using std::sin, std::cos, std::atan2, std::rint;
    T sn = sin(rad), cs = cos(rad);
    return AngleT<T>(sn, cs, rint( (rad - atan2(sn, cs)) / (2 * Math::pi()) ),
                 true);
  }

  template<typename T>
  inline T AngleT<T>::radians() const {
    T r = radians0();
    // Preserve sign of +/-0
    return _n == 0 ? r : r + 2 * Math::pi() * _n;
  }

  template<typename T>
  inline T AngleT<T>::radians0() const {
    using std::atan2;
    return atan2(_s, _c);
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::lam(T psi) {
    using std::sinh;
    return AngleT<T>(sinh(psi), 1, 0);
  }

  template<typename T>
  inline T AngleT<T>::lam() const {
    using std::asinh;
    return asinh(t());
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::NaN() {
    return AngleT<T>(Math::NaN(), Math::NaN(), 0, true);
  }

  template<typename T>
  inline T AngleT<T>::ncardinal() const {
    using std::signbit, std::fabs;
    int iq = (signbit(_s) ? -1 : 1) * (signbit(_c) ?
                                       ( -_c >= fabs(_s) ? 2 : 1 ) :
                                       (  _c >= fabs(_s) ? 0 : 1 ));
    return 4 * _n + iq;
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::eps() {
    return AngleT<T>(std::numeric_limits<T>::epsilon() / (1 << 20), 1, 0, true);
  }

  // AngleT<T> AngleT<T>::operator+() const { return *this; }
  template<typename T>
  inline AngleT<T> AngleT<T>::operator-() const {
    return AngleT<T>(-_s, _c, -_n, true);
  }

  template<typename T>
  inline AngleT<T>& AngleT<T>::operator+=(const AngleT<T>& p) {
    using std::rint;
    T q = ncardinal() + p.ncardinal();
    T c = _c * p._c - _s * p._s;
    _s = _s * p._c + _c * p._s;
    _c = c;
    _n += p._n;
    q -= ncardinal();
    _n += rint(q / 4);
    return *this;
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::operator+(const AngleT<T>& p) const {
    AngleT<T> t = *this; t += p;
    return t;
  }

  template<typename T>
  inline AngleT<T>& AngleT<T>::operator-=(const AngleT<T>& p) {
    *this += -p;
    return *this;
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::operator-(const AngleT<T>& p) const {
    AngleT<T> t = *this; t -= p;
    return t;
  }

  template<typename T>
  inline bool AngleT<T>::zerop(T mult) const {
    using std::fabs;
    return _n == 0 &&_c > 0 &&
      fabs(_s) <= mult * std::numeric_limits<T>::epsilon();
  }

  template<typename T>
  inline bool AngleT<T>::operator==(const AngleT<T>& p) const {
    AngleT<T> t = *this; t -= p;
    return t.zerop();
  }

  template<typename T>
  inline AngleT<T>& AngleT<T>::round() {
    _s = rnd(_s); _c = rnd(_c);
    return *this;
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::base() const {
    return AngleT<T>(_s, _c, 0, true);
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::rebase(const AngleT<T>& c) const {
    // This is exact for c = cardinal direction
    // return (*this - c).base() + c;
    AngleT<T> t = *this;
    return t.setn0(((*this - c).base() + c).n0());
  }

  template<typename T>
  inline AngleT<T>& AngleT<T>::renormalize() {
    using std::hypot;
    T h = hypot(_s, _c); _s /= h; _c /= h;
    return *this;
  }

  template<typename T>
  inline AngleT<T>& AngleT<T>::setn(T n) {
    using std::rint;
    _n = rint(n);
    return *this;
  }

  template<typename T>
  inline T AngleT<T>::n0() const {
      using std::signbit;
      return (_n - (_s == 0 && signbit(_s) && _c < 0 ? 1 : 0)) + 0;
    }

  template<typename T>
  inline AngleT<T>& AngleT<T>::setn0(T n) {
    using std::rint, std::signbit;
    _n = rint(n) + (_s == 0 && signbit(_s) && _c < 0 ? 1 : 0);
    return *this;
  }

  template<typename T>
  inline AngleT<T>& AngleT<T>::setquadrant(unsigned q) {
    using std::copysign;
    _s = copysign(_s, T(             q  & 2U ? -1 : 1 ));
    _c = copysign(_c, T( ((q >> 1) ^ q) & 1U ? -1 : 1 ));
    return *this;
  }

  template<typename T>
  inline unsigned AngleT<T>::quadrant() const {
    using std::signbit;
    return 2U * signbit(_s) + (signbit(_c) ^ signbit(_s));
  }

  template<typename T>
  inline AngleT<T>& AngleT<T>::reflect(bool flips, bool flipc, bool swapp) {
    using std::swap;
    if (flips) _s *= -1;
    if (flipc) _c *= -1;
    if (swapp) swap(_s, _c);
    return *this;
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::flipsign(T mult) const {
    using std::signbit;
    return signbit(mult) ? -*this : *this;
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::modang(T m) const {
    using std::signbit;
    return signbit(m) ? AngleT<T>::NaN() :
      // Avoid nans if m == inf.
      AngleT<T>( _s * (m > 1 ? 1 : m),
             _c / (m > 1 ? m : 1),
             _n );
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::cardinal(T q) {
    using std::isfinite, std::rint, std::remainder;
    if (!isfinite(q)) return AngleT<T>::NaN();
    q = rint(q);
    int iq = int(remainder(q, T(4)));
    // iq is in [-2, 2];
    // We could fold iq = -2 to iq = 2; but this way works too.
    T s, c, z = 0;
    switch (iq) {
    case -2: s = -z; c = -1; break;
    case -1: s = -1; c =  z; break;
    case  1: s =  1; c =  z; break;
    case  2: s =  z; c = -1; break;
    default:
      // iq = 0, but distinguish q = +/-0
      s =  q != 0 ? z : q; c =  1;
      break;
    }
    return AngleT<T>(s, c, rint((q - iq) / 4), true);
  }

  template<typename T>
  inline AngleT<T> AngleT<T>::nearest(unsigned ind) const {
    using std::fabs, std::copysign;
    T s, c;
    if (ind == 0U) {
      if (fabs(_c) >= fabs(_s)) {
        s = copysign(T(0), _s); c = copysign(T(1), _c);
      } else {
        s = copysign(T(1), _s); c = copysign(T(0), _c);
      }
    } else if ((ind & 1U) == 0U) { // ind nonzero and even
      s = copysign(T(0), _s); c = copysign(T(1), _c);
    } else {                    // ind odd
      s = copysign(T(1), _s); c = copysign(T(0), _c);
    }
    return AngleT<T>(s, c, _n, true);
  }

  using Angle = AngleT<Math::real>;

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ANGLE_HPP

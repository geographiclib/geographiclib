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
    real n() const { return _n; }
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
    Angle& rnd();
    Angle rnded() const;
    Angle base() const;
    Angle rebase(const Angle& c) const;
    Angle& renormalize();
    Angle& setn(real n = 0);
    Angle& setquadrant(unsigned q);
    unsigned quadrant() const;
    Angle& reflect(bool flips, bool flipc = false, bool swapp = false);
    Angle flipsign(int mult) const;

    // Backward compatibility
    // static Angle degrees(real deg) { return Angle(deg); }
    // static Angle aux(int id, real s, real c, bool donorm);
    // real degrees0() const { return Math::atan2d(_s, _c); }
    // real radians0() const { using std::atan2; return atan2(_s, _c); }
    // real degrees() const { return degrees0(); }
    // real y() const { return s(); }
    // real x() const { return c(); }
    // Math::real& y() { return _s; }
    // Math::real& x() { return _c; }
    // Angle& normalize() { return *this; }
    // Angle normalized() const { return *this; }
    // Math::real tan() const { return t(); }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ANGLE_HPP

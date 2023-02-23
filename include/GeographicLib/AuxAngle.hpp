/**
 * \file AuxAngle.hpp
 * \brief Header for the GeographicLib::experimental::AuxAngle class.
 *
 * \note This is just sample code.  It is not part of GeographicLib itself.
 *
 * This file is an implementation of the methods described in
 * - C. F. F. Karney,
 *   On auxiliary latitudes,
 *   Technical Report, SRI International, December 2022.
 *   https://arxiv.org/abs/2212.05818
 * .
 * Copyright (c) Charles Karney (2022-2023) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_AUXANGLE_HPP)
#define GEOGRAPHICLIB_AUXANGLE_HPP 1

#include <GeographicLib/Math.hpp>

namespace GeographicLib {

  /**
   * \brief An accurate representation of angles.
   *
   * \note This is just sample code.  It is not part of GeographicLib itself.
   *
   * This class is an implementation of the methods described in
   * - C. F. F. Karney,
   *   On auxiliary latitudes,
   *   Technical Report, SRI International, December 2022.
   *   https://arxiv.org/abs/2212.05818
   *
   * An angle is represented be the \e y and \e x coordinates of a point in the
   * 2d plane.  The two coordinates are proportional to the sine and cosine of
   * the angle.  This allows angles close to the cardinal points to be
   * represented accurately.  It also saves on unnecessary recomputations of
   * trigonometric functions of the angle.  Only angles in [&minus;180&deg;,
   * 180&deg;] can be represented.  (A possible extension would be to keep
   * count of the number of turns.)
   *
   * @tparam T the floating-point type to use for real numbers.
   **********************************************************************/
  template<typename T = Math::real>
  class AuxAngle {
  public:
    /**
     * The floating-point type for real numbers.  This just connects to the
     * template parameters for the class.
     **********************************************************************/
    typedef T real;
    /**
     * The constructor.
     *
     * @param[in] y the \e y coordinate.
     * @param[in] x the \e x coordinate.
     *
     * \note the \e y coordinate is specified \e first.
     * \warning either \e x or \e y can be infinite, but not both.
     *
     * The defaults (\e x = 1 and \e y = 0) are such that
     * + no arguments gives an angle of 0;
     * + 1 argument specifies the tangent of the angle.
     **********************************************************************/
    AuxAngle(T y = 0, T x = 1) : _y(y), _x(x) {}
    /**
     * @return the \e y component.  This is the sine of the angle if the
     *   AuxAngle has been normalized.
     **********************************************************************/
    T y() const { return _y; }
    /**
     * @return the \e x component.  This is the cosine of the angle if the
     *   AuxAngle has been normalized.
     **********************************************************************/
    T x() const { return _x; }
    /**
     * @return a reference to the \e y component.  This allows this component
     *   to be altered.
     **********************************************************************/
    T& y() { return _y; }
    /**
     * @return a reference to the \e x component.  This allows this component
     *   to be altered.
     **********************************************************************/
    T& x() { return _x; }
    /**
     * @return the AuxAngle converted to the conventional angle measured in
     *   degrees.
     **********************************************************************/
    T degrees() const;
    /**
     * @return the AuxAngle converted to the conventional angle measured in
     *   radians.
     **********************************************************************/
    T radians() const;
    /**
     * @return the lambertian of the AuxAngle.
     **********************************************************************/
    T lam() const;
    /**
     * @return the lambertian of the AuxAngle in degrees.
     **********************************************************************/
    T lamd() const;
    /**
     * @return the tangent of the angle.
     **********************************************************************/
    T tan() const { return _y / _x; }
    /**
     * @return a new normalized AuxAngle with the point lying on the unit
     *   circle and the \e y and \e x components are equal to the sine and
     *   cosine of the angle.
     **********************************************************************/
    AuxAngle normalized() const;
    /**
     * Normalize the AuxAngle in place so that the \e y and \e x components are
     *   equal to the sine and cosine of the angle.
     **********************************************************************/
    void normalize() { *this = normalized(); }
    /**
     * Set the quadrant for the AuxAngle.
     *
     * @param[in] p the AuxAngle from which the quadrant information is taken.
     * @return the new AuxAngle in the same quadrant as \e p.
     **********************************************************************/
    AuxAngle copyquadrant(const AuxAngle& p) const;
    /**
     * Add an AuxAngle.
     *
     * @param[in] p the AuxAngle to be added.
     * @return a reference to the new AuxAngle.
     *
     * The addition is done in place, altering the current AuxAngle.
     *
     * \warning Neither *this nor \e p should have an infinite component.  If
     * necessary, invoke AuxAngle::normalize on these angles first.
     **********************************************************************/
    AuxAngle& operator+=(const AuxAngle& p);
    /**
     * Convert degrees to an AuxAngle.
     *
     * @param[in] d the angle measured in degrees.
     * @return the corresponding AuxAngle.
     *
     * This allows a new AuxAngle to be initialized as an angle in degrees with
     * @code
     *   AuxAngle<real> phi = AuxAngle<real>::degrees(d);
     * @endcode
     * This is the so-called "named constructor" idiom.
     **********************************************************************/
    static AuxAngle degrees(T d);
    /**
     * Convert radians to an AuxAngle.
     *
     * @param[in] r the angle measured in radians.
     * @return the corresponding AuxAngle.
     *
     * This allows a new AuxAngle to be initialized as an angle in radians with
     * @code
     *   AuxAngle<real> phi = AuxAngle<real>::radians(r);
     * @endcode
     * This is the so-called "named constructor" idiom.
     **********************************************************************/
    static AuxAngle radians(T r);
    /**
     * Convert lambertian to an AuxAngle.
     *
     * @param[in] psi the lambertian of the angle.
     * @return the corresponding AuxAngle.
     *
     * This allows a new AuxAngle to be initialized given the lambertian with
     * @code
     *   AuxAngle<real> chi = AuxAngle<real>::lam(psi);
     * @endcode
     * This is the so-called "named constructor" idiom.
     **********************************************************************/
    static AuxAngle lam(T psi);
    /**
     * Convert lambertian in degrees to an AuxAngle.
     *
     * @param[in] psid the lambertian of the angle in degrees.
     * @return the corresponding AuxAngle.
     *
     * This allows a new AuxAngle to be initialized given the lambertian with
     * @code
     *   AuxAngle<real> chi = AuxAngle<real>::lamd(psid);
     * @endcode
     * This is the so-called "named constructor" idiom.
     **********************************************************************/
    static AuxAngle lamd(T psid);
    /**
     * @return a "NaN" AuxAngle.
     **********************************************************************/
    static AuxAngle NaN();
    /**
     * Compute the absolute error in another angle.
     *
     * @tparam T1 the floating-point type of the other angle.
     * @param[in] p the other angle
     * @return the absolute error between p and *this considered as angles in
     *   radians.
     **********************************************************************/
    template<typename T1>
    T AbsError(const AuxAngle<T1>& p) const;
    /**
     * Compute the relative error in another angle.
     *
     * @tparam T1 the floating-point type of the other angle.
     * @param[in] p the other angle
     * @return the relative error between p.tan() and this->tan().
     **********************************************************************/
    template<typename T1>
    T RelError(const AuxAngle<T1>& p) const;
  private:
    T _y, _x;
  };

  /// \cond SKIP
  template<typename T>
  inline AuxAngle<T> AuxAngle<T>::degrees(T d) {
    T y, x;
    Math::sincosd(d, y, x);
    return AuxAngle(y, x);
  }

  template<typename T>
  inline AuxAngle<T> AuxAngle<T>::radians(T r) {
    using std::sin; using std::cos;
    return AuxAngle(sin(r), cos(r));
  }

  template<typename T>
  inline AuxAngle<T> AuxAngle<T>::lam(T psi) {
    using std::sinh;
    return AuxAngle(sinh(psi));
  }

  template<typename T>
  inline AuxAngle<T> AuxAngle<T>::lamd(T psid) {
    using std::sinh;
    return AuxAngle(sinh(psid * Math::degree<T>()));
  }

  template<typename T>
  inline T AuxAngle<T>::degrees() const {
    return Math::atan2d(_y, _x);
  }

  template<typename T>
  inline T AuxAngle<T>::radians() const {
    using std::atan2; return atan2(_y, _x);
  }

  template<typename T>
  inline T AuxAngle<T>::lam() const {
    using std::asinh; return asinh( tan() );
  }

  template<typename T>
  inline T AuxAngle<T>::lamd() const {
    using std::asinh; return asinh( tan() ) / Math::degree<T>();
  }

  template<typename T> template<typename T1>
  inline T AuxAngle<T>::AbsError(const AuxAngle<T1>& p) const {
    using std::fabs;
    return fabs((AuxAngle(-T(p.y()), T(p.x())) += *this).radians());
  }

  template<typename T> template<typename T1>
  inline T AuxAngle<T>::RelError(const AuxAngle<T1>& p) const {
    using std::fabs;
    return fabs((T(p.y()) / T(p.x()) - tan()) / tan());
  }
  /// \endcond

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_AUXANGLE_HPP

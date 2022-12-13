/**
 * \file AuxLatitude.hpp
 * \brief Header for the GeographicLib::AuxLatitude and GeographicLib::AuxAngle
 * classes.
 *
 * \note This is just sample code.  It is not part of GeographicLib itself.
 *
 * This file is an implementation of the methods described in
 * - C. F. F. Karney,
 *   On auxiliary latitudes,
 *   Technical Report, SRI International, December 2022.
 *   https://arxiv.org/abs/2212.05818
 * .
 * Copyright (c) Charles Karney (2022) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(AUXLATITUDE_HPP)
#define AUXLATITUDE_HPP 1

#include <GeographicLib/Math.hpp>

#if !defined(GEOGRAPHICLIB_AUXLATITUDE_ORDER)
/**
 * The order of the series approximation used in AuxLatitude.
 * GEOGRAPHICLIB_AUXLATITUDE_ORDER can be set to one of [4, 6, 8].
 **********************************************************************/
#  define GEOGRAPHICLIB_AUXLATITUDE_ORDER 6
#endif

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
   * represented accurately.  Only angles in [&minus;180&deg;, 180&deg;] can be
   * represented.  (A possible extension would be to keep count of the number
   * of turns.)
   *
   * @tparam T the floating-point type to use for real numbers.
   **********************************************************************/
  template<typename T = double>
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
    AuxAngle(real y = 0, real x = 1) : _y(y), _x(x) {}
    /**
     * @return the \e y component.  This is the sine of the angle if the
     *   AuxAngle has been normalized.
     **********************************************************************/
    real y() const { return _y; }
    /**
     * @return the \e x component.  This is the cosine of the angle if the
     *   AuxAngle has been normalized.
     **********************************************************************/
    real x() const { return _x; }
    /**
     * @return a reference to the \e y component.  This allows this component
     *   to be altered.
     **********************************************************************/
    real& y() { return _y; }
    /**
     * @return a reference to the \e x component.  This allows this component
     *   to be altered.
     **********************************************************************/
    real& x() { return _x; }
    /**
     * @return the AuxAngle converted to the conventional angle measured in
     *   degrees.
     **********************************************************************/
    real degrees() const;
    /**
     * @return the AuxAngle converted to the conventional angle measured in
     *   radians.
     **********************************************************************/
    real radians() const;
    /**
     * @return the tangent of the angle.
     **********************************************************************/
    real tan() const { return _y / _x; }
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
    static AuxAngle degrees(real d);
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
    static AuxAngle radians(real r);
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
    real AbsError(const AuxAngle<T1>& p) const;
    /**
     * Compute the relative error in another angle.
     *
     * @tparam T1 the floating-point type of the other angle.
     * @param[in] p the other angle
     * @return the relative error between p.tan() and this->tan().
     **********************************************************************/
    template<typename T1>
    real RelError(const AuxAngle<T1>& p) const;
  private:
    real _y, _x;
  };

  /// \cond SKIP
  template<typename T>
  inline AuxAngle<T> AuxAngle<T>::degrees(real d) {
    real y, x;
    Math::sincosd(d, y, x);
    return AuxAngle(y, x);
  }

  template<typename T>
  inline AuxAngle<T> AuxAngle<T>::radians(real r) {
    using std::sin; using std::cos;
    return AuxAngle(sin(r), cos(r));
  }

  template<typename T>
  inline T AuxAngle<T>::degrees() const {
    return Math::atan2d(_y, _x);
  }

  template<typename T>
  inline T AuxAngle<T>::radians() const {
    using std::atan2; return atan2(_y, _x);
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

  /**
   * \brief Conversions between auxiliary latitudes.
   *
   * \note This is just sample code.  It is not part of GeographicLib itself.
   *
   * This class is an implementation of the methods described in
   * - C. F. F. Karney,
   *   On auxiliary latitudes,
   *   Technical Report, SRI International, December 2022.
   *   https://arxiv.org/abs/2212.05818
   *
   * The provides accurate conversions between geographic (\e phi, &phi;),
   * parametric (\e beta, &beta;), geocentric (\e theta, &theta;), rectifying
   * (\e mu, &mu;), conformal (\e chi, &chi;), and authalic (\e xi, &xi;)
   * latitudes for an ellipsoid of revolution.  A latitude is represented by an
   * AuxAngle in order to maintain precision close to the poles.
   *
   * The class implements two methods for the conversion:
   * - Direct evaluation of the defining equations, the \e exact method.  These
   *   equations are formulated so as to preserve relative accuracy of the
   *   tangent of the latitude, ensuring high accuracy near the equator and the
   *   poles.  Newton's method is used for those conversions that can't be
   *   expressed in closed form.
   * - Expansions in powers of &e n, the third flattening, the \e series
   *   method.  This delivers full accuracy for abs(\e f) &le; 1/150.  Here, \e
   *   f is the flattening of the ellipsoid.
   *
   * The series method is the preferred method of conversion for any conversion
   * involving &mu;, &chi;, or &xi;, with abs(\e f) &le; 1/150.  The equations
   * for the conversions between &phi;, &beta;, and &theta; are sufficiently
   * simple that the exact method should be used for such conversions and also
   * for conversions with with abs(\e f) &gt; 1/150.
   *
   * @tparam T the floating-point type to use.
   *
   * Example of use:
   * \include example-AuxLatitude.cpp
   *
   * For more information on this projection, see \ref auxlat.
   **********************************************************************/
  template<typename T = double>
  class AuxLatitude {
  public:
    /**
     * The floating-point type for real numbers.  This just connects to the
     * template parameters for the class.
     **********************************************************************/
    typedef T real;
    /**
     * The type used to represent angles.
     **********************************************************************/
    typedef AuxAngle<real> angle;
    /**
     * The different auxiliary latitudes.
     **********************************************************************/
    enum aux {
      /**
       * Geographic latitude, \e phi, &phi;
       * @hideinitializer
       **********************************************************************/
      GEOGRAPHIC = 0,
      /**
       * Parametric latitude, \e beta, &beta;
       * @hideinitializer
       **********************************************************************/
      PARAMETRIC = 1,
      /**
       * %Geocentric latitude, \e theta, &theta;
       * @hideinitializer
       **********************************************************************/
      GEOCENTRIC = 2,
      /**
       * Rectifying latitude, \e mu, &mu;
       * @hideinitializer
       **********************************************************************/
      RECTIFYING = 3,
      /**
       * Conformal latitude, \e chi, &chi;
       * @hideinitializer
       **********************************************************************/
      CONFORMAL  = 4,
      /**
       * Authalic latitude, \e xi, &xi;
       * @hideinitializer
       **********************************************************************/
      AUTHALIC   = 5,
      /**
       * The total number of auxiliary latitudes
       * @hideinitializer
       **********************************************************************/
      AUXNUMBER  = 6,
      /**
       * An alias for GEOGRAPHIC
       * @hideinitializer
       **********************************************************************/
      COMMON = GEOGRAPHIC,
      /**
       * An alias for GEOGRAPHIC
       * @hideinitializer
       **********************************************************************/
      GEODETIC = GEOGRAPHIC,
      /**
       * An alias for PARAMETRIC
       * @hideinitializer
       **********************************************************************/
      REDUCED = PARAMETRIC,
    };
    /**
     * Constructor
     *
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.
     *
     * \note the constructor does not precompute the coefficients for the
     * Fourier series for the series conversions.  These are computed and saved
     * when first needed.
     **********************************************************************/
    AuxLatitude(real f);
    /**
     * Constructor
     *
     * @param[in] a equatorial radius.
     * @param[in] b polar semi-axis.
     *
     * \note the constructor does not precompute the coefficients for the
     * Fourier series for the series conversions.  These are computed and saved
     * when first needed.
     **********************************************************************/
    AuxLatitude(real a, real b);
    /**
     * Convert between any two auxiliary latitudes.
     *
     * @param[in] auxin an AuxLatitude::aux indicating the type of
     *   auxiliary latitude \e zeta.
     * @param[in] auxout an AuxLatitude::aux indicating the type of
     *   auxiliary latitude \e eta.
     * @param[in] zeta the input auxiliary latitude.
     * @param[in] series if true use the Taylor series instead of the exact
     *   equations [default false].
     * @return the output auxiliary latitude \e eta.
     *
     * With \e series = true, the Fourier coefficients for a specific \e auxin
     * and \e auxout are computed and saved on the first call; the saved
     * coefficients are used on subsequent calls.  The series method is
     * accurate for abs(\e f) &le; 1/150; for other \e f, the exact method
     * should be used
     **********************************************************************/
    angle Convert(int auxin, int auxout, const angle& zeta,
                  bool series = false) const;
    /**
     * Convert geographic latitude to an auxiliary latitude \e eta.
     *
     * @param[in] auxout an AuxLatitude::aux indicating the auxiliary
     *   latitude returned.
     * @param[in] phi the geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e eta) / d
     *   tan(\e phi).
     * @return the auxiliary latitude \e eta.
     *
     * This uses the exact equations.
     **********************************************************************/
    angle ToAuxiliary(int auxout, const angle& phi,
                      real* diff = nullptr) const;
    /**
     * Convert an auxiliary latitude \e zeta to geographic latitude.
     *
     * @param[in] auxin an AuxLatitude::aux indicating the type of
     *   auxiliary latitude \e zeta.
     * @param[in] zeta the input auxiliary latitude.
     * @param[out] niter optional pointer to the number of iterations.
     * @return the geographic latitude \e phi.
     *
     * This uses the exact equations.
     **********************************************************************/
    angle FromAuxiliary(int auxin, const angle& zeta,
                        int* niter = nullptr) const;
    /**
     * @return \e f, the flattening of the ellipsoid.
     **********************************************************************/
    real Flattening() const { return _f; }
    /**
     * The order of the series expansions.  This is set at compile time to
     * either 4, 6, or 8, by the preprocessor macro
     * GEOGRAPHICLIB_AUXLATITUDE_ORDER.
     * @hideinitializer
     **********************************************************************/
    static const int Lmax = GEOGRAPHICLIB_AUXLATITUDE_ORDER;
  private:
    /**
     * Convert geographic latitude to parametric latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e beta) / d
     *   tan(\e phi).
     * @return \e beta, the parametric latitude
     **********************************************************************/
    angle Parametric(const angle& phi, real* diff = nullptr) const;
    /**
     * Convert geographic latitude to geocentric latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e theta) / d
     *   tan(\e phi).
     * @return \e theta, the geocentric latitude.
     **********************************************************************/
    angle Geocentric(const angle& phi, real* diff = nullptr) const;
    /**
     * Convert geographic latitude to rectifying latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e mu) / d
     *   tan(\e phi).
     * @return \e mu, the rectifying latitude.
     **********************************************************************/
    angle Rectifying(const angle& phi, real* diff = nullptr) const;
    /**
     * Convert geographic latitude to conformal latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e chi) / d
     *   tan(\e phi).
     * @return \e chi, the conformal latitude.
     **********************************************************************/
    angle Conformal(const angle& phi, real* diff = nullptr) const;
    /**
     * Convert geographic latitude to authalic latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e xi) / d
     *   tan(\e phi).
     * @return \e xi, the authalic latitude.
     **********************************************************************/
    angle Authalic(const angle& phi, real* diff = nullptr) const;
    // Maximum number of iterations for Newton's method
    static const int numit_ = 1000;
    real tol_, bmin_, bmax_;         // Static consts for Newton's method
    // Ellipsoid parameters
    real _f, _fm1, _e2, _e2m1, _e12, _e12p1, _n, _e, _e1, _n2, _q;
    // To hold computed Fourier coefficients
    mutable real _c[Lmax * AUXNUMBER * AUXNUMBER];
    // 1d index into AUXNUMBER x AUXNUMBER data
    static int ind(int auxout, int auxin) {
      return (auxout >= 0 && auxout < AUXNUMBER &&
              auxin  >= 0 && auxin  < AUXNUMBER) ?
        AUXNUMBER * auxout + auxin : -1;
    }
    // the function sqrt(1 + tphi^2), convert tan to sec
    static real sc(real tphi)
    { using std::hypot; return hypot(real(1), tphi); }
    // the function tphi / sqrt(1 + tphi^2), convert tan to sin
    static real sn(real tphi) {
      using std::isfinite; using std::isnan; using std::copysign;
      return isfinite(tphi) || isnan(tphi) ? tphi / sc(tphi) :
        copysign(real(1), tphi);
    }
    // The symmetric elliptic integral RD
    static real RD(real x, real y, real z);
    // The symmetric elliptic integral RF
    static real RF(real x, real y, real z);
    // the function atanh(e * sphi)/e; works for e^2 = 0 and e^2 < 0
    real atanhee(real tphi) const;
    // the function atanh(e * sphi)/e + sphi / (1 - (e * sphi)^2);
    real q(real tphi) const;
    // The divided difference of (q(1) - q(sphi)) / (1 - sphi)
    real Dq(real tphi) const;
    // Populate [_c[Lmax * k], _c[Lmax * (k + 1)])
    void fillcoeff(int auxin, int auxout, int k) const;
    // Clenshaw applied to sum(c[k] * sin( (2*k+2) * zeta), i, 0, K-1)
    // if alt, use the Reinsch optimizations
    static real Clenshaw(real szeta, real czeta, const real c[], int K,
                         bool alt = true);
  };

} // namespace GeographicLib

#endif  // AUXLATITUDE_HPP

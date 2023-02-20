/**
 * \file AuxLatitude.hpp
 * \brief Header for the GeographicLib::experimental::AuxLatitude class.
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

#if !defined(GEOGRAPHICLIB_AUXLATITUDE_HPP)
#define GEOGRAPHICLIB_AUXLATITUDE_HPP 1

#include <GeographicLib/Math.hpp>
#include "AuxAngle.hpp"

#if !defined(GEOGRAPHICLIB_AUXLATITUDE_ORDER)
/**
 * The order of the series approximation used in AuxLatitude.
 * GEOGRAPHICLIB_AUXLATITUDE_ORDER can be set to one of [4, 6, 8].  Use order
 * appropriate for double precision, 6, also for GEOGRAPHICLIB_PRECISION == 5
 * to enable truncation errors to be measured easily.
 **********************************************************************/
#  define GEOGRAPHICLIB_AUXLATITUDE_ORDER \
  (GEOGRAPHICLIB_PRECISION == 2 || GEOGRAPHICLIB_PRECISION == 5 ? 6 : \
   (GEOGRAPHICLIB_PRECISION == 1 ? 4 : 8))
#endif

namespace GeographicLib {
namespace experimental {

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
   * latitudes for an ellipsoid of revolution.  A latitude is represented by
   * the class AuxAngle in order to maintain precision close to the poles.
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
   **********************************************************************/
  template<typename T = Math::real>
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
    typedef AuxAngle<T> angle;
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
      PHI = GEOGRAPHIC,
      /**
       * An alias for PARAMETRIC
       * @hideinitializer
       **********************************************************************/
      BETA = PARAMETRIC,
      /**
       * An alias for GEOCENTRIC
       * @hideinitializer
       **********************************************************************/
      THETA = GEOCENTRIC,
      /**
       * An alias for RECTIFYING
       * @hideinitializer
       **********************************************************************/
      MU = RECTIFYING,
      /**
       * An alias for CONFORMAL
       * @hideinitializer
       **********************************************************************/
      CHI = CONFORMAL,
      /**
       * An alias for AUTHALIC
       * @hideinitializer
       **********************************************************************/
      XI = AUTHALIC,
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
    AuxLatitude(T f);
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
    AuxLatitude(T a, T b);
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
     * should be used.
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
                      T* diff = nullptr) const;
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
     * Return the rectifying radius.
     *
     * @param[in] a the equatorial radius.
     * @param[in] series if true use the Taylor series instead of the exact
     *   expression [default false].
     * @return the rectifying radius in the same units as \e a.
     **********************************************************************/
    T RectifyingRadius(T a, bool series = false) const;
    /**
     * Return the authalic radius squared.
     *
     * @param[in] a the equatorial radius.
     * @param[in] series if true use the Taylor series instead of the exact
     *   expression [default false].
     * @return the authalic radius squared in the same units as \e a.
     **********************************************************************/
    T AuthalicRadiusSquared(T a, bool series = false) const;
    /**
     * @return \e f, the flattening of the ellipsoid.
     **********************************************************************/
    T Flattening() const { return _f; }
    /**
     * The order of the series expansions.  This is set at compile time to
     * either 4, 6, or 8, by the preprocessor macro
     * GEOGRAPHICLIB_AUXLATITUDE_ORDER.
     * @hideinitializer
     **********************************************************************/
    static const int Lmax = GEOGRAPHICLIB_AUXLATITUDE_ORDER;
  private:
    // Maximum number of iterations for Newton's method
    static const int numit_ = 1000;
    T tol_, bmin_, bmax_;       // Static consts for Newton's method
    // the function atanh(e * sphi)/e + sphi / (1 - (e * sphi)^2);
  protected:
    /**
     * Convert geographic latitude to parametric latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e beta) / d
     *   tan(\e phi).
     * @return \e beta, the parametric latitude
     **********************************************************************/
    angle Parametric(const angle& phi, T* diff = nullptr) const;
    /**
     * Convert geographic latitude to geocentric latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e theta) / d
     *   tan(\e phi).
     * @return \e theta, the geocentric latitude.
     **********************************************************************/
    angle Geocentric(const angle& phi, T* diff = nullptr) const;
    /**
     * Convert geographic latitude to rectifying latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e mu) / d
     *   tan(\e phi).
     * @return \e mu, the rectifying latitude.
     **********************************************************************/
    angle Rectifying(const angle& phi, T* diff = nullptr) const;
    /**
     * Convert geographic latitude to conformal latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e chi) / d
     *   tan(\e phi).
     * @return \e chi, the conformal latitude.
     **********************************************************************/
    angle Conformal(const angle& phi, T* diff = nullptr) const;
    /**
     * Convert geographic latitude to authalic latitude
     *
     * @param[in] phi geographic latitude.
     * @param[out] diff optional pointer to the derivative d tan(\e xi) / d
     *   tan(\e phi).
     * @return \e xi, the authalic latitude.
     **********************************************************************/
    angle Authalic(const angle& phi, T* diff = nullptr) const;
    // Ellipsoid parameters
    T _f, _fm1, _e2, _e2m1, _e12, _e12p1, _n, _e, _e1, _n2, _q;
    // To hold computed Fourier coefficients
    mutable T _c[Lmax * AUXNUMBER * AUXNUMBER];
    // 1d index into AUXNUMBER x AUXNUMBER data
    static int ind(int auxout, int auxin) {
      return (auxout >= 0 && auxout < AUXNUMBER &&
              auxin  >= 0 && auxin  < AUXNUMBER) ?
        AUXNUMBER * auxout + auxin : -1;
    }
    // the function sqrt(1 + tphi^2), convert tan to sec
    static T sc(T tphi)
    { using std::hypot; return hypot(T(1), tphi); }
    // the function tphi / sqrt(1 + tphi^2), convert tan to sin
    static T sn(T tphi) {
      using std::isinf; using std::copysign;
      return isinf(tphi) ? copysign(T(1), tphi) : tphi / sc(tphi);
    }
    // The symmetric elliptic integral RD
    static T RD(T x, T y, T z);
    // The symmetric elliptic integral RF
    static T RF(T x, T y, T z);
    // The symmetric elliptic integral RF
    static T RG(T x, T y);
    // Populate [_c[Lmax * k], _c[Lmax * (k + 1)])
    void fillcoeff(int auxin, int auxout, int k) const;
    // Clenshaw applied to sum(c[k] * sin( (2*k+2) * zeta), i, 0, K-1);
    // if !sinp then subst sine->cosine.
    static T Clenshaw(bool sinp, T szeta, T czeta, const T c[], int K);
    // the function atanh(e * sphi)/e; works for e^2 = 0 and e^2 < 0
    T atanhee(T tphi) const;
  private:
    T q(T tphi) const;
    // The divided difference of (q(1) - q(sphi)) / (1 - sphi)
    T Dq(T tphi) const;
  };

} // namespace experimental
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_AUXLATITUDE_HPP

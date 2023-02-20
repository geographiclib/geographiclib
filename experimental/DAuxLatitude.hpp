/**
 * \file DAuxLatitude.hpp
 * \brief Header for the GeographicLib::experimental::DAuxLatitude class.
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

#if !defined(GEOGRAPHICLIB_DAUXLATITUDE_HPP)
#define GEOGRAPHICLIB_DAUXLATITUDE_HPP 1

#include "AuxLatitude.hpp"

namespace GeographicLib {
namespace experimental {

  /**
   * \brief Divided differences of auxiliary latitudes.
   *
   * \note This is just sample code.  It is not part of GeographicLib itself.
   *
   * This class is an implementation of the methods described in
   * - C. F. F. Karney,
   *   On auxiliary latitudes,
   *   Technical Report, SRI International, December 2022.
   *   https://arxiv.org/abs/2212.05818
   *
   * @tparam T the floating-point type to use.
   **********************************************************************/
  template<typename T
           /// \cond SKIP
           = Math::real
           /// \endcond
           >
  class DAuxLatitude : public AuxLatitude<T> {
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
     * The base class.
     **********************************************************************/
    typedef AuxLatitude<T> aux;
    /**
     * Constructor
     *
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.
     **********************************************************************/
    DAuxLatitude(T f) : AuxLatitude<T>(f) {}
    /**
     * Constructor
     *
     * @param[in] a equatorial radius.
     * @param[in] b polar semi-axis.
     **********************************************************************/
    DAuxLatitude(T a, T b) : AuxLatitude<T>(a, b) {}
    /**
     * The divided difference of one auxiliary latitude with respect to
     * another.
     *
     * @param[in] auxin an AuxLatitude::aux indicating the type of
     *   auxiliary latitude \e zeta.
     * @param[in] auxout an AuxLatitude::aux indicating the type of
     *   auxiliary latitude \e eta.
     * @param[in] zeta1 the first of the input auxiliary latitudeas.
     * @param[in] zeta2 the second of the input auxiliary latitude.
     * @return the divided difference (\e eta2 - \e eta1) / (\e zeta2 - \e
     *   zeta1).
     *
     * \note This routine uses the series method.
     *
     * In the expression for the divided difference above, the angle quantities
     * should be understood as the conventional measure of angle (either in
     * radians or in degrees).
     *
     * The Fourier coefficients for a specific \e auxin and \e auxout are
     * computed and saved on the first call; the saved coefficients are used on
     * subsequent calls.  The series method is accurate for abs(\e f) &le;
     * 1/150.
     **********************************************************************/
    T DConvert(int auxin, int auxout,
               const angle& zeta1, const angle& zeta2) const;
    T DRectifying(const angle& phi1, const angle& phi2) const;
    // overflow for tphi1, tphi2 >~ sqrt(mx)
    T DIsometric(const angle& phi1, const angle& phi2) const;
    T DParametric(const angle& phi1, const angle& phi2) const;
    // Divided difference: (eta2 - eta1) / Delta.  Delta is EITHER 1, giving
    // the plain difference OR (zeta2 - zeta1) in radians, giving the divided
    // difference.  Other values will give nonsense.
    static T DClenshaw(bool sinp, T Delta,
                       T szet1, T czet1, T szet2, T czet2,
                       const T c[], int K);
    // Dasinh(x, y) / Datan(x, y)
    // overflow for x, y >~ sqrt(mx)
    static T Dlam(T x, T y) {
      using std::isnan; using std::isinf;
      return x == y ? aux::sc(x) :
        (isnan(x) || isnan(y) ? std::numeric_limits<T>::quiet_NaN() :
         (isinf(x) || isinf(y) ? std::numeric_limits<T>::infinity() :
          Dasinh(x, y) / Datan(x, y)));
    }
    // Dp0Dpsi in terms of chi
    static T Dp0Dpsi(T x, T y) {
      using std::isnan; using std::isinf; using std::copysign;
      return x == y ? aux::sn(x) :
        (isnan(x + y) ? x + y : // N.B. nan for inf-inf
         (isinf(x) ? copysign(T(1), x) :
          (isinf(y) ? copysign(T(1), y) :
           Dasinh(h(x), h(y)) * Dh(x, y) / Dasinh(x, y))));
    }
  protected:                    // so TestAux can access these functions
    /* NOT USED!
    // (sc(y) - sc(x)) / (y - x)
    static T Dsc(T x, T y) {
      using std::isnan; using std::isinf; using std::copysign;
      return isnan(x + y) ? x + y : // N.B. nan for inf-inf
        (isinf(x) ? copysign(T(1), x) :
         (isinf(y) ? copysign(T(1), y) :
          (x + y) / (aux::sc(x) + aux::sc(y))));
    }
    */
    // (sn(y) - sn(x)) / (y - x)
    static T Dsn(T x, T y) {
      T sc1 = aux::sc(x);
      if (x == y) return 1 / (sc1 * (1 + x*x));
      T sc2 = aux::sc(y), sn1 = aux::sn(x), sn2 = aux::sn(y);
      return x * y > 0 ?
        (sn1/sc2 + sn2/sc1) / ((sn1 + sn2) * sc1 * sc2) :
        (sn2 - sn1) / (y - x);
    }
    static T Datan(T x, T y) {
      using std::isinf; using std::atan;
      T d = y - x, xy = x*y;
      return x == y ? 1 / (1 + xy) :
        (isinf(xy) && xy > 0 ? 0 :
         (2 * xy > -1 ? atan( d / (1 + xy) ) : atan(y) - atan(x)) / d);
    }
    static T Dasinh(T x, T y) {
      using std::isinf; using std::asinh;
      T d = y - x, xy = x*y, hx = aux::sc(x), hy = aux::sc(y);
      // KF formula for x*y < 0 is asinh(y*hx - x*hy) / (y - x)
      // but this has problem if x*y overflows to -inf
      return x == y ? 1 / hx :
        (isinf(d) ? 0 :
         (xy > 0 ? asinh(d * (x*y < 1 ? (x + y) / (x*hy + y*hx) :
                              (1/x + 1/y) / (hy/y + hx/x))) :
          asinh(y) - asinh(x)) / d);
    }
    T Datanhee(T tphi1, T tphi2) const;
    // h(tan(x)) = tan(x) * sin(x) / 2
    static T h(T x) { return x * aux::sn(x) / 2; }
    static T Dh(T x, T y) {
      using std::isnan; using std::isinf; using std::copysign;
      if (isnan(x + y))
        return x + y;           // N.B. nan for inf-inf
      if (isinf(x))
        return copysign(1/T(2), x);
      if (isinf(y))
        return copysign(1/T(2), y);
      T sx = aux::sn(x), sy = aux::sn(y), d = sx*x + sy*y;
      if (d / 2 == 0)
        return (x + y) / 2;     // Handle underflow
      if (x * y <= 0)
        return (h(y) - h(x)) / (y - x); // Does not include x = y = 0
      T scx = aux::sc(x), scy = aux::sc(y);
      return ((x + y) / (2 * d)) *
        (Math::sq(sx*sy) + Math::sq(sy/scx) + Math::sq(sx/scy));
    }
  private:
    static T Dsin(T x, T y) {
      using std::sin; using std::cos;
      T d = (x - y) / 2;
      return cos((x + y)/2) * (d != 0 ? sin(d) / d : 1);
    }
    // (E(x) - E(y)) / (x - y)
    T DE(const angle& X, const angle& Y) const;
  };

} // namespace experimental
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_DAUXLATITUDE_HPP

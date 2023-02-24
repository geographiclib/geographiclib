/**
 * \file DAuxLatitude.hpp
 * \brief Header for the GeographicLib::DAuxLatitude class.
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

#include <GeographicLib/AuxLatitude.hpp>

namespace GeographicLib {

  /**
   * \brief Divided differences of auxiliary latitudes.
   *
   * This class is an implementation of the methods described in
   * - C. F. F. Karney,
   *   On auxiliary latitudes,
   *   Technical Report, SRI International, December 2022.
   *   https://arxiv.org/abs/2212.05818
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT DAuxLatitude : public AuxLatitude {
  private:
    typedef Math::real real;
    typedef AuxLatitude base;
  public:
    /**
     * Constructor
     *
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.
     **********************************************************************/
    DAuxLatitude(real f) : AuxLatitude(f) {}
    /**
     * Constructor
     *
     * @param[in] a equatorial radius.
     * @param[in] b polar semi-axis.
     **********************************************************************/
    DAuxLatitude(real a, real b) : AuxLatitude(a, b) {}
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
    real DConvert(int auxin, int auxout,
                  const AuxAngle& zeta1, const AuxAngle& zeta2) const;
    real DRectifying(const AuxAngle& phi1, const AuxAngle& phi2) const;
    // overflow for tphi1, tphi2 >~ sqrt(mx)
    real DIsometric(const AuxAngle& phi1, const AuxAngle& phi2) const;
    real DParametric(const AuxAngle& phi1, const AuxAngle& phi2) const;
    // Divided difference: (eta2 - eta1) / Delta.  Delta is EITHER 1, giving
    // the plain difference OR (zeta2 - zeta1) in radians, giving the divided
    // difference.  Other values will give nonsense.
    static real DClenshaw(bool sinp, real Delta,
                          real szet1, real czet1, real szet2, real czet2,
                          const real c[], int K);
    // Dasinh(x, y) / Datan(x, y)
    // overflow for x, y >~ sqrt(mx)
    static real Dlam(real x, real y) {
      using std::isnan; using std::isinf;
      return x == y ? base::sc(x) :
        (isnan(x) || isnan(y) ? std::numeric_limits<real>::quiet_NaN() :
         (isinf(x) || isinf(y) ? std::numeric_limits<real>::infinity() :
          Dasinh(x, y) / Datan(x, y)));
    }
    // Dp0Dpsi in terms of chi
    static real Dp0Dpsi(real x, real y) {
      using std::isnan; using std::isinf; using std::copysign;
      return x == y ? base::sn(x) :
        (isnan(x + y) ? x + y : // N.B. nan for inf-inf
         (isinf(x) ? copysign(real(1), x) :
          (isinf(y) ? copysign(real(1), y) :
           Dasinh(h(x), h(y)) * Dh(x, y) / Dasinh(x, y))));
    }
  protected:                    // so TestAux can access these functions
    /// \cond SKIP
    // (sn(y) - sn(x)) / (y - x)
    static real Dsn(real x, real y) {
      real sc1 = base::sc(x);
      if (x == y) return 1 / (sc1 * (1 + x*x));
      real sc2 = base::sc(y), sn1 = base::sn(x), sn2 = base::sn(y);
      return x * y > 0 ?
        (sn1/sc2 + sn2/sc1) / ((sn1 + sn2) * sc1 * sc2) :
        (sn2 - sn1) / (y - x);
    }
    static real Datan(real x, real y) {
      using std::isinf; using std::atan;
      real d = y - x, xy = x*y;
      return x == y ? 1 / (1 + xy) :
        (isinf(xy) && xy > 0 ? 0 :
         (2 * xy > -1 ? atan( d / (1 + xy) ) : atan(y) - atan(x)) / d);
    }
    static real Dasinh(real x, real y) {
      using std::isinf; using std::asinh;
      real d = y - x, xy = x*y, hx = base::sc(x), hy = base::sc(y);
      // KF formula for x*y < 0 is asinh(y*hx - x*hy) / (y - x)
      // but this has problem if x*y overflows to -inf
      return x == y ? 1 / hx :
        (isinf(d) ? 0 :
         (xy > 0 ? asinh(d * (x*y < 1 ? (x + y) / (x*hy + y*hx) :
                              (1/x + 1/y) / (hy/y + hx/x))) :
          asinh(y) - asinh(x)) / d);
    }
    real Datanhee(real tphi1, real tphi2) const;
    // h(tan(x)) = tan(x) * sin(x) / 2
    static real h(real x) { return x * base::sn(x) / 2; }
    static real Dh(real x, real y) {
      using std::isnan; using std::isinf; using std::copysign;
      if (isnan(x + y))
        return x + y;           // N.B. nan for inf-inf
      if (isinf(x))
        return copysign(1/real(2), x);
      if (isinf(y))
        return copysign(1/real(2), y);
      real sx = base::sn(x), sy = base::sn(y), d = sx*x + sy*y;
      if (d / 2 == 0)
        return (x + y) / 2;     // Handle underflow
      if (x * y <= 0)
        return (h(y) - h(x)) / (y - x); // Does not include x = y = 0
      real scx = base::sc(x), scy = base::sc(y);
      return ((x + y) / (2 * d)) *
        (Math::sq(sx*sy) + Math::sq(sy/scx) + Math::sq(sx/scy));
    }
    /// \endcond
  private:
    static real Dsin(real x, real y) {
      using std::sin; using std::cos;
      real d = (x - y) / 2;
      return cos((x + y)/2) * (d != 0 ? sin(d) / d : 1);
    }
    // (E(x) - E(y)) / (x - y)
    real DE(const AuxAngle& X, const AuxAngle& Y) const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_DAUXLATITUDE_HPP

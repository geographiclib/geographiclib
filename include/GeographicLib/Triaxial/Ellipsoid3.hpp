/**
 * \file Ellipsoid3.hpp
 * \brief Header for GeographicLib::Triaxial::Ellipsoid3 class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ELLIPSOID3_HPP)
#define GEOGRAPHICLIB_ELLIPSOID3_HPP 1

#include <iostream>
#include <array>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Angle.hpp>

namespace GeographicLib {
  /**
   * \brief Namespace for operations on triaxial ellipsoids
   *
   * The algorithms on triaxial ellipsoids are, for the most part, distinct
   * from those in the rest of %GeographicLib, which deal with biaxial
   * ellipsoids; it is therefore convenient to put them in a distinct
   * namespace.  The header files are included via
   * GeographicLib/Triaxial/Class.hpp.
   **********************************************************************/
  namespace Triaxial {
  class GEOGRAPHICLIB_EXPORT Ellipsoid3 {
  public:
    using vec3 = std::array<Math::real, 3>;
  private:
    friend class Cartesian3;  // For access to cart2toellipint normvec
    friend class Geodesic3;   // For Flip
    using real = Math::real;
    using ang = Angle;
    static void normvec(vec3& r) {
      real h = Math::hypot3(r[0], r[1], r[2]);
      // No checking for h = 0.  Result will be NaNs
      r[0] /= h; r[1] /= h; r[2] /= h;
    }
    static void Flip(ang& bet, ang& omg, ang& alp) {
      bet.reflect(false, true);
      omg.reflect(true);
      alp.reflect(true, true);
    }
    real _a, _b, _c;            // semi-axes
    real _e2, _k2, _kp2, _k, _kp;
    bool _oblate, _prolate, _biaxial;
    void cart2toellipint(vec3 r, ang& bet, ang& omg, vec3 axes) const;
  public:
    Ellipsoid3();
    Ellipsoid3(real a, real b, real c);
    Ellipsoid3(real b, real e2, real k2, real kp2);
    void Norm(vec3& r) const;
    void Norm(vec3& r, vec3& v) const;
    real a() const { return _a; }
    real b() const { return _b; }
    real c() const { return _c; }
    real e2() const { return _e2; }
    real k2() const { return _k2; }
    real kp2() const { return _kp2; }
    real k() const { return _k; }
    real kp() const { return _kp; }
    bool oblate() const { return _oblate; }
    bool prolate() const { return _prolate; }
    bool biaxial() const { return _biaxial; }
    static bool AngNorm(Angle& bet, Angle& omg, Angle& alp,
                        bool alt = false) {
      using std::signbit;
      // If !alt, put bet in [-pi/2,pi/2]
      // If  alt, put omg in [0, pi]
      bool flip = alt ? signbit(omg.s()) : signbit(bet.c());
      if (flip)
        Flip(bet, omg, alp);
      if (0) {
        if (bet.c() == 0 && bet.s() * alp.c() > 0)
          alp.reflect(true, true);
        if (bet.c() == 0 && alp.c() == 0)
          alp.reflect(alp.s() * bet.s() > 0); // alp.s() = -bet.s();
      }
      return flip;
    }
    static bool AngNorm(Angle& bet, Angle& omg,
                        bool alt = false) {
      using std::signbit;
      // If !alt, put bet in [-pi/2,pi/2]
      // If  alt, put omg in [0, pi]
      bool flip = alt ? signbit(omg.s()) : signbit(bet.c());
      if (flip) {
        ang alp;
        Flip(bet, omg, alp);
      }
      return flip;
    }
    void cart2toellip(vec3 r, Angle& bet, Angle& omg) const;
    void cart2toellip(vec3 r, vec3 v,
                      Angle& bet, Angle& omg, Angle& alp) const;
    void cart2toellip(Angle bet, Angle omg,
                      vec3 v, Angle& alp) const;
    void elliptocart2(Angle bet, Angle omg, vec3& r) const;
    void elliptocart2(Angle bet, Angle omg, Angle alp,
                      vec3& r, vec3& v) const;
    /**
     * A global instantiation of Ellipsoid3 with the parameters for the
     * Earth.
     **********************************************************************/
    static const Ellipsoid3& Earth();

  };

  } // namespace Triaxial
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIAL_HPP

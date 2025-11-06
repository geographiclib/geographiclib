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
   * namespace.  In order to minimize the change of name clashes, the classes
   * in the namespace include "3" in their names.  The header files are
   * included via GeographicLib/Triaxial/Class.hpp.
   **********************************************************************/
  namespace Triaxial {
  /**
   * \brief A triaxial ellipsoid.
   *
   * The class holds the basic information about a triaxial ellipsoid
   * given by
   * \f[ S(\mathbf R) =
   * \frac{X^2}{a^2} + \frac{Y^2}{b^2} + \frac{Z^2}{c^2} - 1 = 0, \f]
   * where the semiaxes satisfy \f$ a \ge b \ge c > 0\f$.
   * It is useful to characterize the shape of the ellipsoid by
   * \f[
   * \begin{align}
   *   e  &= \frac{\sqrt{a^2-c^2}}b,\\
   *   k  &= \frac{\sqrt{b^2-c^2}}{\sqrt{a^2-c^2}},\\
   *   k' &= \frac{\sqrt{a^2-b^2}}{\sqrt{a^2-c^2}};
   * \end{align}
   * \f]
   * note that \f$k^2 + k'^2 = 1\f$.  The spherical limit \f$ e\rightarrow 0
   * \f$ is nonuniform since the values of \f$k\f$ and \f$k'\f$ depend on how
   * the limit is taken.  In this cases, it's convenient to specify the
   * ellipsoid in terms of these parameters.  The semiaxes are related to these
   * parameters by
   * \f[
   * [a,b,c] = b \bigl[ \sqrt{1 + e^2k'^2}, 1, \sqrt{1 - e^2k^2} \bigr].
   * \f]
   *
   * Positions on the ellipsoid are given in term so the ellipsoidal
   * coordinates given in Section 2 of
   * - C. F. F. Karney,<br>
   *   <a href="https://arxiv.org/abs/2511.01621">
   *   Jacobi's solution for geodesics on a triaxial ellipsoid</a>,<br>
   *   Technical Report, SRI International, Nov. 2025.<br>
   *   <a href="https://arxiv.org/abs/2511.01621">arxiv:2511.01621</a>
   * .
   * Ellipsoidal latitude * \f$\beta\f$ and the ellipsoidal longitude
   * \f$\omega\f$ which are defined * by
   * \f[
   * \mathbf R = \begin{bmatrix}
   * a \cos\omega \sqrt{k^2\cos^2\beta + k'^2} \\
   * b \cos\beta \sin\omega \\
   * c \sin\beta \sqrt{k^2 + k'^2\sin^2\omega}
   * \end{bmatrix}.
   * \f]
   * Headings are given by the direction \f$ \alpha \f$ measured clockwise from
   * a line of constant \f$ \omega \f$.  Conversions between cartesian and
   * ellipsoidal coordinates is provided by cart2toellip() and elliptocart2().
   *
   * The ellipsoid coordinates "cover" the ellipsoid twice; the replacement
   * \f[
   * \begin{align}
   * \omega & \rightarrow -\omega,\\
   * \beta & \rightarrow \pi-\beta,\\
   * \alpha & \rightarrow \pi+\alpha,
   * \end{align}
   * \f]
   * leaves the position and direction unchanged; see AngNorm(),
   *
   * Example of use:
   * \include example-Ellipsoid3.cpp
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT Ellipsoid3 {
  public:
    /**
     * A type to hold three-dimensional positions and directions in cartesian
     * coordinates.
     **********************************************************************/
    using vec3 = std::array<Math::real, 3>;
  private:
    /// \cond SKIP
    friend class Cartesian3;  // For access to cart2toellipint normvec
    friend class Geodesic3;   // For Flip
    /// \endcond
    using real = Math::real;
    using ang = Angle;
    static void normvec(vec3& R) {
      real h = Math::hypot3(R[0], R[1], R[2]);
      // No checking for h = 0.  Result will be NaNs (we rely on this in
      // Cartesian3::cart2rand).
      R[0] /= h; R[1] /= h; R[2] /= h;
    }
    static void Flip(ang& bet, ang& omg, ang& alp) {
      bet.reflect(false, true);
      omg.reflect(true);
      alp.reflect(true, true);
    }
    real _a, _b, _c;            // semiaxes
    real _e2, _k2, _kp2, _k, _kp;
    bool _oblate, _prolate, _biaxial;
    void cart2toellipint(vec3 R, ang& bet, ang& omg, vec3 axes) const;
    /**
     * @return \e k the oblateness parameter.
     **********************************************************************/
    real k() const { return _k; }
    /**
     * @return \e kp the prolateness parameter.
     **********************************************************************/
    real kp() const { return _kp; }
    /**
     * @return whether the ellipsoid is oblate.
     *
     * This is determined by the condition kp2() == 0.
     **********************************************************************/
    bool oblate() const { return _oblate; }
    /**
     * @return whether the ellipsoid is prolate.
     *
     * This is determined by the condition k2() == 0.
     **********************************************************************/
    bool prolate() const { return _prolate; }
    /**
     * @return whether the ellipsoid is oblate or prolate.
     **********************************************************************/
    bool biaxial() const { return _biaxial; }
  public:
    /**
     * The default constructor for a unit sphere in the oblate limit.
     **********************************************************************/
    Ellipsoid3();
    /**
     * An ellipsoid specified by its semiaxes.
     *
     * @param[in] a the major semiaxis.
     * @param[in] b the median semiaxis.
     * @param[in] c the minor semiaxis.
     * @exception GeographicErr if the required ordering is semiaxes is
     *   violated.
     *
     * The semiaxes must satisfy \e a &ge; \e b &ge; \e c &gt; 0.
     * If \e a = \e c (a sphere), then the oblate limit is taken.
     **********************************************************************/
    Ellipsoid3(real a, real b, real c);
    /**
     * An ellipsoid specified by its median semiaxis and shape.
     *
     * @param[in] b the middle semiaxis.
     * @param[in] e2 the eccentricity squared \f$e^2 = (a^2 - c^2)/b^2\f$.
     * @param[in] k2 the oblateness parameter squared \f$k^2 = (b^2 - c^2) /
     *  (a^2 - c^2)\f$.
     * @param[in] kp2 the prolateness parameter squared \f$k'^2= (a^2 - b^2) /
     *   (a^2 - c^2)\f$.
     * @exception GeographicErr if the required ordering is semiaxes is
     *   violated.
     *
     * This form of the constructor is important when the eccentricity is small
     * and giving \e e2 allows for more precision.
     *
     * In the case of a sphere with \e e2 = 0, this constructor distinguishes
     * between and "oblate sphere" (\e k2 = 1), a "prolate sphere" (\e k2 = 0),
     * and a "triaxial sphere" (\e k2 &isin (0,1)).  These distinctions matter
     * when ellipsoidal coordinates are used.
     *
     * \note The constructor normalizes \e k2 and \e kp2 so that \e k2 + \e kp2
     *   = 1.
     **********************************************************************/
    Ellipsoid3(real b, real e2, real k2, real kp2);
    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the major semiaxis.
     **********************************************************************/
    real a() const { return _a; }
    /**
     * @return \e b the median semiaxis.
     **********************************************************************/
    real b() const { return _b; }
    /**
     * @return \e c the minor semiaxis.
     **********************************************************************/
    real c() const { return _c; }
    /**
     * @return \e e2 the eccentricity squared.
     **********************************************************************/
    real e2() const { return _e2; }
    /**
     * @return \e k2 the oblateness parameter squared.
     **********************************************************************/
    real k2() const { return _k2; }
    /**
     * @return \e kp2 the prolateness parameter squared.
     **********************************************************************/
    real kp2() const { return _kp2; }
    ///@}
    /** \name Normalizing functions
     **********************************************************************/
    ///@{
    /**
     * Scale a position to ensure it lies on the ellipsoid
     *
     * @param[inout] R the position.
     *
     * The components of \e R are scaled so that it lies on the ellipsoid.  The
     * resulting position is not in general the closest point on the ellipsoid.
     * Use Cartesian3::carttocart2() for that.
     **********************************************************************/
    void Norm(vec3& R) const;
    /**
     * Scale a position and direction to the ellipsoid
     *
     * @param[inout] R the position.
     * @param[inout] V the position.
     *
     * The components of \e R are scaled so that it lies on the ellipsoid.  Then
     * \e V is projected to be tangent to the surface and is normalized to a
     * unit vector.
     **********************************************************************/
    void Norm(vec3& R, vec3& V) const;
    /**
     * Set the sheet for coordinates.
     *
     * @param[inout] bet the ellipsoidal latitude.
     * @param[inout] omg the ellipsoidal longitude.
     * @param[inout] alp the heading.
     * @param[in] alt if true switch to the alternate sheet.
     * @return whether the coordinates were changed.
     *
     * If alt = false (the default), the conventional sheet is used switching
     * the values of \e bet, \e omg, and \e alp, so that \e bet &isin;
     * [-&pi;/2, &pi/2].
     *
     * If alt = true, the alternate sheet is used switching the values of \e
     * bet, \e omg, and \e alp, so that \e omg &isin; [0, &pi;].
     *
     * This routine does not change \e n (the number of turns) for the
     * coordinates.
     **********************************************************************/
    static bool AngNorm(Angle& bet, Angle& omg, Angle& alp, bool alt = false) {
      using std::signbit;
      // If !alt, put bet in [-pi/2,pi/2]
      // If  alt, put omg in [0, pi]
      bool flip = alt ? signbit(omg.s()) : signbit(bet.c());
      if (flip)
        Flip(bet, omg, alp);
      return flip;
    }
    /**
     * Set the sheet for coordinates.
     *
     * @param[inout] bet the ellipsoidal latitude.
     * @param[inout] omg the ellipsoidal longitude.
     * @param[in] alt if true switch to the alternate sheet.
     * @return whether the coordinated were changed.
     *
     * This acts precisely the same and AngNorm(Angle&, Angle&, Angle&, bool)
     * except that \e alp is omitted.
     **********************************************************************/
    static bool AngNorm(Angle& bet, Angle& omg, bool alt = false) {
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
    ///@}
    /** \name Coordinate conversions.
     **********************************************************************/
    ///@{
    /**
     * Convert a cartesian position to ellipsoidal coordinates.
     *
     * @param[in] R the cartesian position.
     * @param[out] bet the ellipsoidal latitude.
     * @param[out] omg the ellipsoidal longitude.
     *
     * \note \e R must lie on the surface of the ellipsoid.  The "2" in "cart2"
     * is used to emphasize this.
     **********************************************************************/
    void cart2toellip(vec3 R, Angle& bet, Angle& omg) const;
    /**
     * Convert a cartesian position and direction to ellipsoidal coordinates.
     *
     * @param[in] R the cartesian position.
     * @param[in] V the cartesian direction.
     * @param[out] bet the ellipsoidal latitude.
     * @param[out] omg the ellipsoidal longitude.
     * @param[out] alp the azimuth.
     *
     * \note \e R must lie on the surface of the ellipsoid and \e V must be
     *   tangent to the surface at that point.  The "2" in "cart2" is used to
     *   emphasize this.
     **********************************************************************/
    void cart2toellip(vec3 R, vec3 V,
                      Angle& bet, Angle& omg, Angle& alp) const;
    /**
     * Convert an ellipsoid position and cartesian direction to a heading.
     *
     * @param[in] bet the ellipsoidal latitude.
     * @param[in] omg the ellipsoidal longitude.
     * @param[in] V the cartesian direction.
     * @param[out] alp the azimuth.
     *
     * This is a variant of cart2toellip(vec3, vec3, Angle&, Angle&, Angle&)
     * where \e bet and \e omg are used to ensure that the correct sheet is
     * used in determining \e alp.
     *
     * \note \e V must be tangent to the surface of the ellipsoid.  The "2" in
     *   "cart2" is used to emphasize this.
     **********************************************************************/
    void cart2toellip(Angle bet, Angle omg,
                      vec3 V, Angle& alp) const;
    /**
     * Convert ellipsoidal coordinates to a cartesian position.
     *
     * @param[in] bet the ellipsoidal latitude.
     * @param[in] omg the ellipsoidal longitude.
     * @param[out] R the cartesian position.
     **********************************************************************/
    void elliptocart2(Angle bet, Angle omg, vec3& R) const;
    /**
     * Convert coordinates and heading to a cartesian position and direction.
     *
     * @param[in] bet the ellipsoidal latitude.
     * @param[in] omg the ellipsoidal longitude.
     * @param[in] alp the azimuth.
     * @param[out] R the cartesian position.
     * @param[out] V the cartesian direction.
     **********************************************************************/
    void elliptocart2(Angle bet, Angle omg, Angle alp,
                      vec3& R, vec3& V) const;
    ///@}
    /**
     * A global instantiation of Ellipsoid3 with the parameters for the
     * Earth.
     **********************************************************************/
    static const Ellipsoid3& Earth();

  };

  } // namespace Triaxial
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIAL_HPP

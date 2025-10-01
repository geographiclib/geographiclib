/**
 * \file Conformal3.hpp
 * \brief Header for GeographicLib::Triaxial::Conformal3 class
 *
 * Copyright (c) Charles Karney (2014-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CONFORMAL3_HPP)
#define GEOGRAPHICLIB_CONFORMAL3_HPP 1

#include <GeographicLib/Angle.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

namespace GeographicLib {
  namespace Triaxial {

  /**
   * \brief Jacobi's conformal projection of a triaxial ellipsoid
   *
   * This is a conformal projection of the ellipsoid to a plane in which the
   * grid lines are straight; see Jacobi,
   * <a href="https://books.google.com/books?id=ryEOAAAAQAAJ&pg=PA212">
   * Vorlesungen &uuml;ber Dynamik, &sect;28</a>.  The constructor takes the
   * semiaxes of the ellipsoid (which must be in order).  Member functions map
   * the ellipsoidal coordinates &omega; and &beta; separately to \e x and \e
   * y.  Jacobi's coordinates have been multiplied by \f$eb/2\f$ so that the
   * customary results are returned in the cases of a sphere or an ellipsoid of
   * revolution.  In particular the scale satisfies \f$m\ge 1\f$ with \f$m =
   * 1\f$ at \f$\cos\beta = \sin\omega = 1\f$.  (Note: usually %GeographicLib
   * uses \f$k\f$ to denote the scale.  However, in order to avoid confusion
   * with the oblateness parameter \f$k\f$ used to specify the triaxial
   * ellipsoid, the symbol \f$m\f$ is used for the scale in this class.)
   *
   * The ellipsoid is oriented so that the major principal ellipse, \f$Z=0\f$,
   * is the equator, \f$\beta=0\f$, while the minor principal ellipse,
   * \f$X=0\f$, is the meridian, \f$\omega\pm\pi/2\f$.  The four umbilical
   * points, \f$\sin\omega = \cos\beta = 0\f$, lie on median principal ellipse
   * in the plane \f$Y=0\f$.
   *
   * In this projection the easting, \e x, depends only on the longitude
   * \f$\omega\f$ and the northing, \e y depends only on the latitude
   * \f$\beta\f$.  Thus lines of constant latitude and longitude map to
   * straight lines.  For a general ellipsoid, each octant of the ellipsoid
   * maps to a finite rectangle of dimensions x0() x y0().  For an oblate
   * (resp. prolate) ellipsoid, x0() (resp. y0()) diverges.
   *
   * For any ellipsoid, we define a "conformal sphere" which has exactly the
   * same extent for the mapping of an octant.  This allows points on the
   * ellipsoid to be conformally mapped to a sphere.  This in turn allows any
   * triaxial ellipsoid to be conformally mapped to any other triaxial
   * ellipsoid.
   *
   * The units for the easting and the northing are the same as the units uses
   * to specify the size of the ellipsoid.
   *
   * For more information on this projection, see \ref jacobi.
   *
   * Example of use:
   * \include example-Conformal3.cpp
   *
   * <a href="Conformal3Proj.1.html">Conformal3Proj</a> is a command-line
   * utility providing access to the functionality of Conformal3.
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT Conformal3 {
  private:
    using real = Math::real;
    using ang = Angle;
    Ellipsoid3 _t, _s;
    EllipticFunction _ex, _ey, _exs, _eys;
    real _x, _y;

    real a() const { return _t.a(); }
    real b() const { return _t.b(); }
    real c() const { return _t.c(); }
    real e2() const { return _t.e2(); }
    real k2() const { return _t.k2(); }
    real kp2() const { return _t.kp2(); }
    static real Pi(const EllipticFunction& ell, ang phi);
    static ang Piinv(const EllipticFunction& ell, real x);
    static real F(const EllipticFunction& ell, ang phi);
    static ang Finv(const EllipticFunction& ell, real x);
    static ang omegashift(ang omg, int dir) {
      return omg - ang::cardinal(dir);
    }
    // Jacobi has m = 2/sqrt(lambda2 - lambda3)
    // Setting lambda2 = Cayley's k = -(b^2 \sin^2\beta + c^2 \cos^2\beta)
    //                              = -(b^2 - (b^2 - c^2) \cos^2\beta)
    // Setting lambda3 = Cayley's h = -(b^2 \cos^2\omega + a^2 \sin^2\omega)
    //                              = -(b^2 + (a^2 - b^2) \sin^2\omega)
    // lambda2 - lambda3 =
    //      - (b^2 - (b^2 - c^2) \cos^2\beta)
    //      + (b^2 + (a^2 - b^2) \sin^2\omega)
    // (a^2 - c^2) * (k'^2 * \sin^2\omega + k^2 \cos^2\beta)
    // m = 2/sqrt(a^2 - c^2) * 1/sqrt(k'^2 * \sin^2\omega + k^2 \cos^2\beta)
    real invscale(ang bet, ang omg) const {
      using std::sqrt;
      omg = omegashift(omg, +1);
      return sqrt(k2() * Math::sq(bet.c()) + kp2() * Math::sq(omg.c()));
    }
    static Ellipsoid3 EquivSphere(real x, real y,
                                  EllipticFunction& ellx,
                                  EllipticFunction& elly);
    real sphericalscale(real ma, real mb) const;
  public:
    using vec3 = Ellipsoid3::vec3;
    /**
     * Constructor for a triaxial ellipsoid defined by Ellipsoid3 object.
     *
     * @param[in] t the Ellipsoid3 object.
     **********************************************************************/
    Conformal3(const Ellipsoid3& t = Ellipsoid3{});
    /**
     * Constructor for a triaxial ellipsoid with semiaxes.
     *
     * @param[in] a the largest semiaxis.
     * @param[in] b the middle semiaxis.
     * @param[in] c the smallest semiaxis.
     * @exception GeographicErr if the required ordering is semiaxes is
     *   violated.
     *
     * The semiaxes must satisfy \e a &ge; \e b &ge; \e c &gt; 0.
     * If \e a = \e c (a sphere), then the oblate limit is taken.
     **********************************************************************/
    Conformal3(real a, real b, real c);
    /**
     * Alternate constructor for a triaxial ellipsoid.
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
     * \note The constructor normalizes \e k2 and \e kp2 to ensure then \e k2 +
     * \e kp2 = 1.
     **********************************************************************/
    Conformal3(real b, real e2, real k2, real kp2);
    /**
     * @return the quadrant length in the \e x direction.
     **********************************************************************/
    Math::real x0() const { return _x; }
    /**
     * The \e x projection.
     *
     * @param[in] omg the longitude &omega; as an Angle.
     * @return the easting.
     **********************************************************************/
    Math::real x(Angle omg) const;
    /**
     * The \e x projection.
     *
     * @param[in] omg the longitude &omega; (in degrees).
     * @return \e x.
     **********************************************************************/
    Math::real x(real omg) const {
      return x(ang(omg));
    }
    /**
     * The inverse \e x projection.
     *
     * @param[in] x the easting.
     * @return the longitude &omega; as an Angle.
     **********************************************************************/
    Angle omega(real x) const;
    /**
     * The inverse \e x projection.
     *
     * @param[in] x the easting.
     * @return the longitude &omega; in degrees.
     **********************************************************************/
    Math::real omegad(real x) const {
      return real(omega(x));
    }
    /**
     * @return the quadrant length in the \e y direction.
     **********************************************************************/
    Math::real y0() const { return _y; }
    /**
     * The \e y projection.
     *
     * @param[in] bet the latitude &beta; as an Angle
     * @return the northing.
     **********************************************************************/
    Math::real y(Angle bet) const;
    /**
     * The \e y projection.
     *
     * @param[in] bet the latitude &beta; (in degrees).
     * @return the northing.
     **********************************************************************/
    Math::real y(real bet) const {
      return y(ang(bet));
    }
    /**
     * The inverse \e y projection.
     *
     * @param[in] y the northing.
     * @return the latitude &beta; as an Angle.
     **********************************************************************/
    Angle beta(real y) const;
    /**
     * The inverse \e y projection.
     *
     * @param[in] y the northing.
     * @return the latitude &beta; in degrees.
     **********************************************************************/
    Math::real betad(real y) const {
      return real(beta(y));
    }
    /**
     * The forward projection.
     *
     * @param[in] bet the ellipsoidal latitude.
     * @param[in] omg the ellipsoidal longitude.
     * @param[out] x the easting.
     * @param[out] y the easting.
     **********************************************************************/
    void Forward(Angle bet, Angle omg, real& x, real& y) const;
    /**
     * The forward projection in degrees
     *
     * @param[in] bet the ellipsoidal latitude (in degrees).
     * @param[in] omg the ellipsoidal longitude (in degrees).
     * @param[out] x the easting.
     * @param[out] y the easting.
     **********************************************************************/
    void Forward(real bet, real omg, real& x, real& y) const {
      Forward(ang(bet), ang(omg), x, y);
    }
    /**
     * The forward projection with the scale.
     *
     * @param[in] bet the ellipsoidal latitude.
     * @param[in] omg the ellipsoidal longitude.
     * @param[out] x the easting.
     * @param[out] y the easting.
     * @param[out] m the scale.
     *
     * The meridian convergence for this projection is 0.
     **********************************************************************/
    void Forward(Angle bet, Angle omg, real& x, real& y, real& m) const;
    /**
     * The forward projection in degrees with the scale.
     *
     * @param[in] bet the ellipsoidal latitude (in degrees).
     * @param[in] omg the ellipsoidal longitude (in degrees).
     * @param[out] x the easting.
     * @param[out] y the easting.
     * @param[out] m the scale.
     *
     * The meridian convergence for this projection is 0.
     **********************************************************************/
    void Forward(real bet, real omg, real& x, real& y, real& m) const {
      Forward(ang(bet), ang(omg), x, y, m);
    }
    /**
     * The reverse projection.
     *
     * @param[in] x the easting.
     * @param[in] y the easting.
     * @param[out] bet the ellipsoidal latitude.
     * @param[out] omg the ellipsoidal longitude.
     **********************************************************************/
    void Reverse(real x, real y, Angle& bet, Angle& omg) const;
    /**
     * The reverse projection in degrees.
     *
     * @param[in] x the easting.
     * @param[in] y the easting.
     * @param[out] bet the ellipsoidal latitude (in degrees).
     * @param[out] omg the ellipsoidal longitude (in degrees).
     **********************************************************************/
    void Reverse(real x, real y, real& bet, real& omg) const {
      ang beta, omga;
      Reverse(x, y, beta, omga);
      bet = real(beta); omg = real(omga);
    }
    /**
     * The reverse projection with the scale.
     *
     * @param[in] x the easting.
     * @param[in] y the easting.
     * @param[out] bet the ellipsoidal latitude.
     * @param[out] omg the ellipsoidal longitude.
     * @param[out] m the scale.
     *
     * The meridian convergence for this projection is 0.
     **********************************************************************/
    void Reverse(real x, real y, Angle& bet, Angle& omg, real& m) const;
    /**
     * The reverse projection in degrees with the scale.
     *
     * @param[in] x the easting.
     * @param[in] y the easting.
     * @param[out] bet the ellipsoidal latitude (in degrees).
     * @param[out] omg the ellipsoidal longitude (in degrees).
     * @param[out] m the scale.
     *
     * The meridian convergence for this projection is 0.
     **********************************************************************/
    void Reverse(real x, real y, real& bet, real& omg, real& m) const {
      ang beta, omga;
      Reverse(x, y, beta, omga, m);
      bet = real(beta); omg = real(omga);
    }
    /**
     * The forward projection onto the conformal sphere.
     *
     * @param[in] bet the ellipsoidal latitude.
     * @param[in] omg the ellipsoidal longitude.
     * @param[out] r the projected position on the sphere
     * @param[out] v the projected direction for due north on the triaxial
     *   ellipsoid.
     * @param[out] m the scale.
     **********************************************************************/
    void ForwardSphere(Angle bet, Angle omg, vec3& r, vec3& v, real& m) const;
    /**
     * The forward projection in degrees onto the conformal sphere.
     *
     * @param[in] bet the ellipsoidal latitude (in degrees).
     * @param[in] omg the ellipsoidal longitude (in degrees).
     * @param[out] r the projected position on the sphere
     * @param[out] v the projected direction for due north on the triaxial
     *   ellipsoid.
     * @param[out] m the scale.
     **********************************************************************/
    void ForwardSphere(real bet, real omg, vec3& r, vec3& v, real& m) const {
      ForwardSphere(ang(bet), ang(omg), r, v, m);
    }
    /**
     * The reverse projection from the conformal sphere.
     *
     * @param[in] r the position on the sphere.
     * @param[in] v a reference direction.
     * @param[out] bet the ellipsoidal latitude.
     * @param[out] omg the ellipsoidal longitude.
     * @param[out] gam the meridian convergence.
     * @param[out] m the scale.
     *
     * Here \e v is the direction of grid north in some projection and \e gam
     * is the resulting meridian convergence.
     **********************************************************************/
    void ReverseSphere(vec3 r, vec3 v, Angle& bet, Angle& omg,
                       Angle& gam, real& m) const;
    /**
     * The reverse projection in degrees from the conformal sphere.
     *
     * @param[in] r the position on the sphere.
     * @param[in] v a reference direction.
     * @param[out] bet the ellipsoidal latitude (in degrees).
     * @param[out] omg the ellipsoidal longitude (in degrees).
     * @param[out] gam the meridian convergence (in degrees).
     * @param[out] m the scale.
     *
     * Here \e v is the direction of grid north in some projection and \e gam
     * is the resulting meridian convergence.
     **********************************************************************/
    void ReverseSphere(vec3 r, vec3 v, real& bet, real& omg,
                       real& gam, real& m) const {
      ang beta, omga, gama;
      ReverseSphere(r, v, beta, omga, gama, m);
      bet = real(beta); omg = real(omga); gam = real(gama);
    }
    /**
     * The forward conformal projection to some other ellipsoid.
     *
     * @param[in] alt the Conformal3 object for the other ellipsoid.
     * @param[in] bet the ellipsoidal latitude.
     * @param[in] omg the ellipsoidal longitude.
     * @param[out] betalt the ellipsoidal latitude on the other ellipsoid.
     * @param[out] omgalt the ellipsoidal longitude on the other ellipsoid.
     * @param[out] gam the meridian convergence.
     * @param[out] m the scale.
     **********************************************************************/
    void ForwardOther(const Conformal3& alt, Angle bet, Angle omg,
                      Angle& betalt, Angle& omgalt, Angle& gam, real& m)
      const;
    /**
     * The forward conformal projection in degrees to some other ellipsoid.
     *
     * @param[in] alt the Conformal3 object for the other ellipsoid.
     * @param[in] bet the ellipsoidal latitude (in degrees).
     * @param[in] omg the ellipsoidal longitude (in degrees).
     * @param[out] betalt the ellipsoidal latitude on the other ellipsoid (in
     *   degrees).
     * @param[out] omgalt the ellipsoidal longitude on the other ellipsoid (in
     *   degrees).
     * @param[out] gam the meridian convergence (in degrees).
     * @param[out] m the scale.
     **********************************************************************/
    void ForwardOther(const Conformal3& alt, real bet, real omg,
                      real& betalt, real& omgalt, real& gam, real& m)
      const {
      ang betalta, omgalta, gama;
      ForwardOther(alt, ang(bet), ang(omg), betalta, omgalta, gama, m);
      betalt = real(betalta); omgalt = real(omgalta); gam = real(gama);
    }
    /**
     * The reverse conformal projection from some other ellipsoid.
     *
     * @param[in] alt the Conformal3 object for the other ellipsoid.
     * @param[in] betalt the ellipsoidal latitude on the other ellipsoid.
     * @param[in] omgalt the ellipsoidal longitude on the other ellipsoid.
     * @param[out] bet the ellipsoidal latitude.
     * @param[out] omg the ellipsoidal longitude.
     * @param[out] gam the meridian convergence.
     * @param[out] m the scale.
     **********************************************************************/
    void ReverseOther(const Conformal3& alt, Angle betalt, Angle omgalt,
                      Angle& bet, Angle& omg, Angle& gam, real& m)
      const;
    /**
     * The reverse conformal projection in degrees from some other ellipsoid.
     *
     * @param[in] alt the Conformal3 object for the other ellipsoid.
     * @param[in] betalt the ellipsoidal latitude on the other ellipsoid (in
     *   degrees).
     * @param[in] omgalt the ellipsoidal longitude on the other ellipsoid (in
     *   degrees).
     * @param[out] bet the ellipsoidal latitude (in degrees).
     * @param[out] omg the ellipsoidal longitude (in degrees).
     * @param[out] gam the meridian convergence (in degrees).
     * @param[out] m the scale.
     **********************************************************************/
    void ReverseOther(const Conformal3& alt, real betalt, real omgalt,
                      real& bet, real& omg, real& gam, real& m)
      const {
      ang beta, omga, gama;
      ReverseOther(alt, ang(betalt), ang(omgalt), beta, omga, gama, m);
      bet = real(beta); omg = real(omga); gam = real(gama);
    }

    /**
     * @return the Ellipsoid3 object for this projection.
     **********************************************************************/
    const Ellipsoid3& t() { return _t; }
    /**
     * @return the Ellipsoid3 object for the conformal sphere.
     **********************************************************************/
    const Ellipsoid3& s() { return _s; }
  };

  } // namespace Triaxial
} // namespace GeographicLib
#endif  // GEOGRAPHICLIB_CONFORMAL3_HPP

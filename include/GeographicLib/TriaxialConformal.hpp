/**
 * \file TriaxialConformal.hpp
 * \brief Header for GeographicLib::TriaxialConformal class
 *
 * \note This is just sample code.  It is not part of GeographicLib
 * itself.
 *
 * Copyright (c) Charles Karney (2014-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIALCONFORMAL_HPP)
#define GEOGRAPHICLIB_TRIAXIALCONFORMAL_HPP 1

#include <GeographicLib/Angle.hpp>
#include <GeographicLib/Triaxial.hpp>
#include <GeographicLib/EllipticFunction.hpp>

namespace GeographicLib {
  /**
   * \brief Jacobi's conformal projection of a triaxial ellipsoid
   *
   * <b>NOTE:</b> This is just sample code.  It is not part of GeographicLib
   * itself.
   *
   * This is a conformal projection of the ellipsoid to a plane in which the
   * grid lines are straight; see Jacobi,
   *  <a href="https://books.google.com/books?id=ryEOAAAAQAAJ&pg=PA212">
   * Vorlesungen &uuml;ber Dynamik, &sect;28</a>.  The constructor takes the
   * semi-axes of the ellipsoid (which must be in order).  Member functions map
   * the ellipsoidal coordinates &omega; and &beta; separately to \e x and \e
   * y.  Jacobi's coordinates have been multiplied by \f$eb/2\f$ so that the
   * customary results are returned in the cases of a sphere or an ellipsoid of
   * revolution.  In particular the scale satisfies \f$k\ge 1\f$ with \f$k =
   * 1\f$ at \f$\cos\beta = \sin\omega = 1\f$.
   *
   * The ellipsoid is oriented so that the major principal ellipse, \f$Z=0\f$,
   * is the equator, \f$\beta=0\f$, while the minor principal ellipse,
   * \f$X=0\f$, is the meridian, \f$\omega\pm\pi/2\f$.  The four umbilic
   * points, \f$\sin\omega = \cos\beta = 0\f$, lie on median principal ellipse
   * in the plane \f$Y=0\f$.
   *
   * For more information on this projection, see \ref jacobi.
   *
   * Example of use:
   * \include example-TriaxialConformal.cpp
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT TriaxialConformal {
  private:
    typedef Math::real real;
    Triaxial _t;
    EllipticFunction _ex, _ey, _exalt, _eyalt;
    real a() const { return _t.a(); }
    real b() const { return _t.b(); }
    real c() const { return _t.c(); }
    real e2() const { return _t.e2(); }
    real k2() const { return _t.k2(); }
    real kp2() const { return _t.kp2(); }
    static real Pi(const EllipticFunction& ell, Angle phi);
    static Angle Piinv(const EllipticFunction& ell, real x);
    static Angle omegashift(Angle omg, int dir) {
      return omg - Angle::cardinal(dir);
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
    real scale(Angle bet, Angle omg) const {
      using std::sqrt;
      omg = omegashift(omg, +1);
      return 1/sqrt(k2() * Math::sq(bet.c()) + kp2() * Math::sq(omg.c()));
    }
  public:
    TriaxialConformal(const Triaxial& t);
    /**
     * Constructor for a trixial ellipsoid with semi-axes.
     *
     * @param[in] a the largest semi-axis.
     * @param[in] b the middle semi-axis.
     * @param[in] c the smallest semi-axis.
     *
     * The semi-axes must satisfy \e a &ge; \e b &ge; \e c > 0 and \e a >
     * \e c.
     **********************************************************************/
    TriaxialConformal(real a, real b, real c);
    /**
     * Alternate constructor for a triaxial ellipsoid.
     *
     * @param[in] b the middle semi-axis.
     * @param[in] e2 the smallest semi-axis.
     * @param[in] k2 the relative magnitude of \e a &minus; \e b.
     * @param[in] kp2 the relative magnitude of \e b &minus; \e c.
     *
     * This form can be used to specify a sphere.
     **********************************************************************/
    TriaxialConformal(real b, real e2, real k2, real kp2);
    /**
     * @return the quadrant length in the \e x direction (in meters).
     **********************************************************************/
    Math::real x() const;
    Math::real x2() const;
    /**
     * The \e x projection.
     *
     * @param[in] omg the longitude &omega; as an Angle.
     * @return the easting (in meters).
     **********************************************************************/
    Math::real x(Angle omg) const;
    Math::real x2(Angle omg) const;
    /**
     * The \e x projection.
     *
     * @param[in] omg the longitude &omega; (in degrees).
     * @return \e x (in meters).
     **********************************************************************/
    Math::real x(real omg) const {
      return x(Angle(omg));
    }
    /**
     * The inverse \e x projection.
     *
     * @param[in] x the easting (in meters).
     * @return the longitude &omega; as an Angle.
     **********************************************************************/
    Angle omega(real x) const;
    /**
     * The inverse \e x projection.
     *
     * @param[in] x the easting (in meters).
     * @return the longitude &omega; in degrees.
     **********************************************************************/
    Math::real omegad(real x) const {
      Angle omg = omega(x);
      return real(omg);
    }
    /**
     * @return the quadrant length in the \e y direction (in meters).
     **********************************************************************/
    Math::real y() const;
    Math::real y2() const;
    /**
     * The \e y projection.
     *
     * @param[in] bet the latitude &beta; as an Angle
     * @return the northing (in meters).
     **********************************************************************/
    Math::real y(Angle bet) const;
    Math::real y2(Angle bet) const;
    /**
     * The \e y projection.
     *
     * @param[in] bet the latitude &beta; (in degrees).
     * @return the northing (in meters).
     **********************************************************************/
    Math::real y(real bet) const {
      return y(Angle(bet));
    }
    /**
     * The inverse \e y projection.
     *
     * @param[in] y the northing (in meters).
     * @return the latitude &beta; as an Angle.
     **********************************************************************/
    Angle beta(real y) const;
    /**
     * The inverse \e y projection.
     *
     * @param[in] y the northing (in meters).
     * @return the latitude &beta; in degrees.
     **********************************************************************/
    Math::real betad(real y) const {
      Angle bet = beta(y);
      return real(bet);
    }
    void Forward(Angle bet, Angle omg, real& x, real& y) const;
    void Forward(real bet, real omg, real& x, real& y) const {
      Forward(Angle(bet), Angle(omg), x, y);
    }
    void Forward(Angle bet, Angle omg, real& x, real& y, real& m) const;
    void Forward(real bet, real omg, real& x, real& y, real& m) const {
      Forward(Angle(bet), Angle(omg), x, y, m);
    }
    void Reverse(real x, real y, Angle& bet, Angle& omg) const;
    void Reverse(real x, real y, real& bet, real& omg) const {
      Angle beta, omga;
      Reverse(x, y, beta, omga);
      bet = real(beta); omg = real(omga);
    }
    void Reverse(real x, real y, Angle& bet, Angle& omg, real& m) const;
    void Reverse(real x, real y, real& bet, real& omg, real& m) const {
      Angle beta, omga;
      Reverse(x, y, beta, omga, m);
      bet = real(beta); omg = real(omga);
    }
    const Triaxial& t() { return _t; }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALCONFORMAL_HPP

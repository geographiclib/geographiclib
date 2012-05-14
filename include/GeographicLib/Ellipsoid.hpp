/**
 * \file Ellipsoid.hpp
 * \brief Header for GeographicLib::Ellipsoid class
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ELLIPSOID_HPP)
#define GEOGRAPHICLIB_ELLIPSOID_HPP \
  "$Id$"

#include <string>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <GeographicLib/AlbersEqualArea.hpp>

namespace GeographicLib {

  /**
   * \brief Properties of an ellipsoid
   *
   * This class returns various properties of the ellipsoid and converts
   * between various types of latitudes.  The latitude conversions are also
   * possible using the various projections supported by %GeographicLib; but
   * Ellipsoid provides more direct access (sometimes using private functions
   * of the projection classes).  Ellipsoid::RectifyingLatitude,
   * Ellipsoid::InverseRectifyingLatitude, and Ellipsoid::MeridianDistance
   * provide functionality which can be provided by the Geodesic class.
   * However Geodesic uses a series approximation (valid for abs \e f < 1/150),
   * whereas Ellipsoid computes these quantities using EllipticFunction which
   * provides accurate results even when \e f is large.
   *
   * Example of use:
   * \include example-Ellipsoid.cpp
   **********************************************************************/

  class GEOGRAPHIC_EXPORT Ellipsoid {
  private:
    typedef Math::real real;
    static const int numit_ = 10;
    real _a, _f, _f1, _f12, _e2, _e12, _n, _b, _stol;
    const TransverseMercatorExact _tm;
    const AlbersEqualArea _au;
    static real tand(real x) throw() {
      return
        std::abs(x) == real(90) ? (x < 0 ?
                                   - TransverseMercatorExact::overflow_
                                   : TransverseMercatorExact::overflow_) :
        std::tan(x * Math::degree<real>());
    }
    static real atand(real x) throw()
    { return std::atan(x) / Math::degree<real>(); }

  public:

    /** \name Constructor
     **********************************************************************/
    ///@{

    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters).
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.  If \e f > 1, set flattening
     *   to 1/\e f.
     *
     * An exception is thrown if either of the axes of the ellipsoid is
     * non-positive.
     **********************************************************************/
    Ellipsoid(real a, real f);
    ///@}

    /** \name %Ellipsoid dimensions.
     **********************************************************************/
    ///@{

    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _a; }

    /**
     * @return \e b the polar semi-axis (meters).
     **********************************************************************/
    Math::real MinorRadius() const throw() { return _b; }

    /**
     * @return \e L the distance between the equator and a pole along a
     *   meridian (meters).  For a sphere \e L = (\e pi / 2) \e a.  The radius
     *   of a sphere with the same meridian length is \e L / (\e pi / 2).
     **********************************************************************/
    Math::real QuarterMeridian() const throw();

    /**
     * @return \e A the total area of the ellipsoid (meters<sup>2</sup>).  For
     *   a sphere \e A = 4\e pi <i>a</i><sup>2</sup>.  The radius of a sphere
     *   with the same area is sqrt(\e A / (4 \e pi)).
     **********************************************************************/
    Math::real Area() const throw();

    /**
     * @return \e V the total volume of the ellipsoid (meters<sup>3</sup>).
     *   For a sphere \e V = (4\e pi / 3) <i>a</i><sup>3</sup>.  The radius of
     *   a sphere with the same volume is cbrt(\e V / (4 \e pi / 3)).
     **********************************************************************/
    Math::real Volume() const throw()
    { return (4 * Math::pi<real>()) * Math::sq(_a) * _b / 3; }
    ///@}

    /** \name %Ellipsoid shape
     **********************************************************************/
    ///@{

    /**
     * @return \e f = (\e a - \e b) / \e a, the flattening of the ellipsoid.
     *   This is the value used in the constructor.  This is zero, positive, or
     *   negative for a sphere, oblate ellipsoid, or prolate ellipsoid.
     **********************************************************************/
    Math::real Flattening() { return _f; }

    /**
     * @return \e f' = (\e a - \e b) / \e b, the second flattening of the
     *   ellipsoid.  This is zero, positive, or negative for a sphere, oblate
     *   ellipsoid, or prolate ellipsoid.
     **********************************************************************/
    Math::real SecondFlattening() { return _f / (1 - _f); }

    /**
     * @return \e n = (\e a - \e b) / (\e a + \e b), the third flattening of
     *   the ellipsoid.  This is zero, positive, or negative for a sphere,
     *   oblate ellipsoid, or prolate ellipsoid.
     **********************************************************************/
    Math::real ThirdFlattening() { return _n; }

    /**
     * @return <i>e</i><sup>2</sup> = (<i>a</i><sup>2</sup> -
     *   <i>b</i><sup>2</sup>) / <i>a</i><sup>2</sup>, the eccentricity squared
     *   of the the ellipsoid.  This is zero, positive, or negative for a
     *   sphere, oblate ellipsoid, or prolate ellipsoid.
     **********************************************************************/
    Math::real EccentricitySq() { return _e2; }

    /**
     * @return <i>e'</i><sup>2</sup> = (<i>a</i><sup>2</sup> -
     *   <i>b</i><sup>2</sup>) / <i>b</i><sup>2</sup>, the second eccentricity
     *   squared of the the ellipsoid.  This is zero, positive, or negative for
     *   a sphere, oblate ellipsoid, or prolate ellipsoid.
     **********************************************************************/
    Math::real SecondEccentricitySq() { return _e12; }

    /**
     * @return <i>e''</i><sup>2</sup> = (<i>a</i><sup>2</sup> -
     *   <i>b</i><sup>2</sup>) / (<i>a</i><sup>2</sup> + <i>b</i><sup>2</sup>),
     *   the third eccentricity squared of the the ellipsoid.  This is zero,
     *   positive, or negative for a sphere, oblate ellipsoid, or prolate
     *   ellipsoid.
     **********************************************************************/
    Math::real ThirdEccentricitySq() { return _e2 / (2 - _e2); }
    ///@}

    /** \name Latitude conversion.
     **********************************************************************/
    ///@{

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e beta the parametric latitude (degrees).
     *
     * The geographic latitude, \e phi, is the angle beween the equatorial
     * plane and a vector normal to the surface of the ellipsoid.
     *
     * The parametric latitude (also called the reduced latitude), \e beta,
     * allows the cartesian coordinated of a meridian to be expressed
     * conveniently in parametric form as
     * - \e R = \e a cos \e beta
     * - \e Z = \e b sin \e beta
     * .
     * where \e a and \e b are the equatorial radius and the polar semi-axis.
     * For a sphere \e beta = \e phi.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e beta lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real ParametricLatitude(real phi) const throw();

    /**
     * @param[in] beta the parametric latitude (degrees).
     * @return \e phi the geographic latitude (degrees).
     *
     * \e beta must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e phi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real InverseParametricLatitude(real beta) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e theta the geocentric latitude (degrees).
     *
     * The geocentric latitude, \e theta, is the angle beween the equatorial
     * plane and a line between the center of the ellipsoid and a point on the
     * ellipsoid.  For a sphere \e theta = \e phi.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e theta lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real GeocentricLatitude(real phi) const throw();

    /**
     * @param[in] theta the geocentric latitude (degrees).
     * @return \e phi the geographic latitude (degrees).
     *
     * \e theta must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e phi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real InverseGeocentricLatitude(real theta) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e mu the rectifying latitude (degrees).
     *
     * The rectifying latitude, \e mu, has the property that the distance along
     * a meridian of the ellipsoid between two points with rectifying latitudes
     * <i>mu</i><sub>1</sub> and <i>mu</i><sub>2</sub> is equal to
     * (<i>mu</i><sub>2</sub> - <i>mu</i><sub>1</sub>) \e L / 90<sup>o</sup>,
     * where \e L = QuarterMeridian().  For a sphere \e mu = \e phi.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e mu lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real RectifyingLatitude(real phi) const throw();

    /**
     * @param[in] mu the rectifying latitude (degrees).
     * @return \e phi the geographic latitude (degrees).
     *
     * \e mu must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e phi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real InverseRectifyingLatitude(real mu) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e xi the authalic latitude (degrees).
     *
     * The authalic latitude, \e xi, has the property that the area of the
     * ellipsoid between two circles with authalic latitudes
     * <i>xi</i><sub>1</sub> and <i>xi</i><sub>2</sub> is equal to (sin
     * <i>xi</i><sub>2</sub> - sin <i>xi</i><sub>1</sub>) \e A / 2, where \e A
     * = Area().  For a sphere \e xi = \e phi.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e xi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real AuthalicLatitude(real phi) const throw();

    /**
     * @param[in] xi the authalic latitude (degrees).
     * @return \e phi the geographic latitude (degrees).
     *
     * \e xi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e phi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real InverseAuthalicLatitude(real xi) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e chi the conformal latitude (degrees).
     *
     * The conformal latitude, \e chi, gives the mapping of the ellipsoid to a
     * sphere which which is conformal (angles are preserved) and in which the
     * equator of the ellipsoid maps to the equator of the sphere.  For a
     * sphere \e chi = \e phi.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e chi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real ConformalLatitude(real phi) const throw();

    /**
     * @param[in] chi the conformal latitude (degrees).
     * @return \e phi the geographic latitude (degrees).
     *
     * \e chi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.  The returned value
     * \e phi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real InverseConformalLatitude(real chi) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e psi the isometric latitude (degrees).
     *
     * The isometric latitude gives the mapping of the ellipsoid to a plane
     * which which is conformal (angles are preserved) and in which the equator
     * of the ellipsoid maps to a straight line of constant scale; this mapping
     * defines the Mercator projection.  For a sphere \e psi =
     * sinh<sup>-1</sup> tan \e phi.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.
     **********************************************************************/
    Math::real IsometricLatitude(real phi) const throw();

    /**
     * @param[in] psi the isometric latitude (degrees).
     * @return \e phi the geographic latitude (degrees).
     *
     * The returned value \e phi lies in [-90<sup>o</sup>, 90<sup>o</sup>].
     **********************************************************************/
    Math::real InverseIsometricLatitude(real psi) const throw();
    ///@}

    /** \name Other quantities.
     **********************************************************************/
    ///@{

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e R = \e a cos \e beta the radius of a circle of latitude \e
     *   phi (meters).  \e R (\e pi / 180<sup>o</sup>) gives meters per degree
     *   longitude measured along a circle of latitude.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.
     **********************************************************************/
    Math::real CircleRadius(real phi) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e Z = \e b sin \e beta the distance of a circle of latitude \e
     *   phi from the equator measured parallel to the ellipsoid axis (meters).
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.
     **********************************************************************/
    Math::real CircleHeight(real phi) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e s the distance along a meridian
     *   between the equator and a point of latitude \e phi (meters).  \e s is
     *   given by \e s = \e mu \e L / 90<sup>o</sup>, where \e L =
     *   QuarterMeridian()).
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.
     **********************************************************************/
    Math::real MeridianDistance(real phi) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e rho the meridional radius of curvature of the ellipsoid at
     *   latitude \e phi (meters); this is the curvature of the meridian.  \e
     *   rho is given by \e rho = (180<sup>o</sup> / \e pi) d\e s / d\e phi,
     *   where \e s = MeridianDistance(); thus \e rho (\e pi / 180<sup>o</sup>)
     *   gives meters per degree latitude measured along a meridian.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.
     **********************************************************************/
    Math::real MeridionalCurvatureRadius(real phi) const throw();

    /**
     * @param[in] phi the geographic latitude (degrees).
     * @return \e nu the transverse radius of curvature of the ellipsoid at
     *   latitude \e phi (meters); this is the curvature of a curve on the
     *   ellipsoid which also lies in a plane perpendicular to the ellipsoid
     *   and to the meridian.  \e nu is related to \e R = CircleRadius() by \e
     *   R = \e nu cos \e phi.
     *
     * \e phi must lie in the range [-90<sup>o</sup>, 90<sup>o</sup>]; the
     * result is undefined if this condition does not hold.
     **********************************************************************/
    Math::real TransverseCurvatureRadius(real phi) const throw();
    ///@}

    /**
     * A global instantiation of Ellipsoid with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    static const Ellipsoid WGS84;

  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ELLIPSOID_HPP

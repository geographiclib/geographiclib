/**
 * \file Gnomonic.hpp
 * \brief Header for GeographicLib::Gnomonic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GNOMONIC_HPP)
#define GEOGRAPHICLIB_GNOMONIC_HPP "$Id$"

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/Constants.hpp"

namespace GeographicLib {

  /**
   * \brief %Gnomonic Projection.
   *
   * %Gnomonic projection centered at an arbitrary position \e C on the
   * ellipsoid.  The projection of \e P is defined as follows: compute the
   * geodesic line from \e C to \e P; compute the reduced length \e m12,
   * geodesic scale \e M12, and \e rho = \e m12/\e M12; finally \e x = \e rho
   * sin \e azi1; \e y = \e rho cos \e azi1, where \e azi1 is the azimuth of
   * the geodesic at \e C.  The Forward and Reverse methods also return the
   * azimuth \e azi of the geodesic at \e P and reciprocal scale \e rk in the
   * azimuthal direction.  The scale in the radial direction if 1/\e
   * rk<sup>2</sup>.
   *
   * For a sphere, \e rho is reduces to \e a tan(\e s12/\e a), where \e s12 is
   * the length of the geodesic from \e C to \e P, and the gnomonic projection
   * has the property that all geodesics appear as straight lines.  For an
   * ellipsoid, this property holds only for geodesics interesting the centers.
   * However geodesics segments close to the center are approximately straight;
   * the deviation from straightness is of order \e f (\e r/\e a)<sup>3</sup>
   * where \e a and \e f are the major radius and the flattening of the
   * ellipsoid and \e r is the maximum distance of the geodesic from the
   * center.
   *
   * The conversions all take place using a GeographicLib::Geodesic object (by
   * default GeographicLib::Geodesic::WGS84).  For more information on
   * geodesics see \ref geodesic.
   *
   * <b>CAUTION:</b> The definition of this projection for a sphere is
   * standard.  However, there is no standard for how it should be extended to
   * an ellipsoid.  The choices are:
   * - Declare that the projection is undefined for an ellipsoid.
   * - Project to a tangent plane from the center of the ellipsoid.  This
   *   causes great ellipses to appear as straight lines in the projection;
   *   i.e., it generalizes the spherical great circle to a great ellipse.
   *   This was proposed by Roy Williams, Geometry of Navigation (Horwood,
   *   Chichester, 1998).
   * - Project to the conformal sphere with the constant of integration chosen
   *   so that the values of the latitude match for the center point and
   *   perform a central projection onto the plane tangent to the conformal
   *   sphere at the center point.  This causes normal sections through the
   *   center point to appear as straight lines in the projection; i.e., it
   *   generalizes the spherical great circle to a normal section.  This was
   *   proposed by I. G. Letoval'tsev, Generalization of the %Gnomonic
   *   Projection for a Spheroid and the Principal Geodetic Problems Involved
   *   in the Alignment of Surface Routes, Geodesy and Aerophotography (5),
   *   271-274 (1963).
   * - The projection given here.  This causes geodesics close to the center
   *   point to appear as straight lines in the projection; i.e., it
   *   generalizes the spherical great circle to a geodesic.
   **********************************************************************/

  class Gnomonic {
  private:
    typedef Math::real real;
    const Geodesic _earth;
    real _a, _f;
    static const real eps0, eps;
    static const int numit = 10;
  public:

    /**
     * Constructor for Gnomonic setting the Geodesic object to use
     * for geodesic calculations.  By default this uses the WGS84 ellipsoid.
     **********************************************************************/
    explicit Gnomonic(const Geodesic& earth = Geodesic::WGS84)
      throw()
      : _earth(earth)
      , _a(_earth.MajorRadius())
      , _f(_earth.InverseFlattening() ?
           1/std::abs(_earth.InverseFlattening()) : 0)
    {}

    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * gnomonic easting \e x (meters) and northing \e y (meters).  The center
     * of the projection is at latitude \e lat0 (degrees) and longitude \e lon0
     * (degrees).  Also return the azimuth \e azi (degrees) and the reciprocal
     * of the azimuthal scale \e rk.  \e lat0 and \e lat should be in the range
     * [-90, 90] and \e lon0 and \e lon should be in the range [-180, 360].
     * The scale of the projection is 1/\e rk<sup>2</sup> in the "radial"
     * direction, \e azi clockwise from true north, and is 1/\e rk in the
     * direction perpendicular to this.  If the point lies "over the horizon",
     * i.e., if \e rk <= 0, then NaNs are returned for \e x and \e y (the
     * correct values are returned for \e azi and \e rk).  A call to Forward
     * followed by a call to Reverse will return the original (\e lat, \e lon)
     * (to within roundoff) provided the point in not over the horizon.
     **********************************************************************/
    void Forward(real lat0, real lon0, real lat, real lon,
                 real& x, real& y, real& azi, real& rk) const throw();

    /**
     * Convert from gnomonic easting \e x (meters) and northing \e y (meters)
     * to latitude \e lat (degrees) and longitude \e lon (degrees).  The center
     * of the projection is at latitude \e lat0 (degrees) and longitude \e lon0
     * (degrees).  Also return the azimuth \e azi (degrees) and the reciprocal
     * of the azimuthal scale \e rk.  \e lat0 should be in the range [-90, 90]
     * and \e lon0 should be in the range [-180, 360].  \e lat will be in the
     * range [-90, 90] and \e lon will be in the range [-180, 180).  The scale
     * of the projection is 1/\e rk<sup>2</sup> in the "radial" direction, \e
     * azi clockwise from true north, and is 1/\e rk in the direction
     * perpendicular to this.  Even though all inputs should return a valid \e
     * lat and \e lon, it's possible that the procedure fails to converge for
     * very large \e x or \e y; in this case NaNs are returned for all the
     * output arguments.  A call to Reverse followed by a call to Forward will
     * return the original (\e x, \e y) (to roundoff).
     **********************************************************************/
    void Reverse(real lat0, real lon0, real x, real y,
                 real& lat, real& lon, real& azi, real& rk) const throw();

    /**
     * The major radius of the ellipsoid (meters).  This is that value of \e a
     * inherited from the Geodesic object used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _earth.MajorRadius(); }

    /**
     * The inverse flattening of the ellipsoid.  This is that value of \e r
     * inherited from the Geodesic object used in the constructor.  A value of
     * 0 is returned for a sphere (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw()
    { return _earth.InverseFlattening(); }
  };

} // namespace GeographicLib

#endif

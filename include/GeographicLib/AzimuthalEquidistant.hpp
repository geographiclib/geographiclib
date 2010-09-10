/**
 * \file AzimuthalEquidistant.hpp
 * \brief Header for GeographicLib::AzimuthalEquidistant class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_AZIMUTHALEQUIDISTANT_HPP)
#define GEOGRAPHICLIB_AZIMUTHALEQUIDISTANT_HPP "$Id$"

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/Constants.hpp"

namespace GeographicLib {

  /**
   * \brief Azimuthal Equidistant Projection.
   *
   * Azimuthal equidistant projection centered at an arbitrary position on the
   * ellipsoid.  For a point in projected space (\e x, \e y), the geodesic
   * distance from the center position is hypot(\e x, \e y) and the azimuth of
   * the geodesic from the center point is atan2(\e x, \e y).  The Forward and
   * Reverse methods also return the azimuth \e azi of the geodesic at (\e x,
   * \e y) and reciprocal scale \e rk in the azimuthal direction which,
   * together with the basic properties of the projection, serve to specify
   * completely the local affine transformation between geographic and
   * projected coordinates.
   *
   * The conversions all take place using a GeographicLib::Geodesic object (by
   * default GeographicLib::Geodesic::WGS84).  For more information on
   * geodesics see \ref geodesic.
   **********************************************************************/

  class AzimuthalEquidistant {
  private:
    typedef Math::real real;
    const Geodesic _earth;
    static const real eps;
  public:

    /**
     * Constructor for AzimuthalEquidistant setting the Geodesic object to use
     * for geodesic calculations.  By default this uses the WGS84 ellipsoid.
     **********************************************************************/
    explicit AzimuthalEquidistant(const Geodesic& earth = Geodesic::WGS84)
      throw() : _earth(earth) {}

    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * azimuthal equidistant easting \e x (meters) and northing \e y (meters).
     * The center of the projection is at latitude \e lat0 (degrees) and
     * longitude \e lon0 (degrees).  Also return the azimuth \e azi (degrees)
     * and the reciprocal of the azimuthal scale \e rk.  \e lat0 and \e lat
     * should be in the range [-90, 90] and \e lon0 and \e lon should be in the
     * range [-180, 360].  The scale of the projection is 1 in the "radial"
     * direction, \e azi clockwise from true north, and is 1/\e rk in the
     * direction perpendicular to this.  A call to Forward followed by a call
     * to Reverse will return the original (\e lat, \e lon) (to within
     * roundoff).
     **********************************************************************/
    void Forward(real lat0, real lon0, real lat, real lon,
                 real& x, real& y, real& azi, real& rk) const throw();

    /**
     * Convert from azimuthal equidistant easting \e x (meters) and northing \e
     * y (meters) to latitude \e lat (degrees) and longitude \e lon (degrees).
     * The center of the projection is at latitude \e lat0 (degrees) and
     * longitude \e lon0 (degrees).  Also return the azimuth \e azi (degrees)
     * and the reciprocal of the azimuthal scale \e rk.  \e lat0 should be in
     * the range [-90, 90] and \e lon0 should be in the range [-180, 360].  \e
     * lat will be in the range [-90, 90] and \e lon will be in the range
     * [-180, 180).  The scale of the projection is 1 in the "radial"
     * direction, \e azi clockwise from true north, and is 1/\e rk in the
     * direction perpendicular to this.  A call to Reverse followed by a call
     * to Forward will return the original (\e x, \e y) (to roundoff) only if
     * the geodesic to (\e x, \e y) is a shortest path.
     **********************************************************************/
    void Reverse(real lat0, real lon0, real x, real y,
                 real& lat, real& lon, real& azi, real& rk) const throw();

    /**
     * AzimuthalEquidistant::Forward without returning the azimuth and scale.
     **********************************************************************/
    void Forward(real lat0, real lon0, real lat, real lon,
                 real& x, real& y) const throw() {
      real azi, rk;
      Forward(lat0, lon0, lat, lon, x, y, azi, rk);
    }

    /**
     * AzimuthalEquidistant::Reverse without returning the azimuth and scale.
     **********************************************************************/
    void Reverse(real lat0, real lon0, real x, real y,
                 real& lat, real& lon) const throw() {
      real azi, rk;
      Reverse(lat0, lon0, x, y, lat, lon, azi, rk);
    }

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

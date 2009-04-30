/**
 * \file AzimuthalEquidistant.hpp
 * \brief Header for GeographicLib::AzimuthalEquidistant class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(AZIMUTHALEQUIDISTANT_HPP)
#define AZIMUTHALEQUIDISTANT_HPP "$Id$"

#include "GeographicLib/Geodesic.hpp"

namespace GeographicLib {

  /**
   * \brief Azimuthal Equidistant Projection.
   *
   * Azimuthal equidistant projection centered at an arbitrary position on the
   * ellipsoid.  For a point in projected space (\e x, \e y), the geodesic
   * distance from the center position is hypot(\e x, \e y) and the azimuth of
   * the geodesic from the center point is atan2(\e x, \e y).  The Forward and
   * Reverse methods also return the azimuth \e azi of the geodesic at (\e x,
   * \e y) and reduced length \e m of the geodesic which, together with the
   * basic properties of the projection, serve to specify completely the local
   * affine transformation between geographic and projected coordinates.
   *
   * The conversions all take place using a GeographicLib::Geodesic object (by
   * default GeographicLib::Geodesic::WGS84).  For more information on
   * geodesics see \ref geodesic.
   **********************************************************************/

  class AzimuthalEquidistant {
  private:
    const Geodesic& _earth;
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) throw()
    { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) throw()
    { return ::hypot(x, y); }
#endif
  public:

    /**
     * Constructor for AzimuthalEquidistant setting the Geodesic object to use
     * for geodesic calculations.  By default this uses the WGS84 ellipsoid.
     **********************************************************************/
    AzimuthalEquidistant(const Geodesic& earth = Geodesic::WGS84) throw()
      : _earth(earth) {}

    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * azimuthal equidistant easting \e x (meters) and northing \e y (meters).
     * The center of the projection is at latitude \e lat0 (degrees) and
     * longitude \e lon0 (degrees).  Also return the azimuth \e azi (degrees)
     * and the reduced length \e m (meters).  \e lat0 and \e lat should be in
     * the range [-90, 90] and \e lon0 and \e lon should be in the range [-180,
     * 360].  The scale of the projection is 1 in the "radial" direction, \e
     * azi clockwise from true north, and is hypot(\e x, \e y)/\e m in the
     * direction perpendicular to this.  A call to Forward followed by a call
     * to Reverse will return the original (\e lat, \e lon) (to within
     * roundoff).
     **********************************************************************/
    void Forward(double lat0, double lon0, double lat, double lon,
		 double& x, double& y, double& azi, double& m) const throw();

    /**
     * Convert from azimuthal equidistant easting \e x (meters) and northing \e
     * y (meters) to latitude \e lat (degrees) and longitude \e lon (degrees).
     * The center of the projection is at latitude \e lat0 (degrees) and
     * longitude \e lon0 (degrees).  Also return the azimuth \e azi (degrees)
     * and the reduced length \e m (meters).  \e lat0 should be in the range
     * [-90, 90] and \e lon0 should be in the range [-180, 360].  \e lat will
     * be in the range [-90, 90] and \e lon will be in the range [-180, 180).
     * The scale of the projection is 1 in the "radial" direction, \e azi
     * clockwise from true north, and is hypot(\e x, \e y)/\e m in the
     * direction perpendicular to this.  A call to Reverse followed by a call
     * to Forward will return the original (\e x, \e y) (to roundoff) only if
     * the geodesic to (\e x, \e y) is a shortest path.
     **********************************************************************/
    void Reverse(double lat0, double lon0, double x, double y,
		 double& lat, double& lon,
		 double& azi, double& m) const throw();

  };

} // namespace GeographicLib

#endif

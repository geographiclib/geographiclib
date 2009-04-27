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
   * ellipsoid.  For an point in projected space (\e x, \e y), the distance
   * from the center position is hypot(\e x, \e y) and the azimuth from the
   * center point is atan2(\e x, \e y).  The Forward and Reverse methods also
   * return the azimuth and reduced length which, together with the basic
   * properties of the projection, serve to specify completely the local affine
   * transformation between geographic and projected coordinates.
   *
   * The conversions all take place using a GeographicLib::Geodesic object (by
   * default GeographicLib::Geodesic::WGS84).
   **********************************************************************/

  class AzimuthalEquidistant {
  private:
    const Geodesic& _earth;
    double _lat0, _lon0;
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) throw()
    { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) throw()
    { return ::hypot(x, y); }
#endif
  public:

    /**
     * Constructor setting the origin to latitude = \e lat0 and longitude = \e
     * lon0 (degrees).  \e lat0 should be in the range [-90, 90] and \e lon0
     * should be in the range [-180, 360].  The optional \e earth argument
     * (default Geodesic::WGS84) specifies the Geodesic object to use for the
     * transformation.
     **********************************************************************/
    AzimuthalEquidistant(double lat0, double lon0,
			 const Geodesic& earth = Geodesic::WGS84) throw()
      : _earth(earth)
    { Reset(lat0, lon0); }

    /**
     * Default constructor sets the origin to \e lat0 = 0 and \e lon0 = 0.  The
     * optional \e earth argument (default Geodesic::WGS84) specifies the
     * Geodesic object to use for the transformation.
     **********************************************************************/
    AzimuthalEquidistant(const Geodesic& earth = Geodesic::WGS84) throw()
      : _earth(earth)
    { Reset(0.0, 0.0); }

    /**
     * Change the origin to latitude = \e lat0 and longitude = \e lon0
     * (degrees).  \e lat0 should be in the range [-90, 90] and \e lon0 should
     * be in the range [-180, 360].
     **********************************************************************/
    void Reset(double lat0, double lon0) throw();

    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * azimuthal equidistant easting \e x (meters) and northing \e y (meters).
     * Also return the azimuth \e azi (degrees) and the reduced length \e m
     * (meters).  \e lat should be in the range [-90, 90] and \e lon should be
     * in the range [-180, 360].  The scale of the projection is 1 in the
     * "radial" direction, \e azi clockwise from true north, and is hypot(\e x,
     * \e y)/\e m in the direction perpendicular to this.
     **********************************************************************/
    void Forward(double lat, double lon,
		 double& x, double& y,
		 double& azi, double& m) const throw();

    /**
     * Convert from azimuthal equidistant easting \e x (meters) and northing \e
     * y (meters) to latitude \e lat (degrees) and longitude \e lon (degrees).
     * Also return the azimuth \e azi (degrees) and the reduced length \e m
     * (meters).  \e lat will be in the range [-90, 90] and \e lon will be in
     * the range [-180, 180).  The scale of the projection is 1 in the "radial"
     * direction, \e azi clockwise from true north, and is hypot(\e x, \e y)/\e
     * m in the direction perpendicular to this.
     **********************************************************************/
    void Reverse(double x, double y,
		 double& lat, double& lon,
		 double& azi, double& m) const throw();

    /**
     * Return the latitude of the origin (degrees).
     **********************************************************************/
    double LatitudeOrigin() const throw() { return _lat0; }

    /**
     * Return the longitude of the origin (degrees).
     **********************************************************************/
    double LongitudeOrigin() const throw() { return _lon0; }

  };

} // namespace GeographicLib

#endif

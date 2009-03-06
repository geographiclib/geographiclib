/**
 * \file LocalCartesian.hpp
 * \brief Header for GeographicLib::LocalCartesian class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(LOCALCARTESIAN_HPP)
#define LOCALCARTESIAN_HPP "$Id$"

#include "GeographicLib/Geocentric.hpp"

namespace GeographicLib {

  /**
   * \brief Local Cartesian coordinates
   *
   * Convert between geodetic coordinates latitude = \e lat, longitude = \e
   * lon, height = \e h (measured vertically from the surface of the ellipsoid)
   * to local cartesian coordinates (\e x, \e y, \e z).  The origin of local
   * cartesian coordinate system is at \e lat = \e lat0, \e lon = \e lon0, \e h
   * = \e h0. The \e z axis is normal to the ellipsoid; the \e y axis points due
   * north.  The plane \e z = - \e h0 is tangent to the ellipsoid.
   *
   * The conversions all take place via geocentric coordinates using a
   * GeographicLib::GeoCentric object (by default
   * GeographicLib::Geocentric::WGS84).
   **********************************************************************/

  class LocalCartesian {
  private:
    const Geocentric& _earth;
    double _lat0, _lon0, _h0;
    double _x0, _y0, _z0,
      _rxx, _rxy, _rxz,
      _ryx, _ryy, _ryz,
      _rzx, _rzy, _rzz;
  public:

    /**
     * Constructor setting the origin to latitude = \e lat0, longitude = \e
     * lon0 (degrees), height = \e h0 (meters).  The optional \e earth argument
     * (default Geocentric::WGS84) specifies the Geocentric object to use for
     * the transformation.
     **********************************************************************/
    LocalCartesian(double lat0, double lon0, double h0 = 0,
		   const Geocentric& earth = Geocentric::WGS84) throw()
      : _earth(earth)
    { Reset(lat0, lon0, h0); }

    /**
     * Default constructor sets the origin to \e lat0 = 0, \e lon0 = 0, \e h0 =
     * 0.  The optional \e earth argument (default Geocentric::WGS84) specifies
     * the Geocentric object to use for the transformation.
     **********************************************************************/
    LocalCartesian(const Geocentric& earth = Geocentric::WGS84) throw()
      : _earth(earth)
    { Reset(0.0, 0.0, 0.0); }

    /**
     * Change the origin to latitude = \e lat0, longitude = \e lon0 (degrees),
     * height = \e h0 (meters).
     **********************************************************************/
    void Reset(double lat0, double lon0, double h0 = 0) throw();

    /**
     * Convert from geodetic coordinates \e lat, \e lon (degrees), \e h
     * (meters) to local cartesian coordinates \e x, \e y, \e z (meters).  \e
     * lat should be in the range [-90, 90]; \e lon and \e lon0 should be in
     * the range [-180, 360].
     **********************************************************************/
    void Forward(double lat, double lon, double h,
		 double& x, double& y, double& z) const throw();

    /**
     * Convert from local cartesian \e x, \e y, \e z (meters) to geodetic
     * coordinates \e lat, \e lon (degrees), \e h (meters).  The value of \e
     * lon returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(double x, double y, double z,
		 double& lat, double& lon, double& h) const throw();

    /**
     * Return the latitude of the origin (degrees).
     **********************************************************************/
    double LatitudeOrigin() const throw() { return _lat0; }

    /**
     * Return the longitude of the origin (degrees).
     **********************************************************************/
    double LongitudeOrigin() const throw() { return _lon0; }

    /**
     * Return the height of the origin (meters).
     **********************************************************************/
    double HeightOrigin() const throw() { return _h0; }
  };

} // namespace GeographicLib

#endif

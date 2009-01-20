/**
 * \file LocalCartesian.hpp
 * \brief Header for GeographicLib::LocalCartesian class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(LOCALCARTESIAN_HPP)
#define LOCALCARTESIAN_HPP "$Id$"

namespace GeographicLib {

  /**
   * \brief Local Cartesian coordinates
   *
   * Convert between geodetic coordinates latitude = \e lat, longitude = \e
   * lon, height = \e h (measured vertically from the surface of the ellipsoid)
   * to local cartesian coordinates (\e x, \e y, \e z).  The origin of local
   * cartesian coordinate system is at \e lat = \e lat0, \e lon = \e lon0, \e h
   * = 0. The \e z axis is normal to the ellipsoid; the \e y axis points due
   * north.
   *
   * The conversions all take place via ECEF coordinates using
   * GeographicLib::ECEF::WGS84.  (As presently written, there's no provision
   * for changing the ellipsoid.)
   **********************************************************************/

  class LocalCartesian {
  private:
    double _lat0, _lon0;
    double _x0, _y0, _z0,
      _rxx, _rxy, _rxz,
      _ryx, _ryy, _ryz,
      _rzx, _rzy, _rzz;
  public:
    /**
     * Constructor setting the origin to latitude = \e lat0, longitude = \e
     * lon0 (degrees).
     **********************************************************************/
    LocalCartesian(double lat0, double lon0) {
      Reset(lat0, lon0);
    }
    /**
     * Default constructor sets the origin to \e lat0 = 0, \e lon0 = 0.
     **********************************************************************/
    LocalCartesian() { Reset(0.0, 0.0); }
    /**
     * Change the origin to latitude = \e lat0, longitude = \e lon0 (degrees).
     **********************************************************************/
    void Reset(double lat0, double lon0);
    /**
     * Convert from geodetic coordinates \e lat, \e lon (degrees), \e h
     * (meters) to local cartesian \e x, \e y, \e z (meters).
     **********************************************************************/
    void Forward(double lat, double lon, double h,
		 double& x, double& y, double& z) const;
    /**
     * Convert from local cartesian \e x, \e y, \e z (meters) to geodetic
     * coordinates \e lat, \e lon (degrees), \e h (meters).
     **********************************************************************/
    void Reverse(double x, double y, double z,
		 double& lat, double& lon, double& h) const;
    /**
     * Return the latitude of the origin (degrees).
     **********************************************************************/
    double LatitudeOrigin() const { return _lat0; }
    /**
     * Return the longitude of the origin (degrees).
     **********************************************************************/
    double LongitudeOrigin() const { return _lon0; }
  };

} // namespace GeographicLib

#endif

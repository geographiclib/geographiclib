/**
 * \file LocalCartesian.hpp
 * \brief Header for GeographicLib::LocalCartesian class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_LOCALCARTESIAN_HPP)
#define GEOGRAPHICLIB_LOCALCARTESIAN_HPP "$Id: LocalCartesian.hpp 6867 2010-09-11 13:04:26Z karney $"

#include "GeographicLib/Geocentric.hpp"
#include "GeographicLib/Constants.hpp"

namespace GeographicLib {

  /**
   * \brief Local Cartesian coordinates
   *
   * Convert between geodetic coordinates latitude = \e lat, longitude = \e
   * lon, height = \e h (measured vertically from the surface of the ellipsoid)
   * to local cartesian coordinates (\e x, \e y, \e z).  The origin of local
   * cartesian coordinate system is at \e lat = \e lat0, \e lon = \e lon0, \e h
   * = \e h0. The \e z axis is normal to the ellipsoid; the \e y axis points
   * due north.  The plane \e z = - \e h0 is tangent to the ellipsoid.
   *
   * The conversions all take place via geocentric coordinates using a
   * Geocentric object (by default Geocentric::WGS84).
   **********************************************************************/

  class LocalCartesian {
  private:
    typedef Math::real real;
    const Geocentric _earth;
    real _lat0, _lon0, _h0;
    real _x0, _y0, _z0,
      _rxx, _rxy, _rxz,
      _ryx, _ryy, _ryz,
      _rzx, _rzy, _rzz;
  public:

    /**
     * Constructor setting the origin.
     *
     * @param[in] lat0 latitude at origin (degrees).
     * @param[in] lon0 longitude at origin (degrees).
     * @param[in] h0 height above ellipsoid at origin (meters); default 0.
     * @param[in] earth Geocentric object for the transformation; default
     *   Geocentric::WGS84.
     **********************************************************************/
    LocalCartesian(real lat0, real lon0, real h0 = 0,
                   const Geocentric& earth = Geocentric::WGS84) throw()
      : _earth(earth)
    { Reset(lat0, lon0, h0); }

    /**
     * Default constructor.
     *
     * @param[in] earth Geocentric object for the transformation; default
     *   Geocentric::WGS84.
     *
     * Sets \e lat0 = 0, \e lon0 = 0, \e h0 = 0.
     **********************************************************************/
    explicit LocalCartesian(const Geocentric& earth = Geocentric::WGS84)
      throw()
      : _earth(earth)
    { Reset(real(0), real(0), real(0)); }

    /**
     * Reset the origin.
     *
     * @param[in] lat0 latitude at origin (degrees).
     * @param[in] lon0 longitude at origin (degrees).
     * @param[in] h0 height above ellipsoid at origin (meters); default 0.
     **********************************************************************/
    void Reset(real lat0, real lon0, real h0 = 0)
      throw();

    /**
     * Convert from geodetic to local cartesian coordinates.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] h height of point above the ellipsoid (meters).
     * @param[out] x local cartesian coordinate (meters).
     * @param[out] y local cartesian coordinate (meters).
     * @param[out] z local cartesian coordinate (meters).
     *
     * \e lat should be in the range [-90, 90]; \e lon and \e lon0 should be in
     * the range [-180, 360].
     **********************************************************************/
    void Forward(real lat, real lon, real h, real& x, real& y, real& z)
      const throw();

    /**
     * Convert from local cartesian to geodetic to coordinates.
     *
     * @param[in] x local cartesian coordinate (meters).
     * @param[in] y local cartesian coordinate (meters).
     * @param[in] z local cartesian coordinate (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] h height of point above the ellipsoid (meters).
     *
     * The value of \e lon returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(real x, real y, real z, real& lat, real& lon, real& h)
      const throw();

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return latitude of the origin (degrees).
     **********************************************************************/
    Math::real LatitudeOrigin() const throw() { return _lat0; }

    /**
     * @return longitude of the origin (degrees).
     **********************************************************************/
    Math::real LongitudeOrigin() const throw() { return _lon0; }

    /**
     * @return height of the origin (meters).
     **********************************************************************/
    Math::real HeightOrigin() const throw() { return _h0; }

    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value of \e a inherited from the Geocentric object used in the
     *   constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _earth.MajorRadius(); }

    /**
     * @return \e r the inverse flattening of the ellipsoid.  This is the
     *   value of \e r inherited from the Geocentric object used in the
     *   constructor.  A value of 0 is returned for a sphere (infinite inverse
     *   flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw()
    { return _earth.InverseFlattening(); }
    ///@}
  };

} // namespace GeographicLib

#endif

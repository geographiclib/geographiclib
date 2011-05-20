/**
 * \file AzimuthalEquidistant.hpp
 * \brief Header for GeographicLib::AzimuthalEquidistant class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_AZIMUTHALEQUIDISTANT_HPP)
#define GEOGRAPHICLIB_AZIMUTHALEQUIDISTANT_HPP "$Id: AzimuthalEquidistant.hpp 6867 2010-09-11 13:04:26Z karney $"

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
   * The conversions all take place using a Geodesic object (by default
   * Geodesic::WGS84).  For more information on geodesics see \ref geodesic.
   **********************************************************************/

  class AzimuthalEquidistant {
  private:
    typedef Math::real real;
    const Geodesic _earth;
    static const real eps;
  public:

    /**
     * Constructor for AzimuthalEquidistant.
     *
     * @param[in] earth the Geodesic object to use for geodesic calculations.
     *   By default this uses the WGS84 ellipsoid.
     **********************************************************************/
    explicit AzimuthalEquidistant(const Geodesic& earth = Geodesic::WGS84)
      throw() : _earth(earth) {}

    /**
     * Forward projection, from geographic to azimuthal equidistant.
     *
     * @param[in] lat0 latitude of center point of projection (degrees).
     * @param[in] lon0 longitude of center point of projection (degrees).
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] azi azimuth of geodesic at point (degrees).
     * @param[out] rk reciprocal of azimuthal scale at point.
     *
     * \e lat0 and \e lat should be in the range [-90, 90] and \e lon0 and \e
     * lon should be in the range [-180, 360].  The scale of the projection is
     * 1 in the "radial" direction, \e azi clockwise from true north, and is
     * 1/\e rk in the direction perpendicular to this.  A call to Forward
     * followed by a call to Reverse will return the original (\e lat, \e lon)
     * (to within roundoff).
     **********************************************************************/
    void Forward(real lat0, real lon0, real lat, real lon,
                 real& x, real& y, real& azi, real& rk) const throw();

    /**
     * Reverse projection, from azimuthal equidistant to geographic.
     *
     * @param[in] lat0 latitude of center point of projection (degrees).
     * @param[in] lon0 longitude of center point of projection (degrees).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] azi azimuth of geodesic at point (degrees).
     * @param[out] rk reciprocal of azimuthal scale at point.
     *
     * \e lat0 should be in the range [-90, 90] and \e lon0 should be in the
     * range [-180, 360].  \e lat will be in the range [-90, 90] and \e lon
     * will be in the range [-180, 180).  The scale of the projection is 1 in
     * the "radial" direction, \e azi clockwise from true north, and is 1/\e rk
     * in the direction perpendicular to this.  A call to Reverse followed by a
     * call to Forward will return the original (\e x, \e y) (to roundoff) only
     * if the geodesic to (\e x, \e y) is a shortest path.
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

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value inherited from the Geodesic object used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _earth.MajorRadius(); }

    /**
     * @return \e r the inverse flattening of the ellipsoid.  This is the
     *   value inherited from the Geodesic object used in the constructor.  A
     *   value of 0 is returned for a sphere (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw()
    { return _earth.InverseFlattening(); }
    ///@}
  };

} // namespace GeographicLib

#endif

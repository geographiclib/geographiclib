/**
 * \file Geocentric.hpp
 * \brief Header for GeographicLib::Geocentric class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEOCENTRIC_HPP)
#define GEOGRAPHICLIB_GEOCENTRIC_HPP "$Id: Geocentric.hpp 6867 2010-09-11 13:04:26Z karney $"

#include "GeographicLib/Constants.hpp"

namespace GeographicLib {

  /**
   * \brief %Geocentric coordinates
   *
   * Convert between geodetic coordinates latitude = \e lat, longitude = \e
   * lon, height = \e h (measured vertically from the surface of the ellipsoid)
   * to geocentric coordinates (\e x, \e y, \e z).  The origin of geocentric
   * coordinates is at the center of the earth.  The \e z axis goes thru the
   * north pole, \e lat = 90<sup>o</sup>.  The \e x axis goes thru \e lat = 0,
   * \e lon = 0.  Geocentric coordinates are also known as earth centered,
   * earth fixed (ECEF) coordinates.
   *
   * The conversion from geographic to geocentric coordinates is
   * straightforward.  For the reverse transformation we use
   * - H. Vermeille,
   *   <a href="http://dx.doi.org/10.1007/s00190-002-0273-6"> Direct
   *   transformation from geocentric coordinates to geodetic coordinates</a>,
   *   J. Geodesy 76, 451&ndash;454 (2002).
   * .
   * Several changes have been made to ensure that the method returns accurate
   * results for all finite inputs (even if \e h is infinite).  See
   * \ref geocentric for details.
   *
   * The errors in these routines are close to round-off.  Specifically, for
   * points within 5000 km of the surface of the ellipsoid (either inside or
   * outside the ellipsoid), the error is bounded by 7 nm for the WGS84
   * ellipsoid.  See \ref geocentric for further information on the errors.
   **********************************************************************/

  class Geocentric {
  private:
    typedef Math::real real;
    const real _a, _r, _f, _e2, _e2m, _e2a, _e4a, _maxrad;
    static inline real sq(real x) throw() { return x * x; }
  public:

    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters)
     * @param[in] r reciprocal flattening.  Setting \e r = 0 implies \e r = inf
     *   or flattening = 0 (i.e., a sphere).  Negative \e r indicates a prolate
     *   ellipsoid.
     *
     * An exception is thrown if either of the axes of the ellipsoid is
     * non-positive.
     **********************************************************************/
    Geocentric(real a, real r);

    /**
     * Convert from geodetic to geocentric coordinates.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] h height of point above the ellipsoid (meters).
     * @param[out] x geocentric coordinate (meters).
     * @param[out] y geocentric coordinate (meters).
     * @param[out] z geocentric coordinate (meters).
     *
     * \e lat should be in the range [-90, 90]; \e lon and \e lon0 should be in
     * the range [-180, 360].
     **********************************************************************/
    void Forward(real lat, real lon, real h, real& x, real& y, real& z)
      const throw();

    /**
     * Convert from geocentric to geodetic to coordinates.
     *
     * @param[in] x geocentric coordinate (meters).
     * @param[in] y geocentric coordinate (meters).
     * @param[in] z geocentric coordinate (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] h height of point above the ellipsoid (meters).
     *
     * In general there are multiple solutions and the result which maximizes
     * \e h is returned.  If there are still multiple solutions with different
     * latitutes (applies only if \e z = 0), then the solution with \e lat > 0
     * is returned.  If there are still multiple solutions with different
     * longitudes (applies only if \e x = \e y = 0) then \e lon = 0 is
     * returned.  The value of \e h returned satisfies \e h >= - \e a (1 - \e
     * e<sup>2</sup>) / sqrt(1 - \e e<sup>2</sup> sin<sup>2</sup>\e lat).  The
     * value of \e lon returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(real x, real y, real z, real& lat, real& lon, real& h)
      const throw();

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _a; }

    /**
     * @return \e r the inverse flattening of the ellipsoid.  This is the
     *   value used in the constructor.  A value of 0 is returned for a sphere
     *   (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return _r; }
    ///@}

    /**
     * A global instantiation of Geocentric with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static Geocentric WGS84;
  };

} // namespace GeographicLib
#endif

/**
 * \file ECEF.hpp
 * \brief Header for GeographicLib::ECEF class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(ECEF_HPP)
#define ECEF_HPP "$Id$"

namespace GeographicLib {

  /**
   * \brief Earth centered, earth fixed coordinates
   *
   * Convert between geodetic coordinates latitude = \e lat, longitude = \e
   * lon, height = \e h (measured vertically from the surface of the ellipsoid)
   * to earth centered, earth fixed (ECEF) coordinates (\e x, \e y, \e z).  The
   * origin of ECEF coordinates is at the center of the earth.  The \e z axis
   * goes thru the north pole, \e lat = 90<sup>o</sup>.  The \e x axis goes
   * thru \e lat = 0, \e lon = 0.

   * The conversion from geographic to ECEF coordinates is straightforward.
   * For the reverse transformation we use
   * H. Vermeille, <a href="http://dx.doi.org/10.1007/s00190-002-0273-6">
   * Direct transformation from geocentric coordinates to geodetic
   * coordinates</a>, J. Geodesy 76, 451&ndash;454 (2002).  Several changes
   * have been made to ensure that the method returns accurate results for all
   * inputs (provided that \e h is finite).
   * 
   * There's a cube-root singularity in the reverse transformation at sqrt(\e
   * x<sup>2</sup> + \e y<sup>2</sup>) = \e a\e e<sup>2</sup>, \e z = 0.  This
   * point maps to the equator.  However (for the WGS84 ellipsoid), changing \e
   * z = 1 nm, changes the latitude to 7.5" or a distance of 229 m from the
   * equator.  Because of this, when measuring the error in the reverse
   * trasnformation, we distinguish the error in \e h (which is well behaved)
   * from the error in the latitude and longitude.  For the latter, we further
   * distinguish points outside the ellipsoid, in which case we convert the
   * error in the latitude and longitude into a distance on the surface of the
   * ellipsoid, from points inside the ellipsoid, where we instead apply the
   * forward transformation on the result and measure the discrepancy in ECEF
   * coordinates.  The error in the height is bounded by 8 nm max(1, \e h/\e
   * a).  The error in the footprint position (for points outside the
   * ellipsoid) is bounded by 4 nm and the reverse-forward discrepancy (for
   * points inside the ellipsoid) is bounded by 6 nm.
   **********************************************************************/

  class ECEF {
  private:
    const double _a, _f, _e2, _e12, _b, _tol, _maxrad;
    const int _numit;
    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    static inline double cbrt(double x) {
      double y = std::pow(std::abs(x), 1/3.0);
      return x < 0 ? -y : y;
    }
#endif
  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters) and flattening \e f.
     **********************************************************************/
    ECEF(double a, double f);
    /**
     * Convert from geodetic coordinates \e lat, \e lon (degrees), \e h
     * (meters) to ECEF \e x, \e y, \e z (meters).
     **********************************************************************/
    void Forward(double lat, double lon, double h,
		 double& x, double& y, double& z) const;
    /**
     * Convert from ECEF coordinates \e x, \e y, \e z (meters) to geodetic \e
     * lat, \e lon (degrees), \e h (meters).  In general there are multiple
     * solutions and the result which minimizes the absolute value of \e h is
     * returned.  If there are still multiple solutions with different
     * latitutes (applies only if \e z = 0), then the solution with \e lat > 0
     * is returned.  If there are still multiple solutions with different
     * longitudes (applies only if \e x = \e y = 0) then \e lon = 0 is
     * returned.  The value of \e h returned satisfies \e h >= - \e a (1 - \e
     * e<sup>2</sup>) / sqrt(1 - \e e<sup>2</sup> sin<sup>2</sup>\e lat).
     **********************************************************************/
    void Reverse(double x, double y, double z,
		 double& lat, double& lon, double& h) const;

    /**
     * A global instantiation of ECEF with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static ECEF WGS84;
  };

} //namespace GeographicLib
#endif

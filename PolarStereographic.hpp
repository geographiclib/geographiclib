/**
 * \file PolarStereographic.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(POLARSTEREOGRAPHIC_HPP)
#define POLARSTEREOGRAPHIC_HPP "$Id$"

namespace GeographicLib {

  /**
   * \brief Polar Stereographic Projection
   *
   * Implementation taken from the report, J. P.Snyder, <a
   * href="http://pubs.er.usgs.gov/usgspubs/pp/pp1395"> Map Projections: A
   * Working Manual</a>, USGS Professional Paper 1395 (1987), pp. 160&ndash;163.
   **********************************************************************/
  class PolarStereographic {
  private:
    const double _a, _f, _k0, _e, _e2m, _c, _tol;
    const int _numit;
    static inline double sq(double x) { return x * x; }
  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters), flattening \e f,
     * and central scale factor \e k0.
     **********************************************************************/
    PolarStereographic(double a, double f, double k0);
    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * polar stereographic easting \e x (meters) and northing \e y (meters).
     * The projection is about the pole given by \e northp (false means south,
     * true means north).  Also return the meridian convergence \e gamma
     * (degrees) and the scale \e k.  No false easting or northing is added.
     **********************************************************************/
    void Forward(bool northp, double lat, double lon,
		 double& x, double& y, double& gamma, double& k) const;
    /**
     * Convert from polar stereogrphic easting \e x (meters) and northing \e y
     * (meters) to latitude \e lat (degrees) and longitude \e lon (degrees) .
     * The hemisphere is given by \e northp (false means south, true means
     * north).  Also return the meridian convergence \e gamma (degrees) and the
     * scale \e k.  No false easting or northing is added.
     **********************************************************************/
    void Reverse(bool northp, double x, double y,
		 double& lat, double& lon, double& gamma, double& k) const;
    /**
     * A global instantiation of PolarStereographic with thw WGS84 ellipsoid
     * and the UPS scale factor.
     **********************************************************************/
    const static PolarStereographic UPS;
  };

} // namespace GeographicLib

#endif

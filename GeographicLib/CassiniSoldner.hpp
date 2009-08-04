/**
 * \file CassiniSoldner.hpp
 * \brief Header for GeographicLib::CassiniSoldner class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(CASSINISOLDNER_HPP)
#define CASSINISOLDNER_HPP "$Id$"

#include "GeographicLib/Geodesic.hpp"

namespace GeographicLib {

  /**
   * \brief Cassini-Solder Projection.
   *
   * Cassini-Solder projection centered at an arbitrary position on the
   * ellipsoid.
   *
   * The conversions all take place using a GeographicLib::Geodesic object (by
   * default GeographicLib::Geodesic::WGS84).  For more information on
   * geodesics see \ref geodesic.
   **********************************************************************/

  class CassiniSoldner {
  private:
    const Geodesic& _earth;
    GeodesicLine _meridian;
    static inline double sq(double x) throw() { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) throw()
    { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) throw()
    { return ::hypot(x, y); }
#endif
  public:

    /**
     * Constructor for CassiniSoldner setting the Geodesic object to use
     * for geodesic calculations.  By default this uses the WGS84 ellipsoid.
     **********************************************************************/
    CassiniSoldner(const Geodesic& earth = Geodesic::WGS84) throw()
      : _earth(earth) {}

    CassiniSoldner(double lat0, double lon0,
		   const Geodesic& earth = Geodesic::WGS84) throw()
      : _earth(earth) {
      Reset(lat0, lon0);
    }

    void Reset(double lat0, double lon0) throw();

    void Forward(double lat, double lon,
		 double& x, double& y, double& azi, double& m) const throw();

    void Reverse(double x, double y,
		 double& lat, double& lon,
		 double& azi, double& m) const throw();

    /**
     * Has this object been initialized with an origin?
     **********************************************************************/
    bool Init() const throw() { return _meridian.Init(); }

    /**
     * Return the latitude of the origin (degrees).
     **********************************************************************/
    double LatitudeOrigin() const throw()
    { return _meridian.Latitude(); }

    /**
     * Return the longitude of the origin (degrees).
     **********************************************************************/
    double LongitudeOrigin() const throw()
    { return _meridian.Longitude(); }

  };

} // namespace GeographicLib

#endif

/**
 * \file Geodesic.hpp
 * \brief Header for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(GEODESIC_HPP)
#define GEODESIC_HPP "$Id$"

#include <cmath>

namespace GeographicLib {

  /**
   * \brief Geodesic calculations
   *
   **********************************************************************/

  class Geodesic {
  private:
    const double _a, _f, _e2, _ep2, _b, _eps2, _tol;
    const int _numit;
    static const int maxpow = 8;
    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) { return ::hypot(x, y); }
#endif
  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters) and inverse flattening
     * \e invf.
     **********************************************************************/
    Geodesic(double a, double invf);
    /**
     * Perform the direct geodesic calculation.
     **********************************************************************/
    void Direct(double lat1, double lon1, double alpha1, double s,
		 double& lat2, double& lon2, double& alpha2) const;
    /**
     * Perform the inverse geodesic calculation.
     **********************************************************************/
    void Inverse(double lat1, double lon1, double lat2, double lon2,
		 double& s, double& alpha1, double& alpha2) const;

    /**
     * A global instantiation of Geodesic with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static Geodesic WGS84;
  };

  class GeodesicLine {
  private:
    static const int maxpow = 8;
    double _lat0, _lon0, _bearing0;
    double _b, _f, _u2, _sig;
    double _d[maxpow];
  protected:
    GeodesicLine(const Geodesic& g,
		 double lat, double lon, double bearing);
  public:
    // A default public constuctor
    GeodesicLine() : _b(0) {};
    void Position(double s, double& lat, double& lon, double& bearing)
      const throw();
    double Longitude() const throw() { return _lon0; }
    double Latitude() const throw() { return _lat0; }
    double Bearing() const throw() { return _bearing0; }
    
} //namespace GeographicLib
#endif

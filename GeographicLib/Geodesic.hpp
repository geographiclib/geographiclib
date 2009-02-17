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

  class GeodesicLine;

  class Geodesic {
  protected:
    friend class GeodesicLine;
    static const double eps2, tol;
    const double _a, _f, _f1, _e2, _ep2, _b;
    static const int maxpow = 8;
    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) { return ::hypot(x, y); }
#endif
    static double SinSeries(double x, const double c[], int n) throw();
    static inline double AngNormalize(double x) throw() {
      // Place angle in [-180, 180).  Assumes x is in [-540, 540).
      return x >= 180 ? x - 360 : x < -180 ? x + 360 : x;
    }
    static inline double AngRound(double x) throw() {
      // The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
      // for doubles = 0.7 pm on the earth if x is an angle in degrees.  (This
      // is about 1000 times more resolution than we get with angles around 90
      // degrees.)  We use this to avoid having to deal with near singular
      // cases when x is non-zero but tiny (e.g., 1.0e-200).
      const double z = 0.0625;	// 1/16
      double y = std::abs(x);
      // The compiler mustn't "simplify" z - (z - y) to y
      y = y < z ? z - (z - y) : y;
      return x < 0 ? -y : y;
    }

    static double sigmaScale(double u2) throw();
    static void sigmaCoeffSet(double u2, double c[]) throw();
    static void sCoeffSet(double u2, double c[]) throw();
    static double dlambdaScale(double f, double mu) throw();
    static void dlambdaCoeffSet(double f, double mu, double e[]) throw();

  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters) and inverse flattening
     * \e invf.
     **********************************************************************/
    Geodesic(double a, double invf);
    /**
     * Set up to do a series of ranges
     **********************************************************************/
    GeodesicLine Line(double lat1, double lon1, double bearing1) const throw();
    /**
     * Perform the direct geodesic calculation.
     **********************************************************************/
    void Direct(double lat1, double lon1, double bearing1, double s12,
		double& lat2, double& lon2, double& bearing2) const throw();
    /**
     * Perform the inverse geodesic calculation.
     **********************************************************************/
    void Inverse(double lat1, double lon1, double lat2, double lon2,
		 double& s12, double& bearing1, double& bearing2) const throw();

    /**
     * A global instantiation of Geodesic with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static Geodesic WGS84;
  };


  class GeodesicLine {
  private:
    static const int maxpow = 8;

    int _bsign;
    double _lat1, _lon1, _bearing1;
    double _sScale, _S1, _cosalpha0, _sinalpha0, _f1, _dlambdaScale, _chi1;
    double _sigmaCoeff[maxpow], _dlambdaCoeff[maxpow];

  protected:
    friend class Geodesic;
    GeodesicLine(const Geodesic& g,
		 double lat, double lon, double bearing);
  public:
    // A default public constuctor
    GeodesicLine() : _sScale(0) {};
    void Position(double s, double& lat, double& lon, double& bearing)
      const throw();
    double Longitude() const throw() { return _lon1; }
    double Latitude() const throw() { return _lat1; }
    double Bearing() const throw() { return _bsign * _bearing1; }
  };

} //namespace GeographicLib
#endif

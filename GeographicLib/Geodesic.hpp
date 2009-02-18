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
    static double SinSeries(double x, const double c[], int n) throw();
    static double SinSeries(double sinx, double cosx, const double c[], int n)
      throw();
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

    static double sigScale(double u2) throw();
    static void sigCoeffSet(double u2, double c[]) throw();
    static void sCoeffSet(double u2, double c[]) throw();
    static double dlamScale(double f, double mu) throw();
    static void dlamCoeffSet(double f, double mu, double e[]) throw();
  private:
    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) { return ::hypot(x, y); }
#endif
    double ChiDiff(double sbet1, double cbet1,
		   double sbet2, double cbet2,
		   double salp1, double calp1,
		   double& salp2, double& calp2,
		   double& sig12,
		   double& ssig1, double& csig1,
		   double& ssig2, double& csig2,
		   double& u2,
		   double c[]) const throw();

  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters) and inverse flattening
     * \e invf.
     **********************************************************************/
    Geodesic(double a, double invf);
    /**
     * Set up to do a series of ranges
     **********************************************************************/
    GeodesicLine Line(double lat1, double lon1, double head1) const throw();
    /**
     * Perform the direct geodesic calculation.
     **********************************************************************/
    void Direct(double lat1, double lon1, double head1, double s12,
		double& lat2, double& lon2, double& head2) const throw();
    /**
     * Perform the inverse geodesic calculation.
     **********************************************************************/
    void Inverse(double lat1, double lon1, double lat2, double lon2,
		 double& s12, double& head1, double& head2) const throw();

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
    double _lat1, _lon1, _head1;
    double _sScale, _S1, _calp0, _salp0, _f1, _dlamScale, _chi1;
    double _sigCoeff[maxpow], _dlamCoeff[maxpow];
    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) { return ::hypot(x, y); }
#endif

  protected:
    friend class Geodesic;
    GeodesicLine(const Geodesic& g,
		 double lat1, double lon1, double head1);
  public:
    // A default public constuctor
    GeodesicLine() : _sScale(0) {};
    void Position(double s12, double& lat2, double& lon2, double& head2)
      const throw();
    double Longitude() const throw() { return _lon1; }
    double Latitude() const throw() { return _lat1; }
    double Heading() const throw() { return _bsign * _head1; }
  };

} //namespace GeographicLib
#endif

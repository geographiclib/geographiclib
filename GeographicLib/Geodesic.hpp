/**
 * \file Geodesic.hpp
 * \brief Header for GeographicLib::Geodesic and GeographicLib::GeodesicLine classes
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(GEODESIC_HPP)
#define GEODESIC_HPP "$Id$"

#include <cmath>

namespace GeographicLib {

  class GeodesicLine;

  /**
   * \brief Geodesic calculations
   *
   * Direct and inverse geodesic calculations.
   **********************************************************************/

  class Geodesic {
  private:
    friend class GeodesicLine;
    static const int maxpow = 8, head2sense = 1;

    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) { return _hypot(x, y); }
#else
    static inline double hypot(double x, double y) { return ::hypot(x, y); }
#endif
    double Chi12(double sbet1, double cbet1, double sbet2, double cbet2,
		 double salp1, double calp1, double& salp2, double& calp2,
		 double& sig12,
		 double& ssig1, double& csig1, double& ssig2, double& csig2,
		 double& u2, bool diffp, double& dchi12, double c[]) const throw();

    static const double eps2, tol;
    const double _a, _f, _f1, _e2, _ep2, _b;
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
    static inline void SinCosNorm(double& sinx, double& cosx) throw() {
      double r = hypot(sinx, cosx);
      sinx /= r;
      cosx /= r;
    }

    static double tauScale(double u2) throw();
    static void tauCoeff(double u2, double c[]) throw();
    static void sigCoeff(double u2, double c[]) throw();
    static double dlamScale(double f, double mu) throw();
    static void dlamCoeff(double f, double mu, double e[]) throw();
    static double dlamScalemu(double f, double mu) throw();
    static void dlamCoeffmu(double f, double mu, double e[]) throw();

  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters) and inverse flattening
     * \e invf.  Setting \e invf <= 0 implies \e invf = inf or flattening = 0
     * (i.e., a sphere).
     **********************************************************************/
    Geodesic(double a, double invf);
    /**
     * Perform the direct geodesic calculation.  Given a latitude, \e lat1,
     * longitude, \e lon1, and heading \e head1 (in degrees) for point 1 and a
     * range, \e s12 (in meters) from point 1 to point 2, return the latitude,
     * \e lat2, longitude, \e lon2, and forward heading, \e head2 (in degees)
     * for point 2.
     **********************************************************************/
    void Direct(double lat1, double lon1, double head1, double s12,
		double& lat2, double& lon2, double& head2) const throw();
    /**
     * Set up to do a series of ranges.  This returns a GeodesicLine object
     * with point 1 given by latitude, \e lat1, longitude, \e lon1, and heading
     * \e head1 (in degrees).  Calls to GeodesicLine::Position return the
     * position and heading for point 2 a specified distance away.
     **********************************************************************/
    GeodesicLine Line(double lat1, double lon1, double head1) const throw();
    /**
     * Perform the inverse geodesic calculation.  Given a latitude, \e lat1,
     * longitude, \e lon1, for point 1 and a latitude, \e lat2, longitude, \e
     * lon2, for point 2 (all in degrees), return the geodesic distance, \e s12
     * (in meters), and the forward headings, \e head1 and \e head2 (in
     * degrees), at points 1 and 2.
     **********************************************************************/
    void Inverse(double lat1, double lon1, double lat2, double lon2,
		 double& s12, double& head1, double& head2) const throw();

    /**
     * A global instantiation of Geodesic with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static Geodesic WGS84;
  };

  /**
   * \brief A geodesic line.
   *
   * Calculate multiple points on a single geodesic line specified by a point 1
   * and a heading at that point.  Geodesic.Line allows point 1 and the heading
   * to be specified and returns a GeodesicLine object.  GeodesicLine.Position
   * returns the position and head at point 2 a give distance away.  An example
   * of use of this class is:
   \verbatim
   // Print positions on a geodesic going through latitude 30,
   // longitude 10 at heading 80.  Points at intervals of 10km
   // in the range [-1000km, 1000km] are given.
   GeodesicLine line(Geodesic::WGS84.Line(30.0, 10.0, 80.0));
   double step = 10e3;
   for (int s = -100; s <= 100; ++s) {
     double lat2, lon2, head2;
     double s12 = s * step;
     line.Position(s12, lat2, lon2, head2);
     cout << s12 << " " << lat2 << " " << lon2 << " " << head2 << "\n";
   }
   \endverbatim
   **********************************************************************/

  class GeodesicLine {
  private:
    friend class Geodesic;
    static const int maxpow = 8;

    int _bsign;
    double _lat1, _lon1, _head1;
    double  _f1, _salp0, _calp0,
      _ssig1, _csig1, _stau1, _ctau1, _slam1, _clam1,
      _sScale, _dlamScale, _dtau1, _dchi1;
    double _sigCoeff[maxpow], _dlamCoeff[maxpow];

    GeodesicLine(const Geodesic& g,
		 double lat1, double lon1, double head1);
  public:
    /**
     * A default constructor.  If Position is called on the resulting object,
     * it returns immediately (without doing any calculations).  The object
     * should be set with a call to Geodesic::Line.  Use Init() to test whether
     * object is still in this uninitialized state.
     **********************************************************************/
    GeodesicLine() : _sScale(0) {};

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward heading,
     * \e head2 (in degrees) of the point 2 which is a distance, \e s12
     * (meters), from point 1.  \e s12 can be signed.
     **********************************************************************/
    void Position(double s12, double& lat2, double& lon2, double& head2)
      const throw();

    /**
     * Has this object been initialize so that Position can be called?
     **********************************************************************/
    bool Init() const throw() { return _sScale > 0; }
    /**
     * Return the latitude of point 1 (in degrees).
     **********************************************************************/
    double Latitude() const throw() { return _lat1; }
    /**
     * Return the longitude of point 1 (in degrees).
     **********************************************************************/
    double Longitude() const throw() { return _lon1; }
    /**
     * Return the heading of the geodesic line as it passes through point 1.
     **********************************************************************/
    double Heading() const throw() { return _bsign * _head1; }
  };

} //namespace GeographicLib
#endif

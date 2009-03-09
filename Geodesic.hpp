/**
 * \file Geodesic.hpp
 * \brief Header for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEODESIC_HPP)
#define GEODESIC_HPP "$Id$"

#define ITER 1

// These correspond to the max and min orders in Geodesic.cpp
#define MINPOW 1
#define MAXPOW 8
#define MAXPOW_CLAMP(x) ((x) > MAXPOW ? MAXPOW : ((x) < MINPOW ? MINPOW : (x)))

#if !defined(MAXPOW_DEFAULT)
#define MAXPOW_DEFAULT 6
#endif

#define MAXPOW_TAUSC MAXPOW_CLAMP(MAXPOW_DEFAULT)
#define MAXPOW_TAUCOEF MAXPOW_CLAMP(MAXPOW_DEFAULT)
#define MAXPOW_SIGCOEF MAXPOW_CLAMP(MAXPOW_DEFAULT)
#define MAXPOW_LAMSC MAXPOW_CLAMP(MAXPOW_DEFAULT-1)
#define MAXPOW_LAMCOEF MAXPOW_CLAMP(MAXPOW_DEFAULT-2)

#define ALTAZI 0

#include <cmath>

namespace GeographicLib {

  class GeodesicLine;

  /**
   * \brief %Geodesic calculations
   *
   * The shortest path between two points on the ellipsoid at (\e lat1, \e
   * lon1) and (\e lat2, \e lon2) is called the geodesic.  Its length is \e s12
   * and the geodesic from point 1 to point 2 has azimuths \e azi1 and \e azi2
   * at the two end points.  (The azimuth is the heading measured clockwise
   * from north.  \e azi2 is the "forward" azimuth, i.e., the heading that
   * takes you beyond point 2 not back to point 1.)
   *
   * Given \e lat1, \e lon1, \e azi1, and \e s12, we can determine \e lat2, \e
   * lon2, \e azi2.  This is the \e direct geodesic problem.  (If \e s12 is
   * sufficiently large that the geodesic wraps more than halfway around the
   * earth, there will be a true geodesic between the points with a smaller \e
   * s12.)
   *
   * Given \e lat1, \e lon1, \e lat2, and \e lon2, we can determine \e azi1, \e
   * azi2, \e s12.  This is the \e inverse geodesic problem.  Usually, the
   * solution to the inverse problem is unique.  In cases where there are
   * muliple solutions (all with the same \e s12, of course), all the solutions
   * can be easily generated once a particular solution is provided.
   *
   * The calculations are accurate to better than 12 nm.  (See \ref geoderrors
   * for details.)
   **********************************************************************/

  class Geodesic {
  private:
    friend class GeodesicLine;
    static const int azi2sense = 1;
    static const int maxpow_taucoef = MAXPOW_TAUCOEF;
    static const int maxpow_lamcoef = MAXPOW_LAMCOEF;

    static inline double sq(double x) throw() { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) throw()
    { return _hypot(x, y); }
    static inline double cbrt(double x) throw() {
      double y = std::pow(std::abs(x), 1/3.0);
      return x < 0 ? -y : y;
    }
#else
    static inline double hypot(double x, double y) throw()
    { return ::hypot(x, y); }
    static inline double cbrt(double x) throw() { return ::cbrt(x); }
#endif
    double Lambda12(double sbet1, double cbet1, double sbet2, double cbet2,
		    double salp1, double calp1,
		    double& salp2, double& calp2,
		    double& sig12,
		    double& ssig1, double& csig1, double& ssig2, double& csig2,
		    double& u2, bool diffp, double& dlam12, double c[])
      const throw();

    static const double eps2, tol, tol1, xthresh;
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
     * Constructor for a ellipsoid radius \e a (meters) and reciprocal flattening
     * \e r.  Setting \e r <= 0 implies \e r = inf or flattening = 0
     * (i.e., a sphere).
     **********************************************************************/
    Geodesic(double a, double r) throw();

    /**
     * Perform the direct geodesic calculation.  Given a latitude, \e lat1,
     * longitude, \e lon1, and azimuth \e azi1 (in degrees) for point 1 and a
     * range, \e s12 (in meters) from point 1 to point 2, return the latitude,
     * \e lat2, longitude, \e lon2, and forward azimuth, \e azi2 (in degees)
     * for point 2.
     **********************************************************************/
    void Direct(double lat1, double lon1, double azi1, double s12,
		double& lat2, double& lon2, double& azi2) const throw();

    /**
     * Set up to do a series of ranges.  This returns a GeodesicLine object
     * with point 1 given by latitude, \e lat1, longitude, \e lon1, and azimuth
     * \e azi1 (in degrees).  Calls to GeodesicLine::Position return the
     * position and azimuth for point 2 a specified distance away.  Using
     * GeodesicLine::Position is approximately 2.5 faster than calling
     * Geodesic::Direct.
     **********************************************************************/
    GeodesicLine Line(double lat1, double lon1, double azi1) const throw();

    /**
     * Perform the inverse geodesic calculation.  Given a latitude, \e lat1,
     * longitude, \e lon1, for point 1 and a latitude, \e lat2, longitude, \e
     * lon2, for point 2 (all in degrees), return the geodesic distance, \e s12
     * (in meters), and the forward azimuths, \e azi1 and \e azi2 (in
     * degrees), at points 1 and 2.
     **********************************************************************/
    void Inverse(double lat1, double lon1, double lat2, double lon2,
		 double& s12, double& azi1, double& azi2) const throw();


    /**
     * A global instantiation of Geodesic with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static Geodesic WGS84;
#if ITER
    mutable int iter, iterx;
#endif
  };

  /**
   * \brief A geodesic line.
   *
   * GeodesicLine facilitates the determination of a series of points on a
   * single geodesic.  Geodesic.Line returns a GeodesicLine object with the
   * geodesic defined by by \e lat1, \e lon1, and \e azi1.
   * GeodesicLine.Position returns the \e lat2, \e lon2, and \e azi2 given \e
   * s12.  An example of use of this class is:
   \verbatim
   // Print positions on a geodesic going through latitude 30,
   // longitude 10 at azimuth 80.  Points at intervals of 10km
   // in the range [-1000km, 1000km] are given.
   GeodesicLine line(Geodesic::WGS84.Line(30.0, 10.0, 80.0));
   double step = 10e3;
   for (int s = -100; s <= 100; ++s) {
     double lat2, lon2, azi2;
     double s12 = s * step;
     line.Position(s12, lat2, lon2, azi2);
     cout << s12 << " " << lat2 << " " << lon2 << " " << azi2 << "\n";
   }
   \endverbatim
   * The default copy constructor and assignment operators work with this
   * class, so that, for example, the previous example could start
   \verbatim
   GeodesicLine line;
   line = Geodesic::WGS84.Line(30.0, 10.0, 80.0);
   ...
   \endverbatim
   * Similarly, a vector can be used to hold GeodesicLine objects.
   *
   * The calculations are accurate to better than 12 nm.  (See \ref geoderrors
   * for details.)
   **********************************************************************/

  class GeodesicLine {
  private:
    friend class Geodesic;
    static const int maxpow_taucoef = MAXPOW_TAUCOEF;
    static const int maxpow_sigcoef = MAXPOW_SIGCOEF;
    static const int maxpow_lamcoef = MAXPOW_LAMCOEF;

    int _bsign;
    double _lat1, _lon1, _azi1;
    double  _f1, _salp0, _calp0,
      _ssig1, _csig1, _stau1, _ctau1, _schi1, _cchi1,
      _sScale, _dlamScale, _dtau1, _dlam1;
    double _sigCoeff[maxpow_taucoef > maxpow_sigcoef ?
		     maxpow_taucoef : maxpow_sigcoef],
      _dlamCoeff[maxpow_lamcoef];

    GeodesicLine(const Geodesic& g, double lat1, double lon1, double azi1)
      throw();
  public:

    /**
     * A default constructor.  If GeodesicLine::Position is called on the
     * resulting object, it returns immediately (without doing any
     * calculations).  The object should be set with a call to Geodesic::Line.
     * Use Init() to test whether object is still in this uninitialized state.
     **********************************************************************/
    GeodesicLine() throw() : _sScale(0) {};

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward azimuth,
     * \e azi2 (in degrees) of the point 2 which is a distance, \e s12
     * (meters), from point 1.  \e s12 can be signed.
     **********************************************************************/
    void Position(double s12, double& lat2, double& lon2, double& azi2)
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
     * Return the azimuth of the geodesic line as it passes through point 1.
     **********************************************************************/
    double Azimuth() const throw() { return _bsign * _azi1; }
  };

} //namespace GeographicLib
#endif

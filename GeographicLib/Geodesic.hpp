/**
 * \file Geodesic.hpp
 * \brief Header for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEODESIC_HPP)
#define GEOGRAPHICLIB_GEODESIC_HPP "$Id$"

#if !defined(GEOD_TAU_ORD)
/**
 * The order of the series relating \e s and \e sigma when expanded in powers
 * of \e k1.
 **********************************************************************/
#define GEOD_TAU_ORD 5
#endif

#if !defined(GEOD_ETA_ORD)
/**
 * The order of the series relating \e lambda and \e eta when expanded in
 * powers of \e fp.
 **********************************************************************/
#define GEOD_ETA_ORD 5
#endif

#include <cmath>

namespace GeographicLib {

  class GeodesicLine;

  /**
   * \brief %Geodesic calculations
   *
   * The shortest path between two points on a ellipsoid at (\e lat1, \e lon1)
   * and (\e lat2, \e lon2) is called the geodesic.  Its length is \e s12 and
   * the geodesic from point 1 to point 2 has azimuths \e azi1 and \e azi2 at
   * the two end points.  (The azimuth is the heading measured clockwise from
   * north.  \e azi2 is the "forward" azimuth, i.e., the heading that takes you
   * beyond point 2 not back to point 1.)
   *
   * If we fix the first point and increase \e s12 by \e ds12, then the second
   * point is displaced \e ds12 in the direction \e azi2.  Similarly we
   * increase \e azi1 by \e dazi1 (radians), the the second point is displaced
   * \e m12 \e dazi1 in the direction \e azi2 + 90<sup>o</sup>.  The quantity
   * \e m12 is called the "reduced length" and is symmetric under interchange
   * of the two points.  On a flat surface, he have \e m12 = \e s12.  The ratio
   * \e s12/\e m12 gives the azimuthal scale for an azimuthal equidistant
   * projection.
   *
   * Given \e lat1, \e lon1, \e azi1, and \e s12, we can determine \e lat2, \e
   * lon2, \e azi2, \e m12.  This is the \e direct geodesic problem.  (If \e
   * s12 is sufficiently large that the geodesic wraps more than halfway around
   * the earth, there will be another geodesic between the points with a
   * smaller \e s12.)
   *
   * Given \e lat1, \e lon1, \e lat2, and \e lon2, we can determine \e azi1, \e
   * azi2, \e s12, \e m12.  This is the \e inverse geodesic problem.  Usually,
   * the solution to the inverse problem is unique.  In cases where there are
   * muliple solutions (all with the same \e s12, of course), all the solutions
   * can be easily generated once a particular solution is provided.
   *
   * As an alternative to using distance to measure \e s12, the class can also
   * use the arc length \e a12 (in degrees) on the auxiliary sphere.  This is a
   * mathematical construct used in solving the geodesic problems.  However, an
   * arc length in excess of 180<sup>o</sup> indicates that the geodesic is not
   * a shortest path.  In addition, the arc length between an equatorial
   * crossing and the next extremum of latitude for a geodesic is
   * 90<sup>o</sup>.
   *
   * The calculations are accurate to better than 15 nm.  (See \ref geoderrors
   * for details.)
   *
   * For more information on geodesics see \ref geodesic.
   **********************************************************************/

  class Geodesic {
  private:
    friend class GeodesicLine;
    friend class CassiniSoldner;
    static const int tauord = GEOD_TAU_ORD;
    static const int ntau = tauord;
    static const int nsig = tauord;
    static const int zetord = GEOD_TAU_ORD;
    static const int nzet = zetord;
    static const int etaord = GEOD_ETA_ORD;
    // etaCoeff is multiplied by etaFactor which is O(f), so we reduce the
    // order to which etaCoeff is computed by 1.
    static const int neta = etaord > 0 ? etaord - 1 : 0;
    static const unsigned maxit = 50;

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
    void Lengths(double k1, double sig12,
		 double ssig1, double csig1, double ssig2, double csig2,
		 double cbet1, double cbet2,
		 double& s12s, double& m12a, double& m0,
		 double tc[], double zc[]) const throw();
    static void Evolute(double R, double z, double& c, double& s) throw();
    void InverseStart(double sbet1, double cbet1, double sbet2, double cbet2,
		      double lam12, double slam12, double clam12,
		      double& salp1, double& calp1,
		      double tc[], double zc[]) const throw();
    double Lambda12(double sbet1, double cbet1, double sbet2, double cbet2,
		    double salp1, double calp1,
		    double& salp2, double& calp2, double& sig12,
		    double& ssig1, double& csig1, double& ssig2, double& csig2,
		    double& k1, bool diffp, double& dlam12,
		    double tc[], double zc[], double ec[])
      const throw();

    static const double eps2, tol0, tol1, tol2, xthresh;
    const double _a, _f, _f1, _e2, _ep2, _n, _b;
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

    // These are Maxima generated functions to provide series approximations to
    // the integrals for the ellipsoidal geodesic.
    static double tauFactorm1(double k1) throw();
    static void tauCoeff(double k1, double t[]) throw();
    static void sigCoeff(double k1, double tp[]) throw();
    static double zetFactorm1(double k1) throw();
    static void zetCoeff(double k1, double t[]) throw();
    static double etaFactor(double f, double k1) throw();
    static void etaCoeff(double f, double k1, double h[]) throw();

  public:

    /**
     * Constructor for a ellipsoid with equatorial radius \e a (meters) and
     * reciprocal flattening \e r.  Setting \e r = 0 implies \e r = inf or
     * flattening = 0 (i.e., a sphere).  Negative \e r indicates a prolate
     * ellipsoid.
     **********************************************************************/
    Geodesic(double a, double r) throw();

    /**
     * Perform the direct geodesic calculation.  Given a latitude, \e lat1,
     * longitude, \e lon1, and azimuth \e azi1 (degrees) for point 1 and a
     * range, \e s12 (meters) from point 1 to point 2, return the latitude, \e
     * lat2, longitude, \e lon2, and forward azimuth, \e azi2 (degrees) for
     * point 2 and the reduced length \e m12 (meters).  If \e arcmode (default
     * false) is set to true, \e s12 is interpreted as the arc length \e a12
     * (degrees) on the auxiliary sphere.  An arc length greater that 180
     * degrees results in a geodesic which is not a shortest path.  For a
     * prolate ellipsoid, an additional condition is necessary for a shortest
     * path: the longitudinal extent must not exceed of 180 degrees.  Returned
     * value is the arc length \e a12 (degrees) if \e arcmode is false,
     * otherwise it is the distance \e s12 (meters)
     **********************************************************************/
    double Direct(double lat1, double lon1, double azi1, double s12,
		  double& lat2, double& lon2, double& azi2, double& m12,
		  bool arcmode = false) const throw();

    /**
     * Perform the inverse geodesic calculation.  Given a latitude, \e lat1,
     * longitude, \e lon1, for point 1 and a latitude, \e lat2, longitude, \e
     * lon2, for point 2 (all in degrees), return the geodesic distance, \e s12
     * (meters), the forward azimuths, \e azi1 and \e azi2 (degrees) at points
     * 1 and 2, and the reduced length \e m12 (meters).  Returned value is the
     * arc length \e a12 (degrees) on the auxiliary sphere.  The routine uses
     * an iterative method.  If the method fails to converge, the negative of
     * the distances (\e s12, \e m12, and \e a12) and reverse of the azimuths
     * are returned.  This is not expected to happen with ellipsoidal models of
     * the earth.  Please report all cases where this occurs.
     **********************************************************************/
    double Inverse(double lat1, double lon1, double lat2, double lon2,
		   double& s12, double& azi1, double& azi2, double& m12)
      const throw();

    /**
     * Set up to do a series of ranges.  This returns a GeodesicLine object
     * with point 1 given by latitude, \e lat1, longitude, \e lon1, and azimuth
     * \e azi1 (degrees).  Calls to GeodesicLine::Position return the
     * position and azimuth for point 2 a specified distance away.  Using
     * GeodesicLine::Position is approximately 2.1 faster than calling
     * Geodesic::Direct.
     **********************************************************************/
    GeodesicLine Line(double lat1, double lon1, double azi1) const throw();

    /**
     * A global instantiation of Geodesic with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static Geodesic WGS84;
  };

  /**
   * \brief A geodesic line.
   *
   * GeodesicLine facilitates the determination of a series of points on a
   * single geodesic.  Geodesic.Line returns a GeodesicLine object with the
   * geodesic defined by by \e lat1, \e lon1, and \e azi1.
   * GeodesicLine.Position returns the \e lat2, \e lon2, \e azi2, and \e m12
   * given \e s12.  An example of use of this class is:
   \verbatim
   // Print positions on a geodesic going through latitude 30,
   // longitude 10 at azimuth 80.  Points at intervals of 10km
   // in the range [-1000km, 1000km] are given.
   GeodesicLine line(Geodesic::WGS84.Line(30.0, 10.0, 80.0));
   double step = 10e3;
   for (int s = -100; s <= 100; ++s) {
     double lat2, lon2, azi2, m12;
     double s12 = s * step;
     line.Position(s12, lat2, lon2, azi2, m12);
     cout << s12 << " " << lat2 << " " << lon2 << " "
          << azi2 << " " << m12 << "\n";
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
    friend class CassiniSoldner;
    static const int ntau = Geodesic::ntau;
    static const int nsig = Geodesic::nsig;
    static const int nzet = Geodesic::nzet;
    static const int neta = Geodesic::neta;

    double _lat1, _lon1, _azi1;
    double  _b, _f1, _salp0, _calp0, _u2,
      _ssig1, _csig1, _stau1, _ctau1, _somg1, _comg1,
      _taufm1, _zetfm1, _etaf, _dtau1, _dzet1, _dlam1;
    double _tauCoeff[ntau ? ntau : 1], _sigCoeff[nsig ? nsig : 1],
      _zetCoeff[nzet ? nzet : 1], _etaCoeff[neta ? neta : 1];

    GeodesicLine(const Geodesic& g, double lat1, double lon1, double azi1)
      throw();
  public:

    /**
     * A default constructor.  If GeodesicLine::Position is called on the
     * resulting object, it returns immediately (without doing any
     * calculations).  The object should be set with a call to Geodesic::Line.
     * Use Init() to test whether object is still in this uninitialized state.
     **********************************************************************/
    GeodesicLine() throw() : _b(0) {};

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward azimuth,
     * \e azi2 (degrees) of the point 2 which is a distance, \e s12 (in
     * meters), from point 1.  Also return the reduced length \e m12 (meters).
     * \e s12 can be signed.  If \e arcmode (default false) is set to true, \e
     * s12 is interpreted as the arc length \e a12 (in degrees) on the
     * auxiliary sphere.  Returned value is the arc length \e a12 (degrees) if
     * \e arcmode is false, otherwise it is the distance \e s12 (meters).
     **********************************************************************/
    double Position(double s12, double& lat2, double& lon2, double& azi2,
		    double &m12, bool arcmode = false)
      const throw();

    /**
     * Has this object been initialized so that Position can be called?
     **********************************************************************/
    bool Init() const throw() { return _b > 0; }

    /**
     * Return the latitude of point 1 (degrees).
     **********************************************************************/
    double Latitude() const throw() { return _lat1; }

    /**
     * Return the longitude of point 1 (degrees).
     **********************************************************************/
    double Longitude() const throw() { return _lon1; }

    /**
     * Return the azimuth (degrees) of the geodesic line as it passes through
     * point 1.
     **********************************************************************/
    double Azimuth() const throw() { return _azi1; }
  };

} //namespace GeographicLib
#endif

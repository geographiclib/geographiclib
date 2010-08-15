/**
 * \file Geodesic.hpp
 * \brief Header for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEODESIC_HPP)
#define GEOGRAPHICLIB_GEODESIC_HPP "$Id$"

#include "GeographicLib/Constants.hpp"

#if !defined(GEOD_ORD)
/**
 * The order of the expansions used by Geodesic.
 **********************************************************************/
#define GEOD_ORD (GEOGRAPHICLIB_PREC == 1 ? 6 : GEOGRAPHICLIB_PREC == 0 ? 3 : 7)
#endif

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
    typedef Math::real real;
    friend class GeodesicLine;
    static const int nA1 = GEOD_ORD, nC1 = GEOD_ORD, nC1p = GEOD_ORD,
      nA2 = GEOD_ORD, nC2 = GEOD_ORD, nA3 = GEOD_ORD, nC3 = GEOD_ORD,
      nA3x = GEOD_ORD, nC3x = (nA3x * (nA3x - 1)) / 2;
    static const unsigned maxit = 50;

    static inline real sq(real x) throw() { return x * x; }
    void Lengths(real eps, real sig12,
                 real ssig1, real csig1, real ssig2, real csig2,
                 real cbet1, real cbet2,
                 real& s12s, real& m12a, real& m0,
                 real tc[], real zc[]) const throw();
    static real Astroid(real R, real z) throw();
    real InverseStart(real sbet1, real cbet1, real sbet2, real cbet2,
                      real lam12,
                      real& salp1, real& calp1,
                      real& salp2, real& calp2,
                      real C1a[], real C2a[]) const throw();
    real Lambda12(real sbet1, real cbet1, real sbet2, real cbet2,
                  real salp1, real calp1,
                  real& salp2, real& calp2, real& sig12,
                  real& ssig1, real& csig1, real& ssig2, real& csig2,
                  real& eps, bool diffp, real& dlam12,
                  real C1a[], real C2a[], real C3a[])
      const throw();
    void A3coeff() throw();
    void C3coeff() throw();

    static const real eps2, tol0, tol1, tol2, xthresh;
    const real _a, _r, _f, _f1, _e2, _ep2, _n, _b, _etol2;
    real _A3x[nA3x], _C3x[nC3x];
    static real SinSeries(real sinx, real cosx, const real c[], int n)
      throw();
    static inline real AngNormalize(real x) throw() {
      // Place angle in [-180, 180).  Assumes x is in [-540, 540).
      return x >= 180 ? x - 360 : x < -180 ? x + 360 : x;
    }
    static inline real AngRound(real x) throw() {
      // The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
      // for reals = 0.7 pm on the earth if x is an angle in degrees.  (This
      // is about 1000 times more resolution than we get with angles around 90
      // degrees.)  We use this to avoid having to deal with near singular
      // cases when x is non-zero but tiny (e.g., 1.0e-200).
      const real z = real(0.0625); // 1/16
      volatile real y = std::abs(x);
      // The compiler mustn't "simplify" z - (z - y) to y
      y = y < z ? z - (z - y) : y;
      return x < 0 ? -y : y;
    }
    static inline void SinCosNorm(real& sinx, real& cosx) throw() {
      real r = Math::hypot(sinx, cosx);
      sinx /= r;
      cosx /= r;
    }
    // These are Maxima generated functions to provide series approximations to
    // the integrals for the ellipsoidal geodesic.
    static real A1m1f(real eps) throw();
    static void C1f(real eps, real c[]) throw();
    static void C1pf(real eps, real c[]) throw();
    static real A2m1f(real eps) throw();
    static void C2f(real eps, real c[]) throw();
    real A3f(real eps) const throw();
    void C3f(real eps, real c[]) const throw();

  public:

    /**
     * Constructor for a ellipsoid with equatorial radius \e a (meters) and
     * reciprocal flattening \e r.  Setting \e r = 0 implies \e r = inf or
     * flattening = 0 (i.e., a sphere).  Negative \e r indicates a prolate
     * ellipsoid.  An exception is thrown if \e a is not positive.
     **********************************************************************/
    Geodesic(real a, real r);

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
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2, real& azi2, real& m12,
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
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12, real& azi1, real& azi2, real& m12)
      const throw();

    /**
     * Set up to do a series of ranges.  This returns a GeodesicLine object
     * with point 1 given by latitude, \e lat1, longitude, \e lon1, and azimuth
     * \e azi1 (degrees).  Calls to GeodesicLine::Position return the
     * position and azimuth for point 2 a specified distance away.  Using
     * GeodesicLine::Position is approximately 2.1 faster than calling
     * Geodesic::Direct.
     **********************************************************************/
    GeodesicLine Line(real lat1, real lon1, real azi1)
      const throw();

    /**
     * The major radius of the ellipsoid (meters).  This is that value of \e a
     * used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _a; }

    /**
     * The inverse flattening of the ellipsoid.  This is that value of \e r
     * used in the constructor.  A value of 0 is returned for a sphere
     * (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return _r; }

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
   GeographicLib::GeodesicLine
     line(GeographicLib::Geodesic::WGS84.Line(30.0, 10.0, 80.0));
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
   GeographicLib::GeodesicLine line;
   line = GeographicLib::Geodesic::WGS84.Line(30.0, 10.0, 80.0);
   ...
   \endverbatim
   * Similarly, a vector can be used to hold GeodesicLine objects.
   *
   * The calculations are accurate to better than 12 nm.  (See \ref geoderrors
   * for details.)
   **********************************************************************/

  class GeodesicLine {
  private:
    typedef Math::real real;
    friend class Geodesic;
    static const int nC1 = Geodesic::nC1, nC1p = Geodesic::nC1p,
      nC2 = Geodesic::nC2, nC3 = Geodesic::nC3;

    real _lat1, _lon1, _azi1;
    real _a, _r,  _b, _f1, _salp0, _calp0, _k2,
      _ssig1, _csig1, _stau1, _ctau1, _somg1, _comg1,
      _A1m1, _A2m1, _A3c, _B11, _B21, _B31;
    // index zero elements of these arrays are unused
    real _C1a[nC1 + 1], _C1pa[nC1p + 1], _C2a[nC2 + 1], _C3a[nC3];

    static inline real sq(real x) throw() { return x * x; }
    GeodesicLine(const Geodesic& g, real lat1, real lon1, real azi1)
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
    Math::real Position(real s12, real& lat2, real& lon2,
                        real& azi2, real &m12, bool arcmode = false)
      const throw();

    /**
     * Return the scale of the geodesic line extending an arc length \e a12
     * (degrees) from point 1 to point 2.  \e M12 (a number) measures the
     * convergence of initially parallel geodesics.  It is defined by the
     * following construction: starting at point 1 proceed at azimuth \e azi1 +
     * 90<sup>o</sup> a small distance \e dt; turn -90<sup>o</sup> and proceed
     * a distance \e s12 (\e not the arc length \e a12); the distance to point
     * 2 is given by \e M12 \e dt.  \e M21 is defined analogously.
     **********************************************************************/
    void Scale(real a12, real& M12, real& M21) const throw();

    /**
     * Has this object been initialized so that Position can be called?
     **********************************************************************/
    bool Init() const throw() { return _b > 0; }

    /**
     * Return the latitude of point 1 (degrees).
     **********************************************************************/
    Math::real Latitude() const throw() { return Init() ? _lat1 : 0; }

    /**
     * Return the longitude of point 1 (degrees).
     **********************************************************************/
    Math::real Longitude() const throw() { return Init() ? _lon1 : 0; }

    /**
     * Return the azimuth (degrees) of the geodesic line as it passes through
     * point 1.
     **********************************************************************/
    Math::real Azimuth() const throw() { return Init() ? _azi1 : 0; }

    /**
     * Return the azimuth (degrees) of the geodesic line as it crosses the
     * equator in a northward direction.
     **********************************************************************/
    Math::real EquatorialAzimuth() const throw() {
      return Init() ? atan2(_salp0, _calp0) / Constants::degree() : 0;
    }

    /**
     * Return the arc length (degrees) between the northward equatorial
     * crossing and point 1.
     **********************************************************************/
    Math::real EquatorialArc() const throw() {
      return Init() ? atan2(_ssig1, _csig1) / Constants::degree() : 0;
    }

    /**
     * The major radius of the ellipsoid (meters).  This is that value of \e a
     * inherited from the Geodesic object used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return Init() ? _a : 0; }

    /**
     * The inverse flattening of the ellipsoid.  This is that value of \e r
     * inherited from the Geodesic object used in the constructor.  A value of
     * 0 is returned for a sphere (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return Init() ? _r : 0; }
  };

} // namespace GeographicLib
#endif

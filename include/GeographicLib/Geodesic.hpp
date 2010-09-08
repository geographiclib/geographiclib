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
   * The Geodesic class provides the solution of the direct and inverse
   * geodesic problems via Geodesic::Direct and Geodesic::Inverse.  Additional
   * functionality if provided by the GeodesicLine class, which allows a
   * sequence of points along a geodesic to be computed and calculates the
   * geodesic scale and the area under a geodesic via GeodesicLine::Scale and
   * GeodesicLine::Area.
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
      nA2 = GEOD_ORD, nC2 = GEOD_ORD,
      nA3 = GEOD_ORD, nA3x = nA3,
      nC3 = GEOD_ORD, nC3x = (nC3 * (nC3 - 1)) / 2,
      nC4 = GEOD_ORD, nC4x = (nC4 * (nC4 + 1)) / 2;
    static const unsigned maxit = 50;

    static inline real sq(real x) throw() { return x * x; }
    void Lengths(real eps, real sig12,
                 real ssig1, real csig1, real ssig2, real csig2,
                 real cbet1, real cbet2,
                 real& s12s, real& m12a, real& m0,
                 bool scalep, real& M12, real& M21,
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

    static const real eps2, tol0, tol1, tol2, xthresh;
    const real _a, _r, _f, _f1, _e2, _ep2, _n, _b, _c2, _etol2;
    real _A3x[nA3x], _C3x[nC3x], _C4x[nC4x];
    static real SinCosSeries(bool sinp,
                             real sinx, real cosx, const real c[], int n)
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
    void A3coeff() throw();
    real A3f(real eps) const throw();
    void C3coeff() throw();
    void C3f(real eps, real c[]) const throw();
    void C4coeff() throw();
    void C4f(real k2, real c[]) const throw();

    enum captype {
      CAP_NONE = 0U,
      CAP_C1   = 1U<<0,
      CAP_C1p  = 1U<<1,
      CAP_C2   = 1U<<2,
      CAP_C3   = 1U<<3,
      CAP_C4   = 1U<<4,
      CAP_ALL  = 0x1FU,
      OUT_ALL  = 0x7F80U,
    };
  public:

    /**
     * Bit masks for what calculations to do.  These masks do double duty.
     * They signify to the GeographicLib::GeodesicLine::GeodesicLine
     * constructor and to GeographicLib::Geodesic::Line what capabilities
     * should be included in the GeodesicLine object.  They also specify which
     * results to return in the general Geodesic::Direct(real, real, real,
     * bool, real, unsigned, real&, real&, real&, real&, real&, real&, real&,
     * real&) and Geodesic::Inverse routines.
     **********************************************************************/
    enum mask {
      /**
       * No capabilities, no output.
       **********************************************************************/
      NONE          = 0U,
      /**
       * Calculate latitude.  (It's not necessary to include this as a
       * capability to GeodesicLine because this is included by default.)
       **********************************************************************/
      LATITUDE      = 1U<<7  | CAP_NONE,
      /**
       * Calculate longitude.
       **********************************************************************/
      LONGITUDE     = 1U<<8  | CAP_C3,
      /**
       * Calculate azimuth.  (It's not necessary to include this as a
       * capability to GeodesicLine because this is included by default.)
       **********************************************************************/
      AZIMUTH       = 1U<<9  | CAP_NONE,
      /**
       * Calculate distance.
       **********************************************************************/
      DISTANCE      = 1U<<10 | CAP_C1,
      /**
       * Allow distance to be used as input in the direct geodesic problem.
       **********************************************************************/
      DISTANCE_IN   = 1U<<11 | CAP_C1 | CAP_C1p,
      /**
       * Calculate reduced length.
       **********************************************************************/
      REDUCEDLENGTH = 1U<<12 | CAP_C1 | CAP_C2,
      /**
       * Calculate geodesic scale.
       **********************************************************************/
      GEODESICSCALE = 1U<<13 | CAP_C1 | CAP_C2,
      /**
       * Calculate area.
       **********************************************************************/
      AREA          = 1U<<14 | CAP_C4,
      /**
       * All capabilities.  Calculate everything.
       **********************************************************************/
      ALL           = OUT_ALL| CAP_ALL,
    };

    /** \name Constructor
     **********************************************************************/
    ///@{
    /**
     * Constructor for a ellipsoid with equatorial radius \e a (meters) and
     * reciprocal flattening \e r.  Setting \e r = 0 implies \e r = inf or
     * flattening = 0 (i.e., a sphere).  Negative \e r indicates a prolate
     * ellipsoid.  An exception is thrown if \e a is not positive.
     **********************************************************************/
    Geodesic(real a, real r);
    ///@}

    /** \name Direct geodesic problem specified in terms of distance.
     **********************************************************************/
    ///@{
    /**
     * Perform the direct geodesic calculation.  Given a latitude, \e lat1
     * (degrees), longitude, \e lon1 (degrees), and azimuth \e azi1 (degrees)
     * for point 1 and a distance, \e s12 (meters), from point 1 to point 2,
     * return various of
     * - the latitude \e lat2 (degrees) via an argument
     * - the longitude \e lat2 (degrees) via an argument
     * - the azimuth \e azi2 (degrees) via an argument
     * - the reduced length \e m12 (meters) via an argument
     * - the geodesic scale \e M12 (dimensionless) via an argument
     * - the geodesic scale \e M21 (dimensionless) via an argument
     * - the area \e S12 (meters<sup>2</sup>) via an argument
     * - the arc length on the auxiliary sphere \e a12 (degrees), via the
     *   return value of the function
     * .
     * Which quanitities are computed and returned depends on the arguments to
     * the function.  The function return value (the arc length) is always
     * computed.  If either point is at a pole, the azimuth is defined by
     * keeping the longitude fixed and writing \e lat = 90 - \e eps or -90 + \e
     * eps and taking the limit \e eps -> 0 from above.  An arc length greater
     * that 180 degrees results in a geodesic which is not a shortest path.
     * (For a prolate ellipsoid, an additional condition is necessary for a
     * shortest path: the longitudinal extent must not exceed of 180 degrees.)
     **********************************************************************/
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2)
      const throw() {
      real t;
      return GenDirect(lat1, lon1, azi1, false, s12,
                       LATITUDE | LONGITUDE,
                       lat2, lon2, t, t, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Direct.
     **********************************************************************/
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2, real& azi2)
      const throw() {
      real t;
      return GenDirect(lat1, lon1, azi1, false, s12,
                       LATITUDE | LONGITUDE | AZIMUTH,
                       lat2, lon2, azi2, t, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Direct.
     **********************************************************************/
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2, real& azi2, real& m12)
      const throw() {
      real t;
      return GenDirect(lat1, lon1, azi1, false, s12,
                       LATITUDE | LONGITUDE | AZIMUTH | REDUCEDLENGTH,
                       lat2, lon2, azi2, t, m12, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Direct.
     **********************************************************************/
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2, real& azi2,
                      real& M12, real& M21)
      const throw() {
      real t;
      return GenDirect(lat1, lon1, azi1, false, s12,
                       LATITUDE | LONGITUDE | AZIMUTH | GEODESICSCALE,
                       lat2, lon2, azi2, t, t, M12, M21, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Direct.
     **********************************************************************/
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2, real& azi2,
                      real& m12, real& M12, real& M21)
      const throw() {
      real t;
      return GenDirect(lat1, lon1, azi1, false, s12,
                       LATITUDE | LONGITUDE | AZIMUTH |
                       REDUCEDLENGTH | GEODESICSCALE,
                       lat2, lon2, azi2, t, m12, M12, M21, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Direct.
     **********************************************************************/
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2, real& azi2,
                      real& m12, real& M12, real& M21, real& S12)
      const throw() {
      real t;
      return GenDirect(lat1, lon1, azi1, false, s12,
                       LATITUDE | LONGITUDE | AZIMUTH |
                       REDUCEDLENGTH | GEODESICSCALE | AREA,
                       lat2, lon2, azi2, t, m12, M12, M21, S12);
    }
    ///@}

    /** \name Direct geodesic problem specified in terms of arc length.
     **********************************************************************/
    ///@{
    /**
     * Perform the direct geodesic calculation.  Given a latitude, \e lat1
     * (degrees), longitude, \e lon1 (degrees), and azimuth \e azi1 (degrees)
     * for point 1 and an arc length on the auxiliary sphere, \e a12 (degrees),
     * from point 1 to point 2, return various of
     * - the latitude \e lat2 (degrees) via an argument
     * - the longitude \e lat2 (degrees) via an argument
     * - the azimuth \e azi2 (degrees) via an argument
     * - the distance \e s12 (meters) via an argument
     * - the reduced length \e m12 (meters) via an argument
     * - the geodesic scale \e M12 (dimensionless) via an argument
     * - the geodesic scale \e M21 (dimensionless) via an argument
     * - the area \e S12 (meters<sup>2</sup>) via an argument
     * .
     * Which quanitities are computed and returned depends on the arguments to
     * the function.  If either point is at a pole, the azimuth is defined by
     * keeping the longitude fixed and writing \e lat = 90 - \e eps or -90 + \e
     * eps and taking the limit \e eps -> 0 from above.  An arc length greater
     * that 180 degrees results in a geodesic which is not a shortest path.
     * (For a prolate ellipsoid, an additional condition is necessary for a
     * shortest path: the longitudinal extent must not exceed of 180 degrees.)
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2) const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE,
                lat2, lon2, t, t, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2, real& azi2) const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE | AZIMUTH,
                lat2, lon2, azi2, t, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2, real& azi2, real& s12)
      const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE | AZIMUTH | DISTANCE,
                lat2, lon2, azi2, s12, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2, real& azi2,
                   real& s12, real& m12) const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE | AZIMUTH | DISTANCE |
                REDUCEDLENGTH,
                lat2, lon2, azi2, s12, m12, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2, real& azi2, real& s12,
                   real& M12, real& M21) const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE | AZIMUTH | DISTANCE |
                GEODESICSCALE,
                lat2, lon2, azi2, s12, t, M12, M21, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2, real& azi2, real& s12,
                   real& m12, real& M12, real& M21) const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE | AZIMUTH | DISTANCE |
                REDUCEDLENGTH | GEODESICSCALE,
                lat2, lon2, azi2, s12, m12, M12, M21, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2, real& azi2, real& s12,
                   real& m12, real& M12, real& M21, real& S12)
      const throw() {
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE | AZIMUTH | DISTANCE |
                REDUCEDLENGTH | GEODESICSCALE | AREA,
                lat2, lon2, azi2, s12, m12, M12, M21, S12);
    }
    ///@}

    /** \name General version of the direct geodesic solution.
     **********************************************************************/
    ///@{

    /**
     * Generic direct geodesic calculation.
     **********************************************************************/
    Math::real GenDirect(real lat1, real lon1, real azi1,
                         bool arcmode, real a12, unsigned outmask,
                         real& lat2, real& lon2, real& azi2,
                         real& s12, real& m12, real& M12, real& M21,
                         real& S12) const throw();
    ///@}

    /** \name Inverse geodesic problem.
     **********************************************************************/
    ///@{
    /**
     * Perform the inverse geodesic calculation.  Given a latitude, \e lat1,
     * longitude, \e lon1, for point 1 and a latitude, \e lat2, longitude, \e
     * lon2, for point 2 (all in degrees), return the geodesic distance, \e s12
     * (meters), the forward azimuths, \e azi1 and \e azi2 (degrees) at points
     * 1 and 2, and the reduced length \e m12 (meters).  Returned value is the
     * arc length \e a12 (degrees) on the auxiliary sphere.  The routine uses
     * an iterative method.  If the method fails to converge, FIX FIX FIX the
     * negative of the distances (\e s12, \e m12, and \e a12) and reverse of
     * the azimuths are returned.  This is not expected to happen with
     * ellipsoidal models of the earth.  Please report all cases where this
     * occurs.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12) const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE,
                        s12, t, t, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& azi1, real& azi2) const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        AZIMUTH,
                        t, azi1, azi2, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12, real& azi1, real& azi2)
      const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE | AZIMUTH,
                        s12, azi1, azi2, t, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12, real& azi1, real& azi2, real& m12)
      const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE | AZIMUTH | REDUCEDLENGTH,
                        s12, azi1, azi2, m12, t, t, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12, real& azi1, real& azi2,
                       real& M12, real& M21) const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE | AZIMUTH | GEODESICSCALE,
                        s12, azi1, azi2, t, M12, M21, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12, real& azi1, real& azi2, real& m12,
                       real& M12, real& M21) const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE | AZIMUTH |
                        REDUCEDLENGTH | GEODESICSCALE,
                        s12, azi1, azi2, m12, M12, M21, t);
    }

    /**
     * See the documentation for GeographicLib::Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12, real& azi1, real& azi2, real& m12,
                       real& M12, real& M21, real& S12) const throw() {
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE | AZIMUTH |
                        REDUCEDLENGTH | GEODESICSCALE | AREA,
                        s12, azi1, azi2, m12, M12, M21, S12);
    }

    ///@}

    /** \name General version of inverse geodesic solution.
     **********************************************************************/
    ///@{
    /**
     * General version of inverse geodesic calculation.
     **********************************************************************/
    Math::real GenInverse(real lat1, real lon1, real lat2, real lon2,
                       unsigned outmask,
                       real& s12, real& azi1, real& azi2,
                       real& m12, real& M12, real& M21, real& S12)
      const throw();
    ///@}

    /** \name Interface to GeographicLib::GeodesicLine.
     **********************************************************************/
    ///@{

    /**
     * Set up to do a series of ranges.  This returns a GeodesicLine object
     * with point 1 given by latitude, \e lat1, longitude, \e lon1, and azimuth
     * \e azi1 (degrees).  Calls to GeodesicLine::Position return the
     * position and azimuth for point 2 a specified distance away.  Using
     * GeodesicLine::Position is approximately 2.1 faster than calling
     * Geodesic::Direct.
     **********************************************************************/
    GeodesicLine Line(real lat1, real lon1, real azi1, unsigned caps = ALL)
      const throw();

    ///@}

    /** \name Inspector functions.
     **********************************************************************/
    ///@{

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
     * Return total area of ellipsoid in meters<sup>2</sup>.  (Does not require
     * GeodesicLine::AreaEnable to have been called.)  The area of a polygon
     * encircling a pole can be found by adding GeodesicLine::EllipsoidArea()/2
     * to the sum of GeodesicLine::Area for each side of the polygon.
     **********************************************************************/
    Math::real EllipsoidArea() const throw() {
      return 4 * Constants::pi() * _c2;
    }
    ///@}

    /**
     * A global instantiation of Geodesic with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    const static Geodesic WGS84;


    /** \name Deprecated function.
     **********************************************************************/
    ///@{
    /**
     * <b>DEPRECATED</b> Perform the direct geodesic calculation.  Given a
     * latitude, \e lat1, longitude, \e lon1, and azimuth \e azi1 (degrees) for
     * point 1 and a range, \e s12 (meters) from point 1 to point 2, return the
     * latitude, \e lat2, longitude, \e lon2, and forward azimuth, \e azi2
     * (degrees) for point 2 and the reduced length \e m12 (meters).  If either
     * point is at a pole, the azimuth is defined by keeping the longitude
     * fixed and writing \e lat = 90 - \e eps or -90 + \e eps and taking the
     * limit \e eps -> 0 from above.  If \e arcmode (default false) is set to
     * true, \e s12 is interpreted as the arc length \e a12 (degrees) on the
     * auxiliary sphere.  An arc length greater that 180 degrees results in a
     * geodesic which is not a shortest path.  For a prolate ellipsoid, an
     * additional condition is necessary for a shortest path: the longitudinal
     * extent must not exceed of 180 degrees.  Returned value is the arc length
     * \e a12 (degrees) if \e arcmode is false, otherwise it is the distance \e
     * s12 (meters).
     **********************************************************************/
    Math::real Direct(real lat1, real lon1, real azi1, real s12,
                      real& lat2, real& lon2, real& azi2, real& m12,
                      bool arcmode) const throw() {
      if (arcmode) {
        real s12x;
        ArcDirect(lat1, lon1, azi1, s12, lat2, lon2, azi2, s12x, m12);
        return s12x;
      } else
        return Direct(lat1, lon1, azi1, s12, lat2, lon2, azi2, m12);
    }
    ///@}
  };

} // namespace GeographicLib
#endif

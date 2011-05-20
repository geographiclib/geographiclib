/**
 * \file Geodesic.hpp
 * \brief Header for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEODESIC_HPP)
#define GEOGRAPHICLIB_GEODESIC_HPP "$Id: Geodesic.hpp 6867 2010-09-11 13:04:26Z karney $"

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
   * Given \e lat1, \e lon1, \e azi1, and \e s12, we can determine \e lat2, \e
   * lon2, and \e azi2.  This is the \e direct geodesic problem and its
   * solution is given by the function Geodesic::Direct.  (If \e s12 is
   * sufficiently large that the geodesic wraps more than halfway around the
   * earth, there will be another geodesic between the points with a smaller \e
   * s12.)
   *
   * Given \e lat1, \e lon1, \e lat2, and \e lon2, we can determine \e azi1, \e
   * azi2, and \e s12.  This is the \e inverse geodesic problem, whose solution
   * is given by Geodesic::Inverse.  Usually, the solution to the inverse
   * problem is unique.  In cases where there are muliple solutions (all with
   * the same \e s12, of course), all the solutions can be easily generated
   * once a particular solution is provided.
   *
   * The standard way of specifying the direct problem is the specify the
   * distance \e s12 to the second point.  However it is sometimes useful
   * instead to specify the the arc length \e a12 (in degrees) on the auxiliary
   * sphere.  This is a mathematical construct used in solving the geodesic
   * problems.  The solution of the direct problem in this form is provide by
   * Geodesic::ArcDirect.  An arc length in excess of 180<sup>o</sup> indicates
   * that the geodesic is not a shortest path.  In addition, the arc length
   * between an equatorial crossing and the next extremum of latitude for a
   * geodesic is 90<sup>o</sup>.
   *
   * This class can also calculate several other quantities related to
   * geodesics.  These are:
   * - <i>reduced length</i>.  If we fix the first point and increase \e azi1
   *   by \e dazi1 (radians), the the second point is displaced \e m12 \e dazi1
   *   in the direction \e azi2 + 90<sup>o</sup>.  The quantity \e m12 is
   *   called the "reduced length" and is symmetric under interchange of the
   *   two points.  On a flat surface, we have \e m12 = \e s12.  The ratio \e
   *   s12/\e m12 gives the azimuthal scale for an azimuthal equidistant
   *   projection.
   * - <i>geodesic scale</i>.  Consider a reference geodesic and a second
   *   geodesic parallel to this one at point 1 and separated by a small
   *   distance \e dt.  The separation of the two geodesics at point 2 is \e
   *   M12 \e dt where \e M12 is called the "geodesic scale".  \e M21 is
   *   defined similarly (with the geodesics being parallel at point 2).  On a
   *   flat surface, we have \e M12 = \e M21 = 1.  The quantity 1/\e M12 gives
   *   the scale of the Cassini-Soldner projection.
   * - <i>area</i>.  Consider the quadrilateral bounded by the following lines:
   *   the geodesic from point 1 to point 2, the meridian from point 2 to the
   *   equator, the equator from \e lon2 to \e lon1, the meridian from the
   *   equator to point 1.  The area of this quadrilateral is represented by \e
   *   S12 with a clockwise traversal of the perimeter counting as a positive
   *   area and it can be used to compute the area of any simple geodesic
   *   polygon.
   *
   * Overloaded versions of Geodesic::Direct, Geodesic::ArcDirect, and
   * Geodesic::Inverse allow these quantities to be returned.  In addition
   * there are general functions Geodesic::GenDirect, and Geodesic::GenInverse
   * which allow an arbitrary set of results to be computed.
   *
   * Additional functionality if provided by the GeodesicLine class, which
   * allows a sequence of points along a geodesic to be computed.
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
     * They signify to the GeodesicLine::GeodesicLine constructor and to
     * Geodesic::Line what capabilities should be included in the GeodesicLine
     * object.  They also specify which results to return in the general
     * routines Geodesic::GenDirect and Geodesic::GenInverse routines.
     * GeodesicLine::mask is a duplication of this enum.
     **********************************************************************/
    enum mask {
      /**
       * No capabilities, no output.
       * @hideinitializer
       **********************************************************************/
      NONE          = 0U,
      /**
       * Calculate latitude \e lat2.  (It's not necessary to include this as a
       * capability to GeodesicLine because this is included by default.)
       * @hideinitializer
       **********************************************************************/
      LATITUDE      = 1U<<7  | CAP_NONE,
      /**
       * Calculate longitude \e lon2.
       * @hideinitializer
       **********************************************************************/
      LONGITUDE     = 1U<<8  | CAP_C3,
      /**
       * Calculate azimuths \e azi1 and \e azi2.  (It's not necessary to
       * include this as a capability to GeodesicLine because this is included
       * by default.)
       * @hideinitializer
       **********************************************************************/
      AZIMUTH       = 1U<<9  | CAP_NONE,
      /**
       * Calculate distance \e s12.
       * @hideinitializer
       **********************************************************************/
      DISTANCE      = 1U<<10 | CAP_C1,
      /**
       * Allow distance \e s12 to be used as input in the direct geodesic
       * problem.
       * @hideinitializer
       **********************************************************************/
      DISTANCE_IN   = 1U<<11 | CAP_C1 | CAP_C1p,
      /**
       * Calculate reduced length \e m12.
       * @hideinitializer
       **********************************************************************/
      REDUCEDLENGTH = 1U<<12 | CAP_C1 | CAP_C2,
      /**
       * Calculate geodesic scales \e M12 and \e M21.
       * @hideinitializer
       **********************************************************************/
      GEODESICSCALE = 1U<<13 | CAP_C1 | CAP_C2,
      /**
       * Calculate area \e S12.
       * @hideinitializer
       **********************************************************************/
      AREA          = 1U<<14 | CAP_C4,
      /**
       * All capabilities.  Calculate everything.
       * @hideinitializer
       **********************************************************************/
      ALL           = OUT_ALL| CAP_ALL,
    };

    /** \name Constructor
     **********************************************************************/
    ///@{
    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters)
     * @param[in] r reciprocal flattening.  Setting \e r = 0 implies \e r = inf
     *   or flattening = 0 (i.e., a sphere).  Negative \e r indicates a prolate
     *   ellipsoid.
     *
     * An exception is thrown if either of the axes of the ellipsoid is
     * non-positive.
     **********************************************************************/
    Geodesic(real a, real r);
    ///@}

    /** \name Direct geodesic problem specified in terms of distance.
     **********************************************************************/
    ///@{
    /**
     * Perform the direct geodesic calculation where the length of the geodesic
     * is specify in terms of distance.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi1 azimuth at point 1 (degrees).
     * @param[in] s12 distance between point 1 and point 2 (meters); it can be
     *   signed.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees).
     * @param[out] azi2 (forward) azimuth at point 2 (degrees).
     * @param[out] m12 reduced length of geodesic (meters).
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless).
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless).
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>).
     * @return \e a12 arc length of between point 1 and point 2 (degrees).
     *
     * If either point is at a pole, the azimuth is defined by keeping the
     * longitude fixed and writing \e lat = 90 - \e eps or -90 + \e eps and
     * taking the limit \e eps -> 0 from above.  An arc length greater that 180
     * degrees signifies a geodesic which is not a shortest path.  (For a
     * prolate ellipsoid, an additional condition is necessary for a shortest
     * path: the longitudinal extent must not exceed of 180 degrees.)
     *
     * The following functions are overloaded versions of Geodesic::Direct
     * which omit some of the output parameters.  Note, however, that the arc
     * length is always computed and returned as the function value.
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

    /**
     * See the documentation for Geodesic::Direct.
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
     * See the documentation for Geodesic::Direct.
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
     * See the documentation for Geodesic::Direct.
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
     * See the documentation for Geodesic::Direct.
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
     * See the documentation for Geodesic::Direct.
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
    ///@}

    /** \name Direct geodesic problem specified in terms of arc length.
     **********************************************************************/
    ///@{
    /**
     * Perform the direct geodesic calculation where the length of the geodesic
     * is specify in terms of arc length.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi1 azimuth at point 1 (degrees).
     * @param[in] a12 arc length between point 1 and point 2 (degrees); it can
     *   be signed.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees).
     * @param[out] azi2 (forward) azimuth at point 2 (degrees).
     * @param[out] s12 distance between point 1 and point 2 (meters).
     * @param[out] m12 reduced length of geodesic (meters).
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless).
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless).
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>).
     *
     * If either point is at a pole, the azimuth is defined by keeping the
     * longitude fixed and writing \e lat = 90 - \e eps or -90 + \e eps and
     * taking the limit \e eps -> 0 from above.  An arc length greater that 180
     * degrees signifies a geodesic which is not a shortest path.  (For a
     * prolate ellipsoid, an additional condition is necessary for a shortest
     * path: the longitudinal extent must not exceed of 180 degrees.)
     *
     * The following functions are overloaded versions of Geodesic::Direct
     * which omit some of the output parameters.
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

    /**
     * See the documentation for Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2) const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE,
                lat2, lon2, t, t, t, t, t, t);
    }

    /**
     * See the documentation for Geodesic::ArcDirect.
     **********************************************************************/
    void ArcDirect(real lat1, real lon1, real azi1, real a12,
                   real& lat2, real& lon2, real& azi2) const throw() {
      real t;
      GenDirect(lat1, lon1, azi1, true, a12,
                LATITUDE | LONGITUDE | AZIMUTH,
                lat2, lon2, azi2, t, t, t, t, t);
    }

    /**
     * See the documentation for Geodesic::ArcDirect.
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
     * See the documentation for Geodesic::ArcDirect.
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
     * See the documentation for Geodesic::ArcDirect.
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
     * See the documentation for Geodesic::ArcDirect.
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
    ///@}

    /** \name General version of the direct geodesic solution.
     **********************************************************************/
    ///@{

    /**
     * The general direct geodesic calculation.  Geodesic::Direct and
     * Geodesic::ArcDirect are defined in terms of this function.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi1 azimuth at point 1 (degrees).
     * @param[in] arcmode boolean flag determining the meaning of the second
     *   parameter.
     * @param[in] s12_a12 if \e arcmode is false, this is the distance between
     *   point 1 and point 2 (meters); otherwise it is the arc length between
     *   point 1 and point 2 (degrees); it can be signed.
     * @param[in] outmask a bitor'ed combination of Geodesic::mask values
     *   specifying which of the following parameters should be set.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees).
     * @param[out] azi2 (forward) azimuth at point 2 (degrees).
     * @param[out] s12 distance between point 1 and point 2 (meters).
     * @param[out] m12 reduced length of geodesic (meters).
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless).
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless).
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>).
     * @return \e a12 arc length of between point 1 and point 2 (degrees).
     *
     * The Geodesic::mask values possible for \e outmask are
     * - \e outmask |= Geodesic::LATITUDE for the latitude \e lat2.
     * - \e outmask |= Geodesic::LONGITUDE for the latitude \e lon2.
     * - \e outmask |= Geodesic::AZIMUTH for the latitude \e azi2.
     * - \e outmask |= Geodesic::DISTANCE for the distance \e s12.
     * - \e outmask |= Geodesic::REDUCEDLENGTH for the reduced length \e
     *   m12.
     * - \e outmask |= Geodesic::GEODESICSCALE for the geodesic scales \e
     *   M12 and \e M21.
     * - \e outmask |= Geodesic::AREA for the area \e S12.
     * .
     * The function value \e a12 is always computed and returned and this
     * equals \e s12_a12 is \e arcmode is true.  If \e outmask includes
     * Geodesic::DISTANCE and \e arcmode is false, then \e s12 = \e s12_a12.
     * It is not necessary to include Geodesic::DISTANCE_IN in \e outmask; this
     * is automatically included is \e arcmode is false.
     **********************************************************************/
    Math::real GenDirect(real lat1, real lon1, real azi1,
                         bool arcmode, real s12_a12, unsigned outmask,
                         real& lat2, real& lon2, real& azi2,
                         real& s12, real& m12, real& M12, real& M21,
                         real& S12) const throw();
    ///@}

    /** \name Inverse geodesic problem.
     **********************************************************************/
    ///@{
    /**
     * Perform the inverse geodesic calculation.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] lat2 latitude of point 2 (degrees).
     * @param[in] lon2 longitude of point 2 (degrees).
     * @param[out] s12 distance between point 1 and point 2 (meters).
     * @param[out] azi1 azimuth at point 1 (degrees).
     * @param[out] azi2 (forward) azimuth at point 1 (degrees).
     * @param[out] m12 reduced length of geodesic (meters).
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless).
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless).
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>).
     * @return \e a12 arc length of between point 1 and point 2 (degrees).
     *
     * If either point is at a pole, the azimuth is defined by keeping the
     * longitude fixed and writing \e lat = 90 - \e eps or -90 + \e eps and
     * taking the limit \e eps -> 0 from above.  If the routine fails to
     * converge, then all the requested outputs are set to Math::NaN().  This
     * is not expected to happen with ellipsoidal models of the earth; please
     * report all cases where this occurs.
     *
     * The following functions are overloaded versions of Geodesic::Inverse
     * which omit some of the output parameters.  Note, however, that the arc
     * length is always computed and returned as the function value.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12, real& azi1, real& azi2, real& m12,
                       real& M12, real& M21, real& S12) const throw() {
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE | AZIMUTH |
                        REDUCEDLENGTH | GEODESICSCALE | AREA,
                        s12, azi1, azi2, m12, M12, M21, S12);
    }

    /**
     * See the documentation for Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& s12) const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        DISTANCE,
                        s12, t, t, t, t, t, t);
    }

    /**
     * See the documentation for Geodesic::Inverse.
     **********************************************************************/
    Math::real Inverse(real lat1, real lon1, real lat2, real lon2,
                       real& azi1, real& azi2) const throw() {
      real t;
      return GenInverse(lat1, lon1, lat2, lon2,
                        AZIMUTH,
                        t, azi1, azi2, t, t, t, t);
    }

    /**
     * See the documentation for Geodesic::Inverse.
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
     * See the documentation for Geodesic::Inverse.
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
     * See the documentation for Geodesic::Inverse.
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
     * See the documentation for Geodesic::Inverse.
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
    ///@}

    /** \name General version of inverse geodesic solution.
     **********************************************************************/
    ///@{
    /**
     * The general inverse geodesic calculation.  Geodesic::Inverse is defined
     * in terms of this function.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] lat2 latitude of point 2 (degrees).
     * @param[in] lon2 longitude of point 2 (degrees).
     * @param[in] outmask a bitor'ed combination of Geodesic::mask values
     *   specifying which of the following parameters should be set.
     * @param[out] s12 distance between point 1 and point 2 (meters).
     * @param[out] azi1 azimuth at point 1 (degrees).
     * @param[out] azi2 (forward) azimuth at point 1 (degrees).
     * @param[out] m12 reduced length of geodesic (meters).
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless).
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless).
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>).
     * @return \e a12 arc length of between point 1 and point 2 (degrees).
     *
     * The Geodesic::mask values possible for \e outmask are
     * - \e outmask |= Geodesic::DISTANCE for the distance \e s12.
     * - \e outmask |= Geodesic::AZIMUTH for the latitude \e azi2.
     * - \e outmask |= Geodesic::REDUCEDLENGTH for the reduced length \e
     *   m12.
     * - \e outmask |= Geodesic::GEODESICSCALE for the geodesic scales \e
     *   M12 and \e M21.
     * - \e outmask |= Geodesic::AREA for the area \e S12.
     * .
     * The arc length is always computed and returned as the function value.
     **********************************************************************/
    Math::real GenInverse(real lat1, real lon1, real lat2, real lon2,
                       unsigned outmask,
                       real& s12, real& azi1, real& azi2,
                       real& m12, real& M12, real& M21, real& S12)
      const throw();
    ///@}

    /** \name Interface to GeodesicLine.
     **********************************************************************/
    ///@{

    /**
     * Set up to do a series of ranges.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi1 azimuth at point 1 (degrees).
     * @param[in] caps bitor'ed combination of Geodesic::mask values
     *   specifying the capabilities the GeodesicLine object should possess,
     *   i.e., which quantities can be returned in calls to
     *   GeodesicLib::Position.
     *
     * The Geodesic::mask values are
     * - \e caps |= Geodesic::LATITUDE for the latitude \e lat2; this is
     *   added automatically
     * - \e caps |= Geodesic::LONGITUDE for the latitude \e lon2
     * - \e caps |= Geodesic::AZIMUTH for the latitude \e azi2; this is
     *   added automatically
     * - \e caps |= Geodesic::DISTANCE for the distance \e s12
     * - \e caps |= Geodesic::REDUCEDLENGTH for the reduced length \e m12
     * - \e caps |= Geodesic::GEODESICSCALE for the geodesic scales \e M12
     *   and \e M21
     * - \e caps |= Geodesic::AREA for the area \e S12
     * - \e caps |= Geodesic::DISTANCE_IN permits the length of the
     *   geodesic to be given in terms of \e s12; without this capability the
     *   length can only be specified in terms of arc length.
     * .
     * The default value of \e caps is Geodesic::ALL which turns on all the
     * capabilities.
     *
     * If the point is at a pole, the azimuth is defined by keeping the \e lon1
     * fixed and writing \e lat1 = 90 - \e eps or -90 + \e eps and taking the
     * limit \e eps -> 0 from above.
     **********************************************************************/
    GeodesicLine Line(real lat1, real lon1, real azi1, unsigned caps = ALL)
      const throw();

    ///@}

    /** \name Inspector functions.
     **********************************************************************/
    ///@{

    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _a; }

    /**
     * @return \e r the inverse flattening of the ellipsoid.  This is the
     *   value used in the constructor.  A value of 0 is returned for a sphere
     *   (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return _r; }

    /**
     * @return total area of ellipsoid in meters<sup>2</sup>.  The area of a
     *   polygon encircling a pole can be found by adding
     *   Geodesic::EllipsoidArea()/2 to the sum of \e S12 for each side of the
     *   polygon.
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
    Math::real Direct(real lat1, real lon1, real azi1, real s12_a12,
                      real& lat2, real& lon2, real& azi2, real& m12,
                      bool arcmode) const throw() {
      if (arcmode) {
        real a12 = s12_a12, s12;
        ArcDirect(lat1, lon1, azi1, a12, lat2, lon2, azi2, s12, m12);
        return s12;
      } else {
        real s12 = s12_a12;
        return Direct(lat1, lon1, azi1, s12, lat2, lon2, azi2, m12);
      }
    }
    ///@}
  };

} // namespace GeographicLib
#endif

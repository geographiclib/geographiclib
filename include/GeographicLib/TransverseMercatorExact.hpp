/**
 * \file TransverseMercatorExact.hpp
 * \brief Header for GeographicLib::TransverseMercatorExact class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRANSVERSEMERCATOREXACT_HPP)
#define GEOGRAPHICLIB_TRANSVERSEMERCATOREXACT_HPP "$Id: TransverseMercatorExact.hpp 6867 2010-09-11 13:04:26Z karney $"

#include "GeographicLib/Constants.hpp"
#include "GeographicLib/EllipticFunction.hpp"

namespace GeographicLib {

  /**
   * \brief An exact implementation of the Transverse Mercator Projection
   *
   * Implementation of the Transverse Mercator Projection given in
   *  - L. P. Lee,
   *    <a href="http://dx.doi.org/10.3138/X687-1574-4325-WM62"> Conformal
   *    Projections Based On Jacobian Elliptic Functions</a>, Part V of
   *    Conformal Projections Based on Elliptic Functions,
   *    (B. V. Gutsell, Toronto, 1976), 128pp.,
   *    ISBN: 0919870163
   *    (also appeared as:
   *    Monograph 16, Suppl. No. 1 to Canadian Cartographer, Vol 13).
   *
   * This method gives the correct results for forward and reverse
   * transformations subject to the branch cut rules (see the description of
   * the \e extendp argument to the constructor).  The maximum error is about 8
   * nm (ground distance) for the forward and reverse transformations.  The
   * error in the convergence is 2e-15&quot;, the relative error in the scale
   * is 7e-12%%.  (See \ref tmerrors for the weasel words.)  The method is
   * "exact" in the sense that the errors are close to the round-off limit and
   * that no changes are needed in the algorithms for them to be used with
   * reals of a higher precision.  Thus the errors using long double (with a
   * 64-bit fraction) are about 2000 times smaller than using double (with a
   * 53-bit fraction).
   *
   * This algorithm is about 4.5 times slower than the 6th-order Kr&uuml;ger
   * method, TransverseMercator, taking about 11 us for a combined forward and
   * reverse projection on a 2.66 GHz Intel machine (g++, version 4.3.0, -O3).
   *
   * The ellipsoid parameters and the central scale are set in the constructor.
   * The central meridian (which is a trivial shift of the longitude) is
   * specified as the \e lon0 argument of the TransverseMercatorExact::Forward
   * and TransverseMercatorExact::Reverse functions.  The latitude of origin is
   * taken to be the equator.  See the documentation on TransverseMercator for
   * how to include a false easting, false northing, or a latitude of origin.
   *
   * See TransverseMercatorExact.cpp for more information on the
   * implementation.
   *
   * See \ref transversemercator for a discussion of this projection.
   **********************************************************************/

  class TransverseMercatorExact {
  private:
    typedef Math::real real;
    static const real tol, tol1, tol2, taytol, overflow;
    static const int numit = 10;
    const real _a, _r, _f, _k0, _mu, _mv, _e, _ep2;
    const bool _extendp;
    const EllipticFunction _Eu, _Ev;
    static inline real sq(real x) throw() { return x * x; }
    // tan(x) for x in [-pi/2, pi/2] ensuring that the sign is right
    static inline real tanx(real x) throw() {
      real t = std::tan(x);
      return x >= 0 ? (t >= 0 ? t : overflow) : (t < 0 ? t : -overflow);
    }

    real taup(real tau) const throw();
    real taupinv(real taup) const throw();

    void zeta(real u, real snu, real cnu, real dnu,
              real v, real snv, real cnv, real dnv,
              real& taup, real& lam) const throw();

    void dwdzeta(real u, real snu, real cnu, real dnu,
                 real v, real snv, real cnv, real dnv,
                 real& du, real& dv) const throw();

    bool zetainv0(real psi, real lam, real& u, real& v) const throw();
    void zetainv(real taup, real lam, real& u, real& v) const throw();

    void sigma(real u, real snu, real cnu, real dnu,
               real v, real snv, real cnv, real dnv,
               real& xi, real& eta) const throw();

    void dwdsigma(real u, real snu, real cnu, real dnu,
                  real v, real snv, real cnv, real dnv,
                  real& du, real& dv) const throw();

    bool sigmainv0(real xi, real eta, real& u, real& v) const throw();
    void sigmainv(real xi, real eta, real& u, real& v) const throw();

    void Scale(real tau, real lam,
               real snu, real cnu, real dnu,
               real snv, real cnv, real dnv,
               real& gamma, real& k) const throw();

  public:

    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters)
     * @param[in] r reciprocal flattening.
     * @param[in] k0 central scale factor.
     * @param[in] extendp use extended domain.
     *
     * The transverse Mercator projection has a branch point singularity at \e
     * lat = 0 and \e lon - \e lon0 = 90 (1 - \e e) or (for
     * TransverseMercatorExact::UTM) x = 18381 km, y = 0m.  The \e extendp
     * argument governs where the branch cut is placed.  With \e extendp =
     * false, the "standard" convention is followed, namely the cut is placed
     * along x > 18381 km, y = 0m.  Forward can be called with any \e lat and
     * \e lon then produces the transformation shown in Lee, Fig 46.  Reverse
     * analytically continues this in the +/- \e x direction.  As a
     * consequence, Reverse may map multiple points to the same geographic
     * location; for example, for TransverseMercatorExact::UTM, \e x =
     * 22051449.037349 m, \e y = -7131237.022729 m and \e x = 29735142.378357
     * m, \e y = 4235043.607933 m both map to \e lat = -2 deg, \e lon = 88 deg.
     *
     * With \e extendp = true, the branch cut is moved to the lower left
     * quadrant.  The various symmetries of the transverse Mercator projection
     * can be used to explore the projection on any sheet.  In this mode the
     * domains of \e lat, \e lon, \e x, and \e y are restricted to
     * - the union of
     *   - \e lat in [0, 90] and \e lon - \e lon0 in [0, 90]
     *   - \e lat in (-90, 0] and \e lon - \e lon0 in [90 (1 - \e e), 90]
     * - the union of
     *   - \e x/(\e k0 \e a) in [0, inf) and
     *     \e y/(\e k0 \e a) in [0, E(\e e^2)]
     *   - \e x/(\e k0 \e a) in [K(1 - \e e^2) - E(1 - \e e^2), inf) and
     *     \e y/(\e k0 \e a) in (-inf, 0]
     * .
     * See \ref extend for a full discussion of the treatment of the branch
     * cut.
     *
     * The method will work for all ellipsoids used in terrestial geodesy.  The
     * method cannot be applied directly to the case of a sphere (\e r = inf)
     * because some the constants characterizing this method diverge in that
     * limit, and in practise, \e r should be smaller than about
     * 1/numeric_limits< real >::%epsilon().  However, TransverseMercator
     * treats the sphere exactly.  An exception is thrown if either axis of the
     * ellipsoid or \e k0 is not positive or if \e r < 1.
     **********************************************************************/
    TransverseMercatorExact(real a, real r, real k0, bool extendp = false);

    /**
     * Forward projection, from geographic to transverse Mercator.
     *
     * @param[in] lon0 central meridian of the projection (degrees).
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added. \e lat should be in the range
     * [-90, 90]; \e lon and \e lon0 should be in the range [-180, 360].
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y, real& gamma, real& k) const throw();

    /**
     * Reverse projection, from transverse Mercator to geographic.
     *
     * @param[in] lon0 central meridian of the projection (degrees).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added.  \e lon0 should be in the range
     * [-180, 360].  The value of \e lon returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon, real& gamma, real& k) const throw();

    /**
     * TransverseMercatorExact::Forward without returning the convergence and
     * scale.
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y) const throw() {
      real gamma, k;
      Forward(lon0, lat, lon, x, y, gamma, k);
    }

    /**
     * TransverseMercatorExact::Reverse without returning the convergence and
     * scale.
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon) const throw() {
      real gamma, k;
      Reverse(lon0, x, y, lat, lon, gamma, k);
    }

    /** \name Inspector functions
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
     * @return \e k0 central scale for the projection.  This is the value of \e
     *   k0 used in the constructor and is the scale on the central meridian.
     **********************************************************************/
    Math::real CentralScale() const throw() { return _k0; }
    ///@}

    /**
     * A global instantiation of TransverseMercatorExact with the WGS84
     * ellipsoid and the UTM scale factor.  However, unlike UTM, no false
     * easting or northing is added.
     **********************************************************************/
    const static TransverseMercatorExact UTM;
  };

} // namespace GeographicLib

#endif

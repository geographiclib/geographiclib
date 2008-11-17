/**
 * \file TransverseMercatorExact.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(TRANSVERSEMERCATOREXACT_HPP)
#define TRANSVERSEMERCATOREXACT_HPP "$Id$"

#include <cmath>
#include "GeographicLib/EllipticFunction.hpp"

namespace GeographicLib {

  /**
   * \brief An exact implementation of the Transverse Mercator Projection
   *
   * Implementation of the Transverse Mercator Projection given in
   *  - L. P. Lee,
   *    Conformal Projections Based on Elliptic Functions,
   *    (B. V. Gutsell, Toronto, 1976), 128pp.,
   *    ISBN: 0919870163.
   *    [Also appeared as:
   *    Monograph 16, Suppl. No. 1 to Canadian Cartographer, Vol 13.]
   *
   * The method entails using the Thompson Transverse Mercator as an
   * intermediate projection.  The projections from the intermediate
   * coordinates to [phi, lam] and [x, y] are given by elliptic functions.  The
   * inverse of these projections are found by Newton's method with a suitable
   * starting guess.
   *
   * This method gives the correct results for forward and reverse
   * transformations subject to the branch cut rules (see the description of
   * the \e foldp argument to the constructor..  The maximum error is about
   * 10nm (ground distance) for the forward and reverse transformations.  The
   * error in the convergence is 5e-9", the relative error in the scale is
   * 2e-11%%.  The method is "exact" in the sense that the errors are close to
   * the round-off limit and that no changes are needed in the algorithms for
   * them to be used with reals of a higher precision (e.g., long double).
   *
   * This algorithm is about 2.5 times slower than the 6th-order Krueger method
   * taking about 15us for a combined forward and reverse projection on a
   * 2.6GHz Intel machine (g++, version 4.3.0, -O3).
   *
   * This implementation and notation closely follows Lee, with the following
   * exceptions:
   * <center><table>
   * <tr><th>Lee    <th>here    <th>Description
   * <tr><td>x/a    <td>xi      <td>Northing (unit Earth)
   * <tr><td>y/a    <td>eta     <td>Easting (unit Earth)
   * <tr><td>s/a    <td>sigma   <td>xi + i * eta
   * <tr><td>y      <td>x       <td>Easting
   * <tr><td>x      <td>y       <td>Northing
   * <tr><td>k      <td>e       <td>eccentricity
   * <tr><td>k^2    <td>mu      <td>elliptic function parameter
   * <tr><td>k'^2   <td>mv      <td>elliptic function complementary parameter
   * <tr><td>m      <td>k       <td>scale
   * </table></center>
   *
   * Minor alterations have been made in some of Lee's expressions in an
   * attempt to control round-off.  For example atanh(sin(phi)) is replaced by
   * asinh(tan(phi)) which maintains accuracy near phi = pi/2.  Such changes
   * are noted in the code.
   *
   * Loose ends:
   *
   * Testing only done for WGS84 eccentricity.  Things that might go wrong are
   * other eccentricities: (1) Failure to converge most likely for much higher
   * eccentricity.  This will require refining the starting guesses for the
   * Newton's iterations.  (2) Failure with a sphere.  The basic formulation
   * carries over to the sphere with no problem.  Some attention might need to
   * be paid to the treatment of the singularity for phi=0, lam=90.
   *
   * The singularity at phi=90 is handled safely.  There's another singularity
   * (with the intermediate projection) at phi=0, lam=90*(1-e).  This is
   * handled by using a Taylor expansion about the singularity.  This gives a
   * good enough starting guess for Newton's method to converge.  However
   * detailed testing in the immediate neighborhood of the singularity has not
   * been done.  If there is a problem it can be handled easily by treating the
   * singularity specially.
   *
   * The initial guesses for Newton's method are a little ad hoc.  Probably
   * better guesses can be used and so one or more iterations of Newton's
   * method can be skipped.
   *
   * The reverse transformation of points outside the curved portion of the
   * equator, i.e., phi = 0, +/-lam in [90 * (1 - e), 90 * (1 + e)], correctly
   * transform to points in the opposite hemisphere in a narrow strip about lam
   * = +/- 90.  No attempt has been made to gauge the accuracy of the reverse
   * transformation of these points.
   * 
   */

  class TransverseMercatorExact {
  private:
    static const double tol, tol1, tol2, taytol, ahypover;
    static const int numit = 10;
    const double _a, _f, _k0, _mu, _mv, _e, _ep2;
    const bool _foldp;
    const EllipticFunction _Eu, _Ev;
    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    // These have poor relative accuracy near x = 0.  However, for mapping
    // applications, we only need good absolute accuracy.
    // For small arguments we would write
    //
    // asinh(x) = asin(x) -x^3/3-5*x^7/56-63*x^11/1408-143*x^15/5120 ...
    // atanh(x) = atan(x) +2*x^3/3+2*x^7/7+2*x^11/11+2*x^15/15
    //
    // The accuracy of asinh is also bad for large negative arguments.  This is
    // easy to fix in the definition of asinh.  Instead we call these functions
    // with positive arguments and enforce the correct parity separately.
    static inline double asinh(double x) { return log(x + sqrt(1 + sq(x))); }
    static inline double atanh(double x) { return log((1 + x)/(1 - x))/2; }
#endif
    double psi(double phi) const;
    double psiinv(double psi) const;

    void zeta(double u, double snu, double cnu, double dnu,
	      double v, double snv, double cnv, double dnv,
	      double& psi, double& lam) const;

    void dwdzeta(double u, double snu, double cnu, double dnu,
		 double v, double snv, double cnv, double dnv,
		 double& du, double& dv) const;

    bool zetainv0(double psi, double lam, double& u, double& v) const;
    void zetainv(double psi, double lam, double& u, double& v) const;

    void sigma(double u, double snu, double cnu, double dnu,
	       double v, double snv, double cnv, double dnv,
	       double& xi, double& eta) const;

    void dwdsigma(double u, double snu, double cnu, double dnu,
		  double v, double snv, double cnv, double dnv,
		  double& du, double& dv) const;

    bool sigmainv0(double xi, double eta, double& u, double& v) const;
    void sigmainv(double xi, double eta, double& u, double& v) const;

    void Scale(double phi, double lam,
	       double snu, double cnu, double dnu,
	       double snv, double cnv, double dnv,
	       double& gamma, double& k) const;

  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters), flattening \e f, and
     * central scale factor \e k0.  With \e foldp, the "standard" conventions
     * for cutting the transverse Mercator space are followed.  Forward can be
     * called with any \e lat and \e lon then produces the transformation shown
     * in Lee Fig 46.  Reverse analytically continues this in the +/- \e x
     * direction (with a cut at \e y = 0).
     *
     * With !\e foldp, the domains of \e lat, \e lon, \e x, and \e y are
     * restricted to
     * - the union of
     *   - \e lat in [0, 90] and \e lon - \e lon0 in [0, 90]
     *   - \e lat in (-90, 0] and \e lon - \e lon0 in [90 (1 - \e e), 90]
     * - the union of
     *   - \e x/(\e k0 \e a) in [0, inf) and
     *     \e y/(\e k0 \e a) in [0, E(\e e^2)]
     *   - \e x/(\e k0 \e a) in [K(1 - \e e^2) - E(1 - \e e^2), inf) and
     *     \e y/(\e k0 \e a) in (-inf, 0]
     * .
     * This allows the multi-valued nature of the projection to be explored
     * using its symmetries.
     **********************************************************************/
    TransverseMercatorExact(double a, double f, double k0, bool foldp = true);
    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * transverse Mercator easting \e x (meters) and northing \e y (meters).
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added.
     **********************************************************************/
    void Forward(double lon0, double lat, double lon,
		 double& x, double& y, double& gamma, double& k) const;
    /**
     * Convert from transverse Mercator easting \e x (meters) and northing \e y
     * (meters) to latitude \e lat (degrees) and longitude \e lon (degrees) .
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added.
     **********************************************************************/
    void Reverse(double lon0, double x, double y,
		 double& lat, double& lon, double& gamma, double& k) const;
    /**
     * A global instantiation of TransverseMercatorExact with the WGS84
     * ellipsoid and the UTM scale factor.
     **********************************************************************/
    const static TransverseMercatorExact UTM;
  };

} // namespace GeographicLib

#endif

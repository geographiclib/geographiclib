/**
 * \file TransverseMercatorExact.hpp
 * \brief Header for GeographicLib::TransverseMercatorExact class
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
   * the \e extendp argument to the constructor)..  The maximum error is about
   * 8 nm (ground distance) for the forward and reverse transformations.  The
   * error in the convergence is 2e-15", the relative error in the scale is
   * 7e-12%%.  (See \ref errors for the weasel words.)  The method is "exact"
   * in the sense that the errors are close to the round-off limit and that no
   * changes are needed in the algorithms for them to be used with reals of a
   * higher precision.  Thus the errors using long double (with a 64-bit
   * fraction) are about 2000 times smaller than using double (with a 53-bit
   * fraction).
   *
   * This algorithm is about 4.5 times slower than the 6th-order Kr&uuml;ger
   * method, GeographicLib::TransverseMercator, taking about 11 us for a
   * combined forward and reverse projection on a 2.6GHz Intel machine (g++,
   * version 4.3.0, -O3).
   *
   * See TransverseMercatorExact.cpp for more information on the
   * implementation.
   *
   * See \ref transversemercator for a discussion of this projection.
   **********************************************************************/

  class TransverseMercatorExact {
  private:
    static const double tol, tol1, tol2, taytol, ahypover;
    static const int numit = 10;
    const double _a, _f, _k0, _mu, _mv, _e, _ep2;
    const bool _extendp;
    const EllipticFunction _Eu, _Ev;
    static inline double sq(double x) throw() { return x * x; }
#if defined(_MSC_VER)
    static inline double hypot(double x, double y) throw()
    { return _hypot(x, y); }
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
    static inline double asinh(double x) throw() {
      return std::log(x + std::sqrt(1 + sq(x)));
    }
    static inline double atanh(double x) throw() {
      return std::log((1 + x)/(1 - x))/2;
    }
#else
    static inline double hypot(double x, double y) throw()
    { return ::hypot(x, y); }
    static inline double asinh(double x) throw() { return ::asinh(x); }
    static inline double atanh(double x) throw() { return ::atanh(x); }
#endif
    double psi(double phi) const throw();
    double psiinv(double psi) const throw();

    void zeta(double u, double snu, double cnu, double dnu,
	      double v, double snv, double cnv, double dnv,
	      double& psi, double& lam) const throw();

    void dwdzeta(double u, double snu, double cnu, double dnu,
		 double v, double snv, double cnv, double dnv,
		 double& du, double& dv) const throw();

    bool zetainv0(double psi, double lam, double& u, double& v) const throw();
    void zetainv(double psi, double lam, double& u, double& v) const throw();

    void sigma(double u, double snu, double cnu, double dnu,
	       double v, double snv, double cnv, double dnv,
	       double& xi, double& eta) const throw();

    void dwdsigma(double u, double snu, double cnu, double dnu,
		  double v, double snv, double cnv, double dnv,
		  double& du, double& dv) const throw();

    bool sigmainv0(double xi, double eta, double& u, double& v) const throw();
    void sigmainv(double xi, double eta, double& u, double& v) const throw();

    void Scale(double phi, double lam,
	       double snu, double cnu, double dnu,
	       double snv, double cnv, double dnv,
	       double& gamma, double& k) const throw();

  public:

    /**
     * Constructor for a ellipsoid radius \e a (meters), inverse flattening \e
     * invf, and central scale factor \e k0.  The transverse Mercator
     * projection has a branch point singularity at \e lat = 0 and \e lon - \e
     * lon0 = 90 (1 - \e e) or (for TransverseMercatorExact::UTM) x = 18381 km
     * , y = 0m.  The \e extendp argument governs where the branch cut is
     * placed.  With \e extendp = false, the "standard" convention is followed,
     * namely the cut is placed along x > 18381 km, y = 0m.  Forward can be
     * called with any \e lat and \e lon then produces the transformation shown
     * in Lee, Fig 46.  Reverse analytically continues this in the +/- \e x
     * direction.  As a consequence, Reverse may map multiple points to the
     * same geographic location; for example, for TransverseMercatorExact::UTM,
     * \e x = 22051449.037349 m, \e y = -7131237.022729 m and \e x =
     * 29735142.378357 m, \e y = 4235043.607933 m both map to \e lat = -2 deg,
     * \e lon = 88 deg.
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
     * method cannot be applied directly to the case of a sphere (\e invf =
     * inf) because some the constants characterizing this method diverge in
     * that limit.  However, GeographicLib::TransverseMercator treats the
     * sphere exactly.
     **********************************************************************/
    TransverseMercatorExact(double a, double invf, double k0,
			    bool extendp = false) throw();

    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * transverse Mercator easting \e x (meters) and northing \e y (meters).
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added.  \e lat should be in the range
     * [-90, 90]; \e lon and \e lon0 should be in the range [-180, 360].
     **********************************************************************/
    void Forward(double lon0, double lat, double lon,
		 double& x, double& y,
		 double& gamma, double& k) const throw();

    /**
     * Convert from transverse Mercator easting \e x (meters) and northing \e y
     * (meters) to latitude \e lat (degrees) and longitude \e lon (degrees) .
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added.  The value of \e lon returned is
     * in the range [-180, 180).
     **********************************************************************/
    void Reverse(double lon0, double x, double y,
		 double& lat, double& lon,
		 double& gamma, double& k) const throw();

    /**
     * A global instantiation of TransverseMercatorExact with the WGS84
     * ellipsoid and the UTM scale factor.  However, unlike UTM, no false
     * easting or northing is added.
     **********************************************************************/
    const static TransverseMercatorExact UTM;
  };

} // namespace GeographicLib

#endif

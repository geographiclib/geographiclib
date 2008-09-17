/**
 * \file TransverseMercatorExact.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/TransverseMercatorExact.hpp"
#include "GeographicLib/Constants.hpp"
#include <limits>

/*
 * Implementation of the Transverse Mercator Projection given in
 *
 *    L. P. Lee,
 *    Conformal Projections Based on Elliptic Functions,
 *    (B. V. Gutsell, Toronto, 1976), 128pp.,
 *    ISBN: 0919870163.
 *
 *    [Also appeared as:
 *    Monograph 16, Suppl. No. 1 to Canadian Cartographer, Vol 13.]
 *
 * The method entails using the Thompson Transverse Mercator as an intermediate
 * projection.  The projections from the intermediate coordinates to [phi, lam]
 * and [x, y] are given by elliptic functions.  The inverse of these
 * projections are found by Newton's method with a suitable starting guess.
 *
 * The same basic technique was used by
 *
 *   J. Dozier,
 *   Improved Algorithm for Calculation of UTM and Geodetic Coordinates
 *   NOAA Technical Report NESS 81, Sept. 1980
 *
 * The differences here are:
 *
 * * Use of real instead of complex arithmetic.  This uses the separation into
 *   real and imaginary parts given by Lee and allows better control over
 *   branch cuts.
 *
 * * Different (and possibly better?) algorithms for certain elliptic
 *   functions.
 *
 * * Better starting guesses for Newton's method.  Dozier's guesses lead to the
 *   wrong solution is some cases.  This disqualifies his method for routine
 *   use.
 *
 * This method gives the correct results for forward and reverse
 * transformations with the proviso that the reverse transformation may not
 * yield a sensible result if the input [x,y] is too far outside the set of
 * [x,y] obtained from the forward transformation.  The maximum error is about
 * 12nm (ground distance) for the forward and reverse transformations.  The
 * error in the convergence is 3e-9", the relative error in the scale is 1e-11.
 * The method is "exact" in the sense that the errors are close to the
 * round-off limit and that no changes are needed in the algorithms for them to
 * be used with reals of a higher precision (e.g., long double).
 *
 * This algorithm is about 2.5 times slower than the 6th-order series method
 * taking about 12us for a combined forward and reverse projection on a 2.6GHz
 * Intel machine (g++, version 4.3.0, -O3).
 *
 * This implementation and notation closely follows Lee, with the following
 * exceptions:
 *
 *     Lee    here    Description
 *     x/a    xi      Northing (unit Earth)
 *     y/a    eta     Easting (unit Earth)
 *     s/a    sigma   xi + i * eta
 *     y      x       Easting
 *     x      y       Northing
 *     k      e       eccentricity
 *     k^2    mu      elliptic function parameter
 *     k'^2   mv      elliptic function complementary parameter
 *     m      k       scale
 *
 * Minor alterations have been made in some of Lee's expressions in an attempt
 * to control round-off.  For example atanh(sin(phi)) is replaced by
 * asinh(tan(phi)) which maintains accuracy near phi = pi/2.  Such changes are
 * noted in the code.
 *
 * Loose ends:
 *
 * Testing only done for WGS84 eccentricity.  Things that might go wrong are
 * other eccentricities: (1) Failure to converge most likely for much higher
 * eccentricity.  This will require refining the starting guesses for the
 * Newton's iterations.  (2) Failure with a sphere.  The basic formulation
 * carries over to the sphere with no problem.  Some attention might need to be
 * paid to the treatment of the singularity for phi=0, lam=90.
 *
 * The singularity at phi=90 is handled safely.  There's another singularity
 * (with the intermediate projection) at phi=0, lam=90*(1-e).  This is handled
 * by using a Taylor expansion about the singularity.  This gives a good enough
 * starting guess for Newton's method to converge.  However detailed testing in
 * the immediate neighborhood of the singularity has not been done.  If there
 * is a problem it can be handled easily by treating the singularity specially.
 *
 * The initial guesses for Newton's method are a little ad hoc.  Probably
 * better guesses can be used and so one or more iterations of Newton's method
 * can be skipped.
 *
 * An analysis of the reverse transformation of "out-of-bounds" points has not
 * been done.  In particular it would be nice to determine the limits of the
 * correct reverse transformation.  These limits will of course depend
 * sensitively on the initial guess for Newton's method for the reverse
 * transformation.  (The out-of-bounds points are in the "opposite"
 * hemisphere.)
 * 
 */


namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = TRANSVERSEMERCATOREXACT_HPP;
}

namespace GeographicLib {

  TransverseMercatorExact::TransverseMercatorExact(double a, double invf,
						   double k0)
    : _a(a)
    , _f(1 / invf)
    , _k0(k0)
    , _mu(_f * (2 - _f))
    , _mv(1 - _mu)
    , _e(sqrt(_mu))
    , _tol(0.1*std::numeric_limits<double>::epsilon())
    , _tol1(0.1*sqrt(_tol))
    , Eu(_mu)
    , Ev(_mv)
    , _numit(10)
  {}

  const TransverseMercatorExact
  TransverseMercatorExact::UTM(Constants::WGS84_a, Constants::WGS84_invf,
			       Constants::UTM_k0);

  void TransverseMercatorExact::zeta(double u,
				     double snu, double cnu, double dnu,
				     double v,
				     double snv, double cnv, double dnv,
				     double& psi, double& lam) const {
    // Lee 54.17 but write
    // atanh(snu * dnv) = asinh(snu * dnv / sqrt(kv * snv^2 + cnu^2 * dnv^2)
    psi = asinh(snu * dnv / sqrt(_mv * snv * snv + cnu * cnu  * dnv * dnv))
      - _e * atanh(_e * snu / dnv);
    lam = atan2(dnu * snv, cnu * cnv)
      - _e * atan2(_e * cnu * snv, dnu * cnv);
  }

  void TransverseMercatorExact::dwdzeta(double u,
					double snu, double cnu, double dnu,
					double v,
					double snv, double cnv, double dnv,
					double& du, double& dv) const {
    // Lee 54.20
    double d = (1 - dnu * dnu * snv * snv);
    d = _mv * d * d;
    du = cnu * dnu * dnv * (cnv * cnv - _mu * snu * snu * snv * snv) / d;
    dv = -snu * snv * cnv * (dnu * dnu * dnv * dnv + _mu * cnu * cnu) / d;
  }

  void TransverseMercatorExact::sigma(double u,
				      double snu, double cnu, double dnu,
				      double v,
				      double snv, double cnv, double dnv,
				      double& xi, double& eta) const {
    // Lee 55.4 writing
    // dnu^2 + dnv^2 - 1 = _mu * cnu^2 + _mv * cnv^2
    double den = _mu * cnu * cnu + _mv * cnv * cnv;
    xi = Eu.E(snu, cnu, dnu) - _mu * snu * cnu * dnu / den;
    eta = v - Ev.E(snv, cnv, dnv) + _mv * snv * cnv * dnv / den;
  }

  void TransverseMercatorExact::dwdsigma(double u,
					 double snu, double cnu, double dnu,
					 double v,
					 double snv, double cnv, double dnv,
					 double& du, double& dv) const {
    // Reciprocal of 55.9: dw/ds = dn(w)^2/_mv, expanding complex dn(w) using
    // A+S 16.21.4
    // cnv^2 + _mu * snu^2 * snv^2 = 1 - dnu^2 * snv^2

    double den = (1 - dnu * dnu * snv * snv);
    den *= _mv * den;
    double
      dnr = dnu * cnv * dnv,
      dni = - _mu * snu * cnu * snv;
    du = (dnr * dnr - dni * dni) / den;
    dv = 2 * dnr * dni / den;
  }

  void TransverseMercatorExact::Scale(double phi,
				      double snu, double cnu, double dnu,
				      double snv, double cnv, double dnv,
				      double& gamma, double& k) const {
    // Lee 55.12 -- negated for our sign convention.  gamma gives the bearing
    // (clockwise from true north) of grid north
    gamma = atan2(_mv * snu * snv * cnv, cnu * dnu * dnv);
    k = sqrt(1 - _mu * sin(phi) * sin(phi))/cos(phi) *
      sqrt((1 - snu * snu * dnv * dnv) / (_mu * cnu * cnu + _mv * cnv * cnv));
  }

  void TransverseMercatorExact::Forward(double lon0, double lat, double lon,
					double& x, double& y,
					double& gamma, double& k) const {
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon - lon0 > 180)
      lon -= lon0 - 360;
    else if (lon - lon0 <= -180)
      lon -= lon0 + 360;
    else
      lon -= lon0;
    // Now lon in (-180, 180]
    // Explicitly enforce the parity
    int
      latsign = lat < 0 ? -1 : 1,
      lonsign = lon < 0 ? -1 : 1;
    lon *= lonsign;
    lat *= latsign;
    bool backside = lon > 90;
    if (backside) {
      if (lat == 0)
	latsign = -1;
      lon = 180 - lon;
    }
    double
      phi = lat * Constants::degree,
      lam = lon * Constants::degree;
    // u,v = coordinates for the Thompson TM, Lee 54
    double u, snu, cnu, dnu, v, snv, cnv, dnv;
    if (lat < 90) {
      double
	// Lee 9.4, replace atanh(sin(phi)) by more accurate asinh(tan(phi))
	psi = asinh(tan(phi)) - _e * atanh(_e * sin(phi));

      // Starting point for Newton's method to invert zeta(w).  Probably should
      // break this bit of magic into a separate function.
      if (psi < _e * Constants::pi && lam > (1 - 2 * _e) * Constants::pi/2) {
	// Expand about psi = 0, lam = (1 - _e) * pi/2
	// zeta - i*(lam - (1 - _e) *pi/2) = -(_mv * _e) / 3 * (w - i Ev.K())^3
	// mapping arg = [-pi/2, -pi/6] to [-pi/2, pi/2]
	double
	  rad = hypot(psi, lam - (1 - _e) * Constants::pi/2),
	  ang = atan2(psi, - (lam - (1 - _e) * Constants::pi/2));
	rad = std::pow(3 / (_mv * _e) * rad, 1/3.0);
	ang /= 3;
	u = rad * sin(ang);
	v = -rad * cos(ang) + Ev.K();
      } else {
	// Use spherical TM, Lee 12.6
	v = atanh(sin(lam) / cosh(psi));
	u = atan2(sinh(psi), cos(lam));
	// But scale u to put 90,0 on the right place
	u *= Eu.K() / (Constants::pi/2);
      }

      // Min iterations = 2, max iterations = 5; mean = 3.1
      for (int i = 0; i < _numit; ++i) {
	Eu.sncndn(u, snu, cnu, dnu);
	Ev.sncndn(v, snv, cnv, dnv);
	double psi1, lam1, du1, dv1;
	zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, psi1, lam1);
	dwdzeta(u, snu, cnu, dnu, v, snv, cnv, dnv, du1, dv1);
	psi1 -= psi;
	lam1 -= lam;
	double
	  delu = psi1 * du1 - lam1 * dv1,
	  delv = psi1 * dv1 + lam1 * du1;
	u -= delu;
	v -= delv;
	double delw2 = delu * delu + delv * delv;
	if (delw2 < _tol)
	  break;
      }

    } else {
      u = Eu.K();
      v = 0;
    }

    Eu.sncndn(u, snu, cnu, dnu);
    Ev.sncndn(v, snv, cnv, dnv);

    double xi, eta;
    sigma(u, snu, cnu, dnu, v, snv, cnv, dnv, xi, eta);
    if (backside)
      xi = 2 * Eu.E() - xi;
    y = xi * _a * _k0 * latsign;
    x = eta * _a * _k0 * lonsign;
    Scale(phi, snu, cnu, dnu, snv, cnv, dnv, gamma, k);
    gamma /= Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= latsign * lonsign;
    k *= _k0;
  }

  void TransverseMercatorExact::Reverse(double lon0, double x, double y,
					double& lat, double& lon,
					double& gamma, double& k) const {
    // This undoes the steps in Forward.
    double
      xi = y / (_a * _k0),
      eta = x / (_a * _k0);
    // Explicitly enforce the parity
    int
      latsign = y < 0 ? -1 : 1,
      lonsign = x < 0 ? -1 : 1;
    xi *= latsign;
    eta *= lonsign;
    bool backside = xi > Eu.E();
    if (backside)
      xi = 2 * Eu.E()- xi;

    // u,v = coordinates for the Thompson TM, Lee 54
    double u, snu, cnu, dnu, v, snv, cnv, dnv;

    // Starting point for Newton's method to invert sigma(w).  Probably should
    // break this bit of magic into a separate function.
    if ((eta > 0.75 * Ev.KE() && xi < 0.25 * Eu.E())
	|| (eta > Ev.KE())) {
      // Expand about xi = 0, eta = Ev.KE
      // s - i*Ev.KE() = -_mv / 3 * (w - i Ev.K())^3
      // mapping arg = [-pi/2, -pi/6] to [-pi/2, pi/2]
      double
	rad = hypot(xi, eta - Ev.KE()),
	ang = atan2(xi, -(eta - Ev.KE()));
      rad = std::pow(3 / _mv * rad, 1/3.0);
      ang /= 3;
      u = rad * sin(ang);
      v = -rad * cos(ang) + Ev.K();
    } else {
      u = xi * Eu.K()/Eu.E();
      v = eta * Eu.K()/Eu.E();
    }

    // Min iterations = 2, max iterations = 7; mean = 2.7
    for (int i = 0; i < _numit; ++i) {
      Eu.sncndn(u, snu, cnu, dnu);
      Ev.sncndn(v, snv, cnv, dnv);
      double xi1, eta1, du1, dv1;
      sigma(u, snu, cnu, dnu, v, snv, cnv, dnv, xi1, eta1);
      dwdsigma(u, snu, cnu, dnu, v, snv, cnv, dnv, du1, dv1);
      xi1 -= xi;
      eta1 -= eta;
      double
	delu = xi1 * du1 - eta1 * dv1,
	delv = xi1 * dv1 + eta1 * du1;
      u -= delu;
      v -= delv;
      double delw2 = delu * delu + delv * delv;
      if (delw2 < _tol)
	break;
    }

    Eu.sncndn(u, snu, cnu, dnu);
    Ev.sncndn(v, snv, cnv, dnv);
    double psi, lam, phi;
    if (v != 0 || u != Eu.K()) {
      zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, psi, lam);
      // Solve
      // psi = q - e * atanh(e * tanh(q))
      // for q = asinh(tan(phi))
      double q = psi;
      for (int i = 0; i < _numit; ++i) {
	// Max 3 trips through for double precision
	double
	  t = tanh(q),
	  dq = -(q - _e * atanh(_e * t) - psi) *
	  (1 - _mu * t * t) / _mv;
	q += dq;
	if (std::abs(dq) < _tol1)
	  break;
      }
      phi = atan(sinh(q));
    } else {
      phi = Constants::pi/2;
      lam = 0;
    }
    lat = phi / Constants::degree * latsign;
    lon = lam / Constants::degree;
    if (backside)
      lon = 180 - lon;
    lon *= lonsign;
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon + lon0 >= 180)
      lon += lon0 - 360;
    else if (lon + lon0 < -180)
      lon += lon0 + 360;
    else
      lon += lon0;
    Scale(phi, snu, cnu, dnu, snv, cnv, dnv, gamma, k);
    if (backside)
      y = 2 * Eu.E() - y;
    y *= _a * _k0 * latsign;
    x *= _a * _k0 * lonsign;
    gamma /= Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= latsign * lonsign;
    k *= _k0;
  }

} // namespace GeographicLib

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

#include <iostream>
#include <iomanip>

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
 * yield a sensible result if the input [x,y] is far outside the set of [x,y]
 * obtained from the forward transformation.  The maximum error is about 12nm
 * (ground distance) for the forward and reverse transformations.  The error in
 * the convergence is 3e-9", the relative error in the scale is 1e-11.  The
 * method is "exact" in the sense that the errors are close to the round-off
 * limit and that no changes are needed in the algorithms for them to be used
 * with reals of a higher precision (e.g., long double).
 *
 * This algorithm is about 2.5 times slower than the 6th-order Krueger method
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
 * The reverse transformation of points outside the curved portion of the
 * equator, i.e., phi = 0, +/-lam in [90 * (1 - e), 90 * (1 + e)], correctly
 * transform to points in the opposite hemisphere in a narrow strip about lam =
 * +/- 90.  No attempt has been made to gauge the accuracy of the reverse
 * transformation of these points.
 * 
 */


namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = TRANSVERSEMERCATOREXACT_HPP;
}

#if defined(_MSC_VER)
#define hypot _hypot
#endif

namespace GeographicLib {

  const double TransverseMercatorExact::tol =
    std::numeric_limits<double>::epsilon();
  const double TransverseMercatorExact::tol1 = 0.1*sqrt(tol);
  const double TransverseMercatorExact::taytol = std::pow(tol, 0.6);
      // Overflow value for asinh(tan(pi/2)) etc.
  const double TransverseMercatorExact::ahypover =
    double(std::numeric_limits<double>::digits)
    /log(double(std::numeric_limits<double>::radix)) + 2;

  TransverseMercatorExact::TransverseMercatorExact(double a, double invf,
						   double k0, bool fold)
    : _a(a)
    , _f(1 / invf)
    , _k0(k0)
    , _mu(_f * (2 - _f))
    , _mv(1 - _mu)
    , _e(sqrt(_mu))
    , _fold(fold)
    , Eu(_mu)
    , Ev(_mv)
  {}

  const TransverseMercatorExact
  TransverseMercatorExact::UTM(Constants::WGS84_a, Constants::WGS84_invf,
			       Constants::UTM_k0);

  double  TransverseMercatorExact::psi0(double phi) {
    // Rewrite asinh(tan(phi)) = atanh(sin(phi)) which is more accurate.
    // Write tan(phi) this way to ensure that sign(tan(phi)) = sign(phi)
    return asinh(sin(phi) / std::max(cos(phi), 0.1 * tol));
  }

  double  TransverseMercatorExact::psi(double phi) const {
    // Lee 9.4, replace atanh(sin(phi)) by more accurate asinh(tan(phi))
    // Write tan(phi) this way to ensure that sign(tan(phi)) = sign(phi)
    return psi0(phi) - _e * atanh(_e * sin(phi));
  }

  double TransverseMercatorExact::psiinv(double psi) const {
    // This is the inverse of psi.  Use Newton's method to solve for q in
    //
    //   psi = q - e * atanh(e * tanh(q))
    //
    // and then substitute phi = atan(sinh(q)).  Note that
    // dpsi/dq = (1 - e^2)/(1 - e^2 * tanh(q)^2)
    double q = psi;		// Initial guess
    for (int i = 0; i < numit; ++i) {
      // Max 3 trips through for double precision
      double
	t = tanh(q),
	dq = -(q - _e * atanh(_e * t) - psi) * (1 - _mu * sq(t)) / _mv;
      q += dq;
      if (std::abs(dq) < tol1)
	break;
    }
    return atan(sinh(q));
  }

  void TransverseMercatorExact::zeta(double u,
				     double snu, double cnu, double dnu,
				     double v,
				     double snv, double cnv, double dnv,
				     double& psi, double& lam) const {
    // Lee 54.17 but write
    // atanh(snu * dnv) = asinh(snu * dnv / sqrt(cnu^2 + _mv * snu^2 * snv^2))
    // atanh(_e * snu / dnv) = asinh(_e * snu / sqrt(_mu * cnu^2 + _mv * cnv^2))
    double
      d1 = sqrt(sq(cnu) + _mv * sq(snu * snv)),
      d2 = sqrt(_mu * sq(cnu) + _mv * sq(cnv));
    psi =
      // Overflow to values s.t. tanh = 1.
      (d1 ? asinh(snu * dnv / d1) : snu < 0 ? -ahypover : ahypover)
      - (d2 ? _e * asinh(_e * snu / d2) : snu < 0 ? -ahypover : ahypover);
    lam = (d1 != 0 && d2 != 0) ?
      atan2(dnu * snv, cnu * cnv) - _e * atan2(_e * cnu * snv, dnu * cnv) : 0;
  }

  void TransverseMercatorExact::dwdzeta(double u,
					double snu, double cnu, double dnu,
					double v,
					double snv, double cnv, double dnv,
					double& du, double& dv) const {
    // Lee 54.21 but write (1 - dnu^2 * snv^2) = (cnv^2 + _mu * snu^2 * snv^2)
    // (see A+S 16.21.4)
    double d = _mv * sq(sq(cnv) + _mu * sq(snu * snv));
    du =  cnu * dnu * dnv * (sq(cnv) - _mu * sq(snu * snv)) / d;
    dv = -snu * snv * cnv * (sq(dnu * dnv) + _mu * sq(cnu)) / d;
  }


  // Starting point for zetainv
  bool TransverseMercatorExact::zetainv0(double psi, double lam,
					  double& u, double& v) const {
    bool retval = false;
    if (psi < -_e * Constants::pi/4 &&
	lam > (1 - 2 * _e) * Constants::pi/2 &&
	psi < lam - (1 - _e) * Constants::pi/2) {
      // N.B. this branch is normally not taken because psi < 0 is converted
      // psi > 0 by Forward.
      //
      // There's a log singularity at w = w0 = Eu.K() + i * Ev.K(),
      // corresponding to the south pole, where we have, approximately
      //
      //   psi = _e + i * pi/2 - _e * atanh(cos(i * (w - w0)/(1 + _mu/2)))
      //
      // Inverting this gives:
      double
	psix = 1 - psi / _e,
	lamx = (Constants::pi/2 - lam) / _e;
      u = asinh(sin(lamx) / hypot(cos(lamx), sinh(psix))) * (1 + _mu/2);
      v = atan2(cos(lamx), sinh(psix)) * (1 + _mu/2);
      u = Eu.K() - u;
      v = Ev.K() - v;
    } else if (psi < _e * Constants::pi/2 &&
	       lam > (1 - 2 * _e) * Constants::pi/2) {
      // At w = w0 = i * Ev.K(), we have
      //
      //     zeta = zeta0 = i * (1 - _e) * pi/2
      //     zeta' = zeta'' = 0
      //
      // including the next term in the Taylor series gives:
      //
      // zeta = zeta0 - (_mv * _e) / 3 * (w - w0)^3
      //
      // When inverting this, we map arg(w - w0) = [-90, 0] to
      // arg(zeta - zeta0) = [-90, 180]
      double
	dlam = lam - (1 - _e) * Constants::pi/2,
	rad = hypot(psi, dlam),
	// atan2(dlam-psi, psi+dlam) + 45d gives arg(zeta - zeta0) in range
	// [-135, 225).  Subtracting 180 (since multiplier is negative) makes
	// range [-315, 45).  Multiplying by 1/3 (for cube root) gives range
	// [-105, 15).  In particular the range [-90, 180] in zeta space maps
	// to [-90, 0] in w space as required.
	ang = atan2(dlam-psi, psi+dlam) - 0.75 * Constants::pi;
      // Error using this guess is about 0.21 * (rad/e)^(5/3)
      retval = rad < _e * taytol;
      rad = std::pow(3 / (_mv * _e) * rad, 1/3.0);
      ang /= 3;
      u = rad * cos(ang);
      v = rad * sin(ang) + Ev.K();
    } else {
      // Use spherical TM, Lee 12.6 -- writing atanh(sin(lam) / cosh(psi)) =
      // asinh(sin(lam) / hypot(cos(lam), sinh(psi))).  This takes care of the
      // log singularity at zeta = Eu.K() (corresponding to the north pole)
      v = asinh(sin(lam) / hypot(cos(lam), sinh(psi)));
      u = atan2(sinh(psi), cos(lam));
      // But scale to put 90,0 on the right place
      u *= Eu.K() / (Constants::pi/2);
      v *= Eu.K() / (Constants::pi/2);
    }
    return retval;
  }

  // Invert zeta using Newton's method
  void  TransverseMercatorExact::zetainv(double psi, double lam,
					 double& u, double& v) const {
    if (zetainv0(psi, lam, u, v))
      return;

    // Min iterations = 2, max iterations = 4; mean = 3.1
    for (int i = 0; i < numit; ++i) {
      double snu, cnu, dnu, snv, cnv, dnv;
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
      double delw2 = sq(delu) + sq(delv);
      if (delw2 < tol)
	break;
    }
  }

  void TransverseMercatorExact::sigma(double u,
				      double snu, double cnu, double dnu,
				      double v,
				      double snv, double cnv, double dnv,
				      double& xi, double& eta) const {
    // Lee 55.4 writing
    // dnu^2 + dnv^2 - 1 = _mu * cnu^2 + _mv * cnv^2
    double d = _mu * sq(cnu) + _mv * sq(cnv);
    xi = Eu.E(snu, cnu, dnu) - _mu * snu * cnu * dnu / d;
    eta = v - Ev.E(snv, cnv, dnv) + _mv * snv * cnv * dnv / d;
  }

  void TransverseMercatorExact::dwdsigma(double u,
					 double snu, double cnu, double dnu,
					 double v,
					 double snv, double cnv, double dnv,
					 double& du, double& dv) const {
    // Reciprocal of 55.9: dw/ds = dn(w)^2/_mv, expanding complex dn(w) using
    // A+S 16.21.4
    double d = _mv * sq(sq(cnv) + _mu * sq(snu * snv));
    double
      dnr = dnu * cnv * dnv,
      dni = - _mu * snu * cnu * snv;
    du = (sq(dnr) - sq(dni)) / d;
    dv = 2 * dnr * dni / d;
  }

  // Starting point for sigmainv
  bool  TransverseMercatorExact::sigmainv0(double xi, double eta,
					   double& u, double& v) const {
    bool retval = false;
    if (eta > 1.25 * Ev.KE() ||
	(xi < -0.25 * Eu.E() && xi < eta - Ev.KE())) {
      // sigma as a simple pole at w = w0 = Eu.K() + i * Ev.K() and sigma is
      // approximated by
      // 
      // sigma = (Eu.E() + i * Ev.KE()) + 1/(w - w0)
      double
	x = xi - Eu.E(),
	y = eta - Ev.KE(),
	r2 = sq(x) + sq(y);
      u = Eu.K() + x/r2;
      v = Ev.K() - y/r2;      
    } else if ((eta > 0.75 * Ev.KE() && xi < 0.25 * Eu.E())
	       || eta > Ev.KE()) {
      // At w = w0 = i * Ev.K(), we have
      //
      //     sigma = sigma0 = i * Ev.KE()
      //     sigma' = sigma'' = 0
      //
      // including the next term in the Taylor series gives:
      //
      // sigma = sigma0 - _mv / 3 * (w - w0)^3
      //
      // When inverting this, we map arg(w - w0) = [-pi/2, -pi/6] to
      // arg(sigma - sigma0) = [-pi/2, pi/2]
      // mapping arg = [-pi/2, -pi/6] to [-pi/2, pi/2]
      double
	deta = eta - Ev.KE(),
	rad = hypot(xi, deta),
	// Map the range [-90, 180] in sigma space to [-90, 0] in w space.  See
	// discussion in zetainv0 on the cut for ang.
	ang = atan2(deta-xi, xi+deta) - 0.75 * Constants::pi;
      // Error using this guess is about 0.068 * rad^(5/3)
      retval = rad < 2 * taytol;
      rad = std::pow(3 / _mv * rad, 1/3.0);
      ang /= 3;
      u = rad * cos(ang);
      v = rad * sin(ang) + Ev.K();
      /*
      double
	rad = hypot(xi, eta - Ev.KE()),
	//ang = atan2(xi, -(eta - Ev.KE())),
	// Subtract pi because the multiplier of the cube term is negative
	ang = atan2(eta - Ev.KE(), xi) - Constants::pi;
      rad = std::pow(3 / _mv * rad, 1/3.0);
      ang /= 3;
      //ang1 /= 3;
      u = rad * cos(ang);
      v = rad * sin(ang) + Ev.K();
      //u = rad * sin(ang);
      //v = -rad * cos(ang) + Ev.K();
      // std::cerr << std::setprecision(18);
      // std::cerr << sin(ang) << " " << cos(ang1) << " "
      // << -cos(ang) << " " << sin(ang1) << std::endl;
      */
    } else {
      // Else use w = sigma * Eu.K/Eu.E (which is correct in the limit _e -> 0)
      u = xi * Eu.K()/Eu.E();
      v = eta * Eu.K()/Eu.E();
    }
    return retval;
  }

  // Invert sigma using Newton's method
  void  TransverseMercatorExact::sigmainv(double xi, double eta,
					  double& u, double& v) const {
    if (sigmainv0(xi, eta, u, v))
      return;

    // Min iterations = 2, max iterations = 6; mean = 2.7
    for (int i = 0; i < numit; ++i) {
      double snu, cnu, dnu, snv, cnv, dnv;
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
      double delw2 = sq(delu) + sq(delv);
      if (delw2 < tol)
	break;
    }
  }
  

  void TransverseMercatorExact::Scale(double phi, double lam,
				      double snu, double cnu, double dnu,
				      double snv, double cnv, double dnv,
				      double& gamma, double& k) const {
    double c = cos(phi);
    if (c > tol1) {
      // Lee 55.12 -- negated for our sign convention.  gamma gives the bearing
      // (clockwise from true north) of grid north
      gamma = atan2(_mv * snu * snv * cnv, cnu * dnu * dnv);
      // Lee 55.13 with nu given by Lee 9.1 -- in sqrt change the numerator
      // from
      //
      //    (1 - snu^2 * dnv^2) to (_mv * snv^2 + cnu^2 * dnv^2)
      //
      // to maintain accuracy near phi = 90 and change the denomintor from
      //
      //    (dnu^2 + dnv^2 - 1) to (_mu * cnu^2 + _mv * cnv^2)
      //
      // to maintain accuracy near phi = 0, lam = 90 * (1 - e).  Similarly
      // rewrite sqrt term in 9.1 as
      //
      //    _mv + _mu * c^2 instead of 1 - _mu * sin(phi)^2
      //
      // Finally replace 1/cos(phi) by cosh(psi0(phi)) which, near phi = pi/2,
      // behaves in a manner consistent with the last sqrt term.  (At least
      // that's the idea.)
      k = sqrt(_mv + _mu * sq(c)) * cosh(psi0(phi)) *
	sqrt( (_mv * sq(snv) + sq(cnu * dnv)) /
	      (_mu * sq(cnu) + _mv * sq(cnv)) );
    } else {
      // Near the pole
      gamma = lam;
      k = 1;
    }
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
      latsign = _fold && lat < 0 ? -1 : 1,
      lonsign = _fold && lon < 0 ? -1 : 1;
    lon *= lonsign;
    lat *= latsign;
    bool backside = _fold && lon > 90;
    if (backside) {
      if (lat == 0)
	latsign = -1;
      lon = 180 - lon;
    }
    double
      phi = lat * Constants::degree,
      lam = lon * Constants::degree;

    // u,v = coordinates for the Thompson TM, Lee 54
    double u, v;
    if (lat == 90) {
      u = Eu.K();
      v = 0;
    } else if (lat == 0 && lon == 90 * (1 - _e)) {
      u = 0;
      v = Ev.K();
    } else
      zetainv(psi(phi), lam, u, v);

    double snu, cnu, dnu, snv, cnv, dnv;
    Eu.sncndn(u, snu, cnu, dnu);
    Ev.sncndn(v, snv, cnv, dnv);

    double xi, eta;
    sigma(u, snu, cnu, dnu, v, snv, cnv, dnv, xi, eta);
    if (backside)
      xi = 2 * Eu.E() - xi;
    y = xi * _a * _k0 * latsign;
    x = eta * _a * _k0 * lonsign;
    Scale(phi, lam, snu, cnu, dnu, snv, cnv, dnv, gamma, k);
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
      latsign = _fold && y < 0 ? -1 : 1,
      lonsign = _fold && x < 0 ? -1 : 1;
    xi *= latsign;
    eta *= lonsign;
    bool backside = _fold && xi > Eu.E();
    if (backside)
      xi = 2 * Eu.E()- xi;

    // u,v = coordinates for the Thompson TM, Lee 54
    double u, v;
    if (xi == 0 && eta == Ev.KE()) {
      u = 0;
      v = Ev.K();
    } else
      sigmainv(xi, eta, u, v);

    double snu, cnu, dnu, snv, cnv, dnv;
    Eu.sncndn(u, snu, cnu, dnu);
    Ev.sncndn(v, snv, cnv, dnv);
    double phi, lam;
    if (v != 0 || u != Eu.K()) {
      double psi;
      zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, psi, lam);
      phi = psiinv(psi);
      lat = phi / Constants::degree;
      lon = lam / Constants::degree;
    } else {
      phi = Constants::pi/2;
      lat = 90;
      lon = lam = 0;
    }
    Scale(phi, lam, snu, cnu, dnu, snv, cnv, dnv, gamma, k);
    gamma /= Constants::degree;
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
    lat *= latsign;
    if (backside)
      y = 2 * Eu.E() - y;
    y *= _a * _k0 * latsign;
    x *= _a * _k0 * lonsign;
    if (backside)
      gamma = 180 - gamma;
    gamma *= latsign * lonsign;
    k *= _k0;
  }

} // namespace GeographicLib

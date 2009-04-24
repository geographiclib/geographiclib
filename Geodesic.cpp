/**
 * \file Geodesic.cpp
 * \brief Implementation for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 *
 * This is a reformulation of the geodesic problem.  The notation is as
 * follows:
 * - at a general point (no suffix or 1 or 2 as suffix)
 *   - phi = latitude
 *   - beta = latitude on auxilliary sphere
 *   - omega = longitude on auxilliary sphere
 *   - lambda = longitude
 *   - alpha = azimuth of great circle
 *   - sigma = arc length along greate circle
 *   - s = distance
 *   - tau = scaled distance (= sigma at multiples of pi/2)
 * - at northwards equator crossing
 *   - beta = phi = 0
 *   - omega = lambda = 0
 *   - alpha = alpha0
 *   - sigma = s = 0
 * - a 12 suffix means a difference, e.g., s12 = s2 - s1.
 * - s and c prefixes mean sin and cos
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/Constants.hpp"
#include <algorithm>
#include <limits>

#define GEODESIC_CPP "$Id$"

RCSID_DECL(GEODESIC_CPP)
RCSID_DECL(GEODESIC_HPP)

namespace GeographicLib {

  using namespace std;

  // Underflow guard.  We require
  //   eps2 * epsilon() > 0
  //   eps2 + epsilon() == epsilon()
  const double Geodesic::eps2 = sqrt(numeric_limits<double>::min());
  const double Geodesic::tol0 = numeric_limits<double>::epsilon();
  const double Geodesic::tol1 = 100 * tol0;
  const double Geodesic::tol2 = sqrt(numeric_limits<double>::epsilon());
  const double Geodesic::xthresh = 1000 * tol2;

  Geodesic::Geodesic(double a, double r) throw()
    : _a(a)
    , _f(r != 0 ? 1 / r : 0)
    , _f1(1 - _f)
    , _e2(_f * (2 - _f))
    , _ep2(_e2 / sq(_f1))	// e2 / (1 - e2)
    , _b(_a * _f1)
  {}

  const Geodesic Geodesic::WGS84(Constants::WGS84_a(), Constants::WGS84_r());

  double Geodesic::SinSeries(double sinx, double cosx,
			     const double c[], int n) throw() {
    // Evaluate y = sum(c[i - 1] * sin(2 * i * x), i, 1, n) using Clenshaw
    // summation.  (Indices into c offset by 1.)
    // Approx operation count = (n + 5) mult and (2 * n + 2) add
    double
      ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
      y0 = n & 1 ? c[--n] : 0, y1 = 0;	      // Accumulators for sum
    // Now n is even
    while (n) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + c[--n];
      y0 = ar * y1 - y0 + c[--n];
    }
    return 2 * sinx * cosx * y0; // sin(2 * x) * y0
  }

  GeodesicLine Geodesic::Line(double lat1, double lon1, double azi1)
    const throw() {
    return GeodesicLine(*this, lat1, lon1, azi1);
  }

  double Geodesic::Direct(double lat1, double lon1, double azi1, double s12,
			  double& lat2, double& lon2, double& azi2,
			  bool arcmode)
    const throw() {
    GeodesicLine l(*this, lat1, lon1, azi1);
    return l.Position(s12, lat2, lon2, azi2, arcmode);
  }

  double Geodesic::Inverse(double lat1, double lon1, double lat2, double lon2,
			   double& s12, double& azi1, double& azi2, double& m12)
    const throw() {
    lon1 = AngNormalize(lon1);
    double lon12 = AngNormalize(AngNormalize(lon2) - lon1);
    // If very close to being on the same meridian, then make it so.
    // Not sure this is necessary...
    lon12 = AngRound(lon12);
    // Make longitude difference positive.
    int lonsign = lon12 >= 0 ? 1 : -1;
    lon12 *= lonsign;
    if (lon12 == 180)
      lonsign = 1;
    // If really close to the equator, treat as on equator.
    lat1 = AngRound(lat1);
    lat2 = AngRound(lat2);
    // Swap points so that point with higher (abs) latitude is point 1
    int swapp = abs(lat1) >= abs(lat2) ? 1 : -1;
    if (swapp < 0) {
      lonsign *= -1;
      swap(lat1, lat2);
    }
    // Make lat1 <= 0
    int latsign = lat1 < 0 ? 1 : -1;
    lat1 *= latsign;
    lat2 *= latsign;
    // Now we have
    //
    //     0 <= lon12 <= 180
    //     -90 <= lat1 <= 0
    //     lat1 <= lat2 <= -lat1
    //
    // longsign, swapp, latsign register the transformation to bring the
    // coordinates to this canonical form.  In all cases, 1 means no change was
    // made.  We make these transformations so that there are few cases to
    // check, e.g., on verifying quadrants in atan2.  In addition, this
    // enforces some symmetries in the results returned.

    double phi, sbet1, cbet1, sbet2, cbet2, n1;

    phi = lat1 * Constants::degree();
    // Ensure cbet1 = +eps at poles
    sbet1 = _f1 * sin(phi);
    cbet1 = lat1 == -90 ? eps2 : cos(phi);
    // n = 1/sqrt(1 - e2 * sq(sin(phi)))
    n1 = 1/hypot(sbet1, cbet1);
    sbet1 *= n1; cbet1 *= n1;

    phi = lat2 * Constants::degree();
    // Ensure cbet2 = +eps at poles
    sbet2 = _f1 * sin(phi);
    cbet2 = abs(lat2) == 90 ? eps2 : cos(phi);
    SinCosNorm(sbet2, cbet2);

    double
      lam12 = lon12 * Constants::degree(),
      slam12 = lon12 == 180 ? 0 :sin(lam12),
      clam12 = cos(lam12);	// lon12 == 90 isn't interesting

    double sig12, calp1, salp1, calp2, salp2,
      c[ntau > neta ? (ntau ? ntau : 1) : (neta ? neta : 1)];
    // Enumerate all the cases where the geodesic is a meridian.  This includes
    // coincident points.
    bool meridian = lat1 == -90 || (_f >= 0 ? slam12 : lam12) == 0;
    if (!meridian && _f < 0 && lon12 == 180) {
      // For _f < 0 and lam12 = 180, need to check if we're beyond singular
      // point.  If lon12 == 180 then define bet2[ab] with
      //
      // tan(bet2[ab]) + tan(bet1) + H * (eta(bet2) + eta(bet1)) = +/- H * pi
      //
      // if bet2b < bet2 < bet2a, the geodesic is not a meridian
      double h0 = etaFactor(_f, 1.0);
      etaCoeff(_f, 1.0, c);
      double
	sbet12a = sbet2 * cbet1 + cbet2 * sbet1,
	cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
	bet12a = atan2(sbet12a, cbet12a), // bet12a = bet2 + bet1
	x = ( sbet12a / (cbet1 * cbet2)
	      + h0 * (bet12a + (SinSeries(sbet2, cbet2, c, neta) +
				SinSeries(sbet1, cbet1, c, neta))) ) /
	(h0 * Constants::pi());
      meridian = x <= -1;
    }
    if (meridian) {
      // Head to the target longitude
      calp1 = clam12; salp1 = slam12;
      // At the target we're heading north
      calp2 = 1; salp2 = 0;

      double
	// tan(bet) = tan(sig) * cos(alp),
	ssig1 = sbet1, csig1 = calp1 * cbet1,
	ssig2 = sbet2, csig2 = calp2 * cbet2;
      SinCosNorm(ssig1, csig1);
      SinCosNorm(ssig2, csig2);

      // sig12 = sig2 - sig1
      sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, 0.0),
		    csig1 * csig2 + ssig1 * ssig2);

      tauCoeff(_ep2, c);
      double
	dtau = dtauFactor(_ep2),
	et = (1 + dtau) * (SinSeries(ssig2, csig2, c, ntau) -
			   SinSeries(ssig1, csig1, c, ntau));
      zetCoeff(_ep2, c);
      double
	dzet = dzetFactor(_ep2),
	ez = (1 + dzet) * (SinSeries(ssig2, csig2, c, nzet) -
			   SinSeries(ssig1, csig1, c, nzet));

      m12 = _a * (sqrt(1 - _e2 * sq(cbet2)) * csig1 * ssig2 -
		  sqrt(1 - _e2 * sq(cbet1)) * ssig1 * csig2)
	- _b * csig1 * csig2 * ( (dtau - dzet) * sig12 + (et - ez) );

      s12 = _b * ((1 + dtau) * sig12 + et);
      sig12 /= Constants::degree();
    } else if (sbet1 == 0 &&	// and sbet2 == 0
	       (_f <= 0 ||
	       // Mimic the way Lambda12 works with calp1 = 0
		lam12 <= Constants::pi() - _f * Constants::pi())) {
      // Geodesic runs along equator
      calp1 = calp2 = 0; salp1 = salp2 = 1;
      s12 = _a * lam12;
      m12 = _b * sin(lam12 / _f1);
      sig12 = lon12 / _f1;
    } else {

      // Now point1 and point2 belong within a hemisphere bounded by a
      // meridian.

      // Figure a starting point for Newton's method
      InverseStart(sbet1, cbet1, n1, sbet2, cbet2,
		   lam12, slam12, clam12, salp1, calp1, c);

      // Newton's method
      double ssig1, csig1, ssig2, csig2, u2;
      double ov = 0;
      unsigned numit = 0;
      for (unsigned trip = 0; numit < maxit; ++numit) {
	double dv;
	double v = Lambda12(sbet1, cbet1, sbet2, cbet2, salp1, calp1,
			    salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
			    u2, trip < 1, dv, c) - lam12;
	if (abs(v) <= eps2 || !(trip < 1)) {
	  if (abs(v) > max(tol1, ov))
	    numit = maxit;
	  break;
	}
	double
	  dalp1 = -v/dv;
	double
	  sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
	  nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
	calp1 = calp1 * cdalp1 - salp1 * sdalp1;
	salp1 = max(0.0, nsalp1);
	SinCosNorm(salp1, calp1);
	// In some regimes we don't get quadratic convergence because slope ->
	// 0.  So use convergernce conditions based on epsilon instead of
	// sqrt(epsilon).  The first criterion is a test on abs(v) against 100
	// * epsilon.  The second takes credit for an anticipated reduction in
	// abs(v) by v/ov (due to the latest update in alp1) and checks this
	// against epsilon.
	if (abs(v) < tol1 || sq(v) < ov * tol0) ++trip;
	ov = abs(v);
      }

      tauCoeff(u2, c);
      double
	dtau = dtauFactor(u2),
	et = (1 + dtau) * (SinSeries(ssig2, csig2, c, ntau) -
			     SinSeries(ssig1, csig1, c, ntau));
      zetCoeff(u2, c);
      double
	dzet = dzetFactor(u2),
	ez = (1 + dzet) * (SinSeries(ssig2, csig2, c, nzet) -
			   SinSeries(ssig1, csig1, c, nzet));

      m12 = _a * (sqrt(1 - _e2 * sq(cbet2)) * csig1 * ssig2 -
		  sqrt(1 - _e2 * sq(cbet1)) * ssig1 * csig2)
	- _b * csig1 * csig2 * ( (dtau - dzet) * sig12 + (et - ez) );

      s12 = _b * ((1 + dtau) * sig12 + et);

      sig12 /= Constants::degree();
      if (numit >= maxit) {
	// Signal failure to converge by negating the distance and azimuths.
	s12 *= -1; sig12 *= -1;
	salp1 *= -1; calp1 *= -1;
	salp2 *= -1; calp2 *= -1;
      }
    }

    // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
    if (swapp < 0) {
      swap(salp1, salp2);
      swap(calp1, calp2);
    }

    // minus signs give range [-180, 180). 0- converts -0 to +0.
    azi1 = 0 - atan2(- swapp * lonsign * salp1,
		     + swapp * latsign * calp1) / Constants::degree();
    azi2 = 0 - atan2(- swapp * lonsign * salp2,
		     + swapp * latsign * calp2) / Constants::degree();
    // Returned value in [0, 180], unless it's negated to signal convergence
    // failure
    return sig12;
  }

  void Geodesic::InverseStart(double sbet1, double cbet1, double n1,
			      double sbet2, double cbet2,
			      double lam12, double slam12, double clam12,
			      double& salp1, double& calp1,
			      double c[]) const throw() {
    // Figure a starting point for Newton's method
    double
      // How close to antipodal lat?
      // bet12 = bet2 - bet1 in [0, pi);  bet12a = bet2 + bet1 in (-pi, 0]
      sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
      sbet12a = sbet2 * cbet1 + cbet2 * sbet1;

    salp1 = cbet2 * slam12;
    calp1 = clam12 >= 0 ?
      // The factor _f1*n1 applies a spheroidal correction for close points.
      // This saves 1 iteration of Newton's method in the case of short lines.
      sbet12 * _f1*n1 + cbet2 * sbet1 * sq(slam12) / (1 + clam12) :
      sbet12a - cbet2 * sbet1 * sq(slam12) / (1 - clam12);

    double
      ssig12 = hypot(salp1, calp1),
      csig12 = sbet1 * sbet2 + cbet1 * cbet2 * clam12;

    if (csig12 >= 0 || ssig12 >= 3 * abs(_f) * Constants::pi() * sq(cbet1)) {
      // Nothing to do, zeroth order spherical approximation is OK
    } else {
      // Scale lam12 and bet2 to x, y coordinate system where antipodal point
      // is at origin and singular point is at y = 0, x = -1.
      double x, y, lamscale, betscale;
      if (_f >= 0) {		// In fact f == 0 does not get here
	// x = dlong, y = dlat
	lamscale = -cbet1 * etaFactor(_f, sq(sbet1)) * Constants::pi();
	betscale = lamscale * cbet1;
	x = (lam12 - Constants::pi()) / lamscale;
	y = sbet12a / betscale;
      } else {			// _f < 0
	// x = dlat, y = dlong
	double
	  h0 = etaFactor(_f, 1.0),
	  cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
	  bet12a = atan2(sbet12a, cbet12a);
	etaCoeff(_f, 1.0, c);
	// In the case of lon12 = 180, this repeats a calculation made in
	// Inverse.
	x = ( sbet12a / (cbet1 * cbet2)
	      + h0 * (bet12a +  (SinSeries(sbet2, cbet2, c, neta) +
				 SinSeries(sbet1, cbet1, c, neta))) ) /
	  (h0 * Constants::pi());
	betscale = x < -0.01 ? sbet12a / x : -_f * sq(cbet1) * Constants::pi();
	lamscale = betscale / cbet1;
	y = (lam12 - Constants::pi()) / lamscale;
      }

      if (y > -100 * tol1 && x >  -1 - xthresh) {
	// strip near cut
	if (_f >= 0) {
	  salp1 = min(1.0, -x); calp1 = - sqrt(1 - sq(salp1));
	} else {
	  calp1 = max(-1.0, x); salp1 =   sqrt(1 - sq(calp1));
	}
      } else {
	// Estimate alp2, by solving calp2 * (salp2 + x) - y * salp2 = 0.  (For
	// f < 0, we're solving for pi/2 - alp2 and calp2 and salp2 are
	// swapped.)
	double salp2, calp2;
	if (y == 0) {
	  salp2 = 1; calp2 = 0; // This applies only for x < -1
	} else if (y > -0.027 && x > -1.09 && x < -0.91) {
	  // Near singular point we have
	  //     t^3 - 2*a*t - 2 = -t^2 * y^(2/3) approx 0
	  // where a = (x + 1)/|y|^(2/3), t = calp2/|y|^(1/3)
	  double
	    y3 = cbrt(-y),
	    a = (x + 1) / sq(y3),
	    a3 = sq(a) * a,
	    disc = 729 - 216 * a3, // sq(3 * sqrt(3) * sqrt(27 - 8 * a3))
	    b = 4 * a3 - 27,
	    t3 = a;		// t = - 3 / t3
	  if (disc >= 0) {
	    // b < 0 here, so use neg sqrt
	    double s = cbrt( (b - sqrt(disc)) / 4 );
	    t3 += s + sq(a)/s;
	  } else {
	    double ang = atan2(sqrt(-disc), b) + 2 * Constants::pi();
	    t3 += 2 * a * cos(ang/3);
	  }
	  calp2 = - 3 * y3 / t3;
	  // calp2 is small so sqrt is safe
	  salp2 = sqrt(1 - sq(calp2));
	} else {
	  salp2 = 0;
	  calp2 = 1;
	}
	// Now apply Newton's method to solve for salp2, calp2
	for (unsigned i = 0; i < 30; ++i) {
	  double
	    v = calp2 * (salp2 + x) - y * salp2,
	    dv = - calp2 * y - salp2 * x +
	    (calp2 - salp2) * (calp2 + salp2),
	    da = -v/dv,
	    sda = sin(da),
	    cda = cos(da),
	    nsalp2 = salp2 * cda + calp2 * sda;
	  if (v == 0)
	    break;
	  calp2 = max(0.0, calp2 * cda - salp2 * sda);
	  salp2 = max(0.0, nsalp2);
	  SinCosNorm(salp2, calp2);
	  if (abs(da) < tol2)
	    break;
	}
	// estimate omg12a = pi - omg12
	double
	  omg12a = lamscale * ( _f >= 0
				? hypot(y,  salp2 + x) * salp2
				: hypot(x, -calp2 + y) * calp2 ),
	  somg12 = sin(omg12a), comg12 = -cos(omg12a);
	// Update spherical estimate of alp1 using omg12 instead of lam12
	salp1 = cbet2 * somg12;
	calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
      }
    }
    SinCosNorm(salp1, calp1);
  }

  double Geodesic::Lambda12(double sbet1, double cbet1,
			    double sbet2, double cbet2,
			    double salp1, double calp1,
			    double& salp2, double& calp2,
			    double& sig12,
			    double& ssig1, double& csig1,
			    double& ssig2, double& csig2,
			    double& u2,
			    bool diffp, double& dlam12, double c[])
    const throw() {

    if (sbet1 == 0 && calp1 == 0)
      // Break degeneracy of equatorial line.  This cases has already been
      // handled.
      calp1 = -eps2;

    double
      // sin(alp1) * cos(bet1) = sin(alp0),
      salp0 = salp1 * cbet1,
      calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0

    double somg1, comg1, somg2, comg2, omg12, lam12, mu;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
    ssig1 = sbet1; somg1 = salp0 * sbet1;
    csig1 = comg1 = calp1 * cbet1;
    SinCosNorm(ssig1, csig1);
    SinCosNorm(somg1, comg1);

    // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
    // about this case, since this can yield singularities in the Newton
    // iteration.
    // sin(alp2) * cos(bet2) = sin(alp0),
    salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
    // calp2 = sqrt(1 - sq(salp2))
    //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
    // and subst for calp0 and rearrange to give (choose positive sqrt
    // to give alp2 in [0, pi/2]).
    calp2 = cbet2 != cbet1 || abs(sbet2) != -sbet1 ?
      sqrt(sq(calp1 * cbet1) + (cbet1 < -sbet1 ?
				(cbet2 - cbet1) * (cbet1 + cbet2) :
				(sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
      abs(calp1);
    // tan(bet2) = tan(sig2) * cos(alp2)
    // tan(omg2) = sin(alp0) * tan(sig2).
    ssig2 = sbet2; somg2 = salp0 * sbet2;
    csig2 = comg2 = calp2 * cbet2;
    SinCosNorm(ssig2, csig2);
    SinCosNorm(somg2, comg2);

    // sig12 = sig2 - sig1, limit to [0, pi]
    sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, 0.0),
		  csig1 * csig2 + ssig1 * ssig2);

    // omg12 = omg2 - omg1, limit to [0, pi]
    omg12 = atan2(max(comg1 * somg2 - somg1 * comg2, 0.0),
		  comg1 * comg2 + somg1 * somg2);
    double eta12, h0;
    mu = sq(calp0);
    etaCoeff(_f, mu, c);
    eta12 = SinSeries(ssig2, csig2, c, neta) - SinSeries(ssig1, csig1, c, neta);
    h0 = etaFactor(_f, mu),
    lam12 = omg12 + salp0 * h0 * (sig12 + eta12);
    u2 = mu * _ep2;

    if (diffp) {
      if (calp2 == 0)
	dlam12 = - 2 * sqrt(1 - _e2 * sq(cbet1)) / sbet1;
      else {
	tauCoeff(u2, c);
	double
	  dtau = dtauFactor(u2),
	  et = (1 + dtau) * (SinSeries(ssig2, csig2, c, ntau) -
			     SinSeries(ssig1, csig1, c, ntau));
	zetCoeff(u2, c);
	double
	  dzet = dzetFactor(u2),
	  ez = (1 + dzet) * (SinSeries(ssig2, csig2, c, nzet) -
			     SinSeries(ssig1, csig1, c, nzet));

	dlam12 = (sqrt(1 - _e2 * sq(cbet2)) * csig1 * ssig2 -
		   sqrt(1 - _e2 * sq(cbet1)) * ssig1 * csig2)
	  - _f1 * csig1 * csig2 * ( (dtau - dzet) * sig12 + (et - ez) );
	dlam12 /= calp2 * cbet2;
      }
    }

    return lam12;
  }

  GeodesicLine::GeodesicLine(const Geodesic& g,
			     double lat1, double lon1, double azi1) throw() {
    azi1 = Geodesic::AngNormalize(azi1);
    // Normalize azimuth at poles.  Evaluate azimuths at lat = +/- (90 - eps).
    if (lat1 == 90) {
      lon1 -= azi1 - (azi1 >= 0 ? 180 : -180);
      azi1 = -180;
    } else if (lat1 == -90) {
      lon1 += azi1;
      azi1 = 0;
    }
    // Guard against underflow in salp0
    azi1 = Geodesic::AngRound(azi1);
    lon1 = Geodesic::AngNormalize(lon1);
    _lat1 = lat1;
    _lon1 = lon1;
    _azi1 = azi1;
    _f1 = g._f1;
    // alp1 is in [0, pi]
    double
      alp1 = azi1 * Constants::degree(),
      // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
      // problems directly than to skirt them.
      salp1 = azi1 == 180 ? 0 : sin(alp1),
      calp1 = azi1 ==  90 ? 0 : cos(alp1);
    double cbet1, sbet1, phi;
    phi = lat1 * Constants::degree();
    // Ensure cbet1 = +eps at poles
    sbet1 = _f1 * sin(phi);
    cbet1 = abs(lat1) == 90 ? Geodesic::eps2 : cos(phi);
    Geodesic::SinCosNorm(sbet1, cbet1);

    // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    _salp0 = salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    _calp0 = Geodesic::hypot(calp1, salp1 * sbet1);
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles sce cbet1 = +eps.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    _ssig1 = sbet1; _somg1 = _salp0 * sbet1;
    _csig1 = _comg1 = sbet1 != 0 || calp1 != 0 ? cbet1 * calp1 : 1;
    Geodesic::SinCosNorm(_ssig1, _csig1); // sig1 in (-pi, pi]
    Geodesic::SinCosNorm(_somg1, _comg1);

    double mu = Geodesic::sq(_calp0);
    _u2 = mu * g._ep2;

    _sScale = g._b * (1 + Geodesic::dtauFactor(_u2));
    Geodesic::tauCoeff(_u2, _sigCoeff);
    _dtau1 = Geodesic::SinSeries(_ssig1, _csig1, _sigCoeff, ntau);
    {
      double s = sin(_dtau1), c = cos(_dtau1);
      // tau1 = sig1 + dtau1
      _stau1 = _ssig1 * c + _csig1 * s;
      _ctau1 = _csig1 * c - _ssig1 * s;
    }
    Geodesic::sigCoeff(_u2, _sigCoeff);
    // Not necessary because sigCoeff reverts tauCoeff
    //    _dtau1 = -SinSeries(_stau1, _ctau1, _sigCoeff, nsig);

    _etaFactor = _salp0 * Geodesic::etaFactor(g._f, mu);
    Geodesic::etaCoeff(g._f, mu, _etaCoeff);
    _dlam1 = Geodesic::SinSeries(_ssig1, _csig1, _etaCoeff, neta);
  }

  void GeodesicLine::ArcPosition(double sig12, double ssig12, double csig12,
				 double& lat2, double& lon2, double& azi2)
  const throw() {
    double omg12, lam12, lon12;
    double ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
    csig2 = _csig1 * csig12 - _ssig1 * ssig12;
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = _calp0 * ssig2;
    // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
    cbet2 = Geodesic::hypot(_salp0, _calp0 * csig2);
    // tan(omg2) = sin(alp0) * tan(sig2)
    somg2 = _salp0 * ssig2; comg2 = csig2;  // No need to normalize
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = _salp0; calp2 = _calp0 * csig2; // No need to normalize
    // omg12 = omg2 - omg1
    omg12 = atan2(somg2 * _comg1 - comg2 * _somg1,
		  comg2 * _comg1 + somg2 * _somg1);
    lam12 = omg12 + _etaFactor *
      ( sig12 +
	(Geodesic::SinSeries(ssig2, csig2, _etaCoeff, neta)  - _dlam1));
    lon12 = lam12 / Constants::degree();
    // Can't use AngNormalize because longitude might have wrapped multiple
    // times.
    lon12 = lon12 - 360 * floor(lon12/360 + 0.5);
    lat2 = atan2(sbet2, _f1 * cbet2) / Constants::degree();
    lon2 = Geodesic::AngNormalize(_lon1 + lon12);
    // minus signs give range [-180, 180). 0- converts -0 to +0.
    azi2 = 0 - atan2(-salp2, calp2) / Constants::degree();
  }

  double GeodesicLine::Position(double s12,
				double& lat2, double& lon2, double& azi2,
				bool arcmode)
  const throw() {
    if (_sScale == 0)
      // Uninitialized
      return 0;
    double sig12, ssig12, csig12;
    if (arcmode) {
      // Interpret s12 as spherical arc length
      sig12 = s12 * Constants::degree();
      double s12a = abs(s12);
      s12a -= 180 * floor(s12a / 180);
      ssig12 = s12a ==  0 ? 0 : sin(sig12);
      csig12 = s12a == 90 ? 0 : cos(sig12);
    } else {
      // Interpret s12 as distance
      double
	tau12 = s12 / _sScale,
	s = sin(tau12),
	c = cos(tau12);
      sig12 = tau12 + (_dtau1 +
		       // tau2 = tau1 + tau12
		       Geodesic::SinSeries(_stau1 * c + _ctau1 * s,
					   _ctau1 * c - _stau1 * s,
					   _sigCoeff, nsig));
      ssig12 = sin(sig12);
      csig12 = cos(sig12);
    }
    ArcPosition(sig12, ssig12, csig12, lat2, lon2, azi2);
    return arcmode ? s12 : sig12 /  Constants::degree();
  }

  // Generated by Maxima on 06:15:35 Fri, 3/13/2009 (GMT-5)
  //
  // Minor edits have been made: (1) remove trailing spaces, (2) promote
  // integers of 10 digits or more to doubles; (3) change 1/4 to 1/4.0; (4)
  // rudimentary line-breaking.

  // Generated by Maxima on 13:48:48 Thu, 4/23/2009 (GMT-5)

  // The scale factor, T-1, to convert tau to s / b
  double Geodesic::dtauFactor(double u2) throw() {
    double
      eps = u2 / (2 * (1 + sqrt(1 + u2)) + u2),
      eps2 = sq(eps),
      t;
    switch (tauord/2) {
    case 0:
      t = 0;
      break;
    case 1:
      t = eps2/4;
      break;
    case 2:
      t = eps2*(eps2+16)/64;
      break;
    case 3:
      t = eps2*(eps2*(eps2+4)+64)/256;
      break;
    case 4:
      t = eps2*(eps2*(eps2*(25*eps2+64)+256)+4096)/16384;
      break;
    default:
      STATIC_ASSERT(tauord >= 0 && tauord <= 8, "Bad value of tauord");
      t = 0;
    }
    return (t + eps) / (1 - eps);
  }

  // Coefficients, t[k], of sine series to convert sigma to tau
  void Geodesic::tauCoeff(double u2, double t[]) throw() {
    double
      eps = u2 / (2 * (1 + sqrt(1 + u2)) + u2),
      eps2 = sq(eps),
      d = eps;
    switch (ntau) {
    case 0:
      break;
    case 1:
      t[0] = -d/2;
      break;
    case 2:
      t[0] = -d/2;
      d *= eps;
      t[1] = -d/16;
      break;
    case 3:
      t[0] = d*(3*eps2-8)/16;
      d *= eps;
      t[1] = -d/16;
      d *= eps;
      t[2] = -d/48;
      break;
    case 4:
      t[0] = d*(3*eps2-8)/16;
      d *= eps;
      t[1] = d*(eps2-2)/32;
      d *= eps;
      t[2] = -d/48;
      d *= eps;
      t[3] = -5*d/512;
      break;
    case 5:
      t[0] = d*((6-eps2)*eps2-16)/32;
      d *= eps;
      t[1] = d*(eps2-2)/32;
      d *= eps;
      t[2] = d*(9*eps2-16)/768;
      d *= eps;
      t[3] = -5*d/512;
      d *= eps;
      t[4] = -7*d/1280;
      break;
    case 6:
      t[0] = d*((6-eps2)*eps2-16)/32;
      d *= eps;
      t[1] = d*((64-9*eps2)*eps2-128)/2048;
      d *= eps;
      t[2] = d*(9*eps2-16)/768;
      d *= eps;
      t[3] = d*(3*eps2-5)/512;
      d *= eps;
      t[4] = -7*d/1280;
      d *= eps;
      t[5] = -7*d/2048;
      break;
    case 7:
      t[0] = d*(eps2*(eps2*(19*eps2-64)+384)-1024)/2048;
      d *= eps;
      t[1] = d*((64-9*eps2)*eps2-128)/2048;
      d *= eps;
      t[2] = d*((72-9*eps2)*eps2-128)/6144;
      d *= eps;
      t[3] = d*(3*eps2-5)/512;
      d *= eps;
      t[4] = d*(35*eps2-56)/10240;
      d *= eps;
      t[5] = -7*d/2048;
      d *= eps;
      t[6] = -33*d/14336;
      break;
    case 8:
      t[0] = d*(eps2*(eps2*(19*eps2-64)+384)-1024)/2048;
      d *= eps;
      t[1] = d*(eps2*(eps2*(7*eps2-18)+128)-256)/4096;
      d *= eps;
      t[2] = d*((72-9*eps2)*eps2-128)/6144;
      d *= eps;
      t[3] = d*((96-11*eps2)*eps2-160)/16384;
      d *= eps;
      t[4] = d*(35*eps2-56)/10240;
      d *= eps;
      t[5] = d*(9*eps2-14)/4096;
      d *= eps;
      t[6] = -33*d/14336;
      d *= eps;
      t[7] = -429*d/262144;
      break;
    default:
      STATIC_ASSERT(ntau >= 0 && ntau <= 8, "Bad value of ntau");
    }
  }

  // Coefficients, t'[k], of sine series to convert tau to sigma
  void Geodesic::sigCoeff(double u2, double tp[]) throw() {
    double
      eps = u2 / (2 * (1 + sqrt(1 + u2)) + u2),
      eps2 = sq(eps),
      d = eps;
    switch (nsig) {
    case 0:
      break;
    case 1:
      tp[0] = d/2;
      break;
    case 2:
      tp[0] = d/2;
      d *= eps;
      tp[1] = 5*d/16;
      break;
    case 3:
      tp[0] = d*(16-9*eps2)/32;
      d *= eps;
      tp[1] = 5*d/16;
      d *= eps;
      tp[2] = 29*d/96;
      break;
    case 4:
      tp[0] = d*(16-9*eps2)/32;
      d *= eps;
      tp[1] = d*(30-37*eps2)/96;
      d *= eps;
      tp[2] = 29*d/96;
      d *= eps;
      tp[3] = 539*d/1536;
      break;
    case 5:
      tp[0] = d*(eps2*(205*eps2-432)+768)/1536;
      d *= eps;
      tp[1] = d*(30-37*eps2)/96;
      d *= eps;
      tp[2] = d*(116-225*eps2)/384;
      d *= eps;
      tp[3] = 539*d/1536;
      d *= eps;
      tp[4] = 3467*d/7680;
      break;
    case 6:
      tp[0] = d*(eps2*(205*eps2-432)+768)/1536;
      d *= eps;
      tp[1] = d*(eps2*(4005*eps2-4736)+3840)/12288;
      d *= eps;
      tp[2] = d*(116-225*eps2)/384;
      d *= eps;
      tp[3] = d*(2695-7173*eps2)/7680;
      d *= eps;
      tp[4] = 3467*d/7680;
      d *= eps;
      tp[5] = 38081*d/61440;
      break;
    case 7:
      tp[0] = d*(eps2*((9840-4879*eps2)*eps2-20736)+36864)/73728;
      d *= eps;
      tp[1] = d*(eps2*(4005*eps2-4736)+3840)/12288;
      d *= eps;
      tp[2] = d*(eps2*(8703*eps2-7200)+3712)/12288;
      d *= eps;
      tp[3] = d*(2695-7173*eps2)/7680;
      d *= eps;
      tp[4] = d*(41604-141115*eps2)/92160;
      d *= eps;
      tp[5] = 38081*d/61440;
      d *= eps;
      tp[6] = 459485*d/516096;
      break;
    case 8:
      tp[0] = d*(eps2*((9840-4879*eps2)*eps2-20736)+36864)/73728;
      d *= eps;
      tp[1] = d*(eps2*((120150-86171*eps2)*eps2-142080)+115200)/368640;
      d *= eps;
      tp[2] = d*(eps2*(8703*eps2-7200)+3712)/12288;
      d *= eps;
      tp[3] = d*(eps2*(1082857*eps2-688608)+258720)/737280;
      d *= eps;
      tp[4] = d*(41604-141115*eps2)/92160;
      d *= eps;
      tp[5] = d*(533134-2200311*eps2)/860160;
      d *= eps;
      tp[6] = 459485*d/516096;
      d *= eps;
      tp[7] = 109167851*d/82575360;
      break;
    default:
      STATIC_ASSERT(nsig >= 0 && nsig <= 8, "Bad value of nsig");
    }
  }

  // The scale factor, Z-1
  double Geodesic::dzetFactor(double u2) throw() {
    double
      eps = u2 / (2 * (1 + sqrt(1 + u2)) + u2),
      eps2 = sq(eps),
      t;
    switch (zetord/2) {
    case 0:
      t = 0;
      break;
    case 1:
      t = eps2/4;
      break;
    case 2:
      t = eps2*(9*eps2+16)/64;
      break;
    case 3:
      t = eps2*(eps2*(25*eps2+36)+64)/256;
      break;
    case 4:
      t = eps2*(eps2*(eps2*(1225*eps2+1600)+2304)+4096)/16384;
      break;
    default:
      STATIC_ASSERT(zetord >= 0 && zetord <= 8, "Bad value of zetord");
      t = 0;
    }
    return t * (1 - eps) - eps;
  }

  // Coefficients, z[k], of sine series to convert sigma to zeta
  void Geodesic::zetCoeff(double u2, double z[]) throw() {
    double
      eps = u2 / (2 * (1 + sqrt(1 + u2)) + u2),
      eps2 = sq(eps),
      d = eps;
    switch (nzet) {
    case 0:
      break;
    case 1:
      z[0] = d/2;
      break;
    case 2:
      z[0] = d/2;
      d *= eps;
      z[1] = 3*d/16;
      break;
    case 3:
      z[0] = d*(eps2+8)/16;
      d *= eps;
      z[1] = 3*d/16;
      d *= eps;
      z[2] = 5*d/48;
      break;
    case 4:
      z[0] = d*(eps2+8)/16;
      d *= eps;
      z[1] = d*(eps2+6)/32;
      d *= eps;
      z[2] = 5*d/48;
      d *= eps;
      z[3] = 35*d/512;
      break;
    case 5:
      z[0] = d*(eps2*(eps2+2)+16)/32;
      d *= eps;
      z[1] = d*(eps2+6)/32;
      d *= eps;
      z[2] = d*(15*eps2+80)/768;
      d *= eps;
      z[3] = 35*d/512;
      d *= eps;
      z[4] = 63*d/1280;
      break;
    case 6:
      z[0] = d*(eps2*(eps2+2)+16)/32;
      d *= eps;
      z[1] = d*(eps2*(35*eps2+64)+384)/2048;
      d *= eps;
      z[2] = d*(15*eps2+80)/768;
      d *= eps;
      z[3] = d*(7*eps2+35)/512;
      d *= eps;
      z[4] = 63*d/1280;
      d *= eps;
      z[5] = 77*d/2048;
      break;
    case 7:
      z[0] = d*(eps2*(eps2*(41*eps2+64)+128)+1024)/2048;
      d *= eps;
      z[1] = d*(eps2*(35*eps2+64)+384)/2048;
      d *= eps;
      z[2] = d*(eps2*(69*eps2+120)+640)/6144;
      d *= eps;
      z[3] = d*(7*eps2+35)/512;
      d *= eps;
      z[4] = d*(105*eps2+504)/10240;
      d *= eps;
      z[5] = 77*d/2048;
      d *= eps;
      z[6] = 429*d/14336;
      break;
    case 8:
      z[0] = d*(eps2*(eps2*(41*eps2+64)+128)+1024)/2048;
      d *= eps;
      z[1] = d*(eps2*(eps2*(47*eps2+70)+128)+768)/4096;
      d *= eps;
      z[2] = d*(eps2*(69*eps2+120)+640)/6144;
      d *= eps;
      z[3] = d*(eps2*(133*eps2+224)+1120)/16384;
      d *= eps;
      z[4] = d*(105*eps2+504)/10240;
      d *= eps;
      z[5] = d*(33*eps2+154)/4096;
      d *= eps;
      z[6] = 429*d/14336;
      d *= eps;
      z[7] = 6435*d/262144;
      break;
    default:
      STATIC_ASSERT(nzet >= 0 && nzet <= 8, "Bad value of nzet");
    }
  }

  // The scale factor, H, to convert eta to changes in lambda
  double Geodesic::etaFactor(double f, double mu) throw() {
  double
    e2 = f*(2-f),
    ep2 = e2/sq(1-f),
    u2 = ep2*mu,
    eps = u2 / (2 * (1 + sqrt(1 + u2)) + u2),
    fp = (f - eps) / (1 - eps),
    nu = fp !=0 ? 2 * eps / fp : mu / (1 - mu/2),
    nu2 = sq(nu);
  double g;
  switch (etaord) {
    case 0:
      g = 0;
      break;
    case 1:
      g = 1;
      break;
    case 2:
      g = 1;
      break;
    case 3:
      g = 1;
      break;
    case 4:
      g = (64-fp*sq(fp)*nu2)/64;
      break;
    case 5:
      g = (fp*sq(fp)*(-fp*sq(nu2)-16*nu2)+1024)/1024;
      break;
    case 6:
      g = (fp*sq(fp)*(fp*(4*fp*nu2-sq(nu2))-16*nu2)+1024)/1024;
      break;
    case 7:
      g = (fp*sq(fp)*(fp*(fp*(fp*nu2*((8-nu2)*nu2+32)+32*nu2)-8*sq(nu2))-128*nu2)+8192)/8192;
      break;
    case 8:
      g = (fp*sq(fp)*(fp*(fp*(fp*(fp*nu2*(nu2*(15*nu2+56)+384)+nu2*((128-16*nu2)*nu2+512))+512*nu2)-128*sq(nu2))-2048*nu2)+131072)/131072;
      break;
    default:
      STATIC_ASSERT(etaord >= 0 && etaord <= 8, "Bad value of etaord");
      g = 0;
    }
  return -f * (2 - f)/(2 - fp) * g;
  }

  // Coefficients, h[k], of sine series to convert sigma to eta
  void Geodesic::etaCoeff(double f, double mu, double h[]) throw() {
    double
      e2 = f*(2-f),
      ep2 = e2/sq(1-f),
      u2 = ep2*mu,
      eps = u2 / (2 * (1 + sqrt(1 + u2)) + u2),
      fp = (f - eps) / (1 - eps),
      nu = fp !=0 ? 2 * eps / fp : mu / (1 - mu/2),
      nu2 = sq(nu);
    double s = fp * nu, d = s;
    switch (neta) {
    case 0:
      break;
    case 1:
      h[0] = d/8;
      break;
    case 2:
      h[0] = d*(2-fp)/16;
      d *= s;
      h[1] = d/64;
      break;
    case 3:
      h[0] = d*(fp*(fp*(-nu2-16)-32)+64)/512;
      d *= s;
      h[1] = d*(4-3*fp)/256;
      d *= s;
      h[2] = 5*d/1536;
      break;
    case 4:
      h[0] = d*(fp*(fp*(-2*nu2+fp*(-nu2-16)-32)-64)+128)/1024;
      d *= s;
      h[1] = d*(fp*(fp*(-nu2-8)-24)+32)/2048;
      d *= s;
      h[2] = d*(10-9*fp)/3072;
      d *= s;
      h[3] = 7*d/8192;
      break;
    case 5:
      h[0] = d*(fp*(fp*(fp*(fp*((4-nu2)*nu2-32)-4*nu2-64)-8*nu2-128)-256)+512)/4096;
      d *= s;
      h[1] = d*(fp*(fp*(-nu2-2*fp-8)-24)+32)/2048;
      d *= s;
      h[2] = d*(fp*(fp*(-7*nu2-32)-144)+160)/49152;
      d *= s;
      h[3] = d*(7-7*fp)/8192;
      d *= s;
      h[4] = 21*d/81920;
      break;
    case 6:
      h[0] = d*(fp*(fp*(fp*(fp*(fp*(nu2*(3*nu2+32)-128)+(32-8*nu2)*nu2-256)-32*nu2-512)-64*nu2-1024)-2048)+4096)/32768;
      d *= s;
      h[1] = d*(fp*(fp*(fp*(fp*(40-7*nu2)*nu2-128)-64*nu2-512)-1536)+2048)/131072;
      d *= s;
      h[2] = d*(fp*(fp*(2*fp*nu2-7*nu2-32)-144)+160)/49152;
      d *= s;
      h[3] = d*(fp*(fp*(-3*nu2-8)-56)+56)/65536;
      d *= s;
      h[4] = d*(42-45*fp)/163840;
      d *= s;
      h[5] = 11*d/131072;
      break;
    case 7:
      h[0] = d*(fp*(fp*(fp*(fp*(fp*(fp*(nu2*((320-61*nu2)*nu2+1280)-4096)+nu2*(192*nu2+2048)-8192)+(2048-512*nu2)*nu2-16384)-2048*nu2-32768)-4096*nu2-65536)-131072)+262144)/2097152;
      d *= s;
      h[1] = d*(fp*(fp*(fp*(fp*(fp*(nu2*(37*nu2+192)+256)+(320-56*nu2)*nu2)-1024)-512*nu2-4096)-12288)+16384)/1048576;
      d *= s;
      h[2] = d*(fp*(fp*(fp*(fp*((560-91*nu2)*nu2+768)+256*nu2)-896*nu2-4096)-18432)+20480)/6291456;
      d *= s;
      h[3] = d*(fp*(fp*(fp*(23*nu2+40)-48*nu2-128)-896)+896)/1048576;
      d *= s;
      h[4] = d*(fp*(fp*(-165*nu2-240)-2880)+2688)/10485760;
      d *= s;
      h[5] = d*(88-99*fp)/1048576;
      d *= s;
      h[6] = 429*d/14680064;
      break;
    default:
      STATIC_ASSERT(neta >= 0 && neta <= 7, "Bad value of neta");
    }
  }

} // namespace GeographicLib


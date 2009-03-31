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
 *   - chi = longitude on auxilliary sphere
 *   - lambda = longitude
 *   - alpha = azimuth of great circle
 *   - sigma = arc length along greate circle
 *   - s = distance
 *   - tau = scaled distance (= sigma at multiples of pi/2)
 * - at northwards equator crossing
 *   - beta = phi = 0
 *   - chi = lambda = 0
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
			   double& s12, double& azi1, double& azi2)
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
    // n = sqrt(1 - e2 * sq(sin(phi)))
    n1 = hypot(sbet1, cbet1);
    sbet1 /= n1; cbet1 /= n1;

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
      s12 = _b * tauFactor(_ep2) *
	(sig12 + (SinSeries(ssig2, csig2, c, ntau) -
		  SinSeries(ssig1, csig1, c, ntau)));
      sig12 /= Constants::degree();
    } else if (sbet1 == 0 &&	// and sbet2 == 0
	       (_f <= 0 ||
	       // Mimic the way Lambda12 works with calp1 = 0
		lam12 <= Constants::pi() - _f * Constants::pi())) {
      // Geodesic runs along equator
      calp1 = calp2 = 0; salp1 = salp2 = 1;
      s12 = _a * lam12;
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
      s12 = _b * tauFactor(u2) *
	(sig12 + (SinSeries(ssig2, csig2, c, ntau) -
		  SinSeries(ssig1, csig1, c, ntau)));
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
    // Returned value in [0, 180]
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
      // The factor _f1/n1 applies a spheroidal correction for close points.
      // This saves 1 iteration of Newton's method in the case of short lines.
      sbet12 * _f1/n1 + cbet2 * sbet1 * sq(slam12) / (1 + clam12) :
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
	// estimate chi12a = pi - chi12
	double
	  chi12a = lamscale * ( _f >= 0
				? hypot(y,  salp2 + x) * salp2
				: hypot(x, -calp2 + y) * calp2 ),
	  schi12 = sin(chi12a), cchi12 = -cos(chi12a);
	// Update spherical estimate of alp1 using chi12 instead of lam12
	salp1 = cbet2 * schi12;
	calp1 = sbet12a - cbet2 * sbet1 * sq(schi12) / (1 - cchi12);
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

    double schi1, cchi1, schi2, cchi2, chi12, lam12, mu;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(chi1) = sin(alp0) * tan(sig1) = tan(chi1)=tan(alp1)*sin(bet1)
    ssig1 = sbet1; schi1 = salp0 * sbet1;
    csig1 = cchi1 = calp1 * cbet1;
    SinCosNorm(ssig1, csig1);
    SinCosNorm(schi1, cchi1);

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
    // tan(chi2) = sin(alp0) * tan(sig2).
    ssig2 = sbet2; schi2 = salp0 * sbet2;
    csig2 = cchi2 = calp2 * cbet2;
    SinCosNorm(ssig2, csig2);
    SinCosNorm(schi2, cchi2);

    // sig12 = sig2 - sig1, limit to [0, pi]
    sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, 0.0),
		  csig1 * csig2 + ssig1 * ssig2);

    // chi12 = chi2 - chi1, limit to [0, pi]
    chi12 = atan2(max(cchi1 * schi2 - schi1 * cchi2, 0.0),
		  cchi1 * cchi2 + schi1 * schi2);
    double eta12, h0;
    mu = sq(calp0);
    etaCoeff(_f, mu, c);
    eta12 = SinSeries(ssig2, csig2, c, neta) - SinSeries(ssig1, csig1, c, neta);
    h0 = etaFactor(_f, mu),
    lam12 = chi12 + salp0 * h0 * (sig12 + eta12);

    if (diffp) {
      double dalp0, dsig1, dchi1, dalp2, dsig2, dchi2;
      // Differentiate sin(alp) * cos(bet) = sin(alp0),
      dalp0 = cbet1 * calp1 / calp0;
      dalp2 = calp2 != 0 ? calp1 * cbet1/ (calp2 * cbet2) :
	calp1 >= 0 ? 1 : -1;
      // Differentiate tan(bet) = tan(sig) * cos(alp) and clear
      // calp from the denominator with tan(alp0)=cos(sig)*tan(alp),
      dsig1 = ssig1 * salp0 / calp0;
      dsig2 = ssig2 * salp0 / calp0 * dalp2;
      // Differentiate tan(chi) = sin(alp0) * tan(sig).  Substitute
      //   tan(sig) = tan(bet) / cos(alp) = tan(chi) / sin(alp0)
      //   cos(chi) / cos(sig) = 1 / cos(bet)
      // to give
      dchi1 = (sbet1 * sq(cchi1) + schi1 * salp0 / (calp0 * cbet1));
      dchi2 = (sbet2 * sq(cchi2) + schi2 * salp0 / (calp0 * cbet2)) * dalp2;

      double deta12, dmu, dh0, dlamsig;
      etaCoeffmu(_f, mu, c);
      dmu = - 2 * calp0 * salp0 * dalp0;
      deta12 = dmu * (SinSeries(ssig2, csig2, c, neta) -
		      SinSeries(ssig1, csig1, c, neta));
      dh0 = etaFactormu(_f, mu) * dmu;

      // Derivative of salp0 * h0 * (sig + eta) wrt sig.  This
      // is from integral form of this expression.
      dlamsig = - _e2 * salp0 *
	(dsig2 / (sqrt(1 - _e2 * (1 - mu * sq(ssig2))) + 1) -
	 dsig1 / (sqrt(1 - _e2 * (1 - mu * sq(ssig1))) + 1)) ;

      dlam12 =
	(dchi2 - dchi1) + dlamsig +
	// Derivative wrt mu
	(dalp0 * calp0 * h0 + salp0 * dh0) * (sig12 + eta12) +
	salp0 * h0 * deta12;
    }

    u2 = mu * _ep2;
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
    // Evaluate chi1 with tan(chi1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and chi coincide.
    // No atan2(0,0) ambiguity at poles sce cbet1 = +eps.
    // With alp0 = 0, chi1 = 0 for alp1 = 0, chi1 = pi for alp1 = pi.
    _ssig1 = sbet1; _schi1 = _salp0 * sbet1;
    _csig1 = _cchi1 = sbet1 != 0 || calp1 != 0 ? cbet1 * calp1 : 1;
    Geodesic::SinCosNorm(_ssig1, _csig1); // sig1 in (-pi, pi]
    Geodesic::SinCosNorm(_schi1, _cchi1);

    double
      mu = Geodesic::sq(_calp0),
      u2 = mu * g._ep2;

    _sScale = g._b * Geodesic::tauFactor(u2);
    Geodesic::tauCoeff(u2, _sigCoeff);
    _dtau1 = Geodesic::SinSeries(_ssig1, _csig1, _sigCoeff, ntau);
    {
      double s = sin(_dtau1), c = cos(_dtau1);
      // tau1 = sig1 + dtau1
      _stau1 = _ssig1 * c + _csig1 * s;
      _ctau1 = _csig1 * c - _ssig1 * s;
    }
    Geodesic::sigCoeff(u2, _sigCoeff);
    // Not necessary because sigCoeff reverts tauCoeff
    //    _dtau1 = -SinSeries(_stau1, _ctau1, _sigCoeff, nsig);

    _etaFactor = _salp0 * Geodesic::etaFactor(g._f, mu);
    Geodesic::etaCoeff(g._f, mu, _etaCoeff);
    _dlam1 = Geodesic::SinSeries(_ssig1, _csig1, _etaCoeff, neta);
  }

  void GeodesicLine::ArcPosition(double sig12, double ssig12, double csig12,
				 double& lat2, double& lon2, double& azi2)
  const throw() {
    double chi12, lam12, lon12;
    double ssig2, csig2, sbet2, cbet2, schi2, cchi2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
    csig2 = _csig1 * csig12 - _ssig1 * ssig12;
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = _calp0 * ssig2;
    // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
    cbet2 = Geodesic::hypot(_salp0, _calp0 * csig2);
    // tan(chi2) = sin(alp0) * tan(sig2)
    schi2 = _salp0 * ssig2; cchi2 = csig2;  // No need to normalize
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = _salp0; calp2 = _calp0 * csig2; // No need to normalize
    // chi12 = chi2 - chi1
    chi12 = atan2(schi2 * _cchi1 - cchi2 * _schi1,
		  cchi2 * _cchi1 + schi2 * _schi1);
    lam12 = chi12 + _etaFactor *
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

  // The scale factor, T, to convert tau to s / b
  double Geodesic::tauFactor(double u2) throw() {
    switch (tauord) {
    case 0:
      return 1;
      break;
    case 1:
      return (u2+4)/4;
      break;
    case 2:
      return ((16-3*u2)*u2+64)/64;
      break;
    case 3:
      return (u2*(u2*(5*u2-12)+64)+256)/256;
      break;
    case 4:
      return (u2*(u2*((320-175*u2)*u2-768)+4096)+16384)/16384;
      break;
    case 5:
      return (u2*(u2*(u2*(u2*(441*u2-700)+1280)-3072)+16384)+65536)/65536;
      break;
    case 6:
      return (u2*(u2*(u2*(u2*((7056-4851*u2)*u2-11200)+20480)-49152)+262144)+
	1048576)/1048576;
      break;
    case 7:
      return (u2*(u2*(u2*(u2*(u2*(u2*(14157*u2-19404)+28224)-44800)+81920)-
	196608)+1048576)+4194304)/4194304;
      break;
    case 8:
      return (u2*(u2*(u2*(u2*(u2*(u2*((3624192-2760615*u2)*u2-4967424)+
	7225344)-11468800)+20971520)-50331648)+268435456)+1073741824.0)/
	1073741824.0;
      break;
    default:
      STATIC_ASSERT(tauord >= 0 && tauord <= 8, "Bad value of tauord");
      return 0;
    }
  }

  // Coefficients, t[k], of sine series to convert sigma to tau
  void Geodesic::tauCoeff(double u2, double t[]) throw() {
    double d = u2;
    switch (ntau) {
    case 0:
      break;
    case 1:
      t[0] = -d/8;
      break;
    case 2:
      t[0] = d*(u2-2)/16;
      d *= u2;
      t[1] = -d/256;
      break;
    case 3:
      t[0] = d*((64-37*u2)*u2-128)/1024;
      d *= u2;
      t[1] = d*(u2-1)/256;
      d *= u2;
      t[2] = -d/3072;
      break;
    case 4:
      t[0] = d*(u2*(u2*(47*u2-74)+128)-256)/2048;
      d *= u2;
      t[1] = d*((32-27*u2)*u2-32)/8192;
      d *= u2;
      t[2] = d*(3*u2-2)/6144;
      d *= u2;
      t[3] = -5*d/131072;
      break;
    case 5:
      t[0] = d*(u2*(u2*((752-511*u2)*u2-1184)+2048)-4096)/32768;
      d *= u2;
      t[1] = d*(u2*(u2*(22*u2-27)+32)-32)/8192;
      d *= u2;
      t[2] = d*((384-423*u2)*u2-256)/786432;
      d *= u2;
      t[3] = d*(10*u2-5)/131072;
      d *= u2;
      t[4] = -7*d/1310720;
      break;
    case 6:
      t[0] = d*(u2*(u2*(u2*(u2*(731*u2-1022)+1504)-2368)+4096)-8192)/65536;
      d *= u2;
      t[1] = d*(u2*(u2*((22528-18313*u2)*u2-27648)+32768)-32768)/8388608;
      d *= u2;
      t[2] = d*(u2*(u2*(835*u2-846)+768)-512)/1572864;
      d *= u2;
      t[3] = d*((160-217*u2)*u2-80)/2097152;
      d *= u2;
      t[4] = d*(35*u2-14)/2621440;
      d *= u2;
      t[5] = -7*d/8388608;
      break;
    case 7:
      t[0] = d*(u2*(u2*(u2*(u2*((374272-278701*u2)*u2-523264)+770048)-1212416)+
	2097152)-4194304)/33554432;
      d *= u2;
      t[1] = d*(u2*(u2*(u2*(u2*(15003*u2-18313)+22528)-27648)+32768)-32768)/
	8388608;
      d *= u2;
      t[2] = d*(u2*(u2*((53440-50241*u2)*u2-54144)+49152)-32768)/100663296;
      d *= u2;
      t[3] = d*(u2*(u2*(251*u2-217)+160)-80)/2097152;
      d *= u2;
      t[4] = d*((2240-3605*u2)*u2-896)/167772160;
      d *= u2;
      t[5] = d*(21*u2-7)/8388608;
      d *= u2;
      t[6] = -33*d/234881024;
      break;
    case 8:
      t[0] = d*(u2*(u2*(u2*(u2*(u2*(u2*(428731*u2-557402)+748544)-1046528)+
	1540096)-2424832)+4194304)-8388608)/67108864;
      d *= u2;
      t[1] = d*(u2*(u2*(u2*(u2*((480096-397645*u2)*u2-586016)+720896)-884736)+
	1048576)-1048576)/268435456;
      d *= u2;
      t[2] = d*(u2*(u2*(u2*(u2*(92295*u2-100482)+106880)-108288)+98304)-65536)/
	201326592;
      d *= u2;
      t[3] = d*(u2*(u2*((128512-136971*u2)*u2-111104)+81920)-40960)/
	1073741824.0;
      d *= u2;
      t[4] = d*(u2*(u2*(9555*u2-7210)+4480)-1792)/335544320;
      d *= u2;
      t[5] = d*((672-1251*u2)*u2-224)/268435456;
      d *= u2;
      t[6] = d*(231*u2-66)/469762048;
      d *= u2;
      t[7] = -429*d/17179869184.0;
      break;
    default:
      STATIC_ASSERT(ntau >= 0 && ntau <= 8, "Bad value of ntau");
    }
  }

  // Coefficients, t'[k], of sine series to convert tau to sigma
  void Geodesic::sigCoeff(double u2, double tp[]) throw() {
    double d = u2;
    switch (nsig) {
    case 0:
      break;
    case 1:
      tp[0] = d/8;
      break;
    case 2:
      tp[0] = d*(2-u2)/16;
      d *= u2;
      tp[1] = 5*d/256;
      break;
    case 3:
      tp[0] = d*(u2*(71*u2-128)+256)/2048;
      d *= u2;
      tp[1] = d*(5-5*u2)/256;
      d *= u2;
      tp[2] = 29*d/6144;
      break;
    case 4:
      tp[0] = d*(u2*((142-85*u2)*u2-256)+512)/4096;
      d *= u2;
      tp[1] = d*(u2*(383*u2-480)+480)/24576;
      d *= u2;
      tp[2] = d*(58-87*u2)/12288;
      d *= u2;
      tp[3] = 539*d/393216;
      break;
    case 5:
      tp[0] = d*(u2*(u2*(u2*(20797*u2-32640)+54528)-98304)+196608)/1572864;
      d *= u2;
      tp[1] = d*(u2*((383-286*u2)*u2-480)+480)/24576;
      d *= u2;
      tp[2] = d*(u2*(2907*u2-2784)+1856)/393216;
      d *= u2;
      tp[3] = d*(539-1078*u2)/393216;
      d *= u2;
      tp[4] = 3467*d/7864320;
      break;
    case 6:
      tp[0] = d*(u2*(u2*(u2*((41594-27953*u2)*u2-65280)+109056)-196608)+
	393216)/3145728;
      d *= u2;
      tp[1] = d*(u2*(u2*(u2*(429221*u2-585728)+784384)-983040)+983040)/
	50331648;
      d *= u2;
      tp[2] = d*(u2*((5814-5255*u2)*u2-5568)+3712)/786432;
      d *= u2;
      tp[3] = d*(u2*(111407*u2-86240)+43120)/31457280;
      d *= u2;
      tp[4] = d*(6934-17335*u2)/15728640;
      d *= u2;
      tp[5] = 38081*d/251658240;
      break;
    case 7:
      tp[0] = d*(u2*(u2*(u2*(u2*(u2*(7553633*u2-10733952)+15972096)-25067520)+
	41877504)-75497472)+150994944)/1207959552.0;
      d *= u2;
      tp[1] = d*(u2*(u2*(u2*((429221-314863*u2)*u2-585728)+784384)-983040)+
	983040)/50331648;
      d *= u2;
      tp[2] = d*(u2*(u2*(u2*(1133151*u2-1345280)+1488384)-1425408)+950272)/
	201326592;
      d *= u2;
      tp[3] = d*(u2*((111407-118621*u2)*u2-86240)+43120)/31457280;
      d *= u2;
      tp[4] = d*(u2*(2563145*u2-1664160)+665664)/1509949440.0;
      d *= u2;
      tp[5] = d*(38081-114243*u2)/251658240;
      d *= u2;
      tp[6] = 459485*d/8455716864.0;
      break;
    case 8:
      tp[0] = d*(u2*(u2*(u2*(u2*(u2*((15107266-11062823*u2)*u2-21467904)+
	31944192)-50135040)+83755008)-150994944)+301989888)/2415919104.0;
      d *= u2;
      tp[1] = d*(u2*(u2*(u2*(u2*(u2*(112064929*u2-151134240)+206026080)-
	281149440)+376504320)-471859200)+471859200)/24159191040.0;
      d *= u2;
      tp[2] = d*(u2*(u2*(u2*((2266302-1841049*u2)*u2-2690560)+2976768)-
	2850816)+1900544)/402653184;
      d *= u2;
      tp[3] = d*(u2*(u2*(u2*(174543337*u2-182201856)+171121152)-132464640)+
	66232320)/48318382080.0;
      d *= u2;
      tp[4] = d*(u2*((5126290-6292895*u2)*u2-3328320)+1331328)/3019898880.0;
      d *= u2;
      tp[5] = d*(u2*(45781749*u2-25590432)+8530144)/56371445760.0;
      d *= u2;
      tp[6] = d*(918970-3216395*u2)/16911433728.0;
      d *= u2;
      tp[7] = 109167851*d/5411658792960.0;
      break;
    default:
      STATIC_ASSERT(nsig >= 0 && nsig <= 8, "Bad value of nsig");
    }
  }

  // The scale factor, H, to convert eta to changes in lambda
  double Geodesic::etaFactor(double f, double mu) throw() {
  double g;
  switch (etaord) {
    case 0:
      g = 0;
      break;
    case 1:
      g = -1;
      break;
    case 2:
      g = (f*mu-4)/4;
      break;
    case 3:
      g = (f*(f*(4-3*mu)*mu+4*mu)-16)/16;
      break;
    case 4:
      g = (f*(f*(f*mu*(mu*(25*mu-54)+32)+(32-24*mu)*mu)+32*mu)-128)/128;
      break;
    case 5:
      g = (f*(f*(f*(f*mu*(mu*((720-245*mu)*mu-720)+256)+mu*(mu*(200*mu-432)+
	256))+(256-192*mu)*mu)+256*mu)-1024)/1024;
      break;
    case 6:
      g = (f*(f*(f*(f*(f*mu*(mu*(mu*(mu*(1323*mu-4900)+6800)-4224)+1024)+mu*
	(mu*((2880-980*mu)*mu-2880)+1024))+mu*(mu*(800*mu-1728)+1024))+(1024-
	768*mu)*mu)+1024*mu)-4096)/4096;
      break;
    case 7:
      g = (f*(f*(f*(f*(f*(f*mu*(mu*(mu*(mu*((34020-7623*mu)*mu-60200)+52800)-
	23040)+4096)+mu*(mu*(mu*(mu*(5292*mu-19600)+27200)-16896)+4096))+mu*
	(mu*((11520-3920*mu)*mu-11520)+4096))+mu*(mu*(3200*mu-6912)+4096))+
	(4096-3072*mu)*mu)+4096*mu)-16384)/16384;
      break;
    case 8:
      g = (f*(f*(f*(f*(f*(f*(f*mu*(mu*(mu*(mu*(mu*(mu*(184041*mu-960498)+
	2063880)-2332400)+1459200)-479232)+65536)+mu*(mu*(mu*(mu*((544320-
	121968*mu)*mu-963200)+844800)-368640)+65536))+mu*(mu*(mu*(mu*(84672*mu-
	313600)+435200)-270336)+65536))+mu*(mu*((184320-62720*mu)*mu-184320)+
	65536))+mu*(mu*(51200*mu-110592)+65536))+(65536-49152*mu)*mu)+65536*
	mu)-262144)/262144;
      break;
    default:
      STATIC_ASSERT(etaord >= 0 && etaord <= 8, "Bad value of etaord");
      g = 0;
    }
    return f * g;
  }

  // Coefficients, h[k], of sine series to convert sigma to eta
  void Geodesic::etaCoeff(double f, double mu, double h[]) throw() {
    double s = f * mu, d = s;
    switch (neta) {
    case 0:
      break;
    case 1:
      h[0] = d/8;
      break;
    case 2:
      h[0] = d*(f*(4-3*mu)+4)/32;
      d *= s;
      h[1] = d/64;
      break;
    case 3:
      h[0] = d*(f*(f*(mu*(51*mu-112)+64)-48*mu+64)+64)/512;
      d *= s;
      h[1] = d*(f*(18-13*mu)+8)/512;
      d *= s;
      h[2] = 5*d/1536;
      break;
    case 4:
      h[0] = d*(f*(f*(f*(mu*((764-255*mu)*mu-768)+256)+mu*(204*mu-448)+256)-
	192*mu+256)+256)/2048;
      d *= s;
      h[1] = d*(f*(f*(mu*(79*mu-190)+120)-52*mu+72)+32)/2048;
      d *= s;
      h[2] = d*(f*(72-51*mu)+20)/6144;
      d *= s;
      h[3] = 7*d/8192;
      break;
    case 5:
      h[0] = d*(f*(f*(f*(f*(mu*(mu*(mu*(701*mu-2646)+3724)-2304)+512)+mu*
	((1528-510*mu)*mu-1536)+512)+mu*(408*mu-896)+512)-384*mu+512)+512)/
	4096;
      d *= s;
      h[1] = d*(f*(f*(f*(mu*((1610-487*mu)*mu-1816)+704)+mu*(316*mu-760)+480)-
	208*mu+288)+128)/8192;
      d *= s;
      h[2] = d*(f*(f*(mu*(813*mu-2056)+1360)-408*mu+576)+160)/49152;
      d *= s;
      h[3] = d*(f*(70-49*mu)+14)/16384;
      d *= s;
      h[4] = 21*d/81920;
      break;
    case 6:
      h[0] = d*(f*(f*(f*(f*(f*(mu*(mu*(mu*((74558-16411*mu)*mu-134064)+118720)-
	51200)+8192)+mu*(mu*(mu*(11216*mu-42336)+59584)-36864)+8192)+mu*
	((24448-8160*mu)*mu-24576)+8192)+mu*(6528*mu-14336)+8192)-6144*mu+
	8192)+8192)/65536;
      d *= s;
      h[1] = d*(f*(f*(f*(f*(mu*(mu*(mu*(12299*mu-51072)+80360)-56960)+15360)+
	mu*((25760-7792*mu)*mu-29056)+11264)+mu*(5056*mu-12160)+7680)-3328*mu+
	4608)+2048)/131072;
      d *= s;
      h[2] = d*(f*(f*(f*(mu*((10567-3008*mu)*mu-12712)+5280)+mu*(1626*mu-4112)+
	2720)-816*mu+1152)+320)/98304;
      d *= s;
      h[3] = d*(f*(f*(mu*(485*mu-1266)+860)-196*mu+280)+56)/65536;
      d *= s;
      h[4] = d*(f*(540-375*mu)+84)/327680;
      d *= s;
      h[5] = 11*d/131072;
      break;
    case 7:
      h[0] = d*(f*(f*(f*(f*(f*(f*(mu*(mu*(mu*(mu*(mu*(803251*mu-4262272)+
	9306208)-10659328)+6707200)-2162688)+262144)+mu*(mu*(mu*((2385856-
	525152*mu)*mu-4290048)+3799040)-1638400)+262144)+mu*(mu*(mu*(358912*mu-
	1354752)+1906688)-1179648)+262144)+mu*((782336-261120*mu)*mu-786432)+
	262144)+mu*(208896*mu-458752)+262144)-196608*mu+262144)+262144)/
	2097152;
      d *= s;
      h[1] = d*(f*(f*(f*(f*(f*(mu*(mu*(mu*((1579066-317733*mu)*mu-3149568)+
	3154560)-1587200)+319488)+mu*(mu*(mu*(196784*mu-817152)+1285760)-
	911360)+245760)+mu*((412160-124672*mu)*mu-464896)+180224)+mu*(80896*mu-
	194560)+122880)-53248*mu+73728)+32768)/2097152;
      d *= s;
      h[2] = d*(f*(f*(f*(f*(mu*(mu*(mu*(346689*mu-1534256)+2588016)-1980928)+
	583680)+mu*((676288-192512*mu)*mu-813568)+337920)+mu*(104064*mu-
	263168)+174080)-52224*mu+73728)+20480)/6291456;
      d *= s;
      h[3] = d*(f*(f*(f*(mu*((123082-33633*mu)*mu-154232)+66640)+mu*(15520*mu-
	40512)+27520)-6272*mu+8960)+1792)/2097152;
      d *= s;
      h[4] = d*(f*(f*(mu*(35535*mu-94800)+65520)-12000*mu+17280)+2688)/
	10485760;
      d *= s;
      h[5] = d*(f*(1386-957*mu)+176)/2097152;
      d *= s;
      h[6] = 429*d/14680064;
      break;
    default:
      STATIC_ASSERT(neta >= 0 && neta <= 7, "Bad value of neta");
    }
  }

  // The derivative of etaFactor with respect to mu, dH/dmu
  double Geodesic::etaFactormu(double f, double mu) throw() {
  double g;
  switch (etaord) {
    case 0:
      g = 0;
      break;
    case 1:
      g = 0;
      break;
    case 2:
      g = 1/4.0;
      break;
    case 3:
      g = (f*(2-3*mu)+2)/8;
      break;
    case 4:
      g = (f*(f*(mu*(75*mu-108)+32)-48*mu+32)+32)/128;
      break;
    case 5:
      g = (f*(f*(f*(mu*((540-245*mu)*mu-360)+64)+mu*(150*mu-216)+64)-96*mu+64)+
	64)/256;
      break;
    case 6:
      g = (f*(f*(f*(f*(mu*(mu*(mu*(6615*mu-19600)+20400)-8448)+1024)+mu*((8640-
	3920*mu)*mu-5760)+1024)+mu*(2400*mu-3456)+1024)-1536*mu+1024)+1024)/
	4096;
      break;
    case 7:
      g = (f*(f*(f*(f*(f*(mu*(mu*(mu*((85050-22869*mu)*mu-120400)+79200)-
	23040)+2048)+mu*(mu*(mu*(13230*mu-39200)+40800)-16896)+2048)+mu*
	((17280-7840*mu)*mu-11520)+2048)+mu*(4800*mu-6912)+2048)-3072*mu+2048)+
	2048)/8192;
      break;
    case 8:
      g = (f*(f*(f*(f*(f*(f*(mu*(mu*(mu*(mu*(mu*(1288287*mu-5762988)+10319400)-
	9329600)+4377600)-958464)+65536)+mu*(mu*(mu*((2721600-731808*mu)*mu-
	3852800)+2534400)-737280)+65536)+mu*(mu*(mu*(423360*mu-1254400)+
	1305600)-540672)+65536)+mu*((552960-250880*mu)*mu-368640)+65536)+mu*
	(153600*mu-221184)+65536)-98304*mu+65536)+65536)/262144;
      break;
    default:
      STATIC_ASSERT(etaord >= 0 && etaord <= 8, "Bad value of etaord");
      g = 0;
    }
    return sq(f) * g;
  }

  // Coefficients of sine series to convert sigma to the derivative of eta
  // with respect to mu, dh[k]/dmu
  void Geodesic::etaCoeffmu(double f, double mu, double hp[]) throw() {
    double s = f * mu, d = f;
    switch (neta) {
    case 0:
      break;
    case 1:
      hp[0] = d/8;
      break;
    case 2:
      hp[0] = d*(f*(2-3*mu)+2)/16;
      d *= s;
      hp[1] = d/32;
      break;
    case 3:
      hp[0] = d*(f*(f*(mu*(153*mu-224)+64)-96*mu+64)+64)/512;
      d *= s;
      hp[1] = d*(f*(36-39*mu)+16)/512;
      d *= s;
      hp[2] = 5*d/512;
      break;
    case 4:
      hp[0] = d*(f*(f*(f*(mu*((573-255*mu)*mu-384)+64)+mu*(153*mu-224)+64)-96*
	mu+64)+64)/512;
      d *= s;
      hp[1] = d*(f*(f*(mu*(158*mu-285)+120)-78*mu+72)+32)/1024;
      d *= s;
      hp[2] = d*(f*(18-17*mu)+5)/512;
      d *= s;
      hp[3] = 7*d/2048;
      break;
    case 5:
      hp[0] = d*(f*(f*(f*(f*(mu*(mu*(mu*(3505*mu-10584)+11172)-4608)+512)+mu*
	((4584-2040*mu)*mu-3072)+512)+mu*(1224*mu-1792)+512)-768*mu+512)+512)/
	4096;
      d *= s;
      hp[1] = d*(f*(f*(f*(mu*((6440-2435*mu)*mu-5448)+1408)+mu*(1264*mu-2280)+
	960)-624*mu+576)+256)/8192;
      d *= s;
      hp[2] = d*(f*(f*(mu*(4065*mu-8224)+4080)-1632*mu+1728)+480)/49152;
      d *= s;
      hp[3] = d*(f*(280-245*mu)+56)/16384;
      d *= s;
      hp[4] = 21*d/16384;
      break;
    case 6:
      hp[0] = d*(f*(f*(f*(f*(f*(mu*(mu*(mu*((186395-49233*mu)*mu-268128)+
	178080)-51200)+4096)+mu*(mu*(mu*(28040*mu-84672)+89376)-36864)+4096)+
	mu*((36672-16320*mu)*mu-24576)+4096)+mu*(9792*mu-14336)+4096)-6144*mu+
	4096)+4096)/32768;
      d *= s;
      hp[1] = d*(f*(f*(f*(f*(mu*(mu*(mu*(36897*mu-127680)+160720)-85440)+
	15360)+mu*((51520-19480*mu)*mu-43584)+11264)+mu*(10112*mu-18240)+7680)-
	4992*mu+4608)+2048)/65536;
      d *= s;
      hp[2] = d*(f*(f*(f*(mu*((52835-18048*mu)*mu-50848)+15840)+mu*(8130*mu-
	16448)+8160)-3264*mu+3456)+960)/98304;
      d *= s;
      hp[3] = d*(f*(f*(mu*(1455*mu-3165)+1720)-490*mu+560)+112)/32768;
      d *= s;
      hp[4] = d*(f*(270-225*mu)+42)/32768;
      d *= s;
      hp[5] = 33*d/65536;
      break;
    case 7:
      hp[0] = d*(f*(f*(f*(f*(f*(f*(mu*(mu*(mu*(mu*(mu*(5622757*mu-25573632)+
	46531040)-42637312)+20121600)-4325376)+262144)+mu*(mu*(mu*((11929280-
	3150912*mu)*mu-17160192)+11397120)-3276800)+262144)+mu*(mu*(mu*
	(1794560*mu-5419008)+5720064)-2359296)+262144)+mu*((2347008-1044480*
	mu)*mu-1572864)+262144)+mu*(626688*mu-917504)+262144)-393216*mu+
	262144)+262144)/2097152;
      d *= s;
      hp[1] = d*(f*(f*(f*(f*(f*(mu*(mu*(mu*((9474396-2224131*mu)*mu-15747840)+
	12618240)-4761600)+638976)+mu*(mu*(mu*(1180704*mu-4085760)+5143040)-
	2734080)+491520)+mu*((1648640-623360*mu)*mu-1394688)+360448)+mu*
	(323584*mu-583680)+245760)-159744*mu+147456)+65536)/2097152;
      d *= s;
      hp[2] = d*(f*(f*(f*(f*(mu*(mu*(mu*(2426823*mu-9205536)+12940080)-
	7923712)+1751040)+mu*((3381440-1155072*mu)*mu-3254272)+1013760)+mu*
	(520320*mu-1052672)+522240)-208896*mu+221184)+61440)/6291456;
      d *= s;
      hp[3] = d*(f*(f*(f*(mu*((738492-235431*mu)*mu-771160)+266560)+mu*(93120*
	mu-202560)+110080)-31360*mu+35840)+7168)/2097152;
      d *= s;
      hp[4] = d*(f*(f*(mu*(49749*mu-113760)+65520)-14400*mu+17280)+2688)/
	2097152;
      d *= s;
      hp[5] = d*(f*(8316-6699*mu)+1056)/2097152;
      d *= s;
      hp[6] = 429*d/2097152;
      break;
    default:
      STATIC_ASSERT(neta >= 0 && neta <= 7, "Bad value of neta");
    }
  }
} // namespace GeographicLib


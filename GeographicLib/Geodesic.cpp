/**
 * \file Geodesic.cpp
 * \brief Implementation for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 *
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/Constants.hpp"
#include <algorithm>
#include <limits>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = GEODESIC_HPP;
}

namespace GeographicLib {

  using namespace std;

  Geodesic::Geodesic(double a, double invf)
    : _a(a)
    , _f(1 / invf)
    , _e2(_f * (2 - _f))
    , _ep2(_e2 / (1 - _e2))
    , _b(_a * (1 - _f))
    // Underflow guard.  We require
    //   eps2 * epsilon() > 0
    //   eps2 + epsilon() == epsilon()
    , _eps2(sqrt(numeric_limits<double>::min()))
    , _tol(0.1 * sqrt(numeric_limits<double>::epsilon()))
  {}

  const Geodesic Geodesic::WGS84(Constants::WGS84_a, Constants::WGS84_invf);

  double Geodesic::SinSeries(double x, const double c[], int n) throw() {
    // Evaluate y = sum(c[i - 1] * sin(i * x), i, 1, n) using Clenshaw summation
    double
      ar = 2 * cos(x);
    double			// Accumulators for sum
      y0 = c[n - 1], y1 = 0;
    for (int j = n; --j;) {	// j = n-1 .. 1
      double y2 = y1;
      y1 = y0;
      y0  = ar * y1 - y2 + c[j - 1];
    }
    return sin(x) * y0;
  }

  double Geodesic::AngNormalize(double x) {
    // Place angle in [-180, 180).  Assumes x is in [-540, 540).
    return x >= 180 ? x - 360 : x < -180 ? x + 360 : x;
  }

  // s = b * q(u2) * ( sigma + t(sigma, u2) )
  double Geodesic::sigmaCoeff0(double u2) {
    return (u2*(u2*(u2*(u2*(u2*(u2*((3624192-2760615*u2)*u2-4967424)+7225344)-11468800)+20971520)-50331648)+268435456)+1073741824.)/1073741824.;
  }

  double Geodesic::sigmaCoeff(double sigma, double u2) {
    double d[maxpow];
    double t = u2;
    c[0]=t*(u2*(u2*(u2*(u2*(u2*(u2*(428731*u2-557402)+748544)-1046528)+1540096)-2424832)+4194304)-8388608)/67108864;
    t *= u2;
    c[1]=t*(u2*(u2*(u2*(u2*((480096-397645*u2)*u2-586016)+720896)-884736)+1048576)-1048576)/268435456;
    t *= u2;
    c[2]=t*(u2*(u2*(u2*(u2*(92295*u2-100482)+106880)-108288)+98304)-65536)/201326592;
    t *= u2;
    c[3]=t*(u2*(u2*((128512-136971*u2)*u2-111104)+81920)-40960)/1073741824.;
    t *= u2;
    c[4]=t*(u2*(u2*(9555*u2-7210)+4480)-1792)/335544320;
    t *= u2;
    c[5]=t*((672-1251*u2)*u2-224)/268435456;
    t *= u2
    c[6]=t*(231*u2-66)/469762048;
    t *= u2
    c[7]=-429*t/17179869184.;
    return SinSeries(2 * sigma x, c, maxpow);
  }

  void Geodesic::sCoeff(double u2, double d[]) {
    double t = u2
    _d[0]=t*(u2*(u2*(u2*(u2*(u2*((15107266-11062823*u2)*u2-21467904)+31944192)-50135040)+83755008)-150994944)+301989888)/2415919104.;
    t *= u2;
    _d[1]=t*(u2*(u2*(u2*(u2*(u2*(112064929*u2-151134240)+206026080)-281149440)+376504320)-471859200)+471859200)/24159191040.;
    t *= u2;
    _d[2]=t*(u2*(u2*(u2*((2266302-1841049*u2)*u2-2690560)+2976768)-2850816)+1900544)/402653184;
    t *= u2;
    _d[3]=t*(u2*(u2*(u2*(174543337*u2-182201856)+171121152)-132464640)+66232320)/48318382080.;
    t *= u2;
    _d[4]=t*(u2*((5126290-6292895*u2)*u2-3328320)+1331328)/3019898880.;
    t *= u2;
    _d[5]=t*(u2*(45781749*u2-25590432)+8530144)/56371445760.;
    t *= u2
    _d[6]=(918970-3216395*u2)*t/16911433728.;
    t *= u2
    _d[7]=109167851*t/5411658792960.;
  }
    
  double Geodesic::lamCoeff0(double f, double mu) {
    return  f*(f*(f*(f*(f*(f*(f*(f*mu*(mu*(mu*(mu*(mu*(mu*(184041*mu-960498)+2063880)-2332400)+1459200)-479232)+65536)+mu*(mu*(mu*(mu*((544320-121968*mu)*mu-963200)+844800)-368640)+65536))+mu*(mu*(mu*(mu*(84672*mu-313600)+435200)-270336)+65536))+mu*(mu*((184320-62720*mu)*mu-184320)+65536))+mu*(mu*(51200*mu-110592)+65536))+(65536-49152*mu)*mu)+65536*mu)-262144)/262144;
  }

  void Geodesic::lamCoeff(double f, double mu; double e[]) {
    double s = f*mu, t = s;
    e[0] = (f*(f*(f*(f*(f*(f*(f*(mu*(mu*(mu*(mu*(mu*((30816920-5080225*mu)*mu-79065664)+110840000)-91205632)+43638784)-11010048)+1048576)+mu*(mu*(mu*(mu*(mu*(3213004*mu-17049088)+37224832)-42637312)+26828800)-8650752)+1048576)+mu*(mu*(mu*((9543424-2100608*mu)*mu-17160192)+15196160)-6553600)+1048576)+mu*(mu*(mu*(1435648*mu-5419008)+7626752)-4718592)+1048576)+mu*((3129344-1044480*mu)*mu-3145728)+1048576)+mu*(835584*mu-1835008)+1048576)-786432*mu+1048576)+1048576)*t/8388608; 
    t *= s; 
    e[1] = (f*(f*(f*(f*(f*(f*(mu*(mu*(mu*(mu*(mu*(2092939*mu-12074982)+29005488)-37129344)+26700800)-10207232)+1605632)+mu*(mu*(mu*((6316264-1270932*mu)*mu-12598272)+12618240)-6348800)+1277952)+mu*(mu*(mu*(787136*mu-3268608)+5143040)-3645440)+983040)+mu*((1648640-498688*mu)*mu-1859584)+720896)+mu*(323584*mu-778240)+491520)-212992*mu+294912)+131072)*t/8388608; 
    t *= s; 
    e[2] = (f*(f*(f*(f*(f*(mu*(mu*(mu*((13101384-2474307*mu)*mu-28018000)+30323072)-16658432)+3727360)+mu*(mu*(mu*(1386756*mu-6137024)+10352064)-7923712)+2334720)+mu*((2705152-770048*mu)*mu-3254272)+1351680)+mu*(416256*mu-1052672)+696320)-208896*mu+294912)+81920)*t/25165824; 
    t *= s; 
    e[3] = (f*(f*(f*(f*(mu*(mu*(mu*(273437*mu-1265846)+2238200)-1799088)+557760)+mu*((492328-134532*mu)*mu-616928)+266560)+mu*(62080*mu-162048)+110080)-25088*mu+35840)+7168)*t/8388608; 
    t *= s; 
    e[4] = (f*(f*(f*(mu*((1333160-353765*mu)*mu-1718160)+761600)+mu*(142140*mu-379200)+262080)-48000*mu+69120)+10752)*t/41943040; 
    t *= s; 
    e[5] = (f*(f*(mu*(39633*mu-107426)+75152)-11484*mu+16632)+2112)*t/25165824; 
    t *= s; 
    e[6] = (f*(16016-11011*mu)+1716)*t/58720256; 
    t *= s; 
    e[7] = 715*t/67108864;
  }
(
  GeodesicLine::GeodesicLine(const Geodesic& g,
			     double lat1, double lon1, double bearing1) {
    bearing1 = AngNormalize(bearing1);
    // Normalize bearing at poles.  Evaluate bearings at lat = +/- (90 - eps).
    if (lat1 == 90) {
      lon1 -= bearing1 - (bearing1 >= 0 ? 180 : -180);
      bearing1 = -180;
    } else if (lat1 == -90) {
      lon1 += bearing1;
      bearing1 = 0;
    }
    if (abs(bearing1) < _eps2)
      bearing1 = 0;		// Guard against underflow in sinalpha0
    lon1 = AngNormalize(lon1);
    int bsign = bearing1 >= 0 ? 1 : -1;
    bearing1 *= bsign;
    // alpha is in [0, pi]
    double
      alpha1 = bearing1 * Constants::degree,
      // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
      // problems directly than to skirt them.
      sinalpha1 = bearing1 == 180 ? 0 : sin(alpha1),
      cosalpha1 = bearing1 ==  90 ? 0 : cos(alpha1);
    double cosbeta1, sinbeta1;
    {
      double
	phi1 = lat1 * Constants::degree,
	// Ensure cosbeta1 = +eps at poles
	c = abs(lat1) == 90 ? _eps2 : cos(phi1), s = (1 - _f) * sin(phi1),
	r = hypot(s, c);
      cosbeta1 = c/r, sinbeta1 = s/r;
    }
    double			// alpha0 in [0, pi/2 - |beta|]
      sinalpha0 = sinalpha1 * cosbeta1,
      // Alt: cosalpha0 = hypot(sinbeta1, cosalpha1 * cosbeta1).  The following
      // is slightly better (consider the case sinalpha1 = 0).
      cosalpha0 = hypot(cosalpha1, sinalpha1 * sinbeta1),
      // Evaluate sigma with tan(beta1) = tan(sigma1) * cos(alpha1).
      // sigma = 0 is nearest northward crossing of equator.
      // With beta1 = 0, alpha1 = pi/2, we have sigma1 = 0 (equatorial line).
      // With beta1 =  pi/2, alpha1 = -pi, sigma1 =  pi/2
      // With beta1 = -pi/2, alpha1 =  0 , sigma1 = -pi/2
      sigma1 = atan2(sinbeta1, cosbeta1 * cosalpha1), // sigma1 in (-pi, pi]
      // Evaluate lam1 with tan(lam1) = sin(alpha0) * tan(sigma1).
      // With alpha0 in (0, pi/2], quadrants for sigma and lam coincide.
      // No atan2(0,0) ambiguity at poles since cosbeta1 = +eps.
      // With alpha0 = 0, lam1 = 0 for alpha1 = 0, lam1 = pi for alpha1 = pi.
      lam1 = atan2(sinalpha0 * sinbeta1, cosbeta1 * cosalpha1),
      u2 = sq(cosalpha0) * _ep2,
      q = SigmaScale(u2),
      t = SigmaToT(sigma, u2),
      s0 = sigma + t;

  }

  void GeodesicLine::Position(double s,
			      double& lat, double& lon, double& bearing)
  const throw() {
    s /= (_q * _b);
    s += s0;
    sigma = aa;
    sinsigma2 = sin(sigma2);
    cossigma2 = cos(sigma2);
    sinbeta2 = cosalpha0 * sinsigma2;
    // Alt: cosbeta2 = hypot(cossigma2, sinalpha0 * sinsigma2);
    cosbeta2 = hypot(sinalpha0, cosalpha0 * cossigma2);
    // tan(lam2)=sin(alpha0)*tan(sigma2)
    lam2 = atan2(sinalpha0 * sinsigma2, cossigma2);
    // tan(alpha0)=cos(sigma2)*tan(alpha2)
    alpha2 = -atan2(-sinalpha0, cosalpha0 * cossigma2);

    lat = atan2(sinbet2, (1 - _f) * cosbet2) / Constants::degree;

    bearing = alpha2 / Constants::degree;
tan(alpha0)=cos(sigma)*tan(alpha),

    alpha2 = 

  void Geodesic::Forward(double lat, double lon, double h,
		     double& x, double& y, double& z) const {
    double
      phi = lat * Constants::degree,
      lam = lon * Constants::degree,
      sphi = sin(phi),
      n = _a/sqrt(1 - _e2 * sq(sphi));
    z = ( sq(1 - _f) * n + h) * sphi;
    x = (n + h) * cos(phi);
    y = x * sin(lam);
    x *= cos(lam);
  }

  void Geodesic::Reverse(double x, double y, double z,
		     double& lat, double& lon, double& h) const {
    double
      rad = hypot(x, y),
      p = sq(rad / _a),
      q = (1 - _e2) * sq(z / _a),
      e4 = sq(_e2),
      r = (p + q - e4) / 6;
    double phi;
    h = hypot(rad, z);		// Distance to center of earth
    if (h > _maxrad) {
      // We really far away (> 12 million light years); treat the earth as a
      // point and h, above, is an acceptable approximation to the height.
      // This avoids overflow, e.g., in the computation of disc below.  It's
      // possible that h has overflowed to inf; but that's OK.
      //
      // Treat the case x, y finite, but rad overflows to +inf by scaling by 2.
      phi = atan2(z/2, hypot(x/2, y/2));
    }
    else if ( !(e4 * q == 0 && r <= 0) ) {
      double
	// Avoid possible division by zero when r = 0 by multiplying equations
	// for s and t by r^3 and r, resp.
	S = e4 * p * q / 4,	// S = r^3 * s
	r2 = sq(r),
	r3 = r * r2,
	disc =  S * (2 * r3 + S);
      double u = r;
      if (disc >= 0) {
	double T3 = r3 + S;
	// Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
	// of precision due to cancellation.  The result is unchanged because
	// of the way the T is used in definition of u.
	T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
	// N.B. cbrt always returns the real root.  cbrt(-8) = -2.
	double T = cbrt(T3);	// T = r * t
	// T can be zero; but then r2 / T -> 0.
	u += T + (T != 0 ? r2 / T : 0);
      } else {
	// T is complex, but the way u is defined the result is real.
	double ang = atan2(sqrt(-disc), r3 + S);
	// There are three possible real solutions for u depending on the
	// multiple of 2*pi here.  We choose multiplier = 1 which leads to a
	// jump in the solution across the line 2 + s = 0; but this
	// nevertheless leads to a continuous (and accurate) solution for k.
	// Other choices of the multiplier lead to poorly conditioned solutions
	// near s = 0 (i.e., near p = 0 or q = 0).
	u += 2 * abs(r) * cos((2 * Constants::pi + ang) / 3.0);
      }
      double
	v = sqrt(sq(u) + e4 * q), // guaranteed positive
	// Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
	// e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
	uv = u < 0 ? e4 * q / (v - u) : u + v, //  u + v, guaranteed positive
	// Need to guard against w going negative due to roundoff in uv - q.
	w = max(0.0, _e2 * (uv - q) / (2 * v)),
	// Rearrange expression for k to avoid loss of accuracy due to
	// subtraction.  Division by 0 not possible because uv > 0, w >= 0.
	k = uv / (sqrt(uv + sq(w)) + w), // guaranteed positive
	d = k * rad / (k + _e2);
      // Probably atan2 returns the result for phi more accurately than the
      // half-angle formula that Vermeille uses.  It's certainly simpler.
      phi = atan2(z, d);
      h = (k + _e2 - 1) * hypot(d, z) / k;
    } else {			// e4 * q == 0 && r <= 0
      // Very near equatorial plane with rad <= a * e^2.  This leads to k = 0
      // using the general formula and division by 0 in formula for h.  So
      // handle this case directly.  The condition e4 * q == 0 implies abs(z) <
      // 1.e-145 for WGS84 so it's OK to treat these points as though z = 0.
      // (But we do take care that the sign of phi matches the sign of z.)
      phi = atan2(sqrt( -6 * r), sqrt(p * (1 - _e2)));
      if (z < 0) phi = -phi;	// for tiny negative z
      h = _a * (_e2 - 1) / sqrt(1 - _e2 * sq(sin(phi)));
    }
    lat = phi / Constants::degree;
    // Negative signs return lon in [-180, 180)
    lon = (rad != 0 ? -atan2(-y, x) : 0) / Constants::degree;
  }

#define TOL 1e-13
#define MAXLOOP 11

  void Geodesic::ReverseA(double x, double y, double z,
		     double& lat, double& lon, double& h) const {

    double ae2,p,tgla,tglax,phi,lam,en,slat,clat;
	int count = MAXLOOP;

	ae2 = _a * _e2;
	p = hypot(x, y);
	tgla = z / p / (1. - _e2);
	do {
		tglax = tgla;
		tgla = z / (p - ae2 / sqrt((1 - _e2) * tgla * tgla + 1.));
	} while ((abs(tgla - tglax) > TOL) && --count);
	 { /* convergence achieved */
		phi = atan(tgla);
		slat = sin(phi);
		clat = cos(phi);
		lam = atan2(y, x);
		en = _a / sqrt(1. - _e2 * slat * slat);
		if (abs(phi) <= .7854)
			h = p / clat - en;
		else
			h = z / slat - en + _e2 * en;
		lam = atan2(y, x);
	}
	 lat = phi / Constants::degree;
	 lon = lam / Constants::degree;

  }

} // namespace GeographicLib


/**
 * \file Geodesic.cpp
 * \brief Implementation for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
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

  // Underflow guard.  We require
  //   eps2 * epsilon() > 0
  //   eps2 + epsilon() == epsilon()
  const double Geodesic::eps2 = sqrt(numeric_limits<double>::min());
  const double Geodesic::tol = 100 * numeric_limits<double>::epsilon();

  Geodesic::Geodesic(double a, double invf)
    : _a(a)
    , _f(invf > 0 ? 1 / invf : 0)
    , _f1(1 - _f)
    , _e2(_f * (2 - _f))
    , _ep2(_e2 / (1 - _e2))
    , _b(_a * _f1)
  {}

  const Geodesic Geodesic::WGS84(Constants::WGS84_a, Constants::WGS84_invf);

  double Geodesic::SinSeries(double sinx, double cosx,
			     const double c[], int n) throw() {
    // Evaluate y = sum(c[i - 1] * sin(2 * i * x), i, 1, n) using Clenshaw
    // summation.  (Indices into c offset by 1.)
    double
      ar = 2 * (sq(cosx) - sq(sinx)), // 2 * cos(2 * x)
      y0 = c[n - 1], y1 = 0;	      // Accumulators for sum
    for (int j = n; --j;) {	      // j = n-1 .. 1
      double y2 = y1;
      y1 = y0; y0  = ar * y1 - y2 + c[j - 1];
    }
    return 2 * sinx * cosx * y0; // sin(2 * x) * y0
  }

  // The scale factor to convert tau to s / b
  double Geodesic::tauScale(double u2) throw() {
    return (u2 * (u2 * (u2 * (u2 * (u2 * (u2 * ((3624192 - 2760615 * u2) * u2 - 4967424) + 7225344) - 11468800) + 20971520) - 50331648) + 268435456) + 1073741824.0) / 1073741824.0;
  }

  // Coefficients of sine series to convert sigma to tau (a reversion of
  // tauCoeff).
  void Geodesic::tauCoeff(double u2, double c[]) throw() {
    double t = u2;
    c[0] = t * (u2 * (u2 * (u2 * (u2 * (u2 * (u2 * (428731 * u2 - 557402) + 748544) - 1046528) + 1540096) - 2424832) + 4194304) - 8388608) / 67108864;
    t  *= u2;
    c[1] = t * (u2 * (u2 * (u2 * (u2 * ((480096 - 397645 * u2) * u2 - 586016) + 720896) - 884736) + 1048576) - 1048576) / 268435456;
    t  *= u2;
    c[2] = t * (u2 * (u2 * (u2 * (u2 * (92295 * u2 - 100482) + 106880) - 108288) + 98304) - 65536) / 201326592;
    t  *= u2;
    c[3] = t * (u2 * (u2 * ((128512 - 136971 * u2) * u2 - 111104) + 81920) - 40960) / 1073741824.0;
    t  *= u2;
    c[4] = t * (u2 * (u2 * (9555 * u2 - 7210) + 4480) - 1792) / 335544320;
    t  *= u2;
    c[5] = t * ((672 - 1251 * u2) * u2 - 224) / 268435456;
    t  *= u2;
    c[6] = t * (231 * u2 - 66) / 469762048;
    t  *= u2;
    c[7] = -429 * t / 17179869184.0;
  }

  // Coefficients of sine series to convert tau to sigma (a reversion of
  // tauCoeff).
  void Geodesic::sigCoeff(double u2, double d[]) throw() {
    double t = u2;
    d[0] = t * (u2 * (u2 * (u2 * (u2 * (u2 * ((15107266 - 11062823 * u2) * u2 - 21467904) + 31944192) - 50135040) + 83755008) - 150994944) + 301989888) / 2415919104.0;
    t  *= u2;
    d[1] = t * (u2 * (u2 * (u2 * (u2 * (u2 * (112064929 * u2 - 151134240) + 206026080) - 281149440) + 376504320) - 471859200) + 471859200) / 24159191040.0;
    t  *= u2;
    d[2] = t * (u2 * (u2 * (u2 * ((2266302 - 1841049 * u2) * u2 - 2690560) + 2976768) - 2850816) + 1900544) / 402653184;
    t  *= u2;
    d[3] = t * (u2 * (u2 * (u2 * (174543337 * u2 - 182201856) + 171121152) - 132464640) + 66232320) / 48318382080.0;
    t  *= u2;
    d[4] = t * (u2 * ((5126290 - 6292895 * u2) * u2 - 3328320) + 1331328) / 3019898880.0;
    t  *= u2;
    d[5] = t * (u2 * (45781749 * u2 - 25590432) + 8530144) / 56371445760.0;
    t  *= u2;
    d[6] = (918970 - 3216395 * u2) * t / 16911433728.0;
    t  *= u2;
    d[7] = 109167851 * t / 5411658792960.0;
  }

  double Geodesic::dlamScale(double f, double mu) throw() {
    double g =
      (f * (f * (f * (f * (f * (f * (f * mu * (mu * (mu * (mu * (mu * (mu * (184041 * mu - 960498) + 2063880) - 2332400) + 1459200) - 479232) + 65536) + mu * (mu * (mu * (mu * ((544320 - 121968 * mu) * mu - 963200) + 844800) - 368640) + 65536)) + mu * (mu * (mu * (mu * (84672 * mu - 313600) + 435200) - 270336) + 65536)) + mu * (mu * ((184320 - 62720 * mu) * mu - 184320) + 65536)) + mu * (mu * (51200 * mu - 110592) + 65536)) + (65536 - 49152 * mu) * mu) + 65536 * mu) - 262144) / 262144;
    return f * g;
  }

  double Geodesic::dlamScalemu(double f, double mu) throw() {
    double h = (f * (f * (f * (f * (f * (f * (mu * (mu * (mu * (mu * (mu * (1288287 * mu - 5762988) + 10319400) - 9329600) + 4377600) - 958464) + 65536) + mu * (mu * (mu * ((2721600 - 731808 * mu) * mu - 3852800) + 2534400) - 737280) + 65536) + mu * (mu * (mu * (423360 * mu - 1254400) + 1305600) - 540672) + 65536) + mu * ((552960 - 250880 * mu) * mu - 368640) + 65536) + mu * (153600 * mu - 221184) + 65536) - 98304 * mu + 65536) + 65536) / 262144;
    return h * sq(f);
  }

  void Geodesic::dlamCoeff(double f, double mu, double e[]) throw() {
    double s = f * mu, t = s;
    e[0] = (f * (f * (f * (f * (f * (f * (f * (mu * (mu * (mu * (mu * (mu * ((30816920 - 5080225 * mu) * mu - 79065664) + 110840000) - 91205632) + 43638784) - 11010048) + 1048576) + mu * (mu * (mu * (mu * (mu * (3213004 * mu - 17049088) + 37224832) - 42637312) + 26828800) - 8650752) + 1048576) + mu * (mu * (mu * ((9543424 - 2100608 * mu) * mu - 17160192) + 15196160) - 6553600) + 1048576) + mu * (mu * (mu * (1435648 * mu - 5419008) + 7626752) - 4718592) + 1048576) + mu * ((3129344 - 1044480 * mu) * mu - 3145728) + 1048576) + mu * (835584 * mu - 1835008) + 1048576) - 786432 * mu + 1048576) + 1048576) * t / 8388608;
    t *= s;
    e[1] = (f * (f * (f * (f * (f * (f * (mu * (mu * (mu * (mu * (mu * (2092939 * mu - 12074982) + 29005488) - 37129344) + 26700800) - 10207232) + 1605632) + mu * (mu * (mu * ((6316264 - 1270932 * mu) * mu - 12598272) + 12618240) - 6348800) + 1277952) + mu * (mu * (mu * (787136 * mu - 3268608) + 5143040) - 3645440) + 983040) + mu * ((1648640 - 498688 * mu) * mu - 1859584) + 720896) + mu * (323584 * mu - 778240) + 491520) - 212992 * mu + 294912) + 131072) * t / 8388608;
    t *= s;
    e[2] = (f * (f * (f * (f * (f * (mu * (mu * (mu * ((13101384 - 2474307 * mu) * mu - 28018000) + 30323072) - 16658432) + 3727360) + mu * (mu * (mu * (1386756 * mu - 6137024) + 10352064) - 7923712) + 2334720) + mu * ((2705152 - 770048 * mu) * mu - 3254272) + 1351680) + mu * (416256 * mu - 1052672) + 696320) - 208896 * mu + 294912) + 81920) * t / 25165824;
    t *= s;
    e[3] = (f * (f * (f * (f * (mu * (mu * (mu * (273437 * mu - 1265846) + 2238200) - 1799088) + 557760) + mu * ((492328 - 134532 * mu) * mu - 616928) + 266560) + mu * (62080 * mu - 162048) + 110080) - 25088 * mu + 35840) + 7168) * t / 8388608;
    t *= s;
    e[4] = (f * (f * (f * (mu * ((1333160 - 353765 * mu) * mu - 1718160) + 761600) + mu * (142140 * mu - 379200) + 262080) - 48000 * mu + 69120) + 10752) * t / 41943040;
    t *= s;
    e[5] = (f * (f * (mu * (39633 * mu - 107426) + 75152) - 11484 * mu + 16632) + 2112) * t / 25165824;
    t *= s;
    e[6] = (f * (16016 - 11011 * mu) + 1716) * t / 58720256;
    t *= s;
    e[7] = 715 * t / 67108864;
  }

  void Geodesic::dlamCoeffmu(double f, double mu, double h[]) throw() {
    double s = f * mu, t = f;
    h[0] = (f * (f * (f * (f * (f * (f * (f * (mu * (mu * (mu * (mu * (mu * ((53929610 - 10160450 * mu) * mu - 118598496) + 138550000) - 91205632) + 32729088) - 5505024) + 262144) + mu * (mu * (mu * (mu * (mu * (5622757 * mu - 25573632) + 46531040) - 42637312) + 20121600) - 4325376) + 262144) + mu * (mu * (mu * ((11929280 - 3150912 * mu) * mu - 17160192) + 11397120) - 3276800) + 262144) + mu * (mu * (mu * (1794560 * mu - 5419008) + 5720064) - 2359296) + 262144) + mu * ((2347008 - 1044480 * mu) * mu - 1572864) + 262144) + mu * (626688 * mu - 917504) + 262144) - 393216 * mu + 262144) + 262144) * t / 2097152;
    t *= s;
    h[1] = (f * (f * (f * (f * (f * (f * (mu * (mu * (mu * (mu * (mu * (8371756 * mu - 42262437) + 87016464) - 92823360) + 53401600) - 15310848) + 1605632) + mu * (mu * (mu * ((18948792 - 4448262 * mu) * mu - 31495680) + 25236480) - 9523200) + 1277952) + mu * (mu * (mu * (2361408 * mu - 8171520) + 10286080) - 5468160) + 983040) + mu * ((3297280 - 1246720 * mu) * mu - 2789376) + 720896) + mu * (647168 * mu - 1167360) + 491520) - 319488 * mu + 294912) + 131072) * t / 4194304;
    t *= s;
    h[2] = (f * (f * (f * (f * (f * (mu * (mu * (mu * ((22927422 - 4948614 * mu) * mu - 42027000) + 37903840) - 16658432) + 2795520) + mu * (mu * (mu * (2426823 * mu - 9205536) + 12940080) - 7923712) + 1751040) + mu * ((3381440 - 1155072 * mu) * mu - 3254272) + 1013760) + mu * (520320 * mu - 1052672) + 522240) - 208896 * mu + 221184) + 61440) * t / 6291456;
    t *= s;
    h[3] = (f * (f * (f * (f * (mu * (mu * (mu * (1093748 * mu - 4430461) + 6714600) - 4497720) + 1115520) + mu * ((1476984 - 470862 * mu) * mu - 1542320) + 533120) + mu * (186240 * mu - 405120) + 220160) - 62720 * mu + 71680) + 14336) * t / 4194304;
    t *= s;
    h[4] = (f * (f * (f * (mu * ((466606 - 141506 * mu) * mu - 515448) + 190400) + mu * (49749 * mu - 113760) + 65520) - 14400 * mu + 17280) + 2688) * t / 2097152;
    t *= s;
    h[5] = (f * (f * (mu * (158532 * mu - 375991) + 225456) - 40194 * mu + 49896) + 6336) * t / 12582912;
    t *= s;
    h[6] = (f * (4004 - 3146 * mu) + 429) * t / 2097152;
    t *= s;
    h[7] = 715 * t / 8388608;
  }

  GeodesicLine Geodesic::Line(double lat1, double lon1, double azi1)
    const throw() {
    return GeodesicLine(*this, lat1, lon1, azi1);
  }

  void Geodesic::Direct(double lat1, double lon1, double azi1, double s12,
			double& lat2, double& lon2, double& azi2)
    const throw() {
    GeodesicLine l(*this, lat1, lon1, azi1);
    l.Position(s12, lat2, lon2, azi2);
  }

  void Geodesic::Inverse(double lat1, double lon1, double lat2, double lon2,
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

    double phi, sbet1, cbet1, sbet2, cbet2;

    phi = lat1 * Constants::degree;
    // Ensure cbet1 = +eps at poles
    sbet1 = _f1 * sin(phi);
    cbet1 = lat1 == -90 ? eps2 : cos(phi);
    SinCosNorm(sbet1, cbet1);

    phi = lat2 * Constants::degree;
    // Ensure cbet2 = +eps at poles
    sbet2 = _f1 * sin(phi);
    cbet2 = abs(lat2) == 90 ? eps2 : cos(phi);
    SinCosNorm(sbet2, cbet2);

    double
      // How close to antipodal lat?
      phi12a = (lat2 + lat1) * Constants::degree,
      chi12 = lon12 * Constants::degree,
      cchi12 = cos(chi12),	// lon12 == 90 isn't interesting
      schi12 = lon12 == 180 ? 0 :sin(chi12);

    double calp1, salp1, calp2, salp2, c[maxpow];
    // Enumerate all the cases where the geodesic is a meridian.  This includes
    // coincident points.
    if (schi12 == 0 || lat1 == -90) {
      // Head to the target longitude
      calp1 = cchi12; salp1 = schi12;
      // At the target we're heading north
      calp2 = 1; salp2 = 0;

      double
	// tan(bet) = tan(sig) * cos(alp),
	ssig1 = sbet1, csig1 = calp1 * cbet1,
	ssig2 = sbet2, csig2 = calp2 * cbet2;
      SinCosNorm(ssig1, csig1);
      SinCosNorm(ssig2, csig2);
	
      // sig12 = sig2 - sig1
      double sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, 0.0),
			   csig1 * csig2 + ssig1 * ssig2);

      tauCoeff(_ep2, c);
      s12 = _b * tauScale(_ep2) *
	(sig12 + (SinSeries(ssig2, csig2, c, maxpow) -
		  SinSeries(ssig1, csig1, c, maxpow)));
    } else if (sbet1 == 0 &&	// and sbet2 == 0
	       // Mimic the way Chi12 works with calp1 = 0
	       chi12 <= Constants::pi - _f * Constants::pi) {
      // Geodesic runs along equator
      calp1 = calp2 = 0; salp1 = salp2 = 1;
      s12 = _a * chi12;
    } else {

      // Now point1 and point2 belong within a hemisphere bounded by a line of
      // longitude (lon = lon12/2 +/- 90).

      double sig12, ssig1, csig1, ssig2, csig2, u2;

      double
	chicrita = -cbet1 * dlamScale(_f, sq(sbet1)) * Constants::pi,
	chicrit = Constants::pi - chicrita;
      if (chi12 == chicrit && cbet1 == cbet2 && sbet2 == -sbet1) {
	sig12 = Constants::pi;
	ssig1 = -1; salp1 = salp2 =ssig2 = 1;
	calp1 = calp2 = csig1 = csig2 = 0;
	u2 = sq(sbet1) * _ep2;
      } else {
	if (chi12 > chicrit && phi12a > - chicrita) {
	  salp1 = min(1.0, (Constants::pi - chi12) / chicrita);
	  calp1 = - sqrt(1 - sq(salp1));
	} else {
	  salp1 = 1;
	  calp1 = sbet2 <= 0 ? -eps2 : eps2;
	}

	for (unsigned i = 0, trip = 0; i < 100; ++i) {
	  double dv;
	  double v = Chi12(sbet1, cbet1, sbet2, cbet2,
			   salp1, calp1, salp2, calp2,
			   sig12, ssig1, csig1, ssig2, csig2,
			   u2, trip == 0, dv, c) - chi12;
	  if (v == 0 || trip > 0)
	    break;
	  double
	    dalp1 = -v/dv,
	    sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
	    nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
	  calp1 = calp1 * cdalp1 - salp1 * sdalp1;
	  salp1 = max(0.0, nsalp1);
	  SinCosNorm(salp1, calp1);
	  if (abs(v) < tol) ++trip;
	}
      }	
      tauCoeff(u2, c);
      s12 =  _b * tauScale(u2) *
	(sig12 + (SinSeries(ssig2, csig2, c, maxpow) -
		  SinSeries(ssig1, csig1, c, maxpow)));
    }

    // Convert calp, salp to head accounting for
    // lonsign, swapp, latsign.  The minus signs up result in [-180, 180).

    if (swapp < 0) {
      swap(salp1, salp2);
      swap(calp1, calp2);
    }

    azi1 = -atan2(- swapp * lonsign * salp1,
		   + swapp * latsign * calp1) / Constants::degree;
    azi2 = -atan2(- azi2sense * swapp * lonsign * salp2,
		   + azi2sense * swapp * latsign * calp2) / Constants::degree;
    return;
  }

  double Geodesic::Chi12(double sbet1, double cbet1,
			 double sbet2, double cbet2,
			 double salp1, double calp1,
			 double& salp2, double& calp2,
			 double& sig12,
			 double& ssig1, double& csig1,
			 double& ssig2, double& csig2,
			 double& u2,
			 bool diffp, double& dchi12, double c[])
    const throw() {

    if (sbet1 == 0 && calp1 == 0)
      // Break degeneracy of equatorial line.  This cases has already been
      // handled.
      calp1 = -eps2;

    double
      // sin(alp1) * cos(bet1) = sin(alp0),
      salp0 = salp1 * cbet1,
      calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0

    double slam1, clam1, slam2, clam2, lam12, chi12, mu;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(lam1) = sin(alp0) * tan(sig1).
    ssig1 = sbet1; slam1 = salp0 * sbet1;
    csig1 = clam1 = calp1 * cbet1;
    SinCosNorm(ssig1, csig1);
    SinCosNorm(slam1, clam1);

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
    // tan(lam2) = sin(alp0) * tan(sig2).
    ssig2 = sbet2; slam2 = salp0 * sbet2;
    csig2 = clam2 = calp2 * cbet2;
    SinCosNorm(ssig2, csig2);
    SinCosNorm(slam2, clam2);

    // sig12 = sig2 - sig1, limit to [0, pi]
    sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, 0.0),
		  csig1 * csig2 + ssig1 * ssig2);

    // lam12 = lam2 - lam1, limit to [0, pi]
    lam12 = atan2(max(clam1 * slam2 - slam1 * clam2, 0.0),
		  clam1 * clam2 + slam1 * slam2);

    double eta12, lamscale;
    mu = sq(calp0);
    dlamCoeff(_f, mu, c);
    eta12 = SinSeries(ssig2, csig2, c, maxpow) -
      SinSeries(ssig1, csig1, c, maxpow);
    lamscale = dlamScale(_f, mu),
    chi12 = lam12 + salp0 * lamscale * (sig12 + eta12);

    if (diffp) {
      double dalp0, dsig1, dlam1, dalp2, dsig2, dlam2;
      // Differentiate sin(alp) * cos(bet) = sin(alp0),
      dalp0 = cbet1 * calp1 / calp0;
      dalp2 = calp2 != 0 ? calp1 * cbet1/ (calp2 * cbet2) :
	calp1 >= 0 ? 1 : -1;
      // Differentiate tan(bet) = tan(sig) * cos(alp) and clear
      // calp from the denominator with tan(alp0)=cos(sig)*tan(alp),
      dsig1 = ssig1 * salp0 / calp0;
      dsig2 = ssig2 * salp0 / calp0 * dalp2;
      // Differentiate tan(lam) = sin(alp0) * tan(sig).  Substitute
      //   tan(sig) = tan(bet) / cos(alp) = tan(lam) / sin(alp0)
      //   cos(lam) / cos(sig) = 1 / cos(bet)
      // to give
      dlam1 = (sbet1 * sq(clam1) + slam1 * salp0 / (calp0 * cbet1));
      dlam2 = (sbet2 * sq(clam2) + slam2 * salp0 / (calp0 * cbet2)) * dalp2;

      double deta12, dmu, dlamscale, dchisig;
      dlamCoeffmu(_f, mu, c);
      dmu = - 2 * calp0 * salp0 * dalp0;
      deta12 = dmu * (SinSeries(ssig2, csig2, c, maxpow) -
		      SinSeries(ssig1, csig1, c, maxpow));
      dlamscale = dlamScalemu(_f, mu) * dmu;

      // Derivative of salp0 * lamscale * (sig + eta) wrt sig.  This
      // is from integral form of this expression.
      dchisig =  - _e2 * salp0 *
	(dsig2 / (sqrt(1 - _e2 * (1 - mu * sq(ssig2))) + 1) -
	 dsig1 / (sqrt(1 - _e2 * (1 - mu * sq(ssig1))) + 1)) ;

      dchi12 =
	(dlam2 - dlam1) + dchisig +
	// Derivative wrt mu
	(dalp0 * calp0 * lamscale + salp0 * dlamscale) * (sig12 + eta12) +
	salp0 * lamscale * deta12;
    }

    u2 = mu * _ep2;
    return chi12;
  }

  GeodesicLine::GeodesicLine(const Geodesic& g,
			     double lat1, double lon1, double azi1) {
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
    _bsign = azi1 >= 0 ? 1 : -1;
    azi1 *= _bsign;
    _lat1 = lat1;
    _lon1 = lon1;
    _azi1 = azi1;
    _f1 = g._f1;
    // alp1 is in [0, pi]
    double
      alp1 = azi1 * Constants::degree,
      // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
      // problems directly than to skirt them.
      salp1 = azi1 == 180 ? 0 : sin(alp1),
      calp1 = azi1 ==  90 ? 0 : cos(alp1);
    double cbet1, sbet1, phi;
    phi = lat1 * Constants::degree;
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
    // Evaluate lam1 with tan(lam1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and lam coincide.
    // No atan2(0,0) ambiguity at poles sce cbet1 = +eps.
    // With alp0 = 0, lam1 = 0 for alp1 = 0, lam1 = pi for alp1 = pi.
    _ssig1 = sbet1; _slam1 = _salp0 * sbet1;
    _csig1 = _clam1 = sbet1 != 0 || calp1 != 0 ? cbet1 * calp1 : 1;

    Geodesic::SinCosNorm(_ssig1, _csig1); // sig1 in (-pi, pi]
    Geodesic::SinCosNorm(_slam1, _clam1);
    double
      mu = Geodesic::sq(_calp0),
      u2 = mu * g._ep2;

    _sScale =  g._b * Geodesic::tauScale(u2);
    Geodesic::tauCoeff(u2, _sigCoeff);
    _dtau1 = Geodesic::SinSeries(_ssig1, _csig1, _sigCoeff, maxpow);
    {
      double s = sin(_dtau1), c = cos(_dtau1);
      // tau1 = sig1 + dtau1
      _stau1 = _ssig1 * c + _csig1 * s;
      _ctau1 = _csig1 * c - _ssig1 * s;
    }
    Geodesic::sigCoeff(u2, _sigCoeff);
    // Not necessary because sigCoeff reverts tauCoeff
    //    _dtau1 = -SinSeries(_stau1, _ctau1, _sigCoeff, maxpow);

    _dlamScale = _salp0 * Geodesic::dlamScale(g._f, mu);
    Geodesic::dlamCoeff(g._f, mu, _dlamCoeff);
    _dchi1 = Geodesic::SinSeries(_ssig1, _csig1, _dlamCoeff, maxpow);
  }

  void GeodesicLine::Position(double s12,
			      double& lat2, double& lon2, double& azi2)
  const throw() {
    if (_sScale == 0)
      // Uninitialized
      return;
    double tau12, sig12, lam12, chi12, lon12, s, c;
    double ssig2, csig2, sbet2, cbet2, slam2, clam2, salp2, calp2;
    tau12 = s12 / _sScale;
    s = sin(tau12); c = cos(tau12);
    sig12 = tau12 + (_dtau1 +
		     // tau2 = tau1 + tau12
		     Geodesic::SinSeries(_stau1 * c + _ctau1 * s,
					 _ctau1 * c - _stau1 * s,
					 _sigCoeff, maxpow));
    s = sin(sig12); c = cos(sig12);
    // sig2 = sig1 + sig12
    ssig2 = _ssig1 * c + _csig1 * s;
    csig2 = _csig1 * c - _ssig1 * s;
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = _calp0 * ssig2;
    // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
    cbet2 = Geodesic::hypot(_salp0, _calp0 * csig2);
    // tan(lam2) = sin(alp0) * tan(sig2)
    slam2 = _salp0 * ssig2; clam2 = csig2;  // No need to normalize
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = _salp0; calp2 = _calp0 * csig2; // No need to normalize
    // lam12 = lam2 - lam1
    lam12 = atan2(slam2 * _clam1 - clam2 * _slam1,
		  clam2 * _clam1 + slam2 * _slam1);
    chi12 = lam12 + _dlamScale *
      ( sig12 +
	(Geodesic::SinSeries(ssig2, csig2, _dlamCoeff, maxpow)  - _dchi1));
    lon12 = _bsign * chi12 / Constants::degree;
    // Can't use AngNormalize because longitude might have wrapped multiple
    // times.
    lon12 = lon12 - 360 * floor(lon12/360 + 0.5);
    lat2 = atan2(sbet2, _f1 * cbet2) / Constants::degree;
    lon2 = Geodesic::AngNormalize(_lon1 + lon12);
    azi2 = -atan2(- Geodesic::azi2sense * _bsign * salp2,
		   + Geodesic::azi2sense * calp2) / Constants::degree;
  }

} // namespace GeographicLib


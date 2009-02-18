/**
 * \file Geodesic.cpp
 * \brief Implementation for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 *
 **********************************************************************/

#define REVHEADING 1
#define DEBUG 0
#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/Constants.hpp"
#include <algorithm>
#include <limits>
#if 1
#include <iostream>
#include <iomanip>
#endif

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
  const double Geodesic::tol = 0.1 * sqrt(numeric_limits<double>::epsilon());

  Geodesic::Geodesic(double a, double invf)
    : _a(a)
    , _f(invf > 0 ? 1 / invf : 0)
    , _f1(1 - _f)
    , _e2(_f * (2 - _f))
    , _ep2(_e2 / (1 - _e2))
    , _b(_a * _f1)
  {}

  const Geodesic Geodesic::WGS84(Constants::WGS84_a, Constants::WGS84_invf);

  double Geodesic::SinSeries(double x, const double c[], int n) throw() {
    // Evaluate y = sum(c[i - 1] * sin(2 * i * x), i, 1, n) using Clenshaw
    // summation
    double
      ar = 2 * cos(2 * x),
      y0 = c[n - 1], y1 = 0;	// Accumulators for sum
    for (int j = n; --j;) {	// j = n-1 .. 1
      double y2 = y1;
      y1 = y0;
      y0  = ar * y1 - y2 + c[j - 1];
    }
    return sin(2 * x) * y0;
  }

  double Geodesic::SinSeries(double sinx, double cosx,
			     const double c[], int n) throw() {
    // Evaluate y = sum(c[i - 1] * sin(2 * i * x), i, 1, n) using Clenshaw
    // summation
    double
      ar = 2 * (sq(cosx) - sq(sinx)),	// cos(2 * x)
      y0 = c[n - 1], y1 = 0;	// Accumulators for sum
    for (int j = n; --j;) {	// j = n-1 .. 1
      double y2 = y1;
      y1 = y0;
      y0  = ar * y1 - y2 + c[j - 1];
    }
    return 2 * sinx * cosx * y0; // sin(2 * x) * y0
  }

  // s = b * q(u2) * ( sig + t(sig, u2) )
  double Geodesic::sigScale(double u2) throw() {
    return (u2*(u2*(u2*(u2*(u2*(u2*((3624192-2760615*u2)*u2-4967424)+7225344)-11468800)+20971520)-50331648)+268435456)+1073741824.)/1073741824.;
  }

  void Geodesic::sigCoeffSet(double u2, double c[]) throw()  {
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
    t *= u2;
    c[6]=t*(231*u2-66)/469762048;
    t *= u2;
    c[7]=-429*t/17179869184.;
  }

  void Geodesic::sCoeffSet(double u2, double d[]) throw() {
    double t = u2;
    d[0]=t*(u2*(u2*(u2*(u2*(u2*((15107266-11062823*u2)*u2-21467904)+31944192)-50135040)+83755008)-150994944)+301989888)/2415919104.;
    t *= u2;
    d[1]=t*(u2*(u2*(u2*(u2*(u2*(112064929*u2-151134240)+206026080)-281149440)+376504320)-471859200)+471859200)/24159191040.;
    t *= u2;
    d[2]=t*(u2*(u2*(u2*((2266302-1841049*u2)*u2-2690560)+2976768)-2850816)+1900544)/402653184;
    t *= u2;
    d[3]=t*(u2*(u2*(u2*(174543337*u2-182201856)+171121152)-132464640)+66232320)/48318382080.;
    t *= u2;
    d[4]=t*(u2*((5126290-6292895*u2)*u2-3328320)+1331328)/3019898880.;
    t *= u2;
    d[5]=t*(u2*(45781749*u2-25590432)+8530144)/56371445760.;
    t *= u2;
    d[6]=(918970-3216395*u2)*t/16911433728.;
    t *= u2;
    d[7]=109167851*t/5411658792960.;
  }
    
  double Geodesic::dlamScale(double f, double mu) throw() {
    return  f*(f*(f*(f*(f*(f*(f*(f*mu*(mu*(mu*(mu*(mu*(mu*(184041*mu-960498)+2063880)-2332400)+1459200)-479232)+65536)+mu*(mu*(mu*(mu*((544320-121968*mu)*mu-963200)+844800)-368640)+65536))+mu*(mu*(mu*(mu*(84672*mu-313600)+435200)-270336)+65536))+mu*(mu*((184320-62720*mu)*mu-184320)+65536))+mu*(mu*(51200*mu-110592)+65536))+(65536-49152*mu)*mu)+65536*mu)-262144)/262144;
  }

  void Geodesic::dlamCoeffSet(double f, double mu, double e[]) throw() {
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

  GeodesicLine Geodesic::Line(double lat1, double lon1, double head1)
    const throw() {
    return GeodesicLine(*this, lat1, lon1, head1);
  }
 
  void Geodesic::Direct(double lat1, double lon1, double head1, double s12,
			double& lat2, double& lon2, double& head2)
    const throw() {
    GeodesicLine l(*this, lat1, lon1, head1);
    l.Position(s12, lat2, lon2, head2);
  }

  void Geodesic::Inverse(double lat1, double lon1, double lat2, double lon2,
		 double& s12, double& head1, double& head2)
    const throw() {
    lon1 = AngNormalize(lon1);
    double lon12 = AngNormalize(AngNormalize(lon2) - lon1);
    // If very close to being on the same meridian, then make it so
    // Not sure this is necessary...
    lon12 = AngRound(lon12);
    // Make longitude difference positive
    int lonsign = lon12 >= 0 ? 1 : -1;
    lon12 *= lonsign;
    // If really close to the equator, treat as on equator
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

    double cbet1, sbet1;
    {
      double
	phi = lat1 * Constants::degree,
	// Ensure cbet1 = +eps at poles
	c = lat1 == -90 ? eps2 : cos(phi), s = _f1 * sin(phi),
	r = hypot(s, c);
      cbet1 = c/r, sbet1 = s/r;
    }
    double cbet2, sbet2;
    {
      double
	phi = lat2 * Constants::degree,
	// Ensure cbet2 = +eps at poles
	c = abs(lat2) == 90 ? eps2 : cos(phi), s = _f1 * sin(phi),
	r = hypot(s, c);
      cbet2 = c/r, sbet2 = s/r;
    }

    double
      chi12 = lon12 * Constants::degree,
      cchi12 = cos(chi12),	// lon12 == 90 isn't interesting
      schi12 = lon12 == 180 ? 0 :sin(chi12);

    double calp1, salp1, calp2, salp2;
    double c[maxpow];
    // Enumerate all the cases where the geodesic is a meridian.  This includes
    // coincident points.
    if (schi12 == 0 || lat1 == -90) {
      calp1 = lon12 == 0 || lat1 == -90 ? 1 : -1;
      salp1 = 0;
      // cos(alp)=cos(lam)*cos(alp0)
      // cos(alp0) = 1
      // lam1 = calp1 > 0 ? 0 : pi
      // cos(lam2) = calp1 * cos(chi12)
      calp2 = cchi12 * calp1;
      // We assume alp in [0, pi], like chi12, so no sign change here.
      salp2 = schi12;

      double
	// tan(bet) = tan(sig)*cos(alp),
	sig1 = atan2(sbet1, calp1 * cbet1),
	sig2 = atan2(sbet2, calp2 * cbet2);

      sigCoeffSet(_ep2, c);
      s12 = _b * sigScale(_ep2) *
	((sig2 - sig1) +
	 (SinSeries(sig2, c, maxpow) - SinSeries(sig1, c, maxpow)));
    } else if (sbet1 == 0 && lon12 <= _f1 * 180) { // and sbet2 == 0
      // Geodesic runs along equator
      calp1 = calp2 = 0;
      salp1 = salp2 = 1;
      s12 = _a * chi12;
    } else {

      // Now point1 and point2 belong within a hemisphere bounded by a line of
      // longitude (lon = lon12/2 +/- 90).  Possible sgular cases left to
      // worry about are:

      if (sbet1 == 0) { 	// and sbet2 == 0 and lon12 > _f1 * 180
	// Break degeneracy of points on equator.  This gets sig1 into the
	// right quadrant.
	sbet1 = 0 * -eps2;
      }

      double sig12, ssig1, csig1, ssig2, csig2, u2;
#if DEBUG
      std::cerr << setprecision(20);
#endif
      double
	salp1a = 0, calp1a = 1,
	salp1b = 0, calp1b = -1;
      double
	chidiffa = ChiDiff(sbet1, cbet1,
			   sbet2, cbet2,
			   salp1a, calp1a,
			   salp2, calp2,
			   sig12,
			   ssig1, csig1,
			   ssig2, csig2,
			   u2, c) - chi12,
	chidiffc = ChiDiff(sbet1, cbet1,
			   sbet2, cbet2,
			   salp1b, calp1b,
			   salp2, calp2,
			   sig12,
			   ssig1, csig1,
			   ssig2, csig2,
			   u2, c) - chi12;
      double chidiff;

      for (int i = 0; i < 128; ++i) {
	salp1 = 0.5 * (salp1a + salp1b);
	calp1 = 0.5 * (calp1a + calp1b);
	if (salp1 == 0)
	  salp1 = 1;
	double r = hypot(salp1, calp1);
	salp1 /= r;
	calp1 /= r;

	chidiff = ChiDiff(sbet1, cbet1,
			  sbet2, cbet2,
			  salp1, calp1,
			  salp2, calp2,
			  sig12,
			  ssig1, csig1,
			  ssig2, csig2,
			  u2, c) - chi12;
	//      cerr << "x " << calp1 << " " << chidiff << "\n";
	if (chidiff == 0)
	  break;
	if (chidiff < 0) {
	  if (calp1a == calp1 && salp1a == salp1)
	    break;
	  calp1a = calp1;
	  salp1a = salp1;
	  chidiffa = chidiff;
	} else {
	  if (calp1b == calp1 && salp1b == salp1)
	    break;
	  calp1b = calp1;
	  salp1b = salp1;
	  chidiffc = chidiff;
	}
      }
	
      sigCoeffSet(u2, c);      
      s12 =  _b * sigScale(u2) *
	(sig12 + (SinSeries(ssig2, csig2, c, maxpow) -
		    SinSeries(ssig1, csig1, c, maxpow)));
      cerr << chidiff*cbet2*_a << "\n";
    }
#if DEBUG
    cerr << lonsign << " " << latsign << " " << swapp << "\n";
#endif
    // Convert calp[12], salp[12] to head[12] accounting for
    // lonsign, swapp, latsign.  The minus signs up result in [-180, 180).
#if REVHEADING
    head1 = -atan2(- lonsign * salp1,
		      + latsign * calp1) / Constants::degree;
    head2 = -atan2(+ lonsign * salp2,
		      - latsign * calp2) / Constants::degree;
#else
    head1 = -atan2(- swapp * lonsign * salp1,
		      + swapp * latsign * calp1) / Constants::degree;
    head2 = -atan2(- swapp * lonsign * salp2,
		      + swapp * latsign * calp2) / Constants::degree;
#endif
    if (swapp < 0)
      swap(head1, head2);
    return;
  }

  double Geodesic::ChiDiff(double sbet1, double cbet1,
			   double sbet2, double cbet2,
			   double salp1, double calp1,
			   double& salp2, double& calp2,
			   double& sig12,
			   double& ssig1, double& csig1,
			   double& ssig2, double& csig2,
			   double& u2,
			   double c[]) const throw() {
      // Pick salp1 > 0 calp1 in (-1, 1), i.e., alp1 in (0, pi).  If
      // sbet1 == 0, we pick calp1 in (-1, 0), i.e., alp in (pi/2, pi)

	if (calp1 == 0)
	  // If sbet1 = sbet2 = 0, and calp1 = calp2 = 0, we want
	  // to ensure that lam1 = -pi, lam2 = 0.
	  calp1 = 0 * sqrt(eps2);

      double
	// Follow GeodesicLine constructor
	salp0 = salp1 * cbet1,
	calp0 = hypot(calp1, salp1 * sbet1);
      
      double slam1, clam1, slam2, clam2, r;
      ssig1 = sbet1;
      csig1 = clam1 = sbet1 != 0 || calp1 > 0 ? calp1 * cbet1 : -1;
      slam1 = salp0 * sbet1;
      r = hypot(ssig1, csig1);
      ssig1 /= r; csig1 /= r;
      r = hypot(slam1, clam1);
      slam1 /= r; clam1 /= r;

      salp2 = salp0 / cbet2;
      // calp2 = sqrt(1 - sq(salp2))
      //           = sqrt(sq(calp0) - sq(sbet2)) / cbet2
      // and subst for calp0 and rearrange to give (choose positive sqrt
      // to give alp2 in [0, pi/2]).  N.B. parens around
      //    sq(sbet1) - sq(sbet2)
      // are needed to maintain accuracy when calph1 is small.
      calp2 = sqrt(sq(calp1 * cbet1) +
		   (sbet1 - sbet2) * (sbet1 + sbet2)) / cbet2;
      ssig2 = sbet2;
      csig2 = clam2 = sbet2 != 0 || calp2 != 0 ? calp2 * cbet2 : 1;
      slam2 = salp0 * sbet2;
      r = hypot(ssig2, csig2);
      ssig2 /= r; csig2 /= r;
      r = hypot(slam2, clam2);
      slam2 /= r; clam2 /= r;

      sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, 0.0),
		      csig1 * csig2 + ssig1 * ssig2);

      double lam12 =
	atan2(max(clam1 * slam2 - slam1 * clam2, 0.0),
	      clam1 * clam2 + slam1 * slam2);
      double mu = sq(calp0);
      u2 = mu * _ep2;
      dlamCoeffSet(_f, mu, c);
      double 
	xchi12 = lam12 +
	salp0 * dlamScale(_f, mu) *
	(sig12 + (SinSeries(ssig2, csig2, c, maxpow) -
		    SinSeries(ssig1, csig1, c, maxpow)));

#if DEBUG
      /*
      std::cerr << salp0 << "\n" << calp0 << "\n"
		<< salp1 << "\n" << calp1 << "\n"
		<< salp2 << "\n" << calp2 << "\n"
		<< ssig1 << "\n" << csig1 << "\n"
		<< ssig2 << "\n" << csig2 << "\n"
		<< slam1 << "\n" << clam1 << "\n"
		<< slam2 << "\n" << clam2 << "\n";
      */
      std::cerr << salp1 << " " << calp1 << " "
		<< atan2(salp1, calp1)/Constants::degree << " "
		<< xchi12/Constants::degree << " "
		<< sig12/Constants::degree << " "
		<< lam12/Constants::degree << "\n";
#endif
      // Compare xchi12 to chi12
      return xchi12;
  }
  
  GeodesicLine::GeodesicLine(const Geodesic& g,
			     double lat1, double lon1, double head1) {
    head1 = Geodesic::AngNormalize(head1);
    // Normalize head at poles.  Evaluate heads at lat = +/- (90 - eps).
    if (lat1 == 90) {
      lon1 -= head1 - (head1 >= 0 ? 180 : -180);
      head1 = -180;
    } else if (lat1 == -90) {
      lon1 += head1;
      head1 = 0;
    }
    // Guard against underflow in salp0
    head1 = Geodesic::AngRound(head1);
    lon1 = Geodesic::AngNormalize(lon1);
    _bsign = head1 >= 0 ? 1 : -1;
    head1 *= _bsign;
    _lat1 = lat1;
    _lon1 = lon1;
    _head1 = head1;
    _f1 = g._f1;
    // alp1 is in [0, pi]
    double
      alp1 = head1 * Constants::degree,
      // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
      // problems directly than to skirt them.
      salp1 = head1 == 180 ? 0 : sin(alp1),
      calp1 = head1 ==  90 ? 0 : cos(alp1);
    double cbet1, sbet1;
    {
      double
	phi = lat1 * Constants::degree,
	// Ensure cbet1 = +eps at poles
	c = abs(lat1) == 90 ? Geodesic::eps2 : cos(phi), s = _f1 * sin(phi),
	r = hypot(s, c);
      cbet1 = c/r, sbet1 = s/r;
    }
    _salp0 = salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    _calp0 = hypot(calp1, salp1 * sbet1);
    double
      // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
      // sig = 0 is nearest northward crossing of equator.
      // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
      // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
      // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
      sig1 = atan2(sbet1, cbet1 * calp1), // sig1 in (-pi, pi]
      // Evaluate lam1 with tan(lam1) = sin(alp0) * tan(sig1).
      // With alp0 in (0, pi/2], quadrants for sig and lam coincide.
      // No atan2(0,0) ambiguity at poles sce cbet1 = +eps.
      // With alp0 = 0, lam1 = 0 for alp1 = 0, lam1 = pi for alp1 = pi.
      lam1 = atan2(_salp0 * sbet1, cbet1 * calp1),
      mu = Geodesic::sq(_calp0),
      u2 = mu * g._ep2;

    _sScale =  g._b * Geodesic::sigScale(u2);
    Geodesic::sigCoeffSet(u2, _sigCoeff);
    _S1 = sig1 + Geodesic::SinSeries(sig1, _sigCoeff, maxpow);
    Geodesic::sCoeffSet(u2, _sigCoeff);

    _dlamScale = _salp0 * Geodesic::dlamScale(g._f, mu);
    Geodesic::dlamCoeffSet(g._f, mu, _dlamCoeff);
    _chi1 = lam1 +
      _dlamScale * (sig1 +
		       Geodesic::SinSeries(sig1, _dlamCoeff, maxpow));
  }

  void GeodesicLine::Position(double s12,
			      double& lat2, double& lon2, double& head2)
  const throw() {
    if (_sScale == 0)
      // Uninitialized
      return;
    double
      S2 = _S1 + s12 / _sScale,
      sig2 = S2 +  Geodesic::SinSeries(S2, _sigCoeff, maxpow),
      ssig2 = sin(sig2),
      csig2 = cos(sig2),
      sbet2 = _calp0 * ssig2,
      // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
      cbet2 = hypot(_salp0, _calp0 * csig2),
      // tan(lam2)=sin(alp0)*tan(sig2)
      lam2 = atan2(_salp0 * ssig2, csig2),
      // tan(alp0)=cos(sig2)*tan(alp2)
      alp2 = atan2(_salp0, _calp0 * csig2),
      chi2 = lam2 +
      _dlamScale * (sig2 +
		       Geodesic::SinSeries(sig2, _dlamCoeff, maxpow));
    lat2 = atan2(sbet2, _f1 * cbet2) / Constants::degree;
    lon2 = Geodesic::AngNormalize(_lon1 +
				  _bsign * (chi2 - _chi1) / Constants::degree);
#if REVHEADING
    head2 = Geodesic::AngNormalize(_bsign * alp2 / Constants::degree + 180);
#else
    head2 = Geodesic::AngNormalize(_bsign * alp2 / Constants::degree);
#endif
  }
 
} // namespace GeographicLib


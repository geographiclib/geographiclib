/**
 * \file Geodesic.cpp
 * \brief Implementation for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 *
 **********************************************************************/

#define REVBEARING 1
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

  // s = b * q(u2) * ( sigma + t(sigma, u2) )
  double Geodesic::sigmaScale(double u2) throw() {
    return (u2*(u2*(u2*(u2*(u2*(u2*((3624192-2760615*u2)*u2-4967424)+7225344)-11468800)+20971520)-50331648)+268435456)+1073741824.)/1073741824.;
  }

  void Geodesic::sigmaCoeffSet(double u2, double c[]) throw()  {
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
    
  double Geodesic::dlambdaScale(double f, double mu) throw() {
    return  f*(f*(f*(f*(f*(f*(f*(f*mu*(mu*(mu*(mu*(mu*(mu*(184041*mu-960498)+2063880)-2332400)+1459200)-479232)+65536)+mu*(mu*(mu*(mu*((544320-121968*mu)*mu-963200)+844800)-368640)+65536))+mu*(mu*(mu*(mu*(84672*mu-313600)+435200)-270336)+65536))+mu*(mu*((184320-62720*mu)*mu-184320)+65536))+mu*(mu*(51200*mu-110592)+65536))+(65536-49152*mu)*mu)+65536*mu)-262144)/262144;
  }

  void Geodesic::dlambdaCoeffSet(double f, double mu, double e[]) throw() {
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

  GeodesicLine Geodesic::Line(double lat1, double lon1, double bearing1)
    const throw() {
    return GeodesicLine(*this, lat1, lon1, bearing1);
  }
 
  void Geodesic::Direct(double lat1, double lon1, double bearing1, double s12,
			double& lat2, double& lon2, double& bearing2)
    const throw() {
    GeodesicLine l(*this, lat1, lon1, bearing1);
    l.Position(s12, lat2, lon2, bearing2);
  }

  void Geodesic::Inverse(double lat1, double lon1, double lat2, double lon2,
		 double& s12, double& bearing1, double& bearing2)
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

    double cosbeta1, sinbeta1;
    {
      double
	phi = lat1 * Constants::degree,
	// Ensure cosbeta1 = +eps at poles
	c = lat1 == -90 ? eps2 : cos(phi), s = _f1 * sin(phi),
	r = hypot(s, c);
      cosbeta1 = c/r, sinbeta1 = s/r;
    }
    double cosbeta2, sinbeta2;
    {
      double
	phi = lat2 * Constants::degree,
	// Ensure cosbeta2 = +eps at poles
	c = abs(lat2) == 90 ? eps2 : cos(phi), s = _f1 * sin(phi),
	r = hypot(s, c);
      cosbeta2 = c/r, sinbeta2 = s/r;
    }

    double
      chi12 = lon12 * Constants::degree,
      coschi12 = cos(chi12),	// lon12 == 90 isn't interesting
      sinchi12 = lon12 == 180 ? 0 :sin(chi12);

    double cosalpha1, sinalpha1, cosalpha2, sinalpha2;
    double c[maxpow];
    // Enumerate all the cases where the geodesic is a meridian.  This includes
    // coincident points.
    if (sinchi12 == 0 || lat1 == -90) {
      cosalpha1 = lon12 == 0 || lat1 == -90 ? 1 : -1;
      sinalpha1 = 0;
      // cos(alpha)=cos(lam)*cos(alpha0)
      // cos(alpha0) = 1
      // lam1 = cosalpha1 > 0 ? 0 : pi
      // cos(lam2) = cosalpha1 * cos(chi12)
      cosalpha2 = coschi12 * cosalpha1;
      // We assume alpha in [0, pi], like chi12, so no sign change here.
      sinalpha2 = sinchi12;

      double
	// tan(beta) = tan(sigma)*cos(alpha),
	sigma1 = atan2(sinbeta1, cosalpha1 * cosbeta1),
	sigma2 = atan2(sinbeta2, cosalpha2 * cosbeta2);

      sigmaCoeffSet(_ep2, c);
      s12 = _b * sigmaScale(_ep2) *
	((sigma2 - sigma1) +
	 (SinSeries(sigma2, c, maxpow) - SinSeries(sigma1, c, maxpow)));
    } else if (sinbeta1 == 0 && lon12 <= _f1 * 180) { // and sinbeta2 == 0
      // Geodesic runs along equator
      cosalpha1 = cosalpha2 = 0;
      sinalpha1 = sinalpha2 = 1;
      s12 = _a * chi12;
    } else {

      // Now point1 and point2 belong within a hemisphere bounded by a line of
      // longitude (lon = lon12/2 +/- 90).  Possible singular cases left to
      // worry about are:

      if (sinbeta1 == 0) { 	// and sinbeta2 == 0 and lon12 > _f1 * 180
	// Break degeneracy of points on equator.  This gets sigma1 into the
	// right quadrant.
	sinbeta1 = 0 * -eps2;
      }

      double sigma1, sigma2, u2;
#if DEBUG
      std::cerr << setprecision(20);
#endif
      double
	sinalpha1a = 0, cosalpha1a = 1,
	sinalpha1c = 0, cosalpha1c = -1;
      double
	chidiffa = ChiDiff(sinbeta1, cosbeta1,
			   sinbeta2, cosbeta2,
			   sinalpha1a, cosalpha1a,
			   sinalpha2, cosalpha2,
			   sigma1, sigma2, u2, c) - chi12,
	chidiffc = ChiDiff(sinbeta1, cosbeta1,
			   sinbeta2, cosbeta2,
			   sinalpha1c, cosalpha1c,
			   sinalpha2, cosalpha2,
			   sigma1, sigma2, u2, c) - chi12;
      double chidiff;

      for (int i = 0; i < 128; ++i) {
	sinalpha1 = 0.5 * (sinalpha1a + sinalpha1c);
	cosalpha1 = 0.5 * (cosalpha1a + cosalpha1c);
	if (sinalpha1 == 0)
	  sinalpha1 = 1;
	double r = hypot(sinalpha1, cosalpha1);
	//	sinalpha1 /= r;
	//	cosalpha1 /= r;
	// cerr << i << " ";
	chidiff = ChiDiff(sinbeta1, cosbeta1,
			   sinbeta2, cosbeta2,
			   sinalpha1/r, cosalpha1/r,
			   sinalpha2, cosalpha2,
			   sigma1, sigma2, u2, c) - chi12;
	//      cerr << "x " << cosalpha1 << " " << chidiff << "\n";
	if (chidiff == 0)
	  break;
	if (chidiff < 0) {
	  //	  if (cosalpha1a == cosalpha1 && sinalpha1a == sinalpha1)
	  //	    break;
	  cosalpha1a = cosalpha1;
	  sinalpha1a = sinalpha1;
	  chidiffa = chidiff;
	} else {
	  //	  if (cosalpha1c == cosalpha1 && sinalpha1c == sinalpha1)
	  //	    break;
	  cosalpha1c = cosalpha1;
	  sinalpha1c = sinalpha1;
	  chidiffc = chidiff;
	}
      }
	
      {
	double r = hypot(sinalpha1, cosalpha1);
	sinalpha1 /= r;
	cosalpha1 /= r;
      }

      sigmaCoeffSet(u2, c);      
      s12 =  _b * sigmaScale(u2) *
	((sigma2 - sigma1) +
	 (SinSeries(sigma2, c, maxpow) - SinSeries(sigma1, c, maxpow)));
      cerr << chidiff*cosbeta2*_a << "\n";
    }
#if DEBUG
    cerr << lonsign << " " << latsign << " " << swapp << "\n";
#endif
    // Convert cosalpha[12], sinalpha[12] to bearing[12] accounting for
    // lonsign, swapp, latsign.  The minus signs up result in [-180, 180).
#if REVBEARING
    bearing1 = -atan2(- lonsign * sinalpha1,
		      + latsign * cosalpha1) / Constants::degree;
    bearing2 = -atan2(+ lonsign * sinalpha2,
		      - latsign * cosalpha2) / Constants::degree;
#else
    bearing1 = -atan2(- swapp * lonsign * sinalpha1,
		      + swapp * latsign * cosalpha1) / Constants::degree;
    bearing2 = -atan2(- swapp * lonsign * sinalpha2,
		      + swapp * latsign * cosalpha2) / Constants::degree;
#endif
    if (swapp < 0)
      swap(bearing1, bearing2);
    return;
  }

  double Geodesic::ChiDiff(double sinbeta1, double cosbeta1,
			   double sinbeta2, double cosbeta2,
			   double sinalpha1, double cosalpha1,
			   double& sinalpha2, double& cosalpha2,
			   double& sigma1, double& sigma2,
			   double& u2,
			   double c[]) const throw() {
      // Pick sinalpha1 > 0 cosalpha1 in (-1, 1), i.e., alpha1 in (0, pi).  If
      // sinbeta1 == 0, we pick cosalpha1 in (-1, 0), i.e., alpha in (pi/2, pi)

	if (cosalpha1 == 0)
	  // If sinbeta1 = sinbeta2 = 0, and cosalpha1 = cosalpha2 = 0, we want
	  // to ensure that lambda1 = -pi, lambda2 = 0.
	  cosalpha1 = 0 * sqrt(eps2);

      double
	// Follow GeodesicLine constructor
	sinalpha0 = sinalpha1 * cosbeta1,
	cosalpha0 = hypot(cosalpha1, sinalpha1 * sinbeta1);
      
      double lambda1;
      if (sinbeta1 == 0 && cosalpha1 <= 0)
	  sigma1 = lambda1 = -Constants::pi;
      else {
	sigma1 = atan2(sinbeta1, cosbeta1 * cosalpha1);
	lambda1 = atan2(sinalpha0 * sinbeta1, cosbeta1 * cosalpha1);
      }

      sinalpha2 = sinalpha0 / cosbeta2;
      // cosalpha2 = sqrt(1 - sq(sinalpha2))
      //           = sqrt(sq(cosalpha0) - sq(sinbeta2)) / cosbeta2
      // and subst for cosalpha0 and rearrange to give (choose positive sqrt
      // to give alpha2 in [0, pi/2])
      cosalpha2 = sqrt(sq(cosalpha1 * cosbeta1) +
		       sq(sinbeta1) - sq(sinbeta2)) / cosbeta2;
      sigma2 = atan2(sinbeta2, cosbeta2 * cosalpha2);
      double lambda2 = atan2(sinalpha0 * sinbeta2, cosbeta2 * cosalpha2);

      double mu = sq(cosalpha0);
      u2 = mu * _ep2;
      dlambdaCoeffSet(_f, mu, c);
      double 
	xchi12 = lambda2 - lambda1 +
	sinalpha0 * dlambdaScale(_f, mu) *
	((sigma2 - sigma1) +
	 (SinSeries(sigma2, c, maxpow) - SinSeries(sigma1, c, maxpow)));

#if DEBUG
      std::cerr << sinalpha1 << " " << cosalpha1 << " "
		<< atan2(sinalpha1, cosalpha1)/Constants::degree << " "
		<< xchi12/Constants::degree << " "
		<< (sigma2 - sigma1)/Constants::degree << " "
		<< (lambda2 - lambda1)/Constants::degree << "\n";
#endif
      // Compare xchi12 to chi12
      return xchi12;
  }
  
  GeodesicLine::GeodesicLine(const Geodesic& g,
			     double lat1, double lon1, double bearing1) {
    bearing1 = Geodesic::AngNormalize(bearing1);
    // Normalize bearing at poles.  Evaluate bearings at lat = +/- (90 - eps).
    if (lat1 == 90) {
      lon1 -= bearing1 - (bearing1 >= 0 ? 180 : -180);
      bearing1 = -180;
    } else if (lat1 == -90) {
      lon1 += bearing1;
      bearing1 = 0;
    }
    // Guard against underflow in sinalpha0
    bearing1 = Geodesic::AngRound(bearing1);
    lon1 = Geodesic::AngNormalize(lon1);
    _bsign = bearing1 >= 0 ? 1 : -1;
    bearing1 *= _bsign;
    _lat1 = lat1;
    _lon1 = lon1;
    _bearing1 = bearing1;
    _f1 = g._f1;
    // alpha1 is in [0, pi]
    double
      alpha1 = bearing1 * Constants::degree,
      // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
      // problems directly than to skirt them.
      sinalpha1 = bearing1 == 180 ? 0 : sin(alpha1),
      cosalpha1 = bearing1 ==  90 ? 0 : cos(alpha1);
    double cosbeta1, sinbeta1;
    {
      double
	phi = lat1 * Constants::degree,
	// Ensure cosbeta1 = +eps at poles
	c = abs(lat1) == 90 ? Geodesic::eps2 : cos(phi), s = _f1 * sin(phi),
	r = hypot(s, c);
      cosbeta1 = c/r, sinbeta1 = s/r;
    }
    _sinalpha0 = sinalpha1 * cosbeta1; // alpha0 in [0, pi/2 - |beta1|]
    // Alt: cosalpha0 = hypot(sinbeta1, cosalpha1 * cosbeta1).  The following
    // is slightly better (consider the case sinalpha1 = 0).
    _cosalpha0 = hypot(cosalpha1, sinalpha1 * sinbeta1);
    double
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
      lambda1 = atan2(_sinalpha0 * sinbeta1, cosbeta1 * cosalpha1),
      mu = Geodesic::sq(_cosalpha0),
      u2 = mu * g._ep2;

    _sScale =  g._b * Geodesic::sigmaScale(u2);
    Geodesic::sigmaCoeffSet(u2, _sigmaCoeff);
    _S1 = sigma1 + Geodesic::SinSeries(sigma1, _sigmaCoeff, maxpow);
    Geodesic::sCoeffSet(u2, _sigmaCoeff);

    _dlambdaScale = _sinalpha0 * Geodesic::dlambdaScale(g._f, mu);
    Geodesic::dlambdaCoeffSet(g._f, mu, _dlambdaCoeff);
    _chi1 = lambda1 +
      _dlambdaScale * (sigma1 +
		       Geodesic::SinSeries(sigma1, _dlambdaCoeff, maxpow));
  }

  void GeodesicLine::Position(double s12,
			      double& lat2, double& lon2, double& bearing2)
  const throw() {
    if (_sScale == 0)
      // Uninitialized
      return;
    double
      S2 = _S1 + s12 / _sScale,
      sigma2 = S2 +  Geodesic::SinSeries(S2, _sigmaCoeff, maxpow),
      sinsigma2 = sin(sigma2),
      cossigma2 = cos(sigma2),
      sinbeta2 = _cosalpha0 * sinsigma2,
      // Alt: cosbeta2 = hypot(cossigma2, sinalpha0 * sinsigma2);
      cosbeta2 = hypot(_sinalpha0, _cosalpha0 * cossigma2),
      // tan(lambda2)=sin(alpha0)*tan(sigma2)
      lambda2 = atan2(_sinalpha0 * sinsigma2, cossigma2),
      // tan(alpha0)=cos(sigma2)*tan(alpha2)
      alpha2 = atan2(_sinalpha0, _cosalpha0 * cossigma2),
      chi2 = lambda2 +
      _dlambdaScale * (sigma2 +
		       Geodesic::SinSeries(sigma2, _dlambdaCoeff, maxpow));
    lat2 = atan2(sinbeta2, _f1 * cosbeta2) / Constants::degree;
    lon2 = Geodesic::AngNormalize(_lon1 +
				  _bsign * (chi2 - _chi1) / Constants::degree);
#if REVBEARING
    bearing2 = Geodesic::AngNormalize(_bsign * alpha2 / Constants::degree + 180);
#else
    bearing2 = Geodesic::AngNormalize(_bsign * alpha2 / Constants::degree);
#endif
  }
 
} // namespace GeographicLib


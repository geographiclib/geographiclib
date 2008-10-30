/**
 * \file TransverseMercator.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/TransverseMercator.hpp"
#include "GeographicLib/Constants.hpp"
#include <limits>

/*
 * Implementation taken from the report:
 *
 * JHS 154, ETRS89 - j‰rjestelm‰‰n liittyv‰t karttaprojektiot,
 * tasokoordinaatistot ja karttalehtijako (Map projections, plane coordinates,
 * and map sheet index for ETRS89), Published by JUHTA, Finnish Geodetic
 * Institute, and the National Land Survey of Finland, 34 p (2006).
 *
 * http://www.jhs-suositukset.fi/suomi/jhs154
 * http://docs.jhs-suositukset.fi/jhs-suositukset/JHS154/JHS154.pdf
 *
 * This achieves approximately 1um accuracy.  This is a straight transcription
 * of the formulas in this paper with the following exceptions:
 *
 *  * use Newton's method instead of plain iteration to solve for latitude in
 *    terms of isometric latitude in the Reverse method
 *
 *  * use more accurate expressions for the convergence and scale.
 *
 * Other accurate implementations are given in
 * http://www.ign.fr/telechargement/MPro/geodesie/CIRCE/NTG_76.pdf
 * http://www.lantmateriet.se/upload/filer/kartor/geodesi_gps_och_detaljmatning/geodesi/Formelsamling/Gauss_Conformal_Projection.pdf
 */

/*
 *
 * Define mu = asin(sin(lam) * cos(phi))
 *           = angular distance from meridian
 *
 * Errors are primarily a function of mu (or x).
 *
 * For each set, define
 *
 *    dxm  = max(erra, errb)
 *           for mth order method (order n^m or e^(2*m)), where
 *
 *           erra = the error in the forward transformation scaled to
 *                  distance on the ground with the scale factor k
 *
 *           errb = the discrepancy in applying the forward transformation
 *                  followed by the reverse transformation and converting
 *                  the result to a distance
 *
 *    dgam = max error in meridian convergence using the 6th order method
 *
 *    dk   = max relative error in scale using the 6th order method
 *
 * Units:
 *
 *    1um = 1e-6 m
 *    d = degrees, ' = minutes, " = seconds
 *    % = 0.01, %% = 0.001
 *
 *     set         dx4     dx5     dx6     dx7     dx8      dgam        dk
 * x<4e5, y<95e5 200nm   5.0nm   5.0nm   5.0nm   5.0nm    6e-11"   2e-12%%
 * x<5e5, y<96e5 210nm   5.0nm   5.0nm   5.0nm   5.0nm    6e-11"   2e-12%%
 *     mu<10     350nm   5.1nm   5.0nm   5.0nm   5.0nm    1e-10"   2e-12%%
 *     mu<15     700nm   6.5nm   5.0nm   5.0nm   5.0nm    1e-10"   2e-12%%
 *     mu<20     1.5um    11nm   5.0nm   5.0nm   5.0nm    1e-10"   2e-12%%
 *     mu<25     3.3um    23nm   5.0nm   5.0nm   5.0nm    2e-10"   2e-12%%
 *     mu<30     7.6um    62nm   5.0nm   5.0nm   5.0nm    4e-10"   2e-12%%
 *     mu<35      18um   180nm   5.0nm   5.0nm   5.0nm   1.0e-9"   6e-12%%
 *     mu<40      47um   570nm    10nm   5.0nm   5.0nm   4.1e-9"   2e-11%%
 *     mu<45     130um   2.0um    35nm   5.0nm   5.0nm   2.0e-8"   1e-10%%
 *     mu<50     400um   8.0um   170nm   6.3nm   5.0nm   1.1e-7"   6e-10%%
 *     mu<55     1.4mm    37um   1.1um    33nm   5.0nm   8.0e-7"  3.8e-9%%
 *     mu<60     5.8mm   210um   8.4um   350nm    17nm   7.1e-6"  3.5e-8%%
 *     mu<65      31mm   1.6mm    94um   5.7um   360nm   9.8e-5"  4.7e-7%%
 *     mu<70     230mm    20mm   1.8mm   170um    17um   2.2e-3"  1.1e-5%%
 *     mu<72     600mm    62mm   6.9mm   820um   100um    0.010"  4.8e-5%%
 *     mu<74     1.8m    230mm    33mm   4.9mm   750um    0.055"  2.7e-4%%
 *     mu<76     6.2m    1.1m    200mm    39mm   7.9mm     0.39"  1.8e-3%%
 *     mu<78      27m    6.3m    1.6m    430mm   .12m       3.7"   0.017%%
 *     mu<80     160m     55m     20m    7.9m    3.2m        55"    0.28%%
 *     mu<82     1.5km   870m    520m    330m    210m        35'     7.9%%
 *     mu<84      27km    28km    37km    53km    81km       19d      39%
 *
 */

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = TRANSVERSEMERCATOR_HPP;
}

#if defined(_MSC_VER)
#define hypot _hypot
#endif

namespace GeographicLib {

  const double TransverseMercator::tol =
    0.1*sqrt(std::numeric_limits<double>::epsilon());

  TransverseMercator::TransverseMercator(double a, double invf, double k0)
    : _a(a)
    , _f(1 / invf)
    , _k0(k0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(_e2))
    , _e2m(1 - _e2)
    , _ep2(_e2 / _e2m)
    , _n(_f / (2 - _f))
  {
#if TM_TX_MAXPOW <= 4
    _b1 = 1/(1+_n)*(sq(_n)*(sq(_n)+16)+64)/64;
    _h[0] = _n*(_n*((555-4*_n)*_n-960)+720)/1440;
    _h[1] = sq(_n)*((96-437*_n)*_n+30)/1440;
    _h[2] = (119-148*_n)*sq(_n)*_n/3360;
    _h[3] = 4397*sq(sq(_n))/161280;
    _hp[0] = _n*(_n*(_n*(164*_n+225)-480)+360)/720;
    _hp[1] = sq(_n)*(_n*(557*_n-864)+390)/1440;
    _hp[2] = (427-1236*_n)*sq(_n)*_n/1680;
    _hp[3] = 49561*sq(sq(_n))/161280;
#elif TM_TX_MAXPOW == 5
    _b1 = 1/(1+_n)*(sq(_n)*(sq(_n)+16)+64)/64;
    _h[0] = _n*(_n*(_n*((-3645*_n-64)*_n+8880)-15360)+11520)/23040;
    _h[1] = sq(_n)*(_n*(_n*(4416*_n-3059)+672)+210)/10080;
    _h[2] = sq(_n)*_n*((-627*_n-592)*_n+476)/13440;
    _h[3] = (4397-3520*_n)*sq(sq(_n))/161280;
    _h[4] = 4583*sq(sq(_n))*_n/161280;
    _hp[0] = _n*(_n*(_n*((328-635*_n)*_n+450)-960)+720)/1440;
    _hp[1] = sq(_n)*(_n*(_n*(4496*_n+3899)-6048)+2730)/10080;
    _hp[2] = sq(_n)*_n*(_n*(15061*_n-19776)+6832)/26880;
    _hp[3] = (49561-171840*_n)*sq(sq(_n))/161280;
    _hp[4] = 34729*sq(sq(_n))*_n/80640;
#elif TM_TX_MAXPOW == 6
    _b1 = 1/(1+_n)*(sq(_n)*(sq(_n)*(sq(_n)+4)+64)+256)/256;
    _h[0] = _n*(_n*(_n*(_n*(_n*(384796*_n-382725)-6720)+932400)-1612800)+
	       1209600)/2419200;
    _h[1] = sq(_n)*(_n*(_n*((1695744-1118711*_n)*_n-1174656)+258048)+80640)/
      3870720;
    _h[2] = sq(_n)*_n*(_n*(_n*(22276*_n-16929)-15984)+12852)/362880;
    _h[3] = sq(sq(_n))*((-830251*_n-158400)*_n+197865)/7257600;
    _h[4] = (453717-435388*_n)*sq(sq(_n))*_n/15966720;
    _h[5] = 20648693*sq(sq(_n))*sq(_n)/638668800;
    _hp[0] = _n*(_n*(_n*(_n*(_n*(31564*_n-66675)+34440)+47250)-100800)+75600)/
      151200;
    _hp[1] = sq(_n)*(_n*(_n*((863232-1983433*_n)*_n+748608)-1161216)+524160)/
      1935360;
    _hp[2] = sq(_n)*_n*(_n*(_n*(670412*_n+406647)-533952)+184464)/725760;
    _hp[3] = sq(sq(_n))*(_n*(6601661*_n-7732800)+2230245)/7257600;
    _hp[4] = (3438171-13675556*_n)*sq(sq(_n))*_n/7983360;
    _hp[5] = 212378941*sq(sq(_n))*sq(_n)/319334400;
#elif TM_TX_MAXPOW == 7
    _b1 = 1/(1+_n)*(sq(_n)*(sq(_n)*(sq(_n)+4)+64)+256)/256;
    _h[0] = _n*(_n*(_n*(_n*(_n*((6156736-5406467*_n)*_n-6123600)-107520)+
		       14918400)-25804800)+19353600)/38707200;
    _h[1] = sq(_n)*(_n*(_n*(_n*(_n*(829456*_n-5593555)+8478720)-5873280)+
		      1290240)+403200)/19353600;
    _h[2] = sq(_n)*_n*(_n*(_n*(_n*(9261899*_n+3564160)-2708640)-2557440)+
		     2056320)/58060800;
    _h[3] = sq(sq(_n))*(_n*(_n*(14928352*_n-9132761)-1742400)+2176515)/
      79833600;
    _h[4] = sq(sq(_n))*_n*((-8005831*_n-1741552)*_n+1814868)/63866880;
    _h[5] = (268433009-261810608*_n)*sq(sq(_n))*sq(_n)/8302694400.;
    _h[6] = 219941297*sq(sq(_n))*sq(_n)*_n/5535129600.;
    _hp[0] = _n*(_n*(_n*(_n*(_n*(_n*(1804025*_n+2020096)-4267200)+2204160)+
		       3024000)-6451200)+4838400)/9676800;
    _hp[1] = sq(_n)*(_n*(_n*(_n*(_n*(4626384*_n-9917165)+4316160)+3743040)-
		      5806080)+2620800)/9676800;
    _hp[2] = sq(_n)*_n*(_n*(_n*((26816480-67102379*_n)*_n+16265880)-21358080)+
		     7378560)/29030400;
    _hp[3] = sq(sq(_n))*(_n*(_n*(155912000*_n+72618271)-85060800)+24532695)/
      79833600;
    _hp[4] = sq(sq(_n))*_n*(_n*(102508609*_n-109404448)+27505368)/63866880;
    _hp[5] = (2760926233.-12282192400.*_n)*sq(sq(_n))*sq(_n)/4151347200.;
    _hp[6] = 1522256789.*sq(sq(_n))*sq(_n)*_n/1383782400.;
#elif TM_TX_MAXPOW >= 8
    _b1 = 1/(1+_n)*(sq(_n)*(sq(_n)*(sq(_n)*(25*sq(_n)+64)+256)+4096)+16384)/
      16384;
    _h[0] = _n*(_n*(_n*(_n*(_n*(_n*(_n*(31777436*_n-37845269)+43097152)-
			       42865200)-752640)+104428800)-180633600)+
	       135475200)/270950400;
    _h[1] = sq(_n)*(_n*(_n*(_n*(_n*(_n*(24749483*_n+14930208)-100683990)+
			      152616960)-105719040)+23224320)+7257600)/
      348364800;
    _h[2] = sq(_n)*_n*(_n*(_n*(_n*((101880889-232468668*_n)*_n+39205760)-
			     29795040)-28131840)+22619520)/638668800;
    _h[3] = sq(sq(_n))*(_n*(_n*(_n*(324154477*_n+1433121792.)-876745056)-
			    167270400)+208945440)/7664025600.;
    _h[4] = sq(sq(_n))*_n*(_n*(_n*(457888660*_n-312227409)-67920528)+
			   70779852)/2490808320.;
    _h[5] = sq(sq(_n))*sq(_n)*((-19841813847.*_n-3665348512.)*_n+3758062126.)/
      116237721600.;
    _h[6] = (1979471673.-1989295244.*_n)*sq(sq(_n))*sq(_n)*_n/49816166400.;
    _h[7] = 191773887257.*sq(sq(sq(_n)))/3719607091200.;
    _hp[0] = _n*(_n*(_n*(_n*(_n*(_n*((37884525-75900428*_n)*_n+42422016)-
			       89611200)+46287360)+63504000)-135475200)+
	       101606400)/203212800;
    _hp[1] = sq(_n)*(_n*(_n*(_n*(_n*(_n*(148003883*_n+83274912)-178508970)+
			      77690880)+67374720)-104509440)+47174400)/
      174182400;
    _hp[2] = sq(_n)*_n*(_n*(_n*(_n*(_n*(318729724*_n-738126169)+294981280)+
			     178924680)-234938880)+81164160)/319334400;
    _hp[3] = sq(sq(_n))*(_n*(_n*((14967552000.-40176129013.*_n)*_n+
				6971354016.)-8165836800.)+2355138720.)/
      7664025600.;
    _hp[4] = sq(sq(_n))*_n*(_n*(_n*(10421654396.*_n+3997835751.)-4266773472.)+
			   1072709352.)/2490808320.;
    _hp[5] = sq(sq(_n))*sq(_n)*(_n*(175214326799.*_n-171950693600.)+
			      38652967262.)/58118860800.;
    _hp[6] = (13700311101.-67039739596.*_n)*sq(sq(_n))*sq(_n)*_n/12454041600.;
    _hp[7] = 1424729850961.*sq(sq(sq(_n)))/743921418240.;
#endif
    // _a1 is the equivalent radius for computing the circumference of
    // ellipse.  Relative error is f^6/16384 = 8.8e-20 for WGS84.
    _a1 = _b1 * _a;
  }

  const TransverseMercator
  TransverseMercator::UTM(Constants::WGS84_a, Constants::WGS84_invf,
			  Constants::UTM_k0);

  void TransverseMercator::Forward(double lon0, double lat, double lon,
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
    // q is isometric latitude
    // JHS 154 has
    //
    //   beta = atan(sinh(q)) = conformal latitude
    //   [xi', eta'] = spheroidal TM coordinates
    //   eta' = atanh(cos(beta) * sin(lam))
    //   xi' = asin(sin(beta)*cosh(eta')
    //
    // We use
    //
    //   tan(beta) = sinh(q)
    //   sin(beta) = tanh(q)
    //   cos(beta) = sech(q)
    //   denom^2    = 1-cos(beta)^2*sin(lam)^2 = 1-sech(q)^2*sin(lam)^2
    //   sin(xip)   = sin(beta)/denom          = tanh(q)/denom
    //   cos(xip)   = cos(beta)*cos(lam)/denom = sech(q)*cos(lam)/denom
    //   cosh(etap) = 1/denom                  = 1/denom
    //   sinh(etap) = cos(beta)*sin(lam)/denom = sech(q)*sin(lam)/denom
    //
    // to eliminate beta and derive more stable expressions for xi',eta'
    double etap, xip;
    if (lat < 90) {
      double
	qp = asinh(tan(phi)),
	qpp = atanh(_e * sin(phi)),
	q = qp - _e * qpp;
      xip = atan2(sinh(q), cos(lam));
      etap = atanh(sin(lam) / cosh(q));
      // convergence and scale for speroidal TM (xip, etap) -- gamma0 =
      // atan(tan(xip) * tanh(etap)) = atan(tan(lam) * sin(beta))
      gamma = atan(tan(lam) * tanh(q));
      // k0 = sqrt(1 - _e2 * sin(phi)^2) * (cos(beta) / cos(phi)) * cosh(etap)
      // Note 1/cos(phi) = cosh(qp);
      // and cos(beta) * cosh(etap) = 1/hypot(sinh(q), cos(lam))
      k = sqrt(_e2m + _e2 * sq(cos(phi))) * cosh(qp) / hypot(sinh(q), cos(lam));
    } else {
      xip = Constants::pi/2;
      etap = 0;
      gamma = lam;
      // See, for example, Lee (1976), p 100.
      k = sqrt( std::pow(1 + _e, 1 + _e) * std::pow(1 - _e, 1 - _e) );
    }
    // {xi',eta'} is {northing,easting} for spheroidal transverse mercator (for
    // eta' = 0, xi' = bet). {xi,eta} is {northing,easting} for transverse
    // mercator with constant scale on the central meridian (for eta = 0, xip =
    // rectifying latitude).  Define
    //
    //   zeta = xi + i*eta
    //   zeta' = xi' + i*eta'
    //
    // The conversion from conformal to rectifying latitude can be expresses as
    // a series in _n:
    //
    //   zeta = zeta' + sum(h[j-1]' * sin(2 * j * zeta'), j = 1..maxpow)
    //
    // where h[j]' = O(_n^j).  The reversion of this series gives
    //
    //   zeta' = zeta - sum(h[j-1] * sin(2 * j * zeta), j = 1..maxpow)
    //
    // which is used in Reverse.
    //
    // Evaluate sums via Clenshaw method.  See
    //    http://mathworld.wolfram.com/ClenshawRecurrenceFormula.html
    //
    // Let
    //
    //    S = sum(c[k] * F[k](x), k = 0..N)
    //    F[n+1](x) = alpha(n,x) * F[n](x) + beta(n,x) * F[n-1](x)
    //
    // Evaluate S with
    //
    //    y[N+2] = y[N+1] = 0
    //    y[k] = alpha(k,x) * y[k+1] + beta(k+1,x) * y[k+2] + c[k]
    //    S = c[0] * F[0](x) + y[1] * F[1](x) + beta(1,x) * F[0](x) * y[2]
    //
    // Here we have
    //
    //    x = 2 * zeta'
    //    F[n](x) = sin(n * x)
    //    a(n, x) = 2 * cos(x)
    //    b(n, x) = -1
    //    [ sin(A+B) - 2*cos(B)*sin(A) + sin(A-B) = 0, A = n*x, B = x ]
    //    N = maxpow
    //    c[k] = _hp[k-1]
    //    S = y[1] * sin(x)
    //
    // For the derivative we have
    //
    //    x = 2 * zeta'
    //    F[n](x) = cos(n * x)
    //    a(n, x) = 2 * cos(x)
    //    b(n, x) = -1
    //    [ cos(A+B) - 2*cos(B)*cos(A) + cos(A-B) = 0, A = n*x, B = x ]
    //    c[0] = 1; c[k] = 2*k*_hp[k-1]
    //    S = (c[0] - y[2])  + y[1] * cos(x)
    double
      c0 = cos(2 * xip), ch0 = cosh(2 * etap),
      s0 = sin(2 * xip), sh0 = sinh(2 * etap),
      ar = 2 * c0 * ch0, ai = -2 * s0 * sh0; // 2 * cos(2*zeta')
    double			// Accumulators for zeta
      xi0 = _hp[maxpow - 1], eta0 = 0,
      xi1 = 0, eta1 = 0,
      xi2, eta2;
    double			// Accumulators for dzeta/dzeta'
      yr0 = 2 * maxpow * _hp[maxpow - 1], yi0 = 0,
      yr1 = 0, yi1 = 0,
      yr2, yi2;
    for (int j = maxpow; --j;) { // j = maxpow-1 .. 1
      xi2 = xi1; eta2 = eta1; yr2 = yr1; yi2 = yi1;
      xi1 = xi0; eta1 = eta0; yr1 = yr0; yi1 = yi0;
      xi0  = ar * xi1 - ai * eta1 - xi2 + _hp[j - 1];
      eta0 = ai * xi1 + ar * eta1 - eta2;
      yr0 = ar * yr1 - ai * yi1 - yr2 + 2 * j * _hp[j - 1];
      yi0 = ai * yr1 + ar * yi1 - yi2;
    }
    ar /= 2; ai /= 2;		// cos(2*zeta')
    yr2 = 1 - yr1 + ar * yr0 - ai * yi0;
    yi2 =   - yi1 + ai * yr0 + ar * yi0;
    ar = s0 * ch0; ai = c0 * sh0; // sin(2*zeta')
    double
      xi  = xip  + ar * xi0 - ai * eta0,
      eta = etap + ai * xi0 + ar * eta0;
    // Fold in change in convergence and scale for spheroidal TM to UTM.
    gamma -= atan2(yi2, yr2);
    k *= _b1 * hypot(yr2, yi2);
    gamma /= Constants::degree;
    y = _a1 * _k0 * (backside ? Constants::pi - xi : xi) * latsign;
    x = _a1 * _k0 * eta * lonsign;
    if (backside)
      gamma = 180 - gamma;
    gamma *= latsign * lonsign;
    k *= _k0;
  }

  void TransverseMercator::Reverse(double lon0, double x, double y,
				   double& lat, double& lon,
				   double& gamma, double& k) const {
    // This undoes the steps in Forward.  The wrinkles are: (1) Use of the
    // reverted series to express zeta' in terms of zeta. (2) Newton's method
    // to solve for phi in terms of q.
    double
      xi = y / (_a1 * _k0),
      eta = x / (_a1 * _k0);
    // Explicitly enforce the parity
    int
      xisign = xi < 0 ? -1 : 1,
      etasign = eta < 0 ? -1 : 1;
    xi *= xisign;
    eta *= etasign;
    bool backside = xi > Constants::pi/2;
    if (backside)
      xi = Constants::pi - xi;
    double
      c0 = cos(2 * xi), ch0 = cosh(2 * eta),
      s0 = sin(2 * xi), sh0 = sinh(2 * eta),
      ar = 2 * c0 * ch0, ai = -2 * s0 * sh0; // 2 * cos(2*zeta)
    double			// Accumulators for zeta'
      xip0 = -_h[maxpow - 1], etap0 = 0,
      xip1 = 0, etap1 = 0,
      xip2, etap2;
    double			// Accumulators for dzeta'/dzeta
      yr0 = - 2 * maxpow * _h[maxpow - 1], yi0 = 0,
      yr1 = 0, yi1 = 0,
      yr2, yi2;
    for (int j = maxpow; --j;) { // j = maxpow-1 .. 1
      xip2 = xip1; etap2 = etap1; yr2 = yr1; yi2 = yi1;
      xip1 = xip0; etap1 = etap0; yr1 = yr0; yi1 = yi0;
      xip0  = ar * xip1 - ai * etap1 - xip2 - _h[j - 1];
      etap0 = ai * xip1 + ar * etap1 - etap2;
      yr0 = ar * yr1 - ai * yi1 - yr2 - 2 * j * _h[j - 1];
      yi0 = ai * yr1 + ar * yi1 - yi2;
    }
    ar /= 2; ai /= 2;		// cos(2*zeta')
    yr2 = 1 - yr1 + ar * yr0 - ai * yi0;
    yi2 =   - yi1 + ai * yr0 + ar * yi0;
    ar = s0 * ch0; ai = c0 * sh0; // sin(2*zeta)
    double
      xip  = xi  + ar * xip0 - ai * etap0,
      etap = eta + ai * xip0 + ar * etap0;
    // Convergence and scale for spheroidal TM to UTM.
    gamma = atan2(yi2, yr2);
    k = _b1 / hypot(yr2, yi2);
    // JHS 154 has
    //
    // 	 beta = asin(sin(xip) / cosh(etap))
    // 	 lam = asin(tanh(etap) / cos(beta)
    // 	 q = asinh(tan(beta))
    //
    // the following eliminates beta and is more stable
    double lam, phi;
    double
      s = sinh(etap),
      c = cos(xip),
      r = hypot(s, c);
    if (r > 0) {
      lam = atan2(s, c);
      // Solve
      // q = qp - e * atanh(e * tanh(qp))
      // for qp = asinh(tan(phi))
      double
	q = asinh(sin(xip)/r),
	qp = q;
      for (int i = 0; i < numit; ++i) {
	double
	  t = tanh(qp),
	  dqp = -(qp - _e * atanh(_e * t) - q) *
	  (1 - _e2 * sq(t)) / _e2m;
	qp += dqp;
	if (std::abs(dqp) < tol)
	  break;
      }
      phi = atan(sinh(qp));
      gamma += atan(tan(xip) * tanh(etap));
      // Note cos(beta) * cosh(etap) = r
      k *= sqrt(_e2m + _e2 * sq(cos(phi))) * cosh(qp) * r;
    } else {
      phi = Constants::pi/2;
      lam = 0;
      k *= sqrt( std::pow(1 + _e, 1 + _e) * std::pow(1 - _e, 1 - _e) );
    }
    lat = phi / Constants::degree * xisign;
    lon = lam / Constants::degree;
    if (backside)
      lon = 180 - lon;
    lon *= etasign;
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon + lon0 >= 180)
      lon += lon0 - 360;
    else if (lon + lon0 < -180)
      lon += lon0 + 360;
    else
      lon += lon0;
    gamma /= Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= xisign * etasign;
    k *= _k0;
  }

} // namespace GeographicLib

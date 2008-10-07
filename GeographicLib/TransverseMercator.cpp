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
 *           for transformation and 4th order formula for the convergence
 * 
 *    dk   = max relative error in scale using the 6th order method for
 *           transformation and 4th order formula for the scale
 * 
 * Units:
 * 
 *    1um = 1e-6 m
 *    d = degrees, ' = minutes, " = seconds
 *    % = 0.01, %% = 0.001
 * 
 *     set         dx4     dx5     dx6     dx7     dx8   dgam   dk
 * x<4e5, y<95e5 200nm   5.0nm   5.0nm   5.0nm   5.0nm   2e-8" 3e-10%%
 * x<5e5, y<96e5 210nm   5.0nm   5.0nm   5.0nm   5.0nm   2e-8" 3e-10%%
 *     mu<10     350nm   5.1nm   5.0nm   5.0nm   5.0nm   4e-8" 4e-10%%
 *     mu<15     700nm   6.5nm   5.0nm   5.0nm   5.0nm   8e-8" 5e-10%%
 *     mu<20     1.5um    11nm   5.0nm   5.0nm   5.0nm   2e-7" 8e-10%%
 *     mu<25     3.3um    23nm   5.0nm   5.0nm   5.0nm   4e-7" 2e-9%%
 *     mu<30     7.6um    62nm   5.0nm   5.0nm   5.0nm   7e-7" 4e-9%%
 *     mu<35      18um   180nm   5.0nm   5.0nm   5.0nm   2e-6" 8e-9%%
 *     mu<40      47um   570nm    10nm   5.0nm   5.0nm   4e-6" 2e-8%%
 *     mu<45     130um   2.0um    35nm   5.0nm   5.0nm   2e-5" 6e-8%%
 *     mu<50     400um   8.0um   170nm   6.3nm   5.0nm   4e-5" 2e-7%%
 *     mu<55     1.4mm    37um   1.1um    33nm   5.0nm   2e-4" 7e-7%%
 *     mu<60     5.8mm   210um   8.4um   350nm    17nm   7e-4" 3e-6%%
 *     mu<65      31mm   1.6mm    94um   5.7um   360nm   4e-3" 2e-5%%
 *     mu<70     230mm    20mm   1.8mm   170um    17um   0.04" 2e-4%%
 *     mu<72     600mm    62mm   6.9mm   820um   100um   0.2"  6e-4%%
 *     mu<74     1.8m    230mm    33mm   4.9mm   750um   0.5"  2e-3%%
 *     mu<76     6.2m    1.1m    200mm    39mm   7.9mm   2"    9e-3%%
 *     mu<78      27m    6.3m    1.6m    430mm   .12m    10"   0.05%%
 *     mu<80     160m     55m     20m    7.9m    3.2m    84"   0.5%%
 *     mu<82     1.5km   870m    520m    330m    210m    32'   8%%
 *     mu<84      27km    28km    37km    53km    81km   17d   54%
 *
 */

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = TRANSVERSEMERCATOR_HPP;
}

namespace GeographicLib {

  TransverseMercator::TransverseMercator(double a, double invf, double k0)
    : _a(a)
    , _f(1 / invf)
    , _k0(k0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(_e2))
    , _e2m(1 - _e2)
    , _ep2(_e2 / _e2m)
    , _n(_f / (2 - _f))
      // _a1 is the equivalent radius for computing the circumference of
      // ellipse.  Relative error is f^6/16384 = 8.8e-20 for WGS84.
    , _tol(0.1*sqrt(std::numeric_limits<double>::epsilon()))
    , _numit(5)
  {
#if TM_TX_MAXPOW <= 4
    _a1 =_a/(1+_n)*(_n*_n*(_n*_n+16)+64)/64;
    _h[0] =_n*(_n*((555-4*_n)*_n-960)+720)/1440;
    _h[1] =_n*_n*((96-437*_n)*_n+30)/1440;
    _h[2] =(119-148*_n)*_n*_n*_n/3360;
    _h[3] =4397*_n*_n*_n*_n/161280;
    _hp[0]=_n*(_n*(_n*(164*_n+225)-480)+360)/720;
    _hp[1]=_n*_n*(_n*(557*_n-864)+390)/1440;
    _hp[2]=(427-1236*_n)*_n*_n*_n/1680;
    _hp[3]=49561*_n*_n*_n*_n/161280;
#elif TM_TX_MAXPOW == 5
    _a1 =_a/(1+_n)*(_n*_n*(_n*_n+16)+64)/64;
    _h[0] =_n*(_n*(_n*((-3645*_n-64)*_n+8880)-15360)+11520)/23040;
    _h[1] =_n*_n*(_n*(_n*(4416*_n-3059)+672)+210)/10080;
    _h[2] =_n*_n*_n*((-627*_n-592)*_n+476)/13440;
    _h[3] =(4397-3520*_n)*_n*_n*_n*_n/161280;
    _h[4] =4583*_n*_n*_n*_n*_n/161280;
    _hp[0]=_n*(_n*(_n*((328-635*_n)*_n+450)-960)+720)/1440;
    _hp[1]=_n*_n*(_n*(_n*(4496*_n+3899)-6048)+2730)/10080;
    _hp[2]=_n*_n*_n*(_n*(15061*_n-19776)+6832)/26880;
    _hp[3]=(49561-171840*_n)*_n*_n*_n*_n/161280;
    _hp[4]=34729*_n*_n*_n*_n*_n/80640;
#elif TM_TX_MAXPOW == 6
    _a1 =_a/(1+_n)*(_n*_n*(_n*_n*(_n*_n+4)+64)+256)/256;
    _h[0] =_n*(_n*(_n*(_n*(_n*(384796*_n-382725)-6720)+932400)-1612800)+
	       1209600)/2419200;
    _h[1] =_n*_n*(_n*(_n*((1695744-1118711*_n)*_n-1174656)+258048)+80640)/
      3870720;
    _h[2] =_n*_n*_n*(_n*(_n*(22276*_n-16929)-15984)+12852)/362880;
    _h[3] =_n*_n*_n*_n*((-830251*_n-158400)*_n+197865)/7257600;
    _h[4] =(453717-435388*_n)*_n*_n*_n*_n*_n/15966720;
    _h[5] =20648693*_n*_n*_n*_n*_n*_n/638668800;
    _hp[0]=_n*(_n*(_n*(_n*(_n*(31564*_n-66675)+34440)+47250)-100800)+75600)/
      151200;
    _hp[1]=_n*_n*(_n*(_n*((863232-1983433*_n)*_n+748608)-1161216)+524160)/
      1935360;
    _hp[2]=_n*_n*_n*(_n*(_n*(670412*_n+406647)-533952)+184464)/725760;
    _hp[3]=_n*_n*_n*_n*(_n*(6601661*_n-7732800)+2230245)/7257600;
    _hp[4]=(3438171-13675556*_n)*_n*_n*_n*_n*_n/7983360;
    _hp[5]=212378941*_n*_n*_n*_n*_n*_n/319334400;
#elif TM_TX_MAXPOW == 7
    _a1 =_a/(1+_n)*(_n*_n*(_n*_n*(_n*_n+4)+64)+256)/256;
    _h[0] =_n*(_n*(_n*(_n*(_n*((6156736-5406467*_n)*_n-6123600)-107520)+
		       14918400)-25804800)+19353600)/38707200;
    _h[1] =_n*_n*(_n*(_n*(_n*(_n*(829456*_n-5593555)+8478720)-5873280)+
		      1290240)+403200)/19353600;
    _h[2] =_n*_n*_n*(_n*(_n*(_n*(9261899*_n+3564160)-2708640)-2557440)+
		     2056320)/58060800;
    _h[3] =_n*_n*_n*_n*(_n*(_n*(14928352*_n-9132761)-1742400)+2176515)/
      79833600;
    _h[4] =_n*_n*_n*_n*_n*((-8005831*_n-1741552)*_n+1814868)/63866880;
    _h[5] =(268433009-261810608*_n)*_n*_n*_n*_n*_n*_n/8302694400.;
    _h[6] =219941297*_n*_n*_n*_n*_n*_n*_n/5535129600.;
    _hp[0]=_n*(_n*(_n*(_n*(_n*(_n*(1804025*_n+2020096)-4267200)+2204160)+
		       3024000)-6451200)+4838400)/9676800;
    _hp[1]=_n*_n*(_n*(_n*(_n*(_n*(4626384*_n-9917165)+4316160)+3743040)-
		      5806080)+2620800)/9676800;
    _hp[2]=_n*_n*_n*(_n*(_n*((26816480-67102379*_n)*_n+16265880)-21358080)+
		     7378560)/29030400;
    _hp[3]=_n*_n*_n*_n*(_n*(_n*(155912000*_n+72618271)-85060800)+24532695)/
      79833600;
    _hp[4]=_n*_n*_n*_n*_n*(_n*(102508609*_n-109404448)+27505368)/63866880;
    _hp[5]=(2760926233.-12282192400.*_n)*_n*_n*_n*_n*_n*_n/4151347200.;
    _hp[6]=1522256789.*_n*_n*_n*_n*_n*_n*_n/1383782400.;
#elif TM_TX_MAXPOW >= 8
    _a1 =_a/(1+_n)*(_n*_n*(_n*_n*(_n*_n*(25*_n*_n+64)+256)+4096)+16384)/
      16384;
    _h[0] =_n*(_n*(_n*(_n*(_n*(_n*(_n*(31777436*_n-37845269)+43097152)-
			       42865200)-752640)+104428800)-180633600)+
	       135475200)/270950400;
    _h[1] =_n*_n*(_n*(_n*(_n*(_n*(_n*(24749483*_n+14930208)-100683990)+
			      152616960)-105719040)+23224320)+7257600)/
      348364800;
    _h[2] =_n*_n*_n*(_n*(_n*(_n*((101880889-232468668*_n)*_n+39205760)-
			     29795040)-28131840)+22619520)/638668800;
    _h[3] =_n*_n*_n*_n*(_n*(_n*(_n*(324154477*_n+1433121792.)-876745056)-
			    167270400)+208945440)/7664025600.;
    _h[4] =_n*_n*_n*_n*_n*(_n*(_n*(457888660*_n-312227409)-67920528)+
			   70779852)/2490808320.;
    _h[5] =_n*_n*_n*_n*_n*_n*((-19841813847.*_n-3665348512.)*_n+3758062126.)/
      116237721600.;
    _h[6] =(1979471673.-1989295244.*_n)*_n*_n*_n*_n*_n*_n*_n/49816166400.;
    _h[7] =191773887257.*_n*_n*_n*_n*_n*_n*_n*_n/3719607091200.;
    _hp[0]=_n*(_n*(_n*(_n*(_n*(_n*((37884525-75900428*_n)*_n+42422016)-
			       89611200)+46287360)+63504000)-135475200)+
	       101606400)/203212800;
    _hp[1]=_n*_n*(_n*(_n*(_n*(_n*(_n*(148003883*_n+83274912)-178508970)+
			      77690880)+67374720)-104509440)+47174400)/
      174182400;
    _hp[2]=_n*_n*_n*(_n*(_n*(_n*(_n*(318729724*_n-738126169)+294981280)+
			     178924680)-234938880)+81164160)/319334400;
    _hp[3]=_n*_n*_n*_n*(_n*(_n*((14967552000.-40176129013.*_n)*_n+
				6971354016.)-8165836800.)+2355138720.)/
      7664025600.;
    _hp[4]=_n*_n*_n*_n*_n*(_n*(_n*(10421654396.*_n+3997835751.)-4266773472.)+
			   1072709352.)/2490808320.;
    _hp[5]=_n*_n*_n*_n*_n*_n*(_n*(175214326799.*_n-171950693600.)+
			      38652967262.)/58118860800.;
    _hp[6]=(13700311101.-67039739596.*_n)*_n*_n*_n*_n*_n*_n*_n/12454041600.;
    _hp[7]=1424729850961.*_n*_n*_n*_n*_n*_n*_n*_n/743921418240.;
#endif
#if TM_SC_MAXPOW <= 4
    _s[0] = (_e2*(_e2*((-1367*_e2-2320)*_e2-4608)-12288)+98304)/98304;
    _s[1] = _e2*(_e2*((184-43*_e2)*_e2+1152)+6144)/49152;
    _s[2] = _e2*_e2*(_e2*(665*_e2+1360)+2304)/98304;
    _s[3] = _e2*_e2*_e2*(611*_e2+592)/98304;
    _s[4] = 59*_e2*_e2*_e2*_e2/32768;
#elif TM_SC_MAXPOW == 5
    _s[0] = (_e2*(_e2*(_e2*((-35577*_e2-54680)*_e2-92800)-184320)-491520)+
	     3932160)/3932160;
    _s[1] = _e2*(_e2*(_e2*((-3603*_e2-1720)*_e2+7360)+46080)+245760)/1966080;
    _s[2] = _e2*_e2*(_e2*(_e2*(23547*_e2+53200)+108800)+184320)/7864320;
    _s[3] = _e2*_e2*_e2*(_e2*(36243*_e2+48880)+47360)/7864320;
    _s[4] = _e2*_e2*_e2*_e2*(7049*_e2+4720)/2621440;
    _s[5] = 1543*_e2*_e2*_e2*_e2*_e2/2621440;
#elif TM_SC_MAXPOW == 6
    _s[0] = (_e2*(_e2*(_e2*(_e2*((-594943*_e2-853848)*_e2-1312320)-2227200)-
		       4423680)-11796480)+94371840)/94371840;
    _s[1] = _e2*(_e2*(_e2*(_e2*((-114471*_e2-115296)*_e2-55040)+235520)+
		      1474560)+7864320)/62914560;
    _s[2] = _e2*_e2*(_e2*(_e2*(_e2*(71271*_e2+188376)+425600)+870400)+
		     1474560)/62914560;
    _s[3] = _e2*_e2*_e2*(_e2*(_e2*(374403*_e2+579888)+782080)+757760)/125829120;
    _s[4] = _e2*_e2*_e2*_e2*(_e2*(111325*_e2+112784)+75520)/41943040;
    _s[5] = _e2*_e2*_e2*_e2*_e2*(48533*_e2+24688)/41943040;
    _s[6] = 77041*_e2*_e2*_e2*_e2*_e2*_e2/377487360;
#elif TM_SC_MAXPOW == 7
    _s[0] = (_e2*(_e2*(_e2*(_e2*(_e2*((-195279827*_e2-266534464)*_e2-
				      382523904)-587919360)-997785600)-
		       1981808640.)-5284823040.)+42278584320.)/42278584320.;
    _s[1] = _e2*(_e2*(_e2*(_e2*(_e2*((-5538513*_e2-6410376)*_e2-6456576)-
				3082240)+13189120)+82575360)+
		 440401920)/3523215360.;
    _s[2] = _e2*_e2*(_e2*(_e2*(_e2*(_e2*(3546477*_e2+15964704)+42196224)+
			       95334400)+194969600)+330301440)/14092861440.;
    _s[3] = _e2*_e2*_e2*(_e2*(_e2*(_e2*(1550463*_e2+2620821)+4059216)+5474560)+
			 5304320)/880803840;
    _s[4] = _e2*_e2*_e2*_e2*(_e2*(_e2*(10243421*_e2+12468400)+12631808)+
			     8458240)/4697620480.;
    _s[5] = _e2*_e2*_e2*_e2*_e2*(_e2*(6696873*_e2+5435696)+2765056)/4697620480.;
    _s[6] = _e2*_e2*_e2*_e2*_e2*_e2*(21098327*_e2+8628592)/42278584320.;
    _s[7] = 69319*_e2*_e2*_e2*_e2*_e2*_e2*_e2/939524096;
#elif TM_SC_MAXPOW >= 8
    _s[0] = (_e2*(_e2*(_e2*(_e2*(_e2*(_e2*((-2378730735.*_e2-3124477232.)*_e2-
					   4264551424.)-6120382464.)-
				 9406709760.)-15964569600.)-31708938240.)-
		  84557168640.)+676457349120.)/676457349120.;
    _s[1] = _e2*(_e2*(_e2*(_e2*(_e2*(_e2*((-293610175*_e2-354464832)*_e2-
					  410264064)-413220864)-197263360)+
			   844103680)+5284823040.)+28185722880.)/225485783040.;
    _s[2] = _e2*_e2*(_e2*(_e2*(_e2*(_e2*((28371816-16236487*_e2)*_e2+
					 127717632)+337569792)+762675200)+
			  1559756800.)+2642411520.)/112742891520.;
    _s[3] = _e2*_e2*_e2*(_e2*(_e2*(_e2*(_e2*(108089577*_e2+198459264)+
					335465088)+519579648)+700743680)+
			 678952960)/112742891520.;
    _s[4] = _e2*_e2*_e2*_e2*(_e2*(_e2*(_e2*(22691857*_e2+30730263)+37405200)+
				  37895424)+25374720)/14092861440.;
    _s[5] = _e2*_e2*_e2*_e2*_e2*(_e2*(_e2*(318566285*_e2+321449904)+260913408)+
				 132722688)/225485783040.;
    _s[6] = _e2*_e2*_e2*_e2*_e2*_e2*(_e2*(498732939*_e2+337573232)+
				     138057472)/676457349120.;
    _s[7] = _e2*_e2*_e2*_e2*_e2*_e2*_e2*(48729311*_e2+16636560)/225485783040.;
    _s[8] = 6204619*_e2*_e2*_e2*_e2*_e2*_e2*_e2*_e2/225485783040.;
#endif
  }

  const TransverseMercator
  TransverseMercator::UTM(Constants::WGS84_a, Constants::WGS84_invf,
			  Constants::UTM_k0);

  void TransverseMercator::Scale(double phi, double lam, double xi, double eta,
				 double& gamma, double& k) const {
    // This returns approximations to the convergence and scale for the
    // transverse Mercator projection.
    double c = cos(phi);
    if (c < _tol) {
      gamma = lam;
      k = 1;
    } else {
      // Evaluate
      //
      //    S = cos(zeta) * sum(s[0] * cos(2*k*zeta), k = 0..scpow)
      //
      // via Clenshaw summation with
      //
      //    x = 2 * zeta
      //    F[n](x) = cos(n * x)
      //    a(n, x) = 2 * cos(x)
      //    b(n, x) = -1
      //    [ cos(A+B) - 2*cos(B)*cos(A) + cos(A-B) = 0, A = n*x, B = x ]
      //    S/cos(zeta) = (s[0] - y[2])  + y[1] * cos(x)
      double
	c0 = cos(2 * xi), ch0 = cosh(2 * eta),
	s0 = sin(2 * xi), sh0 = sinh(2 * eta),
	ar = 2 * c0 * ch0, ai = -2 * s0 * sh0; // 2 * cos(2*zeta)
      double
	yr0 = _s[scpow], yi0 = 0,
	yr1 = 0, yi1 = 0,
	yr2, yi2;
      for (int j = scpow; --j;) { // j = scpow-1 .. 1
	yr2 = yr1; yi2 = yi1;
	yr1 = yr0; yi1 = yi0;
	yr0 = ar * yr1 - ai * yi1 - yr2 + _s[j];
	yi0 = ai * yr1 + ar * yi1 - yi2;
      }
      ar /= 2; ai /= 2;		// cos(2*zeta)
      yr2 = _s[0] - yr1 + ar * yr0 - ai * yi0;
      yi2 =       - yi1 + ai * yr0 + ar * yi0;
      c0 = cos(xi), ch0 = cosh(eta);
      s0 = sin(xi), sh0 = sinh(eta);
      ar = c0 * ch0, ai = -s0 * sh0; // cos(zeta)
      double
	sr = ar * yr2 - ai * yi2,
	si = ai * yr2 + ar * yi2;
      gamma = -atan2(si, sr);
      k = sqrt(_e2m + _e2 * c * c) / c * hypot(sr, si);
    }
  }

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
      l = lon * Constants::degree;
    // q is isometric latitude
    // JHS 154 has
    //
    //   beta = atan(sinh(q)) = conformal latitude
    //   [xi', eta'] = spheroidal TM coordinates
    //   eta' = atanh(cos(beta) * sin(l))
    //   xi' = asin(sin(beta)*cosh(eta')
    //
    // We use
    //
    //   tan(beta) = sinh(q)
    //   sin(beta) = tanh(q)
    //   cos(beta) = sech(q)
    //   denom^2    = 1-cos(beta)^2*sin(l)^2 = 1-sech(q)^2*sin(l)^2
    //   sin(xip)   = sin(beta)/denom        = tanh(q)/denom
    //   cos(xip)   = cos(beta)*cos(l)/denom = sech(q)*cos(l)/denom
    //   cosh(etap) = 1/denom                = 1/denom
    //   sinh(etap) = cos(beta)*sin(l)/denom = sech(q)*sin(l)/denom
    //
    // to eliminate beta and derive more stable expressions for xi',eta'
    double etap, xip;
    if (lat < 90) {
      double
	qp = asinh(tan(phi)),
	qpp = atanh(_e * sin(phi)),
	q = qp - _e * qpp;
      etap = atanh(sin(l) / cosh(q));
      xip = atan2(sinh(q), cos(l));
    } else {
      xip = Constants::pi/2;
      etap = 0;
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
    double
      c0 = cos(2 * xip), ch0 = cosh(2 * etap),
      s0 = sin(2 * xip), sh0 = sinh(2 * etap),
      ar = 2 * c0 * ch0, ai = -2 * s0 * sh0; // 2 * cos(2*zeta')
    double
      xi0 = _hp[maxpow - 1], eta0 = 0,
      xi1 = 0, eta1 = 0,
      xi2, eta2;
    for (int j = maxpow - 1; j--;) {
      xi2 = xi1; eta2 = eta1;
      xi1 = xi0; eta1 = eta0;
      xi0  = ar * xi1 - ai * eta1 - xi2 + _hp[j];
      eta0 = ai * xi1 + ar * eta1 - eta2;
    }
    ar = s0 * ch0; ai = c0 * sh0; // sin(2*zeta')
    double
      xi  = xip  + ar * xi0 - ai * eta0,
      eta = etap + ai * xi0 + ar * eta0;
    Scale(phi, l, xi, eta, gamma, k);
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
    double
      xip0 = -_h[maxpow - 1], etap0 = 0,
      xip1 = 0, etap1 = 0,
      xip2, etap2;
    for (int j = maxpow - 1; j--;) {
      xip2 = xip1; etap2 = etap1;
      xip1 = xip0; etap1 = etap0;
      xip0  = ar * xip1 - ai * etap1 - xip2 - _h[j];
      etap0 = ai * xip1 + ar * etap1 - etap2;
    }
    ar = s0 * ch0; ai = c0 * sh0; // sin(2*zeta)
    double
      xip  = xi  + ar * xip0 - ai * etap0,
      etap = eta + ai * xip0 + ar * etap0;
    // JHS 154 has
    //
    // 	 beta = asin(sin(xip) / cosh(etap))
    // 	 l = asin(tanh(etap) / cos(beta)
    // 	 q = asinh(tan(beta))
    //
    // the following eliminates beta and is more stable

    double l, phi;
    double
      s = sinh(etap),
      c = cos(xip),
      r = hypot(s, c);
    if (r > 0) {
      l = atan2(s, c);
      // Solve
      // q = qp - e * atanh(e * tanh(qp))
      // for qp = asinh(tan(phi))
      double
	q = asinh(sin(xip)/r),
	qp = q;
      for (int i = 0; i < _numit; ++i) {
	double
	  t = tanh(qp),
	  dqp = -(qp - _e * atanh(_e * t) - q) *
	  (1 - _e2 * t * t) / _e2m;
	qp += dqp;
	if (std::abs(dqp) < _tol)
	  break;
      }
      phi = atan(sinh(qp));
    } else {
      phi = Constants::pi/2;
      l = 0;
    }
    lat = phi / Constants::degree * xisign;
    lon = l / Constants::degree;
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
    Scale(phi, l, xi, eta, gamma, k);
    gamma /= Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= xisign * etasign;
    k *= _k0;
  }

} // namespace GeographicLib

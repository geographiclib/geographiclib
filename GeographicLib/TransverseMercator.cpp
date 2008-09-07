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
 *    dxn  = max(error in forward transformation,
 *               discrepancy in forward and reverse transformations)
 *           for nth order method (order e^(2*n))
 *
 *    dgam = max error in meridian convergence using the O(ep2^2) formula
 *
 *    dk   = max relative error in scale using the O(ep2) formula
 *
 *      set         dx4     dx5     dx6     dx7     dx8   dgam   dk
 * x<4e5, y<95e5 .21um     4nm     4nm     4nm     4nm   3e-8"  4e-7%%
 * x<5e5, y<96e5 .22um     4nm     4nm     4nm     4nm   9e-8"  1e-6%%
 *     mu<10     .35um     5nm     4nm     4nm     4nm   5e-6"  3e-5%%
 *     mu<15     .70um     6nm     4nm     4nm     4nm   4e-5"  2e-4%%
 *     mu<20     1.5um    11nm     4nm     4nm     4nm   2e-4"  5e-4%%
 *     mu<25     3.6um    25nm     4nm     4nm     4nm   7e-4"  2e-3%%
 *     mu<30     8.7um    71nm     4nm     4nm     4nm   3e-3"  3e-3%%
 *     mu<35      22um   .22um     6nm     4nm     4nm   7e-3"  6e-3%%
 *     mu<40      61um   .75um    11nm     4nm     4nm   0.02"  0.02%%
 *     mu<45     .18mm   2.8um    49nm     5nm     5nm   0.05"  0.03%%
 *     mu<50     .62mm    12um   .27um    10nm     6nm   0.2"   0.05%%
 *     mu<55     2.4mm    64um   1.8um    58nm     7nm   0.4"   0.1%%
 *     mu<60      12mm   .42mm    17um   .69um    34nm   1"     0.3%%
 *     mu<65      72mm   3.8mm   .22mm    13um   .83um   4"     0.6%%
 *     mu<70     .67m     56mm   5.1mm   .48mm    48um   18"    1.5%%
 *     mu<72     1.9m    .20m     22mm   2.6mm   .32mm   35"    2.4%%
 *     mu<74     6.6m    .88m    .12m     18mm   2.8mm   76"    3.9%%
 *     mu<76      27m    4.7m    .87m    .17m     35mm   3'     7.1%%
 *     mu<78     .14km    33m    8.4m    2.2m    .62m    8'     1.4%
 *     mu<80     1.0km   .36km   .13km    52m     21m    28'    3.3%
 *     mu<82      14km   8.1km   4.8km   3.0km   2.0km   2.4d   10%
 *     mu<84     390km   380km   390km   410km   440km   14d    15%
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
#if TM_MAXPOW <= 4
    , _a1 (_a/(1+_n)*(_n*_n*(_n*_n+16)+64)/64)
    , _h1 (_n*(_n*((555-4*_n)*_n-960)+720)/1440)
    , _h1p(_n*(_n*(_n*(164*_n+225)-480)+360)/720)
    , _h2 (_n*_n*((96-437*_n)*_n+30)/1440)
    , _h2p(_n*_n*(_n*(557*_n-864)+390)/1440)
    , _h3 ((119-148*_n)*_n*_n*_n/3360)
    , _h3p((427-1236*_n)*_n*_n*_n/1680)
    , _h4 (4397*_n*_n*_n*_n/161280)
    , _h4p(49561*_n*_n*_n*_n/161280)
#elif TM_MAXPOW == 5
    , _a1 (_a/(1+_n)*(_n*_n*(_n*_n+16)+64)/64)
    , _h1 (_n*(_n*(_n*((-3645*_n-64)*_n+8880)-15360)+11520)/23040)
    , _h1p(_n*(_n*(_n*((328-635*_n)*_n+450)-960)+720)/1440)
    , _h2 (_n*_n*(_n*(_n*(4416*_n-3059)+672)+210)/10080)
    , _h2p(_n*_n*(_n*(_n*(4496*_n+3899)-6048)+2730)/10080)
    , _h3 (_n*_n*_n*((-627*_n-592)*_n+476)/13440)
    , _h3p(_n*_n*_n*(_n*(15061*_n-19776)+6832)/26880)
    , _h4 ((4397-3520*_n)*_n*_n*_n*_n/161280)
    , _h4p((49561-171840*_n)*_n*_n*_n*_n/161280)
    , _h5 (4583*_n*_n*_n*_n*_n/161280)
    , _h5p(34729*_n*_n*_n*_n*_n/80640)
#elif TM_MAXPOW == 6
    , _a1 (_a/(1+_n)*(_n*_n*(_n*_n*(_n*_n+4)+64)+256)/256)
    , _h1 (_n*(_n*(_n*(_n*(_n*(384796*_n-382725)-6720)+932400)-1612800)+1209600)/2419200)
    , _h1p(_n*(_n*(_n*(_n*(_n*(31564*_n-66675)+34440)+47250)-100800)+75600)/151200)
    , _h2 (_n*_n*(_n*(_n*((1695744-1118711*_n)*_n-1174656)+258048)+80640)/3870720)
    , _h2p(_n*_n*(_n*(_n*((863232-1983433*_n)*_n+748608)-1161216)+524160)/1935360)
    , _h3 (_n*_n*_n*(_n*(_n*(22276*_n-16929)-15984)+12852)/362880)
    , _h3p(_n*_n*_n*(_n*(_n*(670412*_n+406647)-533952)+184464)/725760)
    , _h4 (_n*_n*_n*_n*((-830251*_n-158400)*_n+197865)/7257600)
    , _h4p(_n*_n*_n*_n*(_n*(6601661*_n-7732800)+2230245)/7257600)
    , _h5 ((453717-435388*_n)*_n*_n*_n*_n*_n/15966720)
    , _h5p((3438171-13675556*_n)*_n*_n*_n*_n*_n/7983360)
    , _h6 (20648693*_n*_n*_n*_n*_n*_n/638668800)
    , _h6p(212378941*_n*_n*_n*_n*_n*_n/319334400)
#elif TM_MAXPOW == 7
    , _a1 (_a/(1+_n)*(_n*_n*(_n*_n*(_n*_n+4)+64)+256)/256)
    , _h1 (_n*(_n*(_n*(_n*(_n*((6156736-5406467*_n)*_n-6123600)-107520)+14918400)-25804800)+19353600)/38707200)
    , _h1p(_n*(_n*(_n*(_n*(_n*(_n*(1804025*_n+2020096)-4267200)+2204160)+3024000)-6451200)+4838400)/9676800)
    , _h2 (_n*_n*(_n*(_n*(_n*(_n*(829456*_n-5593555)+8478720)-5873280)+1290240)+403200)/19353600)
    , _h2p(_n*_n*(_n*(_n*(_n*(_n*(4626384*_n-9917165)+4316160)+3743040)-5806080)+2620800)/9676800)
    , _h3 (_n*_n*_n*(_n*(_n*(_n*(9261899*_n+3564160)-2708640)-2557440)+2056320)/58060800)
    , _h3p(_n*_n*_n*(_n*(_n*((26816480-67102379*_n)*_n+16265880)-21358080)+7378560)/29030400)
    , _h4 (_n*_n*_n*_n*(_n*(_n*(14928352*_n-9132761)-1742400)+2176515)/79833600)
    , _h4p(_n*_n*_n*_n*(_n*(_n*(155912000*_n+72618271)-85060800)+24532695)/79833600)
    , _h5 (_n*_n*_n*_n*_n*((-8005831*_n-1741552)*_n+1814868)/63866880)
    , _h5p(_n*_n*_n*_n*_n*(_n*(102508609*_n-109404448)+27505368)/63866880)
    , _h6 ((268433009-261810608*_n)*_n*_n*_n*_n*_n*_n/8302694400.)
    , _h6p((2760926233.-12282192400.*_n)*_n*_n*_n*_n*_n*_n/4151347200.)
    , _h7 (219941297*_n*_n*_n*_n*_n*_n*_n/5535129600.)
    , _h7p(1522256789.*_n*_n*_n*_n*_n*_n*_n/1383782400.)
#elif TM_MAXPOW >= 8
    , _a1 (_a/(1+_n)*(_n*_n*(_n*_n*(_n*_n*(25*_n*_n+64)+256)+4096)+16384)/16384)
    , _h1 (_n*(_n*(_n*(_n*(_n*(_n*(_n*(31777436*_n-37845269)+43097152)-42865200)-752640)+104428800)-180633600)+135475200)/270950400)
    , _h1p(_n*(_n*(_n*(_n*(_n*(_n*((37884525-75900428*_n)*_n+42422016)-89611200)+46287360)+63504000)-135475200)+101606400)/203212800)
    , _h2 (_n*_n*(_n*(_n*(_n*(_n*(_n*(24749483*_n+14930208)-100683990)+152616960)-105719040)+23224320)+7257600)/348364800)
    , _h2p(_n*_n*(_n*(_n*(_n*(_n*(_n*(148003883*_n+83274912)-178508970)+77690880)+67374720)-104509440)+47174400)/174182400)
    , _h3 (_n*_n*_n*(_n*(_n*(_n*((101880889-232468668*_n)*_n+39205760)-29795040)-28131840)+22619520)/638668800)
    , _h3p(_n*_n*_n*(_n*(_n*(_n*(_n*(318729724*_n-738126169)+294981280)+178924680)-234938880)+81164160)/319334400)
    , _h4 (_n*_n*_n*_n*(_n*(_n*(_n*(324154477*_n+1433121792.)-876745056)-167270400)+208945440)/7664025600.)
    , _h4p(_n*_n*_n*_n*(_n*(_n*((14967552000.-40176129013.*_n)*_n+6971354016.)-8165836800.)+2355138720.)/7664025600.)
    , _h5 (_n*_n*_n*_n*_n*(_n*(_n*(457888660*_n-312227409)-67920528)+70779852)/2490808320.)
    , _h5p(_n*_n*_n*_n*_n*(_n*(_n*(10421654396.*_n+3997835751.)-4266773472.)+1072709352.)/2490808320.)
    , _h6 (_n*_n*_n*_n*_n*_n*((-19841813847.*_n-3665348512.)*_n+3758062126.)/116237721600.)
    , _h6p(_n*_n*_n*_n*_n*_n*(_n*(175214326799.*_n-171950693600.)+38652967262.)/58118860800.)
    , _h7 ((1979471673.-1989295244.*_n)*_n*_n*_n*_n*_n*_n*_n/49816166400.)
    , _h7p((13700311101.-67039739596.*_n)*_n*_n*_n*_n*_n*_n*_n/12454041600.)
    , _h8 (191773887257.*_n*_n*_n*_n*_n*_n*_n*_n/3719607091200.)
    , _h8p(1424729850961.*_n*_n*_n*_n*_n*_n*_n*_n/743921418240.)
#endif
    , _tol(0.1*sqrt(std::numeric_limits<double>::epsilon()))
    , _numit(5)
  {}

  const TransverseMercator
  TransverseMercator::UTM(Constants::WGS84_a, Constants::WGS84_invf,
			  Constants::UTM_k0);

  void TransverseMercator::Scale(double phi, double l,
				 double& gamma, double& k) const {
    // This returns approximations to the convergence and scale for the EXACT
    // transverse Mercator projection.  Thus this routine returns a constant
    // scale on the central meridian even though the approximate transverse
    // Mercator projection has a (slightly) varying scale.
    double
      c = cos(phi),
      s = sin(l),
      c2 = c * c,
      s2 = s * s,
      d = 1 - s2 * c2,
      // Accurate to order _ep2^2
      carg = 1 + c2 * c2 * s2 / d * _ep2 *
      (1 + c2 / (3 * d * d) *
       (2 - s2 * (c2 * ((4 * c2 - 1) * s2 - 9) + 8)) * _ep2),
      // Accurate to order _ep2
      cabs = 1 + c2 * c2 * s2 * ((c2 - 2) * s2 + 1) / (2 * d * d) * _ep2;
#if 0
    // To order _ep2^4
    carg = 1+c2^2*s2*_ep2/d-
      ((4*c2^5-c2^4)*s2^3+(-9*c2^4+8*c2^3)*s2^2-2*c2^3*s2)*_ep2^2/(3*d^3)+
      ((28*c2^8-11*c2^7-2*c2^6)*s2^5+(-120*c2^7+75*c2^6+48*c2^5)*s2^4+
       (165*c2^6-336*c2^5+120*c2^4)*s2^3+
       (93*c2^5-60*c2^4)*s2^2)*_ep2^3/(15*d^5)-
      ((836*c2^11-429*c2^10-75*c2^9-17*c2^8)*s2^7+
       (-5264*c2^10+3472*c2^9+1106*c2^8+1264*c2^7)*s2^6+
       (13097*c2^9-13566*c2^8-14460*c2^7+11536*c2^6)*s2^5+
       (-13727*c2^8+52288*c2^7-43064*c2^6+8064*c2^5)*s2^4+
       (-15054*c2^7+20272*c2^6-6048*c2^5)*s2^3+
       (-735*c2^6+504*c2^5)*s2^2)*_ep2^4/(315*d^7);
    // To order _ep2^4
    cabs = 1+((c2^3-2*c2^2)*s2^2+c2^2*s2)*_ep2/(2*d^2)-
      ((3*c2^6+12*c2^5-28*c2^4)*s2^4+(-34*c2^5+124*c2^4-64*c2^3)*s2^3+
       (-61*c2^4+48*c2^3)*s2^2)*_ep2^2/(24*d^4)+
      ((15*c2^9+30*c2^8+140*c2^7-328*c2^6)*s2^6+
       (-147*c2^8-1044*c2^7+4268*c2^6-2688*c2^5)*s2^5+
       (1405*c2^7-8274*c2^6+8480*c2^5-1920*c2^4)*s2^4+
       (3383*c2^6-5280*c2^5+1920*c2^4)*s2^3+
       (280*c2^5-240*c2^4)*s2^2)*_ep2^3/(240*d^6)-
      ((1575*c2^12+2520*c2^11+5880*c2^10+27552*c2^9-64016*c2^8)*s2^8+
       (-17396*c2^11-54680*c2^10-457680*c2^9+1860384*c2^8-1247744*c2^7)*s2^7+
       (114506*c2^10+1432424*c2^9-8037000*c2^8+9297792*c2^7-2874368*c2^6)*s2^6+
       (-1212372*c2^9+10085240*c2^8-16270464*c2^7+
	8397312*c2^6-1032192*c2^5)*s2^5+
       (-3418969*c2^8+8422624*c2^7-6228096*c2^6+1290240*c2^5)*s2^4+
       (-685664*c2^7+985152*c2^6-322560*c2^5)*s2^3)*_ep2^4/(40320*d^8);
#endif
    gamma = atan2(sin(phi) * s * carg, cos(l));
    k = cabs/sqrt(d);
    // Old expressions from JHS 154.  These are Taylor expansions in l as well
    // as ep2.
    // double v2 = 1 + _ep2 * c2;
    // Accurate to _ep2^2 and l^3
    // gamma = l * sin(phi) * (1 + v2 * (2 * v2 - 1) * c2 * l * l / 3)
    // Accurate to _ep2^0 and l^2
    // k = ( 1 + c2 * l * l / 2);
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
    //   denom^2 = 1-cos(beta)^2*sin(l)^2   = 1-sech(q)^2*sin(l)^2)
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
      if (lon < 90)
	xip = atan2(sinh(q), cos(l));
      else
	xip = Constants::pi/2;
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
    //   zeta = zeta' + sum(h[j]' * sin(2 * j * zeta'), j = 1..4)
    //
    // where h[j]' = O(_n^j).  The reversion of this series gives
    //
    //   zeta' = zeta - sum(h[j] * sin(2 * j * zeta), j = 1..4)
    //
    // which is used in Reverse.  Extensions of these series to order _n^8 are
    // given above.
    double xi = 0, eta = 0;
#if TM_MAXPOW >= 8
    xi  += _h8p * sin(16 * xip) * cosh(16 * etap);
    eta += _h8p * cos(16 * xip) * sinh(16 * etap);
#endif
#if TM_MAXPOW >= 7
    xi  += _h7p * sin(14 * xip) * cosh(14 * etap);
    eta += _h7p * cos(14 * xip) * sinh(14 * etap);
#endif
#if TM_MAXPOW >= 6
    xi  += _h6p * sin(12 * xip) * cosh(12 * etap);
    eta += _h6p * cos(12 * xip) * sinh(12 * etap);
#endif
#if TM_MAXPOW >= 5
    xi  += _h5p * sin(10 * xip) * cosh(10 * etap);
    eta += _h5p * cos(10 * xip) * sinh(10 * etap);
#endif
    xi  += _h4p * sin(8 * xip) * cosh(8 * etap);
    eta += _h4p * cos(8 * xip) * sinh(8 * etap);
    xi  += _h3p * sin(6 * xip) * cosh(6 * etap);
    eta += _h3p * cos(6 * xip) * sinh(6 * etap);
    xi  += _h2p * sin(4 * xip) * cosh(4 * etap);
    eta += _h2p * cos(4 * xip) * sinh(4 * etap);
    xi  += _h1p * sin(2 * xip) * cosh(2 * etap);
    eta += _h1p * cos(2 * xip) * sinh(2 * etap);
    xi  += xip;
    eta += etap;
    y = _a1 * _k0 * (backside ? Constants::pi - xi : xi) * latsign;
    x = _a1 * _k0 * eta * lonsign;
    Scale(phi, l, gamma, k);
    gamma /= Constants::degree;
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
    if (backside) {
      xi = Constants::pi - xi;
    }
    double xip = 0, etap = 0;
#if TM_MAXPOW >= 8
    xip  -= _h8 * sin(16 * xi) * cosh(16 * eta);
    etap -= _h8 * cos(16 * xi) * sinh(16 * eta);
#endif
#if TM_MAXPOW >= 7
    xip  -= _h7 * sin(14 * xi) * cosh(14 * eta);
    etap -= _h7 * cos(14 * xi) * sinh(14 * eta);
#endif
#if TM_MAXPOW >= 6
    xip  -= _h6 * sin(12 * xi) * cosh(12 * eta);
    etap -= _h6 * cos(12 * xi) * sinh(12 * eta);
#endif
#if TM_MAXPOW >= 5
    xip  -= _h5 * sin(10 * xi) * cosh(10 * eta);
    etap -= _h5 * cos(10 * xi) * sinh(10 * eta);
#endif
    xip  -= _h4 * sin(8 * xi) * cosh(8 * eta);
    etap -= _h4 * cos(8 * xi) * sinh(8 * eta);
    xip  -= _h3 * sin(6 * xi) * cosh(6 * eta);
    etap -= _h3 * cos(6 * xi) * sinh(6 * eta);
    xip  -= _h2 * sin(4 * xi) * cosh(4 * eta);
    etap -= _h2 * cos(4 * xi) * sinh(4 * eta);
    xip  -= _h1 * sin(2 * xi) * cosh(2 * eta);
    etap -= _h1 * cos(2 * xi) * sinh(2 * eta);
    xip  += xi;
    etap += eta;
    // JHS has
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
    Scale(phi, l, gamma, k);
    gamma /= Constants::degree;
    if (backside)
      gamma = 180 - gamma;
    gamma *= xisign * etasign;
    k *= _k0;
  }

} // namespace GeographicLib

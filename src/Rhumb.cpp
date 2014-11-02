/**
 * \file Rhumb.cpp
 * \brief Implementation for GeographicLib::Rhumb and GeographicLib::RhumbLine
 * classes
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <algorithm>
#include <GeographicLib/Rhumb.hpp>

namespace GeographicLib {

  using namespace std;

  Rhumb::Rhumb(real a, real f, bool exact)
    : _ell(a, f)
    , _exact(exact)
    , _c2(_ell.Area() / 720)
  {
    real n = _ell._n, nx = n;
    switch (maxpow_) {
    case 4:
      _R[1] = nx*(n*(n*(1772*n-5340)+6930)-4725)/14175;
      nx *= n;
      _R[2] = nx*((1590-1747*n)*n-630)/4725;
      nx *= n;
      _R[3] = nx*(104*n-31)/315;
      nx *= n;
      _R[4] = -41*nx/420;
      break;
    case 5:
      _R[1] = nx*(n*(n*(n*(41662*n+58476)-176220)+228690)-155925)/467775;
      nx *= n;
      _R[2] = nx*(n*(n*(18118*n-57651)+52470)-20790)/155925;
      nx *= n;
      _R[3] = nx*((17160-23011*n)*n-5115)/51975;
      nx *= n;
      _R[4] = nx*(5480*n-1353)/13860;
      nx *= n;
      _R[5] = -668*nx/5775;
      break;
    case 6:
      _R[1] = nx*(n*(n*(n*((56868630-114456994*n)*n+79819740)-240540300)+
                     312161850)-212837625)/638512875;
      nx *= n;
      _R[2] = nx*(n*(n*(n*(51304574*n+24731070)-78693615)+71621550)-28378350)/
        212837625;
      nx *= n;
      _R[3] = nx*(n*(n*(1554472*n-6282003)+4684680)-1396395)/14189175;
      nx *= n;
      _R[4] = nx*((3205800-4913956*n)*n-791505)/8108100;
      nx *= n;
      _R[5] = nx*(1092376*n-234468)/2027025;
      nx *= n;
      _R[6] = -313076*nx/2027025;
      break;
    case 7:
      _R[1] = nx*(n*(n*(n*(n*(n*(258618446*n-343370982)+170605890)+239459220)-
                        721620900)+936485550)-638512875)/1915538625;
      nx *= n;
      _R[2] = nx*(n*(n*(n*((153913722-248174686*n)*n+74193210)-236080845)+
                     214864650)-85135050)/638512875;
      nx *= n;
      _R[3] = nx*(n*(n*(n*(114450437*n+23317080)-94230045)+70270200)-20945925)/
        212837625;
      nx *= n;
      _R[4] = nx*(n*(n*(15445736*n-103193076)+67321800)-16621605)/170270100;
      nx *= n;
      _R[5] = nx*((16385640-27766753*n)*n-3517020)/30405375;
      nx *= n;
      _R[6] = nx*(4892722*n-939228)/6081075;
      nx *= n;
      _R[7] = -3189007*nx/14189175;
      break;
    case 8:
      _R[1] = nx*(n*(n*(n*(n*(n*((65947703730LL-13691187484LL*n)*n-
                                 87559600410LL)+43504501950LL)+61062101100LL)-
                        184013329500LL)+238803815250LL)-162820783125LL)/
        488462349375LL;
      nx *= n;
      _R[2] = nx*(n*(n*(n*(n*(n*(30802104839LL*n-63284544930LL)+39247999110LL)+
                           18919268550LL)-60200615475LL)+54790485750LL)-
                  21709437750LL)/162820783125LL;
      nx *= n;
      _R[3] = nx*(n*(n*(n*((5836972287LL-8934064508LL*n)*n+1189171080)-
                        4805732295LL)+3583780200LL)-1068242175)/10854718875LL;
      nx *= n;
      _R[4] = nx*(n*(n*(n*(50072287748LL*n+3938662680LL)-26314234380LL)+
                     17167059000LL)-4238509275LL)/43418875500LL;
      nx *= n;
      _R[5] = nx*(n*(n*(359094172*n-9912730821LL)+5849673480LL)-1255576140)/
        10854718875LL;
      nx *= n;
      _R[6] = nx*((8733508770LL-16053944387LL*n)*n-1676521980)/10854718875LL;
      nx *= n;
      _R[7] = nx*(930092876*n-162639357)/723647925;
      nx *= n;
      _R[8] = -673429061*nx/1929727800;
      break;
    default:
      GEOGRAPHICLIB_STATIC_ASSERT(maxpow_ >= 4 && maxpow_ <= 8,
                                  "Bad value of maxpow_");
    }
  }

  const Rhumb& Rhumb::WGS84() {
    static const Rhumb wgs84(Constants::WGS84_a(), Constants::WGS84_f(), false);
    return wgs84;
  }

  void Rhumb::GenInverse(real lat1, real lon1, real lat2, real lon2, bool areap,
                         real& s12, real& azi12, real& S12) const {
    real
      lon12 = Math::AngDiff(Math::AngNormalize(lon1), Math::AngNormalize(lon2)),
      psi1 = _ell.IsometricLatitude(lat1),
      psi2 = _ell.IsometricLatitude(lat2),
      psi12 = psi2 - psi1,
      h = Math::hypot(lon12, psi12);
    azi12 = 0 - atan2(-lon12, psi12) / Math::degree();
    psi1 *= Math::degree();
    psi2 *= Math::degree();
    real dmudpsi = DIsometricToRectifying(psi2, psi1);
    s12 = h * dmudpsi * _ell.QuarterMeridian() / 90;
    if (areap)
      S12 = _c2 * lon12 * MeanSinXi(psi2, psi1);
  }

  RhumbLine Rhumb::Line(real lat1, real lon1, real azi12) const
  { return RhumbLine(*this, lat1, lon1, azi12, _exact); }

  void Rhumb::GenDirect(real lat1, real lon1, real azi12, real s12, bool areap,
                        real& lat2, real& lon2, real& S12) const
  { Line(lat1, lon1, azi12).GenPosition(s12, areap, lat2, lon2, S12); }

  Math::real Rhumb::DE(real x, real y) const {
    const EllipticFunction& ei = _ell._ell;
    real d = x - y;
    if (x * y <= 0)
      return d ? (ei.E(x) - ei.E(y)) / d : 1;
    // See DLMF: Eqs (19.11.2) and (19.11.4) letting
    // theta -> x, phi -> -y, psi -> z
    //
    // (E(x) - E(y)) / d = E(z)/d - k2 * sin(x) * sin(y) * sin(z)/d
    //
    // tan(z/2) = (sin(x)*Delta(y) - sin(y)*Delta(x)) / (cos(x) + cos(y))
    //          = d * Dsin(x,y) * (sin(x) + sin(y))/(cos(x) + cos(y)) /
    //             (sin(x)*Delta(y) + sin(y)*Delta(x))
    //          = t = d * Dt
    // sin(z) = 2*t/(1+t^2); cos(z) = (1-t^2)/(1+t^2)
    // Alt (this only works for |z| <= pi/2 -- however, this conditions holds
    // if x*y > 0):
    // sin(z) = d * Dsin(x,y) * (sin(x) + sin(y))/
    //          (sin(x)*cos(y)*Delta(y) + sin(y)*cos(x)*Delta(x))
    // cos(z) = sqrt((1-sin(z))*(1+sin(z)))
    real sx = sin(x), sy = sin(y), cx = cos(x), cy = cos(y);
    real Dt = Dsin(x, y) * (sx + sy) /
      ((cx + cy) * (sx * ei.Delta(sy, cy) + sy * ei.Delta(sx, cx))),
      t = d * Dt, Dsz = 2 * Dt / (1 + t*t),
      sz = d * Dsz, cz = (1 - t) * (1 + t) / (1 + t*t);
    return ((sz ? ei.E(sz, cz, ei.Delta(sz, cz)) / sz : 1)
            - ei.k2() * sx * sy) * Dsz;
  }

  Math::real Rhumb::DRectifying(real latx, real laty) const {
    real
      phix = latx * Math::degree(), tbetx = _ell._f1 * tano(phix),
      phiy = laty * Math::degree(), tbety = _ell._f1 * tano(phiy);
    return (Math::pi()/2) * _ell._b * _ell._f1 * DE(atan(tbetx), atan(tbety))
      * Dtan(phix, phiy) * Datan(tbetx, tbety) / _ell.QuarterMeridian();
  }

  Math::real Rhumb::DIsometric(real latx, real laty) const {
    real
      phix = latx * Math::degree(), tx = tano(phix),
      phiy = laty * Math::degree(), ty = tano(phiy);
    return Dasinh(tx, ty) * Dtan(phix, phiy)
      - Deatanhe(sin(phix), sin(phiy)) * Dsin(phix, phiy);
  }

  Math::real Rhumb::SinCosSeries(bool sinp,
                                 real x, real y, const real c[], int n) {
    // N.B. n >= 0 and c[] has n+1 elements 0..n, of which c[0] is ignored.
    //
    // Use Clenshaw summation to evaluate
    //   m = (g(x) + g(y)) / 2         -- mean value
    //   s = (g(x) - g(y)) / (x - y)   -- average slope
    // where
    //   g(x) = sum(c[j]*SC(2*j*x), j = 1..n)
    //   SC = sinp ? sin : cos
    //   CS = sinp ? cos : sin
    //
    // This function returns only s; m is discarded.
    //
    // Write
    //   t = [m; s]
    //   t = sum(c[j] * f[j](x,y), j = 1..n)
    // where
    //   f[j](x,y) = [ (SC(2*j*x)+SC(2*j*y))/2 ]
    //               [ (SC(2*j*x)-SC(2*j*y))/d ]
    //
    //             = [       cos(j*d)*SC(j*p)    ]
    //               [ +/-(2/d)*sin(j*d)*CS(j*p) ]
    // (+/- = sinp ? + : -) and
    //    p = x+y, d = x-y
    //
    //   f[j+1](x,y) = A * f[j](x,y) - f[j-1](x,y)
    //
    //   A = [  2*cos(p)*cos(d)      -sin(p)*sin(d)*d]
    //       [ -4*sin(p)*sin(d)/d   2*cos(p)*cos(d)  ]
    //
    // Let b[n+1] = b[n+2] = [0 0; 0 0]
    //     b[j] = A * b[j+1] - b[j+2] + c[j] * I for j = n..1
    //    t =  (c[0] * I  - b[2]) * f[0](x,y) + b[1] * f[1](x,y)
    // c[0] is not accessed for s = t[2]
    real p = x + y, d = x - y,
      cp = cos(p), cd =     cos(d),
      sp = sin(p), sd = d ? sin(d)/d : 1,
      m = 2 * cp * cd, s = sp * sd;
    // 2x2 matrices stored in row-major order
    const real a[4] = {m, -s * d * d, -4 * s, m};
    real ba[4] = {0, 0, 0, 0};
    real bb[4] = {0, 0, 0, 0};
    real* b1 = ba;
    real* b2 = bb;
    if (n > 0) b1[0] = b1[3] = c[n];
    for (int j = n - 1; j > 0; --j) { // j = n-1 .. 1
      std::swap(b1, b2);
      // b1 = A * b2 - b1 + c[j] * I
      b1[0] = a[0] * b2[0] + a[1] * b2[2] - b1[0] + c[j];
      b1[1] = a[0] * b2[1] + a[1] * b2[3] - b1[1];
      b1[2] = a[2] * b2[0] + a[3] * b2[2] - b1[2];
      b1[3] = a[2] * b2[1] + a[3] * b2[3] - b1[3] + c[j];
    }
    // Here are the full expressions for m and s
    // m =   (c[0] - b2[0]) * f01 - b2[1] * f02 + b1[0] * f11 + b1[1] * f12;
    // s = - b2[2] * f01 + (c[0] - b2[3]) * f02 + b1[2] * f11 + b1[3] * f12;
    if (sinp) {
      // real f01 = 0, f02 = 0;
      real f11 = cd * sp, f12 = 2 * sd * cp;
      // m = b1[0] * f11 + b1[1] * f12;
      s = b1[2] * f11 + b1[3] * f12;
    } else {
      // real f01 = 1, f02 = 0;
      real f11 = cd * cp, f12 = - 2 * sd * sp;
      // m = c[0] - b2[0] + b1[0] * f11 + b1[1] * f12;
      s = - b2[2] + b1[2] * f11 + b1[3] * f12;
    }
    return s;
  }

  Math::real Rhumb::DConformalToRectifying(real chix, real chiy) const {
    return 1 + SinCosSeries(true, chix, chiy,
                            _ell.ConformalToRectifyingCoeffs(), tm_maxord);
  }

  Math::real Rhumb::DRectifyingToConformal(real mux, real muy) const {
    return 1 - SinCosSeries(true, mux, muy,
                            _ell.RectifyingToConformalCoeffs(), tm_maxord);
  }

  Math::real Rhumb::DIsometricToRectifying(real psix, real psiy) const {
    if (_exact) {
      real
        latx = _ell.InverseIsometricLatitude(psix/Math::degree()),
        laty = _ell.InverseIsometricLatitude(psiy/Math::degree());
      return DRectifying(latx, laty) / DIsometric(latx, laty);
    } else
      return DConformalToRectifying(gd(psix), gd(psiy)) * Dgd(psix, psiy);
  }

  Math::real Rhumb::DRectifyingToIsometric(real mux, real muy) const {
    real
      latx = _ell.InverseRectifyingLatitude(mux/Math::degree()),
      laty = _ell.InverseRectifyingLatitude(muy/Math::degree());
    return _exact ?
      DIsometric(latx, laty) / DRectifying(latx, laty) :
      Dgdinv(_ell.ConformalLatitude(latx) * Math::degree(),
             _ell.ConformalLatitude(laty) * Math::degree()) *
      DRectifyingToConformal(mux, muy);
  }

  Math::real Rhumb::MeanSinXi(real psix, real psiy) const {
    return Dlog(cosh(psix), cosh(psiy)) * Dcosh(psix, psiy)
      + SinCosSeries(false, gd(psix), gd(psiy), _R, maxpow_) * Dgd(psix, psiy);
    /* Checking code ...
    real a = Dlog(cosh(psix), cosh(psiy)) * Dcosh(psix, psiy);
    real b = SinCosSeries(false, gd(psix), gd(psiy), _R, maxpow_);
    real c = Dgd(psix, psiy);

    real chix = gd(psix), chiy = gd(psiy);
    real sx = 0, sy = 0;
    for (int i = 1; i <= maxpow_; ++i) {
      sx += _R[i] * cos(2*i*chix);
      sy += _R[i] * cos(2*i*chiy);
    }
    real aa = -(log(cos(chix)) - log(cos(chiy))) / (psix - psiy);
    real bb = (sx - sy) / (chix - chiy);
    real cc = (chix - chiy) / (psix - psiy);
    cout << aa << " " << a << " " << (a-aa)/aa << "\n";
    cout << bb << " " << b << " " << (b-bb)/bb << "\n";
    cout << cc << " " << c << " " << (c-cc)/cc << "\n";
    sx -= log(cos(chix));
    sy -= log(cos(chiy));
    return (sx - sy) / (psix - psiy);
    */
  }

  RhumbLine::RhumbLine(const Rhumb& rh, real lat1, real lon1, real azi12,
                       bool exact)
    : _rh(rh)
    , _exact(exact)
    , _lat1(lat1)
    , _lon1(Math::AngNormalize(lon1))
    , _azi12(Math::AngNormalize(azi12))
  {
    real alp12 = azi12 * Math::degree();
    _salp =     azi12  == -180 ? 0 : sin(alp12);
    _calp = abs(azi12) ==   90 ? 0 : cos(alp12);
    _mu1 = _rh._ell.RectifyingLatitude(lat1);
    _psi1 = _rh._ell.IsometricLatitude(lat1);
    _r1 = _rh._ell.CircleRadius(lat1);
  }

  void RhumbLine::GenPosition(real s12, bool areap,
                              real& lat2, real& lon2, real& S12) const {
    real
      mu12 = s12 * _calp * 90 / _rh._ell.QuarterMeridian(),
      mu2 = _mu1 + mu12;
    real psi2;
    if (abs(mu2) <= 90) {
      if (_calp) {
        lat2 = _rh._ell.InverseRectifyingLatitude(mu2);
        real psi12 = _rh.DRectifyingToIsometric(  mu2 * Math::degree(),
                                                 _mu1 * Math::degree()) * mu12;
        lon2 = _salp * psi12 / _calp;
        psi2 = _psi1 + psi12;
      } else {
        lat2 = _lat1;
        lon2 = _salp * s12 / (_r1 * Math::degree());
        psi2 = _psi1;
      }
      if (areap)
        S12 = _rh._c2 * lon2 *
          _rh.MeanSinXi(_psi1 * Math::degree(), psi2 * Math::degree());
      lon2 = Math::AngNormalize2(_lon1 + lon2);
    } else {
      // Reduce to the interval [-180, 180)
      mu2 = Math::AngNormalize2(mu2);
      // Deal with points on the anti-meridian
      if (abs(mu2) > 90) mu2 = Math::AngNormalize(180 - mu2);
      lat2 = _rh._ell.InverseRectifyingLatitude(mu2);
      lon2 = Math::NaN();
      if (areap)
        S12 = Math::NaN();
    }
  }

} // namespace GeographicLib

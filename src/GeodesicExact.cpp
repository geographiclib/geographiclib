/**
 * \file GeodesicExact.cpp
 * \brief Implementation for GeographicLib::GeodesicExact class
 *
 * Copyright (c) Charles Karney (2012-2022) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * This is a reformulation of the geodesic problem.  The notation is as
 * follows:
 * - at a general point (no suffix or 1 or 2 as suffix)
 *   - phi = latitude
 *   - beta = latitude on auxiliary sphere
 *   - omega = longitude on auxiliary sphere
 *   - lambda = longitude
 *   - alpha = azimuth of great circle
 *   - sigma = arc length along great circle
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

#include <GeographicLib/GeodesicExact.hpp>
#include <GeographicLib/GeodesicLineExact.hpp>

#if defined(_MSC_VER)
// Squelch warnings about potentially uninitialized local variables,
// constant conditional and enum-float expressions and mixing enums
#  pragma warning (disable: 4701 4127 5055 5054)
#endif

namespace GeographicLib {

  using namespace std;

  GeodesicExact::GeodesicExact(real a, real f)
    : maxit2_(maxit1_ + Math::digits() + 10)
      // Underflow guard.  We require
      //   tiny_ * epsilon() > 0
      //   tiny_ + epsilon() == epsilon()
    , tiny_(sqrt(numeric_limits<real>::min()))
    , tol0_(numeric_limits<real>::epsilon())
      // Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse
      // case 52.784459512564 0 -52.784459512563990912 179.634407464943777557
      // which otherwise failed for Visual Studio 10 (Release and Debug)
    , tol1_(200 * tol0_)
    , tol2_(sqrt(tol0_))
    , tolb_(tol0_ * tol2_)      // Check on bisection interval
    , xthresh_(1000 * tol2_)
    , _a(a)
    , _f(f)
    , _f1(1 - _f)
    , _e2(_f * (2 - _f))
    , _ep2(_e2 / Math::sq(_f1)) // e2 / (1 - e2)
    , _n(_f / ( 2 - _f))
    , _b(_a * _f1)
      // The Geodesic class substitutes atanh(sqrt(e2)) for asinh(sqrt(ep2)) in
      // the definition of _c2.  The latter is more accurate for very oblate
      // ellipsoids (which the Geodesic class does not attempt to handle).  Of
      // course, the area calculation in GeodesicExact is still based on a
      // series and so only holds for moderately oblate (or prolate)
      // ellipsoids.
    , _c2((Math::sq(_a) + Math::sq(_b) *
           (_f == 0 ? 1 :
            (_f > 0 ? asinh(sqrt(_ep2)) : atan(sqrt(-_e2))) /
            sqrt(fabs(_e2))))/2) // authalic radius squared
      // The sig12 threshold for "really short".  Using the auxiliary sphere
      // solution with dnm computed at (bet1 + bet2) / 2, the relative error in
      // the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
      // (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a
      // given f and sig12, the max error occurs for lines near the pole.  If
      // the old rule for computing dnm = (dn1 + dn2)/2 is used, then the error
      // increases by a factor of 2.)  Setting this equal to epsilon gives
      // sig12 = etol2.  Here 0.1 is a safety factor (error decreased by 100)
      // and max(0.001, abs(f)) stops etol2 getting too large in the nearly
      // spherical case.
    , _etol2(real(0.1) * tol2_ /
             sqrt( fmax(real(0.001), fabs(_f)) * fmin(real(1), 1 - _f/2) / 2 ))
  {
    if (!(isfinite(_a) && _a > 0))
      throw GeographicErr("Equatorial radius is not positive");
    if (!(isfinite(_b) && _b > 0))
      throw GeographicErr("Polar semi-axis is not positive");

    // Required number of terms in DST for full accuracy for all precisions as
    // a function of n in [-0.99, 0.99].  For each precision the P- and P+
    // columns list the values for negative and positive n.  Values determined
    // by running develop/AreaEst compiled with GEOGRAPHICLIB_PRECISION = 5.
    // For precision 4 and 5, GEOGRAPHICLIB_DIGITS was set to, resp., 384 and
    // 768.  The error criterion is relative error less than or equal to
    // epsilon/2 = 0.5^digits, with digits = 24, 53, 64, 113, 256.  The first 4
    // are the the "standard" values for float, double, long double, and
    // float128; the last is the default for GeographicLib + mpfr.
    //
    //         float    double     long d     quad        mpfr
    // n       f- f+    d-   d+    l-   l+    q-   q+     m-    m+
    // 0.01     4  4     6    6     8    8    16   16     48    48
    // 0.02     4  4     8    8    12   12    24   24     48    48
    // 0.03     4  4     8    8    12   12    24   24     48    48
    // 0.04     4  4    12   12    12   12    24   24     64    64
    // 0.05     4  4    12   12    12   12    24   24     64    64
    // 0.06     4  4    12   12    12   12    24   24     64    64
    // 0.07     4  4    12   12    16   16    32   32     64    64
    // 0.08     4  4    12   12    16   16    32   32     96    64
    // 0.09     4  4    12   12    16   16    32   32     96    96
    // 0.1      4  4    12   12    16   16    32   32     96    96
    // 0.11     6  4    12   12    16   16    32   32     96    96
    // 0.12     6  6    16   16    16   16    32   32     96    96
    // 0.13     6  6    16   16    24   16    48   32     96    96
    // 0.14     6  6    16   16    24   24    48   48     96    96
    // 0.15     6  6    16   16    24   24    48   48     96    96
    // 0.16     6  6    16   16    24   24    48   48     96    96
    // 0.17     6  6    16   16    24   24    48   48     96    96
    // 0.18     6  6    16   16    24   24    48   48     96    96
    // 0.19     6  6    16   16    24   24    48   48    128   128
    // 0.2      6  6    24   16    24   24    48   48    128   128
    // 0.21     6  6    24   24    24   24    48   48    128   128
    // 0.22     8  6    24   24    24   24    48   48    128   128
    // 0.23     8  6    24   24    24   24    48   48    128   128
    // 0.24     8  6    24   24    24   24    48   48    128   128
    // 0.25     8  6    24   24    32   24    48   48    128   128
    // 0.26     8  6    24   24    32   24    64   48    128   128
    // 0.27     8  8    24   24    32   32    64   64    128   128
    // 0.28     8  8    24   24    32   32    64   64    128   128
    // 0.29     8  8    24   24    32   32    64   64    192   192
    // 0.3      8  8    24   24    32   32    64   64    192   192
    // 0.31    12  8    24   24    32   32    64   64    192   192
    // 0.32    12  8    24   24    32   32    64   64    192   192
    // 0.33    12  8    24   24    32   32    64   64    192   192
    // 0.34    12  8    32   24    32   32    64   64    192   192
    // 0.35    12  8    32   24    32   32    64   64    192   192
    // 0.36    12  8    32   24    48   32    96   64    192   192
    // 0.37    12  8    32   32    48   32    96   64    192   192
    // 0.38    12  8    32   32    48   32    96   96    192   192
    // 0.39    12  8    32   32    48   48    96   96    192   192
    // 0.4     12  8    32   32    48   48    96   96    192   192
    // 0.41    12  8    32   32    48   48    96   96    192   192
    // 0.42    12 12    32   32    48   48    96   96    192   192
    // 0.43    12 12    32   32    48   48    96   96    256   192
    // 0.44    12 12    48   32    48   48    96   96    256   256
    // 0.45    12 12    48   32    48   48    96   96    256   256
    // 0.46    16 12    48   32    48   48    96   96    256   256
    // 0.47    16 12    48   32    48   48    96   96    256   256
    // 0.48    16 12    48   32    48   48    96   96    256   256
    // 0.49    16 12    48   48    48   48    96   96    256   256
    // 0.5     16 12    48   48    64   48    96   96    256   256
    // 0.51    16 12    48   48    64   48   128   96    256   256
    // 0.52    16 12    48   48    64   48   128   96    256   256
    // 0.53    16 12    48   48    64   48   128  128    256   256
    // 0.54    16 12    48   48    64   48   128  128    384   256
    // 0.55    16 12    48   48    64   64   128  128    384   384
    // 0.56    24 12    48   48    64   64   128  128    384   384
    // 0.57    24 12    48   48    64   64   128  128    384   384
    // 0.58    24 12    64   48    64   64   128  128    384   384
    // 0.59    24 12    64   48    64   64   128  128    384   384
    // 0.6     24 12    64   48    96   64   192  128    384   384
    // 0.61    24 12    64   48    96   64   192  128    384   384
    // 0.62    24 12    64   48    96   64   192  128    384   384
    // 0.63    24 16    64   48    96   64   192  192    384   384
    // 0.64    24 16    64   64    96   64   192  192    384   384
    // 0.65    24 16    64   64    96   96   192  192    384   384
    // 0.66    24 16    96   64    96   96   192  192    512   384
    // 0.67    24 16    96   64    96   96   192  192    512   512
    // 0.68    32 16    96   64    96   96   192  192    512   512
    // 0.69    32 16    96   64    96   96   192  192    512   512
    // 0.7     32 16    96   64    96   96   192  192    512   512
    // 0.71    32 16    96   64   128   96   192  192    512   512
    // 0.72    32 16    96   64   128   96   256  192    512   512
    // 0.73    32 16    96   96   128   96   256  192    768   512
    // 0.74    32 16    96   96   128   96   256  256    768   768
    // 0.75    48 16    96   96   128   96   256  256    768   768
    // 0.76    48 16   128   96   128  128   256  256    768   768
    // 0.77    48 24   128   96   192  128   256  256    768   768
    // 0.78    48 24   128   96   192  128   384  256    768   768
    // 0.79    48 24   128   96   192  128   384  256    768   768
    // 0.8     48 24   128   96   192  128   384  384    768   768
    // 0.81    48 24   128   96   192  128   384  384   1024   768
    // 0.82    48 24   192   96   192  192   384  384   1024  1024
    // 0.83    64 24   192  128   192  192   384  384   1024  1024
    // 0.84    64 24   192  128   256  192   384  384   1024  1024
    // 0.85    64 24   192  128   256  192   512  384   1024  1024
    // 0.86    64 24   192  128   256  192   512  384   1536  1024
    // 0.87    96 24   192  128   256  192   512  512   1536  1536
    // 0.88    96 24   256  192   384  192   768  512   1536  1536
    // 0.89    96 24   256  192   384  256   768  512   1536  1536
    // 0.9     96 24   256  192   384  256   768  768   2048  1536
    // 0.91   128 24   384  192   384  256   768  768   2048  2048
    // 0.92   128 24   384  192   512  384  1024  768   2048  2048
    // 0.93   192 24   384  256   512  384  1024  768   3072  3072
    // 0.94   192 24   512  256   768  384  1536 1024   3072  3072
    // 0.95   192 24   768  384   768  512  1536 1024   4096  3072
    // 0.96   256 24   768  384  1024  512  2048 1536   4096  4096
    // 0.97   384 24  1024  512  1536  768  3072 2048   6144  6144
    // 0.98   512 16  1536  768  2048 1024  4096 3072   8192  8192
    // 0.99  1024  8  3072 1024  4096 1536  8192 6144  16384 16384
    real n = 100 * _n;
    int N;
#if GEOGRAPHICLIB_PRECISION == 1
    if      (n >= -10 && n <= 11) N = 4;
    else if (n >= -21 && n <= 26) N = 6;
    else if (n >= -30 && n <= 41) N = 8;
    else if (n >= -45 && n <= 62) N = 12;
    else if (n >= -55 && n <= 76) N = 16;
    // ignore the *reduction* in N for n > 97
    else if (n >= -67 && n <= 99) N = 24;
    else if (n >= -74 && n <= 99) N = 32;
    else if (n >= -82 && n <= 99) N = 48;
    else if (n >= -86 && n <= 99) N = 64;
    else if (n >= -90 && n <= 99) N = 96;
    else if (n >= -92 && n <= 99) N = 128;
    else if (n >= -95 && n <= 99) N = 192;
    else if (n >= -96 && n <= 99) N = 256;
    else if (n >= -97 && n <= 99) N = 384;
    else if (n >= -98 && n <= 99) N = 512;
    else                          N = 1024;
#elif GEOGRAPHICLIB_PRECISION == 2
    if      (n >= - 1 && n <=  1) N = 6;
    else if (n >= - 3 && n <=  3) N = 8;
    else if (n >= -11 && n <= 11) N = 12;
    else if (n >= -19 && n <= 20) N = 16;
    else if (n >= -33 && n <= 36) N = 24;
    else if (n >= -43 && n <= 48) N = 32;
    else if (n >= -57 && n <= 63) N = 48;
    else if (n >= -65 && n <= 72) N = 64;
    else if (n >= -75 && n <= 82) N = 96;
    else if (n >= -81 && n <= 87) N = 128;
    else if (n >= -87 && n <= 92) N = 192;
    else if (n >= -90 && n <= 94) N = 256;
    else if (n >= -93 && n <= 96) N = 384;
    else if (n >= -94 && n <= 97) N = 512;
    else if (n >= -96 && n <= 98) N = 768;
    else if (n >= -97 && n <= 99) N = 1024;
    else if (n >= -98 && n <= 99) N = 1536;
    else                          N = 3072;
#elif GEOGRAPHICLIB_PRECISION == 3
    if      (n >= - 1 && n <=  1) N = 8;
    else if (n >= - 6 && n <=  6) N = 12;
    else if (n >= -12 && n <= 13) N = 16;
    else if (n >= -23 && n <= 26) N = 24;
    else if (n >= -35 && n <= 38) N = 32;
    else if (n >= -49 && n <= 54) N = 48;
    else if (n >= -59 && n <= 64) N = 64;
    else if (n >= -70 && n <= 75) N = 96;
    else if (n >= -76 && n <= 81) N = 128;
    else if (n >= -83 && n <= 88) N = 192;
    else if (n >= -87 && n <= 91) N = 256;
    else if (n >= -91 && n <= 94) N = 384;
    else if (n >= -93 && n <= 96) N = 512;
    else if (n >= -95 && n <= 97) N = 768;
    else if (n >= -96 && n <= 98) N = 1024;
    else if (n >= -97 && n <= 99) N = 1536;
    else if (n >= -98 && n <= 99) N = 2048;
    else                          N = 4096;
#elif GEOGRAPHICLIB_PRECISION == 4
    if      (n >= - 1 && n <=  1) N = 16;
    else if (n >= - 6 && n <=  6) N = 24;
    else if (n >= -12 && n <= 13) N = 32;
    else if (n >= -25 && n <= 26) N = 48;
    else if (n >= -35 && n <= 37) N = 64;
    else if (n >= -50 && n <= 52) N = 96;
    else if (n >= -59 && n <= 62) N = 128;
    else if (n >= -71 && n <= 73) N = 192;
    else if (n >= -77 && n <= 79) N = 256;
    else if (n >= -84 && n <= 86) N = 384;
    else if (n >= -87 && n <= 89) N = 512;
    else if (n >= -91 && n <= 93) N = 768;
    else if (n >= -93 && n <= 95) N = 1024;
    else if (n >= -95 && n <= 96) N = 1536;
    else if (n >= -96 && n <= 97) N = 2048;
    else if (n >= -97 && n <= 98) N = 3072;
    else if (n >= -98 && n <= 98) N = 4096;
    else if (n >= -98 && n <= 99) N = 6144;
    else                          N = 8192;
#elif GEOGRAPHICLIB_PRECISION == 5
    if      (n >= - 3 && n <=  3) N = 48;
    else if (n >= - 7 && n <=  8) N = 64;
    else if (n >= -18 && n <= 18) N = 96;
    else if (n >= -28 && n <= 28) N = 128;
    else if (n >= -42 && n <= 43) N = 192;
    else if (n >= -53 && n <= 54) N = 256;
    else if (n >= -65 && n <= 66) N = 384;
    else if (n >= -72 && n <= 73) N = 512;
    else if (n >= -80 && n <= 81) N = 768;
    else if (n >= -85 && n <= 86) N = 1024;
    else if (n >= -89 && n <= 90) N = 1536;
    else if (n >= -92 && n <= 92) N = 2048;
    else if (n >= -94 && n <= 95) N = 3072;
    else if (n >= -96 && n <= 96) N = 4096;
    else if (n >= -97 && n <= 97) N = 6144;
    else if (n >= -98 && n <= 98) N = 8192;
    else                          N = 16384;
    if (Math::digits() > 256) {
      // Scale up N by the number of digits in the precision relative to
      // the number used for the test = 256.
      int M = (Math::digits() * N) / 256;
      while (N < M) N = N % 3 == 0 ? 4*N/3 : 3*N/2;
    }
#else
#error "Bad value for GEOGRAPHICLIB_PRECISION"
#endif
    _fft.reset(N);
    _nC4 = N;
  }

  const GeodesicExact& GeodesicExact::WGS84() {
    static const GeodesicExact wgs84(Constants::WGS84_a(),
                                     Constants::WGS84_f());
    return wgs84;
  }

  GeodesicLineExact GeodesicExact::Line(real lat1, real lon1, real azi1,
                                        unsigned caps) const {
    return GeodesicLineExact(*this, lat1, lon1, azi1, caps);
  }

  Math::real GeodesicExact::GenDirect(real lat1, real lon1, real azi1,
                                      bool arcmode, real s12_a12,
                                      unsigned outmask,
                                      real& lat2, real& lon2, real& azi2,
                                      real& s12, real& m12,
                                      real& M12, real& M21,
                                      real& S12) const {
    // Automatically supply DISTANCE_IN if necessary
    if (!arcmode) outmask |= DISTANCE_IN;
    return GeodesicLineExact(*this, lat1, lon1, azi1, outmask)
      .                         // Note the dot!
      GenPosition(arcmode, s12_a12, outmask,
                  lat2, lon2, azi2, s12, m12, M12, M21, S12);
  }

  GeodesicLineExact GeodesicExact::GenDirectLine(real lat1, real lon1,
                                                 real azi1,
                                                 bool arcmode, real s12_a12,
                                                 unsigned caps) const {
    azi1 = Math::AngNormalize(azi1);
    real salp1, calp1;
    // Guard against underflow in salp0.  Also -0 is converted to +0.
    Math::sincosd(Math::AngRound(azi1), salp1, calp1);
    // Automatically supply DISTANCE_IN if necessary
    if (!arcmode) caps |= DISTANCE_IN;
    return GeodesicLineExact(*this, lat1, lon1, azi1, salp1, calp1,
                             caps, arcmode, s12_a12);
  }

  GeodesicLineExact GeodesicExact::DirectLine(real lat1, real lon1,
                                              real azi1, real s12,
                                              unsigned caps) const {
    return GenDirectLine(lat1, lon1, azi1, false, s12, caps);
  }

  GeodesicLineExact GeodesicExact::ArcDirectLine(real lat1, real lon1,
                                                 real azi1, real a12,
                                                 unsigned caps) const {
    return GenDirectLine(lat1, lon1, azi1, true, a12, caps);
  }

  Math::real GeodesicExact::GenInverse(real lat1, real lon1,
                                       real lat2, real lon2,
                                       unsigned outmask, real& s12,
                                       real& salp1, real& calp1,
                                       real& salp2, real& calp2,
                                       real& m12, real& M12, real& M21,
                                       real& S12) const {
    // Compute longitude difference (AngDiff does this carefully).  Result is
    // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
    // east-going and meridional geodesics.
    using std::isnan;           // Needed for Centos 7, ubuntu 14
    real lon12s, lon12 = Math::AngDiff(lon1, lon2, lon12s);
    // Make longitude difference positive.
    int lonsign = signbit(lon12) ? -1 : 1;
    lon12 *= lonsign; lon12s *= lonsign;
    real
      lam12 = lon12 * Math::degree(),
      slam12, clam12;
    // Calculate sincos of lon12 + error (this applies AngRound internally).
    Math::sincosde(lon12, lon12s, slam12, clam12);
    // the supplementary longitude difference
    lon12s = (Math::hd - lon12) - lon12s;

    // If really close to the equator, treat as on equator.
    lat1 = Math::AngRound(Math::LatFix(lat1));
    lat2 = Math::AngRound(Math::LatFix(lat2));
    // Swap points so that point with higher (abs) latitude is point 1
    // If one latitude is a nan, then it becomes lat1.
    int swapp = fabs(lat1) < fabs(lat2) || isnan(lat2) ? -1 : 1;
    if (swapp < 0) {
      lonsign *= -1;
      swap(lat1, lat2);
    }
    // Make lat1 <= -0
    int latsign = signbit(lat1) ? 1 : -1;
    lat1 *= latsign;
    lat2 *= latsign;
    // Now we have
    //
    //     0 <= lon12 <= 180
    //     -90 <= lat1 <= -0
    //     lat1 <= lat2 <= -lat1
    //
    // longsign, swapp, latsign register the transformation to bring the
    // coordinates to this canonical form.  In all cases, 1 means no change was
    // made.  We make these transformations so that there are few cases to
    // check, e.g., on verifying quadrants in atan2.  In addition, this
    // enforces some symmetries in the results returned.

    real sbet1, cbet1, sbet2, cbet2, s12x, m12x;
    // Initialize for the meridian.  No longitude calculation is done in this
    // case to let the parameter default to 0.
    EllipticFunction E(-_ep2);

    Math::sincosd(lat1, sbet1, cbet1); sbet1 *= _f1;
    // Ensure cbet1 = +epsilon at poles; doing the fix on beta means that sig12
    // will be <= 2*tiny for two points at the same pole.
    Math::norm(sbet1, cbet1); cbet1 = fmax(tiny_, cbet1);

    Math::sincosd(lat2, sbet2, cbet2); sbet2 *= _f1;
    // Ensure cbet2 = +epsilon at poles
    Math::norm(sbet2, cbet2); cbet2 = fmax(tiny_, cbet2);

    // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
    // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
    // a better measure.  This logic is used in assigning calp2 in Lambda12.
    // Sometimes these quantities vanish and in that case we force bet2 = +/-
    // bet1 exactly.  An example where is is necessary is the inverse problem
    // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
    // which failed with Visual Studio 10 (Release and Debug)

    if (cbet1 < -sbet1) {
      if (cbet2 == cbet1)
        sbet2 = copysign(sbet1, sbet2);
    } else {
      if (fabs(sbet2) == -sbet1)
        cbet2 = cbet1;
    }

    real
      dn1 = (_f >= 0 ? sqrt(1 + _ep2 * Math::sq(sbet1)) :
             sqrt(1 - _e2 * Math::sq(cbet1)) / _f1),
      dn2 = (_f >= 0 ? sqrt(1 + _ep2 * Math::sq(sbet2)) :
             sqrt(1 - _e2 * Math::sq(cbet2)) / _f1);

    real a12, sig12;

    bool meridian = lat1 == -Math::qd || slam12 == 0;

    if (meridian) {

      // Endpoints are on a single full meridian, so the geodesic might lie on
      // a meridian.

      calp1 = clam12; salp1 = slam12; // Head to the target longitude
      calp2 = 1; salp2 = 0;           // At the target we're heading north

      real
        // tan(bet) = tan(sig) * cos(alp)
        ssig1 = sbet1, csig1 = calp1 * cbet1,
        ssig2 = sbet2, csig2 = calp2 * cbet2;

      // sig12 = sig2 - sig1
      sig12 = atan2(fmax(real(0), csig1 * ssig2 - ssig1 * csig2),
                                  csig1 * csig2 + ssig1 * ssig2);
      {
        real dummy;
        Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                cbet1, cbet2, outmask | REDUCEDLENGTH,
                s12x, m12x, dummy, M12, M21);
      }
      // Add the check for sig12 since zero length geodesics might yield m12 <
      // 0.  Test case was
      //
      //    echo 20.001 0 20.001 0 | GeodSolve -i
      //
      // In fact, we will have sig12 > pi/2 for meridional geodesic which is
      // not a shortest path.
      if (sig12 < 1 || m12x >= 0) {
        // Need at least 2, to handle 90 0 90 180
        if (sig12 < 3 * tiny_ ||
            // Prevent negative s12 or m12 for short lines
            (sig12 < tol0_ && (s12x < 0 || m12x < 0)))
          sig12 = m12x = s12x = 0;
        m12x *= _b;
        s12x *= _b;
        a12 = sig12 / Math::degree();
      } else
        // m12 < 0, i.e., prolate and too close to anti-podal
        meridian = false;
    }

    // somg12 == 2 marks that it needs to be calculated
    real omg12 = 0, somg12 = 2, comg12 = 0;
    if (!meridian &&
        sbet1 == 0 &&   // and sbet2 == 0
        (_f <= 0 || lon12s >= _f * Math::hd)) {

      // Geodesic runs along equator
      calp1 = calp2 = 0; salp1 = salp2 = 1;
      s12x = _a * lam12;
      sig12 = omg12 = lam12 / _f1;
      m12x = _b * sin(sig12);
      if (outmask & GEODESICSCALE)
        M12 = M21 = cos(sig12);
      a12 = lon12 / _f1;

    } else if (!meridian) {

      // Now point1 and point2 belong within a hemisphere bounded by a
      // meridian and geodesic is neither meridional or equatorial.

      // Figure a starting point for Newton's method
      real dnm;
      sig12 = InverseStart(E, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                           lam12, slam12, clam12,
                           salp1, calp1, salp2, calp2, dnm);

      if (sig12 >= 0) {
        // Short lines (InverseStart sets salp2, calp2, dnm)
        s12x = sig12 * _b * dnm;
        m12x = Math::sq(dnm) * _b * sin(sig12 / dnm);
        if (outmask & GEODESICSCALE)
          M12 = M21 = cos(sig12 / dnm);
        a12 = sig12 / Math::degree();
        omg12 = lam12 / (_f1 * dnm);
      } else {

        // Newton's method.  This is a straightforward solution of f(alp1) =
        // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
        // root in the interval (0, pi) and its derivative is positive at the
        // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
        // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
        // maintained which brackets the root and with each evaluation of
        // f(alp) the range is shrunk, if possible.  Newton's method is
        // restarted whenever the derivative of f is negative (because the new
        // value of alp1 is then further from the solution) or if the new
        // estimate of alp1 lies outside (0,pi); in this case, the new starting
        // guess is taken to be (alp1a + alp1b) / 2.
        //
        // initial values to suppress warnings (if loop is executed 0 times)
        real ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, domg12 = 0;
        unsigned numit = 0;
        // Bracketing range
        real salp1a = tiny_, calp1a = 1, salp1b = tiny_, calp1b = -1;
        for (bool tripn = false, tripb = false;
             numit < maxit2_ || GEOGRAPHICLIB_PANIC;
             ++numit) {
          // 1/4 meridian = 10e6 m and random input.  max err is estimated max
          // error in nm (checking solution of inverse problem by direct
          // solution).  iter is mean and sd of number of iterations
          //
          //           max   iter
          // log2(b/a) err mean  sd
          //    -7     387 5.33 3.68
          //    -6     345 5.19 3.43
          //    -5     269 5.00 3.05
          //    -4     210 4.76 2.44
          //    -3     115 4.55 1.87
          //    -2      69 4.35 1.38
          //    -1      36 4.05 1.03
          //     0      15 0.01 0.13
          //     1      25 5.10 1.53
          //     2      96 5.61 2.09
          //     3     318 6.02 2.74
          //     4     985 6.24 3.22
          //     5    2352 6.32 3.44
          //     6    6008 6.30 3.45
          //     7   19024 6.19 3.30
          real dv;
          real v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                            slam12, clam12,
                            salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                            E, domg12, numit < maxit1_, dv);
          // Reversed test to allow escape with NaNs
          if (tripb || !(fabs(v) >= (tripn ? 8 : 1) * tol0_)) break;
          // Update bracketing values
          if (v > 0 && (numit > maxit1_ || calp1/salp1 > calp1b/salp1b))
            { salp1b = salp1; calp1b = calp1; }
          else if (v < 0 && (numit > maxit1_ || calp1/salp1 < calp1a/salp1a))
            { salp1a = salp1; calp1a = calp1; }
          if (numit < maxit1_ && dv > 0) {
            real
              dalp1 = -v/dv;
            // |dalp1| < pi test moved earlier because GEOGRAPHICLIB_PRECISION
            // = 5 can result in dalp1 = 10^(10^8).  Then sin(dalp1) takes ages
            // (because of the need to do accurate range reduction).
            if (fabs(dalp1) < Math::pi()) {
              real
                sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
                nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
              if (nsalp1 > 0) {
                calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                salp1 = nsalp1;
                Math::norm(salp1, calp1);
                // In some regimes we don't get quadratic convergence because
                // slope -> 0.  So use convergence conditions based on epsilon
                // instead of sqrt(epsilon).
                tripn = fabs(v) <= 16 * tol0_;
                continue;
              }
            }
          }
          // Either dv was not positive or updated value was outside legal
          // range.  Use the midpoint of the bracket as the next estimate.
          // This mechanism is not needed for the WGS84 ellipsoid, but it does
          // catch problems with more eccentric ellipsoids.  Its efficacy is
          // such for the WGS84 test set with the starting guess set to alp1 =
          // 90deg:
          // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
          // WGS84 and random input: mean = 4.74, sd = 0.99
          salp1 = (salp1a + salp1b)/2;
          calp1 = (calp1a + calp1b)/2;
          Math::norm(salp1, calp1);
          tripn = false;
          tripb = (fabs(salp1a - salp1) + (calp1a - calp1) < tolb_ ||
                   fabs(salp1 - salp1b) + (calp1 - calp1b) < tolb_);
        }
        {
          real dummy;
          Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                  cbet1, cbet2, outmask, s12x, m12x, dummy, M12, M21);
        }
        m12x *= _b;
        s12x *= _b;
        a12 = sig12 / Math::degree();
        if (outmask & AREA) {
          // omg12 = lam12 - domg12
          real sdomg12 = sin(domg12), cdomg12 = cos(domg12);
          somg12 = slam12 * cdomg12 - clam12 * sdomg12;
          comg12 = clam12 * cdomg12 + slam12 * sdomg12;
        }
      }
    }

    if (outmask & DISTANCE)
      s12 = real(0) + s12x;     // Convert -0 to 0

    if (outmask & REDUCEDLENGTH)
      m12 = real(0) + m12x;     // Convert -0 to 0

    if (outmask & AREA) {
      real
        // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
        salp0 = salp1 * cbet1,
        calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0
      real alp12,
        // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
        A4 = Math::sq(_a) * calp0 * salp0 * _e2;
      if (A4 != 0) {
        real
          k2 = Math::sq(calp0) * _ep2,
          // From Lambda12: tan(bet) = tan(sig) * cos(alp)
          ssig1 = sbet1, csig1 = calp1 * cbet1,
          ssig2 = sbet2, csig2 = calp2 * cbet2;
        Math::norm(ssig1, csig1);
        Math::norm(ssig2, csig2);
        I4Integrand i4(_ep2, k2);
        vector<real> C4a(_nC4);
        _fft.transform(i4, C4a.data());
        real
          B41 = DST::integral(ssig1, csig1, C4a.data(), _nC4),
          B42 = DST::integral(ssig2, csig2, C4a.data(), _nC4);
        S12 = A4 * (B42 - B41);
      } else
        // Avoid problems with indeterminate sig1, sig2 on equator
        S12 = 0;

      if (!meridian && somg12 == 2) {
        somg12 = sin(omg12); comg12 = cos(omg12);
      }

      if (!meridian &&
          // omg12 < 3/4 * pi
          comg12 > -real(0.7071) &&     // Long difference not too big
          sbet2 - sbet1 < real(1.75)) { // Lat difference not too big
        // Use tan(Gamma/2) = tan(omg12/2)
        // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
        // with tan(x/2) = sin(x)/(1+cos(x))
        real domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
        alp12 = 2 * atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                           domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
      } else {
        // alp12 = alp2 - alp1, used in atan2 so no need to normalize
        real
          salp12 = salp2 * calp1 - calp2 * salp1,
          calp12 = calp2 * calp1 + salp2 * salp1;
        // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
        // salp12 = -0 and alp12 = -180.  However this depends on the sign
        // being attached to 0 correctly.  The following ensures the correct
        // behavior.
        if (salp12 == 0 && calp12 < 0) {
          salp12 = tiny_ * calp1;
          calp12 = -1;
        }
        alp12 = atan2(salp12, calp12);
      }
      S12 += _c2 * alp12;
      S12 *= swapp * lonsign * latsign;
      // Convert -0 to 0
      S12 += 0;
    }

    // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
    if (swapp < 0) {
      swap(salp1, salp2);
      swap(calp1, calp2);
      if (outmask & GEODESICSCALE)
        swap(M12, M21);
    }

    salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
    salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

    // Returned value in [0, 180]
    return a12;
  }

  Math::real GeodesicExact::GenInverse(real lat1, real lon1,
                                       real lat2, real lon2,
                                       unsigned outmask,
                                       real& s12, real& azi1, real& azi2,
                                       real& m12, real& M12, real& M21,
                                       real& S12) const {
    outmask &= OUT_MASK;
    real salp1, calp1, salp2, calp2,
      a12 =  GenInverse(lat1, lon1, lat2, lon2,
                        outmask, s12, salp1, calp1, salp2, calp2,
                        m12, M12, M21, S12);
    if (outmask & AZIMUTH) {
      azi1 = Math::atan2d(salp1, calp1);
      azi2 = Math::atan2d(salp2, calp2);
    }
    return a12;
  }

  GeodesicLineExact GeodesicExact::InverseLine(real lat1, real lon1,
                                               real lat2, real lon2,
                                               unsigned caps) const {
    real t, salp1, calp1, salp2, calp2,
      a12 = GenInverse(lat1, lon1, lat2, lon2,
                       // No need to specify AZIMUTH here
                       0u, t, salp1, calp1, salp2, calp2,
                       t, t, t, t),
      azi1 = Math::atan2d(salp1, calp1);
    // Ensure that a12 can be converted to a distance
    if (caps & (OUT_MASK & DISTANCE_IN)) caps |= DISTANCE;
    return GeodesicLineExact(*this, lat1, lon1, azi1, salp1, calp1, caps,
                             true, a12);
  }

  void GeodesicExact::Lengths(const EllipticFunction& E,
                              real sig12,
                              real ssig1, real csig1, real dn1,
                              real ssig2, real csig2, real dn2,
                              real cbet1, real cbet2, unsigned outmask,
                              real& s12b, real& m12b, real& m0,
                              real& M12, real& M21) const {
    // Return m12b = (reduced length)/_b; also calculate s12b = distance/_b,
    // and m0 = coefficient of secular term in expression for reduced length.

    outmask &= OUT_ALL;
    // outmask & DISTANCE: set s12b
    // outmask & REDUCEDLENGTH: set m12b & m0
    // outmask & GEODESICSCALE: set M12 & M21

    // It's OK to have repeated dummy arguments,
    // e.g., s12b = m0 = M12 = M21 = dummy

    if (outmask & DISTANCE)
      // Missing a factor of _b
      s12b = E.E() / (Math::pi() / 2) *
        (sig12 + (E.deltaE(ssig2, csig2, dn2) - E.deltaE(ssig1, csig1, dn1)));
    if (outmask & (REDUCEDLENGTH | GEODESICSCALE)) {
      real
        m0x = - E.k2() * E.D() / (Math::pi() / 2),
        J12 = m0x *
        (sig12 + (E.deltaD(ssig2, csig2, dn2) - E.deltaD(ssig1, csig1, dn1)));
      if (outmask & REDUCEDLENGTH) {
        m0 = m0x;
        // Missing a factor of _b.  Add parens around (csig1 * ssig2) and
        // (ssig1 * csig2) to ensure accurate cancellation in the case of
        // coincident points.
        m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) -
          csig1 * csig2 * J12;
      }
      if (outmask & GEODESICSCALE) {
        real csig12 = csig1 * csig2 + ssig1 * ssig2;
        real t = _ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
        M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
        M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
      }
    }
  }

  Math::real GeodesicExact::Astroid(real x, real y) {
    // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
    // This solution is adapted from Geocentric::Reverse.
    real k;
    real
      p = Math::sq(x),
      q = Math::sq(y),
      r = (p + q - 1) / 6;
    if ( !(q == 0 && r <= 0) ) {
      real
        // Avoid possible division by zero when r = 0 by multiplying equations
        // for s and t by r^3 and r, resp.
        S = p * q / 4,            // S = r^3 * s
        r2 = Math::sq(r),
        r3 = r * r2,
        // The discriminant of the quadratic equation for T3.  This is zero on
        // the evolute curve p^(1/3)+q^(1/3) = 1
        disc = S * (S + 2 * r3);
      real u = r;
      if (disc >= 0) {
        real T3 = S + r3;
        // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
        // of precision due to cancellation.  The result is unchanged because
        // of the way the T is used in definition of u.
        T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
        // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
        real T = cbrt(T3); // T = r * t
        // T can be zero; but then r2 / T -> 0.
        u += T + (T != 0 ? r2 / T : 0);
      } else {
        // T is complex, but the way u is defined the result is real.
        real ang = atan2(sqrt(-disc), -(S + r3));
        // There are three possible cube roots.  We choose the root which
        // avoids cancellation.  Note that disc < 0 implies that r < 0.
        u += 2 * r * cos(ang / 3);
      }
      real
        v = sqrt(Math::sq(u) + q),    // guaranteed positive
        // Avoid loss of accuracy when u < 0.
        uv = u < 0 ? q / (v - u) : u + v, // u+v, guaranteed positive
        w = (uv - q) / (2 * v);           // positive?
      // Rearrange expression for k to avoid loss of accuracy due to
      // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
      k = uv / (sqrt(uv + Math::sq(w)) + w);   // guaranteed positive
    } else {               // q == 0 && r <= 0
      // y = 0 with |x| <= 1.  Handle this case directly.
      // for y small, positive root is k = abs(y)/sqrt(1-x^2)
      k = 0;
    }
    return k;
  }

  Math::real GeodesicExact::InverseStart(EllipticFunction& E,
                                         real sbet1, real cbet1, real dn1,
                                         real sbet2, real cbet2, real dn2,
                                         real lam12, real slam12, real clam12,
                                         real& salp1, real& calp1,
                                         // Only updated if return val >= 0
                                         real& salp2, real& calp2,
                                         // Only updated for short lines
                                         real& dnm) const {
    // Return a starting point for Newton's method in salp1 and calp1 (function
    // value is -1).  If Newton's method doesn't need to be used, return also
    // salp2 and calp2 and function value is sig12.
    real
      sig12 = -1,               // Return value
      // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
      sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
      cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
    real sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
    bool shortline = cbet12 >= 0 && sbet12 < real(0.5) &&
      cbet2 * lam12 < real(0.5);
    real somg12, comg12;
    if (shortline) {
      real sbetm2 = Math::sq(sbet1 + sbet2);
      // sin((bet1+bet2)/2)^2
      // =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
      sbetm2 /= sbetm2 + Math::sq(cbet1 + cbet2);
      dnm = sqrt(1 + _ep2 * sbetm2);
      real omg12 = lam12 / (_f1 * dnm);
      somg12 = sin(omg12); comg12 = cos(omg12);
    } else {
      somg12 = slam12; comg12 = clam12;
    }

    salp1 = cbet2 * somg12;
    calp1 = comg12 >= 0 ?
      sbet12 + cbet2 * sbet1 * Math::sq(somg12) / (1 + comg12) :
      sbet12a - cbet2 * sbet1 * Math::sq(somg12) / (1 - comg12);

    real
      ssig12 = hypot(salp1, calp1),
      csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

    if (shortline && ssig12 < _etol2) {
      // really short lines
      salp2 = cbet1 * somg12;
      calp2 = sbet12 - cbet1 * sbet2 *
        (comg12 >= 0 ? Math::sq(somg12) / (1 + comg12) : 1 - comg12);
      Math::norm(salp2, calp2);
      // Set return value
      sig12 = atan2(ssig12, csig12);
    } else if (fabs(_n) > real(0.1) || // Skip astroid calc if too eccentric
               csig12 >= 0 ||
               ssig12 >= 6 * fabs(_n) * Math::pi() * Math::sq(cbet1)) {
      // Nothing to do, zeroth order spherical approximation is OK
    } else {
      // Scale lam12 and bet2 to x, y coordinate system where antipodal point
      // is at origin and singular point is at y = 0, x = -1.
      real x, y, lamscale, betscale;
      real lam12x = atan2(-slam12, -clam12); // lam12 - pi
      if (_f >= 0) {            // In fact f == 0 does not get here
        // x = dlong, y = dlat
        {
          real k2 = Math::sq(sbet1) * _ep2;
          E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
          lamscale = _e2/_f1 * cbet1 * 2 * E.H();
        }
        betscale = lamscale * cbet1;

        x = lam12x / lamscale;
        y = sbet12a / betscale;
      } else {                  // _f < 0
        // x = dlat, y = dlong
        real
          cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
          bet12a = atan2(sbet12a, cbet12a);
        real m12b, m0, dummy;
        // In the case of lon12 = 180, this repeats a calculation made in
        // Inverse.
        Lengths(E, Math::pi() + bet12a,
                sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                cbet1, cbet2, REDUCEDLENGTH, dummy, m12b, m0, dummy, dummy);
        x = -1 + m12b / (cbet1 * cbet2 * m0 * Math::pi());
        betscale = x < -real(0.01) ? sbet12a / x :
          -_f * Math::sq(cbet1) * Math::pi();
        lamscale = betscale / cbet1;
        y = lam12x / lamscale;
      }

      if (y > -tol1_ && x > -1 - xthresh_) {
        // strip near cut
        // Need real(x) here to cast away the volatility of x for min/max
        if (_f >= 0) {
          salp1 = fmin(real(1), -x); calp1 = - sqrt(1 - Math::sq(salp1));
        } else {
          calp1 = fmax(real(x > -tol1_ ? 0 : -1), x);
          salp1 = sqrt(1 - Math::sq(calp1));
        }
      } else {
        // Estimate alp1, by solving the astroid problem.
        //
        // Could estimate alpha1 = theta + pi/2, directly, i.e.,
        //   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
        //   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
        //
        // However, it's better to estimate omg12 from astroid and use
        // spherical formula to compute alp1.  This reduces the mean number of
        // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
        // (min 0 max 5).  The changes in the number of iterations are as
        // follows:
        //
        // change percent
        //    1       5
        //    0      78
        //   -1      16
        //   -2       0.6
        //   -3       0.04
        //   -4       0.002
        //
        // The histogram of iterations is (m = number of iterations estimating
        // alp1 directly, n = number of iterations estimating via omg12, total
        // number of trials = 148605):
        //
        //  iter    m      n
        //    0   148    186
        //    1 13046  13845
        //    2 93315 102225
        //    3 36189  32341
        //    4  5396      7
        //    5   455      1
        //    6    56      0
        //
        // Because omg12 is near pi, estimate work with omg12a = pi - omg12
        real k = Astroid(x, y);
        real
          omg12a = lamscale * ( _f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k );
        somg12 = sin(omg12a); comg12 = -cos(omg12a);
        // Update spherical estimate of alp1 using omg12 instead of lam12
        salp1 = cbet2 * somg12;
        calp1 = sbet12a - cbet2 * sbet1 * Math::sq(somg12) / (1 - comg12);
      }
    }
    // Sanity check on starting guess.  Backwards check allows NaN through.
    if (!(salp1 <= 0))
      Math::norm(salp1, calp1);
    else {
      salp1 = 1; calp1 = 0;
    }
    return sig12;
  }

  Math::real GeodesicExact::Lambda12(real sbet1, real cbet1, real dn1,
                                     real sbet2, real cbet2, real dn2,
                                     real salp1, real calp1,
                                     real slam120, real clam120,
                                     real& salp2, real& calp2,
                                     real& sig12,
                                     real& ssig1, real& csig1,
                                     real& ssig2, real& csig2,
                                     EllipticFunction& E,
                                     real& domg12,
                                     bool diffp, real& dlam12) const
    {

    if (sbet1 == 0 && calp1 == 0)
      // Break degeneracy of equatorial line.  This case has already been
      // handled.
      calp1 = -tiny_;

    real
      // sin(alp1) * cos(bet1) = sin(alp0)
      salp0 = salp1 * cbet1,
      calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0

    real somg1, comg1, somg2, comg2, somg12, comg12, cchi1, cchi2, lam12;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
    ssig1 = sbet1; somg1 = salp0 * sbet1;
    csig1 = comg1 = calp1 * cbet1;
    // Without normalization we have schi1 = somg1.
    cchi1 = _f1 * dn1 * comg1;
    Math::norm(ssig1, csig1);
    // Math::norm(somg1, comg1); -- don't need to normalize!
    // Math::norm(schi1, cchi1); -- don't need to normalize!

    // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
    // about this case, since this can yield singularities in the Newton
    // iteration.
    // sin(alp2) * cos(bet2) = sin(alp0)
    salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
    // calp2 = sqrt(1 - sq(salp2))
    //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
    // and subst for calp0 and rearrange to give (choose positive sqrt
    // to give alp2 in [0, pi/2]).
    calp2 = cbet2 != cbet1 || fabs(sbet2) != -sbet1 ?
      sqrt(Math::sq(calp1 * cbet1) +
           (cbet1 < -sbet1 ?
            (cbet2 - cbet1) * (cbet1 + cbet2) :
            (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
      fabs(calp1);
    // tan(bet2) = tan(sig2) * cos(alp2)
    // tan(omg2) = sin(alp0) * tan(sig2).
    ssig2 = sbet2; somg2 = salp0 * sbet2;
    csig2 = comg2 = calp2 * cbet2;
    // Without normalization we have schi2 = somg2.
    cchi2 = _f1 * dn2 * comg2;
    Math::norm(ssig2, csig2);
    // Math::norm(somg2, comg2); -- don't need to normalize!
    // Math::norm(schi2, cchi2); -- don't need to normalize!

    // sig12 = sig2 - sig1, limit to [0, pi]
    sig12 = atan2(fmax(real(0), csig1 * ssig2 - ssig1 * csig2),
                                csig1 * csig2 + ssig1 * ssig2);

    // omg12 = omg2 - omg1, limit to [0, pi]
    somg12 = fmax(real(0), comg1 * somg2 - somg1 * comg2);
    comg12 =               comg1 * comg2 + somg1 * somg2;
    real k2 = Math::sq(calp0) * _ep2;
    E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
    // chi12 = chi2 - chi1, limit to [0, pi]
    real
      schi12 = fmax(real(0), cchi1 * somg2 - somg1 * cchi2),
      cchi12 =               cchi1 * cchi2 + somg1 * somg2;
    // eta = chi12 - lam120
    real eta = atan2(schi12 * clam120 - cchi12 * slam120,
                     cchi12 * clam120 + schi12 * slam120);
    real deta12 = -_e2/_f1 * salp0 * E.H() / (Math::pi() / 2) *
      (sig12 + (E.deltaH(ssig2, csig2, dn2) - E.deltaH(ssig1, csig1, dn1)));
    lam12 = eta + deta12;
    // domg12 = deta12 + chi12 - omg12
    domg12 = deta12 + atan2(schi12 * comg12 - cchi12 * somg12,
                            cchi12 * comg12 + schi12 * somg12);
    if (diffp) {
      if (calp2 == 0)
        dlam12 = - 2 * _f1 * dn1 / sbet1;
      else {
        real dummy;
        Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                cbet1, cbet2, REDUCEDLENGTH,
                dummy, dlam12, dummy, dummy, dummy);
        dlam12 *= _f1 / (calp2 * cbet2);
      }
    }

    return lam12;
  }

  Math::real GeodesicExact::I4Integrand::asinhsqrt(real x) {
    // return asinh(sqrt(x))/sqrt(x)
    using std::sqrt; using std::asinh; using std::asin;
    return x == 0 ? 1 :
      (x > 0 ? asinh(sqrt(x))/sqrt(x) :
       asin(sqrt(-x))/sqrt(-x)); // NaNs end up here
  }
  Math::real GeodesicExact::I4Integrand::t(real x) {
    // This differs by from t as defined following Eq 61 in Karney (2013) by
    // the final subtraction of 1.  This changes nothing since Eq 61 uses the
    // difference of two evaluations of t and improves the accuracy(?).
    using std::sqrt;
    // Group terms to minimize roundoff
    // with x = ep2, this is the same as
    // e2/(1-e2) + (atanh(e)/e - 1)
    return x + (sqrt(1 + x) * asinhsqrt(x) - 1);
  }
  Math::real GeodesicExact::I4Integrand::td(real x) {
    // d t(x) / dx
    using std::sqrt;
    return x == 0 ? 4/real(3) :
      // Group terms to minimize roundoff
      1 + (1 - asinhsqrt(x) / sqrt(1+x)) / (2*x);
  }
  // Math::real GeodesicExact::I4Integrand::Dt(real x, real y) {
  //   // ( t(x) - t(y) ) / (x - y)
  //   using std::sqrt; using std::fabs; using std::asinh; using std::asin;
  //   if (x == y) return td(x);
  //   if (x * y <= 0) return ( t(x) - t(y) ) / (x - y);
  //   real
  //     sx = sqrt(fabs(x)), sx1 = sqrt(1 + x),
  //     sy = sqrt(fabs(y)), sy1 = sqrt(1 + y),
  //     z = (x - y) / (sx * sy1 + sy * sx1),
  //     d1 = 2 * sx * sy,
  //     d2 = 2 * (x * sy * sy1 + y * sx * sx1);
  //   return x > 0 ?
  //     ( 1 + (asinh(z)/z) / d1 - (asinh(sx) + asinh(sy)) / d2 ) :
  //     // NaNs fall through to here
  //     ( 1 - (asin (z)/z) / d1 - (asin (sx) + asin (sy)) / d2 );
  // }
  Math::real GeodesicExact::I4Integrand::DtX(real y) const {
    // idiot version:
    // return ( tX - t(y) ) / (X - y);
    using std::sqrt; using std::fabs; using std::asinh; using std::asin;
    if (X == y) return tdX;
    if (X * y <= 0) return ( tX - t(y) ) / (X - y);
    real
      sy = sqrt(fabs(y)), sy1 = sqrt(1 + y),
      z = (X - y) / (sX * sy1 + sy * sX1),
      d1 = 2 * sX * sy,
      d2 = 2 * (X * sy * sy1 + y * sXX1);
    return X > 0 ?
      ( 1 + (asinh(z)/z) / d1 - (asinhsX + asinh(sy)) / d2 ) :
      // NaNs fall through to here
      ( 1 - (asin (z)/z) / d1 - (asinhsX + asin (sy)) / d2 );
  }
  GeodesicExact::I4Integrand::I4Integrand(real ep2, real k2)
    : X( ep2 )
    , tX( t(X) )
    , tdX( td(X) )
    , _k2( k2 )
  {
    using std::fabs; using std::sqrt; using std::asinh; using std::asin;
    sX = sqrt(fabs(X));     // ep
    sX1 =  sqrt(1 + X);     // 1/(1-f)
    sXX1 = sX * sX1;
    asinhsX = X > 0 ? asinh(sX) : asin(sX); // atanh(e)
  }
  Math::real GeodesicExact::I4Integrand::operator()(real sig) const {
    using std::sin;
    real ssig = sin(sig);
    return - DtX(_k2 * Math::sq(ssig)) * ssig/2;
  }

} // namespace GeographicLib

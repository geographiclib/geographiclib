// Estimate the required order of the DST for a given GEOGRAPHICLIB_PRECISION
// and n in [-0.91, 0.91].  For each n, vary alp0 to determine the case with
// the maximum error.  Ensure that this is smaller by the machine epsilon than
// the spherical terms.

// This code allows either FFTW or kissfft to be used.  kissfft prefers because
// it allows the use of mpreal and mpreal is needed to reliably measure the
// error for quad precision.

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>

#include <GeographicLib/Math.hpp>
#include <GeographicLib/Utility.hpp>

#define FFTW 0
#if GEOGRAPHICLIB_PRECISION == 5 && FFTW
#undef FFTW
#define FFTW 0
#endif

#if FFTW
#include <fftw3.h>

#if GEOGRAPHICLIB_PRECISION == 1
#define fftw_2r2_kind fftwf_2r2_kind
#define fftw_plan fftwf_plan
#define fftw_plan_r2r_1d fftwf_plan_r2r_1d
#define fftw_execute fftwf_execute
#define fftw_destroy_plan fftwf_destroy_plan
#elif GEOGRAPHICLIB_PRECISION == 3
#define fftw_2r2_kind fftwl_2r2_kind
#define fftw_plan fftwl_plan
#define fftw_plan_r2r_1d fftwl_plan_r2r_1d
#define fftw_execute fftwl_execute
#define fftw_destroy_plan fftwl_destroy_plan
#elif GEOGRAPHICLIB_PRECISION == 4
#define fftw_2r2_kind fftwq_2r2_kind
#define fftw_plan fftwq_plan
#define fftw_plan_r2r_1d fftwq_plan_r2r_1d
#define fftw_execute fftwq_execute
#define fftw_destroy_plan fftwq_destroy_plan
#endif

#else
#include <../src/kissfft.hh>
#endif

using namespace GeographicLib;
using namespace std;

class I4Integrand {
  Math::real X, tX, tdX, sX, sX1, sXX1, asinhsX, _k2;
  // return asinh(sqrt(x))/sqrt(x)
  static Math::real asinhsqrt(Math::real x) {
    using std::sqrt; using std::asinh; using std::asin;
    return x == 0 ? 1 :
      (x > 0 ? asinh(sqrt(x))/sqrt(x) :
       asin(sqrt(-x))/sqrt(-x)); // NaNs end up here
  }
  // This differs by from t as defined following Eq 61 in Karney (2013) by
  // the final subtraction of 1.  This changes nothing since Eq 61 uses the
  // difference of two evaluations of t and improves the accuracy(?).
  static Math::real t(Math::real x) {
    using std::sqrt;
    // Group terms to minimize roundoff
    // with x = ep2, this is the same as
    // e2/(1-e2) + (atanh(e)/e - 1)
    return x + (sqrt(1 + x) * asinhsqrt(x) - 1);
  }
  // d t(x) / dx
  static Math::real td(Math::real x) {
    using std::sqrt;
    return x == 0 ? 4/Math::real(3) :
      // Group terms to minimize roundoff
      1 + (1 - asinhsqrt(x) / sqrt(1+x)) / (2*x);
  }
  // ( t(x) - t(y) ) / (x - y)
  static Math::real Dt(Math::real x, Math::real y) {
    using std::sqrt; using std::fabs; using std::asinh; using std::asin;
    if (x == y) return td(x);
    if (x * y <= 0) return ( t(x) - t(y) ) / (x - y);
    Math::real
      sx = sqrt(fabs(x)), sx1 = sqrt(1 + x),
      sy = sqrt(fabs(y)), sy1 = sqrt(1 + y),
      z = (x - y) / (sx * sy1 + sy * sx1),
      d1 = 2 * sx * sy,
      d2 = 2 * (x * sy * sy1 + y * sx * sx1);
    return x > 0 ?
      ( 1 + (asinh(z)/z) / d1 - (asinh(sx) + asinh(sy)) / d2 ) :
      // NaNs fall through to here
      ( 1 - (asin (z)/z) / d1 - (asin (sx) + asin (sy)) / d2 );
  }
  // ( t(X) - t(y) ) / (X - y)
  Math::real DtX(Math::real y) const {
    using std::sqrt; using std::fabs; using std::asinh; using std::asin;
    if (X == y) return tdX;
    if (X * y <= 0) return ( tX - t(y) ) / (X - y);
    Math::real
      sy = sqrt(fabs(y)), sy1 = sqrt(1 + y),
      z = (X - y) / (sX * sy1 + sy * sX1),
      d1 = 2 * sX * sy,
      d2 = 2 * (X * sy * sy1 + y * sXX1);
    return X > 0 ?
      ( 1 + (asinh(z)/z) / d1 - (asinhsX + asinh(sy)) / d2 ) :
      // NaNs fall through to here
      ( 1 - (asin (z)/z) / d1 - (asinhsX + asin (sy)) / d2 );
  }

public:
  I4Integrand(Math::real ep2, Math::real k2)
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
  Math::real operator()(Math::real sig) const {
    using std::sin;
    Math::real ssig = sin(sig);
    return - DtX(_k2 * Math::sq(ssig)) * ssig/2;
  }
};

Math::real CosSeries(Math::real sinx, Math::real cosx, const Math::real c[], int n) {
  // Evaluate
  // y = sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
  // using Clenshaw summation.
  // Approx operation count = (n + 5) mult and (2 * n + 2) add
  c += n ;                    // Point to one beyond last element
  Math::real
    ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
    y0 = n & 1 ? *--c : 0, y1 = 0;          // accumulators for sum
  // Now n is even
  n /= 2;
  while (n--) {
    // Unroll loop x 2, so accumulators return to their original role
    y1 = ar * y0 - y1 + *--c;
    y0 = ar * y1 - y0 + *--c;
  }
  return cosx * (y0 - y1);    // cos(x) * (y0 - y1)
}

Math::real SinSeries(Math::real sinx, Math::real cosx, const Math::real c[], int n) {
  // Evaluate
  // y = sum(c[i] * sin((2*i+1) * x), i, 0, n-1)
  // using Clenshaw summation.
  // Approx operation count = (n + 5) mult and (2 * n + 2) add
  c += n ;                    // Point to one beyond last element
  Math::real
    ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
    y0 = n & 1 ? *--c : 0, y1 = 0;          // accumulators for sum
  // Now n is even
  n /= 2;
  while (n--) {
    // Unroll loop x 2, so accumulators return to their original role
    y1 = ar * y0 - y1 + *--c;
    y0 = ar * y1 - y0 + *--c;
  }
  return sinx * (y0 + y1);    // sin(x) * (y0 + y1)
}

Math::real fft_check(const vector<Math::real>& vals,
                     const vector<Math::real>& tx,
                     bool centerp = false) {
  Math::real maxerr = 0;
  int N = vals.size();
  for (int i = 0; i < N; ++i) {
    Math::real
      sig = (2*i + (centerp ? 1 : 2))/Math::real(4*N) * Math::pi(),
      err = fabs(vals[i] - SinSeries(sin(sig), cos(sig),
                                     tx.data(), tx.size()));
    maxerr = fmax(err, maxerr);
  }
  return maxerr;
}

// Implement DST-III (centerp = false) or DST-IV (centerp = true)
void fft_transform(const vector<Math::real>& in, vector<Math::real>& out,
                   bool centerp = false, bool check = false) {
  int N = in.size(); out.resize(N);
#if FFTW
  fftw_r2r_kind kind = centerp ? FFTW_RODFT11 : FFTW_RODFT01;
  fftw_plan p;
#if GEOGRAPHICLIB_PRECISION == 4
  vector<__float128> temp(N);
  for (int i = 0; i < N; ++i) temp[i] = __float128(in[i]);
  p = fftw_plan_r2r_1d(N, temp.data(), temp.data(), kind, FFTW_ESTIMATE);
  fftw_execute(p);
  for (int i = 0; i < N; ++i) out[i] = Math::real(temp[i])/N;
#else
  out = in;
  p = fftw_plan_r2r_1d(N, out.data(), out.data(), kind, FFTW_ESTIMATE);
  fftw_execute(p);
  for (int i = 0; i < N; ++i) out[i] /= N;
#endif
  fftw_destroy_plan(p);
#else
  vector<Math::real> tempin(4*N);
  if (centerp) {
    for (int i = 0; i < N; ++i) {
      tempin[i] = in[i];
      tempin[N+i] = in[N-1-i];
      tempin[2*N+i] = -in[i];
      tempin[3*N+i] = -in[N-1-i];
    }
  } else {
    tempin[0] = 0;
    for (int i = 0; i < N; ++i) tempin[i+1] = in[i];
    for (int i = 1; i < N; ++i) tempin[N+i] = tempin[N-i];
    for (int i = 0; i < 2*N; ++i) tempin[2*N+i] = -tempin[i];
  }
  kissfft<Math::real> fft(2*N, false);
  vector<complex<Math::real>> tempout(2*N);
  fft.transform_real(tempin.data(), tempout.data());
  for (int i = 0, j = 1; i < N; ++i, j+=2) {
    if (centerp)
      tempout[j] *= exp(complex<Math::real>(0, -j*Math::pi()/(4*N)));
    out[i] = -tempout[j].imag() / (2*N);
  }
#endif
  if (check)
    cout << "err(" << centerp << ") = " << fft_check(in, out, centerp) << "\n";
}

void fft_transform2(const vector<Math::real>& oldin,
                    const vector<Math::real>& newin,
                    const vector<Math::real>& oldout,
                    vector<Math::real>& newout,
                    bool check = false) {
  // oldin and oldout are a transform pair with centerp = false.
  // newin at the corresponding centerp = true values

  // newout is the centerp = false transform for the combined set of values, so
  // newout.size() = 2*oldout.size()

  // N.B. the combined set of input values are order newin[0], oldin[0],
  // newin[1], oldin[1], newin[2], oldin[2], etc.

  // oldin is only touched with check = true
  int N = newin.size();
  fft_transform(newin, newout, true, false);
  newout.resize(2*N);
  for (int i = N; i < 2*N; ++i)
    newout[i] = (-oldout[2*N-1-i] + newout[2*N-1-i])/2;
  for (int i = 0; i < N; ++i)
    newout[i] = (oldout[i] + newout[i])/2;
  if (check) {
    vector<Math::real> tempin(2*N);
    for (int i = 0; i < N; ++i) {
      tempin[2*i  ] = newin[i];
      tempin[2*i+1] = oldin[i];
    }
    if (check)
      cout << "err(2) = " << fft_check(tempin, newout, false) << "\n";

  }
}

void I4f(Math::real n, Math::real alp0, int N,
         vector<Math::real>& I4, bool centerp = false, bool check = false) {
  Math::real
    ep2 = 4*n/Math::sq(1 - n),
    k2 = ep2 * Math::sq( Math::cosd(alp0) );
  I4Integrand i4(ep2, k2);
  vector<Math::real> in(N);
  for (int i = 0; i < N; ++i) {
    Math::real sig = (2*i + (centerp ? 1 : 2))/Math::real(4*N) * Math::pi();
    in[i] = i4(sig);
  }
  fft_transform(in, I4, centerp, check);
}

// Scale to do the integal -1/(2*i + 1) and to normalize term to estimate error A4/c2;

void C4f(Math::real n, Math::real alp0, int N,
         vector<Math::real>& C4, bool centerp = false) {
  I4f(n, alp0, N, C4, centerp);
  Math::real
    ep2 = 4*n/Math::sq(1 - n),
    e2 = 4*n/Math::sq(1 + n),
    c2 = (1 + Math::sq((1-n)/(1+n)) *
          (n == 0 ? 1 :
           (n > 0 ? asinh(sqrt(ep2)) : atan(sqrt(-e2))) /
           sqrt(fabs(e2))))/2, // authalic radius squared
    A4 = Math::cosd(alp0) * Math::sind(alp0) * e2;
  for (int i = 0; i < N; ++i)
    C4[i] *= - A4/((2*i + 1) * c2);
}

Math::real maxerr(const vector<Math::real> C4a, const vector<Math::real> C4b) {
  int Na = C4a.size(), Nb = C4b.size(), N0 = min(Na, Nb), N1 = max(Na, Nb);
  vector<Math::real> diff = Nb > Na ? C4b : C4a;
  for (int i = 0; i < N0; ++i)
    diff[i] -= Nb > Na ? C4a[i] : C4b[i];
  Math::real err = 0;
  if (true) {
    // Simplest (and most conservative) estimate of error is to sum the abs
    // difference of the coefficients
    for (int i = 0; i < N1; ++i) {
      err += fabs(diff[i]);
    }
  } else {
    // Step through sig looking for the largest error
    for (int i = 0; i < 4*N0; ++i) {
      Math::real sig = Math::pi()/2 * i / Math::real(4*N0);
      err = fmax(err, fabs( CosSeries(sin(sig), cos(sig), diff.data(), N1) ));
    }
  }
  return err;
}

int main(int argc, const char* const argv[]) {
  try {
    Utility::set_digits();
    int Nmax, prec;
    if (argc != 3) { cerr << "AreaEst Nmax prec\n"; return 1; }
    Nmax = Utility::val<int>(string(argv[1]));
    prec = Utility::val<int>(string(argv[2]));
    vector<Math::real> C4ref, C4;
    Math::real eps;
    switch (prec) {
    case 1: eps = numeric_limits<float>::epsilon() / 2; break;
    case 2: eps = numeric_limits<double>::epsilon() / 2; break;
    case 3: eps = numeric_limits<long double>::epsilon() / 2; break;
    case 4: eps = pow(0.5, 112) / 2; break;
    default: eps = numeric_limits<double>::epsilon() / 2; break;
    }
    // The range of n in [-0.91, 0.91] the includes 1/20 < b/a < 20
    for (int in = -91; in <= 91; ++in) {
      vector<Math::real> C4x;
      Math::real n = in/Math::real(100),
        errx = -1,
        maxalp0 = -1;
      int N = 0;
      // Pick N = 2^k and 3*2^k: [4, 6, 8, 12, 16, 24, 32, ...]
      for (N = 4; N <= Nmax; N = N % 3 == 0 ? 4*N/3 : 3*N/2) {
        errx = -1;
        maxalp0 = -1;
        Math::real alp0 = 10, err;
        C4f(n, alp0, 2*N, C4ref);
        C4f(n, alp0, N, C4);
        err = maxerr(C4, C4ref);
        if (err > eps) continue;
        bool ok = true;
        for (int a = 10; a < 90; ++a) {
          alp0 = a/Math::real(10);
          C4f(n, alp0, 2*N, C4ref);
          C4f(n, alp0, N, C4);
          err = maxerr(C4, C4ref);
          if (err > eps) { ok = false; break; }
          if (err > errx) {
            errx = err; C4x = C4;
            maxalp0 = alp0;
          }
        }
        if (!ok) continue;
        Math::real alp00 = maxalp0;
        for (int a = -9; a < 10; ++a) {
          alp0 = alp00 + a;
          C4f(n, alp0, 2*N, C4ref);
          C4f(n, alp0, N, C4);
          err = maxerr(C4, C4ref);
          if (err > eps) { ok = false; break; }
          if (err > errx) {
            errx = err; C4x = C4;
            maxalp0 = alp0;
          }
        }
        if (!ok) continue;
        alp00 = maxalp0;
        for (int a = -9; a < 10; ++a) {
          alp0 = alp00 + a/Math::real(10);
          C4f(n, alp0, 2*N, C4ref);
          C4f(n, alp0, N, C4);
          err = maxerr(C4, C4ref);
          if (err > eps) { ok = false; break; }
          if (err > errx) {
            errx = err; C4x = C4;
            maxalp0 = alp0;
          }
        }
        if (!ok) continue;
        alp00 = maxalp0;
        for (int a = -9; a < 10; ++a) {
          alp0 = alp00 + a/Math::real(100);
          C4f(n, alp0, 2*N, C4ref);
          C4f(n, alp0, N, C4);
          err = maxerr(C4, C4ref);
          if (err > eps) { ok = false; break; }
          if (err > errx) {
            errx = err; C4x = C4;
            maxalp0 = alp0;
          }
        }
        if (!ok) continue;
        alp00 = maxalp0;
        for (int a = -9; a < 10; ++a) {
          alp0 = alp00 + a/Math::real(1000);
          C4f(n, alp0, 2*N, C4ref);
          C4f(n, alp0, N, C4);
          err = maxerr(C4, C4ref);
          if (err > eps) { ok = false; break; }
          if (err > errx) {
            errx = err; C4x = C4;
            maxalp0 = alp0;
          }
        }
        if (!ok) continue;
        alp00 = maxalp0;
        for (int a = -9; a < 10; ++a) {
          alp0 = alp00 + a/Math::real(10000);
          C4f(n, alp0, 2*N, C4ref);
          C4f(n, alp0, N, C4);
          err = maxerr(C4, C4ref);
          if (err > eps) { ok = false; break; }
          if (err > errx) {
            errx = err; C4x = C4;
            maxalp0 = alp0;
          }
        }
        if (ok) break;
      }
      Math::real erry = 0;
      // Assess summing last few scaled coefficients as a less expensive
      // error metric.
      for (int i = (31*N)/32; i < N; ++i)
        erry = erry+fabs(C4x[i]);
      cout << n << " " << N << " " << maxalp0 << " " << errx << " "
           << erry/errx << endl;
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
}

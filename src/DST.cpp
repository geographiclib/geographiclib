/**
 * \file DST.cpp
 * \brief Implementation for GeographicLib::DST class
 *
 * Copyright (c) Charles Karney (2022) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/DST.hpp>
#include <memory>

#include "kissfft.hh"

namespace GeographicLib {

  using namespace std;

  DST::DST(unsigned N)
    : _fft(make_shared<fft_t>(fft_t(2*N, false)))
  {}

  void DST::reserve(unsigned N) {
    _fft->assign(2*N, false);
  }

  void DST::fft_transform(const vector<Math::real>& in,
                          vector<Math::real>& out, bool centerp) const {
    // Implement DST-III (centerp = false) or DST-IV (centerp = true)
    // in and out can be the same array
    int N = in.size(); out.resize(N);
    if (N == 0) return;
    vector<Math::real> tempin(4*N);
    if (centerp) {
      for (int i = 0; i < N; ++i) {
        tempin[i] = in[i];
        tempin[N+i] = in[N-1-i];
        tempin[2*N+i] = -in[i];
        tempin[3*N+i] = -in[N-1-i];
      }
    } else {
      tempin[0] = 0;            // set [0]
      for (int i = 0; i < N; ++i) tempin[i+1] = in[i]; // set [1,N]
      for (int i = 1; i < N; ++i) tempin[N+i] = tempin[N-i]; // set [N+1,2*N-1]
      for (int i = 0; i < 2*N; ++i) tempin[2*N+i] = -tempin[i]; // [2*N, 4*N-1]
    }
    vector<complex<Math::real>> tempout(2*N);
    _fft->assign(2*N, false);
    _fft->transform_real(tempin.data(), tempout.data());
    if (centerp) {
      Math::real d = -Math::pi()/(4*N);
      for (int i = 0, j = 1; i < N; ++i, j+=2)
        tempout[j] *= exp(complex<Math::real>(0, j*d));
    }
    for (int i = 0, j = 1; i < N; ++i, j+=2) {
      out[i] = -tempout[j].imag() / (2*N);
    }
  }

  void DST::fft_transform2(const vector<Math::real>& newin,
                           const vector<Math::real>& oldout,
                           vector<Math::real>& newout) const {
    // oldout is the transform for N points with centerp = false.
    // newin at the corresponding N points with centerp = true values

    // newout is the centerp = false transform for the combined set of values,
    // so newout.size() = 2*oldout.size()

    // newin and newout can be the same array
    // oldout and newout cannot be the same array
    int N = newin.size();
    if (oldout.size() != unsigned(N))
      throw GeographicErr("Mismatch of array sizes in DST::fft_transform2");
    fft_transform(newin, newout, true);
    newout.resize(2*N);
    for (int i = N; i < 2*N; ++i)
      newout[i] = (-oldout[2*N-1-i] + newout[2*N-1-i])/2;
    for (int i = 0; i < N; ++i)
      newout[i] = (oldout[i] + newout[i])/2;
  }

  // void DST::transform(const vector<Math::real>& x,
  //                     vector<Math::real>& F) const {
  //   fft_transform(x, F, false);
  // }

  void DST::transform(function<Math::real(Math::real)> f, int N,
                      vector<Math::real>& F) const {
    F.resize(N);
    Math::real d = Math::pi()/(2 * N);
    for (int i = 0; i < N; ++i)
      F[i] = f( (i + 1) * d );
    fft_transform(F, F, false);
  }

  void DST::refine(function<Math::real(Math::real)> f,
                   const vector<Math::real>& oldF,
                   vector<Math::real>& newF) const {
    // oldF and newF can be the same arrays
    int N = oldF.size();
    vector<Math::real> temp(N);
    Math::real d = Math::pi()/(4 * N);
    for (int i = 0; i < N; ++i)
      temp[i] = f( (2*i + 1) * d );
    fft_transform2(temp, oldF, temp);
    newF.swap(temp);
  }

  Math::real DST::eval(const vector<Math::real>& F,
                       Math::real sinx, Math::real cosx) {
    // Evaluate
    // y = sum(F[i] * sin((2*i+1) * x), i, 0, n-1)
    // using Clenshaw summation.
    // Approx operation count = (n + 5) mult and (2 * n + 2) add
    int n = F.size();
    Math::real
      ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
      y0 = n & 1 ? F[--n] : 0, y1 = 0;          // accumulators for sum
    // Now n is even
    while (n > 0) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + F[--n];
      y0 = ar * y1 - y0 + F[--n];
    }
    return sinx * (y0 + y1);    // sin(x) * (y0 + y1)
  }

  Math::real DST::evalx(const vector<Math::real>& F, Math::real x) {
    return eval(F, sin(x), cos(x));
  }

  Math::real DST::integral(const vector<Math::real>& F,
                           Math::real sinx, Math::real cosx) {
    // Evaluate
    // y = -sum(F[i]/(2*i+1) * cos((2*i+1) * x), i, 0, n-1)
    // using Clenshaw summation.
    // Approx operation count = (n + 5) mult and (2 * n + 2) add
    int n = F.size(), l = n;
    Math::real
      ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
      y0 = n & 1 ? F[--n]/(2*(--l)+1) : 0, y1 = 0; // accumulators for sum
    // Now n is even
    while (n > 0) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + F[--n]/(2*(--l)+1);
      y0 = ar * y1 - y0 + F[--n]/(2*(--l)+1);
    }
    return cosx * (y1 - y0);    // cos(x) * (y1 - y0)
  }

  Math::real DST::integralx(const vector<Math::real>& F, Math::real x) {
    return integral(F, sin(x), cos(x));
  }

} // namespace GeographicLib

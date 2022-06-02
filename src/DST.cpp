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

  DST::DST(int N)
    : _N(N < 0 ? 0 : N)
    , _fft(make_shared<fft_t>(fft_t(2 * _N, false)))
    , _data(4*_N)
    , _temp(4*_N)
  {}

  void DST::reset(int N) {
    N = N < 0 ? 0 : N;
    if (N == _N) return;
    _N = N;
    _fft->assign(2 * _N, false);
    _data.resize(4*_N);
    _temp.resize(4*_N);
  }

  void DST::fft_transform(real F[], bool centerp) const {
    // Implement DST-III (centerp = false) or DST-IV (centerp = true).

    // Elements (0,N], resp. [0,N), of _data should be set on input for centerp
    // = false, resp. true.  F must have a size of at least N and on output
    // elements [0,N) of F contain the transform.
    if (_N == 0) return;
    if (centerp) {
      for (int i = 0; i < _N; ++i) {
        _data[_N+i] = _data[_N-1-i];
        _data[2*_N+i] = -_data[i];
        _data[3*_N+i] = -_data[_N-1-i];
      }
    } else {
      _data[0] = 0;            // set [0]
      for (int i = 1; i < _N; ++i) _data[_N+i] = _data[_N-i]; // set [N+1,2*N-1]
      for (int i = 0; i < 2*_N; ++i) _data[2*_N+i] = -_data[i]; // [2*N, 4*N-1]
    }
    complex<real>* ctemp = reinterpret_cast<complex<real>*>(_temp.data());
    _fft->transform_real(_data.data(), ctemp);
    if (centerp) {
      real d = -Math::pi()/(4*_N);
      for (int i = 0, j = 1; i < _N; ++i, j+=2)
        ctemp[j] *= exp(complex<real>(0, j*d));
    }
    for (int i = 0, j = 1; i < _N; ++i, j+=2) {
      F[i] = -ctemp[j].imag() / (2*_N);
    }
  }

  void DST::fft_transform2(real F[]) const {
    // Elements [0,N), of _data should be set to the N grid center values and F
    // should have size of at least 2*N.  On input elements [0,N) of F contain
    // the size N transform; on output elements [0,2*N) of F contain the size
    // 2*N transform.

    fft_transform(F+_N, true);
    for (int i = 0; i < _N; ++i) _data[i] = F[i+_N];
    for (int i = _N; i < 2*_N; ++i)
      F[i] = (-_data[2*_N-1-i] + F[2*_N-1-i])/2;
    for (int i = 0; i < _N; ++i)
      F[i] = (_data[i] + F[i])/2;
  }

  void DST::transform(function<real(real)> f, real F[]) const {
    real d = Math::pi()/(2 * _N);
    for (int i = 1; i <= _N; ++i)
      _data[i] = f( i * d );
    fft_transform(F, false);
  }

  void DST::refine(function<real(real)> f, real F[]) const {
    real d = Math::pi()/(4 * _N);
    for (int i = 0; i < _N; ++i)
      _data[i] = f( (2*i + 1) * d );
    fft_transform2(F);
  }

  Math::real DST::eval(real sinx, real cosx, const real F[], int N) {
    // Evaluate
    // y = sum(F[i] * sin((2*i+1) * x), i, 0, N-1)
    // using Clenshaw summation.
    // Approx operation count = (N + 5) mult and (2 * N + 2) add
    real
      ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
      y0 = N & 1 ? F[--N] : 0, y1 = 0;          // accumulators for sum
    // Now N is even
    while (N > 0) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + F[--N];
      y0 = ar * y1 - y0 + F[--N];
    }
    return sinx * (y0 + y1);    // sin(x) * (y0 + y1)
  }

  Math::real DST::integral(real sinx, real cosx, const real F[], int N) {
    // Evaluate
    // y = -sum(F[i]/(2*i+1) * cos((2*i+1) * x), i, 0, N-1)
    // using Clenshaw summation.
    // Approx operation count = (N + 5) mult and (2 * N + 2) add
    int l = N;
    real
      ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
      y0 = N & 1 ? F[--N]/(2*(--l)+1) : 0, y1 = 0; // accumulators for sum
    // Now N is even
    while (N > 0) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + F[--N]/(2*(--l)+1);
      y0 = ar * y1 - y0 + F[--N]/(2*(--l)+1);
    }
    return cosx * (y1 - y0);    // cos(x) * (y1 - y0)
  }

} // namespace GeographicLib

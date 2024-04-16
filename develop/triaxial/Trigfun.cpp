/**
 * \file Trigfun.cpp
 * \brief Implementation for GeographicLib::Trigfun class
 *
 * Copyright (c) Charles Karney (2022) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "Trigfun.hpp"
#include <iostream>

#include "kissfft.hh"

namespace GeographicLib {

  using namespace std;

  Trigfun Trigfun::initbycoeffs(std::vector<real> C,
                                bool odd, bool sym,
                                real halfp) {
    return Trigfun(int(C.size()) - (!odd && !sym ? 1 : 0),
                   odd, sym, C, halfp);
  }

  Trigfun Trigfun::initbysamples(std::vector<real> F,
                                 bool odd, bool sym, bool centerp,
                                 Math::real halfp) {
    typedef kissfft<real> fft_t;
    int N = int(F.size()) - (!odd && !sym && !centerp ? 1 : 0),
      M = N * (sym ? 4 : 2);    // The size of the sample array over a period
    vector<real> H(M, Math::NaN());
    if (!centerp) {
      if (odd) H[0] = 0;
      // real slope = (odd & !sym) ? F[N-1] / N : 0;
      for (int i = 0; i < N; ++i)
        // H[i + (odd ? 1 : 0)] = F[i] - slope * i;
        H[i + (odd ? 1 : 0)] = F[i];
      if (!odd) {
        H[N] = sym ? 0 : F[N];
      }
      // Now H[0:N] is populated
      if (sym) {
        for (int i = 0; i < N; ++i)
          H[2*N - i] = (odd ? 1 : -1) * H[i];
      }
      // Now H[0:M/2] is populated
      for (int i =  1; i < M/2; ++i)
        H[M - i] = (odd ? -1 : 1) * H[i];
      // Now H[0:M-1] is populated
    } else {
      for (int i = 0; i < N; ++i)
        H[i] = F[i];
      // Now H[0:N-1] is populated
      if (sym) {
        for (int i = 0; i < N; ++i)
          H[2*N - i - 1] = (odd ? 1 : -1) * H[i];
      }
      // Now H[0:M/2-1] is populated
      for (int i =  0; i < M/2; ++i)
        H[M - i - 1] = (odd ? -1 : 1) * H[i];
      // Now H[0:M-1] is populated
    }
    /*
    cout << "F\n";
    for (int i = 0; i < M; ++i)
      cout << i << " " << H[i] << "\n";
    */
    fft_t fft(M/2, false);
    // Leave an extra slot
    vector<complex<real>> cF(M/2 + 1);
    fft.transform_real(H.data(), cF.data());
    cF[M/2] = cF[0].imag(); cF[0] = cF[0].real();
    if (centerp) {
      for (int i = 1; i <= M/2; ++i)
        cF[i] *= exp(complex<real>(0, i * (-Math::pi() / M)));
    }
    if (!sym) {
      H.resize(N+1);
      if (!odd) {
        for (int i = 0; i <= N; ++i)
          H[i] = cF[i].real() / N;
        H[0] /= 2;
        H[N] = centerp ? 0 : H[N]/2;
        /*
          cout << "H\n";
          for (int i = 0; i <= N; ++i)
          cout << i << " " << H[i] << "\n";
        */
      } else {
        for (int i = 0; i <= N; ++i)
          H[i] = -cF[i].imag() / N;
        H[0] = 0;
        H[N] = !centerp ? 0 : H[N]/2;
        }
        /*
          cout << "cF\n";
          for (int i = 0; i <= M/2; ++i)
          cout << cF[i] << "\n";
          cout << "H\n";
          for (int i = 0; i <= N; ++i)
          cout << i << " " << H[i] << "\n";
        */
    } else {                    // sym
      H.resize(N);
      if (!odd) {
        for (int i = 0; i < N; ++i)
          H[i] = cF[2*i+1].real() / (2*N);
        /*
        cout << "H\n";
        for (int i = 0; i < N; ++i)
          cout << i << " " << H[i] << "\n";
        */
      } else {
        for (int i = 0; i < N; ++i)
          H[i] = -cF[2*i+1].imag() / (2*N);
        /*
        cout << "H\n";
        for (int i = 0; i < N; ++i)
          cout << i << " " << H[i] << "\n";
        */
      }
    }
    //    if (centerp) cout << "SIZE " << F.size() << " " << H.size() << "\n";
    Trigfun t(N, odd, sym, H, halfp);
    real err = t.check(F, centerp);
    cout << err << "\n";
    return t;
  }

  Math::real Trigfun::check(const vector<real>& F, bool centerp) const {
    real err = 0, maxval = 0;
    real d = (_sym ? _h/2 : _h) / _nN;
    for (int i = 0; i < (centerp ? _nN : _nN + 1); ++i) {
      real a = centerp ? F[i] :
        (_odd ? (i == 0 ? 0 : F[i-1]) :
         (_sym && i == _nN ? 0 : F[i])),
        x = d * i + (centerp ? d/2 : 0),
        b = eval(x);
      maxval = fmax(maxval, a);
      err = err + fabs(a - b);
    }
    return err / (numeric_limits<real>::epsilon() *
                  maxval * (centerp ? _nN : _nN + 1));
  }
  void Trigfun::refine(const Trigfun& tb, const Trigfun& tref) {
    bool debug = false;
    if (debug) {
      cout << _C.size() << " " << tb._C.size() << " " << tref._C.size() << "\n";
      for (int i = 0; i < _nN + (_sym ? 0 : 1); ++i)
        cout << i << " " << _C[i] << " " << tb._C[i] << "\n";
      for (int i = 0; i < 2*_nN + (_sym ? 0 : 1); ++i)
        cout << i << " " << tref._C[i] << "\n";
    }
    _C.resize(2 * _nN + (_sym ? 0 : 1));
    for (int i = 0; i < _nN; ++i)
      _C[2*_nN + (_sym ? 0 : 1) - 1 - i] =
        (_odd ? -1 : 1) * (_C[i] - tb._C[i])/2;
    if (_odd && !_sym) _C[_nN] = tb._C[_nN];
    for (int i = 0; i < _nN; ++i)
      _C[i] = (_C[i] + tb._C[i])/2;
    if (debug) {
    for (int i = 0; i < 2*_nN + (_sym ? 0 : 1); ++i)
      cout << i << " " << _C[i] << " " << tref._C[i] << " "
           << _C[i] - tref._C[i] << "\n";
    }
    _nN *= 2;
  }

  Math::real Trigfun::eval(real z) const {
    // Evaluate
    // y = sum(c[k] * sin((k+1/2) * pi/q * z), k, 0, N - 1) if  odd && sym
    // y = sum(c[k] * cos((k+1/2) * pi/q * z), k, 0, N - 1) if !odd && sym
    // y = c[0] * pi/h * z +
    //     sum(c[k] * sin(k * pi/h * z), k, 1, N) if odd && !sym
    // y = c[0] +
    //     sum(c[k] * cos(k * pi/h * z), k, 1, N) if !odd && !sym
    real y = Math::pi()/(_sym ? _q : _h) * z;
    int k = !_sym ? _nN+1 : _nN, k0 = !_sym ? 1 : 0;
    // cout <<"C " << _C[8] << " " << _C.size() << " " << k << " " << k0 << "\n";
    real u0 = 0, u1 = 0,        // accumulators for sum
      x = 2 * cos(y);
    for (; k > k0;) {
      real t = x * u0 - u1 + _C[--k];
      u1 = u0; u0 = t;
    }
    // sym
    //   y = u0*f0(zeta) - u1*fm1(zeta)
    //   f0 = odd ? sin(y/2) : cos(y/2)
    //   fm1 = odd ? -sin(y/2) : cos(y/2)
    //   y = odd ? sin(y/2) * (u0 + u1) : cos(y/2) * (u0 - u1)
    // !sym
    //   y = u0*f1(zeta) - u1*f0(zeta)
    //   f1 = odd ? sin(y) : cos(y)
    //   f0 = odd ? 0 : 1
    return _sym ? (_odd ? sin(y/2) * (u0 + u1) : cos(y/2) * (u0 - u1)) :
      _C[0] * (_odd ? y : 1) + (_odd ? sin(y) : x/2) * u0 - (_odd ? 0 : u1);
  }
#if 0
  Math::real Trigfun::eval(real x);
  Math::real Trigfun::eval(int i);
  std::vector<Math::real> Trigfun::Coeffs();
  std::vector<Math::real> Trigfun::Samples();
#endif

} // namespace GeographicLib

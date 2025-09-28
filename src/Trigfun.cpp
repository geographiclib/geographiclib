/**
 * \file Trigfun.cpp
 * \brief Implementation for GeographicLib::Trigfun class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

/// \cond SKIP
// For to_string
#include <string>
#include <iostream>
#include <iomanip>
#include <GeographicLib/Trigfun.hpp>
#include "kissfft.hh"

#define USE_ANGLE 0

namespace GeographicLib {

  using namespace std;

  Trigfun::Trigfun(int n, const std::function<real(real)>& f,
                   bool odd, bool sym, real halfp, bool centerp) {
    if (n > 0) {
      int M = n + (!(odd || sym || centerp) ? 1 : 0);
      real p = halfp / (sym ? 2 : 1), d = p / n,
        o = centerp ? d/2 : ( odd ? d : 0 );
      vector<real> F(M);
      for (int i = 0; i < M; ++i)
        F[i] = f(o + d * i);
      *this = initbysamples(F, odd, sym, halfp, centerp);
    } else
      *this = Trigfun{vector<real>{1,0}, odd,sym, halfp};
  }

  Trigfun::Trigfun(const function<real(real)>& f, bool odd, bool sym,
                   real halfp, int nmax, real tol, real scale) {
    int n = 16;
    Trigfun t(n, f, odd, sym, halfp, false);
    while (n <= nmax) {
      int K = chop(t._coeff, tol, scale);
      if (K < n) {
        t._m = K;
        t._n = t._sym ? K : K - 1;
        t._coeff.resize(K);
        *this = t;
        return;
      }
      Trigfun tx(n, f, odd, sym, halfp, true);
      t.refine(tx);
      n *= 2;
    }
    *this = t;
  }

  Trigfun::Trigfun(const function<real(real, real)>& f, bool odd, bool sym,
                   real halfp, int nmax, real tol, real scale) {
    // Initialize with 2 samples
    Trigfun t(2,
              [&f] (real x) -> real
              { return f(x, Math::NaN()); },
              odd, sym, halfp, false);
    while (t._n <= nmax) {
      int K = chop(t._coeff, tol, scale);
      if (K < t._n) {
        t._m = K;
        t._n = t._sym ? K : K - 1;
        t._coeff.resize(K);
        *this = t;
        return;
      }
      // bool centerp = true;
      int M = t._n;
      real p = halfp / (sym ? 2 : 1), d = p / t._n, o = d/2;
      vector<real> F(M);
      for (int i = 0; i < M; ++i)
        F[i] = f(o + d * i, t(o + d * i));
      t.refine(initbysamples(F, odd, sym, halfp, true));
      // n *= 2;
    }
    *this = t;
    return;
  }

  Trigfun Trigfun::initbysamples(const vector<real>& F,
                                 bool odd, bool sym, real halfp, bool centerp) {
    if (!(isfinite(halfp) && halfp > 0))
      throw GeographicErr("Trigfun::initbysamples halfp not positive");
    using fft_t = kissfft<real>;
    int n = int(F.size()) - (!(odd || sym || centerp) ? 1 : 0),
      M = n * (sym ? 4 : 2);    // The size of the sample array over a period
    vector<real> H(M, Math::NaN());
    if (!centerp) {
      if (odd) H[0] = 0;
      // real slope = (odd & !sym) ? F[n-1] / n : 0;
      for (int i = 0; i < n; ++i)
        // H[i + (odd ? 1 : 0)] = F[i] - slope * i;
        H[i + (odd ? 1 : 0)] = F[i];
      if (!odd) {
        H[n] = sym ? 0 : F[n];
      }
      // Now H[0:n] is populated
      if (sym) {
        for (int i = 0; i < n; ++i)
          H[2*n - i] = (odd ? 1 : -1) * H[i];
      }
      // Now H[0:M/2] is populated
      for (int i =  1; i < M/2; ++i)
        H[M - i] = (odd ? -1 : 1) * H[i];
      // Now H[0:M-1] is populated
    } else {
      for (int i = 0; i < n; ++i)
        H[i] = F[i];
      // Now H[0:n-1] is populated
      if (sym) {
        for (int i = 0; i < n; ++i)
          H[2*n - i - 1] = (odd ? 1 : -1) * H[i];
      }
      // Now H[0:M/2-1] is populated
      for (int i =  0; i < M/2; ++i)
        H[M - i - 1] = (odd ? -1 : 1) * H[i];
      // Now H[0:M-1] is populated
    }
    //    cout << "FFT size " << M/2 << "\n";
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
      H.resize(n+1);
      if (!odd) {
        for (int i = 0; i <= n; ++i)
          H[i] = cF[i].real() / n;
        H[0] /= 2;
        H[n] = centerp ? 0 : H[n]/2;
      } else {
        for (int i = 0; i <= n; ++i)
          H[i] = -cF[i].imag() / n;
        H[0] = 0;
        H[n] = !centerp ? 0 : H[n]/2;
        }
    } else {                    // sym
      H.resize(n);
      if (!odd) {
        for (int i = 0; i < n; ++i)
          H[i] = cF[2*i+1].real() / (2*n);
      } else {
        for (int i = 0; i < n; ++i)
          H[i] = -cF[2*i+1].imag() / (2*n);
      }
    }
    //    if (centerp) cout << "SIZE " << F.size() << " " << H.size() << "\n";
    Trigfun t(H, odd, sym, halfp);
    if constexpr (debug_) {
      real err = t.check(F, centerp, 0);
      if (err > 100) {
        cout << t._n << " " << err << endl;
        throw GeographicErr("initbysamples error");
      }
    }
    return t;
  }

  Math::real Trigfun::check(const vector<real>& F, bool centerp, real tol)
    const {
    real err = 0, maxval = 0;
    real d = (_sym ? _h/2 : _h) / _n;
    for (int i = 0; i < (centerp ? _n : _n + 1); ++i) {
      real a = centerp ? F[i] :
        (_odd ? (i == 0 ? 0 : F[i-1]) :
         (_sym && i == _n ? 0 : F[i])),
        x = d * i + (centerp ? d/2 : 0),
        b = (*this)(x);
      maxval = fmax(maxval, fabs(a));
      err = err + fabs(a - b);
    }
    //    cout << "Maxval " << maxval << "\n";
    return err / (tolerance(tol) *
                  maxval * (centerp ? _n : _n + 1));
  }

  void Trigfun::refine(const Trigfun& tb) {
    int m = 2 * _n + (_sym ? 0 : 1);
    _coeff.resize(m);
    for (int i = 0; i < _n; ++i)
      _coeff[2*_n + (_sym ? 0 : 1) - 1 - i] =
        (_odd ? -1 : 1) * (_coeff[i] - tb._coeff[i])/2;
    if (_odd && !_sym) _coeff[_n] = tb._coeff[_n];
    for (int i = 0; i < _n; ++i)
      _coeff[i] = (_coeff[i] + tb._coeff[i])/2;
    _max = -1;
    _n *= 2;
    _m = m;
  }

  void Trigfun::setsecular(real f0) {
    if (!(_odd && !_sym))
      throw GeographicErr
        ("Trigfun: cannot set secular term unless odd && !sym");
    _coeff[0] = f0 / Math::pi();
  }

  Math::real Trigfun::Max() const {
    if (_max < 0) {
      _max = 0;
      for (int k = _m; k > (_sym ? 0 : 1);)
        _max += fabs(_coeff[--k]);
    }
    return _max;
  }

#if !USE_ANGLE
  Math::real Trigfun::operator()(real z) const {
    // Evaluate
    // y = sum(c[k] * sin((k+1/2) * pi/q * z), k, 0, n - 1) if  odd && sym
    // y = sum(c[k] * cos((k+1/2) * pi/q * z), k, 0, n - 1) if !odd && sym
    // y = c[0] * pi/h * z +
    //     sum(c[k] * sin(k * pi/h * z), k, 1, n) if odd && !sym
    // y = c[0] +
    //     sum(c[k] * cos(k * pi/h * z), k, 1, n) if !odd && !sym
    if (_coeff.empty()) return 0;
    real y = Math::pi()/(_sym ? _q : _h) * z;
    int k = _m,
      k0 = _odd && !_sym ? 1 : 0; // add secular term at the end
    real u0 = 0, u1 = 0,          // accumulators for sum
      x = 2 * cos(y);
    for (; k > k0;) {
      real t = x * u0 - u1 + _coeff[--k];
      u1 = u0; u0 = t;
    }
    // sym
    //   v = u0*f0(zeta) - u1*fm1(zeta)
    //   f0 = odd ? sin(y/2) : cos(y/2)
    //   fm1 = odd ? -sin(y/2) : cos(y/2)
    //   v = odd ? sin(y/2) * (u0 + u1) : cos(y/2) * (u0 - u1)
    // !sym && !odd
    //   v = u0*f0(zeta) - u1*fm1(zeta)
    //   f0 = 1
    //   fm1 = cos(y)
    //   v = u0 - cos(y) * u1 = u0 - (x/2) * u1
    // !sym && odd
    //   v = u0*f1(zeta) - u1*f0(zeta)
    //   f1 = sin(y)
    //   f0 = 0
    //   v = sin(y) * u0 + secular term
    return _sym ? (_odd ? sin(y/2) * (u0 + u1) : cos(y/2) * (u0 - u1)) :
      !_odd ? u0 - (x/2) * u1 : _coeff[0] * y + sin(y) * u0;
  }

#else
  // Implementation in terms of Angle
  Math::real Trigfun::operator()(real z) const {
    // Evaluate
    // y = sum(c[k] * sin((k+1/2) * pi/q * z), k, 0, n - 1) if  odd && sym
    // y = sum(c[k] * cos((k+1/2) * pi/q * z), k, 0, n - 1) if !odd && sym
    // y = c[0] * pi/h * z +
    //     sum(c[k] * sin(k * pi/h * z), k, 1, n) if odd && !sym
    // y = c[0] +
    //     sum(c[k] * cos(k * pi/h * z), k, 1, n) if !odd && !sym
    if (_coeff.empty()) return 0;
    real y = Math::pi()/(_sym ? _q : _h) * z;
    return (_odd && !_sym ? _coeff[0] * y : 0) + eval(Angle::radians(y/2));
  }

  Math::real Trigfun::eval(Angle phi) const {
    // Evaluate
    // y = sum(c[k] * sin((2*k+1) * phi), k, 0, n - 1) if  odd && sym
    // y = sum(c[k] * cos((2*k+1) * phi), k, 0, n - 1) if !odd && sym
    // y = c[0] * 0 +
    //     sum(c[k] * sin(2*k * phi), k, 1, n) if odd && !sym
    // y = c[0] +
    //     sum(c[k] * cos(2*k * phi), k, 1, n) if !odd && !sym
    // Note that secular term c[0] is ignored for odd && !sym.
    if (_coeff.empty()) return 0;
    int k = _m,
      k0 = _odd && !_sym ? 1 : 0; // secular term excluded
    real u0 = 0, u1 = 0,          // accumulators for sum
      x = 2 * (phi.c() - phi.s()) * (phi.c() + phi.s());
    for (; k > k0;) {
      real t = x * u0 - u1 + _coeff[--k];
      u1 = u0; u0 = t;
    }
    // sym
    //   v = u0*f0(zeta) - u1*fm1(zeta)
    //   f0 = odd ? sin(phi) : cos(phi)
    //   fm1 = odd ? -sin(phi) : cos(phi)
    //   v = odd ? sin(phi) * (u0 + u1) : cos(phi) * (u0 - u1)
    // !sym && !odd
    //   v = u0*f0(zeta) - u1*fm1(zeta)
    //   f0 = 1
    //   fm1 = cos(2*phi)
    //   v = u0 - cos(2*phi) * u1 = u2 - (x/2) * u1
    // !sym && odd
    //   v = u0*f1(zeta) - u1*f0(zeta)
    //   f1 = sin(2*phi)
    //   f0 = 0
    //   v = sin(2*phi) * u0 + no secular term
    return _sym ? (_odd ? phi.s() * (u0 + u1) : phi.c() * (u0 - u1)) :
      !_odd ? u0 - (x/2) * u1 : 2 * phi.s() * phi.c() * u0;
  }
#endif

  Trigfun Trigfun::integral() const {
    if (_odd && !_sym && _coeff[0] != 0)
      throw GeographicErr("Trigfun: cannot integrate a secular term");
    vector <real> c(_coeff);
    real mult = (_odd ? -1 : 1) * (_sym ? _q : _h) / Math::pi();
    for (int i = 0; i < _m; ++i)
      c[i] *= mult / (i + (_sym ? real(0.5) : 0));
    if (!_sym)
      c[0] = _odd ? 0 : _coeff[0] * mult;
    return Trigfun(c, !_odd, _sym, _h);
  }

  // root sig 2
  Math::real Trigfun::root(ind indicator,
                           real z, const function<real(real)>& fp,
                           real x0,
                           int* countn, int* countb, real tol) const {
    if (!(_odd && !_sym))
      throw GeographicErr
        ("Trigfun: cannot take root unless odd && !sym");
    // y = pi/h * x
    // f(x) = c[0] * y + sum(c[k] * sin(k * y), k, 1, n)
    real hr = Math::pi() * _coeff[0], s = _h / hr,
      x00 = s * z, dx = fabs(s) * Max();
    x0 = isfinite(x0) ? fmin(x00 + dx, fmax(x00 - dx, x0)) : x00;
    //    cout << "QQG " << dx << "\n";
    if (dx == 0) return x0;
    auto ffp = [this, &fp]
      (real x) -> pair<real, real>
      { return pair<real, real>(this->operator()(x), fp(x)); };
    return root(indicator,ffp, z, x0, x00 - dx, x00 + dx, _h, fabs(hr),
                s > 0 ? 1 : -1, countn, countb, tol);
  }

  // root sig 4
  Math::real Trigfun::root(ind indicator,
                           const function<pair<real, real>(real)>& ffp,
                           real z,
                           real x0, real xa, real xb,
                           real xscale, real zscale, int s,
                           int* countn, int* countb,
                           real tol) {
    // Solve v = f(x) - z = 0
    if (!(xa <= x0 && x0 <= xb)) return Math::NaN();
    if (x0 == xa && x0 == xb)
      return x0;
    tol = tolerance(tol);
    real vtol = tol * zscale/100,
      xtol = pow(tol, real(0.75)) * xscale,
      x = x0, oldx = Math::infinity(), oldv = oldx, olddx = oldx;
    int k = 0, b = 0;
    real p = Math::pi()/2 * 0;
    if constexpr (debug_) {
      cout << "SCALE " << xscale << " " << zscale << "\n";
      pair<real, real> vala = ffp(xa);
      pair<real, real> val0 = ffp(x0);
      pair<real, real> valb = ffp(xb);
      cout << "DAT " << s << " " << x0-xa << " " << xb-x0 << " " << z << "\n";
      cout << "DAT "
           << xa << " " << vala.first - z << " " << vala.second << "\n";
      cout << "DAT "
           << x0 << " " << val0.first - z << " " << val0.second << "\n";
      cout << "DAT "
           << xb << " " << valb.first - z << " " << valb.second << "\n";
      if ((vala.first - z) * (valb.first - z) > 0)
        cout << "DATBAD\n";
    }
    for (; k < maxit_ ||
           (throw_ && (throw GeographicLib::GeographicErr
                       ("Convergence failure Trigfun::root case=" +
                        to_string(indicator)), false));) {
      // TODO: This inverse problem uses lots of iterations
      //   20 60 -90 180 127.4974 24.6254 2.4377
      // Need to figure out why.  (Probably fixed by now.)
      ++k;
      pair<real, real> val = ffp(x);
      real v = val.first - z,
        vp = val.second,
        dx = - v/vp;
      if constexpr (debug_)
        cout << "XX " << k << " " << xa-p << " " << x-p << " " << xb-p << " "
             << dx << " " << x + dx-p << " " << v << " " << vp << endl;
      if (!(fabs(v) > (k < 2 ? 0 : vtol))) {
        if constexpr (debug_) cout << "break1 " << k << " " << fabs(v) << endl;
        break;
      } else if (s*v > 0)
        xb = fmin(xb, x);
      else
        xa = fmax(xa, x);
      x += dx;
      if (!(xa <= x && x <= xb) || fabs(v) > oldv ||
          (k > 2 && 2 * fabs(dx) > olddx)) {
        if constexpr (debug_)
          cout << "bis " << k << " " << xa-x << " " << x-xb << " ";
        x = (xa + xb)/2;
        ++b;
        if (x == oldx) {
          if constexpr (debug_)
            cout << "break3 " << k << " " << x << " " << dx << "\n";
          break;
        }
      } else if (!(fabs(dx) > xtol)) {
        if constexpr (debug_)
          cout << "break2 " << k << " " << dx << " " << xtol << endl;
        break;
      }
      if constexpr (debug_)
        cout << "GAPS " << k << " " << dx << " " <<  x-xa << " " << xb-x << " "
             << oldx << " " << x << " " << (oldx - x) << "\n";
      oldx = x;
      oldv = fabs(v);
      olddx = fabs(dx);
    }
    if (countn) *countn += k;
    if (countb) *countb += b;
    if constexpr (debug_)
      cout << "return " << x << "\n";
    return x;
  }

  Math::real Trigfun::inversep(real z,
                               const function<real(real)>& fp,
                               real dx0,
                               int* countn, int* countb, real tol) const {
    real hr = Math::pi() * _coeff[0], nslope = _h / hr;
    return root(INVERSEP, z, fp, z * nslope + dx0, countn, countb, tol) -
      nslope * z;
  }

  Trigfun Trigfun::inverse(const function<real(real)>& fp,
                           int* countn, int* countb,
                           int nmax, real tol, real scale) const {
    if (!(_odd && !_sym && isfinite(_coeff[0]) && _coeff[0] != 0))
      throw GeographicErr("Can only invert Trigfun with a secular term");
    int s = _coeff[0] > 0 ? 1 : -1;
    real hp = _h, hr = Math::pi() * _coeff[0],
      nhp = hr * s, nhr = hp * s, c0p = nhr / Math::pi();
     Trigfun t(
               [this, &fp, countn, countb, tol]
               (real z, real dx0) -> real
               { return inversep(z, fp, dx0, countn, countb, tol); },
               _odd, _sym, nhp, nmax, tol, scale);
     t._coeff[0] = c0p;
     return t;
  }

  int Trigfun::chop(const vector<real>& c, real tol, real scale) {
    // This is a clone of Chebfun's standardChop function.  For C++, the return
    // value is number of terms to retain.  Index of last term is one less than
    // this.
    //
    // See J. L. Aurentz and L. N. Trefethen, "Chopping a Chebyshev series",
    // https://doi.org/10.1145/2998442 (2017) and
    // https://arxiv.org/abs/1512.01803 (2015).
    //
    // Input:
    //
    // COEFFS  A nonempty row or column vector of real or complex numbers
    //         which typically will be Chebyshev or Fourier coefficients.
    //
    // TOL     A number in (0,1) representing a target relative accuracy.
    //         TOL will typically will be set to the Chebfun EPS parameter,
    //         sometimes multiplied by a factor such as vglobal/vlocal in
    //         construction of local pieces of global chebfuns.
    //         Default value: machine epsilon (MATLAB EPS).
    //
    // Output:
    //
    // CUTOFF  A positive integer.
    //         If CUTOFF == length(COEFFS), then we are "not happy":
    //         a satisfactory chopping point has not been found.
    //         If CUTOFF < length(COEFFS), we are "happy" and CUTOFF
    //         represents the last index of COEFFS that should be retained.
    //
    // Examples:
    //
    // coeffs = 10.^-(1:50); random = cos((1:50).^2);
    // standardChop(coeffs) // = 18
    // standardChop(coeffs + 1e-16*random) // = 15
    // standardChop(coeffs + 1e-13*random) // = 13
    // standardChop(coeffs + 1e-10*random) // = 50
    // standardChop(coeffs + 1e-10*random, 1e-10) // = 10

    // Jared Aurentz and Nick Trefethen, July 2015.
    //
    // Copyright 2017 by The University of Oxford and The Chebfun Developers.
    // See http://www.chebfun.org/ for Chebfun information.

    // STANDARDCHOP normally chops COEFFS at a point beyond which it is smaller
    // than TOL^(2/3).  COEFFS will never be chopped unless it is of length at
    // least 17 and falls at least below TOL^(1/3).  It will always be chopped
    // if it has a long enough final segment below TOL, and the final entry
    // COEFFS(CUTOFF) will never be smaller than TOL^(7/6).  All these
    // statements are relative to MAX(ABS(COEFFS)) and assume CUTOFF > 1.
    // These parameters result from extensive experimentation involving
    // functions such as those presented in the paper cited above.  They are
    // not derived from first principles and there is no claim that they are
    // optimal.

    // Check magnitude of TOL:
    tol = tolerance(tol);
    if (tol >= 1) return 1;

    // Make sure c has length at least 17:
    int n = int(c.size());
    // Change 17 in original code to 16 to accommodate trig expansions which
    // may only have 2^n terms.
    if (n < 16) return n;

    // Step 1: Convert c to a new monotonically nonincreasing
    //         vector ENVELOPE normalized to begin with the value 1.

    vector<real> m(n);
    int j = n;
    m[--j] = fabs(c[n - 1]);
    for (; j;) {
      --j;
      m[j] = fmax(fabs(c[j]), m[j + 1]);
    }
    if (m[0] == 0) return 1;
    if (scale >= 0) m[0] = fmax(scale, m[0]);
    for (j = n; j;)
      m[--j] /= m[0];

    // Step 2: Scan ENVELOPE for a value PLATEAUPOINT, the first point J-1, if
    // any, that is followed by a plateau.  A plateau is a stretch of
    // coefficients ENVELOPE(J),...,ENVELOPE(J2), J2 = round(1.25*J+5) <= n,
    // with the property that ENVELOPE(J2)/ENVELOPE(J) > R.  The number R
    // ranges from R = 0 if ENVELOPE(J) = TOL up to R = 1 if ENVELOPE(J) =
    // TOL^(2/3).  Thus a potential plateau whose starting value is ENVELOPE(J)
    // ~ TOL^(2/3) has to be perfectly flat to count, whereas with ENVELOPE(J)
    // ~ TOL it doesn't have to be flat at all.  If a plateau point is found,
    // then we know we are going to chop the vector, but the precise chopping
    // point CUTOFF still remains to be determined in Step 3.

    int j2 = 0, plateauPoint = n;
    real logtol = log(tol);
    for (j = 2; j <= n; ++j) {  // j is a MATLAB index (starts at 1)
      j2 = int(round(1.25*j + 5));
      if (j2 > n) return n;
      real e1 = m[j-1],
        e2 = m[j2-1],
        r = 3 * (1 - log(e1)/logtol);
      if ( e1 == 0 || e2/e1 > r ) {
        // a plateau has been found: go to Step 3
        plateauPoint = j - 1;
        break;
      }
    }

    // Step 3: fix CUTOFF at a point where ENVELOPE, plus a linear function
    // included to bias the result towards the left end, is minimal.
    //
    // Some explanation is needed here.  One might imagine that if a plateau is
    // found, then one should simply set CUTOFF = PLATEAUPOINT and be done,
    // without the need for a Step 3. However, sometimes CUTOFF should be
    // smaller or larger than PLATEAUPOINT, and that is what Step 3 achieves.
    //
    // CUTOFF should be smaller than PLATEAUPOINT if the last few coefficients
    // made negligible improvement but just managed to bring the vector
    // ENVELOPE below the level TOL^(2/3), above which no plateau will ever be
    // detected.  This part of the code is important for avoiding situations
    // where a coefficient vector is chopped at a point that looks "obviously
    // wrong" with PLOTCOEFFS.
    //
    // CUTOFF should be larger than PLATEAUPOINT if, although a plateau has
    // been found, one can nevertheless reduce the amplitude of the
    // coefficients a good deal further by taking more of them.  This will
    // happen most often when a plateau is detected at an amplitude close to
    // TOL, because in this case, the "plateau" need not be very flat.  This
    // part of the code is important to getting an extra digit or two beyond
    // the minimal prescribed accuracy when it is easy to do so.

    if ( m[plateauPoint - 1] == 0 ) return plateauPoint;
    real tol76 = tol * sqrt(cbrt(tol)); // tol^(7/6)
    int j3 = 0;
    for (j = 0; j < n; ++j) {
      if (m[j] >= tol76) ++j3;
    }
    if ( j3 < j2 ) {
      j2 = j3 + 1;
      m[j2 - 1] = tol76;
    }
    vector<real> cc(j2);
    // Replace log10 by log.  This involved no change in the logic.
    real tol3 = logtol/(3 * (j2 - 1));
    int d = 0;
    for (j = 0; j < j2; ++j) {
      cc[j] = log(m[j]) - tol3 * j;
      if (j > 0 && cc[j] < cc[d]) d = j;
    }
    return max(d, 1);
  }

  TrigfunExt::TrigfunExt(const function<real(real)>& fp, real halfp,
                         bool sym, real scale)
      : _fp(fp)
      , _sym(sym)
        // N.B. tol defaults to epsilon() here.  We need to compute the
        // integral accurately.
      , _f(Trigfun(_fp, false, _sym, halfp, 1 << 16, 0, scale).integral())
      , _tol(sqrt(numeric_limits<real>::epsilon()))
      , _nmax(int(ceil(real(1.5) * _f.NCoeffs())))
      , _invp(false)
    {}

} // namespace GeographicLib
/// \endcond

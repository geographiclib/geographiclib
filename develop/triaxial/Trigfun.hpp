/**
 * \file Trigfun.hpp
 * \brief Header for GeographicLib::Trigfun class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIGFUN_HPP)
#define GEOGRAPHICLIB_TRIGFUN_HPP 1

#include <GeographicLib/Constants.hpp>

#include <functional>
#include <utility>
#include <vector>
#include <iostream>

namespace GeographicLib {

  /**
   * \brief Representing a function by a Fourier series
   *
   * This class mimic the functionality of Chebfun's 'trig' representation of
   * periodic functions.  Key differences are:
   *
   * - Only odd or even functions are allowed (i.e., only sine of only cosine
   *   terms in the Fourier series).
   * - Can specify that the function has symmetry about the quarter-period
   *   point so that the Fouries series only includes odd terms.
   * - The integral of a trigfun is counted as a trigfun even if it includes a
   *   secular term.
   * - The inverse function of the integral is also a trigfun (only makes sense
   *   if the orginal function is either strictly positive or strictly
   *   negative.
   *
   * Here we compute FFTs using the kissfft package
   * https://github.com/mborgerding/kissfft by Mark Borgerding.
   *
   * Example of use:
   * XX include example-Trigfun.cpp
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT Trigfun {
  private:
    typedef Math::real real;
    int _m,                     // Number of coefficients in series
      _n;                       // Number of samples in half/quarter period
    bool _odd, _sym;
    std::vector<real> _coeff;
    real _h, _q, _p;            // half, quarter, whole period
    mutable real _max;
    static int chop(const std::vector<real>& c, real tol, real scale = -1);
    // Function samples over half/quarter period of !sym/sym
    // odd sym cent  samples            nF  nC
    //  f   f   f    |-|-|-|-|-|-|-|-|  n+1 n+1 (4)
    //  t   f   f    --|-|-|-|-|-|-|-|  n   n+1 (1), (2), (3), (7)
    //  f   t   f    |-|-|-|-|-|-|-|--  n   n   (1)
    //  t   t   f    --|-|-|-|-|-|-|-|  n   n   (1)
    //  f   f   t    -|-|-|-|-|-|-|-|-  n   n+1 (4), (6)
    //  t   f   t    -|-|-|-|-|-|-|-|-  n   n+1 (5)
    //  f   t   t    -|-|-|-|-|-|-|-|-  n   n
    //  t   t   t    -|-|-|-|-|-|-|-|-  n   n
    //
    // (1) missing end terms presumed zero
    // (2) included last term is usually zero, if non zero, gives secular term
    // (3) zeroth coeff used for secular term
    // (4) zeroth coeff gives constant.
    // (5) secular term should have been removed from samples
    // (6) last coeff is zero (but not for centerp)
    // (7) last coeff is zero (but not for !centerp)
    // Function is represented by (y = pi/h * x)
    // sym = false, sample in f_i = f(h * i/n)
    // odd = true (n samples, n+1 coeffs)
    //    f_0 = 0, need f_i for i in (0, n], f_n defines linear contrib
    // f(x) = c[0] * y + sum(c[k] * sin(k * y), k, 1, n)
    // F(x) = 0 + (-h/pi) * sum( c[k]/k * cos(k * y), k, 1, n)
    // odd = false (n+1 samples, n+1 coeffs)
    //   need f_i for i in [0, n]
    // f(x) = c[0] + sum(c[k] * cos(k * y), k, 1, n)
    // F(x) = (h/pi) * (c[0] * y + sum( c[k]/k * sin(k * y), k, 1, n))
    // sym = true, sample in f_i = f(q * i/n)
    // odd = true (n samples, n coeffs)
    //   f_0 = 0, need f_i for i in (0, n] (n samples)
    // f(x) = sum(c[k] * sin(2*(k+1/2) * y), k, 0, n - 1)
    // F(x) = -(q/pi) * sum(c[k]/(k+1/2) * cos(2*(k+1/2) * y), k, 0, n - 1)
    // odd = false (n samples, n coeffs)
    //   f_n = 0, need f_i for i in [0, n) (n samples)
    // f(x) = sum(c[k] * cos(2*(k+1/2) * y), k, 0, n - 1)
    // F(x) = (q/pi) * sum(c[k]/(k+1/2) * sin(2*(k+1/2) * y), k, 0, n - 1)

    Trigfun(const std::vector<Math::real>& c, bool odd, bool sym, real h)
      : _m(int(c.size()))
      , _n(sym ? _m : _m - 1)
      , _odd(odd)
      , _sym(sym)
      , _coeff(c)
      , _h(h)
      , _q(_h/2)
      , _p(2*_h)
      , _max(-1)
    {
      (void) _p;
    }
    /**
     * Constructor given a function.  Specify n = 0 to do auto
     **********************************************************************/
    Trigfun(const std::function<real(real)>& f, bool odd, bool sym,
            bool centerp, real halfp, int n, int nmax = 1 << 16,
            real tol = std::numeric_limits<real>::epsilon(),
            real scale = -1);

  public:
    /**
     * Default constructor specifying the number of points to use.
     **********************************************************************/
    Trigfun() : _m(0), _max(-1) {}
    Trigfun(const std::function<real(real)>& f, bool odd, bool sym,
            real halfp, int nmax = 1 << 16,
            real tol = std::numeric_limits<real>::epsilon(),
            real scale = -1);
    Trigfun(const std::function<real(real, real)>& f, bool odd, bool sym,
            real halfp, int nmax = 1 << 16,
            real tol = std::numeric_limits<real>::epsilon(),
            real scale = -1);
    real check(const std::vector<real>& F, bool centerp,
               real tol = std::numeric_limits<real>::epsilon()) const;
    //    real eval(real x) const;
    real operator()(real x) const;
    static Trigfun initbysamples(const std::vector<real>& F,
                                 bool odd, bool sym, bool centerp, real halfp);
    void refine(const Trigfun& tb);
    Trigfun integral() const;
    enum ind {
      NONE = 0,
      INV1,
      ARCPOS0,
      FFUNROOT,
      GFUNROOT,
      INVERSEP,
      OTHER,
    };

    // Given z, find x, s.t. z = f(x); fp is the derivative f'.  This defines x
    // = finv(z).  x0 is an estimate of x (NaN means no information)
    // root sig 1
    real root(real z, const std::function<real(real)>& fp,
              int* countn = nullptr, int* countb = nullptr,
              real tol = 0, ind indicator = NONE) const;
    // root sig 2
    real root(real z, const std::function<real(real)>& fp,
              real x0,
              int* countn = nullptr, int* countb = nullptr,
              real tol = 0, ind indicator = NONE) const;
    // Solve f(x) = z for x, given x in [xa, xb];
    // fp(x) = df(x)/dx
    // s is sign of fp
    // xscale and zscale of scales for x and f(x).
    // root sig 3
    static real root(const std::function<real(real)>& f,
                     real z, const std::function<real(real)>& fp,
                     real x0, real xa, real xb,
                     real xscale = 1, real zscale = 1, int s = 1,
                     int* countn = nullptr, int* countb = nullptr,
                     real tol = 0,
                     ind indicator = NONE);
    // ffp returns a pair [f(x), fp(x)]
    // root sig 4
    static real root(const std::function<std::pair<real, real>(real)>& ffp,
                     real z,
                     real x0, real xa, real xb,
                     real xscale = 1, real zscale = 1, int s = 1,
                     int* countn = nullptr, int* countb = nullptr,
                     real tol = 0, ind indicator = NONE);
    // Given z, return dx = finv(z) - nslope * z
    // dx0 is an estimate of dx (NaN means no information)
    real inversep(real z, const std::function<Math::real(Math::real)>& fp,
                  real dx0 = Math::NaN(),
                  int* countn = nullptr, int* countb = nullptr,
                  real tol = 0) const;
    Trigfun invert(const std::function<real(real)>& fp,
                   int* countn = nullptr, int* countb = nullptr,
                   int nmax = 1 << 16, real tol = 0, real scale = -1) const;
    int NCoeffs() const { return _m; }
    real Max() const;
    real HalfPeriod() const { return _h; }
    real HalfRange() const {
      return _odd && !_sym ? _coeff[0] * Math::pi() : Max();
    }
    real Slope() const {
      return _odd && !_sym ? HalfRange() / HalfPeriod() : 0;
    }
#if 0
    /**
     * Evaluate the Fourier sum given the sine and cosine of the angle
     *
     * @param[in] sinx sin&sigma;.
     * @param[in] cosx cos&sigma;.
     * @param[in] F the array of Fourier coefficients.
     * @param[in] n the number of Fourier coefficients.
     * @return the value of the Fourier sum.
     **********************************************************************/
    real eval(int i);
    std::vector<real> Coeffs();
    std::vector<real> Samples();
#endif
  };

  /**
   * \brief A function defined by its derivative and its inverse
   *
   * This builds on the Trigfun class.
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT TrigfunExt {
  private:
    typedef Math::real real;
    std::function<real(real)> _fp;
    bool _sym;
    Trigfun _f;
    real _tol;
    int _nmax;
    bool _invp;
    Trigfun _finv;
    int _countn, _countb;
  public:
    TrigfunExt() {}
    /**
     * Constructor specifying the derivative, an even periodic function
     **********************************************************************/
    TrigfunExt(const std::function<real(real)>& fp, real halfp,
               bool sym = false, real scale = -1);
    real operator()(real x) const { return _f(x); }
    real deriv(real x) const { return _fp(x); }
    void ComputeInverse() {
      if (!_invp && !_sym) {
        _countn = _countb = 0;
        _finv = _f.invert(_fp, &_countn, &_countb, _nmax, _tol);
        _invp = true;
      }
    }
    // Approximate inverse using _finv
    real inv0(real z) const {
      if (!_invp) return Math::NaN();
      return _sym ? Math::NaN() : _finv(z);
    }
    // Accurate inverse by direct Newton (not using _finv)
    real inv1(real z, int* countn = nullptr, int* countb = nullptr) const {
      return _sym ? Math::NaN() : _f.root(z, _fp, countn, countb, 0,
                                          Trigfun::INV1);
    }
    // Accurate inverse correcting result from _finv
    real inv2(real z, int* countn = nullptr, int* countb = nullptr) const {
      if (!_invp) return Math::NaN();
      return _sym ? Math::NaN() : _f.root(z, _fp, _finv(z), countn, countb);
    }
    real inv(real z, int* countn = nullptr, int* countb = nullptr) const {
      return _invp ? inv2(z, countn, countb) : inv1(z, countn, countb);
    }
    int NCoeffs() const { return _f.NCoeffs(); }
    int NCoeffsInv() const {
      if (!_invp) return -1;
      return _finv.NCoeffs(); }
    std::pair<int, int> InvCounts() const {
      if (!_invp) return std::pair<int, int>(-1, -1);
      return std::pair<int, int>(_countn, _countb);
    }
    real Max() const { return _f.Max(); }
    real HalfPeriod() const { return _f.HalfPeriod(); }
    real HalfRange() const { return _f.HalfRange(); }
    real Slope() const { return _f.Slope(); }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIGFUN_HPP

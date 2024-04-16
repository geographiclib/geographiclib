/**
 * \file Trigfun.hpp
 * \brief Header for GeographicLib::Trigfun class
 *
 * Copyright (c) Charles Karney (2022-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIGFUN_HPP)
#define GEOGRAPHICLIB_TRIGFUN_HPP 1

#include <GeographicLib/Constants.hpp>

#include <functional>
#include <memory>
#include <vector>

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
    int _nN;                    // Number of samples in half/quarter period
    bool _odd, _sym;
    std::vector<real> _C;
    real _q, _h, _p;              // quarter, half, whole period
    // Function samples over half/quarter period of !sym/sym
    // odd sym cent  samples            nF  nC
    //  f   f   f    |-|-|-|-|-|-|-|-|  N+1 N+1 (4)
    //  t   f   f    --|-|-|-|-|-|-|-|  N   N+1 (1), (2), (3), (7)
    //  f   t   f    |-|-|-|-|-|-|-|--  N   N   (1)
    //  t   t   f    --|-|-|-|-|-|-|-|  N   N   (1)
    //  f   f   t    -|-|-|-|-|-|-|-|-  N   N+1 (4), (6)
    //  t   f   t    -|-|-|-|-|-|-|-|-  N   N+1 (5)
    //  f   t   t    -|-|-|-|-|-|-|-|-  N   N  
    //  t   t   t    -|-|-|-|-|-|-|-|-  N   N  
    //
    // (1) missing end terms presumed zero
    // (2) included last term is usually zero, if non zero, gives secular term
    // (3) zeroth coeff used for secular term
    // (4) zeroth coeff gives constant.
    // (5) secular term should have been removed from samples
    // (6) last coeff is zero (but not for centerp)
    // (7) last coeff is zero (but not for !centerp)
    // Function is represented by (y = pi/h * x)
    // sym = false, sample in f_i = f(h * i/N)
    // odd = true (N samples, N+1 coeffs)
    //    f_0 = 0, need f_i for i in (0, N], f_N defines linear contrib
    // f(x) = c[0] * y + sum(c[k] * sin(k * pi/h * x), k, 1, N)
    // odd = false (N+1 samples, N+1 coeffs)
    //   need f_i for i in [0, N]
    // f(x) = c[0] + sum(c[k] * cos(k * pi/h * x), k, 1, N)
    // sym = true, sample in f_i = f(q * i/N)
    // odd = true (N samples, N coeffs)
    //   f_0 = 0, need f_i for i in (0, N] (N samples)
    // f(x) = sum(c[k] * sin((k+1/2) * pi/q * x), k, 0, N - 1)
    // odd = false (N samples, N coeffs)
    //   f_N = 0, need f_i for i in [0, N) (N samples)
    // f(x) = sum(c[k] * cos((k+1/2) * pi/q * x), k, 0, N - 1)
    
  private:
    Trigfun(int N, bool odd, bool sym, const std::vector<Math::real>& C, real h)
      : _nN(N)
      , _odd(odd)
      , _sym(sym)
      , _C(C)
      , _h(h) {
      _q = _h/2; _p = 2*_h;
    }

  public:
    /**
     * Constructor specifying the number of points to use.
     *
     * @param[in] N the number of points to use.
     **********************************************************************/
    Trigfun(int N = 0);
    static Trigfun initbycoeffs(std::vector<real> C, bool odd, bool sym,
                                real halfp);
    static Trigfun initbysamples(std::vector<real> F, bool odd, bool sym,
                                 bool centerp,
                                 real halfp);
    real check(const std::vector<real>& F, bool centerp) const;
    real eval(real x) const;
    void refine(const Trigfun& tb, const Trigfun& tref);
#if 0
    /**
     * Evaluate the Fourier sum given the sine and cosine of the angle
     *
     * @param[in] sinx sin&sigma;.
     * @param[in] cosx cos&sigma;.
     * @param[in] F the array of Fourier coefficients.
     * @param[in] N the number of Fourier coefficients.
     * @return the value of the Fourier sum.
     **********************************************************************/
    real eval(int i);
    std::vector<real> Coeffs();
    std::vector<real> Samples();
#endif
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIGFUN_HPP

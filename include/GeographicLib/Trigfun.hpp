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

#if defined(_MSC_VER)
// Squelch warnings about dll vs vector
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {
  namespace Triaxial {
  class GeodesicLine3;
  class Conformal3;
  }
  /**
   * \brief Representing a function by a Fourier series
   *
   * This class mimic the functionality of Chebfun's 'trig' representation of
   * periodic functions.  Key differences are:
   *
   * - Only odd or even functions are allowed (i.e., only sine of only cosine
   *   terms in the Fourier series).
   * - Can specify that the function has symmetry about the quarter period
   *   point so that the Fourier series only includes odd harmonics.
   * - The integral of a Trigfun is counted as a Trigfun even if it includes a
   *   secular term.
   * - The inverse function of the integral is also a Trigfun (only makes sense
   *   if the original function is either strictly positive or strictly
   *   negative).
   *
   * Assuming the half period = \e h, the function f(x) is represented as
   * follows, here \f$ = (\pi/h) x\f$:
   * - \e sym = false,
   *   - \e odd = true,
   *     \f$ f(x) = C_0 y + \sum_{m = 1}^{M-1} C_m \sin my \f$;
   *   - \e odd = false,
   *     \f$ f(x) = \sum_{m = 0}^{M-1} C_m \cos my \f$;
   * - \e sym = true,
   *   - \e odd = true,
   *     \f$ f(x) = \sum_{m = 0}^{M-1} C_m \sin 2(m+\frac12) y \f$;
   *   - \e odd = false,
   *     \f$ f(x) = \sum_{m = 0}^{M-1} C_m \cos 2(m+\frac12) y \f$.
   * .
   *
   * Here we compute FFTs using the kissfft package
   * https://github.com/mborgerding/kissfft by Mark Borgerding.
   *
   * Example of use:
   * \include example-Trigfun.cpp
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT Trigfun {
  private:
    /// \cond SKIP
    friend class TrigfunExt;              // For access to root sig 2
    friend class Triaxial::GeodesicLine3; // For access to root sig 2
    friend class Triaxial::Conformal3;    // For access to root sig 2
    /// \endcond
    using real = Math::real;
    static constexpr bool debug_ = false;
    static constexpr bool throw_ = true; // exception on convergence failure
    static constexpr int maxit_ = 300;
    int _m,                     // Number of coefficients in series
      _n;                       // Number of samples in half/quarter period
    bool _odd, _sym;
    std::vector<real> _coeff;
    real _h, _q;                // half, quarter, whole period
    mutable real _max;
    static int chop(const std::vector<real>& c, real tol, real scale = -1);
    static real tolerance(real tol) {
      static const real eps = std::numeric_limits<real>::epsilon();
      return tol <= 0 ? eps : tol;
    }
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

    Trigfun(const std::vector<real>& c, bool odd, bool sym, real h)
      : _m(int(c.size()))
      , _n(sym ? _m : _m - 1)
      , _odd(odd)
      , _sym(sym)
      , _coeff(c)
      , _h(h)
      , _q(_h/2)
      , _max(-1)
    {}
    void refine(const Trigfun& tb);
    real check(const std::vector<real>& F, bool centerp, real tol) const;
    // Given z, return dx = finv(z) - nslope * z
    // dx0 is an estimate of dx (NaN means no information)
    // the "p" in the function name mean periodic (vs secular)
    real inversep(real z, const std::function<real(real)>& fp,
                  real dx0, int* countn, int* countb, real tol) const;
    static Trigfun initbysamples(const std::vector<real>& F,
                                 bool odd, bool sym, real halfp, bool centerp);
    /**
     * Tags to indicate which routine is invoking root().  Because the root
     * functions may be called recursively, each invocation is tagged by an
     * indicator value \e ind.  This is merely an aid to debugging.
     **********************************************************************/
    enum ind {
      NONE = 0,
      INV1,
      INV2,
      ARCPOS0,
      FFUNROOT,
      GFUNROOT,
      INVERSEP,
      PIINV,
      FINV,
      KINV,
      OTHER,
    };

    /**
     * Given \e z, find \e x, such that \e z = \e f(\e x).
     *
     * @param[in] indicator a numeric indicator to track this call (can be
     *   safely set to Trigfun::OTHER).
     * @param[in] z the value of \e f(\e x).
     * @param[in] fp the derivative of \e f(\e x).
     * @param[in] x0 an estimate of the solution, i.e., \e z &asymp; \e f(\e
     *   x0).  Use Math::NaN() to indicate that no estimate is known.
     * @param[in] countn if not nullptr, a pointer to an integer that gets
     *   incremented by the number of iterations.
     * @param[in] countb if not nullptr, a pointer to an integer that gets
     *   incremented by the number of bisection steps (which indicates how well
     *   Newton's method is working).
     * @param[in] tol the tolerance using in terminating the root finding.  \e
     *   tol = 0 (the default) mean to use the machine epsilon.
     * @return the root \e x = \e f <sup>&minus;1</sup>(\e z).
     *
     * Newton's method is used to find the root.  At each step the bounds are
     * adjusted.  If any Newton step gives a result which lies outside the
     * bounds, a bisection step is taken instead.
     *
     * \warning The routine assumes that there's a unique root.  This, in turn,
     *   requires that \e f include a secular term.
     **********************************************************************/
    // root sig 2
    real root(ind indicator, real z, const std::function<real(real)>& fp,
              real x0, int* countn, int* countb, real tol) const;
    /**
     * A general purpose Newton solver for \e z = \e f(\e x).
     *
     * @param[in] indicator a numeric indicator to track this call (can be
     *   safely set to Trigfun::OTHER).
     * @param[in] ffp a function returning \e f(\e x) and \e f'(\e x) as a
     *   pair.
     * @param[in] z the value of \e f(\e x).
     * @param[in] x0 an estimate of the solution, i.e., \e z &cong; \e f(\e
     *   x0).
     * @param[in] xa a lower estimate of the solution.
     * @param[in] xb an upper estimate of the solution.
     * @param[in] xscale a representative scale for \e x.
     * @param[in] zscale a representative scale for \e z.
     * @param[in] s &plusmn;1 depending on whether \e f is an increasing or
     *   decreasing function.
     * @param[in] countn if not nullptr, a pointer to an integer that gets
     *   incremented by the number of iterations.
     * @param[in] countb if not nullptr, a pointer to an integer that gets
     *   incremented by the number of bisection steps (which indicates how well
     *   Newton's method is working).
     * @param[in] tol the tolerance using in terminating the root finding.  \e
     *   tol = 0 (the default) mean to use the machine epsilon.
     * @return the root \e x = \e f <sup>&minus;1</sup>(\e z).
     *
     * This is a static function, so \e f(\e x) need not be a Trigfun.  \e ffp
     * provides both the function an its derivative in one function call to
     * accommodate the (common) situation where the two values can be
     * efficiently computed together.
     *
     * Newton's method is used to find the inverse function.  At each step the
     * bounds are adjusted.  If any Newton step gives a result which lies
     * outside the bounds, a bisection step is taken instead.
     *
     * \warning The routine assumes that there's a unique root lying in the
     *   interval [\e xa, \e xb] and that \e x0 lies in the same interval.
     **********************************************************************/
    // root sig 4
    static real root(ind indicator,
                     const std::function<std::pair<real, real>(real)>& ffp,
                     real z,
                     real x0, real xa, real xb,
                     real xscale = 1, real zscale = 1, int s = 1,
                     int* countn = nullptr, int* countb = nullptr,
                     real tol = 0);
    /**
     * Produce a Trigfun for the inverse of \e f.
     *
     * @param[in] fp the derivative of \e f(\e x).
     * @param[in] countn if not nullptr, a pointer to an integer that gets
     *   incremented by the number of iterations need to create the inverse.
     * @param[in] countb if not nullptr, a pointer to an integer that gets
     *   incremented by the number of bisection steps (which indicates how well
     *   Newton's method is working) needed to create the inverse.
     * @param[in] nmax the maximum number of points in a quarter period
     *   (default 2^16 = 65536).
     * @param[in] tol the tolerance using in terminating the root finding.  \e
     *   tol = 0 (the default) mean to use the machine epsilon.
     * @param[in] scale; if \e scale is negative (the default), \e tol sets the
     *   error relative to the largest Fourier coefficient.  Otherwise, the error
     *   is relative to the maximum of the largest Fourier coefficient and \e
     *   scale.
     * @return the Trigfun representation of \e f <sup>&minus;1</sup>(\e z).
     *
     * As with the normal constructor this routine successively doubles the
     * number of sample points, which are computed using Newton's method.  A
     * good starting guess for Newton's method is provided by the previous
     * Fourier approximation.  As a result the average number of Newton
     * iterations per sample point is about 1 or 2.
     *
     * \note Computing the inverse is only possible with a Trigfun with a
     *   secular term.
     **********************************************************************/
    Trigfun inverse(const std::function<real(real)>& fp,
                    int* countn, int* countb,
                    int nmax, real tol, real scale) const;
  public:
    /**
     * Default constructor specifying with the function \e f(\e x) = 0.
     **********************************************************************/
    Trigfun()
      : _m(1)
      , _n(0)
      , _odd(false)
      , _sym(false)
      , _coeff(1, 0)
      , _h(Math::pi())
      , _q(_h/2)
      , _max(-1)
    {}
    /**
     * Construct a Trigfun with a given number of samples a function.
     *
     * @param[in] n the number of samples.
     * @param[in] f the function.
     * @param[in] odd is the function odd?  If it's not odd, then it is even.
     * @param[in] sym is the function symmetric about the quarter period
     *   point (so it contains only odd Fourier harmonics)?
     * @param[in] halfp the half period.
     * @param[in] centerp whether to sample on a centered grid (default false).
     *
     * For \e sym = false, \e n is the number of samples in a half period, and
     * spacing between the samples is \e halfp/\e n.  (If \e centerp = false
     * and \e oddp = false, the function is, in fact sampled \e n + 1 times.)
     * The number of points given to the FFT routine is 2\e n.
     *
     * For \e sym = true, \e n is the number of samples in a quarter period, and
     * spacing between the samples is \e halfp/(2\e n).  The number of points
     * given to the FFT routine is 4\e n.
     *
     * \note In order for the FFT method to operate efficiently, \e n should be
     * the product of a small factors (typically a power of 2).
     *
     * \warning \e f must be a periodic function and it must be either even or
     *   odd.  With \e odd = true and \e sym = false, the secular term can be
     *   set with setsecular().
     **********************************************************************/
    Trigfun(int n, const std::function<real(real)>& f,
            bool odd, bool sym, real halfp, bool centerp = false);
    /**
     * Construct a Trigfun from a function of one argument.
     *
     * @param[in] f the function.
     * @param[in] odd is the function odd?  If it's not odd, then it is even.
     * @param[in] sym is the function symmetric about the quarter period
     *   point (so it contains only odd Fourier harmonics)?
     * @param[in] halfp the half period.
     * @param[in] nmax the maximum number of points in a half/quarter period
     *   (default 2^16 = 65536).
     * @param[in] tol the tolerance, the default value 0 means use the machine
     *   epsilon.
     * @param[in] scale; if \e scale is negative (the default), \e tol sets the
     *   error relative to the largest Fourier coefficient.  Otherwise, the error
     *   is relative to the maximum of the largest Fourier coefficient and \e
     *   scale.
     *
     * The constructor successively doubles the number of sample points and
     * updating the Fourier coefficients accordingly until the high order
     * coefficients become sufficiently small.  At that point Fourier series is
     * truncated discarding some of the trailing coefficients.  This mimics the
     * method used by Chebfun.  In particular, the method used to truncate the
     * series is taken from Aurentz and L. N. Trefethen,
     * <a href="https://doi.org/10.1145/2998442"> Chopping a Chebyshev
     * series</a> (2017); <a href="https://arxiv.org/abs/1512.01803">
     * preprint</a>.
     *
     * \warning \e f must be a periodic function and it must be either even or
     *   odd.  With \e odd = true and \e sym = false, the secular term can be
     *   set with setsecular().
     **********************************************************************/
    Trigfun(const std::function<real(real)>& f, bool odd, bool sym,
            real halfp, int nmax = 1 << 16,
            real tol = 0,
            real scale = -1);
    /**
     * Construct a Trigfun from a function of two arguments.
     *
     * @param[in] f the function.
     * @param[in] odd is the function odd?  If it's not odd, then it is even.
     * @param[in] sym is the function symmetric about the quarter period
     *   point (so it contains only odd Fourier harmonics)?
     * @param[in] halfp the half period.
     * @param[in] nmax the maximum number of points in a half/quarter period
     *   (default 2^16 = 65536).
     * @param[in] tol the tolerance, the default value 0 means use the machine
     *   epsilon.
     * @param[in] scale; if \e scale is negative (the default), \e tol sets the
     *   error relative to the largest Fourier coefficient.  Otherwise, the error
     *   is relative to the maximum of the largest Fourier coefficient and \e
     *   scale.
     *
     * This accommodates the situation where the inverse of a Trigfun \e g is
     * being computed using inverse().  In this case \e f(\e x, \e y0) returns
     * the value \e y such that \e g(\e y) = \e x.  This is typically found
     * using Newton's method which requires a starting guess \e y0.  In the
     * implementation of inverse(), the Fourier representation is successively
     * refined by doubling the number samples.  At each stage, a good estimate
     * of the function values at the new points is found by using the current
     * Fourier representation.
     *
     * \warning \e f must be a periodic function and it must be either even or
     *   odd.  With \e odd = true and \e sym = false, the secular term can be
     *   set with setsecular().
     **********************************************************************/
    Trigfun(const std::function<real(real, real)>& f, bool odd, bool sym,
            real halfp, int nmax = 1 << 16,
            real tol = 0,
            real scale = -1);
    /**
     * Set the coefficient of the secular term
     *
     * @param[in] f0 the value of \e f(\e halfp).
     *
     * \warning This throws an error unless \e odd = true and \e sym = false.
     **********************************************************************/
    void setsecular(real f0);
    /**
     * Evaluate the Trigfun.
     *
     * @param[in] x the function argument.
     * @return the function value \e f(\e x).
     **********************************************************************/
    real operator()(real x) const;
    // For support of Angle
    // real eval(Angle phi) const;

    /**
     * The integral of a Trigfun.
     *
     * @return the integral.
     *
     * \warning The secular term (only present with \e odd = true and \e sym =
     *   false) is ignored when taking the integral.
     **********************************************************************/
    Trigfun integral() const;
    /**
     * Produce a Trigfun for the inverse of \e f.
     *
     * @param[in] fp the derivative of \e f(\e x).
     * @param[in] nmax the maximum number of points in a quarter period
     *   (default 2^16 = 65536).
     * @param[in] tol the tolerance using in terminating the root finding.  \e
     *   tol = 0 (the default) mean to use the machine epsilon.
     * @param[in] scale; if \e scale is negative (the default), \e tol sets the
     *   error relative to the largest Fourier coefficient.  Otherwise, the error
     *   is relative to the maximum of the largest Fourier coefficient and \e
     *   scale.
     * @return the Trigfun representation of \e f <sup>&minus;1</sup>(\e z).
     *
     * As with the normal constructor this routine successively doubles the
     * number of sample points, which are computed using Newton's method.  A
     * good starting guess for Newton's method is provided by the previous
     * Fourier approximation.  As a result the average number of Newton
     * iterations per sample point is about 1 or 2.
     *
     * \e scale is used when \e f(\e x) is a correction term added to a larger
     * contribution; and it would then be the magnitude of the larger
     * contribution.

     * \note Computing the inverse is only possible with a Trigfun with a
     *   secular term.
     **********************************************************************/
    Trigfun inverse(const std::function<real(real)>& fp,
                    int nmax = 1 << 16, real tol = 0, real scale = -1) const {
      return inverse(fp, nullptr, nullptr, nmax, tol, scale);
    }

    /**
     * @return whether the function is odd or not.  If it's not odd, then it is
     *   even.
     **********************************************************************/
    bool Odd() const { return _odd; }
    /**
     * @return whether the function is symmetric about the quarter period
     *   point.  If it is it, then the Fourier series has only odd terms.
     **********************************************************************/
    bool Symmetric() const { return _sym; }
    /**
     * @return the half period of the function.
     **********************************************************************/
    real HalfPeriod() const { return _h; }
    /**
     * @return the number of terms in the  Fourier series.
     **********************************************************************/
    int NCoeffs() const { return _m; }
    /**
     * @return an estimate of the amplitude of the oscillating component of \e
     *   f.
     *
     * \note This estimate excludes any constant or secular terms in the
     *   series.  The estimate is found by summing the absolute values of the
     *   remaining coefficients (and is thus an overestimate).
     **********************************************************************/
    real Max() const;
    /**
     * @return the (approximate) half range of the function.
     *
     * For a Trigfun containing a secular contribution this is the value of the
     * function at the half perioid.  Otherwise Max() is returned.
     **********************************************************************/
    real HalfRange() const {
      return _odd && !_sym ? _coeff[0] * Math::pi() : Max();
    }
    /**
     * @return the average slope of the function.
     *
     * For a Trigfun containing a secular contribution this is the slope of the
     * secular component.  Otherwise 0 is returned..
     **********************************************************************/
    real Slope() const {
      return _odd && !_sym ? HalfRange() / HalfPeriod() : 0;
    }
  };

  /**
   * \brief A function defined by its derivative and its inverse.
   *
   * This is an extension of Trigfun which allows a function to be defined by
   * its derivative.  In this case the derivative must be even, so that its
   * integral is odd (and taken to be zero at the origin).
   *
   * In addition, this class offers a flexible interface to computing the
   * inverse of the function.  If the inverse is only required at a few points
   *
   * Example of use:
   * \include example-TrigfunExt.cpp
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT TrigfunExt {
  private:
    /// \cond SKIP
    friend class Triaxial::GeodesicLine3; // For access internal inv, inv1
    /// \endcond
    using real = Math::real;
    std::function<real(real)> _fp;
    bool _sym;
    Trigfun _f, _finv;
    real _tol;
    int _nmax;
    bool _invp;
    int _countn, _countb;

    // Approximate inverse using _finv
    real inv0(real z) const {
      if (!_invp) return Math::NaN();
      return _sym ? Math::NaN() : _finv(z);
    }
    // Accurate inverse by direct Newton (not using _finv)
    real inv1(real z, int* countn, int* countb) const {
      return _sym ? Math::NaN() : _f.root(Trigfun::INV1, z, _fp, Math::NaN(),
                                          countn, countb, 0);
    }
    // Accurate inverse correcting result from _finv
    real inv2(real z, int* countn, int* countb) const {
      if (!_invp) return Math::NaN();
      return _sym ? Math::NaN() :
        _f.root(Trigfun::INV2, z, _fp, _finv(z), countn, countb, 0);
    }
    real inv(real z, int* countn, int* countb) const {
      return _invp ? inv2(z, countn, countb) : inv1(z, countn, countb);
    }
  public:
    TrigfunExt() {}
    /**
     * Constructor specifying the derivative, an even periodic function
     *
     * @param[in] fp the derivative function, \e fp(\e x) = \e f'(\e x).
     * @param[in] halfp the half period.
     * @param[in] sym is \e fp symmetric about the quarter period point (so it
     *   contains only odd Fourier harmonics)?
     * @param[in] scale; this is passed to the Trigfun constructor when finding
     *   the Fourier series for \e fp.
     *
     * \warning \e fp must be an even periodic function.  In addition \e fp
     *   must be nonnegative for the inverse of \e f to be computed (in this
     *   case, \e f is a monotonically increasing function).  The inverse is
     *   undefined for \e sym = true.
     **********************************************************************/
    TrigfunExt(const std::function<real(real)>& fp, real halfp,
               bool sym = false, real scale = -1);
    /**
     * Evaluate the TrigfunExt.
     *
     * @param[in] x the function argument.
     * @return the function value \e f(\e x).
     **********************************************************************/
    real operator()(real x) const { return _f(x); }
    /**
     * Evaluate the derivative for TrigfunExt.
     *
     * @param[in] x the function argument.
     * @return the value of the derivative \e fp(\e x).  This uses the function
     *   object passed to the constructor.
     **********************************************************************/
    real deriv(real x) const { return _fp(x); }
    /**
     * Evaluate the inverse of \e f
     *
     * @param[in] z the value of \e f(\e x)
     * @return the value of \e x = \e f <sup>&minus;1</sup>(\e z).
     *
     * This compute the inverse using Newton's method with the derivative
     * function \e fp supplied on construction.  Initially, the starting guess
     * is based on just the secular component of \e f(\e x).  However, if
     * ComputeInverse() is called, a rough Trigfun approximation to the inverse
     * is found and this is used as the starting point for Newton's method.
     **********************************************************************/
    real inv(real z) const {
      return _invp ? inv2(z, nullptr, nullptr) : inv1(z, nullptr, nullptr);
    }
    /**
     * Compute a coarse Fourier series approximation of the inverse.
     *
     * This is used to provide a better starting guess for Newton's method in
     * inv().  Because ComputeInverse() is fairly expensive, this only makes
     * sense if inv() will be called many times.  In order to limit the expense
     * in computing this approximate inverse, the number of Fourier components
     * in the Trigfun for the inverse is limited to 3/2 of the number of
     * components for \e f and the tolerance is set to the square root of the
     * machine epsilon.
     **********************************************************************/
    void ComputeInverse() {
      if (!_invp && !_sym) {
        _countn = _countb = 0;
        _finv = _f.inverse(_fp, &_countn, &_countb, _nmax, _tol, -1);
        _invp = true;
      }
    }
    /**
     * @return whether the function is symmetric about the quarter period
     *   point.  If it is it, then the Fourier series has only odd harmonics.
     **********************************************************************/
    bool Symmetric() const { return _sym; }
    /**
     * @return the number of terms in the  Fourier series for \e f.
     **********************************************************************/
    int NCoeffs() const { return _f.NCoeffs(); }
    /**
     * @return the number of terms in the Fourier series for
     *   \e f<sup>&minus;1</sup>.
     **********************************************************************/
    int NCoeffsInv() const {
      if (!_invp) return -1;
      return _finv.NCoeffs(); }
    /**
     * @return Max() for the underlying Trigfun.
     **********************************************************************/
    real Max() const { return _f.Max(); }
    /**
     * @return HalfPeriod() for the underlying Trigfun.
     **********************************************************************/
    real HalfPeriod() const { return _f.HalfPeriod(); }
    /**
     * @return HalfRange() for the underlying Trigfun.
     **********************************************************************/
    real HalfRange() const { return _f.HalfRange(); }
    /**
     * @return Slope() for the underlying Trigfun.
     **********************************************************************/
    real Slope() const { return _f.Slope(); }
  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_TRIGFUN_HPP

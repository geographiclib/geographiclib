/**
 * \file Constants.hpp
 * \brief Header for GeographicLib::Constants class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010, 2011) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CONSTANTS_HPP)
#define GEOGRAPHICLIB_CONSTANTS_HPP "$Id$"

#include <GeographicLib/Config.h>

/**
 * Are C++0X math functions available?
 **********************************************************************/
#if !defined(GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
#  if defined(__GXX_EXPERIMENTAL_CXX0X__)
#    define GEOGRAPHICLIB_CPLUSPLUS0X_MATH 1
#  else
#    define GEOGRAPHICLIB_CPLUSPLUS0X_MATH 0
#  endif
#endif

/**
 * A compile-time assert.  Use C++0X static_assert, if available.
 **********************************************************************/
#if !defined(STATIC_ASSERT)
#  if defined(__GXX_EXPERIMENTAL_CXX0X__)
#    define STATIC_ASSERT static_assert
#  elif defined(_MSC_VER) && _MSC_VER >= 1600
#    define STATIC_ASSERT static_assert
#  else
#    define STATIC_ASSERT(cond,reason) \
            { enum{ STATIC_ASSERT_ENUM = 1/int(cond) }; }
#  endif
#endif

#if defined(__GNUC__)
// Suppress "defined but not used" warnings
#  define RCSID_DECL(x) namespace \
{ char VAR_ ## x [] __attribute__((used)) = x; }
#else
/**
 * Insertion of RCS Id strings into the object file.
 **********************************************************************/
#  define RCSID_DECL(x) namespace { char VAR_ ## x [] = x; }
#endif

#if defined(_WIN32) && defined(GEOGRAPHIC_SHARED_LIB)
#  if defined(Geographic_EXPORTS)
#    define GEOGRAPHIC_EXPORT __declspec(dllexport)
#  else
#    define GEOGRAPHIC_EXPORT __declspec(dllimport)
#  endif
#else
#  define GEOGRAPHIC_EXPORT
#endif

RCSID_DECL(GEOGRAPHICLIB_CONSTANTS_HPP)

#if !defined(GEOGRAPHICLIB_PREC)
/**
 * The precision of floating point numbers used in %GeographicLib.  0 means
 * float; 1 (default) means double; 2 means long double.  Nearly all the
 * testing has been carried out with doubles and that's the recommended
 * configuration.  Note that with Microsoft Visual Studio, long double is the
 * same as double.
 **********************************************************************/
#define GEOGRAPHICLIB_PREC 1
#endif

#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>

/**
 * \brief Namespace for %GeographicLib
 *
 * All of %GeographicLib is defined within the GeographicLib namespace.  In
 * addtion all the header files are included via %GeographicLib/filename.  This
 * minimizes the likelihood of conflicts with other packages.
 **********************************************************************/
namespace GeographicLib {

  /**
   * \brief Mathematical functions needed by %GeographicLib
   *
   * Define mathematical functions in order to localize system dependencies and
   * to provide generic versions of the functions.  In addition define a real
   * type to be used by %GeographicLib.
   **********************************************************************/
  class GEOGRAPHIC_EXPORT Math {
  private:
    void dummy() {
      STATIC_ASSERT((GEOGRAPHICLIB_PREC) >= 0 && (GEOGRAPHICLIB_PREC) <= 2,
                    "Bad value of precision");
    }
    Math();                     // Disable constructor
  public:

#if defined(HAVE_LONG_DOUBLE)
    /**
     * The extended precision type for real numbers, used for some testing.
     * This is long double on computers with this type; otherwise it is double.
     **********************************************************************/
    typedef long double extended;
#else
    typedef double extended;
#endif

#if GEOGRAPHICLIB_PREC == 1
    /**
     * The real type for %GeographicLib. Nearly all the testing has been done
     * with \e real = double.  However, the algorithms should also work with
     * float and long double (where available).
     **********************************************************************/
    typedef double real;
#elif GEOGRAPHICLIB_PREC == 0
    typedef float real;
#elif GEOGRAPHICLIB_PREC == 2
    typedef extended real;
#else
    typedef double real;
#endif

    /**
     * @return \e pi
     **********************************************************************/
    template<typename T>
    static inline T pi() throw() { return std::atan2(T(0), -T(1)); }
    /**
     * A synonym for pi<real>().
     **********************************************************************/
    static inline real pi() throw() { return pi<real>(); }
    /**
     * <b>DEPRECATED</b> A synonym for pi<extened>().
     **********************************************************************/
    static inline extended epi() throw() { return pi<extended>(); }

    /**
     * @return the number of radians in a degree.
     **********************************************************************/
    template<typename T>
    static inline T degree() throw() { return pi<T>() / T(180); }
    /**
     * A synonym for degree<real>().
     **********************************************************************/
    static inline real degree() throw() { return degree<real>(); }
    /**
     * <b>DEPRECATED</b> A synonym for degree<extened>().
     **********************************************************************/
    static inline extended edegree() throw() { return degree<extended>(); }

    /**
     * Square a number.
     * @param[in] x
     * @return \e x<sup>2</sup>.
     **********************************************************************/
    template<typename T>
    static inline T sq(T x) throw() { return x * x; }

#if defined(DOXYGEN)
    /**
     * The hypotenuse function avoiding underflow and overflow.
     *
     * @param[in] x
     * @param[in] y
     * @return sqrt(\e x<sup>2</sup> + \e y<sup>2</sup>).
     **********************************************************************/
    template<typename T>
    static inline T hypot(T x, T y) throw() {
      x = std::abs(x);
      y = std::abs(y);
      T a = (std::max)(x, y),
        b = (std::min)(x, y) / (a ? a : 1);
      return a * std::sqrt(1 + b * b);
    }
#elif GEOGRAPHICLIB_CPLUSPLUS0X_MATH
    template<typename T>
    static inline T hypot(T x, T y) throw() { return std::hypot(x, y); }
#elif defined(_MSC_VER)
    static inline double hypot(double x, double y) throw()
    { return _hypot(x, y); }
    static inline float hypot(float x, float y) throw()
    { return _hypotf(x, y); }
#if defined(HAVE_LONG_DOUBLE)
    static inline long double hypot(long double x, long double y) throw()
    { return _hypot(x, y); }
#endif
#else
    // Use overloading to define generic versions
    static inline double hypot(double x, double y) throw()
    { return ::hypot(x, y); }
    static inline float hypot(float x, float y) throw()
    { return ::hypotf(x, y); }
#if defined(HAVE_LONG_DOUBLE)
    static inline long double hypot(long double x, long double y) throw()
    { return ::hypotl(x, y); }
#endif
#endif

#if defined(DOXYGEN) || (defined(_MSC_VER) && !GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
    /**
     * exp(\e x) - 1 accurate near \e x = 0.  This is taken from
     * N. J. Higham, Accuracy and Stability of Numerical Algorithms, 2nd
     * Edition (SIAM, 2002), Sec 1.14.1, p 19.
     *
     * @param[in] x
     * @return exp(\e x) - 1.
     **********************************************************************/
    template<typename T>
    static inline T expm1(T x) throw() {
      volatile T
        y = std::exp(x),
        z = y - 1;
      // The reasoning here is similar to that for log1p.  The expression
      // mathematically reduces to exp(x) - 1, and the factor z/log(y) = (y -
      // 1)/log(y) is a slowly varying quantity near y = 1 and is accurately
      // computed.
      return std::abs(x) > 1 ? z : z == 0 ?  x : x * z / std::log(y);
    }
#elif GEOGRAPHICLIB_CPLUSPLUS0X_MATH
    template<typename T>
    static inline T expm1(T x) throw() { return std::expm1(x); }
#else
    static inline double expm1(double x) throw() { return ::expm1(x); }
    static inline float expm1(float x) throw() { return ::expm1f(x); }
#if defined(HAVE_LONG_DOUBLE)
    static inline long double expm1(long double x) throw()
    { return ::expm1l(x); }
#endif
#endif

#if defined(DOXYGEN) || (defined(_MSC_VER) && !GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
    /**
     * log(\e x + 1) accurate near \e x = 0.
     *
     * This is taken from D. Goldberg,
     * <a href="http://dx.doi.org/10.1145/103162.103163">What every computer
     * scientist should know about floating-point arithmetic</a> (1991),
     * Theorem 4.  See also, Higham (op. cit.), Answer to Problem 1.5, p 528.
     *
     * @param[in] x
     * @return log(\e x + 1).
     **********************************************************************/
    template<typename T>
    static inline T log1p(T x) throw() {
      volatile T
        y = 1 + x,
        z = y - 1;
      // Here's the explanation for this magic: y = 1 + z, exactly, and z
      // approx x, thus log(y)/z (which is nearly constant near z = 0) returns
      // a good approximation to the true log(1 + x)/x.  The multiplication x *
      // (log(y)/z) introduces little additional error.
      return z == 0 ? x : x * std::log(y) / z;
    }
#elif GEOGRAPHICLIB_CPLUSPLUS0X_MATH
    template<typename T>
    static inline T log1p(T x) throw() { return std::log1p(x); }
#else
    static inline double log1p(double x) throw() { return ::log1p(x); }
    static inline float log1p(float x) throw() { return ::log1pf(x); }
#if defined(HAVE_LONG_DOUBLE)
    static inline long double log1p(long double x) throw()
    { return ::log1pl(x); }
#endif
#endif

#if defined(DOXYGEN) || (defined(_MSC_VER) && !GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
    /**
     * The inverse hyperbolic sine function.  This is defined in terms of
     * Math::log1p(\e x) in order to maintain accuracy near \e x = 0.  In
     * addition, the odd parity of the function is enforced.
     *
     * @param[in] x
     * @return asinh(\e x).
     **********************************************************************/
    template<typename T>
    static inline T asinh(T x) throw() {
      T y = std::abs(x);     // Enforce odd parity
      y = log1p(y * (1 + y/(hypot(T(1), y) + 1)));
      return x < 0 ? -y : y;
    }
#elif GEOGRAPHICLIB_CPLUSPLUS0X_MATH
    template<typename T>
    static inline T asinh(T x) throw() { return std::asinh(x); }
#else
    static inline double asinh(double x) throw() { return ::asinh(x); }
    static inline float asinh(float x) throw() { return ::asinhf(x); }
#if defined(HAVE_LONG_DOUBLE)
    static inline long double asinh(long double x) throw()
    { return ::asinhl(x); }
#endif
#endif

#if defined(DOXYGEN) || (defined(_MSC_VER) && !GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
    /**
     * The inverse hyperbolic tangent function.  This is defined in terms of
     * Math::log1p(\e x) in order to maintain accuracy near \e x = 0.  In
     * addition, the odd parity of the function is enforced.
     *
     * @param[in] x
     * @return atanh(\e x).
     **********************************************************************/
    template<typename T>
    static inline T atanh(T x) throw() {
      T y = std::abs(x);     // Enforce odd parity
      y = log1p(2 * y/(1 - y))/2;
      return x < 0 ? -y : y;
    }
#elif GEOGRAPHICLIB_CPLUSPLUS0X_MATH
    template<typename T>
    static inline T atanh(T x) throw() { return std::atanh(x); }
#else
    static inline double atanh(double x) throw() { return ::atanh(x); }
    static inline float atanh(float x) throw() { return ::atanhf(x); }
#if defined(HAVE_LONG_DOUBLE)
    static inline long double atanh(long double x) throw()
    { return ::atanhl(x); }
#endif
#endif

#if defined(DOXYGEN) || (defined(_MSC_VER) && !GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
    /**
     * The cube root function.
     *
     * @param[in] x
     * @return the real cube root of \e x.
     **********************************************************************/
    template<typename T>
    static inline T cbrt(T x) throw() {
      T y = std::pow(std::abs(x), 1/T(3)); // Return the real cube root
      return x < 0 ? -y : y;
    }
#elif GEOGRAPHICLIB_CPLUSPLUS0X_MATH
    template<typename T>
    static inline T cbrt(T x) throw() { return std::cbrt(x); }
#else
    static inline double cbrt(double x) throw() { return ::cbrt(x); }
    static inline float cbrt(float x) throw() { return ::cbrtf(x); }
#if defined(HAVE_LONG_DOUBLE)
    static inline long double cbrt(long double x) throw() { return ::cbrtl(x); }
#endif
#endif

    /**
     * Test for finiteness.
     *
     * @param[in] x
     * @return true if number is finite, false if NaN or infinite.
     **********************************************************************/
    template<typename T>
    static inline bool isfinite(T x) throw() {
#if defined(DOXYGEN)
      return std::abs(x) <= (std::numeric_limits<T>::max)();
#elif (defined(_MSC_VER) && !GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
      return _finite(x) != 0;
#else
      return std::isfinite(x);
#endif
    }

    /**
     * The NaN (not a number)
     *
     * @return NaN if available, otherwise return the max real.
     **********************************************************************/
    template<typename T>
    static inline T NaN() throw() {
      return std::numeric_limits<T>::has_quiet_NaN ?
        std::numeric_limits<T>::quiet_NaN() :
        (std::numeric_limits<T>::max)();
    }
    /**
     * A synonym for NaN<real>().
     **********************************************************************/
    static inline real NaN() throw() { return NaN<real>(); }

    /**
     * Test for NaN.
     *
     * @param[in] x
     * @return true if argument is a NaN.
     **********************************************************************/
    template<typename T>
    static inline bool isnan(T x) throw() {
#if defined(DOXYGEN) || (defined(_MSC_VER) && !GEOGRAPHICLIB_CPLUSPLUS0X_MATH)
      return x != x;
#else
      return std::isnan(x);
#endif
    }

    /**
     * Infinity
     *
     * @return infinity if available, otherwise return the max real.
     **********************************************************************/
    template<typename T>
    static inline T infinity() throw() {
      return std::numeric_limits<T>::has_infinity ?
        std::numeric_limits<T>::infinity() :
        (std::numeric_limits<T>::max)();
    }
    /**
     * A synonym for infinity<real>().
     **********************************************************************/
    static inline real infinity() throw() { return infinity<real>(); }
  };

  /**
   * \brief %Constants needed by %GeographicLib
   *
   * Define constants specifying the WGS84 ellipsoid, the UTM and UPS
   * projections, and various unit conversions.
   **********************************************************************/
  class GEOGRAPHIC_EXPORT Constants {
  private:
    typedef Math::real real;
    Constants();                // Disable constructor

  public:
    /**
     * <b>DEPRECATED</b> A synonym for Math::pi<real>().
     **********************************************************************/
    static inline Math::real pi() throw() { return Math::pi<real>(); }
    /**
     * A synonym for Math::degree<real>().
     **********************************************************************/
    static inline Math::real degree() throw() { return Math::degree<real>(); }
    /**
     * @return the number of radians in an arcminute.
     **********************************************************************/
    static inline Math::real arcminute() throw()
    { return Math::degree<real>() / 60; }
    /**
     * @return the number of radians in an arcsecond.
     **********************************************************************/
    static inline Math::real arcsecond() throw()
    { return Math::degree<real>() / 3600; }

    /** \name Ellipsoid parameters
     **********************************************************************/
    ///@{
    /**
     * @return the equatorial radius of WGS84 ellipsoid
     **********************************************************************/
    template<typename T>
    static inline T WGS84_a() throw() { return T(6378137) * meter<T>(); }
    /**
     * A synonym for WGS84_a<real>().
     **********************************************************************/
    static inline Math::real WGS84_a() throw() { return WGS84_a<real>(); }
    /**
     * @return the flattening of WGS84 ellipsoid
     **********************************************************************/
    template<typename T>
    static inline T WGS84_f() throw() {
      // 1/298.257223563
      return T(1) / ( T(298) + T(257223563) / T(1000000000) );
    }
    /**
     * A synonym for WGS84_f<real>().
     **********************************************************************/
    static inline Math::real WGS84_f() throw() { return WGS84_f<real>(); }
    /**
     * @return the reciprocal flattening of WGS84 ellipsoid
     **********************************************************************/
    template<typename T>
    static inline T WGS84_r() throw() { return 1/WGS84_f<T>(); }
    /**
     * A synonym for WGS84_r<real>().
     **********************************************************************/
    static inline Math::real WGS84_r() throw() { return WGS84_r<real>(); }
    /**
     * @return the central scale factor for UTM
     **********************************************************************/
    template<typename T>
    static inline T UTM_k0() throw() {return T(9996) / T(10000); } // 0.9996
    /**
     * A synonym for UTM_k0<real>().
     **********************************************************************/
    static inline Math::real UTM_k0() throw() { return UTM_k0<real>(); }
    /**
     * @return the central scale factor for UPS
     **********************************************************************/
    template<typename T>
    static inline T UPS_k0() throw() { return T(994) / T(1000); } // 0.994
    /**
     * A synonym for UPS_k0<real>().
     **********************************************************************/
    static inline Math::real UPS_k0() throw() { return UPS_k0<real>(); }
    ///@}

    /** \name SI units
     **********************************************************************/
    ///@{
    /**
     * @return the number of meters in a meter.
     *
     * This is unity, but this lets the internal system of units be changed if
     * necessary.
     **********************************************************************/
    template<typename T>
    static inline T meter() throw() { return T(1); }
    /**
     * A synonym for meter<real>().
     **********************************************************************/
    static inline Math::real meter() throw() { return meter<real>(); }
    /**
     * @return the number of meters in a kilometer.
     **********************************************************************/
    static inline Math::real kilometer() throw()
    { return 1000 * meter<real>(); }
    /**
     * @return the number of meters in a nautical mile (approximately 1 arc
     *   minute)
     **********************************************************************/
    static inline Math::real nauticalmile() throw()
    { return 1852 * meter<real>(); }
    ///@}

    /** \name Anachronistic British units
     **********************************************************************/
    ///@{
    /**
     * @return the number of meters in an international foot.
     **********************************************************************/
    static inline Math::real foot() throw()
    { return real(0.0254L) * 12 * meter<real>(); }
    /**
     * @return the number of meters in a yard.
     **********************************************************************/
    static inline Math::real yard() throw() { return 3 * foot(); }
    /**
     * @return the number of meters in a fathom.
     **********************************************************************/
    static inline Math::real fathom() throw() { return 2 * yard(); }
    /**
     * @return the number of meters in a chain.
     **********************************************************************/
    static inline Math::real chain() throw() { return 22 * yard(); }
    /**
     * @return the number of meters in a furlong.
     **********************************************************************/
    static inline Math::real furlong() throw() { return 10 * chain(); }
    /**
     * @return the number of meters in a statute mile.
     **********************************************************************/
    static inline Math::real mile() throw() { return 8 * furlong(); }
    ///@}

    /** \name Anachronistic US units
     **********************************************************************/
    ///@{
    /**
     * @return the number of meters in a US survey foot.
     **********************************************************************/
    static inline Math::real surveyfoot() throw()
    { return real(1200) / real(3937) * meter<real>(); }
    ///@}
  };

  /**
   * \brief An accumulator for sums.
   *
   * This allow many numbers of floating point type \e T to be added together
   * with twice the normal precision.  Thus if \e T is double, the effective
   * precision of the sum is 106 bits or about 32 decimal places.  The core
   * idea is the error free transformation of a sum, D. E. Knuth, TAOCP, Vol 2,
   * 4.2.2, Theorem B.
   *
   * Two implementations are provided.  The "fast" one is an implementation of
   * Algorithm 4.1, of T. Ogita, S. M. Rump, S. Oishi,
   * <a href="http://dx.doi.org/10.1137/030601818"> Accurate sum and dot
   * product</a>, SIAM J. Sci. Comp., 26(6) 1955-1988 (2005).  The accumulator
   * is represented by a two numbers _s + _t where _s is the result of the
   * normal sum and _t accumulates (approximately) the exact roundoff errors in
   * the summation of _s.
   *
   * The slow "non-fast" implementation follows J. R. Shewchuk,
   * <a href="http://dx.doi.org/10.1007/PL00009321"> Adaptive Precision
   * Floating-Point Arithmetic and Fast Robust Geometric Predicates</a>,
   * Discrete & Computational Geometry 18(3) 305-363 (1997).  In this case,
   * with each addition of a number to the accumulator, _s and _t are adjusted
   * so that _s represents the sum to an accuracy of 1 ulp.  This may result in
   * considerably greater accuracy than the fast implementation (depending on
   * the input data).
   *
   * Approximate timings (summing vector<double> per addition)
   * - double:                      2ns = 1
   * - Accumulator<double, true>:   7ns = 3.5
   * - Accumulator<double, true>:  21ns = 10 -- Optimize() on each addition
   * - Accumulator<double, false>: 23ns = 11
   * .
   * Thus the slow method is about 3 times slower than the fast method.
   * However all times are negligible with the typical timings for geodesic
   * calculations which are roughly 2us.  Thus the default mode is taken to be
   * slow, \e fast = false.
   *
   * In the documentation of the member functions, \e sum stands for the value
   * currently held in the accumulator.
   **********************************************************************/
  template<typename T = Math::real, bool fast = false>
  class Accumulator {
  private:
    // _s accumulates for the straight sum
    // _t accumulates the errors.
    T _s, _t;
    // Error free transformation of a sum.  Note that t can be the same as one
    // of the first two arguments.
    static inline T sum(T u, T v, T& t) {
      volatile T s = u + v;
      volatile T up = s - v;
      volatile T vpp = s - up;
      up -= u;
      vpp -= v;
      t = -(up + vpp);
      // u + v =       s      + t
      //       = round(u + v) + t
      return s;
    }
    // Same as sum, but requires abs(u) >= abs(v).  This isn't currently used.
    static inline T fastsum(T u, T v, T& t) {
      volatile T s = u + v;
      volatile T vp = s - u;
      t = v - vp;
      return s;
    }
    void Add(T y) throw() {
      if (fast) {
        _s = sum(_s, y, y);     // Accumulate sum to _s with error in y;
        _t += y;                // accumulate y in _t
      } else {                  // Here's Shewchuk's solution...
        T u;                    // hold exact sum as [s, t, u]
        y  = sum(y, _t,  u);    // Accumulate starting at least significant end
        _s = sum(y, _s, _t);
        // Start is _s, _t decreasing and non-adjacent.  Sum is now (s + t + u)
        // exactly with s, t, u non-adjacent and in decreasing order (except
        // for possible zeros).  The following code tries to normalize the
        // result.  Ideally, we want _s = round(s+t+u) and _u = round(s+t+u -
        // _s).  The follow does an approximate job (and maintains the
        // decreasing non-adjacent property).  Here are two "failures" using
        // 3-bit floats:
        //
        // Case 1: _s is not equal to round(s+t+u) -- off by 1 ulp
        // [12, -1] - 8 -> [4, 0, -1] -> [4, -1] = 3 should be [3, 0] = 3
        //
        // Case 2: _s+_t is not as close to s+t+u as it shold be
        // [64, 5] + 4 -> [64, 8, 1] -> [64,  8] = 72 (off by 1)
        //                    should be [80, -7] = 73 (exact)
        //
        // "Fixing" these problems is probably not worth the expense.  The
        // representation inevitably leads to small errors in the accumulated
        // values.  The additional errors illustrated here amount to 1 ulp of
        // the less significant word during each addition to the Accumulator
        // and an additional possible error of 1 ulp in the reported sum.
        //
        // Incidentally, the "ideal" representation described above is not
        // canonical, because _s = round(_s + _t) may not be true.  For
        // example, with 3-bit floats:
        //
        // [128, 16] + 1 -> [160, -16] -- 160 = round(145).
        // But [160, 0] - 16 -> [128, 16] -- 128 = round(144).
        //
        if (_s == 0)            // This implies t == 0,
          _s = u;               // so result is u
        else
          _t += u;              // otherwise just accumulate u to t.
      }
    }
    // Perhaps these should both be _s + _t?
    T Sum() const throw() { return fast ? _s + _t : _s; }
    T Sum(T y) const throw() {
      Accumulator a(*this);
      a.Add(y);
      return a.Sum();
    }
  public:
    /**
     * Construct from a \e T.  This is not declared explicit, so that you can
     * write <code>Accumulator<double> a = 5;</code>.
     *
     * @param[in] y set \e sum = \e y.
     **********************************************************************/
    Accumulator(T y = T(0)) throw() : _s(y), _t(0) {
      STATIC_ASSERT(!std::numeric_limits<T>::is_integer,
                    "Accumulator type is not floating point");
    };
    /**
     * Set the accumulator to a number.
     *
     * @param[in] y set \e sum = \e y.
     **********************************************************************/
    Accumulator& operator=(T y) throw() { _s = y; _t = 0; return *this; }
    /**
     * Return the value held in the accumulator.
     *
     * @return \e sum.
     **********************************************************************/
    T operator()() const throw() { return Sum(); }
    /**
     * Return the result of adding a number to \e sum (but don't change \e sum).
     *
     * @param[in] y the number to be added to the sum.
     * @return \e sum + \e y.
     **********************************************************************/
    T operator()(T y) const throw() { return Sum(y); }
    /**
     * Add a number to the accumulator.
     *
     * @param[in] y set \e sum += \e y.
     **********************************************************************/
    Accumulator& operator+=(T y) throw() { Add(y); return *this; }
    /**
     * Subtract a number from the accumulator.
     *
     * @param[in] y set \e sum -= \e y.
     **********************************************************************/
    Accumulator& operator-=(T y) throw() { Add(-y); return *this; }
    /**
     * Multiply accumulator by an integer.  To avoid loss of accuracy, use only
     * integers such that \e n * \e T is exactly representable as a \e T (i.e.,
     * +/- powers of two).  Use \e n = -1 to negate \e sum.
     *
     * @param[in] n set \e sum *= \e n.
     **********************************************************************/
    Accumulator& operator*=(int n) throw() { _s *= n; _t *= n; return *this; }
    /**
     * Optimize how \e sum is stored in the accumulator.  This is a no-op for
     * the non-fast implementation.  It is rarely necessary to do this; however
     * the accuracy might be improved if this function is called every time a
     * million numbers (for example) have been added to the accumulator.
     **********************************************************************/
    void Optimize() throw() { if (fast) _s = sum(_s, _t, _t); }
    /**
     * Test equality of an Accumulator with a number.
     **********************************************************************/
    bool operator==(T y) const throw() { return Sum() == y; }
    /**
     * Test inequality of an Accumulator with a number.
     **********************************************************************/
    bool operator!=(T y) const throw() { return Sum() != y; }
    /**
     * Less operator on an Accumulator and a number.
     **********************************************************************/
    bool operator<(T y) const throw() { return Sum() < y; }
    /**
     * Less or equal operator on an Accumulator and a number.
     **********************************************************************/
    bool operator<=(T y) const throw() { return Sum() <= y; }
    /**
     * Greater operator on an Accumulator and a number.
     **********************************************************************/
    bool operator>(T y) const throw() { return Sum() > y; }
    /**
     * Greater or equal operator on an Accumulator and a number.
     **********************************************************************/
    bool operator>=(T y) const throw() { return Sum() >= y; }
  };

  /**
   * \brief Exception handling for %GeographicLib
   *
   * A class to handle exceptions.  It's derived from std::runtime_error so it
   * can be caught by the usual catch clauses.
   **********************************************************************/
  class GeographicErr : public std::runtime_error {
  public:

    /**
     * Constructor
     *
     * @param[in] msg a string message, which is accessible in the catch
     *   clause, via what().
     **********************************************************************/
    GeographicErr(const std::string& msg) : std::runtime_error(msg) {}
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_CONSTANTS_HPP

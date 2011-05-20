/**
 * \file Constants.hpp
 * \brief Header for GeographicLib::Constants class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CONSTANTS_HPP)
#define GEOGRAPHICLIB_CONSTANTS_HPP "$Id: Constants.hpp 6867 2010-09-11 13:04:26Z karney $"

/**
 * A simple compile-time assert.  This is designed to be compatible with the
 * C++0X static_assert.
 **********************************************************************/
#if !defined(STATIC_ASSERT)
#define STATIC_ASSERT(cond,reason) { enum{ STATIC_ASSERT_ENUM=1/int(cond) }; }
#endif

#if defined(__GNUC__)
// Suppress "defined but not used" warnings
#define RCSID_DECL(x) namespace \
{ char VAR_ ## x [] __attribute__((used)) = x; }
#else
/**
 * Insertion of RCS Id strings into the object file.
 **********************************************************************/
#define RCSID_DECL(x) namespace { char VAR_ ## x [] = x; }
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
  class Math {
  private:
    void dummy() {
      STATIC_ASSERT((GEOGRAPHICLIB_PREC) >= 0 && (GEOGRAPHICLIB_PREC) <= 2,
                    "Bad value of precision");
    }
    Math();                     // Disable constructor
  public:

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
    typedef long double real;
#else
    typedef double real;
#endif

#if defined(DOXYGEN)
    /**
     * The hypotenuse function avoiding underflow and overflow.
     *
     * @param[in] x
     * @param[in] y
     * @return sqrt(\e x<sup>2</sup> + \e y<sup>2</sup>).
     **********************************************************************/
    static inline real hypot(real x, real y) throw() {
      x = std::abs(x);
      y = std::abs(y);
      real
        a = std::max(x, y),
        b = std::min(x, y) / a;
      return a * std::sqrt(1 + b * b);
    }
#elif defined(_MSC_VER)
    static inline double hypot(double x, double y) throw()
    { return _hypot(x, y); }
    static inline float hypot(float x, float y) throw()
    { return _hypotf(x, y); }
#if !defined(__NO_LONG_DOUBLE_MATH)
    static inline long double hypot(long double x, long double y) throw()
    { return _hypot(x, y); }
#endif
#else
    // Use overloading to define generic versions
    static inline double hypot(double x, double y) throw()
    { return ::hypot(x, y); }
    static inline float hypot(float x, float y) throw()
    { return ::hypotf(x, y); }
#if !defined(__NO_LONG_DOUBLE_MATH)
    static inline long double hypot(long double x, long double y) throw()
    { return ::hypotl(x, y); }
#endif
#endif

#if defined(DOXYGEN) || defined(_MSC_VER)
    /**
     * exp(\e x) - 1 accurate near \e x = 0.  This is taken from
     * N. J. Higham, Accuracy and Stability of Numerical Algorithms, 2nd
     * Edition (SIAM, 2002), Sec 1.14.1, p 19.
     *
     * @param[in] x
     * @return exp(\e x) - 1.
     **********************************************************************/
    static inline real expm1(real x) throw() {
      volatile real
        y = std::exp(x),
        z = y - 1;
      // The reasoning here is similar to that for log1p.  The expression
      // mathematically reduces to exp(x) - 1, and the factor z/log(y) = (y -
      // 1)/log(y) is a slowly varying quantity near y = 1 and is accurately
      // computed.
      return std::abs(x) > 1 ? z : z == 0 ?  x : x * z / std::log(y);
    }
#else
    static inline double expm1(double x) throw() { return ::expm1(x); }
    static inline float expm1(float x) throw() { return ::expm1f(x); }
#if !defined(__NO_LONG_DOUBLE_MATH)
    static inline long double expm1(long double x) throw()
    { return ::expm1l(x); }
#endif
#endif

#if defined(DOXYGEN) || defined(_MSC_VER)
    /**
     * log(\e x + 1) accurate near \e x = 0.  This is taken See
     * D. Goldberg,
     * <a href="http://docs.sun.com/source/806-3568/ncg_goldberg.html"> What
     * every computer scientist should know about floating-point arithmetic</a>
     * (1991), Theorem 4.  See also, Higham (op. cit.), Answer to Problem 1.5,
     * p 528.
     *
     * @param[in] x
     * @return log(\e x + 1).
     **********************************************************************/
    static inline real log1p(real x) throw() {
      volatile real
        y = 1 + x,
        z = y - 1;
      // Here's the explanation for this magic: y = 1 + z, exactly, and z
      // approx x, thus log(y)/z (which is nearly constant near z = 0) returns
      // a good approximation to the true log(1 + x)/x.  The multiplication x *
      // (log(y)/z) introduces little additional error.
      return z == 0 ? x : x * std::log(y) / z;
    }
#else
    static inline double log1p(double x) throw() { return ::log1p(x); }
    static inline float log1p(float x) throw() { return ::log1pf(x); }
#if !defined(__NO_LONG_DOUBLE_MATH)
    static inline long double log1p(long double x) throw()
    { return ::log1pl(x); }
#endif
#endif

#if defined(DOXYGEN) || defined(_MSC_VER)
    /**
     * The inverse hyperbolic sine function.  This is defined in terms of
     * Math::log1p(\e x) in order to maintain accuracy near \e x = 0.  In
     * addition, the odd parity of the function is enforced.
     *
     * @param[in] x
     * @return asinh(\e x).
     **********************************************************************/
    static inline real asinh(real x) throw() {
      real y = std::abs(x);     // Enforce odd parity
      y = log1p(y * (1 + y/(hypot(real(1), y) + 1)));
      return x < 0 ? -y : y;
    }
#else
    static inline double asinh(double x) throw() { return ::asinh(x); }
    static inline float asinh(float x) throw() { return ::asinhf(x); }
#if !defined(__NO_LONG_DOUBLE_MATH)
    static inline long double asinh(long double x) throw()
    { return ::asinhl(x); }
#endif
#endif

#if defined(DOXYGEN) || defined(_MSC_VER)
    /**
     * The inverse hyperbolic tangent function.  This is defined in terms of
     * Math::log1p(\e x) in order to maintain accuracy near \e x = 0.  In
     * addition, the odd parity of the function is enforced.
     *
     * @param[in] x
     * @return atanh(\e x).
     **********************************************************************/
    static inline real atanh(real x) throw() {
      real y = std::abs(x);     // Enforce odd parity
      y = log1p(2 * y/(1 - y))/2;
      return x < 0 ? -y : y;
    }
#else
    static inline double atanh(double x) throw() { return ::atanh(x); }
    static inline float atanh(float x) throw() { return ::atanhf(x); }
#if !defined(__NO_LONG_DOUBLE_MATH)
    static inline long double atanh(long double x) throw()
    { return ::atanhl(x); }
#endif
#endif

#if defined(DOXYGEN) || defined(_MSC_VER)
    /**
     * The cube root function.
     *
     * @param[in] x
     * @return the real cube root of \e x.
     **********************************************************************/
    static inline real cbrt(real x) throw() {
      real y = std::pow(std::abs(x), 1/real(3)); // Return the real cube root
      return x < 0 ? -y : y;
    }
#else
    static inline double cbrt(double x) throw() { return ::cbrt(x); }
    static inline float cbrt(float x) throw() { return ::cbrtf(x); }
#if !defined(__NO_LONG_DOUBLE_MATH)
    static inline long double cbrt(long double x) throw() { return ::cbrtl(x); }
#endif
#endif

    /**
     * Test for finiteness.
     *
     * @param[in] x
     * @return true if number is finite, false if NaN or infinite.
     **********************************************************************/
    static inline bool isfinite(real x) throw() {
#if defined(DOXYGEN)
      return std::abs(x) <= std::numeric_limits<real>::max();
#elif defined(_MSC_VER)
      return _finite(x) != 0;
#else
      return std::isfinite(x) != 0;
#endif
    }

    /**
     * The NaN (not a number)
     *
     * @return NaN if available, otherwise return the max real.
     **********************************************************************/
    static inline real NaN() throw() {
      return std::numeric_limits<real>::has_quiet_NaN ?
        std::numeric_limits<real>::quiet_NaN() :
        std::numeric_limits<real>::max();
    }
  };

  /**
   * \brief %Constants needed by %GeographicLib
   *
   * Define constants specifying the WGS84 ellipsoid, the UTM and UPS
   * projections, and various unit conversions.
   **********************************************************************/
  class Constants {
  private:
    typedef Math::real real;
    Constants();                // Disable constructor

  public:
    /**
     * @return \e pi
     **********************************************************************/
    static inline Math::real pi() throw()
    // good for about 123-bit accuracy
    { return real(3.141592653589793238462643383279502884L); }
    /**
     * @return the number of radians in a degree.
     **********************************************************************/
    static inline Math::real degree() throw() { return pi() / 180; }
    /**
     * @return the number of radians in an arcminute.
     **********************************************************************/
    static inline Math::real arcminute() throw() { return degree() / 60; }
    /**
     * @return the number of radians in an arcsecond.
     **********************************************************************/
    static inline Math::real arcsecond() throw() { return arcminute() / 60; }

    /** \name Ellipsoid parameters
     **********************************************************************/
    ///@{
    /**
     * @return equatorial radius of WGS84 ellipsoid
     **********************************************************************/
    static inline Math::real WGS84_a() throw() { return 6378137 * meter(); }
    /**
     * @return reciprocal flattening of WGS84 ellipsoid
     **********************************************************************/
    static inline Math::real WGS84_r() throw() { return real(298.257223563L); }
    /**
     * @return central scale factor for UTM
     **********************************************************************/
    static inline Math::real UTM_k0() throw() {return real(0.9996L); }
    /**
     * @return central scale factor for UPS
     **********************************************************************/
    static inline Math::real UPS_k0() throw() { return real(0.994L); }
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
    static inline Math::real meter() throw() { return real(1); }
    /**
     * @return the number of meters in a kilometer.
     **********************************************************************/
    static inline Math::real kilometer() throw() { return 1000 * meter(); }
    ///@}

    /**
     * @return the number of meters in a nautical mile (approximately 1 arc
     *   minute)
     **********************************************************************/
    static inline Math::real nauticalmile() throw() { return 1852 * meter(); }

    /** \name Anachronistic British units
     **********************************************************************/
    ///@{
    /**
     * @return the number of meters in an international foot.
     **********************************************************************/
    static inline Math::real foot() throw()
    { return real(0.0254L) * 12 * meter(); }
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
    { return real(1200) / real(3937) * meter(); }
    ///@}
  };

  /**
   * \brief %Exception handling for %GeographicLib
   *
   * A class to handle exceptions.  It's derived off std::runtime_error so it
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

#endif

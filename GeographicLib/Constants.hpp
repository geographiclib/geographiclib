/**
 * \file Constants.hpp
 * \brief Header for GeographicLib::Constants class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CONSTANTS_HPP)
#define GEOGRAPHICLIB_CONSTANTS_HPP "$Id$"

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
{ char VAR_ ## x [] __attribute__((unused)) = x; }
#else
/**
 * Insertion of RCS Id strings into the object file.
 **********************************************************************/
#define RCSID_DECL(x) namespace { char VAR_ ## x [] = x; }
#endif

RCSID_DECL(GEOGRAPHICLIB_CONSTANTS_HPP)

#include <cmath>
#include <cstdlib>
#if defined(_MSC_VER)
#include <float.h>		// For _finite
#endif

namespace GeographicLib {

  /**
   * \brief Mathematic functions needed by %GeographicLib
   *
   * Define mathematical functions in a way to localize system dependencies.
   * In addition define a real type to be used by %GeographicLib.
   **********************************************************************/
  class Math {
  private:
    Math();			// Disable constructor
  public:
    
    /**
     * The real type for %GeographicLib. Nearly all the testing has been done
     * with \e real_t = double.  However, the algorithms should also work with
     * float and long double (where available).
     **********************************************************************/
    typedef double real_t;

#if !defined(_MSC_VER)
    /**
     * sqrt(\e x<sup>2</sup> + \e y<sup>2</sup>)
     **********************************************************************/
    static inline real_t hypot(real_t x, real_t y) throw()
    { return ::hypot(x, y); }
    /**
     * asinh(\e x)
     **********************************************************************/
    static inline real_t asinh(real_t x) throw() { return ::asinh(x); }
    /**
     * atanh(\e x)
     **********************************************************************/
    static inline real_t atanh(real_t x) throw() { return ::atanh(x); }
    /**
     * \e x<sup>1/3</sup>
     **********************************************************************/
    static inline real_t cbrt(real_t x) throw() { return ::cbrt(x); }
    /**
     * Is \e x a regular number
     **********************************************************************/
    static inline int isfinite(real_t x) throw() { return std::isfinite(x); }
    /**
     * Convert a string to a real number
     **********************************************************************/
    static inline real_t strtod(const char *nptr, char **endptr) {
      return real_t( sizeof(real_t) == sizeof(double) ?
		     std::strtod(nptr, endptr) :
		     sizeof(real_t) == sizeof(float) ?
		     ::strtof(nptr, endptr) :
		     ::strtold(nptr, endptr) );
    }
#else
    static inline real_t hypot(real_t x, real_t y) throw() {
      return real_t( sizeof(real_t) == sizeof(double) ?
		     _hypot(x, y) : _hypot(x, y) );
    }
    // These have poor relative accuracy near x = 0.  However, for mapping
    // applications, we only need good absolute accuracy.
    // For small arguments we would write
    //
    // asinh(x) = asin(x) -x^3/3-5*x^7/56-63*x^11/1408-143*x^15/5120 ...
    // atanh(x) = atan(x) +2*x^3/3+2*x^7/7+2*x^11/11+2*x^15/15
    //
    // The accuracy of asinh is also bad for large negative arguments.  This is
    // easy to fix in the definition of asinh.  Instead we call these functions
    // with positive arguments and enforce the correct parity separately.
    static inline real_t asinh(real_t x) throw() {
      return std::log(x + std::sqrt(1 + x * x));
    }
    static inline real_t atanh(real_t x) throw() {
      return std::log((1 + x)/(1 - x))/2;
    }
    static inline real_t cbrt(real_t x) throw() {
      real_t y = std::pow(std::abs(x), 1/real_t(3));
      return x < 0 ? -y : y;
    }
    static inline int isfinite(real_t x) throw() { return _finite(x); }
    static inline real_t strtod(const char *nptr, char **endptr) {
      return real_t( std::strtod(nptr, endptr) );
    }
#endif
  };

  /**
   * \brief %Constants needed by %GeographicLib
   *
   * Define constants specifying the WGS84 ellipsoid, the UTM and UPS
   * projections, and various unit conversions.
   **********************************************************************/
  class Constants {
  private:
    typedef Math::real_t real_t;
    Constants();		// Disable constructor

  public:
    /**
     * pi
     **********************************************************************/
    static inline Math::real_t pi() throw()
    // good for about 123-bit accuracy
    { return real_t(3.141592653589793238462643383279502884L); }
    /**
     * Factor to convert from degrees to radians
     **********************************************************************/
    static inline Math::real_t degree() throw() { return pi() / 180; }

    /** \name Ellipsoid parameters
     **********************************************************************/
    ///@{
    /**
     * Major radius of WGS84 ellipsoid
     **********************************************************************/
    static inline Math::real_t WGS84_a() throw() { return 6378137 * meter(); }
    /**
     * Reciprocal flattening of WGS84 ellipsoid
     **********************************************************************/
    static inline Math::real_t WGS84_r() throw()
    { return real_t(298.257223563L); }
    /**
     * Central scale factor for UTM
     **********************************************************************/
    static inline Math::real_t UTM_k0() throw() {return real_t(0.9996L); }
    /**
     * Central scale factor for UPS
     **********************************************************************/
    static inline Math::real_t UPS_k0() throw() { return real_t(0.994L); }
    ///@}

    /** \name SI units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from meters to meters (i.e., 1, but this lets the
     * internal system of units be changed if necessary).
     **********************************************************************/
    static inline Math::real_t meter() throw() { return 1; }
    /**
     * Factor to convert from kilometers to meters.
     **********************************************************************/
    static inline Math::real_t kilometer() throw() { return 1000 * meter(); }
    ///@}

    /**
     * Factor to convert from nautical miles (approximately 1 arc minute) to
     * meters.
     **********************************************************************/
    static inline Math::real_t nauticalmile() throw() { return 1852 * meter(); }

    /** \name Anachronistic British units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from international feet to meters.
     **********************************************************************/
    static inline Math::real_t foot() throw()
    { return real_t(0.0254L) * 12 * meter(); }
    /**
     * Factor to convert from yards to meters.
     **********************************************************************/
    static inline Math::real_t yard() throw() { return 3 * foot(); }
    /**
     * Factor to convert from fathoms to meters.
     **********************************************************************/
    static inline Math::real_t fathom() throw() { return 2 * yard(); }
    /**
     * Factor to convert from chains to meters.
     **********************************************************************/
    static inline Math::real_t chain() throw() { return 22 * yard(); }
    /**
     * Factor to convert from furlongs to meters.
     **********************************************************************/
    static inline Math::real_t furlong() throw() { return 10 * chain(); }
    /**
     * Factor to convert from statute miles to meters.
     **********************************************************************/
    static inline Math::real_t mile() throw() { return 8 * furlong(); }
    ///@}

    /** \name Anachronistic US units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from US survery feet to meters.
     **********************************************************************/
    static inline Math::real_t surveyfoot() throw()
    { return real_t(1200) / real_t(3937) * meter(); }
    ///@}
  };

} // namespace GeographicLib

#endif

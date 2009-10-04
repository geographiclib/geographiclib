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
namespace GeographicLib {

  /**
   * \brief Mathematic functions needed by %GeographicLib
   *
   * Define mathematical functions in a way to localize system dependencies.
   * In addition define a real type to be used by %GeographicLib.
   **********************************************************************/
  class Math {
  private:
    Math();                     // Disable constructor
  public:
    
    /**
     * The real type for %GeographicLib. Nearly all the testing has been done
     * with \e real = double.  However, the algorithms should also work with
     * float and long double (where available).
     **********************************************************************/
    typedef double real;

#if !defined(_MSC_VER)
    /**
     * sqrt(\e x<sup>2</sup> + \e y<sup>2</sup>)
     **********************************************************************/
    static inline real hypot(real x, real y) throw()
    { return ::hypot(x, y); }
    /**
     * asinh(\e x)
     **********************************************************************/
    static inline real asinh(real x) throw() { return ::asinh(x); }
    /**
     * atanh(\e x)
     **********************************************************************/
    static inline real atanh(real x) throw() { return ::atanh(x); }
    /**
     * \e x<sup>1/3</sup>
     **********************************************************************/
    static inline real cbrt(real x) throw() { return ::cbrt(x); }
#else
    static inline real hypot(real x, real y) throw() {
      return real( sizeof(real) <= sizeof(double) ?
                   _hypot(x, y) : _hypotf(float(x), float(y)) );
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
    static inline real asinh(real x) throw() {
      return std::log(x + std::sqrt(1 + x * x));
    }
    static inline real atanh(real x) throw() {
      return std::log((1 + x)/(1 - x))/2;
    }
    static inline real cbrt(real x) throw() {
      real y = std::pow(std::abs(x), 1/real(3));
      return x < 0 ? -y : y;
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
    typedef Math::real real;
    Constants();                // Disable constructor

  public:
    /**
     * pi
     **********************************************************************/
    static inline Math::real pi() throw()
    // good for about 123-bit accuracy
    { return real(3.141592653589793238462643383279502884L); }
    /**
     * Factor to convert from degrees to radians
     **********************************************************************/
    static inline Math::real degree() throw() { return pi() / 180; }

    /** \name Ellipsoid parameters
     **********************************************************************/
    ///@{
    /**
     * Major radius of WGS84 ellipsoid
     **********************************************************************/
    static inline Math::real WGS84_a() throw() { return 6378137 * meter(); }
    /**
     * Reciprocal flattening of WGS84 ellipsoid
     **********************************************************************/
    static inline Math::real WGS84_r() throw() { return real(298.257223563L); }
    /**
     * Central scale factor for UTM
     **********************************************************************/
    static inline Math::real UTM_k0() throw() {return real(0.9996L); }
    /**
     * Central scale factor for UPS
     **********************************************************************/
    static inline Math::real UPS_k0() throw() { return real(0.994L); }
    ///@}

    /** \name SI units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from meters to meters (i.e., 1, but this lets the
     * internal system of units be changed if necessary).
     **********************************************************************/
    static inline Math::real meter() throw() { return real(1); }
    /**
     * Factor to convert from kilometers to meters.
     **********************************************************************/
    static inline Math::real kilometer() throw() { return 1000 * meter(); }
    ///@}

    /**
     * Factor to convert from nautical miles (approximately 1 arc minute) to
     * meters.
     **********************************************************************/
    static inline Math::real nauticalmile() throw() { return 1852 * meter(); }

    /** \name Anachronistic British units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from international feet to meters.
     **********************************************************************/
    static inline Math::real foot() throw()
    { return real(0.0254L) * 12 * meter(); }
    /**
     * Factor to convert from yards to meters.
     **********************************************************************/
    static inline Math::real yard() throw() { return 3 * foot(); }
    /**
     * Factor to convert from fathoms to meters.
     **********************************************************************/
    static inline Math::real fathom() throw() { return 2 * yard(); }
    /**
     * Factor to convert from chains to meters.
     **********************************************************************/
    static inline Math::real chain() throw() { return 22 * yard(); }
    /**
     * Factor to convert from furlongs to meters.
     **********************************************************************/
    static inline Math::real furlong() throw() { return 10 * chain(); }
    /**
     * Factor to convert from statute miles to meters.
     **********************************************************************/
    static inline Math::real mile() throw() { return 8 * furlong(); }
    ///@}

    /** \name Anachronistic US units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from US survery feet to meters.
     **********************************************************************/
    static inline Math::real surveyfoot() throw()
    { return real(1200) / real(3937) * meter(); }
    ///@}
  };

} // namespace GeographicLib

#endif

/**
 * \file Constants.hpp
 * \brief Header for GeographicLib::Constants class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(CONSTANTS_HPP)
#define CONSTANTS_HPP "$Id$"

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

RCSID_DECL(CONSTANTS_HPP)

namespace GeographicLib {

  /**
   * \brief %Constants needed by %GeographicLib
   *
   * Define constants specifying the WGS84 ellipsoid, the UTM and UPS
   * projections, and various unit conversions.
   **********************************************************************/
  class Constants {
  public:
    /**
     * pi
     **********************************************************************/
    static inline double pi() throw()
    // good for about 123-bit accuracy
    { return 3.141592653589793238462643383279502884; }
    /**
     * Factor to convert from degrees to radians
     **********************************************************************/
    static inline double degree() throw() { return pi() / 180; }

    /** \name Ellipsoid parameters
     **********************************************************************/
    ///@{
    /**
     * Major radius of WGS84 ellipsoid
     **********************************************************************/
    static inline double WGS84_a() throw() { return 6378137 * meter(); }
    /**
     * reciprocal flattening of WGS84 ellipsoid
     **********************************************************************/
    static inline double WGS84_r() throw() { return 298.257223563; }
    /**
     * Central scale factor for UTM
     **********************************************************************/
    static inline double UTM_k0() throw() {return 0.9996; }
    /**
     * Central scale factor for UPS
     **********************************************************************/
    static inline double UPS_k0() throw() { return 0.994; }
    ///@}

    /** \name SI units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from meters to meters (i.e., 1, but this lets the
     * internal system of units be changed if necessary).
     **********************************************************************/
    static inline double meter() throw() { return 1; }
    /**
     * Factor to convert from kilometers to meters.
     **********************************************************************/
    static inline double kilometer() throw() { return 1000 * meter(); }
    ///@}

    /**
     * Factor to convert from nautical miles (approximately 1 arc minute) to
     * meters.
     **********************************************************************/
    static inline double nauticalmile() throw() { return 1852 * meter(); }

    /** \name Anachronistic British units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from international feet to meters.
     **********************************************************************/
    static inline double foot() throw() { return 0.0254 * 12 * meter(); }
    /**
     * Factor to convert from yards to meters.
     **********************************************************************/
    static inline double yard() throw() { return 3 * foot(); }
    /**
     * Factor to convert from fathoms to meters.
     **********************************************************************/
    static inline double fathom() throw() { return 2 * yard(); }
    /**
     * Factor to convert from chains to meters.
     **********************************************************************/
    static inline double chain() throw() { return 22 * yard(); }
    /**
     * Factor to convert from furlongs to meters.
     **********************************************************************/
    static inline double furlong() throw() { return 10 * chain(); }
    /**
     * Factor to convert from statute miles to meters.
     **********************************************************************/
    static inline double mile() throw() { return 8 * furlong(); }
    ///@}

    /** \name Anachronistic US units
     **********************************************************************/
    ///@{
    /**
     * Factor to convert from US survery feet to meters.
     **********************************************************************/
    static inline double surveyfoot() throw()
    { return 1200.0 / 3937.0 * meter(); }
    ///@}
  };

} // namespace GeographicLib

#endif

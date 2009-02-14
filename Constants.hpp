/**
 * \file Constants.hpp
 * \brief Header for GeographicLib::Constants class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(CONSTANTS_HPP)
#define CONSTANTS_HPP "$Id$"

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
    static const double pi;
    /**
     * Conversion from degrees to radians
     **********************************************************************/
    static const double degree;

    /** \name Ellipsoid parameters
     **********************************************************************/
    ///@{
    /**
     * Major radius of WGS84 ellipsoid
     **********************************************************************/
    static const double WGS84_a;
    /**
     * Inverse flattening of WGS84 ellipsoid
     **********************************************************************/
    static const double WGS84_invf;
    /**
     * Central scale factor for UTM
     **********************************************************************/
    static const double UTM_k0;
    /**
     * Central scale factor for UPS
     **********************************************************************/
    static const double UPS_k0;
    ///@}

    /** \name SI units
     **********************************************************************/
    ///@{
    /**
     * Convert meters to meters (i.e., 1, but this lets the internal
     * system of units be changed if necessary).
     **********************************************************************/
    static const double meter;
    /**
     * Convert kilometers to meters.
     **********************************************************************/
    static const double kilometer;
    ///@}

    /**
     * Convert nautical miles (approximately 1 arc minute) to meters.
     **********************************************************************/
    static const double nauticalmile;

    /** \name Anachronistic British units
     **********************************************************************/
    ///@{
    /**
     * Convert international feet to meters.
     **********************************************************************/
    static const double foot;
    /**
     * Convert yards to meters.
     **********************************************************************/
    static const double yard;
    /**
     * Convert fathoms to meters.
     **********************************************************************/
    static const double fathom;
    /**
     * Convert chains to meters.
     **********************************************************************/
    static const double chain;
    /**
     * Convert statute miles to meters.
     **********************************************************************/
    static const double mile;
    ///@}

    /** \name Anachronistic US units
     **********************************************************************/
    ///@{
    /**
     * Convert US survery feet to meters.
     **********************************************************************/
    static const double surveyfoot;
    ///@}
  };

} // namespace GeographicLib

#endif

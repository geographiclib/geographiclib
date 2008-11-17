/**
 * \file Constants.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(CONSTANTS_HPP)
#define CONSTANTS_HPP "$Id$"

namespace GeographicLib {

  /**
   * \brief Constants needed by GeographicLib
   **********************************************************************/
  class Constants {
  public:
    static const double pi, degree, huge, WGS84_a, WGS84_invf, UPS_k0, UTM_k0;
    static const double meter, kilometer, // SI units
      nauticalmile,			  // Approx 1 arcminute
      foot, yard, fathom, chain, mile,	  // Anachronistic British units
      surveyfoot;			  // Anachronistic US unit
  };

} // namespace GeographicLib

#endif

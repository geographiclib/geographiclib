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
  };

} // namespace GeographicLib

#endif

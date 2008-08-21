/**
 * \file Constants.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/Constants.hpp"
#include <cmath>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = CONSTANTS_HPP;
}

namespace GeographicLib {

#if defined(M_PI)
  const double Constants::pi = M_PI;
#else
  const double Constants::pi = atan2(0.0, -1.0);
#endif
  const double Constants::degree = Constants::pi / 180;

  const double Constants::WGS84_a = 6378137.0;
  const double Constants::WGS84_invf = 298.257223563;
  const double Constants::UPS_k0 = 0.994;
  const double Constants::UTM_k0 = 0.9996;

} // namespace GeographicLib

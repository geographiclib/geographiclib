/**
 * \file Constants.cpp
 * \brief Implementation for GeographicLib::Constants class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/Constants.hpp"
#include <cmath>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = CONSTANTS_HPP;
}

namespace GeographicLib {

  using namespace std;

#if defined(M_PI)
  const double Constants::pi = M_PI;
#else
  const double Constants::pi = atan2(0.0, -1.0);
#endif
  const double Constants::degree = pi / 180;

  // All these constants are exact

  const double Constants::meter = 1.0;
  const double Constants::WGS84_a = 6378137.0 * meter;
  const double Constants::WGS84_invf = 298.257223563;
  const double Constants::UTM_k0 = 0.9996;
  const double Constants::UPS_k0 = 0.994;

  const double Constants::kilometer = 1000.0 * meter;
  const double Constants::nauticalmile = 1.852 * kilometer;
  const double Constants::foot = 0.0254 * 12 * meter;
  const double Constants::yard = 3 * foot;
  const double Constants::fathom = 2 * yard;
  const double Constants::chain = 22 * yard;
  const double Constants::mile = 1760 * yard;
  const double Constants::surveyfoot = 1200.0 / 3937.0 * meter;
} // namespace GeographicLib

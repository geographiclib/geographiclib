/**
 * \file AzimuthalEquidistant.cpp
 * \brief Implementation for GeographicLib::AzimuthalEquidistant class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/AzimuthalEquidistant.hpp"
#include "GeographicLib/Constants.hpp"
#include <cmath>
#include <limits>

#define AZIMUTHALEQUIDISTANT_CPP "$Id$"

RCSID_DECL(AZIMUTHALEQUIDISTANT_CPP)
RCSID_DECL(AZIMUTHALEQUIDISTANT_HPP)

namespace GeographicLib {

  using namespace std;

  const double AzimuthalEquidistant::eps =
    0.01 * sqrt(numeric_limits<double>::min());

  void AzimuthalEquidistant::Forward(double lat0, double lon0,
				     double lat, double lon,
				     double& x, double& y,
				     double& azi, double& rk) const throw() {
    double sig, s, azi0, m;
    sig = _earth.Inverse(lat0, lon0, lat, lon, s, azi0, azi, m);
    azi0 *= Constants::degree();
    x = s * sin(azi0);
    y = s * cos(azi0);
    rk = sig > eps ? m / s : 1;
  }

  void AzimuthalEquidistant::Reverse(double lat0, double lon0,
				     double x, double y,
				     double& lat, double& lon,
				     double& azi, double& rk) const throw() {
    double
      azi0 = atan2(x, y) / Constants::degree(),
      s = hypot(x, y);
    double sig, m;
    sig = _earth.Direct(lat0, lon0, azi0, s, lat, lon, azi, m);
    rk = sig > eps ? m / s : 1;
  }

} // namespace GeographicLib

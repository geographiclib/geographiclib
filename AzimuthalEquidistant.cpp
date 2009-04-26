/**
 * \file AzimuthalEquidistant.cpp
 * \brief Implementation for GeographicLib::AzimuthalEquidistant class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/AzimuthalEquidistant.hpp"
#include "GeographicLib/Constants.hpp"
#include <cmath>
#include <stdexcept>

#define AZIMUTHALEQUIDISTANT_CPP "$Id$"

RCSID_DECL(AZIMUTHALEQUIDISTANT_CPP)
RCSID_DECL(AZIMUTHALEQUIDISTANT_HPP)

namespace GeographicLib {

  using namespace std;

  void AzimuthalEquidistant::Reset(double lat0, double lon0) throw() {
    _lat0 = lat0;
    _lon0 = lon0 >= 180 ? lon0 - 360 : lon0 < -180 ? lon0 + 360 : lon0;
  }

  void AzimuthalEquidistant::Forward(double lat, double lon,
				     double& x, double& y,
				     double& azi, double& m) const throw() {
    double s12, azi1;
    _earth.Inverse(_lat0, _lon0, lat, lon, s12, azi1, azi, m);
    azi1 *= Constants::degree();
    x = s12 * sin(azi1);
    y = s12 * cos(azi1);
  }

  void AzimuthalEquidistant::Reverse(double x, double y,
				     double& lat, double& lon,
				     double& azi, double& m) const throw() {
    double
      azi1 = atan2(x, y) / Constants::degree(),
      s12 = hypot(x, y);
    _earth.Direct(_lat0, _lon0, azi1, s12, lat, lon, azi, m);
  }

} // namespace GeographicLib

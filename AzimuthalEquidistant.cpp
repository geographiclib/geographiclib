/**
 * \file AzimuthalEquidistant.cpp
 * \brief Implementation for GeographicLib::AzimuthalEquidistant class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/AzimuthalEquidistant.hpp"
#include <limits>

#define GEOGRAPHICLIB_AZIMUTHALEQUIDISTANT_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_AZIMUTHALEQUIDISTANT_CPP)
RCSID_DECL(GEOGRAPHICLIB_AZIMUTHALEQUIDISTANT_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real_t AzimuthalEquidistant::eps =
    real_t(0.01L) * sqrt(numeric_limits<real_t>::min());

  void AzimuthalEquidistant::Forward(real_t lat0, real_t lon0,
                                     real_t lat, real_t lon,
                                     real_t& x, real_t& y,
                                     real_t& azi, real_t& rk) const throw() {
    real_t sig, s, azi0, m;
    sig = _earth.Inverse(lat0, lon0, lat, lon, s, azi0, azi, m);
    azi0 *= Constants::degree();
    x = s * sin(azi0);
    y = s * cos(azi0);
    rk = sig > eps ? m / s : 1;
  }

  void AzimuthalEquidistant::Reverse(real_t lat0, real_t lon0,
                                     real_t x, real_t y,
                                     real_t& lat, real_t& lon,
                                     real_t& azi, real_t& rk) const throw() {
    real_t
      azi0 = atan2(x, y) / Constants::degree(),
      s = Math::hypot(x, y);
    real_t sig, m;
    sig = _earth.Direct(lat0, lon0, azi0, s, lat, lon, azi, m);
    rk = sig > eps ? m / s : 1;
  }

} // namespace GeographicLib

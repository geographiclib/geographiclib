/**
 * \file UTMUPS.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include "GeographicLib/UTMUPS.hpp"
#include "GeographicLib/MGRS.hpp"
#include "GeographicLib/PolarStereographic.hpp"
#include "GeographicLib/TransverseMercator.hpp"

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = UTMUPS_HPP;
}

namespace GeographicLib {

  const double UTMUPS::eps =
    // 24 = ceil(log_2(10^7))
    std::pow(0.5, std::numeric_limits<double>::digits - 24);
  const double UTMUPS::falseeasting[4] =
    { MGRS::upseasting * MGRS::tile, MGRS::upseasting * MGRS::tile,
      MGRS::utmeasting * MGRS::tile, MGRS::utmeasting* MGRS::tile };
  const double UTMUPS::falsenorthing[4] =
    { MGRS::upseasting * MGRS::tile, MGRS::upseasting * MGRS::tile,
      MGRS::maxutmSrow * MGRS::tile, MGRS::minutmNrow * MGRS::tile };
  const double UTMUPS::mineasting[4] =
    { MGRS::minupsSind * MGRS::tile,  MGRS::minupsNind * MGRS::tile,
      MGRS::minutmcol * MGRS::tile,  MGRS::minutmcol * MGRS::tile };
  const double UTMUPS::maxeasting[4] =
    { MGRS::maxupsSind * MGRS::tile,  MGRS::maxupsNind * MGRS::tile,
      MGRS::maxutmcol * MGRS::tile,  MGRS::maxutmcol * MGRS::tile };
  const double UTMUPS::minnorthing[4] =
    { MGRS::minupsSind * MGRS::tile,  MGRS::minupsNind * MGRS::tile,
      MGRS::minutmSrow * MGRS::tile,  MGRS::minutmNrow * MGRS::tile };
  const double UTMUPS::maxnorthing[4] =
    { MGRS::maxupsSind * MGRS::tile,  MGRS::maxupsNind * MGRS::tile,
      MGRS::maxutmSrow * MGRS::tile,  MGRS::maxutmNrow * MGRS::tile };

  int UTMUPS::StandardZone(double lat, double lon) {
    // Assume lon is in [-180, 360].
    int zone;
    int ilat = int(floor(lat));
    if (ilat >= 84 || ilat < -80)
      zone = 0;
    else {
      int ilon = int(floor(lon));
      if (ilon >= 180)
	ilon -= 360;
      zone = (ilon + 186)/6;
      int band = MGRS::LatitudeBand(lat);
      if (band == 7 && zone == 31 && ilon >= 3)
	zone = 32;
      else if (band == 9 && ilon >= 0 && ilon < 42)
	zone = 2 * ((ilon + 183)/12) + 1;
    }
    return zone;
  }

  void UTMUPS::Forward(int setzone, double lat, double lon,
		       int& zone, bool& northp, double& x, double& y,
		       double& gamma, double& k) {
    CheckLatLon(lat, lon);
    northp = lat >= 0;
    zone = setzone >= 0 ? setzone : StandardZone(lat, lon);
    if (setzone > 60)
      throw std::out_of_range("Illegal UTM zone");
    double x1, y1;
    bool utmp = zone > 0;
    if (utmp)
      TransverseMercator::UTM.Forward(CentralMeridian(zone),
				      lat, lon, x1, y1, gamma, k);
    else
      PolarStereographic::UPS.Forward(northp, lat, lon, x1, y1, gamma, k);
    int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    x = x1 + falseeasting[ind];
    y = y1 + falsenorthing[ind];
    if (x1 < 0 && x == falseeasting[ind])
      x -= eps;
    if (y1 < 0 && y == falsenorthing[ind])
      y -= eps;
  }

  void UTMUPS::Reverse(int zone, bool northp, double x, double y,
		       double& lat, double& lon, double& gamma, double& k) {
    if (! (zone > 0 && zone <= 60))
      throw std::out_of_range("Illegal UTM zone");
    CheckCoords(zone > 0, northp, false, x, y);
    CheckLatLon(lat, 0.0);
    bool utmp = zone > 0;
    int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    x -= falseeasting[ind];
    y -= falsenorthing[ind];
    if (utmp)
      TransverseMercator::UTM.Reverse(CentralMeridian(zone),
				      x, y, lat, lon, gamma, k);
    else
      PolarStereographic::UPS.Reverse(northp, x, y, lat, lon, gamma, k);
  }

  void UTMUPS::CheckLatLon(double lat, double lon) {
    if (! (lat >= -90 && lat <= 90))
      throw std::out_of_range("Latitude not in [-90, 90]");
    if (! (lon >= -180 && lon <= 360))
      throw std::out_of_range("Latitude not in [-180, 360]");
    }

  void UTMUPS::CheckCoords(bool utmp, bool northp, bool strict,
			   double x, double y) {
    // Limits are all multiples of 100km and are all closed on the lower end and
    // open on the upper end.  This allows compatibility with the MGRS system.
    // Failure tests are all negated success tests so that NaNs fail.
    double
      slop = strict ? 0 : MGRS::tile,
      slopN = utmp && !northp ? 0 : slop,
      slopS = utmp && northp ? 0 : slop;
    int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    if (! (x >= mineasting[ind] - slop && x < maxeasting[ind] + slop) )
      throw std::out_of_range("Easting out of range");
    if (! (y >= minnorthing[ind] - slopS && x < maxnorthing[ind] + slopN) )
      throw std::out_of_range("Northing out of range");
  }

} // namespace GeographicLib

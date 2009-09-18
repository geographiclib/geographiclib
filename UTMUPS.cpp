/**
 * \file UTMUPS.cpp
 * \brief Implementation for GeographicLib::UTMUPS class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/UTMUPS.hpp"
#include "GeographicLib/MGRS.hpp"
#include "GeographicLib/PolarStereographic.hpp"
#include "GeographicLib/TransverseMercator.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <iomanip>

#define GEOGRAPHICLIB_UTMUPS_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_UTMUPS_CPP)
RCSID_DECL(GEOGRAPHICLIB_UTMUPS_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real_t UTMUPS::falseeasting[4] =
    { MGRS::upseasting * MGRS::tile, MGRS::upseasting * MGRS::tile,
      MGRS::utmeasting * MGRS::tile, MGRS::utmeasting* MGRS::tile };
  const Math::real_t UTMUPS::falsenorthing[4] =
    { MGRS::upseasting * MGRS::tile, MGRS::upseasting * MGRS::tile,
      MGRS::maxutmSrow * MGRS::tile, MGRS::minutmNrow * MGRS::tile };
  const Math::real_t UTMUPS::mineasting[4] =
    { MGRS::minupsSind * MGRS::tile,  MGRS::minupsNind * MGRS::tile,
      MGRS::minutmcol * MGRS::tile,  MGRS::minutmcol * MGRS::tile };
  const Math::real_t UTMUPS::maxeasting[4] =
    { MGRS::maxupsSind * MGRS::tile,  MGRS::maxupsNind * MGRS::tile,
      MGRS::maxutmcol * MGRS::tile,  MGRS::maxutmcol * MGRS::tile };
  const Math::real_t UTMUPS::minnorthing[4] =
    { MGRS::minupsSind * MGRS::tile,  MGRS::minupsNind * MGRS::tile,
      MGRS::minutmSrow * MGRS::tile,
      (MGRS::minutmNrow + MGRS::minutmSrow - MGRS::maxutmSrow) * MGRS::tile };
  const Math::real_t UTMUPS::maxnorthing[4] =
    { MGRS::maxupsSind * MGRS::tile,  MGRS::maxupsNind * MGRS::tile,
      (MGRS::maxutmSrow + MGRS::maxutmNrow - MGRS::minutmNrow) * MGRS::tile,
      MGRS::maxutmNrow * MGRS::tile };

  int UTMUPS::StandardZone(real_t lat, real_t lon)  throw() {
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

  void UTMUPS::Forward(real_t lat, real_t lon,
		       int& zone, bool& northp, real_t& x, real_t& y,
		       real_t& gamma, real_t& k,
		       int setzone) {
    CheckLatLon(lat, lon);
    northp = lat >= 0;
    zone = setzone >= 0 ? setzone : StandardZone(lat, lon);
    if (setzone > 60)
      throw out_of_range("Illegal UTM zone requested " + setzone);
    real_t x1, y1;
    bool utmp = zone > 0;
    if (utmp) {
      real_t
	lon0 = CentralMeridian(zone),
	dlon = lon - lon0;
      dlon = abs(dlon - 360 * floor((dlon + 180)/360));
      if (dlon > 60)
	// Check isn't really necessary because CheckCoords catches this case.
	// But this allows a more meaningful error message to be given.
	throw out_of_range("Longitude " + str(lon)
				+ "d more than 60d from center of UTM zone "
				+ str(zone));
      TransverseMercator::UTM.Forward(lon0, lat, lon, x1, y1, gamma, k);
    } else {
      if (abs(lat) < 70)
	// Check isn't really necessary ... (see above).
	throw out_of_range("Latitude " + str(lat)
				+ "d more than 20d from "
				+ (northp ? "N" : "S") + " pole");
      PolarStereographic::UPS.Forward(northp, lat, lon, x1, y1, gamma, k);
    }
    int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    x = x1 + falseeasting[ind];
    y = y1 + falsenorthing[ind];
    if (! CheckCoords(zone > 0, northp, x, y, false) )
      throw out_of_range("Latitude " + str(lat) +
			      ", longitude " + str(lon) +
			      " out of legal range for " +
			      (utmp ? "UTM zone " + str(zone) : "UPS"));
  }

  void UTMUPS::Reverse(int zone, bool northp, real_t x, real_t y,
		       real_t& lat, real_t& lon, real_t& gamma, real_t& k) {
    if (! (zone >= 0 && zone <= 60))
      throw out_of_range("Illegal UTM zone " + str(zone));
    CheckCoords(zone > 0, northp, x, y);
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

  void UTMUPS::CheckLatLon(real_t lat, real_t lon) {
    if (! (lat >= -90 && lat <= 90))
      throw out_of_range("Latitude " + str(lat) +
			      "d not in [-90d, 90d]");
    if (! (lon >= -180 && lon <= 360))
      throw out_of_range("Latitude " + str(lon) +
			      "d not in [-180d, 360d]");
    }

  bool UTMUPS::CheckCoords(bool utmp, bool northp, real_t x, real_t y,
			   bool throwp) {
    // Limits are all multiples of 100km and are all closed on the both ends.
    // Failure tests are all negated success tests so that NaNs fail.
    real_t slop = MGRS::tile;
    int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    if (! (x >= mineasting[ind] - slop && x <= maxeasting[ind] + slop) ) {
      if (!throwp) return false;
      throw out_of_range("Easting " + str(x/1000)
			      + "km not in "
			      + (utmp ? "UTM" : "UPS") + " range for "
			      + (northp ? "N" : "S" )
			      + " hemisphere ["
			      + str((mineasting[ind] - slop)/1000) + "km, "
			      + str((maxeasting[ind] + slop)/1000) + "km]");
    }
    if (! (y >= minnorthing[ind] - slop && y <= maxnorthing[ind] + slop) ) {
      if (!throwp) return false;
      throw out_of_range("Northing " + str(y/1000)
			      + "km not in "
			      + (utmp ? "UTM" : "UPS") + " range for "
			      + (northp ? "N" : "S" )
			      + " hemisphere ["
			      + str((minnorthing[ind] - slop)/1000) + "km, "
			      + str((maxnorthing[ind] + slop)/1000) + "km]");
    }
    return true;
  }

  void UTMUPS::DecodeZone(const std::string& zonestr,
			  int& zone, bool& northp) {
    unsigned zlen = unsigned(zonestr.size());
    if (zlen == 0)
      throw out_of_range("Empty zone specification");
    if (zlen > 3)
      throw out_of_range("More than 3 characters in zone specification " +
			 zonestr);
    char hemi = toupper(zonestr[zlen - 1]);
    northp = hemi == 'N';
    if (! (northp || hemi == 'S'))
      throw out_of_range(string("Illegal hemisphere letter ") + hemi +
			 " in " + zonestr);
    if (zlen == 1)
      zone = 0;
    else {
      const char* c = zonestr.c_str();
      char* q;
      zone = strtol(c, &q, 10);
      if (q - c != int(zlen) - 1)
	throw out_of_range("Extra text in UTM/UPS zone " + zonestr);
      if (q > c && zone == 0)
	// Don't allow 0N as an alternative to N for UPS coordinates
	throw out_of_range("Illegal zone 0 in " + zonestr);
    }
  }

  std::string UTMUPS::EncodeZone(int zone, bool northp) {
    ostringstream os;
    if (! (zone >= 0 && zone <= 60)) {
      os << "Illegal UTM zone " << zone;
      throw out_of_range(os.str());
    }
    if (zone)
      os << setfill('0') << setw(2) << zone;
    os << (northp ? 'N' : 'S');
    return os.str();
  }

  Math::real_t UTMUPS::UTMShift() throw() { return real_t(MGRS::utmNshift); }

} // namespace GeographicLib

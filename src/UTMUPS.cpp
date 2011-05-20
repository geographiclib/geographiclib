/**
 * \file UTMUPS.cpp
 * \brief Implementation for GeographicLib::UTMUPS class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/UTMUPS.hpp"
#include "GeographicLib/MGRS.hpp"
#include "GeographicLib/PolarStereographic.hpp"
#include "GeographicLib/TransverseMercator.hpp"
#include <iomanip>

#define GEOGRAPHICLIB_UTMUPS_CPP "$Id: UTMUPS.cpp 6814 2010-02-04 22:33:41Z karney $"

RCSID_DECL(GEOGRAPHICLIB_UTMUPS_CPP)
RCSID_DECL(GEOGRAPHICLIB_UTMUPS_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real UTMUPS::falseeasting[4] =
    { MGRS::upseasting * MGRS::tile, MGRS::upseasting * MGRS::tile,
      MGRS::utmeasting * MGRS::tile, MGRS::utmeasting* MGRS::tile };
  const Math::real UTMUPS::falsenorthing[4] =
    { MGRS::upseasting * MGRS::tile, MGRS::upseasting * MGRS::tile,
      MGRS::maxutmSrow * MGRS::tile, MGRS::minutmNrow * MGRS::tile };
  const Math::real UTMUPS::mineasting[4] =
    { MGRS::minupsSind * MGRS::tile,  MGRS::minupsNind * MGRS::tile,
      MGRS::minutmcol * MGRS::tile,  MGRS::minutmcol * MGRS::tile };
  const Math::real UTMUPS::maxeasting[4] =
    { MGRS::maxupsSind * MGRS::tile,  MGRS::maxupsNind * MGRS::tile,
      MGRS::maxutmcol * MGRS::tile,  MGRS::maxutmcol * MGRS::tile };
  const Math::real UTMUPS::minnorthing[4] =
    { MGRS::minupsSind * MGRS::tile,  MGRS::minupsNind * MGRS::tile,
      MGRS::minutmSrow * MGRS::tile,
      (MGRS::minutmNrow + MGRS::minutmSrow - MGRS::maxutmSrow) * MGRS::tile };
  const Math::real UTMUPS::maxnorthing[4] =
    { MGRS::maxupsSind * MGRS::tile,  MGRS::maxupsNind * MGRS::tile,
      (MGRS::maxutmSrow + MGRS::maxutmNrow - MGRS::minutmNrow) * MGRS::tile,
      MGRS::maxutmNrow * MGRS::tile };

  int UTMUPS::StandardZone(real lat, real lon, int setzone) {
    if (setzone < MINPSEUDOZONE || setzone > MAXZONE)
      throw GeographicErr("Illegal zone requested " + str(setzone));
    if (setzone >= MINZONE)
      return setzone;
    // Assume lon is in [-180, 360].
    if (setzone == UTM || (lat >= -80 && lat < 84)) {
      // Assume lon is in [-180, 360].
      int ilon = int(floor(lon));
      if (ilon >= 180)
        ilon -= 360;
      int zone = (ilon + 186)/6;
      int band = MGRS::LatitudeBand(lat);
      if (band == 7 && zone == 31 && ilon >= 3)
        zone = 32;
      else if (band == 9 && ilon >= 0 && ilon < 42)
        zone = 2 * ((ilon + 183)/12) + 1;
      return zone;
    } else
      return UPS;
  }

  void UTMUPS::Forward(real lat, real lon,
                       int& zone, bool& northp, real& x, real& y,
                       real& gamma, real& k,
                       int setzone, bool mgrslimits) {
    CheckLatLon(lat, lon);
    bool northp1 = lat >= 0;
    int zone1 = StandardZone(lat, lon, setzone);
    real x1, y1, gamma1, k1;
    bool utmp = zone1 != UPS;
    if (utmp) {
      real
        lon0 = CentralMeridian(zone1),
        dlon = lon - lon0;
      dlon = abs(dlon - 360 * floor((dlon + 180)/360));
      if (dlon > 60)
        // Check isn't really necessary because CheckCoords catches this case.
        // But this allows a more meaningful error message to be given.
        throw GeographicErr("Longitude " + str(lon)
                            + "d more than 60d from center of UTM zone "
                            + str(zone1));
      TransverseMercator::UTM.Forward(lon0, lat, lon, x1, y1, gamma1, k1);
    } else {
      if (abs(lat) < 70)
        // Check isn't really necessary ... (see above).
        throw GeographicErr("Latitude " + str(lat) + "d more than 20d from "
                            + (northp1 ? "N" : "S") + " pole");
      PolarStereographic::UPS.Forward(northp1, lat, lon, x1, y1, gamma1, k1);
    }
    int ind = (utmp ? 2 : 0) + (northp1 ? 1 : 0);
    x1 += falseeasting[ind];
    y1 += falsenorthing[ind];
    if (! CheckCoords(zone1 != UPS, northp1, x1, y1, mgrslimits, false) )
      throw GeographicErr("Latitude " + str(lat) + ", longitude " + str(lon)
                          + " out of legal range for "
                          + (utmp ? "UTM zone " + str(zone1) : "UPS"));
    zone = zone1;
    northp = northp1;
    x = x1;
    y = y1;
    gamma = gamma1;
    k = k1;
  }

  void UTMUPS::Reverse(int zone, bool northp, real x, real y,
                       real& lat, real& lon, real& gamma, real& k,
                       bool mgrslimits) {
    if (! (zone >= MINZONE && zone <= MAXZONE))
      throw GeographicErr("Zone " + str(zone) + " not in range [0, 60]");
    bool utmp = zone != UPS;
    CheckCoords(utmp, northp, x, y, mgrslimits);
    int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    x -= falseeasting[ind];
    y -= falsenorthing[ind];
    if (utmp)
      TransverseMercator::UTM.Reverse(CentralMeridian(zone),
                                      x, y, lat, lon, gamma, k);
    else
      PolarStereographic::UPS.Reverse(northp, x, y, lat, lon, gamma, k);
  }

  void UTMUPS::CheckLatLon(real lat, real lon) {
    if (! (lat >= -90 && lat <= 90))
      throw GeographicErr("Latitude " + str(lat) + "d not in [-90d, 90d]");
    if (! (lon >= -180 && lon <= 360))
      throw GeographicErr("Latitude " + str(lon) + "d not in [-180d, 360d]");
    }

  bool UTMUPS::CheckCoords(bool utmp, bool northp, real x, real y,
                           bool mgrslimits, bool throwp) {
    // Limits are all multiples of 100km and are all closed on the both ends.
    // Failure tests are all negated success tests so that NaNs fail.
    real slop = mgrslimits ? 0 : MGRS::tile;
    int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    if (! (x >= mineasting[ind] - slop && x <= maxeasting[ind] + slop) ) {
      if (!throwp) return false;
      throw GeographicErr("Easting " + str(x/1000) + "km not in "
                          + (mgrslimits ? "MGRS/" : "")
                          + (utmp ? "UTM" : "UPS") + " range for "
                          + (northp ? "N" : "S" ) + " hemisphere ["
                          + str((mineasting[ind] - slop)/1000) + "km, "
                          + str((maxeasting[ind] + slop)/1000) + "km]");
    }
    if (! (y >= minnorthing[ind] - slop && y <= maxnorthing[ind] + slop) ) {
      if (!throwp) return false;
      throw GeographicErr("Northing " + str(y/1000) + "km not in "
                          + (mgrslimits ? "MGRS/" : "")
                          + (utmp ? "UTM" : "UPS") + " range for "
                          + (northp ? "N" : "S" ) + " hemisphere ["
                          + str((minnorthing[ind] - slop)/1000) + "km, "
                          + str((maxnorthing[ind] + slop)/1000) + "km]");
    }
    return true;
  }

  void UTMUPS::DecodeZone(const std::string& zonestr, int& zone, bool& northp) {
    unsigned zlen = unsigned(zonestr.size());
    if (zlen == 0)
      throw GeographicErr("Empty zone specification");
    if (zlen > 3)
      throw GeographicErr("More than 3 characters in zone specification "
                          + zonestr);
    char hemi = toupper(zonestr[zlen - 1]);
    bool northp1 = hemi == 'N';
    if (! (northp1 || hemi == 'S'))
      throw GeographicErr(string("Illegal hemisphere letter ") + hemi + " in "
                          + zonestr + ", specify N or S");
    if (zlen == 1)
      zone = UPS;
    else {
      const char* c = zonestr.c_str();
      char* q;
      int zone1 = strtol(c, &q, 10);
      if (q == c)
        throw GeographicErr("No zone number found in " + zonestr);
      if (q - c != int(zlen) - 1)
        throw GeographicErr("Extra text " +
                            zonestr.substr(q - c, int(zlen) - 1 - (q - c)) +
                            " in UTM/UPS zone " + zonestr);
      if (zone1 == UPS)
        // Don't allow 0N as an alternative to N for UPS coordinates
        throw GeographicErr("Illegal zone 0 in " + zonestr +
                            ", use just " + hemi + " for UPS");
      if (!(zone1 >= MINUTMZONE && zone1 <= MAXUTMZONE))
        throw GeographicErr("Zone " + str(zone1) + " not in range [1, 60]");
      zone = zone1;
    }
    northp = northp1;
  }

  std::string UTMUPS::EncodeZone(int zone, bool northp) {
    if (! (zone >= MINZONE && zone <= MAXZONE))
        throw GeographicErr("Zone " + str(zone) + " not in range [0, 60]");
    ostringstream os;
    if (zone != UPS)
      os << setfill('0') << setw(2) << zone;
    os << (northp ? 'N' : 'S');
    return os.str();
  }

  Math::real UTMUPS::UTMShift() throw() { return real(MGRS::utmNshift); }

} // namespace GeographicLib

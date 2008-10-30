/**
 * \file UTMUPS.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

// Convert between Geographic coordinates and UTM/UPS.
//
// Properties
//
// The conversions are closed, i.e., output from Forward is legal input for
// Reverse and vice versa.  The error is about 5nm in each direction.  However,
// the conversion from legal UTM/UPS coordinates to geographic coordinates and
// back might throw an error if the initial point is within 5nm of the edge of
// the allowed range for the UTM/UPS coordinates.
//
// The simplest way to guarantee the closed property is to define allowed
// ranges for the eastings and northings for UTM and UPS coordinates.  The UTM
// boundaries are the same for all zones.  (The only place the exceptional
// nature of the zone boundaries is evident is when converting to UTM/UPS
// coordinates requesting the standard zone.)  The MGRS lettering scheme
// imposes natural limits on UTM/UPS coordinates which may be converted into
// MGRS coordinates.  For the conversion to/from geographic coordinates these
// ranges have been extended by 100km in order to provide a generous overlap
// between UTM and UPS and between UTM zones.

#if !defined(UTMUPS_HPP)
#define UTMUPS_HPP "$Id$"

#include <string>
#include <sstream>

namespace GeographicLib {

  class UTMUPS {
  private:
    static const double falseeasting[4];
    static const double falsenorthing[4];
    static const double mineasting[4];
    static const double maxeasting[4];
    static const double minnorthing[4];
    static const double maxnorthing[4];
    static double CentralMeridian(int zone) { return 6 * zone - 183.0; }
    template<typename T> static std::string str(T x) {
      std::ostringstream s; s << x; return s.str();
    }
    static void CheckLatLon(double lat, double lon);
    // Throw an error if easting or northing are outside standard ranges.
    static void CheckCoords(bool utmp, bool northp, double x, double y);
  public:

    // Return the standard zone for a given latitude and longitude (degrees).
    // Return 0 if in the standard regions for UPS otherwise return the UTM
    // zone.  This includes the Norway and Svalbard exceptions.  The tests on
    // latitudes and longitudes are all closed on the lower end open on the
    // upper.  Thus for UTM zone 38, latitude is in [-80, 84) and longitude is
    // in [42, 48).  This is exact.
    static int StandardZone(double lat, double lon);

    // Convert latitude and longitude (degrees) to either UTM or UPS
    // coordinates.  If setzone < 0 then use standard zone (which includes the
    // selection of UTM vs UPS), else if setzone == 0 use UPS, else use
    // specified UTM zone.  zone and northp give resulting zone (with zone == 0
    // indicating UPS) and hemisphere, x and y are the returned easting and
    // northing (meters); gamma and k give convergence (degrees) and scale.
    // Throw error if the resulting easting or northing is outside the allowed
    // range (see Reverse).  The accuracy of the conversion is about 5nm.
    static void Forward(int setzone, double lat, double lon,
			int& zone, bool& northp, double& x, double& y,
			double& gamma, double& k);

    // Convert UTM or UPS coordinate to latitude and longitude (degrees).  zone
    // and northp give the zone (with zone == 0 indicating UPS) and hemisphere.
    // x and y are the easting and northing (meters).  Throw error if easting
    // or northing is outside the allowed range (see below).  This also returns
    // gamma and k (see Forward).  The accuracy of the conversion is about 5nm.
    //
    // UTM eastings are allowed to be in the range [0km, 1000km], northings are
    // allowed to be in in [0km, 9600km] for the northern hemisphere and in
    // [900km, 10000km] for the southern hemisphere.
    //
    // UPS eastings and northings are allowed to be in the range [1200km,
    // 2800km] in the northern hemisphere and in [700km, 3100km] in the
    // southern hemisphere.
    //
    // These ranges are 100km larger than allowed for the conversions to MGRS.
    // (100km is the maximum extra padding consistent with eastings remaining
    // non-negative.)  This allows generous overlaps between zones and UTM and
    // UPS.  No checks are performed beyond these (e.g., to limit the distance
    // outside the standard zone boundaries).
    static void Reverse(int zone, bool northp, double x, double y,
			double& lat, double& lon, double& gamma, double& k);

    // Forward without convergence and scale
    static void Forward(int setzone, double lat, double lon,
			int& zone, bool& northp, double& x, double& y) {
      double gamma, k;
      Forward(setzone, lat, lon, zone, northp, x, y, gamma, k);
    }

    // Reverse without convergence and scale
    static void Reverse(int zone, bool northp, double x, double y,
			double& lat, double& lon) {
      double gamma, k;
      Reverse(zone, northp, x, y, lat, lon, gamma, k);
    }
  };

} // namespace GeographicLib
#endif

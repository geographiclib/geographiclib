/**
 * \file UTMUPS.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(UTMUPS_HPP)
#define UTMUPS_HPP "$Id$"

namespace GeographicLib {

  class UTMUPS {
  private:
    // The smallest length s.t., 1.0e7 - eps < 1.07 (approx 1.9 nm)
    static const double eps;
    static const double falseeasting[4];
    static const double falsenorthing[4];
    static const double mineasting[4];
    static const double maxeasting[4];
    static const double minnorthing[4];
    static const double maxnorthing[4];
    static double CentralMeridian(int zone) { return 6 * zone - 183.0; }
  public:
    static int StandardZone(double lat, double lon);

    // Lat/Long -> UTM/UPS.  If setzone < 0 then use standard zone, else if
    // setzone == 0 use UPS else use specified UTM zone.  zone and northp give
    // resulting zone (0 => UPS) and hemisphere, x and y are returned easting
    // and northing; gamma and k give convergence and scale.
    static void Forward(int setzone, double lat, double lon,
			int& zone, bool& northp, double& x, double& y,
			double& gamma, double& k);
    // UTM/UPS -> Lat/Long reversing Forward.
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
    static void CheckLatLon(double lat, double lon);
    // Throw an error if easting or northing are outside standard ranges.  If
    // strict, restrict to MSGS subset.  Otherwise allow an extra 100km.
    static void CheckCoords(bool utmp, bool northp, bool strict,
			    double x, double y);
  };

} // namespace GeographicLib
#endif

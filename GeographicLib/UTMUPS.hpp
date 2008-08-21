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

    static void Forward(int setzone, double lat, double lon,
			int& zone, bool& northp, double& x, double& y);
    static void Reverse(int zone, bool northp, double x, double y,
			double& lat, double& lon);
    static void CheckCoords(bool utmp, bool northp, bool strict,
			    double x, double y);
  };

} // namespace GeographicLib
#endif

/**
 * \file MGRS.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(MGRS_HPP)
#define MGRS_HPP "$Id$"

#include <cmath>
#include <string>

namespace GeographicLib {
  class MGRS {
  private:
    static const std::string utmcols[3];
    static const std::string utmrow;
    static const std::string upscols[4];
    static const std::string upsrows[2];
    static const std::string latband;
    static const std::string upsband;
    static const std::string digits;
    static const int mineasting[4];
    static const int maxeasting[4];
    static const int minnorthing[4];
    static const int maxnorthing[4];
    enum {
      base = 10,
      // Top-level tiles are 10^5 m = 100km on a side
      tilelevel = 5,
      // Period of UTM row letters
      utmrowperiod = 20,
      // Row letters are shifted by 5 for even zones
      utmevenrowshift = 5,
    };
    static int lookup(const std::string& s, char c) {
      std::string::size_type r = s.find(c);
      return r == std::string::npos ? -1 : int(r);
    }
  public:
    static void Forward(int zone, bool northp, double x, double y,
			int prec, std::string& mgrs);
    // In case latitude is already known
    static void Forward(int zone, bool northp, double x, double y, double lat,
			int prec, std::string& mgrs);
    static void Reverse(const std::string& mgrs,
			int& zone, bool& northp, double& x, double& y,
			int& prec);
    static int LatitudeBand(double lat) {
      int ilat = int(floor(lat));
      return std::max(-10, std::min(9, (ilat + 80)/8 - 10));
    }
    static void CheckCoords(bool utmp, bool northp,
			    double x, double y);
    enum {
      tile = 100000,		// Size MGRS blocks
      minutmcol = 1,
      maxutmcol = 9,
      minutmSrow = 10,
      maxutmSrow = 100,		// Also used for UTM S false northing
      minutmNrow = 0,		// Also used for UTM N false northing
      maxutmNrow = 95,
      minupsSind = 8,		// These 4 ind's apply to easting and northing
      maxupsSind = 32,
      minupsNind = 13,
      maxupsNind = 27,
      upseasting = 20,		// Also used for UPS false northing
      utmeasting = 5,		// UTM false easting
    };
  };

} // namespace GeographicLib
#endif

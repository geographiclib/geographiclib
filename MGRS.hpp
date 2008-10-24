/**
 * \file MGRS.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

// Convert between UTM/UPS and MSGS.  Semantics is as follows:
//
// A MGRS coordinate represents a square area with size depending on the
// precision.  Forward returns the MGRS square enclosing the point.  Reverse
// returns the centerp ? (center of the square) : (SW corner of the square).
// By default, centerp is true.
//
// UTM eastings are allowed to be in the range [100km, 900km), northings are
// allowed to be in in [0km, 9500km) for the northern hemisphere and in
// [1000km, 10000km) for the southern hemisphere.
//
// UPS eastings/northings are allowed to be in the range [1300km, 2700km) in
// the northern hemisphere and in [800km, 3200km) in the southern hemisphere.
//
// These restrictions are dictated by the allowed letters in MGRS coordinates
// and allow plenty of overlap between zones and between UTM and UPS.
//
// The band letter for UTM-style MSRS coordinates is dictated by the latitude.
// For reverse lookups a neighboring band letter is permitted provided that the
// band letter is valid for some portion of the top-level 100km square.

#if !defined(MGRS_HPP)
#define MGRS_HPP "$Id$"

#include <cmath>
#include <string>
#include <sstream>

namespace GeographicLib {
  class MGRS {
  private:
    static const double eps;
    static const double angeps;
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
      // Maximum precision is um
      maxprec = 5 + 6,
    };
    static void CheckCoords(bool utmp, bool northp, double& x, double& y);
    static int lookup(const std::string& s, char c) {
      std::string::size_type r = s.find(toupper(c));
      return r == std::string::npos ? -1 : int(r);
    }
    template<typename T> static std::string str(T x) {
      std::ostringstream s; s << x; return s.str();
    }
    static int UTMRow(int iband, int icol, int irow);
  public:
    static void Forward(int zone, bool northp, double x, double y,
			int prec, std::string& mgrs);
    // In case latitude is already known.  Latitude is ignored for zone == 0
    // (UPS), otherwise resulting latitude band is checked for consistency.
    static void Forward(int zone, bool northp, double x, double y, double lat,
			int prec, std::string& mgrs);
    // If centerp (default), return center of square, else return SW corner
    static void Reverse(const std::string& mgrs,
			int& zone, bool& northp, double& x, double& y,
			int& prec, bool centerp = true);
    // Public because UTMUPS::StandardZone needs to use this
    static int LatitudeBand(double lat) {
      int ilat = int(floor(lat));
      return (std::max)(-10, (std::min)(9, (ilat + 80)/8 - 10));
    }
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

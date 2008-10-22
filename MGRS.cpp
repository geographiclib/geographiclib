/**
 * \file MGRS.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/MGRS.hpp"
#include "GeographicLib/UTMUPS.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = MGRS_HPP;
}

namespace GeographicLib {

  const double MGRS::eps =
    // 7 = ceil(log_2(90))
    std::pow(0.5, std::numeric_limits<double>::digits - 7);
  const std::string MGRS::utmcols[3] =
    { "ABCDEFGH", "JKLMNPQR", "STUVWXYZ" };
  const std::string MGRS::utmrow  = "ABCDEFGHJKLMNPQRSTUV";
  const std::string MGRS::upscols[4] =
    { "JKLPQRSTUXYZ", "ABCFGHJKLPQR", "RSTUXYZ", "ABCFGHJ" };
  const std::string MGRS::upsrows[2] =
    { "ABCDEFGHJKLMNPQRSTUVWXYZ", "ABCDEFGHJKLMNP" };
  const std::string MGRS::latband = "CDEFGHJKLMNPQRSTUVWX";
  const std::string MGRS::upsband = "ABYZ";
  const std::string MGRS::digits  = "0123456789";

  const int MGRS::mineasting[4] =
    { MGRS::minupsSind,  MGRS::minupsNind,
      MGRS::minutmcol,  MGRS::minutmcol };
  const int MGRS::maxeasting[4] =
    { MGRS::maxupsSind,  MGRS::maxupsNind,
      MGRS::maxutmcol,  MGRS::maxutmcol };
  const int MGRS::minnorthing[4] =
    { MGRS::minupsSind,  MGRS::minupsNind,
      MGRS::minutmSrow,  MGRS::minutmNrow };
  const int MGRS::maxnorthing[4] =
    { MGRS::maxupsSind,  MGRS::maxupsNind,
      MGRS::maxutmSrow,  MGRS::maxutmNrow };

  void MGRS::Forward(int zone, bool northp, double x, double y, double lat,
		     int prec, std::string& mgrs) {
    bool utmp = zone != 0;
    CheckCoords(utmp, northp, x, y);
    if (!(zone >= 0 || zone <= 60))
      throw std::out_of_range("Zone " + str(zone) + " not in [0,60]");
    if (!(prec >= 0 || prec <= maxprec))
      throw std::out_of_range("MGRS precision " + str(prec) + " not in [0, "
			      + str(int(maxprec)) + "]");
    int
      zone1 = zone - 1,
      z = utmp ? 2 : 0;
    // Space for zone, 3 block letters, easting + northing
    mgrs.resize(z + 3 + 2 * prec);
    if (utmp) {
      mgrs[0] = digits[ zone / base ];
      mgrs[1] = digits[ zone % base ];
      // This isn't necessary...!  Keep y non-neg
      // if (!northp) y -= maxutmSrow * tile;
    }
    int
      xh = int(floor(x)) / tile,
      yh = int(floor(y)) / tile;
    double
      xf = x - tile * xh,
      yf = y - tile * yh;
    if (utmp) {
      int
	// Correct fuzziness in latitude near equator
	iband = std::abs(lat) > eps ? LatitudeBand(lat) : (northp ? 0 : -1),
	icol = xh - minutmcol,
	irow = UTMRow(iband, icol, yh % utmrowperiod);
      if (irow != yh - (northp ? 0 : maxutmSrow))
	throw std::out_of_range("Latitude " + str(lat)
				+ " is inconsistent with UTM coordinates");
      mgrs[z++] = latband[10 + iband];
      mgrs[z++] = utmcols[zone1 % 3][icol];
      mgrs[z++] = utmrow[(yh + (zone1 & 1 ? utmevenrowshift : 0))
			 % utmrowperiod];
    } else {
      bool eastp = xh >= upseasting;
      int iband = (northp ? 2 : 0) + (eastp ? 1 : 0);
      mgrs[z++] = upsband[iband];
      mgrs[z++] = upscols[iband][xh - (eastp ? upseasting :
				       northp ? minupsNind : minupsSind)];
      mgrs[z++] = upsrows[northp][yh - (northp ? minupsNind : minupsSind)];
    }
    double mult = std::pow(double(base), (std::min)(prec - tilelevel, 0));
    int
      ix = int(floor(xf * mult)),
      iy = int(floor(yf * mult));
    for (int c = (std::min)(prec, int(tilelevel)); c--;) {
      mgrs[z + c] = digits[ ix % base ];
      ix /= base;
      mgrs[z + c + prec] = digits[ iy % base ];
      iy /= base;
    }
    if (prec > tilelevel) {
      xf -= floor(xf * mult);
      yf -= floor(yf * mult);
      mult = std::pow(double(base), prec - tilelevel);
      ix = int(floor(xf * mult));
      iy = int(floor(yf * mult));
      for (int c = prec - tilelevel; c--;) {
	mgrs[z + c + tilelevel] = digits[ ix % base ];
	ix /= base;
	mgrs[z + c + tilelevel + prec] = digits[ iy % base ];
	iy /= base;
      }
    }
      
  }

  void MGRS::Forward(int zone, bool northp, double x, double y,
		     int prec, std::string& mgrs) {
    double lat, lon;
    if (zone)
      UTMUPS::Reverse(zone, northp, x, y, lat, lon);
    else
      // Latitude isn't needed for UPS specs.
      lat = 0;
    Forward(zone, northp, x, y, lat, prec, mgrs);
  }

  void MGRS::Reverse(const std::string& mgrs,
		     int& zone, bool& northp, double& x, double& y,
		     int& prec, bool centerp) {
    int
      p = 0,
      len = int(mgrs.size());
    zone = 0;
    while (p < len) {
      int i = lookup(digits, mgrs[p]);
      if (i < 0)
	break;
      zone = 10 * zone + i;
      ++p;
    }
    if (p > 0 && (zone == 0 || zone > 60))
      throw std::out_of_range("Zone " + str(zone) + " not in [1,60]");
    if (p > 2)
      throw std::out_of_range("More than 2 digits at start of MGRS "
			      + mgrs.substr(0, p));
    if (len - p < 3)
      throw std::out_of_range("MGRS string " + mgrs + " too short");
    bool utmp = zone != 0;
    int zone1 = zone - 1;
    const std::string& band = utmp ? latband : upsband;
    int iband = lookup(band, mgrs[p++]);
    if (iband < 0)
      throw std::out_of_range("Band letter " + str(mgrs[p-1])
			      + " not in " + band);
    northp = iband >= (utmp ? 10 : 2);
    const std::string& col = utmp ? utmcols[zone1 % 3] : upscols[iband];
    const std::string& row = utmp ? utmrow : upsrows[northp];
    int icol = lookup(col, mgrs[p++]);
    if (icol < 0)
      throw std::out_of_range("Column letter " + str(mgrs[p-1])
			      + " not in " + col);
    int irow = lookup(row, mgrs[p++]);
    if (irow < 0)
      throw std::out_of_range("Row letter " + str(mgrs[p-1])
			      + " not in " + row);
    if (utmp) {
      if (zone1 & 1)
	irow = (irow + utmrowperiod - utmevenrowshift) % utmrowperiod;
      iband -= 10;
      irow = UTMRow(iband, icol, irow);
      if (irow == maxutmSrow)
	throw std::out_of_range("Block " + mgrs.substr(p-2, 2)
				+ " not in zone/band " + mgrs.substr(0, p-2));

      irow = northp ? irow : irow + 100;
      icol = icol + minutmcol;
    } else {
      bool eastp = iband & 1;
      icol += eastp ? upseasting : northp ? minupsNind : minupsSind;
      irow += northp ? minupsNind : minupsSind;
    }
    prec = (len - p)/2;
    double unit = tile;
    x = unit * icol;
    y = unit * irow;
    for (int i = 0; i < prec; ++i) {
      unit /= base;
      int
	ix = lookup(digits, mgrs[p + i]),
	iy = lookup(digits, mgrs[p + i + prec]);
      if (ix < 0 || iy < 0)
	throw std::out_of_range("Encountered a non-digit in " + mgrs.substr(p));
      x += unit * ix;
      y += unit * iy;
    }
    if ((len - p) % 2) {
      if (lookup(digits, mgrs[len - 1]) < 0)
	throw std::out_of_range("Encountered a non-digit in " + mgrs.substr(p));
      else
	throw std::out_of_range("Not an even number of digits in "
				+ mgrs.substr(p));
    }
    if (prec > maxprec)
      throw std::out_of_range("More than " + str(2*maxprec) + " digits in "
			      + mgrs.substr(p));
    if (centerp) {
      x += unit/2;
      y += unit/2;
    }
  }

  void MGRS::CheckCoords(bool utmp, bool northp, double x, double y) {
    // Limits are all multiples of 100km and are all closed on the lower end and
    // open on the upper end.
    int
      ix = int(floor(x / tile)),
      iy = int(floor(y / tile)),
      ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    if (! (ix >= mineasting[ind] && ix < maxeasting[ind]) )
      throw std::out_of_range("Easting " + str(int(floor(x/1000)))
			      + "km not in ["
			      + str(mineasting[ind]*tile/1000) + "km, "
			      + str(maxeasting[ind]*tile/1000) + "km)");
    if (! (iy >= minnorthing[ind] && iy < maxnorthing[ind]) )
      throw std::out_of_range("Northing " + str(int(floor(y/1000)))
			      + "km not in ["
			      + str(minnorthing[ind]*tile/1000) + "km, "
			      + str(maxnorthing[ind]*tile/1000) + "km)");
  }

  int MGRS::UTMRow(int iband, int icol, int irow) {
    // Input is MGRS (periodic) row index and output is true row index.  Band
    // index is in [-10, 10) (as returned by LatitudeBand).  Column index
    // origin is easting = 100km.  Returns maxutmSrow if irow and iband are
    // incompatible.  Row index origin is equator.

    // Estimate center row number for latitude band
    // 90 deg = 100 tiles; 1 band = 8 deg = 100*8/90 tiles
    double c = 100 * (8 * iband + 4)/90.0;
    bool northp = iband >= 0;
    int
      minrow = iband > -10 ? int(floor(c - 4.3 - 0.1 * northp)) : -90,
      maxrow = iband < 9 ? int(floor(c + 4.4 - 0.1 * northp)) : 94,
      baserow = (minrow + maxrow) / 2 - utmrowperiod / 2;
    // Add maxutmSrow = 5 * utmrowperiod to ensure operand is positive
    irow = (irow - baserow + maxutmSrow) % utmrowperiod + baserow;
    if (irow < minrow || irow > maxrow) {
      // Northing = 71*100km and 80*100km intersect band boundaries
      // The following deals with these special cases.
      int
	// Fold [-10,-1] -> [9,0]
	sband = iband >= 0 ? iband : - iband  - 1,
	// Fold [-90,-1] -> [89,0]
	srow = irow >= 0 ? irow : -irow - 1,
	// Fold [4,7] -> [3,0]
	scol = icol < 4 ? icol : -icol + 7;
      if ( ! ( (srow == 70 && sband == 8 && scol >= 2) ||
	       (srow == 71 && sband == 7 && scol <= 2) ||
	       (srow == 79 && sband == 9 && scol >= 1) ||
	       (srow == 80 && sband == 8 && scol <= 1) ) )
	irow = maxutmSrow;
    }
    return irow;
  }

} // namespace GeographicLib

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

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = MGRS_HPP;
}

namespace GeographicLib {

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
      throw std::out_of_range("Illegal zone");
    if (!(prec >= 0 || zone <= 12))
      throw std::out_of_range("Illegal precision");
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
      mgrs[z++] = latband[10 + LatitudeBand(lat)];
      mgrs[z++] = utmcols[zone1 % 3][xh - minutmcol];
      mgrs[z++] = utmrow[(yh + (zone1 & 1 ? utmevenrowshift : 0))
	% utmrowperiod];
    } else {
      bool eastp = x >= upseasting * tile;
      int iband = (northp ? 2 : 0) + (eastp ? 1 : 0);
      mgrs[z++] = upsband[iband];
      mgrs[z++] = upscols[iband][xh - (eastp ? upseasting :
				       northp ? minupsNind : minupsSind)];
      mgrs[z++] = upsrows[northp][yh - (northp ? minupsNind : minupsSind)];
    }
    double mult = std::pow(double(base), prec - tilelevel);
    int
      ix = int(floor(xf * mult)),
      iy = int(floor(yf * mult));
    for (int c = prec; c--;) {
      mgrs[z + c] = digits[ ix % base ];
      ix /= base;
      mgrs[z + c + prec] = digits[ iy % base ];
      iy /= base;
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
		     int& prec) {
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
      throw std::out_of_range("Illegal zone");
    if (p > 2)
      throw std::out_of_range("No more than 2 digits at start");
    if (len - p < 3)
      throw std::out_of_range("MGRS string too short");
    bool utmp = zone != 0;
    int zone1 = zone - 1;
    const std::string& band = utmp ? latband : upsband;
    int iband = lookup(band, mgrs[p++]);
    if (iband < 0)
      throw std::out_of_range("Illegal band letter");
    northp = iband >= (utmp ? 10 : 2);
    const std::string& col = utmp ? utmcols[zone1 % 3] : upscols[iband];
    const std::string& row = utmp ? utmrow : upsrows[northp];
    int icol = lookup(col, mgrs[p++]);
    if (icol < 0)
      throw std::out_of_range("Illegal column letter");
    int irow = lookup(row, mgrs[p++]);
    if (irow < 0)
      throw std::out_of_range("Illegal row letter");
    if (utmp) {
      if (zone1 & 1)
	irow = (irow + utmrowperiod - utmevenrowshift) % utmrowperiod;
      iband -= 10;
      // Estimate center row number for latitude band
      // 90 deg = 100 tiles; 1 band = 8 deg = 100*8/90 tiles
      double c = 100 * (8 * iband + 4)/90.0;
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
	  throw std::out_of_range("Bad band");
      }
      irow = northp ? irow : irow + 100;
      icol = icol + minutmcol;
    } else {
      bool eastp = iband & 1;
      icol += eastp ? upseasting : northp ? minupsNind : minupsSind;
      irow += northp ? minupsNind : minupsSind;
    }
    if ((len - p) % 2)
      throw std::out_of_range("Need even number of digits");
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
	throw std::out_of_range("Expected a digit");
      x += unit * ix;
      y += unit * iy;
    }
    x += unit/2;
    y += unit/2;
  }

  void MGRS::CheckCoords(bool utmp, bool northp, double x, double y) {
    // Limits are all multiples of 100km and are all closed on the lower end and
    // open on the upper end.  This allows compatibility with the MGRS system.
    int
      ix = int(floor(x / MGRS::tile)),
      iy = int(floor(y / MGRS::tile)),
      ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    if (! (ix >= mineasting[ind] && ix < maxeasting[ind]) )
      throw std::out_of_range("Easting out of range");
    if (! (iy >= minnorthing[ind] && ix < maxnorthing[ind]) )
      throw std::out_of_range("Northing out of range");
  }

} // namespace GeographicLib

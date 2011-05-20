/**
 * \file MGRS.cpp
 * \brief Implementation for GeographicLib::MGRS class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/MGRS.hpp"

#define GEOGRAPHICLIB_MGRS_CPP "$Id: MGRS.cpp 6785 2010-01-05 22:15:42Z karney $"

RCSID_DECL(GEOGRAPHICLIB_MGRS_CPP)
RCSID_DECL(GEOGRAPHICLIB_MGRS_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real MGRS::eps =
    // 25 = ceil(log_2(2e7)) -- use half circumference here because northing
    // 195e5 is a legal in the "southern" hemisphere.
    pow(real(0.5), numeric_limits<real>::digits - 25);
  const Math::real MGRS::angeps =
    // 7 = ceil(log_2(90))
    pow(real(0.5), numeric_limits<real>::digits - 7);
  const string MGRS::hemispheres = "SN";
  const string MGRS::utmcols[3] =
    { "ABCDEFGH", "JKLMNPQR", "STUVWXYZ" };
  const string MGRS::utmrow  = "ABCDEFGHJKLMNPQRSTUV";
  const string MGRS::upscols[4] =
    { "JKLPQRSTUXYZ", "ABCFGHJKLPQR", "RSTUXYZ", "ABCFGHJ" };
  const string MGRS::upsrows[2] =
    { "ABCDEFGHJKLMNPQRSTUVWXYZ", "ABCDEFGHJKLMNP" };
  const string MGRS::latband = "CDEFGHJKLMNPQRSTUVWX";
  const string MGRS::upsband = "ABYZ";
  const string MGRS::digits  = "0123456789";

  const int MGRS::mineasting[4] =
    { minupsSind, minupsNind, minutmcol, minutmcol };
  const int MGRS::maxeasting[4] =
    { maxupsSind, maxupsNind, maxutmcol, maxutmcol };
  const int MGRS::minnorthing[4] =
    { minupsSind, minupsNind,
      minutmSrow,  minutmSrow - (maxutmSrow - minutmNrow) };
  const int MGRS::maxnorthing[4] =
    { maxupsSind, maxupsNind,
      maxutmNrow + (maxutmSrow - minutmNrow), maxutmNrow };

  void MGRS::Forward(int zone, bool northp, real x, real y, real lat,
                     int prec, std::string& mgrs) {
    bool utmp = zone != 0;
    CheckCoords(utmp, northp, x, y);
    if (!(zone >= 0 || zone <= 60))
      throw GeographicErr("Zone " + str(zone) + " not in [0,60]");
    if (!(prec >= 0 || prec <= maxprec))
      throw GeographicErr("MGRS precision " + str(prec) + " not in [0, "
                          + str(int(maxprec)) + "]");
    // Fixed char array for accumulating string.  Allow space for zone, 3 block
    // letters, easting + northing.  Don't need to allow for terminating null.
    char mgrs1[2 + 3 + 2 * maxprec];
    int
      zone1 = zone - 1,
      z = utmp ? 2 : 0,
      mlen = z + 3 + 2 * prec;
    if (utmp) {
      mgrs1[0] = digits[ zone / base ];
      mgrs1[1] = digits[ zone % base ];
      // This isn't necessary...!  Keep y non-neg
      // if (!northp) y -= maxutmSrow * tile;
    }
    int
      xh = int(floor(x)) / tile,
      yh = int(floor(y)) / tile;
    real
      xf = x - tile * xh,
      yf = y - tile * yh;
    if (utmp) {
      int
        // Correct fuzziness in latitude near equator
        iband = abs(lat) > angeps ? LatitudeBand(lat) : (northp ? 0 : -1),
        icol = xh - minutmcol,
        irow = UTMRow(iband, icol, yh % utmrowperiod);
      if (irow != yh - (northp ? minutmNrow : maxutmSrow))
        throw GeographicErr("Latitude " + str(lat)
                            + " is inconsistent with UTM coordinates");
      mgrs1[z++] = latband[10 + iband];
      mgrs1[z++] = utmcols[zone1 % 3][icol];
      mgrs1[z++] = utmrow[(yh + (zone1 & 1 ? utmevenrowshift : 0))
                         % utmrowperiod];
    } else {
      bool eastp = xh >= upseasting;
      int iband = (northp ? 2 : 0) + (eastp ? 1 : 0);
      mgrs1[z++] = upsband[iband];
      mgrs1[z++] = upscols[iband][xh - (eastp ? upseasting :
                                       northp ? minupsNind : minupsSind)];
      mgrs1[z++] = upsrows[northp][yh - (northp ? minupsNind : minupsSind)];
    }
    real mult = pow(real(base), max(tilelevel - prec, 0));
    int
      ix = int(floor(xf / mult)),
      iy = int(floor(yf / mult));
    for (int c = min(prec, int(tilelevel)); c--;) {
      mgrs1[z + c] = digits[ ix % base ];
      ix /= base;
      mgrs1[z + c + prec] = digits[ iy % base ];
      iy /= base;
    }
    if (prec > tilelevel) {
      xf -= floor(xf / mult);
      yf -= floor(yf / mult);
      mult = pow(real(base), prec - tilelevel);
      ix = int(floor(xf * mult));
      iy = int(floor(yf * mult));
      for (int c = prec - tilelevel; c--;) {
        mgrs1[z + c + tilelevel] = digits[ ix % base ];
        ix /= base;
        mgrs1[z + c + tilelevel + prec] = digits[ iy % base ];
        iy /= base;
      }
    }
    mgrs.resize(mlen);
    copy(mgrs1, mgrs1 + mlen, mgrs.begin());
  }

  void MGRS::Forward(int zone, bool northp, real x, real y,
                     int prec, std::string& mgrs) {
    real lat, lon;
    if (zone)
      UTMUPS::Reverse(zone, northp, x, y, lat, lon);
    else
      // Latitude isn't needed for UPS specs.
      lat = 0;
    Forward(zone, northp, x, y, lat, prec, mgrs);
  }

  void MGRS::Reverse(const std::string& mgrs,
                     int& zone, bool& northp, real& x, real& y,
                     int& prec, bool centerp) {
    int
      p = 0,
      len = int(mgrs.size());
    int zone1 = 0;
    while (p < len) {
      int i = lookup(digits, mgrs[p]);
      if (i < 0)
        break;
      zone1 = 10 * zone1 + i;
      ++p;
    }
    if (p > 0 && (zone1 == 0 || zone1 > 60))
      throw GeographicErr("Zone " + str(zone1) + " not in [1,60]");
    if (p > 2)
      throw GeographicErr("More than 2 digits at start of MGRS "
                          + mgrs.substr(0, p));
    if (len - p < 3)
      throw GeographicErr("MGRS string " + mgrs + " too short");
    bool utmp = zone1 != 0;
    int zonem1 = zone1 - 1;
    const string& band = utmp ? latband : upsband;
    int iband = lookup(band, mgrs[p++]);
    if (iband < 0)
      throw GeographicErr("Band letter " + str(mgrs[p-1]) + " not in "
                          + (utmp ? "UTM" : "UPS") + " set " + band);
    bool northp1 = iband >= (utmp ? 10 : 2);
    const string& col = utmp ? utmcols[zonem1 % 3] : upscols[iband];
    const string& row = utmp ? utmrow : upsrows[northp1];
    int icol = lookup(col, mgrs[p++]);
    if (icol < 0)
      throw GeographicErr("Column letter " + str(mgrs[p-1]) + " not in "
                          + (utmp ? "zone " + mgrs.substr(0, p-2) :
                             "UPS band " + str(mgrs[p-2]))
                          + " set " + col );
    int irow = lookup(row, mgrs[p++]);
    if (irow < 0)
      throw GeographicErr("Row letter " + str(mgrs[p-1]) + " not in "
                          + (utmp ? "UTM" :
                             "UPS " + str(hemispheres[northp1]))
                          + " set " + row);
    if (utmp) {
      if (zonem1 & 1)
        irow = (irow + utmrowperiod - utmevenrowshift) % utmrowperiod;
      iband -= 10;
      irow = UTMRow(iband, icol, irow);
      if (irow == maxutmSrow)
        throw GeographicErr("Block " + mgrs.substr(p-2, 2)
                            + " not in zone/band " + mgrs.substr(0, p-2));

      irow = northp1 ? irow : irow + 100;
      icol = icol + minutmcol;
    } else {
      bool eastp = iband & 1;
      icol += eastp ? upseasting : northp1 ? minupsNind : minupsSind;
      irow += northp1 ? minupsNind : minupsSind;
    }
    int prec1 = (len - p)/2;
    real
      unit = tile,
      x1 = unit * icol,
      y1 = unit * irow;
    for (int i = 0; i < prec1; ++i) {
      unit /= base;
      int
        ix = lookup(digits, mgrs[p + i]),
        iy = lookup(digits, mgrs[p + i + prec1]);
      if (ix < 0 || iy < 0)
        throw GeographicErr("Encountered a non-digit in " + mgrs.substr(p));
      x1 += unit * ix;
      y1 += unit * iy;
    }
    if ((len - p) % 2) {
      if (lookup(digits, mgrs[len - 1]) < 0)
        throw GeographicErr("Encountered a non-digit in " + mgrs.substr(p));
      else
        throw GeographicErr("Not an even number of digits in "
                            + mgrs.substr(p));
    }
    if (prec1 > maxprec)
      throw GeographicErr("More than " + str(2*maxprec) + " digits in "
                          + mgrs.substr(p));
    if (centerp) {
      x1 += unit/2;
      y1 += unit/2;
    }
    zone = zone1;
    northp = northp1;
    x = x1;
    y = y1;
    prec = prec1;
  }

  void MGRS::CheckCoords(bool utmp, bool& northp, real& x, real& y) {
    // Limits are all multiples of 100km and are all closed on the lower end
    // and open on the upper end -- and this is reflected in the error
    // messages.  However if a coordinate lies on the excluded upper end (e.g.,
    // after rounding), it is shifted down by eps.  This also folds UTM
    // northings to the correct N/S hemisphere.
    int
      ix = int(floor(x / tile)),
      iy = int(floor(y / tile)),
      ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
    if (! (ix >= mineasting[ind] && ix < maxeasting[ind]) ) {
      if (ix == maxeasting[ind] && x == maxeasting[ind] * tile)
        x -= eps;
      else
        throw GeographicErr("Easting " + str(int(floor(x/1000)))
                            + "km not in MGRS/"
                            + (utmp ? "UTM" : "UPS") + " range for "
                            + (northp ? "N" : "S" ) + " hemisphere ["
                            + str(mineasting[ind]*tile/1000) + "km, "
                            + str(maxeasting[ind]*tile/1000) + "km)");
    }
    if (! (iy >= minnorthing[ind] && iy < maxnorthing[ind]) ) {
      if (iy == maxnorthing[ind] && y == maxnorthing[ind] * tile)
        y -= eps;
      else
        throw GeographicErr("Northing " + str(int(floor(y/1000)))
                            + "km not in MGRS/"
                            + (utmp ? "UTM" : "UPS") + " range for "
                            + (northp ? "N" : "S" ) + " hemisphere ["
                            + str(minnorthing[ind]*tile/1000) + "km, "
                            + str(maxnorthing[ind]*tile/1000) + "km)");
    }

    // Correct the UTM northing and hemisphere if necessary
    if (utmp) {
      if (northp && iy < minutmNrow) {
        northp = false;
        y += utmNshift;
      } else if (!northp && iy >= maxutmSrow) {
        if (y == maxutmSrow * tile)
          // If on equator retain S hemisphere
          y -= eps;
        else {
          northp = true;
          y -= utmNshift;
        }
      }
    }
  }

  int MGRS::UTMRow(int iband, int icol, int irow) throw() {
    // Input is MGRS (periodic) row index and output is true row index.  Band
    // index is in [-10, 10) (as returned by LatitudeBand).  Column index
    // origin is easting = 100km.  Returns maxutmSrow if irow and iband are
    // incompatible.  Row index origin is equator.

    // Estimate center row number for latitude band
    // 90 deg = 100 tiles; 1 band = 8 deg = 100*8/90 tiles
    real c = 100 * (8 * iband + 4)/real(90);
    bool northp = iband >= 0;
    int
      minrow = iband > -10 ?
      int(floor(c - real(4.3) - real(0.1) * northp)) : -90,
      maxrow = iband <   9 ?
      int(floor(c + real(4.4) - real(0.1) * northp)) :  94,
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

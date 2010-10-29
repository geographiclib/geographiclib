/**
 * \file OSGB.cpp
 * \brief Implementation for GeographicLib::OSGB class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <iostream>
#include "GeographicLib/OSGB.hpp"

#define GEOGRAPHICLIB_OSGB_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_OSGB_CPP)
RCSID_DECL(GEOGRAPHICLIB_OSGB_HPP)

namespace GeographicLib {

  using namespace std;

  const string OSGB::letters = "ABCDEFGHJKLMNOPQRSTUVWXYZ";
  const string OSGB::digits = "0123456789";

  const TransverseMercator
  OSGB::OSGBTM(MajorRadius(), InverseFlattening(), CentralScale());

  Math::real OSGB::computenorthoffset() throw() {
    real x, y;
    OSGBTM.Forward(real(0), OriginLatitude(), real(0), x, y);
    return FalseNorthing() - y;
  }

  const Math::real OSGB::northoffset = computenorthoffset();

  void OSGB::GridReference(real x, real y, int prec, std::string& gridref) {
    CheckCoords(x, y);
    if (!(prec >= 0 || prec <= maxprec))
      throw GeographicErr("OSGB precision " + str(prec) + " not in [0, "
                          + str(int(maxprec)) + "]");
    char grid[2 + 2 * maxprec];
    int
      xh = int(floor(x)) / tile,
      yh = int(floor(y)) / tile;
    real
      xf = x - tile * xh,
      yf = y - tile * yh;
    xh += tileoffx;
    yh += tileoffy;
    int z = 0;
    grid[z++] = letters[(tilegrid - (yh / tilegrid) - 1)
                        * tilegrid + (xh / tilegrid)];
    grid[z++] = letters[(tilegrid - (yh % tilegrid) - 1)
                        * tilegrid + (xh % tilegrid)];
    real mult = pow(real(base), max(tilelevel - prec, 0));
    int
      ix = int(floor(xf / mult)),
      iy = int(floor(yf / mult));
    for (int c = min(prec, int(tilelevel)); c--;) {
      grid[z + c] = digits[ ix % base ];
      ix /= base;
      grid[z + c + prec] = digits[ iy % base ];
      iy /= base;
    }
    if (prec > tilelevel) {
      xf -= floor(xf / mult);
      yf -= floor(yf / mult);
      mult = pow(real(base), prec - tilelevel);
      ix = int(floor(xf * mult));
      iy = int(floor(yf * mult));
      for (int c = prec - tilelevel; c--;) {
        grid[z + c + tilelevel] = digits[ ix % base ];
        ix /= base;
        grid[z + c + tilelevel + prec] = digits[ iy % base ];
        iy /= base;
      }
    }
    int mlen = z + 2 * prec;
    gridref.resize(mlen);
    copy(grid, grid + mlen, gridref.begin());
  }

  void OSGB::GridReference(const std::string& gridref,
                           real& x, real& y, int& prec,
                           bool centerp) {
    int
      len = int(gridref.size()),
      p = 0;
    char grid[2 + 2 * maxprec];
    for (int i = 0; i < len; ++i) {
      if (!isspace(gridref[i])) {
        if (p >= 2 + 2 * maxprec)
          throw GeographicErr("OSGB string " + gridref + " too long");
        grid[p++] = gridref[i];
      }
    }
    len = p;
    p = 0;
    if (len < 2)
      throw GeographicErr("OSGB string " + gridref + " too short");
    if (len % 2)
      throw GeographicErr("OSGB string " + gridref +
                          " has odd number of characters");
    int
      xh = 0,
      yh = 0;
    while (p < 2) {
      int i = lookup(letters, grid[p++]);
      if (i < 0)
        throw GeographicErr("Illegal prefix character " + gridref);
      yh = yh * tilegrid + tilegrid - (i / tilegrid) - 1;
      xh = xh * tilegrid + (i % tilegrid);
    }
    xh -= tileoffx;
    yh -= tileoffy;

    int prec1 = (len - p)/2;
    real
      unit = tile,
      x1 = unit * xh,
      y1 = unit * yh;
    for (int i = 0; i < prec1; ++i) {
      unit /= base;
      int
        ix = lookup(digits, grid[p + i]),
        iy = lookup(digits, grid[p + i + prec1]);
      if (ix < 0 || iy < 0)
        throw GeographicErr("Encountered a non-digit in " + gridref);
      x1 += unit * ix;
      y1 += unit * iy;
    }
    if (centerp) {
      x1 += unit/2;
      y1 += unit/2;
    }
    x = x1;
    y = y1;
    prec = prec1;
  }

  void OSGB::CheckCoords(real x, real y) {
    // Limits are all multiples of 100km and are all closed on the lower end
    // and open on the upper end -- and this is reflected in the error
    // messages.
    if (! (x >= minx && x < maxx) )
      throw GeographicErr("Easting " + str(int(floor(x/1000)))
                          + "km not in OSGB range ["
                          + str(minx/1000) + "km, "
                          + str(maxx/1000) + "km)");
    if (! (y >= miny && y < maxy) )
      throw GeographicErr("Northing " + str(int(floor(y/1000)))
                          + "km not in OSGB range ["
                          + str(miny/1000) + "km, "
                          + str(maxy/1000) + "km)");
  }
} // namespace GeographicLib

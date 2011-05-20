/**
 * \file LocalCartesian.cpp
 * \brief Implementation for GeographicLib::LocalCartesian class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/LocalCartesian.hpp"

#define GEOGRAPHICLIB_LOCALCARTESIAN_CPP "$Id: LocalCartesian.cpp 6785 2010-01-05 22:15:42Z karney $"

RCSID_DECL(GEOGRAPHICLIB_LOCALCARTESIAN_CPP)
RCSID_DECL(GEOGRAPHICLIB_LOCALCARTESIAN_HPP)

namespace GeographicLib {

  using namespace std;

  void LocalCartesian::Reset(real lat0, real lon0, real h0) throw() {
    _lat0 = lat0;
    _lon0 = lon0 >= 180 ? lon0 - 360 : lon0 < -180 ? lon0 + 360 : lon0;
    _h0 = h0;
    _earth.Forward(_lat0, _lon0, _h0, _x0, _y0, _z0);
    real
      phi = lat0 * Constants::degree(),
      sphi = sin(phi),
      cphi = abs(_lat0) == 90 ? 0 : cos(phi),
      lam = lon0 * Constants::degree(),
      slam = _lon0 == -180 ? 0 : sin(lam),
      clam = abs(_lon0) == 90 ? 0 : cos(lam);
    // Local x axis in geocentric coords
    _rxx = -slam; _rxy = clam; _rxz = 0;
    // Local y axis in geocentric coords
    _ryx = -clam * sphi; _ryy = -slam * sphi; _ryz = cphi;
    // Local z axis in geocentric coords
    _rzx = clam * cphi; _rzy = slam * cphi; _rzz = sphi;
  }

  void LocalCartesian::Forward(real lat, real lon, real h,
                               real& x, real& y, real& z) const throw() {
    real xc, yc, zc;
    _earth.Forward(lat, lon, h, xc, yc, zc);
    xc -= _x0; yc -= _y0; zc -= _z0;
    x = _rxx * xc + _rxy * yc + _rxz * zc;
    y = _ryx * xc + _ryy * yc + _ryz * zc;
    z = _rzx * xc + _rzy * yc + _rzz * zc;
  }

  void LocalCartesian::Reverse(real x, real y, real z,
                               real& lat, real& lon, real& h) const throw() {
    real
      xc = _x0 + _rxx * x + _ryx * y + _rzx * z,
      yc = _y0 + _rxy * x + _ryy * y + _rzy * z,
      zc = _z0 + _rxz * x + _ryz * y + _rzz * z;
    _earth.Reverse(xc, yc, zc, lat, lon, h);
  }

} // namespace GeographicLib

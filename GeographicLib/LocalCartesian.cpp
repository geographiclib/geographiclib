/**
 * \file LocalCartesian.cpp
 * \brief Implementation for GeographicLib::LocalCartesian class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/LocalCartesian.hpp"
#include "GeographicLib/Constants.hpp"
#include <cmath>
#include <stdexcept>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = LOCALCARTESIAN_HPP;
}

namespace GeographicLib {

  using namespace std;

  void LocalCartesian::Reset(double lat0, double lon0, double h0) throw() {
    _lat0 = lat0;
    _lon0 = lon0;
    _h0 = h0;
    ECEF::WGS84.Forward(lat0, lon0, h0, _x0, _y0, _z0);
    double
      phi = lat0 * Constants::degree,
      sphi = sin(phi),
      cphi = cos(phi),
      lam = lon0 * Constants::degree,
      slam = sin(lam),
      clam = cos(lam);
    // Local x axis in ECEF coords
    _rxx = -slam; _rxy = clam; _rxz = 0;
    // Local y axis in ECEF coords
    _ryx = -clam * sphi; _ryy = -slam * sphi; _ryz = cphi;
    // Local z axis in ECEF coords
    _rzx = clam * cphi; _rzy = slam * cphi; _rzz = sphi;
  }

  void LocalCartesian::Forward(double lat, double lon, double h,
			       double& x, double& y, double& z) const throw() {
    double xc, yc, zc;
    ECEF::WGS84.Forward(lat, lon, h, xc, yc, zc);
    xc -= _x0; yc -= _y0; zc -= _z0;
    x = _rxx * xc + _rxy * yc + _rxz * zc;
    y = _ryx * xc + _ryy * yc + _ryz * zc;
    z = _rzx * xc + _rzy * yc + _rzz * zc;
  }

  void LocalCartesian::Reverse(double x, double y, double z,
			       double& lat, double& lon, double& h)
    const throw() {
    double
      xc = _x0 + _rxx * x + _ryx * y + _rzx * z,
      yc = _y0 + _rxy * x + _ryy * y + _rzy * z,
      zc = _z0 + _rxz * x + _ryz * y + _rzz * z;
    ECEF::WGS84.Reverse(xc, yc, zc, lat, lon, h);
  }

} // namespace GeographicLib

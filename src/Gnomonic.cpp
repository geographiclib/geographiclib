/**
 * \file Gnomonic.cpp
 * \brief Implementation for GeographicLib::Gnomonic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/Gnomonic.hpp"

#define GEOGRAPHICLIB_GNOMONIC_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GNOMONIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_GNOMONIC_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real Gnomonic::eps0 = numeric_limits<real>::epsilon();
  const Math::real Gnomonic::eps = real(0.01) * sqrt(eps0);

  void Gnomonic::Forward(real lat0, real lon0, real lat, real lon,
                         real& x, real& y, real& azi, real& rk)
    const throw() {
    real sig, s, azi0, m;
    sig = _earth.Inverse(lat0, lon0, lat, lon, s, azi0, azi, m);
    const GeodesicLine line(_earth.Line(lat0, lon0, azi0));
    real M, Mx;
    line.Scale(sig, M, Mx);
    rk = M;
    if (M <= 0)
      x = y = Math::NaN();
    else {
      real rho = m/M;
      azi0 *= Constants::degree();
      x = rho * sin(azi0);
      y = rho * cos(azi0);
    }
  }

  void Gnomonic::Reverse(real lat0, real lon0, real x, real y,
                         real& lat, real& lon, real& azi, real& rk)
    const throw() {
    real
      azi0 = atan2(x, y) / Constants::degree(),
      rho = min(Math::hypot(x, y), _a/(2 * eps0));
    GeodesicLine line(_earth.Line(lat0, lon0, azi0));
    real lat1, lon1, azi1, M, s;
    int count = numit;
    if (rho * _f < _a / 2)
      s = _a * std::atan(rho/_a);
    else {
      real m, Mx, ang = 90;
      int trip = _f == 0 ? 1 : -std::log(rho/_a) / std::log(_f);
      while (count--) {
        s = line.Position(ang, lat1, lon1, azi1, m, true);
        line.Scale(ang, M, Mx);
        if (trip < 0 && M > 0)
          break;
        // Estimate new arc length assuming dM/da = -1.
        ang += (M - m/rho)/Constants::degree();
        if (M > 0)
          --trip;
      }
      s -= (m/M - rho) * M * M;
    }
    int trip = 0;
    // Reset count if previous iteration was OK; otherwise skip next iteration
    count = count < 0 ? 0 : numit;
    while (count--) {
      real m, Mx;
      real ang = line.Position(s, lat1, lon1, azi1, m);
      line.Scale(ang, M, Mx);
      if (trip)
        break;
      else if (M <= 0) {
        count = -1;
        break;
      }
      real ds = (m/M - rho) * M * M;
      s -= ds;
      if (std::abs(ds) < eps * _a)
        ++trip;
    }
    if (count >= 0) {
      lat = lat1; lon = lon1; azi = azi1; rk = M;
    } else
      lat = lon = azi = rk = Math::NaN();
    return;
  }

} // namespace GeographicLib


/**
 * \file CassiniSoldner.cpp
 * \brief Implementation for GeographicLib::CassiniSoldner class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/CassiniSoldner.hpp"

#define GEOGRAPHICLIB_CASSINISOLDNER_CPP "$Id: CassiniSoldner.cpp 6863 2010-09-08 21:28:58Z karney $"

RCSID_DECL(GEOGRAPHICLIB_CASSINISOLDNER_CPP)
RCSID_DECL(GEOGRAPHICLIB_CASSINISOLDNER_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real CassiniSoldner::eps1 =
    real(0.01) * sqrt(numeric_limits<real>::epsilon());
  const Math::real CassiniSoldner::eps2 = sqrt(numeric_limits<real>::min());

  void CassiniSoldner::Reset(real lat0, real lon0) throw() {
    _meridian = _earth.Line(lat0, lon0, real(0),
                            Geodesic::LATITUDE | Geodesic::LONGITUDE |
                            Geodesic::DISTANCE | Geodesic::DISTANCE_IN |
                            Geodesic::AZIMUTH);
    real
      phi = LatitudeOrigin() * Constants::degree(),
      f = _earth.InverseFlattening() != 0 ? 1 / _earth.InverseFlattening() : 0;
    _sbet0 = (1 - f) * sin(phi);
    _cbet0 = abs(LatitudeOrigin()) == 90 ? 0 : cos(phi);
    SinCosNorm(_sbet0, _cbet0);
  }

  void CassiniSoldner::Forward(real lat, real lon, real& x, real& y,
                               real& azi, real& rk) const throw() {
    if (!Init())
      return;
    real dlon = AngNormalize(lon - LongitudeOrigin());
    real sig12, s12, azi1, azi2;
    lat = AngRound(lat);
    sig12 = _earth.Inverse(lat, -abs(dlon), lat, abs(dlon), s12, azi1, azi2);
    if (sig12 < 100 * eps2)
      sig12 = s12 = 0;
    sig12 *= real(0.5);
    s12 *= real(0.5);
    if (s12 == 0) {
      real da = (azi2 - azi1)/2;
      if (abs(dlon) <= 90) {
        azi1 = 90 - da;
        azi2 = 90 + da;
      } else {
        azi1 = -90 - da;
        azi2 = -90 + da;
      }
    }
    if (dlon < 0) {
      azi2 = azi1;
      s12 = -s12;
      sig12 = -sig12;
    }
    x = s12;
    azi = AngNormalize(azi2);
    GeodesicLine perp(_earth.Line(lat, dlon, azi2, Geodesic::GEODESICSCALE));
    real t;
    perp.GenPosition(true, -sig12,
                     Geodesic::GEODESICSCALE,
                     t, t, t, t, t, t, rk, t);

    real
      alp0 = perp.EquatorialAzimuth() * Constants::degree(),
      calp0 = cos(alp0), salp0 = sin(alp0),
      sbet1 = lat >=0 ? calp0 : -calp0,
      cbet1 = abs(dlon) <= 90 ? abs(salp0) : -abs(salp0),
      sbet01 = sbet1 * _cbet0 - cbet1 * _sbet0,
      cbet01 = cbet1 * _cbet0 + sbet1 * _sbet0,
      sig01 = atan2(sbet01, cbet01) / Constants::degree();
    _meridian.GenPosition(true, sig01,
                          Geodesic::DISTANCE,
                          t, t, t, y, t, t, t, t);
  }

  void CassiniSoldner::Reverse(real x, real y, real& lat, real& lon,
                               real& azi, real& rk) const throw() {
    if (!Init())
      return;
    real lat1, lon1;
    real azi0, t;
    _meridian.Position(y, lat1, lon1, azi0);
    _earth.Direct(lat1, lon1, azi0 + 90, x, lat, lon, azi, rk, t);
  }

} // namespace GeographicLib

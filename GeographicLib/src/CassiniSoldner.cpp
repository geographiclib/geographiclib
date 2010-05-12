
/**
 * \file CassiniSoldner.cpp
 * \brief Implementation for GeographicLib::CassiniSoldner class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/CassiniSoldner.hpp"

#define GEOGRAPHICLIB_CASSINISOLDNER_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_CASSINISOLDNER_CPP)
RCSID_DECL(GEOGRAPHICLIB_CASSINISOLDNER_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real CassiniSoldner::eps1 =
    real(0.01) * sqrt(numeric_limits<real>::epsilon());
  const Math::real CassiniSoldner::eps2 = sqrt(numeric_limits<real>::min());

  void CassiniSoldner::Reset(real lat0, real lon0) throw() {
    _meridian = _earth.Line(lat0, lon0, real(0));
    real phi = LatitudeOrigin() * Constants::degree();
    _sbet0 = _earth._f1 * sin(phi);
    _cbet0 = abs(LatitudeOrigin()) == 90 ? 0 : cos(phi);
    Geodesic::SinCosNorm(_sbet0, _cbet0);
  }

  void CassiniSoldner::Forward(real lat, real lon, real& x, real& y,
                               real& azi, real& m) const throw() {
    if (!Init())
      return;
    real dlon = Geodesic::AngNormalize(lon - LongitudeOrigin());
    real sig12, s12, azi1, azi2, m12;
    lat = Geodesic::AngRound(lat);
    sig12 = _earth.Inverse(lat, -abs(dlon), lat, abs(dlon),
                           s12, azi1, azi2, m12);
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
    azi = Geodesic::AngNormalize(azi2);
    GeodesicLine perp = _earth.Line(lat, dlon, azi2);
    m = Scale(perp, sig12);

    real
      sbet1 = lat >=0 ? perp._calp0 : -perp._calp0,
      cbet1 = abs(dlon) <= 90 ? abs(perp._salp0) : -abs(perp._salp0),
      sbet01 = sbet1 * _cbet0 - cbet1 * _sbet0,
      cbet01 = cbet1 * _cbet0 + sbet1 * _sbet0,
      sig01 = atan2(sbet01, cbet01) / Constants::degree();
    real latx, lonx, azix, m12x;
    y = _meridian.Position(sig01, latx, lonx, azix, m12x, true);
  }

  void CassiniSoldner::Reverse(real x, real y, real& lat, real& lon,
                               real& azi, real& m) const throw() {
    if (!Init())
      return;
    real lat1, lon1;
    real azi0, m0;
    _meridian.Position(y, lat1, lon1, azi0, m0);
    GeodesicLine perp = _earth.Line(lat1, lon1, azi0 + 90);
    real sig12 = perp.Position(x, lat, lon, azi, m0);
    m = Scale(perp, sig12);
  }

  Math::real CassiniSoldner::Scale(const GeodesicLine& perp, real sig12)
    const throw() {
    // Return
    // M12 = csig1 * csig2 
    // + ( ssig1 * ssig2 * sqrt(1+k2*sq(ssig2))
    //   - ssig1 * csig2 * (J(sig2) - J(sig1)) ) / sqrt(1+k2*sq(ssig1))
    // Simplify by setting ssig1 = 1, csig1 = 0, B11 = B21 = 0
    sig12 *= Constants::degree();
    real
      ssig2 = cos(sig12), csig2 = -sin(sig12), // sig2 = pi/2 + sig12
      et = (1 + perp._A1m1)
      * Geodesic::SinSeries(ssig2, csig2, perp._C1, Geodesic::nC1),
      ez = (1 + perp._A2m1)
      * Geodesic::SinSeries(ssig2, csig2, perp._C2, Geodesic::nC2);
    return ( ssig2 * sqrt(1 + perp._k2 * sq(ssig2)) -
             csig2 * ( (perp._A1m1 - perp._A2m1) * sig12 + (et - ez) ) )
      / sqrt(1 + perp._k2);
  }

} // namespace GeographicLib

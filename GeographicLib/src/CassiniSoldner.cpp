
/**
 * \file CassiniSoldner.cpp
 * \brief Implementation for GeographicLib::CassiniSoldner class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://www.petrel.org/geographic/
 **********************************************************************/

#include "GeographicLib/CassiniSoldner.hpp"
#include <limits>

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
    if (sig12 == 0)
      return 1;
    // Result is symmetric in sig12, however numerical accuracy is better if
    // sigc-sig2 is smaller.
    sig12 = abs(sig12) * Constants::degree();
    // Point on meridian, sigma = pi/2
    real
      // ssig1 = 1, csig1 = 0,
      dtau1 = Geodesic::SinSeries(real(1), real(0), perp._tauCoeff,
                                  Geodesic::ntau),
      dzet1 = Geodesic::SinSeries(real(1), real(0), perp._zetCoeff,
                                  Geodesic::nzet);
    // Find semi-conjugate point -- initial guess, sigma = pi
    real ssigc = 0, csigc = -1;
    real wc, dtauc, dzetc;
    for (unsigned trip = 0, numit = 0; numit < maxit; ++numit) {
      wc = sqrt(1 + perp._u2 * sq(ssigc));
      dtauc = Geodesic::SinSeries(ssigc, csigc, perp._tauCoeff, Geodesic::ntau);
      dzetc = Geodesic::SinSeries(ssigc, csigc, perp._zetCoeff, Geodesic::nzet);
      real
        sig1c = atan2(-csigc, ssigc),
        et = (1 + perp._taufm1) * ( dtauc - dtau1 ),
        ez = (1 + perp._zetfm1) * ( dzetc - dzet1 ),
        j1c = ( (perp._taufm1 - perp._zetfm1) * sig1c + (et - ez) ),
        v = - 2 * wc * csigc * ssigc + 2 * sq(csigc) * j1c;
      if (abs(v) <= eps2 || trip > 1)
        break;
      real dv = - 2 * wc * (1 - 2 * sq(ssigc)) - 4 * csigc * ssigc * j1c;
      real dsig = -v/dv;
      real
        sdsigc = sin(dsig),
        cdsigc = cos(dsig),
        nssigc = ssigc * cdsigc + csigc * sdsigc;
      csigc = csigc * cdsigc - ssigc * sdsigc;
      ssigc = nssigc;
      if (abs(v) < eps1)
        ++trip;
    }
    // Scale meridian to conjugate
    real m1c = - sqrt(1 + perp._u2 ) * csigc;

    // Target point, sigma = sig12 + pi/2
    real
      ssig2 = cos(sig12), csig2 = -sin(sig12),
      dtau2 = Geodesic::SinSeries(ssig2, csig2, perp._tauCoeff, Geodesic::ntau),
      dzet2 = Geodesic::SinSeries(ssig2, csig2, perp._zetCoeff, Geodesic::nzet),
      et = (1 + perp._taufm1) * ( dtauc - dtau2 ),
      ez = (1 + perp._zetfm1) * ( dzetc - dzet2 ),
      sig2c = atan2(ssigc * csig2 - csigc * ssig2,
                    csigc * csig2 + ssigc * ssig2),
      j2c = ( (perp._taufm1 - perp._zetfm1) * sig2c + (et - ez) ),
      // Scale target to conjugate
      m2c = ((wc * csig2 * ssigc -
              sqrt(1 + perp._u2 * sq(ssig2)) * ssig2 * csigc)
             - csig2 * csigc * j2c);
    return m2c/m1c;
  }

} // namespace GeographicLib

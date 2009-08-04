/**
 * \file CassiniSoldner.cpp
 * \brief Implementation for GeographicLib::CassiniSoldner class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/CassiniSoldner.hpp"
#include "GeographicLib/Constants.hpp"
#include <limits>
#include <cmath>
#include <stdexcept>
#include <iostream>

#define CASSINISOLDNER_CPP "$Id$"

RCSID_DECL(CASSINISOLDNER_CPP)
RCSID_DECL(CASSINISOLDNER_HPP)

namespace GeographicLib {

  using namespace std;

  void CassiniSoldner::Reset(double lat0, double lon0) throw() {
    _meridian = _earth.Line(lat0, lon0, 0.0);
  }

  void CassiniSoldner::Forward(double lat, double lon,
			       double& x, double& y,
			       double& azi, double& m) const throw() {
    const double eps1 = 0.01 * sqrt(numeric_limits<double>::epsilon());
    const double eps2 = sqrt(numeric_limits<double>::min());
    const unsigned maxit = 10;
    double
      dlon = Geodesic::AngNormalize(lon) - LongitudeOrigin();
    double sig12, s12, azi1, azi2, m12;
    lat = Geodesic::AngRound(lat);
    sig12 = _earth.Inverse(lat, -abs(dlon), lat, abs(dlon),
			   s12, azi1, azi2, m12);
    if (dlon < 0) {
      azi2 = azi1;
      s12 = -s12;
    }
    if (sig12 == 0) {
      m12 = 0;
      s12 = 0;
      azi2 = Geodesic::AngNormalize(90 + (lat >= 0 ? 1 : -1) * dlon);
    }
    double
      phi = lat * Constants::degree(),
      sbet = _earth._f1 * sin(phi),
      cbet = abs(lat) == 90 ? 0 : cos(phi),
      alp2 = azi2 * Constants::degree(),
      salp2 = abs(azi2) == 180 ? 0 : sin(alp2),
      calp2 = abs(azi2) == 90 ? 0 : cos(alp2);
    Geodesic::SinCosNorm(sbet, cbet);
    double
      salp0 = salp2 * cbet,
      calp0 = Geodesic::hypot(calp2, salp2 * sbet),
      sphi0 = calp0,
      cphi0 = _earth._f1 * abs(salp0),
      lat0 = atan2(sphi0, cphi0) / Constants::degree();
    GeodesicLine perp = _earth.Line(lat0, 0.0, 90.0);
    // Find semi-conjugate point -- initial guess
    double
      ssigc = 0,
      csigc = -1,
      dtau0 = Geodesic::SinSeries(1.0, 0.0, perp._tauCoeff, Geodesic::ntau),
      dzet0 = Geodesic::SinSeries(1.0, 0.0, perp._zetCoeff, Geodesic::nzet);
    double v;
    int i = 0;
    for (unsigned trip = 0, numit = 0; numit < maxit; ++numit) {
      double
	sig12 = atan2(-csigc, ssigc),
	et = (1 + perp._taufm1) *
	( Geodesic::SinSeries(ssigc, csigc, perp._tauCoeff, Geodesic::ntau)
	  - dtau0 ),
	ez = (1 + perp._zetfm1) *
	( Geodesic::SinSeries(ssigc, csigc, perp._zetCoeff, Geodesic::nzet)
	  - dzet0 ),
	wc = sqrt(1 + perp._u2 * sq(ssigc)),
	j12 = ( (perp._taufm1 - perp._zetfm1) * sig12 + (et - ez) );
      v = - 2 * wc * csigc * ssigc + 2 * sq(csigc) * j12;
      std::cout << i << " " << ssigc << " " << csigc
		<< " " << v << "\n";
      if (abs(v) <= eps2 || trip > 1)
	break;
      double dv = - 2 * wc * (1 - 2 * sq(ssigc)) - 4 * csigc * ssigc * j12;
      double dsig = -v/dv;
      double
	sdsigc = sin(dsig),
	cdsigc = cos(dsig),
	nssigc = ssigc * cdsigc + csigc * sdsigc;
      csigc = csigc * cdsigc - ssigc * sdsigc;
      ssigc = nssigc;
      if (abs(v) < eps1)
	++trip;
    }
  }

  void CassiniSoldner::Reverse(double x, double y,
			       double& lat, double& lon,
			       double& azi, double& m) const throw() {
    double lat1, lon1;
    _meridian.Position(y, lat1, lon1, azi, m);
    _earth.Direct(lat1, lon1, 90.0, x, lat, lon, azi, m);
  }

} // namespace GeographicLib

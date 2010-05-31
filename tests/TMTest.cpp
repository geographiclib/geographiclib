/**
 * \file TMTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with ProjExact.o
 * EllipticFunction.o Proj.o
 *
 * See \ref transversemercatortest for usage information.
 **********************************************************************/

#include "GeographicLib/TransverseMercator.hpp"
#include "GeographicLib/TransverseMercatorExact.hpp"
#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/Constants.hpp"
#include <vector>
#include <algorithm>

#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

GeographicLib::Math::real
dist(GeographicLib::Math::real a, GeographicLib::Math::real r,
     long double lat0, long double lon0,
     GeographicLib::Math::real lat1, GeographicLib::Math::real lon1) {
  using namespace GeographicLib;
  typedef Math::real real;
  real
    phi = real(lat0) * Constants::degree(),
    f = r != 0 ? 1/r : 0,
    e2 = f * (2 - f),
    sinphi = sin(phi),
    n = 1/sqrt(1 - e2 * sinphi * sinphi),
      // See Wikipedia article on latitude
    hlon = std::cos(phi) * n,
    hlat = (1 - e2) * n * n * n;
  long double dlon = (long double)(lon1) - lon0;
  if (dlon >= 180) dlon -= 360;
  else if (dlon < -180) dlon += 360;
  return a * Constants::degree() *
    Math::hypot(real((long double)(lat1) - lat0) * hlat, real(dlon) * hlon);
}

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"TMTest [-s] [-d]\n\
$Id$\n\
\n\
Read in TMcoords.dat on standard input and test TransverseMercatorExact\n\
or (if -s is given) TransverseMercator.  If -d dump the error for wach\n\
point; otherwise summarise errors.\n";
  return retval;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;
  bool series = false;
  bool dump = false;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-s")
      series = true;
    else if (arg == "-d")
      dump = true;
    else
      return usage(arg != "-h");
  }

  try {
    real minlat = series ? 0 : (dump ? -100 : -15);
    const unsigned nbins = 101;
    std::vector<real> d(nbins);
    std::vector<real> errv(nbins, 0);
    std::vector<real> errvg(nbins, 0);
    std::vector<real> errvk(nbins, 0);
    real esterr = sizeof(real) == sizeof(double) ? (series ? 3e-9 : 8e-9) :
      (series ? 4e-12 : 4e-12);
    for (unsigned i = 0; i < nbins; ++i)
      d[i] = 100e3 * i;
    d[0] = 10e3;
    d[nbins - 1] = 10001966;
    const TransverseMercator& tm = TransverseMercator::UTM;
    const TransverseMercatorExact tme(Constants::WGS84_a(),
                                      Constants::WGS84_r(),
                                      Constants::UTM_k0(),
                                      true);
    real
      a = series ? tm.MajorRadius() : tme.MajorRadius(),
      r = series ? tm.InverseFlattening() : tme.InverseFlattening();
    const Geodesic geod(a, r);
    long double lat0l, lon0l, x0l, y0l, gam0l, k0l;
    while (std::cin >> lat0l >> lon0l >> x0l >> y0l >> gam0l >> k0l) {
      real
        lat0 = lat0l,
        lon0 = lon0l,
        x0 = x0l,
        y0 = y0l,
        gam0 = gam0l,
        k0 = k0l;
      if (lat0 < minlat)
        continue;
      real azi1, azi2, s12, m12;
      real errf, errgf, errkf, errr, errgr, errkr;
      geod.Inverse(std::max(lat0,real(0)), lon0, std::max(lat0,real(0)), -lon0,
                   s12, azi1, azi2, m12);
      s12 /= 2;
      real lat, lon, x, y, gam, k;
      if (series)
        tm.Forward(0, lat0, lon0, x, y, gam, k);
      else
        tme.Forward(0, lat0, lon0, x, y, gam, k);
      errf = real(Math::hypot((long double)(x) - x0l,
                              (long double)(y) - y0l)) / k0;
      errgf = real(std::abs((long double)(gam) - gam0));
      errkf = real(std::abs((long double)(k) - k0));
      if (series)
        tm.Reverse(0, x0, y0, lat, lon, gam, k);
      else
        tme.Reverse(0, x0, y0, lat, lon, gam, k);
      errr = dist(a, r, lat0l, lon0l, lat, lon);
      errgr = real(std::abs((long double)(gam) - gam0));
      errkr = real(std::abs((long double)(k) - k0));

      real
        err = std::max(errf, errr),
        errg = std::max(errgf, errgr)
        - esterr/(a * std::sin((90 - lat0) * Constants::degree())
                  * Constants::degree()),
        errk = std::max(errkf, errkr) / k0;
      if (dump)
        std::cout << std::fixed << std::setprecision(12)
                  << lat0 << " " << lon0 << " "
                  << std::scientific << std::setprecision(4)
                  << errf << " " << errr << " "
                  << errgf << " " << errgr << " "
                  << errkf << " " << errkr << "\n";
      else
        for (unsigned i = 0; i < nbins; ++i) {
          if (s12 <= d[i]) {
            errv[i] = std::max(err, errv[i]);
            errvg[i] = std::max(errg, errvg[i]);
            errvk[i] = std::max(errk, errvk[i]);
          }
        }
    }
    if (!dump)
      for (unsigned i = 0; i < nbins; ++i)
        std::cout << int(d[i]/1000) << " "
                  << errv[i] << " "
                  << errvg[i] << " "
                  << errvk[i] << "\n";
  }
  catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

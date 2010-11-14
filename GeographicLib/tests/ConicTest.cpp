/**
 * \file ProjTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/LambertConformalConic.hpp"
#include "GeographicLib/Constants.hpp"
#include "GeographicLib/DMS.hpp"
#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

GeographicLib::Math::real
dist(GeographicLib::Math::real a, GeographicLib::Math::real r,
     GeographicLib::Math::real lat0, GeographicLib::Math::real lon0,
     GeographicLib::Math::real lat1, GeographicLib::Math::real lon1) {
  using namespace GeographicLib;
  typedef Math::real real;
  real
    phi = lat0 * Constants::degree(),
    f = r != 0 ? 1/r : 0,
    e2 = f * (2 - f),
    sinphi = sin(phi),
    n = 1/sqrt(1 - e2 * sinphi * sinphi),
      // See Wikipedia article on latitude
    hlon = std::cos(phi) * n,
    hlat = (1 - e2) * n * n * n,
    dlon = lon1 - lon0;
  if (dlon >= 180) dlon -= 360;
  else if (dlon < -180) dlon += 360;
  return a * Constants::degree() *
    Math::hypot((lat1 - lat0) * hlat, dlon * hlon);
}

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"ConicTest -l -s\n\
$Id$\n\
\n\
Checks conic projections\n";
  return retval;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  using namespace std;
  typedef Math::real real;
  if (false) {
    for (int i = -12; i <= 0; ++i) {
      Math::extended
        colat = std::pow(Math::extended(10), i),
        lat = 90 - colat,
        coslat1 = cos(lat*Math::edegree()),
        coslat2 = sin(colat*Math::edegree());
      std::cout << i << " " << coslat1 << " " << coslat2 << " " <<  (coslat1-coslat2)/coslat2 << "\n";
      std::cout << i << " " << asin(coslat1)/Math::edegree() - colat << "\n";
    }
    return 0;
  }
  bool lambert = true;
  bool albert = false;
  bool checkstdlats = false;
  real a = Constants::WGS84_a(), r = Constants::WGS84_r();
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-l") {
      lambert = true;
      albert = false;
    } else if (arg == "-s") {
      checkstdlats = true;
    } else if (arg == "-e") {
      if (m + 2 >= argc) return usage(1);
      try {
        a = DMS::Decode(std::string(argv[m + 1]));
        r = DMS::Decode(std::string(argv[m + 2]));
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
        return 1;
      }
      m += 2;
    } else
      return usage(arg != "-h");
  }

  try {
    std::cout << std::setprecision(14);
    Math::extended quant = 1e12L;
    while (true) {
      Math::extended lat1, lat2, lat0, k0;
      if (!(cin >> lat1 >> lat2 >> lat0 >> k0))
        break;
      int
        sign1 = lat1 < 0 ? -1 : 1,
        sign2 = lat2 < 0 ? -1 : 1;
      lat1 = floor(sign1 * lat1 * quant + 0.5L);
      lat2 = floor(sign2 * lat2 * quant + 0.5L);
      Math::extended
        colat1 = (90 * quant - lat1) / quant,
        colat2 = (90 * quant - lat2) / quant;
      lat1 /= quant;
      lat2 /= quant;
      Math::extended
        sinlat1 = sign1 * (lat1 < 45 ? sin(lat1 * Math::edegree())
                           : cos(colat1 * Math::edegree())),
        sinlat2 = sign2 * (lat2 < 45 ? sin(lat2 * Math::edegree())
                           : cos(colat2 * Math::edegree())),
        coslat1 = (lat1 < 45 ? cos(lat1 * Math::edegree())
                   : sin(colat1 * Math::edegree())),
        coslat2 = (lat2 < 45 ? cos(lat2 * Math::edegree())
                   : sin(colat2 * Math::edegree()));
      lat1 *= sign1;
      lat2 *= sign2;
      const LambertConformalConic lam(a, r, /* real(lat1), real(lat2), */
                                      real(sinlat1), real(coslat1),
                                      real(sinlat2), real(coslat2),
                                      real(1));
      Math::extended lat0a = lam.OriginLatitude(), k0a = lam.CentralScale();
      if (abs(lat0a-lat0) > 0.4e-13L || abs(k0a - k0) > 7e-15L * k0 )
        cout << lat1 << " " << lat2 << " " << lat0 << " " << lat0a << " " << lat0a - lat0 << " " << (k0a - k0)/k0 << "\n";
      /*
      const LambertConformalConic lamb(a, r, real(sin(lat1*degree)), real(cos(lat1*degree)), real(sin(lat2*degree)), real(cos(lat2*degree)), real(1));
      lat0a = lamb.OriginLatitude(); k0a = lamb.CentralScale();
      cout << lat1 << " " << lat2 << " " << lat0 << " " << lat0a << " " << lat0a - lat0 << " " << k0a - k0 << "\n";
      */
    }
  }
  catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

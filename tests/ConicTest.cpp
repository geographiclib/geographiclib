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

GeographicLib::Math::extended
dist(GeographicLib::Math::extended a, GeographicLib::Math::extended r,
     GeographicLib::Math::extended lat0, GeographicLib::Math::extended lon0,
     GeographicLib::Math::extended lat1, GeographicLib::Math::extended lon1) {
  using namespace GeographicLib;
  typedef Math::extended extended;
  extended
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

class TestData {
  // Read test data with one line of buffering
public:
  typedef GeographicLib::Math::extended extended;
private:
  std::istream& _is;
  bool _usesaved;               // Should GetNext use saved values?
  bool _valid;                  // Are there saved values?
  extended _lat0, _lat, _lon, _x, _y, _k;
public:
  TestData(std::istream& is) : _is(is), _usesaved(false), _valid(false) {}
  bool GetNext(extended& lat0, extended& lat, extended& lon,
               extended& x, extended& y, extended& k) {
    if (_usesaved)
      _usesaved = false;
    else {
      _valid = (_is >> _lat0 >> _lat >> _lon >> _x >> _y >> _k);
      if (!_valid)
        return false;
    }
    lat0 = _lat0; lat = _lat; lon = _lon; x = _x; y = _y; k = _k;
    return true;
  }
  bool BackUp() {
    if (!_valid || _usesaved)
      return false;             // Can't backup up again
    else
      return _usesaved = true;  // Set flag for GetNext
  }
};

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
  typedef Math::extended extended;
  if (false) {
    for (int i = -12; i <= 0; ++i) {
      extended
        colat = std::pow(extended(10), i),
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
    if (checkstdlats) {         // stdin contains lat1 lat2 lat0 k0
      std::cout << std::setprecision(14);
      extended quant = 1e12L;
      while (true) {
        extended lat1, lat2, lat0, k0;
        if (!(cin >> lat1 >> lat2 >> lat0 >> k0))
          break;
        int
          sign1 = lat1 < 0 ? -1 : 1,
          sign2 = lat2 < 0 ? -1 : 1;
        lat1 = floor(sign1 * lat1 * quant + 0.5L);
        lat2 = floor(sign2 * lat2 * quant + 0.5L);
        extended
          colat1 = (90 * quant - lat1) / quant,
          colat2 = (90 * quant - lat2) / quant;
        lat1 /= quant;
        lat2 /= quant;
        extended
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
        extended lat0a = lam.OriginLatitude(), k0a = lam.CentralScale();
        if (abs(lat0a-lat0) > 0.4e-13L || abs(k0a - k0) > 7e-15L * k0 )
          cout << lat1 << " " << lat2 << " " << lat0 << " " << lat0a << " " << lat0a - lat0 << " " << (k0a - k0)/k0 << "\n";
      }
    } else { // Check projection
      // stdin contains lat0 lat lon x y k
      TestData testset(std::cin);
      cout << setprecision(8);
      while (true) {
        extended lat0, lat, lon, x, y, k;
        if (!testset.GetNext(lat0, lat, lon, x, y, k))
          break;
        if (!testset.BackUp())
          break;
        int
          sign0 = lat0 < 0 ? -1 : 1;
        extended quant = 1e12L;
        extended
          lat00 = floor(sign0 * lat0 * quant + 0.5L),
          colat00 = (90 * quant - lat00) / quant;
        lat00 /= quant;
        real
          sinlat0 = real(sign0 * (lat00 < 45 ? sin(lat00 * Math::edegree())
                                  : cos(colat00 * Math::edegree()))),
          coslat0 = real(lat00 < 45 ? cos(lat00 * Math::edegree())
                         : sin(colat00 * Math::edegree()));
        const LambertConformalConic lcc(a, r,
                                        sinlat0, coslat0, sinlat0, coslat0,
                                        real(1));
        real maxerrf = 0, maxerrr = 0, maxerrk = 0;
        real latf = 0, lonf = 0, latr = 0, lonr = 0, latk = 0, lonk = 0;
        // std::cout << "New lat0: " << lat0 << "\n";
        while (true) {
          extended lat0x;
          if (!testset.GetNext(lat0x, lat, lon, x, y, k))
            break;
          if (lat0 != lat0x) {
            testset.BackUp();
            break;
          }
          // std::cout << "Process: " << lat0 << " " << lat << " " << lon << " " << k << "\n";
          real lata, lona, xa, ya, gammaa, gammab, ka, kb;
          lcc.Forward(real(0), real(lat), real(lon), xa, ya, gammaa, ka);
          real errf = Math::hypot(extended(xa) - x, extended(ya) - y) / k;
          real errk = abs(extended(ka) - k);
          lcc.Reverse(real(0), real(x), real(y), lata, lona, gammab, kb);
          real errr = dist(extended(a), extended(r),
                           lat, lon, extended(lata), extended(lona));
          std::cout << lata << " " << lona << " " << xa << " " << ya << " "
                    << ka << " " << kb << " "
                    << gammaa << " " << gammab << "\n";
          errk = max(errk, real(abs(extended(kb) - k)));
          errk /= k;
          if (errf > maxerrf) {
            maxerrf = errf;
            latf = lat;
            lonf = lon;
          }
          if (errr > maxerrr) {
            maxerrr = errr;
            latr = lat;
            lonr = lon;
          }
          if (errk > maxerrk && abs(lat) < 89) {
            maxerrk = errk;
            latk = lat;
            lonk = lon;
          }
          std::cout << lat0 << " " << lat << " " << lon << " "
                    << errf << " " << errr << " " << errk << "\n";
        }
        std::cout << "Max errf: " << lat0 << " "
                  << maxerrf << " " << latf << " " << lonf << "\n";
        std::cout << "Max errr: " << lat0 << " "
                  << maxerrr << " " << latr << " " << lonr << "\n";
        std::cout << "Max errk: " << lat0 << " "
                  << maxerrk << " " << latk << " " << lonk << "\n";
      }
    }
  }
  catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

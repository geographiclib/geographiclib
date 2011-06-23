/**
 * \file ProjTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/LambertConformalConic.hpp"
#include "GeographicLib/AlbersEqualArea.hpp"
#include "GeographicLib/Constants.hpp"
#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/DMS.hpp"
#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

GeographicLib::Math::extended
dist(GeographicLib::Math::extended a, GeographicLib::Math::extended f,
     GeographicLib::Math::extended lat0, GeographicLib::Math::extended lon0,
     GeographicLib::Math::extended lat1, GeographicLib::Math::extended lon1) {
  using namespace GeographicLib;
  typedef Math::extended extended;
  extended
    phi = lat0 * Math::degree<extended>(),
    e2 = f * (2 - f),
    sinphi = sin(phi),
    n = 1/sqrt(1 - e2 * sinphi * sinphi),
      // See Wikipedia article on latitude
    hlon = std::cos(phi) * n,
    hlat = (1 - e2) * n * n * n,
    dlon = lon1 - lon0;
  if (dlon >= 180) dlon -= 360;
  else if (dlon < -180) dlon += 360;
  return a * Math::degree<extended>() *
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
  bool lambert = true;
  bool albers = false;
  bool checkstdlats = false;
  real a = Constants::WGS84_a(), f = Constants::WGS84_f();
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-l") {
      lambert = true;
      albers = false;
    } else if (arg == "-a") {
      lambert = false;
      albers = true;
    } else if (arg == "-s") {
      checkstdlats = true;
    } else if (arg == "-e") {
      if (m + 2 >= argc) return usage(1);
      try {
        a = DMS::Decode(std::string(argv[m + 1]));
        f = DMS::DecodeFraction(std::string(argv[m + 2]));
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
        return 1;
      }
      m += 2;
      if (f > 1) f = 1/f;
    } else
      return usage(arg != "-h");
  }

  try {
    if (checkstdlats) {         // stdin contains lat1 lat2 lat0 k0
      std::cout << std::setprecision(17);
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
          sinlat1 = sign1 * (lat1 < 45 ? sin(lat1 * Math::degree<extended>())
                             : cos(colat1 * Math::degree<extended>())),
          sinlat2 = sign2 * (lat2 < 45 ? sin(lat2 * Math::degree<extended>())
                             : cos(colat2 * Math::degree<extended>())),
          coslat1 = (lat1 < 45 ? cos(lat1 * Math::degree<extended>())
                     : sin(colat1 * Math::degree<extended>())),
          coslat2 = (lat2 < 45 ? cos(lat2 * Math::degree<extended>())
                     : sin(colat2 * Math::degree<extended>()));
        lat1 *= sign1;
        lat2 *= sign2;
        const LambertConformalConic lam(a, f, /* real(lat1), real(lat2), */
                                        real(sinlat1), real(coslat1),
                                        real(sinlat2), real(coslat2),
                                        real(1));
        const AlbersEqualArea alb(a, f, /* real(lat1), real(lat2), */
                                        real(sinlat1), real(coslat1),
                                        real(sinlat2), real(coslat2),
                                        real(1));
        extended
          lat0a = albers ? alb.OriginLatitude() : lam.OriginLatitude();
          //k0a = albers ? alb.CentralScale() : lam.CentralScale();
        if (!(abs(lat0a-lat0) <= 4.5e-14))
          cout << lat1 << " " << lat2 << " " << lat0
               << " " << lat0a << " " << lat0a - lat0 << "\n";
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
          sinlat0 = real(sign0 * (lat00 < 45 ?
                                  sin(lat00 * Math::degree<extended>()) :
                                  cos(colat00 * Math::degree<extended>()))),
          coslat0 = real(lat00 < 45 ? cos(lat00 * Math::degree<extended>())
                         : sin(colat00 * Math::degree<extended>()));
        const LambertConformalConic lcc(a, f,
                                        sinlat0, coslat0, sinlat0, coslat0,
                                        real(1));
        const AlbersEqualArea alb(a, f,
                                  sinlat0, coslat0, sinlat0, coslat0,
                                  real(1));
        real maxerrf = 0, maxerrr = 0, maxerrkf = 0, maxerrkr = 0;
        real latf = 0, lonf = 0, latr = 0, lonr = 0,
          latkf = 0, lonkf = 0, latkr = 0, lonkr = 0;
        // std::cout << "New lat0: " << lat0 << "\n";
        while (true) {
          extended lat0x;
          if (!testset.GetNext(lat0x, lat, lon, x, y, k))
            break;
          if (lat0 != lat0x) {
            testset.BackUp();
            break;
          }
          real latb, lonb, xa, ya, gammaa, gammab, ka, kb;
          if (albers)
            alb.Forward(real(0), real(lat), real(lon), xa, ya, gammaa, ka);
          else
            lcc.Forward(real(0), real(lat), real(lon), xa, ya, gammaa, ka);
          real errf = Math::hypot(extended(xa) - x, extended(ya) - y);
          if (lambert)
            errf /= k;
          real errkf = abs(extended(ka) - k)/k;
          if (albers)
            alb.Reverse(real(0), real(x), real(y), latb, lonb, gammab, kb);
          else
            lcc.Reverse(real(0), real(x), real(y), latb, lonb, gammab, kb);
          real errr = dist(extended(a), extended(f),
                           lat, lon, extended(latb), extended(lonb));
          /*
          std::cout << latb << " " << lonb << " " << xa << " " << ya << " "
                    << ka << " " << kb << " "
                    << gammaa << " " << gammab << "\n";
          */
          real errkr = real(abs(extended(kb) - k))/k;
          if (!(errf <= maxerrf)) {
            maxerrf = errf;
            latf = lat;
            lonf = lon;
          }
          if (!(errr <= maxerrr)) {
            maxerrr = errr;
            latr = lat;
            lonr = lon;
          }
          if (!(errkf <= maxerrkf || abs(lat) >= 89)) {
            maxerrkf = errkf;
            latkf = lat;
            lonkf = lon;
          }
          if (!(errkr <= maxerrkr || abs(lat) >= 89)) {
            maxerrkr = errkr;
            latkr = lat;
            lonkr = lon;
          }
          std::cout << lat0 << " " << lat << " " << lon << " "
                    << errf << " " << errr << " "
                    << errkf << " " << errkr << "\n";
        }
        std::cout << "Max errf: " << lat0 << " "
                  << maxerrf << " " << latf << " " << lonf << "\n";
        std::cout << "Max errr: " << lat0 << " "
                  << maxerrr << " " << latr << " " << lonr << "\n";
        std::cout << "Max errkf: " << lat0 << " "
                  << maxerrkf << " " << latkf << " " << lonkf << "\n";
        std::cout << "Max errkr: " << lat0 << " "
                  << maxerrkr << " " << latkr << " " << lonkr << "\n";
      }
    }
  }
  catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

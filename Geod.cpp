/**
 * \file Geod.cpp
 * \brief Command line utility for geodesic calculations
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o Geod Geod.cpp Geodesic.cpp Constants.cpp
 *
 * See \ref geod for usage information.
 **********************************************************************/

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "GeographicLib/Geodesic.hpp"

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: Geod [-r] [-l lat0 lon0 h0] [-h]\n\
$Id$\n\
\n\
TODO\n\
-h prints this help\n";
  return retval;
}

int main(int argc, char* argv[]) {
  bool linecalc = false, reverse = false;
  double latlonbearing[3] = {0, 0, 0};
  for (int m = 1; m < argc; ++m) {
    std::string arg = std::string(argv[m]);
    if (arg == "-r") {
      reverse = true;
      linecalc = false;
    } else if (arg == "-l") {
      reverse = false;
      linecalc = true;
      for (unsigned i = 0; i < 3; ++i) {
	if (++m == argc) return usage(1);
	std::string a = std::string(argv[m]);
	std::istringstream str(a);
	if (!(str >> latlonbearing[i])) return usage(1);
      }
    } else
      return usage(arg != "-h");
  }

  const GeographicLib::Geodesic& geod = GeographicLib::Geodesic::WGS84;
  GeographicLib::GeodesicLine l;

  std::cout << std::setprecision(16);

  if (linecalc) {
    l = geod.Line(latlonbearing[0], latlonbearing[1], latlonbearing[2]);
    while (true) {
      double s12, lat2, lon2, bearing2;
      std::cin >> s12;
      if (!std::cin.good())
	break;
      l.Position(s12, lat2, lon2, bearing2);
      std::cout << lat2 << " " << lon2 << " " << bearing2 << "\n";
    }
  } else if (reverse) {
    while (true) {
      double lat1, lon1, bearing1, lat2, lon2, bearing2, s12;
      std::cin >> lat1 >> lon1 >> lat2 >> lon2;
      if (!std::cin.good())
	break;
      geod.Inverse(lat1, lon1, lat2, lon2, s12, bearing1, bearing2);
      std::cout << bearing1 << " " << bearing2 << " " << s12 << "\n";
    }
  } else {
    while (true) {
      double lat1, lon1, bearing1, lat2, lon2, bearing2, s12;
      std::cin >> lat1 >> lon1 >> bearing1 >> s12;
      if (!std::cin.good())
      break;
      geod.Direct(lat1, lon1, bearing1, s12, lat2, lon2, bearing2);
      std::cout << lat2 << " " << lon2 << " " << bearing2 << "\n";
    }
  }
  return 0;
}

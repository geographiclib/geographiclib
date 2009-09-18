/**
 * \file EquidistantTest.cpp
 * \brief Command line utility for azimuthal equidistant and Cassini-Soldner
 * projections
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o EquidistantTest EquidistantTest.cpp Geodesic.cpp AzimuthalEquidistant.cpp CassiniSoldner.cpp
 *
 * See \ref equidistanttest for usage information.
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/AzimuthalEquidistant.hpp"
#include "GeographicLib/CassiniSoldner.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: EquidistantTest [-c lat0 lon0] [-z lat0 lon0] [-r] [-h]\n\
$Id$\n\
\n\
Convert geodetic coordinates to either azimuthal equidistant or\n\
Cassini-Soldner coordinates.  The center of the projection (lat0, lon0)\n\
is specified by either the -c option (for Cassini-Soldner) or -z option\n\
(for azimuthal equidistant).  At least one of these options must be\n\
given (the last one given is used).  The WGS84 model of the earth is\n\
used.\n\
\n\
Geodetic coordinates are provided on standard input as a set of lines\n\
containing (blank separated) latitude and longitude (decimal degrees).\n\
For each set of geodetic coordinates, the corresponding projected\n\
coordinates x, y (meters) are printed on standard output together with\n\
the azimuth azi (degrees) and reciprocal scale rk.  For Cassini-Soldner,\n\
azi is the bearing of the easting direction and the scale in the\n\
northing direction is 1/rk.  For azimuthal equidistant, azi is the\n\
bearing of the radial direction and the scale in the azimuthal direction\n\
is 1/rk.\n\
\n\
If -r is given the reverse transformation is performed.  x and y are\n\
given on standard input and each line of the standard output gives\n\
latitude, longitude, azi, and rk.\n\
\n\
-h prints this help.\n";
  return retval;
}

int main(int argc, char* argv[]) {
  bool azimuthal = false, cassini = false, reverse = false;
  double lat0 = 0, lon0 = 0;
  for (int m = 1; m < argc; ++m) {
    std::string arg = std::string(argv[m]);
    if (arg == "-r")
      reverse = true;
    else if (arg == "-c" || arg == "-z") {
      cassini = arg == "-c";
      azimuthal = arg != "-c";
      for (unsigned i = 0; i < 2; ++i) {
	if (++m == argc) return usage(1);
	std::string a = std::string(argv[m]);
	std::istringstream str(a);
	if (!(str >> (i ? lon0 : lat0))) return usage(1);
      }
    } else
      return usage(arg != "-h");
  }

  if (!(azimuthal || cassini)) {
    std::cerr << "Must specify \"-c lat lon\" or \"-z lat lon\"\n";
    return 1;
  }

  if ( !(-90 <= lat0 && lat0 <= 90) ) {
    std::cerr << "Latitude not in range [-90, 90]\n";
    return 1;
  } else if ( !(-180 <= lon0 && lon0 <= 360) ) {
    std::cerr << "Longitude not in range [-180, 360]\n";
    return 1;
  }

  const GeographicLib::CassiniSoldner cs = cassini ?
    GeographicLib::CassiniSoldner(lat0, lon0, GeographicLib::Geodesic::WGS84) :
    GeographicLib::CassiniSoldner(GeographicLib::Geodesic::WGS84);
  const GeographicLib::AzimuthalEquidistant az(GeographicLib::Geodesic::WGS84);

  std::string s;
  int retval = 0;
  std::cout << std::setprecision(16);
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      double lat, lon, x, y, a, m;
      if (!(reverse ?
	    (str >> x >> y) :
	    (str >> lat >> lon)))
	throw std::out_of_range("Incomplete input: " + s);
      if (reverse) {
	if (cassini)
	  cs.Reverse(x, y, lat, lon, a, m);
	else
	  az.Reverse(lat0, lon0, x, y, lat, lon, a, m);
	std::cout << lat << " " << lon << " " << a << " " << m << "\n";
      } else {
	if ( !(-90 <= lat && lat <= 90) )
	  throw std::out_of_range("Latitude not in range [-90, 90]\n");
	if ( !(-180 <= lon && lat <= 360) )
	  throw std::out_of_range("Longitude not in range [-180, 360]");
	if (cassini)
	  cs.Forward(lat, lon, x, y, a, m);
	else
	  az.Forward(lat0, lon0, lat, lon, x, y, a, m);
	std::cout << x << " " << y << " " << a << " " << m << "\n";
      }
    }
    catch (std::out_of_range& e) {
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }

  return retval;
}

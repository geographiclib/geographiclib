/**
 * \file EquidistantTest.cpp
 * \brief Command line utility for azimuthal equidistant and Cassini-Soldner
 * projections
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geodesic.o AzimuthalEquidistant.o
 * CassiniSoldner.o
 *
 * See \ref equidistanttest for usage information.
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/AzimuthalEquidistant.hpp"
#include "GeographicLib/CassiniSoldner.hpp"
#include "GeographicLib/DMS.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>

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
containing (blank separated) latitude and longitude (degrees or DMS).\n\
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
  using namespace GeographicLib;
  typedef Math::real real;
  bool azimuthal = false, cassini = false, reverse = false;
  real lat0 = 0, lon0 = 0;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-r")
      reverse = true;
    else if (arg == "-c" || arg == "-z") {
      cassini = arg == "-c";
      azimuthal = arg != "-c";
      if (m + 2 >= argc) return usage(1);
      try {
        DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                          lat0, lon0);
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding arguments of " << arg << ": "
                  << e.what() << "\n";
        return 1;
      }
      m += 2;
    } else
      return usage(arg != "-h");
  }

  if (!(azimuthal || cassini)) {
    std::cerr << "Must specify \"-c lat lon\" or \"-z lat lon\"\n";
    return 1;
  }

  const CassiniSoldner cs = cassini ?
    CassiniSoldner(lat0, lon0, Geodesic::WGS84) :
    CassiniSoldner(Geodesic::WGS84);
  const AzimuthalEquidistant az(Geodesic::WGS84);

  std::string s;
  int retval = 0;
  std::cout << std::setprecision(16);
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      real lat, lon, x, y, a, m;
      std::string stra, strb;
      if (!(str >> stra >> strb))
        throw GeographicErr("Incomplete input: " + s);
      if (reverse) {
        x = DMS::Decode(stra);
        y = DMS::Decode(strb);
      } else
        DMS::DecodeLatLon(stra, strb, lat, lon);
      std::string strc;
      if (str >> strc)
        throw GeographicErr("Extraneous input: " + strc);
      if (reverse) {
        if (cassini)
          cs.Reverse(x, y, lat, lon, a, m);
        else
          az.Reverse(lat0, lon0, x, y, lat, lon, a, m);
        std::cout << lat << " " << lon << " " << a << " " << m << "\n";
      } else {
        if ( !(-90 <= lat && lat <= 90) )
          throw GeographicErr("Latitude not in range [-90, 90]");
        if ( !(-180 <= lon && lat <= 360) )
          throw GeographicErr("Longitude not in range [-180, 360]");
        if (cassini)
          cs.Forward(lat, lon, x, y, a, m);
        else
          az.Forward(lat0, lon0, lat, lon, x, y, a, m);
        std::cout << x << " " << y << " " << a << " " << m << "\n";
      }
    }
    catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }

  return retval;
}

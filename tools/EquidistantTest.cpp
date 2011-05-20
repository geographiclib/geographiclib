/**
 * \file EquidistantTest.cpp
 * \brief Command line utility for azimuthal equidistant and Cassini-Soldner
 * projections
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geodesic.o GeodesicLine.o
 * AzimuthalEquidistant.o Gnomonic.o CassiniSoldner.o
 *
 * See \ref equidistanttest for usage information.
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/AzimuthalEquidistant.hpp"
#include "GeographicLib/CassiniSoldner.hpp"
#include "GeographicLib/Gnomonic.hpp"
#include "GeographicLib/DMS.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: EquidistantTest -(c|z|g) lat0 lon0 [-r] [-e a r] [-h]\n\
$Id: EquidistantTest.cpp 6859 2010-09-06 14:45:33Z karney $\n\
\n\
Perform projections based on geodesics.  Convert geodetic coordinates to\n\
either azimuthal equidistant, Cassini-Soldner, or gnomonic coordinates.\n\
The center of the projection (lat0, lon0) is specified by either the -c\n\
option (for Cassini-Soldner), the -z option (for azimuthal equidistant),\n\
or the -g option (for gnomonic).  At least one of these options must be\n\
given (the last one given is used).\n\
\n\
By default, the WGS84 ellipsoid is used.  Specifying \"-e a r\" sets the\n\
equatorial radius of the ellipsoid to \"a\" and the reciprocal flattening\n\
to r.  Setting r = 0 results in a sphere.  Specify r < 0 for a prolate\n\
ellipsoid.\n\
\n\
Geodetic coordinates are provided on standard input as a set of lines\n\
containing (blank separated) latitude and longitude (degrees or DMS).\n\
For each set of geodetic coordinates, the corresponding projected\n\
coordinates x, y (meters) are printed on standard output together with\n\
the azimuth azi (degrees) and reciprocal scale rk.  For Cassini-Soldner,\n\
azi is the bearing of the easting direction and the scale in the easting\n\
direction is 1 and the scale in the northing direction is 1/rk.  For\n\
azimuthal equidistant and gnomonic, azi is the bearing of the radial\n\
direction and the scale in the azimuthal direction is 1/rk.  For\n\
azimuthal equidistant and gnomonic, the scales in the radial direction\n\
are 1 and 1/rk^2, respectively.\n\
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
  bool azimuthal = false, cassini = false, gnomonic = false, reverse = false;
  real lat0 = 0, lon0 = 0;
  real
    a = Constants::WGS84_a(),
    r = Constants::WGS84_r();
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-r")
      reverse = true;
    else if (arg == "-c" || arg == "-z" || arg == "-g") {
      cassini = arg == "-c";
      azimuthal = arg == "-z";
      gnomonic = arg == "-g";
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

  if (!(azimuthal || cassini || gnomonic)) {
    std::cerr
      << "Must specify \"-c lat lon\" or \"-z lat lon\" or \"-g lat lon\"\n";
    return 1;
  }

  const Geodesic geod(a, r);
  const CassiniSoldner cs = cassini ?
    CassiniSoldner(lat0, lon0, geod) : CassiniSoldner(geod);
  const AzimuthalEquidistant az(geod);
  const Gnomonic gn(geod);

  std::string s;
  int retval = 0;
  std::cout << std::fixed;
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      real lat, lon, x, y, azi, rk;
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
          cs.Reverse(x, y, lat, lon, azi, rk);
        else if (azimuthal)
          az.Reverse(lat0, lon0, x, y, lat, lon, azi, rk);
        else
          gn.Reverse(lat0, lon0, x, y, lat, lon, azi, rk);
        std::cout << std::setprecision(15)
                  << lat << " " << lon << " " << azi << " "
                  << std::setprecision(16)
                  << rk << "\n";
      } else {
        if (cassini)
          cs.Forward(lat, lon, x, y, azi, rk);
        else if (azimuthal)
          az.Forward(lat0, lon0, lat, lon, x, y, azi, rk);
        else
          gn.Forward(lat0, lon0, lat, lon, x, y, azi, rk);
        std::cout << std::setprecision(10)
                  << x << " " << y << " "
                  << std::setprecision(15)
                  << azi << " "
                  << std::setprecision(16)
                  << rk << "\n";
      }
    }
    catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }

  return retval;
}

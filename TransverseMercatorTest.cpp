/**
 * \file TransverseMercatorTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o TransverseMercatorTest TransverseMercatorTest.cpp TransverseMercatorExact.cpp EllipticFunction.cpp TransverseMercator.cpp
 *
 * See \ref transversemercatortest for usage information.
 **********************************************************************/

#include "GeographicLib/EllipticFunction.hpp"
#include "GeographicLib/TransverseMercatorExact.hpp"
#include "GeographicLib/TransverseMercator.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"TransverseMercatorTest [-r] [-t|-s]\n\
$Id$\n\
\n\
Convert between geographic coordinates and transverse Mercator coordinates.\n\
\n\
Read lines with latitude and longitude (or easting and northing if -r is\n\
specified) from standard input and print latitude, longitude, easting,\n\
northing, convergence, and scale.  Units are degrees and meters.\n\
\n\
By default, the WGS84 is ellipsoid is used, central meridian = 0, UTM\n\
central scale (0.9996), and false easting and false northing are zero.\n\
\n\
If -r is given, the reverse projection is performed (the inputs are easting\n\
and northing).\n\
\n\
If -s is given, the sixth-order Krueger series approximation to the\n\
transverse Mercator projection is used instead of the exact projection.\n\
\n\
If -t is specified, an ellipsoid of eccentricity 0.1 is used, central scale\n\
= 1, 1/4 meridian distance = 1.  In addition, the cut in the exact\n\
transverse Mercator projection at northing = 0 is removed.  The domain of\n\
latitude (lat) and longitude (lon) is the union of\n\
    lat in [0, 90]  and lon in [0, 90]\n\
    lat in (-90, 0] and lon in [81, 90]\n\
The domain of easting (x) and northing (x) is the union of\n\
    x in [0, inf)       and y in [0, 1]\n\
    x in [1.71..., inf) and y in (-inf, 0]\n\
\n\
-s and -t are mutually exclusive (the last flag specified is the operative\n\
one).\n\
\n\
-h prints this help.\n";
  return retval;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real_t real_t;
  bool reverse = false, testing = false, series = false;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-r")
      reverse = true;
    else if (arg == "-t") {
      testing = true;
      series = false;
    } else if (arg == "-s") {
      testing = false;
      series = true;
    } else
      return usage(arg != "-h");
  }

  real_t e, a;
  if (testing) {
    e = real_t(0.1L);
    GeographicLib::EllipticFunction temp(e * e);
    a = 1/temp.E();
  }
  const GeographicLib::TransverseMercatorExact& TME = testing ?
    GeographicLib::TransverseMercatorExact
    (a, (std::sqrt(1 - e * e) + 1) / (e * e), real_t(1), true) :
    GeographicLib::TransverseMercatorExact::UTM;

  const GeographicLib::TransverseMercator& TMS =
    GeographicLib::TransverseMercator::UTM;

  std::string s;
  int retval = 0;
  std::cout << std::setprecision(16);
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      real_t lat, lon, x, y;
      if (!(reverse ?
            (str >> x >> y) :
            (str >> lat >> lon)))
        throw  std::out_of_range("Incomplete input: " + s);
      real_t gamma, k;
      if (reverse) {
        if (series)
          TMS.Reverse(real_t(0), x, y, lat, lon, gamma, k);
        else
          TME.Reverse(real_t(0), x, y, lat, lon, gamma, k);
      } else {
        if (series)
          TMS.Forward(real_t(0), lat, lon, x, y, gamma, k);
        else
          TME.Forward(real_t(0), lat, lon, x, y, gamma, k);
      }
      std::cout << lat << " " << lon << " "
                << x << " " << y << " "
                << gamma << " " << k << "\n";
    }
    catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }

  return retval;
}

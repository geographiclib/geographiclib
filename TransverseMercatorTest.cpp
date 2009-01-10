/**
 * \file TransverseMercatorTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o TransverseMercatorTest TransverseMercatorTest.cpp TransverseMercatorExact.cpp Constants.cpp EllipticFunction.cpp TransverseMercator.cpp
 *
 * Here is the usage (obtained from "TranverseMercatorTest -h")
\verbatim
TransverseMercatorTest [-r] [-t|-s]

Convert between geographic coordinates and transverse Mercator coordinates.

Read lines with latitude and longitude (or easting and northing if -r is
specified) from standard input and print latitude, longitude, easting,
northing, convergence, and scale.  Units are degrees and meters.

By default, the WGS84 is ellipsoid is used, central meridian = 0, UTM
central scale, and false easting and false northing are zero.

If -r is given, the reverse projection is performed (the inputs are easting
and northing).

If -s is given, the sixth-order Krueger series approximation to the
transverse Mercator projection is used instead of the exact projection.

If -t is specified, an ellipsoid of eccentricity 0.1 is used, central scale
= 1, 1/4 meridian distance = 1.  In addition, the cut in the exact
transverse Mercator projection at northing = 0 is removed.  The domain of
latitude (lat) and longitude (lon) is the union of
    lat in [0, 90]  and lon in [0, 90]
    lat in (-90, 0] and lon in [81, 90]
The domain of easting (x) and northing (x) is the union of
    x in [0, inf)       and y in [0, 1]
    x in [1.71..., inf) and y in (-inf, 0]

-s and -t are mutually exclusive (the last flag specified is the operative
one).

-h prints this help
\endverbatim
 **********************************************************************/

#include <string>
#include <iostream>
#include <iomanip>
#include "GeographicLib/EllipticFunction.hpp"
#include "GeographicLib/TransverseMercatorExact.hpp"
#include "GeographicLib/TransverseMercator.hpp"

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
central scale, and false easting and false northing are zero.\n\
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
-h prints this help\n";
  return retval;
}

int main(int argc, char* argv[]) {
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

  double e, a;
  if (testing) {
    e = 0.1;
    GeographicLib::EllipticFunction temp(e * e);
    a = 1/temp.E();
  }
  const GeographicLib::TransverseMercatorExact& TME = testing ?
    GeographicLib::TransverseMercatorExact
    (a, (sqrt(1 - e * e) + 1) / (e * e), 1.0, true) :
    GeographicLib::TransverseMercatorExact::UTM;

  const GeographicLib::TransverseMercator& TMS =
    GeographicLib::TransverseMercator::UTM;

  std::cout << std::setprecision(16);
  while (true) {
    double lat, lon, x, y;
    if (reverse)
      std::cin >> x >> y;
    else
      std::cin >> lat >> lon;
    if (!std::cin.good())
      break;
    double gamma, k;
    if (reverse) {
      if (series)
	TMS.Reverse(0.0, x, y, lat, lon, gamma, k);
      else
	TME.Reverse(0.0, x, y, lat, lon, gamma, k);
    } else {
      if (series)
	TMS.Forward(0.0, lat, lon, x, y, gamma, k);
      else
	TME.Forward(0.0, lat, lon, x, y, gamma, k);
    }
    std::cout << lat << " " << lon << " "
	      << x << " " << y << " "
	      << gamma << " " << k << "\n";
  }
  return 0;
}

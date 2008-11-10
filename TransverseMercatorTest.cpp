/**
 * \file TransverseMercatorExact.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o TransverseMercatorTest TransverseMercatorTest.cpp TransverseMercatorExact.cpp Constants.cpp EllipticFunction.cpp
 *
 **********************************************************************/

#include <string>
#include <iostream>
#include <iomanip>
#include "GeographicLib/EllipticFunction.hpp"
#include "GeographicLib/TransverseMercatorExact.hpp"

namespace {
  char RCSID[] = "$Id$";
}

int main(int argc, char* argv[]) {
  bool reverse = false, testing = false;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-r")
      reverse = true;
    else if (arg == "-t")
      testing = true;
    else {
    std::cerr << "TransverseMercatorTest [-r] [-t]\n\
\n\
Read lines with latitutde and longitude (or easting and northing if -r\n\
is specified) from stdin and print latitude, longitude, easting,\n\
northing, convergence and scale.  By default, WGS84 ellipsoid assumed, central\n\
meridian = 0, UTM central scale, and false easting and false northing\n\
are zero.  With [-t] ellipsoid of eccentricity 0.1 is using, central scale = 1,\n\
1/4 meridian distance = 1.\n";
    return 1;
    }
  }

  double e, a;
  if (testing) {
    e = 0.1;
    GeographicLib::EllipticFunction temp(e * e);
    a = 1/temp.E();
  }
  const GeographicLib::TransverseMercatorExact& TM = testing ?
    GeographicLib::TransverseMercatorExact
    (a, (sqrt(1 - e * e) + 1) / (e * e), 1.0, false) :
    GeographicLib::TransverseMercatorExact::UTM;

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
    if (reverse)
      TM.Reverse(0.0, x, y, lat, lon, gamma, k);
    else
      TM.Forward(0.0, lat, lon, x, y, gamma, k);
    std::cout << lat << " " << lon << " "
	      << x << " " << y << " "
	      << gamma << " " << k << "\n";
  }
  return 0;
}

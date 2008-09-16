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
#include "GeographicLib/TransverseMercatorExact.hpp"

namespace {
  char RCSID[] = "$Id$";
}

int main(int argc, char* argv[]) {
  bool reverse = false;
  if (argc == 2 && std::string(argv[1]) == "-r")
    reverse = true;
  else if (argc != 1) {
    std::cerr << "TransverseMercatorTest [-r]\n\
\n\
Read lines with latitutde and longitude (or easting and northing if -r\n\
is specified) from stdin and print latitude, longitude, easting,\n\
northing, convergence and scale.  WGS84 ellipsoid assumed, central\n\
meridian = 0, UTM central scale, and false easting and false northing\n\
are zero.\n";
    return 1;
  }
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
      GeographicLib::TransverseMercatorExact::UTM.Reverse(0.0, x, y,
							  lat, lon, gamma, k);
    else
      GeographicLib::TransverseMercatorExact::UTM.Forward(0.0, lat, lon,
							  x, y, gamma, k);
    std::cout << lat << " " << lon << " "
	      << x << " " << y << " "
	      << gamma << " " << k << "\n";
  }
  return 0;
}

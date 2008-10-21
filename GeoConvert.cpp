/**
 * \file GeoConvert.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o GeoConvert GeoConvert.cpp GeoCoords.cpp MGRS.cpp UTMUPS.cpp DMS.cpp Constants.cpp TransverseMercator.cpp PolarStereographic.cpp
 *
 **********************************************************************/

#include <iostream>
#include "GeoCoords.hpp"
#include <sstream>
#include <string>
#include <stdexcept>
#include <iomanip>

int main(int argc, char* argv[]) {
  int outputmode = 0;		// Lat/Lon; 1 = DMS; 2 = UTM/UPS; 3 = MGRS
  int prec = 0;
  int zone = -2;		// -2 = track input, -1 = standard

  for (int m = 1; m < argc; ++m) {
    std::string arg = std::string(argv[m]);
    if (arg == "-g")
      outputmode = 0;
    else if (arg == "-d")
      outputmode = 1;
    else if (arg == "-u")
      outputmode = 2;
    else if (arg == "-m")
      outputmode = 3;
    else if (arg == "-c")
      outputmode = 4;
    else if (arg == "-p") {
      std::string a = std::string(argv[++m]);
      std::istringstream str(a);
      str >> prec;
    } else if (arg == "-z") {
      std::string a = std::string(argv[++m]);
      std::istringstream str(a);
      str >> zone;
    } else if (arg == "-s")
      zone = -1;
    else {
      ( arg == "-h" ? std::cout : std::cerr )	<<
"Usage: GeoConvert [-g|-d|-u|-m|-c] [-p prec] [-z zone] [-s] [-h]\n\
$Id$\n\
\n\
Convert geographic coordinates to\n\
\n\
    -g latitude and longitude (decimal degrees), default output\n\
    -d latitude and longitude (degrees mins secs)\n\
    -u UTM or UPS\n\
    -m MGRS\n\
    -c meridian convergence and scale\n\
\n\
Geographic coordinates are given on standard input as\n\
\n\
latitude and longitude (decimal degrees or degrees minutes seconds).\n\
Latitude is given first unless hemisphere is specified, e.g.,\n\
\n\
    33.3 44.4\n\
    E44.4 N33.3\n\
    33d18'N 44d24'E\n\
\n\
UTM or UPS given as zone+hemisphere easting northing or easting northing\n\
zone+hemisphere.  The zone is absent for a UPS specification.  E.g.,\n\
\n\
    38N 444140.54 3684706.36\n\
    444140.54 3684706.36 38N\n\
    S 2173854.98 2985980.58\n\
    2173854.98 2985980.58 S\n\
\n\
MRGS is used to specify the center of a grid square, e.g.,\n\
\n\
    38SMB4484\n\
    38SMB44140847064\n\
\n\
-p prec (default 0) sets the precision relative to 1m.  This gives the\n\
number of digits after the decimal point for UTM/UPS.  The number of digits\n\
per coordinate for MGRS is 5 + prec.  For decimal degrees, the number\n\
of digits after the decimal point is 5 + prec.  For DMS (degree,\n\
minute, seconds) output, the number of digits after the decimal point\n\
in the seconds components is 1 + prec.  The minimum value of prec is -5\n\
and the maximum is 9 for UTM/UPS, 10 for decimal degrees and DMS, and 6 for\n\
MGRS.\n\
\n\
UTM/UPS and MGS are given in zone of the input if applicable, otherwise\n\
in the standard zone.\n\
\n\
-z zone sets the zone for output.  Use zone = 0 to specify a UPS (polar)\n\
zone.\n\
\n\
-s uses the standard zone\n\
\n\
For example, the point\n\
\n\
    79.9S 6.1E\n\
\n\
corresponds to 3 possible MGRS coordinates\n\
\n\
    32CMS4324728161 (standard UTM zone = 32)\n\
    31CEM6066227959 (neighboring UTM zone = 31)\n\
      BBZ1945517770 (neighboring UPS zone)\n\
\n\
then\n\
    echo 31CEM6066227959 | GeoConvert -p -3 -m       ==> 31CEM6027\n\
    echo 31CEM6066227959 | GeoConvert -p -3 -m -s    ==> 32CMS4328\n\
    echo 31CEM6066227959 | GeoConvert -p -3 -m -z 0  ==>   BBZ1917\n\
\n\
-h prints this help\n";
      return arg == "-h" ? 0 : 1;
    }
  }

  GeographicLib::GeoCoords p;
  std::string s;
  std::string os;
  int retval = 0;
  if (!(zone >= -2 && zone <= 60)) {
    std::cerr << "Illegal zone " << zone << "\n";
    return 1;
  }
  while (true) {
    std::getline(std::cin, s);
    if (!std::cin.good())
      break;
    try {
      p.Reset(s);
      if (zone != -2)
	p.SetAltZone(zone);
      switch (outputmode) {
      case 0:
	os = p.GeoRepresentation(prec);
	break;
      case 1:
	os = p.DMSRepresentation(prec);
	break;
      case 2:
	os = p.AltUTMUPSRepresentation(prec);
	break;
      case 3:
	os = p.AltMGRSRepresentation(prec);
	break;
      case 4:
	{
	  double
	    gamma = p.AltConvergence(),
	    k = p.AltScale();
	  std::ostringstream ss;
	  ss << std::fixed
	     << std::setprecision(std::max(0, std::min(8, prec) + 5)) << gamma
	     << " "
	     << std::setprecision(std::max(0, std::min(8, prec) + 7)) << k;
	  os = ss.str();
	}
      }
    }
    catch (std::out_of_range& e) {
      os = std::string("ERROR: ") + e.what();
      retval = 1;
    }
    std::cout << os << std::endl;
  }
  return retval;
}

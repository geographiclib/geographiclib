/**
 * \file GeoConvert.cpp
 * \brief Command line utility for geographic coordinate conversions
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with GeoCoords.o MGRS.o UTMUPS.o DMS.o
 * TransverseMercator.o PolarStereographic.o
 *
 * See \ref geoconvert for usage information.
 **********************************************************************/

#include "GeographicLib/GeoCoords.hpp"
#include "GeographicLib/DMS.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: GeoConvert [-g|-d|-u|-m|-c] [-p prec] [-z zone] [-s] [-t] [-n] [-h]\n\
$Id: GeoConvert.cpp 6827 2010-05-20 19:56:18Z karney $\n\
\n\
Convert geographic coordinates to\n\
\n\
    -g latitude and longitude (decimal degrees), default output\n\
    -d latitude and longitude (degrees mins secs)\n\
    -u UTM or UPS\n\
    -m MGRS\n\
    -c meridian convergence and scale\n\
\n\
The WGS84 model of the earth is used.  Geographic coordinates are given on\n\
standard input as:\n\
\n\
Latitude and longitude (decimal degrees or degrees minutes seconds).  d,\n\
', and \" are used to denote degrees, minutes, and seconds, with the least\n\
significant designator optional.  Latitude is given first unless a\n\
hemisphere is specified, e.g., the following are all equivalent\n\
\n\
    33.3 44.4\n\
    E44.4 N33.3\n\
    33d18'N 44d24'E\n\
    44d24 33d18N\n\
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
per coordinate for MGRS is 5 + prec.  For decimal degrees, the number of\n\
digits after the decimal point is 5 + prec.  For DMS (degree, minute,\n\
seconds) output, the number of digits after the decimal point in the\n\
seconds components is 1 + prec; if this is negative then use minutes (prec\n\
= -2 or -3) or degrees (prec <= -4) as the least significant component.\n\
Print convergence, resp. scale, with 5 + prec, resp. 7 + prec, digits after\n\
the decimal point.  The minimum value of prec is -5 and the maximum is 9\n\
for UTM/UPS, 9 for decimal degrees, 10 for DMS, 6 for MGRS, and 8 for\n\
convergence and scale.\n\
\n\
MGRS coordinates are given by truncating (instead of rounding) the\n\
coordinates to the requested precision.  For example is prec = -3, the\n\
result is the 1km square enclosing the position, for example,\n\
\n\
    echo 38N 444800 3684700 | GeoConvert -m -p -3   ==> 38SMB4484\n\
\n\
If the -n option is given, then, on input, an MSGS coordinate refers to\n\
the south-west corner of the MGRS square instead of the center.  Thus:\n\
\n\
    echo 38SMB4484 | GeoConvert -u         ==> 38N 444500 3684500\n\
    echo 38SMB4484 | GeoConvert -u -n      ==> 38N 444000 3684000\n\
\n\
Convergence is the bearing of grid north given as degrees clockwise from\n\
true north.\n\
\n\
UTM/UPS and MGRS are given in zone of the input if applicable, otherwise in\n\
the standard zone.\n\
\n\
-z zone sets the zone for output.  Use either a positive number for a UTM\n\
zone or zone = 0 to specify UPS.  Alternatively use a zone+hemisphere\n\
designation (and the hemisphere is ignored).\n\
\n\
-s uses the standard UPS and UTM zone boundaries.\n\
\n\
-t is similar to -s but forces UPS regions to the closest UTM zone.\n\
\n\
For example, the point\n\
\n\
    79.9S 6.1E\n\
\n\
corresponds to possible MGRS coordinates\n\
\n\
    32CMS4324728161 (standard UTM zone = 32)\n\
    31CEM6066227959 (neighboring UTM zone = 31)\n\
      BBZ1945517770 (neighboring UPS zone)\n\
\n\
then\n\
    echo 79.9S 6.1E      | GeoConvert -p -3 -m       ==> 32CMS4328\n\
    echo 31CEM6066227959 | GeoConvert -p -3 -m       ==> 31CEM6027\n\
    echo 31CEM6066227959 | GeoConvert -p -3 -m -s    ==> 32CMS4328\n\
    echo 31CEM6066227959 | GeoConvert -p -3 -m -z 0  ==>   BBZ1917\n\
\n\
-h prints this help.\n";
  return retval;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;
  enum { GEOGRAPHIC, DMS, UTMUPS, MGRS, CONVERGENCE };
  int outputmode = GEOGRAPHIC;
  int prec = 0;
  int zone = UTMUPS::MATCH;
  bool centerp = true;

  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-g")
      outputmode = GEOGRAPHIC;
    else if (arg == "-d")
      outputmode = DMS;
    else if (arg == "-u")
      outputmode = UTMUPS;
    else if (arg == "-m")
      outputmode = MGRS;
    else if (arg == "-c")
      outputmode = CONVERGENCE;
    else if (arg == "-n")
      centerp = false;
    else if (arg == "-p") {
      if (++m == argc) return usage(1);
      std::istringstream str(argv[m]);
      char c;
      if (!(str >> prec) || (str >> c)) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
      }
    } else if (arg == "-z") {
      if (++m == argc) return usage(1);
      std::string zonestr(argv[m]);
      try {
        bool northp;
        UTMUPS::DecodeZone(zonestr, zone, northp);
      }
      catch (const std::exception&) {
        std::istringstream str(zonestr);
        char c;
        if (!(str >> zone) || (str >> c)) {
          std::cerr << "Zone " << zonestr
                    << " is not a number or zone+hemisphere\n";
          return 1;
        }
        if (!(zone >= UTMUPS::MINZONE && zone <= UTMUPS::MAXZONE)) {
          std::cerr << "Zone " << zone << " not in [0, 60]\n";
          return 1;
        }
      }
    } else if (arg == "-s")
      zone = UTMUPS::STANDARD;
    else if (arg == "-t")
      zone = UTMUPS::UTM;
    else
      return usage(arg != "-h");
  }

  GeoCoords p;
  std::string s;
  std::string os;
  int retval = 0;

  while (std::getline(std::cin, s)) {
    try {
      p.Reset(s, centerp);
      p.SetAltZone(zone);
      switch (outputmode) {
      case GEOGRAPHIC:
        os = p.GeoRepresentation(prec);
        break;
      case DMS:
        os = p.DMSRepresentation(prec);
        break;
      case UTMUPS:
        os = p.AltUTMUPSRepresentation(prec);
        break;
      case MGRS:
        os = p.AltMGRSRepresentation(prec);
        break;
      case CONVERGENCE:
        {
          real
            gamma = p.AltConvergence(),
            k = p.AltScale();
          os =
            DMS::Encode(gamma, std::max(-5, std::min(8, prec)) + 5, DMS::NUMBER)
            + " " +
            DMS::Encode(k, std::max(-5, std::min(8, prec)) + 7, DMS::NUMBER);
        }
      }
    }
    catch (const std::exception& e) {
      // Write error message cout so output lines match input lines
      os = std::string("ERROR: ") + e.what();
      retval = 1;
    }
    std::cout << os << "\n";
  }
  return retval;
}

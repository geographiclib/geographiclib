/**
 * \file GeoConvert.cpp
 * \brief Command line utility for geographic coordinate conversions
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010, 2011) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with GeoCoords.o MGRS.o UTMUPS.o DMS.o
 * TransverseMercator.o PolarStereographic.o
 *
 * See the <a href="GeoConvert.1.html">man page</a> for usage
 * information.
 *
 * $Id$
 **********************************************************************/

#include "GeographicLib/GeoCoords.hpp"
#include "GeographicLib/DMS.hpp"
#include <sstream>
#include <iostream>

#include "GeoConvert.usage"

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
      if (++m == argc) return usage(1, true);
      std::istringstream str(argv[m]);
      char c;
      if (!(str >> prec) || (str >> c)) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
      }
    } else if (arg == "-z") {
      if (++m == argc) return usage(1, true);
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
    else if (arg == "--version") {
      std::cout << PROGRAM_NAME << ": $Id$\n"
                << "GeographicLib version " << GEOGRAPHICLIB_VERSION << "\n";
      return 0;
    } else
      return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
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

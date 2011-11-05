/**
 * \file GeoConvert.cpp
 * \brief Command line utility for geographic coordinate conversions
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010, 2011) <charles@karney.com>
 * and licensed under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with GeoCoords.o MGRS.o UTMUPS.o DMS.o
 * TransverseMercator.o PolarStereographic.o
 *
 * See the <a href="GeoConvert.1.html">man page</a> for usage
 * information.
 **********************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/GeoCoords.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

#include "GeoConvert.usage"

int main(int argc, char* argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    enum { GEOGRAPHIC, DMS, UTMUPS, MGRS, CONVERGENCE };
    int outputmode = GEOGRAPHIC;
    int prec = 0;
    int zone = UTMUPS::MATCH;
    bool centerp = true, swaplatlong = false;
    std::string istring, ifile, ofile;

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
      else if (arg == "-w")
        swaplatlong = true;
      else if (arg == "-p") {
        if (++m == argc) return usage(1, true);
        try {
          prec = Utility::num<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
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
      else if (arg == "--input-string") {
        if (++m == argc) return usage(1, true);
        istring = argv[m];
      } else if (arg == "--input-file") {
        if (++m == argc) return usage(1, true);
        ifile = argv[m];
      } else if (arg == "--output-file") {
        if (++m == argc) return usage(1, true);
        ofile = argv[m];
      } else if (arg == "--version") {
        std::cout
          << argv[0]
          << ": $Id$\n"
          << "GeographicLib version " << GEOGRAPHICLIB_VERSION_STRING << "\n";
        return 0;
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
    }

    if (!ifile.empty() && !istring.empty()) {
      std::cerr << "Cannot specify --input-string and --input-file together\n";
      return 1;
    }
    if (ifile == "-") ifile.clear();
    std::ifstream infile;
    std::istringstream instring;
    if (!ifile.empty()) {
      infile.open(ifile.c_str());
      if (!infile.is_open()) {
        std::cerr << "Cannot open " << ifile << " for reading\n";
        return 1;
      }
    } else if (!istring.empty()) {
      std::string::size_type m = 0;
      while (true) {
        m = istring.find(';', m);
        if (m == std::string::npos)
          break;
        istring[m] = '\n';
      }
      instring.str(istring);
    }
    std::istream* input = !ifile.empty() ? &infile :
      (!istring.empty() ? &instring : &std::cin);

    std::ofstream outfile;
    if (ofile == "-") ofile.clear();
    if (!ofile.empty()) {
      outfile.open(ofile.c_str());
      if (!outfile.is_open()) {
        std::cerr << "Cannot open " << ofile << " for writing\n";
        return 1;
      }
    }
    std::ostream* output = !ofile.empty() ? &outfile : &std::cout;

    GeoCoords p;
    std::string s;
    std::string os;
    int retval = 0;

    while (std::getline(*input, s)) {
      try {
        p.Reset(s, centerp, swaplatlong);
        p.SetAltZone(zone);
        switch (outputmode) {
        case GEOGRAPHIC:
          os = p.GeoRepresentation(prec, swaplatlong);
          break;
        case DMS:
          os = p.DMSRepresentation(prec, swaplatlong);
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
              Utility::str<real>(gamma, std::max(-5,std::min(8,prec))+5)
              + " " +
              Utility::str<real>(k, std::max(-5,std::min(8,prec))+7);
          }
        }
      }
      catch (const std::exception& e) {
        // Write error message to cout so output lines match input lines
        os = std::string("ERROR: ") + e.what();
        retval = 1;
      }
      *output << os << "\n";
    }
    return retval;
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
}

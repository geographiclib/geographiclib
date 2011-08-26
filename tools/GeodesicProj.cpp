/**
 * \file GeodesicProj.cpp
 * \brief Command line utility for geodesic projections
 *
 * Copyright (c) Charles Karney (2009, 2010, 2011) <charles@karney.com>
 * and licensed under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geodesic.o GeodesicLine.o
 * AzimuthalEquidistant.o Gnomonic.o CassiniSoldner.o
 *
 * See the <a href="GeodesicProj.1.html">man page</a> for usage
 * information.
 **********************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/AzimuthalEquidistant.hpp>
#include <GeographicLib/CassiniSoldner.hpp>
#include <GeographicLib/Gnomonic.hpp>
#include <GeographicLib/DMS.hpp>

#include "GeodesicProj.usage"

int main(int argc, char* argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    bool azimuthal = false, cassini = false, gnomonic = false, reverse = false;
    real lat0 = 0, lon0 = 0;
    real
      a = Constants::WGS84_a<real>(),
      f = Constants::WGS84_f<real>();
    std::string istring, ifile, ofile;

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-r")
        reverse = true;
      else if (arg == "-c" || arg == "-z" || arg == "-g") {
        cassini = arg == "-c";
        azimuthal = arg == "-z";
        gnomonic = arg == "-g";
        if (m + 2 >= argc) return usage(1, true);
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
        if (m + 2 >= argc) return usage(1, true);
        try {
          a = DMS::Decode(std::string(argv[m + 1]));
          f = DMS::DecodeFraction(std::string(argv[m + 2]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        m += 2;
      } else if (arg == "--input-string") {
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

    if (!(azimuthal || cassini || gnomonic)) {
      std::cerr << "Must specify \"-z lat0 lon0\" or "
                << "\"-c lat0 lon0\" or \"-g lat0 lon0\"\n";
      return 1;
    }

    const Geodesic geod(a, f);
    const CassiniSoldner cs = cassini ?
      CassiniSoldner(lat0, lon0, geod) : CassiniSoldner(geod);
    const AzimuthalEquidistant az(geod);
    const Gnomonic gn(geod);

    std::string s;
    int retval = 0;
    std::cout << std::fixed;
    while (std::getline(*input, s)) {
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
          *output << DMS::Encode(lat, 15, DMS::NUMBER) << " "
                  << DMS::Encode(lon, 15, DMS::NUMBER) << " "
                  << DMS::Encode(azi, 15, DMS::NUMBER) << " "
                  << DMS::Encode(rk, 16, DMS::NUMBER) << "\n";
        } else {
          if (cassini)
            cs.Forward(lat, lon, x, y, azi, rk);
          else if (azimuthal)
            az.Forward(lat0, lon0, lat, lon, x, y, azi, rk);
          else
            gn.Forward(lat0, lon0, lat, lon, x, y, azi, rk);
          *output << DMS::Encode(x, 10, DMS::NUMBER) << " "
                  << DMS::Encode(y, 10, DMS::NUMBER) << " "
                  << DMS::Encode(azi, 15, DMS::NUMBER) << " "
                  << DMS::Encode(rk, 16, DMS::NUMBER) << "\n";
        }
      }
      catch (const std::exception& e) {
        *output << "ERROR: " << e.what() << "\n";
        retval = 1;
      }
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

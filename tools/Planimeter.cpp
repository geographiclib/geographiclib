/**
 * \file Planimeter.cpp
 * \brief Command line utility for measuring the area of geodesic polygons
 *
 * Copyright (c) Charles Karney (2010, 2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geodesic.o GeodesicLine.o DMS.o
 *
 * See the <a href="Planimeter.1.html">man page</a> for usage
 * information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/PolygonArea.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/GeoCoords.hpp>

#include "Planimeter.usage"

int main(int argc, char* argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    real
      a = Constants::WGS84_a<real>(),
      f = Constants::WGS84_f<real>();
    bool reverse = false, sign = true, polyline = false;
    std::string istring, ifile, ofile;

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-r")
        reverse = !reverse;
      else if (arg == "-s")
        sign = !sign;
      else if (arg == "-l")
        polyline = !polyline;
      else if (arg == "-e") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          a = Utility::num<real>(std::string(argv[m + 1]));
          f = Utility::fract<real>(std::string(argv[m + 2]));
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

    const Geodesic geod(a, f);
    PolygonArea poly(geod, polyline);
    GeoCoords p;

    std::string s;
    real perimeter, area;
    unsigned num;
    while (std::getline(*input, s)) {
      bool endpoly = s.empty();
      if (!endpoly) {
        try {
          p.Reset(s);
          if (Math::isnan(p.Latitude()) || Math::isnan(p.Longitude()))
            endpoly = true;
        }
        catch (const GeographicErr&) {
          endpoly = true;
        }
      }
      if (endpoly) {
        num = poly.Compute(reverse, sign, perimeter, area);
        if (num > 0) {
          *output << num << " "
                  << Utility::str<real>(perimeter, 8);
          if (!polyline)
            *output << " " << Utility::str<real>(area, 3);
          *output << "\n";
        }
        poly.Clear();
      } else
        poly.AddPoint(p.Latitude(), p.Longitude());
    }
    num = poly.Compute(reverse, sign, perimeter, area);
    if (num > 0) {
      *output << num << " "
              << Utility::str<real>(perimeter, 8);
      if (!polyline)
        *output << " " << Utility::str<real>(area, 3);
      *output << "\n";
    }
    poly.Clear();
    return 0;
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

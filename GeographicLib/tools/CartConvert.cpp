/**
 * \file CartConvert.cpp
 * \brief Command line utility for geodetic to cartesian coordinate conversions
 *
 * Copyright (c) Charles Karney (2009, 2010, 2011) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geocentric.o LocalCartesian.o
 *
 * See the <a href="CartConvert.1.html">man page</a> for usage
 * information.
 **********************************************************************/

#include "GeographicLib/Geocentric.hpp"
#include "GeographicLib/LocalCartesian.hpp"
#include "GeographicLib/DMS.hpp"
#include <iostream>
#include <sstream>

#include "CartConvert.usage"

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;
  bool localcartesian = false, reverse = false;
  real
    a = Constants::WGS84_a<real>(),
    r = Constants::WGS84_r<real>();
  real lat0 = 0, lon0 = 0, h0 = 0;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-r")
      reverse = true;
    else if (arg == "-l") {
      localcartesian = true;
      if (m + 3 >= argc) return usage(1, true);
      try {
        DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                          lat0, lon0);
        h0 = DMS::Decode(std::string(argv[m + 3]));
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding arguments of -l: " << e.what() << "\n";
        return 1;
      }
      m += 3;
    } else if (arg == "-e") {
      if (m + 2 >= argc) return usage(1, true);
      try {
        a = DMS::Decode(std::string(argv[m + 1]));
        r = DMS::Decode(std::string(argv[m + 2]));
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
        return 1;
      }
      m += 2;
    } else if (arg == "--version") {
      std::cout
        << PROGRAM_NAME
        << ": $Id$\n"
        << "GeographicLib version " << GEOGRAPHICLIB_VERSION << "\n";
      return 0;
    } else
      return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
  }

  const Geocentric ec(a, r);
  const LocalCartesian lc(lat0, lon0, h0, ec);

  std::string s;
  int retval = 0;
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      real lat, lon, h, x, y, z;
      std::string stra, strb, strc;
      if (!(str >> stra >> strb >> strc))
        throw  GeographicErr("Incomplete input: " + s);
      if (reverse) {
        x = DMS::Decode(stra);
        y = DMS::Decode(strb);
        z = DMS::Decode(strc);
      } else {
        DMS::DecodeLatLon(stra, strb, lat, lon);
        h = DMS::Decode(strc);
      }
      std::string strd;
      if (str >> strd)
        throw GeographicErr("Extraneous input: " + strd);
      if (reverse) {
        if (localcartesian)
          lc.Reverse(x, y, z, lat, lon, h);
        else
          ec.Reverse(x, y, z, lat, lon, h);
        std::cout << DMS::Encode(lat, 15, DMS::NUMBER) << " "
                  << DMS::Encode(lon, 15, DMS::NUMBER) << " "
                  << DMS::Encode(h, 12, DMS::NUMBER) << "\n";
      } else {
        if (localcartesian)
          lc.Forward(lat, lon, h, x, y, z);
        else
          ec.Forward(lat, lon, h, x, y, z);
        std::cout << DMS::Encode(x, 10, DMS::NUMBER) << " "
                  << DMS::Encode(y, 10, DMS::NUMBER) << " "
                  << DMS::Encode(z, 10, DMS::NUMBER) << "\n";
      }
    }
    catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }

  return retval;
}

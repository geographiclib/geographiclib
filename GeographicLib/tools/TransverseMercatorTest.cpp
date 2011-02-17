/**
 * \file TransverseMercatorTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010, 2011) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with TransverseMercatorExact.o
 * EllipticFunction.o TransverseMercator.o
 *
 * See the <a href="TransverseMercatorTest.1.html">man page</a> for usage
 * information.
 *
 * $Id$
 **********************************************************************/

#include "GeographicLib/EllipticFunction.hpp"
#include "GeographicLib/TransverseMercatorExact.hpp"
#include "GeographicLib/TransverseMercator.hpp"
#include "GeographicLib/DMS.hpp"
#include <iostream>
#include <sstream>

#include "TransverseMercatorTest.usage"

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;
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
      return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
  }

  real e, a;
  if (testing) {
    e = 1/real(10);
    EllipticFunction temp(e * e);
    a = 1/temp.E();
  }
  const TransverseMercatorExact& TME = testing ?
    TransverseMercatorExact(a, (std::sqrt(1 - e * e) + 1) / (e * e),
                            real(1), true) :
    TransverseMercatorExact::UTM;

  const TransverseMercator& TMS = TransverseMercator::UTM;

  std::string s;
  int retval = 0;
  std::cout << std::fixed;
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      real lat, lon, x, y;
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
      real gamma, k;
      if (reverse) {
        if (series)
          TMS.Reverse(real(0), x, y, lat, lon, gamma, k);
        else
          TME.Reverse(real(0), x, y, lat, lon, gamma, k);
      } else {
        if (series)
          TMS.Forward(real(0), lat, lon, x, y, gamma, k);
        else
          TME.Forward(real(0), lat, lon, x, y, gamma, k);
      }
      std::cout << DMS::Encode(lat, 15, DMS::NUMBER) << " "
                << DMS::Encode(lon, 15, DMS::NUMBER) << " "
                << DMS::Encode(x, 10, DMS::NUMBER) << " "
                << DMS::Encode(y, 10, DMS::NUMBER) << " "
                << DMS::Encode(gamma, 16, DMS::NUMBER) << " "
                << DMS::Encode(k, 16, DMS::NUMBER) << "\n";
    }
    catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }

  return retval;
}

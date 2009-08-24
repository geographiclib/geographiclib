/**
 * \file Geoid.cpp
 * \brief Command line utility for evaluation geoid heights
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o GeoidEval GeoidEval.cpp Geoid.cpp DMS.cpp
 *
 * See \ref geod for usage information.
 **********************************************************************/

#include "GeographicLib/Geoid.hpp"
#include "GeographicLib/DMS.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: Geoid filename\n\
$Id$\n\
\n";
  return retval;
}

int main(int argc, char* argv[]) {
  if (argc != 2)
    return usage(1);
  std::string geoid = std::string(argv[1]);
  Geoid gx(geoid);
  gx.CacheAll();
  std::cout << std::fixed << std::setprecision(4);
  std::string s;
  int retval = 0;
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream  str(s);
      std::string stra, strb;
      if (!(str >> stra >> strb))
	throw std::out_of_range("Incomplete input: " + s);
      double lat, lon;
      GeographicLib::DMS::DecodeLatLon(stra, strb, lat, lon);
      double h = gx(lat, lon);
      std::cout << h << "\n";
    }
    catch (std::exception& e) {
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }
  return retval;
}

/*
  Errors:
  EGM84 MAX RMS
  30' 1.546 0.070
  15' 0.470 0.018

  EGM96
  15' grid (vs 3.75') 1.152 0.040
  5' grid (vs 2.5')   0.140 0.005

  EMG2009
  15' grid (vs 2.5') 2.670 0.079*
  5' grid (vs 2.5')  0.478 0.012*
  2.5' grid (vs 1')  0.143 0.003*
  1' grid (vs 2.5')  0.021 0.001*


  
0.1314 => 1.085 * 0.1314 = 


degree:ev(%pi/180,0)$
rnd():=bfloat(random(1.0))+0.5^53$
quant(x):=entier(1e3*x+0.5b0)*1e-3$
rndlat():=quant(asin(2*rnd()-1)/degree);
rndlon():=quant(180*(2*rnd()-1))$
set_random_state(make_random_state(true))$

for i:1 thru 100 do print(rndlat(), rndlon())$
*/

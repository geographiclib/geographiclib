/**
 * \file Geod.cpp
 * \brief Command line utility for geodesic calculations
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 *
 * Compile with
 *
 *   g++ -g -O3 -I.. -o Geod Geod.cpp Geodesic.cpp Constants.cpp
 *
 * See \ref geod for usage information.
 **********************************************************************/

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/DMS.hpp"

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: Geod [-i | -l lat1 lon1 head1] [-n] [-d] [-f] [-p prec] [-h]\n\
$Id$\n\
\n\
-i: inverse instead of direct\n\
-l: line calc\n\
-n: international ellipsoid\n\
-d: DMS\n\
-f: full geoid line\n\
-p: prec relative to 1m\n\
-h prints this help\n";
  return retval;
}


std::string LatLonString(double lat, double lon, int prec, bool dms) {
  using namespace GeographicLib;
  if (dms)
    return
      DMS::Encode(lat, prec + 5, DMS::LATITUDE) + " " +
      DMS::Encode(lon, prec + 5, DMS::LONGITUDE);
  else {
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec + 5)
       << lat << " " << lon;
    return os.str();
  }
}

std::string AzimuthString(double head, int prec, bool dms) {
  using namespace GeographicLib;
  if (dms)
    return DMS::Encode(head, prec + 5, DMS::LONGITUDE);
  else {
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec + 5)
       << head;
    return os.str();
  }
}

double ReadAzimuth(const std::string& s) {
  using namespace GeographicLib;
  DMS::flag ind;
  double head = DMS::Decode(s, ind);
  if (!(head >= -180 && head <= 360))
    throw std::out_of_range("Azimuth " + s + " not in range [-180,360]");
  if (head >= 180) head -= 360;
  if (ind == DMS::LATITUDE)
    throw std::out_of_range("Azimuth " + s +
			    " has a latitude hemisphere, N/S");
  return head;
}

int main(int argc, char* argv[]) {
  bool linecalc = false, inverse = false, international = false,
    dms = false, full = false;
  double lat1, lon1, head1, lat2, lon2, head2, s12;
  int prec = 3;
  for (int m = 1; m < argc; ++m) {
    std::string arg = std::string(argv[m]);
    if (arg == "-i") {
      inverse = true;
      linecalc = false;
    } else if (arg == "-n")
      international = true;
    else if (arg == "-d")
      dms = true;
    else if (arg == "-f")
      full = true;
    else if (arg == "-p") {
      if (++m == argc) return usage(1);
      std::string a = std::string(argv[m]);
      std::istringstream str(a);
      if (!(str >> prec)) return usage(1);
    } else if (arg == "-l") {
      inverse = false;
      linecalc = true;
      if (m + 3 >= argc) return usage(1);
      try {
	GeographicLib::DMS::DecodeLatLon(std::string(argv[m + 1]),
					 std::string(argv[m + 2]),
					 lat1, lon1);
	head1 = ReadAzimuth(std::string(argv[m + 3]));
	m += 3;
      }
      catch (std::out_of_range& e) {
	std::cerr << "ERROR: " << e.what() << "\n";
	return usage(1);
      }
    } else
      return usage(arg != "-h");
  }

  /* FIX ME */
  const GeographicLib::Geodesic internat(6378388.0, 297.0);
  const GeographicLib::Geodesic& geod = international ? internat :
    GeographicLib::Geodesic::WGS84;
  GeographicLib::GeodesicLine l;
  if (linecalc)
    l = geod.Line(lat1, lon1, head1);

  // Max precision = 9: 1 nm in distance, 10^-14 deg (= 1.1 nm),
  // 10^-10 sec (= 3 nm).
  prec = std::min(9, std::max(0, prec));
  std::cout << std::fixed << std::setprecision(prec);
  std::string s;
  int retval = 0;
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      if (linecalc) {
	if (!(str >> s12))
	  throw std::out_of_range("Incomplete input: " + s);
	l.Position(s12, lat2, lon2, head2);
	if (full)
	  std::cout << LatLonString(lat1, lon1, prec, dms) << " " <<
	    AzimuthString(head1, prec, dms) << " ";
	std::cout << LatLonString(lat2, lon2, prec, dms) << " " <<
	  AzimuthString(head2, prec, dms);
	if (full)
	  std::cout << " " << s12;
	std::cout << "\n";
      } else if (inverse) {
	std::string slat1, slon1, slat2, slon2;
	if (!(str >> slat1 >> slon1 >> slat2 >> slon2))
	  throw std::out_of_range("Incomplete input: " + s);
	GeographicLib::DMS::DecodeLatLon(slat1, slon1, lat1, lon1);
	GeographicLib::DMS::DecodeLatLon(slat2, slon2, lat2, lon2);
	geod.Inverse(lat1, lon1, lat2, lon2, s12, head1, head2);
	std::cerr << geod.itera << " " << geod.iterb << " ";
	if (full)
	  std::cout << LatLonString(lat1, lon1, prec, dms) << " ";
	std::cout << AzimuthString(head1, prec, dms) << " ";
	if (full)
	  std::cout << LatLonString(lat2, lon2, prec, dms) << " ";
	std::cout << AzimuthString(head2, prec, dms) << " " << s12 << "\n";
      } else {
	std::string slat1, slon1, shead1;
	if (!(str >> slat1 >> slon1 >> shead1 >> s12))
	  throw std::out_of_range("Incomplete input: " + s);
	GeographicLib::DMS::DecodeLatLon(slat1, slon1, lat1, lon1);
	head1 = ReadAzimuth(shead1);
	geod.Direct(lat1, lon1, head1, s12, lat2, lon2, head2);
	if (full)
	  std::cout << LatLonString(lat1, lon1, prec, dms) << " " <<
	    AzimuthString(head1, prec, dms) << " ";
	std::cout << LatLonString(lat2, lon2, prec, dms) << " " <<
	  AzimuthString(head2, prec, dms);
	if (full)
	  std::cout << " " << s12;
	std::cout << "\n";
      }
    }
    catch (std::out_of_range& e) {
      // Write error message cout so output lines match input lines
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }
  return retval;
}

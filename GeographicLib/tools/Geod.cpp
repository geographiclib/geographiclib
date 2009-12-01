/**
 * \file Geod.cpp
 * \brief Command line utility for geodesic calculations
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geodesic.o DMS.o
 *
 * See \ref geod for usage information.
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/DMS.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: Geod [-l lat1 lon1 azi1 | -i] [-a] [-n | -e a r]\n\
            [-d] [-b] [-f] [-p prec] [-h]\n\
$Id$\n\
\n\
Perform geodesic calculations.\n\
\n\
The shortest path between two points on the ellipsoid at (lat1, lon1) and\n\
(lat2, lon2) is called the geodesic.  Its length is s12 and the geodesic\n\
from point 1 to point 2 has azimuths azi1 and azi2 at the two end\n\
points.  The reduced length of the geodesic, m12, is defined such that\n\
if the initial azimuth is perturbed by dazi1 (radians) then the second\n\
point is displaced by m12*dazi1 in the direction perpendicular to the\n\
geodesic.  On a flat surface, we have m12 = s12.\n\
\n\
Geod operates in one of three modes:\n\
\n\
(1) It accepts lines on the standard input containing \"lat1 lon1 azi1\n\
    s12\" and prints \"lat2 lon2 azi2 m12\" on standard output.  This is\n\
    the direct geodesic calculation.\n\
\n\
(2) Command line arguments \"-l lat1 lon1 azi1\" specify a geodesic line.\n\
    Geod then accepts a sequence of s12 values (one per line) on\n\
    standard input and prints \"lat2 lon2 azi2 m12\" for each.  This\n\
    generates a sequence of points on a single geodesic.\n\
\n\
(3) With the -i command line argument, Geod performs the inverse\n\
    geodesic calculation.  It reads lines containing \"lat1 lon1 lat2\n\
    lon2\" and prints the corresponding values of \"azi1 azi2 s12 m12\".\n\
\n\
By default, the WGS84 ellipsoid is used.  Specifying \"-e a r\" sets the\n\
equatorial radius of the ellipsoid to \"a\" and the reciprocal flattening\n\
to r.  Setting r = 0 results in a sphere.  Specify r < 0 for a prolate\n\
ellipsoid.  The -n option uses the international ellipsoid (equivalent to\n\
\"-e 6378388 297\").\n\
\n\
Output of angles is as decimal degrees.  If -d is specified the output\n\
is as degrees, minutes, seconds.  Input can be in either style.  d, ',\n\
and \" are used to denote degrees, minutes, and seconds, with the least\n\
significant designator optional.  By default, latitude precedes\n\
longitude for each point; however on input either may be given first by\n\
appending N or S to the latitude and E or W to the longitude.  Azimuths\n\
(measured clockwise from north) give the heading of the geodesic.  The\n\
azimuth azi2 is the forward azimuth (the heading beyond point 2).  If\n\
the -b flag is given, azi2 is converted to a back azimuth (the direction\n\
back to point 1) for output.\n\
\n\
s12 is given in meters, unless the -a flag is given.  In that case, s12\n\
(on both input and output) are given as the arc length on the auxiliary\n\
sphere a12 (measured in degrees).  In these terms, 180 degrees is the\n\
distance from one equator crossing to the next or from the minimum\n\
latitude to the maximum latitude.  Distances greater than 180 degrees do\n\
not correspond to shortest paths.  m12 is always given in meters.\n\
\n\
The output lines consist of the four quantities needed to complete the\n\
specification of the geodesic.  With the -f option, each line of output\n\
is a complete geodesic specification consisting of nine quantities\n\
\n\
    lat1 lon1 azi1 lat2 lon2 azi2 s12 a12 m12\n\
\n\
where here s12 is the distance and a12 the arc length.\n\
\n\
-p prec (default 3) gives the precision of the output relative to 1m.\n\
The minimum value of prec is 0 (1 m accuracy) and the maximum value is\n\
10 (0.1 nm accuracy, but then the last digits are unreliable).\n\
\n\
-h prints this help.\n";
  return retval;
}

typedef GeographicLib::Math::real real;

std::string LatLonString(real lat, real lon, int prec, bool dms) {
  using namespace GeographicLib;
  if (dms)
    return
      DMS::Encode(lat, prec + 5, DMS::LATITUDE) + " " +
      DMS::Encode(lon, prec + 5, DMS::LONGITUDE);
  else {
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec + 5) << lat << " " << lon;
    return os.str();
  }
}

std::string AzimuthString(real azi, int prec, bool dms) {
  using namespace GeographicLib;
  if (dms)
    return DMS::Encode(azi, prec + 5, DMS::AZIMUTH);
  else {
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec + 5)
       << (azi >= 180 ? azi - 360 : azi);
    return os.str();
  }
}

real ReadAzimuth(const std::string& s) {
  using namespace GeographicLib;
  DMS::flag ind;
  real azi = DMS::Decode(s, ind);
  if (!(azi >= -180 && azi <= 360))
    throw std::out_of_range("Azimuth " + s + " not in range [-180,360]");
  if (azi >= 180) azi -= 360;
  if (ind == DMS::LATITUDE)
    throw std::out_of_range("Azimuth " + s
                            + " has a latitude hemisphere, N/S");
  return azi;
}

std::string DistanceStrings(real s12, real a12,
                            bool full, bool arcmode, int prec, bool dms) {
  using namespace GeographicLib;
  std::ostringstream os;
  if (full || !arcmode)
    os << std::fixed << std::setprecision(prec) << s12;
  if (full)
    os << " ";
  if (full || arcmode) {
    if (dms)
      os << DMS::Encode(a12, prec + 5, DMS::NONE);
    else
      os << std::fixed << std::setprecision(prec + 5) << a12;
  }
  return os.str();
}

real ReadDistance(const std::string& s, bool arcmode) {
  using namespace GeographicLib;
  real s12;
  if (arcmode) {
    DMS::flag ind;
    s12 = DMS::Decode(s, ind);
    if (ind != DMS::NONE)
      throw std::out_of_range("Arc angle " + s
                              + " includes a hemisphere, N/E/W/S");
  } else {
    std::istringstream is(s);
    if (!(is >> s12))
      throw std::out_of_range("Could not read distance: " + s);
    // is >> s12 gobbles final E in 1234E, so look for last character which is
    // legal as the final character in a number (digit or period).
    int pos = std::min(int(is.tellg()), int(s.find_last_of("0123456789.")) + 1);
    if (pos != int(s.size()))
      throw std::out_of_range("Extra text "
                         + s.substr(pos) + " in distance " + s);
  }
  return s12;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  bool linecalc = false, inverse = false, arcmode = false,
    dms = false, full = false;
  real
    a = GeographicLib::Constants::WGS84_a(),
    r = GeographicLib::Constants::WGS84_r();
  real lat1, lon1, azi1, lat2, lon2, azi2, s12, m12, a12;
  real azi2sense = 0;
  int prec = 3;

  for (int m = 1; m < argc; ++m) {
    std::string arg = std::string(argv[m]);
    if (arg == "-i") {
      inverse = true;
      linecalc = false;
    } else if (arg == "-a")
      arcmode = true;
    else if (arg == "-l") {
      inverse = false;
      linecalc = true;
      if (m + 3 >= argc) return usage(1);
      try {
        GeographicLib::DMS::DecodeLatLon(std::string(argv[m + 1]),
                                         std::string(argv[m + 2]),
                                         lat1, lon1);
        azi1 = ReadAzimuth(std::string(argv[m + 3]));
        m += 3;
      }
      catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return usage(1);
      }
    } else if (arg == "-n") {
      a = 6378388;
      r = 297;
    } else if (arg == "-e") {
      for (unsigned i = 0; i < 2; ++i) {
        if (++m == argc) return usage(1);
        std::string s = std::string(argv[m]);
        std::istringstream str(s);
        if (!(str >> (i ? r : a))) return usage(1);
      }
    }
    else if (arg == "-d")
      dms = true;
    else if (arg == "-b")
      azi2sense = 180;
    else if (arg == "-f")
      full = true;
    else if (arg == "-p") {
      if (++m == argc) return usage(1);
      std::string s = std::string(argv[m]);
      std::istringstream str(s);
      if (!(str >> prec)) return usage(1);
    } else
      return usage(arg != "-h");
  }

  const GeographicLib::Geodesic geod(a, r);
  GeographicLib::GeodesicLine l;
  if (linecalc)
    l = geod.Line(lat1, lon1, azi1);

  // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
  // 10^-11 sec (= 0.3 nm).
  prec = std::min(10, std::max(0, prec));
  std::cout << std::fixed << std::setprecision(prec);
  std::string s;
  int retval = 0;
  while (std::getline(std::cin, s)) {
    try {
      std::istringstream str(s);
      if (inverse) {
        std::string slat1, slon1, slat2, slon2;
        if (!(str >> slat1 >> slon1 >> slat2 >> slon2))
          throw std::out_of_range("Incomplete input: " + s);
        std::string strc;
        if (str >> strc)
          throw std::out_of_range("Extraneous input: " + strc);
        GeographicLib::DMS::DecodeLatLon(slat1, slon1, lat1, lon1);
        GeographicLib::DMS::DecodeLatLon(slat2, slon2, lat2, lon2);
        a12 = geod.Inverse(lat1, lon1, lat2, lon2, s12, azi1, azi2, m12);
        if (full)
          std::cout << LatLonString(lat1, lon1, prec, dms) << " ";
        std::cout << AzimuthString(azi1, prec, dms) << " ";
        if (full)
          std::cout << LatLonString(lat2, lon2, prec, dms) << " ";
        std::cout << AzimuthString(azi2 + azi2sense, prec, dms) << " "
                  << DistanceStrings(s12, a12, full, arcmode, prec, dms)
                  << " " << m12 << "\n";
      } else {
        if (linecalc) {
          std::string ss12;
          if (!(str >> ss12))
            throw std::out_of_range("Incomplete input: " + s);
          std::string strc;
          if (str >> strc)
            throw std::out_of_range("Extraneous input: " + strc);
          s12 = ReadDistance(ss12, arcmode);
          a12 = l.Position(s12, lat2, lon2, azi2, m12, arcmode);
        } else {
          std::string slat1, slon1, sazi1, ss12;
          if (!(str >> slat1 >> slon1 >> sazi1 >> ss12))
            throw std::out_of_range("Incomplete input: " + s);
          std::string strc;
          if (str >> strc)
            throw std::out_of_range("Extraneous input: " + strc);
          GeographicLib::DMS::DecodeLatLon(slat1, slon1, lat1, lon1);
          azi1 = ReadAzimuth(sazi1);
          s12 = ReadDistance(ss12, arcmode);
          a12 =
            geod.Direct(lat1, lon1, azi1, s12, lat2, lon2, azi2, m12, arcmode);
        }
        if (arcmode)
          std::swap(s12, a12);
        if (full)
          std::cout << LatLonString(lat1, lon1, prec, dms) << " "
                    << AzimuthString(azi1, prec, dms) << " ";
        std::cout << LatLonString(lat2, lon2, prec, dms) << " "
                  << AzimuthString(azi2 + azi2sense, prec, dms);
        if (full)
          std::cout << " "
                    << DistanceStrings(s12, a12, full, arcmode, prec, dms);
        std::cout << " " << m12 << "\n";
      }
    }
    catch (const std::exception& e) {
      // Write error message cout so output lines match input lines
      std::cout << "ERROR: " << e.what() << "\n";
      retval = 1;
    }
  }
  return retval;
}

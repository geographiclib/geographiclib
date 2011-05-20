/**
 * \file GeoidEval.cpp
 * \brief Command line utility for evaluation geoid heights
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geoid.o DMS.o
 *
 * See \ref geoideval for usage information.
 **********************************************************************/

#include "GeographicLib/Geoid.hpp"
#include "GeographicLib/DMS.hpp"
#include "GeographicLib/GeoCoords.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>

int usage(int retval) {
  using namespace GeographicLib;
  std::string
    geoidpath = Geoid::GeoidPath(),
    defaultpath = Geoid::DefaultPath();
  if (geoidpath.empty())
    geoidpath = "UNDEFINED";
  ( retval ? std::cerr : std::cout ) <<
"Usage:\n\
  GeoidEval [-n name] [-d dir] [-l] [-a] [-c south west north east] \\\n\
            [-z zone] [-msltohae] [-haetomsl] [-v] [-h]\n\
$Id: GeoidEval.cpp 6867 2010-09-11 13:04:26Z karney $\n\
\n\
Read in positions on standard input and print out the corresponding\n\
geoid heights on standard output.  In addition print the northerly and\n\
easterly gradients of the geoid height (i.e., the rate at which the\n\
geoid height changes per unit distance along the WGS84 ellipsoid in\n\
the specified directions).\n\
\n\
Positions are given as latitude and longitude, UTM/UPS, or MGRS, in\n\
any of the formats accepted by GeoConvert.  (MGRS coordinates signify\n\
the center of the corresponing MGRS square.)  If the -z option is\n\
specified then the specified zone is prepended to each line of input\n\
(which must be in UTM/UPS coordinates).  This allows a file with UTM\n\
eastings and northings in a single zone to be used as standard input.\n\
\n\
With the -msltohae or -haetomsl options, each line of input should\n\
also include a height (in meters) as the last item.  In this case the\n\
height is converted in the indicated direction and the input line is\n\
echoed to standard output with the height adjusted, thereby allowing\n\
vertical datum of 3d data to be changed.\n\
\n\
By default the EGM96 geoid is used with a 5\' grid.  This may be\n\
overriden with the -n option.  The name specified should be one of\n\
\n\
                                  bilinear error    cubic error\n\
   name         geoid    grid     max     rms       max     rms\n\
   egm84-30     EGM84    30\'      1.546m  70mm      0.274m  14mm\n\
   egm84-15     EGM84    15\'      0.413m  18mm      0.020m   1mm\n\
   egm96-15     EGM96    15\'      1.152m  40mm      0.169m   7mm\n\
   egm96-5      EGM96     5\'      0.140m   5mm      0.003m   1mm\n\
   egm2008-5    EGM2008   5\'      0.478m  12mm      0.294m   5mm\n\
   egm2008-2_5  EGM2008   2.5\'    0.135m   3mm      0.031m   1mm\n\
   egm2008-1    EGM2008   1\'      0.025m   1mm      0.003m   1mm\n\
\n\
(Some of the geoids may not be available.)  The errors listed here\n\
are estimates of the quantization and interpolation errors in the\n\
reported heights compared to the specified geoid.\n\
\n\
Cubic interpolation is used to compute the geoid height unless\n\
-l is specified in which case bilinear interpolation is used.\n\
Cubic interpolation is more accurate; however it results in\n\
small discontinuities in the returned height on cell boundaries.\n\
The gradients are computed by differentiating the interpolated\n\
results.\n\
\n\
GeoidEval will load the geoid data from the directory specified by\n\
the -d option.  If this is not provided, it will look up the value of\n\
GEOID_PATH (currently " << geoidpath << ") in the environment.\n\
If this is not defined, it will use the compile-time value of\n"
<< defaultpath << ".\n\
\n\
By default, the data file is randomly read to compute the geoid\n\
heights at the input positions.  Usually this is sufficient for\n\
interactive use.  If many heights are to be computed, GeoidEval\n\
allows a block of data to be read into memory and heights within the\n\
corresponding rectangle can then be computed without any disk access.\n\
If -a is specified all the geoid data is read; in the case of\n\
egm2008-1, this requires about 0.5 GB of RAM.  The -c option allows\n\
a rectangle of data to be cached.  The evaluation of heights\n\
outside the cached rectangle causes the necessary data to be read\n\
from disk.\n\
\n\
Regardless of whether any cache is requested (with the -a or -c\n\
options), the data for the last grid cell in cached.  This allows\n\
the geoid height along a continuous path to be returned with little\n\
disk overhead.\n\
\n\
The -v option causes the data about the current geoid to be printed\n\
to standard error.\n\
\n\
-h prints this help.\n";
  return retval;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;
  bool cacheall = false, cachearea = false, verbose = false, cubic = true;
  real caches, cachew, cachen, cachee;
  std::string dir;
  std::string geoid = "egm96-5";
  std::string zone;
  Geoid::convertflag heightmult = Geoid::NONE;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-a") {
      cacheall = true;
      cachearea = false;
    }
    else if (arg == "-c") {
      if (m + 4 >= argc) return usage(1);
      cacheall = false;
      cachearea = true;
      try {
        DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                          caches, cachew);
        DMS::DecodeLatLon(std::string(argv[m + 3]), std::string(argv[m + 4]),
                          cachen, cachee);
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding argument of -c: " << e.what() << "\n";
        return 1;
      }
      m += 4;
    } else if (arg == "-msltohae")
      heightmult = Geoid::GEOIDTOELLIPSOID;
    else if (arg == "-haetomsl")
      heightmult = Geoid::ELLIPSOIDTOGEOID;
    else if (arg == "-z") {
      if (++m == argc) return usage(1);
      zone = argv[m];
      bool northp;
      int zonenum;
      try {
        UTMUPS::DecodeZone(zone, zonenum, northp);
        zone += ' ';
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding zone: " << e.what() << "\n";
        return 1;
      }
    } else if (arg == "-n") {
      if (++m == argc) return usage(1);
      geoid = argv[m];
    } else if (arg == "-d") {
      if (++m == argc) return usage(1);
      dir = argv[m];
    } else if (arg == "-l") {
      cubic = false;
    } else if (arg == "-v")
      verbose = true;
    else
      return usage(arg != "-h");
  }

  int retval = 0;
  try {
    Geoid g(geoid, dir, cubic);
    try {
      if (cacheall)
        g.CacheAll();
      else if (cachearea)
        g.CacheArea(caches, cachew, cachen, cachee);
    }
    catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << "\nProceeding without a cache\n";
    }
    if (verbose) {
      std::cerr << "Geoid file: "    << g.GeoidFile()     << "\n"
                << "Description: "   << g.Description()   << "\n"
                << "Interpolation: " << g.Interpolation() << "\n"
                << "Date & Time: "   << g.DateTime()      << "\n"
                << "Offset (m): "    << g.Offset()        << "\n"
                << "Scale (m): "     << g.Scale()         << "\n"
                << "Max error (m): " << g.MaxError()      << "\n"
                << "RMS error (m): " << g.RMSError()      << "\n";
      if (g.Cache())
        std::cerr << "Caching:"
                  << "\n SW Corner: " << g.CacheSouth() << " " << g.CacheWest()
                  << "\n NE Corner: " << g.CacheNorth() << " " << g.CacheEast()
                  << "\n";
    }

    std::cout << std::fixed;
    GeoCoords p;
    std::string s;
    const char* spaces = " \t\n\v\f\r,"; // Include comma as space
    if (heightmult)
      std::cout << std::setprecision(3);
    while (std::getline(std::cin, s)) {
      try {
        real height = 0;
        if (heightmult) {
          std::string::size_type pb = s.find_last_not_of(spaces);
          std::string::size_type pa = s.find_last_of(spaces, pb);
          std::string::size_type px = s.find_last_not_of(spaces, pa);
          if (pa == std::string::npos || pb == std::string::npos ||
              px == std::string::npos)
            throw GeographicErr("Incomplete input: " + s);
          height = DMS::Decode(s.substr(pa + 1, pb - pa));
          s = s.substr(0, px + 1);
        }
        p.Reset(zone + s);
        real gradn, grade;
        real h = g(p.Latitude(), p.Longitude(), gradn, grade);
        if (heightmult)
          std::cout << s << " " << height + real(heightmult) * h << "\n";
        else
          std::cout << std::setprecision(4) << h << " " << std::setprecision(2)
                    << gradn * 1e6 << "e-6 " << grade * 1e6 << "e-6\n";
      }
      catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << "\n";
        retval = 1;
      }
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Error reading " << geoid << ": " << e.what() << "\n";
    retval = 1;
  }
  return retval;
}

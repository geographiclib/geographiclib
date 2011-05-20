/**
 * \file Planimeter.cpp
 * \brief Command line utility for measuring the area of geodesic polygons
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geodesic.o GeodesicLine.o DMS.o
 *
 * See \ref geod for usage information.
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/DMS.hpp"
#include "GeographicLib/GeoCoords.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"Usage: Planimeter [-e a r] [-h]\n\
$Id: Planimeter.cpp 6867 2010-09-11 13:04:26Z karney $\n\
\n\
Measure the area of a geodesic polygon.  Reads polygon vertices from\n\
standard input, one per line.  Vertices may be given as latitude and\n\
longitude, UTM/UPS, or MGRS coordinates (interpreted in the same way as\n\
GeoConvert).  (MGRS coordinates signify the center of the corresponing\n\
MGRS square.)  The end of input, a blank line, or a line which can't be\n\
interpreted as a vertex signals the end of one polygon and the start of\n\
the next.  For each polygon print a summary line with the number of\n\
points, the perimeter (in meters), and the area (in meters^2).  Areas\n\
are traversed in a clockwise sense (in other words, the included area\n\
is to the right of the perimeter).  By this rule, a polygon tranversed\n\
in a counter-clockwise sense will return the area of ellipsoid excluded\n\
by the polygon.  Only simple polygons are supported for the area\n\
computation.  Polygons may include one or both poles.\n\
\n\
By default, the WGS84 ellipsoid is used.  Specifying \"-e a r\" sets the\n\
equatorial radius of the ellipsoid to \"a\" and the reciprocal flattening\n\
to r.  Setting r = 0 results in a sphere.  Specify r < 0 for a prolate\n\
ellipsoid.\n\
\n\
-h prints this help.\n";
  return retval;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;

  class GeodesicPolygon {
  private:
    const Geodesic& _g;
    const real _area0;          // Full ellipsoid area
    unsigned _num;
    int _crossings;
    real _area, _perimeter, _lat0, _lon0, _lat1, _lon1;
    // Copied from Geodesic class
    static inline real AngNormalize(real x) throw() {
      // Place angle in [-180, 180).  Assumes x is in [-540, 540).
      return x >= 180 ? x - 360 : x < -180 ? x + 360 : x;
    }
    static inline int transit(real lon1, real lon2) {
      // Return 1 or -1 if crossing prime meridian in east or west direction.
      // Otherwise return zero.
      lon1 = AngNormalize(lon1);
      lon2 = AngNormalize(lon2);
      // treat lon12 = -180 as an eastward geodesic, so convert to 180.
      real lon12 = -AngNormalize(lon1 - lon2); // In (-180, 180]
      int cross =
        lon1 < 0 && lon2 >= 0 && lon12 > 0 ? 1 :
        lon2 < 0 && lon1 >= 0 && lon12 < 0 ? -1 : 0;
      return cross;
    }
  public:
    GeodesicPolygon(const Geodesic& g) throw()
      : _g(g)
      , _area0(_g.EllipsoidArea())
    {
      Clear();
    }
    void Clear() throw() {
      _num = 0;
      _crossings = 0;
      _area = _perimeter = 0;
      _lat0 = _lon0 = _lat1 = _lon1 = 0;
    }
    void AddPoint(real lat, real lon) throw() {
      if (_num == 0) {
        _lat0 = _lat1 = lat;
        _lon0 = _lon1 = lon;
      } else {
        real s12, S12, t;
        _g.GenInverse(_lat1, _lon1, lat, lon,
                      Geodesic::DISTANCE | Geodesic::AREA,
                      s12, t, t, t, t, t, S12);
        _perimeter += s12;
        _area += S12;
        _crossings += transit(_lon1, lon);
        if (_area > _area0/2)
          _area -= _area0;
        if (_area <= -_area0/2)
          _area += _area0;
        _lat1 = lat;
        _lon1 = lon;
      }
      ++_num;
    }
    unsigned Compute(real& perimeter, real& area) const throw() {
      real s12, S12, t;
      _g.GenInverse(_lat1, _lon1, _lat0, _lon0,
                    Geodesic::DISTANCE | Geodesic::AREA,
                    s12, t, t, t, t, t, S12);
      perimeter = _perimeter + s12;
      area = _area + S12;
      int crossings = _crossings + transit(_lon1, _lon0);
      if (crossings & 1) {
        if (area < 0)
          area += _area0/2;
        else
          area -= _area0/2;
      }
      if (area > _area0)
        area -= _area0;
      if (area < 0)
        area += _area0;
      return _num;
    }
  };

  real
    a = Constants::WGS84_a(),
    r = Constants::WGS84_r();
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-e") {
      if (m + 2 >= argc) return usage(1);
      try {
        a = DMS::Decode(std::string(argv[m + 1]));
        r = DMS::Decode(std::string(argv[m + 2]));
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
        return 1;
      }
      m += 2;
    } else
      return usage(arg != "-h");
  }

  const Geodesic geod(a, r);
  GeodesicPolygon poly(geod);
  GeoCoords p;
  std::cout << std::fixed;

  std::string s;
  real perimeter, area;
  unsigned num;
  while (std::getline(std::cin, s)) {
    try {
      p.Reset(s);
    }
    catch (const GeographicErr&) {
      num = poly.Compute(perimeter, area);
      if (num > 0)
        std::cout << num << " "
                  << std::setprecision(8) << perimeter << " "
                  << std::setprecision(3) << area << "\n";
      poly.Clear();
      continue;
    }
    poly.AddPoint(p.Latitude(), p.Longitude());
  }
  num = poly.Compute(perimeter, area);
  if (num > 0)
    std::cout << num << " "
              << std::setprecision(8) << perimeter << " "
              << std::setprecision(3) << area << "\n";
  poly.Clear();
  return 0;
}

/**
 * \file Planimeter.cpp
 * \brief Command line utility for measuring the area of geodesic polygons
 *
 * Copyright (c) Charles Karney (2010, 2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile with -I../include and link with Geodesic.o GeodesicLine.o DMS.o
 *
 * See the <a href="Planimeter.1.html">man page</a> for usage
 * information.
 *
 * $Id$
 **********************************************************************/

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/DMS.hpp"
#include "GeographicLib/GeoCoords.hpp"
#include <iostream>

#include "Planimeter.usage"

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;

  class Accumulator {
    // Compute a sum following W. M. Kahan, CACM 8(1), 40 (1965).
  private:
    real _s, _s2;
  public:
    Accumulator() throw() : _s(0), _s2(0) {};
    void Clear() throw() { _s = 0; _s2 = 0; }
    // Accumulate y
    void Add(real y) throw() {
      _s2 += y;
      volatile real t = _s + _s2;
      _s2 += _s - t;
      _s = t;
    }
    void Negate() throw() { _s *= -1; _s2 *= -1; }
    // Return sum +  y (don't accumulate)
    real Sum(real y) const throw() { return _s + (_s2 + y); }
    // Return sum
    real Sum() const throw() { return _s; }
  };

  class GeodesicPolygon {
  private:
    const Geodesic& _g;
    const real _area0;          // Full ellipsoid area
    unsigned _num;
    int _crossings;
    Accumulator _area, _perimeter;
    real _lat0, _lon0, _lat1, _lon1;
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
      _area.Clear();
      _perimeter.Clear();
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
        _perimeter.Add(s12);
        _area.Add(S12);
        if (_area.Sum() > _area0/2)
          _area.Add(-_area0);
        else if (_area.Sum() <= -_area0/2)
          _area.Add(_area0);
        _crossings += transit(_lon1, lon);
        _lat1 = lat;
        _lon1 = lon;
      }
      ++_num;
    }
    unsigned Compute(bool reverse, bool sign,
                     real& perimeter, real& area) const throw() {
      real s12, S12, t;
      if (_num < 2) {
        perimeter = area = 0;
        return _num;
      }
      _g.GenInverse(_lat1, _lon1, _lat0, _lon0,
                    Geodesic::DISTANCE | Geodesic::AREA,
                    s12, t, t, t, t, t, S12);
      perimeter = _perimeter.Sum(s12);
      Accumulator area1(_area);
      area1.Add(S12);
      if (area1.Sum() > _area0/2)
        area1.Add(-_area0);
      else if (area1.Sum() <= -_area0/2)
        area1.Add(_area0);
      int crossings = _crossings + transit(_lon1, _lon0);
      if (crossings & 1) {
        if (area1.Sum() < 0)
          area1.Sum(_area0/2);
        else
          area1.Sum(-_area0/2);
      }
      // area is with the clockwise sense.  If !reverse convert to
      // counter-clockwise convention.
      if (!reverse)
        area1.Negate();
      // If sign put area in (-area0/2, area0/2], else put area in [0, area0)
      if (sign) {
        if (area1.Sum() > _area0/2)
          area1.Add(-_area0);
        else if (area1.Sum() <= -_area0/2)
          area1.Add(_area0);
      } else {
        if (area1.Sum() >= _area0)
          area1.Add(-_area0);
        else if (area1.Sum() < 0)
          area1.Add(_area0);
      }
      area = area1.Sum();
      return _num;
    }
  };

  real
    a = Constants::WGS84_a<real>(),
    r = Constants::WGS84_r<real>();
  bool reverse = false, sign = false;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-r")
      reverse = !reverse;
    else if (arg == "-s")
      sign = !sign;
    else if (arg == "-e") {
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
    } else
      return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
  }

  const Geodesic geod(a, r);
  GeodesicPolygon poly(geod);
  GeoCoords p;

  std::string s;
  real perimeter, area;
  unsigned num;
  while (std::getline(std::cin, s)) {
    try {
      p.Reset(s);
      if (p.Latitude() != p.Latitude() || p.Longitude() != p.Longitude())
        throw GeographicErr("NAN");
    }
    catch (const GeographicErr&) {
      num = poly.Compute(reverse, sign, perimeter, area);
      if (num > 0)
        std::cout << num << " "
                  << DMS::Encode(perimeter, 8, DMS::NUMBER) << " "
                  << DMS::Encode(area, 4, DMS::NUMBER) << "\n";
      poly.Clear();
      continue;
    }
    poly.AddPoint(p.Latitude(), p.Longitude());
  }
  num = poly.Compute(reverse, sign, perimeter, area);
  if (num > 0)
        std::cout << num << " "
                  << DMS::Encode(perimeter, 8, DMS::NUMBER) << " "
                  << DMS::Encode(area, 4, DMS::NUMBER) << "\n";
  poly.Clear();
  return 0;
}

// Example of using the GeographicLib::PolygonArea class

#include <iostream>
#include <exception>
#include <GeographicLib/PolygonArea.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicExact.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/GeodesicLineExact.hpp>
#include <GeographicLib/Constants.hpp>

using namespace std;
using namespace GeographicLib;

void test(double a, double f,
          double lat0, double lon0, double azi0, double s0) {
  double ds0 = 100e3, w = 1;
  Geodesic geod(a, f);
  PolygonArea poly(geod);
  poly.AddPoint(lat0, lon0);
  poly.AddEdge(azi0, s0);
  int ndiv = max(int(floor(s0/ds0+0.5)), 1);
  double ds = s0 / ndiv;
  GeodesicLine l(geod, lat0, lon0, azi0);
  for (int i = ndiv; i >= 0; --i) {
    double azi, lat, lon;
    l.Position(i * ds, lat, lon, azi);
    azi -= 90;
    geod.Direct(lat, lon, azi, w, lat, lon);
    poly.AddPoint(lat, lon);
  }
  double perimeter, area;
  poly.Compute(false, true, perimeter, area);
  area /= s0 * w;
  if (area < 0.9 || area > 1.1)
    std::cerr << "Polygon " << area << " "
              << a << " "
              << f << " "
              << lat0 << " "
              << lon0 << " "
              << azi0 << " "
              << s0 << "\n";
}

void testExact(double a, double f,
               double lat0, double lon0, double azi0, double s0) {
  double ds0 = 100e3, w = 1;
  GeodesicExact geod(a, f);
  PolygonAreaExact poly(geod);
  poly.AddPoint(lat0, lon0);
  poly.AddEdge(azi0, s0);
  int ndiv = max(int(floor(s0/ds0+0.5)), 1);
  double ds = s0 / ndiv;
  GeodesicLineExact l(geod, lat0, lon0, azi0);
  for (int i = ndiv; i >= 0; --i) {
    double azi, lat, lon;
    l.Position(i * ds, lat, lon, azi);
    azi -= 90;
    geod.Direct(lat, lon, azi, w, lat, lon);
    poly.AddPoint(lat, lon);
  }
  double perimeter, area;
  poly.Compute(false, true, perimeter, area);
  area /= s0 * w;
  if (area < 0.9 || area > 1.1)
    std::cerr << "PolygonExact " << area << " "
              << a << " "
              << f << " "
              << lat0 << " "
              << lon0 << " "
              << azi0 << " "
              << s0 << "\n";
}

void testBoth(double a, double f,
              double lat0, double lon0, double azi0, double s0) {
  test(a, f, lat0, lon0, azi0, s0);
  testExact(a, f, lat0, lon0, azi0, s0);
}

int main() {
  try {
    {
      int azi0 = 0;
      int s0 = 0;
      int lon0 = 0;
      //      for (int azi0 = -180; azi0 <= 180; azi0 += 10) {
      for (s0 = 10; s0 <= 100; s0 += 10) {
        for (lon0 = -400; lon0 <= 400; lon0 += 10) {
          testBoth(20e6/Math::pi(), 0.2, 30, lon0, azi0-1e-10, s0*1e6);
        }
      }
      return 0;
      GeodesicExact geod(20e6/Math::pi(), 0);
      PolygonAreaExact poly(geod);
      int l0 = 136;
      poly.AddPoint(0,-135+l0);
      poly.AddEdge(90,30e6);
      for (int i = 135; i >= -135; --i)
        poly.AddPoint(90/10e6, i+l0);
      double perimeter, area;
      unsigned n = poly.Compute(false, true, perimeter, area);
      cout << n << " " << perimeter << " " << area << "\n";
    }
    Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
    // Alternatively: const Geodesic& geod = Geodesic::WGS84();
    PolygonArea poly(geod);
    poly.AddPoint( 52,  0);     // London
    poly.AddPoint( 41,-74);     // New York
    poly.AddPoint(-23,-43);     // Rio de Janeiro
    poly.AddPoint(-26, 28);     // Johannesburg
    double perimeter, area;
    unsigned n = poly.Compute(false, true, perimeter, area);
    cout << n << " " << perimeter << " " << area << "\n";
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

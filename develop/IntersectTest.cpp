#include <iostream>
#include <iomanip>
#include <string>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    Utility::set_digits();
    bool exact = false, arcdist = false;
    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-E") {
        exact = true;
      } else if (arg == "-a")
        arcdist = true;
      else {
        int retval = arg == "-h" || arg == "--help" ? 0 : 1;
        (retval ? std::cerr : std::cout) << "Usage: "
                                         << argv[0] << " [-E] [-a]\n";
        return retval;
      }
    }
    int num = 6;
    real a = Constants::WGS84_a(), f = Constants::WGS84_f();
    Geodesic g(a, f, exact);
    // Intersection at
    //   -1.44147956008236583 27.97257917717199337
    // all digits correct
    GeodesicLine lx(g, 4, 20, -56);
    GeodesicLine ly(g, -30, -40, 80);
    if (arcdist) {
      std::cout << "Arc distance calc "
                << (exact ? "exact" : "series") << "\n";
      std::cout << std::setprecision(3) << std::scientific;
      real ax = 0, ay = 0;
      for (int i = 0; i < num; ++i) {
        real latx, lonx, azix, laty, lony, aziy, zeta, gamx, gamy;
        lx.ArcPosition(ax, latx, lonx, azix);
        ly.ArcPosition(ay, laty, lony, aziy);
        zeta = g.Inverse(latx, lonx, laty, lony, gamx, gamy);
        real mux = Math::AngDiff(azix, gamx), muy = Math::AngDiff(aziy, gamy);
        if (Math::AngDiff(mux, muy) < 0) {
          mux = -mux; muy = -muy;
        }
        real smux = Math::sind(mux), cmux = Math::cosd(mux),
          smuy = Math::sind(muy), cmuy = Math::cosd(muy),
          szeta = Math::sind(zeta), czeta = Math::cosd(zeta);
        real delx = Math::atan2d(smuy * szeta,
                                 smuy * cmux * czeta - cmuy * smux),
          dely = Math::atan2d(smux * szeta,
                              -smux * cmuy * czeta + cmux * smuy);
        std::cout << i << " "
                  << delx * Math::degree() << " "
                  << dely * Math::degree() << "\n";
        ax += delx; ay += dely;
      }
      std::cout << std::setprecision(17) << std::fixed;
      real lat, lon;
      lx.ArcPosition(ax, lat, lon);
      std::cout << lat << " " << lon << "\n";
    } else {
      std::cout << "Distance calc "
                << (exact ? "exact" : "series") << "\n";
      std::cout << std::setprecision(3) << std::scientific;
      using std::atanh, std::sin, std::cos, std::atan2, std::sqrt;
      real ax = 0, ay = 0;
      real b = a * (1 - f), e = sqrt(f * (2 - f)),
        R = sqrt(a*a/2 + b*b/2 * atanh(e) / e);
      for (int i = 0; i < num; ++i) {
        real latx, lonx, azix, laty, lony, aziy, zeta, z, gamx, gamy;
        lx.Position(ax, latx, lonx, azix);
        ly.Position(ay, laty, lony, aziy);
        g.Inverse(latx, lonx, laty, lony, z, gamx, gamy);
        zeta = z/R;
        real mux = Math::AngDiff(azix, gamx), muy = Math::AngDiff(aziy, gamy);
        if (Math::AngDiff(mux, muy) < 0) {
          mux = -mux; muy = -muy;
        }
        real smux = Math::sind(mux), cmux = Math::cosd(mux),
          smuy = Math::sind(muy), cmuy = Math::cosd(muy),
          szeta = sin(zeta), czeta = cos(zeta);
        real delx = atan2(smuy * szeta,
                          smuy * cmux * czeta - cmuy * smux),
          dely = atan2(smux * szeta,
                       -smux * cmuy * czeta + cmux * smuy);
        std::cout << i << " " << delx << " " << dely << "\n";
        ax += delx*R; ay += dely*R;
      }
      std::cout << std::setprecision(17) << std::fixed;
      real lat, lon;
      lx.Position(ax, lat, lon);
      std::cout << lat << " " << lon << "\n";
    }
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

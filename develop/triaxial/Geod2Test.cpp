/**
 * \file GeodSolve2.cpp
 **********************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Geodesic.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"
#include "TriaxialODE.hpp"

int main(int argc, const char* const argv[]) {
  try {
    using namespace std;
    using namespace GeographicLib;
    using std::fmax; using std::pow; using std::sqrt;
    typedef Math::real real;
    typedef Angle ang;
    Utility::set_digits();
    bool skip = false;
    bool ode = false, extended = false, dense = false, normp = false,
      alt = false;
    real lg2eps = 20;
    for (int m = 1; m < argc; ++m) {
      string arg(argv[m]);
      if (arg == "--skip")
        skip = true;
      else if (arg == "--alt")
        alt = true;
      else if (arg == "--ode")
        ode = true;
      else if (arg == "--dense")
        dense = true;
      else if (arg == "--normp")
        normp = true;
      else if (arg == "-p") {
        if (++m == argc) {
          std::cerr << "epspow missing\n";
          return 1;
        }
        try {
          lg2eps = Utility::fract<real>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "epspow " << argv[m] << " is not a number\n";
          return 1;
        }
      } else {
        std::cerr << "Bad input\n";
        return 1;
      }
    }
    real eps = pow(1/real(2), lg2eps);
    real f = alt ? 0 : Constants::WGS84_f(),
      a = alt ? sqrt(real(2)) : Constants::WGS84_a(),
      b = alt ? 1 : a,
      c = alt ? sqrt(1/real(2)) : (1 - f) * a,
      lat1, lon1, azi1, lat2, lon2, azi2, s12,
      angnorm = 1e-9 * Math::degree() / 3600,
      distnorm = alt ? 1e-9 * a / 6.4e6 : 1.e-9;
    Triaxial::vec3 r1, v1, r2, v2, r2a, v2a;
    Triaxial t(a, b, c);
    Geodesic g(a, f);
    real err1a = 0, err2a = 0, err1b = 0, err2b = 0,
      nstepa = 0, nstepb = 0;
    int cnt = 0;
    while (cin >> lat1 >> lon1 >> azi1 >> lat2 >> lon2 >> azi2 >> s12) {
      ++cnt;
      ang phi1(lat1), bet1((1 - f) * phi1.s(), phi1.c()),
        phi2(lat2), bet2((1 - f) * phi2.s(), phi2.c()),
        omg1(lon1), alp1(azi1), omg2(lon2), alp2(azi2);
      t.elliptocart2(bet1, omg1, alp1, r1, v1);
      t.elliptocart2(bet2, omg2, alp2, r2, v2);
      real lat2a = 0, lon2a = 0, azi2a = 0;
      if (!skip) {
        if (ode) {
          TriaxialODE ts(t, r1, v1, extended, dense, normp, eps);
          (void) ts.Position(s12, r2a, v2a);
          nstepa += ts.NSteps();
          nstepb = fmax(nstepb, real(ts.NSteps()));
        } else {
          if (alt) {
            ang bet2a, omg2a, alp2a;
            t.Direct(bet1, omg1, alp1, s12, bet2a, omg2a, alp2a);
            t.elliptocart2(bet2a, omg2a, alp2a, r2a, v2a);
          } else {
            g.Direct(lat1, lon1, azi1, s12, lat2a, lon2a, azi2a);
            ang phi2a(lat2a), bet2a((1 - f) * phi2a.s(), phi2a.c()),
              omg2a(lon2a), alp2a(azi2a);
            t.elliptocart2(bet2a, omg2a, alp2a, r2a, v2a);
          }
        }
      }
      pair<real, real> d = Triaxial::EuclideanDiff(r2, v2, r2a, v2a);
      err1a += d.first;
      err2a += d.second;
      err1b = fmax(err1b, d.first);
      err2b = fmax(err2b, d.second);
    }
    err1a /= cnt; err2a /= cnt; nstepa /= cnt;
    err1a /= distnorm; err1b /= distnorm;
    err2a /= angnorm; err2b /= angnorm;
    cout << lg2eps << " "
         << err1a << " " << err2a << " "
         << err1b << " " << err2b << " "
         << nstepa << " " << int(nstepb) << "\n";
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

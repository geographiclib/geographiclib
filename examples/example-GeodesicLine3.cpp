// Example of using the Triaxial::GeodesicLine3 class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Triaxial/Geodesic3.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    Triaxial::Geodesic3 g(Triaxial::Ellipsoid3::Earth());
    Triaxial::GeodesicLine3 l = g.Line(90, 0, 135);
    // The approximate distance between opposite umbilical points
    double s0 = 20003986;
    int num = 20;
    for (int i = 0; i <= num; ++i) {
      double bet2, omg2, alp2, s12 = i*s0/num;
      l.Position(s12, bet2, omg2, alp2);
      cout << s12 << " " << bet2 << " " << omg2 << " " << alp2 << "\n";
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

// Example of using the Triaxial::Conformal3 class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Triaxial/Conformal3.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    Triaxial::Conformal3 proj(Triaxial::Ellipsoid3::Earth());
    cout  << fixed << setprecision(1)
          << "Ellipsoid parameters: a = " << proj.t().a()
          << ", b = " << proj.t().b() << ", c = " << proj.t().c() << "\n"
          << setprecision(3)
          << "Quadrants: x = " << proj.x0() << ", y = " << proj.y0() << "\n";
    cout << "Coordinates angle (deg) x (m) y (m):\n";
    cout << setprecision(2);
    for (int i = -180; i <= 180; i += 15) {
      double omg = i + 90, bet = i,
        x = proj.x(Angle(omg)), y = proj.y(Angle(bet));
      cout << i << " " << x << " " << y << "\n";
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

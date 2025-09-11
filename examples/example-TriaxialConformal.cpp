// Example of using the GeographicLib::JacobiConformal class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/TriaxialConformal.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    TriaxialConformal proj(Triaxial::Earth());
    cout  << fixed << setprecision(1)
          << "Ellipsoid parameters: a = " << proj.t().a()
          << ", b = " << proj.t().b() << ", c = " << proj.t().c() << "\n"
          << setprecision(3)
          << "Quadrants: x = " << proj.x() << ", y = " << proj.y() << "\n";
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

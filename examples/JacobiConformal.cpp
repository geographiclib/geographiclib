#include <iostream>
#include <iomanip>
#include <GeographicLib/Utility.hpp>
#include "JacobiConformal.hpp"

using namespace std;
using namespace GeographicLib;

int main() {
  Utility::set_digits();
  Math::real a = 6378137+35, b = 6378137-35, c = 6356752;
  JacobiConformal jc(a, b, c);
  cout  << fixed << setprecision(1)
        << "Ellipsoid parameters: a = "
        << a << ", b = " << b << ", c = " << c << "\n"
        << setprecision(10)
        << "Quadrants: x = " << jc.x() << ", y = " << jc.y() << "\n"
        << "Quadrants: x = " << jc.x_0() << ", y = " << jc.y_0() << "\n"
        << "Quadrants: x = " << jc.x_b() << ", y = " << jc.y_b() << "\n"
        << "Quadrants: x = " << jc.x_bx() << ", y = " << jc.y_bx() << "\n"
        << "Quadrants: x = " << jc.x_c() << ", y = " << jc.y_c() << "\n"
        << "Quadrants: x = " << jc.x_bx() << ", y = " << jc.y_bx() << "\n"
        << "Quadrants: x = " << jc.x_cx() << ", y = " << jc.y_cx() << "\n";
  Math::real f = 1/Math::degree();
  cout << "Scaled coordinates (angle x y):\n";
  for (int i = 0; i <= 90; i += 5) {
    Math::real omg = i, bet = i;
    cout << i << " " << (jc.x()-jc.x(90-omg))*f << " " << jc.y(bet)*f << "\n";
    cout << i << " " << jc.x_0(omg)*f << " " << jc.y_0(bet)*f << "\n";
    cout << i << " " << jc.x_b(omg)*f << " " << jc.y_b(bet)*f << "\n";
    cout << i << " " << jc.x_c(omg)*f << " " << jc.y_c(bet)*f << "\n";
  }
}

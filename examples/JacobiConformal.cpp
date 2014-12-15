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
  cout  << fixed << setprecision(10)
        << "Quadrants " << jc.x() << " " << jc.y() << "\n";
  Math::real q = Math::pi() * sqrt((a*a + b*b) / (a*a + b*b - 2*c*c));
  for (int i = 0; i <= 90; i += 1) {
    if (i == -180) continue;
    Math::real omg = i, bet = i;
    cout << i << " " << jc.x(omg)*90/q << " " << jc.y(bet)*90/q << "\n";
  }
}

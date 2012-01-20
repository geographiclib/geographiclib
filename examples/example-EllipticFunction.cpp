// Example of using GeographicLib::EllipticFunction class
// $Id$

#include <iostream>
#include <cmath>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/EllipticFunction.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    EllipticFunction ell(0.1);  // parameter m = 0.1
    // See Abramowitz and Stegun, table 17.1
    cout << ell.K() << " " << ell.E() << "\n";
    double phi = 20 * Math::degree();
    // See Abramowitz and Stegun, table 17.6 with
    // alpha = asin(sqrt(m)) = 18.43 deg and phi = 20 deg
    cout << ell.E(phi) << " "
         << ell.E(sin(phi), cos(phi), sqrt(1 - ell.m() * Math::sq(sin(phi))))
         << "\n";
  }
  catch (const GeographicErr& e) {
    cout << "Caught exception: " << e.what() << "\n";
  }
  return 0;
}

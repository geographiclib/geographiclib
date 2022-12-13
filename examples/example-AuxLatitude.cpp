// Example of using the GeographicLib::AuxLatitude class.  See the paper
//
// - C. F. F. Karney,
//   On auxiliary latitudes,
//   Technical Report, SRI International, December 2022.
//   https://arxiv.org/abs/2212.05818

#include <iostream>
#include <iomanip>
#include <exception>
#include "AuxLatitude.hpp"

using namespace std;

int main() {
  try {
    typedef GeographicLib::AuxLatitude<double> latitude;
    typedef latitude::angle angle;
    double a = 2, b = 1;        // Equatorial radius and polar semi-axis
    latitude aux(a, b);
    bool series = false;        // Don't use series method
    int auxin = latitude::GEOGRAPHIC;
    cout << setprecision(9) << fixed;
    for (int l = 0; l <= 90; ++l) {
      angle phi(angle::degrees(l));
      for (int auxout = 0; auxout < latitude::AUXNUMBER; ++auxout) {
        angle eta = aux.Convert(auxin, auxout, phi, series);
        cout << (auxout ? " " : "") << eta.degrees();
      }
      cout << "\n";
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

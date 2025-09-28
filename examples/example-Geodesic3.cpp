// Example of using the Triaxial::Geodesic3 class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Triaxial/Geodesic3.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    Triaxial::Geodesic3 g(Triaxial::Ellipsoid3::Earth());
    double bet1 = 10, omg1 = 30, bet2 = 80, omg2 = -50,
      alp1, alp2, s12;
    g.Inverse(bet1, omg1, bet2, omg2, s12, alp1, alp2);
    cout << "geodesic between [" << bet1 << "," << omg1 << "] and ["
         << bet2 << "," << omg2 << "]\n";
    cout << fixed << setprecision(3)
         << "s12 = " << s12 << ", alp1 = "
         << alp1 << ", alp2 = " << alp2 << "\n";
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

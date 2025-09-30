// Example of using the Triaxial::Cartesian3 class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <random>
#include <GeographicLib/Triaxial/Cartesian3.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    using Triaxial::Cartesian3;
    Cartesian3 cart(3, 2, 1);
    unsigned long long
      s1 = std::random_device()(),
      s2 = std::random_device()();
    std::seed_seq seq{s1, s2};
    std::mt19937 g(seq);
    cout << "10 random points on the ellipsoid with semiaxes [3, 2, 1]\n"
         << "X Y Z beta omega\n" << fixed << setprecision(3);
    Cartesian3::vec3 r;
    Angle bet, omg;
    for (int i = 0; i < 10; ++i) {
      cart.cart2rand(g, r);
      cart.cart2toany(r, Cartesian3::ELLIPSOIDAL, bet, omg);
      cout << r[0] << " " << r[1] << " " << r[2] << " "
           << double(bet) << " " << double(omg) << "\n";
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

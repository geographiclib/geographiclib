// Example of using the Triaxial::Ellipsoid3 class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    using Triaxial::Ellipsoid3;
    using ang = Angle;
    using vec3 = Ellipsoid3::vec3;
    Ellipsoid3 t(3, 2, 1);
    ang bet{40}, omg{20}, alp{80};
    vec3 r, v;
    t.elliptocart2(bet, omg, alp, r, v);
    cout  << fixed << setprecision(3)
          << "[bet, omg, alp] = ["
          << double(bet) << " " << double(omg) << " " << double(alp)
          << "]\n => r = [" << r[0] << " " << r[1] << " " << r[2]
          << "], v = [" << v[0] << " " << v[1] << " " << v[2] << "]\n";
    ang betn, omgn, alpn;
    t.cart2toellip(r, v, betn, omgn, alpn);
    cout << " => [bet, omg, alp] = ["
         << double(betn) << " " << double(omgn) << " " << double(alpn) << "]\n";
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

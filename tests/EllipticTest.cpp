#include <iostream>
#include <iomanip>
#include <GeographicLib/EllipticFunction.hpp>

using namespace GeographicLib;
int main() {
  typedef GeographicLib::Math::real real;
  try {
    EllipticFunction ell(0.99, 0.8);
    real dphi = Math::pi<real>()/(2*100);
    dphi = 0.05;
    std::cout << std::fixed << std::setprecision(10);
    for (int i = -100; i <= 100; ++i) {
      real phi = i * dphi;
      std::cout << phi << " "
                << ell.F(phi) << " "
                << ell.E(phi) << " "
                << ell.D(phi) << " "
                << ell.Pi(phi) << " "
                << ell.G(phi) << "\n";
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
  return 0;
}

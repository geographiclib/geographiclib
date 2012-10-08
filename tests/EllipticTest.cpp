#include <iostream>
#include <iomanip>
#include <GeographicLib/EllipticFunction.hpp>

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions
#  pragma warning (disable: 4127)
#endif

using namespace GeographicLib;
int main() {
  typedef GeographicLib::Math::real real;
  try {
    // For tabulated values in A+S
    real
      ASalpha = 30*Math::degree<double>(),
      k2 = Math::sq(std::sin(ASalpha)),
      alpha2 = 0.3;

    EllipticFunction ell(k2, alpha2);
    real dphi = Math::degree<real>();
    std::cout << std::fixed << std::setprecision(10);
    for (int i = 0; i <= 90; i += 15) {
      real phi = i * dphi;
      std::cout << i << " "
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

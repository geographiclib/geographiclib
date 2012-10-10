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
    {
      real alpha2 = 0.8, k2 = -0.4;
      EllipticFunction ellG(k2,alpha2);
      EllipticFunction ellH(k2,k2/alpha2);
      
      std::cout << std::setprecision(10);
      for (int i = -179; i <= 180; i += 10) {
        real
          phi = i * Math::degree<real>(),
          sn = sin(phi), cn = cos(phi), dn = ellG.Delta(sn, cn),
          g = ellG.G(phi),
          h = (k2/alpha2)*ellH.H(phi) + sqrt(1-k2/alpha2)/sqrt(1-alpha2)*
          atan2(sqrt(1-alpha2)*sqrt(1-k2/alpha2)*sn, dn*cn);
        
        std::cout << i << " " << g << " " << h << " " << h-g << "\n";
      }
      return 0;
    }
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
                << ell.G(phi) << " "
                << ell.H(phi) << "\n";
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

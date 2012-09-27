#include <iostream>
#include <iomanip>
#include <GeographicLib/EllipticFunction.hpp>

using namespace GeographicLib;
int main() {
  typedef GeographicLib::Math::real real;
  try {
    EllipticFunction ell(-99.0);
    real dphi = Math::pi<real>()/(2*100);
    std::cout << std::fixed << std::setprecision(16);
    for (int i = 0; i <= 100; ++i) {
      real
        phi = i ? i * dphi : 0.000001,
        e1 = ell.E(phi),
        c = 1/Math::sq(sin(phi)),
        e2 = EllipticFunction::RF(c - 1, c - ell.m(), c) -
        (ell.m()/3) * EllipticFunction::RD(c - 1, c - ell.m(), c),
        e3 = ell.m1() * EllipticFunction::RF(c - 1, c - ell.m(), c) +
        (ell.m() * ell.m1()/3) * EllipticFunction::RD(c - 1, c, c - ell.m()) +
        ell.m() * std::sqrt( (c - 1)/(c * (c - ell.m())) ),
        e4 = - (ell.m1()/3) * EllipticFunction::RD(c - ell.m(), c, c - 1) + 
        std::sqrt( (c - ell.m())/(c * (c - 1)) );
      std::cout << phi << " "
                << e1 << " "
                << e2 << " "
                << e3 << " "
                << e4 << "\n";
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

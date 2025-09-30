// Example of using the GeographicLib::TrigfunExt class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Trigfun.hpp>
#include <GeographicLib/EllipticFunction.hpp>

using namespace std;
using namespace GeographicLib;

int main(int argc, const char* const argv[]) {
  try {
    // Integrate 1/sqrt(1 - k2*sin(x)^2) and compare this to the elliptic
    // integral of the first kine.
    double k2 = Math::sq(0.8);
    auto f = [k2] (double x) -> double
    { return 1/sqrt(1 - k2 * Math::sq(sin(x))); };
    TrigfunExt tfe(f, Math::pi()/2);
    cout << "Number of coefficients: " << tfe.NCoeffs() << "\n";
    EllipticFunction ell(k2);
    cout << "x tf(x) tf(x)-F(x)\n";
    for (int i = 0; i <= 20; ++i) {
      double x = i/10.0;
      cout << x << " " << tfe(x) << " " << tfe(x) - ell.F(x) << "\n";
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

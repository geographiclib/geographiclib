// Example of using the GeographicLib::Trigfun class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Trigfun.hpp>

using namespace std;
using namespace GeographicLib;

int main(int argc, const char* const argv[]) {
  try {
    // Approximate 1/sqrt(1 - k2*sin(x)^2) by a Fourier series
    double k2 = Math::sq(0.8);
    auto f = [k2] (double x) -> double
    { return 1/sqrt(1 - k2 * Math::sq(sin(x))); };
    Trigfun tf(f, false, false, Math::pi()/2);
    cout << "Number of coefficients: " << tf.NCoeffs() << "\n";
    cout << "x tf(x) tf(x)-f(x)\n";
    for (int i = 0; i <= 20; ++i) {
      double x = i/10.0;
      cout << x << " " << tf(x) << " " << tf(x) - f(x) << "\n";
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

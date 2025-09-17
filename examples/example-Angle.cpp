// Example of using the GeographicLib::Angle class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Angle.hpp>

int main(int argc, const char* const argv[]) {
  try {
    using ang = GeographicLib::Angle;
    // Print table of parametric latitudes for f = 0.5
    double f = 0.5;
    std::cout << std::fixed << std::setprecision(4);
    for (double d = 0; d <= 90; d+=10) {
      ang phi{d};
      std::cout << double(d) << " " << double(phi.modang(1-f)) << "\n";
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

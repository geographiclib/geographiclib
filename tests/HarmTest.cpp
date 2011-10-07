#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <GeographicLib/SphericalHarmonic.hpp>

using namespace GeographicLib;
int main() {
  try {
    int type = 0;
    int
      N = type ? 2190 : 360,
      k = ((N + 1) * (N + 2)) / 2;
    std::vector<double> C(k);
    std::vector<double> S(k);
    {
      std::ifstream f(type ?
                      "/scratch/egm2008/harm/data_EGM2008_2190.dat" :
                      "/scratch/egm2008/harm/data_EGM96_360.dat",
                      std::ios::binary);
      if (!f.good())
        throw GeographicErr("Cannot open coefficient file");
      f.read(reinterpret_cast<char *>(&C[0]), k * sizeof(double));
      f.read(reinterpret_cast<char *>(&S[0]), k * sizeof(double));
    }
    //    for (int i = 0; i < k; ++i)
    //      std::cout << i << " " << C[i] << " " << S[i] << "\n";
    double lat, lon;
    std::cout << std::setprecision(17);
    double a(1.0L), r(1.2L);
    while (std::cin >> lat >> lon) {
      double
        phi = Math::degree<double>() * lat,
        lam = Math::degree<double>() * lon,
        z = r * sin(phi),
        x = r * cos(phi) * cos(lam),
        y = r * cos(phi) * sin(lam);
      std::cout << SphericalHarmonic::Value(N, C, S, x, y, z, a)
                << "\n";
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

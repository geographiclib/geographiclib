#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <GeographicLib/SphericalHarmonic.hpp>

typedef GeographicLib::SphericalHarmonic::work real;

using namespace GeographicLib;
int main() {
  try {
    int type = 2;
    int N;
    std::string name;
    switch (type) {
    case 0:
      N = 360;
      name = "data_EGM96_360.dat";
      break;
    case 1:
      N = 2190;
      name = "data_EGM2008_2190.dat";
      break;
    case 2:
    default:
      N = 5;
      name = "harmtest.dat";
      break;
    }
    int k = ((N + 1) * (N + 2)) / 2;
    std::vector<double> C(k);
    std::vector<double> S(k);
    name = "/scratch/egm2008/harm/" + name;
    {
      std::ifstream f(name.c_str(), std::ios::binary);
      if (!f.good())
        throw GeographicErr("Cannot open coefficient file");
      f.read(reinterpret_cast<char *>(&C[0]), k * sizeof(double));
      f.read(reinterpret_cast<char *>(&S[0]), k * sizeof(double));
    }
    //    for (int i = 0; i < k; ++i)
    //      std::cout << i << " " << C[i] << " " << S[i] << "\n";
    real lat, lon;
    std::cout << std::setprecision(17);
    real a(1.0L), r(1.201L);
    while (std::cin >> lat >> lon) {
      real
        phi = Math::degree<real>() * lat,
        lam = Math::degree<real>() * lon,
        x = r * cos(phi) * cos(lam),
        y = r * cos(phi) * sin(lam),
        z = r * sin(phi);
      std::cout << SphericalHarmonic::Value(N, C, S, x, y, z, a)
                << "\n";
      real
        d = 1e-6,
        dr1 = (SphericalHarmonic::Value(N, C, S,
                                        (r+d) * cos(phi) * cos(lam),
                                        (r+d) * cos(phi) * sin(lam),
                                        (r+d) * sin(phi),
                                        a) -
               SphericalHarmonic::Value(N, C, S,
                                        (r-d) * cos(phi) * cos(lam),
                                        (r-d) * cos(phi) * sin(lam),
                                        (r-d) * sin(phi),
                                        a)) / (2 * d),
        dr2;
      real
        dl1 = (SphericalHarmonic::Value(N, C, S,
                                        r * cos(phi) * cos(lam+d),
                                        r * cos(phi) * sin(lam+d),
                                        r * sin(phi),
                                        a) -
               SphericalHarmonic::Value(N, C, S,
                                        r * cos(phi) * cos(lam-d),
                                        r * cos(phi) * sin(lam-d),
                                        r * sin(phi),
                                        a)) / (2 * d),
        dl2;
      real
        dt1 = (SphericalHarmonic::Value(N, C, S,
                                        r * cos(phi+d) * cos(lam),
                                        r * cos(phi+d) * sin(lam),
                                        r * sin(phi+d),
                                        a) -         
               SphericalHarmonic::Value(N, C, S,     
                                        r * cos(phi-d) * cos(lam),
                                        r * cos(phi-d) * sin(lam),
                                        r * sin(phi-d),
                                        a)) / (- 2 * d),
        dt2;
      SphericalHarmonic::Value(N, C, S, x, y, z, a, dr2, dl2, dt2);
      std::cout << dr1 << " " << dr2 << "\n"
                << dl1 << " " << dl2 << "\n"
                << dt1 << " " << dt2 << "\n";
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

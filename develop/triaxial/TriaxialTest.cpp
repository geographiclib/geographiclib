#include <iostream>
#include <limits>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include "Triaxial.hpp"

using namespace GeographicLib;
using namespace std;

int main() {
  try {
    Utility::set_digits();
    using std::sqrt;
    using std::pow;
    typedef Math::real real;
    typedef Triaxial::vec3 vec3;
    real a = 6378172, b = 6378103, c = 6356753;
    Triaxial t(a, b, c);
    vec3 r1{t.axes[0] * t.kp, 0, t.axes[2] * t.k},
      v1{0,1,0};
    real s12 = 4 * EllipticFunction::RG(Math::sq(t.axes[0]),
                                        Math::sq(t.axes[2]));
    vec3 r2, v2;
    int kmax = Math::digits()-4;
    if (0) {
      for (int k = 10; k <= kmax; ++k) {
        real eps = pow(real(2), -k);
        int n = t.Direct(r1, v1, s12, r2, v2, eps);
        cout << k << " " << n << " "
             << sqrt(Math::sq(r1[0]+r2[0]) +
                     Math::sq(r1[1]+r2[1]) +
                     Math::sq(r1[2]+r2[2])) / t.b << endl;
      }
    } else {
      int n = 0, imax = 1000;
#if GEOGRAPHICLIB_PRECISION >= 4
      real eps = pow(numeric_limits<real>::epsilon(), real(0.875));
#else
      real eps = numeric_limits<real>::epsilon() * 16;
#endif
      for (int i = 0; i < imax; ++i) {
        n += t.Direct(r1, v1, s12, r2, v2, eps);
        r1 = r2;
        v1 = v2;
      }
      cout << n/real(imax) << "\n";
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
}

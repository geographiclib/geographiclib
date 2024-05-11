#include <iostream>
#include <limits>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/AuxAngle.hpp>
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"

using namespace GeographicLib;
using namespace std;

void ODEtest(Math::real a, Math::real b, Math::real c) {
    typedef Math::real real;
    typedef Triaxial::vec3 vec3;
    Triaxial t(a, b, c);
    vec3 r1{t.axes[0] * t.kp, 0, t.axes[2] * t.k},
      v1{0,1,0};
    real s12 = 4 * EllipticFunction::RG(Math::sq(t.axes[0]),
                                        Math::sq(t.axes[2]));
    vec3 r2, v2;
    int kmax = Math::digits()-4;
    if (1) {
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

/*
void TriaxialTest0() {
  typedef Math::real real;
  // Triaxial t(sqrt(real(2)), 1, sqrt(real(0.5)));
  Triaxial t(1.01,1,0.8);
  // Triaxial t(sqrt(real(3)), 1, 1/sqrt(real(3)));
  real mu = real(0.0000001);
  //  cout << t.k2 << " " << t.kp2 << " " << t.e2 << "\n";
  geod_fun fa(t.k2, t.kp2, -t.e2, mu, true);
  fa.NCoeffsInv();
  cerr << fa.NCoeffs() << " "
            << fa.NCoeffsInv() << " "
            << fa.InvCounts().first << "\n";
  geod_fun fb(t.k2, t.kp2, -t.e2, mu, false);
  fb.NCoeffsInv();
  cerr << fb.NCoeffs() << " "
            << fb.NCoeffsInv() << " "
            << fb.InvCounts().first << "\n";
  for (int k = 0; k <= 360; k += 3) {
    real x = k * Math::degree(),
      u = fa.ell.F(x),
      fu = fa.fun(u),
      fx = fb.fun(x),
      uu = fa.fun.inv(fu),
      xx = fb.fun.inv(fx);
      cout << k << " " << fu << " " << fx << " " << fu - fx << " "
                << u - uu << " " << x - xx << "\n";
  }
}
*/
void TriaxialTest1(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  Triaxial t(a, b, c);
  real k2 = t.k2, kp2 = t.kp2, e2 = t.e2;
  if (0) {
    int num = 100, numk = int(round(num*k2)), numkp = num-numk;
    for (int k = -numkp; k <= numk; ++k) {
      if (k == 0) continue;
      real gam = k < 0 ?
        - kp2 * Math::sq(k / real(numkp)) : k2 * Math::sq(k / real(numk));
      if (0) {
        geod_fun fphia(k2, kp2, e2, -gam, false); fphia.NCoeffsInv();
        geod_fun fphib(k2, kp2, e2, -gam, true); fphib.NCoeffsInv();
        geod_fun fomga(kp2, k2, -e2, gam, false); fomga.NCoeffsInv();
        geod_fun fomgb(kp2, k2, -e2, gam, true); fomgb.NCoeffsInv();
         cout << k << " " << gam << " "
                   << fphia.NCoeffs() << " " << fphia.NCoeffsInv() << " "
                   << fomga.NCoeffs() << " " << fomga.NCoeffsInv() << " "
                   << fphib.NCoeffs() << " " << fphib.NCoeffsInv() << " "
                   << fomgb.NCoeffs() << " " << fomgb.NCoeffsInv() << "\n";
      } else {
        geod_fun fphia(k2, kp2, e2, -gam); fphia.NCoeffsInv();
        geod_fun fomga(kp2, k2, -e2, gam); fomga.NCoeffsInv();
         cout << k << " " << gam << " "
                   << fphia.txp() << " "
                   << fphia.NCoeffs() << " " << fphia.NCoeffsInv() << " "
                   << fomga.txp() << " "
                   << fomga.NCoeffs() << " " << fomga.NCoeffsInv() << "\n";
      }
    }
  }
  for (int k = 3; k <= 16; ++k) {
    for (int s = -1; s <= 1; s += 2) {
      real gam = pow(real(10), -k) * s;
      geod_fun fphia(k2, kp2, e2, -gam); fphia.NCoeffsInv();
      geod_fun fomga(kp2, k2, -e2, gam); fomga.NCoeffsInv();
       cout << k << " " << gam << " "
                 << fphia.HalfPeriod() << " "
                 << fphia.NCoeffs() << " " << fphia.NCoeffsInv() << " "
                 << fomga.HalfPeriod() << " "
                 << fomga.NCoeffs() << " " << fomga.NCoeffsInv() << "\n";
    }
  }
}

void UmbCheck(const geodu_fun& f) {
  typedef Math::real real;
  int num = 25; real maxval = 5, dx = maxval/num;
  for (int i = -num; i <= num; ++i) {
    real x = i * dx;
    cout << "c " << x << " " << x - f.inv(f(x)) << "\n";
  }
}

void UmbilicTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  Triaxial t(a, b, c);
  real k2 = t.k2, kp2 = t.kp2, e2 = t.e2;
  geodu_fun fbet(k2, kp2, e2);
  cout << fbet.NCoeffs() << " " << fbet.NCoeffsInv() << " "
       << fbet.InvCounts().first << " " << fbet.InvCounts().second << "\n";
  geodu_fun fomg(kp2, k2, -e2);
  cout << fomg.NCoeffs() << " " << fomg.NCoeffsInv() << " "
       << fomg.InvCounts().first << " " << fomg.InvCounts().second << "\n";
  cout << "\n";
  UmbCheck(fbet);
  cout << "\n";
  UmbCheck(fomg);
}

void DirectfunTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  Triaxial t(a, b, c);
  // Circumpolar
  TriaxialLine lca(t, AuxAngle::degrees(real(1)), AuxAngle::degrees(real(0)),
                  AuxAngle::degrees(real(90)));
  if (1) {
  TriaxialLine lcb(t, AuxAngle::degrees(real(89)), AuxAngle::degrees(real(0)),
                  AuxAngle::degrees(real(90)));
  }
  // Umbilic
  TriaxialLine lu(t, AuxAngle::degrees(real(90)), AuxAngle::degrees(real(0)),
                 AuxAngle::degrees(real(135)));
  if (1) {
  // Circumpolar
  TriaxialLine ltb(t, AuxAngle::degrees(real(90)), AuxAngle::degrees(real(1)),
                  AuxAngle::degrees(real(180)));
  }
  TriaxialLine lta(t, AuxAngle::degrees(real(90)), AuxAngle::degrees(real(89)),
                  AuxAngle::degrees(real(180)));
}

int main() {
  try {
    Utility::set_digits();
    typedef Math::real real;
    real a = 6378172, b = 6378103, c = 6356753;
    if (0)
      ODEtest(a, b, c);
    /*
    if (0)
      TriaxialTest0();
    */
    if (0) {
      real a = 1.01, b = 1, c = 0.8;
      //    a = sqrt(2.0); c = 1/a;
      //    a = 1.2; c = 0.99;
      a = 6378172; b = 6378103; c = 6356753;
      a -= 34; b += 34;
      TriaxialTest1(a, b, c);
    }
    if (0) {
      real a = 1.01, b = 1, c = 0.8;
      a = 6378172; b = 6378103; c = 6356753;
      // a -= 34; b += 34;
      UmbilicTest(a, b, c);
    }
    if (1) {
      real a = 1.01, b = 1, c = 0.8;
      a = 6378172; b = 6378103; c = 6356753;
      //      a -= 34; b += 34;
      DirectfunTest(a, b, c);
      a = sqrt(real(2)); b = 1; c = 1/a;
      DirectfunTest(a, b, c);
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

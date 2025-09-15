// Example of using the Triaxial::Conformal3 class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>
#include <GeographicLib/Triaxial/Conformal3.hpp>
#include <GeographicLib/Triaxial/Geodesic3.hpp>

using namespace std;
using namespace GeographicLib;
using namespace Triaxial;

int main() {
  try {
    Utility::set_digits();
    typedef Math::real real;
    if (0) {
      real kp2 = Math::sq(numeric_limits<real>::epsilon());
      EllipticFunction ell(1-kp2, 0, kp2, 1);
      cout << "VALS " << ell.K() << " " << ell.E() << "\n";
      real nx = 0.75, ny = 0.25;
      auto ksolve = [nx, ny]
        (real k2) -> pair<real, real>
        {
          real kp2 = 1 - k2;
          EllipticFunction ell(k2);
          EllipticFunction ellp(kp2, 0, k2, 1);
          real f = nx * ell.K() - ny * ellp.K(),
          fp = (nx * (ell .E() - kp2 * ell .K()) +
                ny * (ellp.E() - k2  * ellp.K())) / (2 * k2 * kp2);
          return pair<real, real>(f, fp);
        };
      auto x = ksolve(kp2);
      cout << "KK " << x.first << " " << x.second << "\n";
      return 0;
    }
    real a = Constants::Triaxial_Earth_a(),
       b = Constants::Triaxial_Earth_b(),
      c = Constants::Triaxial_Earth_c();
    Ellipsoid3 t(a, b, c);
    Conformal3 p3(t);
    Geodesic3 g3(t);
    Geodesic g2(p3.SphereRadius(), 0);
    cout  << fixed << setprecision(1)
          << "Ellipsoid parameters: a = " << t.a()
          << ", b = " << t.b() << ", c = " << t.c() << "\n"
          << setprecision(8)
          << "Quadrants: x = " << p3.x0()/t.b()
          << ", y = " << p3.y0()/t.b() << "\n";
    cout << "PARAMS\n"
      << t.b() << " " << t.e2() << " "
      << t.k2() << " " << t.kp2() << "\n"
      << p3.s().b() << " " << p3.s().e2() << " "
      << p3.s().k2() << " " << p3.s().kp2() << "\n"
      << p3.s0().b() << " " << p3.s0().e2() << " "
      << p3.s0().k2() << " " << p3.s0().kp2() << "\n";
    if (0) {
      if (t.kp2() == 0) {
        real e = sqrt(t.e2()), sg = sinh(e * atanh(e));
        cout << "SCALEO " << (t.c()/t.b()) * (hypot(sg, real(1)) + sg) << "\n";
      } else if (t.k2() == 0) {
        real e = sqrt(t.e2()), sg = sinh(e * atan(e));
        cout << "SCALEP " << (t.a()/t.b()) / (hypot(sg, real(1)) + sg) << "\n";
      } else
        cout << "SCALEX " << sqrt(p3.s().k2() * p3.s().kp2()/(t.k2() * t.kp2()))
          * (t.b()/p3.s().b()) << "\n";
    }
    Angle bet1{90-real(0)/10000}, omg1{0+real(0)/10000}, phi1, lam1, gam;
    real s12 = 0.1, m;
    p3.ForwardSphere(bet1, omg1, phi1, lam1, gam, m);
    Angle bet1q, omg1q, gamq; real mq;
    p3.ReverseSphere(phi1, lam1, bet1q, omg1q, gamq, mq);
    cout << "CONVSCALE " << real(gam) << " " << m << "\n";
    cout << "CONVSCALEQ " << real(gamq) << " " << mq << " "
         << real(bet1q) << " " << real(omg1q) << "\n";
    if (0) cout << "CC0 " << real(bet1) << " " << real(omg1) << " "
                << real(phi1) << " " << real(lam1) << "\n";
    for (int i = -180; i < 180; i +=45) {
      Angle alp1{real(i)}, bet2, omg2, alp2;
      g3.Direct(bet1, omg1, alp1, s12, bet2, omg2, alp2);
      Angle alp1x, alp2x;
      real s12x;
      g3.Inverse(bet1, omg1, bet2, omg2, s12x, alp1x, alp2x);
      Angle phi2, lam2, gam2; real m2;
      p3.ForwardSphere(bet2, omg2, phi2, lam2, gam2, m2);
      if (0) cout << "CC " << real(bet2) << " " << real(omg2) << " "
                  << real(phi2) << " " << real(lam2) << "\n";
      real alp1s, alp2s, s12s;
      g2.Inverse(real(phi1), real(lam1), real(phi2), real(lam2),
                 s12s, alp1s, alp2s);
      Angle gamx{alp1s}; gamx -= alp1;
      cout << i << " " << real(gamx.base()) << " " << (s12s/m)/s12 << "\n";
    }
    return 0;
  }  
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

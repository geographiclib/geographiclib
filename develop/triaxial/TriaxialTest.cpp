#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"
#include "TriaxialODE.hpp"

using namespace GeographicLib;
using namespace std;

void ODEtest(Math::real a, Math::real b, Math::real c) {
    typedef Math::real real;
    typedef Triaxial::vec3 vec3;
    Triaxial t(a, b, c);
    vec3 r1{t.axes()[0] * sqrt(t.kp2()), 0, t.axes()[2] * sqrt(t.k2())},
      v1{0,1,0};
    real s12 = 4 * EllipticFunction::RG(Math::sq(t.axes()[0]),
                                        Math::sq(t.axes()[2]));
    vec3 r2, v2;
    int kmax = Math::digits()-4;
    TriaxialODE direct(t, r1, v1);
    if (1) {
      for (int k = 10; k <= kmax; ++k) {
        real eps = pow(real(2), -k);
        int n = direct.Position(s12, r2, v2, eps);
        cout << k << " " << n << " "
             << sqrt(Math::sq(r1[0]+r2[0]) +
                     Math::sq(r1[1]+r2[1]) +
                     Math::sq(r1[2]+r2[2])) / t.b() << endl;
      }
    } else {
      int n = 0, imax = 1000;
#if GEOGRAPHICLIB_PRECISION >= 4
      real eps = pow(numeric_limits<real>::epsilon(), real(0.875));
#else
      real eps = numeric_limits<real>::epsilon() * 16;
#endif
      for (int i = 0; i < imax; ++i) {
        n += direct.Position(s12, r2, v2, eps);
        r1 = r2;
        v1 = v2;
      }
      cout << n/real(imax) << "\n";
    }
}

void DirectfunTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  typedef Angle ang;
  Triaxial t(a, b, c);
  // Circumpolar
  TriaxialLine lca(t, ang(1), ang(0), ang(90));
  if (1) {
  TriaxialLine lcb(t, ang(real(89.999)), ang(0), ang(90));
  }
  // Umbilic
  TriaxialLine lu(t, ang(90), ang(0), ang(135));
  if (1) {
  // Circumpolar
  TriaxialLine ltb(t, ang(90), ang(real(0.001)), ang(180));
  }
  TriaxialLine lta(t, ang(90), ang(89), ang(180));
}

void PositionTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  typedef Angle ang;
  Triaxial t(a, b, c);
  ang bet1, omg1, alp1, bet2, omg2, alp2;
  bet1 = ang(-0.0); omg1 = ang(105); alp1 = ang(90);
  //  bet1 = ang(1); omg1 = ang(2); alp1 = ang(3);
  TriaxialLine l(t, bet1, omg1, alp1);
  cout << fixed << setprecision(6);
  real ds = 1/real(10);
  //  for (int s12 = 0; s12 <= 30; ++s12) {
  {real s12 = 0.996504/ds;
    l.Position(s12*ds, bet2, omg2, alp2);
    (void) Triaxial::AngNorm(bet2, omg2, alp2);
    cout << s12 << " "
         << real(omg2) << "\n";
  }
}

int main() {
  try {
    Utility::set_digits();
    typedef Math::real real;
    typedef Angle ang;
   if (0) {
      Triaxial t(sqrt(2.0), 1.0, 1/sqrt(2.0));
      TriaxialLine l = t.Inverse(
                     ang(-45.0),ang(-30.0),
                     ang(-90.0),ang(0.0));
      real s12 = l.Distance();
      /*      real s12 = 0.361729;
      TriaxialLine l(t, ang(-45.0),ang(-30.0),
                      ang(135.0));
      l.SetDistance(s12);
      */
      cout << s12 << "\n";
      ang bet2, omg2, alp2;
      l.pos1(bet2, omg2, alp2);
        cout << "POS1 " << s12 << " "
             << real(bet2) << " "
             << real(omg2) << " "
             << real(alp2) << "\n";
      cout << fixed << setprecision(6);
      for (int i = -3; i <= 3; ++i) {
        real ss = 0 + i * 1e-6;
        l.Position(ss, bet2, omg2, alp2);
        cout << i << " " << ss << " "
             << real(bet2) << " "
             << real(omg2) << " "
             << real(alp2) << "\n";
      }
      for (int i = -3; i <= 3; ++i) {
        real ss = s12 + i * 1e-6;
        l.Position(ss, bet2, omg2, alp2);
        cout << i << " " << ss << " "
             << real(bet2) << " "
             << real(omg2) << " "
             << real(alp2) << "\n";
      }
      return 0;
    }
    if (0) {
    {
      TriaxialLine ll(Triaxial(sqrt(2.0), 1.0, 1/sqrt(2.0)),
                      ang(90), ang(180),
                      ang(158.253574));
      ll.SetDistance(0.898324);
      ang bet1x, omg1x, alp1x, bet2x, omg2x, alp2x;
      ll.pos1(bet1x, omg1x, alp1x);
      real s12x = ll.Distance();
      ll.Position(s12x, bet2x, omg2x, alp2x);
      std::cout << "OUT\n"
                << real(bet1x) << " " << real(omg1x) << " "
                << real(alp1x) << " "
                << signbit(bet1x.c()) << signbit(omg1x.s()) << "\n"
                << real(bet2x) << " " << real(omg2x) << " "
                << real(alp2x) << "\n"
                << s12x << "\n";
    }
    return 0;
        {
      TriaxialLine ll(Triaxial(sqrt(2.0), 1.0, 1/sqrt(2.0)),
                      ang(-90), ang(0),
                      ang(-87.089900));
      ll.SetDistance(2.736330);
      ang bet1x, omg1x, alp1x, bet2x, omg2x, alp2x;
      ll.pos1(bet1x, omg1x, alp1x);
      real s12x = ll.Distance();
      ll.Position(s12x, bet2x, omg2x, alp2x);
      std::cout << "OUT\n"
                << real(bet1x) << " " << real(omg1x) << " "
                << real(alp1x) << " "
                << signbit(bet1x.c()) << signbit(omg1x.s()) << "\n"
                << real(bet2x) << " " << real(omg2x) << " "
                << real(alp2x) << "\n"
                << s12x << "\n";
    }
    return 0;
    }

    real a = 6378172, b = 6378103, c = 6356753;
    if (0) {
      using std::sqrt;
      a = sqrt(real(2)); b = 1; c = 1/a;
      //      real a = 1.01, b = 1, c = 0.8;
      //      a = 6378172; b = 6378103; c = 6356753;
      //      a -= 34; b += 34;
      PositionTest(a, b, c);
      }
    if (1)
      ODEtest(a, b, c);
    /*
    if (0)
      TriaxialTest0();
    */
    if (0) {
      real a = 1.01, b = 1, c = 0.8;
      a = 6378172; b = 6378103; c = 6356753;
      a -= 34; b += 34;
      DirectfunTest(a, b, c);
      //      a = sqrt(real(2)); b = 1; c = 1/a;
      //      DirectfunTest(a, b, c);
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

/**
 * \file elliptictest.cpp
 * \brief Test EllipticFunction
 *
 * Copyright (c) Charles Karney (2026) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <iostream>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>

using namespace std;
using namespace GeographicLib;

static int check(Math::real x, Math::real y) {
  static const Math::real m = 0.5e-13;
  Math::real d = y == 0 ? 0 : (fabs(y) < 1 ? m/10 : m);
  if (fabs(x - y) <= d)
    return 0;
  cout << "check fails: " << x << " != " << y << " +/- " << d << "\n";
  return 1;
}

static int check(Math::cmplx x, Math::cmplx y) {
  static const Math::real m = 0.5e-13;
  Math::cmplx d(y.real() == 0 ? 0 : (fabs(y.real()) < 1 ? m/10 : m),
                y.imag() == 0 ? 0 : (fabs(y.imag()) < 1 ? m/10 : m));
  if (fabs(x.real() - y.real()) <= d.real() &&
      fabs(x.imag() - y.imag()) <= d.imag())
    return 0;
  cout << "check fails: " << x << " != " << y << " +/- " << d << "\n";
  return 1;
}

int dotests() {
  // The numerical checks in Carlson 1995, Sec 3.

  //   RF(1, 2, 0) = 1.3110287771461
  //   RF(i, -i, 0) = RF(0.5, 1, 0) = 1.8540746773014
  //   RF(i-1, i, 0) = 0.79612586584234 - i*1.2138566698365
  //   RF(2, 3, 4) = 0.58408284167715
  //   RF(i, -i, 2) = 1.0441445654064
  //   RF(i-1, i, 1-i) = 0.93912050218619 - i*0.53296252018635

  //   RC(0, 1/4) = pi = 3.1415926535898
  //   RC(9/4, 2) = ln2 = 0.69314718055995
  //   RC(0, i) = (1-i) * 1.1107207345396
  //   RC(-i, i) = 1.2260849569072 - i*0.34471136988768
  //   RC(1/4, -2) = ln2/3 = 0.23104906018665
  //   RC(i, -1) = 0.77778596920447 + i*0.19832484993429

  //   RJ(0, 1, 2, 3) = 0.77688623778582
  //   RJ(2, 3, 4, 5) = 0.14297579667157
  //   RJ(2, 3, 4, -1+i) = 0.13613945827771 - i*0.38207561624427
  //   RJ(i, -i, 0, 2) = 1.6490011662711
  //   RJ(-1+i, -1-i, 1, 2) = 0.94148358841220
  //   RJ(i, -i, 0, 1-i) = 1.8260115229009 + i*1.2290661908643
  //   RJ(-1+i, -1-i, 1, -3+i) = -0.61127970812028 - i*1.0684038390007
  //   RJ(-1+i, -2-i, -i, -1+i) = 1.8249027393704 - i*1.2218475784827
  //   RJ(2, 3, 4, -0.5) = 0.24723819703052
  //   RJ(2, 3, 4, -5) = -0.12711230042964

  //   RD(0, 2, 1) = 1.7972103521034
  //   RD(2, 3, 4) = 0.16510527294261
  //   RD(i, -i, 2) = 0.65933854154220
  //   RD(0, i, -i) = 1.2708196271910 + i*2.7811120159521
  //   RD(0, i-1, i) = -1.8577235439239 - i*0.96193450888839
  //   RD(-2-i, -i, -1+i) = 1.8249027393704 - i*1.2218475784827

  //   RG(0, 16, 16) = 2E(0) = pi = 3.1415926535898
  //   RG(2, 3, 4) = 1.7255030280692
  //   RG(0, i, -i) = 0.42360654239699
  //   RG(i-1, i, 0) = 0.44660591677018 + i*0.70768352357515
  //   RG(-i, i-1, i) = 0.36023392184473 + i*0.40348623401722
  //   RG(0, 0.0796, 4) = E(0.99) = 1.0284758090288
  static const int ncases = 35;
  static const Math::real testcases[ncases][11] = {
    //c    ==x==    ==y==   ==z==   ==p==              ==result==
    // RF tests
    {0,   1  , 0,  2   , 0,  0, 0,  0 , 0,   1.3110287771461, 0              },
    {0,   0.5, 0,  1   , 0,  0, 0,  0 , 0,   1.8540746773014, 0              },
    {0,   0  , 1,  0   ,-1,  0, 0,  0 , 0,   1.8540746773014, 0              },
    {0,  -1  , 1,  0   , 1,  0, 0,  0 , 0,   .79612586584234,-1.2138566698365},
    {0,   2  , 0,  3   , 0,  4, 0,  0 , 0,   .58408284167715, 0              },
    {0,   0  , 1,  0   ,-1,  2, 0,  0 , 0,   1.0441445654064, 0              },
    {0,  -1  , 1,  0   , 1,  1,-1,  0 , 0,   .93912050218619,-.53296252018635},
    // RC tests
    {1,   0  , 0,  0.25, 0,  0, 0,  0 , 0,   3.1415926535898, 0              },
    {1,  2.25, 0,  2   , 0,  0, 0,  0 , 0,   .69314718055995, 0              },
    {1,   0  , 0,  0   , 1,  0, 0,  0 , 0,   1.1107207345396,-1.1107207345396},
    {1,   0  ,-1,  0   , 1,  0, 0,  0 , 0,   1.2260849569072,-.34471136988768},
    {1,  0.25, 0, -2   , 0,  0, 0,  0 , 0,   .23104906018665, 0              },
    {1,   0  , 1, -1   , 0,  0, 0,  0 , 0,   .77778596920447, .19832484993429},
    // RJ tests
    {2,   0  , 0,  1   , 0,  2, 0,  3 , 0,   .77688623778582, 0              },
    {2,   2  , 0,  3   , 0,  4, 0,  5 , 0,   .14297579667157, 0              },
    {2,   2  , 0,  3   , 0,  4, 0, -1 , 1,   .13613945827771,-.38207561624427},
    {2,   0  , 1,  0   ,-1,  0, 0,  2 , 0,   1.6490011662711, 0              },
    {2,  -1  , 1, -1   ,-1,  1, 0,  2 , 0,   .94148358841220, 0              },
    {2,   0  , 1,  0   ,-1,  0, 0,  1 ,-1,   1.8260115229009, 1.2290661908643},
    {2,  -1  , 1, -1   ,-1,  1, 0, -3 , 1,  -.61127970812028,-1.0684038390007},
    {2,  -1  , 1, -2   ,-1,  0,-1, -1 , 1,   1.8249027393704,-1.2218475784827},
    {2,   2  , 0,  3   , 0,  4, 0, -.5, 0,   .24723819703052, 0              },
    {2,   2  , 0,  3   , 0,  4, 0, -5 , 0,  -.12711230042964, 0              },
    // RD tests
    {3,   0  , 0,  2   , 0,  1, 0,  0 , 0,   1.7972103521034, 0              },
    {3,   2  , 0,  3   , 0,  4, 0,  0 , 0,   .16510527294261, 0              },
    {3,   0  , 1,  0   ,-1,  2, 0,  0 , 0,   .65933854154220, 0              },
    {3,   0  , 0,  0   , 1,  0,-1,  0 , 0,   1.2708196271910, 2.7811120159521},
    {3,   0  , 0, -1   , 1,  0, 1,  0 , 0,  -1.8577235439239,-.96193450888839},
    {3,  -2  ,-1,  0   ,-1, -1, 1,  0 , 0,   1.8249027393704,-1.2218475784827},
    // RG tests
    {4,   0  , 0, 16   , 0, 16, 0,  0 , 0,   3.1415926535898, 0              },
    {4,   2  , 0,  3   , 0,  4, 0,  0 , 0,   1.7255030280692, 0              },
    {4,   0  , 0,  0   , 1,  0,-1,  0 , 0,   .42360654239699, 0              },
    {4,  -1  , 1,  0   , 1,  0, 0,  0 , 0,   .44660591677018, .70768352357515},
    {4,   0  ,-1, -1   , 1,  0, 1,  0 , 0,   .36023392184473, .40348623401722},
    {4,   0  , 0, .0796, 0,  4, 0,  0 , 0,   1.0284758090288, 0              },
  };
  int n = 0;
  for (int j = 0; j < ncases; ++j) {
    int c = int(testcases[j][0]);
    Math::cmplx
      x(testcases[j][1], testcases[j][2]),
      y(testcases[j][3], testcases[j][4]),
      z(testcases[j][5], testcases[j][6]),
      p(testcases[j][7], testcases[j][8]),
      e(testcases[j][9], testcases[j][10]);
    switch (c) {
    case 0:
      if (x.imag() == 0 && y.imag() == 0 && z.imag() == 0)
        n += check(EllipticFunction::RF(x.real(), y.real(), z.real()),
                   e.real());
      n += check(EllipticFunction::RF(x, y, z), e);
      break;
    case 1:
      if (x.imag() == 0 && y.imag() == 0)
        n += check(EllipticFunction::RC(x.real(), y.real()), e.real());
      n += check(EllipticFunction::RC(x, y), e);
      break;
    case 2:
      if (x.imag() == 0 && y.imag() == 0 && z.imag() == 0 && p.imag() == 0)
        n += check(EllipticFunction::RJ(x.real(), y.real(), z.real(), p.real()),
                   e.real());
      n += check(EllipticFunction::RJ(x, y, z, p), e);
      break;
    case 3:
      if (x.imag() == 0 && y.imag() == 0 && z.imag() == 0)
        n += check(EllipticFunction::RD(x.real(), y.real(), z.real()),
                   e.real());
      n += check(EllipticFunction::RD(x, y, z), e);
      break;
    case 4:
    default:
      if (x.imag() == 0 && y.imag() == 0 && z.imag() == 0)
        n += check(EllipticFunction::RG(x.real(), y.real(), z.real()),
                   e.real());
      n += check(EllipticFunction::RG(x, y, z), e);
      break;
    }
  }
  return n;
}

int main() {
  Utility::set_digits();
  {
    EllipticFunction ell(0.5);
#if GEOGRAPHICLIB_COMPLEX_JACOBI_AM
    Math::real K = ell.K(), d = K/45;
    for (int ii = -45; ii <= 90; ++ii) {
      for (int ir = -45; ir <= 45; ++ir) {
        Math::cmplx phi = {ir*d, ii*d}, z = ell.am(phi);
        cout << " " << z.real() << " " << z.imag();
      }
      cout << "\n";
    }
    return 0;
#endif
    for (int i = -730; i <= 730; i+= 10) {
    // { int i = 90;
      Math::real phi = i * Math::degree(),
        u = ell.F(phi),
        sn, cn, dn,
        phia = ell.am(u, sn, cn, dn);
      cout << i << " " << phia-phi << "\n";
    }
    return 0;
  }
  int n = dotests();
  if (n) {
    cout << n << " failure" << (n > 1 ? "s" : "") << "\n";
      return 1;
  }
}

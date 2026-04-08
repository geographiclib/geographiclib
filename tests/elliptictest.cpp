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
  static const int ncases = 35;
  static const Math::real testcases[ncases][11] = {
    //c,   xr,xi,    yr,yi, zr,zi,  pr,pi,                er,ei
    {0,   1  , 0,  2   , 0,  0, 0,  0 , 0,   1.3110287771461, 0              },
    {0,   0.5, 0,  1   , 0,  0, 0,  0 , 0,   1.8540746773014, 0              },
    {0,   0  , 1,  0   ,-1,  0, 0,  0 , 0,   1.8540746773014, 0              },
    {0,  -1  , 1,  0   , 1,  0, 0,  0 , 0,   .79612586584234,-1.2138566698365},
    {0,   2  , 0,  3   , 0,  4, 0,  0 , 0,   .58408284167715, 0              },
    {0,   0  , 1,  0   ,-1,  2, 0,  0 , 0,   1.0441445654064, 0              },
    {0,  -1  , 1,  0   , 1,  1,-1,  0 , 0,   .93912050218619,-.53296252018635},
    {1,   0  , 0,  0.25, 0,  0, 0,  0 , 0,   3.1415926535898, 0              },
    {1,  2.25, 0,  2   , 0,  0, 0,  0 , 0,   .69314718055995, 0              },
    {1,   0  , 0,  0   , 1,  0, 0,  0 , 0,   1.1107207345396,-1.1107207345396},
    {1,   0  ,-1,  0   , 1,  0, 0,  0 , 0,   1.2260849569072,-.34471136988768},
    {1,  0.25, 0, -2   , 0,  0, 0,  0 , 0,   .23104906018665, 0              },
    {1,   0  , 1, -1   , 0,  0, 0,  0 , 0,   .77778596920447, .19832484993429},
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
    {3,   0  , 0,  2   , 0,  1, 0,  0 , 0,   1.7972103521034, 0              },
    {3,   2  , 0,  3   , 0,  4, 0,  0 , 0,   .16510527294261, 0              },
    {3,   0  , 1,  0   ,-1,  2, 0,  0 , 0,   .65933854154220, 0              },
    {3,   0  , 0,  0   , 1,  0,-1,  0 , 0,   1.2708196271910, 2.7811120159521},
    {3,   0  , 0, -1   , 1,  0, 1,  0 , 0,  -1.8577235439239,-.96193450888839},
    {3,  -2  ,-1,  0   ,-1, -1, 1,  0 , 0,   1.8249027393704,-1.2218475784827},
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
  int n = dotests();
  if (n) {
    cout << n << " failure" << (n > 1 ? "s" : "") << "\n";
      return 1;
  }
}

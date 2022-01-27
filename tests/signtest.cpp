#include <iostream>
#include <GeographicLib/Math.hpp>

using namespace std;
using namespace GeographicLib;

int atantest(Math::real y, Math::real x, Math::real r) {
  r *= 180;
  Math::real r1 = Math::atan2d(y, x);
  int e = 0;
  if (isnan(r)) {
    if (!isnan(r1)) ++e;
  } else {
    if (!( r == r1 && signbit(r) == signbit(r1) )) ++e;
  }
  if (e)
    cout << "atan2d(" << y << ", " << x << ") != " << r << " (" << r1 << ")\n";
  return e;
}

int main() {
  Math::real inf = Math::infinity(), nan = Math::NaN();
  int n = 0;
  n += atantest(+0.0, -0.0, +1.0);
  n += atantest(-0.0, -0.0, -1.0);
  n += atantest(+0.0, +0.0, +0.0);
  n += atantest(-0.0, +0.0, -0.0);
  n += atantest(+0.0, -1.0, +1.0);
  n += atantest(-0.0, -1.0, -1.0);
  n += atantest(+0.0, +1.0, +0.0);
  n += atantest(-0.0, +1.0, -0.0);
  n += atantest(-1.0, +0.0, -0.5);
  n += atantest(-1.0, -0.0, -0.5);
  n += atantest(+1.0, +0.0, +0.5);
  n += atantest(+1.0, -0.0, +0.5);
  n += atantest(+1.0, -inf, +1.0);
  n += atantest(-1.0, -inf, -1.0);
  n += atantest(+1.0, +inf, +0.0);
  n += atantest(-1.0, +inf, -0.0);
  n += atantest(+inf, +1.0, +0.5);
  n += atantest(+inf, -1.0, +0.5);
  n += atantest(-inf, +1.0, -0.5);
  n += atantest(-inf, -1.0, -0.5);
  n += atantest(+inf, -inf, +0.75);
  n += atantest(-inf, -inf, -0.75);
  n += atantest(+inf, +inf, +0.25);
  n += atantest(-inf, +inf, -0.25);
  n += atantest( nan, +1.0,  nan);
  n += atantest(+1.0,  nan,  nan);
  if (n) {
    cout << n << " failure" << (n > 1 ? "s" : "") << "\n";
    return 1;
  }
}

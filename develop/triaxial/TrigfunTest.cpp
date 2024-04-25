#include <iostream>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include "Trigfun.hpp"

using namespace GeographicLib;
using namespace std;

void BasicTest() {
  typedef Math::real real;
  int n = 8;
  real pi = Math::pi<real>();
  vector<real> f1(n+1);       // n+1 samples
  for (int i = 0; i <= n; ++i)
    f1[i] = sqrt(real(i)+1);
  vector<real> f0(f1); f0.resize(n); // n samples
  vector<real> f2(f0); f2[n-1]=0;      // n samples, last one is 0
  if (0) {
    Trigfun::initbysamples(f1, false, false, false, pi);
    Trigfun::initbysamples(f0, false, false, true, pi);
    Trigfun::initbysamples(f2, true, false, false, pi);
    Trigfun::initbysamples(f0, true, false, true, pi);
    Trigfun::initbysamples(f0, false, true, false, pi);
    Trigfun::initbysamples(f0, false, true, true, pi);
    Trigfun::initbysamples(f0, true, true, false, pi);
    Trigfun::initbysamples(f0, true, true, true, pi);
    cout << "++++++\n";
  }
  vector<real> f1a(n/2 + 1);
  vector<real> f1b(n/2);
  for (int i = 0; i <= n; i += 2)
    f1a[i/2] = f1[i];
  for (int i = 1; i <= n; i += 2)
    f1b[i/2] = f1[i];
  if (1) {
    {
      Trigfun
        t = Trigfun::initbysamples(f1, false, false, false, pi),
        ta = Trigfun::initbysamples(f1a, false, false, false, pi),
        tb = Trigfun::initbysamples(f1b, false, false, true, pi);
      ta.refine(tb);
      cout << ta.check(f1, false) << "\n";
      cout << "++++++\n";
    }
    vector<real> f2a(f1b); f2a.resize(n/2); f2a[n/2-1] = 0;
    vector<real> f2b(f1a); f2b.resize(n/2);
    {
      Trigfun
        t = Trigfun::initbysamples(f2, true, false, false, pi),
        ta = Trigfun::initbysamples(f2a, true, false, false, pi),
        tb = Trigfun::initbysamples(f2b, true, false, true, pi);
      ta.refine(tb);
      cout << ta.check(f2, false) << "\n";
      cout << "++++++\n";
    }
    vector<real> f0a(f1a); f0a.resize(n/2);
    vector<real> f0b(f1b);
    {
      Trigfun
        t = Trigfun::initbysamples(f0, false, true, false, pi),
        ta = Trigfun::initbysamples(f0a, false, true, false, pi),
        tb = Trigfun::initbysamples(f0b, false, true, true, pi);
      ta.refine(tb);
      cout << ta.check(f0, false) << "\n";
      cout << "++++++\n";
    }
    {
      Trigfun
        t = Trigfun::initbysamples(f0, true, true, false, pi),
        ta = Trigfun::initbysamples(f0b, true, true, false, pi),
        tb = Trigfun::initbysamples(f0a, true, true, true, pi);
      ta.refine(tb);
      cout << ta.check(f0, false) << "\n";
      cout << "++++++\n";
    }
  }
}

void EllipTest(Math::real k2, Math::real kp2) {
  typedef Math::real real;
  EllipticFunction ell(k2, 0, kp2, 1);
  auto f = [&ell] (real x) -> real {
    return 1/ell.Delta(sin(x), cos(x)); };
  Trigfun t(f, false, false, false, Math::pi()/2, 0);
  Trigfun F(t.integral());
  Trigfun Finv(F.invert(f));
  int countn = 0, countb = 0;
  for (int i = 0; i <= 100; ++i) {
    real z = i/real(10),
      x = Finv(z),
      xa = F.root(z, f, countn, countb),
      xc = ell.am(z);
    cout << "roots " << z << " " << x << " "
         << z-F(x) << " " << x  - xa <<  " "
         << x - xc << "\n";
  }
  return;
  if (0) {
    for (int i = 0; i <= 100; ++i) {
      real x = i/real(10);
      cout << i << " " << t(x) << " " << F(x) << " " << F(x) - ell.F(x) << "\n";
    }
  }
  if (0) {
    real z = 6,
      x = F.root(z, f, countn, countb);
    cout << "roots " << z << " " << x << " " << z - F(x) << "\n";
    cout << "counts " << countn << " " << countb << "\n";
    return;
  }
  for (int i = 0; i <= 100; ++i) {
    real z = i/real(10),
      x = F.root(z, f, countn, countb);
    cout << "roots " << z << " " << x << " " << z-F(x) << "\n";
  }
  cout << "counts " << countn << " " << countb << "\n";
}

int main() {
  try {
    Utility::set_digits();
    typedef Math::real real;
    if (0) BasicTest();
    real kp2 = 0.00001, k2 = 1 - kp2;
    //    k2 = -0.1; kp2 = 1 - k2;
    EllipTest(k2, kp2);
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

#include <iostream>
#include <limits>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include "Trigfun.hpp"
#include "Triaxial.hpp"

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

void TriaxialTest(Math::real a, Math::real b, Math::real c) {
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
                 << fphia.ell.K() << " "
                 << fphia.NCoeffs() << " " << fphia.NCoeffsInv() << " "
                 << fomga.ell.K() << " "
                 << fomga.NCoeffs() << " " << fomga.NCoeffsInv() << "\n";
    }
  }
}
int main() {
  try {
    Utility::set_digits();
    typedef Math::real real;
    if (0) BasicTest();
    if (0) {
      real kp2 = 0.00001, k2 = 1 - kp2;
      //    k2 = -0.1; kp2 = 1 - k2;
      EllipTest(k2, kp2);
    }
    if (0) {
      TriaxialTest0();
    }
    real a = 1.01, b = 1, c = 0.8;
    //    a = sqrt(2.0); c = 1/a;
    //    a = 1.2; c = 0.99;
    a = 6378172; b = 6378103; c = 6356753;
    TriaxialTest(a-34, b+34, c);
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    cerr << "Caught unknown exception\n";
    return 1;
  }
}

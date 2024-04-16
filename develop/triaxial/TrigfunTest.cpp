#include <iostream>
#include <GeographicLib/Utility.hpp>
#include "Trigfun.hpp"

using namespace GeographicLib;
using namespace std;
int main() {
  try {
    typedef Math::real real;
    Utility::set_digits();
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
        ta.refine(tb, t);
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
        ta.refine(tb, t);
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
        ta.refine(tb, t);
        cout << ta.check(f0, false) << "\n";
        cout << "++++++\n";
      }
      {
        Trigfun
          t = Trigfun::initbysamples(f0, true, true, false, pi),
          ta = Trigfun::initbysamples(f0b, true, true, false, pi),
          tb = Trigfun::initbysamples(f0a, true, true, true, pi);
        ta.refine(tb, t);
        cout << ta.check(f0, false) << "\n";
        cout << "++++++\n";
      }
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

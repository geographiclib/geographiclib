// Example of using the GeographicLib::DST class

#include <iostream>
#include <exception>
#include <vector>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/DST.hpp>

using namespace std;
using namespace GeographicLib;

class sawtooth {
private:
  double _a, _b;
public:
  sawtooth(double a) : _a(a) {}
  // only called for x in (0, pi/2].  DST assumes function is periodic, period
  // 2*pi, is odd about 0, and is even about pi/2.
  double operator()(double x) const { return _a * x; }
};

int main() {
  try {
    sawtooth f(Math::pi()/4);
    DST dst;
    vector<double> tx, txa;
    int N = 8;
    dst.transform(f, N, tx);
    cout << "Transform of sawtooth based on " << N << " points\n"
         << "approx 1, -1/9, 1/25, -1/49, ...\n";
    for (int i = 0; i < min(10,N); ++i)
      cout << tx[i] << "\n";
    dst.refine(f, tx, tx);
    cout << "Add another " << N << " points\n";
    for (int i = 0; i < min(10,N); ++i)
      cout << tx[i] << "\n";
    dst.transform(f, 2*N, txa);
    cout << "Retransform of sawtooth based on " << 2*N << " points\n";
    for (int i = 0; i < min(10,N); ++i)
      cout << txa[i] << "\n";
    int M = 10;
    cout << "Table of values and integral\n";
    for (int i = 0; i <= M; ++i) {
      double x = i*Math::pi()/(2*M), sinx = sin(x), cosx = cos(x);
      cout << x << " "
           << DST::eval(txa, sinx, cosx) << " "
           << DST::integral(txa, sinx, cosx) << "\n";
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

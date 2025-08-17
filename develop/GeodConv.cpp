// Generate test data for triaxial ellipsoid from WGS84 test data.

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/AuxLatitude.hpp>

using namespace GeographicLib;
using namespace std;

string nicestr(Math::real x, int prec, bool azi = false) {
  typedef Math::real real;
  // Extra real case because pow(float, int) returns a double
  static const real eps = real(pow(real(10), -20));
  Math::real y = round(x);
  if (fabs(x - y) <= eps)
    x = azi && y == -180 ? 180 : y + real(0);
  ostringstream os;
  os << fixed << setprecision(prec) << x;
  string s = os.str();
  string::size_type p = s.find_last_not_of('0');
  if (p == string::npos)
    s = "0";
  else {
    if (s[p] == '.')
      --p;
    s = s.substr(0, p+1);
  }
  return s;
}

Math::real redlat(Math::real lat, Math::real n) {
  Math::real s, c;
  Math::sincosd(lat, s, c);
  return Math::atan2d(s * n, c);
}

Math::real rnd(Math::real ang) {
  typedef Math::real real;
  // Real casts because 1000000000000 isn't exactly representable as a float.
  return round(real(1000000000000) * ang) / real(1000000000000);
}

int main() {
  try {
    typedef Math::real real;
    Utility::set_digits();
    real a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    // f = 1/298.257223563
    a = 1;
    f = 124 / ( sqrt(real(340804677)) + 18523 );
    // f = 1/298.2572249059396272
    // relative error = 1 part in 2e8
    // e2 = f * (2-f) = 124/18523
    real lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, S12;
    real bet1, bet2, M12, M21;
    int prec = 18;
    Geodesic g(a, f, true);
    while (cin
           >> lat1 >> lon1 >> azi1
           >> lat2 >> lon2 >> azi2
           >> s12 >> a12 >> m12 >> S12) {
      lon1 = rnd(lon1); lon2 = rnd(lon2);
      bet1 = rnd(redlat(lat1, 1-f)); bet2 = rnd(redlat(lat2, 1-f));
      g.Inverse(redlat(bet1,1/(1-f)), lon1, redlat(bet2,1/(1-f)), lon2,
                s12, azi1, azi2, m12, M12, M21);
      cout << nicestr(bet1, prec) << " " << nicestr(lon1, prec) << " "
           << nicestr(azi1, prec, true) << " "
           << nicestr(bet2, prec) << " " << nicestr(lon2, prec) << " "
           << nicestr(azi2, prec, true) << " "
           << nicestr(s12, prec+2) << " " << nicestr(m12, prec+2) << " "
           << nicestr(M12, prec+2) << " " << nicestr(M21, prec+2) << "\n";
    }
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

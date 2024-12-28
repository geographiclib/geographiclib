#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
#include <GeographicLib/Utility.hpp>

using namespace std;
using namespace GeographicLib;

int main(int argc, char* argv[]) {
  Utility::set_digits();
  int prec = numeric_limits<double>::digits10 - 1;
  if (argc == 2) {
    string arg(argv[1]);
    prec = Utility::val<int>(arg);
  }
  string s;
  // readarray is used to read in coefficient data rapidly.  Thus 8.3n is
  // stored in its IEEE double representation.  This is fine is the working
  // precision is double.  However, when working at higher precision, how
  // should be interpret the constant 8.3 appearing in a published table?
  // Possibilities are
  //
  // (a) treat this as an exact decimal number 83/10;
  //
  // (b) treat this as the approximate decimal representation of an exact
  // double precision number 2336242306698445/2^48 =
  // 8.300000000000000710542735760100185871124267578125
  //
  // Here use (a) if the number of significant digits in the number is 15 or
  // less.  Otherwise, we use (b).
  //
  // We implement this as follows.  Any double which can be represented as a
  // decimal number with precision 14 = digis10 - 1 (= 15 sig figs) is treated
  // as an approximation to that decimal number.  The high precision number is
  // then obtained by reading the decimal number at that precision.  Otherwise
  // the double is treated as exact.  The high precision number is obtained by
  // adding zeros in the binary fraction.
  //
  // N.B. printing with precision 14 = digis10 - 1 allows short numbers to be
  // represended with trailing zeros.  This isn't necessarily the case with
  // precision = digits10, e.g., 8.3 becomes 8.300000000000001e+00
  //
  // This prescription doesn't exactly implement the method proposed.  If the
  // published table of numbers includes 8.300000000000001, this will be
  // interpreted as 8.3.

  while (getline(cin, s)) {
    double x = Utility::val<double>(s);
    Math::real X = Utility::val<Math::real>(s);
    ostringstream ostr;
    ostr << scientific << setprecision(prec) << x;
    double y = Utility::val<double>(ostr.str());
    Math::real Y = Utility::val<Math::real>(ostr.str());
    if (X != Y || x != y) cout << s << " " << ostr.str() << "\n";
  }
}

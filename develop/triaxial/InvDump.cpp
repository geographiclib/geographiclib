#include <iostream>
#include <iomanip>
#include <GeographicLib/Utility.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"

using namespace GeographicLib;
using namespace std;

/*
Run with
1 1.5 0.3333333 0.6666667 0          0         90
1 1.5 0.3333333 0.6666667 42.70330   0         90
1 1.5 0.3333333 0.6666667 87.52250   0         90
1 1.5 0.3333333 0.6666667 90         0        135
1 1.5 0.3333333 0.6666667 90        10.15216  180
1 1.5 0.3333333 0.6666667 90        39.25531  180
1 1.5 0.3333333 0.6666667 90        90        180
1 0.75 1 0 
1 3    0 1

*/

int main() {
  typedef Math::real real;
  //  cout << "Enter b e2 k2 kp2 bet omg alp\n";
  real b, e2, k2, kp2, bet, omg, alp;
  cin >> b >> e2 >> k2 >> kp2 >> bet >> omg >> alp;
  Triaxial t(b, e2, k2, kp2);
  TriaxialLine l(t, bet, omg, alp);
  l.Optimize();
  cout << setprecision(8);
  l.inversedump(cout);
}

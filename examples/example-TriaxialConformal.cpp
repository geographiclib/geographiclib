// Example of using the GeographicLib::JacobiConformal class.

#include <iostream>
#include <iomanip>
#include <exception>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/TriaxialConformal.hpp>
#include <GeographicLib/TriaxialGeodesic.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    typedef Math::real real;
    Utility::set_digits();
    // These parameters were derived from the EGM2008 geoid; see 2011-07-04
    // E-mail to PROJ.4 list, "Analyzing the bumps in the EGM2008 geoid".  The
    // longitude of the major axis is -15.  These are close to the values given
    // by Milan Bursa, Vladimira Fialova, "Parameters of the Earth's tri-axial
    // level ellipsoid", Studia Geophysica et Geodaetica 37(1), 1-13 (1993):
    //
    //    longitude of major axis = -14.93 +/- 0.05
    //    a = 6378171.36 +/- 0.30
    //    a/(a-c) = 297.7738 +/- 0.0003
    //    a/(a-b) = 91449 +/- 60
    // which gives: a = 6378171.36, b = 6378101.61, c = 6356751.84
    Math::real a = 6378137+35, b = 6378137-35, c = 6356752;
    TriaxialConformal proj(a, b, c);
    if (1) {
      Triaxial t(2e6, 1e6, 0.5e6);
      t = proj.t();
      TriaxialConformal p2(t);
      TriaxialGeodesic g2(t);
      real bet1 = 85, omg1 = 5, s12 = 0.1, x1, y1, k1;
      p2.Forward(bet1, omg1, x1, y1, k1);
      cout << "SCALE " << k1 << "\n";
      for (int i = -180; i < 180; i +=10) {
        real alp1 = real(i), bet2, omg2, alp2;
        g2.Direct(bet1, omg1, alp1, s12, bet2, omg2, alp2);
        real x2, y2;
        p2.Forward(bet2, omg2, x2, y2);
        x2 -= x1; y2 -= y1;
        cout << fixed << setprecision(5);
        real s12a = hypot(x2, y2) / k1,
          alp1a = Math::atan2d(x2, y2);
        cout << i << " " << alp1a << " " << s12a << "\n";
      }
      return 0;
    }
    cout  << fixed << setprecision(1)
          << "Ellipsoid parameters: a = "
          << a << ", b = " << b << ", c = " << c << "\n"
          << setprecision(3)
          << "Quadrants: x = " << proj.x() << ", y = " << proj.y() << "\n";
    cout << "Coordinates angle (deg) x (m) y (m):\n";
    for (int i = -540; i <= 540; i += 15) {
      Math::real omg = i, bet = i;
      Math::real x = proj.x(Angle(omg)), y = proj.y(Angle(bet));
      Angle omg1 = proj.omega(x), bet1 = proj.beta(y);
      cout << i << " " << x << " " << y << " "
           << Math::real(omg1) << " " << Math::real(bet1) << "\n";

    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <utility>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/Utility.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"
#include "Angle.hpp"
#include "TriaxialODE.hpp"

using namespace GeographicLib;
using namespace std;

int usage(int retval, bool /*brief*/) {
  ( retval ? cerr : cout ) << "Bad input\n";
  return retval;
}

string nicestr(Math::real x, int prec, bool azi = false) {
  typedef Math::real real;
  static const real eps = pow(real(10), -20);
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

void report(const Triaxial& t, int bet1, int omg1, int bet2, int omg2) {
#if GEOGRAPHICLIB_PRECISION <= 2
  int prec = 12;
#elif GEOGRAPHICLIB_PRECISION == 3
  int prec = 14;
#else
  int prec = 18;
#endif
  typedef Math::real real;
  typedef Angle ang;
  bool odep = false;
  ang bet1x(bet1), omg1x(omg1), bet2x(bet2), omg2x(omg2);
  real s12;
  ang alp1, alp2;
  TriaxialLine l =
    t.Inverse(bet1x, omg1x, bet2x, omg2x, alp1, alp2, s12);
  Triaxial::vec3 r2, v2;
  real m12 = 0, M12 = 1, M21 = 1;
  if (odep) {
    // ang bet1a, omg1a;
    // l.pos1(bet1a, omg1a, alp1);
    TriaxialODE direct(t, bet1, omg1, real(alp1));
    direct.Position(s12, r2, v2, m12, M12, M21);
    // t.cart2toellip(bet2x, omg2x, v2, alp2);
  }
  cout << bet1 << " " << omg1 << " "
       << nicestr(real(alp1), prec, true) << " "
       << bet2 << " " << omg2 << " "
       << nicestr(real(alp2), prec, true) << " "
       << nicestr(s12, prec+2);
  if (odep)
    cout << " " << nicestr(m12, prec+2) << " "
         << nicestr(M12, prec+2) << " " << nicestr(M21, prec+2) << endl;
  else
    cout << endl;
}

Math::real vecdiff(const Triaxial::vec3& a, const Triaxial::vec3& b) {
  return Math::hypot3(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

void print3(Triaxial::vec3 v) {
  cout << v[0] << " " << v[1] << " " << v[2] << "\n";
}

void errreport(const Triaxial& t,
               Math::real bet1, Math::real omg1, Math::real alp1,
               Math::real bet2, Math::real omg2, Math::real alp2,
               Math::real s12,
               Math::real /*m12*/, Math::real /*M12*/, Math::real /*M21*/) {
  typedef Math::real real;
  typedef Angle ang;
  bool debug = false;
#if GEOGRAPHICLIB_PRECISION > 3
  static const real eps = real(1e-20);
#else
  static const real eps = numeric_limits<real>::epsilon()/2;
#endif
  ang
    bet1x(bet1), omg1x(omg1), alp1x(alp1),
    bet2x(bet2), omg2x(omg2), alp2x(alp2), alp1a, alp2a;
  real s12a, errs = 0;
  TriaxialLine l0 =
    t.Inverse(bet1x, omg1x, bet2x, omg2x, alp1a, alp2a, s12a);
  errs = fabs(s12 - s12a);

  Triaxial::vec3 r1, v1, r2, v2;
  Triaxial::vec3 r1a, v1a, r2a, v2a;
  ang bet1a, omg1a, bet2a, omg2a;
  // direct checks for inverse calculation using alp1a, alp2a, s12a
  t.elliptocart2(bet1x, omg1x, alp1a, r1, v1);
  t.elliptocart2(bet2x, omg2x, alp2a, r2, v2);
  if (debug) {
    cout << "2X "
         << real(bet2x) << " " << real(omg2x) << " " << real(alp2a) << "\n";
    print3(r2); print3(v2);
  }
  TriaxialLine l1i(t, bet1x, omg1x, alp1a);
  TriaxialLine l2i(t, bet2x, omg2x, alp2a);
  real errr1i = 0, errv1i = 0,
    errr2i = 0, errv2i = 0;
  l2i.Position(-s12a, bet1a, omg1a, alp1a);
  t.elliptocart2(bet1a, omg1a, alp1a, r1a, v1a);
  errr1i = vecdiff(r1, r1a); errv1i = vecdiff(v1, v1a);
  l1i.Position(s12a, bet2a, omg2a, alp2a);
  if (debug)
    cout << "1X "
         << real(bet1x) << " " << real(omg1x) << " " << real(alp1a) << "\n";
  t.elliptocart2(bet2a, omg2a, alp2a, r2a, v2a);
  if (debug) {
    cout << "2A "
         << real(bet2a) << " " << real(omg2a) << " " << real(alp2a) << "\n";
    print3(r2a); print3(v2a);
  }
  errr2i = vecdiff(r2, r2a); errv2i = vecdiff(v2, v2a);

  // direct checks for test sets using alp1x, alp2x, s12
  t.elliptocart2(bet1x, omg1x, alp1x, r1, v1);
  t.elliptocart2(bet2x, omg2x, alp2x, r2, v2);
  TriaxialLine l1(t, bet1x, omg1x, alp1x);
  TriaxialLine l2(t, bet2x, omg2x, alp2x);
  real errr1 = 0, errv1 = 0,
    errr2 = 0, errv2 = 0;
  l2.Position(-s12, bet1a, omg1a, alp1a);
  t.elliptocart2(bet1a, omg1a, alp1a, r1a, v1a);
  errr1 = vecdiff(r1, r1a); errv1 = vecdiff(v1, v1a);
  l1.Position(s12, bet2a, omg2a, alp2a);
  t.elliptocart2(bet2a, omg2a, alp2a, r2a, v2a);
  errr2 = vecdiff(r2, r2a); errv2 = vecdiff(v2, v2a);

  cout << fixed << setprecision(0)
       << ceil(errs/eps) << " "
       << ceil(errr2i/eps) << " " << ceil(errv2i/eps) << " "
       << ceil(errr1i/eps) << " " << ceil(errv1i/eps) << " "
       << ceil(errr2/eps) << " " << ceil(errv2/eps) << " "
       << ceil(errr1/eps) << " " << ceil(errv1/eps) << endl;
}

void angletest() {
  typedef Math::real real;
  typedef Angle ang;
  ang o1 = ang::cardinal(2);
  ang o2 = ang::cardinal(-2);
  cout << o1.s() << " " << o2.s() << " "
       << real(o1) << " " << real(o2) << "\n";
  ang a1(-180);
  ang a2(-180);
  ang a3 = a2 - a1;
  cout << real(a1) << " " << real(a2) << " " << real(a3) << "\n";
  for (int i = -360; i <= 720; i += 30) {
    ang a1(i);
    for (int j = -360; j <= 720; j += 30) {
      ang a2(j);
      ang a3 = a2 - a1;
      int q = int(round(real(a3)));
      if (q != j - i)
        cout << "ERROR " << q - (j-i) << " " << q << " "
             << j << " " << i << "\n";
    }
  }
}

void dfinvtest() {
  Triaxial t(sqrt(real(2)), 1, sqrt(1/real(2)));
  TriaxialLine l(t, Angle(90), Angle(0), Angle(45));
  l.Optimize();
}

int main(int argc, const char* const argv[]) {
  try {
    typedef Math::real real;
    Utility::set_digits();
    if (0) {
      dfinvtest();
      return 0;
    }
    if (0) {
      angletest();
      return 0;
    }
    int num = 1000;
    int skew = 10;
    int div = 1;
    {
      real a = 1, b = 1, c = 1, e2 = -1, k2 = -1, kp2 = -1;
      for (int m = 1; m < argc; ++m) {
        string arg(argv[m]);
        if (arg == "-t") {
          if (m + 3 >= argc) return usage(1, true);
          try {
            a = Utility::val<real>(string(argv[m + 1]));
            b = Utility::val<real>(string(argv[m + 2]));
            c = Utility::val<real>(string(argv[m + 2]));
          }
          catch (const exception& e) {
            cerr << "Error decoding arguments of -e: " << e.what() << "\n";
            return 1;
          }
          m += 2;
        } else if (arg == "-e") {
          if (m + 4 >= argc) return usage(1, true);
          try {
            b = Utility::val<real>(string(argv[m + 1]));
            e2 = Utility::fract<real>(string(argv[m + 2]));
            k2 = Utility::fract<real>(string(argv[m + 3]));
            kp2 = Utility::fract<real>(string(argv[m + 4]));
          }
          catch (const exception& e) {
            cerr << "Error decoding arguments of -e: " << e.what() << "\n";
            return 1;
          }
          m += 4;
        } else
          return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
      }
      // testobl.txt  -e 1 3/4 3 0
      // testsetb.txt -e 1 1   2 1
      // testset.txt  -e 1 3/2 1 2
      // testpro.txt  -e 1 3   0 3
      // testspha.txt -e 1 0   3 0
      // testsphb.txt -e 1 0   2 1
      // testsphc.txt -e 1 0   1 2
      // testsphd.txt -e 1 0   0 3
      // testhu.txt   -e 1 7577279780/1130005142289 1942065235 6378137
      // equiv to     -t 6378172/6378102 1 6356752/6378102
      Triaxial t = e2 < 0 ? Triaxial(a, b, c) : Triaxial(b, e2, k2, kp2);
      // Triaxial t(1, 1, 1/real(2));
      // Triaxial t(2, 1, 1);
      int bet1, omg1, bet2, omg2;
      real alp1, alp2, s12, m12, M12, M21;
      while (cin >> bet1 >> omg1 >> alp1 >> bet2 >> omg2 >> alp2 >> s12
             >> m12 >> M12 >> M21) {
        if (0)
          report(t, bet1, omg1, bet2, omg2);
        else
          errreport(t, bet1, omg1, alp1, bet2, omg2, alp2, s12, m12, M12, M21);
      }
      return 0;
    }
    for (int m = 1; m < argc; ++m) {
      string arg(argv[m]);
      if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        try {
          num = Utility::val<int>(string(argv[m]));
        }
        catch (const exception&) {
          cerr << "num " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-d") {
        if (++m == argc) return usage(1, true);
        try {
          div = Utility::val<int>(string(argv[m]));
        }
        catch (const exception&) {
          cerr << "div " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-k") {
        if (++m == argc) return usage(1, true);
        try {
          skew = Utility::val<int>(string(argv[m]));
        }
        catch (const exception&) {
          cerr << "skew " << argv[m] << " is not a number\n";
          return 1;
        }
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
    }
    if (1) {
      Triaxial t(sqrt(real(2)), 1, sqrt(1/real(2)));
      for (int bet1 = -90; bet1 <= 90; bet1 += div)
        for (int omg1 = -180+div; omg1 <= 180; omg1 += div)
          for (int bet2 = -90; bet2 <= 90; bet2 += div)
            for (int omg2 = -180+div; omg2 <= 180; omg2 += div) {
              /*
              bool umb1 = fabs(bet1) > 85 && fabs( fabs(omg1) - 90 ) > 87,
                umb2 = fabs(bet2) > 85 && fabs( fabs(omg2) - 90 ) > 87;
              if (!(umb1 ^ umb2)) continue;
              if (rnd(gen) != 0) continue;
              */
              if (bet1 == 90 && omg1 == 125 &&
                  bet2 == 90 && omg2 == -140)
                report(t, bet1, omg1, bet2, omg2);
            }
    }
    div = div + skew + num;
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

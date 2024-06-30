#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <utility>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/Utility.hpp>
#include "Angle.hpp"
typedef GeographicLib::Angle AuxAngle;
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"
#include "Angle.hpp"
#include "TriaxialODE.hpp"

using namespace GeographicLib;
using namespace std;

int usage(int retval, bool /*brief*/) {
  ( retval ? std::cerr : std::cout ) << "Bad input\n";
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
  AuxAngle bet1x = AuxAngle::degrees(bet1), omg1x = AuxAngle::degrees(omg1), 
    bet2x = AuxAngle::degrees(bet2), omg2x = AuxAngle::degrees(omg2);
  TriaxialLine l =
    t.Inverse(bet1x, omg1x, bet2x, omg2x);
  real s12 = l.Distance();
  AuxAngle bet1a, omg1a, alp1, bet2a, omg2a, alp2;
  Triaxial::vec3 r2, v2;
  real m12, M12, M21;
  l.pos1(bet1a, omg1a, alp1);
  TriaxialODE direct(t, bet1, omg1, alp1.degrees());
  direct.Position(s12, r2, v2, m12, M12, M21);
  t.cart2toellip(bet2x, omg2x, v2, alp2);
  cout << bet1 << " " << omg1 << " "
       << nicestr(alp1.degrees(), prec, true) << " "
       << bet2 << " " << omg2 << " "
       << nicestr(alp2.degrees(), prec, true) << " "
       << nicestr(s12, prec+2) << " " << nicestr(m12, prec+2) << " "
       << nicestr(M12, prec+2) << " " << nicestr(M21, prec+2) << endl;
}

Math::real vecdiff(const Triaxial::vec3& a, const Triaxial::vec3& b) {
  return Triaxial::hypot3(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

void errreport(const Triaxial& t,
               Math::real bet1, Math::real omg1, Math::real alp1,
               Math::real bet2, Math::real omg2, Math::real alp2,
               Math::real s12,
               Math::real /*m12*/, Math::real /*M12*/, Math::real /*M21*/) {
  typedef Math::real real;
#if GEOGRAPHICLIB_PRECISION > 3
  static const real eps = real(1e-20);
#else
  static const real eps = numeric_limits<real>::epsilon()/2;
#endif
  AuxAngle
    bet1x = AuxAngle::degrees(bet1), omg1x = AuxAngle::degrees(omg1),
    alp1x = AuxAngle::degrees(alp1),
    bet2x = AuxAngle::degrees(bet2), omg2x = AuxAngle::degrees(omg2),
    alp2x = AuxAngle::degrees(alp2);
  TriaxialLine l0 =
    t.Inverse(bet1x, omg1x, bet2x, omg2x);
  real s12a = l0.Distance(), errs = fabs(s12 - s12a);
  Triaxial::vec3 r1, v1, r2, v2;
  t.elliptocart2(bet1x, omg1x, alp1x, r1, v1);
  t.elliptocart2(bet2x, omg2x, alp2x, r2, v2);
  AuxAngle bet1a, omg1a, alp1a, bet2a, omg2a, alp2a;
  Triaxial::vec3 r1a, v1a, r2a, v2a;
  TriaxialLine l1(t, bet1x, omg1x, alp1x);
  TriaxialLine l2(t, bet2x, omg2x, alp2x);
  l2.Position(-s12, bet1a, omg1a, alp1a);
  t.elliptocart2(bet1a, omg1a, alp1a, r1a, v1a);
  real errr1 = vecdiff(r1, r1a), errv1 = vecdiff(v1, v1a);
  l1.Position(s12, bet2a, omg2a, alp2a);
  t.elliptocart2(bet2a, omg2a, alp2a, r2a, v2a);
  real errr2 = vecdiff(r2, r2a), errv2 = vecdiff(v2, v2a);
  cout << fixed << setprecision(0)
       << ceil(errs/eps) << " "
       << ceil(errr2/eps) << " " << ceil(errv2/eps) << " "
       << ceil(errr1/eps) << " " << ceil(errv1/eps) << endl;
}

void angletest() {
  typedef Math::real real;
  typedef Angle ang;
  ang o1 = ang::cardinal(2);
  ang o2 = ang::cardinal(-2);
  cout << o1.s() << " " << o2.s() << " " << real(o1) << " " << real(o2) << "\n";
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
        cout << "ERROR " << q - (j-i) << " " << q << " " << j << " " << i << "\n";
    }
  }
}

int main(int argc, const char* const argv[]) {
  try {
    typedef Math::real real;
    Utility::set_digits();
    if (0) {
      angletest();
      return 0;
    }
    int num = 1000;
    int skew = 10;
    int div = 1;
    Triaxial::vec3 r, v;
    if (0) {
      Triaxial t(sqrt(real(2)), 1, sqrt(1/real(2)));
      real z = 0;
      AuxAngle bet0 = AuxAngle::degrees(90),
        omg0 = AuxAngle::degrees(180),
        alp0 = AuxAngle::degrees(-90);
      for (int sb = 0; sb < 2; ++sb) {
        AuxAngle bet(bet0); bet.x() = sb ? -z : z;
        for (int so = 0; so < 2; ++so) {
          AuxAngle omg(omg0); omg.y() = so ? -z : z;
          for (int sa1 = 0; sa1 < 2; ++sa1) {
            AuxAngle alp(alp0); alp.x() = sa1 ? -z : z;
            for (int sa2 = 0; sa2 < 2; ++sa2) {
              alp.y() = sa2 ? -1 : 1;
                   t.elliptocart2(bet, omg, alp, r, v);
              cout << bet.x() << " " << bet.y() << " "
                   << omg.x() << " " << omg.y() << " "
                   << alp.x() << " " << alp.y() << " "
                   << r[0] << " " << r[1] << " " << r[2] << " "
                   << v[0] << " " << v[1] << " " << v[2] << "\n";
            }
          }
        }
      }
      // ry = +/-0 = cos(bet) * sin(omg)
      // vy = +/-0 = sin(alp) * cos(alp)
      cout << "======\n";
      Triaxial::vec3 r0(r), v0(v);
      r0[1] = v0[1] = z;
      for (int sr = 0; sr < 2; ++sr) {
        r = r0; r[1] = sr ? -z : z;
        for (int sv = 0; sv < 2; ++sv) {
          v = v0; v[1] = sv ? -z : z;
          AuxAngle bet, omg, alp;
          t.cart2toellip(r, v, bet, omg, alp);
          cout << r[1] << " "
               << v[1] << " "
               << bet.x() << " " << bet.y() << " "
               << omg.x() << " " << omg.y() << " "
               << alp.x() << " " << alp.y() << "\n";
        }
      }
      // cos(alp) = -0, sin(alp) = copysign(1, -vy)
      cout << "======\n";
      for (int sb = 0; sb < 2; ++sb) {
        AuxAngle bet(bet0); bet.x() = sb ? -z : z;
        for (int so = 0; so < 2; ++so) {
          AuxAngle omg(omg0); omg.y() = so ? -z : z;
          for (int sv = 0; sv < 2; ++sv) {
            v = v0; v[1] = sv ? -z : z;
            AuxAngle alp;
            t.cart2toellip(bet, omg, v, alp);
            cout << bet.x() << " " << bet.y() << " "
                 << omg.x() << " " << omg.y() << " "
                 << v[1] << " "
                 << alp.x() << " " << alp.y() << "\n";
          }
        }
      }
      // cos(alp) = -0, sin(alp) = copysign(1, -vy)
      return 0;
    }
    {
      Triaxial t(sqrt(real(2)), 1, sqrt(1/real(2)));
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
      std::string arg(argv[m]);
      if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        try {
          num = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "num " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-d") {
        if (++m == argc) return usage(1, true);
        try {
          div = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "div " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-k") {
        if (++m == argc) return usage(1, true);
        try {
          skew = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "skew " << argv[m] << " is not a number\n";
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
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
}

#include <iostream>
#include <iomanip>
#include <GeographicLib/EllipticFunction.hpp>

namespace GeographicLib {

  class JacobiConformal {
    Math::real _a, _b, _c;
    const EllipticFunction _ex0, _exa, _ey0, _eya;
    void norm(Math::real& x, Math::real& y) {
      Math::real z = Math::hypot(x, y);
      x /= z; y /= z;
    }
  public:
    JacobiConformal(Math::real a, Math::real b, Math::real c)
      : _a(a)
      , _b(b)
      , _c(c)
      , _ex0(+(a-b) * (a+b) / ((a-c) * (a+c)) * Math::sq(c/b),
             +(a-b) * (a+b) / ((a-c) * (a+c)),
             +(b-c) * (b+c) / ((a-c) * (a+c)) * Math::sq(a/b),
             +(b-c) * (b+c) / ((a-c) * (a+c)) )
      , _exa(+(a-b) * (a+b) / ((a-c) * (a+c)) * Math::sq(c/b),
             -(a-b) * (a+b) / Math::sq(b),
             +(b-c) * (b+c) / ((a-c) * (a+c)) * Math::sq(a/b),
             +Math::sq(a/b) )
      , _ey0(-(b-c) * (b+c) / ((a-b) * (a+b)) * Math::sq(a/c),
             -(b-c) * (b+c) / ((a-b) * (a+b)),
             +(a-c) * (a+c) / ((a-b) * (a+b)) * Math::sq(b/c),
             +(a-c) * (a+c) / ((a-b) * (a+b)))
      , _eya(-(b-c) * (b+c) / ((a-b) * (a+b)) * Math::sq(a/c),
             -(b-c) * (b+c) / Math::sq(c),
             +(a-c) * (a+c) / ((a-b) * (a+b)) * Math::sq(b/c),
             +Math::sq(b/c))
    {}
    Math::real x(Math::real somg, Math::real comg) {
      using std::sqrt;
      Math::real somg1 = somg, comg1 = sqrt(_ex0.alphap2()) * comg;
      norm(somg1, comg1);
      Math::real domg1 = _ex0.Delta(somg1, comg1);
      return 2*_b/sqrt((_a-_c) * (_a+_c)) * _ex0.G(somg1, comg1, domg1);
    }
    Math::real x(Math::real omg) {
      Math::real
        a = omg * Math::degree(),
        somg = abs(omg) == 180 ? 0 : sin(a),
        comg = abs(omg) ==  90 ? 0 : cos(a);
      return x(somg, comg);
    }
    Math::real y(Math::real sbet, Math::real cbet) {
      using std::sqrt;
      Math::real sbet1 = sbet, cbet1 = sqrt(_ey0.alphap2()) * cbet;
      norm(sbet1, cbet1);
      Math::real dbet1 = _ey0.Delta(sbet1, cbet1);
      return 2*_c/sqrt((_a-_b) * (_a+_b)) * _ey0.G(sbet1, cbet1, dbet1);
    }
    Math::real y(Math::real bet) {
      Math::real
        a = bet * Math::degree(),
        sbet = abs(bet) == 180 ? 0 : sin(a),
        cbet = abs(bet) ==  90 ? 0 : cos(a);
      return y(sbet, cbet);
    }

    Math::real x() {
      return 2*Math::sq(_a)/(_b*sqrt((_a-_c) * (_a+_c))) * _exa.Pi();
    }
    Math::real x2(Math::real somg, Math::real comg) {
      using std::sqrt; using std::atan;
      Math::real somg1 = somg, comg1 = sqrt(_ex0.alphap2()) * comg;
      norm(somg1, comg1);
      Math::real domg1 = _ex0.Delta(somg1, comg1);
      return 2*Math::sq(_a)/(_b*sqrt((_a-_c) * (_a+_c))) *
        _exa.Pi(somg1, comg1, domg1)
        - 2 * atan( (_a-_b)*(_a+_b)*somg1*comg1 /
                    sqrt( Math::sq(_a*somg1) * (_b-_c)*(_b+_c) +
                          Math::sq(_b*comg1) * (_a-_c)*(_a+_c) ) );
    }
    Math::real x2(Math::real omg) {
      Math::real
        a = omg * Math::degree(),
        somg = abs(omg) == 180 ? 0 : sin(a),
        comg = abs(omg) ==  90 ? 0 : cos(a);
      return x2(somg, comg);
    }
    Math::real y() {
      return  2*Math::sq(_b)/(_c*sqrt((_a-_b) * (_a+_b))) * _eya.Pi();
    }
    Math::real y2(Math::real sbet, Math::real cbet) {
      using std:: sqrt;
      Math::real sbet1 = sbet, cbet1 = sqrt(_ey0.alphap2()) * cbet;
      norm(sbet1, cbet1);
      Math::real dbet1 = _ey0.Delta(sbet1, cbet1);
      return 2*Math::sq(_b)/(_c*sqrt((_a-_b) * (_a+_b))) *
        _eya.Pi(sbet1, cbet1, dbet1)
        - 2*(Math::atanh((_b-_c)*(_b+_c)*sbet1*cbet1/
                         sqrt( Math::sq(_b*sbet1) * (_a-_c)*(_a+_c) +
                               Math::sq(_c*cbet1) * (_a-_b)*(_a+_b) )) *0+
             Math::asinh((_b-_c)*(_b+_c)*sbet*cbet/
                         sqrt( Math::sq(_b*sbet) * (_a-_b)*(_a+_b) +
                               Math::sq(_c*cbet) * (_a-_c)*(_a+_c) )));
    }
    Math::real y2(Math::real bet) {
      Math::real
        a = bet * Math::degree(),
        sbet = abs(bet) == 180 ? 0 : sin(a),
        cbet = abs(bet) ==  90 ? 0 : cos(a);
      return y2(sbet, cbet);
    }

  };
}
using namespace std;
using namespace GeographicLib;

int main() {
  Math::real a = 6378137+35, b = /*6378137-35*/6356752+70, c = 6356752;
  JacobiConformal jc(a, b, c);
  cout  << fixed << setprecision(30)
        << "Quadrants " << jc.x() << " " << jc.y() << "\n";
  for (int i = -180; i <= 180; i += 5) {
    if (i == -180) continue;
    Math::real omg = i, bet = i;
    cout //         << i << " " << jc.x(omg) << " " << jc.y(bet) << "\n"
         << i << " " << jc.x2(omg) << " " << jc.y2(bet) << "\n";
  }
}

#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Utility.hpp>

using namespace GeographicLib;
using namespace std;

Math::real area0(Math::real a, Math::real b, Math::real c) {
  return (a*b*c)*
  EllipticFunction::RG(1/Math::sq(a),1/Math::sq(b),1/Math::sq(c));
}

Math::real area1(Math::real a, Math::real b, Math::real c) {
  Math::real p = 1.6075; p = 1.6;
  return pow( (pow(a*b, p) + pow(b*c, p) + pow(c*a, p)) / 3, 1/p );
}

class RandLoc {
private:
  typedef Math::real real;
  real _a, _b, _c, _d;
  seed_seq _seq;
  mt19937 _gen;
  uniform_real_distribution<double> _rnd;
  normal_distribution<double> _norm;
public:
  RandLoc(Math::real a, Math::real b, Math::real c,
          unsigned seed1, unsigned seed2)
    : _a(a)
    , _b(b)
    , _c(c)
    , _seq({seed1, seed2})
    , _gen(_seq)
    , _rnd(0.0, 1.0)
    , _norm(0.0, 1.0)
  {
    _d = fmin(_a, fmin(_b, _c));
  }
  // Return random position as parametric latitude and longitude
  void RandPos(Math::real& beta, Math::real& omega) {
    real Xn, Yn, Zn, g;
    do {
      Xn = _norm(_gen);
      Yn = _norm(_gen);
      Zn = _norm(_gen);
      real R = hypot(Xn, hypot(Yn, Zn));
      Xn /= R; Yn /= R; Zn /= R;
      real Xu = Xn / _a, Yu = Yn / _b, Zu = Zn / _c;
      g = _d * hypot(Xu, hypot(Yu, Zu));
    } while (_rnd(_gen) > g);
    omega = Math::atan2d(Yn, Xn);
    beta = Math::atan2d(Zn, hypot(Xn, Yn));
  }
  Math::real RandAzi() {
    return (2*real(_rnd(_gen)) - 1) * real(Math::hd);
  }
};

string trim0(string s) {
  // Trims fixed point number of trailing 0s + decimal point
  string::size_type r = s.find('.');
  int l = int(s.length());
  if (r != string::npos) {
    int p = int(r);
    if (l == 0) return("0");
    int i = l - 1;
    for (; i > p; --i) {
      if (s[i] != '0') { break; }
    }
    // i points to last non 0
    if (i == p)
      s = s.substr(0, i);         // Trim off period
    else
      s = s.substr(0, i+1);
  }
  l = int(s.length());
  if (l) {
    if (s[0] == '+') { s = s.substr(1); --l; }
    if (s[0] == '0') { s = s.substr(1); --l; }
    else if (l > 1 && s[0] == '-' && s[1] == '0') {
      s = s.substr(1); s[0] = '-'; --l;
    }
  }
  if (l == 0) s = "0";
  return s;
}

string fmt(Math::real x, int prec) {
  return trim0(Utility::str(x, prec));
}

int usage(int retval, bool /*brief*/) {
  ( retval ? cerr : cout ) << "Bad input\n";
  return retval;
}

void rotate(Math::real& beta, Math::real& omega, Math::real& azi, int dir) {
  //  cerr << beta << " " << omega << " " << azi << endl;
  Math::real omegax = omega;
  if (dir > 0)  swap(beta,omega);
  beta = remainder(beta-dir*90, Math::real(360));
  omega = remainder(omega+dir*90, Math::real(360));
  azi = remainder(90-azi, Math::real(360));
  if (dir < 0)  swap(beta,omega);
  if (fabs(beta) > 90 ||
      (fabs(beta) == 90 &&
       (omegax == -180 || (omegax == 0 && signbit(omegax))))) {
    beta = remainder(180-beta, Math::real(360));
    omega = -omega;
    //if (omega == -180) omega = 180;
    azi = remainder(180+azi, Math::real(360));
  }
  // cerr << beta << " " << omega << " " << azi << endl;
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

int main(int argc, const char* const argv[]) {
  try {
    typedef Math::real real;
    Utility::set_digits();
    if (0) {
      int num = 100;
      for (int i = 0; i <= num; ++i) {
        real a = pow(10.0, 2*i/real(num));
        a = 1 + 0.01 * (i-num)/real(num);
        for (int j = i; j <= num; ++j) {
          real b = pow(10.0, 2*j/real(num));
          b = 1 + 0.01 * (j-num)/real(num);
          for (int k = j; k <= num; ++k) {
            real c = pow(10.0, 2*k/real(num));
            c = 1 + 0.01 * (k-num)/real(num);
            real a0 = area0(a, b, c);
            real a1 = area1(a, b, c);
            cout << a << " " << b << " " << c << " "
                      << (a1 - a0) / a0 << "\n";
          }
        }
      }
      return 0;
    }
    real a = 1, f = 0;
    unsigned seed = 0;
    int prec0 = 0,              // Precisions relative to 1m
      prec1 = 9,
      num = 1000;
    for (int m = 1; m < argc; ++m) {
      string arg(argv[m]);
      if (arg == "-e") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          a = Utility::val<real>(string(argv[m + 1]));
          f = Utility::fract<real>(string(argv[m + 2]));
        }
        catch (const exception& e) {
          cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        m += 2;
      } else if (arg == "-s") {
        if (++m == argc) return usage(1, true);
        try {
          seed = Utility::val<unsigned>(string(argv[m]));
        }
        catch (const exception&) {
          cerr << "seed " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        try {
          num = Utility::val<int>(string(argv[m]));
        }
        catch (const exception&) {
          cerr << "num " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-p0") {
        if (++m == argc) return usage(1, true);
        try {
          prec0 = Utility::val<int>(string(argv[m]));
        }
        catch (const exception&) {
          cerr << "prec0 " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-p1") {
        if (++m == argc) return usage(1, true);
        try {
          prec1 = Utility::val<int>(string(argv[m]));
        }
        catch (const exception&) {
          cerr << "prec1 " << argv[m] << " is not a number\n";
          return 1;
        }
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
    }
    if (seed == 0)
      seed = random_device()();
    cerr << "-s " << seed << "\n";
    if (num + prec0 + prec1 + seed == 0)
      cout << "foo\n";
    RandLoc rnd(a, a, a*(1 - f), seed, 0);
    Geodesic geod(a, f, true);
    int preci = prec0 + 5;
    real rndfac = pow(real(10), preci);
    if (0) {
      for (int i = 0; i < num; ++i) {
        real beta1, omega1;
        rnd.RandPos(beta1, omega1);
        beta1 = round(beta1 * rndfac) / rndfac;
        omega1 = round(omega1 * rndfac) / rndfac;
        if (omega1 == -180) omega1 = 180;
        real phi1 = Math::atand(Math::tand(beta1) / (1 - f));
        real beta2, omega2;
        rnd.RandPos(beta2, omega2);
        beta2 = round(beta2 * rndfac) / rndfac;
        omega2 = round(omega2 * rndfac) / rndfac;
        if (omega2 == -180) omega2 = 180;
        real phi2 = Math::atand(Math::tand(beta2) / (1 - f));
        real azi1, azi2, s12, m12;
        geod.Inverse(phi1, omega1, phi2, omega2, s12, azi1, azi2, m12);
        if (f < 0) {
          swap(beta1,omega1);
          beta1 = remainder(beta1-90, 360);
          omega1 = remainder(omega1+90, 360);
          azi1 = remainder(90-azi1, 360);
          if (fabs(beta1) > 90) {
            beta1 = remainder(180-beta1, 360);
            omega1 = -omega1;
            if (omega1 == -180) omega1 = 180;
            azi1 = remainder(180+azi1, 360);
          }
          swap(beta2,omega2);
          beta2 = remainder(beta2-90, 360);
          omega2 = remainder(omega2+90, 360);
          azi2 = remainder(90-azi2, 360);
          if (fabs(beta2) > 90) {
            beta2 = remainder(180-beta2, 360);
            omega2 = -omega2;
            if (omega2 == -180) omega2 = 180;
            azi2 = remainder(180+azi2, 360);
          }
        }
        cout << fmt(beta1, preci) << " "
                  << fmt(omega1, preci) << " "
                  << fmt(azi1, prec1 + 5) << " "
                  << fmt(beta2, preci) << " "
                  << fmt(omega2, preci) << " "
                  << fmt(azi2, prec1 + 5) << " "
                  << fmt(s12, prec1 + 7)<< " "
                  << fmt(m12, prec1 + 7) << "\n";
      }
    } else {
      real bet1, omg1, bet2, omg2;
      real alp1, alp2, s12, m12, M12, M21;
      real phi1, phi2;
#if GEOGRAPHICLIB_PRECISION <= 2
      int prec = 12;
#elif GEOGRAPHICLIB_PRECISION == 3
      int prec = 14;
#else
      int prec = 18;
#endif
      while (cin >> bet1 >> omg1 >> alp1 >> bet2 >> omg2 >> alp2 >> s12
             >> m12 >> M12 >> M21) {
        if (f < 0) {
          rotate(bet1, omg1, alp1, 1);
          rotate(bet2, omg2, alp2, 1);
        }
        phi1 = Math::atand(Math::tand(bet1) / (1 - f));
        phi2 = Math::atand(Math::tand(bet2) / (1 - f));
        geod.Inverse(phi1, omg1, phi2, omg2,
                     s12, alp1, alp2, m12, M12, M21);
        if (f < 0) {
          rotate(bet1, omg1, alp1, -1);
          rotate(bet2, omg2, alp2, -1);
        }
        if (false && f < 0 &&
            bet1 + bet2 == 0 && Math::cosd(alp1) * Math::cosd(alp2) < 0) {
          // Attempt to standardize the azimuth when there's more than one
          // shortest geodesic.  SKIP FOR NOW.
          alp1 = 180 - alp1;
          alp2 = 180 - alp2;
        }
        cout << int(bet1) << " " << int(omg1) << " "
             << nicestr(real(alp1), prec, true) << " "
             << int(bet2) << " " << int(omg2) << " "
             << nicestr(real(alp2), prec, true) << " "
             << nicestr(s12, prec+2) << " " << nicestr(m12, prec+2) << " "
             << nicestr(M12, prec+2) << " " << nicestr(M21, prec+2) << endl;
      }
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

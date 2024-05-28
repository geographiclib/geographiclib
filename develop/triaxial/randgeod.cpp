#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Utility.hpp>

using namespace GeographicLib;

Math::real area0(Math::real a, Math::real b, Math::real c) {
  return (a*b*c)*
  EllipticFunction::RG(1/Math::sq(a),1/Math::sq(b),1/Math::sq(c));
}

Math::real area1(Math::real a, Math::real b, Math::real c) {
  using std::pow;
  Math::real p = 1.6075; p = 1.6;
  return pow( (pow(a*b, p) + pow(b*c, p) + pow(c*a, p)) / 3, 1/p );
}

class RandLoc {
private:
  typedef Math::real real;
  real _a, _b, _c, _d;
  std::seed_seq _seq;
  std::mt19937 _gen;
  std::uniform_real_distribution<double> _rnd;
  std::normal_distribution<double> _norm;
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
    using std::fmin;
    _d = fmin(_a, fmin(_b, _c));
  }
  // Return random position as parametric latitude and longitude
  void RandPos(Math::real& beta, Math::real& omega) {
    using std::hypot; using std::pow;
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

std::string trim0(std::string s) {
  // Trims fixed point number of trailing 0s + decimal point
  std::string::size_type r = s.find('.');
  int l = int(s.length());
  if (r != std::string::npos) {
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

std::string fmt(Math::real x, int prec) {
  return trim0(Utility::str(x, prec));
}

int usage(int retval, bool /*brief*/) {
  ( retval ? std::cerr : std::cout ) << "Bad input\n";
  return retval;
}

int main(int argc, const char* const argv[]) {
  try {
    using std::swap;
    using std::remainder;
    using std::round;
    using std::fabs;
    typedef Math::real real;
    Utility::set_digits();
    {
      using std::pow;
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
            std::cout << a << " " << b << " " << c << " "
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
      std::string arg(argv[m]);
      if (arg == "-e") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          a = Utility::val<real>(std::string(argv[m + 1]));
          f = Utility::fract<real>(std::string(argv[m + 2]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        m += 2;
      } else if (arg == "-s") {
        if (++m == argc) return usage(1, true);
        try {
          seed = Utility::val<unsigned>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "seed " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        try {
          num = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "num " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-p0") {
        if (++m == argc) return usage(1, true);
        try {
          prec0 = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "prec0 " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-p1") {
        if (++m == argc) return usage(1, true);
        try {
          prec1 = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "prec1 " << argv[m] << " is not a number\n";
          return 1;
        }
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
    }
    if (seed == 0)
      seed = std::random_device()();
    std::cerr << "-s " << seed << "\n";
    if (num + prec0 + prec1 + seed == 0)
      std::cout << "foo\n";
    RandLoc rnd(a, a, a*(1 - f), seed, 0);
    Geodesic geod(a, f, true);
    int preci = prec0 + 5;
    real rndfac = pow(real(10), preci);
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
      std::cout << fmt(beta1, preci) << " "
                << fmt(omega1, preci) << " "
                << fmt(azi1, prec1 + 5) << " "
                << fmt(beta2, preci) << " "
                << fmt(omega2, preci) << " "
                << fmt(azi2, prec1 + 5) << " "
                << fmt(s12, prec1 + 7)<< " "
                << fmt(m12, prec1 + 7) << "\n";
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

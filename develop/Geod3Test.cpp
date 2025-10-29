#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <utility>
#include <vector>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Angle.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>
#include <GeographicLib/Triaxial/Geodesic3.hpp>
#include <GeographicLib/Triaxial/Cartesian3.hpp>

#define DIAGS 0
#define ODEMOD 0

#if !defined(HAVE_BOOST)
#  define HAVE_BOOST 0
#endif

#if HAVE_BOOST
#if defined(_MSC_VER)
// Squelch warning triggered by boost:
//   4127: conditional expression is constant
#  pragma warning (disable: 4127)
#endif
#include "../experimental/TriaxialGeodesicODE.hpp"
#endif

using namespace std;
using namespace GeographicLib;
using namespace Triaxial;
#if HAVE_BOOST
using experimental::TriaxialGeodesicODE;
#endif

int usage(int retval, bool /*brief*/) {
  ( retval ? cerr : cout ) << "Bad input\n";
  return retval;
}

string nicestr(Math::real x, int prec, bool azi = false) {
  using real = Math::real;
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

void report(const Geodesic3& tg, int bet1, int omg1, int bet2, int omg2,
            bool odep) {
#if GEOGRAPHICLIB_PRECISION <= 2
  int prec = 12;
#elif GEOGRAPHICLIB_PRECISION == 3
  int prec = 14;
#else
  int prec = 18;
#endif
  using real = Math::real;
  using ang = Angle;
  ang bet1x(bet1), omg1x(omg1), bet2x(bet2), omg2x(omg2);
  real s12;
  ang alp1, alp2;
  GeodesicLine3 l =
    tg.Inverse(bet1x, omg1x, bet2x, omg2x, s12, alp1, alp2);
  real m12 = 0, M12 = 1, M21 = 1;
  if (odep) {
#if HAVE_BOOST
    Ellipsoid3::vec3 r2, v2;
    // ang bet1a, omg1a;
    // l.pos1(bet1a, omg1a, alp1);
    TriaxialGeodesicODE direct(tg.t(), ang(bet1), ang(omg1), alp1, true, true);
    direct.Position(s12, r2, v2, m12, M12, M21);
    // t.cart2toellip(bet2x, omg2x, v2, alp2);
#endif
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

Math::real vecdiff(const Ellipsoid3::vec3& a, const Ellipsoid3::vec3& b) {
  return Math::hypot3(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

void print3(Ellipsoid3::vec3 v) {
  cout << v[0] << " " << v[1] << " " << v[2] << "\n";
}

void errreport(const Geodesic3& tg,
               Math::real bet1, Math::real omg1, Math::real alp1,
               Math::real bet2, Math::real omg2, Math::real alp2,
               Math::real s12,
               Math::real /*m12*/, Math::real /*M12*/, Math::real /*M21*/) {
  using real = Math::real;
  using ang = Angle;
  bool debug = false, invp = true, invdirp = true,
    dirp = true, swapp = false, dirdiag = false;
  // invp = false;
  // invdirp = false;
  // dirp = false;
  // swapp = true;
  int num = 0;
  if (swapp) swap(omg1, omg2);
#if GEOGRAPHICLIB_PRECISION > 3
  static const real eps = real(1e-20);
#else
  static const real eps = numeric_limits<real>::epsilon()/2;
#endif
  ang
    bet1x(bet1), omg1x(omg1), alp1x(alp1),
    bet2x(bet2), omg2x(omg2), alp2x(alp2), alp1a, alp2a;
  real errs = 0, errr1i = 0, errv1i = 0, errr2i = 0, errv2i = 0,
    errr1 = 0, errv1 = 0, errr2 = 0, errv2 = 0;
  Ellipsoid3::vec3 r1, v1, r2, v2;
  Ellipsoid3::vec3 r1a, v1a, r2a, v2a;
  ang bet1a, omg1a, bet2a, omg2a;
  int invcountn = 0, invcountb = 0;
  if (invp) {
    real s12a;
#if DIAG
    GeodesicLine3 l0 =
      tg.Inverse(bet1x, omg1x, bet2x, omg2x, s12a, alp1a, alp2a,
                 &invcountn, &invcountb);
#else
    GeodesicLine3 l0 =
      tg.Inverse(bet1x, omg1x, bet2x, omg2x, s12a, alp1a, alp2a);
    (void) invcountn;
    (void) invcountb;
#endif
    errs = fabs(s12 - s12a);
    // direct checks for inverse calculation using alp1a, alp2a, s12a
    tg.t().elliptocart2(bet1x, omg1x, alp1a, r1, v1);
    tg.t().elliptocart2(bet2x, omg2x, alp2a, r2, v2);
    if (debug) {
      cout << "2X "
           << real(bet2x) << " " << real(omg2x) << " " << real(alp2a) << "\n";
      print3(r2); print3(v2);
    }
    if (invdirp) {
      GeodesicLine3 l1i(tg, bet1x, omg1x, alp1a);
      GeodesicLine3 l2i(tg, bet2x, omg2x, alp2a);
#if DIAG
      if (0) {
        int nfbet, nfomg, ngbet, ngomg;
        l1i.ncoeffs(nfbet, nfomg, ngbet, ngomg);
        cout <<  "LICOEFFS "
             << nfbet << " " << nfomg << " " << ngbet << " " << ngomg << "\n";
        l2i.ncoeffs(nfbet, nfomg, ngbet, ngomg);
        cout <<  "L2COEFFS "
             << nfbet << " " << nfomg << " " << ngbet << " " << ngomg << "\n";
      }
#endif
      l2i.Position(-s12a, bet1a, omg1a, alp1a);
      tg.t().elliptocart2(bet1a, omg1a, alp1a, r1a, v1a);
      errr1i = vecdiff(r1, r1a); errv1i = vecdiff(v1, v1a);
      l1i.Position(s12a, bet2a, omg2a, alp2a);
      if (debug)
        cout << "1X "
             << real(bet1x) << " " << real(omg1x) << " " << real(alp1a) << "\n";
      tg.t().elliptocart2(bet2a, omg2a, alp2a, r2a, v2a);
      if (debug) {
        cout << "2A "
             << real(bet2a) << " " << real(omg2a) << " " << real(alp2a) << "\n";
        print3(r2a); print3(v2a);
      }
      errr2i = vecdiff(r2, r2a); errv2i = vecdiff(v2, v2a);
    }
  }
  int nfbet, nfomg, ngbet, ngomg, countn = 0, countb = 0;
  if (dirp) {
    // direct checks for test sets using alp1x, alp2x, s12
    tg.t().elliptocart2(bet1x, omg1x, alp1x, r1, v1);
    tg.t().elliptocart2(bet2x, omg2x, alp2x, r2, v2);
    GeodesicLine3 l1(tg, bet1x, omg1x, alp1x);
    GeodesicLine3 l2(tg, bet2x, omg2x, alp2x);
    l2.Position(-s12, bet1a, omg1a, alp1a);
    tg.t().elliptocart2(bet1a, omg1a, alp1a, r1a, v1a);
    errr1 = vecdiff(r1, r1a); errv1 = vecdiff(v1, v1a);
#if DIAG
    l1.Position(s12, bet2a, omg2a, alp2a, &countn, &countb);
    l1.ncoeffs(nfbet, nfomg, ngbet, ngomg);
#else
    l1.Position(s12, bet2a, omg2a, alp2a);
    (void) nfbet;
    (void) nfomg;
    (void) ngbet;
    (void) ngomg;
    (void) countn;
    (void) countb;
#endif
    tg.t().elliptocart2(bet2a, omg2a, alp2a, r2a, v2a);
    errr2 = vecdiff(r2, r2a); errv2 = vecdiff(v2, v2a);
    real ds = real(10)/max(1,num);
    for (int i = 1; i < num; ++i) {
      // Extra 2 x (num - 1) calls
      l1.Position(i*ds, bet2a, omg2a, alp2a);
      l2.Position(-i*ds, bet1a, omg1a, alp1a);
    }
  }
  cout << fixed << setprecision(0)
       << ceil(errs/eps) << " "
       << ceil(errr2i/eps) << " " << ceil(errv2i/eps) << " "
       << ceil(errr1i/eps) << " " << ceil(errv1i/eps) << " "
       << ceil(errr2/eps) << " " << ceil(errv2/eps) << " "
       << ceil(errr1/eps) << " " << ceil(errv1/eps) << endl;
  if (dirdiag)
    cout << "DIAG "
         << nfbet << " " << nfomg << " " << ngbet << " " << ngomg << " "
         << countn  << " " << countb << " "
         << invcountn  << " " << invcountb << "\n";
}

#if HAVE_BOOST
void errODE(TriaxialGeodesicODE& l,
            Math::real bet1, Math::real omg1, Math::real alp1,
            Math::real bet2, Math::real omg2, Math::real alp2,
            Math::real s12,
            Math::real m12, Math::real M12, Math::real M21,
            bool steps) {
  using real = Math::real;
  using ang = Angle;
#if GEOGRAPHICLIB_PRECISION > 3
  static const real eps = real(1e-20);
#else
  static const real eps = numeric_limits<real>::epsilon()/2;
#endif
  ang
    bet1x(bet1), omg1x(omg1), alp1x(alp1),
    bet2x(bet2), omg2x(omg2), alp2x(alp2);
  s12 *= 1;
  Ellipsoid3::vec3 r1, v1, r2, v2;
  l.t().elliptocart2(bet1x, omg1x, alp1x, r1, v1);
  l.t().elliptocart2(bet2x, omg2x, alp2x, r2, v2);
  Ellipsoid3::vec3 r1a, v1a, r2a, v2a;
  real m12a = 0, M12a = 0, M21a = 0;
  if (1) {
    l.Reset(r1, v1);
    l.Position(s12, r2a, v2a, m12a, M12a, M21a);
    int n1 = l.NSteps(), i1 = l.IntSteps();
    l.Reset(r2, v2);
    l.Position(-s12, r1a, v1a);
    int n2 = l.NSteps(), i2 = l.IntSteps();
    real errr1 = vecdiff(r1, r1a), errv1 = vecdiff(v1, v1a),
      errr2 = vecdiff(r2, r2a), errv2 = vecdiff(v2, v2a),
      errm12 = fabs(m12a - m12),
      errM12 = fabs(M12a - M12), errM21 = fabs(M21a - M21);
    (void) n1; (void) i1; (void) n2; (void) i2;
    cout << fixed << setprecision(0)
         << ceil(errr2/eps) << " " << ceil(errv2/eps) << " "
         << ceil(errr1/eps) << " " << ceil(errv1/eps);
    if (l.Extended())
      cout << " " << ceil(errm12/eps) << " "
           << ceil(errM12/eps) << " "  << ceil(errM21/eps);
    cout << endl;
    if (steps)
      cout << "STEPS " << n1 << " " << n2 << " " << i1 << " " << i2 << "\n";
  } else {
    int num = 10;
    vector<real> s12v(num);
    vector<Ellipsoid3::vec3> rv, vv;
    for (int i = 1; i <= num; ++i)
      s12v[i-1] = i*s12/num;
    l.Reset(r1, v1);
    l.Position(s12v, rv, vv);
    r2a = rv[num-1]; v2a = vv[num-1];
    for (int i = 1; i <= num; ++i)
      s12v[i-1] = -i*s12/num;
    l.Reset(r2, v2);
    l.Position(s12v, rv, vv);
    r1a = rv[num-1]; v1a = vv[num-1];
    real errr1 = vecdiff(r1, r1a), errv1 = vecdiff(v1, v1a),
      errr2 = vecdiff(r2, r2a), errv2 = vecdiff(v2, v2a);
    cout << fixed << setprecision(0)
         << ceil(errr2/eps) << " " << ceil(errv2/eps) << " "
         << ceil(errr1/eps) << " " << ceil(errv1/eps);
    cout << endl;
  }
}
#endif

void angletest() {
  using real = Math::real;
  using ang = Angle;
  ang o1 = ang::cardinal(2);
  ang o2 = ang::cardinal(-2);
  cout << o1.s() << " " << o2.s() << " "
       << real(o1) << " " << real(o2) << "\n";
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

#if HAVE_BOOST
void msols(TriaxialGeodesicODE& l, Math::real s,
           Math::real& m12, Math::real& M12,
           Math::real& m12p, Math::real& M12p) {
  Ellipsoid3::vec3 R, V;
#if ODEMOD
  // Return real M12p
  l.Position(s, R, V, m12, M12, m12p, M12p);
#else
  // Compute M12p, ill-conditioned when m12 is small
  Math::real M21;
  l.Position(s, R, V, m12, M12, M21);
  m12p = M21; M12p = - (1 - M12 * M21) / m12;
#endif
}

void redlenreport(TriaxialGeodesicODE& l) {
  using real = Math::real;
  const Ellipsoid3& t = l.t();
  real a = t.a(), b = t.b(), c = t.c(), m12, M12, m12p, M12p;
  cout << fixed << setprecision(12);
  // Minor ellipse
  real s = 4 * EllipticFunction::RG(Math::sq(b), Math::sq(c));
  l.Reset({0, b, 0}, {0, 0, 1});
  msols(l, s, m12, M12, m12p, M12p);
  cout << M12 << " " << m12 << " " << M12p << " " << m12p << "\n";
  l.Reset({0, 0, c}, {0, 1, 0});
  msols(l, s, m12, M12, m12p, M12p);
  cout << M12 << " " << m12 << " " << M12p << " " << m12p << "\n";
  // Median ellipse
  s = 4 * EllipticFunction::RG(Math::sq(a), Math::sq(c));
  l.Reset({a, 0, 0}, {0, 0, 1});
  msols(l, s, m12, M12, m12p, M12p);
  cout << M12 << " " << m12 << " " << M12p << " " << m12p << "\n";
  l.Reset({0, 0, c}, {1, 0, 0});
  msols(l, s, m12, M12, m12p, M12p);
  cout << M12 << " " << m12 << " " << M12p << " " << m12p << "\n";
  Ellipsoid3::vec3 R0, V0;
  l.t().elliptocart2(Angle(90),Angle(0),Angle(0), R0, V0);
  l.Reset(R0, V0);
  msols(l, s, m12, M12, m12p, M12p);
  cout << M12 << " " << m12 << " " << M12p << " " << m12p << "\n";
  // Major ellipse
  s = 4 * EllipticFunction::RG(Math::sq(a), Math::sq(b));
  l.Reset({a, 0, 0}, {0, 1, 0});
  msols(l, s, m12, M12, m12p, M12p);
  cout << M12 << " " << m12 << " " << M12p << " " << m12p << "\n";
  l.Reset({0, b, 0}, {1, 0, 0});
  msols(l, s, m12, M12, m12p, M12p);
  cout << M12 << " " << m12 << " " << M12p << " " << m12p << "\n";
}
#endif

int main(int argc, const char* const argv[]) {
  try {
    using real = Math::real;
    Utility::set_digits();
    bool odep = false, reportp = false,
      odetest = false, extended = false, dense = false, normp = false,
      steps = false,
      inverse = false,          // just for timing
      redlen = false;
    real a = 1, b = 1, c = 1, e2 = -1, k2 = -1, kp2 = -1, eps = 0;
    int num = 0, numl = 1;
    unsigned long long seed = 0;
    for (int m = 1; m < argc; ++m) {
      string arg(argv[m]);
      if (arg == "-t") {
        if (m + 3 >= argc) return usage(1, true);
        try {
          a = Utility::val<real>(string(argv[m + 1]));
          b = Utility::val<real>(string(argv[m + 2]));
          c = Utility::val<real>(string(argv[m + 3]));
          e2 = -1;
        }
        catch (const exception& e) {
          cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        m += 3;
      } else if (arg == "-e") {
        if (m + 4 >= argc) return usage(1, true);
        try {
          b = Utility::val<real>(string(argv[m + 1]));
          e2 = Utility::fract<real>(string(argv[m + 2]));
          k2 = Utility::fract<real>(string(argv[m + 3]));
          kp2 = Utility::fract<real>(string(argv[m + 4]));
          a = -1;
        }
        catch (const exception& e) {
          cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        m += 4;
      } else if (arg == "--ode")
        odep = true;
      else if (arg == "--report")
        reportp = true;
      else if (arg == "--odetest")
        odetest = true;
      else if (arg == "--steps")
        steps = true;
      else if (arg == "-x")
        extended = true;
      else if (arg == "--redlen")
        redlen = extended = true;
      else if (arg == "--eps") {
        if (m + 1 >= argc) return usage(1, true);
        try {
          using std::pow;
          eps = pow(std::numeric_limits<real>::epsilon(),
                    Utility::fract<real>(std::string(argv[m + 1])));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of --eps: "
                    << e.what() << "\n";
          return 1;
        }
        m += 1;
      } else if (arg == "--dense")
        dense = true;
      else if (arg == "--normp")
        normp = true;
      else if (arg == "-i")
        inverse = true;
      else if (arg == "--num") {
        if (m + 1 >= argc) return usage(1, true);
        try {
          num = Utility::val<int>(string(argv[m + 1]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of --eps: "
                    << e.what() << "\n";
          return 1;
        }
        m += 1;
      } else if (arg == "--numl") {
        if (m + 1 >= argc) return usage(1, true);
        try {
          numl = Utility::val<int>(string(argv[m + 1]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of --eps: "
                    << e.what() << "\n";
          return 1;
        }
        m += 1;
      } else if (arg == "--seed") {
        if (++m == argc) return usage(1, true);
        try {
          seed = Utility::val<unsigned long long>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
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
    Geodesic3 tg = e2 < 0 ? Geodesic3(a, b, c) :
      Geodesic3(b, e2, k2, kp2);
    // Triaxial tg(1, 1, 1/real(2));
    // Triaxial tg(2, 1, 1);
    real bet1, omg1, bet2, omg2;
    real alp1, alp2, s12, m12, M12, M21;
#if HAVE_BOOST
    TriaxialGeodesicODE l(tg.t(), extended, dense, normp, eps);
#endif
    if (num <= 0) {
      if (redlen) {
#if HAVE_BOOST
        redlenreport(l);
#endif
        return 0;
      }
      while (cin >> bet1 >> omg1 >> alp1 >> bet2 >> omg2 >> alp2 >> s12
             >> m12 >> M12 >> M21) {
        if (odetest) {
#if HAVE_BOOST
          errODE(l, bet1, omg1, alp1, bet2, omg2, alp2, s12, m12, M12, M21,
                 steps);
#endif
          (void) steps;
        } else if (reportp)
          report(tg, int(bet1), int(omg1), int(bet2), int(omg2), odep);
        else
          errreport(tg, bet1, omg1, alp1, bet2, omg2, alp2, s12,
                    m12, M12, M21);
      }
    } else {
      unsigned long long s1 = seed;
      if (seed == 0)
        s1 = std::random_device()();
      std::seed_seq seq{s1};
      std::mt19937 g(seq);
      Cartesian3 tc(tg.t());
      real sum = 0;
      vector<real> s12v(odetest ? numl : 0);
      vector<Ellipsoid3::vec3> R2v, V2v;
      if (!inverse) {
        real s0, bet1, omg1, alp1, bet2, omg2, alp2;
        tg.Inverse(-90, 0, 90, 180, s0, alp1, alp2);
        std::uniform_real_distribution<long double> uni(0, (long double)(s0));
        for (int i = 0; i < num; ++i) {
          Ellipsoid3::vec3 R1, V1;
          tc.cart2rand(g, R1, V1);
          if (!odetest) {
            tc.cart2toany(R1, V1, Cartesian3::ELLIPSOIDAL, bet1, omg1, alp1);
            GeodesicLine3 lx(tg, bet1, omg1, alp1);
            for (int j = 0; j < numl; ++j) {
              real s12 = real(uni(g));
              if (numl == 1) s12 = fmax(s12, real(uni(g)));
              lx.Position(s12, bet2, omg2, alp2);
              sum += fabs(alp2);
            }
          } else {
#if HAVE_BOOST
            l.Reset(R1, V1);
            for (int j = 0; j < numl; ++j) {
              s12v[j] = real(uni(g));
              if (numl == 1) s12v[j] = fmax(s12v[j], real(uni(g)));
            }
            l.Position(s12v, R2v, V2v);
            if (numl > 0) sum += fabs(R2v[0][0]);
#endif
          }
        }
      } else {
        for (int i = 0; i < num; ++i) {
          real s12, bet1, omg1, alp1, bet2, omg2, alp2;
          Ellipsoid3::vec3 R1, R2;
          tc.cart2rand(g, R1);
          tc.cart2toany(R1, Cartesian3::ELLIPSOIDAL, bet1, omg1);
          tc.cart2rand(g, R2);
          tc.cart2toany(R2, Cartesian3::ELLIPSOIDAL, bet2, omg2);
          (void) tg.Inverse(bet1, omg1, bet2, omg2, s12, alp1, alp2);
          sum += s12;
        }
      }
      if (sum < 0) cout << sum << "\n";
    }
    (void) odep;
    (void) reportp;
    (void) odetest;
    (void) extended;
    (void) dense;
    (void) normp;
    (void) eps;
    return 0;
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

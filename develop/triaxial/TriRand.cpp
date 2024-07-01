#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/Utility.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"

using namespace GeographicLib;
using namespace std;

class RandLoc {
private:
  typedef Math::real real;
  std::seed_seq _seq;
  std::mt19937 _gen;
  int _div;
  std::uniform_int_distribution<int> _rndlat;
  std::uniform_int_distribution<int> _rndlon;
  std::uniform_int_distribution<int> _rndskew;
  std::uniform_int_distribution<int> _rndlat0;
  std::uniform_int_distribution<int> _rndlon0;
public:
  RandLoc(int div, int skew, unsigned seed)
    : _seq({seed})
    , _gen(_seq)
    , _div(div)
    , _rndlat(-90/div, 90/div)
    , _rndlon(-180/div+1, 180/div)
    , _rndskew(0, skew - 1)
    , _rndlat0(-1, 1)
    , _rndlon0(-1, 2)
  {}
  // Return random position
  void RandPos(int& bet, int& omg) {
    bet = _rndskew(_gen) == 0 ? 90 * _rndlat0(_gen) :
      _div * _rndlat(_gen);
    omg = _rndskew(_gen) == 0 ? 90 * _rndlon0(_gen) :
      _div * _rndlon(_gen);
  }
};

int usage(int retval, bool /*brief*/) {
  ( retval ? std::cerr : std::cout ) << "Bad input\n";
  return retval;
}

void report(const Triaxial& t, int bet1, int omg1, int bet2, int omg2) {
  typedef Math::real real;
  typedef Angle ang;
  cout << bet1 << " " << omg1 << " " << bet2 << " " << omg2 << " " << flush;
  TriaxialLine l =
    t.Inverse(ang(bet1), ang(omg1),
              ang(bet2), ang(omg2));
  real s12 = l.Distance();
  ang bet1a, omg1a, alp1a, bet2a, omg2a, alp2a;
  l.pos1(bet1a, omg1a, alp1a);
  l.Position(s12, bet2a, omg2a, alp2a);
  Triaxial::vec3 r2, v2;
  t.elliptocart2(bet2a, omg2a, alp2a, r2, v2);
  t.cart2toellip(ang(bet2), ang(omg2), v2,
                 alp2a);
#if GEOGRAPHICLIB_PRECISION <= 2
  int prec = 12;
#elif GEOGRAPHICLIB_PRECISION == 3
  int prec = 14;
#else
  int prec = 18;
#endif
  cout << fixed << setprecision(prec) << alp1a.degrees0() << " "
       << alp2a.degrees0() << " "
       << setprecision(prec+2) <<l.Distance() << endl;
}
int main(int argc, const char* const argv[]) {
  try {
    typedef Math::real real;
    Utility::set_digits();
    unsigned seed = 0;
    int num = 1000;
    int skew = 10;
    int div = 1;
    {
      Triaxial t(sqrt(real(2)), 1, sqrt(1/real(2)));
      int bet1, omg1, bet2, omg2;
      real alp1, alp2, s12;
      while (cin >> bet1 >> omg1 >> bet2 >> omg2 >> alp1 >> alp2 >> s12) {
        report(t, bet1, omg1, bet2, omg2);
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
      } else if (arg == "-s") {
        if (++m == argc) return usage(1, true);
        try {
          seed = Utility::val<unsigned>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "ssed " << argv[m] << " is not a number\n";
          return 1;
        }
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
    }
    if (seed == 0)
      seed = std::random_device()();
    std::cerr << "-s " << seed << "\n";
    RandLoc r(div, skew, seed);
    Triaxial t(sqrt(real(2)), 1, sqrt(1/real(2)));
    cout << fixed;
    if (1) {
      std::seed_seq seq({seed});
      std::mt19937 gen(seq);
      std::uniform_int_distribution<int> rnd(0, skew - 1);
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
    } else {
      for (int i = 0; i < num; ++i) {
      int bet1, omg1, bet2, omg2;
      r.RandPos(bet1, omg1);
      r.RandPos(bet2, omg2);
      report(t, bet1, omg1, bet2, omg2);
    }
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

/**
 * \file Cart3Convert.cpp
 * \brief Command line utility for computing geodesics on a triaxial ellipsoid
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="Cart3Convert.1.html">man page</a> for usage information.
 **********************************************************************/

// Counterpart of CartConvert

// Usual flags
//   -r (x, y, z to ellipsoidal to geographic)
//   -e (supplement with -t)
//   -w longfirst
//   -p prec
//   -h help
//   -d dms
//   -: colon
//   --help
//   --version
//   --comment-delimiter
//   + other I/O related flags

// New flags
//   -E ellipsoidal coords (default)
//   -G geodetic
//   -P parametric
//   -C geocentric
//   -3 include height (only for -E and -G)
//   -D include direction (only for -E without -3)

// Can't specify -3 and -D together
// Height for -G has its usual definition
// Height for -E is the increase in the minor semiaxis for the confocal
// ellipsoid (H = u - c).

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Angle.hpp>
#include <GeographicLib/Triaxial/Cartesian3.hpp>

#include "Cart3Convert.usage"

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    using namespace Triaxial;
    Utility::set_digits();
    using real = Math::real;
    using ang = Angle;
    using vec3 = Ellipsoid3::vec3;
    using coord = Cartesian3::coord;
    enum { GEODETIC, PARAMETRIC, GEOCENTRIC, ELLIPSOIDAL };
    coord mode = coord::ELLIPSOIDAL;
    bool threed = false, direction = false, reverse = false, dms = false,
      longfirst = false, randompts = false;
    real
      a = Constants::Triaxial_Earth_a(),
      b = Constants::Triaxial_Earth_b(),
      c = Constants::Triaxial_Earth_c(),
      e2 = -1, k2 = 0, kp2 = 0;
    int prec = 3, nrand = 0;
    unsigned long long seed = 0;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-E")
        mode = coord::ELLIPSOIDAL;
      else if (arg == "-G")
        mode = coord::GEODETIC;
      else if (arg == "-P")
        mode = coord::PARAMETRIC;
      else if (arg == "-C")
        mode = coord::GEOCENTRIC;
      else if (arg == "-GX")
        mode = coord::GEODETIC_X;
      else if (arg == "-PX")
        mode = coord::PARAMETRIC_X;
      else if (arg == "-CX")
        mode = coord::GEOCENTRIC_X;
      else if (arg == "-R") {
        if (++m == argc) return usage(1, true);
        try {
          nrand = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Number of randoms " << argv[m] << " is not a number\n";
          return 1;
        }
        randompts = true;
      } else if (arg == "-3")
        threed = true;
      else if (arg == "-D")
        direction = true;
      else if (arg == "-r")
        reverse = true;
      else if (arg == "--seed") {
        if (++m == argc) return usage(1, true);
        try {
          seed = Utility::val<unsigned long long>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-t") {
        if (m + 3 >= argc) return usage(1, true);
        try {
          a = Utility::val<real>(std::string(argv[m + 1]));
          b = Utility::val<real>(std::string(argv[m + 2]));
          c = Utility::val<real>(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -t: " << e.what() << "\n";
          return 1;
        }
        e2 = -1;
        m += 3;
      } else if (arg == "-e") {
        // Cayley ellipsoid sqrt([2,1,1/2]) is
        // -e 1 3/2 1/3 2/3
        if (m + 4 >= argc) return usage(1, true);
        try {
          b = Utility::val<real>(std::string(argv[m + 1]));
          e2 = Utility::fract<real>(std::string(argv[m + 2]));
          k2 = Utility::fract<real>(std::string(argv[m + 3]));
          kp2 = Utility::fract<real>(std::string(argv[m + 4]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        a = -1;
        m += 4;
      } else if (arg == "-d") {
        dms = true;
        dmssep = '\0';
      } else if (arg == "-:") {
        dms = true;
        dmssep = ':';
      } else if (arg == "-w")
        longfirst = !longfirst;
      else if (arg == "-p") {
        if (++m == argc) return usage(1, true);
        try {
          prec = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "--input-string") {
        if (++m == argc) return usage(1, true);
        istring = argv[m];
      } else if (arg == "--input-file") {
        if (++m == argc) return usage(1, true);
        ifile = argv[m];
      } else if (arg == "--output-file") {
        if (++m == argc) return usage(1, true);
        ofile = argv[m];
      } else if (arg == "--line-separator") {
        if (++m == argc) return usage(1, true);
        if (std::string(argv[m]).size() != 1) {
          std::cerr << "Line separator must be a single character\n";
          return 1;
        }
        lsep = argv[m][0];
      } else if (arg == "--comment-delimiter") {
        if (++m == argc) return usage(1, true);
        cdelim = argv[m];
      } else if (arg == "--version") {
        std::cout << argv[0] << ": GeographicLib version "
                  << GEOGRAPHICLIB_VERSION_STRING << "\n";
        return 0;
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
    }

    Cartesian3 tc(e2 >= 0 ? Ellipsoid3(b, e2, k2, kp2) : Ellipsoid3(a, b, c));

    if (randompts) {
      if (threed) {
        std::cerr << "Cannot specify -3 for random points\n";
        return 1;
      }
      if (!ifile.empty() || !istring.empty()) {
        std::cerr
          << "Cannot specify --input-{string,file} with random points\n";
      }
    } else {
      if (direction && threed) {
        std::cerr << "Cannot specify both -3 and -D\n";
        return 1;
      }
    }

    if (!ifile.empty() && !istring.empty()) {
      std::cerr << "Cannot specify --input-string and --input-file together\n";
      return 1;
    }
    if (ifile == "-") ifile.clear();
    std::ifstream infile;
    std::istringstream instring;
    if (!ifile.empty()) {
      infile.open(ifile.c_str());
      if (!infile.is_open()) {
        std::cerr << "Cannot open " << ifile << " for reading\n";
        return 1;
      }
    } else if (!istring.empty()) {
      std::string::size_type m = 0;
      while (true) {
        m = istring.find(lsep, m);
        if (m == std::string::npos)
          break;
        istring[m] = '\n';
      }
      instring.str(istring);
    }
    std::istream* input = !ifile.empty() ? &infile :
      (!istring.empty() ? &instring : &std::cin);

    std::ofstream outfile;
    if (ofile == "-") ofile.clear();
    if (!ofile.empty()) {
      outfile.open(ofile.c_str());
      if (!outfile.is_open()) {
        std::cerr << "Cannot open " << ofile << " for writing\n";
        return 1;
      }
    }
    std::ostream* output = !ofile.empty() ? &outfile : &std::cout;

    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    using std::round, std::log10, std::ceil, std::signbit;
    int disprec = std::max(0, prec + int(round(log10(6400000/b)))),
      angprec = prec + 5, vecprec = prec + 7;
    if (randompts) {
      //                         C a r t 3 C o n
      unsigned long long s1 = 0x4361727433436f6eULL, s2 = seed;
      if (seed == 0) {
        s1 = std::random_device()();
        s2 = std::random_device()();
      }
      std::seed_seq seq{s1, s2};
      std::mt19937 g(seq);
      for (int i = 0; i < nrand; ++i) {
        vec3 r, v = {0,0,0};
        if (direction)
          tc.cart2rand(g, r, v);
        else
          tc.cart2rand(g, r);
        if (!reverse) {
          *output << Utility::str(r[0], disprec) << " "
                  << Utility::str(r[1], disprec) << " "
                  << Utility::str(r[2], disprec);
          if (direction)
            *output << " "
                    << Utility::str(v[0], vecprec) << " "
                    << Utility::str(v[1], vecprec) << " "
                    << Utility::str(v[2], vecprec);
        } else {
          ang bet, omg, alp;
          tc.cart2toany(r, v, mode, bet, omg, alp);
          *output << ang::LatLonString(bet, omg,
                                       angprec, dms, dmssep, longfirst);
          if (direction)
            *output << " " << ang::AzimuthString(alp, angprec, dms, dmssep);
        }
        *output << "\n";
      }
      return 0;
    }
    std::string s, eol, sbet, somg, salp, sh, sx, sy, sz, svx, svy, svz, strc;
    std::istringstream str;
    int retval = 0;
    vec3 r = {0,0,0}, v = {0,0,0};
    ang bet, omg, alp;
    real h = 0;
    while (std::getline(*input, s)) {
      try {
        eol = "\n";
        if (!cdelim.empty()) {
          std::string::size_type m = s.find(cdelim);
          if (m != std::string::npos) {
            eol = " " + s.substr(m) + "\n";
            s = s.substr(0, m);
          }
        }
        // READ
        str.clear(); str.str(s);
        if (reverse) {
          if (!(str >> sx >> sy >> sz))
            throw GeographicErr("Incomplete input: " + s);
          r[0] = Utility::val<real>(sx);
          r[1] = Utility::val<real>(sy);
          r[2] = Utility::val<real>(sz);
          if (direction) {
            if (!(str >> svx >> svy >> svz))
              throw GeographicErr("Incomplete input: " + s);
            v[0] = Utility::val<real>(svx);
            v[1] = Utility::val<real>(svy);
            v[2] = Utility::val<real>(svz);
          }
        } else {
          if (!(str >> sbet >> somg))
            throw GeographicErr("Incomplete input: " + s);
          ang::DecodeLatLon(sbet, somg, bet, omg, longfirst);
          if (mode != coord::ELLIPSOIDAL && (bet.n() != 0 || signbit(bet.c())))
            throw GeographicErr("Latitude outside range [-90,90]: " + s);
          if (threed) {
            if (!(str >> sh))
              throw GeographicErr("Incomplete input: " + s);
            h = Utility::val<real>(sh);
          } else if (direction) {
            if (!(str >> salp))
              throw GeographicErr("Incomplete input: " + s);
            alp = ang::DecodeAzimuth(salp);
          }
        }
        if (str >> strc)
          throw GeographicErr("Extraneous input: " + strc);
        // PROCESS
        if (reverse) {
          if (threed)
            tc.carttoany(r, mode, bet, omg, h);
          else if (direction)
            tc.cart2toany(r, v, mode, bet, omg, alp);
          else
            tc.cart2toany(r, mode, bet, omg);
        } else {
          if (threed)
            tc.anytocart(mode, bet, omg, h, r);
          else if (direction)
            tc.anytocart2(mode, bet, omg, alp, r, v);
          else
            tc.anytocart2(mode, bet, omg, r);
        }
        // WRITE
        if (reverse) {
          *output << ang::LatLonString(bet, omg,
                                       angprec, dms, dmssep, longfirst);
          if (threed)
            *output << " " << Utility::str(h, disprec);
          else if (direction)
            *output << " " << ang::AzimuthString(alp, angprec, dms, dmssep);
        } else {
          *output << Utility::str(r[0], disprec) << " "
                  << Utility::str(r[1], disprec) << " "
                  << Utility::str(r[2], disprec);
          if (direction)
            *output << " "
                    << Utility::str(v[0], vecprec) << " "
                    << Utility::str(v[1], vecprec) << " "
                    << Utility::str(v[2], vecprec);
        }
        *output << eol;
      }
      catch (const std::exception& e) {
        // Write error message cout so output lines match input lines
        *output << "ERROR: " << e.what() << "\n";
        retval = 1;
      }
    }
    return retval;
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

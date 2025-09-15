/**
 * \file Conformal3Proj.cpp
 * \brief Command line utility for computing geodesics on a triaxial ellipsoid
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="Conformal3Proj.1.html">man page</a> for usage information.
 **********************************************************************/

// Usual flags
//   -r reverse
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
//   -S project to sphere
//   -Q report sphere radius on stderr

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Angle.hpp>
#include <GeographicLib/Triaxial/Conformal3.hpp>

#include "Conformal3Proj.usage"

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    using namespace Triaxial;
    typedef Math::real real;
    typedef Angle ang;
    Utility::set_digits();
    bool sphere = false, reverse = false, dms = false,
      longfirst = false, report = false;
    real
      a = Constants::Triaxial_Earth_a(),
      b = Constants::Triaxial_Earth_b(),
      c = Constants::Triaxial_Earth_c(),
      e2 = -1, k2 = 0, kp2 = 0;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-r")
        reverse = true;
      else if (arg == "-S")
        sphere = true;
      else if (arg == "-Q")
        report = true;
      else if (arg == "-t") {
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

    Conformal3 tc(e2 >= 0 ? Ellipsoid3(b, e2, k2, kp2) : Ellipsoid3(a, b, c));

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
      angprec = prec + 5, scalprec = prec + 7;
    if (report)
      std::cerr << "sphere-radius "
                << Utility::str(tc.SphereRadius(), disprec) << "\n";
    std::string s, eol, stra, strb, strc;
    std::istringstream str;
    int retval = 0;
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
        if (!(str >> stra >> strb))
          throw GeographicErr("Incomplete input: " + s);
        if (str >> strc)
          throw GeographicErr("Extraneous input: " + strc);
        ang gam{}; real k = 1;
        if (reverse) {
          ang bet, omg;
          if (sphere) {
            real lat, lon;
            DMS::DecodeLatLon(stra, strb, lat, lon, longfirst);
            ang phi{lat}, lam{lon};
            tc.ReverseSphere(phi, lam, bet, omg, gam, k);
          } else {
            real x = Utility::val<real>(stra), y = Utility::val<real>(strb);
            tc.Reverse(x, y, bet, omg, k);
          }
          *output << ang::LatLonString(bet, omg,
                                       angprec, dms, dmssep, longfirst);
        } else {
          ang bet, omg;
          ang::DecodeLatLon(stra, strb, bet, omg, longfirst);
          if (sphere) {
            ang phi, lam;
            tc.ForwardSphere(bet, omg, phi, lam, gam, k);
            *output << ang::LatLonString(phi, lam,
                                         angprec, dms, dmssep, longfirst);
          } else {
            real x, y;
            tc.Forward(bet, omg, x, y, k);
            *output << Utility::str(x, disprec) << " "
                    << Utility::str(y, disprec);
          }
        }
        if (sphere)
          *output << " " << ang::AzimuthString(gam, angprec, dms, dmssep);
        *output << " " << Utility::str(k, scalprec) << eol;
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

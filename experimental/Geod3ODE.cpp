/**
 * \file Geod3ODE.cpp
 * \brief Command line utility for using an ODE solver for geodesics on a
 * triaxial ellipsoid
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * Use "Geod3ODE --help" for brief documentation.
 **********************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Angle.hpp>

#if defined(_MSC_VER)
// Squelch warning triggered by boost:
//   4127: conditional expression is constant
#  pragma warning (disable: 4127)
#endif
#include "TriaxialGeodesicODE.hpp"

// #include "GeodSolve.usage"

using real = GeographicLib::Math::real;
using ang = GeographicLib::Angle;

std::string ErrorString(real err, int prec) {
  std::ostringstream s;
  s << std::scientific << std::setprecision(prec) << err;
  return s.str();
}

int usage(int retval, bool /*brief*/) {
  (retval ? std::cerr : std:: cout ) << "Usage:\n"
"\n"
"  This implements some of the functionality of Geod3Solve(1) by integrating the\n"
"  ordinary differential equations for the geodesic.  Only direct geodesic\n"
"  calculations are supported.\n"
"\n"
"  The following options of Geod3Solve(1) are supported\n"
"    -t a b c | -e b e2 k2 kp2 \n"
"    -L bet1 omg1 alp1\n"
"    -u\n"
"    -d | -:\n"
"    -w \n"
"    -f\n"
"    -p prec\n"
"\n"
"  The following options of Geod3Solve(1) are not supported\n"
"    -i\n"
"    -e2\n"
"    -u\n"
"\n"
"  The following are new options\n"
"\n"
"    -b\n"
"      bufferd mode (only useful with the -L option).  Causes all the s12 values\n"
"      to be buffered and fed into the integrator at the end.  This sorts the\n"
"      entries so that the integrator doesn't have to the continually restarted.\n"
"\n"
"    -x\n"
"      extended mode.  Computes and prints the reduced length, m12, and the\n"
"      geodesic scales, M12, M21.\n"
"\n"
"    --eps eps\n"
"      sets the eps parameter in the constructor for TriaxialGeodesicODE.\n"
"\n"
"    --dense\n"
"      use the dense solver allowing interpolated way points to be computed\n"
"      inexpensively.\n"
"\n"
"    --normp\n"
"      force the solution vector onto the ellipsoid when computing the\n"
"      acceleration.\n"
"\n"
"    --errors\n"
"      print error estimates, the distance from the ellipsoid (in meters) and\n"
"      the deviation of the velocity from a unit tangential vector.\n"
"\n"
"    --steps\n"
"      print the number of integration steps and the number of times the\n"
"      acceleration was computed.\n";
  return retval;
}

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    using namespace Triaxial;
    using namespace experimental;
    Utility::set_digits();
    bool dms = false, longfirst = false,
      linecalc = false, extended = false, dense = false, normp = false,
      buffered = false, full = false, errors = false, steps = false;
    real
      a = Constants::Triaxial_Earth_a(),
      b = Constants::Triaxial_Earth_b(),
      c = Constants::Triaxial_Earth_c(),
      e2 = -1, k2 = 0, kp2 = 0, eps = 0;
    ang bet1 = ang(0), omg1 = ang(0), alp1 = ang(0), bet2, omg2, alp2;
    real s12, m12, M12, M21;
    std::vector<real> s12v;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-t") {
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
        // -e 1 3/2 1 2
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
      } else if (arg == "-x")
        extended = true;
      else if (arg == "-L") {
        linecalc = true;
        if (m + 3 >= argc) return usage(1, true);
        try {
          ang::DecodeLatLon(std::string(argv[m + 1]),
                              std::string(argv[m + 2]),
                              bet1, omg1, longfirst);
          alp1 = ang::DecodeAzimuth(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -L: " << e.what() << "\n";
          return 1;
        }
        m += 3;
      } else if (arg == "--eps") {
        if (m + 1 >= argc) return usage(1, true);
        try {
          using std::pow;
          eps = pow(std::numeric_limits<real>::epsilon(),
                    Utility::fract<real>(std::string(argv[m + 1])));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of --eps: " << e.what() << "\n";
          return 1;
        }
        m += 1;
      } else if (arg == "--dense")
        dense = true;
      else if (arg == "--normp")
        normp = true;
      else if (arg == "--errors")
        errors = true;
      else if (arg == "--steps")
        steps = true;
      else if (arg == "-b")
        buffered = true;
      else if (arg == "-f")
        full = true;
      else if (arg == "-d") {
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

    Ellipsoid3 t(e2 >= 0 ? Ellipsoid3(b, e2, k2, kp2) : Ellipsoid3(a, b, c));
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

    using std::round, std::log10;
    int disprec = int(round(log10(6400000/b)));
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol, sbet1, somg1, salp1, sbet2, somg2, salp2, ss12, strc;
    std::istringstream str;
    int retval = 0;
    buffered = buffered && linecalc;
    errors = errors && !buffered;
    TriaxialGeodesicODE l = linecalc ?
      TriaxialGeodesicODE(t, bet1, omg1, alp1, extended, dense, normp, eps) :
      TriaxialGeodesicODE(t, extended, dense, normp, eps);

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
        str.clear(); str.str(s);
        if (!(linecalc ? (str >> ss12) :
              (str >> sbet1 >> somg1 >> salp1 >> ss12)))
          throw GeographicErr("Incomplete input: " + s);
        if (str >> strc)
          throw GeographicErr("Extraneious input: " + s);
        s12 = Utility::val<real>(ss12);
        if (linecalc) {
          if (buffered) s12v.push_back(s12);
        } else {
          ang::DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
          alp1 = ang::DecodeAzimuth(salp1);
          l.Reset(bet1, omg1, alp1);
        }
        if (!buffered) {
          auto errs = l.Position(s12, bet2, omg2, alp2, m12, M12, M21);
          if (full)
            *output << ang::LatLonString(bet1, omg1, prec, dms, dmssep,
                                         longfirst) << " "
                    << ang::AzimuthString(alp1, prec, dms, dmssep) << " ";
          *output << ang::LatLonString(bet2, omg2, prec, dms, dmssep,
                                       longfirst) << " "
                  << ang::AzimuthString(alp2, prec, dms, dmssep);
          if (full)
            *output << " " << Utility::str(s12, prec + disprec);
          if (extended)
            *output << " " << Utility::str(m12, prec + disprec)
                    << " " << Utility::str(M12, prec+7)
                    << " " << Utility::str(M21, prec+7);
          if (steps)
            *output << " " << l.NSteps() << " " << l.IntSteps();
          if (errors)
            *output << " " << ErrorString(errs.first, 2)
                    << " " << ErrorString(errs.second, 2);
          *output << eol;
        }
      }
      catch (const std::exception& e) {
        if (buffered)
          s12v.push_back(Math::NaN());
        else
          // Write error message cout so output lines match input lines
          *output << "ERROR: " << e.what() << " " << s << "\n";
        retval = 1;
      }
    }

    if (buffered) {
      std::vector<ang> bet2v, omg2v, alp2v;
      std::vector<real> m12v, M12v, M21v;
      l.Position(s12v, bet2v, omg2v, alp2v, m12v, M12v, M21v);
      for (size_t i = 0; i < s12v.size(); ++i) {
          if (full)
            *output << ang::LatLonString(bet1, omg1, prec, dms, dmssep,
                                         longfirst) << " "
                    << ang::AzimuthString(alp1, prec, dms, dmssep) << " ";
          *output << ang::LatLonString(bet2v[i], omg2v[i], prec, dms, dmssep,
                                  longfirst) << " "
                  << ang::AzimuthString(alp2v[i], prec, dms, dmssep);
          if (full)
            *output << " " << Utility::str(s12v[i], prec + disprec);
          if (extended)
            *output << " " << Utility::str(m12v[i], prec + disprec)
                    << " " << Utility::str(M12v[i], prec+7)
                    << " " << Utility::str(M21v[i], prec+7);
          *output << eol;
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

/**
 * \file Geod3Solve.cpp
 * \brief Command line utility for computing geodesics on a triaxial ellipsoid
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="Geod3Solve.1.html">man page</a> for usage information.
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
#include <GeographicLib/Triaxial/Geodesic3.hpp>

#include "Geod3Solve.usage"

using real = GeographicLib::Math::real;
using ang = GeographicLib::Angle;

void BiaxialCoords(bool fwd, real f, ang& bet, ang& omg) {
  using std::isnan, std::signbit;
  if (isnan(f)) return;
  // fwd: biaxial -> triaxial
  // !fwd: triaxial -> biaxial
  if (fwd)
    bet = bet.modang(1 - f);
  if (signbit(f)) {
    omg -= ang::cardinal(1);
    std::swap(bet, omg);
    omg += ang::cardinal(1);
  }
  if (!fwd)
    bet = bet.modang(1 / (1 - f));
}

void BiaxialCoords(bool fwd, real f, ang& bet, ang& omg, ang& alp) {
  using std::isnan, std::signbit;
  if (isnan(f)) return;
  // fwd: biaxial -> triaxial
  // !fwd: triaxial -> biaxial
  BiaxialCoords(fwd, f, bet, omg);
  if (signbit(f)) alp.reflect(false, false, true);
}

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    using namespace Triaxial;
    using std::signbit, std::isnan, std::fabs;
    Utility::set_digits();
    bool inverse = false,
      dms = false, full = false, unroll = false,
      longfirst = false,
      linecalc = false;
    real
      a = Constants::Triaxial_Earth_a(),
      b = Constants::Triaxial_Earth_b(),
      c = Constants::Triaxial_Earth_c(),
      // Markers to determine how the ellipsoid is specified
      e2 = -1, f = Math::NaN(), k2 = 0, kp2 = 0;
    ang bet1, omg1, alp1, bet2, omg2, alp2;
    real s12;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-i") {
        inverse = true;
        linecalc = false;
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
        f = Math::NaN(); e2 = -1;
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
        a = -1; f = Math::NaN();
        m += 4;
      } else if (arg == "-e2") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          b = Utility::val<real>(std::string(argv[m + 1]));
          f = Utility::fract<real>(std::string(argv[m + 2]));
        }
        // f > 0, k2 = 1, kp2 = 0
        // a = b; c = b * (1 - f)
        // e2 = factor(subst([a = A, b = A, c = A*(1-f)],(a^2 - c^2)/b^2))
        //    = f * (2 - f)
        // f < 0, k2 = 0, kp2 = 1
        // a = A * (1 - f); c = b
        // e2 = factor(subst([a = A*(1-f), b = A, c = A],(a^2 - c^2)/b^2))
        //    = -f * (2 - f)
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -e2: " << e.what() << "\n";
          return 1;
        }
        a = e2 = -1;
        m += 2;
      } else if (arg == "-L") {
        linecalc = true;
        inverse = false;
        if (m + 3 >= argc) return usage(1, true);
        try {
          ang::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            bet1, omg1, longfirst);
          alp1 = ang::DecodeAzimuth(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -L: " << e.what() << "\n";
          return 1;
        }
        m += 3;
      } else if (arg == "-u")
        unroll = true;
      else if (arg == "-d") {
        dms = true;
        dmssep = '\0';
      } else if (arg == "-:") {
        dms = true;
        dmssep = ':';
      } else if (arg == "-w")
        longfirst = !longfirst;
      else if (arg == "-f")
        full = true;
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

    Geodesic3 t = e2 >= 0 ? Geodesic3(b, e2, k2, kp2) :
      !isnan(f) ? Geodesic3(b, fabs(f) * (2 - f),
                                   signbit(f) ? 0 : 1, signbit(f) ? 1 : 0) :
      Geodesic3(a, b, c);

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

    if (linecalc) {
      BiaxialCoords(true, f, bet1, omg1, alp1);
    }
    std::unique_ptr<GeodesicLine3> lp = linecalc ?
      std::make_unique<GeodesicLine3>(t, bet1, omg1, alp1) : nullptr;
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    using std::round, std::log10, std::ceil;
    int disprec = std::max(0, prec + int(round(log10(6400000/b)))),
      angprec = prec + 5;
    std::string s, eol, sbet1, somg1, salp1, sbet2, somg2, salp2, ss12, strc;
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
        str.clear(); str.str(s);

        if (inverse) {
          if (!(str >> sbet1 >> somg1 >> sbet2 >> somg2))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          ang::DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
          ang::DecodeLatLon(sbet2, somg2, bet2, omg2, longfirst);
          BiaxialCoords(true, f, bet1, omg1);
          BiaxialCoords(true, f, bet2, omg2);
          GeodesicLine3 l = t.Inverse(bet1, omg1, bet2, omg2,
                                      s12, alp1, alp2);
          if (unroll && full) {
            l.Position(s12, bet2, omg2, alp2);
          }
          BiaxialCoords(false, f, bet1, omg1, alp1);
          BiaxialCoords(false, f, bet2, omg2, alp2);
          if (full)
            *output << ang::LatLonString(bet1, omg1,
                                         angprec, dms, dmssep, longfirst)
                    << " " << ang::AzimuthString(alp1, angprec, dms, dmssep)
                    << " "
                    << ang::LatLonString(bet2, omg2,
                                         angprec, dms, dmssep, longfirst)
                    << " " << ang::AzimuthString(alp2, angprec, dms, dmssep)
                    << " "
                    << Utility::str(s12, disprec) << eol;
          else
            *output << ang::AzimuthString(alp1, angprec, dms, dmssep) << " "
                    << ang::AzimuthString(alp2, angprec, dms, dmssep) << " "
                    << Utility::str(s12, disprec) << eol;
        } else {
          if (linecalc) {
            if (!(str >> ss12))
              throw GeographicErr("Incomplete input: " + s);
          } else {
            if (!(str >> sbet1 >> somg1 >> salp1 >> ss12))
              throw GeographicErr("Incomplete input: " + s);
          }
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          s12 = Utility::val<real>(ss12);
          if (linecalc)
            lp->Position(s12, bet2, omg2, alp2);
          else {
            ang::DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
            alp1 = ang::DecodeAzimuth(salp1);
            BiaxialCoords(true, f, bet1, omg1, alp1);
            t.Direct(bet1, omg1, alp1, s12, bet2, omg2, alp2);
          }
          if (!unroll) {
            Ellipsoid3::AngNorm(bet2, omg2, alp2, !isnan(f) && signbit(f));
            bet2 = bet2.base();
            omg2 = omg2.base();
            alp2 = alp2.base();
          }
          BiaxialCoords(false, f, bet1, omg1, alp1);
          BiaxialCoords(false, f, bet2, omg2, alp2);
          if (full)
            *output << ang::LatLonString(bet1, omg1,
                                         angprec, dms, dmssep, longfirst)
                    << " " << ang::AzimuthString(alp1, angprec, dms, dmssep)
                    << " "
                    << ang::LatLonString(bet2, omg2,
                                         angprec, dms, dmssep, longfirst)
                    << " " << ang::AzimuthString(alp2, angprec, dms, dmssep)
                    << " "
                    << Utility::str(s12, disprec) << eol;
          else
            *output << ang::LatLonString(bet2, omg2,
                                         angprec, dms, dmssep, longfirst)
                    << " "
                    << ang::AzimuthString(alp2, angprec, dms, dmssep) << eol;
        }
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

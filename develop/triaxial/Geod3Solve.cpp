/**
 * \file GeodSolve3.cpp
 * \brief Command line utility for computing geodesics on a triaxial ellipsoid
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="GeodSolve.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include "TriaxialLine.hpp"
#include "TriaxialODE.hpp"

// #include "GeodSolve.usage"

typedef GeographicLib::Math::real real;

void DecodeLatLon(const std::string& stra, const std::string& strb,
                  real& lat, real& lon,
                  bool longfirst) {
  using namespace GeographicLib;
  real a, b;
  DMS::flag ia, ib;
  a = DMS::Decode(stra, ia);
  b = DMS::Decode(strb, ib);
  if (ia == DMS::NONE && ib == DMS::NONE) {
    // Default to lat, long unless longfirst
    ia = longfirst ? DMS::LONGITUDE : DMS::LATITUDE;
    ib = longfirst ? DMS::LATITUDE : DMS::LONGITUDE;
  } else if (ia == DMS::NONE)
    ia = DMS::flag(DMS::LATITUDE + DMS::LONGITUDE - ib);
  else if (ib == DMS::NONE)
    ib = DMS::flag(DMS::LATITUDE + DMS::LONGITUDE - ia);
  if (ia == ib)
    throw GeographicErr("Both " + stra + " and "
                        + strb + " interpreted as "
                        + (ia == DMS::LATITUDE ? "latitudes" : "longitudes"));
  lat = ia == DMS::LATITUDE ? a : b;
  lon = ia == DMS::LATITUDE ? b : a;
}

real DecodeAzimuth(const std::string& azistr) {
  using namespace GeographicLib;
  DMS::flag ind;
  real azi = DMS::Decode(azistr, ind);
  if (ind == DMS::LATITUDE)
    throw GeographicErr("Azimuth " + azistr
                        + " has a latitude hemisphere, N/S");
  return azi;
}

std::string BetOmgString(real bet, real omg, int prec, bool dms, char dmssep,
                         bool longfirst) {
  using namespace GeographicLib;
  std::string
    betstr = dms ? DMS::Encode(bet, prec + 5, DMS::LATITUDE, dmssep) :
    DMS::Encode(bet, prec + 5, DMS::NUMBER),
    omgstr = dms ? DMS::Encode(omg, prec + 5, DMS::LONGITUDE, dmssep) :
    DMS::Encode(omg, prec + 5, DMS::NUMBER);
  return
    (longfirst ? omgstr : betstr) + " " + (longfirst ? betstr : omgstr);
}

std::string AzimuthString(real alp, int prec, bool dms, char dmssep) {
  using namespace GeographicLib;
  return dms ? DMS::Encode(alp, prec + 5, DMS::AZIMUTH, dmssep) :
    DMS::Encode(alp, prec + 5, DMS::NUMBER);
}

std::string DistanceString(real s12, int prec) {
  using namespace GeographicLib;
  return  Utility::str(s12, prec);
}

real ReadDistance(const std::string& s, bool arcmode, bool fraction = false) {
  using namespace GeographicLib;
  return fraction ? Utility::fract<real>(s) :
    (arcmode ? DMS::DecodeAngle(s) : Utility::val<real>(s));
}

int usage(int retval, bool /*brief*/) { return retval; }

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    enum { NONE = 0, LINE, DIRECT, INVERSE };
    Utility::set_digits();
    bool inverse = false, arcmode = false,
      dms = false, full = false, unroll = false,
      longfirst = false, azi2back = false, fraction = false,
      lineseq = false, cart = false;
    real
      a = 6378172, b = 6378102, c = 6356752;
    real bet1, omg1, alp1, bet2, omg2, alp2, s12,
      mult = 1;
    int linecalc = NONE, prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);
    real ds = 0;
    long nmin = 0, nmax = 0;

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-i") {
        inverse = true;
        linecalc = NONE;
      } else if (arg == "-a")
        arcmode = !arcmode;
      else if (arg == "-F")
        fraction = true;
      else if (arg == "-L") {
        inverse = false;
        linecalc = LINE;
        if (m + 3 >= argc) return usage(1, true);
        try {
          DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                       bet1, omg1, longfirst);
          alp1 = DecodeAzimuth(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -L: " << e.what() << "\n";
          return 1;
        }
        m += 3;
      } else if (arg == "-D") {
        inverse = false;
        linecalc = DIRECT;
        if (m + 4 >= argc) return usage(1, true);
        try {
          DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                       bet1, omg1, longfirst);
          alp1 = DecodeAzimuth(std::string(argv[m + 3]));
          s12 = ReadDistance(std::string(argv[m + 4]), arcmode);
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -D: " << e.what() << "\n";
          return 1;
        }
        m += 4;
      } else if (arg == "-I") {
        inverse = false;
        linecalc = INVERSE;
        if (m + 4 >= argc) return usage(1, true);
        try {
          DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                       bet1, omg1, longfirst);
          DecodeLatLon(std::string(argv[m + 3]), std::string(argv[m + 4]),
                       bet2, omg2, longfirst);
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -I: " << e.what() << "\n";
          return 1;
        }
        m += 4;
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
        m += 3;
      } else if (arg == "-s") {
        if (m + 3 >= argc) return usage(1, true);
        try {
          ds = Utility::val<real>(std::string(argv[m + 1]));
          nmin = Utility::val<long>(std::string(argv[m + 2]));
          nmax = Utility::val<long>(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -t: " << e.what() << "\n";
          return 1;
        }
        m += 3;
        lineseq = true;
      } else if (arg == "--cart")
        cart = true;
      else if (arg == "-u")
        unroll = true;
      else if (arg == "-d") {
        dms = true;
        dmssep = '\0';
      } else if (arg == "-:") {
        dms = true;
        dmssep = ':';
      } else if (arg == "-w")
        longfirst = !longfirst;
      else if (arg == "-b")
        azi2back = true;
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

    const Triaxial t(a, b, c);
    TriaxialLine ls(t);
    if (linecalc)
      ls = TriaxialLine(t, bet1, omg1, alp1);
    if (lineseq) {
      if (cart) {
        Triaxial::vec3 r1, v1;
        t.elliptocart2(AuxAngle::degrees(bet1), AuxAngle::degrees(omg1),
                       AuxAngle::degrees(alp1), r1, v1);
        std::vector<Triaxial::vec3> r2, v2;
        TriaxialODE direct(t, r1, v1);
        direct.Position(ds, nmin, nmax, r2, v2);
        for (size_t i = 0; i < r2.size(); ++i) {
          AuxAngle bet2, omg2, alp2;
          t.cart2toellip(r2[i], v2[i], bet2, omg2, alp2);
          *output << BetOmgString(bet2.degrees(), omg2.degrees(),
                                  prec, dms, dmssep, longfirst)
                  << " " << AzimuthString(alp2.degrees(), prec, dms, dmssep)
                  << "\n";
        }
      } else {
        for (long n = nmin; n <= nmax; ++n) {
          real s12 = n * ds;
          ls.Position(s12, bet2, omg2, alp2, unroll);
          *output << BetOmgString(bet2, omg2, prec, dms, dmssep, longfirst)
                  << " " << AzimuthString(alp2, prec, dms, dmssep)
                  << "\n";
        }
      }
      return 0;
    }
    using std::round; using std::log10;
    int disprec = int(round(log10(b/6400000)));
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol, sbet1, somg1, sbet2, somg2, salp1, ss12, strc;
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
          DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
          DecodeLatLon(sbet2, somg2, bet2, omg2, longfirst);
          if (full) {
            if (unroll) {
              real e;
              omg2 = omg1 + Math::AngDiff(omg1, omg2, e);
              omg2 += e;
            } else {
              omg1 = Math::AngNormalize(omg1);
              omg2 = Math::AngNormalize(omg2);
            }
            *output << BetOmgString(bet1, omg1, prec, dms, dmssep, longfirst)
                    << " ";
          }
          *output << AzimuthString(alp1, prec, dms, dmssep) << " ";
          if (full)
            *output << BetOmgString(bet2, omg2, prec, dms, dmssep, longfirst)
                    << " ";
          if (azi2back) {
            using std::copysign;
            // map +/-0 -> -/+180; +/-180 -> -/+0
            // this depends on abs(alp2) <= 180
            alp2 = copysign(alp2 + copysign(real(Math::hd), -alp2), -alp2);
          }
          *output << AzimuthString(alp2, prec, dms, dmssep) << " "
                  << DistanceString(s12, prec + disprec);
          *output << eol;
        } else {
          if (linecalc) {
            if (!(str >> ss12))
              throw GeographicErr("Incomplete input: " + s);
            if (str >> strc)
              throw GeographicErr("Extraneous input: " + strc);
            // In fraction mode input is read as a distance
            s12 = ReadDistance(ss12, !fraction && arcmode, fraction) * mult;
            ls.Position(s12, bet2, omg2, alp2);
          } else {
            if (!(str >> sbet1 >> somg1 >> salp1 >> ss12))
              throw GeographicErr("Incomplete input: " + s);
            if (str >> strc)
              throw GeographicErr("Extraneous input: " + strc);
            DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
            alp1 = DecodeAzimuth(salp1);
            s12 = ReadDistance(ss12, arcmode);
            ls = TriaxialLine(t, bet1, omg1, alp1);
            ls.Position(s12, bet2, omg2, alp2);
          }
          if (azi2back) {
            using std::copysign;
            // map +/-0 -> -/+180; +/-180 -> -/+0
            // this depends on abs(alp2) <= 180
            alp2 = copysign(alp2 + copysign(real(Math::hd), -alp2), -alp2);
          }
          *output << BetOmgString(bet2, omg2, prec, dms, dmssep, longfirst)
                  << " " << AzimuthString(alp2, prec, dms, dmssep);
          if (full)
            *output << " "
                    << DistanceString(s12, prec + disprec);
          *output << eol;
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

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
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"
#include "TriaxialODE.hpp"

// #include "GeodSolve.usage"

typedef GeographicLib::Math::real real;
typedef GeographicLib::Angle ang;

void DecodeLatLon(const std::string& stra, const std::string& strb,
                  ang& lat, ang& lon,
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
  lat = ang(ia == DMS::LATITUDE ? a : b);
  lon = ang(ia == DMS::LATITUDE ? b : a);
}

ang DecodeAzimuth(const std::string& azistr) {
  using namespace GeographicLib;
  DMS::flag ind;
  real azi = DMS::Decode(azistr, ind);
  if (ind == DMS::LATITUDE)
    throw GeographicErr("Azimuth " + azistr
                        + " has a latitude hemisphere, N/S");
  return ang(azi);
}

std::string BetOmgString(ang bet, ang omg, int prec, bool dms, char dmssep,
                         bool longfirst) {
  using namespace GeographicLib;
  std::string
    betstr = dms ? DMS::Encode(real(bet), prec + 5, DMS::LATITUDE, dmssep) :
    DMS::Encode(real(bet), prec + 5, DMS::NUMBER),
    omgstr = dms ? DMS::Encode(real(omg), prec + 5, DMS::LONGITUDE, dmssep) :
    DMS::Encode(real(omg), prec + 5, DMS::NUMBER);
  return
    (longfirst ? omgstr : betstr) + " " + (longfirst ? betstr : omgstr);
}

std::string AzimuthString(ang alp, int prec, bool dms, char dmssep) {
  using namespace GeographicLib;
  return dms ? DMS::Encode(real(alp), prec + 5, DMS::AZIMUTH, dmssep) :
    DMS::Encode(real(alp), prec + 5, DMS::NUMBER);
}

int usage(int retval, bool /*brief*/) { return retval; }

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    typedef GeographicLib::Angle ang;
    Utility::set_digits();
    bool inverse = false,
      dms = false, full = false, unroll = true,
      longfirst = false,
      lineseq = false, cart = false, bench = false, reverse = false,
      debug = false;
    real
      a = 6378172, b = 6378102, c = 6356752;
    ang bet1, omg1, alp1, bet2, omg2, alp2;
    real s12;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);
    real ds = 0;
    long nmin = 0, nmax = 0;

    Triaxial t(a, b, c);
    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-i") {
        inverse = true;
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
        t = Triaxial(a, b, c);
        m += 3;
      } else if (arg == "-e") {
        // Cayley ellipsoid sqrt([2,1,1/2]) is
        // -e 1 3/2 1/3 2/3
        if (m + 4 >= argc) return usage(1, true);
        real e2, k2, kp2;
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
        t = Triaxial(b, e2, k2, kp2);
        a = t.a(); c = t.c();
        m += 4;
      } else if (arg == "-s") {
        if (m + 3 >= argc) return usage(1, true);
        try {
          ds = Utility::val<real>(std::string(argv[m + 1]));
          nmin = Utility::val<long>(std::string(argv[m + 2]));
          nmax = Utility::val<long>(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -s: " << e.what() << "\n";
          return 1;
        }
        if (nmin > 0 || nmax < 0) {
          std::cerr << "Bad values of nmin or nmax\n";
          return 1;
        }
        m += 3;
        lineseq = true;
      } else if (arg == "--cart")
        cart = true;
      else if (arg == "-u")
        unroll = !unroll;
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
      else if (arg == "--bench")
        bench = true;
      else if (arg == "-r")
        reverse = true;
      else if (arg == "--debug")
        debug = true;
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

    real errmult = 1000000;
    if (!ifile.empty() && !istring.empty()) {
      std::cerr << "Cannot specify --input-string and --input-file together\n";
      return 1;
    }
    if (inverse && cart) {
      std::cerr << "Cannot specify -i and --cart together\n";
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

    t.debug(debug);
    using std::round; using std::log10; using std::fabs;
    int disprec = int(round(log10(6400000/b)));
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
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
// Need to specify modes of operation here!
// bet1 omg1 bet2 omg2 -> alp1 alp2 s12 (Jac)  -i
// bet1 omg1 bet2 omg2 -> bet1 omg1 alp1 bet1 omg1 alp2 s12 m12 M12 M21
//                        (Jac and ODE)  -f
// bet1 omg1 alp1 s12 -> bet2 omg2 alp2 (Jac or cart)
// bet1 omg1 alp1 ds nmin nmax -> bet2 omg2 alp2 ... (Jac or ODE) -s, -s --cart
// bet1 omg1 alp1 ds nmin nmax -> bet2 omg2 alp2 dr dv (Jac and ODE + bench)
//                         -s --bench
// Full line -> inverse error dr Jac -i --bench
// Full line -> forward direct error dr dv Jac or ODE --bench, --bench --cart
// Full line -> reverse direct error dr dv Jav or ODE
//                         --bench -r, --bench --cart -r


        if (!bench && (inverse || full)) {
          if (!(str >> sbet1 >> somg1 >> sbet2 >> somg2))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
          DecodeLatLon(sbet2, somg2, bet2, omg2, longfirst);
          (void) t.Inverse(bet1, omg1, bet2, omg2, alp1, alp2, s12);
          if (full) {
            TriaxialODE l(t, bet1, omg1, alp1);
            ang bet2a, omg2a, alp2a;
            real m12, M12, M21;
            (void) l.Position(s12, bet2a, omg2a, alp2a, m12, M12, M21);
            // Strip trailing 0's and convert -180 to 180 with
            // sed -e 's/\.\([0-9]*[1-9]\)0*\b/.\1/g'
            //     -e 's/\.00*\b//g'
            //     -e 's/ -180 / 180 /g'
            *output << BetOmgString(bet1, omg1, prec, dms, dmssep, longfirst)
                    << " " << AzimuthString(alp1, prec, dms, dmssep)
                    << " "
                    << BetOmgString(bet2, omg2, prec, dms, dmssep, longfirst)
                    << " " << AzimuthString(alp2, prec, dms, dmssep)
                    << " "
                    << Utility::str(s12, prec + disprec) << " "
                    << Utility::str(m12, prec + disprec) << " "
                    << Utility::str(M12, prec+7) << " "
                    << Utility::str(M21, prec+7) << eol;
          } else
            *output << AzimuthString(alp1, prec, dms, dmssep) << " "
                    << AzimuthString(alp2, prec, dms, dmssep) << " "
                    << Utility::str(s12, prec + disprec) << eol;
        } else if (lineseq) {
          if (!(str >> sbet1 >> somg1 >> salp1))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
          alp1 = DecodeAzimuth(salp1);
          std::vector<ang> bet2v, omg2v, alp2v, bet2w, omg2w, alp2w;
          if (bench || !cart) {
            int m = nmax - nmin + 1;
            bet2v.resize(m); omg2v.resize(m); alp2v.resize(m);
            TriaxialLine l(t, bet1, omg1, alp1);
            for (int i = nmin, k = 0; i <= nmax; ++i, ++k)
              l.Position(i * ds, bet2v[k], omg2v[k], alp2v[k]);
            if (!unroll) {
              for (int k = 0; k <= nmax - nmin; ++k)
                Triaxial::AngNorm(bet2v[k], omg2v[k], alp2v[k]);
            }
          }
          if (bench || cart) {
            TriaxialODE l(t, bet1, omg1, alp1);
            l.Position(ds, nmin, nmax, bet2w, omg2w, alp2w);
          }
          for (int k = 0; k <= nmax - nmin; ++k) {
            if (!bench) {
              if (cart)
                *output << BetOmgString(bet2w[k], omg2w[k], prec, dms,
                                        dmssep, longfirst) << " "
                        << AzimuthString(alp2w[k], prec, dms, dmssep) << eol;
              else
                *output << BetOmgString(bet2v[k], omg2v[k], prec, dms,
                                        dmssep, longfirst) << " "
                        << AzimuthString(alp2v[k], prec, dms, dmssep) << eol;
            } else {
              std::pair<real, real> diff =
                t.EuclideanDiff(bet2v[k], omg2v[k], alp2v[k],
                                bet2w[k], omg2w[k], alp2w[k]);
              *output << BetOmgString(bet2v[k], omg2v[k], prec, dms,
                                      dmssep, longfirst) << " "
                      << AzimuthString(alp2v[k], prec, dms, dmssep) << " "
                      << Utility::str(diff.first * errmult) << " "
                      << Utility::str(diff.second * errmult) << eol;
            }
          }
        } else if (bench) {
          if (!(str >> sbet1 >> somg1 >> salp1
                >> sbet2 >> somg2 >> salp2 >> ss12))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
          alp1 = DecodeAzimuth(salp1);
          DecodeLatLon(sbet2, somg2, bet2, omg2, longfirst);
          alp2 = DecodeAzimuth(salp2);
          s12 = Utility::val<real>(ss12);
          if (inverse) {
            real s12a;
            t.Inverse(bet1, omg1, bet2, omg2, alp1, alp2, s12a);
            *output << Utility::str(fabs(s12 - s12a)) << eol;
          } else if (reverse) {
            ang bet1a, omg1a, alp1a;
            if (cart) {
              TriaxialODE l(t, bet2, omg2, alp2);
              l.Position(-s12, bet1a, omg1a, alp1a);
            } else
              t.Direct(bet2, omg2, alp2, -s12, bet1a, omg1a, alp1a);
            std::pair<real, real> diff =
              t.EuclideanDiff(bet1a, omg1a, alp1a, bet1, omg1, alp1);
            *output << Utility::str(diff.first * errmult) << " "
                    << Utility::str(diff.second * errmult) << eol;
          } else {
            ang bet2a, omg2a, alp2a;
            if (cart) {
              TriaxialODE l(t, bet1, omg1, alp1);
              l.Position(s12, bet2a, omg2a, alp2a);
            } else
              t.Direct(bet1, omg1, alp1, s12, bet2a, omg2a, alp2a);
            std::pair<real, real> diff =
              t.EuclideanDiff(bet2a, omg2a, alp2a, bet2, omg2, alp2);
            *output << Utility::str(diff.first * errmult) << " "
                    << Utility::str(diff.second * errmult) << eol;
          }
        } else {
          if (!(str >> sbet1 >> somg1 >> salp1 >> ss12))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DecodeLatLon(sbet1, somg1, bet1, omg1, longfirst);
          alp1 = DecodeAzimuth(salp1);
          s12 = Utility::val<real>(ss12);
          if (cart) {
              TriaxialODE l(t, bet1, omg1, alp1);
              l.Position(s12, bet2, omg2, alp2);
          } else
            t.Direct(bet1, omg1, alp1, s12, bet2, omg2, alp2);
          *output << BetOmgString(bet2, omg2, prec, dms,
                                  dmssep, longfirst) << " "
                  << AzimuthString(alp2, prec, dms, dmssep) << eol;
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

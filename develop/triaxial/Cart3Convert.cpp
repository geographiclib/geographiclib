// Counterpart of CartConvert

// Usual flags
//   -r (x, y, z to ellipsoidal to geographic)
//   -e (supplement with -t)
//   -w longfirt
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
#include "Angle.hpp"
#include "Triaxial.hpp"

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
    Utility::set_digits();
    typedef Triaxial::vec3 vec3;
    enum { GEODETIC, PARAMETRIC, GEOCENTRIC, ELLIPSOIDAL };
    int mode = ELLIPSOIDAL;
    bool threed = false, direction = false, reverse = false, dms = false,
      longfirst = false;
    real
      a = 6378172, b = 6378102, c = 6356752;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    Triaxial t(a, b, c);
    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-E")
        mode = ELLIPSOIDAL;
      else if (arg == "-G")
        mode = GEODETIC;
      else if (arg == "-P")
        mode = PARAMETRIC;
      else if (arg == "-C")
        mode = GEOCENTRIC;
      else if (arg == "-3")
        threed = true;
      else if (arg == "-D")
        direction = true;
      else if (arg == "-r")
        reverse = true;
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

    if (direction && mode != ELLIPSOIDAL) {
      std::cerr << "Can only specify -D with ellipsoidal conversions\n";
      return 1;
    }
    if (threed && !(mode == ELLIPSOIDAL || mode == GEODETIC)) {
      std::cerr
        << "Can only specify -3 with ellipsoidal or geodetic conversions\n";
      return 1;
    }
    if (direction && threed) {
      std::cerr
        << "Cannot specify both -3 and -D\n";
      return 1;
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

    using std::round, std::log10, std::ceil, std::signbit;
    int disprec = int(round(log10(6400000/b)));
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
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
          DecodeLatLon(sbet, somg, bet, omg, longfirst);
          if (mode != ELLIPSOIDAL && (bet.n() != 0 || signbit(bet.c())))
            throw GeographicErr("Latitude outside range [-90,90]: " + s);
          if (threed) {
            if (!(str >> sh))
              throw GeographicErr("Incomplete input: " + s);
            h = Utility::val<real>(sh);
          } else if (direction) {
            if (!(str >> salp))
              throw GeographicErr("Incomplete input: " + s);
            alp = DecodeAzimuth(salp);
          }
        }
        if (str >> strc)
          throw GeographicErr("Extraneous input: " + strc);
        // PROCESS
        switch (mode) {
        case ELLIPSOIDAL:
          if (reverse) {
            if (threed)
              t.carttoellip(r, bet, omg, h);
            else if (direction)
              t.cart2toellip(r, v, bet, omg, alp);
            else
              t.cart2toellip(r, bet, omg);
          } else {
            if (threed)
              t.elliptocart(bet, omg, h, r);
            else if (direction)
              t.elliptocart2(bet, omg, alp, r, v);
            else
              t.elliptocart2(bet, omg, r);
          }
          break;
        case GEODETIC:
          if (reverse) {
            if (threed)
              t.carttogeod(r, bet, omg, h);
            else
              t.cart2togeod(r, bet, omg);
          } else {
            if (threed)
              t.geodtocart(bet, omg, h, r);
            else
              t.geodtocart2(bet, omg, r);
          }
          break;
        case PARAMETRIC:
          if (reverse)
            t.cart2toparam(r, bet, omg);
          else
            t.paramtocart2(bet, omg, r);
          break;
        case GEOCENTRIC:
        default:
          if (reverse)
            t.cart2togeocen(r, bet, omg);
          else
            t.geocentocart2(bet, omg, r);
          break;
        }
        // WRITE
        if (reverse) {
          *output << BetOmgString(bet, omg, prec, dms, dmssep, longfirst);
          if (threed)
            *output << " " << Utility::str(h, prec + disprec);
          else if (direction)
            *output << " " << AzimuthString(alp, prec, dms, dmssep);
        } else {
          *output << Utility::str(r[0], prec + disprec) << " "
                  << Utility::str(r[1], prec + disprec) << " "
                  << Utility::str(r[2], prec + disprec);
          if (direction)
          *output << " "
                  << Utility::str(v[0], prec + disprec) << " "
                  << Utility::str(v[1], prec + disprec) << " "
                  << Utility::str(v[2], prec + disprec);
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

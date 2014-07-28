/**
 * \file RhumbSolve.cpp
 * \brief Command line utility for rhumb line calculations
 *
 * Copyright (c) Charles Karney (2014) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
#include <GeographicLib/Ellipsoid.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions and potentially
// uninitialized local variables
#  pragma warning (disable: 4127 4701)
#endif

#include "RhumbSolve.usage"

using namespace GeographicLib;
typedef Math::real real;

class Rhumb;

class RhumbLine {
private:
  const Ellipsoid& _ell;
  bool _exact;
  real _lat1, _lon1, _salp, _calp, _mu1, _psi1, _r1;
  RhumbLine& operator=(const RhumbLine&); // copy assignment not allowed
  friend Rhumb;
  // Threshold for using dmu/dpsi instead of mu12/psi12.
  static inline real tol() {
    static const real
      tol = 90 * Math::cbrt(std::numeric_limits<real>::epsilon());
    return tol;
  }
  static const int tm_maxord = GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER;
  static inline real overflow() {
    // Overflow value s.t. atan(overflow_) = pi/2
    static const real
      overflow = 1 / Math::sq(std::numeric_limits<real>::epsilon());
    return overflow;
  }
  static inline real tano(real x) {
    using std::abs; using std::tan;
    return
      2 * abs(x) == Math::pi() ? (x < 0 ? - overflow() : overflow()) :
      tan(x);
  }
  static inline real sinc(real x) { using std::sin; return x ? sin(x) / x : 1; }

  // Use divided differences to determine (psi2 - psi1) / (mu2 - mu1)
  // accurately
  //
  // Definition: Df(x,y,d) = (f(x)-f(y))/d with d=x-y
  // See:
  //   W. M. Kahan and R. J. Fateman,
  //   Symbolic computation of divided differences,
  //   SIGSAM Bull. 33(3), 7-28 (1999)
  //   http://dx.doi.org/10.1145/334714.334716
  //   http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf

  static inline real Dtan(real x, real y, real d) {
    real tx = tano(x), ty = tano(y), txy = tx * ty;
    return d ? (2 * txy > -1 ? (1 + txy) * tano(d) : tx - ty) / d :
      1 + txy;
  }
  static inline real Dasinh(real x, real y, real d) {
    real hx = Math::hypot(real(1), x), hy = Math::hypot(real(1), y);
    return d ? Math::asinh(x*y > 0 ? d * (x + y) / (x*hy + y*hx) :
                           x*hy - y*hx) / d :
      1 / hx;
  }
  static inline real Dgdinv(real x, real y, real d) {
    real tx = tano(x), ty = tano(y), dt = Dtan(x, y, d);
    return Dasinh(tx, ty, dt * d) * dt;
  }
  real DRectifyingToConformal(real x, real y, real d) const {
    real s = 0, p = x + y;
    for (int j = tm_maxord; j; --j)
      s += j * _ell.RectifyingToConformalCoeffs()[j] * cos(j * p) * sinc(j * d);
    return 1 - 2 * s;
  }
  real DRectifyingToIsometric(real x, real y, real d) const {
    real
      chix = _ell.ConformalLatitude
      (_ell.InverseRectifyingLatitude(x/Math::degree())) * Math::degree(),
      chiy = _ell.ConformalLatitude
      (_ell.InverseRectifyingLatitude(y/Math::degree())) * Math::degree(),
      dchi = DRectifyingToConformal(x, y, d);
    return Dgdinv(chix, chiy, dchi * d) * dchi;
  }
public:
  RhumbLine(const Ellipsoid& ell, real lat1, real lon1, real azi,
            bool exact = false)
    : _ell(ell)
    , _exact(exact)
    , _lat1(lat1)
    , _lon1(Math::AngNormalize(lon1))
  {
    using std::abs; using std::sin; using std::cos;
    azi = Math::AngNormalize(azi);
    real alp = azi * Math::degree();
    _salp =     azi  == -180 ? 0 : sin(alp);
    _calp = abs(azi) ==   90 ? 0 : cos(alp);
    _mu1 = _ell.RectifyingLatitude(lat1);
    _psi1 = _ell.IsometricLatitude(lat1);
    _r1 = _ell.CircleRadius(lat1);
  }
  void Position(real s12, real& lat2, real& lon2) const {
    using std::abs;
    real
      mu12 = s12 * _calp * 90 / _ell.QuarterMeridian(),
      mu2 = _mu1 + mu12;
    if (abs(mu2) <= 90) {
      if (_calp) {
        lat2 = _ell.InverseRectifyingLatitude(mu2);
        real psi12;
        if (_exact)
          psi12 = !(abs(mu12) <= tol()) ?
            _ell.IsometricLatitude(lat2) - _psi1 :
            // use dpsi/dmu = (2/pi)*Q/R
            mu12 * 4 * _ell.QuarterMeridian() /
            (Math::pi() * (_r1 + _ell.CircleRadius(lat2)));
        else
          psi12 = DRectifyingToIsometric(mu2 * Math::degree(),
                                         _mu1 * Math::degree(),
                                         mu12 * Math::degree()) * mu12;
        lon2 = _salp * psi12 / _calp;
      } else {
        lat2 = _lat1;
        lon2 = _salp * s12 / (_r1 * Math::degree());
      }
      lon2 = Math::AngNormalize2(_lon1 + lon2);
    } else {
      // Reduce to the interval [-180, 180)
      mu2 = Math::AngNormalize2(mu2);
      // Deal with points on the anti-meridian
      if (abs(mu2) > 90) mu2 = Math::AngNormalize(180 - mu2);
      lat2 = _ell.InverseRectifyingLatitude(mu2);
      lon2 = Math::NaN();
    }
  }
};

class Rhumb {
private:
  Ellipsoid _ell;
  bool _exact;
  static inline real gd(real x)
  { using std::atan; using std::sinh; return atan(sinh(x)); }

  // Use divided differences to determine (mu2 - mu1) / (psi2 - psi1)
  // accurately
  //
  // Definition: Df(x,y,d) = (f(x)-f(y))/d with d=x-y
  // See:
  //   W. M. Kahan and R. J. Fateman,
  //   Symbolic computation of divided differences,
  //   SIGSAM Bull. 33(3), 7-28 (1999)
  //   http://dx.doi.org/10.1145/334714.334716
  //   http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf

  static inline real Datan(real x, real y, real d) {
    using std::atan;
    real xy = x * y;
    return d ? (2 * xy > -1 ? atan( d / (1 + xy) ) : atan(x) - atan(y)) / d :
      1 / (1 + xy);
  }
  static inline real Dsinh(real x, real y, real d) {
    using std::sinh; using std::cosh;
    return cosh((x + y)/2) * (d ? 2 * sinh(d/2) / d : 1);
  }
  static inline real Dgd(real x, real y, real d) {
    using std::sinh;
    real sx = sinh(x), sy = sinh(y), ds = Dsinh(x, y, d);
    return Datan(sx, sy, ds * d) * ds;
  }
  real DConformalToRectifying(real x, real y, real d) const {
    real s = 0, p = x + y;
    for (int j = RhumbLine::tm_maxord; j; --j)
      s += j * _ell.ConformalToRectifyingCoeffs()[j] * cos(j * p) *
        RhumbLine::sinc(j * d);
    return 1 + 2 * s;
  }
  real DIsometricToRectifying(real x, real y, real d) const {
    real chix = gd(x), chiy = gd(y), dchi = Dgd(x, y, d);
    return DConformalToRectifying(chix, chiy, dchi * d) * dchi;
  }
public:
  Rhumb(real a, real f, bool exact = false) : _ell(a, f), _exact(exact) {}
  void Inverse(real lat1, real lon1, real lat2, real lon2,
               real& s12, real& azi) const {
    using std::atan2; using std::abs;
    real
      lon12 = Math::AngDiff(Math::AngNormalize(lon1), Math::AngNormalize(lon2)),
      psi1 = _ell.IsometricLatitude(lat1),
      psi2 = _ell.IsometricLatitude(lat2),
      psi12 = psi2 - psi1,
      h = Math::hypot(lon12, psi12);
    azi = 0 - atan2(-lon12, psi12) / Math::degree();
    lon12 = abs(lon12);
    real dmudpsi;
    if (_exact) {
      real mu12 = _ell.RectifyingLatitude(lat2) - _ell.RectifyingLatitude(lat1);
      dmudpsi = !(abs(mu12) <= RhumbLine::tol()) ? mu12 / psi12 :
        // use dmu/dpsi = (pi/2)*R/Q
        Math::pi() * (_ell.CircleRadius(lat1) + _ell.CircleRadius(lat2)) /
        (4 * _ell.QuarterMeridian());
    } else
      dmudpsi = DIsometricToRectifying(psi2 * Math::degree(),
                                       psi1 * Math::degree(),
                                       psi12 * Math::degree());
    s12 = h * dmudpsi * _ell.QuarterMeridian() / 90;
  }

  RhumbLine Line(real lat1, real lon1, real azi) const
  { return RhumbLine(_ell, lat1, lon1, azi, _exact); }

  void Direct(real lat1, real lon1, real azi, real s12,
              real& lat2, real& lon2) const
  { Line(lat1, lon1, azi).Position(s12, lat2, lon2); }
};

std::string LatLonString(real lat, real lon, int prec, bool dms, char dmssep) {
  return dms ?
    DMS::Encode(lat, prec + 5, DMS::LATITUDE, dmssep) + " " +
    DMS::Encode(lon, prec + 5, DMS::LONGITUDE, dmssep) :
    DMS::Encode(lat, prec + 5, DMS::NUMBER) + " " +
    DMS::Encode(lon, prec + 5, DMS::NUMBER);
}

std::string AzimuthString(real azi, int prec, bool dms, char dmssep) {
  return dms ? DMS::Encode(azi, prec + 5, DMS::AZIMUTH, dmssep) :
    DMS::Encode(azi >= 180 ? azi - 360 : azi, prec + 5, DMS::NUMBER);
}

int main(int argc, char* argv[]) {
  try {
    Utility::set_digits();
    bool linecalc = false, inverse = false, dms = false, exact = false;
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    real lat1, lon1, azi = Math::NaN(), lat2, lon2, s12;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-i") {
        inverse = true;
        linecalc = false;
      } else if (arg == "-l") {
        inverse = false;
        linecalc = true;
        if (m + 3 >= argc) return usage(1, true);
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            lat1, lon1);
          azi = DMS::DecodeAzimuth(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -l: " << e.what() << "\n";
          return 1;
        }
        m += 3;
      } else if (arg == "-e") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          a = Utility::num<real>(std::string(argv[m + 1]));
          f = Utility::fract<real>(std::string(argv[m + 2]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        m += 2;
      }
      else if (arg == "-d") {
        dms = true;
        dmssep = '\0';
      } else if (arg == "-:") {
        dms = true;
        dmssep = ':';
      } else if (arg == "-p") {
        if (++m == argc) return usage(1, true);
        try {
          prec = Utility::num<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-E")
        exact = true;
      else if (arg == "--input-string") {
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
        std::cout
          << argv[0] << ": GeographicLib version "
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

    const Rhumb rh(a, f, exact);
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    int retval = 0;
    std::string s;
    if (linecalc) {
      const RhumbLine rhl(rh.Line(lat1, lon1, azi));
      while (std::getline(*input, s)) {
        try {
          std::istringstream str(s);
          if (!(str >> s12))
            throw GeographicErr("Incomplete input: " + s);
          std::string strc;
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          rhl.Position(s12, lat2, lon2);
          *output << LatLonString(lat2, lon2, prec, dms, dmssep) << "\n";
        }
        catch (const std::exception& e) {
          // Write error message cout so output lines match input lines
          *output << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    } else if (inverse) {
      while (std::getline(*input, s)) {
        try {
          std::istringstream str(s);
          std::string slat1, slon1, slat2, slon2;
          if (!(str >> slat1 >> slon1 >> slat2 >> slon2))
            throw GeographicErr("Incomplete input: " + s);
          std::string strc;
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DMS::DecodeLatLon(slat1, slon1, lat1, lon1);
          DMS::DecodeLatLon(slat2, slon2, lat2, lon2);
          rh.Inverse(lat1, lon1, lat2, lon2, s12, azi);
          *output << AzimuthString(azi, prec, dms, dmssep) << " "
                  << Utility::str(s12, prec) << "\n";
        }
        catch (const std::exception& e) {
          // Write error message cout so output lines match input lines
          *output << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    } else {
      while (std::getline(*input, s)) {
        try {
          std::istringstream str(s);
          std::string slat1, slon1, sazi;
          if (!(str >> slat1 >> slon1 >> sazi >> s12))
            throw GeographicErr("Incomplete input: " + s);
          std::string strc;
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DMS::DecodeLatLon(slat1, slon1, lat1, lon1);
          azi = DMS::DecodeAzimuth(sazi);
          rh.Direct(lat1, lon1, azi, s12, lat2, lon2);
          *output << LatLonString(lat2, lon2, prec, dms, dmssep) << "\n";
        }
        catch (const std::exception& e) {
          // Write error message cout so output lines match input lines
          *output << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
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

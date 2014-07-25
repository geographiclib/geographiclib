/**
 * \file RhumbSolve.cpp
 * \brief Command line utility for rhumb line calculations
 *
 * Copyright (c) Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/Ellipsoid.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions and potentially
// uninitialized local variables
#  pragma warning (disable: 4127 4701)
#endif

int usage(int retval, bool) {
  ( retval ? std::cerr : std::cout ) << "Usage\n"
"    RhumbSolve [ -i | -l lat1 lon1 azi ] [ -e a f ] [ -d | -: ]\n"
"    [ -p prec ] [ -h | --help ]\n"
"\n";
  return retval;
}

using namespace GeographicLib;
typedef Math::real real;

class RhumbLine {
private:
  const Ellipsoid& _ell;
  real _lat1, _lon1, _salp, _calp, _mu1, _psi1, _r1;
public:
  RhumbLine(const Ellipsoid& ell, real lat1, real lon1, real azi)
    : _ell(ell)
    , _lat1(lat1)
    , _lon1(Math::AngNormalize(lon1))
  {
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
    real mu2 = _mu1 + s12 * _calp * 90 / _ell.QuarterMeridian();
    if (abs(mu2) <= 90) {
      if (_calp) {
        lat2 = _ell.InverseRectifyingLatitude(mu2);
        real psi12 = _ell.IsometricLatitude(lat2) - _psi1;
        lon2 = Math::AngNormalize(_lon1 + _salp * psi12 / _calp);
      } else {
        lat2 = _lat1;
        lon2 = Math::AngNormalize(_lon1 + s12 / (_r1 * Math::degree()));
      }
    } else
      lat2 = lon2 = Math::NaN();
  }
};    

class Rhumb {
private:
  const Ellipsoid _ell;
  friend class RhumbLine;
public:
  Rhumb(real a, real f)
    : _ell(a, f)
  {}
  void Inverse(real lat1, real lon1, real lat2, real lon2,
               real& s12, real& azi) const {
    using std::atan2;
    real
      lon12 = Math::AngDiff(Math::AngNormalize(lon1),
                            Math::AngNormalize(lon2)),
      psi12 = _ell.IsometricLatitude(lat2) - _ell.IsometricLatitude(lat1),
      mu12 =  _ell.RectifyingLatitude(lat2) - _ell.RectifyingLatitude(lat1);
    azi = 0 - atan2(-lon12, psi12) / Math::degree();
    s12 = psi12 ?
      mu12 * (Math::hypot(lon12, psi12) / psi12) * _ell.QuarterMeridian() / 90 :
      (lon12 * Math::degree()) * _ell.CircleRadius(lat1);
  }

  RhumbLine Line(real lat1, real lon1, real azi) const
  { return RhumbLine(_ell, lat1, lon1, azi); }

  void Direct(real lat1, real lon1, real azi, real s12,
              real& lat2, real& lon2) const
  { Line(lat1, lon1, azi).Position(s12, lat2, lon2); }
};

    
std::string LatLonString(real lat, real lon, int prec, bool dms, char dmssep) {
  using namespace GeographicLib;
  return dms ?
    DMS::Encode(lat, prec + 5, DMS::LATITUDE, dmssep) + " " +
    DMS::Encode(lon, prec + 5, DMS::LONGITUDE, dmssep) :
    DMS::Encode(lat, prec + 5, DMS::NUMBER) + " " +
    DMS::Encode(lon, prec + 5, DMS::NUMBER);
}

std::string AzimuthString(real azi, int prec, bool dms, char dmssep) {
  using namespace GeographicLib;
  return dms ? DMS::Encode(azi, prec + 5, DMS::AZIMUTH, dmssep) :
    DMS::Encode(azi >= 180 ? azi - 360 : azi, prec + 5, DMS::NUMBER);
}

int main(int argc, char* argv[]) {
  try {
    using namespace GeographicLib;
    Utility::set_digits();
    bool linecalc = false, inverse = false, dms = false;
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    real lat1, lon1, azi, lat2, lon2, s12;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char dmssep = char(0);

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
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
    }

    const Rhumb rh(a, f);
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    int retval = 0;
    std::string s;
    if (linecalc) {
      const RhumbLine rhl(rh.Line(lat1, lon1, azi));
      while (std::getline(std::cin, s)) {
        try {
          std::istringstream str(s);
          if (!(str >> s12))
            throw GeographicErr("Incomplete input: " + s);
          std::string strc;
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          rhl.Position(s12, lat2, lon2);
          std::cout << LatLonString(lat2, lon2, prec, dms, dmssep) << "\n";
        }
        catch (const std::exception& e) {
          // Write error message cout so output lines match input lines
          std::cout << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    } else if (inverse) {
      while (std::getline(std::cin, s)) {
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
          std::cout << AzimuthString(azi, prec, dms, dmssep) << " "
                    << Utility::str(s12, prec) << "\n";
        }
        catch (const std::exception& e) {
          // Write error message cout so output lines match input lines
          std::cout << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    } else {
      while (std::getline(std::cin, s)) {
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
          std::cout << LatLonString(lat2, lon2, prec, dms, dmssep) << "\n";
        }
        catch (const std::exception& e) {
          // Write error message cout so output lines match input lines
          std::cout << "ERROR: " << e.what() << "\n";
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

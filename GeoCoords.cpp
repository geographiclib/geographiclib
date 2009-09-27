/**
 * \file GeoCoords.cpp
 * \brief Implementation for GeographicLib::GeoCoords class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/GeoCoords.hpp"
#include "GeographicLib/MGRS.hpp"
#include "GeographicLib/DMS.hpp"
#include <vector>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cerrno>

#define GEOGRAPHICLIB_GEOCOORDS_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GEOCOORDS_CPP)
RCSID_DECL(GEOGRAPHICLIB_GEOCOORDS_HPP)

namespace GeographicLib {

  using namespace std;

  void GeoCoords::Reset(const std::string& s, bool centerp) {
    vector<string> sa;
    bool in = false;
    for (unsigned i = 0; i < s.size(); ++i) {
      if (isspace(s[i]) || s[i] == ',') {
        in = false;
        continue;
      }
      if (!in)
        sa.push_back("");
      in = true;
      sa.back().push_back(s[i]);
    }
    if (sa.size() == 1) {
      int prec;
      MGRS::Reverse(sa[0], _zone, _northp, _easting, _northing, prec, centerp);
      UTMUPS::Reverse(_zone, _northp, _easting, _northing,
                      _lat, _long, _gamma, _k);
    } else if (sa.size() == 2) {
      DMS::DecodeLatLon(sa[0], sa[1], _lat, _long);
      UTMUPS::Forward( _lat, _long,
                       _zone, _northp, _easting, _northing, _gamma, _k);
    } else if (sa.size() == 3) {
      unsigned zoneind, coordind;
      if (sa[0].size() > 0 && isalpha(sa[0][sa[0].size() - 1])) {
        zoneind = 0;
        coordind = 1;
      } else if (sa[2].size() > 0 && isalpha(sa[2][sa[2].size() - 1])) {
        zoneind = 2;
        coordind = 0;
      } else
        throw out_of_range("Neither " + sa[0] + " nor " + sa[2] +
                           " of the form UTM/UPS Zone + Hemisphere" +
                           " (ex: 38N, 09S, N)");
      UTMUPS::DecodeZone(sa[zoneind], _zone, _northp);
      for (unsigned i = 0; i < 2; ++i) {
        istringstream str(sa[coordind + i]);
        real_t x;
        if (!(str >> x))
          throw out_of_range("Bad number " + sa[coordind + i]);
        if (int(str.tellg()) != int(sa[coordind + i].size()))
          throw out_of_range(string("Extra text in UTM/UPS ") +
                             (i == 0 ? "easting " : "northing ") +
                             sa[coordind + i]);
        (i ? _northing : _easting) = x;
      }
      UTMUPS::Reverse(_zone, _northp, _easting, _northing,
                      _lat, _long, _gamma, _k);
      FixHemisphere();
    } else
      throw out_of_range("Coordinate requires 1, 2, or 3 elements");
    CopyToAlt();
  }


  string GeoCoords::GeoRepresentation(int prec) const {
    prec = max(0, min(9, prec) + 5);
    ostringstream os;
    os << fixed << setprecision(prec)
       << _lat << " " << _long;
    return os.str();
  }

  string GeoCoords::DMSRepresentation(int prec) const {
    prec = max(0, min(10, prec) + 5);
    return DMS::Encode(_lat, unsigned(prec), DMS::LATITUDE) +
      " " + DMS::Encode(_long, unsigned(prec), DMS::LONGITUDE);
  }

  string GeoCoords::MGRSRepresentation(int prec) const {
    // Max precision is um
    prec = max(0, min(6, prec) + 5);
    string mgrs;
    MGRS::Forward(_zone, _northp, _easting, _northing, _lat, prec, mgrs);
    return mgrs;
  }

  string GeoCoords::AltMGRSRepresentation(int prec) const {
    // Max precision is um
    prec = max(0, min(6, prec) + 5);
    string mgrs;
    MGRS::Forward(_alt_zone, _northp, _alt_easting, _alt_northing, _lat, prec,
                  mgrs);
    return mgrs;
  }

  void GeoCoords::UTMUPSString(int zone, real_t easting, real_t northing,
                               int prec, std::string& utm) const {
    ostringstream os;
    prec = max(-5, min(9, prec));
    real_t scale = prec < 0 ? pow(real_t(10), -prec) : real_t(1);
    os << UTMUPS::EncodeZone(zone, _northp) << fixed << setfill('0');
    os << " "
       << setprecision(max(0, prec))
       << easting / scale;
    if (prec < 0 && abs(easting / scale) > real_t(0.5L))
      os << setw(-prec) << 0;
    os << " "
       << setprecision(max(0, prec))
       << northing / scale;
    if (prec < 0 && abs(northing / scale) > real_t(0.5L))
      os << setw(-prec) << 0;
    utm = os.str();
  }

  string GeoCoords::UTMUPSRepresentation(int prec) const {
    string utm;
    UTMUPSString(_zone, _easting, _northing, prec, utm);
    return utm;
  }

  string GeoCoords::AltUTMUPSRepresentation(int prec) const {
    string utm;
    UTMUPSString(_alt_zone, _alt_easting, _alt_northing, prec, utm);
    return utm;
  }

  void GeoCoords::FixHemisphere() {
    if (_lat == 0 || (_northp && _lat >= 0) || (!_northp && _lat < 0))
      // Allow either hemisphere for equator
      return;
    if (_zone > 0) {
      _northing += (_northp ? 1 : -1) * UTMUPS::UTMShift();
      _northp = !_northp;
    } else
      throw out_of_range("Hemisphere mixup");
  }

} // namespace GeographicLib

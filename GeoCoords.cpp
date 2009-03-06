/**
 * \file GeoCoords.cpp
 * \brief Implementation for GeographicLib::GeoCoords class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include <vector>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cerrno>
#include "GeographicLib/GeoCoords.hpp"
#include "GeographicLib/MGRS.hpp"
#include "GeographicLib/DMS.hpp"

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = GEOCOORDS_HPP;
}

namespace GeographicLib {

  using namespace std;

  void GeoCoords::Reset(const std::string& s) {
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
      MGRS::Reverse(sa[0], _zone, _northp, _easting, _northing, prec);
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
      char hemi = toupper(sa[zoneind][sa[zoneind].size() - 1]);
      _northp = hemi == 'N';
      if (! (_northp || hemi == 'S'))
	throw out_of_range(string("Illegal hemisphere letter ") + hemi
			   + " in " + sa[zoneind]);
      const char* c = sa[zoneind].c_str();
      char* q;
      _zone = strtol(c, &q, 10);
      if (q - c != int(sa[zoneind].size()) - 1)
	throw out_of_range("Extra text in UTM/UPS zone " + sa[zoneind]);
      if (q > c && _zone == 0)
	// Don't allow 0N as an alternative to N for UPS coordinates
	throw out_of_range("Illegal zone 0 in " + sa[zoneind]);
      if (q == c)
	_zone = 0;
      for (unsigned i = 0; i < 2; ++i) {
	const char* c = sa[coordind + i].c_str();
	errno = 0;
	double x = strtod(c, &q);
	if (errno ==  ERANGE || !isfinite(x))
	  throw out_of_range("Number " + sa[coordind + i] + " out of range");
	if (q - c != int(sa[coordind + i].size()))
	  throw out_of_range(string("Extra text in UTM/UPS ") +
			     (i == 0 ? "easting " : "northing ") +
			     sa[coordind + i]);
	if (i == 0)
	  _easting = x;
	else
	  _northing = x;
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

  void GeoCoords::UTMUPSString(int zone, double easting, double northing,
			       int prec, std::string& utm) const {
    ostringstream os;
    os << fixed << setfill('0');
    if (zone)
      os << setw(2) << zone;
    prec = max(-5, min(9, prec));
    double scale = prec < 0 ? pow(10.0, -prec) : 1.0;
    os << (_northp ? 'N' : 'S') << " "
       << setprecision(max(0, prec))
       << easting / scale;
    if (prec < 0 && abs(easting / scale) > 0.5)
      os << setw(-prec) << 0;
    os << " "
       << setprecision(max(0, prec))
       << northing / scale;
    if (prec < 0 && abs(northing / scale) > 0.5)
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
    if (_lat == 0 || (_northp && _lat > 0) || (!_northp && _lat < 0))
      // Allow either hemisphere for equator
      return;
    if (_zone > 0) {
      _northing += (_northp ? 1 : -1) * MGRS::utmNshift;
      _northp = !_northp;
    } else
      throw out_of_range("Hemisphere mixup");
  }

} // namespace GeographicLib

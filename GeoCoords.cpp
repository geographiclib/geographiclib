/**
 * \file GeoCoords.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#include <vector>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cerrno>
#include "GeoCoords.hpp"
#include "MGRS.hpp"
#include "DMS.hpp"

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = GEOCOORDS_HPP;
}

#if defined(_MSC_VER)
#define isnan(x) false
#define isinf(x) false
#else
#define isnan(x) std::isnan(x)
#define isinf(x) std::isinf(x)
#endif

namespace GeographicLib {

  void GeoCoords::Reset(const std::string& s) {
    std::vector<std::string> sa;
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
      double a, b;
      DMS::flag ia, ib;
      a = DMS::Decode(sa[0], ia);
      b = DMS::Decode(sa[1], ib);
      if (ia == DMS::NONE && ib == DMS::NONE) {
	// Default to lat, long
	ia = DMS::LATITUDE;
	ib = DMS::LONGITUDE;
      } else if (ia == DMS::NONE)
	ia = DMS::flag(DMS::LATITUDE + DMS::LONGITUDE - ib);
      else if (ib == DMS::NONE)
	ib = DMS::flag(DMS::LATITUDE + DMS::LONGITUDE - ia);
      if (ia == ib)
	throw std::out_of_range("Both " + sa[0] + " and " + sa[1] +
				" interpreted as "
				+ (ia == DMS::LATITUDE ? "latitudes"
				   : "longitudes"));
      if (ia == DMS::LATITUDE) {
	_lat = a;
	_long = b;
      } else {
	_lat = b;
	_long = a;
      }
      UTMUPS::Forward(-1, _lat, _long,
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
	throw std::out_of_range("Neither " + sa[0] + " nor " + sa[2] +
				" of the form UTM/UPS Zone + Hemisphere" +
				" (ex: 38N, 09S, N)");
      char hemi = toupper(sa[zoneind][sa[zoneind].size() - 1]);
      _northp = hemi == 'N';
      if (! (_northp || hemi != 'S'))
	throw std::out_of_range(std::string("Illegal hemisphere letter ") + hemi
				+ " in " + sa[zoneind]);
      const char* c = sa[zoneind].c_str();
      char* q;
      _zone = std::strtol(c, &q, 10);
      if (q - c != int(sa[zoneind].size()) - 1)
	throw std::out_of_range("Extra text in UTM/UPS zone " + sa[zoneind]);
      if (q > c && _zone == 0)
	// Don't allow 0N as an alternative to N for UPS coordinates
	throw std::out_of_range("Illegal zone 0 in " + sa[zoneind]);
      if (q == c)
	_zone = 0;
      for (unsigned i = 0; i < 2; ++i) {
	const char* c = sa[coordind + i].c_str();
	errno = 0;
	double x = std::strtod(c, &q);
	if (errno ==  ERANGE || isnan(x) || isinf(x))
	  throw std::out_of_range("Number " + sa[coordind + i]
				  + " out of range");
	if (q - c != int(sa[coordind + i].size()))
	  throw std::out_of_range(std::string("Extra text in UTM/UPS ") +
				  (i == 0 ? "easting " : "northing ") +
				  sa[coordind + i]);
	if (i == 0)
	  _easting = x;
	else
	  _northing = x;
      }
      UTMUPS::Reverse(_zone, _northp, _easting, _northing,
		      _lat, _long, _gamma, _k);
    } else
      throw std::out_of_range("Coordinate requires 1, 2, or 3 elements");
    CopyToAlt();
  }


  std::string GeoCoords::GeoRepresentation(int prec) const {
    prec = (std::max)(0, (std::min)(10, prec) + 5);
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec)
       << _lat << " " << _long;
    return os.str();
  }

  std::string GeoCoords::DMSRepresentation(int prec) const {
    prec = (std::max)(0, (std::min)(10, prec) + 5);
    std::ostringstream os;
    DMS::component trailing = DMS::component((std::min)(prec/2, 2));
    prec -= 2 * trailing;
    os << DMS::Encode(_lat, trailing, unsigned(prec), DMS::LATITUDE) << " "
       << DMS::Encode(_long, trailing, unsigned(prec), DMS::LONGITUDE);
    return os.str();
  }

  std::string GeoCoords::MGRSRepresentation(int prec) const {
    // Max precision is um
    prec = (std::max)(0, (std::min)(6, prec) + 5);
    std::string mgrs;
    MGRS::Forward(_zone, _northp, _easting, _northing, _lat, prec, mgrs);
    return mgrs;
  }

  std::string GeoCoords::AltMGRSRepresentation(int prec) const {
    // Max precision is um
    prec = (std::max)(0, (std::min)(6, prec) + 5);
    std::string mgrs;
    MGRS::Forward(_alt_zone, _northp, _alt_easting, _alt_northing, _lat, prec,
		  mgrs);
    return mgrs;
  }

  void GeoCoords::UTMUPSString(int zone, double easting, double northing,
			    int prec, std::string& utm) const {
    std::ostringstream os;
    os << std::fixed << std::setfill('0');
    if (zone)
      os << std::setw(2) << zone;
    prec = (std::max)(-5, (std::min)(9, prec));
    double scale = prec < 0 ? std::pow(10.0, -prec) : 1.0;
    os << (_northp ? 'N' : 'S') << " "
       << std::setprecision((std::max)(0, prec))
       << easting / scale;
    if (prec < 0)
      os << std::setw(-prec) << 0;
    os << " "
       << std::setprecision((std::max)(0, prec))
       << northing / scale;
    if (prec < 0)
      os << std::setw(-prec) << 0;
    utm = os.str();
  }

  std::string GeoCoords::UTMUPSRepresentation(int prec) const {
    std::string utm;
    UTMUPSString(_zone, _easting, _northing, prec, utm);
    return utm;
  }

  std::string GeoCoords::AltUTMUPSRepresentation(int prec) const {
    std::string utm;
    UTMUPSString(_alt_zone, _alt_easting, _alt_northing, prec, utm);
    return utm;
  }

} // namespace GeographicLib

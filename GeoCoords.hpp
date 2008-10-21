/**
 * \file GeoCoords.cpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#ifndef GEOCOORDS_HPP
#define GEOCOORDS_HPP "$Id$"

#include <cmath>
#include <string>
#include <cstdlib>
#include "UTMUPS.hpp"

// Class to parse a coordinate position as a string and interpret it as
//
// MGRS reference
// latitude longitude
// UTM/UPS position

// The input string is broken into space (or comma) separated pieces and Basic
// decision on which format is based on number of components
//
// 1 = MGRS
// 2 = "Lat Long" or "Long Lat"
// 3 = "Zone Easting Northing" or "Easting Northing Zone"

// The following inputs are approximately the same (Ar Ramadi Bridge, Iraq)
//
// Latitude and Longitude
//     33.44      43.27
//     N33d26.4'  E43d16.2'
//     43d16'12"E 33d26'24"N
// MGRS
//     38SLC301
//     38SLC391014
//     38SLC3918701405
//     37SHT9708
// UTM
//     38N 339188 3701405
//     897039 3708229 37N

// Once the input string has been parsed, you can print the result out in any of
// the formats, decimal degrees, degrees minutes seconds, MGRS, UTM/UPS.

// MGRS parsing interprets the grid references as square area at the specified
// precision (1m, 10m, 100m, etc.).  The center of this square is then taken to
// be the precise position.  Thus:
//
// 38SMB           = 38N 450000 3650000
// 38SMB4484       = 38N 444500 3684500
// 38SMB44148470   = 38N 444145 3684705
//
// Similarly when a MGRS grid reference is reported the enclosing grid square
// is returned.  Thus 38N 444180 3684790 converted to a MGRS reference at
// precision -2 (100m) is 38SMB441847 and not 38SMB442848

// Latitude and Longitude parsing.  Latitude precedes longitude, unless a N, S,
// E, W hemisphere designator is used on one or both coordinates.  Thus
//
// 40 -75
// N40 W75
// -75 N40
// 75W 40N
// E-75 -40S
//
// are all the same position.  The coodinates may be given in decimal degrees,
// degrees and decimal minutes, degrees, minutes, seconds, etc.  Use d, ', and
// " to make off the degrees, minutes and seconds.  Thus
//
// 40d30'30" 40d30'30 40d30.5' 40d30.5 40.508333333
//
// all specify the same angle.  The leading sign applies to all components so
// -1d30 is -(1+30/60) = -1.5.

// UTM/UPS details.  For UTM zones (-80 <= Lat <= 84), the zone designator is
// made up of a zone number (for 1 to 60) and a hemisphere letter (N or S),
// e.g., 38N.  The latitude zone designer (C-M in the southern hemispher and
// N-X in the northern) should NOT be used.  (This is the third character of
// the MGRS reference.)  The zone for the poles (where UPS is employed) is a
// hemisphere letter by itself, i.e., N or S.

// Straddling zones (UTM/UTM or UTM/UPS).  If the input string is in UTM, UPS,
// or MGRS format then the zone, easting, and northing are saved (in addition
// to the latitude and longitude).  The same zone is then used for output.
// This allows the use of the same grid system either side of a zone border
// (UTM to UTM or UTM or UPS).  Use SetAltZone(zone) or SetAltZone() (to
// indicate the standard zone) followed by Alt{MGRS,UTM}Representation to
// report the result in another zone.  Use zone = 0 to specify the use of UPS.

namespace GeographicLib {

  class GeoCoords {
  private:
    double _lat, _long, _easting, _northing, _gamma, _k;
    bool _northp;
    int _zone;			// 0 = poles, -1 = undefined
    mutable double _alt_easting, _alt_northing, _alt_gamma, _alt_k;
    mutable int _alt_zone;

    void CopyToAlt() const {
      _alt_easting = _easting;
      _alt_northing = _northing;
      _alt_gamma = _gamma;
      _alt_k = _k;
      _alt_zone = _zone;
    }      
    void UTMUPSString(int zone, double easting, double northing,
		   int prec, std::string& utm) const;
  public:
    GeoCoords()
      // This is the N pole
      : _lat(90.0)
      , _long(0.0)
      , _easting(2000000.0)
      , _northing(2000000.0)
      , _northp(true)
      , _zone(0)
    { CopyToAlt();}
    // Specify a location as a 1-element, 2-element, or 3-element string.
    GeoCoords(const std::string& s) { Reset(s); }
    // Specify the location in terms of latitude and longitude.  Use zone
    // argument to force UTM/UPS representation to use a specified zone (zone =
    // 0 means UPS).  Omitted the third argument causes the standard zone to be
    // used.
    GeoCoords(double latitude, double longitude, int zone = -1) {
      Reset(latitude, longitude, zone);
    }
    // Specify the location in terms of UPS/UPS easting and northing.  Use zone
    // = 0 to indicate UPS.
    GeoCoords(int zone, bool northp, double easting, double northing) {
      Reset(zone, northp, easting, northing);
    }
    // Specify a location as a 1-element, 2-element, or 3-element string.
    void Reset(const std::string& s);
    // Specify the location in terms of latitude and longitude.  Use zone
    // argument to force UTM/UPS representation to use a specified zone (zone =
    // 0 means UPS).  Omitted the third argument causes the standard zone to be
    // used.
    void Reset(double latitude, double longitude, int zone = -1) {
      _lat = latitude;
      _long = longitude;
      UTMUPS::Forward(zone, _lat, _long,
		      _zone, _northp, _easting, _northing, _gamma, _k);
      CopyToAlt();
    }
    // Specify the location in terms of UPS/UPS easting and northing.  Use zone =
    // 0 to indicate UPS.
    void Reset(int zone, bool northp, double easting, double northing) {
      _zone = zone;
      _northp = northp;
      _easting = easting;
      _northing = northing;
      UTMUPS::Reverse(_zone, _northp, _easting, _northing,
		      _lat, _long, _gamma, _k);
      CopyToAlt();
    }
    double Latitude() const { return _lat; }
    double Longitude() const { return _long; }
    double Easting() const { return _easting; }
    double Northing() const { return _northing; }
    double Convergence() const { return _gamma; }
    double Scale() const { return _k; }
    bool Northp() const { return _northp; }
    char Hemisphere() const { return _northp ? 'N' : 'S'; }
    // This returns the zone corresponding to the input (return zone = 0 for
    // UPS)
    int Zone() const { return _zone; }
    // Set the zone (default is the standard zone) for the alternate
    // representation.
    void SetAltZone(int zone = -1) const {
      if (zone == _zone)
	CopyToAlt();
      else {
	bool northp;
	UTMUPS::Forward(zone, _lat, _long,
			_alt_zone, northp,
			_alt_easting, _alt_northing, _alt_gamma, _alt_k);
	if (_alt_zone == _zone)
	  CopyToAlt();
      }
    }
    int AltZone() const { return _alt_zone; }
    double AltEasting() const { return _alt_easting; }
    double AltNorthing() const { return _alt_northing; }
    double AltConvergence() const { return _alt_gamma; }
    double AltScale() const { return _alt_k; }

    // Latitude and longitude signed decimal degrees
    std::string GeoRepresentation(int prec = 0) const;
    // Latitude and longitude degrees, minutes, seconds + hemisphere
    std::string DMSRepresentation(int prec = 0) const;
    // Return MGRS string
    std::string MGRSRepresentation(int prec = 0) const;
    // Easting, northing, zone+hemisphere
    std::string UTMUPSRepresentation(int prec = 0) const; 
    // Report location using alternative zone
    std::string AltMGRSRepresentation(int prec = 0) const;
    std::string AltUTMUPSRepresentation(int prec = 0) const;
  };

} // namespace GeographicLib
#endif	// GEOCOORDS_HPP

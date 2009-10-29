/**
 * \file GeoCoords.hpp
 * \brief Header for GeographicLib::GeoCoords class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#ifndef GEOGRAPHICLIB_GEOCOORDS_HPP
#define GEOGRAPHICLIB_GEOCOORDS_HPP "$Id$"

#include "GeographicLib/UTMUPS.hpp"
#include "GeographicLib/Constants.hpp"
#include <string>

namespace GeographicLib {

  /**
   * \brief Conversion between geographic coordinates
   *
   * This class stores a geographic position which may be set via the
   * constructors or Reset via
   * - latitude and longitude
   * - UTM or UPS coordinates
   * - a string represention of these or an MGRS coordinate string
   *
   * The state consists of the latitude and longitude and the supplied UTM or
   * UPS coordinates (possibly derived from the MGRS coordinates).  If latitude
   * and longitude were given then the UTM/UPS coordinates follows the standard
   * conventions.
   *
   * The mutable state consists of the UTM or UPS coordinates for a alternate
   * zone.  A method SetAltZone is provided to set the alternate UPS/UTM zone.
   *
   * Methods are provided to return the geographic coordinates, the input UTM
   * or UPS coordinates (and associated meridian convergence and scale), or
   * alternate UTM or UPS coordinates (and their associated meridian
   * convergence and scale).
   *
   * Once the input string has been parsed, you can print the result out in any
   * of the formats, decimal degrees, degrees minutes seconds, MGRS, UTM/UPS.
   **********************************************************************/
  class GeoCoords {
  private:
    typedef Math::real real;
    real _lat, _long, _easting, _northing, _gamma, _k;
    bool _northp;
    int _zone;                  // See UTMUPS::zone
    mutable real _alt_easting, _alt_northing, _alt_gamma, _alt_k;
    mutable int _alt_zone;

    void CopyToAlt() const throw() {
      _alt_easting = _easting;
      _alt_northing = _northing;
      _alt_gamma = _gamma;
      _alt_k = _k;
      _alt_zone = _zone;
    }
    void UTMUPSString(int zone, real easting, real northing,
                      int prec, std::string& utm) const;
    void FixHemisphere();
  public:

    /**
     * The default contructor is equivalent to \e latitude = 90<sup>o</sup>, \e
     * longitude = 0<sup>o</sup>.
     **********************************************************************/
    GeoCoords() throw()
      // This is the N pole
      : _lat(90)
      , _long(0)
      , _easting(2000000)
      , _northing(2000000)
      , _northp(true)
      , _zone(0)
    { CopyToAlt();}

    /**
     * Parse as a string and interpret it as a geographic position.  The input
     * string is broken into space (or comma) separated pieces and Basic
     * decision on which format is based on number of components
     * -# MGRS
     * -# "Lat Long" or "Long Lat"
     * -# "Zone Easting Northing" or "Easting Northing Zone"
     *
     * The following inputs are approximately the same (Ar Ramadi Bridge, Iraq)
     * - Latitude and Longitude
     *   -  33.44      43.27
     *   -  N33d26.4'  E43d16.2'
     *   -  43d16'12"E 33d26'24"N
     * - MGRS
     *   -  38SLC301
     *   -  38SLC391014
     *   -  38SLC3918701405
     *   -  37SHT9708
     * - UTM
     *   -  38N 339188 3701405
     *   -  897039 3708229 37N
     *
     * Latitude and Longitude parsing.  Latitude precedes longitude, unless a
     * N, S, E, W hemisphere designator is used on one or both coordinates.
     * Thus
     * - 40 -75
     * - N40 W75
     * - -75 N40
     * - 75W 40N
     * - E-75 -40S
     * .
     * are all the same position.  The coodinates may be given in decimal
     * degrees, degrees and decimal minutes, degrees, minutes, seconds, etc.
     * Use d, ', and " to make off the degrees, minutes and seconds.  Thus
     * - 40d30'30"
     * - 40d30'30
     * - 40d30.5'
     * - 40d30.5
     * - 40.508333333
     * .
     * all specify the same angle.  The leading sign applies to all components
     * so -1d30 is -(1+30/60) = -1.5.  Latitudes must be in the range [-90, 90]
     * and longitudes in the range [-180, 360].  Internally longitudes are
     * reduced to the range [-180, 180).
     *
     * UTM/UPS parsing.  For UTM zones (-80 <= Lat <= 84), the zone designator
     * is made up of a zone number (for 1 to 60) and a hemisphere letter (N or
     * S), e.g., 38N.  The latitude zone designer ([C&ndash;M] in the southern
     * hemisphere and [N&ndash;X] in the northern) should NOT be used.  (This
     * is part of the MGRS coordinate.)  The zone designator for the poles
     * (where UPS is employed) is a hemisphere letter by itself, i.e., N or S.
     *
     * MGRS parsing interprets the grid references as square area at the
     * specified precision (1m, 10m, 100m, etc.).  If \e centerp = true (the
     * default), the center of this square is then taken to be the precise
     * position; thus:
     * - 38SMB           = 38N 450000 3650000
     * - 38SMB4484       = 38N 444500 3684500
     * - 38SMB44148470   = 38N 444145 3684705
     * .
     * Otherwise, the "south-west" corner of the square is used, i.e.,
     * - 38SMB           = 38N 400000 3600000
     * - 38SMB4484       = 38N 444000 3684000
     * - 38SMB44148470   = 38N 444140 3684700
     **********************************************************************/
    explicit GeoCoords(const std::string& s, bool centerp = true)
    { Reset(s, centerp); }

    /**
     * Specify the location in terms of \e latitude (degrees) and \e longitude
     * (degrees).  Use \e zone to force the UTM/UPS representation to use a
     * specified zone using the rules given in UTMUPS::zone.
     **********************************************************************/
    GeoCoords(real latitude, real longitude, int zone = UTMUPS::STANDARD) {
      Reset(latitude, longitude, zone);
    }

    /**
     * Specify the location in terms of UPS/UPS \e zone (zero means UPS),
     * hemisphere \e northp (false means south, true means north), \e easting
     * (meters) and \e northing (meters).
     **********************************************************************/
    GeoCoords(int zone, bool northp, real easting, real northing) {
      Reset(zone, northp, easting, northing);
    }

    /**
     * Reset the location as a 1-element, 2-element, or 3-element string.  See
     * GeoCoords(const string& s, bool centerp).
     **********************************************************************/
    void Reset(const std::string& s, bool centerp = true);

    /**
     * Reset the location in terms of \e latitude and \e longitude.  See
     * GeoCoords(real latitude, real longitude, int zone).
     **********************************************************************/
    void Reset(real latitude, real longitude, int zone = UTMUPS::STANDARD) {
      UTMUPS::Forward(latitude, longitude,
                      _zone, _northp, _easting, _northing, _gamma, _k,
                      zone);
      _lat = latitude;
      _long = longitude;
      if (_long >= 180)
        _long -= 360;
      CopyToAlt();
    }

    /**
     * Reset the location in terms of UPS/UPS \e zone, hemisphere \e northp, \e
     * easting, and \e northing.  See GeoCoords(int zone, bool northp,
     * real easting, real northing).
     **********************************************************************/
    void Reset(int zone, bool northp, real easting, real northing) {
      UTMUPS::Reverse(zone, northp, easting, northing,
                      _lat, _long, _gamma, _k);
      _zone = zone;
      _northp = northp;
      _easting = easting;
      _northing = northing;
      FixHemisphere();
      CopyToAlt();
    }

    /**
     * Return latitude (degrees)
     **********************************************************************/
    Math::real Latitude() const throw() { return _lat; }

    /**
     * Return longitude (degrees)
     **********************************************************************/
    Math::real Longitude() const throw() { return _long; }

    /**
     * Return easting (meters)
     **********************************************************************/
    Math::real Easting() const throw() { return _easting; }

    /**
     * Return northing (meters)
     **********************************************************************/
    Math::real Northing() const throw() { return _northing; }

    /**
     * Return meridian convergence (degrees) for the UTM/UPS projection.
     **********************************************************************/
    Math::real Convergence() const throw() { return _gamma; }

    /**
     * Return scale for the UTM/UPS projection.
     **********************************************************************/
    Math::real Scale() const throw() { return _k; }

    /**
     * Return hemisphere (false means south, true means north).
     **********************************************************************/
    bool Northp() const throw() { return _northp; }

    /**
     * Return hemisphere letter N or S.
     **********************************************************************/
    char Hemisphere() const throw() { return _northp ? 'N' : 'S'; }

    /**
     * Return the zone corresponding to the input (return 0 for UPS).
     **********************************************************************/
    int Zone() const throw() { return _zone; }

    /**
     * Use zone number, \e zone, for the alternate representation.  See
     * UTMUPS::zone for more information on the interpretation of \e zone.
     * Note that \e zone == UTMUPS::STANDARD (the default) use the standard UPS
     * or UTM zone, UTMUPS::MATCH does nothing retaining the existing alternate
     * representation.  Before this is called the alternate zone is the input
     * zone.
     **********************************************************************/
    void SetAltZone(int zone = UTMUPS::STANDARD) const {
      if (zone == UTMUPS::MATCH)
        return;
      zone = UTMUPS::StandardZone(_lat, _long, zone);
      if (zone == _zone)
        CopyToAlt();
      else {
        bool northp;
        UTMUPS::Forward(_lat, _long,
                        _alt_zone, northp,
                        _alt_easting, _alt_northing, _alt_gamma, _alt_k,
                        zone);
      }
    }

    /**
     * Returns the current alternate zone (return 0 for UPS).
     **********************************************************************/
    int AltZone() const throw() { return _alt_zone; }

    /**
     * Return easting (meters) for alternate zone.
     **********************************************************************/
    Math::real AltEasting() const throw() { return _alt_easting; }

    /**
     * Return northing (meters) for alternate zone.
     **********************************************************************/
    Math::real AltNorthing() const throw() { return _alt_northing; }

    /**
     * Return meridian convergence (degrees) for altermate zone.
     **********************************************************************/
    Math::real AltConvergence() const throw() { return _alt_gamma; }

    /**
     * Return scale for altermate zone.
     **********************************************************************/
    Math::real AltScale() const throw() { return _alt_k; }

    /**
     * Return string with latitude and longitude as signed decimal degrees.
     * Precision \e prec specifies accuracy of representation as follows:
     * - prec = -5 (min), 1d
     * - prec = 0, 10<sup>-5</sup>d (about 1m)
     * - prec = 3, 10<sup>-8</sup>d
     * - prec = 9 (max), 10<sup>-14</sup>d
     **********************************************************************/
    std::string GeoRepresentation(int prec = 0) const;

    /**
     * Return string with latitude and longitude as degrees, minutes, seconds,
     * and hemisphere.  Precision \e prec specifies accuracy of representation
     * as follows:
     * - prec = -5 (min), 1d
     * - prec = -4, 0.1d
     * - prec = -3, 1'
     * - prec = -2, 0.1'
     * - prec = -1, 1"
     * - prec = 0, 0.1" (about 3m)
     * - prec = 1, 0.01"
     * - prec = 10 (max), 10<sup>-11</sup>"
     **********************************************************************/
    std::string DMSRepresentation(int prec = 0) const;

    /**
     * Return MGRS string.  This gives the coordinates of the enclosing grid
     * square with size given by the precision \e prec.  Thus 38N 444180
     * 3684790 converted to a MGRS coordinate at precision -2 (100m) is
     * 38SMB441847 and not 38SMB442848.  Precision \e prec specifies the
     * precision of the MSGRS string as follows:
     * - prec = -5 (min), 100km
     * - prec = -4, 10km
     * - prec = -3, 1km
     * - prec = -2, 100m
     * - prec = -1, 10m
     * - prec = 0, 1m
     * - prec = 1, 0.1m
     * - prec = 6 (max), 1um
     **********************************************************************/
    std::string MGRSRepresentation(int prec = 0) const;

    /**
     * Return string consisting of UTM/UPS zone designator, easting, and
     * northing,  Precision \e prec specifies accuracy of representation
     * as follows:
     * - prec = -5 (min), 100km
     * - prec = -3, 1km
     * - prec = 0, 1m
     * - prec = 3, 1mm
     * - prec = 6, 1um
     * - prec = 9 (max), 1nm
     **********************************************************************/
    std::string UTMUPSRepresentation(int prec = 0) const;

    /**
     * Return MGRS string using alternative zone.  See MGRSRepresentation for
     * the interpretation of \e prec.
     **********************************************************************/
    std::string AltMGRSRepresentation(int prec = 0) const;

    /**
     * Return string consisting of alternate UTM/UPS zone designator, easting,
     * and northing.  See UTMUPSRepresentation for the interpretation of \e
     * prec.
     **********************************************************************/
    std::string AltUTMUPSRepresentation(int prec = 0) const;
  };

} // namespace GeographicLib
#endif  // GEOGRAPHICLIB_GEOCOORDS_HPP

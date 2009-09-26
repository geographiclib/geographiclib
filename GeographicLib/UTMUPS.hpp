/**
 * \file UTMUPS.hpp
 * \brief Header for GeographicLib::UTMUPS class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/


#if !defined(GEOGRAPHICLIB_UTMUPS_HPP)
#define GEOGRAPHICLIB_UTMUPS_HPP "$Id$"

#include "GeographicLib/Constants.hpp"
#include <string>
#include <sstream>

namespace GeographicLib {

  /**
   * \brief Convert between Geographic coordinates and UTM/UPS
   *
   * UTM and UPS are defined
   * - <a href="http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf">
   *   The Universal Grids: Universal Transverse Mercator (UTM) and Universal
   *   Polar Stereographic (UPS)</a>, Defense Mapping Agency, Technical Manual
   *   TM8358.2 (1989).
   * .
   * Section 2-3 defines UTM and section 3-2.4 defines UPS.  This document also
   * includes approximate algorithms for the computation of the underlying
   * transverse Mercator and polar stereographic projections.  Here we
   * substitute much more accurate algorithms given by
   * GeographicLib:TransverseMercator and GeographicLib:PolarStereographic.
   *
   * In this implementation, the conversions are closed, i.e., output from
   * Forward is legal input for Reverse and vice versa.  The error is about 5nm
   * in each direction.  However, the conversion from legal UTM/UPS coordinates
   * to geographic coordinates and back might throw an error if the initial
   * point is within 5nm of the edge of the allowed range for the UTM/UPS
   * coordinates.
   *
   * The simplest way to guarantee the closed property is to define allowed
   * ranges for the eastings and northings for UTM and UPS coordinates.  The
   * UTM boundaries are the same for all zones.  (The only place the
   * exceptional nature of the zone boundaries is evident is when converting to
   * UTM/UPS coordinates requesting the standard zone.)  The MGRS lettering
   * scheme imposes natural limits on UTM/UPS coordinates which may be
   * converted into MGRS coordinates.  For the conversion to/from geographic
   * coordinates these ranges have been extended by 100km in order to provide a
   * generous overlap between UTM and UPS and between UTM zones.
   *
   * The <a href="http://www.nga.mil">NGA</a> software package
   * <a href="http://earth-info.nga.mil/GandG/geotrans/index.html">geotrans</a>
   * also provides conversions to and from UTM and UPS.  Version 2.4.2 (and
   * earlier) suffers from some drawbacks:
   * - Inconsistent rules are used to determine the whether a particular UTM or
   *   UPS coordinate is legal.  A more systematic approach is taken here.
   * - The underlying projections are not very accurately implemented.
   **********************************************************************/
  class UTMUPS {
  private:
    typedef Math::real_t real_t;
    static const real_t falseeasting[4];
    static const real_t falsenorthing[4];
    static const real_t mineasting[4];
    static const real_t maxeasting[4];
    static const real_t minnorthing[4];
    static const real_t maxnorthing[4];
    static real_t CentralMeridian(int zone) throw()
    { return real_t(6 * zone - 183); }
    template<typename T> static std::string str(T x) {
      std::ostringstream s; s << x; return s.str();
    }
    static void CheckLatLon(real_t lat, real_t lon);
    // Throw an error if easting or northing are outside standard ranges.  If
    // throwp = false, return bool instead.
    static bool CheckCoords(bool utmp, bool northp, real_t x, real_t y,
			    bool throwp = true);
    UTMUPS();			// Disable constructor
  public:

    /**
     * Return the standard zone for latitude \e lat (degrees) and longitude \e
     * lon (degrees).  Return 0 if in the standard regions for UPS otherwise
     * return the UTM zone.  This includes the Norway and Svalbard exceptions.
     * The tests on latitudes and longitudes are all closed on the lower end
     * open on the upper.  Thus for UTM zone 38, latitude is in [-80, 84) and
     * longitude is in [42, 48).  This is exact.
     **********************************************************************/
    static int StandardZone(real_t lat, real_t lon) throw();

    /**
     * Convert geographic coordinates to UTM or UPS coordinate.  Given latitude
     * \e lat (degrees), and longitude \e lon (degrees), return \e zone (zero
     * indicates UPS), hemisphere \e northp (false means south, true means
     * north), easting \e x (meters), and northing \e y (meters).  The prefered
     * zone for the result can be specified with \e setzone (negative means
     * result of UTMUPS::StandardZone, zero means UPS, positive means a
     * particular UTM zone), Throw error if the resulting easting or northing
     * is outside the allowed range (see Reverse). This also returns meridian
     * convergence \e gamma (degrees) and scale \e k.  The accuracy of the
     * conversion is about 5nm.
     *
     * To extent the standard UTM zones into the UPS regions use \e setzone =
     * UTMUPS::StandardZone(max(-80.0, min(80.0, \e lat))).
     **********************************************************************/
    static void Forward(real_t lat, real_t lon,
			int& zone, bool& northp, real_t& x, real_t& y,
			real_t& gamma, real_t& k,
			int setzone = -1);

    /**
     * Convert UTM or UPS coordinate to geographic coordinates .  Given zone \e
     * zone (\e zone == 0 indicates UPS), hemisphere \e northp (false means
     * south, true means north), easting \e x (meters), and northing \e y
     * (meters), return latitude \e lat (degrees) and longitude \e lon
     * (degrees).  Throw error if easting or northing is outside the allowed
     * range (see below).  This also returns meridian convergence \e gamma
     * (degrees) and scale \e k.  The accuracy of the conversion is about 5nm.
     *
     * UTM eastings are allowed to be in the range [0km, 1000km], northings are
     * allowed to be in in [0km, 9600km] for the northern hemisphere and in
     * [900km, 10000km] for the southern hemisphere.  (However UTM northings
     * can be continued across the equator.  So the actual limits on the
     * northings are [-9100km, 9600km] for the "northern" hemisphere and
     * [900km, 19600km] for the "southern" hemisphere.)
     *
     * UPS eastings and northings are allowed to be in the range [1200km,
     * 2800km] in the northern hemisphere and in [700km, 3100km] in the
     * southern hemisphere.
     *
     * These ranges are 100km larger than allowed for the conversions to MGRS.
     * (100km is the maximum extra padding consistent with eastings remaining
     * non-negative.)  This allows generous overlaps between zones and UTM and
     * UPS.  No checks are performed beyond these (e.g., to limit the distance
     * outside the standard zone boundaries).
     **********************************************************************/
    static void Reverse(int zone, bool northp, real_t x, real_t y,
			real_t& lat, real_t& lon, real_t& gamma, real_t& k);

    /**
     * Forward without returning convergence and scale.
     **********************************************************************/
    static void Forward(real_t lat, real_t lon,
			int& zone, bool& northp, real_t& x, real_t& y,
			int setzone = -1) {
      real_t gamma, k;
      Forward(lat, lon, zone, northp, x, y, gamma, k, setzone);
    }

    /**
     * Reverse without returning convergence and scale.
     **********************************************************************/
    static void Reverse(int zone, bool northp, real_t x, real_t y,
			real_t& lat, real_t& lon) {
      real_t gamma, k;
      Reverse(zone, northp, x, y, lat, lon, gamma, k);
    }

    /**
     * Decode a UTM/UPS zone string, \e zonestr, returning the resulting \e
     * zone and hemisphere thru \e northp (true for northern and false for
     * southern hemispheres).  For UTM, \e zonestr has the form of a zone
     * number in the range [1,60] followed by a hemisphere letter, N or S.  For
     * UPS, it consists just of the hemisphere letter.  The returned value of
     * \e zone is 0 for UPS.  Note well that "38S" indicates the southern
     * hemisphere of zone 38 and not latitude band S, [32,40].  N, 01S, 2N, 38S
     * are legal.  0N, 001S, 61N, 38P are illegal.
     **********************************************************************/
    static void DecodeZone(const std::string& zonestr,
			   int& zone, bool& northp);

    /**
     * Encode a UTM/UPS zone string given the \e zone and hemisphere \e northp.
     * \e zone must be in the range [0,60] with \e zone = 0 indicating UPS (but
     * the resulting string does not contain "0").  This reverses DecodeZone.
     **********************************************************************/
    static std::string EncodeZone(int zone, bool northp);

    /**
     * The shift necessary to align N and S halves of a UTM zone
     * (10<sup>7</sup>).
     **********************************************************************/
    static real_t UTMShift() throw();

  };

} // namespace GeographicLib
#endif

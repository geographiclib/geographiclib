/**
 * \file UTMUPS.hpp
 * \brief Header for GeographicLib::UTMUPS class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/


#if !defined(GEOGRAPHICLIB_UTMUPS_HPP)
#define GEOGRAPHICLIB_UTMUPS_HPP "$Id: UTMUPS.hpp 6867 2010-09-11 13:04:26Z karney $"

#include "GeographicLib/Constants.hpp"
#include <sstream>

namespace GeographicLib {

  /**
   * \brief Convert between Geographic coordinates and UTM/UPS
   *
   * UTM and UPS are defined
   * - J. W. Hager, J. F. Behensky, and B. W. Drew,
   *   <a href="http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf">
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
    typedef Math::real real;
    static const real falseeasting[4];
    static const real falsenorthing[4];
    static const real mineasting[4];
    static const real maxeasting[4];
    static const real minnorthing[4];
    static const real maxnorthing[4];
    static real CentralMeridian(int zone) throw()
    { return real(6 * zone - 183); }
    template<typename T> static std::string str(T x) {
      std::ostringstream s; s << x; return s.str();
    }
    static void CheckLatLon(real lat, real lon);
    // Throw an error if easting or northing are outside standard ranges.  If
    // throwp = false, return bool instead.
    static bool CheckCoords(bool utmp, bool northp, real x, real y,
                            bool msgrlimits = false, bool throwp = true);
    UTMUPS();                   // Disable constructor

  public:

    /**
     * In this class we bring together the UTM and UPS coordinates systems.
     * The UTM divides the earth between latitudes -80 and 84 into 60 zones
     * numbered 1 thru 60.  Zone assign zone number 0 to the UPS regions,
     * covering the two poles.  Within UTMUPS, non-negative zone numbers refer
     * to one of the "physical" zones, 0 for UPS and [1, 60] for UTM.  Negative
     * "pseudo-zone" numbers are used to select one of the physical zones.
     **********************************************************************/
    enum zonespec {
      /**
       * The smallest pseudo-zone number.
       **********************************************************************/
      MINPSEUDOZONE = -3,
      /**
       * If a coordinate already include zone information (e.g., it is an MGRS
       * coordinate), use that, otherwise apply the UTMUPS::STANDARD rules.
       **********************************************************************/
      MATCH = -3,
      /**
       * Apply the standard rules for UTM zone assigment extending the UTM zone
       * to each pole to give a zone number in [1, 60].  For example, use UTM
       * zone 38 for longitude in [42, 48).  The rules include the Norway and
       * Svalbard exceptions.
       **********************************************************************/
      UTM = -2,
      /**
       * Apply the standard rules for zone assignment to give a zone number in
       * [0, 60].  If the latitude is not in [-80, 84), then use UTMUPS::UPS =
       * 0, otherwise apply the rules for UTMUPS::UTM.  The tests on latitudes
       * and longitudes are all closed on the lower end open on the upper.
       * Thus for UTM zone 38, latitude is in [-80, 84) and longitude is in
       * [42, 48).
       **********************************************************************/
      STANDARD = -1,
      /**
       * The largest pseudo-zone number.
       **********************************************************************/
      MAXPSEUDOZONE = -1,
      /**
       * The smallest physical zone number.
       **********************************************************************/
      MINZONE = 0,
      /**
       * The zone number used for UPS
       **********************************************************************/
      UPS = 0,
      /**
       * The smallest UTM zone number.
       **********************************************************************/
      MINUTMZONE = 1,
      /**
       * The largest UTM zone number.
       **********************************************************************/
      MAXUTMZONE = 60,
      /**
       * The largest physical zone number.
       **********************************************************************/
      MAXZONE = 60,
    };

    /**
     * The standard zone.
     *
     * @param[in] lat latitude (degrees).
     * @param[in] lon longitude (degrees).
     * @param[in] setzone zone override (optional)
     *
     * This is exact.  If the optional argument \e setzone is given then use
     * that zone if it is non-negative, otherwise apply the rules given in
     * UTMUPS::zonespec.  Throws an error if \e setzone is outsize the range
     * [UTMUPS::MINPSEUDOZONE, UTMUPS::MAXZONE] = [-3, 60].
     **********************************************************************/
    static int StandardZone(real lat, real lon, int setzone = STANDARD);

    /**
     * Forward projection, from geographic to UTM/UPS.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] zone the UTM zone (zero means UPS).
     * @param[out] northp hemisphere of location (true means northern, false
     *   means southern).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     * @param[in] setzone zone override.
     * @param[in] mgrslimits if true enforce the stricted MGRS limits on the
     *   coordinates (default = false).
     *
     * The prefered zone for the result can be specified with \e setzone, see
     * UTMUPS::StandardZone.  Throw error if the resulting easting or northing
     * is outside the allowed range (see Reverse), in which case the arguments
     * are unchanged.  This also returns meridian convergence \e gamma
     * (degrees) and scale \e k.  The accuracy of the conversion is about 5nm.
     **********************************************************************/
    static void Forward(real lat, real lon,
                        int& zone, bool& northp, real& x, real& y,
                        real& gamma, real& k,
                        int setzone = STANDARD, bool mgrslimits = false);

    /**
     * Reverse projection, from  UTM/UPS to geographic.
     *
     * @param[in] zone the UTM zone (zero means UPS).
     * @param[in] northp hemisphere of location (true means northern, false
     *   means southern).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     * @param[in] mgrslimits if true enforce the stricted MGRS limits on the
     *   coordinates (default = false).
     *
     * Throw error if easting or northing is outside the allowed range (see
     * below), in which case the arguments are unchanged.  The accuracy of the
     * conversion is about 5nm.
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
     * UPS.  If \e mgrslimits = true, then all the ranges are shrunk by 100km
     * so that they agree with the stricter MGRS ranges.  No checks are
     * performed besides these (e.g., to limit the distance outside the
     * standard zone boundaries).
     **********************************************************************/
    static void Reverse(int zone, bool northp, real x, real y,
                        real& lat, real& lon, real& gamma, real& k,
                        bool mgrslimits = false);

    /**
     * UTMUPS::Forward without returning convergence and scale.
     **********************************************************************/
    static void Forward(real lat, real lon,
                        int& zone, bool& northp, real& x, real& y,
                        int setzone = STANDARD, bool mgrslimits = false) {
      real gamma, k;
      Forward(lat, lon, zone, northp, x, y, gamma, k, setzone, mgrslimits);
    }

    /**
     * UTMUPS::Reverse without returning convergence and scale.
     **********************************************************************/
    static void Reverse(int zone, bool northp, real x, real y,
                        real& lat, real& lon, bool mgrslimits = false) {
      real gamma, k;
      Reverse(zone, northp, x, y, lat, lon, gamma, k, mgrslimits);
    }

    /**
     * Decode a UTM/UPS zone string.
     *
     * @param[in] zonestr string represention of zone and hemisphere.
     * @param[out] zone the UTM zone (zero means UPS).
     * @param[out] northp the hemisphere (true means northern, false
     *   means southern).
     *
     * For UTM, \e zonestr has the form of a zone number in the range
     * [UTMUPS::MINUTMZONE, UTMUPS::MAXUTMZONE] = [1, 60] followed by a
     * hemisphere letter, N or S.  For UPS, it consists just of the hemisphere
     * letter.  The returned value of \e zone is UTMUPS::UPS = 0 for UPS.  Note
     * well that "38S" indicates the southern hemisphere of zone 38 and not
     * latitude band S, [32, 40].  N, 01S, 2N, 38S are legal.  0N, 001S, 61N,
     * 38P are illegal.  Throws an error is the zone string is malformed.
     **********************************************************************/
    static void DecodeZone(const std::string& zonestr, int& zone, bool& northp);

    /**
     * Encode a UTM/UPS zone string.
     *
     * @param[out] zone the UTM zone (zero means UPS).
     * @param[out] northp the hemisphere (true means northern, false
     *   means southern).
     * @return string represention of zone and hemisphere.
     *
     * \e zone must be in the range [UTMUPS::MINZONE, UTMUPS::MAXZONE] = [0,
     * 60] with \e zone = UTMUPS::UPS, 0, indicating UPS (but the resulting
     * string does not contain "0").  This reverses UTMUPS::DecodeZone.
     **********************************************************************/
    static std::string EncodeZone(int zone, bool northp);

    /**
     * @return shift (meters) necessary to align N and S halves of a UTM zone
     * (10<sup>7</sup>).
     **********************************************************************/
    static Math::real UTMShift() throw();

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the WGS84 ellipsoid (meters).
     *
     * (The WGS84 values is returned because the UTM and UPS projections are
     * based on this ellipsoid.)
     **********************************************************************/
    static Math::real MajorRadius() throw() { return Constants::WGS84_a(); }

    /**
     * @return \e r the inverse flattening of the WGS84 ellipsoid.
     *
     * (The WGS84 values is returned because the UTM and UPS projections are
     * based on this ellipsoid.)
     **********************************************************************/
    static Math::real InverseFlattening() throw()
    { return Constants::WGS84_r(); }
    ///@}
  };

} // namespace GeographicLib
#endif

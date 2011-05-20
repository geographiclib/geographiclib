/**
 * \file DMS.hpp
 * \brief Header for GeographicLib::DMS class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_DMS_HPP)
#define GEOGRAPHICLIB_DMS_HPP "$Id: DMS.hpp 6866 2010-09-11 02:15:29Z karney $"

#include "GeographicLib/Constants.hpp"
#include <sstream>
#include <iomanip>

namespace GeographicLib {

  /**
   * \brief Convert between degrees and %DMS representation.
   *
   * Parse a string representing degree, minutes, and seconds and return the
   * angle in degrees and format an angle in degrees as degree, minutes, and
   * seconds.
   **********************************************************************/
  class DMS {
  private:
    typedef Math::real real;
    static int lookup(const std::string& s, char c) throw() {
      std::string::size_type r = s.find(toupper(c));
      return r == std::string::npos ? -1 : int(r);
    }
    template<typename T> static std::string str(T x) {
      std::ostringstream s; s << x; return s.str();
    }
    static const std::string hemispheres;
    static const std::string signs;
    static const std::string digits;
    static const std::string dmsindicators;
    static const std::string components[3];
    DMS();                      // Disable constructor

  public:

    /**
     * Indicator for presence of hemisphere indicator (N/S/E/W) on latitudes
     * and longitudes.
     **********************************************************************/
    enum flag { 
      /**
       * No indicator present.
       * @hideinitializer
       **********************************************************************/
      NONE = 0,
      /**
       * Latitude indicator (N/S) present.
       * @hideinitializer
       **********************************************************************/
      LATITUDE = 1,
      /**
       * Longitude indicator (E/W) present.
       * @hideinitializer
       **********************************************************************/
      LONGITUDE = 2,
      /**
       * Used in Encode to indicate output of an azimuth in [000, 360) with no
       * letter indicator.
       * @hideinitializer
       **********************************************************************/
      AZIMUTH = 3,
      /**
       * Used in Encode to indicate output of a plain number.
       * @hideinitializer
       **********************************************************************/
      NUMBER = 4,
    };

    /**
     * Indicator for trailing units on an angle.
     **********************************************************************/
    enum component {
      /**
       * Trailing unit is degrees.
       * @hideinitializer
       **********************************************************************/
      DEGREE = 0,
      /**
       * Trailing unit is arc minutes.
       * @hideinitializer
       **********************************************************************/
      MINUTE = 1,
      /**
       * Trailing unit is arc seconds.
       * @hideinitializer
       **********************************************************************/
      SECOND = 2,
    };

    /**
     * Convert a string in DMS to an angle.
     *
     * @param[in] dms string input.
     * @param[out] ind a DMS::flag value indicating the presence of a
     *   hemisphere indicator.
     * @return angle (degrees).
     *
     * Degrees, minutes, and seconds are indicated by the letters d, ', &quot;,
     * and these components may only be given in this order.  Any (but not all)
     * components may be omitted.  The last component indicator may be omitted
     * and is assumed to be tbe next smallest unit (thus 33d10 is interpreted
     * as 33d10').  The final component may be a decimal fraction but the
     * non-final components must be integers.  The integer parts of the minutes
     * and seconds components must be less than 60.  A single leading sign is
     * permitted.  A hemisphere designator (N, E, W, S) may be added to tbe
     * beginning or end of the string.  The result is multiplied by the implied
     * signed of the hemisphere designator (negative for S and W).  In addition
     * \e ind is set to DMS::LATITUDE if N or S is present, to DMS::LONGITUDE
     * if E or W is present, and to DMS::NONE otherwise.  Throws an error on a
     * malformed string.  No check is performed on the range of the result.
     **********************************************************************/
    static Math::real Decode(const std::string& dms, flag& ind);

    /**
     * Convert DMS to an angle.
     *
     * @param[in] d degrees.
     * @param[in] m arc minutes.
     * @param[in] s arc seconds.
     * @return angle (degress)
     *
     * This does not propagate the sign on \e d to the other components, so
     * -3d20' would need to be represented as - DMS::Decode(3.0, 20.0) or
     * DMS::Decode(-3.0, -20.0).
     **********************************************************************/
    static Math::real Decode(real d, real m = 0, real s = 0) throw()
    { return d + (m + s/real(60))/real(60); }

    /**
     * Convert a string to a real number.
     *
     * @param[in] str string input.
     * @return decoded number.
     **********************************************************************/
    static Math::real Decode(const std::string& str);

    /**
     * Convert a pait of strings to latitude and longitude.
     *
     * @param[in] dmsa first string.
     * @param[in] dmsb second string.
     * @param[out] lat latitude.
     * @param[out] lon longitude.
     *
     * By default, the \e lat (resp., \e lon) is assigned to the results of
     * decoding \e dmsa (resp., \e dmsb).  However this is overridden if either
     * \e dmsa or \e dmsb contain a latitude or longitude hemisphere designator
     * (N, S, E, W).  Throws an error if the decoded numbers are out of the
     * ranges [-90<sup>o</sup>, 90<sup>o</sup>] for latitude and
     * [-180<sup>o</sup>, 360<sup>o</sup>] for longitude and, in which case \e
     * lat and \e lon are unchanged.  Finally the longitude is reduced to the
     * range [-180<sup>o</sup>, 180<sup>o</sup>).
     **********************************************************************/
    static void DecodeLatLon(const std::string& dmsa, const std::string& dmsb,
                             real& lat, real& lon);

    /**
     * Convert a string to an angle in degrees.
     *
     * @param[in] angstr input string.
     * @return angle (degrees)
     *
     * No hemisphere designator is allowed and no check is done on the range of
     * the result.
     **********************************************************************/
    static Math::real DecodeAngle(const std::string& angstr);

    /**
     * Convert a string to an azimuth in degrees.
     *
     * @param[in] azistr input string.
     * @return azimuth (degrees)
     *
     * A hemisphere designator E/W can be used; the result is multiplied by -1
     * if W is present.  Throws an error if the result is out of the range
     * [-180<sup>o</sup>, 360<sup>o</sup>].  Finally the azimuth is reduced to
     * the range [-180<sup>o</sup>, 180<sup>o</sup>).
     **********************************************************************/
    static Math::real DecodeAzimuth(const std::string& azistr);

    /**
     * Convert angle (in degrees) into a DMS string.
     *
     * @param[in] angle input angle (degrees)
     * @param[in] trailing DMS::component value indicating the trailing units
     *   on the string and this is given as a decimal number if necessary.
     * @param[in] prec the number of digits after the decimal point for the
     *   trailing component.
     * @param[in] ind DMS::flag value indicated additional formatting.
     * @return formatting string
     *
     * The interpretation of \e ind is as follows:
     * - ind == DMS::NONE, signed result no leading zeros on degrees except in
     *   the units place, e.g., -8d03'.
     * - ind == DMS::LATITUDE, trailing N or S hemisphere designator, no sign,
     *   pad degress to 2 digits, e.g., 08d03'S.
     * - ind == DMS::LONGITUDE, trailing E or W hemisphere designator, no
     *   sign, pad degress to 3 digits, e.g., 008d03'W.
     * - ind == DMS::AZIMUTH, convert to the range [0, 360<sup>o</sup>), no
     *   sign, pad degrees to 3 digits, , e.g., 351d57'.
     * .
     * The integer parts of the minutes and seconds components are always given
     * with 2 digits.
     **********************************************************************/
    static std::string Encode(real angle, component trailing, unsigned prec,
                              flag ind = NONE);

    /**
     * Convert angle into a DMS string selecting the trailing component
     * based on the precision.
     *
     * @param[in] angle input angle (degrees)
     * @param[in] prec the precision relative to 1 degree.
     * @param[in] ind DMS::flag value indicated additional formatting.
     * @return formatting string
     *
     * \e prec indicates the precision relative to 1 degree, e.g., \e prec = 3
     * gives a result accurate to 0.1' and \e prec = 4 gives a result accurate
     * to 1&quot;.  \e ind is interpreted as in DMS::Encode with the additional
     * facility at DMS::NUMBER treats \e angle a number in fixed format with
     * precision \e prec.
     **********************************************************************/
    static std::string Encode(real angle, unsigned prec, flag ind = NONE) {
      if (ind == NUMBER) {
        std::ostringstream s;
        s << std::fixed << std::setprecision(prec) << angle;
        return s.str();
      } else
        return Encode(angle,
                      prec < 2 ? DEGREE : (prec < 4 ? MINUTE : SECOND),
                      prec < 2 ? prec : (prec < 4 ? prec - 2 : prec - 4),
                      ind);
    }

    /**
     * Split angle into degrees and mintues
     *
     * @param[in] ang angle (degrees)
     * @param[out] d degrees (an integer returned as a real)
     * @param[out] m arc minutes.
     **********************************************************************/
    static void Encode(real ang, real& d, real& m) throw() {
      d = int(ang); m = 60 * (ang - d);
    }

    /**
     * Split angle into degrees and mintues and seconds.
     *
     * @param[in] ang angle (degrees)
     * @param[out] d degrees (an integer returned as a real)
     * @param[out] m arc minutes (an integer returned as a real)
     * @param[out] s arc seconds.
     **********************************************************************/
    static void Encode(real ang, real& d, real& m, real& s) throw() {
      d = int(ang); ang = 60 * (ang - d);
      m = int(ang); s = 60 * (ang - m);
    }

  };

} // namespace GeographicLib

#endif

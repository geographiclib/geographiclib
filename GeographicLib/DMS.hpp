/**
 * \file DMS.hpp
 * \brief Header for GeographicLib::DMS class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_DMS_HPP)
#define GEOGRAPHICLIB_DMS_HPP "$Id$"

#include "GeographicLib/Constants.hpp"
#include <string>
#include <sstream>

namespace GeographicLib {

  /**
   * \brief Convert between degrees to %DMS representation.
   *
   * Parse a string representing degree, minutes, and seconds and return the
   * angle in degrees and format an angle in degrees as degree, minutes, and
   * seconds.
   **********************************************************************/
  class DMS {
  private:
    typedef Math::real_t real_t;
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
    DMS();			// Disable constructor

  public:

    /**
     * Indicator for presence of hemisphere indicator (N/S/E/W) on latitudes
     * and longitudes.  AZIMUTH is used in Encode to indicate output in [000,
     * 360) with no letter indicator.
     **********************************************************************/
    enum flag { NONE = 0, LATITUDE = 1, LONGITUDE = 2, AZIMUTH = 3 };

    /**
     * Indicator for trailing units on an angle.
     **********************************************************************/
    enum component { DEGREE = 0, MINUTE = 1, SECOND = 2 };

    /**
     * Read a string \e dms in DMS format and return the resulting angle in
     * degrees.  Degrees, minutes, and seconds are indicated by the letters d,
     * ', ", and these components may only be given in this order.  Any (but
     * not all) components may be omitted.  The last component indicator may be
     * omitted and is assumed to be tbe next smallest unit (thus 33d10 is
     * interpreted as 33d10').  The final component may be a decimal fraction
     * but the non-final components must be integers.  The integer parts of the
     * minutes and seconds components must be less than 60.  A single leading
     * sign is permitted.  A hemisphere designator (N, E, W, S) may be added to
     * tbe beginning or end of the string.  The result is mulitplied by the
     * implied signed of the hemisphere designator (negative for S and W).  In
     * addition \e flag is used to indicate whether such a designator was found
     * and whether it implies that the angle is a latitude (N or S) or
     * longitude (E or W).
     **********************************************************************/
    static Math::real_t Decode(const std::string& dms, flag& ind);

    /**
     * Convert two strings \e dmsa and \e dmsb to a latitude, \e lat, and
     * longitude, \e lon.  By default, the \e lat (resp., \e lon) is assigned
     * to the results of decoding \e dmsa (resp., \e dmsb).  However this is
     * overridden if either \e dmsa or \e dmsb contain a latitude or longitude
     * hemisphere designator (N, S, E, W).
     **********************************************************************/
    static void DecodeLatLon(const std::string& dmsa, const std::string& dmsb,
			     Math::real_t& lat, Math::real_t& lon);

    /**
     * Convert \e degree into a DMS string.  \e trailing indicates the least
     * significant component of the string (and this component is given as a
     * decimal number if necessary).  \e prec indicates the number of digits
     * after the decimal point for the trailing component.  \e flag indicates
     * whether the result should be include a sign (if negative) or a trailing
     * latitude or longitude hemisphere indicator.  In the latter two cases,
     * the integer part of the degrees component is given with 2 (latitude) or
     * 3 (longitude) digits (with leading zeros if necessary).  The integer
     * parts of the minutes and seconds components are always given with 2
     * digits.
     **********************************************************************/
    static std::string Encode(Math::real_t degree,
			      component trailing,
			      unsigned prec,
			      flag ind = NONE);

    /**
     * Convert \e degree into a DMS string selecting the trailing component
     * based on \e prec.  \e prec indicates the precision relative to 1 degree,
     * e.g., \e prec = 3 gives a result accurate to 0.1' and \e prec = 4 gives
     * a result accurate to 1".
     **********************************************************************/
    static std::string Encode(Math::real_t degree,
			      unsigned prec,
			      flag ind = NONE) {
      return Encode(degree,
		    prec < 2 ? DEGREE : (prec < 4 ? MINUTE : SECOND),
		    prec < 2 ? prec : (prec < 4 ? prec - 2 : prec - 4),
		    ind);
    }
  };

} // namespace GeographicLib

#endif

/**
 * \file Georef.hpp
 * \brief Header for GeographicLib::Georef class
 *
 * Copyright (c) Charles Karney (2015) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEOREF_HPP)
#define GEOGRAPHICLIB_GEOREF_HPP 1

#include <GeographicLib/Constants.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs string
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {

  /**
   * \brief Conversions for the World Geographic Reference System (georef)
   *
   * The World Geographic Reference System is described in
   * - https://en.wikipedia.org/wiki/Georef
   * - http://earth-info.nga.mil/GandG/coordsys/grids/georef.pdf
   * .
   * It provide a compact string representation of a geographic area (expressed
   * as latitude and longitude).
   *
   * Example of use:
   * \include example-Georef.cpp
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT Georef {
  private:
    typedef Math::real real;
    static const std::string digits_;
    static const std::string lontile_;
    static const std::string lattile_;
    static const std::string degrees_;
    enum {
      tile_ = 15,               // The size of tile in degrees
      maxlat_ = 89,             // Max latitude
      lonorig_ = -180 / tile_,  // Origin for longitude tiles
      latorig_ = -90 / tile_,   // Origin for latitude tiles
      base_ = 10,               // Base for minutes
      baselen_ = 4,
      maxprec_ = 5 + 6,         // approximately equivalent to MGRS
      mult2_ = 60,              // number of minuites in a degree
      maxlen_ = baselen_ + 2 * maxprec_,
    };
    Georef();                     // Disable constructor

  public:

    /**
     * Convert from geographic coordinates to georef.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] prec the precision of the resulting georef.
     * @param[out] georef the georef string.
     * @exception GeographicErr if \e la is not in [&minus;90&deg;,
     *   90&deg;].
     * @exception std::bad_alloc if memory for \e georef can't be allocated.
     *
     * \e prec specifies the precision of \e georef as follows:
     * - \e prec = &minus;1 (min), 15&deg;
     * - \e prec = 0, 1&deg;
     * - \e prec = 1, converted to \e prec = 2
     * - \e prec = 2, 1'
     * - \e prec = 3, 0.1'
     * - \e prec = 4, 0.01'
     * - \e prec = 5, 0.001'
     * - &hellip;
     * - \e prec = 11 (max), 10<sup>&minus;9</sup>'
     *
     * If \e lat or \e lon is NaN, then \e georef is set to "INVALID".
     **********************************************************************/
    static void Forward(real lat, real lon, int prec, std::string& georef);

    /**
     * Convert from Georef to geographic coordinates.
     *
     * @param[in] georef the Georef.
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] prec the precision of \e georef.
     * @param[in] centerp if true (the default) return the center
     *   \e georef, otherwise return the south-west corner.
     * @exception GeographicErr if \e georef is illegal.
     *
     * The case of the letters in \e georef is ignored.  \e prec specifies the
     * precision of \e georef as follows:
     * - \e prec = &minus;1 (min), 15&deg;
     * - \e prec = 0, 1&deg;
     * - \e prec = 1, not returned
     * - \e prec = 2, 1'
     * - \e prec = 3, 0.1'
     * - \e prec = 4, 0.01'
     * - \e prec = 5, 0.001'
     * - &hellip;
     * - \e prec = 11 (max), 10<sup>&minus;9</sup>'
     *
     * If the first 3 characters of \e georef are "INV", then \e lat and \e lon
     * are set to NaN and \e prec is unchanged.
     **********************************************************************/
    static void Reverse(const std::string& georef, real& lat, real& lon,
                        int& prec, bool centerp = true);

  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_GEOREF_HPP

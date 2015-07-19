/**
 * \file GARS.hpp
 * \brief Header for GeographicLib::GARS class
 *
 * Copyright (c) Charles Karney (2015) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GARS_HPP)
#define GEOGRAPHICLIB_GARS_HPP 1

#include <GeographicLib/Constants.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs string
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {

  /**
   * \brief Conversions for the Global Area Reference System (GARS)
   *
   * The Global Area Reference System is described in
   * - https://en.wikipedia.org/wiki/Global_Area_Reference_System
   * - http://earth-info.nga.mil/GandG/coordsys/grids/gars.html
   * .
   * It provide a compact string representation of a geographic area (expressed
   * as latitude and longitude).
   *
   * Example of use:
   * \include example-GARS.cpp
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT GARS {
  private:
    typedef Math::real real;
    static const std::string digits_;
    static const std::string letters_;
    enum {
      tile_ = 30,                       // The tile in minutes
      maxlat_ = 90 * 60 / tile_ - 1,    // Max latitude tile
      lonorig_ = -180 * 60 / tile_ - 1, // Origin for longitude tiles
      latorig_ = -90 * 60 / tile_,      // Origin for latitude tiles
      baselon_ = 10,                    // Base for longitude tiles
      baselat_ = 24,                    // Base for latitude tiles
      lonlen_ = 3,
      latlen_ = 2,
      baselen_ = lonlen_ + latlen_,
      maxprec_ = 2,
      mult2_ = 2,               // 6th char gives 2x more precision
      mult3_ = 3,               // 7th char gives 3x more precision
      maxlen_ = baselen_ + maxprec_,
    };
    GARS();                     // Disable constructor

  public:

    /**
     * Convert from geographic coordinates to GARS.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] prec the precision of the resulting GARS.
     * @param[out] gars the GARS string.
     * @exception GeographicErr if \e la is not in [&minus;90&deg;,
     *   90&deg;].
     * @exception std::bad_alloc if memory for \e gars can't be allocated.
     *
     * Internally, \e prec is first put in the range [0, 2].  The meaning of \e
     * prec is as follows:
     * - 0, 30' precision, e.g., 006AG;
     * - 1, 15' precision, e.g., 006AG3;
     * - 2, 5' precision, e.g., 006AG39.
     *
     * If \e lat or \e lon is NaN, then \e gars is set to "INVALID".
     **********************************************************************/
    static void Forward(real lat, real lon, int prec, std::string& gars);

    /**
     * Convert from GARS to geographic coordinates.
     *
     * @param[in] gars the GARS.
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] prec the precision of \e gars.
     * @param[in] centerp if true (the default) return the center of the
     *   \e gars, otherwise return the south-west corner.
     * @exception GeographicErr if \e gars is illegal.
     *
     * The case of the letters in \e gars is ignored.  The value of \e prec is
     * in [0, 2].  The meaning of \e prec is as follows:
     * - 0, 30' precision, e.g., 006AG;
     * - 1, 15' precision, e.g., 006AG3;
     * - 2, 5' precision, e.g., 006AG39.
     *
     * If the first 3 characters of \e gars are "INV", then \e lat and \e lon
     * are set to NaN and \e prec is unchanged.
     **********************************************************************/
    static void Reverse(const std::string& gars, real& lat, real& lon,
                        int& prec, bool centerp = true);

  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_GARS_HPP

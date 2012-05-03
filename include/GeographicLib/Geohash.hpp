/**
 * \file Geohash.hpp
 * \brief Header for GeographicLib::Geohash class
 *
 * Copyright (c) Charles Karney (2008-2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEOHASH_HPP)
#define GEOGRAPHICLIB_GEOHASH_HPP \
  "$Id$"

#include <string>
#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief %Geohash coordinates
   *
   * Example of use:
   * \include example-Geohash.cpp
   **********************************************************************/

  class GEOGRAPHIC_EXPORT Geohash {
  private:
    typedef Math::real real;
    static const int maxprec_ = 18;
    static const unsigned long long mask_ = 1ULL << 45;
    static const int decprec_[];
    static const real loneps_;
    static const real lateps_;
    static const real shift_;
    static const std::string lcdigits_;
    static const std::string ucdigits_;
    Geohash();                     // Disable constructor

  public:

    /**
     * Convert from geodetic coordinates to a geohash.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] prec the length of the resulting geohash.
     * @param[out] geohash the geohash.
     *
     * \e lat should be in the range [-90, 90]; \e lon and \e lon0 should be in
     * the range [-180, 360].  Internally, \e prec is first put in the range
     * [0, 18].
     **********************************************************************/
    static void Forward(real lat, real lon, int prec, std::string& geohash);

    /**
     * Convert from a geohash to geodetic coordinates.
     *
     * @param[in] geohash the geohash.
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] prec the length of the geohash.
     * @param[in] centerp if true (the default) return the center of the
     *   geohash, otherwise return the south-west corner.
     *
     * Only the first 18 characters for \e geohash are considered.
     **********************************************************************/
    static void Reverse(const std::string& geohash, real& lat, real& lon,
                        int& prec, bool centerp = true);

    /**
     * The latitude resolution of a geohash.
     *
     * @param[in] prec the length of the geohash.
     * @return the latitude resolution (degrees).
     *
     * Internally, \e prec is first put in the range [0, 18].
     **********************************************************************/
    static real LatitudeResolution(int prec) throw() {
      prec = std::max(0, std::min(int(maxprec_), prec));
      return 180 * std::pow(0.5, 5 * prec / 2);
    }

    /**
     * The longitude resolution of a geohash.
     *
     * @param[in] prec the length of the geohash.
     * @return the longitude resolution (degrees).
     *
     * Internally, \e prec is first put in the range [0, 18].
     **********************************************************************/
    static real LongitudeResolution(int prec) throw() {
      prec = std::max(0, std::min(int(maxprec_), prec));
      return 360 * std::pow(0.5, 5 * prec - 5 * prec / 2);
    }

    /**
     * The geohash precision required to meet a given geodetic resolution.
     *
     * @param[in] res the minimum of resolution in latitude and longitude
     *   (degrees).
     * @return geohash precision.
     *
     * The return precision is in the range [0, 18].
     **********************************************************************/
    static int GeohashPrecision(real res) throw() {
      res = std::abs(res);
      for (int prec = 0; prec < maxprec_; ++prec)
        if (LongitudeResolution(prec) <= res)
          return prec;
      return maxprec_;
    }

    /**
     * The geohash precision required to meet a given geodetic resolution.
     *
     * @param[in] latres the resolution in latitude (degrees).
     * @param[in] lonres the resolution in longitude (degrees).
     * @return geohash precision.
     *
     * The return precision is in the range [0, 18].
     **********************************************************************/
    static int GeohashPrecision(real latres, real lonres) throw() {
      latres = std::abs(latres);
      lonres = std::abs(lonres);
      for (int prec = 0; prec < maxprec_; ++prec)
        if (LatitudeResolution(prec) <= latres &&
            LongitudeResolution(prec) <= lonres)
          return prec;
      return maxprec_;
    }

    /**
     * The decimal geodetic precision required to match a given geohash
     * precision.  This is the number of digits needed after decimal point in a
     * decimal degrees representation.
     *
     * @param[in] prec the length of the geohash.
     * @return the decimal precision (may be negative).
     *
     * Internally, \e prec is first put in the range [0, 18].  The returned
     * decimal precision is in the range [-2, 12].
     **********************************************************************/
    static int DecimalPrecision(int prec) throw() {
      return -int(std::floor(std::log(LatitudeResolution(prec))/
                             std::log(Math::real(10))));
    }

  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_GEOHASH_HPP

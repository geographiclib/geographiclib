/**
 * \file MagneticModel.hpp
 * \brief Header for GeographicLib::MagneticModel class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_MAGNETICMODEL_HPP)
#define GEOGRAPHICLIB_MAGNETICMODEL_HPP "$Id$"

#include <string>
#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geocentric.hpp>

namespace GeographicLib {

  /**
   * \brief Magnetic model
   *
   * Evaluate the earth's magnetic field according to a model.
   * See
   * - http://geomag.org/models/index.html
   * - http://ngdc.noaa.gov/geomag/EMM/emm.shtml
   * - http://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
   **********************************************************************/

  class GEOGRAPHIC_EXPORT MagneticModel {
  private:
    typedef Math::real real;
    std::string _datafile, _name, _date;
    real _t0, _tmin, _tmax, _a, _minh, _maxh;
    int _N;
    Geocentric _earth;
    std::vector<real> _G, _H, _Gt, _Ht;
    void Field(real lat, real lon, real h, real t, bool diffp,
               real& Bx, real& By, real& Bz,
               real& Bxt, real& Byt, real& Bzt) const;
  public:
    /**
     * Construct a magnetic model.
     *
     * @param[in] datafile the name of the datafile with the parameters of the
     *   model.
     * @param[in] earth Geocentric object for converting coordinates; default
     *   Geocentric::WGS84.
     **********************************************************************/
    MagneticModel(const std::string& datafile,
                  const Geocentric& earth = Geocentric::WGS84);

    /**
     * Evaluate the components of the magnetic field.
     *
     * @param[in] lat latitude of the point (degrees).
     * @param[in] lon longitude of the point (degrees).
     * @param[in] h the height of the point above the ellipsoid (meters).
     * @param[in] t the time (years).
     * @param[out] Bx the easterly component of the magnetic field (nanotesla).
     * @param[out] By the northerly component of the magnetic field (nanotesla).
     * @param[out] Bz the vertical (up) component of the magnetic field
     *   (nanotesla).
     **********************************************************************/
    void operator()(real lat, real lon, real h, real t,
                    real& Bx, real& By, real& Bz) const {
      real dummy;
      Field(lat, lon, h, t, false, Bx, By, Bz, dummy, dummy, dummy);
    }

    /**
     * Evaluate the components of the magnetic field and their time derivatives
     *
     * @param[in] lat latitude of the point (degrees).
     * @param[in] lon longitude of the point (degrees).
     * @param[in] h the height of the point above the ellipsoid (meters).
     * @param[in] t the time (years).
     * @param[out] Bx the easterly component of the magnetic field (nanotesla).
     * @param[out] By the northerly component of the magnetic field (nanotesla).
     * @param[out] Bz the vertical (up) component of the magnetic field
     *   (nanotesla).
     * @param[out] Bxt the rate of change of \e Bx (nT/yr).
     * @param[out] Byt the rate of change of \e Bx (nT/yr).
     * @param[out] Bzt the rate of change of \e Bx (nT/yr).
     **********************************************************************/
    void operator()(real lat, real lon, real h, real t,
                    real& Bx, real& By, real& Bz,
                    real& Bxt, real& Byt, real& Bzt) const {
      Field(lat, lon, h, t, true, Bx, By, Bz, Bxt, Byt, Bzt);
    }

  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_MAGNETICMODEL_HPP

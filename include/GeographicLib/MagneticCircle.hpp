/**
 * \file MagneticCircle.hpp
 * \brief Header for GeographicLib::MagneticCircle class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_MAGNETICCIRCLE_HPP)
#define GEOGRAPHICLIB_MAGNETICCIRCLE_HPP "$Id$"

#include <string>
#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/CircularEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Magnetic circle
   *
   * Evaluate the earth's magnetic field according to a circle.
   * See
   * - http://geomag.org/circles/index.html
   * - http://ngdc.noaa.gov/geomag/EMM/emm.shtml
   * - http://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
   **********************************************************************/

  class GEOGRAPHIC_EXPORT MagneticCircle {
  private:
    typedef Math::real real;

    real _a, _cphi, _sphi;
    CircularEngine _circa, _circb;
    MagneticCircle(real a, real cphi, real sphi,
                   const CircularEngine& circa, const CircularEngine& circb)
      : _a(a)
      , _cphi(cphi)
      , _sphi(sphi)
      , _circa(circa)
      , _circb(circb)
    {}

    void Field(real lon, bool diffp,
               real& Bx, real& By, real& Bz,
               real& Bxt, real& Byt, real& Bzt) const;

    friend class MagneticModel; // MagneticModel calls the private constructor

  public:

    /**
     * Evaluate the components of the magnetic field.
     *
     * @param[in] lon longitude of the point (degrees).
     * @param[out] Bx the easterly component of the magnetic field (nanotesla).
     * @param[out] By the northerly component of the magnetic field (nanotesla).
     * @param[out] Bz the vertical (up) component of the magnetic field
     *   (nanotesla).
     **********************************************************************/
    void operator()(real lon, real& Bx, real& By, real& Bz) const {
      real dummy;
      Field(lon, false, Bx, By, Bz, dummy, dummy, dummy);
    }

    /**
     * Evaluate the components of the magnetic field and their time derivatives
     *
     * @param[in] lon longitude of the point (degrees).
     * @param[out] Bx the easterly component of the magnetic field (nanotesla).
     * @param[out] By the northerly component of the magnetic field (nanotesla).
     * @param[out] Bz the vertical (up) component of the magnetic field
     *   (nanotesla).
     * @param[out] Bxt the rate of change of \e Bx (nT/yr).
     * @param[out] Byt the rate of change of \e Bx (nT/yr).
     * @param[out] Bzt the rate of change of \e Bx (nT/yr).
     **********************************************************************/
    void operator()(real lon, real& Bx, real& By, real& Bz,
                    real& Bxt, real& Byt, real& Bzt) const {
      Field(lon, true, Bx, By, Bz, Bxt, Byt, Bzt);
    }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_MAGNETICCIRCLE_HPP

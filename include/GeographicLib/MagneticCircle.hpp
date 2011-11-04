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
   * \brief Magnetic field on a circle of latitude
   *
   * Evaluate the earth's magnetic field according on a circle of constant
   * height and latitude.  This uses a CircleEngine to pre-evaluate the inner
   * sum of the spherical harmonic sum, allowing the values of the field at
   * different latitude to be evaluated rapidly.
   *
   * Use MagneticModel::Circle to create a MagneticCircle object.  (The
   * constructor for this class is private.)
   **********************************************************************/

  class GEOGRAPHIC_EXPORT MagneticCircle {
  private:
    typedef Math::real real;

    real _a, _cphi, _sphi, _t, _dt0;
    bool _interpolate;
    CircularEngine _circ0, _circ1;
    
    MagneticCircle(real a, real cphi, real sphi, real t, real dt0,
                   bool interpolate,
                   const CircularEngine& circ0, const CircularEngine& circ1)
      : _a(a)
      , _cphi(cphi)
      , _sphi(sphi)
      , _t(t)
      , _dt0(dt0)
      , _interpolate(interpolate)
      , _circ0(circ0)
      , _circ1(circ1)
    {}

    void Field(real lon, bool diffp,
               real& Bx, real& By, real& Bz,
               real& Bxt, real& Byt, real& Bzt) const;

    friend class MagneticModel; // MagneticModel calls the private constructor

  public:

    /**
     * Evaluate the components of the magnetic field at a particular longitude.
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
     * at a particular longitude.
     *
     * @param[in] lon longitude of the point (degrees).
     * @param[out] Bx the easterly component of the magnetic field (nanotesla).
     * @param[out] By the northerly component of the magnetic field (nanotesla).
     * @param[out] Bz the vertical (up) component of the magnetic field
     *   (nanotesla).
     * @param[out] Bxt the rate of change of \e Bx (nT/yr).
     * @param[out] Byt the rate of change of \e By (nT/yr).
     * @param[out] Bzt the rate of change of \e Bz (nT/yr).
     **********************************************************************/
    void operator()(real lon, real& Bx, real& By, real& Bz,
                    real& Bxt, real& Byt, real& Bzt) const {
      Field(lon, true, Bx, By, Bz, Bxt, Byt, Bzt);
    }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_MAGNETICCIRCLE_HPP

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
#include <sstream>
#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/SphericalHarmonic.hpp>
#include <GeographicLib/SphericalHarmonic1.hpp>

namespace GeographicLib {

  class MagneticCircle;

  /**
   * \brief Model of the earth's magnetic field
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
    std::string _name, _description, _date, _coeff, _meta;
    real _t0, _tmin, _tmax, _a, _hmin, _hmax;
    int _N, _M, _N1, _M1;
    Geocentric _earth;
    std::vector<real> _G, _H, _G1, _H1;
    void Field(real lat, real lon, real h, real t, bool diffp,
               real& Bx, real& By, real& Bz,
               real& Bxt, real& Byt, real& Bzt) const;
    SphericalHarmonic1 _harma;
    SphericalHarmonic _harmb;
    void ReadMetadata(const std::string& name);
    static bool ParseLine(const std::string& line,
                          std::string& key, std::string& val);
  public:
    /**
     * Construct a magnetic model.
     *
     * @param[in] name the name of the model.
     * @param[in] earth Geocentric object for converting coordinates; default
     *   Geocentric::WGS84.
     **********************************************************************/
    MagneticModel(const std::string& name,
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
     * @param[out] Byt the rate of change of \e By (nT/yr).
     * @param[out] Bzt the rate of change of \e Bz (nT/yr).
     **********************************************************************/
    void operator()(real lat, real lon, real h, real t,
                    real& Bx, real& By, real& Bz,
                    real& Bxt, real& Byt, real& Bzt) const {
      Field(lat, lon, h, t, true, Bx, By, Bz, Bxt, Byt, Bzt);
    }

    /**
     * Create a MagneticCircle object to allow the magnetic field at many
     * points with constant \e lat, \e h, and \e t and varying \e lon to be
     * computed efficienty.
     *
     * @param[in] lat latitude of the point (degrees).
     * @param[in] h the height of the point above the ellipsoid (meters).
     * @param[in] t the time (years).
     * @return a MagneticCircle object whose MagneticCircle::operator()(real
     *   lon) member function computes the field at a particular \e lon.
     *
     * If the field at several points on a circle of latitude need to be
     * calculated then instead of
     \code
  SphericalModel m(...);     // Create a magnetic model
  double lat = 33, lon0 = 44, dlon = 0.01, h = 1000, t = 2012;
  for (int i = 0; i <= 100; ++i) {
    real
      lon = lon0 + i * dlon, Bx, By, Bz;
    m(lat, lon, h, t, Bx, By, Bz);
    std::cout << lon << " " << Bx << " " << By << " " << Bz << "\n";
  }
     \endcode
     * use a MagneticCircle as in
     \code
  SphericalModel m(...);     // Create a magnetic model
  double lat = 33, lon0 = 44, dlon = 0.01, h = 1000, t = 2012;
  MagneticCircle c(m.Circle(lat, h, t)); // the MagneticCircle object
  for (int i = 0; i <= 100; ++i) {
    real
      lon = lon0 + i * dlon, Bx, By, Bz;
    c(lon, Bx, By, Bz);
    std::cout << lon << " " << Bx << " " << By << " " << Bz << "\n";
  }
     \endcode
     * For high-degree models, this will be substantially faster.
     **********************************************************************/
    MagneticCircle Circle(real lat, real h, real t) const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_MAGNETICMODEL_HPP

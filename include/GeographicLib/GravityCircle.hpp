/**
 * \file GravityCircle.hpp
 * \brief Header for GeographicLib::GravityCircle class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GRAVITYCIRCLE_HPP)
#define GEOGRAPHICLIB_GRAVITYCIRCLE_HPP "$Id$"

#include <string>
#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/CircularEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Gravity on a circle of latitude
   *
   * Evaluate the earth's gravity field on a circle of constant height and
   * latitude.  This uses a CircleEngine to pre-evaluate the inner sum of the
   * spherical harmonic sum, allowing the values of the field at different
   * latitude to be evaluated rapidly.
   *
   * Use GravityModel::Circle to create a GravityCircle object.  (The
   * constructor for this class is private.)
   **********************************************************************/

  class GEOGRAPHIC_EXPORT GravityCircle {
  private:
    typedef Math::real real;

    real _a, _f, _invR, _lat, _h, _cphi, _sphi, _amodel, _GMmodel, _dzonal0,
      _corrmult, _gamma0, _gamma, _frot;
    CircularEngine _gravitation, _disturbing, _correction;

    GravityCircle(real a, real f, real invR, real lat, real h,
                  real cphi, real sphi,
                  real amodel, real GMmodel, real dzonal0, real corrmult,
                  real gamma0, real gamma, real frot,
                  const CircularEngine& gravitation,
                  const CircularEngine& disturbing,
                  const CircularEngine& correction)
      : _a(a)
      , _f(f)
      , _invR(invR)
      , _lat(lat)
      , _h(h)
      , _cphi(cphi)
      , _sphi(sphi)
      , _amodel(amodel)
      , _GMmodel(GMmodel)
      , _dzonal0(dzonal0)
      , _corrmult(corrmult)
      , _gamma0(gamma0)
      , _gamma(gamma)
      , _frot(frot)
      , _gravitation(gravitation)
      , _disturbing(disturbing)
      , _correction(correction)
    {}

    friend class GravityModel; // GravityModel calls the private constructor

  public:

    /** \name Compute the gravitational field
     **********************************************************************/
    ///@{
    /**
     * Evaluate the components of the geogravity field at a particular
     * longitude.
     *
     * @param[in] lon longitude of the point (degrees).
     * @return the geoid height (meters).
     **********************************************************************/
    Math::real GeoidHeight(real lon) const throw();
    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value inherited from the GravityModel object used in the
     *   constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _a; }
    /**
     * @return \e f the flattening of the ellipsoid.  This is the value
     *   inherited from the GravityModel object used in the constructor.
     **********************************************************************/
    Math::real Flattening() const throw() { return _f; }
    /**
     * @return the latitude of the circle (degrees).
     **********************************************************************/
    Math::real Latitude() const throw() { return _lat; }
    /**
     * @return the height of the circle (meters).
     **********************************************************************/
    Math::real Height() const throw() { return _h; }
    ///@}
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_GRAVITYCIRCLE_HPP

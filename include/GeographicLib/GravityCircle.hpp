/**
 * \file GravityCircle.hpp
 * \brief Header for GeographicLib::GravityCircle class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GRAVITYCIRCLE_HPP)
#define GEOGRAPHICLIB_GRAVITYCIRCLE_HPP \
  "$Id$"

#include <string>
#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/CircularEngine.hpp>
#include <GeographicLib/GravityModel.hpp>

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
    /*
    enum captype {
      CAP_NONE   = GravityModel::CAP_NONE,
      CAP_G      = GravityModel::CAP_G,
      CAP_T      = GravityModel::CAP_T,
      CAP_DELTA  = GravityModel::CAP_DELTA,
      CAP_C      = GravityModel::CAP_C,
      CAP_GAMMA0 = GravityModel::CAP_GAMMA0,
      CAP_GAMMA  = GravityModel::CAP_GAMMA,
      CAP_ALL    = GravityModel::CAP_ALL,
    };
    */
    enum mask {
      NONE          = GravityModel::NONE,
      GRAVITY       = GravityModel::GRAVITY,
      DISTURBANCE   = GravityModel::DISTURBANCE,
      DISTPOTENTIAL = GravityModel::DISTPOTENTIAL,
      GEOIDHEIGHT   = GravityModel::GEOIDHEIGHT,
      ANOMALY       = GravityModel::ANOMALY,
      ALL           = GravityModel::ALL,
    };

    mask _caps;
    real _a, _f, _lat, _h, _Z, _P, _invR, _cpsi, _spsi,
      _cphi, _sphi, _amodel, _GMmodel, _dzonal0,
      _corrmult, _gamma0, _gamma, _frot;
    CircularEngine _gravitational, _disturbing, _correction;

    GravityCircle(mask caps, real a, real f, real lat, real h,
                  real Z, real P, real cphi, real sphi,
                  real amodel, real GMmodel, real dzonal0, real corrmult,
                  real gamma0, real gamma, real frot,
                  const CircularEngine& gravitational,
                  const CircularEngine& disturbing,
                  const CircularEngine& correction)
      : _caps(caps)
      , _a(a)
      , _f(f)
      , _lat(lat)
      , _h(h)
      , _Z(Z)
      , _P(P)
      , _invR(Math::hypot(_P, _Z))
      , _cpsi(_P * _invR)
      , _spsi(_Z * _invR)
      , _cphi(cphi)
      , _sphi(sphi)
      , _amodel(amodel)
      , _GMmodel(GMmodel)
      , _dzonal0(dzonal0)
      , _corrmult(corrmult)
      , _gamma0(gamma0)
      , _gamma(gamma)
      , _frot(frot)
      , _gravitational(gravitational)
      , _disturbing(disturbing)
      , _correction(correction)
    {}

    friend class GravityModel; // GravityModel calls the private constructor
    Math::real W(real clam, real slam,
                 real& gX, real& gY, real& gZ) const throw();
    Math::real V(real clam, real slam,
                 real& gX, real& gY, real& gZ) const throw();
    Math::real InternalT(real clam, real slam,
                         real& deltaX, real& deltaY, real& deltaZ,
                         bool gradp, bool correct) const throw();
  public:
    /** \name Compute the gravitational field
     **********************************************************************/
    ///@{
    /**
     * Evaluate the gravity.
     **********************************************************************/
    Math::real Gravity(real lon, real& gx, real& gy, real& gz) const throw();
    /**
     * Evaluate the gravity disturbance vector.
     **********************************************************************/
    Math::real Disturbance(real lon, real& deltax, real& deltay, real& deltaz)
      const throw();
    /**
     * Evaluate the geoid height.
     *
     * @param[in] lon longitude of the point (degrees).
     * @return the geoid height (meters).
     **********************************************************************/
    Math::real GeoidHeight(real lon) const throw();
    /**
     * Evaluate the components of the  gravity anomaly vector.
     **********************************************************************/
    void Anomaly(real lon, real& Dg01, real& xi, real& eta) const throw();
    Math::real W(real lon, real& gX, real& gY, real& gZ) const throw();
    Math::real V(real lon, real& GX, real& GY, real& GZ) const throw();
    Math::real T(real lon, real& deltaX, real& deltaY, real& deltaZ)
      const throw();
    Math::real T(real lon) const throw();

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

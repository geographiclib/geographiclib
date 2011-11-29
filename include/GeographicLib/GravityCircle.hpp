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
   * spherical harmonic sum, allowing the values of the field at several
   * different longitudes to be evaluated rapidly.
   *
   * Use GravityModel::Circle to create a GravityCircle object.  (The
   * constructor for this class is private.)
   **********************************************************************/

  class GEOGRAPHIC_EXPORT GravityCircle {
  private:
    typedef Math::real real;
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
    /**
     * A default constructor for the normal gravity.  This sets up an
     * uninitialized object which can be later replaced by the
     * GravityModel::Circle.
     **********************************************************************/
    GravityCircle() : _a(-1) {};

    /** \name Compute the gravitational field
     **********************************************************************/
    ///@{
    /**
     * Evaluate the gravity.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @param[out] gx the easterly component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gy the northerly component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gz the upward component of the acceleration
     *   (m s<sup>-2</sup>).
     * @return \e W the corresponding disturbing potential.
     *
     * The function includes the effects of the earth's rotation.
     **********************************************************************/
    Math::real Gravity(real lon, real& gx, real& gy, real& gz) const throw();
    /**
     * Evaluate the gravity disturbance vector.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @param[out] deltax the easterly component of the disturbance vector
     *   (m s<sup>-2</sup>).
     * @param[out] deltay the northerly component of the disturbance vector
     *   (m s<sup>-2</sup>).
     * @param[out] deltaz the upward component of the disturbance vector
     *   (m s<sup>-2</sup>).
     * @return \e T the corresponding disturbing potential.
     **********************************************************************/
    Math::real Disturbance(real lon, real& deltax, real& deltay, real& deltaz)
      const throw();
    /**
     * Evaluate the geoid height.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @return \e N the height of the geoid above the reference ellipsoid
     *   (meters).
     **********************************************************************/
    Math::real GeoidHeight(real lon) const throw();
    /**
     * Evaluate the components of the gravity anomaly vector.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @param[out] Dg01 the gravity anomaly (degrees).
     * @param[out] xi the northerly component of the deflection of the vertical
     *  (degrees).
     * @param[out] eta the easterly component of the deflection of the vertical
     *  (degrees).
     **********************************************************************/
    void Anomaly(real lon, real& Dg01, real& xi, real& eta) const throw();

    /**
     * Evaluate the components of the acceleration due to gravity and the
     * centrifugal acceleration in geocentric coordinates.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @param[out] gX the \e X component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gY the \e Y component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gZ the \e Z component of the acceleration
     *   (m s<sup>-2</sup>).
     * @return \e W = \e V + \e Phi the sum of the gravitational and
     *   centrifugal potentials (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real W(real lon, real& gX, real& gY, real& gZ) const throw() {
      real clam, slam;
      CircularEngine::cossin(lon, clam, slam);
      return W(clam, slam, gX, gY, gZ);
    }

    /**
     * Evaluate the components of the acceleration due to gravity in geocentric
     * coordinates.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @param[out] GX the \e X component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] GY the \e Y component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] GZ the \e Z component of the acceleration
     *   (m s<sup>-2</sup>).
     * @return \e V = \e W - \e Phi the gravitational potential
     *   (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real V(real lon, real& GX, real& GY, real& GZ) const throw() {
      real clam, slam;
      CircularEngine::cossin(lon, clam, slam);
      return V(clam, slam, GX, GY, GZ);
    }


    /**
     * Evaluate the components of the gravity disturbance in geocentric
     * coordinates.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @param[out] deltaX the \e X component of the gravity disturbance
     *   (m s<sup>-2</sup>).
     * @param[out] deltaY the \e Y component of the gravity disturbance
     *   (m s<sup>-2</sup>).
     * @param[out] deltaZ the \e Z component of the gravity disturbance
     *   (m s<sup>-2</sup>).
     * @return \e T = \e W - \e U the disturbing potentional (also called the
     *   anomalous potential) (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real T(real lon, real& deltaX, real& deltaY, real& deltaZ)
      const throw() {
      real clam, slam;
      CircularEngine::cossin(lon, clam, slam);
      return InternalT(clam, slam, deltaX, deltaY, deltaZ, true, true);
    }

    /**
     * Evaluate disturbing potentional in geocentric coordinates.
     *
     * @param[in] lon the geographic longitude (degrees).
     * @return \e T = \e W - \e U the disturbing potentional (also called the
     *   anomalous potential) (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real T(real lon) const throw() {
      real clam, slam, dummy;
      CircularEngine::cossin(lon, clam, slam);
      return InternalT(clam, slam, dummy, dummy, dummy, false, true);
    }

    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return true if the object has been initialized.
     **********************************************************************/
    bool Init() const throw() { return _a > 0; }
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value inherited from the GravityModel object used in the
     *   constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw()
    { return Init() ? _a : Math::NaN<real>(); }
    /**
     * @return \e f the flattening of the ellipsoid.  This is the value
     *   inherited from the GravityModel object used in the constructor.
     **********************************************************************/
    Math::real Flattening() const throw()
    { return Init() ? _f : Math::NaN<real>(); }
    /**
     * @return the latitude of the circle (degrees).
     **********************************************************************/
    Math::real Latitude() const throw()
    { return Init() ? _lat : Math::NaN<real>(); }
    /**
     * @return the height of the circle (meters).
     **********************************************************************/
    Math::real Height() const throw()
    { return Init() ? _h : Math::NaN<real>(); }
    ///@}
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_GRAVITYCIRCLE_HPP

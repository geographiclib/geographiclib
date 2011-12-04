/**
 * \file GravityModel.hpp
 * \brief Header for GeographicLib::GravityModel class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GRAVITYMODEL_HPP)
#define GEOGRAPHICLIB_GRAVITYMODEL_HPP \
  "$Id$"

#include <string>
#include <sstream>
#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/NormalGravity.hpp>
#include <GeographicLib/SphericalHarmonic.hpp>
#include <GeographicLib/SphericalHarmonic1.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs vector
#pragma warning (push)
#pragma warning (disable: 4251)
#endif

namespace GeographicLib {

  class GravityCircle;

  /**
   * \brief Model of the earth's gravity field
   *
   * Evaluate the earth's gravity field according to a model.
   *
   * Definitions and terminology:
   * See H+M, Sec 2-13.
   * - \e V = gravitational potential;
   * - \e Phi = rotational potential;
   * - \e W = \e V + \e Phi = \e T + \e U = total potential;
   * - <i>V</i><sub>0</sub> = normal gravitation potential;
   * - \e U = <i>V</i><sub>0</sub> + \e Phi = total normal potential;
   * - \e T = \e W - \e U = \e V - <i>V</i><sub>0</sub> = anomalous or
   *   disturbing potential;
   * - <b>g</b> = <b>grad</b> \e W = <b>gamma</b> + <b>delta</b>;
   * - <b>f</b> = <b>grad</b> \e Phi;
   * - <b>Gamma</b> = <b>grad</b> <i>V</i><sub>0</sub>;
   * - <b>gamma</b> = <b>grad</b> \e U;
   * - <b>delta</b> = <b>grad</b> \e T = gravity disturbance vector
   *   = <b>g</b><sub><i>P</i></sub> - <b>gamma</b><sub><i>P</i></sub>;
   * - delta \e g = gravity disturbance = \e g<sub><i>P</i></sub> - \e
   *   gamma<sub><i>P</i></sub>;
   * - Delta <b>g</b> = gravity anomaly vector =
   *   <b>g</b><sub><i>P</i></sub> - <b>gamma</b><sub><i>Q</i></sub>; here the
   *   line \e PQ is perpendicular to ellipsoid and the potential at \e P
   *   equals the normal potential at \e Q;
   * - Delta \e g = gravity anomaly = \e g<sub><i>P</i></sub> - \e
   *   gamma<sub><i>Q</i></sub>;
   * - (\e xi, \e eta) deflection of the vertical, the difference in
   *   directions of <b>g</b><sub><i>P</i></sub> and
   *   <b>gamma</b><sub><i>Q</i></sub>, \e xi = NS, \e eta = EW.
   * - \e X, \e Y, \e Z, geocentric coordinates;
   * - \e x, \e y, \e z, local cartesian coordinates used to denote the east,
   *   north and up directions.
   *
   * See \ref gravity for details of how to install the gravity model and the
   * data format.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT GravityModel {
  private:
    typedef Math::real real;
    friend class GravityCircle;
    static const int idlength_ = 8;
    std::string _name, _dir, _description, _date, _filename, _id;
    real _amodel, _GMmodel, _zeta0, _corrmult;
    SphericalHarmonic::normalization _norm;
    NormalGravity _earth;
    std::vector<real> _C, _S, _CC, _CS, _zonal;
    real _dzonal0;              // A left over contribution to _zonal.
    SphericalHarmonic _gravitational;
    SphericalHarmonic1 _disturbing;
    SphericalHarmonic _correction;
    void ReadMetadata(const std::string& name);
    Math::real InternalT(real X, real Y, real Z,
                         real& deltaX, real& deltaY, real& deltaZ,
                         bool gradp, bool correct) const throw();
    enum captype {
      CAP_NONE   = 0U,
      CAP_G      = 1U<<0,       // implies potentials W and V
      CAP_T      = 1U<<1,
      CAP_DELTA  = 1U<<2 | CAP_T, // delta implies T?
      CAP_C      = 1U<<3,
      CAP_GAMMA0 = 1U<<4,
      CAP_GAMMA  = 1U<<5,
      CAP_ALL    = 0x3FU,
    };

  public:

    /**
     * Bit masks for the capabilities to be given to the GravityCircle object
     * produced by Circle.
     **********************************************************************/
    enum mask {
      /**
       * No capabilities.
       * @hideinitializer
       **********************************************************************/
      NONE = 0U,
      /**
       * Allow calls to GravityCircle::Gravity, GravityCircle::W, and
       * GravityCircle::V.
       * @hideinitializer
       **********************************************************************/
      GRAVITY = CAP_G,
      /**
       * Allow calls to GravityCircle::Disturbance and GravityCircle::T.
       * @hideinitializer
       **********************************************************************/
      DISTURBANCE = CAP_DELTA,
      /**
       * Allow calls to GravityCircle::T(real lon) (i.e., computing the
       * disturbing potential and not the gravity disturbance vector).
       * @hideinitializer
       **********************************************************************/
      DISTURBING_POTENTIAL = CAP_T,
      /**
       * Allow calls to GravityCircle::SphericalAnomaly.
       * @hideinitializer
       **********************************************************************/
      SPHERICAL_ANOMALY = CAP_DELTA | CAP_GAMMA,
      /**
       * Allow calls to GravityCircle::GeoidHeight.
       * @hideinitializer
       **********************************************************************/
      GEOID_HEIGHT = CAP_T | CAP_C | CAP_GAMMA0,
      /**
       * All capabilities.
       * @hideinitializer
       **********************************************************************/
      ALL = CAP_ALL,
    };
    /** \name Setting up the gravity model
     **********************************************************************/
    ///@{
    /**
     * Construct a gravity model.
     *
     * @param[in] name the name of the model.
     * @param[in] path (optional) directory for data file.
     *
     * A filename is formed by appending ".egm" (World Gravity Model) to the
     * name.  If \e path is specified (and is non-empty), then the file is
     * loaded from directory, \e path.  Otherwise the path is given by
     * DefaultGravityPath().  This may throw an exception because the file does
     * not exist, is unreadable, or is corrupt.
     *
     * This file contains the metadata which specifies the properties of the
     * model.  The coefficients for the spherical harmonic sums are obtained
     * from a file obtained by appending ".cof" to metadata file (so the
     * filename ends in ".egm.cof").
     **********************************************************************/
    GravityModel(const std::string& name, const std::string& path = "");
    ///@}

    /** \name Compute gravity in geodetic coordinates
     **********************************************************************/
    ///@{
    /**
     * Evaluate the gravity at an arbitrary point above (or below) the
     * ellipsoid.
     *
     * @param[in] lat the geographic latitude (degrees).
     * @param[in] lon the geographic longitude (degrees).
     * @param[in] h the height above the ellipsoid (meters).
     * @param[out] gx the easterly component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gy the northerly component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gz the upward component of the acceleration
     *   (m s<sup>-2</sup>); this is usually negative.
     * @return \e W the sum of the gravitational and centrifugal potentials.
     *
     * The function includes the effects of the earth's rotation.
     **********************************************************************/
    Math::real Gravity(real lat, real lon, real h,
                       real& gx, real& gy, real& gz) const throw();

    /**
     * Evaluate the gravity disturbance vector at an arbitrary point above (or
     * below) the ellipsoid.
     *
     * @param[in] lat the geographic latitude (degrees).
     * @param[in] lon the geographic longitude (degrees).
     * @param[in] h the height above the ellipsoid (meters).
     * @param[out] deltax the easterly component of the disturbance vector
     *   (m s<sup>-2</sup>).
     * @param[out] deltay the northerly component of the disturbance vector
     *   (m s<sup>-2</sup>).
     * @param[out] deltaz the upward component of the disturbance vector
     *   (m s<sup>-2</sup>).
     * @return \e T the corresponding disturbing potential.
     **********************************************************************/
    Math::real Disturbance(real lat, real lon, real h,
                           real& deltax, real& deltay, real& deltaz)
      const throw();

    /**
     * Evaluate the geoid height.
     *
     * @param[in] lat the geographic latitude (degrees).
     * @param[in] lon the geographic longitude (degrees).
     * @return \e N the height of the geoid above the ReferenceEllipsoid()
     *   (meters).
     *
     * This calls NormalGravity::U for ReferenceEllipsoid().
     **********************************************************************/
    Math::real GeoidHeight(real lat, real lon) const throw();

    /**
     * Evaluate the components of the gravity anomaly vector using the
     * spherical approximation.
     *
     * @param[in] lat the geographic latitude (degrees).
     * @param[in] lon the geographic longitude (degrees).
     * @param[in] h the height above the ellipsoid (meters).
     * @param[out] Dg01 the gravity anomaly (m s<sup>-2</sup>).
     * @param[out] xi the northerly component of the deflection of the vertical
     *  (degrees).
     * @param[out] eta the easterly component of the deflection of the vertical
     *  (degrees).
     **********************************************************************/
    void SphericalAnomaly(real lat, real lon, real h,
                          real& Dg01, real& xi, real& eta) const throw();
    ///@}

    /** \name Compute gravity in geocentric coordinates
     **********************************************************************/
    ///@{
    /**
     * Evaluate the components of the acceleration due to gravity and the
     * centrifugal acceleration in geocentric coordinates.
     *
     * @param[in] X geocentric coordinate of point (meters).
     * @param[in] Y geocentric coordinate of point (meters).
     * @param[in] Z geocentric coordinate of point (meters).
     * @param[out] gX the \e X component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gY the \e Y component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gZ the \e Z component of the acceleration
     *   (m s<sup>-2</sup>).
     * @return \e W = \e V + \e Phi the sum of the gravitational and
     *   centrifugal potentials (m<sup>2</sup> s<sup>-2</sup>).
     *
     * This calls NormalGravity::U for  ReferenceEllipsoid().
     **********************************************************************/
    Math::real W(real X, real Y, real Z,
                 real& gX, real& gY, real& gZ) const throw();

    /**
     * Evaluate the components of the acceleration due to gravity in geocentric
     * coordinates.
     *
     * @param[in] X geocentric coordinate of point (meters).
     * @param[in] Y geocentric coordinate of point (meters).
     * @param[in] Z geocentric coordinate of point (meters).
     * @param[out] GX the \e X component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] GY the \e Y component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] GZ the \e Z component of the acceleration
     *   (m s<sup>-2</sup>).
     * @return \e V = \e W - \e Phi the gravitational potential
     *   (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real V(real X, real Y, real Z,
                 real& GX, real& GY, real& GZ) const throw();

    /**
     * Evaluate the components of the gravity disturbance in geocentric
     * coordinates.
     *
     * @param[in] X geocentric coordinate of point (meters).
     * @param[in] Y geocentric coordinate of point (meters).
     * @param[in] Z geocentric coordinate of point (meters).
     * @param[out] deltaX the \e X component of the gravity disturbance
     *   (m s<sup>-2</sup>).
     * @param[out] deltaY the \e Y component of the gravity disturbance
     *   (m s<sup>-2</sup>).
     * @param[out] deltaZ the \e Z component of the gravity disturbance
     *   (m s<sup>-2</sup>).
     * @return \e T = \e W - \e U the disturbing potential (also called the
     *   anomalous potential) (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real T(real X, real Y, real Z,
                 real& deltaX, real& deltaY, real& deltaZ) const throw()
    { return InternalT(X, Y, Z, deltaX, deltaY, deltaZ, true, true); }

    /**
     * Evaluate disturbing potential in geocentric coordinates.
     *
     * @param[in] X geocentric coordinate of point (meters).
     * @param[in] Y geocentric coordinate of point (meters).
     * @param[in] Z geocentric coordinate of point (meters).
     * @return \e T = \e W - \e U the disturbing potential (also called the
     *   anomalous potential) (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real T(real X, real Y, real Z) const throw() {
      real dummy;
      return InternalT(X, Y, Z, dummy, dummy, dummy, false, true);
    }

    /**
     * Evaluate the components of the acceleration due to normal gravity and the
     * centrifugal acceleration in geocentric coordinates.
     *
     * @param[in] X geocentric coordinate of point (meters).
     * @param[in] Y geocentric coordinate of point (meters).
     * @param[in] Z geocentric coordinate of point (meters).
     * @param[out] gammaX the \e X component of the normal acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gammaY the \e Y component of the normal acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gammaZ the \e Z component of the normal acceleration
     *   (m s<sup>-2</sup>).
     * @return \e U = <i>V</i><sub>0</sub> + \e Phi the sum of the
     *   normal gravitational and centrifugal potentials
     *   (m<sup>2</sup> s<sup>-2</sup>).
     *
     * This calls NormalGravity::U for  ReferenceEllipsoid().
     **********************************************************************/
    Math::real U(real X, real Y, real Z,
                 real& gammaX, real& gammaY, real& gammaZ) const throw()
    { return _earth.U(X, Y, Z, gammaX, gammaY, gammaZ); }

    /**
     * Evaluate the centrifugal acceleration in geocentric coordinates.
     *
     * @param[in] X geocentric coordinate of point (meters).
     * @param[in] Y geocentric coordinate of point (meters).
     * @param[out] fX the \e X component of the centrifugal acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] fY the \e Y component of the centrifugal acceleration
     *   (m s<sup>-2</sup>).
     * @return \e Phi the centrifugal potential (m<sup>2</sup> s<sup>-2</sup>).
     *
     * This calls NormalGravity::Phi for  ReferenceEllipsoid().
     **********************************************************************/
    Math::real Phi(real X, real Y, real& fX, real& fY) const throw()
    { return _earth.Phi(X, Y, fX, fY); }
    ///@}

    /** \name Compute gravity on a circle of constant latitude
     **********************************************************************/
    ///@{
    /**
     * Create a GravityCircle object to allow the gravity field at many points
     * with constant \e lat and \e h and varying \e lon to be computed
     * efficiently.
     *
     * @param[in] lat latitude of the point (degrees).
     * @param[in] h the height of the point above the ellipsoid (meters).
     * @param[in] caps bitor'ed combination of GravityModel::mask values
     *   specifying the capabilities of the resulting GravityCircle object.
     * @return a GravityCircle object whose member functions computes the
     *   gravitational field at a particular values of \e lon.
     *
     * The GravityModel::mask values are
     * - \e caps |= GravityModel::GRAVITY
     * - \e caps |= GravityModel::DISTURBANCE
     * - \e caps |= GravityModel::DISTURBING_POTENTIAL
     * - \e caps |= GravityModel::SPHERICAL_ANOMALY
     * - \e caps |= GravityModel::GEOID_HEIGHT
     * .
     * The default value of \e caps is GravityModel::ALL which turns on all the
     * capabilities.  If an unsupported function is invoked, it will return
     * NaNs.  Note that GravityModel::GEOID_HEIGHT will only be honored if \e h
     * = 0.
     *
     * If the field at several points on a circle of latitude need to be
     * calculated then instead of
     \code
  GravityModel g(...);          // Create a gravity model
  double lat = 33, lon0 = 44, dlon = 0.01;
  for (int i = 0; i <= 100; ++i) {
    real lon = lon0 + i * dlon;
    std::cout << lon << " " << g.GeoidHeight(lat, lon) << "\n";
  }
     \endcode
     * use a GravityCircle as in
     \code
  GravityModel g(...);          // Create a gravity model
  double lat = 33, lon0 = 44, dlon = 0.01;
  GravityCircle c(g.Circle(lat, 0)); // the GravityCircle object
  for (int i = 0; i <= 100; ++i) {
    real lon = lon0 + i * dlon;
    std::cout << lon << " " << c.GeoidHeight(lon) << "\n";
  }
     \endcode
     * For high-degree models, this will be substantially faster.
     **********************************************************************/
    GravityCircle Circle(real lat, real h, unsigned caps = ALL) const;
    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{

    /**
     * @return the NormalGravity object for the reference ellipsoid.
     **********************************************************************/
    const NormalGravity& ReferenceEllipsoid() const throw() { return _earth; }

    /**
     * @return the description of the gravity model, if available, in the data
     *   file; if absent, return "NONE".
     **********************************************************************/
    const std::string& Description() const throw() { return _description; }

    /**
     * @return date of the model; if absent, return "UNKNOWN".
     **********************************************************************/
    const std::string& DateTime() const throw() { return _date; }

    /**
     * @return full file name used to load the gravity model.
     **********************************************************************/
    const std::string& GravityFile() const throw() { return _filename; }

    /**
     * @return "name" used to load the gravity model (from the first argument
     *   of the constructor, but this may be overridden by the model file).
     **********************************************************************/
    const std::string& GravityModelName() const throw() { return _name; }

    /**
     * @return directory used to load the gravity model.
     **********************************************************************/
    const std::string& GravityModelDirectory() const throw() { return _dir; }

    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _earth.MajorRadius(); }

    /**
     * @return \e GM the mass constant of the model
     *   (m<sup>3</sup> s<sup>-2</sup>); this is the product of \e G the
     *   gravitational constant and \e M the mass of the earth (usually
     *   including the mass of the earth's atmosphere).
     **********************************************************************/
    Math::real MassConstant() const throw() { return _GMmodel; }

    /**
     * @return \e GM the mass constant of the ReferenceEllipsoid()
     *   (m<sup>3</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real ReferenceMassConstant() const throw()
    { return _earth.MassConstant(); }

    /**
     * @return \e omega the angular velocity of the model and the
     *   ReferenceEllipsoid() (rad s<sup>-1</sup>).
     **********************************************************************/
    Math::real AngularVelocity() const throw()
    { return _earth.AngularVelocity(); }

    /**
     * @return \e f the flattening of the ellipsoid.
     **********************************************************************/
    Math::real Flattening() const throw() { return _earth.Flattening(); }
    ///@}

    /**
     * @return the default path for gravity model data files.
     *
     * This is the value of the environment variable GRAVITY_PATH, if set;
     * otherwise, it is $GEOGRAPHICLIB_DATA/gravity if the environment variable
     * GEOGRAPHICLIB_DATA is set; otherwise, it is a compile-time default
     * (/usr/local/share/GeographicLib/gravity on non-Windows systems and
     * C:/Documents and Settings/All Users/Application
     * Data/GeographicLib/gravity on Windows systems).
     **********************************************************************/
    static std::string DefaultGravityPath();

    /**
     * @return the default name for the gravity model.
     *
     * This is the value of the environment variable GRAVITY_NAME, if set,
     * otherwise, it is "egm96".  The GravityModel class does not use
     * this function; it is just provided as a convenience for a calling
     * program when constructing a GravityModel object.
     **********************************************************************/
    static std::string DefaultGravityName();
  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_GRAVITYMODEL_HPP

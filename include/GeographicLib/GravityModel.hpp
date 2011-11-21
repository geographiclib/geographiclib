/**
 * \file GravityModel.hpp
 * \brief Header for GeographicLib::GravityModel class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GRAVITYMODEL_HPP)
#define GEOGRAPHICLIB_GRAVITYMODEL_HPP "$Id$"

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

  // class GravityCircle;

  /**
   * \brief Model of the earth's gravity field
   *
   * Evaluate the earth's gravity field according to a model.
   *
   * See xxref gravity for details of how to install the gravity model and the
   * data format.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT GravityModel {
  private:
    typedef Math::real real;
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
  public:

    /** \name Setting up the gravity model
     **********************************************************************/
    ///@{
    /**
     * Construct a gravity model.
     *
     * @param[in] name the name of the model.
     * @param[in] path (optional) directory for data file.
     *
     * A filename is formed by appending ".egm" (World Gravity Model) to
     * the name.  If \e path is specified (and is non-empty), then the file is
     * loaded from directory, \e path.  Otherwise the path is given by the
     * GRAVITY_PATH environment variable.  If that is undefined, a
     * compile-time default path is used
     * (/usr/local/share/GeographicLib/gravity on non-Windows systems and
     * C:/Documents and Settings/All Users/Application
     * Data/GeographicLib/gravity on Windows systems).  This may throw an
     * exception because the file does not exist, is unreadable, or is corrupt.
     *
     * This file contains the metadata which specifies the properties of the
     * model.  The coefficients for the spherical harmonic sums are obtained
     * from a file obtained by appending ".cof" to metadata file (so the
     * filename ends in ".egm.cof").
     **********************************************************************/
    GravityModel(const std::string& name, const std::string& path = "");
    ///@}

    /** \name Compute the gravity field
     **********************************************************************/
    ///@{
    /**
     * Evaluate the components of the gravity.
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
     *   <b>g</b><sub><i>P</i></sub> - <b>gamma</b><sub><i>Q</i></sub>, (where
     *   \e Q is on ellispoid and \e PQ is perpendicular to ellipsoid);
     * - Delta \e g = gravity anomaly = \e g<sub><i>P</i></sub> - \e
     *   gamma<sub><i>Q</i></sub>;
     * - (\e xi, \e eta) deflection of the vertical, the difference in
     *   directions of <b>g</b><sub><i>P</i></sub> and
     *   <b>gamma</b><sub><i>Q</i></sub>, \e xi = NS, \e eta = EW.
     **********************************************************************/
    Math::real Geoid(real lat, real lon) const throw();
    Math::real Disturbing(real lat, real lon, real h) const throw();
    Math::real Disturbing(real lat, real lon, real h,
                          real& gx, real& gy, real& gz) const throw();
    Math::real Gravitational(real lat, real lon, real h,
                             real& gx, real& gy, real& gz) const throw();
    Math::real Normal(real lat, real lon, real h,
                      real& gx, real& gy, real& gz) const throw();
    Math::real Rotational(real lat, real h,
                          real& gy, real& gz) const throw();
    Math::real Total(real lat, real lon, real h,
                     real& gx, real& gy, real& gz) const throw();
    Math::real T(real X, real Y, real Z) const throw() {
      real dummy;
      return InternalT(X, Y, Z, dummy, dummy, dummy, false, true);
    }
    Math::real T(real X, real Y, real Z,
                 real& deltaX, real& deltaY, real& deltaZ) const throw() {
      return InternalT(X, Y, Z, deltaX, deltaY, deltaZ, true, true);
    }
    void Anomaly(real lat, real lon, real h,
                 real& Dg01, real& xi, real& eta) const throw();
    Math::real V(real X, real Y, real Z,
                 real& GX, real& GY, real& GZ) const throw();
    Math::real W(real X, real Y, real Z,
                 real& gX, real& gY, real& gZ) const throw();
    /**
     * Evaluate the components of the acceleration due to normal gravity and the
     * centrifugal acceleration in geocentric coordinates.
     *
     * @param[in] X geocentric coordinate of point (meters).
     * @param[in] Y geocentric coordinate of point (meters).
     * @param[in] Z geocentric coordinate of point (meters).
     * @param[out] gammaX the \e X component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gammaY the \e Y component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gammaZ the \e Z component of the acceleration
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
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value of \e a inherited from the Geocentric object used in the
     *   constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _earth.MajorRadius(); }

    /**
     * @return \e f the flattening of the ellipsoid.  This is the value
     *   inherited from the Geocentric object used in the constructor.
     **********************************************************************/
    Math::real Flattening() const throw() { return _earth.Flattening(); }
    ///@}

    /**
     * @return the default path for gravity model data files.
     *
     * This is the value of the environment variable GRAVITY_PATH, if set,
     * otherwise, it is a compile-time default.
     **********************************************************************/

    static std::string DefaultGravityPath();

    /**
     * @return the default name for the gravity model.
     *
     * This is the value of the environment variable GRAVITY_NAME, if set,
     * otherwise, it is "egm2008".  The GravityModel class does not use
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

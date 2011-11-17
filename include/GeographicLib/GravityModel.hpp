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
                         bool gradp) const throw();
    Math::real InternalV(real X, real Y, real Z,
                         real& gX, real& gY, real& gZ,
                         bool gradp) const throw();
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
      return InternalT(X, Y, Z, dummy, dummy, dummy, false);
    }
    Math::real T(real X, real Y, real Z,
                 real& deltaX, real& deltaY, real& deltaZ) const throw() {
      return InternalT(X, Y, Z, deltaX, deltaY, deltaZ, true);
    }
    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{
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

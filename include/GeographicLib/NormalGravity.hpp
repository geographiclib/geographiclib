/**
 * \file NormalGravity.hpp
 * \brief Header for GeographicLib::NormalGravity class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_NORMALGRAVITY_HPP)
#define GEOGRAPHICLIB_NORMALGRAVITY_HPP "$Id$"

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/SphericalHarmonic.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <vector>

namespace GeographicLib {

  /**
   * \brief The normal gravity of the earth
   *
   * "Normal" gravity refers to an idealization of the earth which is modeled
   * as an rotating ellipsoid.  The eccentricity of the ellipsoid, the rotation
   * speed, and the distribution of mass within the ellipsoid are such that the
   * surface of the ellipsoid is a surface of constant potential (gravitational
   * plus centrifugal).  The acceleration due to gravity is therefore normal to
   * the surface of the ellipsoid.
   *
   * There is a closed solution to this problem which is implemented here.
   * Series "approximations" are only used to evaluate certain combinations of
   * elementary functions where use of the closed expression results in a loss
   * of accuracy for small arguments due to cancellation of the two leading
   * terms.  However these series include sufficient terms to give full machine
   * precision.
   *
   * References:
   * - W. A. Heiskanen and H. Moritz, Physical Geodesy (Freeman, San
   *   Fransisco, 1967), Secs. 2-7, 2-8 (2-9, 2-10), 6-2 (6-3).
   * - H. Moritz, Geodetic Reference System 1980, J. Geod. 54(3), 395-405
   *   (1980) http://dx.doi.org/10.1007/BF02521480
   **********************************************************************/

  class GEOGRAPHIC_EXPORT NormalGravity {
  private:
    static const int maxit_ = 10;
    static const int N_ = 14;
    typedef Math::real real;
    real _a, _GM, _omega, _J2, _omega2, _aomega2;
    real _e2, _ep2, _f, _b, _E, _U0, _gammae, _gammap, _q0, _m, _k, _fstar;
    std::vector<real> _C;
    SphericalHarmonic _harm;
    GeographicLib::Geocentric _earth;
    static Math::real qf(real ep2) throw();
    static Math::real qpf(real ep2) throw();
    Math::real Jn(int n) const throw();
  public:

    /** \name Setting up the normal gravity
     **********************************************************************/
    ///@{
    /**
     * Constructor for the normal gravity.
     *
     * @param[in] a equatorial radius (meters).
     * @param[in] GM gravitational constant radius
     *   (meters<sup>3</sup>/seconds<sup>2</sup>).
     * @param[in] omega the angular velociry (rad s<sup>-1</sup>).
     * @param[in] J2 dynamical form factor or the flattening \e f (depending on
     *   the value of \e flatp).
     * @param[in] flatp if true, then the flattening \e f is specified as the
     *   4th argument instead of <i>J</i><sub>2</sub> (default = false).
     *
     * The shape of the ellipsoid can be given in one of two ways:
     * - geometrically (\e flatp = true), the 4th argument represents the
     *   flattening \e f = (\e a - \e b) / \e a, where \e a and \e b are the
     *   equatorial radius and the polar semi-axis.
     * - physically (\e flatp = false, the default), the 4th argument
     *   represents the dynamical form factor <i>J</i><sub>2</sub> = (\e C - \e
     *   A) / <i>Ma</i><sup>2</sup>, where \e A and \e C are the equatorial and
     *   polar moments of inertia and \e M is the mass of the earth.
     **********************************************************************/
    NormalGravity(real a, real GM, real omega, real J2, bool flatp = false);

    /**
     * A default constructor for the normal gravity.
     **********************************************************************/
    NormalGravity() : _a(-1) {};
    ///@}

    /** \name Compute the gravity
     **********************************************************************/
    ///@{
    /**
     * Evaluate the gravity on the surface of the ellipsoid.
     *
     * @param[in] lat the geographic latitude (degrees).
     * @return \e g the acceleration due to gravity, positive downwards
     *   (m s<sup>-2</sup>).
     *
     * Due to the axial symmetry of the ellipsoid, the result is independent of
     * the value of the longitude.  This acceleration is normal to the surface
     * of the ellipsoid.  It includes the effects of the earth's rotation.
     **********************************************************************/
    Math::real SurfaceGravity(real lat) const throw();

    /**
     * Evaluate the gravity at an arbitrary point above (or below) the
     * ellipsoid.
     *
     * @param[in] lat the geographic latitude (degrees).
     * @param[in] h the height above the ellipsoid (meters).
     * @param[out] gy the northerly component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gz the upward component of the acceleration
     *   (m s<sup>-2</sup>).  (This is usually negative.)
     * @return \e g the magnitude acceleration due to gravity
     *
     * Due to the axial symmetry of the ellipsoid, the result is independent of
     * the value of the longitude and the easterly component of the
     * acceleration vanishes.  The function includes the effects of the earth's
     * rotation.  When \e h = 0, this function gives \e gy = 0 and the returned
     * value matches that of NormalGravity::SurfaceGravity.
     **********************************************************************/
    Math::real Gravity(real lat, real h, real& gy, real& gz) const throw();

    /**
     * Evaluate the components of the acceleration due to gravity and the
     * centrifugal acceleration in geocentric coordinates.
     *
     * @param[in] x geocentric coordinate of point (meters).
     * @param[in] y geocentric coordinate of point (meters).
     * @param[in] z geocentric coordinate of point (meters).
     * @param[out] gx the \e x component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gy the \e y component of the acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gz the \e z component of the acceleration
     *   (m s<sup>-2</sup>).
     * @return \e U the sum of the gravitational and centrifugal potentials
     *   (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real U(real x, real y, real z,
                 real& gx, real& gy, real& gz) const throw();
    /**
     * The same as NormalGravity::U, but evaluated with a finite set of
     * spherical harmonics.
     **********************************************************************/
    Math::real Useries(real x, real y, real z,
                       real& gx, real& gy, real& gz) const throw();

    /**
     * Evaluate the components of the acceleration due to gravity alone in
     * geocentric coordinates.
     *
     * @param[in] x geocentric coordinate of point (meters).
     * @param[in] y geocentric coordinate of point (meters).
     * @param[in] z geocentric coordinate of point (meters).
     * @param[out] gx the \e x component of the acceleration due to gravity
     *   (m s<sup>-2</sup>).
     * @param[out] gy the \e y component of the acceleration due to gravity
     *   (m s<sup>-2</sup>).
     * @param[out] gz the \e z component of the acceleration due to gravity
     *   (m s<sup>-2</sup>).
     * @return the gravitational potential (m<sup>2</sup> s<sup>-2</sup>).
     *
     * This function excludes the centrifugal acceleration and is appropriate
     * to use for space applications.  In terrestrial applications, the
     * function NormalGravity::U (which includes this effect) should usually be
     * used.
     **********************************************************************/
    Math::real V(real x, real y, real z,
                 real& gx, real& gy, real& gz) const throw();

    /**
     * The same as NormalGravity::V, but evaluated with a finite set of
     * spherical harmonics.
     **********************************************************************/
    Math::real Vseries(real x, real y, real z,
                       real& gx, real& gy, real& gz) const throw();
    /**
     * Evaluate the centrifugal acceleration in geocentric coordinates.
     *
     * @param[in] x geocentric coordinate of point (meters).
     * @param[in] y geocentric coordinate of point (meters).
     * @param[out] gx the \e x component of the centrifugal acceleration
     *   (m s<sup>-2</sup>).
     * @param[out] gy the \e y component of the centrifugal acceleration
     *   (m s<sup>-2</sup>).
     * @return \e Phi the centrifugal potential (m<sup>2</sup> s<sup>-2</sup>).
     *
     * \e Phi is independent of \e z, thus \e gz = 0.  This function
     * NormalGravity::U sums the results of NormalGravity::V and
     * NormalGravity::Phi.
     **********************************************************************/
    Math::real Phi(real x, real y, real& gx, real& gy) const throw();
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
     *   the value used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw()
    { return Init() ? _a : Math::NaN<real>(); }

    /**
     * @return \e GM the gravitational constation of the ellipsoid
     *   (m<sup>3</sup> s<sup>-2</sup>).  This is the value used in the
     *   constructor.
     **********************************************************************/
    Math::real GravitationalConstant() const throw()
    { return Init() ? _GM : Math::NaN<real>(); }

    /**
     * @return \e J<sub>n</sub> the dynamical form factors of the ellipsoid.
     *
     * If \e n = 2 (the default), this is the value of <i>J</i><sub>2</sub>
     * used in the constructor.  Otherwise it is the zonal coefficient of the
     * Legendre harmonic sum of the normal gravitational potential.  Note that
     * \e J<sub>n</sub> = 0 if \e is odd.  In most gravity applications, fully
     * normalized Legendre functions are used and the corresponding coefficient
     * is <i>C</i><sub><i>n</i>0</sub> = -\e J<sub>n</sub> / sqrt(2 \e n + 1).
     **********************************************************************/
    Math::real DynamicalFormFactor(int n = 2) const throw()
    { return Init() ? ( n == 2 ? _J2 : Jn(n)) : Math::NaN<real>(); }

    /**
     * @return \e omega the angular velocity of the ellipsoid
     *   (rad s<sup>-1</sup>).  This is the value used in the constructor.
     **********************************************************************/
    Math::real AngularVelocity() const throw()
    { return Init() ? _omega : Math::NaN<real>(); }

    /**
     * @return <i>f</i> the flattening of the ellipsoid (\e a - \e b)/\e a.
     **********************************************************************/
    Math::real Flattening() const throw() { return Init() ? _f : Math::NaN<real>(); }

    /**
     * @return <i>gamma</i><sub>e</sub> the normal gravity at equator
     *   (m s<sup>-2</sup>).
     **********************************************************************/
    Math::real EquatorialGravity() const throw()
    { return Init() ? _gammae : Math::NaN<real>(); }

    /**
     * @return <i>gamma</i><sub>p</sub> the normal gravity at poles
     *   (m s<sup>-2</sup>).
     **********************************************************************/
    Math::real PolarGravity() const throw()
    { return Init() ? _gammap : Math::NaN<real>(); }

    /**
     * @return <i>f*</i> the gravity flattening
     *   (<i>gamma</i><sub>p</sub> - <i>gamma</i><sub>e</sub>) /
     *   <i>gamma</i><sub>e</sub>.
     **********************************************************************/
    Math::real GravityFlattening() const throw()
    { return Init() ? _fstar : Math::NaN<real>(); }

    /**
     * @return <i>U</i><sub>0</sub> the constant normal potential for the
     *   surface of the ellipsoid (m<sup>2</sup> s<sup>-2</sup>).
     **********************************************************************/
    Math::real SurfacePotential() const throw()
    { return Init() ? _U0 : Math::NaN<real>(); }
    
    /**
     * @return the GeographicLib::Geocentric object used by this instance.
     **********************************************************************/
    const Geocentric& Earth() const throw() { return _earth; }
    ///@}

    /**
     * A global instantiation of NormalGravity for the WGS84 ellipsoid.
     **********************************************************************/
    static const NormalGravity WGS84;

    /**
     * A global instantiation of NormalGravity for the GRS80 ellipsoid.
     **********************************************************************/
    static const NormalGravity GRS80;
    void DumpConstants() const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_NORMALGRAVITY_HPP

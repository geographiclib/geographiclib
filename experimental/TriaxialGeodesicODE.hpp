/**
 * \file TriaxialGeodesicODE.hpp
 * \brief Header for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIALGEODESICODE_HPP)
#define GEOGRAPHICLIB_TRIAXIALGEODESICODE_HPP 1

#include <vector>
#include <array>
#include <utility>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Angle.hpp>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

// Boost's dense output expects numeric_limits<real>::digits to be a constant
// and not a function.  So we can't use GEOGRAPHICLIB_PRECISION == 5.  High
// precision results can be obtained with GEOGRAPHICLIB_PRECISION > 5, e.g.,
// GEOGRAPHICLIB_PRECISION == 256.
#if GEOGRAPHICLIB_PRECISION == 5
#  define GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT 0
#else
#  define GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT 1
#endif

// Boost bugs when using high precision:
//    https://github.com/boostorg/odeint/issues/40
//    https://github.com/boostorg/odeint/issues/75 (duplicate)
// fixed in
//    https://github.com/boostorg/odeint/pull/63
//    Commit 68950d8
//
// This will be included in Boost 1.85.  (Fedora 42 uses Boost 1.83.)
//
// In the meantime, put the patch for commit 68950d8 in
// /usr/include/boost/odeint.patch and applied it with
//    patch -p3 -b < odeint.patch
// -> patching file numeric/odeint/algebra/detail/extract_value_type.hpp
//
// Removed my temporary fix to
//    numeric/odeint/stepper/controlled_runge_kutta.hpp

#if __clang__
// Ignore clang warnings for boost headers
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wreorder-ctor"
#endif

#include <boost/numeric/odeint.hpp>
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>
#endif
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

namespace GeographicLib {
  /**
   * \brief Namespace for experimental components of %GeographicLib
   *
   * These routines are distributed as source code with %GeographicLib but are
   * not incorporated into the library itself.
   **********************************************************************/
  namespace experimental {

  /**
   * \brief The ODE solution of direct geodesic problem for triaxial ellipsoids.
   *
   * This determines the course of a geodesic by solving the equations of
   * motion for a particle sliding without friction on the surface of the
   * ellipsoid.  The solution is carried out in cartesian coordinates.  The
   * same approach was used by
   * <a href="https://doi.org/10.1515/jogs-2019-0001"> Panou and Korakitis
   * (2019)</a>.  Significant differences are:
   * * The code is provided.
   * * This method uses a high order method provided by Boost.  This allows
   *   reasonably high accuracy to be achieved using double precision.
   * * It solver optionally offers "dense" output.  It takes as large steps as
   *   possible while meeting the accuracy requirements.  The results at
   *   specific distances are then found by interpolation.  There is little
   *   penalty in requesting many waypoints.
   * * Because the ODE only "works" in one direction, the solver has to restart
   *   if the direction is reversed (but this is hidden from the user).  To
   *   simplify usage, you can provide a vector of distances to the solver and
   *   this is sorted appropriately before calling the underlying Boost
   *   integrator.  This vector can include negative as well as positive
   *   distances.
   * * This class can also, optionally, solve for the reduced length \e m12,
   *   and the geodesic scales \e M12 and \e M21.
   *
   * These routines are distributed as source code with %GeographicLib but,
   * because of the dependency on the Boost library, are not incorporated into
   * the library itself.
   *
   * The recommended way to solve the direct and indirect geodesic problems on
   * a triaxial ellipsoid is with the class Triaxial::Geodesic3.
   *
   * Geod3ODE.cpp is a utility which uses this class to solve direct geodesic
   * problems.  Use `Geod3ODE --help` for brief documentation.
   **********************************************************************/
  class TriaxialGeodesicODE {
  public:
    /**
     * A type to hold three-dimensional positions and velocities in cartesian
     * coordinates.
     **********************************************************************/
    using vec3 = Triaxial::Ellipsoid3::vec3;
  private:
    using real = Math::real;
    using vec6 = std::array<real, 6>;
    using vec10 = std::array<real, 10>;
    using ang = Angle;
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
    using dstep6 =
      boost::numeric::odeint::bulirsch_stoer_dense_out<vec6, real>;
    using dstep10 =
      boost::numeric::odeint::bulirsch_stoer_dense_out<vec10, real>;
#endif
    using step6 = boost::numeric::odeint::bulirsch_stoer<vec6, real>;
    using step10 = boost::numeric::odeint::bulirsch_stoer<vec10, real>;
    const Triaxial::Ellipsoid3 _t;
    const real _b, _eps;
    const vec3 _axesn, _axes2n;
    vec3  _r1, _v1;
    Angle _bet1, _omg1, _alp1;
    bool _extended, _dense, _normp;
    int _dir;
    long _nsteps;
    mutable long _intsteps;     // This is advanced by AccelXX which is const
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
    dstep6 _dstep6;
    dstep10 _dstep10;
#endif
    step6 _step6;
    step10 _step10;
    real _s;
    vec6 _y6;
    vec10 _y10;
    // These private versions of Accel and Norm assume normalized ellipsoid
    // with axes = axesn
    void Norm6(vec6& y) const;
    void Accel6(const vec6& y, vec6& yp) const;
    void Accel6N(const vec6& y, vec6& yp) const;
    void Norm10(vec10& y) const;
    void Accel10(const vec10& y, vec10& yp) const;
    void Accel10N(const vec10& y, vec10& yp) const;
    static std::vector<size_t> sort_indices(const std::vector<real>& v);
    void Reset();

  public:
    /**
     * Basic constructor specifying ellipsoid and starting point.
     *
     * @param[in] t the Ellipsoid3 object.
     * @param[in] R1 the starting position on the ellipsoid.
     * @param[in] V1 the starting velocity tangent to the ellipsoid.
     * @param[in] extended (default false), if true solve for reduced length
     *   and geodesic scale.
     * @param[in] dense (default false), if true use a dense solver allowing
     *   interpolated way points to be computed inexpensively.
     * @param[in] normp (default false), if true force the solution vector onto
     *   the ellipsoid when computing the acceleration.
     * @param[in] eps (default 0), if positive the error threshold for the
     *   integrator; otherwise use a good default value related to the machine
     *   precision.
     *
     * The values \e R1 and \e V1 are normalized to place \e R1 on the
     * ellipsoid and \e V1 tangent to the ellipsoid with unit speed.
     *
     * Internally, the integration scales the ellipsoid so that the median
     * semiaxis \b = 1.  The \e eps parameter is a measure of the absolution
     * error on this scaled ellipsoid.
     **********************************************************************/
    TriaxialGeodesicODE(const Triaxial::Ellipsoid3& t, vec3 R1, vec3 V1,
                        bool extended = false, bool dense = false,
                        bool normp = false, real eps = 0);
    /**
     * Constructor specifying just the ellipsoid.
     *
     * @param[in] t the Ellipsoid3 object.
     * @param[in] extended (default false), if true solve for reduced length
     *   and geodesic scale.
     * @param[in] dense (default false), if true use a dense solver allowing
     *   interpolated way points to be computed inexpensively.
     * @param[in] normp (default false), if true force the solution vector onto
     *   the ellipsoid when computing the acceleration.
     * @param[in] eps (default 0), if positive the error threshold for the
     *   integrator; otherwise use a good default value related to the machine
     *   precision.
     *
     * This form starts the geodesic at \e R1 = [\e a, 0, 0], \e V1 = [0, 0,
     * 1].
     **********************************************************************/
    TriaxialGeodesicODE(const Triaxial::Ellipsoid3& t,
                        bool extended = false, bool dense = false,
                        bool normp = false, real eps = 0);
    /**
     * Basic constructor specifying ellipsoid and starting point in ellipsoidal
     * coordinates.
     *
     * @param[in] t the Ellipsoid3 object.
     * @param[in] bet1 the starting latitude.
     * @param[in] omg1 the starting longitude.
     * @param[in] alp1 the starting azimuth.
     * @param[in] extended (default false), if true solve for reduced length
     *   and geodesic scale.
     * @param[in] dense (default false), if true use a dense solver allowing
     *   interpolated way points to be computed inexpensively.
     * @param[in] normp (default false), if true force the solution vector onto
     *   the ellipsoid when computing the acceleration.
     * @param[in] eps (default 0), if positive the error threshold for the
     *   integrator; otherwise use a good default value related to the machine
     *   precision.
     **********************************************************************/
    TriaxialGeodesicODE(const Triaxial::Ellipsoid3& t,
                        Angle bet1, Angle omg1, Angle alp1,
                        bool extended = false, bool dense = false,
                        bool normp = false, real eps = 0);
    /**
     * Find the position a given distance from the starting point.
     *
     * @param[in] s12 the distance between point 1 and point 2.
     * @param[out] R2 the position of point 2.
     * @param[out] V2 the velocity at point 2.
     * @return a pair of error estimates, the distance from the ellipsoid (in
     *   meters) and the deviation of the velocity from a unit tangential
     *   vector.
     *
     * The returned values \e R2 and \e V2 are normalized to place \e R2 on the
     * ellipsoid and \e V2 tangent to the ellipsoid with unit speed.
     **********************************************************************/
    std::pair<real, real> Position(real s12, vec3& R2, vec3& V2);
    /**
     * Find the position and differential quantities a given distance from the
     * starting point.
     *
     * @param[in] s12 the distance between point 1 and point 2.
     * @param[out] R2 the position of point 2.
     * @param[out] V2 the velocity at point 2.
     * @param[out] m12 the reduced length
     * @param[out] M12 the geodesic scale at point 2 relative to point 1.
     * @param[out] M21 the geodesic scale at point 1 relative to point 2.
     * @return a pair of error estimates, the distance from the ellipsoid (in
     *   meters) and the deviation of the velocity from a unit tangential
     *   vector.
     *
     * The returned values \e R2 and \e V2 are normalized to place \e R2 on the
     * ellipsoid and \e V2 tangent to the ellipsoid with unit speed.
     *
     * If the object was constructed with \e extended = false, NaNs are
     * returned for the differential quantities.
     **********************************************************************/
    std::pair<real, real> Position(real s12, vec3& R2, vec3& V2,
                                   real& m12, real& M12, real& M21);
    /**
     * Find the ellipsoidal coordinates a given distance from the starting
     * point.
     *
     * @param[in] s12 the distance between point 1 and point 2.
     * @param[out] bet2 the latitude at point 2.
     * @param[out] omg2 the longitude at point 2.
     * @param[out] alp2 the azimuth at point 2.
     * @return a pair of error estimates, the distance from the ellipsoid (in
     *   meters) and the deviation of the velocity from a unit tangential
     *   vector.
     **********************************************************************/
    std::pair<real, real> Position(real s12,
                                   Angle& bet2, Angle& omg2, Angle& alp2);
    /**
     * Find the ellipsoidal coordinates and differential quantities a given
     * distance from the starting point.
     *
     * @param[in] s12 the distance between point 1 and point 2.
     * @param[out] bet2 the latitude at point 2.
     * @param[out] omg2 the longitude at point 2.
     * @param[out] alp2 the azimuth at point 2.
     * @param[out] m12 the reduced length
     * @param[out] M12 the geodesic scale at point 2 relative to point 1.
     * @param[out] M21 the geodesic scale at point 1 relative to point 2.
     * @return a pair of error estimates, the distance from the ellipsoid (in
     *   meters) and the deviation of the velocity from a unit tangential
     *   vector.
     *
     * If the object was constructed with \e extended = false, NaNs are
     * returned for the differential quantities.
     **********************************************************************/
    std::pair<real, real> Position(real s12,
                                   Angle& bet2, Angle& omg2, Angle& alp2,
                                   real& m12, real& M12, real& M21);

    /**
     * Find the positions for a series of distances from the starting point.
     *
     * @param[in] s12 a vector of distances between point 1 and points 2.
     * @param[out] R2 a vector of positions of points 2.
     * @param[out] V2 a vector of velocities at points 2.
     *
     * Before starting the integration, the positive and negative vaules in \e
     * s12 are separated and then sorted in order of increasing magnitude.  The
     * results are placed back in the correct positions in the output vectors.
     * \e s12 can include NaNs; this can be used to "punctuate" the results.
     *
     * The returned values \e R2 and \e V2 are normalized to place \e R2 on the
     * ellipsoid and \e V2 tangent to the ellipsoid with unit speed.
     **********************************************************************/
    void Position(const std::vector<real>& s12,
                  std::vector<vec3>& R2, std::vector<vec3>& V2);
    /**
     * Find the positions and differential quantities for a series of distances
     * from the starting point.
     *
     * @param[in] s12 a vector of distances between point 1 and points 2.
     * @param[out] R2 a vector of positions of points 2.
     * @param[out] V2 a vector of velocities at points 2.
     * @param[out] m12 a vector of the reduced lengths.
     * @param[out] M12 a vector of the geodesic scales at points 2 relative to
     *   point 1.
     * @param[out] M21 a vector of the geodesic scales at point 1 relative to
     *   points 2.
     *
     * Before starting the integration, the positive and negative vaules in \e
     * s12 are separated and then sorted in order of increasing magnitude.  The
     * results are placed back in the correct positions in the output vectors.
     * \e s12 can include NaNs; this can be used to "punctuate" the results.
     *
     * The returned values \e R2 and \e V2 are normalized to place \e R2 on the
     * ellipsoid and \e V2 tangent to the ellipsoid with unit speed.
     *
     * If the object was constructed with \e extended = false, NaNs are
     * returned for the differential quantities.
     **********************************************************************/
    void Position(const std::vector<real>& s12,
                  std::vector<vec3>& R2, std::vector<vec3>& V2,
                  std::vector<real>& m12,
                  std::vector<real>& M12, std::vector<real>& M21);
    /**
     * Find the ellipsoidal coordinates for a series of distances from the
     * starting point.
     *
     * @param[in] s12 a vector of distances between point 1 and points 2.
     * @param[out] bet2 a vector of latitudes at points 2.
     * @param[out] omg2 a vector of longitudes at points 2.
     * @param[out] alp2 a vector of azimuths at points 2.
     *
     * Before starting the integration, the positive and negative vaules in \e
     * s12 are separated and then sorted in order of increasing magnitude.  The
     * results are placed back in the correct positions in the output vectors.
     * \e s12 can include NaNs; this can be used to "punctuate" the results.
     **********************************************************************/
    void Position(const std::vector<real>& s12,
                  std::vector<Angle>& bet2, std::vector<Angle>& omg2,
                  std::vector<Angle>& alp2);
    /**
     * Find the ellipsoidal coordinates and differential quantities for a
     * series of distances from the starting point.
     *
     * @param[in] s12 a vector of distances between point 1 and points 2.
     * @param[out] bet2 a vector of latitudes at points 2.
     * @param[out] omg2 a vector of longitudes at points 2.
     * @param[out] alp2 a vector of azimuths at points 2.
     * @param[out] m12 a vector of the reduced lengths.
     * @param[out] M12 a vector of the geodesic scales at points 2 relative to
     *   point 1.
     * @param[out] M21 a vector of the geodesic scales at point 1 relative to
     *   points 2.
     *
     * Before starting the integration, the positive and negative vaules in \e
     * s12 are separated and then sorted in order of increasing magnitude.  The
     * results are placed back in the correct positions in the output vectors.
     * \e s12 can include NaNs; this can be used to "punctuate" the results.
     **********************************************************************/
    void Position(const std::vector<real>& s12,
                  std::vector<Angle>& bet2, std::vector<Angle>& omg2,
                  std::vector<Angle>& alp2,
                  std::vector<real>& m12,
                  std::vector<real>& M12, std::vector<real>& M21);
    /**
     * Reset the starting point.
     *
     * @param[in] R1 the starting position on the ellipsoid.
     * @param[in] V1 the starting velocity tangent to the ellipsoid.
     *
     * The values \e R1 and \e V1 are normalized to place \e R1 on the
     * ellipsoid and \e V1 tangent to the ellipsoid with unit speed.
     **********************************************************************/
    void Reset(vec3 R1, vec3 V1);
    /**
     * Reset the starting point in ellipsoidal coordinates.
     *
     * @param[in] bet1 the starting latitude.
     * @param[in] omg1 the starting longitude.
     * @param[in] alp1 the starting azimuth.
     **********************************************************************/
    void Reset(Angle bet1, Angle omg1, Angle alp1);
    /**
     * @return the number of integration steps since the last Reset()
     **********************************************************************/
    long NSteps() const { return _nsteps; }
    /**
     * @return the number of calls to the acceleration routine since the last
     *   Reset()
     **********************************************************************/
    long IntSteps() const { return _intsteps; }
    /**
     * @return a pair for the current distance bracket.
     *
     * If the object was constructed with \e dense = true this gives the
     * current interval over which an interpolated waypoint can be found.
     * Otherwise the two distances are equal to the last distance calculation.
     **********************************************************************/
    std::pair<real, real> CurrentDistance() const;
    /**
     * @return whether the object was constructed with \e extended = true.
     **********************************************************************/
    bool Extended() const { return _extended; }
    /**
     * The current starting point.
     *
     * @param[out] R1 the starting position.
     * @param[out] V1 the starting velocity.
     **********************************************************************/
    void Position1(vec3& R1, vec3& V1) const { R1 = _r1; V1 = _v1; }
    /**
     * The current starting point in ellipsoidal coordinates.
     *
     * @param[out] bet1 the starting latitude.
     * @param[out] omg1 the starting longitude.
     * @param[out] alp1 the starting azimuth.
     **********************************************************************/
    void Position1(Angle& bet1, Angle& omg1, Angle& alp1) const {
      bet1 = _bet1; omg1 = _omg1; alp1 = _alp1;
    }
    /**
     * @return the Ellipsoid3 object used in the constructor.
     **********************************************************************/
    const Triaxial::Ellipsoid3& t() const { return _t; }
  };

} // namespace experimental
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALGEODESICODE_HPP

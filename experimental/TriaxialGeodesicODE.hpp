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

  class TriaxialGeodesicODE {
  public:
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
    mutable long _nsteps, _intsteps;
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

  public:
    TriaxialGeodesicODE(const Triaxial::Ellipsoid3& t,
                        bool extended = false, bool dense = false,
                        bool normp = false, real eps = 0);
    TriaxialGeodesicODE(const Triaxial::Ellipsoid3& t, vec3 r1, vec3 v1,
                        bool extended = false, bool dense = false,
                        bool normp = false, real eps = 0);
    TriaxialGeodesicODE(const Triaxial::Ellipsoid3& t, Angle bet1, Angle omg1, Angle alp1,
                        bool extended = false, bool dense = false,
                        bool normp = false, real eps = 0);
    std::pair<real, real> Position(real s12, vec3& r2, vec3& v2);
    std::pair<real, real> Position(real s12, vec3& r2, vec3& v2,
                                   real& m12, real& M12, real& M21);
    std::pair<real, real> Position(real s12,
                                   Angle& bet2, Angle& omg2, Angle& alp2);
    std::pair<real, real> Position(real s12,
                                   Angle& bet2, Angle& omg2, Angle& alp2,
                                   real& m12, real& M12, real& M21);

    void Position(const std::vector<real>& s12,
                  std::vector<vec3>& r2, std::vector<vec3>& v2);
    void Position(const std::vector<real>& s12,
                  std::vector<vec3>& r2, std::vector<vec3>& v2,
                  std::vector<real>& m12,
                  std::vector<real>& M12, std::vector<real>& M21);
    void Position(const std::vector<real>& s12,
                  std::vector<Angle>& bet2, std::vector<Angle>& omg2,
                  std::vector<Angle>& alp2);
    void Position(const std::vector<real>& s12,
                  std::vector<Angle>& bet2, std::vector<Angle>& omg2,
                  std::vector<Angle>& alp2,
                  std::vector<real>& m12,
                  std::vector<real>& M12, std::vector<real>& M21);
    void Reset();
    void Reset(vec3 r1, vec3 v1);
    void Reset(Angle bet1, Angle omg1, Angle alp1);
    long NSteps() const { return _nsteps; }
    void NSteps(long nsteps) const { _nsteps = nsteps; }
    long IntSteps() const { return _intsteps; }
    void IntSteps(long intsteps) const { _intsteps = intsteps; }
    std::pair<real, real> CurrentDistance() const;
    bool Extended() const { return _extended; }
    void Position1(vec3& r1, vec3& v1) const { r1 = _r1; v1 = _v1; }
    void Position1(Angle& bet1, Angle& omg1, Angle& alp1) const {
      bet1 = _bet1; omg1 = _omg1; alp1 = _alp1;
    }
    const Triaxial::Ellipsoid3& t() const { return _t; }
  };

} // namespace experimental
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALGEODESICODE_HPP

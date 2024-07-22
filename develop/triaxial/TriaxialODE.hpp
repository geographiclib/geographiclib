/**
 * \file TriaxialODE.hpp
 * \brief Header for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIALODE_HPP)
#define GEOGRAPHICLIB_TRIAXIALODE_HPP 1

#include <array>
#include <GeographicLib/Constants.hpp>

#if !defined(GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT)
#  if GEOGRAPHICLIB_PRECISION < 5
#    define GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT 1
#  else
#    define GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT 0
#  endif
#endif

#include "Angle.hpp"
#include "Triaxial.hpp"
#include <boost/numeric/odeint.hpp>
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>
#endif
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

namespace GeographicLib {

  class GEOGRAPHICLIB_EXPORT TriaxialODE {
  private:
    typedef Math::real real;
    typedef Triaxial::vec3 vec3;
    typedef std::array<real, 6> vec6;
    typedef std::array<real, 10> vec10;
    typedef Angle ang;
#if GEOGRAPHICLIB_BOOST_ODE_DENSE_OUT
    typedef
    boost::numeric::odeint::bulirsch_stoer_dense_out<vec6, real> dstep6;
    typedef
    boost::numeric::odeint::bulirsch_stoer_dense_out<vec10, real> dstep10;
#endif
    typedef
    boost::numeric::odeint::bulirsch_stoer<vec6, real> step6;
    typedef
    boost::numeric::odeint::bulirsch_stoer<vec10, real> step10;
    Triaxial _t;
    real _b, _eps;
    vec3 _axesn, _axes2n, _r1, _v1;
    bool _extended, _interp, _dense;
    int _dir;
    long _nsteps;
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
    vec6 Accel(const vec6& y) const;
    void Norm(vec6& y) const;
    // These private versions of Accel and Norm assume normalized ellipsoid
    // with axes = axesn
    vec10 Accel(const vec10& y) const;
    void Norm(vec10& y) const;

  public:
    TriaxialODE(const Triaxial& t, vec3 r1, vec3 v1,
                bool extended = true, bool interp = true,
                bool dense = true, real eps = 0);
    TriaxialODE(const Triaxial& t, Angle bet1, Angle omg1, Angle alp1,
                bool extended = true, bool interp = true,
                bool dense = true, real eps = 0);
    TriaxialODE(const Triaxial& t, real bet1, real omg1, real alp1,
                bool extended = true, bool interp = true,
                bool dense = true, real eps = 0);
    bool Position(real s12, vec3& r2, vec3& v2);
    bool Position(real s12, vec3& r2, vec3& v2,
                  real& m12, real& M12, real& M21);
    bool Position(real s12, Angle& bet2, Angle& omg2, Angle& alp2);
    bool Position(real s12, Angle& bet2, Angle& omg2, Angle& alp2,
                  real& m12, real& M12, real& M21);
    void Reset();
    long NSteps() const { return _nsteps; }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALODE_HPP

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

#include <iostream>
#include <array>
#include <vector>
#include <limits>
#include <functional>
#include <GeographicLib/Constants.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"

namespace GeographicLib {

  class GEOGRAPHICLIB_EXPORT TriaxialODE {
  private:
    typedef Math::real real;
    typedef Triaxial::vec3 vec3;
    typedef std::array<real, 6> vec6;
    typedef std::array<real, 10> vec10;
    typedef Angle ang;
    Triaxial _t;
    real _b;
    vec3 _axesn, _axes2n, _r1, _v1;
    // These private versions of Accel and Norm assume normalized ellipsoid
    // with axes = axesn
    vec6 Accel(const vec6& y) const;
    void Norm(vec6& y) const;
    // These private versions of Accel and Norm assume normalized ellipsoid
    // with axes = axesn
    vec10 Accel(const vec10& y) const;
    void Norm(vec10& y) const;
  public:
    TriaxialODE(const Triaxial& t, vec3 r1, vec3 v1);
    TriaxialODE(const Triaxial& t, const Angle& bet1, const Angle& omg1,
                const Angle& alp1);
    TriaxialODE(const Triaxial& t, real bet1, real omg1, real alp1);
    int Position(real s12, vec3& r2, vec3& v2, real eps = 0) const;
    int Position(real s12, vec3& r2, vec3& v2, real& m12, real& M12, real& M21,
                 real eps = 0) const;
    void Position(real ds, long nmin, long nmax,
                  std::vector<vec3>& r2, std::vector<vec3>& v2,
                  real eps = 0) const;
    void Position(real ds, long nmin, long nmax,
                  std::vector<vec3>& r2, std::vector<vec3>& v2,
                  std::vector<real>& m12,
                  std::vector<real>& M12, std::vector<real>& M21,
                  real eps = 0) const;
    int Position(real s12, Angle& bet2, Angle& omg2, Angle& alp2,
                 real eps = 0) const;
    int Position(real s12, Angle& bet2, Angle& omg2, Angle& alp2,
                 real& m12, real& M12, real& M21,
                 real eps = 0) const;
    void Position(real ds, long nmin, long nmax,
                  std::vector<Angle>& bet2, std::vector<Angle>& omg2,
                  std::vector<Angle>& alp2, real eps = 0) const;
    void Position(real ds, long nmin, long nmax,
                  std::vector<Angle>& bet2, std::vector<Angle>& omg2,
                  std::vector<Angle>& alp2,
                  std::vector<real>& m12,
                  std::vector<real>& M12, std::vector<real>& M21,
                  real eps = 0) const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALODE_HPP

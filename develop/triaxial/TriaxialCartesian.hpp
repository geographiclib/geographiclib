/**
 * \file TriaxialCartesian.hpp
 * \brief Header for GeographicLib::TriaxialCartesian class
 *
 * Copyright (c) Charles Karney (2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIALCARTESIAN_HPP)
#define GEOGRAPHICLIB_TRIAXIALCARTESIAN_HPP 1

#include <array>
#include <utility>
#include <functional>
#include "Triaxial.hpp"

namespace GeographicLib {

  class /*GEOGRAPHICLIB_EXPORT*/ TriaxialCartesian {
  public:
    typedef Triaxial::vec3 vec3;
  private:
    typedef Math::real real;
    typedef Angle ang;
    const Triaxial _t;
    const real _b;
    const vec3 _axes, _axes2, _linecc2;

    template<int n>
    void cart2togeneric(vec3 r, Angle& phi, Angle& lam) const;
    template<int n>
    void generictocart2(Angle phi, Angle lam, vec3& r) const;
    real cubic(vec3 r2) const;
    template<int n>
    class funp {
    private:
      // Evaluate
      //   f(p) = sum( (r[0]/(p + l[0]))^n, k = 0..2) - 1
      // and it derivative.
      const real _d;
      const vec3 _r, _l;
    public:
      funp(const vec3& r, const vec3& l)
        : _d(std::numeric_limits<real>::epsilon()/2)
        , _r(r)
        , _l(l)
      {
        static_assert(n >= 1 && n <= 2, "Bad power in funp");
      }
      std::pair<real, real> operator()(real p) const;
    };
    static
    real cartsolve(const std::function<std::pair<real, real>(real)>& f,
                   real p0, real pscale);
  public:
    TriaxialCartesian(const Triaxial& t);
    const Triaxial& t() const { return _t; }
    void cart2toellip(vec3 r, Angle& bet, Angle& omg) const {
      _t.cart2toellip(r, bet, omg);
    }
    void cart2toellip(vec3 r, vec3 v,
                      Angle& bet, Angle& omg, Angle& alp) const {
      _t.cart2toellip(r, v, bet, omg, alp);
    }
    void cart2toellip(Angle bet, Angle omg,
                      vec3 v, Angle& alp) const {
      _t.cart2toellip(bet, omg, v, alp);
    }
    void elliptocart2(Angle bet, Angle omg, vec3& r) const {
      _t.elliptocart2(bet, omg, r);
    }
    void elliptocart2(Angle bet, Angle omg, Angle alp,
                      vec3& r, vec3& v) const {
      _t.elliptocart2(bet, omg, alp, r, v);
    }

    void carttoellip(vec3 r, Angle& bet, Angle& omg, real& H) const;
    void elliptocart(Angle bet, Angle omg, real H, vec3& r) const;

    void cart2togeod(vec3 r, Angle& phi, Angle& lam) const;
    void geodtocart2(Angle phi, Angle lam, vec3& r) const;
    void cart2toparam(vec3 r, Angle& phip, Angle& lamp) const;
    void paramtocart2(Angle phip, Angle lamp, vec3& r) const;
    void cart2togeocen(vec3 r, Angle& phipp, Angle& lampp) const;
    void geocentocart2(Angle phipp, Angle lampp, vec3& r) const;

    void cart2tocart(vec3 r2, real h, vec3& r) const;
    void carttocart2(vec3 r, vec3& r2, real& h) const;
    void carttogeod(vec3 r, Angle& phi, Angle& lam, real& h) const;
    void geodtocart(Angle phi, Angle lam, real h, vec3& r) const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALCARTESIAN_HPP

/**
 * \file Cartesian3.hpp
 * \brief Header for GeographicLib::Triaxial::Cartesian3 class
 *
 * Copyright (c) Charles Karney (2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CARTESIAN3_HPP)
#define GEOGRAPHICLIB_CARTESIAN3_HPP 1

#include <utility>
#include <functional>
#include <random>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs random
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {
  namespace Triaxial {

  class GEOGRAPHICLIB_EXPORT Cartesian3 {
  public:
    using vec3 = Ellipsoid3::vec3;
  private:
    using real = Math::real;
#if GEOGRAPHICLIB_PRECISION > 3
    // <random> only supports "standard" floating point types
    using random_prec = Math::extended;
#else
    using random_prec = Math::real;
#endif
    using ang = Angle;
    static const int maxit_ = 20;
    const Ellipsoid3 _t;
    const vec3 _axes, _axes2, _linecc2;
    // mutable because using these objects in a non-const operation
    mutable std::normal_distribution<random_prec> _norm;
    mutable std::uniform_real_distribution<random_prec> _uni;

    template<int n>
    void cart2togeneric(vec3 r, ang& phi, ang& lam) const;
    template<int n>
    void generictocart2(ang phi, ang lam, vec3& r) const;
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
    // real a() const { return t().a(); } // not needed
    real b() const { return t().b(); }
    real c() const { return t().c(); }
  public:
    Cartesian3(const Ellipsoid3& t);
    Cartesian3(real a, real b, real c);
    Cartesian3(real b, real e2, real k2, real kp2);
    const Ellipsoid3& t() const { return _t; }
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

    template <class G> void cart2rand(G& g, vec3& r) const;
    template <class G> void cart2rand(G& g, vec3& r, vec3& v) const;
  };

  template<class G> inline void Cartesian3::cart2rand(G& g, vec3& r) const {
    // This uses the simple rejection technique given by Marples and Williams,
    // Num. Alg. (2023), Algorithm 1 based on the general method of Williamson,
    // Phys. Med. Biol. (1987).
    using std::isfinite;
    while (true) {
      while (true) {
        // guaranteed evaluated left to right
        r = {real(_norm(g)), real(_norm(g)), real(_norm(g))};
        Ellipsoid3::normvec(r); // But catch rare cases where |r| = 0
        if (isfinite(r[0])) break;
      }
      r[0] *= _axes[0]; r[1] *= _axes[1]; r[2] *= _axes[2];
      vec3 up{ r[0] / _axes2[0],  r[1] / _axes2[1],  r[2] / _axes2[2] };
      real q = c() * Math::hypot3(up[0], up[1], up[2]);
      if (real(_uni(g)) < q) break;
    }
  }
  template<class G> inline void Cartesian3::cart2rand(G& g, vec3& r, vec3& v)
  const {
    using std::isfinite;
    cart2rand<G>(g, r);
    while (true) {
      // guaranteed evaluated left to right
      v = {real(_norm(g)), real(_norm(g)), real(_norm(g))};
      vec3 up{ r[0] / _axes2[0],  r[1] / _axes2[1],  r[2] / _axes2[2] };
      real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]), // |up|^2
        // (up . v) / |up|^2
        uv = (v[0] * up[0] + v[1] * up[1] + v[2] * up[2])/u2;
      // v - up * (up . v) / |up|^2
      v[0] -= uv * up[0]; v[1] -= uv * up[1]; v[2] -= uv * up[2];
      Ellipsoid3::normvec(v);   // But catch rare cases where |v| = 0
      if (isfinite(v[0])) break;
    }
  }

  } // namespace Triaxial
} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_CARTESIAN3_HPP

/**
 * \file Triaxial.hpp
 * \brief Header for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIAL_HPP)
#define GEOGRAPHICLIB_TRIAXIAL_HPP 1

#include <iostream>
#include <array>
#include <vector>
#include <limits>
#include <functional>
#include <utility>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include "Angle.hpp"
#include "Trigfun.hpp"

namespace GeographicLib {
  class TriaxialLine;

  class GEOGRAPHICLIB_EXPORT Triaxial {
  public:
    typedef std::array<Math::real, 3> vec3;
  private:
    friend class TriaxialLine;
    typedef Math::real real;
    typedef Angle ang;
    static void normvec(vec3& r) {
      real h = Math::hypot3(r[0], r[1], r[2]);
      r[0] /= h; r[1] /= h; r[2] /= h;
    }
    static void Flip(Angle& bet, Angle& omg, Angle& alp) {
      bet.reflect(false, true);
      omg.reflect(true);
      alp.reflect(true, true);
    }

    static real HybridA(const Triaxial& t,
                        Angle bet1, Angle omg1,
                        Angle alp1,
                        Angle bet2, Angle omg2);
    static Angle findroot(const std::function<Math::real(const Angle&)>& f,
                          Angle xa,  Angle xb,
                          Math::real fa, Math::real fb,
                          int* countn = nullptr, int* countb = nullptr);
    real _a, _b, _c;            // semi-axes
    vec3 _axes;
    real _e2, _k2, _kp2, _k, _kp;
    bool _umbalt,               // how coordinates wrap with umbilical lines
      _debug;                   // print out diagnostics
    static real BigValue() {
      using std::log;
      static real bigval = -3*log(std::numeric_limits<real>::epsilon());
      return bigval;
    }
  public:
    Triaxial();
    Triaxial(real a, real b, real c);
    Triaxial(real b, real e2, real k2, real kp2);
    void Norm(vec3& r) const;
    void Norm(vec3& r, vec3& v) const;
    TriaxialLine Inverse(Angle bet1, Angle omg1, Angle bet2, Angle omg2,
                         Angle& alp1, Angle& alp2, real& s12) const;
    TriaxialLine Inverse(real bet1, real omg1, real bet2, real omg2,
                         real& alp1, real& alp2, real& s12) const;
    TriaxialLine Line(Angle bet1, Angle omg1, Angle alp1) const;
    TriaxialLine Line(real bet1, real omg1, real alp1) const;
    TriaxialLine Direct(Angle bet1, Angle omg1, Angle alp1, real s12,
                        Angle& bet2, Angle& omg2, Angle& alp2) const;
    TriaxialLine Direct(real bet1, real omg1, real alp1, real s12,
                        real& bet2, real& omg2, real& alp2) const;
    real a() const { return _a; }
    real b() const { return _b; }
    real c() const { return _c; }
    real e2() const { return _e2; }
    real k2() const { return _k2; }
    real kp2() const { return _kp2; }
    const vec3& axes() const { return _axes; }
    bool umbalt() const { return _umbalt; }
    void umbalt(bool numbalt) { _umbalt = numbalt; }
    void debug(bool ndebug) { _debug = ndebug; }
    static bool AngNorm(Angle& bet, Angle& omg, Angle& alp,
                        bool alt = false) {
      using std::signbit;
      // If !alt, put bet in [-pi/2,pi/2]
      // If  alt, put omg in [0, pi]
      bool flip = alt ? signbit(omg.s()) : signbit(bet.c());
      if (flip)
        Flip(bet, omg, alp);
      if (0) {
        if (bet.c() == 0 && bet.s() * alp.c() > 0)
          alp.reflect(true, true);
        if (bet.c() == 0 && alp.c() == 0)
          alp.reflect(alp.s() * bet.s() > 0); // alp.s() = -bet.s();
      }
      return flip;
    }
    static bool AngNorm(Angle& bet, Angle& omg,
                        bool alt = false) {
      using std::signbit;
      // If !alt, put bet in [-pi/2,pi/2]
      // If  alt, put omg in [0, pi]
      bool flip = alt ? signbit(omg.s()) : signbit(bet.c());
      if (flip) {
        ang alp;
        Flip(bet, omg, alp);
      }
      return flip;
    }
    void cart2toellip(vec3 r, Angle& bet, Angle& omg) const;
    void cart2toellip(vec3 r, vec3 v,
                      Angle& bet, Angle& omg, Angle& alp) const;
    void cart2toellip(Angle bet, Angle omg,
                      vec3 v, Angle& alp) const;
    void elliptocart2(Angle bet, Angle omg, vec3& r) const;
    void elliptocart2(Angle bet, Angle omg,
                      Angle alp,
                      vec3& r, vec3& v) const;
    real EuclideanInverse(Angle bet1, Angle omg1, Angle bet2, Angle omg2,
                          Angle& alp1, Angle& alp2) const;
    real EuclideanInverse(vec3 r1, vec3 r2,
                          vec3& v1, vec3& v2) const;
    std::pair<real, real> EuclideanDiff(Angle bet1, Angle omg1, Angle alp1,
                                        Angle bet2, Angle omg2, Angle alp2)
      const;
    std::pair<real, real> EuclideanDiff(vec3 r1, vec3 v1,
                                        vec3 r2, vec3 v2)
      const;
    class gamblk {
    public:
      // gamma = (k * cbet * salp)^2 - (kp * somg * calp)^2
      //       = k2*cb2*sa2 - kp2*so2*ca2
      // Need accurate expressions for
      //   k2  - gamma = k2*(sb2+ca2*cb2) + kp2*so2*ca2
      //   kp2 + gamma = k2*cb2*sa2 + kp2*(co2+sa2*so2)
      // If gamma is given, eval new alp given new bet and new omg
      // gamma < 0
      //   ca2 = (k2*cb2-gamma) / (k2*cb2+kp2*so2)
      //   sa2 = (kp2+gamma - kp2*co2) / (k2*cb2+kp2*so2)
      // gamma > 0
      //   ca2 = (k2-gamma - k2*sb2) / (k2*cb2+kp2*so2)
      //   sa2 = (kp2*so2+gamma) / (k2*cb2+kp2*so2)
      // gamma > 0
      //   k2*sb2 = spsi2 * (k2-gamma)
      //   (k2-gamma - k2*sb2) = (k2-gamma)*(1-spsi2) = (k2-gamma)*cpsi2
      //   spsi2 = k2*sb2/(k2-gamma)
      //   cpsi2 = (k2*cb2-gamma)/(k2-gamma)
      real gamma,
      // [nu, nup]
      //   = [sqrt(gam)/k, sqrt(1 - gam/k2)] for gam > 0,
      //   = [sqrt(-gam)/kp, sqrt(1 + gam/kp2)] for gam < 0
      //   unused for gam == 0
        nu, nup;
    };
    gamblk gamma(Angle bet, Angle omg, Angle alp)
      const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIAL_HPP

/**
 * \file Triaxial.hpp
 * \brief Header for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIAL_HPP)
#define GEOGRAPHICLIB_TRIAXIAL_HPP 1

#include <iostream>
#include <array>
#include <functional>
#include <memory>
#include <GeographicLib/Constants.hpp>
#include "Angle.hpp"
#include "Trigfun.hpp"

namespace GeographicLib {
  class TriaxialLine;

  class /*GEOGRAPHICLIB_EXPORT*/ Triaxial {
  public:
    typedef std::array<Math::real, 3> vec3;
  private:
    friend class TriaxialLine;
    friend class TriaxialCartesian;  // For access to cart2toellipint
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

    // Run geodesic from bet1, omg1, alp1, find its first intersection with bet
    // = bet2a and return omg2a - omg2b
    real HybridA(Angle bet1, Angle omg1, Angle alp1,
                 Angle bet2a, Angle omg2b, bool betp) const;
    static Angle findroot(const std::function<real(const Angle&)>& f,
                          Angle xa,  Angle xb,
                          real fa, real fb,
                          int* countn = nullptr, int* countb = nullptr);
    real _a, _b, _c;            // semi-axes
    real _e2, _k2, _kp2, _k, _kp;
    bool _oblate, _prolate, _biaxial;
    bool _umbalt,               // how coordinates wrap with umbilical lines
      _biaxp,                   // special treatment for biaxial non-meridional
      _debug,                   // print out diagnostics
      _hybridalt,               // favor hybrid solution in terms of omg
      _swapomg;                 // allow swapping of omega{1,2}
    // If k'^2 < ellipthresh transform phi -> F(phi, k^2)
    real _ellipthresh;
    mutable std::shared_ptr<TriaxialLine> _umbline;
    static real BigValue() {
      using std::log;
      static real bigval = -3*log(std::numeric_limits<real>::epsilon());
      return bigval;
    }
    void cart2toellipint(vec3 r, Angle& bet, Angle& omg, vec3 axes) const;
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
    real k() const { return _k; }
    real kp() const { return _kp; }
    bool umbalt() const { return _umbalt; }
    void umbalt(bool numbalt) { if (_k2 > 0 && _kp2 > 0) _umbalt = numbalt; }
    bool biaxp() const { return _biaxp; }
    void biaxp(bool biaxp) { _biaxp = biaxp; }
    bool debug() const { return _debug; }
    void debug(bool debug) { _debug = debug; }
    bool hybridalt() const { return _hybridalt; }
    void hybridalt(bool hybridalt) { _hybridalt = hybridalt; }
    bool swapomg() const { return _swapomg; }
    void swapomg(bool swapomg) { _swapomg = swapomg; }
    real ellipthresh() const { return _ellipthresh; }
    void ellipthresh(real ellipthresh) { _ellipthresh = ellipthresh; }
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
    void elliptocart2(Angle bet, Angle omg, Angle alp,
                      vec3& r, vec3& v) const;
    class gamblk {
    public:
      bool transpolar;
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
      //   = [sqrt(gam)/k, sqrt(1 - gam/k2)] for !signbit(gam)
      //   = [sqrt(-gam)/kp, sqrt(1 + gam/kp2)] for signbit(gam)
      //   unused for umbilics
        nu, nup,
        gammax, kx2, kxp2, kx, kxp;
      // Default values for gamma = +/-0
      gamblk() {}
      gamblk(const Triaxial& t, bool neg = false);
      gamblk(const Triaxial& t, Angle bet, Angle omg, Angle alp);
      //       gamblk(real gammax, real nux, real nupx)
      //        : gamma(gammax), nu(nux), nup(nupx) {}
    };
    gamblk gamma(Angle bet, Angle omg, Angle alp)
      const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIAL_HPP

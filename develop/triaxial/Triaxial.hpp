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
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/AuxAngle.hpp>
#include "Trigfun.hpp"

namespace GeographicLib {
  class GEOGRAPHICLIB_EXPORT Triaxial {
  public:
    typedef std::array<Math::real, 3> vec3;
  private:
    typedef Math::real real;
    typedef std::array<real, 6> vec6;
    // These private versions of Accel and Norm assume normalized ellipsoid
    // with axes = axesn
    vec6 Accel(const vec6& y) const;
    void Norm(vec6& y) const;
    static real hypot3(real x, real y, real z) {
#if __cplusplus < 201703L || GEOGRAPHICLIB_PRECISION == 4
      using std::sqrt;
      return sqrt(x*x + y*y + z*z);
#else
      using std::hypot;
      return hypot(x, y, z);
#endif
    }
    static void normvec(vec3& r) {
      real h = hypot3(r[0], r[1], r[2]);
      r[0] /= h; r[1] /= h; r[2] /= h;
    }
  public:
    real a, b, c;               // semi-axes
    vec3 axes, axesn, axes2n;
    real e2, k2, kp2, k, kp;
    Triaxial(real a, real b, real c);
    void Norm(vec3& r) const;
    void Norm(vec3& r, vec3& v) const;
    int Direct(const vec3& r1, const vec3& v1, real s12, vec3& r2, vec3& v2,
               real eps = 0) const;
    void Direct(const vec3& r1, const vec3& v1, real ds,
                long nmin, long nmax,
                std::vector<vec3>& r2, std::vector<vec3>& v2,
                real eps = 0) const;
    static real EllipticThresh() {
      static real thresh = 1/real(8);
      return thresh;
    }
    real gamma(AuxAngle& bet, AuxAngle& omg, AuxAngle& alp) const {
      // gamma = (k * cbet * salp)^2 - (kp * somg * calp)^2
      //       = k2*cb2*sa2 - kp2*so2*ca2
      // Maybe need accurate expressions for
      //   k2  - gamma = k2*(sb2+ca2*cb2) + kp2*so2*ca2
      //   kp2 + gamma = k2*cb2*sa2 + kp2*(co2+sa2*so2)
      // If gamma is given eval new alp given new bet and new omg
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
      using std::fabs; using std::copysign;
      bet.normalize(); omg.normalize(); alp.normalize();
      real a = k * bet.x() * alp.y(), b = kp * omg.y() * alp.x(),
        gam = (a - b) * (a + b);
      // Factor of 8 allows umbilical geodesics with alp in [-180,180] degrees
      // to return gam = 0.  (If in radians the factor could be reduced to 4.)
      if (fabs(gam) < 8 * std::numeric_limits<real>::epsilon() && gam != 0) {
        gam = 0 * gam;
        // fix alp to be consistent with gamma = 0
        alp = AuxAngle(copysign( kp * omg.y(), alp.y() ),
                       copysign( k  * bet.x(), alp.x() )).normalized();
      } else if (fabs(k2 - gam) <= std::numeric_limits<real>::epsilon()) {
        gam = k2;
        // This requires cbet = +/- 1 and salp = +/- 1
        bet = AuxAngle(copysign(real(0), bet.y()), copysign(real(1), bet.x()));
        alp = AuxAngle(copysign(real(1), alp.y()), copysign(real(0), alp.x()));
      } else if (fabs(kp2 + gam) <= std::numeric_limits<real>::epsilon()) {
        gam = -kp2;
        // This requires somg = +/- 1 and calp = +/- 1
        omg = AuxAngle(copysign(real(1), omg.y()), copysign(real(0), omg.x()));
        alp = AuxAngle(copysign(real(0), alp.y()), copysign(real(1), alp.x()));
      }
      return gam;
    }
    static void AngNorm(AuxAngle& bet, AuxAngle& omg, AuxAngle& alp,
                        bool alt = false) {
      using std::signbit;
      // If !alt, put bet in [-pi/2,pi/2]
      // If  alt, put omg in [0, pi]
      if (alt ? signbit(omg.y()) : signbit(bet.x())) {
        bet.x() *= -1;
        alp.y() *= -1; alp.x() *= -1;
        omg.y() *= -1;
      }
    }
    void cart2toellip(const vec3& r, AuxAngle& bet, AuxAngle& omg) const;
    void cart2toellip(const vec3& r, const vec3& v,
                      AuxAngle& bet, AuxAngle& omg, AuxAngle& alp) const;
    void elliptocart2(const AuxAngle& bet, const AuxAngle& omg, vec3& r) const;
    void elliptocart2(const AuxAngle& bet, const AuxAngle& omg,
                      const AuxAngle& alp,
                      vec3& r, vec3& v) const;
  };

  class GEOGRAPHICLIB_EXPORT geod_fun {
  private:
    typedef Math::real real;
    real _kap, _kapp, _eps, _mu;
    bool _tx;
    EllipticFunction _ell;
    TrigfunExt _fun;
    real _tol;
    int _nmax;
    Trigfun _chiinv;
    int _countn, _countb;
    real _max;
    bool _invp;
    // mu > 0
    static real fphip(real c, real kap, real kapp, real eps, real mu) {
      using std::sqrt;
      real c2 = kap * Math::sq(c);
      return sqrt( (1 - eps * c2) / (( kapp + c2) * (c2 + mu)) );
    }
    static real fup(real cn, real kap, real kapp, real eps, real mu) {
      using std::sqrt;
      real c2 = kap * Math::sq(cn);
      return sqrt( (1 - eps * c2) / (( kapp + c2) * (kap + mu)) );
    }
    // mu == 0
    static real dfp(real c, real kap, real kapp, real eps) {
      // function dfp = dfpf(phi, kappa, epsilon)
      // return derivative of Delta f
      using std::cos; using std::sqrt;
      // s = sqrt(1 - kap * sin(phi)^2)
      real c2 = kap * Math::sq(c), s = sqrt(kapp + c2);
      return (1 + eps*kapp) * kap * c / (s * (sqrt(kapp * (1 - eps*c2)) + s));
    }
    static real dfvp(real cn, real dn, real kap, real kapp, real eps) {
      // function dfvp = dfvpf(v, kap, eps)
      // return derivative of Delta f_v
      using std::sqrt;
      return (1 + eps*kapp) * kap *
        (cn / (sqrt(kapp * (1 - (eps*kap) * Math::sq(cn))) + dn));
    }
    // mu < 0
    static real fpsip(real s, real c, real kap, real kapp, real eps, real mu) {
      using std::sqrt;
      real c2 = kap * Math::sq(c) - mu * Math::sq(s);
      return sqrt( (1 - eps * c2) / (( kapp + c2) * c2) ) ;
    }
    static real fvp(real dn, real kap, real kapp, real eps, real /* mu */) {
      using std::sqrt;
      real c2 = kap * Math::sq(dn);
      return sqrt( (1 - eps * c2) / ((kapp + c2) * kap) );
    }
    static real lam(real x) {
      using std::tan; using std::asinh; using std::fabs;
      // A consistent large value for x near pi/2.  Also deals with the issue
      // that tan(pi/2) may be negative, e.g., for long doubles.
      // static real bigval = 10 - log(std::numeric_limits<real>::epsilon());
      static real bigval = Math::infinity();
      return fabs(x) < Math::pi()/2 ? asinh(tan(x)) :
        (x < 0 ? -bigval : bigval);
    }
    static real gd(real x) {
      using std::atan; using std::sinh;
      return atan(sinh(x));
    }
    real root(real z, real x0, int* countn, int* countb) const;
  public:
    geod_fun() {}
    geod_fun(real kap, real kapp, real eps, real mu, real epspow = 1,
             real nmaxmult = 0);
    geod_fun(real kap, real kapp, real eps, real mu, bool tx,
             real epsow, real nmaxmult);
    real operator()(real u) const {
      if (_mu == 0) {
        real phi = gd(u);
        return u - _fun(_tx ? _ell.F(phi) : phi);
      } else
        return _fun(u);
    }
    real deriv(real u) const {
      using std::cosh; using std::sqrt;
      if (_mu == 0) {
        real phi = gd(u), sch = 1/cosh(u);
        return 1 - _fun.deriv(_tx ? _ell.F(phi) : phi) * sch /
          ( _tx ? sqrt(_ell.kp2() + _ell.k2() * Math::sq(sch)) : 1);
      } else
        return _fun.deriv(u);
    }

    // Approximate inverse using _finv
    real inv0(real z) const;
    // Accurate inverse by direct Newton (not using _finv)
    real inv1(real z, int* countn = nullptr, int* countb = nullptr) const;
    // Accurate inverse correcting result from _finv
    real inv(real z, int* countn = nullptr, int* countb = nullptr) const;
    void ComputeInverse();
    real fwd(real phi) const {
      return _mu == 0 ? lam(phi) : (_tx ? _ell.F(phi) : phi);
    }
    real rev(real u) const {
      return _mu == 0 ? gd(u) : (_tx ? _ell.am(u) : u);
    }
    int NCoeffs() const { return _fun.NCoeffs(); }
    int NCoeffsInv() const {
      return _mu == 0 ? _chiinv.NCoeffs() : _fun.NCoeffsInv();
    }
    std::pair<int, int> InvCounts() const {
      return _mu == 0 ? std::pair<int, int>(_countn, _countb) :
        _fun.InvCounts();
    }
    bool txp() const { return _tx; }
    real HalfPeriod() const {
      return _mu == 0 ? Math::infinity() : (_tx ? _ell.K() : Math::pi()/2);
    }
    real Slope() const {
      return _mu == 0 ? 1 : _fun.Slope();
    }
    real Max() const {
      return _max;
    }
  };

  class GEOGRAPHICLIB_EXPORT dist_fun {
  private:
    typedef Math::real real;
    real _kap, _kapp, _eps, _mu;
    bool _tx;
    EllipticFunction _ell;
    TrigfunExt _fun;
    real _max;
    // _mu > 0
    static real gphip(real c, real kap, real kapp, real eps, real mu) {
      using std::sqrt;
      real c2 = kap * Math::sq(c);
      return sqrt((c2 + mu) * (1 - eps * c2) / ( kapp + c2) );
    }
    static real gfphip(real c, real kap, real mu) {
      real c2 = kap * Math::sq(c);
      return c2 + mu;
    }
    static real gup(real cn, real dn, real kap, real kapp, real eps, real mu) {
      using std::sqrt;
      real c2 = kap * Math::sq(cn);
      return sqrt( (kap + mu) * (1 - eps * c2) / ( kapp + c2) ) * Math::sq(dn);
    }
    static real gfup(real cn, real kap, real mu) {
      real c2 = kap * Math::sq(cn);
      return c2 + mu;           // or (kap + mu) * Math::sq(dn);
    }
    // _mu == 0
    static real g0p(real c, real kap, real kapp, real eps) {
      using std::sqrt;
      real c2 = kap * Math::sq(c);
      return sqrt( kap * (1 - eps * c2) / ( kapp + c2) ) * c;
    }
    /*
    static real gf0p(real c, real kap, real kapp) {
      using std::sqrt;
      real c2 = sqrt(kap * kapp) * kap * Math::sq(c);
      return c2;
    }
    */
    static real gf0up(real u, real kap, real kapp) {
      using std::cosh; using std::sqrt;
      // Adjust by sqrt(kap * kappp) to account of factor removed from f
      // functions.
      real c2 = sqrt(kap / kapp) / Math::sq(cosh(u));
      return c2;
    }
    static real g0vp(real cn, real kap, real /*kapp*/, real eps) {
      using std::sqrt;
      real c2 = kap * Math::sq(cn);
      return sqrt( kap * (1 - eps * c2) ) * cn;
    }
    // _mu < 0
    static real gpsip(real s, real c, real kap, real kapp, real eps, real mu) {
      using std::sqrt;
      // kap * cos(phi)^2
      real c2 = kap * Math::sq(c) - mu * Math::sq(s);
      return (kap + mu) *
        sqrt( (1 - eps * c2) / (( kapp + c2) * c2) ) * Math::sq(c);
    }
    static real gfpsip(real c, real kap, real mu) {
      using std::sqrt;
      return (kap + mu) * Math::sq(c);
    }
    static real gvp(real cn, real dn, real kap, real kapp, real eps, real mu) {
      using std::sqrt;
      real c2 = kap * Math::sq(dn);
      return (kap + mu) *
        sqrt( (1 - eps * c2) / (kap * (kapp + c2))) * Math::sq(cn);
    }
    static real gfvp(real cn, real kap, real mu) {
      // alternatively mu + kap * Math::sq(dn)
      return (kap + mu) * Math::sq(cn);
    }
    static real lam(real x) {
      using std::tan; using std::asinh; using std::fabs;
      // A consistent large value for x near pi/2.  Also deals with the issue
      // that tan(pi/2) may be negative, e.g., for long doubles.
      // static real bigval = 10 - log(std::numeric_limits<real>::epsilon());
      static real bigval = Math::infinity();
      return fabs(x) < Math::pi()/2 ? asinh(tan(x)) :
        (x < 0 ? -bigval : bigval);
    }
    static real gd(real x) {
      using std::atan; using std::sinh;
      return atan(sinh(x));
    }
  public:
    dist_fun() {}
    dist_fun(real kap, real kapp, real eps, real mu);
    dist_fun(real kap, real kapp, real eps, real mu, bool tx);
    real operator()(real u) const {
      if (_mu == 0) {
        real phi = gd(u);
        return _fun(_tx ? _ell.F(phi) : phi);
      } else
        return _fun(u);
    }
    real deriv(real u) const {
      using std::cosh; using std::sqrt;
      if (_mu == 0) {
        real phi = gd(u), sch = 1/cosh(u);
        return _fun.deriv(_tx ? _ell.F(phi) : phi) * sch /
          ( _tx ? sqrt(_ell.kp2() + _ell.k2() * Math::sq(sch)) : 1);
      } else
        return _fun.deriv(u);
    }
    real gfderiv(real u) const;
    // Don't need these
    // real inv(real y) const { return _fun.inv(y); }
    // real inv1(real y) const { return _fun.inv1(y); }
    // Use geod_fun versions of these
    // real fwd(real phi) const {
    //   return _mu == 0 ? lam(phi) : (_tx ? _ell.F(phi) : phi);
    // }
    // real rev(real u) const {
    //   return _mu == 0 ? gd(u) : (_tx ? _ell.am(u) : u);
    // }
    int NCoeffs() const { return _fun.NCoeffs(); }
    //    int NCoeffsInv() const { return _fun.NCoeffsInv(); }
    //    std::pair<int, int> InvCounts() const { return _fun.InvCounts(); }
    bool txp() const { return _tx; }
    real HalfPeriod() const {
      return _mu == 0 ? Math::infinity() : (_tx ? _ell.K() : Math::pi()/2);
    }
    real Slope() const {
      return _mu == 0 ? 1 : _fun.Slope();
    }
    real Max() const {
      return _max;
    }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIAL_HPP

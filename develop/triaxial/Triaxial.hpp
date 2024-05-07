/**
 * \file Triaxial.hpp
 * \brief Header for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2022-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIAL_HPP)
#define GEOGRAPHICLIB_TRIAXIAL_HPP 1

#include <iostream>
#include <array>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/EllipticFunction.hpp>
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
  public:
    real a, b, c;               // semi-axes
    vec3 axes, axesn, axes2n;
    real e2, k2, kp2, k, kp;
    Triaxial(real a, real b, real c);
    void Norm(vec3& r) const;
    void Norm(vec3& r, vec3& v) const;
    int Direct(const vec3& r1, const vec3& v1, real s12, vec3& r2, vec3& v2,
               real eps = 1e-10) const;
  };

  class GEOGRAPHICLIB_EXPORT geod_fun {
  private:
    typedef Math::real real;
    real _kap, _kapp, _eps, _mu;
    bool _tx;
    static real fphipf(real phi, real kap, real kapp, real eps, real mu) {
      using std::cos; using std::sqrt;
      real c2 = kap * Math::sq(cos(phi));
      return sqrt( (1 - eps * c2) / (( kapp + c2) * (c2 + mu)) );
    }
    static real fupf(real u, real kap, real kapp, real eps, real mu,
                     EllipticFunction ell) {
      using std::sqrt;
      real sn, cn, dn;
      ell.am(u, sn, cn, dn);
      real c2 = kap * Math::sq(cn);
      return sqrt( (1 - eps * c2) / (( kapp + c2) * (kap + mu)) );
    }
    static real fpsipf(real psi, real kap, real kapp, real eps, real mu) {
      using std::cos; using std::sin; using std::sqrt;
      real c2 = kap * Math::sq(cos(psi)) - mu * Math::sq(sin(psi));
      return sqrt( (1 - eps * c2) / (( kapp + c2) * c2) ) ;
    }
    static real fvpf(real v, real kap, real kapp, real eps, real /* mu */,
                     EllipticFunction ell) {
      using std::sqrt;
      real sn, cn, dn;
      ell.am(v, sn, cn, dn);
      real c2 = kap * Math::sq(dn);
      return sqrt( (1 - eps * c2) / (kapp + c2) * kap );
    }
  public:
    EllipticFunction ell;
    TrigfunExt fun;
    geod_fun(real kap, real kapp, real eps, real mu);
    geod_fun(real kap, real kapp, real eps, real mu, bool tx);
    real val(real u) const { return fun(u); }
    real inv(real y) const { return fun.inv(y); }
    real inv1(real y) const { return fun.inv1(y); }
    real tx(real phi) const { return _tx ? ell.F(phi) : phi; }
    real txinv(real u) const { return _tx ? ell.am(u) : u; }
    int NCoeffs() const { return fun.NCoeffs(); }
    int NCoeffsInv() const { return fun.NCoeffsInv(); }
    std::pair<int, int> InvCounts() const { return fun.InvCounts(); }
    bool txp() const { return _tx; }
  };

  class GEOGRAPHICLIB_EXPORT geodu_fun {
    // geod_fun for umbilic geodesics mu = 0
  private:
    typedef Math::real real;
    real _kap, _kapp, _eps;
    bool _tx;
    static real lam(real x) {
      using std::tan; using std::asinh; using std::fabs;
      // A consistent large value for x near pi/2.  Also deals with the issue
      // that tan(pi/2) may be negative, e.g., for long doubles.
      static real bigval = 10 - log(std::numeric_limits<real>::epsilon());
      return fabs(x) < Math::pi()/2 ? asinh(tan(x)) :
        (x < 0 ? -bigval : bigval);
    }
    static real gd(real x) {
      using std::atan; using std::sinh;
      return atan(sinh(x));
    }
    static real dfpf(real phi, real kap, real kapp, real eps) {
      // function dfp = dfpf(phi, kappa, epsilon)
      // return derivative of Delta f
      using std::cos; using std::sqrt;
      // s = sqrt(1 - kap * sin(phi)^2)
      real c = cos(phi), c2 = kap * Math::sq(c), s = sqrt(kapp + c2);
      return (1 + eps*kapp) * kap * c / (s * (sqrt(kapp * (1 - eps*c2)) + s));
    }
    static real dfvpf(real v, real kap, real kapp, real eps,
                      EllipticFunction ell) {
      // function dfvp = dfvpf(v, kap, eps)
      // return derivative of Delta f_v
      using std::sqrt;
      real sn, cn, dn;
      ell.am(v, sn, cn, dn);
      return (1 + eps*kapp) * kap *
        (cn / (sqrt(kapp * (1 - (eps*kap) * Math::sq(cn))) + dn));
    }
  public:
    // v = F(phi, kap), phi = am(v, kap)
    // u = lam(phi), phi = gd(u)
    // Delta_phi f(phi) is even sym periodic period 2*pi
    // deriv = dfpf
    // Delta_v f(v) is even sym periodic period 4*K
    // deriv = dfvpf
    // ignoring factor of sqrt(kap*kapp)
    // f(u) = u - Delta_phi f(gd(u))
    //      = u - Delta_v f(F(gd(u), kap))

    EllipticFunction ell;
    TrigfunExt fun;
    geodu_fun(real kap, real kapp, real eps);
    geodu_fun(real kap, real kapp, real eps, bool tx);
    real val(real u) const { return fun(u); }
    real inv(real y) const { return fun.inv(y); }
    real inv1(real y) const { return fun.inv1(y); }
    real tx(real phi) const { return _tx ? ell.F(phi) : phi; }
    real txinv(real u) const { return _tx ? ell.am(u) : u; }
    int NCoeffs() const { return fun.NCoeffs(); }
    int NCoeffsInv() const { return fun.NCoeffsInv(); }
    std::pair<int, int> InvCounts() const { return fun.InvCounts(); }
    bool txp() const { return _tx; }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIAL_HPP

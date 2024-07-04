/**
 * \file TriaxialLine.hpp
 * \brief Header for GeographicLib::TriaxialLin class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIALLINE_HPP)
#define GEOGRAPHICLIB_TRIAXIALLINE_HPP 1

#include <utility>
#include <GeographicLib/Constants.hpp>
#include "Angle.hpp"
#include "Trigfun.hpp"
#include "Triaxial.hpp"

namespace GeographicLib {


  class GEOGRAPHICLIB_EXPORT ffun {
  private:
    typedef Math::real real;
    typedef Angle ang;
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
      using std::sqrt;
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
    real root(real z, real x0, int* countn, int* countb) const;
  public:
    ffun() {}
    ffun(real kap, real kapp, real eps, real mu, real epspow = 1,
             real nmaxmult = 0);
    ffun(real kap, real kapp, real eps, real mu, bool tx,
             real epsow, real nmaxmult);
    static real lam(real x) {
      using std::tan; using std::asinh; using std::fabs;
      // A consistent large value for x near pi/2.  Also deals with the issue
      // that tan(pi/2) may be negative, e.g., for long doubles.
      return fabs(x) < Math::pi()/2 ? asinh(tan(x)) :
        (x < 0 ? -1 : 1) * Triaxial::BigValue();
    }
    static real gd(real x) {
      using std::atan; using std::sinh;
      return atan(sinh(x));
    }
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

  class GEOGRAPHICLIB_EXPORT gfun {
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
  public:
    gfun() {}
    gfun(real kap, real kapp, real eps, real mu);
    gfun(real kap, real kapp, real eps, real mu, bool tx);
    real operator()(real u) const {
      if (_mu == 0) {
        real phi = ffun::gd(u);
        return _fun(_tx ? _ell.F(phi) : phi);
      } else
        return _fun(u);
    }
    real deriv(real u) const {
      using std::cosh; using std::sqrt;
      if (_mu == 0) {
        real phi = ffun::gd(u), sch = 1/cosh(u);
        return _fun.deriv(_tx ? _ell.F(phi) : phi) * sch /
          ( _tx ? sqrt(_ell.kp2() + _ell.k2() * Math::sq(sch)) : 1);
      } else
        return _fun.deriv(u);
    }
    real gfderiv(real u) const;
    // Don't need these
    // real inv(real y) const { return _fun.inv(y); }
    // real inv1(real y) const { return _fun.inv1(y); }
    // Use ffun versions of these
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

  class GEOGRAPHICLIB_EXPORT fline {
  private:
    friend class TriaxialLine;
    typedef Math::real real;
    typedef Angle ang;
    Triaxial _t;
    Triaxial::gamblk _gm;
    ffun _fbet, _fomg;
  public:
    real df, deltashift;
    class fics {
      // bundle of data setting the initial conditions for a geodesic
    public:
      Angle bet1, omg1, alp1,   // starting point (omg displaced by -pi/2)
        psi1,                   // nonumbilic angle psi
      // Angles about which quantities oscillate
      // circumpolar:
      //   omg0 not used, bet0 = cardinal(even), alp0 = cardinal(odd)
      // transpolar:
      //   bet0 not used, omg0 = cardinal(even), alp0 = cardinal(even)
      // umbilical
      //   bet0, omg0 = cardinal(even) = middle of starting segment
      //   !umbalt: alp0 = cardinal(odd)
      //   umbalt: alp0 = cardinal(even)
        bet0, omg0, alp0;
      real u0, v0, delta;              // starting point geodesic
      int nN, eE;                  // Northgoing / eastgoing
      fics();
      fics(const fline& f,
           const Angle& bet1, const Angle& omg1, const Angle& alp1);
      void setquadrant(const fline& f, unsigned q);
    };
    class disttx {
      // bundle of data to pass along for distance
    public:
      real betw2, omgw2;
      int ind2;
    };
    fline() {}
    fline(const Triaxial& t, Triaxial::gamblk gm,
                  real epspow = 1, real nmaxmult = 0);
    const ffun& fbet() const { return _fbet; }
    const ffun& fomg() const { return _fomg; }
    const Triaxial& t() const { return _t; }
    real gamma() const { return _gm.gamma; }
    const Triaxial::gamblk& gm() const { return _gm; }
    real Hybrid0(const fics& ic,
                 const Angle& bet2, const Angle& omg2) const;
    disttx Hybrid(const fics& fic, const Angle& bet2,
                Angle& bet2a, Angle& omg2a, Angle& alp2a) const;
    disttx ArcPos0(const fics& fic, const Angle& tau12,
                   Angle& bet2a, Angle& omg2a, Angle& alp2a,
                   bool betp = true) const;
  };

  class GEOGRAPHICLIB_EXPORT gline {
  private:
    typedef Math::real real;
    typedef Angle ang;
    Triaxial _t;
    Triaxial::gamblk _gm;
    gfun _gbet, _gomg;
  public:
    real s0;
    class gics {
      // bundle of data setting the initial conditions for a distance calc
    public:
      real sig1, s13;           // starting point
      gics();
      gics(const gline& g,
          const fline::fics& fic);
    };
    gline() {}
    gline(const Triaxial& t, const Triaxial::gamblk& gam);
    const gfun& gbet() const { return _gbet; }
    const gfun& gomg() const { return _gomg; }
    const Triaxial& t() const { return _t; }
    real gamma() const { return _gm.gamma; }
    const Triaxial::gamblk& gm() const { return _gm; }
    real dist(gics ic, fline::disttx d) const;
  };

  class GEOGRAPHICLIB_EXPORT TriaxialLine {
  private:
    friend class fline;
    friend class gline;
    friend class ffun;
    friend class gfun;
    typedef Math::real real;
    typedef Angle ang;
    Triaxial _t;
    fline _f;
    fline::fics _fic;
    gline _g;
    gline::gics _gic;
    static void solve2(real f0, real g0,
                       const ffun& fx, const ffun& fy,
                       const gfun& gx, const gfun& gy,
                       real& x, real& y,
                       int* countn = nullptr, int* countb = nullptr);
    static void solve2u(real f0, real g0,
                        const ffun& fx, const ffun& fy,
                        const gfun& gx, const gfun& gy,
                        real& x, real& y,
                        int* countn = nullptr, int* countb = nullptr);
    static void newt2(real f0, real g0,
                      const ffun& fx, const ffun& fy,
                      const gfun& gx, const gfun& gy,
                      real x0, real xa, real xb,
                      real xscale, real zscale,
                      real& x, real& y,
                      int* countn = nullptr, int* countb = nullptr);
    static real clamp(real x, real mult = 1) {
      using std::fmax; using std::fmin;
      real z = mult * Triaxial::BigValue();
      return fmax(-z, fmin(z, x));
    }
    static real lamang0(Angle x) {
      // lam(x) when x is an ang -- no clamping
      using std::asinh; using std::fabs;
      return asinh(x.s()/fabs(x.c()));
    }
    static real lamang(Angle x) {
      // lam(x) when x is an ang -- with clamping
      // A consistent large value for x near pi/2.
      return clamp(lamang0(x));
    }
    static real EllipticThresh() {
      static real thresh = 1/real(8);
      return thresh;
    }

  public:
    // remainder with result in
    //    [-y/2, y/2) is alt = false
    //    (-y/2, y/2] is alt = true
    static std::pair<real, real> remx(real x, real y, bool alt = false) {
      using std::remainder; using std::round;
      real z = remainder(x, y);
      if (alt) {
        if (z == -y/2) z = y/2;
      } else {
        if (z == y/2) z = -y/2;
      }
      return std::pair<real, real>(z, round((x - z) / y));
    }
    static real roundx(real x, real y) {
      return remx(x, y).second;
    }
    TriaxialLine(const Triaxial& t) : _t(t) {}
    TriaxialLine(const Triaxial& t,
                 Angle bet1, Angle omg1, Angle alp1);
    TriaxialLine(const Triaxial& t, real bet1, real omg1, real alp1);
    TriaxialLine(fline f, fline::fics fic,
                 gline g, gline::gics gic);
    const ffun& fbet() const { return _f.fbet(); }
    const ffun& fomg() const { return _f.fomg(); }
    const gfun& gbet() const { return _g.gbet(); }
    const gfun& gomg() const { return _g.gomg(); }
    void Position(real s12, Angle& bet2, Angle& omg2, Angle& alp2,
                  int* countn = nullptr, int* countb = nullptr) const;
    void Position(real s12, real& bet2, real& omg2, real& alp2,
                  bool unroll = true,
                  int* countn = nullptr, int* countb = nullptr) const;
    // Find first crossing of bet = bet2 in the direction dir.  NaN returned if
    // no crossing.  Assume bet1 in [-90,0] (or maybe (-90,0]) and bet2 in
    // [bet2, -bet2].   Special cases (dir = +1 assumed), bet2 = bet1...
    // bet1 < 0, alp1 in [-90,90], omg2 = 0
    // bet1 == 0, alp1 in (-90,90), omg2 = 0,
    //                   alp1 = +/-90 omg2 = conj pt
    void Hybrid(const Angle& bet2,
                Angle& bet2a, Angle& omg2a, Angle& alp2a,
                real& s12) const;
    real gamma() const { return _f.gamma(); }
    real Distance() const { return _gic.s13; }
    void SetDistance(real s13) { _gic.s13 = s13; }
    void pos1(Angle& bet1, Angle& omg1, Angle& alp1) const;
    void pos1(real& bet1, real& omg1, real& alp1, bool unroll = true) const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALLINE_HPP

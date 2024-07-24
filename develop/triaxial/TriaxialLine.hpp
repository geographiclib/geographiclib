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

  class GEOGRAPHICLIB_EXPORT TriaxialLine {
  private:
    friend class Triaxial;
    typedef Math::real real;
    typedef Angle ang;

    class ffun {
    private:
      real _kap, _kapp, _eps, _mu, _sqrtkapp;
      bool _tx, _newumb;
      EllipticFunction _ell;
      TrigfunExt _fun;
      real _tol;
      int _nmax;
      Trigfun _dfinv;
      int _countn, _countb;
      real _max;
      bool _invp;
      // mu > 0
      static real fphip(real c, real kap, real kapp, real eps, real mu);
      static real fup(real cn, real kap, real kapp, real eps, real mu);
      // mu == 0
      static real dfp(real c, real kap, real kapp, real eps);
      static real dfvp(real cn, real dn, real kap, real kapp, real eps);
      static real newdfp(real c, real kap, real kapp, real eps);
      static real newdfvp(real cn, real dn, real kap, real kapp, real eps);
      // mu < 0
      static real fpsip(real s, real c, real kap, real kapp,
                        real eps, real mu);
      static real fvp(real dn, real kap, real kapp, real eps, real mu);
      real root(real z, real x0, int* countn, int* countb) const;
    public:
      ffun() {}
      ffun(real kap, real kapp, real eps, real mu, const Triaxial& t,
           real epsow = 1, real nmaxmult = 0);
      real operator()(real u) const {
        if (_mu == 0) {
          real phi = gd(u, _newumb ? _sqrtkapp : 1);
          return u - _fun(_tx ? _ell.F(phi) : phi);
        } else
          return _fun(u);
      }
      real deriv(real u) const {
        using std::cosh; using std::sinh; using std::sqrt;
        if (_mu == 0) {
          real phi = gd(u, _newumb ? _sqrtkapp : 1),
            // sch = dphi/du
            sch = _newumb ?
            _sqrtkapp * cosh(u) / (_kapp + Math::sq(sinh(u))) :
            1/cosh(u);
          // for _tx and _newumb dv/du = 1/sqrt(kapp+sinh(u)^2)
          return 1 - _fun.deriv(_tx ? _ell.F(phi) : phi) * sch /
            ( _tx ? (_newumb ?
                     // Use _ell.k2() == _kap here
                     _sqrtkapp * cosh(u) / sqrt(_kapp + Math::sq(sinh(u))) :
                     sqrt(_ell.kp2() + _ell.k2() * Math::sq(sch))) :
              1);
        } else
          return _fun.deriv(u);
      }

      // Approximate inverse using _finv
      real inv0(real z) const;
      // Accurate inverse by direct Newton (not using _finv)
      real inv1(real z, int* countn = nullptr, int* countb = nullptr) const;
      // Accurate inverse correcting result from _finv
      real inv2(real z, int* countn = nullptr, int* countb = nullptr) const;
      real inv(real z, int* countn = nullptr, int* countb = nullptr) const {
        return _invp ? inv2(z, countn, countb) : inv1(z, countn, countb);
      }
      void ComputeInverse();
      real fwd(real phi) const {
        return _mu == 0 ? lam(phi, _newumb ? _sqrtkapp : 1) :
          (_tx ? _ell.F(phi) : phi);
      }
      real rev(real u) const {
        return _mu == 0 ? gd(u, _newumb ? _sqrtkapp : 1) :
          (_tx ? _ell.am(u) : u);
      }
      int NCoeffs() const { return _fun.NCoeffs(); }
      int NCoeffsInv() const {
        return _mu == 0 ? _dfinv.NCoeffs() : _fun.NCoeffsInv();
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
      void DumpTable() const {
        int ndiv = 5;
        for (int i = -3*ndiv; i <= 3*ndiv; ++i) {
          real u = i/real(ndiv),
            f = (*this)(u),
            u0 = inv0(f),
            u1 = inv1(f),
            u2 = inv2(f);
          std::cout << i << " " << u << " " << f << " 0: "
                    << u0-u << " 1: " << u1-u << " 2: "
                    << u2-u << "\n";
        }
      }
    };

    class gfun {
    private:
      real _kap, _kapp, _eps, _mu, _sqrtkapp;
      bool _tx, _newumb, _gdag;
      EllipticFunction _ell;
      TrigfunExt _fun;
      real _max;
      // _mu > 0
      static real gphip(real c, real kap, real kapp, real eps, real mu);
      static real gfphip(real c, real kap, real mu);
      static real gup(real cn, real dn, real kap, real kapp,
                      real eps, real mu);
      static real gfup(real cn, real kap, real mu);
      // gdagger = g - mu * f variants
      static real gdagphip(real c, real kap, real kapp, real eps, real mu);
      static real gfdagphip(real c, real kap, real mu);
      static real gdagup(real cn, real dn, real kap, real kapp,
                         real eps, real mu);
      static real gfdagup(real cn, real kap, real mu);
      // _mu == 0
      static real g0p(real c, real kap, real kapp, real eps);
      // static real gf0p(real c, real kap, real kapp);
      static real gf0up(real u, real kap, real kapp);
      static real gf0upalt(real u, real kap, real kapp);
      static real g0vp(real cn, real kap, real kapp, real eps);
      // _mu < 0
      static real gpsip(real s, real c, real kap, real kapp,
                        real eps, real mu);
      static real gfpsip(real c, real kap, real mu);
      static real gvp(real cn, real dn, real kap, real kapp,
                      real eps, real mu);
      static real gfvp(real cn, real kap, real mu);
      // gdagger = g - mu * f variants
      static real gdagpsip(real s, real c, real kap, real kapp,
                           real eps, real mu);
      static real gfdagpsip(real s, real c, real kap, real mu);
      static real gdagvp(real cn, real dn, real kap, real kapp,
                         real eps, real mu);
      static real gfdagvp(real dn, real kap, real mu);
    public:
      gfun() {}
      gfun(real kap, real kapp, real eps, real mu, const Triaxial& t);
      real operator()(real u) const {
        if (_mu == 0) {
          real phi = gd(u, _newumb ? _sqrtkapp : 1);
          return _fun(_tx ? _ell.F(phi) : phi);
        } else
          return _fun(u);
      }

      real deriv(real u) const {
        using std::cosh; using std::sinh; using std::sqrt;
        if (_mu == 0) {
          real phi = gd(u, _newumb ? _sqrtkapp : 1),
            // sch = dphi/du
            sch = _newumb ?
            _sqrtkapp * cosh(u) / (_kapp + Math::sq(sinh(u))) :
            1/cosh(u);
          return _fun.deriv(_tx ? _ell.F(phi) : phi) * sch /
            ( _tx ? (_newumb ?
                     // Use _ell.k2() == _kap here
                     _sqrtkapp * cosh(u) / sqrt(_kapp + Math::sq(sinh(u))) :
                     sqrt(_ell.kp2() + _ell.k2() * Math::sq(sch))) :
              1);
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

    class fline {
    private:
      Triaxial _t;
      Triaxial::gamblk _gm;
      ffun _fbet, _fomg;
      bool _invp;
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
        real u0, v0, delta;     // starting point geodesic
        int nN, eE;             // Northgoing / eastgoing
        fics() {}
        fics(const fline& f,
             Angle bet1, Angle omg1, Angle alp1);

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
                   Angle bet2, Angle omg2) const;
      disttx Hybrid(const fics& fic, Angle bet2,
                    Angle& bet2a, Angle& omg2a, Angle& alp2a) const;
      disttx ArcPos0(const fics& fic, Angle tau12,
                     Angle& bet2a, Angle& omg2a, Angle& alp2a,
                     bool betp = true) const;
      void ComputeInverse();
    };

    class gline {
    private:
      Triaxial _t;
      Triaxial::gamblk _gm;
      gfun _gbet, _gomg;
    public:
      real s0;
      class gics {
        // bundle of data setting the initial conditions for a distance calc
      public:
        real sig1, s13;         // starting point
        gics() {}
        gics(const gline& g, const fline::fics& fic);
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
    static real lamang0(Angle x, real mult = 1) {
      // lam(x) when x is an ang -- no clamping
      using std::asinh; using std::fabs;
      return asinh(mult * x.t());
    }
    static real lamang(Angle x, real mult = 1) {
      // lam(x) when x is an ang -- with clamping
      // A consistent large value for x near pi/2.
      return clamp(lamang0(x, mult));
    }
    static real lam(real x, real mult = 1) {
      using std::tan; using std::asinh; using std::fabs;
      // A consistent large value for x near pi/2.  Also deals with the issue
      // that tan(pi/2) may be negative, e.g., for long doubles.
      return fabs(x) < Math::pi()/2 ? asinh(mult * tan(x)) :
        (x < 0 ? -1 : 1) * Triaxial::BigValue();
    }
    static real gd(real x, real mult = 1) {
      using std::atan; using std::sinh;
      return atan(sinh(x) / mult);
    }
    static ang anglam(real u, real mult = 1) {
      using std::sinh;
      return Angle(sinh(u), mult, 0);
    }
    static real mcosh(real u, real mult = 1) {
      using std::cosh; using std::sinh; using std::hypot;
      return mult == 1 ? cosh(u) : hypot(sinh(u), mult) / mult;
    }
  public:
    const ffun& fbet() const { return _f.fbet(); }
    const ffun& fomg() const { return _f.fomg(); }
  private:
    const gfun& gbet() const { return _g.gbet(); }
    const gfun& gomg() const { return _g.gomg(); }

    // remainder with result in
    //    [-y/2, y/2) is alt = false
    //    (-y/2, y/2] is alt = true
    static std::pair<real, real> remx(real x, real y, bool alt = false) {
      using std::remainder; using std::rint;
      real z = remainder(x, y);
      if (alt) {
        if (z == -y/2) z = y/2;
      } else {
        if (z == y/2) z = -y/2;
      }
      return std::pair<real, real>(z, rint((x - z) / y));
    }
    TriaxialLine(fline f, fline::fics fic,
                 gline g, gline::gics gic);

  public:
    TriaxialLine(const Triaxial& t) : _t(t) {}
    TriaxialLine(const Triaxial& t,
                 Angle bet1, Angle omg1, Angle alp1);
    TriaxialLine(const Triaxial& t, real bet1, real omg1, real alp1);
    void Position(real s12, Angle& bet2, Angle& omg2, Angle& alp2,
                  int* countn = nullptr, int* countb = nullptr) const;
    void Position(real s12, real& bet2, real& omg2, real& alp2,
                  bool unroll = true,
                  int* countn = nullptr, int* countb = nullptr) const;
    // Find first crossing of bet = bet2.  NaN returned if
    // no crossing.  Assume bet1 in [-90,0] (or maybe (-90,0]) and bet2 in
    // [bet2, -bet2].   Special cases, bet2 = bet1...
    // bet1 < 0, alp1 in [-90,90], omg2 = 0
    // bet1 == 0, alp1 in (-90,90), omg2 = 0,
    //                   alp1 = +/-90 omg2 = conj pt
    void Hybrid(Angle bet2,
                Angle& bet2a, Angle& omg2a, Angle& alp2a,
                real& s12) const;
    real gamma() const { return _f.gamma(); }
    real Distance() const { return _gic.s13; }
    void SetDistance(real s13) { _gic.s13 = s13; }
    void Offset(real s13, bool reverse = false);
    void pos1(Angle& bet1, Angle& omg1, Angle& alp1) const;
    void pos1(real& bet1, real& omg1, real& alp1, bool unroll = true) const;
    void Optimize();

    real fbetm() const { return gamma() != 0 ?
        fbet().HalfPeriod() * fbet().Slope() :
        fbet().Max();
    }
    real fomgm() const { return gamma() != 0 ?
        fomg().HalfPeriod() * fomg().Slope() :
        fomg().Max();
    }

    real gbetm() const { return gamma() != 0 ?
        gbet().HalfPeriod() * gbet().Slope() :
        gbet().Max();
    }
    real gomgm() const { return gamma() != 0 ?
        gomg().HalfPeriod() * gomg().Slope() :
        gomg().Max();
    }
    real df() const { return _f.df; }
    real deltashift() const { return _f.deltashift; }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALLINE_HPP

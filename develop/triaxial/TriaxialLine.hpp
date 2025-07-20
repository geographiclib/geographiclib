/**
 * \file TriaxialLine.hpp
 * \brief Header for GeographicLib::TriaxialLin class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIALLINE_HPP)
#define GEOGRAPHICLIB_TRIAXIALLINE_HPP 1

#include <utility>
#include <ostream>
#include <set>
#include <algorithm>
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

    class hfun {
      // This combines ffun abd gfun in order to minimize the duplication of
      // code.
      // Establish consistent notation for coordinates:
      // theta = rotating coord, omega-90 or beta for circum- or transpolar
      // phi = librating coord, beta or omega-90 for circum- or transpolar
      // psi = rotating equivalent of phi
      // zeta = a generic version of theta or psi
      // u = possible change of vars for theta
      // v = possible change of vars for psi
      // w = a generic version of u or v
    private:
      real _kap, _kapp, _eps, _mu, _sqrtkap, _sqrtkapp;
      bool _distp, _tx;
      EllipticFunction _ell;
      TrigfunExt _fun;
      // real _tol;
      Trigfun _dfinv;
      int _countn, _countb;
      real _max;
      bool _umb, _meridr, _meridl, _biaxr, _biaxl, _invp;
      // The f functions
      // mu > 0
      static real fthtp(real c, real kap, real kapp, real eps, real mu);
      static real fup(real cn, real kap, real kapp, real eps, real mu);
      // mu == 0
      static real dfp(real c, real kap, real kapp, real eps);
      static real dfvp(real cn, real dn, real kap, real kapp, real eps);
      // mu < 0
      static real fpsip(real s, real c, real kap, real kapp,
                        real eps, real mu);
      static real fvp(real dn, real kap, real kapp, real eps, real mu);
      // biaxial variant for kap = 0, kapp = 1
      static real fthtbiax(real tht, real eps, real mu);
      // biaxial variant for kap = 1, kapp = 0, mu <= 0, !_tx
      static real dfpsibiax(real s, real c, real eps, real mu);
      // NOT USED
      // biaxial variant for kap = 1, kapp = 0, mu <= 0, _tx
      // static real dfvbiax(real dn, real eps, real mu);

      // The g functions
      // _mu > 0
      static real gthtp(real c, real kap, real kapp, real eps, real mu);
      static real gup(real cn, real dn, real kap, real kapp,
                      real eps, real mu);
      // _mu == 0
      static real g0p(real c, real kap, real kapp, real eps);
      static real g0vp(real cn, real kap, real kapp, real eps);
      // _mu < 0
      static real gpsip(real s, real c, real kap, real kapp,
                        real eps, real mu);
      static real gvp(real cn, real dn, real kap, real kapp,
                      real eps, real mu);
      // biaxial variants for kap = 0, kapp = 1, mu >= 0
      static real gthtbiax(real tht, real eps, real mu);
      // biaxial variants for kap = 1, kapp = 0, mu <= 0
      static real gpsibiax(real s, real c, real eps, real mu);
      // NOT USED
      // static real gvbiax(real cn, real dn, real eps, real mu);

      // Return atan(m * tan(x)) keeping result continuous.  Only defined for
      // !signbit(m).
      static real modang(real x, real m) {
        return Angle::radians(x).modang(m).radians();
      }
      static real modang(real x, real m, real& deriv) {
        using std::fabs;
        Angle xa = Angle::radians(x), ya = xa.modang(m);
        // y = modang(x, m)
        // tan(y) = m * tan(x)
        // dy/dx = m * sec(x)^2 / (1 + m^2 * tan(x)^2)
        //       = m / (cos(x)^2 + m^2 * sin(x)^2)
        //       = m * csc(x)^2 / (m^2 + cot(x)^2)
        //       = m * cos(y)^2/cos(x)^2
        //       = sin(y)^2 / (m * sin(x)^2)
        //       -> m for x -> 0
        //       -> 1/m for x -> +/- pi/2
        deriv = fabs(2*xa.s() > 1) ? // just need to avoid s() close to zero
          Math::sq(ya.s()/xa.s()) / m :
          Math::sq(ya.c()/xa.c()) * m;
        return ya.radians();
      }
      real root(real z, real u0, int* countn, int* countb,
                real tol = std::numeric_limits<real>::epsilon()) const;

      // Approximate inverse using _finv
      real inv0(real z) const;
      // Accurate (to tol) inverse by direct Newton (not using _finv)
      real inv1(real z, int* countn = nullptr, int* countb = nullptr) const;
      // Accurate inverse correcting result from _finv
      real inv2(real z, int* countn = nullptr, int* countb = nullptr) const;

    public:
      // Summary of f and g functions
      // _meridr = kap == 0 && mu == 0 (biaxial meridian rotating coordinate)
      //  f = fthtbiax / sqrt(mu), fthtbiax = 1 inversion trivial
      //  g = gthtbiax, gthtbiax = 0
      // _meridl = kapp == 0 && mu -= 0 (biaxial meridian librating coordinate)
      //  f = (atan(sqrt(-mu)*tan(psi)) - dfpsiobl) / sqrt(-mu)
      //    mu == 0 **
      //     atan(sqrt(-mu)*tan(psi)) -> round(psi/pi)*pi
      //  g = gpsibiax
      // mu > 0 (rotating coordinate)
      //  f = ftht or fu
      //  g = gtht or gu **
      // mu < 0 (librating coordinate)
      //  f = fpsi or fv
      //  g = gpsi or gv **
      // _umb = !biaxial && mu == 0 (triaxial umbilic)
      //  f = (u + df)/sqrt(kap*kapp) or (u + dfv)/sqrt(kap*kapp)
      //  g = g0 or g0v
      //
      // Handing of f funtions needs to be handled specially for
      // _merid[lr]
      //   leading order behavior is seperated out
      //   N.B. inverse of f is discontinuous for mu == 0
      //   functions which need special treatment
      //     operator(), deriv, inv, ComputeInverse, Slope, Max
      // _umb
      //   leading order behavior is seperated out
      // otherwise everything is handled by generic routines
      //
      // Handling of g function is always generic
      // _umb
      //   functions which need special treatment
      //    Slope, fwd, rev, Max
      //
      hfun() {}
      hfun(bool distp, real kap, real kapp, real eps, real mu,
           const Triaxial& t);
      real operator()(real u) const;
      // THIS ISN"T USED
      // Angle operator()(const Angle& ang) const;
      real deriv(real u) const;
      real df(real u) const { return _fun(u); }
      real dfp(real u) const { return _fun.deriv(u); }

      real inv(real z, int* countn = nullptr, int* countb = nullptr) const {
        return _invp ? inv2(z, countn, countb) : inv1(z, countn, countb);
      }
      // Angle inv(const Angle& z,
      //           int* countn = nullptr, int* countb = nullptr) const;
      void ComputeInverse();
      real fwd(real zeta) const {
        return _umb ? lam(zeta, _sqrtkapp) :
          (_tx ? _ell.F(zeta) : zeta);
      }
      real rev(real w) const {
        return _umb ? gd(w, _sqrtkapp) :
          (_tx ? _ell.am(w) : w);
      }
      int NCoeffs() const { return _fun.NCoeffs(); }
      int NCoeffsInv() const {
        return _umb ? _dfinv.NCoeffs() : _fun.NCoeffsInv();
      }
      std::pair<int, int> InvCounts() const {
        return _umb ? std::pair<int, int>(_countn, _countb) :
          _fun.InvCounts();
      }
      bool txp() const { return _tx; }
      real HalfPeriod() const {
        return _umb ? Math::infinity() : (_tx ? _ell.K() : Math::pi()/2);
      }
      real Slope() const {
        using std::sqrt;
        return _umb ? 1 :
          !_distp && _meridl ? 0 :
          !_distp && _biaxl ? 1 - sqrt(-_mu) * _fun.Slope() :
          _fun.Slope();
      }
      real Max() const { return _max; }
      real MaxPlus() const {
        using std::fmax;
        return fmax(_max, HalfPeriod() * (Slope() == 0 ? 1 : Slope()) / 1000);
      }
      void inversedump(std::ostream& os, const std::string& name) const;
    };

    class fline {
    private:
      Triaxial _t;
      Triaxial::gamblk _gm;
      real _deltashift;
      hfun _fpsi, _ftht;
      bool _invp;
    public:
      class fics {
        // bundle of data setting the initial conditions for a geodesic
      public:
        //        bool transpolar;
        // alp1 is angle measured from line of const rotating coording
        Angle tht1, phi1, alp1, // rotating, librating starting point
          psi1,                 // phi1 transformed to rotating angle psi
        // Angles about which quantities oscillate
        // circumpolar:
        //   omg0 not used, bet0 = cardinal(even), alp0 = cardinal(odd)
        // transpolar:
        //   bet0 not used, omg0 = cardinal(even), alp0 = cardinal(even)
        // umbilical
        //   bet0, omg0 = cardinal(even) = middle of starting segment
        //   !umbalt: alp0 = cardinal(odd)
        //   umbalt: alp0 = cardinal(even)
          tht0, phi0, alp0;
        // Angle versions of u0, v0, delta, defer for now
        // Angle u0a, v0a, deltaa;
        real u0, v0, delta;     // starting point geodesic
        int Nx, Ex;             // Northgoing / eastgoing relative to tht
        fics() {}
        fics(const fline& f,
             Angle bet1, Angle omg1, Angle alp1);
        void setquadrant(const fline& f, unsigned q);
        void pos1(bool transpolar,
                  Angle& bet10, Angle& omg10, Angle& alp10) const;
      };
      class disttx {
        // bundle of data to pass along for distance
      public:
        real phiw2, thtw2;
        int ind2;
      };
      fline() {}
      fline(const Triaxial&t, bool neg = false);
      fline(const Triaxial& t, Triaxial::gamblk gm);
      const hfun& fpsi() const { return _fpsi; }
      const hfun& ftht() const { return _ftht; }
      const hfun& fbet() const {
        return !transpolar() ? fpsi() : ftht();
      }
      const hfun& fomg() const {
        return transpolar() ? fpsi() : ftht();
      }
      const Triaxial& t() const { return _t; }
      real gamma() const { return _gm.gamma; }
      real gammax() const { return _gm.gammax; }
      real kx2() const { return _gm.kx2; }
      real kxp2() const { return _gm.kxp2; }
      real kx() const { return _gm.kx; }
      real kxp() const { return _gm.kxp; }
      real nu() const { return _gm.nu; }
      real nup() const { return _gm.nup; }
      real deltashift() const { return _deltashift; }
      bool transpolar() const { return _gm.transpolar; }
      const Triaxial::gamblk& gm() const { return _gm; }
      // Run fline to its first intersection with
      // (for betp) bet2 and return omg2
      // (for !betp) omg2 and return bet2
      real Hybrid0(const fics& fic, Angle bet2, Angle omg2,
                   bool betp = true) const;
      // Run fline to its first intersection with bet and return resulting
      // bet2a, omg2a, alp2a (without angle normalization) and distance
      // calculation object
      disttx Hybrid(const fics& fic, Angle betomg2,
                    Angle& bet2a, Angle& omg2a, Angle& alp2a,
                    bool betp = true) const;
      disttx ArcPos0(const fics& fic, Angle tau12,
                     Angle& bet2a, Angle& omg2a, Angle& alp2a,
                     bool betp = true) const;
      void ComputeInverse();
      void inversedump(std::ostream& os) const;
    };

    class gline {
    private:
      Triaxial _t;
      Triaxial::gamblk _gm;
      hfun _gpsi, _gtht;
      bool _invp;
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
      gline(const Triaxial& t, bool neg = false);
      gline(const Triaxial& t, const Triaxial::gamblk& gm);
      real gamma() const { return _gm.gamma; }
      real gammax() const { return _gm.gammax; }
      real kx2() const { return _gm.kx2; }
      real kxp2() const { return _gm.kxp2; }
      real kx() const { return _gm.kx; }
      real kxp() const { return _gm.kxp; }
      bool transpolar() const { return _gm.transpolar; }
      const hfun& gpsi() const { return _gpsi; }
      const hfun& gtht() const { return _gtht; }
      const hfun& gbet() const {
        return !transpolar() ? gpsi() : gtht();
      }
      const hfun& gomg() const {
        return transpolar() ? gpsi() : gtht();
      }
      const Triaxial& t() const { return _t; }
      const Triaxial::gamblk& gm() const { return _gm; }
      real dist(gics ic, fline::disttx d) const;
      void ComputeInverse();
      void inversedump(std::ostream& os) const;
    };

    class zvals {
    public:
      real z, fz, gz;
      zvals(real z0 = 0, real fz0 = 0, real gz0 = 0)
        : z(z0), fz(fz0), gz(gz0) {}
      bool operator<(const zvals& t) const { return z < t.z; }
      bool operator==(const zvals& t) const { return z == t.z; }
    };

    class zset {
    private:
      std::vector<zvals> _s;
    public:
      zset(const zvals& a, const zvals& b)
        : _s({a, b})
      {
        if (a == b)
          // Allow coincident start and end values
          _s.resize(1);
        else if (!(a < b && a.fz <= b.fz && a.gz <= b.gz))
          throw GeographicLib::GeographicErr("bad zset initializer");
      }
      int num() const { return _s.size(); }
      const zvals& val(int i) const { return _s[i]; }
      const zvals& min() const { return _s[0]; }
      const zvals& max() const { return _s.back(); }
      int insert(zvals& t, int flag = 0);
      real bisect() const {
        // return z in the middle of biggest gap
        if (num() == 1)
          return min().z;
        real maxgap = -1; int maxind = 0;
        for (int i = 0; i < num() - 1; ++i) {
          real gap = _s[i+1].z - _s[i].z;
          if (gap > maxgap) {
            maxgap = gap;
            maxind = i;
          }
        }
        return (_s[maxind].z + _s[maxind+1].z) / 2;
      }
    };

    Triaxial _t;
    fline _f;
    fline::fics _fic;
    gline _g;
    gline::gics _gic;

    static void solve2(real f0, real g0,
                       const hfun& fx, const hfun& fy,
                       const hfun& gx, const hfun& gy,
                       real& x, real& y,
                       int* countn = nullptr, int* countb = nullptr);
    static void solve2u(real f0, real g0,
                        const hfun& fx, const hfun& fy,
                        const hfun& gx, const hfun& gy,
                        real& x, real& y,
                        int* countn = nullptr, int* countb = nullptr);
    static void newt2(real f0, real g0,
                      const hfun& fx, const hfun& fy,
                      const hfun& gx, const hfun& gy,
                      real xa, real xb, real xscale,
                      real ya, real yb, real yscale,
                      real fscale, real gscale,
                      real& x, real& y,
                      int* countn = nullptr, int* countb = nullptr);
    static void zsetsinsert(zset& xset, zset& yset,
                            zvals& xfg, zvals& yfg,
                            real f0, real g0);
    static void zsetsdiag(const zset& xset, const zset& yset,
                          real f0, real g0);
    static std::pair<real, real> zsetsbisect(const zset& xset, const zset& yset,
                                             real f0, real g0);
    static real clamp(real x, real mult = 1) {
      using std::fmax; using std::fmin;
      real z = mult * Triaxial::BigValue();
      return Math::clamp(x, -z, z);
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
      using std::tan; using std::asinh; using std::fabs; using std::copysign;
      // A consistent large value for x near pi/2.  Also deals with the issue
      // that tan(pi/2) may be negative, e.g., for long doubles.
      return fabs(x) >= Math::pi()/2 ? copysign(Triaxial::BigValue(), x) :
        asinh(mult * tan(x));
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

    const hfun& fbet() const { return _f.fbet(); }
    const hfun& fomg() const { return _f.fomg(); }
    const hfun& fpsi() const { return _f.fpsi(); }
    const hfun& ftht() const { return _f.ftht(); }
    const hfun& gbet() const { return _g.gbet(); }
    const hfun& gomg() const { return _g.gomg(); }
    const hfun& gpsi() const { return _g.gpsi(); }
    const hfun& gtht() const { return _g.gtht(); }

    // remainder with result in
    //    [-y/2, y/2) if alt = false (default)
    //    (-y/2, y/2] if alt = true
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
    // remainder with x in
    //    [-pi/2, pi/2) if alt = false (default)
    //    (-pi/2, pi/2] if alt = true
    // equivalent to remx(x.radians(), Math::pi(), alt)
    static std::pair<real, real> remx(Angle x, bool alt = false) {
      using std::signbit;
      real m = 0;
      if (signbit(x.c())) {
        x -= ang::cardinal(2);
        ++m;
      }
      if (x.c() == 0) {
        if (alt && signbit(x.s())) {
          x += ang::cardinal(2);
          --m;
        } else if (!alt && !signbit(x.s())) {
          x -= ang::cardinal(2);
          ++m;
        }
      }
      return  std::pair<real, real>(x.radians0(), m + 2*x.n());
    }
    static bool biaxspecial(const Triaxial& t, real gamma) {
      using std::fabs;
      return t._biaxp && (t._k2 == 0 || t._kp2 == 0) &&
        gamma != 0 && fabs(gamma) < t._ellipthresh;
    }
    // Private constructor to assemble the pieces of the class on exiting
    // Triaxial::Inverse.
    TriaxialLine(fline f, fline::fics fic, gline g, gline::gics gic);
    // Private constructor to provide the umbilical fline object.
    TriaxialLine(const Triaxial& t);
  public:
    // TriaxialLine(const Triaxial& t) : _t(t), _f(_t), _g(_t) {}
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
    void Hybrid(Angle betomg2,
                Angle& bet2a, Angle& omg2a, Angle& alp2a,
                real& s12, bool betp = true) const;
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
    //    real df() const { return _f.df; }
    real deltashift() const { return _f.deltashift(); }
    void inversedump(std::ostream& os) const;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALLINE_HPP

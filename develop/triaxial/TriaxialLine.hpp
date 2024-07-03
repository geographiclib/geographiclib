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

  class GEOGRAPHICLIB_EXPORT TriaxialLineF {
  private:
    typedef Math::real real;
    typedef Angle ang;
    Triaxial _t;
    Triaxial::gamblk _gm;
    geod_fun _fbet, _fomg;
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
      real u0, v0,              // starting point in u,v space
        delta;                   //  starting point for umbilic
      int nN, eE,                  // Northgoing / eastgoing
        flip;                    // Is bet or omg on the backside (non-umb)
      fics();
      fics(const TriaxialLineF& f,
           const Angle& bet1, const Angle& omg1, const Angle& alp1);
      void setquadrant(const TriaxialLineF& f, unsigned q);
    };
    class disttx {
      // bundle of data to pass along for distance
    public:
      real betw2, omgw2;
      int ind2;
    };
    TriaxialLineF() {}
    TriaxialLineF(const Triaxial& t, Triaxial::gamblk gm,
                  real epspow = 1, real nmaxmult = 0);
    const geod_fun& fbet() const { return _fbet; }
    const geod_fun& fomg() const { return _fomg; }
    const Triaxial& t() const { return _t; }
    real gamma() const { return _gm.gam; }
    const Triaxial::gamblk& gm() const { return _gm; }
    real Hybrid0(const fics& ic,
                 const Angle& bet2, const Angle& omg2) const;
    disttx Hybrid(const fics& fic, const Angle& bet2,
                Angle& bet2a, Angle& omg2a, Angle& alp2a) const;
    disttx ArcPos0(const fics& fic, const Angle& tau12,
                   Angle& bet2a, Angle& omg2a, Angle& alp2a,
                   bool betp = true) const;
  };

  class GEOGRAPHICLIB_EXPORT TriaxialLineG {
  private:
    typedef Math::real real;
    typedef Angle ang;
    Triaxial _t;
    Triaxial::gamblk _gm;
    dist_fun _gbet, _gomg;
  public:
    real s0;
    class gics {
      // bundle of data setting the initial conditions for a distance calc
    public:
      real sig1, s13;           // starting point
      gics();
      gics(const TriaxialLineG& g,
          const TriaxialLineF::fics& fic);
    };
    TriaxialLineG() {}
    TriaxialLineG(const Triaxial& t, const Triaxial::gamblk& gam);
    const dist_fun& gbet() const { return _gbet; }
    const dist_fun& gomg() const { return _gomg; }
    const Triaxial& t() const { return _t; }
    real gamma() const { return _gm.gam; }
    const Triaxial::gamblk& gm() const { return _gm; }
    real dist(gics ic, TriaxialLineF::disttx d) const;
  };

  class GEOGRAPHICLIB_EXPORT TriaxialLine {
  private:
    typedef Math::real real;
    typedef Angle ang;
    Triaxial _t;
    TriaxialLineF _f;
    TriaxialLineF::fics _fic;
    TriaxialLineG _g;
    TriaxialLineG::gics _gic;
    static void solve2(real f0, real g0,
                       const geod_fun& fx, const geod_fun& fy,
                       const dist_fun& gx, const dist_fun& gy,
                       real& x, real& y,
                       int* countn = nullptr, int* countb = nullptr);
    static void solve2u(real f0, real g0,
                        const geod_fun& fx, const geod_fun& fy,
                        const dist_fun& gx, const dist_fun& gy,
                        real& x, real& y,
                        int* countn = nullptr, int* countb = nullptr);
    static void newt2(real f0, real g0,
                      const geod_fun& fx, const geod_fun& fy,
                      const dist_fun& gx, const dist_fun& gy,
                      real x0, real xa, real xb,
                      real xscale, real zscale,
                      real& x, real& y,
                      int* countn = nullptr, int* countb = nullptr);

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
    TriaxialLine(TriaxialLineF f, TriaxialLineF::fics fic,
                 TriaxialLineG g, TriaxialLineG::gics gic);
    const geod_fun& fbet() const { return _f.fbet(); }
    const geod_fun& fomg() const { return _f.fomg(); }
    const dist_fun& gbet() const { return _g.gbet(); }
    const dist_fun& gomg() const { return _g.gomg(); }
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

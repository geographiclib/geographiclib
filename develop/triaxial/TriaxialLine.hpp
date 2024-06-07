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
#include <GeographicLib/AuxAngle.hpp>
#include "Trigfun.hpp"
#include "Triaxial.hpp"

namespace GeographicLib {

  class GEOGRAPHICLIB_EXPORT TriaxialLineF {
  private:
    typedef Math::real real;
    Triaxial _t;
    Triaxial::gamblk _gm;
    geod_fun _fbet, _fomg;
  public:
    class ics {
      // bundle of data setting the initial conditions for a geodesic
    public:
      AuxAngle bet1, omg1, alp1, // starting point
        psi1;                    // nonumbilic angle psi
      real u0, v0,               // starting point in u,v space
        df, deltashift,          // umbilic constants
        deltamax,                // max value of delta for umbilic
        delta;                   //  starting point for umbilic
      int ibet, iomg, ialp,      // wrapping quantities for bet, omg, alp
        nN, eE,                  // Northgoing / eastgoing
        flip,                    // Is bet or omg on the backside
        bet0, omg0, alp0;        // Reference vals (1 = 0, -1 = 180)
      bool umbalt;               // how coordinates wrap with umbilical lines
      ics();
      ics(const TriaxialLineF& f,
          const AuxAngle& bet1, const AuxAngle& omg1, const AuxAngle& alp1);
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
    disttx Hybrid(const ics& fi, const AuxAngle& bet2,
                AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a) const;
  };

  class GEOGRAPHICLIB_EXPORT TriaxialLineG {
  private:
    typedef Math::real real;
    Triaxial _t;
    Triaxial::gamblk _gm;
    dist_fun _gbet, _gomg;
  public:
    class ics {
      // bundle of data setting the initial conditions for a distance calc
    public:
      real s0, sig1;           // starting point
      ics();
      ics(const TriaxialLineG& g,
          const TriaxialLineF::ics& fic);
    };
    TriaxialLineG() {}
    TriaxialLineG(const Triaxial& t, const Triaxial::gamblk& gam);
    const dist_fun& gbet() const { return _gbet; }
    const dist_fun& gomg() const { return _gomg; }
    const Triaxial& t() const { return _t; }
    real gamma() const { return _gm.gam; }
    const Triaxial::gamblk& gm() const { return _gm; }
    real dist(ics ic, TriaxialLineF::disttx d) const;
  };

  class GEOGRAPHICLIB_EXPORT TriaxialLine {
  private:
    typedef Math::real real;
    Triaxial _t;
    Triaxial::gamblk _gm;
    //    real _gammao, _rtgam, _gamp, _rtgamp;
    AuxAngle _bet1, _omg1, _alp1;
    int _ibet, _iomg, _ialp;
    TriaxialLineF _f;
    TriaxialLineG _g;
    int _nN, _eE,                 // Northgoing / eastgoing
      _flip,                    // Is bet or omg on the backside
      _bet0, _omg0, _alp0;      // Reference vals (1 = 0, -1 = 180)
    AuxAngle _psi1;             // nonumbilic constants
    real _u0, _v0;
    real _df, _deltashift, _s0; // umbilic constants
    real _deltamax;             // max value of delta for umbilic
    real _delta, _sig1;         //  starting point
    bool _umbalt, _distinit;
    real _s13;                  // Distance to the reference point 3
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
                 AuxAngle bet1, AuxAngle omg1, AuxAngle alp1);
    TriaxialLine(const Triaxial& t, real bet1, real omg1, real alp1);
    const geod_fun& fbet() const { return _f.fbet(); }
    const geod_fun& fomg() const { return _f.fomg(); }
    const dist_fun& gbet() const { return _g.gbet(); }
    const dist_fun& gomg() const { return _g.gomg(); }
    void distinit();
    void Position(real s12, AuxAngle& bet2, AuxAngle& omg2, AuxAngle& alp2,
                  int* ibet2 = nullptr, int* iomg2 = nullptr,
                  int* ialp2 = nullptr,
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
    void Hybrid(const AuxAngle& bet2, int dir,
                AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a,
                real& s12) const;
    real gamma() const { return _gm.gam; }
    real Distance() const { return _s13; }
    void SetDistance(real s13) { _s13 = s13; }
  };
  class GEOGRAPHICLIB_EXPORT BareLine {
    // Method for computing the course of the line only (not distance).
    // Constructor sets starting point.  The Hybrid method returns tau12, bet2,
    // omg2, alp2, betw2, omgw2, ind2.  Reinit resets starting point from the
    // results of Hybrid, possibly (a) reversing the direction, (b) changing
    // sign of bet, (c) flipping beta / omega to another sheet.  Promote
    // returns a TriaxialLine with distance method initialized and s13 computed
    // according to tau12.
  public:
    TriaxialLineF f;
    TriaxialLineF::ics fi;
  };
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALLINE_HPP

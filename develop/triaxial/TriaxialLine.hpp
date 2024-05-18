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
#include "Trigfun.hpp"
#include "Triaxial.hpp"
#include <GeographicLib/AuxAngle.hpp>

namespace GeographicLib {

  class GEOGRAPHICLIB_EXPORT TriaxialLineF {
  private:
    typedef Math::real real;
    real _k2, _kp2, _e2, _gam;
    geod_fun _fbet, _fomg;
  public:
    TriaxialLineF() {}
    TriaxialLineF(real k2, real kp2, real e2, real gam,
                  real epspow = 1, real nmaxmult = 0);
    const geod_fun& fbet() const { return _fbet; }
    const geod_fun& fomg() const { return _fomg; }
    geod_fun& fbet() { return _fbet; }
    geod_fun& fomg() { return _fomg; }

  };

  class GEOGRAPHICLIB_EXPORT TriaxialLineG {
  private:
    typedef Math::real real;
    real _k2, _kp2, _e2, _gam;
    dist_fun _gbet, _gomg;
  public:
    TriaxialLineG() {}
    TriaxialLineG(real k2, real kp2, real e2, real gam);
    const dist_fun& gbet() const { return _gbet; }
    const dist_fun& gomg() const { return _gomg; }

  };

  class GEOGRAPHICLIB_EXPORT TriaxialLine {
  private:
    typedef Math::real real;
    Triaxial _t;
    AuxAngle _bet1, _omg1, _alp1;
    long _ibet, _iomg, _ialp;
    real _gam;
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
    TriaxialLine(const Triaxial& t) : _t(t) {}
    TriaxialLine(const Triaxial& t,
                 const AuxAngle& bet1,
                 const AuxAngle& omg1,
                 const AuxAngle& alp1);
    TriaxialLine(const Triaxial& t, real bet1, real omg1, real alp1)
      : TriaxialLine(t,
                     AuxAngle::degrees(bet1),
                     AuxAngle::degrees(omg1),
                     AuxAngle::degrees(alp1))
    {
      using std::round;
      _ibet = long(round((bet1 - _bet1.degrees()) / Math::td));
      _iomg = long(round((omg1 - Math::qd - _omg1.degrees()) / Math::td));
      _ialp = long(round((alp1 - _alp1.degrees()) / Math::td));
    }
    const geod_fun& fbet() const { return _f.fbet(); }
    const geod_fun& fomg() const { return _f.fomg(); }
    const dist_fun& gbet() const { return _g.gbet(); }
    const dist_fun& gomg() const { return _g.gomg(); }
    void distinit();
    std::pair<long, long> Position(real s12, AuxAngle& bet2, AuxAngle& omg2, AuxAngle& alp2,
                  int* countn = nullptr, int* countb = nullptr) const;
    void Position(real s12, real& bet2, real& omg2, real& alp2,
                  int* countn = nullptr, int* countb = nullptr) const;
  };
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALLINE_HPP

/**
 * \file Geodesic3.hpp
 * \brief Header for GeographicLib::Triaxial::Geodesic3 class
 *
 * Copyright (c) Charles Karney (2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEODESIC3_HPP)
#define GEOGRAPHICLIB_GEODESIC3_HPP 1

#include <functional>
#include <memory>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs vector
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {
  namespace Triaxial {
  class GeodesicLine3;

  class GEOGRAPHICLIB_EXPORT Geodesic3 {
  private:
    // For access to BigValue, _ellipthresh, _biaxp
    friend class GeodesicLine3;
    typedef Math::real real;
    typedef Angle ang;
    Ellipsoid3 _t;

    // Run geodesic from bet1, omg1, alp1, find its first intersection with bet
    // = bet2a and return omg2a - omg2b
    real HybridA(ang bet1, ang omg1, ang alp1,
                 ang bet2a, ang omg2b, bool betp) const;
    static ang findroot(const std::function<real(const ang&)>& f,
                          ang xa,  ang xb,
                          real fa, real fb,
                          int* countn = nullptr, int* countb = nullptr);
    bool _umbalt,               // how coordinates wrap with umbilical lines
      _biaxp,                   // special treatment for biaxial non-meridional
      _debug,                   // print out diagnostics
      _hybridalt,               // favor hybrid solution in terms of omg
      _swapomg;                 // allow swapping of omega{1,2}
    // If k'^2 < ellipthresh transform phi -> F(phi, k^2)
    real _ellipthresh;
    mutable std::shared_ptr<GeodesicLine3> _umbline;
    static real BigValue() {
      using std::log;
      static real bigval = -3*log(std::numeric_limits<real>::epsilon());
      return bigval;
    }
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
      gamblk(const Geodesic3& tg, bool neg = false);
      gamblk(const Geodesic3& tg, ang bet, ang omg, ang alp);
      //       gamblk(real gammax, real nux, real nupx)
      //        : gamma(gammax), nu(nux), nup(nupx) {}
    };
    gamblk gamma(ang bet, ang omg, ang alp)
      const;
    // real a() const { return t().a(); } // not needed
    real b() const { return t().b(); }
    // real c() const { return t().c(); } // not needed
    real e2() const { return t().e2(); }
    real k2() const { return t().k2(); }
    real kp2() const { return t().kp2(); }
    real k() const { return t().k(); }
    real kp() const { return t().kp(); }
    bool oblate() const { return t().oblate(); }
    bool prolate() const { return t().prolate(); }
    bool biaxial() const { return t().biaxial(); }
  public:
    Geodesic3(const Ellipsoid3& t = Ellipsoid3{});
    Geodesic3(real a, real b, real c);
    Geodesic3(real b, real e2, real k2, real kp2);
    const Ellipsoid3& t() const { return _t; }
    GeodesicLine3 Inverse(Angle bet1, Angle omg1, Angle bet2, Angle omg2,
                          real& s12, Angle& alp1, Angle& alp2) const;
    GeodesicLine3 Inverse(real bet1, real omg1, real bet2, real omg2,
                          real& s12, real& alp1, real& alp2) const;
    GeodesicLine3 Line(Angle bet1, Angle omg1, Angle alp1) const;
    GeodesicLine3 Line(real bet1, real omg1, real alp1) const;
    GeodesicLine3 Direct(Angle bet1, Angle omg1, Angle alp1, real s12,
                         Angle& bet2, Angle& omg2, Angle& alp2) const;
    GeodesicLine3 Direct(real bet1, real omg1, real alp1, real s12,
                         real& bet2, real& omg2, real& alp2) const;
    bool umbalt() const { return _umbalt; }
    void umbalt(bool numbalt) {
      if (_t.k2() > 0 && _t.kp2() > 0) _umbalt = numbalt;
    }
    bool biaxp() const { return _biaxp; }
    void biaxp(bool biaxp) { _biaxp = biaxp; }
    bool debug() const { return _debug; }
    void debug(bool debug) { _debug = debug; }
    bool hybridalt() const { return _hybridalt; }
    void hybridalt(bool hybridalt) { _hybridalt = hybridalt; }
    bool swapomg() const { return _swapomg; }
    void swapomg(bool swapomg) { _swapomg = swapomg; }
    // real ellipthresh() const { return _ellipthresh; }
    // void ellipthresh(real ellipthresh) { _ellipthresh = ellipthresh; }
  };

  } // namespace Triaxial
} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

// Include this because all the Geodesic3 methods return a GeodesicLine3.
#include <GeographicLib/Triaxial/GeodesicLine3.hpp>

#endif  // GEOGRAPHICLIB_GEODESIC3_HPP

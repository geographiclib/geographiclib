/**
 * \file Geodesic3.hpp
 * \brief Header for GeographicLib::Triaxial::Geodesic3 class
 *
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEODESIC3_HPP)
#define GEOGRAPHICLIB_GEODESIC3_HPP 1

#include <functional>
#include <memory>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs shared_ptr
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {
  namespace Triaxial {
  class GeodesicLine3;

  /**
   * \brief The solution of the geodesic problem for a triaxial ellipsoid
   *
   * This is an implementation of
   * <a href="https://books.google.com/books?id=RbwGAAAAYAAJ&pg=PA309">
   * Jacobi's method for finding geodesics on a triaxial ellipsoid</a>
   * (1839).  This class offers an interface to the solutions to the direct
   * geodesic problem via the GeodesicLine3 class which contains the meat of
   * Jacobi's direct solution.  In addition it provides a solution to the
   * inverse problem which closely parallels the solution for the biaxial
   * problem given by GeodesicExact.  For more details see \ref triaxial
   * and
   * - C. F. F. Karney,<br>
   *   <a href="https://arxiv.org/abs/2511.01621">
   *   Jacobi's solution for geodesics on a triaxial ellipsoid</a>,<br>
   *   Technical Report, SRI International, Nov. 2025.<br>
   *   <a href="https://arxiv.org/abs/2511.01621">arxiv:2511.01621</a>
   *
   * Data for testing the geodesic routines is available at
   * <a href="https://doi.org/10.5281/zenodo.12510796"> Test set of geodesics
   * on a triaxial ellipsoid (2024)</a>.
   *
   * \note There's a lot of new code here and testing is two orders of
   * magnitude more difficult compared with the biaxial case (an extra
   * parameter to fix the shape of the ellipsoid and geodesics now depend on
   * the longitude of the two end points separately).  I've limited by testing
   * to ellipsoids with \e a/\e c &le; 2.  I don't expect any problems if \e
   * a/\e c &le; 10; but you might run into problems with more eccentric
   * ellipsoids.  The code treats oblate and prolate (biaxial) ellipsoids
   * correctly; but, again, there may be problems with triaxial ellipsoids
   * which are \e extremely close to oblate or prolate i.e., with either \e k
   * or \e k' very small.  (However the triaxial model of the Earth where the
   * difference in the equatorial semiaxes is 70 m, \f$k' = 0.057\f$, is
   * treated just fine.)  While I have made every effort to ensure that the
   * code is error free, it's likely that some bugs remain.  Please use caution
   * with the results and report any problems (via email or with a Github
   * issue).
   *
   * Example of use:
   * \include example-Geodesic3.cpp
   *
   * <a href="Geod3Solve.1.html">Geod3Solve</a> is a command-line utility
   * providing access to the functionality of Geodesic3 and GeodesicLine3.
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT Geodesic3 {
  private:
    /// \cond SKIP
    // For access to BigValue, _ellipthresh, _biaxp
    friend class GeodesicLine3;
    /// \endcond
    using real = Math::real;
    using ang = Angle;
    static constexpr bool debug_ = false; // print out diagnostics
    static constexpr bool throw_ = true; // exception on convergence failure
    // special treatment for biaxial non-meridional
    static constexpr bool biaxp_ = true;
    // favor hybrid solution in terms of omg
    static constexpr bool  hybridalt_ = true;
    // allow swapping of omega{1,2}
    static constexpr bool swapomg_ = false;
    static constexpr int maxit_ = 200;
    Ellipsoid3 _t;

    // Run geodesic from bet1, omg1, alp1, find its first intersection with bet
    // = bet2a and return omg2a - omg2b
    real HybridA(ang bet1, ang omg1, ang alp1,
                 ang bet2a, ang omg2b, bool betp) const;
    static ang findroot(const std::function<real(const ang&)>& f,
                          ang xa,  ang xb,
                          real fa, real fb,
                          int* countn = nullptr, int* countb = nullptr);
    bool _umbalt;               // how coordinates wrap with umbilical lines
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
    /**
     * Constructor for a triaxial ellipsoid defined by Ellipsoid3 object.
     *
     * @param[in] t the Ellipsoid3 object.
     **********************************************************************/
    Geodesic3(const Ellipsoid3& t = Ellipsoid3{});
    /**
     * Constructor for a triaxial ellipsoid with semiaxes.
     *
     * @param[in] a the largest semiaxis.
     * @param[in] b the middle semiaxis.
     * @param[in] c the smallest semiaxis.
     * @exception GeographicErr if the required ordering is semiaxes is
     *   violated.
     *
     * The semiaxes must satisfy \e a &ge; \e b &ge; \e c &gt; 0.
     * If \e a = \e c (a sphere), then the oblate limit is taken.
     **********************************************************************/
    Geodesic3(real a, real b, real c);
    /**
     * Alternate constructor for a triaxial ellipsoid.
     *
     * @param[in] b the middle semiaxis.
     * @param[in] e2 the eccentricity squared \f$e^2 = (a^2 - c^2)/b^2\f$.
     * @param[in] k2 the oblateness parameter squared \f$k^2 = (b^2 - c^2) /
     *  (a^2 - c^2)\f$.
     * @param[in] kp2 the prolateness parameter squared \f$k'^2= (a^2 - b^2) /
     *   (a^2 - c^2)\f$.
     * @exception GeographicErr if the required ordering is semiaxes is
     *   violated.
     *
     * \note The constructor normalizes \e k2 and \e kp2 to ensure then \e k2 +
     * \e kp2 = 1.
     **********************************************************************/
    Geodesic3(real b, real e2, real k2, real kp2);
    /**
     * @return the Ellipsoid3 object for this projection.
     **********************************************************************/
    const Ellipsoid3& t() const { return _t; }
    /**
     * Solve the inverse problem
     *
     * @param[in] bet1 the ellipsoidal latitude of point 1.
     * @param[in] omg1 the ellipsoidal longitude of point 1.
     * @param[in] bet2 the ellipsoidal latitude of point 2.
     * @param[in] omg2 the ellipsoidal longitude of point 2.
     * @param[out] s12 the shortest distance between the points.
     * @param[out] alp1 the forward azimuth of the geodesic at point 1.
     * @param[out] alp2 the forward azimuth of the geodesic at point 2.
     * @return a GeodesicLine3 object defining the geodesic.
     **********************************************************************/
    GeodesicLine3 Inverse(Angle bet1, Angle omg1, Angle bet2, Angle omg2,
                          real& s12, Angle& alp1, Angle& alp2) const;
    /**
     * Solve the inverse problem in degrees
     *
     * @param[in] bet1 the ellipsoidal latitude of point 1 (in degrees).
     * @param[in] omg1 the ellipsoidal longitude of point 1 (in degrees).
     * @param[in] bet2 the ellipsoidal latitude of point 2 (in degrees).
     * @param[in] omg2 the ellipsoidal longitude of point 2 (in degrees).
     * @param[out] s12 the shortest distance between the points.
     * @param[out] alp1 the forward azimuth of the geodesic at point 1 (in
     *   degrees).
     * @param[out] alp2 the forward azimuth of the geodesic at point 2 (in
     *   degrees).
     * @return a GeodesicLine3 object defining the geodesic.
     **********************************************************************/
    GeodesicLine3 Inverse(real bet1, real omg1, real bet2, real omg2,
                          real& s12, real& alp1, real& alp2) const;
    /**
     * Solve the direct problem
     *
     * @param[in] bet1 the ellipsoidal latitude of point 1.
     * @param[in] omg1 the ellipsoidal longitude of point 1.
     * @param[in] alp1 the forward azimuth of the geodesic at point 1.
     * @param[in] s12 the distance from point 1 to point 2.
     * @param[out] bet2 the ellipsoidal latitude of point 2.
     * @param[out] omg2 the ellipsoidal longitude of point 2.
     * @param[out] alp2 the forward azimuth of the geodesic at point 2.
     * @return a GeodesicLine3 object defining the geodesic.
     **********************************************************************/
    GeodesicLine3 Direct(Angle bet1, Angle omg1, Angle alp1, real s12,
                         Angle& bet2, Angle& omg2, Angle& alp2) const;
    /**
     * Solve the direct problem in degrees
     *
     * @param[in] bet1 the ellipsoidal latitude of point 1 (in degrees).
     * @param[in] omg1 the ellipsoidal longitude of point 1 (in degrees).
     * @param[in] alp1 the forward azimuth of the geodesic at point 1 (in
     *   degrees).
     * @param[in] s12 the distance from point 1 to point 2.
     * @param[out] bet2 the ellipsoidal latitude of point 2 (in degrees).
     * @param[out] omg2 the ellipsoidal longitude of point 2 (in degrees).
     * @param[out] alp2 the forward azimuth of the geodesic at point 2 (in
     *   degrees).
     * @return a GeodesicLine3 object defining the geodesic.
     **********************************************************************/
    GeodesicLine3 Direct(real bet1, real omg1, real alp1, real s12,
                         real& bet2, real& omg2, real& alp2) const;
    /**
     * A geodesic line defined at a starting point
     *
     * @param[in] bet1 the ellipsoidal latitude of point 1.
     * @param[in] omg1 the ellipsoidal longitude of point 1.
     * @param[in] alp1 the forward azimuth of the geodesic at point 1.
     * @return a GeodesicLine3 object defining the geodesic.
     **********************************************************************/
    GeodesicLine3 Line(Angle bet1, Angle omg1, Angle alp1) const;
    /**
     * A geodesic line defined at a starting point specified in degrees
     *
     * @param[in] bet1 the ellipsoidal latitude of point 1 (in degrees).
     * @param[in] omg1 the ellipsoidal longitude of point 1 (in degrees).
     * @param[in] alp1 the forward azimuth of the geodesic at point 1 (in
     *   degrees).
     * @return a GeodesicLine3 object defining the geodesic.
     **********************************************************************/
    GeodesicLine3 Line(real bet1, real omg1, real alp1) const;
    /**
     * Specify the behavior of umbilical geodesics
     *
     * @param[in] numbalt the new value of the \e umbalt flag.
     *
     * If \e umbalt is true (resp. false) then the latitude (resp. longitude)
     * of an umbilical geodesic is the librating coordinate and the longitude
     * (resp. latitude) is the rotating coordinate.  This has no effect for
     * biaxial ellipsoids.  In this case the latitude (resp. longitude) is
     * always the librating coordinate for oblate (resp. prolate) ellipsoids.
     **********************************************************************/
    void umbalt(bool numbalt) {
      if (_t.k2() > 0 && _t.kp2() > 0) _umbalt = numbalt;
    }
    /**
     * @return the value of the \e umbalt flag; see umbalt(bool)).
     **********************************************************************/
    bool umbalt() const { return _umbalt; }
  };

  } // namespace Triaxial
} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

// Include this because all the Geodesic3 methods return a GeodesicLine3.
#include <GeographicLib/Triaxial/GeodesicLine3.hpp>

#endif  // GEOGRAPHICLIB_GEODESIC3_HPP

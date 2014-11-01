/**
 * \file Rhumb.hpp
 * \brief Header for GeographicLib::Rhumb and GeographicLib::RhumbLine classes
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_RHUMB_HPP)
#define GEOGRAPHICLIB_RHUMB_HPP 1

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Ellipsoid.hpp>

#if !defined(GEOGRAPHICLIB_RHUMBAREA_ORDER)
/**
 * The order of the series approximation used in rhumb area calculations.
 * GEOGRAPHICLIB_RHUMBAREA_ORDER can be set to any integer in [4, 8].
 **********************************************************************/
#  define GEOGRAPHICLIB_RHUMBAREA_ORDER \
  (GEOGRAPHICLIB_PRECISION == 2 ? 6 : \
   (GEOGRAPHICLIB_PRECISION == 1 ? 4 : 8))
#endif

namespace GeographicLib {

  class RhumbLine;

  /**
   * \brief Solve of the direct and inverse rhumb problems.
   *
   * The path of constant azimuth between two points on a ellipsoid at (\e
   * lat1, \e lon1) and (\e lat2, \e lon2) is called the rhumb line (also
   * called the loxodrome).  Its length is \e s12 and its azimuth is \e azi12
   * and \e azi2.  (The azimuth is the heading measured clockwise from north.)
   *
   * Given \e lat1, \e lon1, \e azi12, and \e s12, we can determine \e lat2,
   * and \e lon2.  This is the \e direct rhumb problem and its solution is
   * given by the function Rhumb::Direct.
   *
   * Given \e lat1, \e lon1, \e lat2, and \e lon2, we can determine \e azi12
   * and \e s12.  This is the \e inverse rhumb problem, whose solution is
   * given by Rhumb::Inverse.  This finds the shortest such rhumb line, i.e.,
   * the one that wraps no more than half way around the earth.
   *
   * These routines also optionally calculate the area under the rhumb line, \e
   * S12.  This is the area, measured counter-clockwise, of the rhumb line
   * quadrilateral with corners (<i>lat1</i>,<i>lon1</i>), (0,<i>lon1</i>),
   * (0,<i>lon2</i>), and (<i>lat2</i>,<i>lon2</i>).
   *
   * Note that rhumb lines may be appreciably longer (up to 50%) than the
   * corresponding Geodesic.  For example the distance between London Heathrow
   * and Tokyo Narita via the rhumb line is 11400 km which is 18% longer than
   * the geodesic distance 9600 km.
   *
   * For more information on rhumb lines see \ref rhumb.
   *
   * Example of use:
   * \include example-Rhumb.cpp
   **********************************************************************/

  class  GEOGRAPHICLIB_EXPORT Rhumb {
  private:
    typedef Math::real real;
    friend class RhumbLine;
    Ellipsoid _ell;
    bool _exact;
    real _c2;
    static const int tm_maxord = GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER;
    static const int maxpow_ = GEOGRAPHICLIB_RHUMBAREA_ORDER;
    // _R[0] unused
    real _R[maxpow_ + 1];
    static inline real overflow() {
      // Overflow value s.t. atan(overflow_) = pi/2
      static const real
        overflow = 1 / Math::sq(std::numeric_limits<real>::epsilon());
      return overflow;
    }
    static inline real tano(real x) {
      using std::abs; using std::tan;
      return
        2 * abs(x) == Math::pi() ? (x < 0 ? - overflow() : overflow()) :
        tan(x);
    }
    static inline real gd(real x)
    { using std::atan; using std::sinh; return atan(sinh(x)); }

    // Use divided differences to determine (mu2 - mu1) / (psi2 - psi1)
    // accurately
    //
    // Definition: Df(x,y,d) = (f(x) - f(y)) / (x - y)
    // See:
    //   W. M. Kahan and R. J. Fateman,
    //   Symbolic computation of divided differences,
    //   SIGSAM Bull. 33(3), 7-28 (1999)
    //   http://dx.doi.org/10.1145/334714.334716
    //   http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf

    static inline real Dlog(real x, real y) {
      real t = x - y;
      return t ? 2 * Math::atanh(t / (x + y)) / t : 1 / x;
    }
    static inline real Dtan(real x, real y) {
      real d = x - y, tx = tano(x), ty = tano(y), txy = tx * ty;
      return d ? (2 * txy > -1 ? (1 + txy) * tano(d) : tx - ty) / d :
        1 + txy;
    }
    static inline real Datan(real x, real y) {
      using std::atan;
      real d = x - y, xy = x * y;
      return d ? (2 * xy > -1 ? atan( d / (1 + xy) ) : atan(x) - atan(y)) / d :
        1 / (1 + xy);
    }
    static inline real Dsin(real x, real y) {
      using std::sin; using std::cos;
      real d = (x - y) / 2;
      return cos((x + y)/2) * (d ? sin(d) / d : 1);
    }
    static inline real Dcos(real x, real y) {
      using std::sin;
      real d = (x - y) / 2;
      return - sin((x + y) / 2) * (d ? sin(d) / d : 1);
    }
    static inline real Dsinh(real x, real y) {
      using std::sinh; using std::cosh;
      real d = (x - y) / 2;
      return cosh((x + y) / 2) * (d ? sinh(d) / d : 1);
    }
    // Dcosh(x,y) = sinh((x+y)/2) * sinh((x-y)/2)/((x-y)/2)
    static inline real Dcosh(real x, real y) {
      using std::sinh;
      real d = (x - y) / 2;
      return sinh((x + y) / 2) * (d ? sinh(d) / d : 1);
    }
    static inline real Dasinh(real x, real y) {
      real d = x - y,
        hx = Math::hypot(real(1), x), hy = Math::hypot(real(1), y);
      return d ? Math::asinh(x*y > 0 ? d * (x + y) / (x*hy + y*hx) :
                             x*hy - y*hx) / d :
        1 / hx;
    }
    static inline real Dgd(real x, real y) {
      using std::sinh;
      return Datan(sinh(x), sinh(y)) * Dsinh(x, y);
    }
    static inline real Dgdinv(real x, real y) {
      return Dasinh(tano(x), tano(y)) * Dtan(x, y);
    }
    // Copied from LambertConformalConic...
    // e * atanh(e * x) = log( ((1 + e*x)/(1 - e*x))^(e/2) ) if f >= 0
    // - sqrt(-e2) * atan( sqrt(-e2) * x)                    if f < 0
    inline real eatanhe(real x) const {
      using std::atan;
      return _ell._f >= 0 ? _ell._e * Math::atanh(_ell._e * x) :
        - _ell._e * atan(_ell._e * x);
    }
    // Copied from LambertConformalConic...
    // Deatanhe(x,y) = eatanhe((x-y)/(1-e^2*x*y))/(x-y)
    inline real Deatanhe(real x, real y) const {
      real t = x - y, d = 1 - _ell._e2 * x * y;
      return t ? eatanhe(t / d) / t : _ell._e2 / d;
    }
    // (E(x) - E(y)) / (x - y) -- E = incomplete elliptic integral of 2nd kind
    real DE(real x, real y) const;
    // (mux - muy) / (phix - phiy) using elliptic integrals
    real DRectifying(real latx, real laty) const;
    // (psix - psiy) / (phix - phiy)
    real DIsometric(real latx, real laty) const;

    // (sum(c[j]*sin(2*j*x),j=1..n) - sum(c[j]*sin(2*j*x),j=1..n)) / (x - y)
    static real SinCosSeries(bool sinp,
                             real x, real y, const real c[], int n);
    // (mux - muy) / (chix - chiy) using Krueger's series
    real DConformalToRectifying(real chix, real chiy) const;
    // (chix - chiy) / (mux - muy) using Krueger's series
    real DRectifyingToConformal(real mux, real muy) const;

    // (mux - muy) / (psix - psiy)
    real DIsometricToRectifying(real psix, real psiy) const;
    // (psix - psiy) / (mux - muy)
    real DRectifyingToIsometric(real mux, real muy) const;

    real MeanSinXi(real psi1, real psi2) const;
  public:

    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters).
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.  If \e f &gt; 1, set
     *   flattening to 1/\e f.
     * @param[in] exact if true (the default) use an addition theorem for
     *   elliptic integrals to compute divided differences; otherwise use
     *   series expansion (accurate for |<i>f</i>| < 0.01).
     * @exception GeographicErr if \e a or (1 &minus; \e f) \e a is not
     *   positive.
     *
     * See \ref rhumb, for a detailed description of the \e exact parameter.
     **********************************************************************/
    Rhumb(real a, real f, bool exact = true);

    /**
     * Solve the direct rhumb problem returning also the area.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi12 azimuth of the rhumb line (degrees).
     * @param[in] s12 distance between point 1 and point 2 (meters); it can be
     *   negative.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees).
     * @param[out] S12 the area under the rhumb line (meters<sup>2</sup>).
     *
     * \e lat1 should be in the range [&minus;90&deg;, 90&deg;]; \e lon1 and \e
     * azi1 should be in the range [&minus;540&deg;, 540&deg;).  The values of
     * \e lon2 and \e azi2 returned are in the range [&minus;180&deg;,
     * 180&deg;).
     *
     * If point 1 is a pole, the cosine of its latitude is taken to be
     * 1/&epsilon;<sup>2</sup> (where &epsilon; is 2<sup>-52</sup>).  This
     * position, which is extremely close to the actual pole, allows the
     * calculation to be carried out in finite terms.  If \e s12 is large
     * enough that the rhumb line crosses a pole, the longitude of point 2
     * is indeterminate (a NaN is returned for \e lon2).
     **********************************************************************/
    void Direct(real lat1, real lon1, real azi12, real s12,
                real& lat2, real& lon2, real& S12) const {
      GenDirect(lat1, lon1, azi12, s12, true, lat2, lon2, S12);
    }

    /**
     * Solve the direct rhumb problem without the area.
     **********************************************************************/
    void Direct(real lat1, real lon1, real azi12, real s12,
                real& lat2, real& lon2) const {
      real dummy;
      GenDirect(lat1, lon1, azi12, s12, false, lat2, lon2, dummy);
    }

    /**
     * The general direct rhumb problem.  Rhumb::Direct is defined in terms
     * of this function.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi12 azimuth of the rhumb line (degrees).
     * @param[in] s12 distance between point 1 and point 2 (meters); it can be
     *   negative.
     * @param[in] areap should S12 be computed?
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees).
     * @param[out] S12 if \e areap, then the area under the rhumb line
     *   (meters<sup>2</sup>); otherwise its value is unchanged.
     **********************************************************************/
    void GenDirect(real lat1, real lon1, real azi12, real s12, bool areap,
                   real& lat2, real& lon2, real& S12) const;

    /**
     * Solve the inverse rhumb problem returning also the area.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] lat2 latitude of point 2 (degrees).
     * @param[in] lon2 longitude of point 2 (degrees).
     * @param[out] s12 rhumb distance between point 1 and point 2 (meters).
     * @param[out] azi12 azimuth of the rhumb line (degrees).
     * @param[out] S12 the area under the rhumb line (meters<sup>2</sup>).
     *
     * The shortest rhumb line is found.  \e lat1 and \e lat2 should be in the
     * range [&minus;90&deg;, 90&deg;]; \e lon1 and \e lon2 should be in the
     * range [&minus;540&deg;, 540&deg;).  The value of \e azi12 returned is in
     * the range [&minus;180&deg;, 180&deg;).
     *
     * If either point is a pole, the cosine of its latitude is taken to be
     * 1/&epsilon;<sup>2</sup> (where &epsilon; is 2<sup>-52</sup>).  This
     * position, which is extremely close to the actual pole, allows the
     * calculation to be carried out in finite terms.
     **********************************************************************/
    void Inverse(real lat1, real lon1, real lat2, real lon2,
                 real& s12, real& azi12, real& S12) const {
      GenInverse(lat1, lon1, lat2, lon2, true, s12, azi12, S12);
    }

    /**
     * Solve the inverse rhumb problem without the area.
     **********************************************************************/
    void Inverse(real lat1, real lon1, real lat2, real lon2,
                 real& s12, real& azi12) const {
      real dummy;
      GenInverse(lat1, lon1, lat2, lon2, false, s12, azi12, dummy);
    }

    /**
     * The general inverse rhumb problem.  Rhumb::Inverse is defined in terms
     * of this function.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] lat2 latitude of point 2 (degrees).
     * @param[in] lon2 longitude of point 2 (degrees).
     * @param[in] areap should S12 be computed?
     * @param[out] s12 rhumb distance between point 1 and point 2 (meters).
     * @param[out] azi12 azimuth of the rhumb line (degrees).
     * @param[out] S12 if \e areap, then the area under the rhumb line
     *   (meters<sup>2</sup>); otherwise its value is unchanged.
     **********************************************************************/
    void GenInverse(real lat1, real lon1, real lat2, real lon2, bool areap,
                    real& s12, real& azi12, real& S12) const;

    /**
     * Set up to compute several points on a single rhumb line.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi12 azimuth of the rhumb line (degrees).
     * @return a RhumbLine object.
     *
     * \e lat1 should be in the range [&minus;90&deg;, 90&deg;]; \e lon1 and \e
     * azi12 should be in the range [&minus;540&deg;, 540&deg;).
     *
     * If point 1 is a pole, the cosine of its latitude is taken to be
     * 1/&epsilon;<sup>2</sup> (where &epsilon; is 2<sup>-52</sup>).  This
     * position, which is extremely close to the actual pole, allows the
     * calculation to be carried out in finite terms.
     **********************************************************************/
    RhumbLine Line(real lat1, real lon1, real azi12) const;

    /** \name Inspector functions.
     **********************************************************************/
    ///@{

    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const { return _ell.MajorRadius(); }

    /**
     * @return \e f the  flattening of the ellipsoid.  This is the
     *   value used in the constructor.
     **********************************************************************/
    Math::real Flattening() const { return _ell.Flattening(); }

    /**
     * A global instantiation of Rhumb with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    static const Rhumb& WGS84();
  };

  /**
   * \brief Find a sequence of points on a single rhumb line.
   *
   * RhumbLine facilitates the determination of a series of points on a single
   * rhumb line.  The starting point (\e lat1, \e lon1) and the azimuth \e
   * azi12 are specified in the call to Rhumb::Line which returns a RhumbLine
   * object.  RhumbLine.Position returns the location of point 2 (and,
   * optionally, the corresponding area, \e S12) a distance \e s12 along the
   * rhumb line.
   *
   * There is no public constructor for this class.  (Use Rhumb::Line to create
   * an instance.)  The Rhumb object used to create a RhumbLine must stay in
   * scope as long as the RhumbLine.
   *
   * Example of use:
   * \include example-RhumbLine.cpp
   **********************************************************************/

  class  GEOGRAPHICLIB_EXPORT RhumbLine {
  private:
    typedef Math::real real;
    friend class Rhumb;
    const Rhumb& _rh;
    bool _exact;
    real _lat1, _lon1, _azi12, _salp, _calp, _mu1, _psi1, _r1;
    RhumbLine& operator=(const RhumbLine&); // copy assignment not allowed
    RhumbLine(const Rhumb& rh, real lat1, real lon1, real azi12,
              bool exact);
  public:
    /**
     * Compute the position of point 2 which is a distance \e s12 (meters) from
     * point 1.  The area is also computed.
     *
     * @param[in] s12 distance between point 1 and point 2 (meters); it can be
     *   negative.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees).
     * @param[out] S12 the area under the rhumb line (meters<sup>2</sup>).
     *
     * The values of \e lon2 and \e azi2 returned are in the range
     * [&minus;180&deg;, 180&deg;).
     *
     * If \e s12 is large enough that the rhumb line crosses a pole, the
     * longitude of point 2 is indeterminate (a NaN is returned for \e lon2).
     **********************************************************************/
    void Position(real s12, real& lat2, real& lon2, real& S12) const {
      GenPosition(s12, true, lat2, lon2, S12);
    }

    /**
     * Compute the position of point 2 which is a distance \e s12 (meters) from
     * point 1.  The area is not computed.
     **********************************************************************/
    void Position(real s12, real& lat2, real& lon2) const {
      real dummy;
      GenPosition(s12, false, lat2, lon2, dummy);
    }

    /**
     * The general position routine.  RhumbLine::Position is defined in term so
     * this function.
     *
     * @param[in] s12 distance between point 1 and point 2 (meters); it can be
     *   negative.
     * @param[in] areap should S12 be computed?
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees).
     * @param[out] S12 if \e areap, then the area under the rhumb line
     *   (meters<sup>2</sup>); otherwise its value is unchanged.
     *
     * The values of \e lon2 and \e azi2 returned are in the range
     * [&minus;180&deg;, 180&deg;).
     *
     * If \e s12 is large enough that the rhumb line crosses a pole, the
     * longitude of point 2 is indeterminate (a NaN is returned for \e lon2).
     **********************************************************************/
    void GenPosition(real s12, bool areap,
                     real& lat2, real& lon2, real& S12) const;

    /** \name Inspector functions
     **********************************************************************/
    ///@{

    /**
     * @return \e lat1 the latitude of point 1 (degrees).
     **********************************************************************/
    Math::real Latitude() const { return _lat1; }

    /**
     * @return \e lon1 the longitude of point 1 (degrees).
     **********************************************************************/
    Math::real Longitude() const { return _lon1; }

    /**
     * @return \e azi12 the azimuth of the rhumb line (degrees).
     **********************************************************************/
    Math::real Azimuth() const { return  _azi12; }

    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value inherited from the Rhumb object used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const { return _rh.MajorRadius(); }

    /**
     * @return \e f the flattening of the ellipsoid.  This is the value
     *   inherited from the Rhumb object used in the constructor.
     **********************************************************************/
    Math::real Flattening() const { return _rh.Flattening(); }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_RHUMB_HPP

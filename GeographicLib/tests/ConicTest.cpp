/**
 * \file ProjTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/LambertConformalConic.hpp"
#include "GeographicLib/Constants.hpp"
#include "GeographicLib/DMS.hpp"
#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace GeographicLib {

  /**
   * \brief Lambert Conformal Conic Projection
   *
   * Implementation taken from the report,
   * - J. P. Snyder,
   *   <a href="http://pubs.er.usgs.gov/usgspubs/pp/pp1395"> Map Projections: A
   *   Working Manual</a>, USGS Professional Paper 1395 (1987),
   *   pp. 107&ndash;109.
   *
   * This is a straightforward implementation of the equations in Snyder except
   * that Newton's method is used to invert the projection and that several of
   * the formulas are modified so that the projection correctly degenerates to
   * the Mercator projection or the polar sterographic projection.
   *
   * The ellipsoid parameters, the standard parallels, and the scale on the
   * standard parallels are set in the constructor.  If the standard parallels
   * are both at a single pole the projection becomes the polar stereographic
   * projection (compare with the PolarStereographic class).  If the standard
   * parallels are symmetric around equator, the projection becomes the
   * Mercator projection.  The central meridian (which is a trivial shift of
   * the longitude) is specified as the \e lon0 argument of the
   * LambertConformalConic::Forward and LambertConformalConic::Reverse
   * functions.  The latitude of origin is taken to be the latitude of tangency
   * which lies between the standard parallels which is given by
   * LambertConformalConic::OriginLatitude.  There is no provision in this
   * class for specifying a false easting or false northing or a different
   * latitude of origin.  However these are can be simply included by the
   * calling function.  For example the Pennsylvania South state coordinate
   * system (<a href="http://www.spatialreference.org/ref/epsg/3364/">
   * EPSG:3364</a>) is obtained by:
   \code
   const double
     a = GeographicLib::Constants::WGS84_a(), r = 298.257222101, // GRS80
     lat1 = 39 + 56/60.0, lat1 = 40 + 58/60.0, // standard parallels
     k1 = 1,                                   // scale
     lat0 = 39 + 20/60.0, lon0 = 75 + 45/60.0, // origin
     fe = 600000, fn = 0;                      // false easting and northing
   // Set up basic projection
   const GeographicLib::LambertConformalConic PASouth(a, r, lat1, lat2, k1);
   double x0, y0;
   {
     // Transform origin point
     PASouth.Forward(lon0, lat0, lon0, x0, y0);
     x0 -= fe; y0 -= fn;         // Combine result with false origin
   }
   double lat, lon, x, y;
   // Sample conversion from geodetic to PASouth grid
   std::cin >> lat >> lon;
   PASouth.Forward(lon0, lat, lon, x, y);
   x -= x0; y -= y0;
   std::cout << x << " " << y << "\n";
   // Sample conversion from PASouth grid to geodetic
   std::cin >> x >> y;
   x += x0; y += y0;
   PASouth.Reverse(lon0, x, y, lat, lon);
   std::cout << lat << " " << lon << "\n";
   \endcode
   **********************************************************************/
  class OLambertConformalConic {
  private:
    typedef Math::real real;
    const real _a, _r, _f, _e2, _e, _e2m;
    real _sign, _n, _nc, _lt0, _t0n, _t0nm1, _scale, _lat0, _k0;
    static const real eps, eps2, epsx, tol, ahypover;
    static const int numit = 5;
    static inline real sq(real x) throw() { return x * x; }
    // e * atanh(e * x) = log( ((1 + e*x)/(1 - e*x))^(e/2) ) if f >= 0
    // - sqrt(-e2) * atan( sqrt(-e2) * x)                    if f < 0
    inline real eatanhe(real x) const throw() {
      return _f >= 0 ? _e * Math::atanh(_e * x) : - _e * std::atan(_e * x);
    }
    inline real mf(real sphi, real cphi) const throw() {
      return cphi/std::sqrt(1 - _e2 * sq(sphi)); // Snyder's m, p 108, eq 14-15
    }
    inline real tf(real sphi, real cphi) const throw() {
      // Snyder's t, p 108, eq 15-9a
      // First factor is sqrt((1 - sphi) / (1 + sphi))
      // Second factor is ((1 + e * sphi)/(1 - e * sphi)) ^ (e/2)
      return (sphi >= 0 ? cphi / (1 + sphi) : (1 - sphi) / cphi) *
        std::exp(eatanhe(sphi));
    }
    inline real logtf(real sphi, real cphi) const throw() {
      // Return log(t).  This retains relative accuracy near the equator where
      // t -> 1.  This is the negative of the standard Mercator northing.  Note
      // that log( sqrt((1 - sin(phi)) / (1 + sin(phi))) ) = -asinh(tan(phi))
      return -Math::asinh(sphi/std::max(epsx, cphi)) + eatanhe(sphi);
    }
    inline real logmtf(real sphi) const throw() {
      // Return log(m/t).  This allows the cancellation of cphi in m and t.
      return std::log((1 + sphi)/std::sqrt(1 - _e2 * sq(sphi))) -
        eatanhe(sphi);
    }
    void Init(real sphi1, real cphi1, real sphi2, real cphi2, real k1) throw();
  public:

    /**
     * Constructor with a single standard parallel.
     *
     * @param[in] a equatorial radius of ellipsoid (meters)
     * @param[in] r reciprocal flattening of ellipsoid.  Setting \e r = 0
     *   implies \e r = inf or flattening = 0 (i.e., a sphere).  Negative \e r
     *   indicates a prolate ellipsoid.
     * @param[in] stdlat standard parallel (degrees), the circle of tangency.
     * @param[in] k0 scale on standard parallel.
     *
     * An exception is thrown if \e a or \e k0 is not positive or if \e stdlat
     * is not in the range [-90, 90].
     **********************************************************************/
    OLambertConformalConic(real a, real r, real stdlat, real k0);

    /**
     * Constructor with two standard parallels.
     *
     * @param[in] a equatorial radius of ellipsoid (meters)
     * @param[in] r reciprocal flattening of ellipsoid.  Setting \e r = 0
     *   implies \e r = inf or flattening = 0 (i.e., a sphere).  Negative \e r
     *   indicates a prolate ellipsoid.
     * @param[in] stdlat1 first standard parallel (degrees).
     * @param[in] stdlat2 second standard parallel (degrees).
     * @param[in] k1 scale on the standard parallels.
     *
     * An exception is thrown if \e a or \e k0 is not positive or if \e stdlat1
     * or \e stdlat2 is not in the range [-90, 90].  In addition, if either \e
     * stdlat1 or \e stdlat2 is a pole, then an exception is thrown if \e
     * stdlat1 is not equal \e stdlat2
     **********************************************************************/
    OLambertConformalConic(real a, real r, real stdlat1, real stdlat2, real k1);

    /**
     * Constructor with two standard parallels specified by sines and cosines.
     *
     * @param[in] a equatorial radius of ellipsoid (meters)
     * @param[in] r reciprocal flattening of ellipsoid.  Setting \e r = 0
     *   implies \e r = inf or flattening = 0 (i.e., a sphere).  Negative \e r
     *   indicates a prolate ellipsoid.
     * @param[in] sinlat1 sine of first standard parallel.
     * @param[in] coslat1 cosine of first standard parallel.
     * @param[in] sinlat2 sine of second standard parallel.
     * @param[in] coslat2 cosine of second standard parallel.
     * @param[in] k1 scale on the standard parallels.
     *
     * This allows parallels close to the poles to be specified accurately.
     **********************************************************************/
    OLambertConformalConic(real a, real r,
                          real sinlat1, real coslat1,
                          real sinlat2, real coslat2,
                          real k1);

    /**
     * Set the scale for the projection.
     *
     * @param[in] lat (degrees).
     * @param[in] k scale at latitude \e lat (default 1).
     *
     * This allows a "latitude of true scale" to be specified.  An exception is
     * thrown if \e k is not positive.
     **********************************************************************/
    void SetScale(real lat, real k = real(1));

    /**
     * Forward projection, from geographic to Lambert conformal conic.
     *
     * @param[in] lon0 central meridian longitude (degrees).
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * The latitude origin is given by LambertConformalConic::LatitudeOrigin().
     * No false easting or northing is added and \e lat should be in the range
     * [-90, 90]; \e lon and \e lon0 should be in the range [-180, 360].  The
     * values of \e x and \e y returned for points which project to infinity
     * (i.e., one or both of the poles) will be large but finite.  The value of
     * \e k returned for the poles in infinite (unless \e lat equals the
     * latitude origin).
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y, real& gamma, real& k) const throw();

    /**
     * Reverse projection, from Lambert conformal conic to geographic.
     *
     * @param[in] lon0 central meridian longitude (degrees).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * The latitude origin is given by LambertConformalConic::LatitudeOrigin().
     * No false easting or northing is added.  \e lon0 should be in the range
     * [-180, 360].  The value of \e lon returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon, real& gamma, real& k) const throw();

    /**
     * LambertConformalConic::Forward without returning the convergence and
     * scale.
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y) const throw() {
      real gamma, k;
      Forward(lon0, lat, lon, x, y, gamma, k);
    }

    /**
     * LambertConformalConic::Reverse without returning the convergence and
     * scale.
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon) const throw() {
      real gamma, k;
      Reverse(lon0, x, y, lat, lon, gamma, k);
    }

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _a; }

    /**
     * @return \e r the inverse flattening of the ellipsoid.  This is the
     *   value used in the constructor.  A value of 0 is returned for a sphere
     *   (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return _r; }

    /**
     * @return latitude of the origin for the projection (degrees).
     *
     * This is the latitude of minimum scale and equals the \e stdlat in the
     * 1-parallel constructor and lies between \e stdlat1 and \e stdlat2 in the
     * 2-parallel constructors.
     **********************************************************************/
    Math::real OriginLatitude() const throw() { return _lat0; }

    /**
     * @return central scale for the projection.  This is the scale on the
     *   latitude of origin..
     **********************************************************************/
    Math::real CentralScale() const throw() { return _k0; }
    ///@}

    /**
     * A global instantiation of LambertConformalConic with the WGS84 ellipsoid
     * and the UPS scale factor and \e stdlat = 0.  This degenerates to the
     * Mercator projection.
     **********************************************************************/
    const static OLambertConformalConic Mercator;
  };

} // namespace GeographicLib

/**
 * \file LambertConformalConic.cpp
 * \brief Implementation for GeographicLib::LambertConformalConic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

namespace GeographicLib {

  using namespace std;

  const Math::real OLambertConformalConic::eps = numeric_limits<real>::epsilon();
  const Math::real OLambertConformalConic::eps2 = sqrt(eps);
  const Math::real OLambertConformalConic::epsx = sq(eps);
  const Math::real OLambertConformalConic::tol = real(0.1) * eps2;
  // Large enough value so that atan(sinh(ahypover)) = pi/2
  const Math::real OLambertConformalConic::ahypover =
    real(numeric_limits<real>::digits) * log(real(numeric_limits<real>::radix))
    + 2;

  OLambertConformalConic::OLambertConformalConic(real a, real r,
                                               real stdlat, real k0)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_f < 1))
      throw GeographicErr("Minor radius is not positive");
    if (!(k0 > 0))
      throw GeographicErr("Scale is not positive");
    if (!(abs(stdlat) <= 90))
      throw GeographicErr("Standard latitude not in [-90, 90]");
    real
      phi = stdlat * Constants::degree(),
      sphi = sin(phi),
      cphi = abs(stdlat) != 90 ? cos(phi) : 0;
    Init(sphi, cphi, sphi, cphi, k0);
  }

  OLambertConformalConic::OLambertConformalConic(real a, real r,
                                               real stdlat1, real stdlat2,
                                               real k1)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_f < 1))
      throw GeographicErr("Minor radius is not positive");
    if (!(k1 > 0))
      throw GeographicErr("Scale is not positive");
    if (!(abs(stdlat1) <= 90))
      throw GeographicErr("Standard latitude 1 not in [-90, 90]");
    if (!(abs(stdlat2) <= 90))
      throw GeographicErr("Standard latitude 2 not in [-90, 90]");
    if (abs(stdlat1) == 90 || abs(stdlat2) == 90)
      if (!(stdlat1 == stdlat2))
        throw GeographicErr
          ("Standard latitudes must be equal is either is a pole");
    real
      phi1 = stdlat1 * Constants::degree(),
      phi2 = stdlat2 * Constants::degree();
    Init(sin(phi1), abs(stdlat1) != 90 ? cos(phi1) : 0,
         sin(phi2), abs(stdlat2) != 90 ? cos(phi2) : 0, k1);
  }

  OLambertConformalConic::OLambertConformalConic(real a, real r,
                                               real sinlat1, real coslat1,
                                               real sinlat2, real coslat2,
                                               real k1)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_f < 1))
      throw GeographicErr("Minor radius is not positive");
    if (!(k1 > 0))
      throw GeographicErr("Scale is not positive");
    if (coslat1 == 0 || coslat2 == 0)
      if (!(coslat1 == coslat2 && sinlat1 == sinlat2))
        throw GeographicErr
          ("Standard latitudes must be equal is either is a pole");
    Init(sinlat1, coslat1, sinlat2, coslat2, k1);
  }

  void OLambertConformalConic::Init(real sphi1, real cphi1,
                                   real sphi2, real cphi2, real k1) throw() {
    // Determine hemisphere of tangent latitude
    _sign = sphi1 + sphi2 > 0 ? 1 : sphi1 + sphi2 < 0  ? -1 :
      atan2(sphi1, cphi1) + atan2(sphi2, cphi2) >= 0 ? 1 : -1;
    // Internally work with tangent latitude positive
    sphi1 *= _sign;
    sphi2 *= _sign;
    real
      m1 = mf(sphi1, cphi1), lt1 = logtf(sphi1, cphi1),
      m2 = mf(sphi2, cphi2), lt2 = logtf(sphi2, cphi2),
      sindiff = abs(sphi1 * cphi2 - cphi1 * sphi2);
    real sphi0, cphi0;          // Use phi0 = tangent latitude
    if (sindiff > eps2) {
      // sphi0 = Snyder's n, p 108, eq 15-8
      sphi0 = (log(m1) - log (m2)) / (lt1 - lt2);
      cphi0 = sindiff < real(0.7) ? sqrt(1 - sq(sphi0)) :
        // cos(phi0) = sqrt( - (sin(phi0) - 1) * (sin(phi0) + 1) )
        sqrt( -(logmtf(sphi1) - logmtf(sphi2)) / (lt1 - lt2)
              * (sphi0 + 1));
    } else {
      // Set phi0 = (phi1 + phi2)/2
      sphi0 = (sphi1 * sqrt((1 + cphi2)/(1 + cphi1)) +
               sphi2 * sqrt((1 + cphi1)/(1 + cphi2)))/2;
      cphi0 = (cphi1 * sqrt((1 + sphi2)/(1 + sphi1)) +
               cphi2 * sqrt((1 + sphi1)/(1 + sphi2)))/2;
    }
    _n = sphi0;                 // Snyder's n
    _nc = cphi0;                // sqrt(1 - sq(n))
    _lat0 = atan2(_sign * sphi0, cphi0) / Constants::degree();
    _lt0 = logtf(sphi0, cphi0); // Snyder's log(t0)
    _t0n = exp(_n * _lt0);      // Snyder's t0^n
    _t0nm1 = Math::expm1(_n * _lt0);      // Snyder's t0^n - 1
    // k1 * m1/t1^n = k1 * m2/t2^n = k1 * n * (Snyder's F)
    _scale = k1 *
      (_nc == 0 ?               // if phi0 = pi/2, n = t = m = 0, so take limit
       2/( sqrt(_e2m) * exp(eatanhe(real(1))) ) :
       (lt1 > lt2 ? m1 / exp(_n * lt1) : m2 / exp(_n * lt2)));
    // Scale at phi0
    _k0 = _scale * _t0n / mf(sphi0, cphi0);
  }

  const OLambertConformalConic
  OLambertConformalConic::Mercator(Constants::WGS84_a(), Constants::WGS84_r(),
                                  real(0), real(1));

  void OLambertConformalConic::Forward(real lon0, real lat, real lon,
                                      real& x, real& y, real& gamma, real& k)
    const throw() {
    if (lon - lon0 > 180)
      lon -= lon0 - 360;
    else if (lon - lon0 <= -180)
      lon -= lon0 + 360;
    else
      lon -= lon0;
    lat *= _sign;
    real
      phi = lat * Constants::degree(),
      lam = lon * Constants::degree(),
      sphi = sin(phi), cphi = abs(lat) != 90 ? cos(phi) : 0,
      m = mf(sphi, cphi),
      lt = logtf(sphi, cphi),
      tn = exp(_n * lt),
      theta = _n * lam, stheta = sin(theta), ctheta = cos(theta);
    x = _a * _scale * tn * (_n > eps2 ? stheta/_n : lam);
    if (_n > real(0.5))
      y = (_t0n - tn * ctheta)/_n;
    else {
      // write as
      // _t0n * ( (1 - ctheta)/_n - (tn/_t0n - 1)/_n * ctheta )
      // _t0n * ( stheta^2/(1 + ctheta) / _n
      //          - ((t/_t0)^n - 1)/_n * ctheta )
      // _t0n * ( stheta^2/(1 + ctheta) / _n
      //          - (exp(_n * log(t/_t0)) - 1)/_n * ctheta )
      // _t0n * (ax - bx)
      real
        ax = (_n > eps2 ? sq(stheta) / _n : lam * theta) / (1 + ctheta),
        fx = lt - _lt0,
        bx = ctheta * (abs(_n * fx) > eps ? Math::expm1(_n * fx) / _n : fx);
      y = _t0n * (ax - bx);
    }
    y *= _a * _scale * _sign;
    k = _scale * (m == 0 && _nc == 0 && sphi > 0 ?
                  sqrt(_e2m) * exp(eatanhe(real(1))) / 2 :
                  tn/m);        // infinite if pole and _n < 1
    gamma = theta * _sign;
  }

  void OLambertConformalConic::Reverse(real lon0, real x, real y,
                                      real& lat, real& lon,
                                      real& gamma, real& k)
    const throw() {
    y *= _sign;
    x /= _a * _scale;
    y /= _a * _scale;
    real
      /*
        rho = hypot(x, rho0 - y)
        rho0 = _a * _scale * _t0n / _n
        _n * rho = hypot(x * _n, _a * _scale * _t0n - y * _n)
        tn    = _n * rho / (_a * _scale)
              = hypot(_n * x/(_a * _scale), _t0n - _n * y/(_a * _scale))
        theta = atan2(_n * x/(_a * _scale), _t0n - _n * y/(_a * _scale))
        lam = theta/_n
      */
      x1 = _n * x, y1 = _n * y,
      tn = Math::hypot(x1, _t0n - y1), theta = atan2(x1, _t0n - y1),
      lam = _n != 0 ? theta/_n : x/_t0n;
    real q;                     // -log(t)
    if (_n != 0 && tn < real(0.5))
      q = -log(tn) / _n;
    else {
      real
        b1 = _t0nm1 - _n * y,
        tnm1 = (sq(x1) +  2 * b1 + sq(b1)) / (tn + 1); // tn - 1
      q = _n != 0 ? - Math::log1p(tnm1)/_n :
        -2 * (_lt0  - y) / (tn + 1); // _t0nm1 -> _n * _lt0
    }
    // Clip to [-ahypover, ahypover] to avoid overflow later
    q = q < ahypover ? (q > -ahypover ? q : -ahypover) : ahypover;
    // Recast Snyder's 15-9 as
    // q = q' - eatanh(tanh(q'))
    // where q = -log(t) and q' = asinh(tan(phi))
    // q is known and we solve for q' by Newton's method.
    // Write f(q') = q' - eatanh(tanh(q')) - q
    // f'(q') = (1 - e^2)/(1 - e^2 * tanh(q')^2)
    // Starting guess is q' = q.
    real qp = q;
    for (int i = 0; i < numit; ++i) {
      real
        t = tanh(qp),
        dqp = -(qp - eatanhe(t) - q) * (1 - _e2 * sq(t)) / _e2m;
      qp += dqp;
      if (abs(dqp) < tol)
        break;
    }
    double
      phi = _sign * atan(sinh(qp)),
      m = mf(tanh(qp), 1/cosh(qp));
    lat = phi / Constants::degree();
    lon = lam / Constants::degree();
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon + lon0 >= 180)
      lon += lon0 - 360;
    else if (lon + lon0 < -180)
      lon += lon0 + 360;
    else
      lon += lon0;
    gamma = _sign * theta / Constants::degree();
    k = _scale * (m == 0 && _nc == 0 && qp > 0 ?
                  sqrt(_e2m) * exp(eatanhe(real(1))) / 2 :
                  tn/m);        // infinite if pole and _n < 1
  }

  void OLambertConformalConic::SetScale(real lat, real k) {
    if (!(k > 0))
      throw GeographicErr("Scale is not positive");
    real x, y, gamma, kold;
    Forward(0, lat, 0, x, y, gamma, kold);
    k /= kold;
    _scale *= k;
    _k0 *= k;
  }

} // namespace GeographicLib





GeographicLib::Math::extended
dist(GeographicLib::Math::extended a, GeographicLib::Math::extended r,
     GeographicLib::Math::extended lat0, GeographicLib::Math::extended lon0,
     GeographicLib::Math::extended lat1, GeographicLib::Math::extended lon1) {
  using namespace GeographicLib;
  typedef Math::extended extended;
  extended
    phi = lat0 * Constants::degree(),
    f = r != 0 ? 1/r : 0,
    e2 = f * (2 - f),
    sinphi = sin(phi),
    n = 1/sqrt(1 - e2 * sinphi * sinphi),
      // See Wikipedia article on latitude
    hlon = std::cos(phi) * n,
    hlat = (1 - e2) * n * n * n,
    dlon = lon1 - lon0;
  if (dlon >= 180) dlon -= 360;
  else if (dlon < -180) dlon += 360;
  return a * Constants::degree() *
    Math::hypot((lat1 - lat0) * hlat, dlon * hlon);
}

class TestData {
  // Read test data with one line of buffering
public:
  typedef GeographicLib::Math::extended extended;
private:
  std::istream& _is;
  bool _usesaved;               // Should GetNext use saved values?
  bool _valid;                  // Are there saved values?
  extended _lat0, _lat, _lon, _x, _y, _k;
public:
  TestData(std::istream& is) : _is(is), _usesaved(false), _valid(false) {}
  bool GetNext(extended& lat0, extended& lat, extended& lon,
               extended& x, extended& y, extended& k) {
    if (_usesaved)
      _usesaved = false;
    else {
      _valid = (_is >> _lat0 >> _lat >> _lon >> _x >> _y >> _k);
      if (!_valid)
        return false;
    }
    lat0 = _lat0; lat = _lat; lon = _lon; x = _x; y = _y; k = _k;
    return true;
  }
  bool BackUp() {
    if (!_valid || _usesaved)
      return false;             // Can't backup up again
    else
      return _usesaved = true;  // Set flag for GetNext
  }
};

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"ConicTest -l -s\n\
$Id$\n\
\n\
Checks conic projections\n";
  return retval;
}

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  using namespace std;
  typedef Math::real real;
  typedef Math::extended extended;
  if (false) {
    for (int i = -12; i <= 0; ++i) {
      extended
        colat = std::pow(extended(10), i),
        lat = 90 - colat,
        coslat1 = cos(lat*Math::edegree()),
        coslat2 = sin(colat*Math::edegree());
      std::cout << i << " " << coslat1 << " " << coslat2 << " " <<  (coslat1-coslat2)/coslat2 << "\n";
      std::cout << i << " " << asin(coslat1)/Math::edegree() - colat << "\n";
    }
    return 0;
  }
  bool lambert = true;
  bool albert = false;
  bool checkstdlats = false;
  real a = Constants::WGS84_a(), r = Constants::WGS84_r();
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-l") {
      lambert = true;
      albert = false;
    } else if (arg == "-s") {
      checkstdlats = true;
    } else if (arg == "-e") {
      if (m + 2 >= argc) return usage(1);
      try {
        a = DMS::Decode(std::string(argv[m + 1]));
        r = DMS::Decode(std::string(argv[m + 2]));
      }
      catch (const std::exception& e) {
        std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
        return 1;
      }
      m += 2;
    } else
      return usage(arg != "-h");
  }

  try {
    if (checkstdlats) {         // stdin contains lat1 lat2 lat0 k0
      std::cout << std::setprecision(14);
      extended quant = 1e12L;
      while (true) {
        extended lat1, lat2, lat0, k0;
        if (!(cin >> lat1 >> lat2 >> lat0 >> k0))
          break;
        int
          sign1 = lat1 < 0 ? -1 : 1,
          sign2 = lat2 < 0 ? -1 : 1;
        lat1 = floor(sign1 * lat1 * quant + 0.5L);
        lat2 = floor(sign2 * lat2 * quant + 0.5L);
        extended
          colat1 = (90 * quant - lat1) / quant,
          colat2 = (90 * quant - lat2) / quant;
        lat1 /= quant;
        lat2 /= quant;
        extended
          sinlat1 = sign1 * (lat1 < 45 ? sin(lat1 * Math::edegree())
                             : cos(colat1 * Math::edegree())),
          sinlat2 = sign2 * (lat2 < 45 ? sin(lat2 * Math::edegree())
                             : cos(colat2 * Math::edegree())),
          coslat1 = (lat1 < 45 ? cos(lat1 * Math::edegree())
                     : sin(colat1 * Math::edegree())),
          coslat2 = (lat2 < 45 ? cos(lat2 * Math::edegree())
                     : sin(colat2 * Math::edegree()));
        lat1 *= sign1;
        lat2 *= sign2;
        const LambertConformalConic lam(a, r, /* real(lat1), real(lat2), */
                                        real(sinlat1), real(coslat1),
                                        real(sinlat2), real(coslat2),
                                        real(1));
        extended lat0a = lam.OriginLatitude(), k0a = lam.CentralScale();
        if (abs(lat0a-lat0) > 0.4e-13L || abs(k0a - k0) > 7e-15L * k0 )
          cout << lat1 << " " << lat2 << " " << lat0 << " " << lat0a << " " << lat0a - lat0 << " " << (k0a - k0)/k0 << "\n";
      }
    } else { // Check projection
      // stdin contains lat0 lat lon x y k
      TestData testset(std::cin);
      cout << setprecision(8);
      while (true) {
        extended lat0, lat, lon, x, y, k;
        if (!testset.GetNext(lat0, lat, lon, x, y, k))
          break;
        if (!testset.BackUp())
          break;
        int
          sign0 = lat0 < 0 ? -1 : 1;
        extended quant = 1e12L;
        extended
          lat00 = floor(sign0 * lat0 * quant + 0.5L),
          colat00 = (90 * quant - lat00) / quant;
        lat00 /= quant;
        real
          sinlat0 = real(sign0 * (lat00 < 45 ? sin(lat00 * Math::edegree())
                                  : cos(colat00 * Math::edegree()))),
          coslat0 = real(lat00 < 45 ? cos(lat00 * Math::edegree())
                         : sin(colat00 * Math::edegree()));
        const LambertConformalConic lcc(a, r,
                                        sinlat0, coslat0, sinlat0, coslat0,
                                        real(1));
        real maxerrf = 0, maxerrr = 0, maxerrk = 0;
        real latf = 0, lonf = 0, latr = 0, lonr = 0, latk = 0, lonk = 0;
        // std::cout << "New lat0: " << lat0 << "\n";
        while (true) {
          extended lat0x;
          if (!testset.GetNext(lat0x, lat, lon, x, y, k))
            break;
          if (lat0 != lat0x) {
            testset.BackUp();
            break;
          }
          // std::cout << "Process: " << lat0 << " " << lat << " " << lon << " " << k << "\n";
          real lata, lona, xa, ya, gamma, ka, kb;
          lcc.Forward(real(0), real(lat), real(lon), xa, ya, gamma, ka);
          real errf = Math::hypot(extended(xa) - x, extended(ya) - y) / k;
          real errk = abs(extended(ka) - k);
          lcc.Reverse(real(0), real(x), real(y), lata, lona, gamma, kb);
          real errr = dist(extended(a), extended(r),
                           lat, lon, extended(lata), extended(lona));
          std::cout << lata << " " << lona << " " << xa << " " << ya << " "
                    << ka << " " << kb << "\n";
          errk = max(errk, real(abs(extended(kb) - k)));
          errk /= k;
          if (errf > maxerrf) {
            maxerrf = errf;
            latf = lat;
            lonf = lon;
          }
          if (errr > maxerrr) {
            maxerrr = errr;
            latr = lat;
            lonr = lon;
          }
          if (errk > maxerrk && abs(lat) < 89) {
            maxerrk = errk;
            latk = lat;
            lonk = lon;
          }
          std::cout << lat0 << " " << lat << " " << lon << " "
                    << errf << " " << errr << " " << errk << "\n";
        }
        std::cout << "Max errf: " << lat0 << " "
                  << maxerrf << " " << latf << " " << lonf << "\n";
        std::cout << "Max errr: " << lat0 << " "
                  << maxerrr << " " << latr << " " << lonr << "\n";
        std::cout << "Max errk: " << lat0 << " "
                  << maxerrk << " " << latk << " " << lonk << "\n";
      }
    }
  }
  catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

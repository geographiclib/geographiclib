/**
 * \file LambertConformalConic.hpp
 * \brief Header for GeographicLib::LambertConformalConic class
 *
 * Copyright (c) Charles Karney (2010, 2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_HPP)
#define GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_HPP "$Id$"

#include "GeographicLib/Constants.hpp"
#include <algorithm>

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
   * This is a implementation of the equations in Snyder except that divided
   * differences have been used to transform the expressions into ones which
   * may be evaluated accurately and that Newton's method is used to invert the
   * projection.  In this implementation, the projection correctly becomes the
   * Mercator projection or the polar sterographic projection when the standard
   * latitude is the equator or a pole.  The accuracy of the projections is
   * about 10 nm.
   *
   * The ellipsoid parameters, the standard parallels, and the scale on the
   * standard parallels are set in the constructor.  Internally, the case with
   * two standard parallels is converted into a single standard parallel, the
   * latitude of tangency (also the latitude of minimum scale), with a scale
   * specified on this parallel.  This latitude is also used as the latitude of
   * origin which is returned by LambertConformalConic::OriginLatitude.  The
   * scale on the latitude of origin is given by
   * LambertConformalConic::CentralScale.  The case with two distinct standard
   * parallels where one is a pole is singular and is disallowed.  The central
   * meridian (which is a trivial shift of the longitude) is specified as the
   * \e lon0 argument of the LambertConformalConic::Forward and
   * LambertConformalConic::Reverse functions.  There is no provision in this
   * class for specifying a false easting or false northing or a different
   * latitude of origin.  However these are can be simply included by the
   * calling function.  For example the Pennsylvania South state coordinate
   * system (<a href="http://www.spatialreference.org/ref/epsg/3364/">
   * EPSG:3364</a>) is obtained by:
   \code
   const double
     a = GeographicLib::Constants::WGS84_a<double>(),
     r = 298.257222101,                        // GRS80
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
  class LambertConformalConic {
  private:
    typedef Math::real real;
    const real _a, _r, _f, _fm, _e2, _e, _e2m;
    real _sign, _n, _nc, _t0nm1, _scale, _lat0, _k0;
    real _scbet0, _tchi0, _scchi0, _psi0, _nrho0;
    static const real eps, epsx, tol, ahypover;
    static const int numit = 5;
    static inline real sq(real x) throw() { return x * x; }
    static inline real hyp(real x) throw() { return Math::hypot(real(1), x); }
    // e * atanh(e * x) = log( ((1 + e*x)/(1 - e*x))^(e/2) ) if f >= 0
    // - sqrt(-e2) * atan( sqrt(-e2) * x)                    if f < 0
    inline real eatanhe(real x) const throw() {
      return _f >= 0 ? _e * Math::atanh(_e * x) : - _e * std::atan(_e * x);
    }
    // Divided differences
    // Definition: Df(x,y) = (f(x)-f(y))/(x-y)
    // See: W. M. Kahan and R. J. Fateman,
    // Symbolic computation of divided differences,
    // SIGSAM Bull. 33(3), 7-28 (1999)
    // http://doi.acm.org/10.1145/334714.334716
    // http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf
    //
    // General rules
    // h(x) = f(g(x)): Dh(x,y) = Df(g(x),g(y))*Dg(x,y)
    // h(x) = f(x)*g(x):
    //        Dh(x,y) = Df(x,y)*g(x) + Dg(x,y)*f(y)
    //                = Df(x,y)*g(y) + Dg(x,y)*f(x)
    //                = Df(x,y)*(g(x)+g(y))/2 + Dg(x,y)*(f(x)+f(y))/2
    //
    // hyp(x) = sqrt(1+x^2): Dhyp(x,y) = (x+y)/(hyp(x)+hyp(y))
    static inline real Dhyp(real x, real y, real hx, real hy) throw()
    // hx = hyp(x)
    { return (x + y) / (hx + hy); }
    // sn(x) = x/sqrt(1+x^2): Dsn(x,y) = (x+y)/((sn(x)+sn(y))*(1+x^2)*(1+y^2))
    static inline real Dsn(real x, real y, real sx, real sy) throw() {
      // sx = x/hyp(x)
      real t = x * y;
      return t > 0 ? (x + y) * sq( (sx * sy)/t ) / (sx + sy) :
        (x - y != 0 ? (sx - sy) / (x - y) : 1);
    }
    // Dlog1p(x,y) = log1p((x-y)/(1+y)/(x-y)
    static inline real Dlog1p(real x, real y) throw() {
      real t = x - y; if (t < 0) { t = -t; y = x; }
      return t != 0 ? Math::log1p(t / (1 + y)) / t : 1 / (1 + x);
    }
    // Dexp(x,y) = exp((x+y)/2) * 2*sinh((x-y)/2)/(x-y)
    static inline real Dexp(real x, real y) throw() {
      real t = (x - y)/2;
      return (t != 0 ? sinh(t)/t : real(1)) * exp((x + y)/2);
    }
    // Dsinh(x,y) = 2*sinh((x-y)/2)/(x-y) * cosh((x+y)/2)
    //   cosh((x+y)/2) = (c+sinh(x)*sinh(y)/c)/2
    //   c=sqrt((1+cosh(x))*(1+cosh(y)))
    //   cosh((x+y)/2) = sqrt( (sinh(x)*sinh(y) + cosh(x)*cosh(y) + 1)/2 )
    static inline real Dsinh(real x, real y, real sx, real sy, real cx, real cy)
      // sx = sinh(x), cx = cosh(x)
      throw() {
      // real t = (x  - y)/2, c = sqrt((1 + cx) * (1 + cy));
      // return (t != 0 ? sinh(t)/t : real(1)) * (c + sx * sy / c) /2;
      real t = (x  - y)/2;
      return (t != 0 ? sinh(t)/t : real(1)) * sqrt((sx * sy + cx * cy + 1) /2);
    }
    // Dasinh(x,y) = asinh((x-y)*(x+y)/(x*sqrt(1+y^2)+y*sqrt(1+x^2)))/(x-y)
    //             = asinh((x*sqrt(1+y^2)-y*sqrt(1+x^2)))/(x-y)
    static inline real Dasinh(real x, real y, real hx, real hy) throw() {
      // hx = hyp(x)
      real t = x - y;
      return t != 0 ?
        Math::asinh(x*y > 0 ? t * (x+y) / (x*hy + y*hx) : x*hy - y*hx) / t :
        1/hx;
    }
    // Deatanhe(x,y) = eatanhe((x-y)/(1-e^2*x*y))/(x-y)
    inline real Deatanhe(real x, real y) const throw() {
      real t = x - y, d = 1 - _e2 * x * y;
      return t != 0 ? eatanhe(t / d) / t : _e2 / d;
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
     * @param[in] k0 scale on the standard parallel.
     *
     * An exception is thrown if \e a or \e k0 is not positive or if \e stdlat
     * is not in the range [-90, 90].
     **********************************************************************/
    LambertConformalConic(real a, real r, real stdlat, real k0);

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
     * stdlat1 is not equal \e stdlat2.
     **********************************************************************/
    LambertConformalConic(real a, real r, real stdlat1, real stdlat2, real k1);

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
     * This routine computes the latitude of origin and the scale at this
     * latitude.  In the case where \e lat1 and \e lat2 are different, the
     * errors in this routines are as follows: if \e dlat = abs(\e lat2 - \e
     * lat1) <= 160<sup>o</sup> and max(abs(\e lat1), abs(\e lat2)) <= 90 -
     * min(0.0002, 2.2e-6(180 - \e dlat), 6e-8\e dlat<sup>2</sup>) (in
     * degrees), then the error in the latitude of origin is less than
     * 4.5e-14<sup>o</sup> and the relative error in the scale is less than
     * 7e-15.
     **********************************************************************/
    LambertConformalConic(real a, real r,
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
     * thrown if \e k is not positive or if \e stdlat is not in the range [-90,
     * 90]
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
     * error in the projection is less than about 10 nm (true distance) and the
     * errors in the meridian convergence and scale are consistent with this.
     * The values of \e x and \e y returned for points which project to
     * infinity (i.e., one or both of the poles) will be large but finite.
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
     * The error in the projection is less than about 10 nm (true distance) and
     * the errors in the meridian convergence and scale are consistent with
     * this.
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
     *   latitude of origin.
     **********************************************************************/
    Math::real CentralScale() const throw() { return _k0; }
    ///@}

    /**
     * A global instantiation of LambertConformalConic with the WGS84
     * ellipsoid, \e stdlat = 0, and \e k0 = 1.  This degenerates to the
     * Mercator projection.
     **********************************************************************/
    static const LambertConformalConic Mercator;
  };

} // namespace GeographicLib

#endif

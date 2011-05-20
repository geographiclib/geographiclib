/**
 * \file LambertConformalConic.hpp
 * \brief Header for GeographicLib::LambertConformalConic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_HPP)
#define GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_HPP "$Id: LambertConformalConic.hpp 6867 2010-09-11 13:04:26Z karney $"

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
  class LambertConformalConic {
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
     * stdlat1 is not equal \e stdlat2
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
    const static LambertConformalConic Mercator;
  };

} // namespace GeographicLib

#endif

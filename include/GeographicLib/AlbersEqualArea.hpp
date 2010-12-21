/**
 * \file AlbersEqualArea.hpp
 * \brief Header for GeographicLib::AlbersEqualArea class
 *
 * Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ALBERSEQUALAREA_HPP)
#define GEOGRAPHICLIB_ALBERSEQUALAREA_HPP "$Id$"

#include "GeographicLib/Constants.hpp"
#include <algorithm>

namespace GeographicLib {

  /**
   * \brief Albers Equal Area Conic Projection
   *
   * Implementation taken from the report,
   * - J. P. Snyder,
   *   <a href="http://pubs.er.usgs.gov/usgspubs/pp/pp1395"> Map Projections: A
   *   Working Manual</a>, USGS Professional Paper 1395 (1987),
   *   pp. 101&ndash;102.
   *
   * This is a implementation of the equations in Snyder except that divided
   * differences will be [have been] used to transform the expressions into
   * ones which may be evaluated accurately.  [In this implementation, the
   * projection correctly becomes the cylindrical equal area or the azimuthal
   * equal area projection when the standard latitude is the equator or a
   * pole.]
   *
   * The ellipsoid parameters, the standard parallels, and the scale on the
   * standard parallels are set in the constructor.  Internally, the case with
   * two standard parallels is converted into a single standard parallel, the
   * latitude of minimum azimuthal scale, with an azimuthal scale specified on
   * this parallel.  This latitude is also used as the latitude of origin which
   * is returned by AlbersEqualArea::OriginLatitude.  The azimuthal scale on
   * the latitude of origin is given by AlbersEqualArea::CentralScale.  The
   * case with two standard parallels at opposite poles is singular and is
   * disallowed.  The central meridian (which is a trivial shift of the
   * longitude) is specified as the \e lon0 argument of the
   * AlbersEqualArea::Forward and AlbersEqualArea::Reverse functions.
   * AlbersEqualArea::Forward and AlbersEqualArea::Reverse also return the
   * meridian convergence, \e gamma, and azimuthal scale, \e k.  A small square
   * aligned with the cardinal directions is projected to a rectangle with
   * dimensions \e k (in the E-W direction) and 1/\e k (in the N-S direction).
   * The E-W sides of the rectangle are oriented \e gamma degrees
   * counter-clockwise from the \e x axis.  There is no provision in this class
   * for specifying a false easting or false northing or a different latitude
   * of origin.
   **********************************************************************/
  class AlbersEqualArea {
  private:
    typedef Math::real real;
    const real _a, _r, _f, _fm, _e2, _e, _e2m, _qZ, _qx;
    real _sign, _lat0, _k0;
    real _n0, _m02, _nrho0, _k2, _txi0, _scxi0, _sxi0;
    static const real eps, epsx, epsx2, tol, tol0, ahypover;
    static const int numit = 5;   // Newton iterations in Reverse
    static const int numit0 = 20; // Newton iterations in Init
    static inline real sq(real x) throw() { return x * x; }
    static inline real hyp(real x) throw() { return Math::hypot(real(1), x); }
    // atanh(      e   * x)/      e   if f > 0
    // atan (sqrt(-e2) * x)/sqrt(-e2) if f < 0
    // x                              if f = 0
    inline real atanhee(real x) const throw() {
      return _f > 0 ? Math::atanh(_e * x)/_e :
        (_f < 0 ? std::atan(_e * x)/_e : x);
    }
    // return atanh(sqrt(x))/sqrt(x) - 1, accurate for small x
    static real atanhxm1(real x) throw();

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
    //        Dh(x,y) = Df(x,y)*(g(x)+g(y))/2 + Dg(x,y)*(f(x)+f(y))/2
    //
    // sn(x) = x/sqrt(1+x^2): Dsn(x,y) = (x+y)/((sn(x)+sn(y))*(1+x^2)*(1+y^2))
    static inline real Dsn(real x, real y, real sx, real sy) throw() {
      // sx = x/hyp(x)
      real t = x * y;
      return t > 0 ? (x + y) * sq( (sx * sy)/t ) / (sx + sy) :
        (x - y != 0 ? (sx - sy) / (x - y) : 1);
    }
    // Datanhee(x,y) = atanhee((x-y)/(1-e^2*x*y))/(x-y)
    inline real Datanhee(real x, real y) const throw() {
      real t = x - y, d = 1 - _e2 * x * y;
      return t != 0 ? atanhee(t / d) / t : 1 / d;
    }
    // DDatanhee(x,y) = (Datanhee(1,y) - Datanhee(1,x))/(y-x)
    real DDatanhee(real x, real y) const throw();
    void Init(real sphi1, real cphi1, real sphi2, real cphi2, real k1) throw();
    real txif(real tphi) const throw();
    real tphif(real txi) const throw();
  public:

    /**
     * Constructor with a single standard parallel.
     *
     * @param[in] a equatorial radius of ellipsoid (meters)
     * @param[in] r reciprocal flattening of ellipsoid.  Setting \e r = 0
     *   implies \e r = inf or flattening = 0 (i.e., a sphere).  Negative \e r
     *   indicates a prolate ellipsoid.
     * @param[in] stdlat standard parallel (degrees), the circle of tangency.
     * @param[in] k0 azimuthal scale on the standard parallel.
     *
     * An exception is thrown if \e a or \e k0 is not positive or if \e stdlat
     * is not in the range [-90, 90].
     **********************************************************************/
    AlbersEqualArea(real a, real r, real stdlat, real k0);

    /**
     * Constructor with two standard parallels.
     *
     * @param[in] a equatorial radius of ellipsoid (meters)
     * @param[in] r reciprocal flattening of ellipsoid.  Setting \e r = 0
     *   implies \e r = inf or flattening = 0 (i.e., a sphere).  Negative \e r
     *   indicates a prolate ellipsoid.
     * @param[in] stdlat1 first standard parallel (degrees).
     * @param[in] stdlat2 second standard parallel (degrees).
     * @param[in] k1 azimuthal scale on the standard parallels.
     *
     * An exception is thrown if \e a or \e k0 is not positive or if \e stdlat1
     * or \e stdlat2 is not in the range [-90, 90].  In addition, an exception
     * is thrown if \e stdlat1 and \e stdlat2 are opposite poles.
     **********************************************************************/
    AlbersEqualArea(real a, real r, real stdlat1, real stdlat2, real k1);

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
     * @param[in] k1 azimuthal scale on the standard parallels.
     *
     * This allows parallels close to the poles to be specified accurately.
     * This routine computes the latitude of origin and the azimuthal scale at
     * this latitude.  If \e dlat = abs(\e lat2 - \e lat1) <= 160<sup>o</sup>,
     * then the error in the latitude of origin is less than
     * 4.5e-14<sup>o</sup>.
     **********************************************************************/
    AlbersEqualArea(real a, real r,
                    real sinlat1, real coslat1,
                    real sinlat2, real coslat2,
                    real k1);

    /**
     * Set the azimuthal scale for the projection.
     *
     * @param[in] lat (degrees).
     * @param[in] k azimuthal scale at latitude \e lat (default 1).
     *
     * This allows a "latitude of conformality" to be specified.  An exception
     * is thrown if \e k is not positive or if \e lat is not in the range (-90,
     * 90).
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
     * @param[out] k azimuthal scale of projection at point; the radial
     *   scale is the 1/\e k.
     *
     * The latitude origin is given by AlbersEqualArea::LatitudeOrigin().  No
     * false easting or northing is added and \e lat should be in the range
     * [-90, 90]; \e lon and \e lon0 should be in the range [-180, 360].  The
     * values of \e x and \e y returned for points which project to infinity
     * (i.e., one or both of the poles) will be large but finite.
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
     * @param[out] k azimuthal scale of projection at point; the radial
     *   scale is the 1/\e k.
     *
     * The latitude origin is given by AlbersEqualArea::LatitudeOrigin().  No
     * false easting or northing is added.  \e lon0 should be in the range
     * [-180, 360].  The value of \e lon returned is in the range [-180, 180).
     * The value of \e lat returned is in the range [-90,90].  If the input
     * point is outside the legal projected space the nearest pole is returned.
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon, real& gamma, real& k) const throw();

    /**
     * AlbersEqualArea::Forward without returning the convergence and
     * scale.
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y) const throw() {
      real gamma, k;
      Forward(lon0, lat, lon, x, y, gamma, k);
    }

    /**
     * AlbersEqualArea::Reverse without returning the convergence and
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
     * This is the latitude of minimum azimuthal scale and equals the \e stdlat
     * in the 1-parallel constructor and lies between \e stdlat1 and \e stdlat2
     * in the 2-parallel constructors.
     **********************************************************************/
    Math::real OriginLatitude() const throw() { return _lat0; }

    /**
     * @return central scale for the projection.  This is the azimuthal scale
     *   on the latitude of origin.
     **********************************************************************/
    Math::real CentralScale() const throw() { return _k0; }
    ///@}

    /**
     * A global instantiation of AlbersEqualArea with the WGS84 ellipsoid, \e
     * stdlat = 0, and \e k0 = 1.  This degenerates to the cylindrical equal
     * area projection.
     **********************************************************************/
    static const AlbersEqualArea CylindricalEqualArea;

    /**
     * A global instantiation of AlbersEqualArea with the WGS84 ellipsoid, \e
     * stdlat = 90<sup>o</sup>, and \e k0 = 1.  This degenerates to the
     * Lambert azimuthal equal area projection.
     **********************************************************************/
    static const AlbersEqualArea AzimuthalEqualAreaNorth;

    /**
     * A global instantiation of AlbersEqualArea with the WGS84 ellipsoid, \e
     * stdlat = -90<sup>o</sup>, and \e k0 = 1.  This degenerates to the
     * Lambert azimuthal equal area projection.
     **********************************************************************/
    static const AlbersEqualArea AzimuthalEqualAreaSouth;
  };

} // namespace GeographicLib

#endif

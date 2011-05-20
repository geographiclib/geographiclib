/**
 * \file TransverseMercator.hpp
 * \brief Header for GeographicLib::TransverseMercator class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRANSVERSEMERCATOR_HPP)
#define GEOGRAPHICLIB_TRANSVERSEMERCATOR_HPP "$Id: TransverseMercator.hpp 6867 2010-09-11 13:04:26Z karney $"

#include "GeographicLib/Constants.hpp"

#if !defined(TM_TX_MAXPOW)
/**
 * The order of the series approximation used in TransverseMercator.
 * TM_TX_MAXPOW can be set to any integer in [4, 8].
 **********************************************************************/
#define TM_TX_MAXPOW \
(GEOGRAPHICLIB_PREC == 1 ? 6 : GEOGRAPHICLIB_PREC == 0 ? 4 : 8)
#endif

namespace GeographicLib {

  /**
   * \brief Transverse Mercator Projection
   *
   * This uses Kr&uuml;ger's method which evaluates the projection and its
   * inverse in terms of a series.  See
   *  - L. Kr&uuml;ger,
   *    <a href="http://dx.doi.org/10.2312/GFZ.b103-krueger28"> Konforme
   *    Abbildung des Erdellipsoids in der Ebene</a> (Conformal mapping of the
   *    ellipsoidal earth to the plane), Royal Prussian Geodetic Institute, New
   *    Series 52, 172 pp. (1912).
   *
   * Kr&uuml;ger's method has been extended from 4th to 6th order.  The maximum
   * errors is 5 nm (ground distance) for all positions within 35 degrees of
   * the central meridian.  The error in the convergence is 2e-15&quot; and the
   * relative error in the scale is 6e-12%%.  (See \ref tmerrors for the weasel
   * words.)  The speed penalty in going to 6th order is only about 1%.
   * TransverseMercatorExact is an alternative implementation of the projection
   * using exact formulas which yield accurate (to 8 nm) results over the
   * entire ellipsoid.
   *
   * The ellipsoid parameters and the central scale are set in the constructor.
   * The central meridian (which is a trivial shift of the longitude) is
   * specified as the \e lon0 argument of the TransverseMercator::Forward and
   * TransverseMercator::Reverse functions.  The latitude of origin is taken to
   * be the equator.  There is no provision in this class for specifying a
   * false easting or false northing or a different latitude of origin.
   * However these are can be simply included by the calling funtcion.  For
   * example, the UTMUPS class applies the false easting and false northing for
   * the UTM projections.  A more complicated example is the British National
   * Grid (<a href="http://www.spatialreference.org/ref/epsg/7405/">
   * EPSG:7405</a>) which requires the use of a latitude of origin.  This is
   * accommodated by (constants from
   * <a href="http://www.ordnancesurvey.co.uk/oswebsite/gps/information/coordinatesystemsinfo/guidecontents/guidea.html">
   * A guide to coordinate systems in Great Britain</a>):
   \code
   const double
     a = 6377563.396, b = 6356256.910, r = a/(a - b), // Airy 1830 ellipsoid
     k0 = 0.9996012717, lat0 = 49, lon0 = -2, // central scale and origin
     fe = 400000, fn = -100000;               // false easting and northing
   // Set up basic projection
   const GeographicLib::TransverseMercator OSGB(a, r, k0);
   double x0, y0;
   {
     // Transform origin point
     OSGB.Forward(lon0, lat0, lon0, x0, y0);
     x0 -= fe; y0 -= fn;         // Combine result with false origin
   }
   double lat, lon, x, y;
   // Sample conversion from geodetic to OSGB grid
   std::cin >> lat >> lon;
   OSGB.Forward(lon0, lat, lon, x, y);
   x -= x0; y -= y0;
   std::cout << x << " " << y << "\n";
   // Sample conversion from OSGB grid to geodetic
   std::cin >> x >> y;
   x += x0; y += y0;
   OSGB.Reverse(lon0, x, y, lat, lon);
   std::cout << lat << " " << lon << "\n";
   \endcode
   *
   * See TransverseMercator.cpp for more information on the implementation.
   *
   * See \ref transversemercator for a discussion of this projection.
   **********************************************************************/

  class TransverseMercator {
  private:
    typedef Math::real real;
    static const int maxpow = TM_TX_MAXPOW;
    static const real tol, overflow;
    static const int numit = 5;
    const real _a, _r, _f, _k0, _e2, _e, _e2m,  _c, _n;
    // _alp[0] and _bet[0] unused
    real _a1, _b1, _alp[maxpow + 1], _bet[maxpow + 1];
    static inline real sq(real x) throw() { return x * x; }
    // tan(x) for x in [-pi/2, pi/2] ensuring that the sign is right
    static inline real tanx(real x) throw() {
      real t = std::tan(x);
      return x >= 0 ? (t >= 0 ? t : overflow) : (t < 0 ? t : -overflow);
    }
    // Return e * atanh(e * x) for f >= 0, else return
    // - sqrt(-e2) * atan( sqrt(-e2) * x) for f < 0
    inline real eatanhe(real x) const throw() {
      return _f >= 0 ? _e * Math::atanh(_e * x) : - _e * std::atan(_e * x);
    }
  public:

    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters)
     * @param[in] r reciprocal flattening.  Setting \e r = 0 implies \e r = inf
     *   or flattening = 0 (i.e., a sphere).  Negative \e r indicates a prolate
     *   ellipsoid.
     * @param[in] k0 central scale factor.
     *
     * An exception is thrown if either of the axes of the ellipsoid or \e k0
     * is not positive.
     **********************************************************************/
    TransverseMercator(real a, real r, real k0);

    /**
     * Forward projection, from geographic to transverse Mercator.
     *
     * @param[in] lon0 central meridian of the projection (degrees).
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added. \e lat should be in the range
     * [-90, 90]; \e lon and \e lon0 should be in the range [-180, 360].
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y, real& gamma, real& k) const throw();

    /**
     * Reverse projection, from transverse Mercator to geographic.
     *
     * @param[in] lon0 central meridian of the projection (degrees).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added.  \e lon0 should be in the range
     * [-180, 360].  The value of \e lon returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon, real& gamma, real& k) const throw();

    /**
     * TransverseMercator::Forward without returning the convergence and scale.
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y) const throw() {
      real gamma, k;
      Forward(lon0, lat, lon, x, y, gamma, k);
    }

    /**
     * TransverseMercator::Reverse without returning the convergence and scale.
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
     * @return \e k0 central scale for the projection.  This is the value of \e
     *   k0 used in the constructor and is the scale on the central meridian.
     **********************************************************************/
    Math::real CentralScale() const throw() { return _k0; }
    ///@}

    /**
     * A global instantiation of TransverseMercator with the WGS84 ellipsoid
     * and the UTM scale factor.  However, unlike UTM, no false easting or
     * northing is added.
     **********************************************************************/
    const static TransverseMercator UTM;
  };

} // namespace GeographicLib

#endif

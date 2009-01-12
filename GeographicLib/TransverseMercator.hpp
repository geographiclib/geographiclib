/**
 * \file TransverseMercator.hpp
 * \brief Header for GeographicLib::TransverseMercator class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(TRANSVERSEMERCATOR_HPP)
#define TRANSVERSEMERCATOR_HPP "$Id$"

#include <cmath>

#if !defined(TM_TX_MAXPOW)
/**
 * The order of the series approximation used in
 * GeographicLib::TransverseMercator.  TM_TX_MAXPOW can be set to any integer
 * in [4, 8].
 **********************************************************************/
#define TM_TX_MAXPOW 6
#endif

namespace GeographicLib {
  /**
   * \brief Transverse Mercator Projection
   *
   * This uses Kr&uuml;ger's method which evaluates the projection and its
   * inverse in terms of a series.  See L. Kr&uuml;ger, <a
   * href="http://dx.doi.org/10.2312/GFZ.b103-krueger28"> Konforme
   * Abbildung des Erdellipsoids in der Ebene</a> (Conformal mapping of the
   * ellipsoidal earth to the plane), Royal Prussian Geodetic Institute, New
   * Series 52, 172 pp. (1912).
   *
   * Kr&uuml;ger's method has been extended from 4th to 6th order.  The maximum
   * errors is 5 nm (ground distance) for all positions within 35 degrees of
   * the central meridian.  The error in the convergence is 2e-15" and the
   * relative error in the scale is 6e-12%%.  (See \ref errors for the weasel
   * words.)  The speed penalty in going to 6th order is only about 1%.
   * GeographicLib::TransverseMercatorExact is an alternative implementation of
   * the projection using exact formulas which yield accurate (to 8 nm)
   * results over the entire ellipsoid.
   *
   * See TransverseMercator.cpp for more information on the implementation.
   *
   * See \ref transversemercator for a discussion of this projection.
   **********************************************************************/

  class TransverseMercator {
  private:
    static const int maxpow =
      TM_TX_MAXPOW > 8 ? 8 : (TM_TX_MAXPOW < 4 ? 4 : TM_TX_MAXPOW);
    static const double tol;
    static const int numit = 5;
    const double _a, _f, _k0, _e2, _e, _e2m, _ep2,  _n;
    double _a1, _b1, _h[maxpow], _hp[maxpow];
    static inline double sq(double x) { return x * x; }
#if defined(_MSC_VER)
    // These have poor relative accuracy near x = 0.  However, for mapping
    // applications, we only need good absolute accuracy.
    // For small arguments we would write
    //
    // asinh(x) = asin(x) -x^3/3-5*x^7/56-63*x^11/1408-143*x^15/5120 ...
    // atanh(x) = atan(x) +2*x^3/3+2*x^7/7+2*x^11/11+2*x^15/15
    //
    // The accuracy of asinh is also bad for large negative arguments.  This is
    // easy to fix in the definition of asinh.  Instead we call these functions
    // with positive arguments and enforce the correct parity separately.
    static inline double asinh(double x) { return log(x + sqrt(1 + sq(x))); }
    static inline double atanh(double x) { return log((1 + x)/(1 - x))/2; }
#endif
  public:
    /**
     * Constructor for a ellipsoid radius \e a (meters), flattening \e f,
     * and central scale factor \e k0.
     **********************************************************************/
    TransverseMercator(double a, double f, double k0);
    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * transverse Mercator easting \e x (meters) and northing \e y (meters).
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added.
     **********************************************************************/
    void Forward(double lon0, double lat, double lon,
		 double& x, double& y, double& gamma, double& k) const;
    /**
     * Convert from transverse Mercator easting \e x (meters) and northing \e y
     * (meters) to latitude \e lat (degrees) and longitude \e lon (degrees) .
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added.
     **********************************************************************/
    void Reverse(double lon0, double x, double y,
		 double& lat, double& lon, double& gamma, double& k) const;
    /**
     * A global instantiation of TransverseMercator with the WGS84 ellipsoid
     * and the UTM scale factor.  However, unlike UTM, no false easting or
     * northing is added.
     **********************************************************************/
    const static TransverseMercator UTM;
  };

} // namespace GeographicLib

#endif

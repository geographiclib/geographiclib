/**
 * \file TransverseMercator.hpp
 * \brief Header for GeographicLib::TransverseMercator class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRANSVERSEMERCATOR_HPP)
#define GEOGRAPHICLIB_TRANSVERSEMERCATOR_HPP "$Id$"

#include "GeographicLib/Constants.hpp"

#if !defined(TM_TX_MAXPOW)
/**
 * The order of the series approximation used in
 * GeographicLib::TransverseMercator.  TM_TX_MAXPOW can be set to any integer
 * in [4, 8].
 **********************************************************************/
#define TM_TX_MAXPOW \
(sizeof(real) == sizeof(double) ? 6 : sizeof(real) == sizeof(float) ? 5 : 7)
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
   * the central meridian.  The error in the convergence is 2e-15" and the
   * relative error in the scale is 6e-12%%.  (See \ref tmerrors for the weasel
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
    typedef Math::real real;
    static const int maxpow = TM_TX_MAXPOW;
    static const real tol;
    static const int numit = 5;
    const real _a, _f, _k0, _e2, _e, _e2m,  _c, _n;
    real _a1, _b1, _h[maxpow], _hp[maxpow];
    static inline real sq(real x) throw() { return x * x; }
    // Return e * atanh(e * x) for f >= 0, else return
    // - sqrt(-e2) * atan( sqrt(-e2) * x) for f < 0
    inline real eatanhe(real x) const throw() {
      return _f >= 0 ? _e * Math::atanh(_e * x) : - _e * atan(_e * x);
    }
  public:

    /**
     * Constructor for a ellipsoid radius \e a (meters), reciprocal flattening
     * \e r, and central scale factor \e k0.  Setting \e r = 0 implies \e r =
     * inf or flattening = 0 (i.e., a sphere).  Negative \e r indicates a
     * prolate spheroid.
     **********************************************************************/
    TransverseMercator(Math::real a, Math::real r, Math::real k0) throw();

    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * transverse Mercator easting \e x (meters) and northing \e y (meters).
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added. \e lat should be in the range
     * [-90, 90]; \e lon and \e lon0 should be in the range [-180, 360].
     **********************************************************************/
    void Forward(Math::real lon0, Math::real lat, Math::real lon,
                 Math::real& x, Math::real& y,
                 Math::real& gamma, Math::real& k) const throw();

    /**
     * Convert from transverse Mercator easting \e x (meters) and northing \e y
     * (meters) to latitude \e lat (degrees) and longitude \e lon (degrees) .
     * The central meridian of the transformation is \e lon0 (degrees).  Also
     * return the meridian convergence \e gamma (degrees) and the scale \e k.
     * No false easting or northing is added.  The value of \e lon returned is
     * in the range [-180, 180).
     **********************************************************************/
    void Reverse(Math::real lon0, Math::real x, Math::real y,
                 Math::real& lat, Math::real& lon,
                 Math::real& gamma, Math::real& k) const throw();

    /**
     * A global instantiation of TransverseMercator with the WGS84 ellipsoid
     * and the UTM scale factor.  However, unlike UTM, no false easting or
     * northing is added.
     **********************************************************************/
    const static TransverseMercator UTM;
  };

} // namespace GeographicLib

#endif

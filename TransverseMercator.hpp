/**
 * \file TransverseMercator.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(TRANSVERSEMERCATOR_HPP)
#define TRANSVERSEMERCATOR_HPP "$Id$"

#include <cmath>

#if !defined(TM_TX_MAXPOW)
#define TM_TX_MAXPOW 6
#endif

namespace GeographicLib {
  /**
   * \brief Transverse Mercator Projection
   *
   * This uses Kr&uuml;ger's method which evaluates the projection and its
   * inverse in terms of a series.  See L. Kr&uuml;ger, <a
   * href="http://dx.doi.org/10.2312/GFZ.b103-krueger28"> Konforme
   * Abbildung des Erdellipsoids in der Ebene</a> (Conformal mpping of the
   * Earth ellipsoid to the plane), Royal Prussian Geodetic Institute, New
   * Series 52, 172 p (1912).
   *
   * This implementation follows closely <a
   * href="http://www.jhs-suositukset.fi/suomi/jhs154"> JHS 154, ETRS89 -
   * j&auml;&auml;rjestelm&auml;&auml;n liittyv&auml;t karttaprojektiot,
   * tasokoordinaatistot ja karttalehtijako</a> (Map projections, plane
   * coordinates, and map sheet index for ETRS89), Published by JUHTA, Finnish
   * Geodetic Institute, and the National Land Survey of Finland, 34 p (2006).
   *
   * http://www.jhs-suositukset.fi/suomi/jhs154
   * http://docs.jhs-suositukset.fi/jhs-suositukset/JHS154/JHS154.pdf
   *
   * This is a straight transcription
   * of the formulas in this paper with the following exceptions:
   *  - use of 6th order series instead of 4th order series.  This reduces the
   *    error to about 5nm for the UTM range of coordinates (instead of 200nm);
   *  - use Newton's method instead of plain iteration to solve for latitude in
   *    terms of isometric latitude in the Reverse method;
   *  - use of Horner's representation for evaluating polynomials and Clenshaw's
   *    method for summing trigonometric series;
   *  - several modifications of the formulas to improve the numerical accuracy;
   *  - evaluating the convergence and scale using the expressions for the
   *    projection or its inverse.
   *
   * Other equivatlent implementations are given in
   *  - http://www.ign.fr/telechargement/MPro/geodesie/CIRCE/NTG_76.pdf
   *  - http://www.lantmateriet.se/upload/filer/kartor/geodesi_gps_och_detaljmatning/geodesi/Formelsamling/Gauss_Conformal_Projection.pdf
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
     * A global instantiation of TransverseMercator with thw WGS84 ellipsoid
     * and the UTM scale factor.
     **********************************************************************/
    const static TransverseMercator UTM;
  };

} // namespace GeographicLib

#endif

/**
 * @file GeodesicLine.java
 * @brief Implementation of the net.sf.geographiclib.GeodesicLine class
 *
 * Copyright (c) Charles Karney (2013) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/
package net.sf.geographiclib;
/**
 * @brief A geodesic line
 *
 * GeodesicLine facilitates the determination of a series of points on a
 * single geodesic.  The starting point (\e lat1, \e lon1) and the azimuth \e
 * azi1 are specified in the constructor.  GeodesicLine.Position returns the
 * location of point 2 a distance \e s12 along the geodesic.  Alternatively
 * GeodesicLine.ArcPosition gives the position of point 2 an arc length \e
 * a12 along the geodesic.
 *
 * The calculations are accurate to better than 15 nm (15 nanometers).  See
 * Sec. 9 of
 * <a href="http://arxiv.org/abs/1102.1215v1">arXiv:1102.1215v1</a> for
 * details.  The algorithms used by this class are based on series expansions
 * using the flattening \e f as a small parameter.  These are only accurate
 * for |<i>f</i>| &lt; 0.02; however reasonably accurate results will be
 * obtained for |<i>f</i>| &lt; 0.2.
 *
 * The algorithms are described in
 * - C. F. F. Karney,
 *   <a href="http://dx.doi.org/10.1007/s00190-012-0578-z">
 *   Algorithms for geodesics</a>,
 *   J. Geodesy <b>87</b>, 43--55 (2013);
 *   DOI: <a href="http://dx.doi.org/10.1007/s00190-012-0578-z">
 *   10.1007/s00190-012-0578-z</a>;
 *   addenda: <a href="http://geographiclib.sf.net/geod-addenda.html">
 *   geod-addenda.html</a>.
 * .
 *
 * Here's an example of using this class
 * @code
 * import net.sf.geographiclib.*;
 * public class GeodesicLineTest {
 *   public static void main(String[] args) {
 *     // Print waypoints between JFK and SIN
 *     Geodesic geod = new Geodesic(Constants.WGS84_a, Constants.WGS84_f);
 *     // Alternatively: Geodesic geod = Geodesic.WGS84;
 *     double
 *       lat1 = 40.640, lon1 = -73.779, // JFK
 *       lat2 =  1.359, lon2 = 103.989; // SIN
 *     GeodesicData g = geod.Inverse(lat1, lon1, lat2, lon2,
 *                                   GeodesicMask.DISTANCE |
 *                                   GeodesicMask.AZIMUTH);
 *     // GeodesicMask.LATITUDE and GeodesicMask.AZIMUTH added automatically
 *     GeodesicLine line = new GeodesicLine(geod, lat1, lon1, g.azi1,
 *                                          GeodesicMask.DISTANCE_IN |
 *                                          GeodesicMask.LONGITUDE);
 *     // Alternatively
 *     // GeodesicLine line =  geod.Line(lat1, lon1, g.azi1,
 *     //                                  GeodesicMask.DISTANCE_IN |
 *     //                                  GeodesicMask.LONGITUDE);
 *     double
 *       s12 = g.s12,
 *       a12 = g.a12,
 *       ds0 = 500e3;            // Nominal distance between points = 500 km
 *     int num = (int)(Math.ceil(s12 / ds0)); // The number of intervals
 *     {
 *       // Use intervals of equal length
 *       double ds = s12 / num;
 *       for (int i = 0; i <= num; ++i) {
 *         g = line.Position(i * ds,
 *                           GeodesicMask.LATITUDE |
 *                           GeodesicMask.LONGITUDE );
 *         System.out.println(i + " " + g.lat2 + " " + g.lon2);
 *       }
 *     }
 *     {
 *       // Slightly faster, use intervals of equal arc length
 *       double da = a12 / num;
 *       for (int i = 0; i <= num; ++i) {
 *         g = line.ArcPosition(i * da,
 *                              GeodesicMask.LATITUDE |
 *                              GeodesicMask.LONGITUDE );
 *         System.out.println(i + " " + g.lat2 + " " + g.lon2);
 *       }
 *     }
 *   }
 * }
 * @endcode
 **********************************************************************/

public class GeodesicLine {

  private static final int nC1_ = Geodesic.nC1_;
  private static final int nC1p_ = Geodesic.nC1p_;
  private static final int nC2_ = Geodesic.nC2_;
  private static final int nC3_ = Geodesic.nC3_;
  private static final int nC4_ = Geodesic.nC4_;

  private double _lat1, _lon1, _azi1;
  private double _a, _f, _b, _c2, _f1, _salp0, _calp0, _k2,
    _salp1, _calp1, _ssig1, _csig1, _dn1, _stau1, _ctau1, _somg1, _comg1,
    _A1m1, _A2m1, _A3c, _B11, _B21, _B31, _A4, _B41;
  // index zero elements of _C1a, _C1pa, _C2a, _C3a are unused
  private double _C1a[], _C1pa[], _C2a[], _C3a[],
    _C4a[];    // all the elements of _C4a are used
  private int _caps;

  /**
   * \name Constructors
   **********************************************************************/
  ///@{

  /**
   * Constructor for a geodesic line staring at latitude \e lat1, longitude
   * \e lon1, and azimuth \e azi1 (all in degrees).
   *
   * @param g A Geodesic object used to compute the necessary information
   *   about the GeodesicLine.
   * @param lat1 latitude of point 1 (degrees).
   * @param lon1 longitude of point 1 (degrees).
   * @param azi1 azimuth at point 1 (degrees).
   *
   * \e lat1 should be in the range [&minus;90&deg;, 90&deg;]; \e lon1 and \e
   * azi1 should be in the range [&minus;540&deg;, 540&deg;).
   *
   * If the point is at a pole, the azimuth is defined by keeping the \e lon1
   * fixed and writing \e lat1 = &plusmn;(90&deg; &minus; &epsilon;) and
   * taking the limit &epsilon; &rarr; 0+.
   **********************************************************************/
  public GeodesicLine(Geodesic g,
                      double lat1, double lon1, double azi1) {
    this(g, lat1, lon1, azi1, GeodesicMask.ALL);
  }

  /**
   * Constructor for a geodesic line staring at latitude \e lat1, longitude
   * \e lon1, and azimuth \e azi1 (all in degrees) with a subset of the
   * capabilities included.
   *
   * @param g A Geodesic object used to compute the necessary information
   *   about the GeodesicLine.
   * @param lat1 latitude of point 1 (degrees).
   * @param lon1 longitude of point 1 (degrees).
   * @param azi1 azimuth at point 1 (degrees).
   * @param caps bitor'ed combination of GeodesicLine mask values
   *   specifying the capabilities the GeodesicLine object should possess,
   *   i.e., which quantities can be returned in calls to
   *   GeodesicLine.Position.
   *
   * The net.sf.geographiclib.GeodesicMask values are
   * - \e caps |= GeodesicMask.LATITUDE for the latitude \e lat2; this is
   *   added automatically;
   * - \e caps |= GeodesicMask.LONGITUDE for the latitude \e lon2;
   * - \e caps |= GeodesicMask.AZIMUTH for the latitude \e azi2; this is
   *   added automatically;
   * - \e caps |= GeodesicMask.DISTANCE for the distance \e s12;
   * - \e caps |= GeodesicMask.REDUCEDLENGTH for the reduced length \e m12;
   * - \e caps |= GeodesicMask.GEODESICSCALE for the geodesic scales \e M12
   *   and \e M21;
   * - \e caps |= GeodesicMask.AREA for the area \e S12;
   * - \e caps |= GeodesicMask.DISTANCE_IN permits the length of the
   *   geodesic to be given in terms of \e s12; without this capability the
   *   length can only be specified in terms of arc length;
   * - \e caps |= GeodesicMask.ALL for all of the above;
   **********************************************************************/
  public GeodesicLine(Geodesic g,
                      double lat1, double lon1, double azi1,
                      int caps) {
    _a = g._a;
    _f = g._f;
    _b = g._b;
    _c2 = g._c2;
    _f1 = g._f1;
    // Always allow latitude and azimuth
    _caps = caps | GeodesicMask.LATITUDE | GeodesicMask.AZIMUTH;

    // Guard against underflow in salp0
    azi1 = Geodesic.AngRound(GeoMath.AngNormalize(azi1));
    lon1 = GeoMath.AngNormalize(lon1);
    _lat1 = lat1;
    _lon1 = lon1;
    _azi1 = azi1;
    // alp1 is in [0, pi]
    double alp1 = azi1 * GeoMath.degree;
    // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
    // problems directly than to skirt them.
    _salp1 =          azi1  == -180 ? 0 : Math.sin(alp1);
    _calp1 = Math.abs(azi1) ==   90 ? 0 : Math.cos(alp1);
    double cbet1, sbet1, phi;
    phi = lat1 * GeoMath.degree;
    // Ensure cbet1 = +epsilon at poles
    sbet1 = _f1 * Math.sin(phi);
    cbet1 = Math.abs(lat1) == 90 ? Geodesic.tiny_ : Math.cos(phi);
    { Pair p = Geodesic.SinCosNorm(sbet1, cbet1);
      sbet1 = p.first; cbet1 = p.second; }
    _dn1 = Math.sqrt(1 + g._ep2 * GeoMath.sq(sbet1));

    // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    _salp0 = _salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    _calp0 = GeoMath.hypot(_calp1, _salp1 * sbet1);
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    _ssig1 = sbet1; _somg1 = _salp0 * sbet1;
    _csig1 = _comg1 = sbet1 != 0 || _calp1 != 0 ? cbet1 * _calp1 : 1;
    { Pair p = Geodesic.SinCosNorm(_ssig1, _csig1);
      _ssig1 = p.first; _csig1 = p.second; } // sig1 in (-pi, pi]
    // Geodesic.SinCosNorm(_somg1, _comg1); -- don't need to normalize!

    _k2 = GeoMath.sq(_calp0) * g._ep2;
    double eps = _k2 / (2 * (1 + Math.sqrt(1 + _k2)) + _k2);

    if ((_caps & GeodesicMask.CAP_C1) != 0) {
      _A1m1 = Geodesic.A1m1f(eps);
      _C1a = new double[nC1_ + 1];
      Geodesic.C1f(eps, _C1a);
      _B11 = Geodesic.SinCosSeries(true, _ssig1, _csig1, _C1a);
      double s = Math.sin(_B11), c = Math.cos(_B11);
      // tau1 = sig1 + B11
      _stau1 = _ssig1 * c + _csig1 * s;
      _ctau1 = _csig1 * c - _ssig1 * s;
      // Not necessary because C1pa reverts C1a
      //    _B11 = -SinCosSeries(true, _stau1, _ctau1, _C1pa, nC1p_);
    }

    if ((_caps & GeodesicMask.CAP_C1p) != 0) {
      _C1pa = new double[nC1p_ + 1];
      Geodesic.C1pf(eps, _C1pa);
    }

    if ((_caps & GeodesicMask.CAP_C2) != 0) {
      _C2a = new double[nC2_ + 1];
      _A2m1 = Geodesic.A2m1f(eps);
      Geodesic.C2f(eps, _C2a);
      _B21 = Geodesic.SinCosSeries(true, _ssig1, _csig1, _C2a);
    }

    if ((_caps & GeodesicMask.CAP_C3) != 0) {
      _C3a = new double[nC3_];
      g.C3f(eps, _C3a);
      _A3c = -_f * _salp0 * g.A3f(eps);
      _B31 = Geodesic.SinCosSeries(true, _ssig1, _csig1, _C3a);
    }

    if ((_caps & GeodesicMask.CAP_C4) != 0) {
      _C4a = new double[nC4_];
      g.C4f(eps, _C4a);
      // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
      _A4 = GeoMath.sq(_a) * _calp0 * _salp0 * g._e2;
      _B41 = Geodesic.SinCosSeries(false, _ssig1, _csig1, _C4a);
    }
  }

  /**
   * A default constructor.  If GeodesicLine.Position is called on the
   * resulting object, it returns immediately (without doing any
   * calculations).  The object can be set with a call to Geodesic.Line.
   * Use Init() to test whether object is still in this uninitialized state.
   **********************************************************************/
  public GeodesicLine() { _caps = 0; }
  ///@}

  /**
   * \name Position in terms of distance
   **********************************************************************/
  ///@{

  /**
   * Compute the position of point 2 which is a distance \e s12 (meters) from
   * point 1.
   *
   * @param s12 distance between point 1 and point 2 (meters); it can be
   *   negative.
   * @return a GeodesicData object with the following fields: \e lat1, \e lon1,
   *   \e azi1, \e lat2, \e lon2, \e azi2, \e s12, \e a12.  Some of these
   *   results may be missing if the GeodesicLine did not include the relevant
   *   capability.
   *
   * The values of \e lon2 and \e azi2 returned are in the range
   * [&minus;180&deg;, 180&deg;).
   *
   * The GeodesicLine object \e must have been constructed with \e caps |=
   * GeodesicLine.DISTANCE_IN; otherwise Double.NaN is returned and no
   * parameters are set.
   **********************************************************************/
  public GeodesicData Position(double s12) {
    return Position(false, s12,
                    GeodesicMask.LATITUDE | GeodesicMask.LONGITUDE |
                    GeodesicMask.AZIMUTH);
  }
  /**
   * Compute the position of point 2 which is a distance \e s12 (meters) from
   * point 1 and with a subset of the geodesic results returned.
   *
   * @param s12 distance between point 1 and point 2 (meters); it can be
   *   negative.
   * @param outmask a bitor'ed combination of GeodesicMask values
   *   specifying which results should be returned.
   * @return a GeodesicData object including the requested results.
   *
   * The GeodesicLine object \e must have been constructed with \e caps |=
   * GeodesicLine.DISTANCE_IN; otherwise Double.NaN is returned and no
   * parameters are set.  Requesting a value which the GeodesicLine object is
   * not capable of computing is not an error; Double.NaN is returned instead.
   **********************************************************************/
  public GeodesicData Position(double s12, int outmask) {
    return Position(false, s12, outmask);
  }
  ///@}

  /**
   * \name Position in terms of arc length.
   **********************************************************************/
  ///@{

  /**
   * Compute the position of point 2 which is an arc length \e a12 (degrees)
   * from point 1.
   *
   * @param a12 arc length between point 1 and point 2 (degrees); it can
   *   be negative.
   * @return a GeodesicData object with the following fields: \e lat1, \e lon1,
   *   \e azi1, \e lat2, \e lon2, \e azi2, \e s12, \e a12.  Some of these
   *   results may be missing if the GeodesicLine did not include the relevant
   *   capability.
   *
   * The values of \e lon2 and \e azi2 returned are in the range
   * [&minus;180&deg;, 180&deg;).
   *
   * The GeodesicLine object \e must have been constructed with \e caps |=
   * GeodesicLine.DISTANCE_IN; otherwise Double.NaN is returned and no
   * parameters are set.
   **********************************************************************/
  public GeodesicData ArcPosition(double a12) {
    return Position(true, a12,
                    GeodesicMask.LATITUDE | GeodesicMask.LONGITUDE |
                    GeodesicMask.AZIMUTH | GeodesicMask.DISTANCE);
  }
  /**
   * Compute the position of point 2 which is an arc length \e a12 (degrees)
   * from point 1 and with a subset of the geodesic results returned.
   *
   * @param a12 arc length between point 1 and point 2 (degrees); it can
   *   be negative.
   * @param outmask a bitor'ed combination of GeodesicMask values
   *   specifying which results should be returned.
   * @return a GeodesicData object giving \e lat1, \e lon2, \e azi2, and \e
   *   a12.
   *
   * The GeodesicLine object \e must have been constructed with \e caps |=
   * GeodesicLine.DISTANCE_IN; otherwise Double.NaN is returned and no
   * parameters are set.  Requesting a value which the GeodesicLine object is
   * not capable of computing is not an error; Double.NaN is returned instead.
   **********************************************************************/
  public GeodesicData ArcPosition(double a12, int outmask) {
    return Position(true, a12, outmask);
  }
  ///@}

  /**
   * \name The general position function.
   **********************************************************************/
  ///@{

  /**
   * The general position function.  GeodesicLine.Position and
   * GeodesicLine.ArcPosition are defined in terms of this function.
   *
   * @param arcmode boolean flag determining the meaning of the second
   *   parameter; if arcmode is false, then the GeodesicLine object must have
   *   been constructed with \e caps |= GeodesicLine.DISTANCE_IN.
   * @param s12_a12 if \e arcmode is false, this is the distance between
   *   point 1 and point 2 (meters); otherwise it is the arc length between
   *   point 1 and point 2 (degrees); it can be negative.
   * @param outmask a bitor'ed combination of GeodesicMask values
   *   specifying which results should be returned.
   * @return a GeodesicData object with the requested results.
   *
   * The GeodesicLine.mask values possible for \e outmask are
   * - \e outmask |= GeodesicLine.LATITUDE for the latitude \e lat2.
   * - \e outmask |= GeodesicLine.LONGITUDE for the latitude \e lon2.
   * - \e outmask |= GeodesicLine.AZIMUTH for the latitude \e azi2.
   * - \e outmask |= GeodesicLine.DISTANCE for the distance \e s12.
   * - \e outmask |= GeodesicLine.REDUCEDLENGTH for the reduced length \e
   *   m12.
   * - \e outmask |= GeodesicLine.GEODESICSCALE for the geodesic scales \e
   *   M12 and \e M21.
   * - \e outmask |= GeodesicLine.AREA for the area \e S12.
   * .
   * Requesting a value which the GeodesicLine object is not capable of
   * computing is not an error; Double.NaN is returned instead.
   **********************************************************************/
  public GeodesicData Position(boolean arcmode, double s12_a12,
                               int outmask) {
    outmask &= _caps & GeodesicMask.OUT_ALL;
    GeodesicData r = new GeodesicData();
    if (!( Init() &&
           (arcmode ||
            (_caps & GeodesicMask.DISTANCE_IN & GeodesicMask.OUT_ALL) != 0) ))
      // Uninitialized or impossible distance calculation requested
      return r;
    r.lat1 = _lat1; r.lon1 = _lon1; r.azi1 = _azi1;

    // Avoid warning about uninitialized B12.
    double sig12, ssig12, csig12, B12 = 0, AB1 = 0;
    if (arcmode) {
      // Interpret s12_a12 as spherical arc length
      r.a12 = s12_a12;
      sig12 = s12_a12 * GeoMath.degree;
      double s12a = Math.abs(s12_a12);
      s12a -= 180 * Math.floor(s12a / 180);
      ssig12 = s12a ==  0 ? 0 : Math.sin(sig12);
      csig12 = s12a == 90 ? 0 : Math.cos(sig12);
    } else {
      // Interpret s12_a12 as distance
      r.s12 = s12_a12;
      double
        tau12 = s12_a12 / (_b * (1 + _A1m1)),
        s = Math.sin(tau12),
        c = Math.cos(tau12);
      // tau2 = tau1 + tau12
      B12 = - Geodesic.SinCosSeries(true,
                                    _stau1 * c + _ctau1 * s,
                                    _ctau1 * c - _stau1 * s,
                                    _C1pa);
      sig12 = tau12 - (B12 - _B11);
      r.a12 = sig12 / GeoMath.degree;
      ssig12 = Math.sin(sig12); csig12 = Math.cos(sig12);
      if (Math.abs(_f) > 0.01) {
        // Reverted distance series is inaccurate for |f| > 1/100, so correct
        // sig12 with 1 Newton iteration.  The following table shows the
        // approximate maximum error for a = WGS_a() and various f relative to
        // GeodesicExact.
        //     erri = the error in the inverse solution (nm)
        //     errd = the error in the direct solution (series only) (nm)
        //     errda = the error in the direct solution (series + 1 Newton) (nm)
        //
        //       f     erri  errd errda
        //     -1/5    12e6 1.2e9  69e6
        //     -1/10  123e3  12e6 765e3
        //     -1/20   1110 108e3  7155
        //     -1/50  18.63 200.9 27.12
        //     -1/100 18.63 23.78 23.37
        //     -1/150 18.63 21.05 20.26
        //      1/150 22.35 24.73 25.83
        //      1/100 22.35 25.03 25.31
        //      1/50  29.80 231.9 30.44
        //      1/20   5376 146e3  10e3
        //      1/10  829e3  22e6 1.5e6
        //      1/5   157e6 3.8e9 280e6
        double
          ssig2 = _ssig1 * csig12 + _csig1 * ssig12,
          csig2 = _csig1 * csig12 - _ssig1 * ssig12;
        B12 = Geodesic.SinCosSeries(true, ssig2, csig2, _C1a);
        double serr = (1 + _A1m1) * (sig12 + (B12 - _B11)) - s12_a12 / _b;
        sig12 = sig12 - serr / Math.sqrt(1 + _k2 * GeoMath.sq(ssig2));
        ssig12 = Math.sin(sig12); csig12 = Math.cos(sig12);
        // Update B12 below
      }
    }

    double omg12, lam12, lon12;
    double ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
    csig2 = _csig1 * csig12 - _ssig1 * ssig12;
    double dn2 = Math.sqrt(1 + _k2 * GeoMath.sq(ssig2));
    if ((outmask & (GeodesicMask.DISTANCE | GeodesicMask.REDUCEDLENGTH |
                    GeodesicMask.GEODESICSCALE)) != 0) {
      if (arcmode || Math.abs(_f) > 0.01)
        B12 = Geodesic.SinCosSeries(true, ssig2, csig2, _C1a);
      AB1 = (1 + _A1m1) * (B12 - _B11);
    }
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = _calp0 * ssig2;
    // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
    cbet2 = GeoMath.hypot(_salp0, _calp0 * csig2);
    if (cbet2 == 0)
      // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
      cbet2 = csig2 = Geodesic.tiny_;
    // tan(omg2) = sin(alp0) * tan(sig2)
    somg2 = _salp0 * ssig2; comg2 = csig2;  // No need to normalize
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = _salp0; calp2 = _calp0 * csig2; // No need to normalize
    // omg12 = omg2 - omg1
    omg12 = Math.atan2(somg2 * _comg1 - comg2 * _somg1,
                  comg2 * _comg1 + somg2 * _somg1);

    if ((outmask & GeodesicMask.DISTANCE) != 0 && arcmode)
      r.s12 = _b * ((1 + _A1m1) * sig12 + AB1);

    if ((outmask & GeodesicMask.LONGITUDE) != 0) {
      lam12 = omg12 + _A3c *
        ( sig12 + (Geodesic.SinCosSeries(true, ssig2, csig2, _C3a)
                   - _B31));
      lon12 = lam12 / GeoMath.degree;
      // Use GeoMath.AngNormalize2 because longitude might have wrapped multiple
      // times.
      lon12 = GeoMath.AngNormalize2(lon12);
      r.lon2 = GeoMath.AngNormalize(_lon1 + lon12);
    }

    if ((outmask & GeodesicMask.LATITUDE) != 0)
      r.lat2 = Math.atan2(sbet2, _f1 * cbet2) / GeoMath.degree;

    if ((outmask & GeodesicMask.AZIMUTH) != 0)
      // minus signs give range [-180, 180). 0- converts -0 to +0.
      r.azi2 = 0 - Math.atan2(-salp2, calp2) / GeoMath.degree;

    if ((outmask &
         (GeodesicMask.REDUCEDLENGTH | GeodesicMask.GEODESICSCALE)) != 0) {
      double
        B22 = Geodesic.SinCosSeries(true, ssig2, csig2, _C2a),
        AB2 = (1 + _A2m1) * (B22 - _B21),
        J12 = (_A1m1 - _A2m1) * sig12 + (AB1 - AB2);
      if ((outmask & GeodesicMask.REDUCEDLENGTH) != 0)
        // Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
        // accurate cancellation in the case of coincident points.
        r.m12 = _b * ((dn2 * (_csig1 * ssig2) - _dn1 * (_ssig1 * csig2))
                    - _csig1 * csig2 * J12);
      if ((outmask & GeodesicMask.GEODESICSCALE) != 0) {
        double t = _k2 * (ssig2 - _ssig1) * (ssig2 + _ssig1) / (_dn1 + dn2);
        r.M12 = csig12 + (t *  ssig2 -  csig2 * J12) * _ssig1 / _dn1;
        r.M21 = csig12 - (t * _ssig1 - _csig1 * J12) *  ssig2 /  dn2;
      }
    }

    if ((outmask & GeodesicMask.AREA) != 0) {
      double
        B42 = Geodesic.SinCosSeries(false, ssig2, csig2, _C4a);
      double salp12, calp12;
      if (_calp0 == 0 || _salp0 == 0) {
        // alp12 = alp2 - alp1, used in atan2 so no need to normalized
        salp12 = salp2 * _calp1 - calp2 * _salp1;
        calp12 = calp2 * _calp1 + salp2 * _salp1;
        // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
        // salp12 = -0 and alp12 = -180.  However this depends on the sign being
        // attached to 0 correctly.  The following ensures the correct behavior.
        if (salp12 == 0 && calp12 < 0) {
          salp12 = Geodesic.tiny_ * _calp1;
          calp12 = -1;
        }
      } else {
        // tan(alp) = tan(alp0) * sec(sig)
        // tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
        // = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
        // If csig12 > 0, write
        //   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
        // else
        //   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
        // No need to normalize
        salp12 = _calp0 * _salp0 *
          (csig12 <= 0 ? _csig1 * (1 - csig12) + ssig12 * _ssig1 :
           ssig12 * (_csig1 * ssig12 / (1 + csig12) + _ssig1));
        calp12 = GeoMath.sq(_salp0) + GeoMath.sq(_calp0) * _csig1 * csig2;
      }
      r.S12 = _c2 * Math.atan2(salp12, calp12) + _A4 * (B42 - _B41);
    }

    return r;
  }

  ///@}

  /**
   * \name Inspector functions
   **********************************************************************/
  ///@{

  /**
   * @return true if the object has been initialized.
   **********************************************************************/
  public boolean Init() { return _caps != 0; }

  /**
   * @return \e lat1 the latitude of point 1 (degrees).
   **********************************************************************/
  public double Latitude()
  { return Init() ? _lat1 : Double.NaN; }

  /**
   * @return \e lon1 the longitude of point 1 (degrees).
   **********************************************************************/
  public double Longitude()
  { return Init() ? _lon1 : Double.NaN; }

  /**
   * @return \e azi1 the azimuth (degrees) of the geodesic line at point 1.
   **********************************************************************/
  public double Azimuth()
  { return Init() ? _azi1 : Double.NaN; }

  /**
   * @return \e azi0 the azimuth (degrees) of the geodesic line as it crosses
   * the equator in a northward direction.
   **********************************************************************/
  public double EquatorialAzimuth() {
    return Init() ?
      Math.atan2(_salp0, _calp0) / GeoMath.degree : Double.NaN;
  }

  /**
   * @return \e a1 the arc length (degrees) between the northward equatorial
   * crossing and point 1.
   **********************************************************************/
  public double EquatorialArc() {
    return Init() ?
      Math.atan2(_ssig1, _csig1) / GeoMath.degree : Double.NaN;
  }

  /**
   * @return \e a the equatorial radius of the ellipsoid (meters).  This is
   *   the value inherited from the Geodesic object used in the constructor.
   **********************************************************************/
  public double MajorRadius()
  { return Init() ? _a : Double.NaN; }

  /**
   * @return \e f the flattening of the ellipsoid.  This is the value
   *   inherited from the Geodesic object used in the constructor.
   **********************************************************************/
  public double Flattening()
  { return Init() ? _f : Double.NaN; }

  /**
   * @return \e caps the computational capabilities that this object was
   *   constructed with.  LATITUDE and AZIMUTH are always included.
   **********************************************************************/
  public int Capabilities() { return _caps; }

  /**
   * @param testcaps a set of bitor'ed GeodesicMask values.
   * @return true if the GeodesicLine object has all these capabilities.
   **********************************************************************/
  public boolean Capabilities(int testcaps) {
    testcaps &= GeodesicMask.OUT_ALL;
    return (_caps & testcaps) == testcaps;
  }
  ///@}

}

/* This is a reformulation of the geodesic problem.  The notation is as
 * follows:
 * - at a general point (no suffix or 1 or 2 as suffix)
 *   - phi = latitude
 *   - beta = latitude on auxiliary sphere
 *   - omega = longitude on auxiliary sphere
 *   - lambda = longitude
 *   - alpha = azimuth of great circle
 *   - sigma = arc length along great circle
 *   - s = distance
 *   - tau = scaled distance (= sigma at multiples of pi/2)
 * - at northwards equator crossing
 *   - beta = phi = 0
 *   - omega = lambda = 0
 *   - alpha = alpha0
 *   - sigma = s = 0
 * - a 12 suffix means a difference, e.g., s12 = s2 - s1.
 * - s and c prefixes mean sin and cos
 **********************************************************************/

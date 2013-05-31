/**
 * @file GeodesicMask.java
 * @brief Implementation of the GeographicLib.GeodesicMask class
 *
 * Copyright (c) Charles Karney (2013) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/
package GeographicLib;
/**
 * @brief Bit masks for what geodesic calculations to do
 *
 * These masks do double duty.  They specify (via the \e outmask parameter)
 * which results to return in the GeodesicData object returned by the general
 * routines Geodesic.Direct and Geodesic.Inverse routines.  They also signify
 * (via the \e caps paramter) to the GeodesicLine.GeodesicLine constructor and
 * to Geodesic.Line what capabilities should be included in the GeodesicLine
 * object.  The GeodesicData object always includes the parameters provided to
 * Geodesic.Direct and Geodesic.Inverse and it always includes the field \e
 * a12.
 **********************************************************************/
public class GeodesicMask {
  /// @cond SKIP
  protected static final int CAP_NONE = 0;
  protected static final int CAP_C1   = 1<<0;
  protected static final int CAP_C1p  = 1<<1;
  protected static final int CAP_C2   = 1<<2;
  protected static final int CAP_C3   = 1<<3;
  protected static final int CAP_C4   = 1<<4;
  protected static final int CAP_ALL  = 0x1F;
  protected static final int OUT_ALL  = 0x7F80;
  /// @endcond

  /**
   * No capabilities, no output.
   * @hideinitializer
   **********************************************************************/
  public static final int NONE          = 0;
  /**
   * Calculate latitude \e lat2.  (It's not necessary to include this as a
   * capability to GeodesicLine because this is included by default.)
   * @hideinitializer
   **********************************************************************/
  public static final int LATITUDE      = 1<<7  | CAP_NONE;
  /**
   * Calculate longitude \e lon2.
   * @hideinitializer
   **********************************************************************/
  public static final int LONGITUDE     = 1<<8  | CAP_C3;
  /**
   * Calculate azimuths \e azi1 and \e azi2.  (It's not necessary to
   * include this as a capability to GeodesicLine because this is included
   * by default.)
   * @hideinitializer
   **********************************************************************/
  public static final int AZIMUTH       = 1<<9  | CAP_NONE;
  /**
   * Calculate distance \e s12.
   * @hideinitializer
   **********************************************************************/
  public static final int DISTANCE      = 1<<10 | CAP_C1;
  /**
   * Allow distance \e s12 to be used as input in the direct geodesic
   * problem.
   * @hideinitializer
   **********************************************************************/
  public static final int DISTANCE_IN   = 1<<11 | CAP_C1 | CAP_C1p;
  /**
   * Calculate reduced length \e m12.
   * @hideinitializer
   **********************************************************************/
  public static final int REDUCEDLENGTH = 1<<12 | CAP_C1 | CAP_C2;
  /**
   * Calculate geodesic scales \e M12 and \e M21.
   * @hideinitializer
   **********************************************************************/
  public static final int GEODESICSCALE = 1<<13 | CAP_C1 | CAP_C2;
  /**
   * Calculate area \e S12.
   * @hideinitializer
   **********************************************************************/
  public static final int AREA          = 1<<14 | CAP_C4;
  /**
   * All capabilities.  Calculate everything.
   * @hideinitializer
   **********************************************************************/
  public static final int ALL           = OUT_ALL| CAP_ALL;
}

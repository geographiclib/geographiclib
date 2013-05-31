/**
 * @file Constants.java
 * @brief Implementation of the GeographicLib.Constants class
 *
 * Copyright (c) Charles Karney (2013) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/
package GeographicLib;

/**
 * @brief %Constants needed by %GeographicLib
 *
 * Define constants specifying the WGS84 ellipsoid.
 ***********************************************************************/
public class Constants {
  /**
   * The equatorial radius of WGS84 ellipsoid (6378137 m).
   **********************************************************************/
  public static final double WGS84_a = 6378137;
  /**
   * The flattening of WGS84 ellipsoid (1/298.257223563).
   **********************************************************************/
  public static final double WGS84_f = 1/298.257223563;
}

/**
 * @file GeographicErr.java
 * @brief Implementation of the net.sf.geographiclib.GeographicErr class
 *
 * Copyright (c) Charles Karney (2013) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/
package net.sf.geographiclib;
/**
 * @brief Exception handling for GeographicLib
 *
 * A class to handle exceptions.  It's derived from RuntimeException so it
 * can be caught by the usual catch clauses.
 **********************************************************************/
public class GeographicErr extends RuntimeException {
  /**
   * Constructor
   *
   * @param msg a string message, which is accessible in the catch
   *   clause via getMessage().
   **********************************************************************/
  public GeographicErr(String msg) { super(msg); }
}

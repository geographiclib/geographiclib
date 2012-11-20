/**
 * \file PolygonArea.hpp
 * \brief Header for GeographicLib::PolygonArea class
 *
 * Copyright (c) Charles Karney (2010-2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_POLYGONAREA_HPP)
#define GEOGRAPHICLIB_POLYGONAREA_HPP 1

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Accumulator.hpp>

namespace GeographicLib {

  /**
   * \brief Polygon areas
   *
   * This computes the area of a geodesic polygon using the method given
   * Section 6 of
   * - C. F. F. Karney,
   *   <a href="http://dx.doi.org/10.1007/s00190-012-0578-z">
   *   Algorithms for geodesics</a>,
   *   J. Geodesy, 2012;
   *   DOI: <a href="http://dx.doi.org/10.1007/s00190-012-0578-z">
   *   10.1007/s00190-012-0578-z</a>;
   *   addenda <a href="http://geographiclib.sf.net/geod-addenda.html">
   *   geod-addenda.html</a>.
   *
   * This class lets you add vertices one at a time to the polygon.  The area
   * and perimeter are accumulated in two times the standard floating point
   * precision to guard against the loss of accuracy with many-sided polygons.
   * At any point you can ask for the perimeter and area so far.  There's an
   * option to treat the points as defining a polyline instead of a polygon; in
   * that case, only the perimeter is computed.
   *
   * Example of use:
   * \include example-PolygonArea.cpp
   *
   * <a href="Planimeter.1.html">Planimeter</a> is a command-line utility
   * providing access to the functionality of PolygonArea.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT PolygonArea {
  private:
    typedef Math::real real;
    Geodesic _earth;
    real _area0;                // Full ellipsoid area
    bool _polyline;             // Assume polyline (don't close and skip area)
    unsigned _mask;
    unsigned _num;
    int _crossings;
    Accumulator<real> _areasum, _perimetersum;
    real _lat0, _lon0, _lat1, _lon1;
    // Copied from Geodesic class (now the Math class)
    static inline real AngNormalize(real x) throw() {
      // Place angle in [-180, 180).  Assumes x is in [-540, 540).
      //
      // g++ 4.4.4 holds a temporary in an extended register causing an error
      // with the triangle 89,0.1;89,90.1;89,-179.9.  The volatile declaration
      // fixes this.  (The bug probably triggered because transit and
      // AngNormalize are inline functions.  So don't port this change over to
      // Geodesic.hpp.)
      volatile real y = x;
      return y >= 180 ? y - 360 : (y < -180 ? y + 360 : y);
    }
    static inline int transit(real lon1, real lon2) {
      // Return 1 or -1 if crossing prime meridian in east or west direction.
      // Otherwise return zero.
      lon1 = AngNormalize(lon1);
      lon2 = AngNormalize(lon2);
      // treat lon12 = -180 as an eastward geodesic, so convert to 180.
      real lon12 = -AngNormalize(lon1 - lon2); // In (-180, 180]
      int cross =
        lon1 < 0 && lon2 >= 0 && lon12 > 0 ? 1 :
        (lon2 < 0 && lon1 >= 0 && lon12 < 0 ? -1 : 0);
      return cross;
    }
  public:

    /**
     * Constructor for PolygonArea.
     *
     * @param[in] earth the Geodesic object to use for geodesic calculations.
     *   By default this uses the WGS84 ellipsoid.
     * @param[in] polyline if true that treat the points as defining a polyline
     *   instead of a polygon (default = false).
     **********************************************************************/
    PolygonArea(const Geodesic& earth, bool polyline = false) throw()
      : _earth(earth)
      , _area0(_earth.EllipsoidArea())
      , _polyline(polyline)
      , _mask(Geodesic::DISTANCE | (_polyline ? 0 : Geodesic::AREA))
    {
      Clear();
    }

    /**
     * Clear PolygonArea, allowing a new polygon to be started.
     **********************************************************************/
    void Clear() throw() {
      _num = 0;
      _crossings = 0;
      _areasum = 0;
      _perimetersum = 0;
      _lat0 = _lon0 = _lat1 = _lon1 = 0;
    }

    /**
     * Add a point to the polygon or polyline.
     *
     * @param[in] lat the latitude of the point (degrees).
     * @param[in] lon the latitude of the point (degrees).
     *
     * \e lat should be in the range [&minus;90&deg;, 90&deg;] and \e
     * lon should be in the range [&minus;540&deg;, 540&deg;).
     **********************************************************************/
    void AddPoint(real lat, real lon) throw();

    /**
     * Return the results so far.
     *
     * @param[in] reverse if true then clockwise (instead of counter-clockwise)
     *   traversal counts as a positive area.
     * @param[in] sign if true then return a signed result for the area if
     *   the polygon is traversed in the "wrong" direction instead of returning
     *   the area for the rest of the earth.
     * @param[out] perimeter the perimeter of the polygon or length of the
     *   polyline (meters).
     * @param[out] area the area of the polygon (meters<sup>2</sup>); only set
     *   if polyline is false in the constructor.
     * @return the number of points.
     **********************************************************************/
    unsigned Compute(bool reverse, bool sign,
                     real& perimeter, real& area) const throw();

    /**
     * Return the results assuming a tentative final test point is added;
     * however, the data for the test point is not saved.  This lets you report
     * a running result for the perimeter and area as the user moves the mouse
     * cursor.  Ordinary floating point arithmetic is used to accumulate the
     * data for the test point; thus the area and perimeter returned are less
     * accurate than if AddPoint and Compute are used.
     *
     * @param[in] lat the latitude of the test point (degrees).
     * @param[in] lon the longitude of the test point (degrees).
     * @param[in] reverse if true then clockwise (instead of counter-clockwise)
     *   traversal counts as a positive area.
     * @param[in] sign if true then return a signed result for the area if
     *   the polygon is traversed in the "wrong" direction instead of returning
     *   the area for the rest of the earth.
     * @param[out] perimeter the approximate perimeter of the polygon or length
     *   of the polyline (meters).
     * @param[out] area the approximate area of the polygon
     *   (meters<sup>2</sup>); only set if polyline is false in the
     *   constructor.
     * @return the number of points.
     *
     * \e lat should be in the range [&minus;90&deg;, 90&deg;] and \e
     * lon should be in the range [&minus;540&deg;, 540&deg;).
     **********************************************************************/
    unsigned TestCompute(real lat, real lon, bool reverse, bool sign,
                         real& perimeter, real& area) const throw();

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value inherited from the Geodesic object used in the constructor.
     **********************************************************************/

    Math::real MajorRadius() const throw() { return _earth.MajorRadius(); }

    /**
     * @return \e f the flattening of the ellipsoid.  This is the value
     *   inherited from the Geodesic object used in the constructor.
     **********************************************************************/
    Math::real Flattening() const throw() { return _earth.Flattening(); }
    ///@}
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_POLYGONAREA_HPP

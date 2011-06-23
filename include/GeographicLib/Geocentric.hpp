/**
 * \file Geocentric.hpp
 * \brief Header for GeographicLib::Geocentric class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010, 2011) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEOCENTRIC_HPP)
#define GEOGRAPHICLIB_GEOCENTRIC_HPP "$Id$"

#include <vector>
#include <algorithm>
#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief %Geocentric coordinates
   *
   * Convert between geodetic coordinates latitude = \e lat, longitude = \e
   * lon, height = \e h (measured vertically from the surface of the ellipsoid)
   * to geocentric coordinates (\e x, \e y, \e z).  The origin of geocentric
   * coordinates is at the center of the earth.  The \e z axis goes thru the
   * north pole, \e lat = 90<sup>o</sup>.  The \e x axis goes thru \e lat = 0,
   * \e lon = 0.  %Geocentric coordinates are also known as earth centered,
   * earth fixed (ECEF) coordinates.
   *
   * The conversion from geographic to geocentric coordinates is
   * straightforward.  For the reverse transformation we use
   * - H. Vermeille,
   *   <a href="http://dx.doi.org/10.1007/s00190-002-0273-6"> Direct
   *   transformation from geocentric coordinates to geodetic coordinates</a>,
   *   J. Geodesy 76, 451&ndash;454 (2002).
   * .
   * Several changes have been made to ensure that the method returns accurate
   * results for all finite inputs (even if \e h is infinite).  The changes are
   * described in Appendix B of
   * - C. F. F. Karney,
   *   <a href="http://arxiv.org/abs/1102.1215v1">Geodesics
   *   on an ellipsoid of revolution</a>,
   *   Feb. 2011;
   *   preprint
   *   <a href="http://arxiv.org/abs/1102.1215v1">arxiv:1102.1215v1</a>.
   * .
   * See \ref geocentric for more information.
   *
   * The errors in these routines are close to round-off.  Specifically, for
   * points within 5000 km of the surface of the ellipsoid (either inside or
   * outside the ellipsoid), the error is bounded by 7 nm for the WGS84
   * ellipsoid.  See \ref geocentric for further information on the errors.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT Geocentric {
  private:
    typedef Math::real real;
    friend class LocalCartesian;
    static const size_t dim_ = 3;
    static const size_t dim2_ = dim_ * dim_;
    const real _a, _f, _r, _e2, _e2m, _e2a, _e4a, _maxrad;
    // Actually this can be static because it doesn't depend on the ellipsoid.
    // But let's be more general than that.
    void Rotation(real sphi, real cphi, real slam, real clam,
                  real M[dim2_]) const throw();
    void IntForward(real lat, real lon, real h, real& x, real& y, real& z,
                    real M[dim2_]) const throw();
    void IntReverse(real x, real y, real z, real& lat, real& lon, real& h,
                    real M[dim2_]) const throw();
  public:

    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters)
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.  If \e f > 1, set flattening
     *   to 1/\e f.
     *
     * An exception is thrown if either of the axes of the ellipsoid is
     * non-positive.
     **********************************************************************/
    Geocentric(real a, real f);

    /**
     * Convert from geodetic to geocentric coordinates.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] h height of point above the ellipsoid (meters).
     * @param[out] x geocentric coordinate (meters).
     * @param[out] y geocentric coordinate (meters).
     * @param[out] z geocentric coordinate (meters).
     *
     * \e lat should be in the range [-90, 90]; \e lon and \e lon0 should be in
     * the range [-180, 360].
     **********************************************************************/
    void Forward(real lat, real lon, real h, real& x, real& y, real& z)
      const throw() {
      IntForward(lat, lon, h, x, y, z, NULL);
    }

    /**
     * Convert from geodetic to geocentric coordinates and return rotation
     * matrix.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] h height of point above the ellipsoid (meters).
     * @param[out] x geocentric coordinate (meters).
     * @param[out] y geocentric coordinate (meters).
     * @param[out] z geocentric coordinate (meters).
     * @param[out] M if the length of the vector is 9, fill with the rotation
     *   matrix in row-major order.
     *
     * Pre-multiplying a unit vector in local cartesian coordinates (east,
     * north, up) by \e M transforms the vector to geocentric coordinates.
     **********************************************************************/
    void Forward(real lat, real lon, real h, real& x, real& y, real& z,
                 std::vector<real>& M)
      const throw() {
      if (M.end() == M.begin() + dim2_) {
        real t[dim2_];
        IntForward(lat, lon, h, x, y, z, t);
        copy(t, t + dim2_, M.begin());
      } else
        IntForward(lat, lon, h, x, y, z, NULL);
    }

    /**
     * Convert from geocentric to geodetic to coordinates.
     *
     * @param[in] x geocentric coordinate (meters).
     * @param[in] y geocentric coordinate (meters).
     * @param[in] z geocentric coordinate (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] h height of point above the ellipsoid (meters).
     *
     * In general there are multiple solutions and the result which maximizes
     * \e h is returned.  If there are still multiple solutions with different
     * latitutes (applies only if \e z = 0), then the solution with \e lat > 0
     * is returned.  If there are still multiple solutions with different
     * longitudes (applies only if \e x = \e y = 0) then \e lon = 0 is
     * returned.  The value of \e h returned satisfies \e h >= - \e a (1 - \e
     * e<sup>2</sup>) / sqrt(1 - \e e<sup>2</sup> sin<sup>2</sup>\e lat).  The
     * value of \e lon returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(real x, real y, real z, real& lat, real& lon, real& h)
      const throw() {
      IntReverse(x, y, z, lat, lon, h, NULL);
    }

    /**
     * Convert from geocentric to geodetic to coordinates.
     *
     * @param[in] x geocentric coordinate (meters).
     * @param[in] y geocentric coordinate (meters).
     * @param[in] z geocentric coordinate (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] h height of point above the ellipsoid (meters).
     * @param[out] M if the length of the vector is 9, fill with the rotation
     *   matrix in row-major order.
     *
     * Pre-multiplying a unit vector in geocentric coordinates by the transpose
     * of \e M transforms the vector to local cartesian coordinates (east,
     * north, up).
     **********************************************************************/
    void Reverse(real x, real y, real z, real& lat, real& lon, real& h,
                 std::vector<real>& M)
      const throw() {
      if (M.end() == M.begin() + dim2_) {
        real t[dim2_];
        IntReverse(x, y, z, lat, lon, h, t);
        copy(t, t + dim2_, M.begin());
      } else
        IntReverse(x, y, z, lat, lon, h, NULL);
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
     * @return \e f the  flattening of the ellipsoid.  This is the
     *   value used in the constructor. 
     **********************************************************************/
    Math::real Flattening() const throw() { return _f; }

    /**
     * <b>DEPRECATED</b>
     * @return \e r the inverse flattening of the ellipsoid.
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return _r; }
    ///@}

    /**
     * A global instantiation of Geocentric with the parameters for the WGS84
     * ellipsoid.
     **********************************************************************/
    static const Geocentric WGS84;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_GEOCENTRIC_HPP

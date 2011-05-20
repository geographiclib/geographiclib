/**
 * \file PolarStereographic.hpp
 * \brief Header for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP)
#define GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP "$Id: PolarStereographic.hpp 6867 2010-09-11 13:04:26Z karney $"

#include "GeographicLib/Constants.hpp"

namespace GeographicLib {

  /**
   * \brief Polar Stereographic Projection
   *
   * Implementation taken from the report,
   * - J. P. Snyder,
   *   <a href="http://pubs.er.usgs.gov/usgspubs/pp/pp1395"> Map Projections: A
   *   Working Manual</a>, USGS Professional Paper 1395 (1987),
   *   pp. 160&ndash;163.
   *
   * This is a straightforward implementation of the equations in Snyder except
   * that Newton's method is used to invert the projection.
   **********************************************************************/
  class PolarStereographic {
  private:
    typedef Math::real real;
    // _Cx used to be _C but g++ 3.4 has a macro of that name
    const real _a, _r, _f, _e2, _e, _e2m, _Cx, _c;
    real _k0;
    static const real tol, overflow;
    static const int numit = 5;
    static inline real sq(real x) throw() { return x * x; }
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
     * An exception is thrown if either of the axes of the ellipsoid is
     * not positive \e a or if \e k0 is not positive.
     **********************************************************************/
    PolarStereographic(real a, real r, real k0);

    /**
     * Set the scale for the projection.
     *
     * @param[in] lat (degrees).
     * @param[in] k scale at latitude \e lat (default 1).
     *
     * This allows a "latitude of true scale" to be specified.  An exception is
     * thrown if \e k is not positive.
     **********************************************************************/
    void SetScale(real lat, real k = real(1));

    /**
     * Forward projection, from geographic to polar stereographic.
     *
     * @param[in] northp the pole which is the center of projection (true means
     *   north, false means south).
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added.  \e lat should be in the range
     * (-90, 90] for \e northp = true and in the range [-90, 90) for \e northp
     * = false; \e lon should be in the range [-180, 360].
     **********************************************************************/
    void Forward(bool northp, real lat, real lon,
                 real& x, real& y, real& gamma, real& k) const throw();

    /**
     * Reverse projection, from polar stereographic to geographic.
     *
     * @param[in] northp the pole which is the center of projection (true means
     *   north, false means south).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added.  The value of \e lon returned is
     * in the range [-180, 180).
     **********************************************************************/
    void Reverse(bool northp, real x, real y,
                 real& lat, real& lon, real& gamma, real& k) const throw();

    /**
     * PolarStereographic::Forward without returning the convergence and scale.
     **********************************************************************/
    void Forward(bool northp, real lat, real lon,
                 real& x, real& y) const throw() {
      real gamma, k;
      Forward(northp, lat, lon, x, y, gamma, k);
    }

    /**
     * PolarStereographic::Reverse without returning the convergence and scale.
     **********************************************************************/
    void Reverse(bool northp, real x, real y,
                 real& lat, real& lon) const throw() {
      real gamma, k;
      Reverse(northp, x, y, lat, lon, gamma, k);
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
     * The central scale for the projection.  This is the value of \e k0 used
     * in the constructor and is the scale at the pole unless overridden by
     * PolarStereographic::SetScale.
     **********************************************************************/
    Math::real CentralScale() const throw() { return _k0; }
    ///@}

    /**
     * A global instantiation of PolarStereographic with the WGS84 ellipsoid
     * and the UPS scale factor.  However, unlike UPS, no false easting or
     * northing is added.
     **********************************************************************/
    const static PolarStereographic UPS;
  };

} // namespace GeographicLib

#endif

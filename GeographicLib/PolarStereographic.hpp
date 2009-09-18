/**
 * \file PolarStereographic.hpp
 * \brief Header for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP)
#define GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP "$Id$"

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
    typedef Math::real_t real_t;
    const real_t _a, _f, _k0, _e2, _e, _e2m, _c;
    static const real_t tol;
    static const int numit = 5;
    static inline real_t sq(real_t x) throw() { return x * x; }
    static inline real_t hypot(real_t x, real_t y) throw()
    { return Math::hypot(x, y); }
    static inline real_t atanh(real_t x) throw() { return Math::atanh(x); }
    // Return e * atanh(e * x) for f >= 0, else return
    // - sqrt(-e2) * atan( sqrt(-e2) * x) for f < 0
    inline real_t eatanhe(real_t x) const throw() {
      return _f >= 0 ? _e * atanh(_e * x) : - _e * atan(_e * x);
    }
  public:

    /**
     * Constructor for a ellipsoid radius \e a (meters), reciprocal flattening
     * \e r, and central scale factor \e k0.  Setting \e r <= 0 implies \e r =
     * inf or flattening = 0 (i.e., a sphere).
     **********************************************************************/
    PolarStereographic(Math::real_t a, Math::real_t r, Math::real_t k0) throw();

    /**
     * Convert from latitude \e lat (degrees) and longitude \e lon (degrees) to
     * polar stereographic easting \e x (meters) and northing \e y (meters).
     * The projection is about the pole given by \e northp (false means south,
     * true means north).  Also return the meridian convergence \e gamma
     * (degrees) and the scale \e k.  No false easting or northing is added.
     * \e lat should be in the range (-90, 90] for \e northp = true and in the
     * range [-90, 90) for \e northp = false; \e lon should be in the range
     * [-180, 360].
     **********************************************************************/
    void Forward(bool northp, Math::real_t lat, Math::real_t lon,
		 Math::real_t& x, Math::real_t& y,
		 Math::real_t& gamma, Math::real_t& k) const throw();

    /**
     * Convert from polar stereogrphic easting \e x (meters) and northing \e y
     * (meters) to latitude \e lat (degrees) and longitude \e lon (degrees) .
     * The hemisphere is given by \e northp (false means south, true means
     * north).  Also return the meridian convergence \e gamma (degrees) and the
     * scale \e k.  No false easting or northing is added.  The value of \e lon
     * returned is in the range [-180, 180).
     **********************************************************************/
    void Reverse(bool northp, Math::real_t x, Math::real_t y,
		 Math::real_t& lat, Math::real_t& lon,
		 Math::real_t& gamma, Math::real_t& k) const throw();

    /**
     * A global instantiation of PolarStereographic with the WGS84 ellipsoid
     * and the UPS scale factor.  However, unlike UPS, no false easting or
     * northing is added.
     **********************************************************************/
    const static PolarStereographic UPS;
  };

} // namespace GeographicLib

#endif

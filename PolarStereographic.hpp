/**
 * \file PolarStereographic.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(POLARSTEREOGRAPHIC_HPP)
#define POLARSTEREOGRAPHIC_HPP "$Id$"

namespace GeographicLib {

  class PolarStereographic {
  private:
    const double _a, _f, _k0, _e, _e2m, _c, _tol;
    const int _numit;
    static inline double sq(double x) { return x * x; }
  public:
    PolarStereographic(double a, double f, double k0);
    void Forward(bool northp, double lat, double lon,
		 double& x, double& y, double& gamma, double& k) const;
    void Reverse(bool northp, double x, double y,
		 double& lat, double& lon, double& gamma, double& k) const;
    // Specialization for WGS84 ellipsoid and UPS scale factor
    const static PolarStereographic UPS;
  };

} // namespace GeographicLib

#endif

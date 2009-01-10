/**
 * \file ECEF.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(ECEF_HPP)
#define ECEF_HPP "$Id$"

namespace GeographicLib {

    // Compute reverse using Pollard's latitude first method, J. Pollard,
    // Iterative vector methods for computing geodetic latitude and height from
    // rectangular coordinates, J. Geodesy 76, 36&ndash;40 (2002).  The advantage of
    // this method is that it's simple and convergence (via Newton's method) is
    // fast.
  class ECEF {
  private:
    const double _a, _f, _e2, _e12, _b, _tol;
    const int _numit;
    static inline double sq(double x) { return x * x; }
  public:
    ECEF(double a, double f);
    void Forward(double lat, double lon, double h,
		 double& x, double& y, double& z) const;
    void Reverse(double x, double y, double z,
		 double& lat, double& lon, double& h) const;
    // Specialization for WGS84 ellipsoid
    const static ECEF WGS84;
  };

} //namespace GeographicLib
#endif

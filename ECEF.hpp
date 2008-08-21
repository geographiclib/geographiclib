/**
 * \file ECEF.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(ECEF_HPP)
#define ECEF_HPP "$Id$"

namespace GeographicLib {

  class ECEF {
  private:
    const double _a, _f, _e2, _e12, _b, _tol;
    const int _numit;
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

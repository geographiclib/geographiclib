/**
 * \file LocalCartesian.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(LOCALCARTESIAN_HPP)
#define LOCALCARTESIAN_HPP "$Id$"

namespace GeographicLib {

  class LocalCartesian {
  private:
    double _x0, _y0, _z0,
      _rxx, _rxy, _rxz,
      _ryx, _ryy, _ryz,
      _rzx, _rzy, _rzz;
  public:
    LocalCartesian() { Reset(0.0, 0.0); }
    LocalCartesian(double lat0, double lon0) {
      Reset(lat0, lon0);
    }
    void Reset(double lat0, double lon0);
    void Forward(double lat, double lon, double h,
		 double& x, double& y, double& z) const;
    void Reverse(double x, double y, double z,
		 double& lat, double& lon, double& h) const;
  };

} // namespace GeographicLib

#endif

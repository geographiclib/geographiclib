/**
 * \file TransverseMercator.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(TRANSVERSEMERCATOR_HPP)
#define TRANSVERSEMERCATOR_HPP "$Id$"

#include <cmath>

namespace GeographicLib {

  class TransverseMercator {
  private:
    const double _a, _f, _k0, _e2, _e, _e2m, _e1,  _n, _a1,
      _h1, _h2, _h3, _h4,
      _h1p, _h2p, _h3p, _h4p,
      _tol;
    const int _numit;
#if defined(_MSC_VER)
    // These have poor relative accuracy near x = 0.  However, for mapping
    // applications, we only need good absolute accuracy.  
    // For small arguments we would write
    //
    // asinh(x) = asin(x) -x^3/3-5*x^7/56-63*x^11/1408-143*x^15/5120 ...
    // atanh(x) = atan(x) +2*x^3/3+2*x^7/7+2*x^11/11+2*x^15/15
    //
    // The accuracy of asinh is also bad for large negative arguments.  This is
    // easy to fix in the definition of asinh.  Instead we call these functions
    // with positive arguments and enforce the correct parity separately.
    static inline double asinh(double x) { return log(x + sqrt(1 + x * x)); }
    static inline double atanh(double x) { return log((1 + x)/(1 - x))/2; }
#endif
    double Convergence(double phi, double l) const;
    double Scale(double phi, double l) const;
  public:
    TransverseMercator(double a, double f, double k0);
    void Forward(double lon0, double lat, double lon,
		 double& x, double& y, double& gamma, double& k) const;
    void Reverse(double lon0, double x, double y,
		 double& lat, double& lon, double& gamma, double& k) const;
    // Specialization for WGS84 ellipsoid and UTM scale factor
    const static TransverseMercator UTM;
  };

} // namespace GeographicLib

#endif

/**
 * \file TransverseMercatorExact.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(TRANSVERSEMERCATOREXACT_HPP)
#define TRANSVERSEMERCATOREXACT_HPP "$Id$"

#include <cmath>
#include "EllipticFunction.hpp"

#if !defined(TM_MAXPOW)
#define TM_MAXPOW 6
#endif

namespace GeographicLib {

  class TransverseMercatorExact {
  private:
    static const double tol, tol1, taytol, ahypover;
    static const int numit = 10;
    const double _a, _f, _k0, _mu, _mv, _e;
    const bool _foldp;
    const EllipticFunction Eu, Ev;
    static inline double sq(double x) { return x * x; }
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
    static inline double asinh(double x) { return log(x + sqrt(1 + sq(x))); }
    static inline double atanh(double x) { return log((1 + x)/(1 - x))/2; }
#endif
    static double psi0(double phi);
    double psi(double phi) const;
    double psiinv(double psi) const;

    void zeta(double u, double snu, double cnu, double dnu,
	      double v, double snv, double cnv, double dnv,
	      double& psi, double& lam) const;

    void dwdzeta(double u, double snu, double cnu, double dnu,
		 double v, double snv, double cnv, double dnv,
		 double& du, double& dv) const;

    bool zetainv0(double psi, double lam, double& u, double& v) const;
    void zetainv(double psi, double lam, double& u, double& v) const;

    void sigma(double u, double snu, double cnu, double dnu,
	       double v, double snv, double cnv, double dnv,
	       double& xi, double& eta) const;

    void dwdsigma(double u, double snu, double cnu, double dnu,
		  double v, double snv, double cnv, double dnv,
		  double& du, double& dv) const;

    bool sigmainv0(double xi, double eta, double& u, double& v) const;
    void sigmainv(double xi, double eta, double& u, double& v) const;

    void Scale(double phi, double lam,
	       double snu, double cnu, double dnu,
	       double snv, double cnv, double dnv,
	       double& gamma, double& k) const;

  public:
    // With foldp, the "standard" conventions for cutting the TM space are
    // followed.  Forward can be called with any latitude and longitude then
    // produces the transformation shown in Lee Fig 46.  Reverse analytically
    // continues this in the +/- x direction (with a cut at y = 0).
    //
    // With not(foldp), the domains of latitude, longitude, x, and y are
    // restricted to
    //
    //   latitude in [0,90] and longitude-lon0 in [0,90]
    //   latitude in (-90,0] and longitude-lon0 in [90*(1-e),90]
    //
    //   x/(k0*a) in [0,inf) and y/(k0*a) in [0,E(e^2)]
    //   x/(k0*a) in [K(1-e^2)-E(1-e^2),inf) and y/(k0*a) in (-inf,0]
    //
    // This allows the multi-valued nature of the projection to be explored
    // using its symmetries.
    TransverseMercatorExact(double a, double f, double k0, bool foldp = true);
    void Forward(double lon0, double lat, double lon,
		 double& x, double& y, double& gamma, double& k) const;
    void Reverse(double lon0, double x, double y,
		 double& lat, double& lon, double& gamma, double& k) const;
    // Specialization for WGS84 ellipsoid and UTM scale factor
    const static TransverseMercatorExact UTM;
  };

} // namespace GeographicLib

#endif

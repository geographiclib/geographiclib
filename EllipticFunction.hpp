/**
 * \file EllipticFunction.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(ELLIPTICFUNCTION_HPP)
#define ELLIPTICFUNCTION_HPP "$Id$"


namespace GeographicLib {

  class EllipticFunction {
  private:
    static const double tol, tolRF, tolRD, tolRG0, tolJAC, tolJAC1;
    enum { num = 10 };		// Max depth required for sncndn.  Probably 5
				// is enough.
    static double RF(double x, double y, double z);
    static double RD(double x, double y, double z);
    static double RG0(double x, double y);
    const double _m, _m1;
    mutable bool _init;
    mutable double _kc, _ec, _kec;
    bool Init() const;
  public:
    EllipticFunction(double m);
    double m() const { return _m; } // Parameter
    double m1() const { return _m1; } // Complementary parameter
    // Complete integral of first kind, K(m)
    double K() const { _init || Init(); return _kc; }
    // Complete integral of second kind, E(m)
    double E() const { _init || Init(); return _ec; }
    // K(m) - E(m)
    double KE() const { _init || Init(); return _kec; }
    // Return sn(x|m), cn(x|m) dn(x|m)
    void sncndn(double x, double& sn, double& cn, double& dn) const;
    // Return incomplete integral of the second kind int dn(w)^2 (A+S 17.2.10)
    double E(double sn, double cn, double dn) const;
  };


} // namespace GeographicLib

#endif

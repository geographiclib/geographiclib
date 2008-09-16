/**
 * \file EllipticFunction.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(ELLIPTICFUNCTION_HPP)
#define ELLIPTICFUNCTION_HPP "$Id: TransverseMercator.hpp 6448 2008-09-08 02:09:51Z ckarney $"


namespace GeographicLib {

  class EllipticFunction {
  private:
    static const double tol, tolRF, tolRD, tolRG0, tolJAC, tolJAC1;
    enum { num = 10 };		// Max depth required for sncndn.  Probably 5
				// is enough.
    static double RF(double x, double y, double z);
    static double RD(double x, double y, double z);
    static double RG0(double x, double y);
    const double _m, _m1, _kc, _ec, _kec;
  public:
    EllipticFunction(double m);
    double m() const { return _m; } // Parameter
    double m1() const { return _m1; } // Complementary parameter
    double K() const { return _kc; }  // Complete integral of first kind, K(m)
    double E() const { return _ec; }  // Complete integral of second kind, E(m)
    double KE() const { return _kec; } // K(m) - E(m)
    // Return sn(x|m), cn(x|m) dn(x|m)
    void sncndn(double x, double& sn, double& cn, double& dn) const;
    // Return incomplete integral of the second kind int dn(w)^2 (A+S 17.2.10)
    double E(double sn, double cn, double dn) const;
  };


} // namespace GeographicLib

#endif

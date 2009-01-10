/**
 * \file EllipticFunction.hpp
 * \brief Header for GeographicLib::EllipticFunction class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(ELLIPTICFUNCTION_HPP)
#define ELLIPTICFUNCTION_HPP "$Id$"


namespace GeographicLib {

  /**
   * \brief Elliptic functions needed for TransverseMercatorExact
   *
   * This provides the subset of elliptic functions needed for
   * TransverseMercatorExact.  For a given ellipsoid, only parameters \e
   * e<sup>2</sup> and 1 - \e e<sup>2</sup> are needed.  This class taken the
   * parameter as a constructor parameters and caches the values of the
   * required complete integrals.  A method is provided for Jacobi elliptic
   * functions and for the incomplete elliptic integral of the second kind in
   * terms of the amplitude.
   *
   * The computation of the elliptic integrals uses the algorithms given in
   * - B. C. Carlson,
   *   Computation of elliptic integrals,
   *   Numerical Algorithms 10, 13&ndash;26 (1995).
   * .
   * The computation of the Jacobi elliptic functions uses the algorithm given
   * in
   * - Roland Bulirsch,
   *   Numerical Calculation of Elliptic Integrals and Elliptic Functions,
   *   Numericshe Mathematik 7, 78&ndash;90 (1965).
   * .
   * The notation follows Abramowitz and Stegun, Chapters 16 and 17.
   **********************************************************************/
  class EllipticFunction {
  private:
    static const double tol, tolRF, tolRD, tolRG0, tolJAC, tolJAC1;
    enum { num = 10 }; // Max depth required for sncndn.  Probably 5 is enough.
    static double RF(double x, double y, double z);
    static double RD(double x, double y, double z);
    static double RG0(double x, double y);
    const double _m, _m1;
    mutable bool _init;
    mutable double _kc, _ec, _kec;
    bool Init() const;
  public:
    /**
     * Constructor with parameter \e m.
     **********************************************************************/
    EllipticFunction(double m);
    /**
     * The parameter.
     **********************************************************************/
    double m() const { return _m; }
    /**
     * The complementary parameter.
     **********************************************************************/
    double m1() const { return _m1; }
    /**
     * The complete integral of first kind, K(m).
     **********************************************************************/
    double K() const { _init || Init(); return _kc; }
    /**
     * The complete integral of second kind, E(m).
     **********************************************************************/
    double E() const { _init || Init(); return _ec; }
    /**
     * The difference K(m) - E(m) (which can be computed directly).
     **********************************************************************/
    double KE() const { _init || Init(); return _kec; }
    /**
     * The Jacobi elliptic functions sn(x|m), cn(x|m), and dn(x|m) with
     * argument \e x.  The results are returned in \e sn, \e cn, and \e dn.
     **********************************************************************/
    void sncndn(double x, double& sn, double& cn, double& dn) const;
    /**
     * The incomplete integral of the second kind = int dn(w)^2 (A+S 17.2.10).
     * Instead of specifying the ampltiude phi, we provide \e sn = sin(phi), \e
     * cn = cos(phi), \e dn = sqrt(1 - m sin(phi)^2).
     **********************************************************************/
    double E(double sn, double cn, double dn) const;
  };


} // namespace GeographicLib

#endif

/**
 * \file EllipticFunction.hpp
 * \brief Header for GeographicLib::EllipticFunction class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP)
#define GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP "$Id$"

#include "GeographicLib/Constants.hpp"

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
   *   <a href="http://dx.doi.org/10.1007/BF02198293"> Computation of elliptic
   *   integrals</a>, Numerical Algorithms 10, 13&ndash;26 (1995).
   * .
   * The computation of the Jacobi elliptic functions uses the algorithm given
   * in
   * - R. Bulirsch,
   *   <a href="http://dx.doi.org/10.1007/BF01397975"> Numerical Calculation of
   *   Elliptic Integrals and Elliptic Functions</a>, Numericshe Mathematik 7,
   *   78&ndash;90 (1965).
   * .
   * The notation follows Abramowitz and Stegun, Chapters 16 and 17.
   **********************************************************************/
  class EllipticFunction {
  private:
    typedef Math::real_t real_t;
    static const real_t tol, tolRF, tolRD, tolRG0, tolJAC, tolJAC1;
    enum { num = 10 }; // Max depth required for sncndn.  Probably 5 is enough.
    static real_t RF(real_t x, real_t y, real_t z) throw();
    static real_t RD(real_t x, real_t y, real_t z) throw();
    static real_t RG0(real_t x, real_t y) throw();
    const real_t _m, _m1;
    mutable bool _init;
    mutable real_t _kc, _ec, _kec;
    bool Init() const throw();
  public:

    /**
     * Constructor with parameter \e m.
     **********************************************************************/
    explicit EllipticFunction(Math::real_t m) throw();

    /**
     * The parameter \e m.
     **********************************************************************/
    Math::real_t m() const throw() { return _m; }

    /**
     * The complementary parameter \e m' = (1 - \e m).
     **********************************************************************/
    Math::real_t m1() const throw() { return _m1; }

    /**
     * The complete integral of first kind, \e K(\e m).
     **********************************************************************/
    Math::real_t K() const throw() { _init || Init(); return _kc; }

    /**
     * The complete integral of second kind, \e E(\e m).
     **********************************************************************/
    Math::real_t E() const throw() { _init || Init(); return _ec; }

    /**
     * The difference \e K(\e m) - \e E(\e m) (which can be computed directly).
     **********************************************************************/
    Math::real_t KE() const throw() { _init || Init(); return _kec; }

    /**
     * The Jacobi elliptic functions sn(<i>x</i>|<i>m</i>),
     * cn(<i>x</i>|<i>m</i>), and dn(<i>x</i>|<i>m</i>) with argument \e x.
     * The results are returned in \e sn, \e cn, and \e dn.
     **********************************************************************/
    void sncndn(Math::real_t x,
                Math::real_t& sn, Math::real_t& cn, Math::real_t& dn)
      const throw();

    /**
     * The incomplete integral of the second kind = int dn(\e w)<sup>2</sup> \e
     * dw (A+S 17.2.10).  Instead of specifying the ampltiude \e phi, we
     * provide \e sn = sin(\e phi), \e cn = cos(\e phi), \e dn = sqrt(1 - \e m
     * sin<sup>2</sup>(\e phi)).
     **********************************************************************/
    Math::real_t E(Math::real_t sn, Math::real_t cn, Math::real_t dn)
      const throw();
  };


} // namespace GeographicLib

#endif

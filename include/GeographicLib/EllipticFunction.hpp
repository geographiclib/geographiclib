/**
 * \file EllipticFunction.hpp
 * \brief Header for GeographicLib::EllipticFunction class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP)
#define GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP "$Id: EllipticFunction.hpp 6866 2010-09-11 02:15:29Z karney $"

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
    typedef Math::real real;
    static const real tol, tolRF, tolRD, tolRG0, tolJAC, tolJAC1;
    enum { num = 10 }; // Max depth required for sncndn.  Probably 5 is enough.
    static real RF(real x, real y, real z) throw();
    static real RD(real x, real y, real z) throw();
    static real RG0(real x, real y) throw();
    const real _m, _m1;
    mutable bool _init;
    mutable real _kc, _ec, _kec;
    bool Init() const throw();
  public:

    /**
     * Constructor.
     *
     * @param[in] m the parameter which must lie in [0, 1].  (No checking
     *   is done.)
     **********************************************************************/
    explicit EllipticFunction(real m) throw();

    /**
     * @return the parameter \e m.
     **********************************************************************/
    Math::real m() const throw() { return _m; }

    /**
     * @return the complementary parameter \e m' = (1 - \e m).
     **********************************************************************/
    Math::real m1() const throw() { return _m1; }

    /**
     * @return the complete integral of first kind, \e K(\e m).
     **********************************************************************/
    Math::real K() const throw() { _init || Init(); return _kc; }

    /**
     * @return the complete integral of second kind, \e E(\e m).
     **********************************************************************/
    Math::real E() const throw() { _init || Init(); return _ec; }

    /**
     * @return the difference \e K(\e m) - \e E(\e m) (which can be computed
     *   directly).
     **********************************************************************/
    Math::real KE() const throw() { _init || Init(); return _kec; }

    /**
     * The Jacobi elliptic functions.
     *
     * @param[in] x the argument.
     * @param[out] sn sn(<i>x</i>|<i>m</i>).
     * @param[out] cn cn(<i>x</i>|<i>m</i>).
     * @param[out] dn dn(<i>x</i>|<i>m</i>).
     **********************************************************************/
    void sncndn(real x, real& sn, real& cn, real& dn) const throw();

    /**
     * The incomplete integral of the second kind.
     *
     * @param[in] phi
     * @return int sqrt(1 -  \e m sin<sup>2</sup>(\e phi)) \e dphi.
     **********************************************************************/
    Math::real E(real phi) const throw();

    /**
     * The incomplete integral of the second kind in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn
     * @param[in] cn
     * @param[in] dn
     * @return int dn(\e w)<sup>2</sup> \e dw (A+S 17.2.10).
     *
     * Instead of specifying the ampltiude \e phi, we provide \e sn = sin(\e
     * phi), \e cn = cos(\e phi), \e dn = sqrt(1 - \e m sin<sup>2</sup>(\e
     * phi)).
     **********************************************************************/
    Math::real E(real sn, real cn, real dn) const throw();
  };


} // namespace GeographicLib

#endif

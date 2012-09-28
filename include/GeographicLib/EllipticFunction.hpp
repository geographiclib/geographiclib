/**
 * \file EllipticFunction.hpp
 * \brief Header for GeographicLib::EllipticFunction class
 *
 * Copyright (c) Charles Karney (2008-2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP)
#define GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP 1

#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief Elliptic functions needed for TransverseMercatorExact
   *
   * This provides the subset of elliptic functions needed for
   * TransverseMercatorExact.  For a given ellipsoid, only parameters
   * <i>e</i><sup>2</sup> and 1 &minus; <i>e</i><sup>2</sup> are needed.  This
   * class taken the parameter as a constructor parameters and caches the
   * values of the required complete integrals.  A method is provided for
   * Jacobi elliptic functions and for the incomplete elliptic integral of the
   * second kind in terms of the amplitude.
   *
   * The computation of the elliptic integrals uses the algorithms given in
   * - B. C. Carlson,
   *   <a href="http://dx.doi.org/10.1007/BF02198293"> Computation of elliptic
   *   integrals</a>, Numerical Algorithms 10, 13--26 (1995).
   * .
   * The computation of the Jacobi elliptic functions uses the algorithm given
   * in
   * - R. Bulirsch,
   *   <a href="http://dx.doi.org/10.1007/BF01397975"> Numerical Calculation of
   *   Elliptic Integrals and Elliptic Functions</a>, Numericshe Mathematik 7,
   *   78--90 (1965).
   * .
   * The notation follows Abramowitz and Stegun, Chapters 16 and 17.
   *
   * Example of use:
   * \include example-EllipticFunction.cpp
   **********************************************************************/
  class GEOGRAPHIC_EXPORT EllipticFunction {
  private:
    typedef Math::real real;
    static const real tol_;
    static const real tolRF_;
    static const real tolRD_;
    static const real tolRG0_;
    static const real tolJAC_;
    enum { num_ = 13 }; // Max depth required for sncndn.  Probably 5 is enough.
    real _m, _m1;
    mutable bool _init;
    mutable real _kc, _ec, _dc;
    bool Init() const throw();
  public:
    /** \name Constructor
     **********************************************************************/
    ///@{
    /**
     * Constructor.
     *
     * @param[in] m the parameter which must lie in [0, 1].  (No checking
     *   is done.)
     **********************************************************************/
    explicit EllipticFunction(real m) throw();
    ///@}

    /** \name Inspector functions.
     **********************************************************************/
    ///@{
    /**
     * @return the parameter \e m.
     **********************************************************************/
    Math::real m() const throw() { return _m; }

    /**
     * @return the complementary parameter \e m' = (1 &minus; \e m).
     **********************************************************************/
    Math::real m1() const throw() { return _m1; }
    ///@}

    /** \name Complete elliptic integrals.
     **********************************************************************/
    ///@{
    /**
     * @return the complete integral of first kind, \e K(\e m)
     *
     * \e K(\e m) is defined in http://dlmf.nist.gov/19.2.E4
     * \f[
     *   K(m) = \int_0^{\pi/2} \frac1{\sqrt{1-m\sin^2\phi}}\,d\phi.
     * \f]
     **********************************************************************/
    Math::real K() const throw() { _init || Init(); return _kc; }

    /**
     * @return the complete integral of second kind, \e E(\e m)
     *
     * \e E(\e m) is defined in http://dlmf.nist.gov/19.2.E5
     * \f[
     *   E(m) = \int_0^{\pi/2} \sqrt{1-m\sin^2\phi}\,d\phi.
     * \f]
     *
     * This function returns the correct result even when \e m is negative.
     **********************************************************************/
    Math::real E() const throw() { _init || Init(); return _ec; }

    /**
     * @return the complete integral \e D(\e m); from
     *
     * \e D(\e m) is defined in http://dlmf.nist.gov/19.2.E6
     * \f[
     *   D(m) = \int_0^{\pi/2} \frac{\sin^2\phi}{\sqrt{1-m\sin^2\phi}}\,d\phi.
     * \f]
     **********************************************************************/
    Math::real D() const throw() { _init || Init(); return _dc; }

    /**
     * @return the difference \e K(\e m) - \e E(\e m).
     **********************************************************************/
    Math::real KE() const throw() { _init || Init(); return _m * _dc; }
    ///@}

    /** \name Incomplete elliptic integrals.
     **********************************************************************/
    ///@{
    /**
     * The incomplete integral of the second kind.
     *
     * @param[in] phi
     * @return &int; sqrt(1 &minus; \e m sin<sup>2</sup>&phi;) d&phi;.
     **********************************************************************/
    Math::real E(real phi) const throw();

    /**
     * The incomplete integral of the second kind with the argument given in
     * degrees.
     *
     * @param[in] ang in <i>degrees</i>.
     * @return &int; sqrt(1 &minus; \e m sin<sup>2</sup>&phi;) d&phi;.
     *
     * \e ang must lie in [&minus;90&deg;, 90&deg;].  This function
     * returns the correct result even when \e m is negative.
     **********************************************************************/
    Math::real Ed(real ang) const throw();

    /**
     * The incomplete integral of the second kind in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn
     * @param[in] cn
     * @param[in] dn
     * @return &int; dn(\e w)<sup>2</sup> \e dw (A+S 17.2.10).
     *
     * Instead of specifying the amplitude &phi;, we provide \e sn = sin&phi;,
     * \e cn = cos&phi;, \e dn = sqrt(1 &minus; \e m sin<sup>2</sup>&phi;).
     **********************************************************************/
    Math::real E(real sn, real cn, real dn) const throw();
    ///@}

    /** \name Elliptic functions.
     **********************************************************************/
    ///@{
    /**
     * The Jacobi elliptic functions.
     *
     * @param[in] x the argument.
     * @param[out] sn sn(<i>x</i>|<i>m</i>).
     * @param[out] cn cn(<i>x</i>|<i>m</i>).
     * @param[out] dn dn(<i>x</i>|<i>m</i>).
     **********************************************************************/
    void sncndn(real x, real& sn, real& cn, real& dn) const throw();
    ///@}

    /** \name Symmetric elliptic integrals.
     **********************************************************************/
    ///@{
    /**
     * Symmetric integral of the first kind <i>R<sub>F</sub></i>.
     *
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @return <i>R<sub>F</sub></i>(\e x, \e y, \e z)
     *
     * <i>R<sub>F</sub></i> is defined in http://dlmf.nist.gov/19.16.E1
     * \f[ R_F(x, y, z) = \frac12
     *       \int_0^\infty[(t + x) (t + y) (t + z)]^{-1/2}\, dt \f]
     * If one of the arguments is zero, it is more efficient to call the
     * two-argument version of this function with the non-zero arguments.
     **********************************************************************/
    static real RF(real x, real y, real z) throw();

    /**
     * Complete symmetric integral of the first kind, <i>R<sub>F</sub></i> with
     * one argument zero.
     *
     * @param[in] x
     * @param[in] y
     * @return <i>R<sub>F</sub></i>(\e x, \e y, 0)
     **********************************************************************/
    static real RF(real x, real y) throw();

    /**
     * Degenerate symmetric integral of the first kind <i>R<sub>C</sub></i>.
     *
     * @param[in] x
     * @param[in] y
     * @return <i>R<sub>C</sub></i>(\e x, \e y) = <i>R<sub>F</sub></i>(\e x, \e
     *   y, \e y)
     *
     * <i>R<sub>C</sub></i> is defined in http://dlmf.nist.gov/19.2.E17
     * \f[ R_C(x, y) = \frac12
     *       \int_0^\infty(t + x)^{-1/2}(t + y)^{-1}\,dt \f]
     **********************************************************************/
    static real RC(real x, real y) throw();

    /**
     * Symmetric integral of the second kind <i>R<sub>G</sub></i>.
     *
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @return <i>R<sub>G</sub></i>(\e x, \e y, \e z)
     *
     * <i>R<sub>G</sub></i> is defined in Carlson, eq 1.5
     * \f[ R_G(x, y, z) = \frac14
     *       \int_0^\infty[(t + x) (t + y) (t + z)]^{-1/2}
     *        \biggl(
     *             \frac x{t + x} + \frac y{t + y} + \frac z{t + z}
     *        \biggr)t\,dt \f]
     * See also http://dlmf.nist.gov/19.16.E3.
     * If one of the arguments is zero, it is more efficient to call the
     * two-argument version of this function with the non-zero arguments.
     **********************************************************************/
    static real RG(real x, real y, real z) throw();

    /**
     * Complete symmetric integral of the second kind, <i>R<sub>G</sub></i>
     * with one argument zero.
     *
     * @param[in] x
     * @param[in] y
     * @return  <i>R<sub>G</sub></i>(\e x, \e y, 0)
     **********************************************************************/
    static real RG(real x, real y) throw();

    /**
     * Symmetric integral of the third kind <i>R<sub>J</sub></i>.
     *
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @param[in] p
     * @return <i>R<sub>J</sub></i>(\e x, \e y, \e z, \e p)
     *
     * <i>R<sub>J</sub></i> is defined in http://dlmf.nist.gov/19.16.E2
     * \f[ R_J(x, y, z, p) = \frac32
     *       \int_0^\infty[(t + x) (t + y) (t + z)]^{-1/2} (t + p)^{-1}\, dt \f]
     **********************************************************************/
    static real RJ(real x, real y, real z, real p) throw();

    /**
     * Degenerate symmetric integral of the third kind <i>R<sub>D</sub></i>.
     *
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @return <i>R<sub>D</sub></i>(\e x, \e y, \e z) = <i>R<sub>J</sub></i>(\e
     *   x, \e y, \e z, \e z)
     *
     * <i>R<sub>D</sub></i> is defined in http://dlmf.nist.gov/19.16.E5
     * \f[ R_D(x, y, z) = \frac32
     *       \int_0^\infty[(t + x) (t + y)]^{-1/2} (t + z)^{-3/2}\, dt \f]
     **********************************************************************/
    static real RD(real x, real y, real z) throw();
    ///@}

  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ELLIPTICFUNCTION_HPP

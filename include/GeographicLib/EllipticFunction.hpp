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
   * The notation follows http://dlmf.nist.gov/19 and http://dlmf.nist.gov/22
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
    real _k2, _kp2, _alpha2;
    mutable bool _init;
    mutable real _Kc, _Ec, _Dc, _Pic, _Gc;
    bool Init() const throw();
  public:
    /** \name Constructor
     **********************************************************************/
    ///@{
    /**
     * Constructor.
     *
     * @param[in] k2 the square of the modulus <i>k</i><sup>2</sup>.
     *   <i>k</i><sup>2</sup> must lie in (-&infin;, 1).  (No checking is
     *   done.)
     **********************************************************************/
    explicit EllipticFunction(real k2 = 0) throw();

    /**
     * Constructor specifying the parameter
     *
     * @param[in] k2 the square of the modulus <i>k</i><sup>2</sup>.
     *   <i>k</i><sup>2</sup> must lie in (-&infin;, 1).  (No checking is
     *   done.)
     * @param[in] alpha2 the parameter &alpha;<sup>2</sup>.
     *   &alpha;<sup>2</sup> must lie in (-&infin;, 1).  (No checking is done.)
     **********************************************************************/
    explicit EllipticFunction(real k2, real alpha2) throw();

    /**
     * Reset the modulus and parameter
     *
     * @param[in] k2 the new value of square of the modulus
     *   <i>k</i><sup>2</sup> which must lie in (-&infin;, 1).  (No checking is
     *   done.)
     * @param[in] alpha2 the new value of parameter &alpha;<sup>2</sup>.
     *   &alpha;<sup>2</sup> must lie in (-&infin;, 1).  (No checking is done.)
     **********************************************************************/
    void Reset(real k2, real alpha2 = 0) throw();
    ///@}

    /** \name Inspector functions.
     **********************************************************************/
    ///@{
    /**
     * @return the square of the modulus <i>k</i><sup>2</sup>.
     **********************************************************************/
    Math::real k2() const throw() { return _k2; }

    /**
     * @return the square of the complementary modulus <i>k'</i><sup>2</sup> =
     *   1 &minus; <i>k</i><sup>2</sup>.
     **********************************************************************/
    Math::real kp2() const throw() { return _kp2; }

    /**
     * @return the square of the modulus <i>k</i><sup>2</sup>.
     *
     * <b>DEPRECATED</b>, use k2() instead.
     **********************************************************************/
    Math::real mx() const throw() { return _k2; }

    /**
     * @return the square of the complementary modulus <i>k'</i><sup>2</sup> =
     *   1 &minus; <i>k</i><sup>2</sup>.
     *
     * <b>DEPRECATED</b>, use kp2() instead.
     **********************************************************************/
    Math::real m1x() const throw() { return _kp2; }
    ///@}

    /** \name Complete elliptic integrals.
     **********************************************************************/
    ///@{
    /**
     * The complete integral of the first kind.
     *
     * @return \e K(\e k).
     *
     * \e K(\e k) is defined in http://dlmf.nist.gov/19.2.E4
     * \f[
     *   K(k) = \int_0^{\pi/2} \frac1{\sqrt{1-k^2\sin^2\phi}}\,d\phi.
     * \f]
     **********************************************************************/
    Math::real K() const throw() { _init || Init(); return _Kc; }

    /**
     * The complete integral of the second kind.
     *
     * @return \e E(\e k)
     *
     * \e E(\e k) is defined in http://dlmf.nist.gov/19.2.E5
     * \f[
     *   E(k) = \int_0^{\pi/2} \sqrt{1-k^2\sin^2\phi}\,d\phi.
     * \f]
     **********************************************************************/
    Math::real E() const throw() { _init || Init(); return _Ec; }

    /**
     * Jahnke's complete integral.
     *
     * @return \e D(\e k).
     *
     * \e D(\e k) is defined in http://dlmf.nist.gov/19.2.E6
     * \f[
     *   D(k) = \int_0^{\pi/2} \frac{\sin^2\phi}{\sqrt{1-k^2\sin^2\phi}}\,d\phi.
     * \f]
     **********************************************************************/
    Math::real D() const throw() { _init || Init(); return _Dc; }

    /**
     * The difference between complete integrals of the first and second kinds.
     *
     * @return \e K(\e k) - \e E(\e k).
     **********************************************************************/
    Math::real KE() const throw() { _init || Init(); return _k2 * _Dc; }

    /**
     * The complete integral of the third kind.
     *
     * @return &Pi(&alpha;<sup>2</sup>, \e k)
     *
     * &Pi;(&alpha;<sup>2</sup>, \e k) is defined in
     * http://dlmf.nist.gov/19.2.E7
     * \f[
     *   \Pi(\alpha^2, k) = \int_0^{\pi/2}
     *     \frac1{\sqrt{1-k^2\sin^2\phi}(1 - \alpha^2\sin^2\phi_)}\,d\phi.
     * \f]
     **********************************************************************/
    Math::real Pi() const throw() { _init || Init(); return _Pic; }

    /**
     * The complete geodesic longitude integral.
     *
     * @return \e G(&alpha;<sup>2</sup>, \e k)
     *
     * \e G(&alpha;<sup>2</sup>, \e k) is defined in
     * http://dlmf.nist.gov/19.2.E7
     * \f[
     *   G(\alpha^2, k) = \int_0^{\pi/2}
     *     \frac{\sqrt{1-k^2\sin^2\phi}}{1 - \alpha^2\sin^2\phi}\,d\phi.
     * \f]
     **********************************************************************/
    Math::real G() const throw() { _init || Init(); return _Gc; }
    ///@}

    /** \name Incomplete elliptic integrals.
     **********************************************************************/
    ///@{
    /**
     * The incomplete integral of the first kind.
     *
     * @param[in] phi
     * @return \e F(&phi;, \e k).
     *
     * \e F(&phi;, \e k) is defined in http://dlmf.nist.gov/19.2.E4
     * \f[
     *   F(\phi, k) = \int_0^\phi \frac1{\sqrt{1-k^2\sin^2\theta}}\,d\theta.
     * \f]
     **********************************************************************/
    Math::real F(real phi) const throw();

    /**
     * The incomplete integral of the first kind in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn = sin(&phi;)
     * @param[in] cn = cos(&phi;)
     * @param[in] dn = sqrt(1 &minus; <i>k</i><sup>2</sup>
     *   sin<sup>2</sup>&phi;)
     * @param[in] n the period number such that &phi; = 2<i>n</i>&pi; +
     *   atan2(\e sn, \e cn)
     * @return \e F(&phi;, \e k).
     **********************************************************************/
    Math::real F(real sn, real cn, real dn, int n) const throw();

    /**
     * The incomplete integral of the second kind.
     *
     * @param[in] phi
     * @return \e E(&phi;, \e k).
     *
     * \e E(&phi;, \e k) is defined in http://dlmf.nist.gov/19.2.E5
     * \f[
     *   E(\phi, k) = \int_0^\phi \sqrt{1-k^2\sin^2\theta}\,d\theta.
     * \f]
     **********************************************************************/
    Math::real E(real phi) const throw();

    /**
     * The incomplete integral of the second kind with the argument given in
     * degrees.
     *
     * @param[in] ang in <i>degrees</i>.
     * @return \e E(&pi; <i>ang</i>/180, \e k).
     **********************************************************************/
    Math::real Ed(real ang) const throw();

    /**
     * The inverse of the incomplete integral of the second kind.
     *
     * @param[in] x
     * @return &phi; such that \e E(&phi;, \e k) = \e x.
     **********************************************************************/
    Math::real Einv(real x) const throw();

    /**
     * The incomplete integral of the second kind in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn = sin(&phi;)
     * @param[in] cn = cos(&phi;)
     * @param[in] dn = sqrt(1 &minus; <i>k</i><sup>2</sup>
     *   sin<sup>2</sup>&phi;)
     * @param[in] n the period number such that &phi; = 2<i>n</i>&pi; +
     *   atan2(\e sn, \e cn)
     * @return \e E(&phi;, \e k).
     **********************************************************************/
    Math::real E(real sn, real cn, real dn, int n) const throw();

    /**
     * The incomplete integral of the second kind in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn = sin(&phi;)
     * @param[in] cn = cos(&phi;)
     * @param[in] dn = sqrt(1 &minus; <i>k</i><sup>2</sup>
     *   sin<sup>2</sup>&phi;)
     * @return \e E(&phi;, \e k).
     *
     * <b>COMPATIBILITY NOTE:</b> This function calls E(real, real, real, int)
     * with a 4th argument of 0.  At some point, E(real, real, real) and will
     * be withdrawn and the interface to E(real, real, real, int) changed so
     * that its 4th defaults to 0.  This will preserve source-level
     * compatibility.
     **********************************************************************/
    Math::real E(real sn, real cn, real dn) const throw();

    /**
     * The incomplete integral of the third kind.
     *
     * @param[in] phi
     * @return &Pi;(&phi;, &alpha;<sup>2</sup>, \e k).
     *
     * \e &Pi;(&phi;, &alpha;<sup>2</sup>, \e k) is defined in
     * http://dlmf.nist.gov/19.2.E7
     * \f[
     *   \Pi(\phi, \alpha^2, k) = \int_0^\phi
     *     \frac1{\sqrt{1-k^2\sin^2\theta}(1 - \alpha^2\sin^2\theta_)}\,d\theta.
     * \f]
     **********************************************************************/
    Math::real Pi(real phi) const throw();

    /**
     * The incomplete integral of the third kind in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn = sin(&phi;)
     * @param[in] cn = cos(&phi;)
     * @param[in] dn = sqrt(1 &minus; <i>k</i><sup>2</sup>
     *   sin<sup>2</sup>&phi;)
     * @param[in] n the period number such that &phi; = 2<i>n</i>&pi; +
     *   atan2(\e sn, \e cn)
     * @return &Pi;(&phi;, &alpha;<sup>2</sup>, \e k).
     **********************************************************************/
    Math::real Pi(real sn, real cn, real dn, int n) const throw();

    /**
     * Jahnke's incomplete elliptic integral.
     *
     * @param[in] phi
     * @return \e D(&phi;, \e k).
     *
     * \e D(&phi;, \e k) is defined in http://dlmf.nist.gov/19.2.E4
     * \f[
     *   D(\phi, k) = \int_0^\phi
     *    \frac{\sin^2\theta}{\sqrt{1-k^2\sin^2\theta}}\,d\theta.
     * \f]
     **********************************************************************/
    Math::real D(real phi) const throw();

    /**
     * Jahnke's incomplete elliptic integral in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn = sin(&phi;)
     * @param[in] cn = cos(&phi;)
     * @param[in] dn = sqrt(1 &minus; <i>k</i><sup>2</sup>
     *   sin<sup>2</sup>&phi;)
     * @param[in] n the period number such that &phi; = 2<i>n</i>&pi; +
     *   atan2(\e sn, \e cn)
     * @return \e D(&phi;, \e k).
     **********************************************************************/
    Math::real D(real sn, real cn, real dn, int n) const throw();

    /**
     * The geodesic longitude integral.
     *
     * @param[in] phi
     * @return \e G(&phi;, &alpha;<sup>2</sup>, \e k).
     *
     * \e G(&phi;, &alpha;<sup>2</sup>, \e k) is defined by
     * \f[
     *   G(\phi, \alpha^2, k) =
     *   \biggl(1 - \frac{k^2}{\alpha^2}\biggr) \Pi(\phi, \alpha^2, k)
     *   + \frac{k^2}{\alpha^2} F(\phi, k) = \int_0^\phi
     *     \frac{\sqrt{1-k^2\sin^2\theta}}{1 - \alpha^2\sin^2\theta}\,d\theta.
     * \f]
     **********************************************************************/
    Math::real G(real phi) const throw();

    /**
     * The geodesic longitude integral in terms of Jacobi elliptic
     * functions
     *
     * @param[in] sn = sin(&phi;)
     * @param[in] cn = cos(&phi;)
     * @param[in] dn = sqrt(1 &minus; <i>k</i><sup>2</sup>
     *   sin<sup>2</sup>&phi;)
     * @param[in] n the period number such that &phi; = 2<i>n</i>&pi; +
     *   atan2(\e sn, \e cn)
     * @return \e G(&phi;, &alpha;<sup>2</sup>, \e k).
     **********************************************************************/
    Math::real G(real sn, real cn, real dn, int n) const throw();

    ///@}

    /** \name Elliptic functions.
     **********************************************************************/
    ///@{
    /**
     * The Jacobi elliptic functions.
     *
     * @param[in] x the argument.
     * @param[out] sn sn(\e x, \e k).
     * @param[out] cn cn(\e x, \e k).
     * @param[out] dn dn(\e x, \e k).
     **********************************************************************/
    void sncndn(real x, real& sn, real& cn, real& dn) const throw();

    /**
     * The &Delta; function.
     *
     * @param[in] sn sin&phi;
     * @param[in] cn cos&phi;
     * @return &Delta; = sqrt(1 &minus; <i>k</i><sup>2</sup>
     *   sin<sup>2</sup>&phi;)
     **********************************************************************/
    Math::real Delta(real sn, real cn) const throw()
    { return sqrt(_k2 < 0 ? 1 - _k2 * sn*sn : _kp2 + _k2 * cn*cn); }
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
     *       \int_0^\infty\frac1{\sqrt{(t + x) (t + y) (t + z)}}\, dt \f]
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
     *       \int_0^\infty\frac1{\sqrt{t + x}(t + y)}\,dt \f]
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
     * @return <i>R<sub>G</sub></i>(\e x, \e y, 0)
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

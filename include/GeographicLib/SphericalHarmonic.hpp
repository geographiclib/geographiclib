/**
 * \file SphericalHarmonic.hpp
 * \brief Header for GeographicLib::SphericalHarmonic class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_SPHERICALHARMONIC_HPP)
#define GEOGRAPHICLIB_SPHERICALHARMONIC_HPP "$Id$"

#include <vector>
#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief Spherical Harmonic series
   **********************************************************************/

  class GEOGRAPHIC_EXPORT SphericalHarmonic {
  public:
    typedef Math::real work;
    typedef Math::real T;
  private:
    // An internal scaling of the coefficients to avoid overflow in
    // intermediate calculations.
    static const T scale_;
    static const T eps_;
  public:
    /**
     * Compute a spherical harmonic sum.
     *
     * @param[in] N the maximum order and degree of the sum.
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[in] a the scaling radius for the coordinates.
     * @return \e V the spherical harmonic sum.
     *
     * Evaluate the spherical harmonic sum \verbatim
 V(x, y, z) = sum(n=0..N)[ q^(n+1) * sum(m=0..n)[
   (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) *
   Pbar[n,m](cos(theta)) ] ]
\endverbatim
     * where
     * - <i>p</i><sup>2</sup> = <i>x</i><sup>2</sup> + <i>y</i><sup>2</sup>.
     * - <i>r</i><sup>2</sup> = <i>p</i><sup>2</sup> + <i>z</i><sup>2</sup>.
     * - \e q = <i>a</i>/<i>r</i>.
     * - \e theta = atan2(\e p, \e z).
     * - \e lambda = atan2(\e y, \e x).
     * - Pbar<sub>\e nm</sub>(\e t) is the fully normalized associate
     *   Legendre function of degree \e n and order \e m.
     *
     * The coefficients \e C<sub>\e nm</sub> and \e S<sub>\e nm</sub>
     * are stored in the one-dimensional vectors \e C and \e S which must
     * contain (at least) (\e N + 1)(\e N + 2)/2 elements, stored in
     * "column-major" order.  Thus for \e N = 3, the order would be:
     * <i>C</i><sub>00</sub>,
     * <i>C</i><sub>10</sub>,
     * <i>C</i><sub>20</sub>,
     * <i>C</i><sub>30</sub>,
     * <i>C</i><sub>11</sub>,
     * <i>C</i><sub>21</sub>,
     * <i>C</i><sub>31</sub>,
     * <i>C</i><sub>22</sub>,
     * <i>C</i><sub>32</sub>,
     * <i>C</i><sub>33</sub>.
     * The first (\e N + 1) elements of \e S should be 0.
     **********************************************************************/
    static Math::real Value(int N,
                            const std::vector<double>& C,
                            const std::vector<double>& S,
                            work x, work y, work z, work a);
    /**
     * Compute a spherical harmonic sum and its gradient.
     *
     * @param[in] N the maximum order and degree of the sum.
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[in] a the scaling radius for the coordinates.
     * @param[out] gradx \e x component of the gradient
     * @param[out] grady \e y component of the gradient
     * @param[out] gradz \e z component of the gradient
     * @return \e V the spherical harmonic sum.
     *
     **********************************************************************/
    static Math::real Value(int N,
                            const std::vector<double>& C,
                            const std::vector<double>& S,
                            work x, work y, work z, work a,
                            work& gradx, work& grady, work& gradz);
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALHARMONIC_HPP

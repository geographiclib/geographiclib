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
   *
   * Sum a spherical harmonic series.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT SphericalHarmonic {
  private:
    typedef Math::real real;
    // Max number of contributions to coefficients (accommodates potential +
    // its 1st, 2nd, 3rd time derivatives in magnetic field models).
    static const int lmax = 4;
    // An internal scaling of the coefficients to avoid overflow in
    // intermediate calculations.
    static const real scale_;
    // Move latitudes near the pole off the axis by this amount.
    static const real eps_;

    // The 1-d index of column major vector for max degree N, degree n, and
    // order m
    static inline int index(int N, int n, int m) throw()
    { return m * N - m * (m - 1) / 2 + n; }

    SphericalHarmonic();        // Disable constructor
  public:
    enum normalization {
      full = 0,
      schmidt = 1,
    };

    /**
     * Compute a spherical harmonic sum.
     *
     * @param[in] N the maximum order and degree of the sum.
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] Cp vector of correction coefficients for cosine terms.
     * @param[in] Sp vector of correction coefficients for sine terms.
     * @param[in] tau multiplier for correction coefficients.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[in] a the scaling radius for the coordinates.
     * @return \e V the spherical harmonic sum.
     *
     * Evaluate the spherical harmonic sum \verbatim
 V(x, y, z) = sum(n = 0..N)[ q^(n+1) * sum(m = 0..n)[
   (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) *
   Pbar[n,m](cos(theta)) ] ]
\endverbatim
     * where
     * - <i>p</i><sup>2</sup> = <i>x</i><sup>2</sup> + <i>y</i><sup>2</sup>,
     * - <i>r</i><sup>2</sup> = <i>p</i><sup>2</sup> + <i>z</i><sup>2</sup>,
     * - \e q = <i>a</i>/<i>r</i>,
     * - \e theta = atan2(\e p, \e z) = the spherical \e colatitude,
     * - \e lambda = atan2(\e y, \e x) = the longitude.
     * - Pbar<sub>\e nm</sub>(\e t) is the fully normalized associated
     *   Legendre function of degree \e n and order \e m.
     *
     * The coefficients \e C<sub>\e nm</sub> and \e S<sub>\e nm</sub> are
     * stored in the one-dimensional vectors \e C and \e S which must contain
     * (\e N + 1)(\e N + 2)/2 elements, stored in "column-major" order.  Thus
     * for \e N = 3, the order would be:
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
     * In general the (\e n,\e m) element is at index \e m*\e N - \e m*(\e m -
     * 1)/2 + \e n.  The first (\e N + 1) elements of \e S should be 0.
     *
     * If \e tau is non-zero, then \e tau \e Cp and \e tau \e Sp are added to
     * the coefficients \e C<sub>\e nm</sub> and \e S<sub>\e nm</sub> in the
     * definition of \e V.  \e Cp and \e Sp are stored in the same order as \e
     * C and \e S, except that the vectors may be shorter (missing elements are
     * assumed to be zero).  Typical usage of \e Cp and \e Sp:
     * - They are both empty, and \e V is computed with \e C and \e S.
     * - The first few even numbered elements of \e Cp (corresponding to \e n
     *   even and \e m = 0) are defined, \e Sp is empty, and \e tau = -1.  This
     *   allows the "normal" potential (the potential of the ellipsoid) to be
     *   subtracted from \e V.
     * - \e Cp and \e Sp represent the secular variation of the potential and
     *   \e tau represents the time.  This allows a simple time varying field
     *   to be modeled (e.g., the world magnetic model).
     *
     * References:
     * - C. W. Clenshaw, A note on the summation of Chebyshev series,
     *   %Math. Tables Aids Comput. 9(51), 118-120 (1955).
     * - R. E. Deakin, Derivatives of the earth's potentials, Geomatics
     *   Research Australasia 68, 31-60, (June 1998).
     * - W. A. Heiskanen and H. Moritz, Physical Geodesy, (Freeman, San
     *   Fransisco, 1967).  (See Sec. 1-14, for a definition of Pbar.)
     * - S. A. Holmes and W. E. Featherstone, A unified approach to the
     *   Clenshaw summation and the recursive computation of very high degree
     *   and order normalised associated Legendre functions, J. Geod. 76(5),
     *   279-299 (2002).
     * - C. C. Tscherning and K. Poder, Some geodetic applications of Clenshaw
     *   summation, Boll. Geod. Sci. Aff. 41(4), 349-375 (1982).
     **********************************************************************/
    static Math::real Value(int N,
                            const std::vector<double>& C,
                            const std::vector<double>& S,
                            const std::vector<double>& Cp,
                            const std::vector<double>& Sp,
                            real tau, real x, real y, real z, real a);
    /**
     * Compute a spherical harmonic sum and its gradient.
     *
     * @param[in] N the maximum order and degree of the sum.
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] Cp vector of correction coefficients for cosine terms.
     * @param[in] Sp vector of correction coefficients for sine terms.
     * @param[in] tau multiplier for correction coefficients.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[in] a the scaling radius for the coordinates.
     * @param[out] gradx \e x component of the gradient
     * @param[out] grady \e y component of the gradient
     * @param[out] gradz \e z component of the gradient
     * @return \e V the spherical harmonic sum.
     *
     * This is the same as the previous function, except that the components of
     * the gradients of the sum in the \e x, \e y, and \e z directions are
     * computed.
     **********************************************************************/
    static Math::real Value(int N,
                            const std::vector<double>& C,
                            const std::vector<double>& S,
                            const std::vector<double>& Cp,
                            const std::vector<double>& Sp,
                            real tau, real x, real y, real z, real a,
                            real& gradx, real& grady, real& gradz);
    template<bool gradp, normalization norm>
      static Math::real TValue(int N,
                               const std::vector<double>& C,
                               const std::vector<double>& S,
                               const std::vector<double>& Cp,
                               const std::vector<double>& Sp,
                               real tau, real x, real y, real z, real a,
                               real& gradx, real& grady, real& gradz);
  };


} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALHARMONIC_HPP

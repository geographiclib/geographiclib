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
#include <GeographicLib/SphericalEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Spherical Harmonic series
   *
   * Sum a spherical harmonic series.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT SphericalHarmonic {
  public:
    enum normalization {
      full = SphericalEngine::full,
      schmidt = SphericalEngine::schmidt,
    };

  private:
    typedef Math::real real;
    SphericalEngine::coeff _c[1];
    real _a;
    normalization _norm;

  public:
    SphericalHarmonic() {}
    SphericalHarmonic(const std::vector<double>& C,
                      const std::vector<double>& S,
                      int N, real a, normalization norm = full)
      : _a(a)
      , _norm(norm)
    { _c[0] = SphericalEngine::coeff(C, S, N); }

    SphericalHarmonic(const std::vector<double>& C,
                      const std::vector<double>& S,
                      int N, int nmx, int mmx,
                      real a, normalization norm = full)
      : _a(a)
      , _norm(norm)
    { _c[0] = SphericalEngine::coeff(C, S, N, nmx, mmx); }

    /**
     * Compute a spherical harmonic sum.
     *
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
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
    Math::real operator()(real x, real y, real z)  const {
      real f[] = {1};
      real v = 0;
      real dummy;
      switch (_norm) {
      case full:
        v = SphericalEngine::GenValue<false, SphericalEngine::full, 1>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      case schmidt:
        v = SphericalEngine::GenValue<false, SphericalEngine::schmidt, 1>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      }
      return v;
    }
    /**
     * Compute a spherical harmonic sum and its gradient.
     *
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[out] gradx \e x component of the gradient
     * @param[out] grady \e y component of the gradient
     * @param[out] gradz \e z component of the gradient
     * @return \e V the spherical harmonic sum.
     *
     * This is the same as the previous function, except that the components of
     * the gradients of the sum in the \e x, \e y, and \e z directions are
     * computed.
     **********************************************************************/
    Math::real operator()(real x, real y, real z,
                          real& gradx, real& grady, real& gradz) const {
      real f[] = {1};
      real v = 0;
      switch (_norm) {
      case full:
        v = SphericalEngine::GenValue<true, SphericalEngine::full, 1>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      case schmidt:
        v = SphericalEngine::GenValue<true, SphericalEngine::schmidt, 1>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      }
      return v;
    }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALHARMONIC_HPP

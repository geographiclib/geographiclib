/**
 * \file SphericalEngine.hpp
 * \brief Header for GeographicLib::SphericalEngine class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_SPHERICALENGINE_HPP)
#define GEOGRAPHICLIB_SPHERICALENGINE_HPP "$Id$"

#include <vector>
#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief Spherical Harmonic series
   *
   * Sum a spherical harmonic series.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT SphericalEngine {
  private:
    typedef Math::real real;
    // An internal scaling of the coefficients to avoid overflow in
    // intermediate calculations.
    static const real scale_;
    // Move latitudes near the pole off the axis by this amount.
    static const real eps_;
    static const std::vector<real> Z_;
    SphericalEngine();        // Disable constructor

  public:
    enum normalization {
      full = 0,
      schmidt = 1,
    };
    class coeff {
    public:
      std::vector<real>::const_iterator Cnm;
      std::vector<real>::const_iterator Snm;
      int N, nmx, mmx;
      // The 1-d index of column major vector for max degree N, degree n, and
      // order m
      inline int index(int n, int m) const throw()
      { return m * N - m * (m - 1) / 2 + n; }
      // Index  of element after m'th column
      inline int rowind(int N1, int m) const throw() {
        // Normally use nmx, however, it may be used in a loop where the
        // starting degree is N1
        return index(std::min(nmx, N1) + 1, m);
      }
      coeff()
        : Cnm(Z_.begin())
        , Snm(Z_.begin())
        , N(-1)
        , nmx(-1)
        , mmx(-1) {}
      coeff(const std::vector<real>& C,
            const std::vector<real>& S,
            int N1, int nmx1, int mmx1)
        : Cnm(C.begin())
        , Snm(S.begin())
        , N(N1)
        , nmx(nmx1)
        , mmx(mmx1) {
        if (!(N >= nmx && nmx >= mmx && mmx >= 0))
          throw GeographicErr("Bad indices for coeff");
        if (!(index(nmx, mmx) < int(C.size()) &&
              index(nmx, mmx) < int(S.size())))
          throw GeographicErr("Arrays too small in coeff");
      }
      coeff(const std::vector<real>& C,
            const std::vector<real>& S,
            int N1)
        : Cnm(C.begin())
        , Snm(S.begin())
        , N(N1)
        , nmx(N1)
        , mmx(N1) {
        if (!(N >= nmx && nmx >= mmx && mmx >= 0))
          throw GeographicErr("Bad indices for coeff");
        if (!(index(nmx, mmx) < int(C.size()) &&
              index(nmx, mmx) < int(S.size())))
          throw GeographicErr("Arrays too small in coeff");
      }
      coeff(const std::vector<real>& C,
            const std::vector<real>& S,
            int N1, int mmx1)
        : Cnm(C.begin())
        , Snm(S.begin())
        , N(N1)
        , nmx(N1)
        , mmx(mmx1) {
        if (!(N >= nmx && nmx >= mmx && mmx >= 0))
          throw GeographicErr("Bad indices for coeff");
        if (!(index(nmx, mmx) < int(C.size()) &&
              index(nmx, mmx) < int(S.size())))
          throw GeographicErr("Arrays too small in coeff");
      }
    };

    template<bool gradp, normalization norm, int L>
      static Math::real GenValue(const coeff c[L], const real f[L],
                                 real x, real y, real z, real a,
                                 real& gradx, real& grady, real& gradz);

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
     * @param[in] norm the normalization of the Legendre functions.
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
    static Math::real Value(int N,
                            const std::vector<double>& C,
                            const std::vector<double>& S,
                            real x, real y, real z, real a,
                            normalization norm = full) {
      coeff c[] = {coeff(C, S, N)};
      real f[] = {1};
      real v = 0;
      real dummy;
      switch (norm) {
      case full:
        v = GenValue<false, full, 1>(c, f, x, y, z, a, dummy, dummy, dummy);
        break;
      case schmidt:
        v = GenValue<false, schmidt, 1>(c, f, x, y, z, a, dummy, dummy, dummy);
        break;
      }
      return v;
    }
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
     * @param[in] norm the normalization of the Legendre functions.
     * @return \e V the spherical harmonic sum.
     *
     * This is the same as the previous function, except that the components of
     * the gradients of the sum in the \e x, \e y, and \e z directions are
     * computed.
     **********************************************************************/
    static Math::real Gradient(int N,
                               const std::vector<double>& C,
                               const std::vector<double>& S,
                               real x, real y, real z, real a,
                               real& gradx, real& grady, real& gradz,
                               normalization norm = full) {
      coeff c[] = {coeff(C, S, N)};
      real f[] = {1};
      real v = 0;
      switch (norm) {
      case full:
        v = GenValue<true, full, 1>(c, f, x, y, z, a, gradx, grady, gradz);
        break;
      case schmidt:
        v = GenValue<true, schmidt, 1>(c, f, x, y, z, a, gradx, grady, gradz);
        break;
      }
      return v;
    }
    /**
     * Compute a spherical harmonic sum with a correction term.
     *
     * @param[in] N the maximum order and degree of the sum.
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] Np vector of correction coefficients for cosine terms.
     * @param[in] Cp vector of correction coefficients for cosine terms.
     * @param[in] Sp vector of correction coefficients for sine terms.
     * @param[in] tau multiplier for correction coefficients.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[in] a the scaling radius for the coordinates.
     * @param[in] norm the normalization of the Legendre functions.
     * @return \e V the spherical harmonic sum.
     *
     * If \e tau is non-zero, then \e tau \e Cp and \e tau \e Sp are added to
     * the coefficients \e C<sub>\e nm</sub> and \e S<sub>\e nm</sub> in the
     * definition of \e V if the degree and order are less than or equal to \e
     * Np.  \e Cp and \e Sp are stored in the same order as \e C and \e S,
     * except that the vectors are shorter if \e Np < \e N.  Typical usage of
     * \e Cp and \e Sp:
     * - They are both empty, and \e V is computed with \e C and \e S.
     * - The first few even numbered elements of \e Cp (corresponding to \e n
     *   even and \e m = 0) are defined, \e Sp is empty, and \e tau = -1.  This
     *   allows the "normal" potential (the potential of the ellipsoid) to be
     *   subtracted from \e V.
     * - \e Cp and \e Sp represent the secular variation of the potential and
     *   \e tau represents the time.  This allows a simple time varying field
     *   to be modeled (e.g., the world magnetic model).
     **********************************************************************/
    static Math::real Value1(int N,
                             const std::vector<double>& C,
                             const std::vector<double>& S,
                             int Np, real tau,
                             const std::vector<double>& Cp,
                             const std::vector<double>& Sp,
                             real x, real y, real z, real a,
                             normalization norm = full) {
      coeff c[] = {coeff(C, S, N), coeff(Cp, Sp, Np)};
      real f[] = {1, tau};
      real v = 0;
      real dummy;
      switch (norm) {
      case full:
        v =  GenValue<false, full, 2>(c, f, x, y, z, a, dummy, dummy, dummy);
        break;
      case schmidt:
        v = GenValue<false, schmidt, 2>(c, f, x, y, z, a, dummy, dummy, dummy);
        break;
      }
      return v;
    }
    /**
     * Compute a spherical harmonic sum with a correction and its gradient.
     *
     * @param[in] N the maximum order and degree of the sum.
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] Np vector of correction coefficients for cosine terms.
     * @param[in] tau multiplier for correction coefficients.
     * @param[in] Cp vector of correction coefficients for cosine terms.
     * @param[in] Sp vector of correction coefficients for sine terms.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[in] a the scaling radius for the coordinates.
     * @param[out] gradx \e x component of the gradient
     * @param[out] grady \e y component of the gradient
     * @param[out] gradz \e z component of the gradient
     * @param[in] norm the normalization of the Legendre functions.
     * @return \e V the spherical harmonic sum.
     *
     * This is the same as the previous function, except that the components of
     * the gradients of the sum in the \e x, \e y, and \e z directions are
     * computed.
     **********************************************************************/
    static Math::real Gradient1(int N,
                                const std::vector<double>& C,
                                const std::vector<double>& S,
                                int Np, real tau,
                                const std::vector<double>& Cp,
                                const std::vector<double>& Sp,
                                real x, real y, real z, real a,
                                real& gradx, real& grady, real& gradz,
                                normalization norm = full) {
      coeff c[] = {coeff(C, S, N), coeff(Cp, Sp, Np)};
      real f[] = {1, tau};
      real v = 0;
      switch (norm) {
      case full:
        v = GenValue<true, full, 2>(c, f, x, y, z, a, gradx, grady, gradz);
        break;
      case schmidt:
        v = GenValue<true, schmidt, 2>(c, f, x, y, z, a, gradx, grady, gradz);
        break;
      }
      return v;
    }

    /**
     * Compute a spherical harmonic sum with a zonal correction term.
     *
     * @param[in] N the maximum order and degree of the sum.
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] Np vector of correction coefficients for cosine terms.
     * @param[in] tau multiplier for correction coefficients.
     * @param[in] Cp vector of correction coefficients for cosine terms.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[in] a the scaling radius for the coordinates.
     * @param[in] norm the normalization of the Legendre functions.
     * @return \e V the spherical harmonic sum.
     **********************************************************************/
    static Math::real Value1Z(int N,
                              const std::vector<double>& C,
                              const std::vector<double>& S,
                              int Np, real tau,
                              const std::vector<double>& Cp,
                              real x, real y, real z, real a,
                              normalization norm = full) {
      // There's no Sp because it doesn't contribute to the m = 0 term.
      coeff c[] = {coeff(C, S, N), coeff(Cp, Cp, Np, 0)};
      real f[] = {1, tau};
      real v = 0;
      real dummy;
      switch (norm) {
      case full:
        v =  GenValue<false, full, 2>(c, f, x, y, z, a, dummy, dummy, dummy);
        break;
      case schmidt:
        v = GenValue<false, schmidt, 2>(c, f, x, y, z, a, dummy, dummy, dummy);
        break;
      }
      return v;
    }

  };


} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALENGINE_HPP

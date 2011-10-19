/**
 * \file SphericalHarmonic2.hpp
 * \brief Header for GeographicLib::SphericalHarmonic2 class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_SPHERICALHARMONIC2_HPP)
#define GEOGRAPHICLIB_SPHERICALHARMONIC2_HPP "$Id$"

#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/SphericalEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Spherical Harmonic series
   *
   * Sum a spherical harmonic series.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT SphericalHarmonic2 {
  public:
    enum normalization {
      full = SphericalEngine::full,
      schmidt = SphericalEngine::schmidt,
    };

  private:
    typedef Math::real real;
    SphericalEngine::coeff _c[3];
    real _a;
    normalization _norm;

  public:
    SphericalHarmonic2() {}
    SphericalHarmonic2(const std::vector<double>& C,
                       const std::vector<double>& S,
                       int N,
                       const std::vector<double>& C1,
                       const std::vector<double>& S1,
                       int N1,
                       const std::vector<double>& C2,
                       const std::vector<double>& S2,
                       int N2,
                       real a, normalization norm = full)
      : _a(a)
      , _norm(norm) {
      _c[0] = SphericalEngine::coeff(C, S, N);
      _c[1] = SphericalEngine::coeff(C1, S1, N1);
      _c[2] = SphericalEngine::coeff(C2, S2, N2);
    }

    SphericalHarmonic2(const std::vector<double>& C,
                       const std::vector<double>& S,
                       int N, int nmx, int mmx,
                       const std::vector<double>& C1,
                       const std::vector<double>& S1,
                       int N1, int nmx1, int mmx1,
                       const std::vector<double>& C2,
                       const std::vector<double>& S2,
                       int N2, int nmx2, int mmx2,
                       real a, normalization norm = full)
      : _a(a)
      , _norm(norm) {
      _c[0] = SphericalEngine::coeff(C, S, N, nmx, mmx);
      _c[1] = SphericalEngine::coeff(C1, S1, N1, nmx1, mmx1);
      _c[2] = SphericalEngine::coeff(C2, S2, N2, nmx2, mmx2);
    }

    /**
     * Compute a spherical harmonic sum with two correction terms.
     *
     * @param[in] tau1 multiplier for correction coefficients.
     * @param[in] tau2 multiplier for correction coefficients.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
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
    Math::real operator()(real tau1, real tau2, real x, real y, real z) const {
      real f[] = {1, tau1, tau2};
      real v = 0;
      real dummy;
      switch (_norm) {
      case full:
        v = SphericalEngine::Value<false, SphericalEngine::full, 3>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      case schmidt:
        v = SphericalEngine::Value<false, SphericalEngine::schmidt, 3>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      }
      return v;
    }

    /**
     * Compute a spherical harmonic sum with two correction terms and its
     * gradient.
     *
     * @param[in] tau1 multiplier for correction coefficients.
     * @param[in] tau2 multiplier for correction coefficients.
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
    Math::real operator()(real tau1, real tau2, real x, real y, real z,
                          real& gradx, real& grady, real& gradz) const {
      real f[] = {1, tau1, tau2};
      real v = 0;
      switch (_norm) {
      case full:
        v = SphericalEngine::Value<true, SphericalEngine::full, 3>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      case schmidt:
        v = SphericalEngine::Value<true, SphericalEngine::schmidt, 3>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      }
      return v;
    }
    CircularEngine Circle(real tau1, real tau2, real p, real z, bool gradp)
      const {
      real f[] = {1, tau1, tau2};
      switch (_norm) {
      case full:
        return gradp ?
          SphericalEngine::Circle<true, SphericalEngine::full, 3>
          (_c, f, p, z, _a) :
          SphericalEngine::Circle<false, SphericalEngine::full, 3>
          (_c, f, p, z, _a);
        break;
      case schmidt:
      default:                  // To avoid compiler warnings
        return gradp ?
          SphericalEngine::Circle<true, SphericalEngine::schmidt, 3>
          (_c, f, p, z, _a) :
          SphericalEngine::Circle<false, SphericalEngine::schmidt, 3>
          (_c, f, p, z, _a);
        break;
      }
    }  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALHARMONIC2_HPP

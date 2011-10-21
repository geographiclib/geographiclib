/**
 * \file SphericalHarmonic1.hpp
 * \brief Header for GeographicLib::SphericalHarmonic1 class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_SPHERICALHARMONIC1_HPP)
#define GEOGRAPHICLIB_SPHERICALHARMONIC1_HPP "$Id$"

#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/SphericalEngine.hpp>
#include <GeographicLib/CircularEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Spherical Harmonic series with a correction to the coefficients.
   *
   * This classes is similar to SphericalHarmonic, except that the coefficients
   * \e C<sub>\e nm</sub> are replaced by \e C<sub>\e nm</sub> + \e tau
   * C'<sub>\e nm</sub> (and similarly for \e S<sub>\e nm</sub>).
   **********************************************************************/

  class GEOGRAPHIC_EXPORT SphericalHarmonic1 {
  public:
    /**
     * Supported normalizations for associate Legendre polynomials.
     **********************************************************************/
    enum normalization {
      /**
       * Fully normalized associated Legendre polynomials.  See
       * SphericalHarmonic::full for documentation.
       *
       * @hideinitializer
       **********************************************************************/
      full = SphericalEngine::full,
      /**
       * Schmidt semi-normalized associated Legendre polynomials.  See
       * SphericalHarmonic::schmidt for documentation.
       *
       * @hideinitializer
       **********************************************************************/
      schmidt = SphericalEngine::schmidt,
    };

  private:
    typedef Math::real real;
    SphericalEngine::coeff _c[2];
    real _a;
    normalization _norm;

  public:
    /**
     * Constructor with a full set of coefficients specified.
     *
     * @param[in] C the coefficients \e C<sub>\e nm</sub>.
     * @param[in] S the coefficients \e S<sub>\e nm</sub>.
     * @param[in] N the maximum degree and order of the sum
     * @param[in] C1 the coefficients \e C'<sub>\e nm</sub>.
     * @param[in] S1 the coefficients \e S'<sub>\e nm</sub>.
     * @param[in] N1 the maximum degree and order of the correction
     *   coefficients \e C'<sub>\e nm</sub> and \e S'<sub>\e nm</sub>.
     * @param[in] a the reference radius appearing in the definition of the
     *   sum.
     * @param[in] norm the normalization for the associated Legendre
     *   polynomials, either SphericalHarmonic1::full (the default) or
     *   SphericalHarmonic1::schmidt.
     *
     * See SphericalHarmonic for the way the coefficients should be stored.  \e
     * N1 should satisfy \e N1 <= \e N.
     *
     * The class stores <i>pointers</i> to the first elements of \e C, \e S, \e
     * C', and \e S'.  These arrays should not be altered or destroyed during
     * the lifetime of a SphericalHarmonic object.
     **********************************************************************/
    SphericalHarmonic1(const std::vector<double>& C,
                       const std::vector<double>& S,
                       int N,
                       const std::vector<double>& C1,
                       const std::vector<double>& S1,
                       int N1,
                       real a, normalization norm = full)
      : _a(a)
      , _norm(norm) {
      if (!(N1 <= N))
        throw GeographicErr("N1 cannot be larger that N");
      _c[0] = SphericalEngine::coeff(C, S, N);
      _c[1] = SphericalEngine::coeff(C1, S1, N1);
    }

    /**
     * Constructor with a subset of coefficients specified.
     *
     * @param[in] C the coefficients \e C<sub>\e nm</sub>.
     * @param[in] S the coefficients \e S<sub>\e nm</sub>.
     * @param[in] N the degree used to determine the layout of \e C and \e S.
     * @param[in] nmx the maximum degree used in the sum.  The sum over \e n is
     *   from 0 thru \e nmx.
     * @param[in] mmx the maximum order used in the sum.  The sum over \e m is
     *   from 0 thru min(\e n, \e mmx).
     * @param[in] C1 the coefficients \e C'<sub>\e nm</sub>.
     * @param[in] S1 the coefficients \e S'<sub>\e nm</sub>.
     * @param[in] N1 the degree used to determine the layout of \e C' and \e
     *   S'.
     * @param[in] nmx1 the maximum degree used for \e C' and \e S'.
     * @param[in] mmx1 the maximum order used for \e C' and \e S'.
     * @param[in] a the reference radius appearing in the definition of the
     *   sum.
     * @param[in] norm the normalization for the associated Legendre
     *   polynomials, either SphericalHarmonic1::full (the default) or
     *   SphericalHarmonic1::schmidt.
     *
     * The class stores <i>pointers</i> to the first elements of \e C, \e S, \e
     * C', and \e S'.  These arrays should not be altered or destroyed during
     * the lifetime of a SphericalHarmonic object.
     **********************************************************************/
    SphericalHarmonic1(const std::vector<double>& C,
                       const std::vector<double>& S,
                       int N, int nmx, int mmx,
                       const std::vector<double>& C1,
                       const std::vector<double>& S1,
                       int N1, int nmx1, int mmx1,
                       real a, normalization norm = full)
      : _a(a)
      , _norm(norm) {
      if (!(nmx1 <= nmx))
        throw GeographicErr("nmx1 cannot be larger that nmx");
      if (!(mmx1 <= mmx))
        throw GeographicErr("mmx1 cannot be larger that mmx");
      _c[0] = SphericalEngine::coeff(C, S, N, nmx, mmx);
      _c[1] = SphericalEngine::coeff(C1, S1, N1, nmx1, mmx1);
    }

    /**
     * A default constructor so that the object can be created when the
     * constructor for another object is initialized.  This default object can
     * then be reset with the default copy assignment operator.
     **********************************************************************/
    SphericalHarmonic1() {}

    /**
     * Compute a spherical harmonic sum with a correction term.
     *
     * @param[in] tau multiplier for correction coefficients \e C' and \e S'.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @return \e V the spherical harmonic sum.
     **********************************************************************/
    Math::real operator()(real tau, real x, real y, real z) const {
      real f[] = {1, tau};
      real v = 0;
      real dummy;
      switch (_norm) {
      case full:
        v = SphericalEngine::Value<false, SphericalEngine::full, 2>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      case schmidt:
        v = SphericalEngine::Value<false, SphericalEngine::schmidt, 2>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      }
      return v;
    }

    /**
     * Compute a spherical harmonic sum with a correction term and its
     * gradient.
     *
     * @param[in] tau multiplier for correction coefficients \e C' and \e S'.
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
    Math::real operator()(real tau, real x, real y, real z,
                          real& gradx, real& grady, real& gradz) const {
      real f[] = {1, tau};
      real v = 0;
      switch (_norm) {
      case full:
        v = SphericalEngine::Value<true, SphericalEngine::full, 2>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      case schmidt:
        v = SphericalEngine::Value<true, SphericalEngine::schmidt, 2>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      }
      return v;
    }

    /**
     * Create a CircularEngine to allow the efficient evaluation of several
     * points on a circle of latitude at a fixed value of \e tau.
     *
     * @param[in] tau the multiplier for the correction coefficients.
     * @param[in] p the radius of the circle.
     * @param[in] z the height of the circle above the equatorial plane.
     * @param[in] gradp if true the returned object will be able to compute the
     *   gradient of the sum.
     * @return the CircularEngine object.
     *
     * SphericalHarmonic1::operator()() exchanges the order of the sums in the
     * definition, i.e., sum(n = 0..N)[sum(m = 0..n)[...]] becomes sum(m =
     * 0..N)[sum(n = m..N)[...]].  SphericalHarmonic1::Circle performs the
     * inner sum over degree \e n (which entails about \e N<sup>2</sup>
     * operations).  Calling CircularEngine::operator()() on the returned
     * object performs the outer sum over the order \e m (about \e N
     * operations).
     *
     * See SphericalHarmonic::Circle for an example of its use.
     **********************************************************************/
    CircularEngine Circle(real tau, real p, real z, bool gradp) const {
      real f[] = {1, tau};
      switch (_norm) {
      case full:
        return gradp ?
          SphericalEngine::Circle<true, SphericalEngine::full, 2>
          (_c, f, p, z, _a) :
          SphericalEngine::Circle<false, SphericalEngine::full, 2>
          (_c, f, p, z, _a);
        break;
      case schmidt:
      default:                  // To avoid compiler warnings
        return gradp ?
          SphericalEngine::Circle<true, SphericalEngine::schmidt, 2>
          (_c, f, p, z, _a) :
          SphericalEngine::Circle<false, SphericalEngine::schmidt, 2>
          (_c, f, p, z, _a);
        break;
      }
    }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALHARMONIC1_HPP

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

  class CircularEngine;

  /**
   * \brief The evaluation engine for SphericalHarmonic
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
      full = 0,
      /**
       * Schmidt semi-normalized associated Legendre polynomials.  See
       * SphericalHarmonic::schmidt for documentation.
       *
       * @hideinitializer
       **********************************************************************/
      schmidt = 1,
    };
    /**
     * \brief Package up coefficients for SphericalEngine
     *
     * This packages up the \e C, \e S coefficients and information about how
     * the coefficients are stored into a single structure.  This allows a
     * vector of coeffs to be passed to SphericalEngine::Value.  This class alo
     * includes functions to aid indexing into \e C and \e S.
     **********************************************************************/
    class coeff {
    private:
      int _N, _nmx, _mmx;
      std::vector<real>::const_iterator _Cnm;
      std::vector<real>::const_iterator _Snm;
    public:
      /**
       * A default constructor
       **********************************************************************/
      coeff()
        : _N(-1)
        , _nmx(-1)
        , _mmx(-1) 
        , _Cnm(Z_.begin())
        , _Snm(Z_.begin()) {}
      /**
       * The general constructor.
       *
       * @param[in] C a vector of coefficients for the cosine terms.
       * @param[in] S a vector of coefficients for the sine terms.
       * @param[in] N the degree giving storage layout for \e C and \e S.
       * @param[in] nmx the maximum degree to be used.
       * @param[in] mmx the maximum order to be used.
       *
       * This requires \e N >= \e nmx >= \e mmx >= -1.  \e C and \e S must also
       * be large enough to hold the coefficients.  Otherwise an exception is
       * thrown.
       **********************************************************************/
      coeff(const std::vector<real>& C,
            const std::vector<real>& S,
            int N, int nmx, int mmx)
        : _N(N)
        , _nmx(nmx)
        , _mmx(mmx)
        , _Cnm(C.begin())
        , _Snm(S.begin() - (_N + 1))
      {
        if (!(_N >= _nmx && _nmx >= _mmx && _mmx >= -1))
          throw GeographicErr("Bad indices for coeff");
        if (!(index(_nmx, _mmx) < int(C.size()) &&
              index(_nmx, _mmx) < int(S.size()) + (_N + 1)))
          throw GeographicErr("Arrays too small in coeff");
      }
      /**
       * The constructor for full coefficient vectors.
       *
       * @param[in] C a vector of coefficients for the cosine terms.
       * @param[in] S a vector of coefficients for the sine terms.
       * @param[in] N the maximum degree and order.
       *
       * This requires \e N >= -1.  \e C and \e S must also be large enough to
       * hold the coefficients.  Otherwise an exception is thrown.
       **********************************************************************/
      coeff(const std::vector<real>& C,
            const std::vector<real>& S,
            int N)
        : _N(N)
        , _nmx(N)
        , _mmx(N)
        , _Cnm(C.begin())
        , _Snm(S.begin() - (_N + 1))
      {
        if (!(_N >= -1))
          throw GeographicErr("Bad indices for coeff");
        if (!(index(_nmx, _mmx) < int(C.size()) &&
              index(_nmx, _mmx) < int(S.size()) + (_N + 1)))
          throw GeographicErr("Arrays too small in coeff");
      }
      /**
       * @return \e N the degree giving storage layout for \e C and \e S.
       **********************************************************************/
      inline int N() const throw() { return _N; }
      /**
       * @return \e nmx the maximum degree to be used.
       **********************************************************************/
      inline int nmx() const throw() { return _nmx; }
      /**
       * @return \e mmx the maximum order to be used.
       **********************************************************************/
      inline int mmx() const throw() { return _mmx; }
      /**
       * The one-dimensional index into \e C and \e S.
       *
       * @param[in] n the degree.
       * @param[in] m the order.
       * @return the one-dimensional index.
       **********************************************************************/
      inline int index(int n, int m) const throw()
      { return m * _N - m * (m - 1) / 2 + n; }
      /**
       * An element of \e C.
       *
       * @param[in] k the one-dimensional index.
       * @return the value of the \e C coefficient.
       **********************************************************************/
      inline Math::real Cv(int k) const { return *(_Cnm + k); }
      /**
       * An element of \e S.
       *
       * @param[in] k the one-dimensional index.
       * @return the value of the \e S coefficient.
       **********************************************************************/
      inline Math::real Sv(int k) const { return *(_Snm + k); }
      /**
       * An element of \e C with checking.
       *
       * @param[in] k the one-dimensional index.
       * @param[in] n the requested degree.
       * @param[in] m the requested order.
       * @param[in] f a multiplier.
       * @return the value of the \e C coefficient multiplied by \e f in \e n
       *   and \e m are in range else 0.
       **********************************************************************/
      inline Math::real Cv(int k, int n, int m, real f) const
      { return m > _mmx || n > _nmx ? 0 : *(_Cnm + k) * f; }
      /**
       * Aan element of \e S with checking.
       *
       * @param[in] k the one-dimensional index.
       * @param[in] n the requested degree.
       * @param[in] m the requested order.
       * @param[in] f a multiplier.
       * @return the value of the \e S coefficient multiplied by \e f in \e n
       *   and \e m are in range else 0.
       **********************************************************************/
      inline Math::real Sv(int k, int n, int m, real f) const
      { return m > _mmx || n > _nmx ? 0 : *(_Snm + k) * f; }
    };

    /**
     * Evaluate a spherical harmonic sum and its gradient.
     *
     * @tparam gradp should the gradient be calculated.
     * @tparam norm the normalization for the associated Legendre polynomials.
     * @tparam L the number of terms in the coefficents.
     * @param[in] c an array of coeff objects.
     * @param[in] f array of coefficient multipliers.  f[0] should be 1.
     * @param[in] x the \e x component of the cartesian position.
     * @param[in] y the \e y component of the cartesian position.
     * @param[in] z the \e z component of the cartesian position.
     * @param[in] a the normalizing radius.
     * @param[out] gradx the \e x component of the gradient.
     * @param[out] grady the \e y component of the gradient.
     * @param[out] gradz the \e z component of the gradient.
     * @result the spherical harmonic sum.
     *
     * See the SphericalHarmonic class for the definition of the sum.
     **********************************************************************/
    template<bool gradp, normalization norm, int L>
      static Math::real Value(const coeff c[], const real f[],
                              real x, real y, real z, real a,
                              real& gradx, real& grady, real& gradz);

    /**
     * Create a CircularEngine object
     *
     * @tparam gradp should the gradient be calculated.
     * @tparam norm the normalization for the associated Legendre polynomials.
     * @tparam L the number of terms in the coefficents.
     * @param[in] c an array of coeff objects.
     * @param[in] f array of coefficient multipliers.  f[0] should be 1.
     * @param[in] p the radius of the circle = sqrt(\e x<sup>2</sup> + \e
     *   y<sup>2</sup>).
     * @param[in] z the height of the circle.
     * @param[in] a the normalizing radius.
     * @result the CircularEngine object.
     **********************************************************************/
    template<bool gradp, SphericalEngine::normalization norm, int L>
      static CircularEngine Circle(const coeff c[], const real f[],
                                   real p, real z, real a);
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALENGINE_HPP

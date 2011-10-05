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
  private:
  public:
    typedef Math::real work;
    typedef Math::real T;
    /**
     * Compute a spherical harmonic sum.
     *
     * @param[in] C vector of coefficients for cosine terms.
     * @param[in] S vector of coefficients for sine terms.
     * @param[in] N the maximum order of the sum
     * @param[in] t the cosine of the geocentric co-latitude
     * @param[in] u the cosine of the geocentric co-latitude
     * @param[in] clam the cosine of the longitude
     * @param[in] slam the cosine of the longitude
     * @param[in] q the ratio \e a/\e r.
     * @return \e U the harmonic sum
     *
     * \e C and \e S must contain (at least) (\e N + 1)(\e N + 2)/2 elements.
     * These are stored in column major order.
     **********************************************************************/
    static Math::real Value(const std::vector<double>& C,
                            const std::vector<double>& S,
                            int N,
                            work t,    // cos(theta)
                            work u,    // sin(theta),
                            work clam, // cos(lambda)
                            work slam, // sin(lambda)
                            work q     // a/r
                            );
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALHARMONIC_HPP

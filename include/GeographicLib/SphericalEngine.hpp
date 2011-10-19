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
      static Math::real Value(const coeff c[], const real f[],
                              real x, real y, real z, real a,
                              real& gradx, real& grady, real& gradz);

    template<bool gradp, SphericalEngine::normalization norm, int L>
      static CircularEngine Circle(const coeff c[], const real f[],
                                   real p, real z, real a);
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALENGINE_HPP

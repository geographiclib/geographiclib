/**
 * \file DST.hpp
 * \brief Header for GeographicLib::DST class
 *
 * Copyright (c) Charles Karney (2022) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_DST_HPP)
#define GEOGRAPHICLIB_DST_HPP 1

#include <GeographicLib/Constants.hpp>

#include <functional>
#include <vector>
#include <memory>

template <typename scalar_t>
class kissfft;

namespace GeographicLib {

  /**
   * \brief Discrete sine transforms
   *
   * This decomposes periodic functions \f$ f(\sigma) \f$ (period \f$ 2\pi \f$)
   * which are odd about \f$ \sigma = 0 \f$ and even about \f$ \sigma = \frac12
   * \pi \f$ into a Fourier series
   * \f[
   *   f(\sigma) = \sum_{l=0}^{N-1} F_l \sin\bigl((2l+1)\sigma\bigr)
   * \f]
   *
   * Example of use:
   * \include example-DST.cpp
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT DST {
  private:
    typedef Math::real real;
    int _N;
    typedef kissfft<real> fft_t;
    std::shared_ptr<fft_t> _fft;
    mutable std::vector<real> _data;
    mutable std::vector<real> _temp;
    // Implement DST-III (centerp = false) or DST-IV (centerp = true)
    void fft_transform(real F[], bool centerp) const;
    // Add another N terms to F
    void fft_transform2(real F[]) const;
  public:
    /**
     * Constructor specifying the number of points to use.
     * @param[in] N the number of points to use.
     **********************************************************************/
    DST(int N = 0);

    /**
     * Reset the given number of points.
     * @param[in] N the number of points to use.
     **********************************************************************/
    void reset(int N);

    /**
     * Return the number of points.
     * @return the number of points to use.
     **********************************************************************/
    int N() const { return _N; }

    /**
     * Determine first \e N terms in the Fourier series
     * @param[in] f the function used for evaluation.
     * @param[out] F the first N coefficients of the Fourier series.
     *
     * The evaluates \f$ f(\sigma) \f$ at \f$ \sigma = i \pi / (2 N) \f$ for
     * integer \f$ i \in (0, N] \f$.
     **********************************************************************/
    void transform(std::function<real(real)> f, real F[]) const;

    /**
     * Refine the Fourier series by doubling the number of points sampled
     * @param[in] f the function used for evaluation.
     * @param[inout] F the refined transform based on 2\e N points.
     *
     * The evaluates \f$ f(\sigma) \f$ at additional points \f$ \sigma = (i -
     * \frac12) \pi / (2 N) \f$ for integer \f$ i \in (0, N] \f$ and combines
     * this with \e oldF to compute \e newF.  This is equivalent to calling
     * transform with twice the value of \e N but is more efficient.
     **********************************************************************/
    void refine(std::function<real(real)> f, real F[]) const;

    /**
     * Evaluate the Fourier sum given the sine and cosine of the angle
     * @param[in] sinx sin&sigma;.
     * @param[in] cosx cos&sigma;.
     * @param[in] F the vector of Fourier coefficients.
     * @param[in] N the number of Fourier coefficients.
     * @return the value of the Fourier sum.
     **********************************************************************/
    static real eval(real sinx, real cosx, const real F[], int N);

    /**
     * Evaluate the integral of Fourier sum given the sine and cosine of the
     * angle
     * @param[in] sinx sin&sigma;.
     * @param[in] cosx cos&sigma;.
     * @param[in] F the vector of Fourier coefficients.
     * @param[in] N the number of Fourier coefficients.
     * @return the value of the integral.
     *
     * The constant of integration is chosen so that the integral is zero at
     * \f$ \sigma = \frac12\pi \f$.
     **********************************************************************/
    static real integral(real sinx, real cosx, const real F[], int N);
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_DST_HPP

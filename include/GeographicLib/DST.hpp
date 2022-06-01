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
    typedef kissfft<GeographicLib::Math::real> fft_t;
    mutable std::shared_ptr<fft_t> _fft;
    // Implement DST-III (centerp = false) or DST-IV (centerp = true)
    void fft_transform(const std::vector<Math::real>& in,
                       std::vector<Math::real>& out,
                       bool centerp) const;
    void fft_transform2(const std::vector<Math::real>& newin,
                        const std::vector<Math::real>& oldout,
                        std::vector<Math::real>& newout) const;
  public:
    /**
     * Constructor specifying the expected number of points to use.
     * @param[in] N the expected number of points to use.
     **********************************************************************/
    DST(unsigned N = 0);

    /**
     * Reserve space for a given number of points.
     * @param[in] N the expected number of points to use.
     **********************************************************************/
    void reserve(unsigned N);
    // void transform(const std::vector<Math::real>& x,
    //                std::vector<Math::real>& tx) const;
    /**
     * Determine first \e N terms in the Fourier series
     * @param[in] f the function used for evaluation.
     * @param[in] N the number of points to use.
     * @param[out] F the first N coefficients of the Fourier series.
     *
     * The evaluates \f$ f(\sigma) \f$ at \f$ \sigma = i \pi / (2 N) \f$ for
     * integer \f$ i \in (0, N] \f$.
     **********************************************************************/
    void transform(std::function<Math::real(Math::real)> f, int N,
                   std::vector<Math::real>& F) const;
    /**
     * Refine the Fourier series by doubling the number of points sampled
     * @param[in] f the function used for evaluation.
     * @param[in] oldF the transform based on  \e N points.
     * @param[out] newF the refined transform based on 2\e N points.
     *
     * The evaluates \f$ f(\sigma) \f$ at additional points \f$ \sigma = (i -
     * \frac12) \pi / (2 N) \f$ for integer \f$ i \in (0, N] \f$ and combines
     * this with \e oldF to compute \e newF.  This is equivalent to calling
     * transform with twice the value of \e N but is more efficient.
     **********************************************************************/
    void refine(std::function<Math::real(Math::real)> f,
                const std::vector<Math::real>& oldF,
                std::vector<Math::real>& newF) const;
    /**
     * Evaluate the Fourier sum given the sine and cosine of the angle
     * @param[in] F the vector of Fourier coefficients.
     * @param[in] sinx sin&sigma;.
     * @param[in] cosx cos&sigma;.
     * @return the value of the Fourier sum.
     **********************************************************************/
    static Math::real eval(const std::vector<Math::real>& F,
                           Math::real sinx, Math::real cosx);
    /**
     * Evaluate the Fourier sum given the angle
     * @param[in] F the vector of Fourier coefficients.
     * @param[in] x &sigma;.
     * @return the value of the Fourier sum.
     **********************************************************************/
    static Math::real evalx(const std::vector<Math::real>& F, Math::real x);
    /**
     * Evaluate the integral of Fourier sum given the sine and cosine of the
     * angle
     * @param[in] F the vector of Fourier coefficients.
     * @param[in] sinx sin&sigma;.
     * @param[in] cosx cos&sigma;.
     * @return the value of the integral.
     *
     * The constant of integration is chosen so that the integral is zero at
     * \f$ \sigma = \frac12\pi \f$.
     **********************************************************************/
    static Math::real integral(const std::vector<Math::real>& F,
                               Math::real sinx, Math::real cosx);
    /**
     * Evaluate the integral of the Fourier sum given the angle
     * @param[in] F the vector of Fourier coefficients.
     * @param[in] x &sigma;.
     * @return the value of the integral.
     *
     * The constant of integration is chosen so that the integral is zero at
     * \f$ \sigma = \frac12\pi \f$.
     **********************************************************************/
    static Math::real integralx(const std::vector<Math::real>& F, Math::real x);
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_DST_HPP

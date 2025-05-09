/**
 * \file Angle.cpp
 * \brief Implementation for the GeographicLib::Angle class.
 *
 * This file is an implementation of the methods described in
 * - C. F. F. Karney,
 *   <a href="https://doi.org/10.1080/00396265.2023.2217604">
 *   On auxiliary latitudes,</a>
 *   Survey Review 56(395), 165--180 (2024);
 *   preprint
 *   <a href="https://arxiv.org/abs/2212.05818">arXiv:2212.05818</a>.
 * .
 * Copyright (c) Charles Karney (2022-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "Angle.hpp"
#include <iostream>

namespace GeographicLib {

  using namespace std;

  Math::real Angle::rnd(real x) {
    // This value of z more-or-less matches the value z = 1/16 in
    // Math::AngRound (where the argument is in degrees).
    static const real z = 1/real(1024);
    GEOGRAPHICLIB_VOLATILE real y = fabs(x);
    GEOGRAPHICLIB_VOLATILE real w = z - y;
    // The compiler mustn't "simplify" z - (z - y) to y
    y = w > 0 ? z - w : y;
    return copysign(y, x);
  }

} // namespace GeographicLib

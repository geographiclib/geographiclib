/**
 * \file Triaxial.hpp
 * \brief Header for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2022-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIAL_HPP)
#define GEOGRAPHICLIB_TRIAXIAL_HPP 1

#include <GeographicLib/Constants.hpp>

namespace GeographicLib {
  class GEOGRAPHICLIB_EXPORT Triaxial {
  private:
    typedef Math::real real;
    real _a, _b, _c;              // semi-axes
    real _e2, _k2, _kp2, _k, _kp;
  public:
  Triaxial(real a, real b, real c);
  };
} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIAL_HPP


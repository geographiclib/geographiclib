/**
 * \file Triaxial.cpp
 * \brief Implementation for GeographicLib::Triaxial class
 *
 * Copyright (c) Charles Karney (2022) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "Triaxial.hpp"
#include <iostream>

#include "kissfft.hh"

namespace GeographicLib {

  using namespace std;

  Triaxial::Triaxial(Math::real a, Math::real b, Math::real c)
    : _a(a)
    , _b(b)
    , _c(c)
  {
    real s = (_a - _c) * (_a + _c);
    _e2 = s / Math::sq(_b);
    if (s == 0)
      // The sphere is a nonuniform limit, we can pick any values in [0,1]
      // s.t. k2 + kp2 = 1.  Here we choose to treat the sphere as an
      // oblate ellipsoid.
        _kp2 = 0; _k2 = 1 - _kp2;
  } else {
    _kp2 = (_a - _b) * (_a + _b) / s;
    _k2  = (_b - _c) * (_b + _c) / s;
  }
  _k = sqrt(_k2); _kp = sqrt(_kp2);

} // namespace GeographicLib

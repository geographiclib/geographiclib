/**
 * \file AuxAngle.cpp
 * \brief Implementation for the GeographicLib::experimental::AuxAngle class.
 *
 * \note This is just sample code.  It is not part of GeographicLib itself.
 *
 * This file is an implementation of the methods described in
 * - C. F. F. Karney,
 *   On auxiliary latitudes,
 *   Technical Report, SRI International, December 2022.
 *   https://arxiv.org/abs/2212.05818
 * .
 * Copyright (c) Charles Karney (2022-2023) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "AuxAngle.hpp"

namespace GeographicLib {
namespace experimental {

  using namespace std;

  template<typename T>
  AuxAngle<T> AuxAngle<T>::NaN() {
    return AuxAngle(std::numeric_limits<T>::quiet_NaN(),
                    std::numeric_limits<T>::quiet_NaN());
  }

  template<typename T>
  AuxAngle<T> AuxAngle<T>::normalized() const {
    if ( isnan( tan() ) ||
         (fabs(_y) > std::numeric_limits<T>::max()/2 &&
          fabs(_x) > std::numeric_limits<T>::max()/2) )
      // deal with
      // (0,0), (inf,inf), (nan,nan), (nan,x), (y,nan), (toobig,toobig)
      return NaN();
    T r = hypot(_y, _x),
      y = _y/r, x = _x/r;
    // deal with r = inf, then one of y,x becomes 1
    if (isnan(y)) y = copysign(T(1), _y);
    if (isnan(x)) x = copysign(T(1), _x);
    return AuxAngle(y, x);
  }

  template<typename T>
  AuxAngle<T> AuxAngle<T>::copyquadrant(const AuxAngle& p) const {
    return AuxAngle(copysign(y(), p.y()), copysign(x(), p.x()));
  }

  template<typename T>
  AuxAngle<T>& AuxAngle<T>::operator+=(const AuxAngle& p) {
    // Do nothing if p.tan() == 0 to preserve signs of y() and x()
    if (p.tan() != 0) {
      T x = _x * p._x - _y * p._y;
      _y = _y * p._x + _x * p._y;
      _x = x;
    }
    return *this;
  }

  /// \cond SKIP
  // Instantiate
  template class AuxAngle<Math::real>;
#if GEOGRAPHICLIB_PRECISION != 2
  template class AuxAngle<double>;
#endif
  /// \endcond

} // namespace experimental
} // namespace GeographicLib

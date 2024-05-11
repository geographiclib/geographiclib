/**
 * \file TriaxialLin.hpp
 * \brief Header for GeographicLib::TriaxialLin class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_TRIAXIALLINE_HPP)
#define GEOGRAPHICLIB_TRIAXIALLINE_HPP 1

#include "Trigfun.hpp"
#include "Triaxial.hpp"
#include <GeographicLib/AuxAngle.hpp>

namespace GeographicLib {
  class GEOGRAPHICLIB_EXPORT TriaxialLine {
  private:
    typedef Math::real real;
    Triaxial _t;
    AuxAngle _bet1, _omg1, _alp1;
    real _gam;
    geod_fun _fbet, _fomg;
    dist_fun _gbet, _gomg;
  public:
    TriaxialLine(const Triaxial& t,
                 const AuxAngle& bet1,
                 const AuxAngle& omg1,
                 const AuxAngle& alp1);
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_TRIAXIALLINE_HPP

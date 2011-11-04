/**
 * \file MagneticCircle.cpp
 * \brief Implementation for GeographicLib::MagneticCircle class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/MagneticCircle.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <GeographicLib/Geocentric.hpp>

#define GEOGRAPHICLIB_MAGNETICCIRCLE_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_MAGNETICCIRCLE_CPP)
RCSID_DECL(GEOGRAPHICLIB_MAGNETICCIRCLE_HPP)

#define MAGNETIC_DEFAULT_PATH "/home/ckarney/geographiclib/magnetic"

namespace GeographicLib {

  using namespace std;

  void MagneticCircle::Field(real lon, bool diffp,
                             real& Bx, real& By, real& Bz,
                             real& Bxt, real& Byt, real& Bzt) const {
    lon = lon >= 180 ? lon - 360 : (lon < -180 ? lon + 360 : lon);
    real
      lam = lon * Math::degree<real>(),
      clam = std::abs(lam) ==   90 ? 0 : cos(lam),
      slam =          lam  == -180 ? 0 : sin(lam);
    real M[Geocentric::dim2_];
    Geocentric::Rotation(_sphi, _cphi, slam, clam, M);
    real BX0, BY0, BZ0, BX1, BY1, BZ1; // Components in geocentric basis
    _circ0(clam, slam, BX0, BY0, BZ0);
    _circ1(clam, slam, BX1, BY1, BZ1);
    if (_interpolate) {
      BX1 = (BX1 - BX0) / _dt0;
      BY1 = (BY1 - BY0) / _dt0;
      BZ1 = (BZ1 - BZ0) / _dt0;
    }
    BX0 += _t * BX1;
    BY0 += _t * BY1;
    BZ0 += _t * BZ1;
    if (diffp) {
      Bxt = - _a * (M[0] * BX1 + M[3] * BY1 + M[6] * BZ1);
      Byt = - _a * (M[1] * BX1 + M[4] * BY1 + M[7] * BZ1);
      Bzt = - _a * (M[2] * BX1 + M[5] * BY1 + M[8] * BZ1);
    }
    Bx = - _a * (M[0] * BX0 + M[3] * BY0 + M[6] * BZ0);
    By = - _a * (M[1] * BX0 + M[4] * BY0 + M[7] * BZ0);
    Bz = - _a * (M[2] * BX0 + M[5] * BY0 + M[8] * BZ0);
  }

} // namespace GeographicLib

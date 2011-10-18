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
    real BX, BY, BZ;            // Components in geocentric basis
    _circa(clam, slam, BX, BY, BZ);
    if (diffp) {
      real BXt, BYt, BZt;
      _circb(clam, slam, BXt, BYt, BZt);
      Bxt = - _a * (M[0] * BXt + M[3] * BYt + M[6] * BZt);
      Byt = - _a * (M[1] * BXt + M[4] * BYt + M[7] * BZt);
      Bzt = - _a * (M[2] * BXt + M[5] * BYt + M[8] * BZt);
    }
    Bx = - _a * (M[0] * BX + M[3] * BY + M[6] * BZ);
    By = - _a * (M[1] * BX + M[4] * BY + M[7] * BZ);
    Bz = - _a * (M[2] * BX + M[5] * BY + M[8] * BZ);
  }

} // namespace GeographicLib

/**
 * \file TriaxialLine.cpp
 * \brief Implementation for GeographicLib::TriaxialLine class
 *
 * Copyright (c) Charles Karney (2024) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include "TriaxialLine.hpp"
#include <iostream>

namespace GeographicLib {

  using namespace std;

  TriaxialLine::TriaxialLine(const Triaxial& t,
                             const AuxAngle& bet1,
                             const AuxAngle& omg1,
                             const AuxAngle& alp1)
    : _t(t)
    , _bet1(bet1.normalized())
    , _omg1(omg1.normalized())
    , _alp1(alp1.normalized())
    , _gam(_t.gamma(_bet1.x(), _omg1.y(), _alp1.y(), _alp1.x()))
    , _fbet(_t.k2, _t.kp2, _t.e2, -_gam)
    , _fomg(_t.kp2, _t.k2, -_t.e2, _gam)
    , _gbet(_t.k2, _t.kp2, _t.e2, -_gam)
    , _gomg(_t.kp2, _t.k2, -_t.e2, _gam)
  {
  
    bool debug = true;
    if (debug) {
      real u = 1, du = sqrt(sqrt(numeric_limits<real>::epsilon()));
      if (1) {
      cout << "gam " << _gam << "\n";
      cout << "fbet " << _fbet.txp() << " " << _fbet(u) << " "
           << _fbet.deriv(u) - (_fbet(u+du/2) - _fbet(u-du/2))/du << " "
           << _fbet.inv(_fbet(u)) - u << "\n";
      cout << "fomg " << _fomg.txp() << " " << _fomg(u) << " "
           << _fomg.deriv(u) - (_fomg(u+du/2) - _fomg(u-du/2))/du << " "
           << _fomg.inv(_fomg(u)) - u << "\n";
      cout << "gbet " << _gbet.txp() << " " << _gbet(u) << " "
           << _gbet.deriv(u) - (_gbet(u+du/2) - _gbet(u-du/2))/du << "\n";
      cout << "gomg " << _gomg.txp() << " " << _gomg(u) << " "
           << _gomg.deriv(u) - (_gomg(u+du/2) - _gomg(u-du/2))/du << "\n";
      cout << "cnts "
           << _fbet.NCoeffs() << " " << _fbet.NCoeffsInv() << " "
           << _fomg.NCoeffs() << " " << _fomg.NCoeffsInv() << " "
           << _gbet.NCoeffs() << " " << _gomg.NCoeffs() << "\n";
      }
      if (0) {
      cout << "gam " << _gam << "\n";
      cout << "fbet " << " " << _fbet.txp() << _gbet.txp() << " "
           << _gbet.deriv(u) / _fbet.deriv(u) - _gbet.gfderiv(u)
           << "\n";
      cout << "fomg " << " " << _fomg.txp() << _gomg.txp() << " "
           << _gomg.deriv(u) / _fomg.deriv(u) - _gomg.gfderiv(u)
           << "\n";
      }
      if (0) {
        cout << "gam " << _gam << "\n";
        geod_fun fbeta(_t.k2, _t.kp2, _t.e2, -_gam, false);
        geod_fun fomga(_t.kp2, _t.k2, -_t.e2, _gam, false);
        dist_fun gbeta(_t.k2, _t.kp2, _t.e2, -_gam, false);
        dist_fun gomga(_t.kp2, _t.k2, -_t.e2, _gam, false);
        geod_fun fbetb(_t.k2, _t.kp2, _t.e2, -_gam, true);
        geod_fun fomgb(_t.kp2, _t.k2, -_t.e2, _gam, true);
        dist_fun gbetb(_t.k2, _t.kp2, _t.e2, -_gam, true);
        dist_fun gomgb(_t.kp2, _t.k2, -_t.e2, _gam, true);
        real x = fbeta(fbeta.fwd(u));
        cout << "fbet "
             << fbeta(fbeta.fwd(u)) -  fbetb(fbetb.fwd(u)) << " "
             << fbeta.rev(fbeta.inv(x)) -  fbetb.rev(fbetb.inv(x)) << "\n";
        x = fomga(fomga.fwd(u));
        cout << "fomg "
             << fomga(fomga.fwd(u)) -  fomgb(fomgb.fwd(u)) << " "
             << fomga.rev(fomga.inv(x)) -  fomgb.rev(fomgb.inv(x)) << "\n";
        cout << "gbet "
             << gbeta(gbeta.fwd(u)) -  gbetb(gbetb.fwd(u)) << "\n";
        cout << "gomg "
             << gomga(gomga.fwd(u)) -  gomgb(gomgb.fwd(u)) << "\n";
      }
    }
  }
} // namespace GeographicLib

/**
 * \file CircularEngine.cpp
 * \brief Implementation for GeographicLib::CircularEngine class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/CircularEngine.hpp>
#include <limits>
#include <iostream>

#define GEOGRAPHICLIB_CIRCULARENGINE_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_CIRCULARENGINE_CPP)
RCSID_DECL(GEOGRAPHICLIB_CIRCULARENGINE_HPP)

namespace GeographicLib {

  using namespace std;

  Math::real CircularEngine::Value(bool gradp, real coslam, real sinlam,
                                   real& gradx, real& grady, real& gradz) const
  {
    gradp = _gradp && gradp;

    // Initialize outer sum
    real vc  = 0, vc2  = 0, vs  = 0, vs2  = 0;   // v [N + 1], v [N + 2]
    // vr, vt, vl and similar w variable accumulate the sums for the
    // derivatives wrt r, theta, and lambda, respectively.
    real vrc = 0, vrc2 = 0, vrs = 0, vrs2 = 0;   // vr[N + 1], vr[N + 2]
    real vtc = 0, vtc2 = 0, vts = 0, vts2 = 0;   // vt[N + 1], vt[N + 2]
    real vlc = 0, vlc2 = 0, vls = 0, vls2 = 0;   // vl[N + 1], vl[N + 2]
    for (int m = _M; m >= 0; --m) {   // m = M .. 0
      // Now Sc[m] = wc, Ss[m] = ws
      // Sc'[m] = wtc, Ss'[m] = wtc
      if (m) {
        real v, A, B;           // alpha[m], beta[m + 1]
        switch (_norm) {
        case full:
          v = 2 * real(2 * m + 3) / (m + 1);
          A = coslam * sqrt(v) * _uq;
          B = - sqrt((v * (2 * m + 5)) / (8 * (m + 2))) * _uq2;
          break;
        case schmidt:
          v = 2 * real(2 * m + 1) / (m + 1);
          A = coslam * sqrt(v) * _uq;
          B = - sqrt((v * (2 * m + 3)) / (8 * (m + 2))) * _uq2;
          break;
        default:
          A = B = 0;
        }
        v = A * vc  + B * vc2  +  _wc[m] ; vc2  = vc ; vc  = v;
        v = A * vs  + B * vs2  +  _ws[m] ; vs2  = vs ; vs  = v;
        if (gradp) {
          v = A * vrc + B * vrc2 +  _wrc[m]; vrc2 = vrc; vrc = v;
          v = A * vrs + B * vrs2 +  _wrs[m]; vrs2 = vrs; vrs = v;
          v = A * vtc + B * vtc2 +  _wtc[m]; vtc2 = vtc; vtc = v;
          v = A * vts + B * vts2 +  _wts[m]; vts2 = vts; vts = v;
          v = A * vlc + B * vlc2 + m*_ws[m]; vlc2 = vlc; vlc = v;
          v = A * vls + B * vls2 - m*_wc[m]; vls2 = vls; vls = v;
        }
      } else {
        real A, B, qs;
        switch (_norm) {
        case full:
          A = sqrt(real(3)) * _uq;       // F[1]/(q*coslam) or F[1]/(q*sinlam)
          B = - sqrt(real(15)/4) * _uq2; // beta[1]/q
          break;
        case schmidt:
          A = _uq;
          B = - sqrt(real(3)/4) * _uq2;
          break;
        default:
          A = B = 0;
        }
        qs = _q / _scale;
        vc = qs * (_wc[m] + A * (coslam * vc + sinlam * vs ) + B * vc2);
        if (gradp) {
          qs /= _r;
          // The components of the gradient in circular coordinates are
          // r: dV/dr
          // theta: 1/r * dV/dtheta
          // lambda: 1/(r*u) * dV/dlambda
          vrc =    -qs*(_wrc[m] + A * (coslam * vrc + sinlam * vrs) + B * vrc2);
          vtc = -_u*qs*(_wtc[m] + A * (coslam * vtc + sinlam * vts) + B * vtc2);
          vlc =  qs/_u*(          A * (coslam * vlc + sinlam * vls) + B * vlc2);
        }
      }
    }

    if (gradp) {
      // Rotate into cartesian (geocentric) coordinates
      gradx = coslam * (_u * vrc + _t * vtc) - sinlam * vlc;
      grady = sinlam * (_u * vrc + _t * vtc) + coslam * vlc;
      gradz =           _t * vrc - _u * vtc                ;
    }
    return vc;
  }

} // namespace GeographicLib

/**
 * \file DAuxLatitude.cpp
 * \brief Implementation for the GeographicLib::DAuxLatitude class.
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

#include <GeographicLib/DAuxLatitude.hpp>
#include <GeographicLib/EllipticFunction.hpp>

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions
#  pragma warning (disable: 4127)
#endif

namespace GeographicLib {

  using namespace std;

  Math::real DAuxLatitude::DRectifying(const AuxAngle& phi1,
                                       const AuxAngle& phi2)
    const {
    // Stipulate that phi1 and phi2 are in [-90d, 90d]
    real x = phi1.radians(), y = phi2.radians();
    if (x == y) {
      real d;
      AuxAngle mu1(base::Rectifying(phi1, &d));
      real tphi1 = phi1.tan(), tmu1 = mu1.tan();
      return
        isfinite(tphi1) ? d * Math::sq(base::sc(tphi1)/base::sc(tmu1)) : 1/d;
    } else if (x * y < 0)
      return (base::Rectifying(phi2).radians() -
              base::Rectifying(phi1).radians()) / (y - x);
    else {
      AuxAngle bet1(base::Parametric(phi1)), bet2(base::Parametric(phi2));
      real dEdbet = DE(bet1, bet2), dbetdphi = DParametric(phi1, phi2);
      return base::_b * dEdbet / base::RectifyingRadius(true) * dbetdphi;
    }
  }

  Math::real DAuxLatitude::DParametric(const AuxAngle& phi1,
                                       const AuxAngle& phi2)
    const {
    real tx = phi1.tan(), ty = phi2.tan(), r;
    // DbetaDphi = Datan(fm1*tx, fm1*ty) * fm1 / Datan(tx, ty)
    // Datan(x, y) = 1/(1 + x^2),                       for x = y
    //             = (atan(y) - atan(x)) / (y-x),       for x*y < 0
    //             = atan( (y-x) / (1 + x*y) ) / (y-x), for x*y > 0
    if (!(tx * ty >= 0))        // This includes, e.g., tx = 0, ty = inf
      r = (atan(base::_fm1 * ty) - atan(base::_fm1 * tx)) /
        (atan(ty) - atan(tx));
    else if (tx == ty) {        // This includes the case tx = ty = inf
      tx *= tx;
      if (tx <= 1)
        r = base::_fm1 * (1 + tx) / (1 + base::_e2m1 * tx);
      else {
        tx = 1/tx;
        r = base::_fm1 * (1 + tx) / (base::_e2m1 + tx);
      }
    } else {
      if (tx * ty <= 1)
        r = atan2(base::_fm1 * (ty - tx), 1 + base::_e2m1 * tx * ty)
          / atan2(        ty - tx , 1 +         tx * ty);
      else {
        tx = 1/tx; ty = 1/ty;
        r = atan2(base::_fm1 * (ty - tx), base::_e2m1 + tx * ty)
          / atan2(        ty - tx ,   1   + tx * ty);
      }
    }
    return r;
  }

  Math::real DAuxLatitude::DE(const AuxAngle& X, const AuxAngle& Y) const {
    AuxAngle Xn(X.normalized()), Yn(Y.normalized());
    // We assume that X and Y are in [-90d, 90d] and have the same sign
    // If not we would include
    //    if (Xn.y() * Yn.y() < 0)
    //      return d != 0 ? (E(X) - E(Y)) / d : 1;

    // The general formula fails for x = y = 0d and x = y = 90d.  Probably this
    // is fixable (the formula works for other x = y.  But let's also stipulate
    // that x != y .

    // Make both positive, so we can do the swap a <-> b trick
    Xn.y() = fabs(Xn.y()); Yn.y() = fabs(Yn.y());
    real x = Xn.radians(), y = Yn.radians(), d = y - x,
      sx = Xn.y(), sy = Yn.y(), cx = Xn.x(), cy = Yn.x(),
      k2;
    // Switch prolate to oblate; we then can use the formulas for k2 < 0
    if (false && base::_f < 0) {
      d = -d; swap(sx, cx); swap(sy, cy);
      k2 = base::_e2;
    } else {
      k2 = -base::_e12;
    }
    // See DLMF: Eqs (19.11.2) and (19.11.4) letting
    // theta -> x, phi -> -y, psi -> z
    //
    // (E(y) - E(x)) / d = E(z)/d - k2 * sin(x) * sin(y) * sin(z)/d
    //                   = (E(z)/sin(z) - k2 * sin(x) * sin(y)) * sin(z)/d
    // tan(z/2) = (sin(x)*Delta(y) - sin(y)*Delta(x)) / (cos(x) + cos(y))
    //          = d * Dsin(x,y) * (sin(x) + sin(y))/(cos(x) + cos(y)) /
    //             (sin(x)*Delta(y) + sin(y)*Delta(x))
    //          = t = d * Dt
    // Delta(x) = sqrt(1 - k2 * sin(x)^2)
    // sin(z) = 2*t/(1+t^2); cos(z) = (1-t^2)/(1+t^2)
    real Dt = Dsin(x, y) * (sx + sy) /
      ((cx + cy) * (sx * sqrt(1 - k2 * sy*sy) + sy * sqrt(1 - k2 * sx*sx))),
      t = d * Dt, Dsz = 2 * Dt / (1 + t*t),
      sz = d * Dsz, cz = (1 - t) * (1 + t) / (1 + t*t),
      sz2 = sz*sz, cz2 = cz*cz, dz2 = 1 - k2 * sz2,
      // E(z)/sin(z)
      Ezbsz = (EllipticFunction::RF(cz2, dz2, 1)
               - k2 * sz2 * EllipticFunction::RD(cz2, dz2, 1) / 3);
    return (Ezbsz - k2 * sx * sy) * Dsz;
  }

  /// \cond SKIP
  Math::real DAuxLatitude::Datanhee(real x, real y) const {
    // atan(e*sn(tphi))/e:
    //  Datan(e*sn(x),e*sn(y))*Dsn(x,y)/Datan(x,y)
    // asinh(e1*sn(fm1*tphi)):
    //  Dasinh(e1*sn(fm1*x)), e1*sn(fm1*y)) *
    // e1 * Dsn(fm1*x, fm1*y) *fm1 / (e * Datan(x,y))
    // = Dasinh(e1*sn(fm1*x)), e1*sn(fm1*y)) *
    //  Dsn(fm1*x, fm1*y) / Datan(x,y)
    return base::_f < 0 ?
      Datan(base::_e * base::sn(x), base::_e * base::sn(y)) * Dsn(x, y) :
      Dasinh(base::_e1 * base::sn(base::_fm1 * x),
             base::_e1 * base::sn(base::_fm1 * y)) *
      Dsn(base::_fm1 * x, base::_fm1 * y);
  }
  /// \endcond

  Math::real DAuxLatitude::DIsometric(const AuxAngle& phi1,
                                      const AuxAngle& phi2)
    const {
    // psi = asinh(tan(phi)) - e^2 * atanhee(tan(phi))
    real tphi1 = phi1.tan(), tphi2 = phi2.tan();
    return isnan(tphi1) || isnan(tphi2) ? numeric_limits<real>::quiet_NaN() :
      (isinf(tphi1) || isinf(tphi2) ? numeric_limits<real>::infinity() :
       (Dasinh(tphi1, tphi2) - base::_e2 * Datanhee(tphi1, tphi2)) /
       Datan(tphi1, tphi2));
  }

  Math::real DAuxLatitude::DConvert(int auxin, int auxout,
                                    const AuxAngle& zeta1,
                                    const AuxAngle& zeta2)
    const {
    int k = base::ind(auxout, auxin);
    if (k < 0) return numeric_limits<real>::quiet_NaN();
    if (auxin == auxout) return 1;
    if ( isnan(base::_c[base::Lmax * (k + 1) - 1]) )
      base::fillcoeff(auxin, auxout, k);
    AuxAngle zeta1n(zeta1.normalized()), zeta2n(zeta2.normalized());
    return 1 + DClenshaw(true, zeta2n.radians() - zeta1n.radians(),
                         zeta1n.y(), zeta1n.x(), zeta2n.y(), zeta2n.x(),
                         base::_c + base::Lmax * k, base::Lmax);
  }

  Math::real DAuxLatitude::DClenshaw(bool sinp, real Delta,
                                     real szet1, real czet1,
                                     real szet2, real czet2,
                                     const real c[], int K) {
    // Evaluate
    // (Clenshaw(sinp, szet2, czet2, c, K) -
    //  Clenshaw(sinp, szet1, czet1, c, K)) / Delta
    // or
    // sum(c[k] * (sin( (2*k+2) * zet2) - sin( (2*k+2) * zet2)), i, 0, K-1)
    //   / Delta
    // (if !sinp, then change sin->cos here.)
    //
    // Delta is EITHER 1, giving the plain difference OR (zeta2 - zeta1) in
    // radians, giving the divided difference.  Other values will give
    // nonsense.
    //
    int k = K;
    // suffices a b denote [1,1], [2,1] elements of matrix/vector
    real D2 = Delta * Delta,
      czetp = czet2 * czet1 - szet2 * szet1,
      szetp = szet2 * czet1 + czet2 * szet1,
      czetm = czet2 * czet1 + szet2 * szet1,
      // sin(zetm) / Delta
      szetmd =  (Delta == 1 ? szet2 * czet1 - czet2 * szet1 :
                 (Delta != 0 ? sin(Delta) / Delta : 1)),
      Xa =  2 * czetp * czetm,
      Xb = -2 * szetp * szetmd,
      u0a = 0, u0b = 0, u1a = 0, u1b = 0; // accumulators for sum
    for (--k; k >= 0; --k) {
      // temporary real = X . U0 - U1 + c[k] * I
      real ta = Xa * u0a + D2 * Xb * u0b - u1a + c[k],
        tb = Xb * u0a +      Xa * u0b - u1b;
      // U1 = U0; U0 = real
      u1a = u0a; u0a = ta;
      u1b = u0b; u0b = tb;
    }
    // P = U0 . F[0] - U1 . F[-1]
    // if sinp:
    //   F[0] = [ sin(2*zet2) + sin(2*zet1),
    //           (sin(2*zet2) - sin(2*zet1)) / Delta]
    //        = 2 * [ szetp * czetm, czetp * szetmd ]
    //   F[-1] = [0, 0]
    // else:
    //   F[0] = [ cos(2*zet2) + cos(2*zet1),
    //           (cos(2*zet2) - cos(2*zet1)) / Delta]
    //        = 2 * [ czetp * czetm, -szetp * szetmd ]
    //   F[-1] = [2, 0]
    real F0a = (sinp ? szetp :  czetp) * czetm,
      F0b = (sinp ? czetp : -szetp) * szetmd,
      Fm1a = sinp ? 0 : 1;  // Fm1b = 0;
    // Don't both to compute sum...
    // divided difference (or difference if Delta == 1)
    return 2 * (F0a * u0b + F0b * u0a  - Fm1a * u1b);
  }

} // namespace GeographicLib

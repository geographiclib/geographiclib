/**
 * \file PolarStereographic.cpp
 * \brief Implementation for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/PolarStereographic.hpp"

#define GEOGRAPHICLIB_POLARSTEREOGRAPHIC_CPP "$Id: PolarStereographic.cpp 6858 2010-09-06 13:52:31Z karney $"

RCSID_DECL(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real PolarStereographic::tol =
    real(0.1)*sqrt(numeric_limits<real>::epsilon());
  // Overflow value s.t. atan(overflow) = pi/2
  const Math::real PolarStereographic::overflow =
    1 / sq(numeric_limits<real>::epsilon());

  PolarStereographic::PolarStereographic(real a, real r, real k0)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
    , _Cx(exp(eatanhe(real(1))))
    , _c( (1 - _f) * _Cx )
    , _k0(k0)
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_k0 > 0))
      throw GeographicErr("Scale is not positive");
  }

  const PolarStereographic
  PolarStereographic::UPS(Constants::WGS84_a(), Constants::WGS84_r(),
                          Constants::UPS_k0());

  // This formulation converts to conformal coordinates by tau = tan(phi) and
  // tau' = tan(phi') where phi' is the conformal latitude.  The formulas are:
  //    tau = tan(phi)
  //    secphi = hypot(1, tau)
  //    sig = sinh(e * atanh(e * tau / secphi))
  //    taup = tan(phip) = tau * hypot(1, sig) - sig * hypot(1, tau)
  //    c = (1 - f) * exp(e * atanh(e))
  //
  // Forward:
  //   rho = (2*k0*a/c) / (hypot(1, taup) + taup)  (taup >= 0)
  //       = (2*k0*a/c) * (hypot(1, taup) - taup)  (taup <  0)
  //
  // Reverse:
  //   taup = ((2*k0*a/c) / rho - rho / (2*k0*a/c))/2
  //
  // Scale:
  //   k = (rho/a) * secphi * sqrt((1-e2) + e2 / secphi^2)
  //
  // In limit rho -> 0, tau -> inf, taup -> inf, secphi -> inf, secphip -> inf
  //   secphip = taup = exp(-e * atanh(e)) * tau = exp(-e * atanh(e)) * secphi

  void PolarStereographic::Forward(bool northp, real lat, real lon,
                                   real& x, real& y, real& gamma, real& k)
    const throw() {
    lat *= northp ? 1 : -1;
    real
      phi = lat * Constants::degree(),
      tau = lat > -90 ? tan(phi) : -overflow,
      secphi = Math::hypot(real(1), tau),
      sig = sinh( eatanhe(tau / secphi) ),
      taup = Math::hypot(real(1), sig) * tau - sig * secphi,
      rho =  Math::hypot(real(1), taup) + abs(taup);
    rho = taup >= 0 ? (lat < 90 ? 1/rho : 0) : rho;
    rho *= 2 * _k0 * _a / _c;
    k = lat < 90 ? (rho / _a) * secphi * sqrt(_e2m + _e2 / sq(secphi)) : _k0;
    lon = lon >= 180 ? lon - 360 : lon < -180 ? lon + 360 : lon;
    real
      lam = lon * Constants::degree();
    x = rho * (lon == -180 ? 0 : sin(lam));
    y = (northp ? -rho : rho) * (abs(lon) == 90 ? 0 : cos(lam));
    gamma = northp ? lon : -lon;
  }

  void PolarStereographic::Reverse(bool northp, real x, real y,
                                   real& lat, real& lon, real& gamma, real& k)
    const throw() {
    real
      rho = Math::hypot(x, y),
      t = rho / (2 * _k0 * _a / _c),
      taup = (1 / t - t) / 2,
      tau = taup * _Cx,
      stol = tol * max(real(1), abs(taup));
    // min iterations = 1, max iterations = 2; mean = 1.99
    for (int i = 0; i < numit; ++i) {
      real
        tau1 = Math::hypot(real(1), tau),
        sig = sinh( eatanhe( tau / tau1 ) ),
        sig1 =  Math::hypot(real(1), sig),
        dtau = - (sig1 * tau - sig * tau1 - taup) * (1 + _e2m * sq(tau)) /
        ( (sig1 * tau1 - sig * tau) * _e2m * tau1 );
      tau += dtau;
      if (abs(dtau) < stol)
        break;
    }
    real
      phi = atan(tau),
      secphi = Math::hypot(real(1), tau);
    k = rho != 0 ? (rho / _a) * secphi * sqrt(_e2m + _e2 / sq(secphi)) : _k0;
    lat = (northp ? 1 : -1) * (rho != 0 ? phi / Constants::degree() : 90);
    lon = -atan2( -x, northp ? -y : y ) / Constants::degree();
    gamma = northp ? lon : -lon;
  }

  void PolarStereographic::SetScale(real lat, real k) {
    if (!(k > 0))
      throw GeographicErr("Scale is not positive");
    if (!(abs(lat) <= 90))
      throw GeographicErr("Latitude must be in [-90d, 90d]");
    real x, y, gamma, kold;
    _k0 = 1;
    Forward(true, lat, 0, x, y, gamma, kold);
    _k0 *= k/kold;
  }

} // namespace GeographicLib

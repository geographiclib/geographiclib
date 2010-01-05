/**
 * \file PolarStereographic.cpp
 * \brief Implementation for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/PolarStereographic.hpp"

#define GEOGRAPHICLIB_POLARSTEREOGRAPHIC_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real PolarStereographic::tol =
    real(0.1)*sqrt(numeric_limits<real>::epsilon());

  PolarStereographic::PolarStereographic(real a, real r, real k0)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
      // _c = sqrt( pow(1 + _e, 1 + _e) * pow(1 - _e, 1 - _e) )
    , _c( sqrt(_e2m) * exp(eatanhe(real(1))) )
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

  void PolarStereographic::Forward(bool northp, real lat, real lon,
                                   real& x, real& y, real& gamma, real& k)
    const throw() {
    real theta = 90 - (northp ? lat : -lat); //  the colatitude
    real rho;
    theta *= Constants::degree();
    real
      // f = ((1 + e cos(theta))/(1 - e cos(theta)))^(e/2),
      ctheta = cos(theta),
      f = exp(eatanhe(ctheta)),
      t2 = 2 * tan(theta/2) * f, // Snyder (15-9) (t2 = 2 * t)
      m = sin(theta) / sqrt(1 - _e2 * sq(ctheta)); // Snyder (14-15)
    rho = _a * _k0 * t2 / _c;                      // Snyder (21-33)
    k = m < numeric_limits<real>::epsilon() ? _k0 : rho / (_a * m);
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
      t2 = rho * _c / (_a * _k0),
      theta = 2 * atan(t2/2);   // initial (spherical) estimate of colatitude
    // Solve from theta using Newton's method on Snyder (15-9) which converges
    // more rapidly than the iteration procedure given by Snyder (7-9).  First
    // rewrite as
    // v(theta) = 2 * tan(theta/2) * f - t2 = 0
    for (int i = 0; i < numit; ++i) {
      real
        ctheta = cos(theta),
        f = exp(eatanhe(ctheta)),
        v = 2 * tan(theta/2) * f - t2,
        // dv/dtheta
        dv = _e2m * f / ((1 - _e2 * sq(ctheta)) * sq(cos(theta/2))),
        dtheta = -v/dv;
      theta += dtheta;
      if (abs(dtheta) < tol)
        break;
    }
    lat = (northp ? 1 : -1) * (90 - theta / Constants::degree());
    // Result is in [-180, 180).  Assume atan2(0,0) = 0.
    lon = -atan2( -x, northp ? -y : y ) / Constants::degree();
    real m = sin(theta) / sqrt(1 - _e2 * sq(cos(theta)));
    k = m == 0 ? _k0 : rho / (_a * m);
    gamma = northp ? lon : -lon;
  }

  void PolarStereographic::SetScale(real lat, real k) {
    if (!(k > 0))
      throw GeographicErr("Scale is not positive");
    real x, y, gamma, kold;
    Forward(true, lat, 0, x, y, gamma, kold);
    _k0 *= k/kold;
  }

} // namespace GeographicLib

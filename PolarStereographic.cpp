/**
 * \file PolarStereographic.cpp
 * \brief Implementation for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008, 2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/PolarStereographic.hpp"
#include <limits>

#define GEOGRAPHICLIB_POLARSTEREOGRAPHIC_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real_t PolarStereographic::tol =
    real_t(0.1L)*sqrt(numeric_limits<real_t>::epsilon());

  PolarStereographic::PolarStereographic(real_t a, real_t r, real_t k0)
    throw()
    : _a(a)
    , _f(r != 0 ? 1 / r : 0)
    , _k0(k0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
      // _c = sqrt( pow(1 + _e, 1 + _e) * pow(1 - _e, 1 - _e) )
    , _c( sqrt(_e2m) * exp(eatanhe(real_t(1))) )
  {}

  const PolarStereographic
  PolarStereographic::UPS(Constants::WGS84_a(), Constants::WGS84_r(),
			  Constants::UPS_k0());

  void PolarStereographic::Forward(bool northp, real_t lat, real_t lon,
				   real_t& x, real_t& y,
				   real_t& gamma, real_t& k) const throw() {
    real_t theta = 90 - (northp ? lat : -lat); //  the colatitude
    real_t rho;
    theta *= Constants::degree();
    real_t
      // f = ((1 + e cos(theta))/(1 - e cos(theta)))^(e/2),
      ctheta = cos(theta),
      f = exp(eatanhe(ctheta)),
      t2 = 2 * tan(theta/2) * f, // Snyder (15-9) (t2 = 2 * t)
      m = sin(theta) / sqrt(1 - _e2 * sq(ctheta)); // Snyder (14-15)
    rho = _a * _k0 * t2 / _c;			   // Snyder (21-33)
    k = m < numeric_limits<real_t>::epsilon() ? _k0 : rho / (_a * m);
    lon = lon >= 180 ? lon - 360 : lon < -180 ? lon + 360 : lon;
    real_t
      lam = lon * Constants::degree();
    x = rho * (lon == -180 ? 0 : sin(lam));
    y = (northp ? -rho : rho) * (abs(lon) == 90 ? 0 : cos(lam));
    gamma = northp ? lon : -lon;
  }

  void PolarStereographic::Reverse(bool northp, real_t x, real_t y,
				   real_t& lat, real_t& lon,
				   real_t& gamma, real_t& k) const throw() {
    real_t
      rho = Math::hypot(x, y),
      t2 = rho * _c / (_a * _k0),
      theta = Constants::pi()/2;	// initial estimate of colatitude
    // Solve from theta using Newton's method on Snyder (15-9) which converges
    // more rapidly than the iteration procedure given by Snyder (7-9).  First
    // rewrite as
    // v(theta) = 2 * tan(theta/2) * f - t2 = 0
    for (int i = 0; i < numit; ++i) {
      real_t
	ctheta = cos(theta),
	f = exp(eatanhe(ctheta)),
	c2 = cos(theta/2),
	v = 2 * tan(theta/2) * f - t2,
	dv = _e2m * f / ((1 - _e2 * sq(ctheta)) * sq(c2)), // dv/dtheta
	dtheta = -v/dv;
      theta += dtheta;
      if (abs(dtheta) < tol)
	break;
    }
    lat = (northp ? 1 : -1) * (90 - theta / Constants::degree());
    // Result is in [-180, 180).  Assume atan2(0,0) = 0.
    lon = -atan2( -x, northp ? -y : y ) / Constants::degree();
    real_t m = sin(theta) / sqrt(1 - _e2 * sq(cos(theta)));
    k = m == 0 ? _k0 : rho / (_a * m);
    gamma = northp ? lon : -lon;
  }

} // namespace GeographicLib

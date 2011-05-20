/**
 * \file LambertConformalConic.cpp
 * \brief Implementation for GeographicLib::LambertConformalConic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/LambertConformalConic.hpp"

#define GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_CPP "$Id: LambertConformalConic.cpp 6856 2010-08-23 12:10:55Z karney $"

RCSID_DECL(GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real LambertConformalConic::eps = numeric_limits<real>::epsilon();
  const Math::real LambertConformalConic::eps2 = sqrt(eps);
  const Math::real LambertConformalConic::epsx = sq(eps);
  const Math::real LambertConformalConic::tol = real(0.1) * eps2;
  // Large enough value so that atan(sinh(ahypover)) = pi/2
  const Math::real LambertConformalConic::ahypover =
    real(numeric_limits<real>::digits) * log(real(numeric_limits<real>::radix))
    + 2;

  LambertConformalConic::LambertConformalConic(real a, real r,
                                               real stdlat, real k0)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_f < 1))
      throw GeographicErr("Minor radius is not positive");
    if (!(k0 > 0))
      throw GeographicErr("Scale is not positive");
    if (!(abs(stdlat) <= 90))
      throw GeographicErr("Standard latitude not in [-90, 90]");
    real
      phi = stdlat * Constants::degree(),
      sphi = sin(phi),
      cphi = abs(stdlat) != 90 ? cos(phi) : 0;
    Init(sphi, cphi, sphi, cphi, k0);
  }

  LambertConformalConic::LambertConformalConic(real a, real r,
                                               real stdlat1, real stdlat2,
                                               real k1)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_f < 1))
      throw GeographicErr("Minor radius is not positive");
    if (!(k1 > 0))
      throw GeographicErr("Scale is not positive");
    if (!(abs(stdlat1) <= 90))
      throw GeographicErr("Standard latitude 1 not in [-90, 90]");
    if (!(abs(stdlat2) <= 90))
      throw GeographicErr("Standard latitude 2 not in [-90, 90]");
    if (abs(stdlat1) == 90 || abs(stdlat2) == 90)
      if (!(stdlat1 == stdlat2))
        throw GeographicErr
          ("Standard latitudes must be equal is either is a pole");
    real
      phi1 = stdlat1 * Constants::degree(),
      phi2 = stdlat2 * Constants::degree();
    Init(sin(phi1), abs(stdlat1) != 90 ? cos(phi1) : 0,
         sin(phi2), abs(stdlat2) != 90 ? cos(phi2) : 0, k1);
  }

  LambertConformalConic::LambertConformalConic(real a, real r,
                                               real sinlat1, real coslat1,
                                               real sinlat2, real coslat2,
                                               real k1)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_f < 1))
      throw GeographicErr("Minor radius is not positive");
    if (!(k1 > 0))
      throw GeographicErr("Scale is not positive");
    if (coslat1 == 0 || coslat2 == 0)
      if (!(coslat1 == coslat2 && sinlat1 == sinlat2))
        throw GeographicErr
          ("Standard latitudes must be equal is either is a pole");
    Init(sinlat1, coslat1, sinlat2, coslat2, k1);
  }

  void LambertConformalConic::Init(real sphi1, real cphi1,
                                   real sphi2, real cphi2, real k1) throw() {
    // Determine hemisphere of tangent latitude
    _sign = sphi1 + sphi2 > 0 ? 1 : sphi1 + sphi2 < 0  ? -1 :
      atan2(sphi1, cphi1) + atan2(sphi2, cphi2) >= 0 ? 1 : -1;
    // Internally work with tangent latitude positive
    sphi1 *= _sign;
    sphi2 *= _sign;
    real
      m1 = mf(sphi1, cphi1), lt1 = logtf(sphi1, cphi1),
      m2 = mf(sphi2, cphi2), lt2 = logtf(sphi2, cphi2),
      sindiff = abs(sphi1 * cphi2 - cphi1 * sphi2);
    real sphi0, cphi0;          // Use phi0 = tangent latitude
    if (sindiff > eps2) {
      // sphi0 = Snyder's n, p 108, eq 15-8
      sphi0 = (log(m1) - log (m2)) / (lt1 - lt2);
      cphi0 = sindiff < real(0.7) ? sqrt(1 - sq(sphi0)) :
        // cos(phi0) = sqrt( - (sin(phi0) - 1) * (sin(phi0) + 1) )
        sqrt( -(logmtf(sphi1) - logmtf(sphi2)) / (lt1 - lt2)
              * (sphi0 + 1));
    } else {
      // Set phi0 = (phi1 + phi2)/2
      sphi0 = (sphi1 * sqrt((1 + cphi2)/(1 + cphi1)) +
               sphi2 * sqrt((1 + cphi1)/(1 + cphi2)))/2;
      cphi0 = (cphi1 * sqrt((1 + sphi2)/(1 + sphi1)) +
               cphi2 * sqrt((1 + sphi1)/(1 + sphi2)))/2;
    }
    _n = sphi0;                 // Snyder's n
    _nc = cphi0;                // sqrt(1 - sq(n))
    _lat0 = atan2(_sign * sphi0, cphi0) / Constants::degree();
    _lt0 = logtf(sphi0, cphi0); // Snyder's log(t0)
    _t0n = exp(_n * _lt0);      // Snyder's t0^n
    _t0nm1 = Math::expm1(_n * _lt0);      // Snyder's t0^n - 1
    // k1 * m1/t1^n = k1 * m2/t2^n = k1 * n * (Snyder's F)
    _scale = k1 *
      (_nc == 0 ?               // if phi0 = pi/2, n = t = m = 0, so take limit
       2/( sqrt(_e2m) * exp(eatanhe(real(1))) ) :
       (lt1 > lt2 ? m1 / exp(_n * lt1) : m2 / exp(_n * lt2)));
    // Scale at phi0
    _k0 = _scale * _t0n / mf(sphi0, cphi0);
  }

  const LambertConformalConic
  LambertConformalConic::Mercator(Constants::WGS84_a(), Constants::WGS84_r(),
                                  real(0), real(1));

  void LambertConformalConic::Forward(real lon0, real lat, real lon,
                                      real& x, real& y, real& gamma, real& k)
    const throw() {
    if (lon - lon0 > 180)
      lon -= lon0 - 360;
    else if (lon - lon0 <= -180)
      lon -= lon0 + 360;
    else
      lon -= lon0;
    lat *= _sign;
    real
      phi = lat * Constants::degree(),
      lam = lon * Constants::degree(),
      sphi = sin(phi), cphi = abs(lat) != 90 ? cos(phi) : 0,
      m = mf(sphi, cphi),
      lt = logtf(sphi, cphi),
      tn = exp(_n * lt),
      theta = _n * lam, stheta = sin(theta), ctheta = cos(theta);
    x = _a * _scale * tn * (_n > eps2 ? stheta/_n : lam);
    if (_n > real(0.5))
      y = (_t0n - tn * ctheta)/_n;
    else {
      // write as
      // _t0n * ( (1 - ctheta)/_n - (tn/_t0n - 1)/_n * ctheta )
      // _t0n * ( stheta^2/(1 + ctheta) / _n
      //          - ((t/_t0)^n - 1)/_n * ctheta )
      // _t0n * ( stheta^2/(1 + ctheta) / _n
      //          - (exp(_n * log(t/_t0)) - 1)/_n * ctheta )
      // _t0n * (ax - bx)
      real
        ax = (_n > eps2 ? sq(stheta) / _n : lam * theta) / (1 + ctheta),
        fx = lt - _lt0,
        bx = ctheta * (abs(_n * fx) > eps ? Math::expm1(_n * fx) / _n : fx);
      y = _t0n * (ax - bx);
    }
    y *= _a * _scale * _sign;
    k = _scale * (m == 0 && _nc == 0 && sphi > 0 ?
                  sqrt(_e2m) * exp(eatanhe(real(1))) / 2 :
                  tn/m);        // infinite if pole and _n < 1
    gamma = theta * _sign;
  }

  void LambertConformalConic::Reverse(real lon0, real x, real y,
                                      real& lat, real& lon,
                                      real& gamma, real& k)
    const throw() {
    y *= _sign;
    x /= _a * _scale;
    y /= _a * _scale;
    real
      /*
        rho = hypot(x, rho0 - y)
        rho0 = _a * _scale * _t0n / _n
        _n * rho = hypot(x * _n, _a * _scale * _t0n - y * _n)
        tn    = _n * rho / (_a * _scale)
              = hypot(_n * x/(_a * _scale), _t0n - _n * y/(_a * _scale))
        theta = atan2(_n * x/(_a * _scale), _t0n - _n * y/(_a * _scale))
        lam = theta/_n
      */
      x1 = _n * x, y1 = _n * y,
      tn = Math::hypot(x1, _t0n - y1), theta = atan2(x1, _t0n - y1),
      lam = _n != 0 ? theta/_n : x/_t0n;
    real q;                     // -log(t)
    if (_n != 0 && tn < real(0.5))
      q = -log(tn) / _n;
    else {
      real
        b1 = _t0nm1 - _n * y,
        tnm1 = (sq(x1) +  2 * b1 + sq(b1)) / (tn + 1); // tn - 1
      q = _n != 0 ? - Math::log1p(tnm1)/_n :
        -2 * (_lt0  - y) / (tn + 1); // _t0nm1 -> _n * _lt0
    }
    // Clip to [-ahypover, ahypover] to avoid overflow later
    q = q < ahypover ? (q > -ahypover ? q : -ahypover) : ahypover;
    // Recast Snyder's 15-9 as
    // q = q' - eatanh(tanh(q'))
    // where q = -log(t) and q' = asinh(tan(phi))
    // q is known and we solve for q' by Newton's method.
    // Write f(q') = q' - eatanh(tanh(q')) - q
    // f'(q') = (1 - e^2)/(1 - e^2 * tanh(q')^2)
    // Starting guess is q' = q.
    real qp = q;
    for (int i = 0; i < numit; ++i) {
      real
        t = tanh(qp),
        dqp = -(qp - eatanhe(t) - q) * (1 - _e2 * sq(t)) / _e2m;
      qp += dqp;
      if (abs(dqp) < tol)
        break;
    }
    double
      phi = _sign * atan(sinh(qp)),
      m = mf(tanh(qp), 1/cosh(qp));
    lat = phi / Constants::degree();
    lon = lam / Constants::degree();
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon + lon0 >= 180)
      lon += lon0 - 360;
    else if (lon + lon0 < -180)
      lon += lon0 + 360;
    else
      lon += lon0;
    gamma = _sign * theta / Constants::degree();
    k = _scale * (m == 0 && _nc == 0 && qp > 0 ?
                  sqrt(_e2m) * exp(eatanhe(real(1))) / 2 :
                  tn/m);        // infinite if pole and _n < 1
  }

  void LambertConformalConic::SetScale(real lat, real k) {
    if (!(k > 0))
      throw GeographicErr("Scale is not positive");
    real x, y, gamma, kold;
    Forward(0, lat, 0, x, y, gamma, kold);
    k /= kold;
    _scale *= k;
    _k0 *= k;
  }

} // namespace GeographicLib

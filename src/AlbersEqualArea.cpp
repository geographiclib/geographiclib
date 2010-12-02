/**
 * \file AlbersEqualArea.cpp
 * \brief Implementation for GeographicLib::AlbersEqualArea class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/AlbersEqualArea.hpp"

#define GEOGRAPHICLIB_ALBERSEQUALAREA_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_ALBERSEQUALAREA_CPP)
RCSID_DECL(GEOGRAPHICLIB_ALBERSEQUALAREA_HPP)

#include <iostream>
#include <iomanip>

namespace GeographicLib {

  using namespace std;

  const Math::real AlbersEqualArea::eps = numeric_limits<real>::epsilon();
  const Math::real AlbersEqualArea::epsx = sq(eps);
  const Math::real AlbersEqualArea::tol = real(0.1) * sqrt(eps);
  const Math::real AlbersEqualArea::ahypover =
    real(numeric_limits<real>::digits) * log(real(numeric_limits<real>::radix))
    + 2;

  AlbersEqualArea::AlbersEqualArea(real a, real r,
                                   real stdlat, real k0)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _fm(1 - _f)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
    , _qp(1 + _e2m * atanhee(real(1)))
    , _qx(_qp / ( 2 * _e2m ))
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
      phi = stdlat * Math::degree(),
      sphi = sin(phi),
      cphi = abs(stdlat) != 90 ? cos(phi) : 0;
    Init(sphi, cphi, sphi, cphi, k0);
  }

  AlbersEqualArea::AlbersEqualArea(real a, real r,
                                   real stdlat1, real stdlat2,
                                   real k1)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _fm(1 - _f)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
    , _qp(1 + _e2m * atanhee(real(1)))
    , _qx(_qp / ( 2 * _e2m ))
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
    real
      phi1 = stdlat1 * Math::degree(),
      phi2 = stdlat2 * Math::degree();
    Init(sin(phi1), abs(stdlat1) != 90 ? cos(phi1) : 0,
         sin(phi2), abs(stdlat2) != 90 ? cos(phi2) : 0, k1);
  }

  AlbersEqualArea::AlbersEqualArea(real a, real r,
                                               real sinlat1, real coslat1,
                                               real sinlat2, real coslat2,
                                               real k1)
    : _a(a)
    , _r(r)
    , _f(_r != 0 ? 1 / _r : 0)
    , _fm(1 - _f)
    , _e2(_f * (2 - _f))
    , _e(sqrt(abs(_e2)))
    , _e2m(1 - _e2)
    , _qp(1 + _e2m * atanhee(real(1)))
    , _qx(_qp / ( 2 * _e2m ))
  {
    if (!(_a > 0))
      throw GeographicErr("Major radius is not positive");
    if (!(_f < 1))
      throw GeographicErr("Minor radius is not positive");
    if (!(k1 > 0))
      throw GeographicErr("Scale is not positive");
    if (coslat1 == 0 && coslat2 == 0 && sinlat1 * sinlat2 <= 0)
      throw GeographicErr
        ("Standard latitudes cannot be opposite poles");
    Init(sinlat1, coslat1, sinlat2, coslat2, k1);
  }

  void AlbersEqualArea::Init(real sphi1, real cphi1,
                                   real sphi2, real cphi2, real k1) throw() {
    {
      real r;
      r = Math::hypot(sphi1, cphi1);
      sphi1 /= r; cphi1 /= r;
      r = Math::hypot(sphi2, cphi2);
      sphi2 /= r; cphi2 /= r;
    }
    // bool polar = (cphi1 == 0);
    cphi1 = max(epsx, cphi1);   // Avoid singularities at poles
    cphi2 = max(epsx, cphi2);
    // Determine hemisphere of tangent latitude
    _sign = sphi1 + sphi2 >= 0 ? 1 : -1;
    // Internally work with tangent latitude positive
    sphi1 *= _sign; sphi2 *= _sign;
    if (sphi1 > sphi2) {
      swap(sphi1, sphi2); swap(cphi1, cphi2); // Make phi1 < phi2
    }
    real
      tphi1 = sphi1/cphi1, tphi2 = sphi2/cphi2;

    // q = (1-e^2)*(sphi/(1-e^2*sphi^2) - atanhee(sphi))
    // qp = q(pi/2) = (1 + (1-e^2)*atanhee(1))
    // atanhee(x) = atanh(e*x)/e
    // q = sxi * qp
    // dq/dphi = 2*(1-e^2)*cphi/(1-e^2*sphi^2)^2
    //
    // n = (m1^2-m2^2)/(q2-q1) -> sin(phi0) for phi1, phi2 -> phi0
    // C = m1^2 + n*q1 = (m1^2*q2-m2^2*q1)/(q2-q1)
    // let
    //   rho(pi/2)/rho(-pi/2) = (1-s)/(1+s)
    //   s = n*qp/C
    //     = qp * (m1^2-m2^2)/(m1^2*q2-m2^2*q1)
    //     = qp * (scbet2^2 - scbet1^2)/(scbet2^2*q2 - scbet1^2*q1)
    //     = (scbet2^2 - scbet1^2)/(scbet2^2*sxi2 - scbet1^2*sxi1)
    //     = (tbet2^2 - tbet1^2)/(scbet2^2*sxi2 - scbet1^2*sxi1)
    // 1-s = -((1-sxi2)*scbet2^2 - (1-sxi1)*scbet1^2)/
    //         (scbet2^2*sxi2 - scbet1^2*sxi1)
    //
    // Define phi0 to give same value of s, i.e.,
    //  s = sphi0 * qp / (m0^2 + sphi0*q0)
    //    = sphi0 * scbet0^2 / (1/qp + sphi0 * scbet0^2 * sxi0)

    real tphi0, n;
    if (tphi1 == tphi2) {
      tphi0 = tphi1;
      n = tphi1/hyp(tphi1);
    } else {
      real
        tbet1 = _fm * tphi1, scbet12 = 1 + sq(tbet1),
        tbet2 = _fm * tphi2, scbet22 = 1 + sq(tbet2),
        txi1 = txif(tphi1), sxi1 = txi1/hyp(txi1),
        txi2 = txif(tphi2), sxi2 = txi2/hyp(txi2),
        snum = (sq(tbet2) - sq(tbet1))/(tphi2 - tphi1),
        sden = (scbet22*sxi2 - scbet12*sxi1)/(tphi2 - tphi1),
        tphi0a = tphi1, tphi0b = tphi2;
      n = (1/scbet12 - 1/scbet22)/(_qp * (sxi2 - sxi1));
      for (unsigned i = 0; i < 100; ++i) {
        real
          tphi0c = (tphi0a + tphi0b)/2,
          sphi0c = tphi0c/hyp(tphi0c),
          scbet0c2 = 1 + sq(_fm * tphi0c),
          txi0c = txif(tphi0c),
          sxi0c = txi0c/hyp(txi0c),
          num = sphi0c * scbet0c2,
          den = 1/_qp + num * sxi0c,
          diff = num * sden - den * snum;
        if (diff > 0) {
          if (tphi0b == tphi0c)
            break;
          else
            tphi0b = tphi0c;
        } else {
          if (tphi0a == tphi0c)
            break;
          else
            tphi0a = tphi0c;
        }
      }
      tphi0 = (tphi0a + tphi0b)/2;
    }
    real txi0 = txif(tphi0);
    _q0 = _qp * txi0/hyp(txi0);
    _n0 = tphi0/hyp(tphi0);
    _m02 = 1/(1 + sq(_fm * tphi0));
    _C0 = _m02 + _n0 * _q0;
    _rho0 = _a * sqrt(_m02) / _n0;
    _k2 = n / _n0;
    _k0 = sqrt(_k2);
    _lat0 = _sign * atan(tphi0)/Constants::degree();
  }

  const AlbersEqualArea
  AlbersEqualArea::CylindricalEqualArea(Constants::WGS84_a(),
                                        Constants::WGS84_r(),
                                        real(0), real(1));

  Math::real AlbersEqualArea::txif(real tphi) const throw() {
    // sxi =
    // ( sphi/(1-e2*sphi^2) + atanhee(sphi) ) /
    // ( 1/(1-e2) + atanhee(1) )
    //
    // txi =
    // ( sphi/(1-e2*sphi^2) + atanhee(sphi) ) /
    // sqrt( ( (1+e2*sphi)*(1-sphi)/( (1-e2*sphi^2) * (1-e2) ) +
    //         atanhee((1-sphi)/(1-e2*sphi)) ) *
    //       ( (1-e2*sphi)*(1+sphi)/( (1-e2*sphi^2) * (1-e2) ) +
    //         atanhee((1+sphi)/(1+e2*sphi)) ) )
    //
    // subst 1-sphi = cphi^2/(1+sphi)
    int s = tphi < 0 ? -1 : 1;  // Enforce odd parity
    tphi *= s;
    real
      cphi2 = 1 / (1 + sq(tphi)),
      sphi = tphi * sqrt(cphi2),
      es1 = _e2 * sphi,
      es2m1 = 1 - es1 * sphi,
      sp1 = 1 + sphi,
      es1m1 = (1 - es1) * sp1,
      es2m1a = _e2m * es2m1,
      es1p1 = sp1 / (1 + es1);
    return s * ( sphi / es2m1 + atanhee(sphi) ) /
      sqrt( ( cphi2 / (es1p1 * es2m1a) + atanhee(cphi2 / es1m1) ) *
            ( es1m1 / es2m1a + atanhee(es1p1) ) );
  }

  Math::real AlbersEqualArea::tphif(real txi) const throw() {
    real
      tphi = txi,
      stol = tol * max(real(1), abs(txi));
    // CHECK: min iterations = 1, max iterations = 2; mean = 1.99
    for (int i = 0; i < numit; ++i) {
      // d(txi)/d(tphi) = (scxi/scphi)^3 * 2*(1-e^2)/(qp*(1-e^2*sphi^2)^2)
      real
        txia = txif(tphi),
        tphi2 = sq(tphi),
        scphi2 = 1 + tphi2,
        scterm = scphi2/(1 + sq(txia)),
        dtphi = (txi - txia) * scterm * sqrt(scterm) *
        _qx * sq(1 - _e2 * tphi2 / scphi2);
      tphi += dtphi;
      if (abs(dtphi) < stol)
        break;
    }
    return tphi;
  }

  void AlbersEqualArea::Forward(real lon0, real lat, real lon,
                                real& x, real& y, real& gamma, real& k)
    const throw() {
    if (lon - lon0 >= 180)
      lon -= lon0 + 360;
    else if (lon - lon0 < -180)
      lon -= lon0 - 360;
    else
      lon -= lon0;
    lat *= _sign;
    // From Snyder, we have
    real
      lam = lon * Math::degree(),
      phi = lat * Math::degree(),
      sphi = sin(phi), cphi = abs(lat) != 90 ? cos(phi) : epsx,
      tphi = sphi/cphi, tbet = _fm * tphi, scbet = hyp(tbet),
      txi = txif(tphi), q = _qp * txi/hyp(txi),
      rho = _a * sqrt(_C0 - _n0 * q)/_n0,
      theta = _k2 * _n0 * lam;
    x = (rho * sin(theta)) / _k0;
    y = (_rho0 - rho * cos(theta)) / _k0;
    k = _k0 * rho * _n0 * scbet / _a;
    y *= _sign;
    gamma = _sign * theta / Math::degree();
  }

  void AlbersEqualArea::Reverse(real lon0, real x, real y,
                                real& lat, real& lon,
                                real& gamma, real& k)
    const throw() {
    y *= _sign;
    real
      rho = Math::hypot(_k0 * x, _rho0 - _k0 * y),
      theta = atan2(_k0 * x, _rho0 - _k0 * y),
      sxi = (_C0 - sq(rho * _n0 / _a))/ (_n0 * _qp),
      txi = sxi/sqrt(1 - sq(sxi)),
      tphi = tphif(txi),
      phi = _sign * atan(tphi),
      scbet = hyp(_fm * tphi),
      lam = theta / (_k2 * _n0);
    gamma = _sign * theta / Math::degree();
    lat = phi / Math::degree();
    lon = lam / Math::degree();
    k = _k0 * rho * _n0 * scbet / _a;
  }

  void AlbersEqualArea::SetScale(real lat, real k) {
    if (!(k > 0))
      throw GeographicErr("Scale is not positive");
    if (!(abs(lat) < 90))
      throw GeographicErr("Latitude for SetScale not in (-90, 90)");
    real x, y, gamma, kold;
    Forward(0, lat, 0, x, y, gamma, kold);
    k /= kold;
    _k0 *= k;
    _k2 = sq(_k0);
  }

} // namespace GeographicLib

/**
 * \file LambertConformalConic.cpp
 * \brief Implementation for GeographicLib::LambertConformalConic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/LambertConformalConic.hpp"

#define GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_LAMBERTCONFORMALCONIC_HPP)

#include <iostream>
#include <iomanip>

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
    , _fm(1 - _f)
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
      phi = stdlat * Math::degree(),
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
    , _fm(1 - _f)
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
      phi1 = stdlat1 * Math::degree(),
      phi2 = stdlat2 * Math::degree();
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
    , _fm(1 - _f)
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
    // Snyder: 15-8: n = (log(m1) - log(m2))/(log(t1)-log(t2))
    //
    // m = cos(beta) = 1/sec(beta) = 1/sqrt(1+tan(beta)^2)
    // beta = parametric lat, tan(beta) = (1-f)*tan(phi)
    //
    // t = tan(pi/4-chi/2) = 1/(sec(chi) + tan(chi)) = sec(chi) - tan(chi)
    // log(t) = -asinh(tan(chi))
    // chi = conformal lat
    // tan(chi) = tan(phi)*cosh(xi) - sinh(xi)*sec(phi)
    // xi = eatanhe(sin(phi)), eatanhe(x) = e * atanh(e*x)
    //
    // n = (log(sec(beta2))-log(sec(beta1)))/(asinh(tan(chi2))-asinh(tan(chi1)))
    //
    //
    // Let log(sec(beta)) = b(tphi), asinh(tan(chi)) = c(tphi)
    // Then n = Db(tphi2, tphi1)/Dc(tphi2, tphi1)
    // In limit tphi2 -> tphi1, n -> sphi1
    //
    bool polar = (cphi1 == 0);
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
    _tphi0 = (sphi1 + sphi2)/(cphi1 + cphi2);
    real
      tbet1 = _fm * tphi1, scbet1 = hyp(tbet1),
      tbet2 = _fm * tphi2, scbet2 = hyp(tbet2);
    real
      scphi1 = 1/cphi1,
      xi1 = eatanhe(sphi1), shxi1 = sinh(xi1), chxi1 = hyp(shxi1),
      tchi1 = chxi1 * tphi1 - shxi1 * scphi1, scchi1 = hyp(tchi1),
      scphi2 = 1/cphi2,
      xi2 = eatanhe(sphi2), shxi2 = sinh(xi2), chxi2 = hyp(shxi2),
      tchi2 = chxi2 * tphi2 - shxi2 * scphi2, scchi2 = hyp(tchi2),
      dshxi = ( Dsinh(xi2, xi1, shxi2, shxi1, chxi2, chxi1) *
                Deatanhe(sphi2, sphi1) * Dsn(tphi2, tphi1, sphi2, sphi1) ),
      psi1 = Math::asinh(tchi1);
    if (tphi2 - tphi1 != 0) {
      real num, numcheck;

      // Db(tphi2, tphi1)
      num = Dlog(scbet2, scbet1) * Dhyp(tbet2, tbet1, scbet2, scbet1)
        * _fm;
      numcheck = (log(scbet2) - log(scbet1))/(tphi2 - tphi1);

      real den, dencheck;

      den = Dasinh(tchi2, tchi1, scchi2, scchi1)
        * // Dchi(tphi2, tphi1)
        ( ( (chxi1  + chxi2 )/2 -
            (shxi1  + shxi2 )/2 * (tphi1 + tphi2)/(scphi1 + scphi2) ) +
          ( (tphi1  + tphi2 )/2 * (shxi1 + shxi2)/(chxi1  + chxi2 ) -
            (scphi1 + scphi2)/2 ) * dshxi );
      dencheck = (Math::asinh(tchi2) - Math::asinh(tchi1))/(tphi2 - tphi1);

      _n = numcheck/dencheck;
      _n = num/den;
      _n = max(-real(1), min(real(1), _n));
      // _n = sin(phi0), _nc = cos(phi0)
      // compute _nc biasing the result towards the mean value
      real t = (1 - _n) * (1 + _n);
      _nc = 1/hyp(_tphi0);
      _nc += (t - sq(_nc))/(sqrt(t) + _nc);
      _nc = max(epsx, _nc);
      _tphi0 = _n / _nc;
    } else {
      _nc = 1/hyp(_tphi0);
      _n = _tphi0 * _nc;
      if (polar)
        _nc = 0;
    }

    _tbet0 = _fm * _tphi0; _scbet0 = hyp(_tbet0);
    _scphi0 = hyp(_tphi0); _sphi0 = _tphi0/_scphi0;
    _xi0 = eatanhe(_sphi0); _shxi0 = sinh(_xi0); _chxi0 = hyp(_shxi0);
    _tchi0 = _chxi0 * _tphi0 - _shxi0 * _scphi0; _scchi0 = hyp(_tchi0);
    _psi0 = Math::asinh(_tchi0);

    _lat0 = atan(_sign * _tphi0) / Math::degree();
    _lt0 = - _psi0; // Snyder's log(t0)
    _t0n = exp(- _n * _psi0);      // Snyder's t0^n
    _t0nm1 = Math::expm1(- _n * _psi0);      // Snyder's t0^n - 1
    // a * k1 * m1/t1^n = a * k1 * m2/t2^n = a * k1 * n * (Snyder's F)
    // = a * k1 / (scbet1 * exp(-n * psi1))
    _scale = _a * k1 /
      (scbet1 * (2 * _n <= 1 ? exp(- _n * psi1) :
                 // exp((1-n)* psi1) * exp(-psi1)
                 // with (1-n) = nc^2/(1+n) and exp(-psi1) = 1/(tchi1+scchi1)
                 exp( (sq(_nc)/(1 + _n)) * psi1 ) / (tchi1 + scchi1)));
    // Scale at phi0 = k0 = k1 * (scbet0*exp(-n*psi0))/(scbet1*exp(-n*psi1))
    //                    = k1 * scbet0/scbet1 * exp(n * (psi1 - psi0))
    // psi1 - psi0 = Dasinh(tchi1, tchi0) * (tchi1 - tchi0)
    _k0 = k1 * (_scbet0/scbet1) *
      exp(_n * Dasinh(tchi1, _tchi0, scchi1, _scchi0) * (tchi1 - _tchi0));
    _nrho0 = _a * _k0 / _scbet0;
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
      lam = lon * Math::degree(),
      phi = lat * Math::degree(),
      sphi = sin(phi), cphi = abs(lat) != 90 ? cos(phi) : epsx,
      tphi = sphi/cphi, tbet = _fm * tphi, scbet = hyp(tbet),
      scphi = 1/cphi,
      xi = eatanhe(sphi), shxi = sinh(xi), chxi = hyp(shxi),
      tchi = chxi * tphi - shxi * scphi, scchi = hyp(tchi),
      psi = Math::asinh(tchi),
      // m = 1/scbet, lt = - psi, tn = exp(_n * lt),
      theta = _n * lam, stheta = sin(theta), ctheta = cos(theta),
      dpsi = Dasinh(tchi, _tchi0, scchi, _scchi0) * (tchi - _tchi0),
      drho = - _scale * Dexp(-_n * psi, -_n * _psi0) * dpsi;
    x = (_nrho0 + _n * drho) * (_n != 0 ? stheta / _n : lam);
    y = _nrho0 * (_n != 0 ? sq(stheta)/((1 + ctheta) * _n) : 0) - drho * ctheta;
    k = _k0 * (scbet/_scbet0) * exp( -_n * dpsi );
    y *= _sign;
    gamma = theta * _sign;
  }

  void LambertConformalConic::Reverse(real lon0, real x, real y,
                                      real& lat, real& lon,
                                      real& gamma, real& k)
    const throw() {
    y *= _sign;
    real
      nx = _n * x, ny = _n * y, y1 = _nrho0 - ny,
      drho = (x*nx - 2*y*_nrho0 + y*ny) / ( Math::hypot(nx, y1) + _nrho0 ),
      dpsi = - Dlog(_t0n + _n * drho/_scale, _t0n) * drho / _scale,
      psi = _psi0 + dpsi, tchi = sinh(psi), scchi = hyp(tchi),
      dtchi = Dsinh(psi, _psi0, tchi, _tchi0, scchi, _scchi0) * dpsi,
      lam = _n != 0 ? atan2( nx, y1 ) / _n : x / y1;
    tchi = _tchi0 + dtchi;      // Update tchi using divided difference

    std::cout << std::fixed << std::setprecision(14);
    std::cout << _nrho0/_n << " " << drho << " " << _nrho0/_n + drho << "\n";
    std::cout << _psi0 << " " << dpsi << " " << psi << "\n";
    std::cout << _tchi0 << " " << dtchi << " " << tchi << "\n";

    // Use Newton's method to solve for tphi
    real
      tphi = tchi,
      stol = tol * max(real(1), abs(tchi));
    // min iterations = 1, max iterations = 2; mean = 1.99
    for (int i = 0; i < numit; ++i) {
      real
        scphi = hyp(tphi),
        shxi = sinh( eatanhe( tphi / scphi ) ),
        tchia = hyp(shxi) * tphi - shxi * scphi,
        dtphi = (tchi - tchia) * (1 + _e2m * sq(tphi)) /
        ( _e2m * scphi * hyp(tchia) );
      tphi += dtphi;
      if (abs(dtphi) < stol)
        break;
    }
    std::cout << _tphi0 << " " << tphi-_tphi0 << " " << tphi << "\n";
    double
      phi = _sign * atan(tphi),
      tbet = _fm * tphi, scbet = hyp(tbet);
    lat = phi / Math::degree();
    lon = lam / Math::degree();
    gamma = _sign * _n * lon;
    // Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    if (lon + lon0 >= 180)
      lon += lon0 - 360;
    else if (lon + lon0 < -180)
      lon += lon0 + 360;
    else
      lon += lon0;
    k = _k0 * (scbet/_scbet0) * exp( -_n * dpsi );
  }

  void LambertConformalConic::SetScale(real lat, real k) {
    if (!(k > 0))
      throw GeographicErr("Scale is not positive");
    if (!(abs(lat) <= 90))
      throw GeographicErr("Latitude for SetScale not in [-90, 90]");
    real x, y, gamma, kold;
    Forward(0, lat, 0, x, y, gamma, kold);
    k /= kold;
    _scale *= k;
    _k0 *= k;
  }

} // namespace GeographicLib

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
    // m = cos(bet) = 1/sec(bet) = 1/sqrt(1+tan(bet)^2)
    // bet = parametric lat, tan(bet) = (1-f)*tan(phi)
    //
    // t = tan(pi/4-chi/2) = 1/(sec(chi) + tan(chi)) = sec(chi) - tan(chi)
    // log(t) = -asinh(tan(chi)) = -psi
    // chi = conformal lat
    // tan(chi) = tan(phi)*cosh(xi) - sinh(xi)*sec(phi)
    // xi = eatanhe(sin(phi)), eatanhe(x) = e * atanh(e*x)
    //
    // n = (log(sec(bet2))-log(sec(bet1)))/(asinh(tan(chi2))-asinh(tan(chi1)))
    //
    // Let log(sec(bet)) = b(tphi), asinh(tan(chi)) = c(tphi)
    // Then n = Db(tphi2, tphi1)/Dc(tphi2, tphi1)
    // In limit tphi2 -> tphi1, n -> sphi1
    //
    {
      real r;
      r = Math::hypot(sphi1, cphi1);
      sphi1 /= r; cphi1 /= r;
      r = Math::hypot(sphi2, cphi2);
      sphi2 /= r; cphi2 /= r;
    }
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
      /*
      dshxi = ( Dsinh(xi2, xi1, shxi2, shxi1, chxi2, chxi1) *
                Deatanhe(sphi2, sphi1) * Dsn(tphi2, tphi1, sphi2, sphi1) ),
      */
      psi1 = Math::asinh(tchi1);
    if (tphi2 - tphi1 != 0) {
      real num, numcheck;

      // Db(tphi2, tphi1)
      num = Dlog(scbet2, scbet1) * Dhyp(tbet2, tbet1, scbet2, scbet1)
        * _fm;
      numcheck = (log(scbet2) - log(scbet1))/(tphi2 - tphi1);

      real den, dencheck;

      den =  Dasinh(tphi2, tphi1, scphi2, scphi1) - Deatanhe(sphi2, sphi1) * Dsn(tphi2, tphi1, sphi2, sphi1);
      //cout << den << "\n";

      dencheck = (Math::asinh(tchi2) - Math::asinh(tchi1))/(tphi2 - tphi1);
      _n = numcheck/dencheck;
      _n = num/den;
      _n = max(-real(1), min(real(1), _n));
      {
        // scbet - scchi
        real s1 = (tphi1 * (2 * shxi1 * chxi1 * scphi1 - _e2 * tphi1) -
                   sq(shxi1) * (1 + 2 * sq(tphi1)));
        real s2 = (tphi2 * (2 * shxi2 * chxi2 * scphi2 - _e2 * tphi2) -
                   sq(shxi2) * (1 + 2 * sq(tphi2)));
        real t1 = tchi1 < 0 ? scbet1 - tchi1 : (s1 + 1)/(scbet1 + tchi1);
        real t2 = tchi2 < 0 ? scbet2 - tchi2 : (s2 + 1)/(scbet2 + tchi2);
        s1 /= scbet1 + scchi1;
        s2 /= scbet2 + scchi2;
        real a2 = -(s2 + t2)/(2*scbet2);
        real a1 = -(s1 + t1)/(2*scbet1);
        real
          dtchi = den/Dasinh(tchi2, tchi1, scchi2, scchi1),
          tbm = ((tbet1 > 0 ? 1/(scbet1+tbet1) : scbet1 - tbet1)+
                 (tbet2 > 0 ? 1/(scbet2+tbet2) : scbet2 - tbet2))/
          (scbet1+scbet2),
          dbet = _e2 * (1/(scbet2+_fm*scphi2)+1/(scbet1+_fm*scphi1))/_fm,
          xi0 = eatanhe(real(1)), shxi0 = sinh(xi0), chxi0 = hyp(shxi0),
          dxi = Deatanhe(sphi1, sphi2) * Dsn(tphi2, tphi1, sphi2, sphi1),
          dxi1 = 2*sinh(Deatanhe(real(1), sphi1)/(scphi1*(tphi1+scphi1))/2),
          dxi2 = 2*sinh(Deatanhe(real(1), sphi2)/(scphi2*(tphi2+scphi2))/2),
          chxi01 = cosh((xi0+xi1)/2), shxi01 = sinh((xi0+xi1)/2),
          chxi02 = cosh((xi0+xi2)/2), shxi02 = sinh((xi0+xi2)/2),
          dshxi1 = dxi1*chxi01, dchxi1 = dxi1*shxi01,
          dshxi2 = dxi2*chxi02, dchxi2 = dxi2*shxi02,
          dxi01 = Deatanhe(real(1), sphi1)/(scphi1*(tphi1+scphi1)),
          dxi02 = Deatanhe(real(1), sphi2)/(scphi2*(tphi2+scphi2)),
          dshxi01 = Dsinh(xi0, xi1, shxi0, shxi1, chxi0, chxi1) * dxi01,
          dshxi02 = Dsinh(xi0, xi2, shxi0, shxi2, chxi0, chxi2) * dxi02,
          /*
            dchxi01 = Dcosh(xi0, xi1, shxi0, shxi1, chxi0, chxi1) * dxi01,
            dchxi02 = Dcosh(xi0, xi2, shxi0, shxi2, chxi0, chxi2) * dxi02,
          */
          dchxi01 = Dhyp(shxi0, shxi1, chxi0, chxi1) * dshxi01,
          dchxi02 = Dhyp(shxi0, shxi2, chxi0, chxi2) * dshxi02,
          /* mu12 = (- scphi1 * dchxi1 + tphi1 * dshxi1
             - scphi2 * dchxi2 + tphi2 * dshxi2), */
          mu12 = (- scphi1 * dchxi01 + tphi1 * dshxi01
                  - scphi2 * dchxi02 + tphi2 * dshxi02),
          ddel = (Dhyp(tphi1, tphi2, scphi1, scphi2) * (dshxi1+dshxi2)/2
                  - (dchxi1+dchxi2)/2
                  - ((scphi1+scphi2)/2  -
                     (tphi1+tphi2)/2 * Dhyp(shxi1,shxi2,chxi1,chxi2))
                  *Dsinh(xi1,xi2,shxi1,shxi2,chxi1,chxi2)*dxi),
          dchi = (mu12-ddel*(scphi2+scphi1))/dtchi,
          tam = dtchi*(dchi-dbet)/(scchi1+scchi2),
          nca = Dlog1p(a2, a1)*
          ( ( (tchi2 >= 0 ? scchi2 + tchi2 : 1/(scchi2 - tchi2)) +
              (tchi1 >= 0 ? scchi1 + tchi1 : 1/(scchi1 - tchi1)) ) /
            (4 * scbet1 * scbet2) )
          * _fm * (tbm - tam) /
          den;
        _nc = sqrt(nca * (1 + _n));
      }
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
      lam = _n != 0 ? atan2( nx, y1 ) / _n : x / y1;
    real tchi;
    if (2 * _n <= 1) {
      real
	psi = _psi0 + dpsi, tchia = sinh(psi), scchi = hyp(tchia),
	dtchi = Dsinh(psi, _psi0, tchia, _tchi0, scchi, _scchi0) * dpsi;
      tchi = _tchi0 + dtchi;	// Update tchi using divided difference
    } else {
      // tchi = sinh(-1/n * log(tn))
      // = sinh((1-1/n) * log(tn) - log(tn))
      // = + sinh((1-1/n) * log(tn)) * cosh(log(tn))
      //   - cosh((1-1/n) * log(tn)) * sinh(log(tn))
      // (1-1/n) = - nc^2/(n*(1+n))
      // cosh(log(tn)) = (tn + 1/tn)/2; sinh(log(tn)) = (tn - 1/tn)/2
      real
	tn = _t0n + _n * drho/_scale,
	sh = sinh( -sq(_nc)/(_n * (1 + _n)) * log(tn) );
      tchi = sh * (tn + 1/tn)/2 - hyp(sh) * (tn - 1/tn)/2;
    }

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

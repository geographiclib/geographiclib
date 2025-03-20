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
#include <iomanip>

namespace GeographicLib {

  using namespace std;

  TriaxialLine::hfun::hfun(bool distp, real kap, real kapp, real eps, real mu,
                           const Triaxial& t)
    : _kap(kap)
    , _kapp(kapp)
    , _eps(eps)
    , _mu(mu)
    , _sqrtkap(sqrt(_kap))
    , _sqrtkapp(sqrt(_kapp))
    , _distp(distp)
    // If oblpro extend special treatment of oblate/prolate cases to mu != 0.
    , _oblpro(t._oblpro)
    , _merid(t._merid)
    , _invp(false)
    , _biaxr(_kap == 0 && (_mu == 0 || (_oblpro && _mu > 0)))
    , _biaxl(_kapp == 0 && (_mu == 0 || (_oblpro && _mu < 0)))
    , _umb(_kap != 0 && _kapp != 0 && _mu == 0)
  {
    (void) _merid;
    // mu in [-kap, kapp], eps in (-inf, 1/kap)
    if (!_distp) {
      if (_biaxr) {
        // oblate/prolate rotating coordinate
        // _kapp == 1, mu < 0 not allowed
        _tx = false;
        // f multiplied by sqrt(mu)
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real tht) -> real
                          // This is a trivial case f' = 1
                          { return fthtoblp(tht, eps, mu); },
                          Math::pi()/2, false);
      } else if (_biaxl) {
        // oblate/prolate librating coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < t._ellipthresh;
        _tx = false;
        // f multiplied by sqrt(-mu)
        /*
          if (_tx) {
          _ell = EllipticFunction(1 + _mu, 0, -_mu, 1);
          _fun = TrigfunExt(
          [kap = _kap, kapp = _kapp,
          eps = _eps, mu = _mu, ell = _ell]
          (real v) -> real
          { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
          return dfvoblp(dn, eps, mu); },
          _ell.K(), false, 1);
          } else
        */
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return dfpsioblp(sin(psi), cos(psi), eps, mu); },
                          Math::pi()/2, false);
      } else if (_mu > 0) {
        _tx = _mu / (_kap + _mu) < t._ellipthresh;
        if (_tx) {
          _ell = EllipticFunction(_kap / (_kap + _mu), 0, _mu / (_kap + _mu), 1);
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp,
                             eps = _eps, mu = _mu, ell = _ell]
                            (real u) -> real
                            { real sn, cn, dn; (void) ell.am(u, sn, cn, dn);
                              return fup(cn, kap, kapp, eps, mu); },
                            _ell.K(), false);
        } else
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                            (real tht) -> real
                            { return fthtp(cos(tht), kap, kapp, eps, mu); },
                            Math::pi()/2, false);
      } else if (_mu < 0) {
        _tx = -_mu / _kap < t._ellipthresh;
        if (_tx) {
          _ell = EllipticFunction((_kap + _mu) / _kap, 0, -_mu / _kap, 1);
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp,
                             eps = _eps, mu = _mu, ell = _ell]
                            (real v) -> real
                            { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                              return fvp(dn, kap, kapp, eps, mu); },
                            _ell.K(), false);
        } else
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                            (real psi) -> real
                            { return fpsip(sin(psi), cos(psi),
                                           kap, kapp, eps, mu); },
                            Math::pi()/2, false);
      } else if (_umb) {
        _tx = _kapp < t._ellipthresh;
        // f multiplied by sqrt(kap*kapp)
        // Include scale = 1 in TrigfunExt constructor because this function gets
        // added to u.
        if (_tx) {
          _ell = EllipticFunction(_kap, 0, _kapp, 1);
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp,
                             eps = _eps, ell = _ell]
                            (real v) -> real
                            { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                              return dfvp(cn, dn, kap, kapp, eps); },
                            2 * _ell.K(), true, 1);
        } else
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp, eps = _eps]
                            (real tht) -> real
                            { return dfp(cos(tht), kap, kapp, eps); },
                            Math::pi(), true, 1);
      } else {
        // _mu == NaN
        _tx = false;
      }
    } else {
      if (_biaxr) {
        // oblate/prolate symmetry coordinate
        // _kapp == 1, mu < 0 not allowed
        _tx = false;
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real tht) -> real
                          // degenerate f' = 0
                          { return gthtoblp(tht, eps, mu); },
                          Math::pi()/2, false);
      } else if (_biaxl) {
        // oblate/prolate non-symmetry coordinate
        // _kap == 1, mu > 0 not allowed
        // DON'T USE tx: _tx = _mu < 0 &&  -_mu < t._ellipthresh;
        _tx = false;
        /*
          if (_tx) {
          _ell = EllipticFunction(1 + _mu, 0, -_mu, 1);
          // Never completed this
          } else
        */
        _fun = TrigfunExt(
                          [eps = _eps, mu = _mu]
                          (real psi) -> real
                          { return gpsioblp(sin(psi), cos(psi), eps, mu); },
                          Math::pi()/2, false);
      } else if (_mu > 0) {
        _tx = _mu / (_kap + _mu) < t._ellipthresh;
        if (_tx) {
          _ell = EllipticFunction(_kap / (_kap + _mu), 0, _mu / (_kap + _mu), 1);
          _fun =TrigfunExt(
                           [kap = _kap, kapp = _kapp,
                            eps = _eps, mu = _mu, ell = _ell]
                           (real u) -> real
                           { real sn, cn, dn; (void) ell.am(u, sn, cn, dn);
                             return gup(cn, dn, kap, kapp, eps, mu); },
                           _ell.K());
        } else
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                            (real tht) -> real
                            { return gthtp(cos(tht), kap, kapp, eps, mu); },
                            Math::pi()/2);
      } else if (_mu < 0) {
        _tx = -_mu / _kap < t._ellipthresh;
        if (_tx) {
          _ell = EllipticFunction((_kap + _mu) / _kap, 0, -_mu / _kap, 1);
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp,
                             eps = _eps, mu = _mu, ell = _ell]
                            (real v) -> real
                            { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                              return gvp(cn, dn, kap, kapp, eps, mu); },
                            _ell.K());
        } else
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp, eps = _eps, mu = _mu]
                            (real psi) -> real
                            { return gpsip(sin(psi), cos(psi),
                                           kap, kapp, eps, mu); },
                            Math::pi()/2);
      } else if (_umb) {
        _tx = _kapp < t._ellipthresh;
        if (_tx) {
          _ell = EllipticFunction(_kap, 0, _kapp, 1);
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp, eps = _eps, ell = _ell]
                            (real v) -> real
                            { real sn, cn, dn; (void) ell.am(v, sn, cn, dn);
                              return g0vp(cn, kap, kapp, eps); },
                            2*_ell.K(), true);
        } else
          _fun = TrigfunExt(
                            [kap = _kap, kapp = _kapp, eps = _eps]
                            (real tht) -> real
                            { return g0p(cos(tht), kap, kapp, eps); },
                            Math::pi(), true);
      } else {
        // _mu == NaN
        _tx = false;
      }
    }
    // N.B. _max < 0 for _umb && eps < 0
    _max = !_distp ?
      ( _umb ? _fun(_tx ? _ell.K() : Math::pi()/2) :
        (_biaxl ? Math::pi()/2 + sqrt(-_mu) * _fun.Max() : _fun.Max()) ) :
      ( _umb ? _fun(_tx ? _ell.K() : Math::pi()/2) : _fun.Max() );
  }

  Math::real TriaxialLine::hfun::operator()(real u) const {
    if (!_distp) {
      if (_biaxl) {
        // This is sqrt(-mu) * f(u)
        return modang(u, sqrt(-_mu)) - sqrt(-_mu) * _fun(u);
      } else if (_umb) {
        // This is sqrt(kap * kapp) * f(u)
        real phi = gd(u, _sqrtkapp);
        return u - _fun(_tx ? _ell.F(phi) : phi);
      } else
        return _fun(u);
    } else {
      if (_umb) {
        real phi = gd(u, _sqrtkapp);
        // cout << "BB " << u << " " << phi << " " << _fun(phi) << "\n";
        return _fun(_tx ? _ell.F(phi) : phi);
      } else
        return _fun(u);
    }
  }

  // THIS ISN"T USED ?
  // Should implement an Angle equivalant of _ell.F(phi)
  Angle TriaxialLine::hfun::operator()(const Angle& ang) const {
    if (_distp) return Angle::NaN();
    if (_biaxr)
      return ang;
    else if (_biaxl && _mu == 0)
      return ang.modang(sqrt(-_mu));
    real u = ang.radians();
    if (_biaxl)
      // This is sqrt(-mu) * f(u)
      return ang.modang(sqrt(-_mu)) - ang::radians(sqrt(-_mu) * _fun(u));
    else if (_umb) {
      // This is sqrt(kap * kapp) * f(u)
      real phi = gd(u, _sqrtkapp);
      return ang::radians(u - _fun(_tx ? _ell.F(phi) : phi));
    } else
      return ang::radians(_fun(u));
  }

  Math::real TriaxialLine::hfun::deriv(real u) const {
    if (!_distp) {
      if (_biaxl) {
        // This is sqrt(-mu) * f'(u)
        /*
          if (_tx) {
          // DLMF (22.6.1): sn^2 + cn^2 =  k^2*sn^2 + dn^2 = 1
          // dn^2 = cn^2 + k'^2*sn^2 = cn^2 - mu*sn^2
          // k'^2 = -mu
          real sn, cn, dn;
          (void) _ell.am(u, sn, cn, dn);
          return sqrt(-_mu) / Math::sq(dn) - _fun.deriv(u);
          } else
        */
        // f0(x) = atan(sqrt(-mu) * tanx(u))
        // f0'(x) = sqrt(-mu) / (cos(u)^2 - mu * sin(u)^2)
        // **HERE**
        return sqrt(-_mu) / (Math::sq(cos(u)) - _mu * Math::sq(sin(u)))
          - sqrt(-_mu) * _fun.deriv(u);
      } else if (_umb) {
        // This is sqrt(kap * kapp) * f'(u)
        real phi = gd(u, _sqrtkapp),
          t = _kapp + Math::sq(sinh(u));
        // dphi/du = _sqrtkapp * cosh(u) / t;
        // for tx w = F(phi, sqrtkap)
        //   dw/dphi = 1/sqrt(kapp + kap*cos(phi)^2)
        //           = sqrt(t)/(sqrtkapp *cosh(u))
        //   N.B. cos(phi)^2 = kapp/t
        //   dw/du = dw/dphi*dphi/du = 1/sqrt(t)
        return 1 - _fun.deriv(_tx ? _ell.F(phi) : phi) /
          ( _tx ? sqrt(t) : t / (_sqrtkapp * cosh(u)) );
      } else
        return _fun.deriv(u);
    } else {
      if (_umb) {
        real phi = gd(u, _sqrtkapp),
          t = _kapp + Math::sq(sinh(u));
        // See comments in ffun::deriv
        return _fun.deriv(_tx ? _ell.F(phi) : phi) /
          ( _tx ? sqrt(t) : t / (_sqrtkapp * cosh(u)) );
      } else
        return _fun.deriv(u);
    }
  }

  Math::real TriaxialLine::hfun::gfderiv(real u) const {
    // return g'(u)/f'(u)
    real sn = 0, cn = 0, dn = 0;
    if (_biaxr)
      return gfthtoblp(u, _mu);
    else if (_biaxl) // mu > 0 not allowed
      return gfpsioblp(sin(u), cos(u), _mu);
    else if (_umb)
      // This includes factor of sqrt(kap * kapp) because of adjustment of
      // definition of f for umbilical geodesics.
      return gf0up(u, _kap, _kapp);
    else {
      if (_tx)
        (void) _ell.am(u, sn, cn, dn);
      if (_mu > 0)
        return _tx ? gfup(cn, _kap, _mu) : gfthtp(cos(u), _kap, _mu);
      else if (_mu < 0)
        return _tx ? gfvp(dn, _kap, _mu) :
          gfpsip(sin(u), cos(u), _kap, _mu);
      else
        return Math::NaN();
    }
  }

  void TriaxialLine::hfun::ComputeInverse() {
    if (!_distp) {
      if (!_invp) {
        if (_biaxl) {
          if (_mu == 0) return; // _fun == 0 and there's an analytic inverse
          // now _mu < 0
          _countn = _countb = 0;
          // Include scale = 1 in TrigfunExt constructor because _dfinv gets
          // added to u.
          // Ars are fun, odd, sym, halfp, nmax, tol, scale
          // **HERE**
          _dfinv = Trigfun(
                           [this]
                           (real z, real u1) -> real
                           {
                             real u0 = modang(z/Slope(), 1/sqrt(-_mu));
                             return root(z, u0 + u1, &_countn, &_countb,
                                         sqrt(numeric_limits<real>::epsilon()))
                               - u0;
                           },
                           true, false, Slope() * HalfPeriod(),
                           int(ceil(real(1.5) * NCoeffs())),
                           sqrt(numeric_limits<real>::epsilon()), HalfPeriod());
        } else if (_umb) {
          _countn = _countb = 0;
          // Include scale = 1 in TrigfunExt constructor because _dfinv gets
          // added to z.
          _dfinv = Trigfun(
                           [this]
                           (real phi, real u1) -> real
                           {
                             real z = lam(phi, _sqrtkapp);
                             return root(z, z + u1, &_countn, &_countb,
                                         sqrt(numeric_limits<real>::epsilon()))
                               - z;
                           },
                           true, true, Math::pi(),
                           int(ceil(real(1.5) * NCoeffs())),
                           sqrt(numeric_limits<real>::epsilon()), 1);
        } else
          _fun.ComputeInverse();
      }
      _invp = true;
    } else {
      if (!(_invp || _biaxr || _umb)) {
        // If _umb, the inverse isn't periodic
        _fun.ComputeInverse();
        /*
          real u = 1, z = _fun(u);
          cout << "HERE " << u << " " << z << " " << _fun.inv0(z) << "\n";
          for (int i = -100; i <= 100; ++i) {
          real u = real(i)/10;
          cout << "DD " << u << " " << _fun(u) << "\n";
          }
        */
        _invp = true;
      }
    }
  }

  Math::real TriaxialLine::hfun::root(real z, real u0,
                                      int* countn, int* countb,
                                      real tol) const {
    if (!_distp) {
      if (!isfinite(z)) return z; // Deals with +/-inf and nan
      if (_biaxl) {
        // z = 1/sqrt(-mu)
        // here f(u) = atan(sqrt(-mu)*tan(u))/sqrt(-mu)-Deltaf(u)
        // let fx(u) = sqrt(-mu) * f(u); fun(u) = sqrt(-mu)*Deltaf(u)
        // z = fx(u) = modang(u, sqrt(-mu)) - fun(u)
        // fun(u) is monotonically increasing/decreasing quasilinear function
        // period pi.  Combined function is quasilinear
        // z = s * u +/- m; period = pi
        // Inverting:
        // u = z/s +/- m/s; period s*pi
        // u = fxinv(z) = modang(z/x, 1/sqrt(-mu)) + funinv(z)
        // funinv(z) periodic function of z period (1-s)*pi
        real d = Max(),
          ua = (z - d) / Slope(),
          ub = (z + d) / Slope();
        u0 = fmin(ub, fmax(ua, u0));
        return Trigfun::root(
                             [this]
                             (real u) -> pair<real, real>
                             { return pair<real, real>((*this)(u), deriv(u)); },
                             z,
                             u0, ua, ub,
                             HalfPeriod(), HalfPeriod()/Slope(), 1,
                             countn, countb, tol, Trigfun::FFUNROOT);
      } else if (_umb) {
        real d = fabs(Max())
          + 2 * numeric_limits<real>::epsilon() * fmax(real(1), fabs(z)),
          ua = z - d,
          ub = z + d;
        u0 = fmin(ub, fmax(ua, u0));
        return Trigfun::root(
                             [this]
                             (real u) -> pair<real, real>
                             { return pair<real, real>((*this)(u), deriv(u)); },
                             z,
                             u0, ua, ub,
                             Math::pi()/2, Math::pi()/2, 1, countn, countb, tol,
                             Trigfun::FFUNROOT);
      } else
        return Math::NaN();
    } else {
      // This function isn't neeed.  General inversion mechanisms in Trigfun
      // suffice.  NO, the trigfun for _umb is not invertible.
      if (!(isfinite(z) && _umb))
        return Math::NaN();       // Deals with +/-inf and nan
      // Now we're dealing with _umb.
      if (fabs(z) >= Max())
        return copysign(Triaxial::BigValue(), z);
      real ua = -Triaxial::BigValue(), ub = -ua;
      u0 = fmin(ub, fmax(ua, u0));
      // Solve z = _fun(_tx ? _ell.F(gd(u)) : gd(u)) for u
      return Trigfun::root(
                           [this]
                           (real u) -> pair<real, real>
                           { return pair<real, real>((*this)(u), deriv(u)); },
                           z,
                           u0, ua, ub,
                           Math::pi()/2, Math::pi()/2, 1, countn, countb, tol,
                           Trigfun::GFUNROOT);
    }
  }

  // Approximate inverse using _dfinv _fun.inv0
  Math::real TriaxialLine::hfun::inv0(real z) const {
    if (_distp) {
      if (!_invp) return Math::NaN();
      // For the inverse in the umbilical case, just use gd(z, _sqrtkapp) and
      // not F(gd(z, _sqrt(kapp)))
      return _umb ? z + _dfinv(gd(z, _sqrtkapp)) :
        (_biaxl ? modang(z/Slope(), 1/sqrt(-_mu)) + (_mu == 0 ? 0 : _dfinv(z)) :
         _fun.inv0(z));
    } else {
      return _invp ? _fun.inv0(z) :
        (_umb ?
         // In limit _eps -> 0
         //   g(u) = atan(_sqrtkap/_sqrtkapp * tanh(u))
         // at u = 0, dg/du = _sqrtkap/_sqrtkapp
         //    u = inf, g = atan(_sqrtkap/_sqrtkapp)
         //
         // For _eps finite
         // at u = 0, dg/du = _sqrtkap/_sqrtkapp * sqrt(1 - _eps*_kap)
         //    u = inf, g = _max
         // Note: at u = 0
         //   du/dphi = _sqrtkapp, so
         //   dg/dphi = _sqrtkap * sqrt(1 - _eps*_kap) -- OK
         //
         // Approximate g(u) for _eps finite by
         // g(u) = _max / atan(_sqrtkap/_sqrtkapp) *
         //   atan(_sqrtkap/_sqrtkapp *
         //        tanh(sqrt(1 - _eps*_kap) * atan(_sqrtkap/_sqrtkapp) / _max
         //             * u))
         // Values at +/- inf and slope at origin match.
         //
         // Solve z = g(u) gives u = u0:
         atanh(_sqrtkapp/_sqrtkap * tan(atan(_sqrtkap/_sqrtkapp) / _max * z)) /
         (sqrt(1 - _eps*_kap) * atan(_sqrtkap/_sqrtkapp) / _max)
         : Math::NaN());
    }
  }

  // Accurate inverse by direct Newton (not using _finv)
  Math::real TriaxialLine::hfun::inv1(real z, int* countn, int* countb) const {
    if (!_distp)
      return _umb ? root(z, z, countn, countb) :
        (_biaxl ? (_mu == 0 ?
                   // In this case _fun.Slope() = 0 and Slope() = 1
                   modang(z/Slope(), 1/sqrt(-_mu)) :
                   root(z, modang(z/Slope(), 1/sqrt(-_mu)), countn, countb)) :
         _fun.inv1(z, countn, countb));
    else {
      if (_biaxr) return Math::NaN();
      return _umb ? root(z, inv0(z), countn, countb) :
        _fun.inv1(z, countn, countb);
    }
  }

  // Accurate inverse correcting result from _finv
  Math::real TriaxialLine::hfun::inv2(real z, int* countn, int* countb) const {
    if (!_invp) return Math::NaN();
    if (!_distp)
      return _biaxl && _mu == 0 ? inv1(z) :
        (_umb || _biaxl ? root(z, inv0(z), countn, countb) :
         _fun.inv2(z, countn, countb));
    else
      return _umb ? root(z, inv0(z), countn, countb) :
        _fun.inv1(z, countn, countb);
  }

  Angle TriaxialLine::hfun::inv(const Angle& z, int* countn, int* countb)
    const {
    if (!_distp) return ang::NaN();
    if (_biaxr)
      return z;
    else if (_biaxl && _mu == 0)
      return z.modang(1/sqrt(-_mu));
    else
      return ang::radians(inv(z.radians(), countn, countb));
  }

  // _mu > 0 && !_tx
  Math::real TriaxialLine::hfun::fthtp(real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }
  // This is non-negative
  Math::real TriaxialLine::hfun::gthtp(real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c);
    return c2 * sqrt((1 - eps * c2) / ((kapp + c2) * (c2 + mu)) );
  }
  Math::real TriaxialLine::hfun::gfthtp(real c, real kap, real /* mu */) {
    real c2 = kap * Math::sq(c);
    return c2;
  }

  // _mu > 0 && _tx
  Math::real TriaxialLine::hfun::fup(real cn, real kap, real kapp,
                                     real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }
  // This is non-negative
  Math::real TriaxialLine::hfun::gup(real cn, real /* dn */,
                                     real kap, real kapp,
                                     real eps, real mu) {
    real c2 = kap * Math::sq(cn);
    return c2 * sqrt( (1 - eps * c2) / ((kapp + c2) * (kap + mu)) );
  }
  Math::real TriaxialLine::hfun::gfup(real cn, real kap, real /* mu */) {
    real c2 = kap * Math::sq(cn);
    return c2;
  }

  // _mu == 0 && !_tx
  Math::real TriaxialLine::hfun::dfp(real c,
                                     real kap, real kapp, real eps) {
    // function dfp = dfpf(phi, kappa, epsilon)
    // return derivative of sqrt(kap * kapp) * Delta f
    // s = sqrt(1 - kap * sin(phi)^2)
    real c2 = kap * Math::sq(c), s = sqrt(kapp + c2);
    return eps*kap * sqrt(kapp) * c / (s * (1 + sqrt(1 - eps*c2)));
  }
  Math::real TriaxialLine::hfun::g0p(real c, real kap, real kapp, real eps) {
    real c2 = kap * Math::sq(c);
    return sqrt( kap * (1 - eps * c2) / (kapp + c2) ) * c;
  }

  // _mu == 0 && _tx
  Math::real TriaxialLine::hfun::dfvp(real cn, real /* dn */,
                                      real kap, real kapp, real eps) {
    // function dfvp = dfvpf(v, kap, eps)
    // return derivative of sqrt(kap * kapp) * Delta f_v
    return eps*kap * sqrt(kapp) * cn /
      (1  + sqrt(1 - eps*kap * Math::sq(cn)));
  }
  Math::real TriaxialLine::hfun::g0vp(real cn, real kap, real /* kapp */,
                                      real eps) {
    real c2 = kap * Math::sq(cn);
    return sqrt( kap * (1 - eps * c2) ) * cn;
  }

  // _mu == 0 (_tx ignored)
  Math::real TriaxialLine::hfun::gf0up(real u, real kap, real kapp) {
    // Subst tan(phi) = sinh(u) /sqrt(kapp) in
    // kap * cos(phi)^2 gives kap*kapp/(kapp + sinh(u)^2)
    // Divide by sqrt(kap * kappp) to account of factor removed from f
    // functions.
    return sqrt(kap * kapp) / ( kapp + Math::sq(sinh(u)) );
  }

  // _mu < 0 && !_tx
  Math::real TriaxialLine::hfun::fpsip(real s, real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * c2) ) ;
  }
  // This is positive
  Math::real TriaxialLine::hfun::gpsip(real s, real c, real kap, real kapp,
                                       real eps, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return sqrt(c2 * (1 - eps * c2) / (kapp + c2)) ;
  }
  Math::real TriaxialLine::hfun::gfpsip(real s, real c, real kap, real mu) {
    real c2 = kap * Math::sq(c) - mu * Math::sq(s);
    return c2;
  }

  // _mu < 0 && _tx
  Math::real TriaxialLine::hfun::fvp(real dn, real kap, real kapp,
                                     real eps, real /* mu */) {
    real c2 = kap * Math::sq(dn);
    return sqrt( (1 - eps * c2) / ((kapp + c2) * kap) );
  }
  // This is positive
  Math::real TriaxialLine::hfun::gvp(real /* cn */, real dn,
                                     real kap, real kapp,
                                     real eps, real /* mu */) {
    real dn2 = Math::sq(dn), c2 = kap * dn2;
    return dn2 * sqrt( kap * (1 - eps * c2) / (kapp + c2) );
  }
  Math::real TriaxialLine::hfun::gfvp(real dn, real kap, real /* mu */) {
    real dn2 = Math::sq(dn), c2 = kap * dn2;
    return c2;
  }

  // oblate/prolate variants for kap = 0, kapp = 1, mu > 0
  Math::real TriaxialLine::hfun::fthtoblp(real /* tht */, real /* eps */,
                                          real /* mu */) {
    // Multiply by f functions by sqrt(abs(mu))
    // return 1 / sqrt(mu);
    return 1;
  }
  Math::real TriaxialLine::hfun::gthtoblp(real /* tht */, real /* eps */,
                                          real /* mu */) {
    return 0;
  }
  Math::real TriaxialLine::hfun::gfthtoblp(real /* tht */, real /* mu */) {
    return 0;
  }

  // oblate/prolate variants for kap = 1, kapp = 0, mu <= 0, !_tx
  Math::real TriaxialLine::hfun::dfpsioblp(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // f functions are multiplied by sqrt(abs(mu)) but don't include this
    // factor here; instead include it in operator()(). etc.  This was we can
    // still use this function in the limit mu -> 0 to determine the conjugate
    // point on a meridian.
    return eps / (1 + sqrt(1 - eps * c2));
  }
  Math::real TriaxialLine::hfun::gpsioblp(real s, real c, real eps, real mu) {
    real c2 = Math::sq(c) - mu * Math::sq(s);
    // return gpsip(s, c, 1, 0, eps, mu) but with the factor c2 canceled
    return sqrt(1 - eps * c2);
  }
  Math::real TriaxialLine::hfun::gfpsioblp(real s, real c, real mu) {
    // cos(phi)^2 = cos(psi)^2 - mu *sin(psi)^2
    // f' = sqrt(1-eps*cos(phi)^2)/cos(phi)^2
    // g' = sqrt(1-eps*cos(phi)^2)
    // g'/f' = cos(phi)^2
    // Adjust by sqrt(-mu) to accommodate this factor in dfpsioblp
    return gfpsip(s, c, 1, mu) / sqrt(-mu);
  }

#if 0
  // oblate/prolate variants for kap = 1, kapp = 0, mu <= 0, _tx
  Math::real TriaxialLine::hfun::dfvoblp(real dn, real eps, real mu) {
    real c2 = Math::sq(dn);
    // Multiply by f functions by sqrt(abs(mu))
    return sqrt(-mu) * eps * dn / (1 + sqrt(1 - eps * c2));
  }
  // This is positive
  Math::real TriaxialLine::hfun::gvoblp(real /* cn */, real dn,
                                        real eps, real /* mu */) {
    // return gvp(dn, cn, 1, 0, epd, mu) but with cancelation
    real c2 = Math::sq(dn);
    return dn * sqrt(1 - eps * c2);
  }
  Math::real TriaxialLine::hfun::gfvoblp(real dn, real mu) {
    return gfvp(dn, 1, mu) / sqrt(-mu);
  }
#endif

} // namespace GeographicLib

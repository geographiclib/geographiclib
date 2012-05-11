/**
 * \file Ellipsoid.cpp
 * \brief Implementation for GeographicLib::Ellipsoid class
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/Ellipsoid.hpp>

#define GEOGRAPHICLIB_ELLIPSOID_CPP \
  "$Id$"

RCSID_DECL(GEOGRAPHICLIB_ELLIPSOID_CPP)
RCSID_DECL(GEOGRAPHICLIB_ELLIPSOID_HPP)

namespace GeographicLib {

  using namespace std;

  Ellipsoid::Ellipsoid(real a, real f)
    : _a(a)
    , _f(f <= 1 ? f : 1/f)
    , _f1(1 - _f)
    , _f12(Math::sq(_f1))
    , _e2(_f * (2 - _f))
    , _e12(_e2 / (1 - _e2))
    , _n(_f / (2  - _f))
    , _b(_a * _f1)
    , _stol(max(_a, _b) * 0.01 * sqrt(numeric_limits<real>::epsilon()))
      // TransverseMercatorExact only handles oblate ellipsoid; if prolate,
      // interchange major and minor radii.
    , _tm(_f >= 0 ? _a : _a * _f1,
          _f >= 0 ? _f : -_f / _f1,
          real(1), false)
    , _au(_a, _f, real(0), real(1), real(0), real(1), real(1))
  {}


  Math::real Ellipsoid::ConformalLatitude(real phi) const throw()
  { return atand(_tm.taup(tand(phi))); }

  Math::real Ellipsoid::InverseConformalLatitude(real chi) const throw()
  { return atand(_tm.taupinv(tand(chi))); }
  
  Math::real Ellipsoid::ParametricLatitude(real phi) const throw()
  { return atand(_f1 * tand(phi)); }

  Math::real Ellipsoid::InverseParametricLatitude(real beta) const throw()
  { return atand(tand(beta) / _f1); }

  Math::real Ellipsoid::GeocentricLatitude(real phi) const throw()
  { return atand(_f12 * tand(phi)); }

  Math::real Ellipsoid::InverseGeocentricLatitude(real theta) const throw()
  { return atand(tand(theta) / _f12); }

  Math::real Ellipsoid::AuthalicLatitude(real phi) const throw()
  { return atand(_au.txif(tand(phi))); }

  Math::real Ellipsoid::InverseAuthalicLatitude(real xi) const throw()
  { return atand(_au.tphif(tand(xi))); }

  Math::real Ellipsoid::QuarterMeridian() const throw()
  { return (_f >= 0 ? _a : _b) * _tm._Eu.E(); }
  
  Math::real Ellipsoid::MeridianDistance(real phi) const throw() {
    if (phi == 0) return 0;
    real
      tbet = _f1 * tand(abs(phi)),
      cbet = 1 / Math::hypot(real(1), tbet),
      sbet = tbet * cbet,
      E = _f >= 0 ?
      _tm._Eu.E(cbet, sbet, sqrt(1 - _tm._Eu.m() * Math::sq(cbet))) :
      _tm._Eu.E(sbet, cbet, sqrt(1 - _tm._Eu.m() * Math::sq(sbet)));
    return (phi < 0 ? -1 : 1) * (_f >= 0 ? _a * (_tm._Eu.E() - E) : _b * E);
  }

  Math::real Ellipsoid::RectifyingLatitude(real phi) const throw() {
    return abs(phi) == 90 ? phi:
      90 * MeridianDistance(phi) / QuarterMeridian();
  }

  Math::real Ellipsoid::InverseRectifyingLatitude(real mu) const throw() {
    if (abs(mu) == 90)
      return mu;
    real
      mdist0 = mu * QuarterMeridian() / 90,
      // Include first order correction in initial guess
      phi = mu + (real(1.5) * _n * sin(2 * mu * Math::degree<real>())
                       / Math::degree<real>());
    for (int i = 0; i < numit_; ++i) {
      // Solve by Newton's method
      real err = MeridianDistance(phi) - mdist0;
      phi = phi - err / (MeridionalRadius(phi) * Math::degree<real>());
      if (abs(err) < _stol)
        break;
    }
    return phi;
  }

  Math::real Ellipsoid::IsometricLatitude(real phi) const throw()
  { return Math::asinh(_tm.taup(tand(phi))) / Math::degree<real>(); }

  Math::real Ellipsoid::InverseIsometricLatitude(real psi) const throw()
  { return atand(_tm.taupinv (sinh(psi * Math::degree<real>()))); }

  Math::real Ellipsoid::CircleRadius(real phi) const throw() {
    return abs(phi) == 90 ? 0 :
      // a * cos(beta)
      _a / Math::hypot(real(1), _f1 * tand(phi));
  }

  Math::real Ellipsoid::TransverseRadius(real phi) const throw() {
    real v = 1 - _e2 * Math::sq(sin(phi * Math::degree<real>()));
    return _a / sqrt(v);
  }

  Math::real Ellipsoid::MeridionalRadius(real phi) const throw() {
    real v = 1 - _e2 * Math::sq(sin(phi * Math::degree<real>()));
    return _a * (1 - _e2) / (v * sqrt(v));
  }
  

} // namespace GeographicLib

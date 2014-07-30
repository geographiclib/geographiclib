/**
 * \file Rhumb.cpp
 * \brief Implementation for GeographicLib::Rhumb and GeographicLib::RhumbLine
 * classes
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/Rhumb.hpp>

namespace GeographicLib {

  using namespace std;

  RhumbLine::RhumbLine(const Ellipsoid& ell, real lat1, real lon1, real azi12,
                       bool exact)
    : _ell(ell)
    , _exact(exact)
    , _lat1(lat1)
    , _lon1(Math::AngNormalize(lon1))
    , _azi12(Math::AngNormalize(azi12))
  {
    real alp12 = azi12 * Math::degree();
    _salp =     azi12  == -180 ? 0 : sin(alp12);
    _calp = abs(azi12) ==   90 ? 0 : cos(alp12);
    _mu1 = _ell.RectifyingLatitude(lat1);
    _psi1 = _ell.IsometricLatitude(lat1);
    _r1 = _ell.CircleRadius(lat1);
  }

  void RhumbLine::Position(real s12, real& lat2, real& lon2) const {
    real
      mu12 = s12 * _calp * 90 / _ell.QuarterMeridian(),
      mu2 = _mu1 + mu12;
    if (abs(mu2) <= 90) {
      if (_calp) {
        lat2 = _ell.InverseRectifyingLatitude(mu2);
        real psi12;
        if (_exact)
          psi12 = !(abs(mu12) <= tol()) ?
            _ell.IsometricLatitude(lat2) - _psi1 :
            // use dpsi/dmu = (2/pi)*Q/R
            mu12 * 4 * _ell.QuarterMeridian() /
            (Math::pi() * (_r1 + _ell.CircleRadius(lat2)));
        else
          psi12 = DRectifyingToIsometric( mu2 * Math::degree(),
                                          _mu1 * Math::degree()) * mu12;
        lon2 = _salp * psi12 / _calp;
      } else {
        lat2 = _lat1;
        lon2 = _salp * s12 / (_r1 * Math::degree());
      }
      lon2 = Math::AngNormalize2(_lon1 + lon2);
    } else {
      // Reduce to the interval [-180, 180)
      mu2 = Math::AngNormalize2(mu2);
      // Deal with points on the anti-meridian
      if (abs(mu2) > 90) mu2 = Math::AngNormalize(180 - mu2);
      lat2 = _ell.InverseRectifyingLatitude(mu2);
      lon2 = Math::NaN();
    }
  }

  Math::real RhumbLine::DRectifyingToConformal(real x, real y) const {
    real d = x - y, p = x + y, s = 0;
    for (int j = tm_maxord; j; --j)
      s += j * _ell.RectifyingToConformalCoeffs()[j] * cos(j * p) * sinc(j * d);
    return 1 - 2 * s;
  }

  Math::real RhumbLine::DRectifyingToIsometric(real x, real y) const {
    real
      chix = _ell.ConformalLatitude
      (_ell.InverseRectifyingLatitude(x/Math::degree())) * Math::degree(),
      chiy = _ell.ConformalLatitude
      (_ell.InverseRectifyingLatitude(y/Math::degree())) * Math::degree();
    return Dgdinv(chix, chiy) * DRectifyingToConformal(x, y);
  }

  const Rhumb& Rhumb::WGS84() {
    static const Rhumb wgs84(Constants::WGS84_a(), Constants::WGS84_f(), false);
    return wgs84;
  }

  void Rhumb::Inverse(real lat1, real lon1, real lat2, real lon2,
                      real& s12, real& azi12) const {
    real
      lon12 = Math::AngDiff(Math::AngNormalize(lon1), Math::AngNormalize(lon2)),
      psi1 = _ell.IsometricLatitude(lat1),
      psi2 = _ell.IsometricLatitude(lat2),
      psi12 = psi2 - psi1,
      h = Math::hypot(lon12, psi12);
    azi12 = 0 - atan2(-lon12, psi12) / Math::degree();
    real dmudpsi;
    if (_exact) {
      real mu12 = _ell.RectifyingLatitude(lat2) - _ell.RectifyingLatitude(lat1);
      dmudpsi = !(abs(mu12) <= RhumbLine::tol()) ? mu12 / psi12 :
        // use dmu/dpsi = (pi/2)*R/Q
        Math::pi() * (_ell.CircleRadius(lat1) + _ell.CircleRadius(lat2)) /
        (4 * _ell.QuarterMeridian());
    } else
      dmudpsi = DIsometricToRectifying(psi2 * Math::degree(),
                                       psi1 * Math::degree());
    s12 = h * dmudpsi * _ell.QuarterMeridian() / 90;
  }

  RhumbLine Rhumb::Line(real lat1, real lon1, real azi12) const
  { return RhumbLine(_ell, lat1, lon1, azi12, _exact); }

  void Rhumb::Direct(real lat1, real lon1, real azi12, real s12,
                     real& lat2, real& lon2) const
  { Line(lat1, lon1, azi12).Position(s12, lat2, lon2); }

  Math::real Rhumb::DConformalToRectifying(real x, real y) const {
    real d = x - y, p = x + y, s = 0;
    for (int j = RhumbLine::tm_maxord; j; --j)
      s += j * _ell.ConformalToRectifyingCoeffs()[j] * cos(j * p) *
        RhumbLine::sinc(j * d);
    return 1 + 2 * s;
  }

  Math::real Rhumb::DIsometricToRectifying(real x, real y) const {
    return DConformalToRectifying(gd(x), gd(y)) * Dgd(x, y);
  }

} // namespace GeographicLib

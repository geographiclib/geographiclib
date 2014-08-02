/**
 * \file Rhumb.cpp
 * \brief Implementation for GeographicLib::Rhumb and GeographicLib::RhumbLine
 * classes
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <algorithm>
#include <GeographicLib/Rhumb.hpp>

namespace GeographicLib {

  using namespace std;

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
      dmudpsi = !(abs(mu12) <= tol()) ? mu12 / psi12 :
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
    return 1 + SinSeries(x, y, _ell.ConformalToRectifyingCoeffs(), tm_maxord);
  }

  Math::real Rhumb::SinSeries(real x, real y, const real c[], int n) {
    // N.B. n >= 0 and c[] has n+1 elements 0..n, of which c[0] is ignored.
    //
    // Use Clenshaw summation to evaluate
    //   m = (g(x) + g(y)) / 2         -- mean value
    //   s = (g(x) - g(y)) / (x - y)   -- average slope
    // where
    //   g(x) = sum(c[j]*sin(2*j*x), j = 1..n)
    //
    // This function returns only s and m is discarded.
    //
    // Write
    //   t = [m; s]
    //   t = sum(c[j] * f[j](x,y), j = 1..n)
    // where
    //   f[j](x,y) = [ (sin(2*j*x)+sin(2*j*y))/2 ]
    //               [ (sin(2*j*x)-sin(2*j*y))/d ]
    //
    //             = [       sin(j*p)*cos(j*d) ]
    //               [ (2/d)*sin(j*d)*cos(j*p) ]
    // and
    //    p = x+y, d = x-y
    //
    //   f[j+1](x,y) = A * f[j](x,y) - f[j-1](x,y)
    //
    //   A = [  2*cos(p)*cos(d)      -sin(p)*sin(d)*d]
    //       [ -4*sin(p)*sin(d)/d   2*cos(p)*cos(d)  ]
    //
    // Let b[n+1] = b[n+2] = [0 0; 0 0]
    //     b[j] = A * b[j+1] - b[j+2] + c[j] * I for j = n..1
    //    t =  b[1] * f[1](x,y)

    real p = x + y, d = x - y,
      cp = cos(p), cd =     cos(d),
      sp = sin(p), sd = d ? sin(d)/d : 1,
      m = 2 * cp * cd, s = sp * sd;
    // 2x2 matrices stored in row-major order
    const real a[4] = {m, -s * d * d, -4 * s, m};
    real ba[4] = {0, 0, 0, 0};
    real bb[4] = {0, 0, 0, 0};
    real* b0 = ba;
    real* b1 = bb;
    if (n > 0) b0[0] = b0[3] = c[n];
    for (int j = n - 1; j > 0; --j) { // j = n-1 .. 1
      std::swap(b0, b1);
      // b0 = A * b1 - b0 + c[j] * I
      b0[0] = a[0] * b1[0] + a[1] * b1[2] - b0[0] + c[j];
      b0[1] = a[0] * b1[1] + a[1] * b1[3] - b0[1];
      b0[2] = a[2] * b1[0] + a[3] * b1[2] - b0[2];
      b0[3] = a[2] * b1[1] + a[3] * b1[3] - b0[3] + c[j];
    }
    b1[0] = sp * cd; b1[2] = 2 * sd * cp;
    // Here is the (unused) expression for m
    // m = b0[0] * b1[0] + b0[1] * b1[2];
    s = b0[2] * b1[0] + b0[3] * b1[2];
    return s;
  }

  Math::real Rhumb::DIsometricToRectifying(real x, real y) const {
    return DConformalToRectifying(gd(x), gd(y)) * Dgd(x, y);
  }

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
          psi12 = !(abs(mu12) <= Rhumb::tol()) ?
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
    return 1 - Rhumb::SinSeries(x, y, _ell.RectifyingToConformalCoeffs(),
                                Rhumb::tm_maxord);
  }

  Math::real RhumbLine::DRectifyingToIsometric(real x, real y) const {
    real
      chix = _ell.ConformalLatitude
      (_ell.InverseRectifyingLatitude(x/Math::degree())) * Math::degree(),
      chiy = _ell.ConformalLatitude
      (_ell.InverseRectifyingLatitude(y/Math::degree())) * Math::degree();
    return Dgdinv(chix, chiy) * DRectifyingToConformal(x, y);
  }

} // namespace GeographicLib

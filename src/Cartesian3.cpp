/**
 * \file Cartesian3.cpp
 * \brief Implementation for GeographicLib::Triaxial::Cartesian3 class
 *
 * Copyright (c) Charles Karney (2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/Triaxial/Cartesian3.hpp>

namespace GeographicLib {
  namespace Triaxial {

  using namespace std;

  Cartesian3::Cartesian3(const Ellipsoid3& t)
    : _t(t)
    , _axes{_t.a(), _t.b(), _t.c()}
    , _axes2{Math::sq(_t.a()), Math::sq(_t.b()), Math::sq(_t.c())}
    , _linecc2{(_t.a() - _t.c()) * (_t.a() + _t.c()),
               (_t.b() - _t.c()) * (_t.b() + _t.c()), 0}
    , _norm(0, 1)
    , _uni(0, 1)
  {}

  Cartesian3::Cartesian3(real a, real b, real c)
    : Cartesian3(Ellipsoid3(a, b, c))
  {}

  Cartesian3::Cartesian3(real b, real e2, real k2, real kp2)
    : Cartesian3(Ellipsoid3(b, e2, k2, kp2))
  {}

  template<int n>
  pair<Math::real, Math::real> Cartesian3::funp<n>::operator()(real p)
    const {
    real fp = 0, fv = -1, fcorr = 0;
    for (int k = 0; k < 3; ++k) {
      if (_r[k] == 0) continue;
      real g = _r[k] / (p + _l[k]);
      if constexpr (n == 2) g *= g;
      real ga = round(g/_d) * _d, gb = g - ga;
      fv = fv + ga; fcorr = fcorr + gb;
      fp = fp - n * g / (p + _l[k]);
    }
    return pair<real, real>(fv + fcorr, fp);
  }

  Math::real Cartesian3::cartsolve(const function<pair<real, real>(real)>& f,
                                   real p0, real pscale) {
    // Solve
    //   f(p) = 0
    // Initial guess is p0; pscale is a scale factor for p; scale factor for f
    // = 1.  This assumes that there's a single solution with p >= 0 and that
    // for p > 0, f' < 0 and f'' > 0
    const real eps = numeric_limits<real>::epsilon(),
      tol = Math::sq(cbrt(eps)), tol2 = Math::sq(tol),
      ptol = pscale * sqrt(eps);
    real p = p0;
    real od = -1;
    for (int i = 0;
         i < maxit_ ||
           (throw_ && (throw GeographicLib::GeographicErr
                       ("Convergence failure Cartesian3::cartsolve"), false));
         ++i) {
      pair<real, real> fx = f(p);
      real fv = fx.first, fp = fx.second;
      // We're done if f(p) <= 0 on initial guess; this can happens when z = 0.
      // However, since Newton converges from below, any negative f(p)
      // indicates convergence.
      if (!(fv > tol2)) break;
      real d = -fv/fp;          // d is positive
      p = p + d;
      // converged if fv <= 8*eps (after first iteration) or
      // d <= max(eps, |p|) * tol and d <= od.
      // N.B. d is always positive.
      if ( (fv <= 8 * eps || d <= fmax(ptol, p) * tol) && d <= od )
        // The condition d <= od means that this won't trip on the first
        // iteration
        break;
      od = d;
    }
    return p;
  }

  Math::real Cartesian3::cubic(vec3 R2) const {
    // Solve sum(R2[i]/(z + lineq2[i]), i,0,2) - 1 = 0 with lineq2[2] = 0.
    // This has three real roots with just one satisifying q >= 0.
    // Express as a cubic equation z^3 + a*z^2 + b*z + c = 0.
    real c = - _linecc2[0]*_linecc2[1] * R2[2],
      b = _linecc2[0]*_linecc2[1]
      - (_linecc2[1] * R2[0] + _linecc2[0] * R2[1] +
         (_linecc2[0] + _linecc2[1]) * R2[2]),
      a = _linecc2[0] + _linecc2[1] - (R2[0] + R2[1] + R2[2]);
    bool recip = b > 0;
    if (recip) {
      // If b positive there a cancellation in p = (3*b - a^2) / 3, so
      // transform to a polynomial in 1/t.  The resulting coefficients are
      real ax = b/c, bx = a/c, cx = 1/c;
      a = ax; b = bx; c = cx;
    }
    // Reduce cubic to w^3 + p*w + q = 0, where z = w - a/3.
    // See https://dlmf.nist.gov/1.11#iii
    real p = (3*b - Math::sq(a)) / 3,
      q = (2*a*Math::sq(a) - 9*a*b + 27*c) / 27;

    // Now switch to https://dlmf.nist.gov/4.43
    // We have 3 real roots, so 4*p^3 + 27*q^2 <= 0
    real A = sqrt(fmax(real(0), -4*p/3)),
      alp = atan2(q, sqrt(fmax(real(0),
                               -(4*p*Math::sq(p)/27 + Math::sq(q)))))/3;
    // alp is in [-pi/3, pi/3]
    // z = A*sin(alp + 2*pi/3 * k) - a/3 for k = -1, 0, 1
    // for the single positive solution we pick k = 1 which gives the
    // algebraically largest result
    real t = A/2 * (cos(alp) * sqrt(real(3)) - sin(alp)) - a/3;

    return recip ? 1/t : t;
  }

  void Cartesian3::carttoellip(vec3 R, Angle& bet, Angle& omg, real& H) const {
    // tol2 = eps^(4/3)
    real tol2 = Math::sq(Math::sq(cbrt(numeric_limits<real>::epsilon())));
    vec3 R2 = {Math::sq(R[0]), Math::sq(R[1]), Math::sq(R[2])};
    real qmax = R2[0] + R2[1] + R2[2],
      qmin = fmax(fmax(R2[2], R2[1] + R2[2] - _linecc2[1]),
                  R2[0] + R2[1] + R2[2] - _linecc2[0]),
      q = qmin;
    do {                       // Executed once (provides the ability to break)
      const funp<1> f(R2, _linecc2);
      pair<real, real> fx = f(q);
      if (!( fx.first > tol2 ))
        break;                  // negative means converged
      q = fmax(qmin, fmin(qmax, cubic(R2)));
      fx = f(q);
      if (!( fabs(fx.first) > tol2 ))
        break;                  // test abs(fv) here
      q = fmax(qmin, q - fx.first/fx.second);
      q = cartsolve(f, q, Math::sq(b()));
    } while (false);
    vec3 axes = {sqrt(_linecc2[0] + q), sqrt(_linecc2[1] + q), sqrt(q)};
    _t.cart2toellipint(R, bet, omg, axes);
    H = axes[2] - c();
  }

  void Cartesian3::elliptocart(Angle bet, Angle omg, real H, vec3& R) const {
    vec3 ax;
    real shift = H * (2*c() + H);
    for (int k = 0; k < 2; ++k)
      ax[k] = sqrt(_axes2[k] + shift);
    ax[2] = c() + H;
    real tx = hypot(_t.k() * bet.c(), _t.kp()),
      tz = hypot(_t.k(), _t.kp() * omg.s());
    R = { ax[0] * omg.c() * tx,
          ax[1] * bet.c() * omg.s(),
          ax[2] * bet.s() * tz };
  }

  template<int n>
  void Cartesian3::cart2togeneric(vec3 R, ang& phi, ang& lam, bool alt) const {
    static_assert(n >= 0 && n <= 2, "Bad coordinate conversion");
    if constexpr (n == 2) {
      R[0] /= _axes2[0];
      R[1] /= _axes2[1];
      R[2] /= _axes2[2];
    } else if constexpr (n == 1) {
      R[0] /= _axes[0];
      R[1] /= _axes[1];
      R[2] /= _axes[2];
    } // else n == 0, R is unchanged
    roty(R, alt ? 1 : 0);
    // nr = [-R[2], R[1], R[0]]
    phi = ang(R[2], hypot(R[0], R[1]));
    lam = ang(R[1], R[0]);      // ang{0, 0} -> 0
  }

  template<int n>
  void Cartesian3::generictocart2(ang phi, ang lam, vec3& R, bool alt) const {
    static_assert(n >= 0 && n <= 2, "Bad coordinate conversion");
    R = {phi.c() * lam.c(), phi.c() * lam.s(), phi.s()};
    roty(R, alt ? -1 : 0);
    if constexpr (n == 2) {
      R[0] *= _axes2[0];
      R[1] *= _axes2[1];
      R[2] *= _axes2[2];
    } else if constexpr (n == 1) {
      R[0] *= _axes[0];
      R[1] *= _axes[1];
      R[2] *= _axes[2];
    } // else n == 0, R is unchanged
    if constexpr (n != 1) {
      real d = Math::hypot3(R[0] / _axes[0], R[1] / _axes[1], R[2] / _axes[2]);
      R[0] /= d; R[1] /= d; R[2] /= d;
    } // else n == 1, d = 1 and R is already on the surface of the ellipsoid
  }

  template<int n>
  Angle Cartesian3::meridianplane(ang lam, bool alt) const {
    if constexpr (n == 2)
      return lam.modang(_axes2[alt ? 2 : 0]/_axes2[1]);
    else if constexpr (n == 1)
      return lam.modang(_axes[alt ? 2 : 0]/_axes[1]);
    else {                      // n == 0
      (void) alt;               // Visual Studio 15 complains about used alt
      return lam;
    }
  }

  void Cartesian3::cardinaldir(vec3 R, ang merid, vec3& N, vec3& E,
                               bool alt) const {
    roty(R, alt ? 1 : 0);
    int i0 = alt ? 2 : 0, i1 = 1, i2 = alt ? 0 : 2;
    vec3 up = {R[0] / _axes2[i0], R[1] / _axes2[i1],
      R[2] / _axes2[i2]};
    Ellipsoid3::normvec(up);
    N = { -R[0]*R[2] / _axes2[i2], -R[1]*R[2] / _axes2[i2],
      Math::sq(R[0]) / _axes2[i0] + Math::sq(R[1]) / _axes2[i1]};
    if (R[0] == 0 && R[1] == 0) {
      real s = copysign(real(1), -R[2]);
      N = {s*merid.c(), s*merid.s(), 0};
    } else
      Ellipsoid3::normvec(N);
    // E = N x up
    E = { N[1]*up[2] - N[2]*up[1], N[2]*up[0] - N[0]*up[2],
      N[0]*up[1] - N[1]*up[0]};
    roty(E, alt ? -1 : 0);
    roty(N, alt ? -1 : 0);
  }

  template<int n>
  void Cartesian3::cart2togeneric(vec3 R, vec3 V,
                                  ang& phi, ang& lam, ang& zet,
                                  bool alt) const {
    cart2togeneric<n>(R, phi, lam, alt);
    vec3 N, E;
    cardinaldir(R, meridianplane<n>(lam, alt), N, E, alt);
    zet = ang(V[0]*E[0] + V[1]*E[1] + V[2]*E[2],
              V[0]*N[0] + V[1]*N[1] + V[2]*N[2]);
  }

  template<int n>
  void Cartesian3::generictocart2(ang phi, ang lam, ang zet,
                                  vec3& R, vec3&V, bool alt) const {
    generictocart2<n>(phi, lam, R, alt);
    vec3 N, E;
    cardinaldir(R, meridianplane<n>(lam, alt), N, E, alt);
    V = {zet.c() * N[0] + zet.s() * E[0],
      zet.c() * N[1] + zet.s() * E[1],
      zet.c() * N[2] + zet.s() * E[2]};
  }

  void Cartesian3::cart2toany(vec3 R,
                              coord coordout, Angle& lat, Angle& lon) const {
    bool alt = coordout > ELLIPSOIDAL;
    switch (coordout) {
    case GEODETIC:
    case GEODETIC_X:
      cart2togeneric<2>(R, lat, lon, alt); break;
    case PARAMETRIC:
    case PARAMETRIC_X:
      cart2togeneric<1>(R, lat, lon,  alt); break;
    case GEOCENTRIC:
    case GEOCENTRIC_X:
      cart2togeneric<0>(R, lat, lon, alt); break;
    case ELLIPSOIDAL:
      _t.cart2toellip(R, lat, lon); break;
    default:
      throw GeographicErr("Bad Cartesian3::coord value " + to_string(coordout));
    }
  }

  void Cartesian3::anytocart2(coord coordin, Angle lat, Angle lon,
                              vec3& R) const {
    bool alt = coordin > ELLIPSOIDAL;
    switch (coordin) {
    case GEODETIC:
    case GEODETIC_X:
      generictocart2<2>(lat, lon, R, alt); break;
    case PARAMETRIC:
    case PARAMETRIC_X:
      generictocart2<1>(lat, lon, R, alt); break;
    case GEOCENTRIC:
    case GEOCENTRIC_X:
      generictocart2<0>(lat, lon, R, alt); break;
    case ELLIPSOIDAL:
      _t.elliptocart2(lat, lon, R); break;
    default:
      throw GeographicErr("Bad Cartesian3::coord value " + to_string(coordin));
    }
  }

  void Cartesian3::anytoany(coord coordin, Angle lat1, Angle lon1,
                            coord coordout, Angle& lat2, Angle& lon2) const {
    vec3 R;
    anytocart2(coordin, lat1, lon1, R);
    cart2toany(R, coordout, lat2, lon2);
  }

  void Cartesian3::cart2toany(vec3 R, vec3 V, coord coordout,
                              Angle& lat, Angle& lon, Angle& azi) const {
    bool alt = coordout > ELLIPSOIDAL;
    switch (coordout) {
    case GEODETIC:
    case GEODETIC_X:
      cart2togeneric<2>(R, V, lat, lon, azi, alt); break;
    case PARAMETRIC:
    case PARAMETRIC_X:
      cart2togeneric<1>(R, V, lat, lon, azi, alt); break;
    case GEOCENTRIC:
    case GEOCENTRIC_X:
      cart2togeneric<0>(R, V, lat, lon, azi, alt); break;
    case ELLIPSOIDAL:
      _t.cart2toellip(R, V, lat, lon, azi); break;
    default:
      throw GeographicErr("Bad Cartesian3::coord value " + to_string(coordout));
    }
  }

  void Cartesian3::anytocart2(coord coordin, Angle lat, Angle lon, Angle azi,
                              vec3& R, vec3& V) const {
    bool alt = coordin > ELLIPSOIDAL;
    switch (coordin) {
    case GEODETIC:
    case GEODETIC_X:
      generictocart2<2>(lat, lon, azi, R, V, alt); break;
    case PARAMETRIC:
    case PARAMETRIC_X:
      generictocart2<1>(lat, lon, azi, R, V, alt); break;
    case GEOCENTRIC:
    case GEOCENTRIC_X:
      generictocart2<0>(lat, lon, azi, R, V, alt); break;
    case ELLIPSOIDAL:
      _t.elliptocart2(lat, lon, azi, R, V); break;
    default:
      throw GeographicErr("Bad Cartesian3::coord value " + to_string(coordin));
    }
  }

  void Cartesian3::carttoany(vec3 R, coord coordout,
                             Angle& lat, Angle& lon, real& h) const {
    switch (coordout) {
    case ELLIPSOIDAL: carttoellip(R, lat, lon, h); break;
    default:
      vec3 R2;
      carttocart2(R, R2, h);
      cart2toany(R2, coordout, lat, lon);
    }
  }

  void Cartesian3::anytocart(coord coordin, Angle lat, Angle lon, real h,
                             vec3& R) const {
    switch (coordin) {
    case ELLIPSOIDAL: elliptocart(lat, lon, h, R); break;
    default:
      vec3 R2;
      anytocart2(coordin, lat, lon, R2);
      cart2tocart(R2, h, R);
    }
  }

  void Cartesian3::cart2tocart(vec3 R2, real h, vec3& R) const {
    vec3 Rn = {R2[0] / _axes2[0], R2[1] / _axes2[1], R2[2] / _axes2[2]};
    real d = h / Math::hypot3(Rn[0], Rn[1], Rn[2]);
    R = R2;
    for (int k = 0; k < 3; ++k)
      R[k] += Rn[k] * d;
  }

  void Cartesian3::carttocart2(vec3 R, vec3& R2, real& h) const {
    const real eps = numeric_limits<real>::epsilon(), ztol = b() * eps/8;
    for (int k = 0; k < 3; ++k)
      if (fabs(R[k]) <= ztol) R[k] = copysign(real(0), R[k]);
    vec3 s = {R[0] * _axes[0], R[1] * _axes[1], R[2] * _axes[2]};
    real p = fmax(fmax(fabs(s[2]), hypot(s[1], s[2]) - _linecc2[1]),
                  Math::hypot3(s[0], s[1], s[2]) - _linecc2[0]);
    const funp<2> f(s, _linecc2);
    p = cartsolve(f, p, Math::sq(b()));
    R2 = R;
    for (int k = 0; k < 3; ++k)
      R2[k] *= _axes2[k] / (p + _linecc2[k]);
    // Deal with case p == 0 (when R2[2] is indeterminate).
    if (p == 0) {
      if (_linecc2[0] == 0)     // sphere
        R2[0] = R[0];
      if (_linecc2[1] == 0)     // sphere or prolate
        R2[1] = R[1];
      R2[2] = copysign(_axes[2], R[2]) *
        sqrt(1 - Math::sq(R2[0]/_axes[0]) - Math::sq(R2[1]/_axes[1]));
    }
    // Easily shown that U.(R-R2)  = (p - c^2) * U.U
    // => h = Uhat.(R-R2) = (p - c^2) * |U|
    h = (p - _axes2[2]) * Math::hypot3(R2[0] / _axes2[0],
                                       R2[1] / _axes2[1],
                                       R2[2] / _axes2[2]);
  }

  } // namespace Triaxial
} // namespace GeographicLib

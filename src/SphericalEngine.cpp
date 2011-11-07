/**
 * \file SphericalEngine.cpp
 * \brief Implementation for GeographicLib::SphericalEngine class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/SphericalEngine.hpp>
#include <limits>
#include <iostream>
#include <GeographicLib/CircularEngine.hpp>

#define GEOGRAPHICLIB_SPHERICALENGINE_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_SPHERICALENGINE_CPP)
RCSID_DECL(GEOGRAPHICLIB_SPHERICALENGINE_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real SphericalEngine::scale_ =
    pow(real(numeric_limits<real>::radix),
        -numeric_limits<real>::max_exponent/2);
  const Math::real SphericalEngine::eps_ =
    Math::sq(numeric_limits<real>::epsilon());

  const std::vector<Math::real> SphericalEngine::Z_(0);

  template<bool gradp, SphericalEngine::normalization norm, int L>
  Math::real SphericalEngine::Value(const coeff c[], const real f[],
                                    real x, real y, real z, real a,
                                    real& gradx, real& grady, real& gradz)
    throw() {
    // General sum
    // V(r, theta, lambda) = sum(n = 0..N) sum(m = 0..n)
    //   q^(n+1) * (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) * P[n,m](t)
    //
    // write t = cos(theta), u = sin(theta), q = a/r.
    //
    // P[n,m] is the fully normalized associated Legendre function (usually
    // denoted Pbar) of degree n and order m.
    //
    // Rewrite outer sum
    // V(r, theta, lambda) = sum(m = 0..N) * P[m,m](t) * q^(m+1) *
    //    [Sc[m] * cos(m*lambda) + Ss[m] * sin(m*lambda)]
    // = "outer sum"
    //
    // where the inner sums are
    //   Sc[m] = sum(n = m..N) q^(n-m) * C[n,m] * P[n,m](t)/P[m,m](t)
    //   Ss[m] = sum(n = m..N) q^(n-m) * S[n,m] * P[n,m](t)/P[m,m](t)
    //
    // Evaluate sums via Clenshaw method.
    //
    // The overall framework is similar to Deakin with the following changes:
    // * use fully normalized associated Legendre functions (instead of the
    //   quasi-normalized ones)
    // * Clenshaw summation is used to roll the computation of cos(m*lambda)
    //   and sin(m*lambda) into the evaluation of the outer sum (rather than
    //   independently computing an array of these trigonometric terms).
    // * Scale the coefficients to guard against overflow when N is large.
    //
    // General framework of Clenshaw;   see
    //    http://mathworld.wolfram.com/ClenshawRecurrenceFormula.html
    //
    // Let
    //    S = sum(k = 0..N) c[k] * F[k](x)
    //    F[n+1](x) = alpha[n](x) * F[n](x) + beta[n](x) * F[n-1](x)
    //
    // Evaluate S with
    //    y[N+2] = y[N+1] = 0
    //    y[k] = alpha[k] * y[k+1] + beta[k+1] * y[k+2] + c[k]
    //    S = c[0] * F[0] + y[1] * F[1] + beta[1] * F[0] * y[2]
    //
    // IF F[0](x) = 1 and beta(0,x) = 0, then F[1](x) = alpha(0,x) and
    // we can continue the recursion for y[k] until y[0]:
    //    S = y[0]
    //
    // Inner sum...
    //
    // let l = n-m; n = l+m
    // Sc[m] = sum(l = 0..N-m) C[l+m,m] * q^l * P[l+m,m](t)/P[m,m](t)
    // F[l] = q^l * P[l+m,m](t)/P[m,m](t)
    //
    // Holmes + Featherstone, Eq. (11):
    //   P[n,m] = sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m))) * t * P[n-1,m] -
    //            sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3))) * P[n-2,m]
    // thus
    //   alpha[l] = t * q * sqrt(((2*n+1)*(2*n+3))/
    //                           ((n-m+1)*(n+m+1)))
    //   beta[l+1] = - q^2 * sqrt(((n-m+1)*(n+m+1)*(2*n+5))/
    //                            ((n-m+2)*(n+m+2)*(2*n+1)))
    //
    // In this case, F[0] = 1 and beta[0] = 0 so the Sc[m] = y[0].
    //
    // Outer sum...
    //
    // V = sum(m = 0..N) Sc[m] * q^(m+1) * cos(m*lambda) * P[m,m](t)
    //   + sum(m = 0..N) Ss[m] * q^(m+1) * cos(m*lambda) * P[m,m](t)
    // F[m] = q^(m+1) * cos(m*lambda) * P[m,m](t) [or sin(m*lambda)]
    //
    // Holmes + Featherstone, Eq. (13):
    //   P[m,m] = u * sqrt((2*m+1)/((m>1?2:1)*m)) * P[m-1,m-1]
    // and
    //   cos((m+1)*lambda) = 2*cos(lambda)*cos(m*lambda) - cos((m-1)*lambda)
    // thus
    //   alpha[m] = 2*cos(lambda) * sqrt((2*m+3)/(2*(m+1))) * u * q
    //            =   cos(lambda) * sqrt( 2*(2*m+3)/(m+1) ) * u * q
    //   beta[m+1] = -sqrt((2*m+3)*(2*m+5)/(4*(m+1)*(m+2))) * u^2 * q^2
    //               * (m == 0 ? sqrt(2) : 1)
    //
    // F[0] = q                                [or 0]
    // F[1] = cos(lambda) * sqrt(3) * u * q^2  [or sin(lambda)]
    // beta[1] = - sqrt(15/4) * u^2 * q^2
    //
    // Here is how the various components of the gradient are computed
    //
    // differentiate wrt r:
    //   d q^(n+1) / dr = (-1/r) * (n+1) * q^(n+1)
    //
    // so multiply C[n,m] by n+1 in inner sum and multiply the sum by -1/r.
    //
    // differentiate wrt lambda
    //   d cos(m*lambda) = -m * sin(m*lambda)
    //   d sin(m*lambda) =  m * cos(m*lambda)
    //
    // so multiply terms by m in outer sum and swap sin and cos variables.
    //
    // differentiate wrt theta
    //  dV/dtheta = -u * dV/dt = -u * V'
    // here ' denotes differentiation wrt to t.
    //   d/dt (Sc[m] * P[m,m](t)) = Sc'[m] * P[m,m](t) + Sc[m] * P'[m,m](t)
    //
    // Now P[m,m](t) = const * u^m, so P'[m,m](t) = -m * t/u^2 * P[m,m](t),
    // thus
    //   d/dt (Sc[m] * P[m,m](t)) = (Sc'[m] - m * t/u^2 Sc[m]) * P'[m,m](t)
    //
    // Clenshaw recursion for Sc[m] reads
    //    y[k] = alpha[k] * y[k+1] + beta[k+1] * y[k+2] + c[k]
    // where alpha'[k] = alpha[k]/t, beta'[k] = c'[k] = 0.  Thus
    //    y'[k] = alpha[k] * y'[k+1] + beta[k+1] * y'[k+2] + alpha[k]/t * y[k+1]
    //
    // Finally, given the derivatives of V, we can compute the components of
    // the gradient in sphierical coordinates and transform the result into
    // cartesian coordinates.

    STATIC_ASSERT(L > 0, "L must be positive");
    STATIC_ASSERT(norm == full || norm == schmidt, "Unknown normalization");
    int N = c[0].nmx(), M = c[0].mmx();

    real
      p = Math::hypot(x, y),
      cl = p ? x / p : 1,       // cos(lambda); at pole, pick lambda = 0
      sl = p ? y / p : 0,       // sin(lambda)
      r = Math::hypot(z, p),
      t = r ? z / r : 0,            // cos(theta); at origin, pick theta = pi/2
      u = r ? max(p / r, eps_) : 1, // sin(theta); but avoid the pole
      q = a / r;
    real
      q2 = Math::sq(q),
      uq = u * q,
      uq2 = Math::sq(uq),
      tu2 = t / Math::sq(u);
    // Initialize outer sum
    real vc  = 0, vc2  = 0, vs  = 0, vs2  = 0;   // v [N + 1], v [N + 2]
    // vr, vt, vl and similar w variable accumulate the sums for the
    // derivatives wrt r, theta, and lambda, respectively.
    real vrc = 0, vrc2 = 0, vrs = 0, vrs2 = 0;   // vr[N + 1], vr[N + 2]
    real vtc = 0, vtc2 = 0, vts = 0, vts2 = 0;   // vt[N + 1], vt[N + 2]
    real vlc = 0, vlc2 = 0, vls = 0, vls2 = 0;   // vl[N + 1], vl[N + 2]
    int k[L];
    for (int m = M; m >= 0; --m) {   // m = M .. 0
      // Initialize inner sum
      real wc  = 0, wc2  = 0, ws  = 0, ws2  = 0; // w [N - m + 1], w [N - m + 2]
      real wrc = 0, wrc2 = 0, wrs = 0, wrs2 = 0; // wr[N - m + 1], wr[N - m + 2]
      real wtc = 0, wtc2 = 0, wts = 0, wts2 = 0; // wt[N - m + 1], wt[N - m + 2]
      for (int l = 0; l < L; ++l)
        k[l] = c[l].index(N, m) + 1;
      for (int n = N; n >= m; --n) {             // n = N .. m; l = N - m .. 0
        real w, A, Ax, B, R;    // alpha[l], beta[l + 1]
        switch (norm) {
        case full:
          w = real(2 * n + 1) / (real(n - m + 1) * (n + m + 1));
          Ax = q * sqrt(w * (2 * n + 3));
          A = t * Ax;
          B = - q2 * sqrt(real(2 * n + 5) / (w * (n - m + 2) * (n + m + 2)));
          break;
        case schmidt:
          w = real(n - m + 1) * (n + m + 1);
          Ax = q * (2 * n + 1) / sqrt(w);
          A = t * Ax;
          B = - q2 * sqrt(w / (real(n - m + 2) * (n + m + 2)));
          break;
        }
        R = c[0].Cv(--k[0]);
        for (int l = 1; l < L; ++l)
          R += c[l].Cv(--k[l], n, m, f[l]);
        R *= scale_;
        w = A * wc + B * wc2 + R; wc2 = wc; wc = w;
        if (gradp) {
          w = A * wrc + B * wrc2 + (n + 1) * R; wrc2 = wrc; wrc = w;
          w = A * wtc + B * wtc2 +    Ax * wc2; wtc2 = wtc; wtc = w;
        }
        if (m) {
          R = c[0].Sv(k[0]);
          for (int l = 1; l < L; ++l)
            R += c[l].Sv(k[l], n, m, f[l]);
          R *= scale_;
          w = A * ws + B * ws2 + R; ws2 = ws; ws = w;
          if (gradp) {
            w = A * wrs + B * wrs2 + (n + 1) * R; wrs2 = wrs; wrs = w;
            w = A * wts + B * wts2 +    Ax * ws2; wts2 = wts; wts = w;
          }
        }
      }
      // Now Sc[m] = wc, Ss[m] = ws
      // Sc'[m] = wtc, Ss'[m] = wtc
      if (m) {
        real v, A, B;           // alpha[m], beta[m + 1]
        switch (norm) {
        case full:
          v = 2 * real(2 * m + 3) / (m + 1);
          A = cl * sqrt(v) * uq;
          B = - sqrt((v * (2 * m + 5)) / (8 * (m + 2))) * uq2;
          break;
        case schmidt:
          v = 2 * real(2 * m + 1) / (m + 1);
          A = cl * sqrt(v) * uq;
          B = - sqrt((v * (2 * m + 3)) / (8 * (m + 2))) * uq2;
          break;
        }
        v = A * vc  + B * vc2  +  wc ; vc2  = vc ; vc  = v;
        v = A * vs  + B * vs2  +  ws ; vs2  = vs ; vs  = v;
        if (gradp) {
          // Include the terms Sc[m] * P'[m,m](t) and  Ss[m] * P'[m,m](t)
          wtc -= m * tu2 * wc; wts -= m * tu2 * ws;
          v = A * vrc + B * vrc2 +  wrc; vrc2 = vrc; vrc = v;
          v = A * vrs + B * vrs2 +  wrs; vrs2 = vrs; vrs = v;
          v = A * vtc + B * vtc2 +  wtc; vtc2 = vtc; vtc = v;
          v = A * vts + B * vts2 +  wts; vts2 = vts; vts = v;
          v = A * vlc + B * vlc2 + m*ws; vlc2 = vlc; vlc = v;
          v = A * vls + B * vls2 - m*wc; vls2 = vls; vls = v;
        }
      } else {
        real A, B, qs;
        switch (norm) {
        case full:
          A = sqrt(real(3)) * uq;       // F[1]/(q*cl) or F[1]/(q*sl)
          B = - sqrt(real(15)/4) * uq2; // beta[1]/q
          break;
        case schmidt:
          A = uq;
          B = - sqrt(real(3)/4) * uq2;
          break;
        }
        qs = q / scale_;
        vc = qs * (wc + A * (cl * vc + sl * vs ) + B * vc2);
        if (gradp) {
          qs /= r;
          // The components of the gradient in spherical coordinates are
          // r: dV/dr
          // theta: 1/r * dV/dtheta
          // lambda: 1/(r*u) * dV/dlambda
          vrc =     - qs * (wrc + A * (cl * vrc + sl * vrs) + B * vrc2);
          vtc = - u * qs * (wtc + A * (cl * vtc + sl * vts) + B * vtc2);
          vlc =   qs / u * (      A * (cl * vlc + sl * vls) + B * vlc2);
        }
      }
    }

    if (gradp) {
      // Rotate into cartesian (geocentric) coordinates
      gradx = cl * (u * vrc + t * vtc) - sl * vlc;
      grady = sl * (u * vrc + t * vtc) + cl * vlc;
      gradz =       t * vrc - u * vtc            ;
    }
    return vc;
  }

  template<bool gradp, SphericalEngine::normalization norm, int L>
  CircularEngine SphericalEngine::Circle(const coeff c[], const real f[],
                                         real p, real z, real a) {

    STATIC_ASSERT(L > 0, "L must be positive");
    STATIC_ASSERT(norm == full || norm == schmidt, "Unknown normalization");
    int N = c[0].nmx(), M = c[0].mmx();

    real
      r = Math::hypot(z, p),
      t = r ? z / r : 0,            // cos(theta); at origin, pick theta = pi/2
      u = r ? max(p / r, eps_) : 1, // sin(theta); but avoid the pole
      q = a / r;
    real
      q2 = Math::sq(q),
      tu2 = t / Math::sq(u);
    CircularEngine circ(M, gradp, norm, scale_, a, r, u, t);
    int k[L];
    for (int m = M; m >= 0; --m) {   // m = M .. 0
      // Initialize inner sum
      real wc  = 0, wc2  = 0, ws  = 0, ws2  = 0; // w [N - m + 1], w [N - m + 2]
      real wrc = 0, wrc2 = 0, wrs = 0, wrs2 = 0; // wr[N - m + 1], wr[N - m + 2]
      real wtc = 0, wtc2 = 0, wts = 0, wts2 = 0; // wt[N - m + 1], wt[N - m + 2]
      for (int l = 0; l < L; ++l)
        k[l] = c[l].index(N, m) + 1;
      for (int n = N; n >= m; --n) {             // n = N .. m; l = N - m .. 0
        real w, A, Ax, B, R;    // alpha[l], beta[l + 1]
        switch (norm) {
        case full:
          w = real(2 * n + 1) / (real(n - m + 1) * (n + m + 1));
          Ax = q * sqrt(w * (2 * n + 3));
          A = t * Ax;
          B = - q2 * sqrt(real(2 * n + 5) / (w * (n - m + 2) * (n + m + 2)));
          break;
        case schmidt:
          w = real(n - m + 1) * (n + m + 1);
          Ax = q * (2 * n + 1) / sqrt(w);
          A = t * Ax;
          B = - q2 * sqrt(w / (real(n - m + 2) * (n + m + 2)));
          break;
        }
        R = c[0].Cv(--k[0]);
        for (int l = 1; l < L; ++l)
          R += c[l].Cv(--k[l], n, m, f[l]);
        R *= scale_;
        w = A * wc + B * wc2 + R; wc2 = wc; wc = w;
        if (gradp) {
          w = A * wrc + B * wrc2 + (n + 1) * R; wrc2 = wrc; wrc = w;
          w = A * wtc + B * wtc2 +    Ax * wc2; wtc2 = wtc; wtc = w;
        }
        if (m) {
          R = c[0].Sv(k[0]);
          for (int l = 1; l < L; ++l)
            R += c[l].Sv(k[l], n, m, f[l]);
          R *= scale_;
          w = A * ws + B * ws2 + R; ws2 = ws; ws = w;
          if (gradp) {
            w = A * wrs + B * wrs2 + (n + 1) * R; wrs2 = wrs; wrs = w;
            w = A * wts + B * wts2 +    Ax * ws2; wts2 = wts; wts = w;
          }
        }
      }
      if (!gradp)
        circ.SetCoeff(m, wc, ws);
      else {
        // Include the terms Sc[m] * P'[m,m](t) and  Ss[m] * P'[m,m](t)
        wtc -= m * tu2 * wc; wts -= m * tu2 * ws;
        circ.SetCoeff(m, wc, ws, wrc, wrs, wtc, wts);
      }
    }

    return circ;
  }

  template
  Math::real SphericalEngine::Value<true, SphericalEngine::full, 1>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<false, SphericalEngine::full, 1>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<true, SphericalEngine::schmidt, 1>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<false, SphericalEngine::schmidt, 1>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);

  template
  Math::real SphericalEngine::Value<true, SphericalEngine::full, 2>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<false, SphericalEngine::full, 2>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<true, SphericalEngine::schmidt, 2>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<false, SphericalEngine::schmidt, 2>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);

  template
  Math::real SphericalEngine::Value<true, SphericalEngine::full, 3>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<false, SphericalEngine::full, 3>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<true, SphericalEngine::schmidt, 3>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);
  template
  Math::real SphericalEngine::Value<false, SphericalEngine::schmidt, 3>
  (const coeff[], const real[], real, real, real, real, real&, real&, real&);

  template
  CircularEngine SphericalEngine::Circle<true, SphericalEngine::full, 1>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<false, SphericalEngine::full, 1>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<true, SphericalEngine::schmidt, 1>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<false, SphericalEngine::schmidt, 1>
  (const coeff[], const real[], real, real, real);

  template
  CircularEngine SphericalEngine::Circle<true, SphericalEngine::full, 2>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<false, SphericalEngine::full, 2>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<true, SphericalEngine::schmidt, 2>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<false, SphericalEngine::schmidt, 2>
  (const coeff[], const real[], real, real, real);

  template
  CircularEngine SphericalEngine::Circle<true, SphericalEngine::full, 3>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<false, SphericalEngine::full, 3>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<true, SphericalEngine::schmidt, 3>
  (const coeff[], const real[], real, real, real);
  template
  CircularEngine SphericalEngine::Circle<false, SphericalEngine::schmidt, 3>
  (const coeff[], const real[], real, real, real);

} // namespace GeographicLib

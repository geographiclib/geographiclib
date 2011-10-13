/**
 * \file SphericalHarmonic.cpp
 * \brief Implementation for GeographicLib::SphericalHarmonic class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/SphericalHarmonic.hpp>
#include <iostream>
#include <iomanip>
#include <limits>

#define GEOGRAPHICLIB_SPHERICALHARMONIC_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_SPHERICALHARMONIC_CPP)
RCSID_DECL(GEOGRAPHICLIB_SPHERICALHARMONIC_HPP)

namespace GeographicLib {

  using namespace std;

  const Math::real SphericalHarmonic::scale_ =
    pow(real(numeric_limits<real>::radix),
        -numeric_limits<real>::max_exponent/2);
  const Math::real SphericalHarmonic::eps_ =
    Math::sq(numeric_limits<real>::epsilon());

  Math::real SphericalHarmonic::Value(int N,
                                      const std::vector<double>& C,
                                      const std::vector<double>& S,
                                      const std::vector<real>& Cp,
                                      const std::vector<real>& Sp,
                                      real x, real y, real z, real a) {
    // General sum
    // V(r, theta, lambda) = sum(n,0,N) sum(m,0,n)
    //   q^(n+1) * (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) * P[n,m](t)
    //
    // write t = cos(theta), u = sin(theta), q = a/r.
    //
    // P[n,m] is the fully normalized associated Legendre function (usually
    // denoted Pbar)
    //
    // Rewrite outer sum
    // V(r, theta, lambda) = sum(m,0,N) * P[m,m](t) * q^(m+1) *
    //    [Sc[m] * cos(m*lambda) + Ss[m] * sin(m*lambda)]
    // = "outer sum"
    //
    // where the inner sums are
    //   Sc[m] = sum(n,m,N) q^(n-m) * C[n,m] * P[n,m](t)/P[m,m](t)
    //   Ss[m] = sum(n,m,N) q^(n-m) * S[n,m] * P[n,m](t)/P[m,m](t)
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
    //    S = sum(c[k] * F[k](x), k = 0..N)
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
    // Sc[m] = sum(l,0,N-m) C[l+m,m] * q^l * P[l+m,m](t)/P[m,m](t)
    // F[l] = q^l * P[l+m,m](t)/P[m,m](t)
    //
    // Holmes + Featherstone, Eq. (11):
    //   P[n,m] = sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m))) * t * P[n-1,m] -
    //            sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3))) * P[n-2,m]
    // thus
    //   alpha[l] = t * q * sqrt(((2*n+1)*(2*n+3))/
    //                                    ((n-m+1)*(n+m+1)))
    //   beta[l+1] = - q^2 * sqrt(((n-m+1)*(n+m+1)*(2*n+5))/
    //                            ((n-m+2)*(n+m+2)*(2*n+1)))
    //
    // In this case, F[0] = 1 and beta[0] = 0 so the Sc[m] = y[0].
    // 
    // Outer sum...
    //
    // V = sum(m,0,N) Sc[m] * q^(m+1) * cos(m*lambda) * P[m,m](t)
    //   + sum(m,0,N) Ss[m] * q^(m+1) * cos(m*lambda) * P[m,m](t)
    // F[m] = q^(m+1) * cos(m*lambda) * P[m,m](t) [or sin(m*lambda)]
    //
    // Holmes + Featherstone, Eq. (13):
    //   P[m,m] = u * sqrt((2*m+1)/((m>1?2:1)*m)) * P[m-1,m-1]
    // and
    //   cos((m+1)*lambda) = 2*cos(lambda)*cos(m*lambda) - cos((m-1)*lambda)
    // thus
    //   alpha[m] = 2*cos(lambda) * sqrt((2*m+3)/(2*(m+1))) * sin(theta) * q
    //            =   cos(lambda) * sqrt( 2*(2*m+3)/(m+1) ) * sin(theta) * q
    //   beta[m+1] = -sqrt((2*m+3)*(2*m+5)/(4*(m+1)*(m+2))) * sin(theta)^2 * q^2
    //               * (m == 0 ? sqrt(2) : 1)
    //
    // F[0] = q                                         [or 0]
    // F[1] = cos(lambda) * sqrt(3) * sin(theta) * q^2  [or sin(lambda)]
    // beta[1] = - sqrt(15/4) * sin(theta)^2 * q^2

    // Check that N is plausible
    if (N < 0)
      throw GeographicErr("N is negative");
    if (sizeof(double) * (N + 1.0) * (N + 2.0) / 2 >
        double(numeric_limits<size_t>::max()))
      throw GeographicErr("N is too large");
    size_t k = (size_t(N + 1) * size_t(N + 2)) / 2;
    if (! (C.size() == k && S.size() == k) )
      throw GeographicErr("Vector C or S is the wrong size");
    size_t kc = Cp.size(), ks = Sp.size();
    if (kc >= k || ks >= k)
      throw GeographicErr("Vector Cp or Sp is too large");

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
      uq2 = Math::sq(uq);
    // Initialize outer sum
    real vc  = 0, vc2  = 0, vs  = 0, vs2  = 0;   // v[N + 1], v[N + 2]
    for (int m = N; m >= 0; --m) {               // m = N .. 0
      // Initialize inner sum
      real wc  = 0, wc2  = 0, ws  = 0, ws2  = 0; // w[N-m+1], w[N-m+2]
      for (int n = N; n >= m; --n) {             // n = N .. m; l = N - m .. 0
        --k;
        // alpha[l], beta[l + 1]
        real w = real(2 * n + 1) / (real(n - m + 1) * (n + m + 1)),
          Ax = q * sqrt(w * (2 * n + 3)), A = t * Ax,
          B = - q2 * sqrt(real(2 * n + 5) / (w * (n - m + 2) * (n + m + 2))),
          R = scale_ * (real(C[k]) - (k < kc ? Cp[k] : 0));
        w = A * wc + B * wc2 + R; wc2  = wc; wc  = w;
        if (m) {
          R = scale_ * (real(S[k]) - (k < ks ? Sp[k] : 0));
          w = A * ws + B * ws2 + R; ws2  = ws; ws  = w;
        }
      }
      if (m) {
        // alpha[m], beta[m + 1]
        real v = 2 * real(2 * m + 3) / (m + 1),
          A = cl * sqrt(v) * uq,
          B = - sqrt((v * (2 * m + 5)) / (8 * (m + 2))) * uq2;
        v = A * vc + B * vc2 + wc; vc2 = vc; vc  = v;
        v = A * vs + B * vs2 + ws; vs2 = vs; vs  = v;
      } else {
        real
          A = sqrt(real(3)) * uq,       // F[1]/(q*cl) or F[1]/(q*sl)
          B = - sqrt(real(15)/4) * uq2, // beta[1]/q
          qs = q / scale_;
        vc = qs * (wc +  A * (cl * vc  + sl * vs ) + B * vc2 );
      }
    }

    return vc;
  }

  Math::real SphericalHarmonic::Value(int N,
                                      const std::vector<double>& C,
                                      const std::vector<double>& S,
                                      const std::vector<real>& Cp,
                                      const std::vector<real>& Sp,
                                      real x, real y, real z,
                                      real a,
                                      real& gradx, real& grady, real& gradz) {
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

    // Check that N is plausible
    if (N < 0)
      throw GeographicErr("N is negative");
    if (sizeof(double) * (N + 1.0) * (N + 2.0) / 2 >
        double(numeric_limits<size_t>::max()))
      throw GeographicErr("N is too large");
    size_t k = (size_t(N + 1) * size_t(N + 2)) / 2;
    if (! (C.size() == k && S.size() == k) )
      throw GeographicErr("Vector C or S is the wrong size");
    size_t kc = Cp.size(), ks = Sp.size();
    if (kc >= k || ks >= k)
      throw GeographicErr("Vector Cp or Sp is too large");

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
    for (int m = N; m >= 0; --m) {               // m = N .. 0
      // Initialize inner sum
      real wc  = 0, wc2  = 0, ws  = 0, ws2  = 0; // w [N - m + 1], w [N - m + 2]
      real wrc = 0, wrc2 = 0, wrs = 0, wrs2 = 0; // wr[N - m + 1], wr[N - m + 2]
      real wtc = 0, wtc2 = 0, wts = 0, wts2 = 0; // wt[N - m + 1], wt[N - m + 2]
      for (int n = N; n >= m; --n) {             // n = N .. m; l = N - m .. 0
        --k;
        // alpha[l], beta[l + 1]
        real w = real(2 * n + 1) / (real(n - m + 1) * (n + m + 1)),
          Ax = q * sqrt(w * (2 * n + 3)), A = t * Ax,
          B = - q2 * sqrt(real(2 * n + 5) / (w * (n - m + 2) * (n + m + 2))),
          R = scale_ * (real(C[k]) - (k < kc ? Cp[k] : 0));
        w = A * wc  + B * wc2  +           R; wc2  = wc ; wc  = w;
        w = A * wrc + B * wrc2 + (n + 1) * R; wrc2 = wrc; wrc = w;
        w = A * wtc + B * wtc2 +    Ax * wc2; wtc2 = wtc; wtc = w;
        if (m) {
          R = scale_ * (real(S[k]) - (k < ks ? Sp[k] : 0));
          w = A * ws  + B * ws2  +           R; ws2  = ws ; ws  = w;
          w = A * wrs + B * wrs2 + (n + 1) * R; wrs2 = wrs; wrs = w;
          w = A * wts + B * wts2 +    Ax * ws2; wts2 = wts; wts = w;
        }
      }
      if (m) {
        // alpha[m], beta[m + 1]
        real v = 2 * real(2 * m + 3) / (m + 1),
          A = cl * sqrt(v) * uq,
          B = - sqrt((v * (2 * m + 5)) / (8 * (m + 2))) * uq2;
        wtc -= m * tu2 * wc; wts -= m * tu2 * ws;
        v = A * vc  + B * vc2  +  wc ; vc2  = vc ; vc  = v;
        v = A * vs  + B * vs2  +  ws ; vs2  = vs ; vs  = v;
        v = A * vrc + B * vrc2 +  wrc; vrc2 = vrc; vrc = v;
        v = A * vrs + B * vrs2 +  wrs; vrs2 = vrs; vrs = v;
        v = A * vtc + B * vtc2 +  wtc; vtc2 = vtc; vtc = v;
        v = A * vts + B * vts2 +  wts; vts2 = vts; vts = v;
        v = A * vlc + B * vlc2 + m*ws; vlc2 = vlc; vlc = v;
        v = A * vls + B * vls2 - m*wc; vls2 = vls; vls = v;
      } else {
        real
          A = sqrt(real(3)) * uq,       // F[1]/(q*cl) or F[1]/(q*sl)
          B = - sqrt(real(15)/4) * uq2, // beta[1]/q
          qs = q / scale_;
        vc  =       qs * (wc +  A * (cl * vc  + sl * vs ) + B * vc2 );
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

    // Rotate into cartesian (geocentric) coordinates
    gradx = cl * (u * vrc + t * vtc) - sl * vlc;
    grady = sl * (u * vrc + t * vtc) + cl * vlc;
    gradz =       t * vrc - u * vtc            ;
    return vc;
  }

} // namespace GeographicLib

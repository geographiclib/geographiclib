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
                                      real x, real y, real z, real a) {
    //
    // Evaluate sums via Clenshaw method.  See
    //    http://mathworld.wolfram.com/ClenshawRecurrenceFormula.html
    //
    // R. E. Deakin, Derivatives of the earth's potentials, Geomatics Research
    // Australasia (68), 31-60, (June 1998)
    //
    // Let
    //
    //    S = sum(c[k] * F[k](x), k = 0..N)
    //    F[n+1](x) = alpha(n,x) * F[n](x) + beta(n,x) * F[n-1](x)
    //
    // Evaluate S with
    //
    //    y[N+2] = y[N+1] = 0
    //    y[k] = alpha(k,x) * y[k+1] + beta(k+1,x) * y[k+2] + c[k]
    //    S = c[0] * F[0](x) + y[1] * F[1](x) + beta(1,x) * F[0](x) * y[2]
    //
    // IF F[0](x) = 1 and beta(0,x) = 0, then F[1](x) = alpha(0,x) and
    // we can continue the recursion for y[k] until y[0]:
    //    S = y[0]
    //
    // General sum
    // V(r, theta, lambda) = sum(n,0,N) sum(m,0,n)
    //   q^(n+1) * (Cnm * cos(m*lambda) + Snm * sin(m*lambda)) * P[n,m](theta)
    //
    // P[n,m] is the fully normalized associated Legendre function (usually
    // denoted Pbar)
    //
    // Switch order of sums:
    // V(r, theta, lambda)
    // = sum(m,0,N) sum(n,m,N)
    //   q^(n+1) * (Cnm * cos(m*lambda) + Snm * sin(m*lambda)) * P[n,m](theta)
    //
    // = sum(m,0,N) * [cos(m*lambda) (sum(n,m,N) q^(n+1) * Cnm * P[n,m](theta))
    //               + sin(m*lambda) (sum(n,m,N) q^(n+1) * Snm * P[n,m](theta))]
    //
    // = sum(m,0,N) * P[m,m](theta) * q^(m+1) *
    //  [cos(m*lambda) (sum(n,m,N) q^(n-m) * Cnm * P[n,m](theta)/P[m,m](theta))
    //  +sin(m*lambda) (sum(n,m,N) q^(n-m) * Snm * P[n,m](theta)/P[m,m](theta))]
    //
    // Inner sum...
    //
    // Let
    //   Sc[m] = (sum(n,m,N) q^(n-m) * Cnm * P[n,m](theta)/P[m,m](theta))
    //   Ss[m] = (sum(n,m,N) q^(n-m) * Snm * P[n,m](theta)/P[m,m](theta))
    //
    // let l = n-m; n = l+m
    // Sc[m] = sum(l,0,N-m) C[l+m,m] * q^l * P[l+m,m](theta)/P[m,m](theta)
    // F[l] = q^l * P[l+m,m](theta)/P[m,m](theta)
    //
    // (See Holmes + Featherstone, Eq. (11))
    // alpha[l] = cos(theta) * q * sqrt(((2*n+1)*(2*n+3))/
    //                                  ((n-m+1)*(n+m+1)))
    // beta[l+1] = - q^2 * sqrt(((n-m+1)*(n+m+1)*(2*n+5))/
    //                          ((n-m+2)*(n+m+2)*(2*n+1)))
    //
    // In this F[0] = 1 and beta[0] = 0 so the sum is given by y[0].
    // 
    // F[1] = q * P[m+1,m]/P[m,m] = cos(theta) * q * sqrt(2*m+3)
    // beta[1] = - q^2 * sqrt((2*m+5)/(4*(m+1)))
    //
    // Outer sum...
    //
    // V(r, theta, lambda) = sum(m,0,N) * P[m,m](theta) * q^(m+1) *
    //  [cos(m*lambda) * Sc[m] + sin(m*lambda) * Ss[m]]
    //
    // = sum(m,0,N) Sc[m] * q^(m+1) * cos(m*lambda) * P[m,m](theta)
    // + sum(m,0,N) Ss[m] * q^(m+1) * cos(m*lambda) * P[m,m](theta)
    //
    // Let F[m] = q^(m+1) * cos(m*lambda) * P[m,m](theta) [or sin(m*lambda)]
    //
    // (See Holmes + Featherstone, Eq. (13) and
    // cos((m+1)*lambda) = 2*cos(lambda)*cos(m*lambda) - cos((m-1)*lambda)
    // alpha[m] = 2*cos(lambda) * sqrt((2*m+3)/(2*(m+1))) * sin(theta) * q
    //          =   cos(lambda) * sqrt( 2*(2*m+3)/(m+1) ) * sin(theta) * q
    // beta[m+1] = - sqrt((2*m+3)*(2*m+5)/(4*(m+1)*(m+2))) * sin(theta)^2 * q^2
    //             * (m == 0 ? sqrt(2) : 1)
    //
    // F[0] = q                                           [or 0]
    // F[1] = cos(lambda) * sqrt(3) * sin(theta) * q^2  [or sin(lambda)]
    // beta[1] = - sqrt(15/4) * sin(theta)^2 * q^2

    // N is limited to 32766 by the requirement that the vectors C and S fit in
    // the address space of 32-bit machines.
    if (sizeof(double) * (N + 1.0) * (N + 2.0) / 2 >
        double(numeric_limits<size_t>::max()))
      throw GeographicErr("N is too large");
    size_t k = (size_t(N + 1) * size_t(N + 2)) / 2;
    if (! (C.size() == k && S.size() == k) )
      throw GeographicErr("Vectors of coefficients  are the wrong size");

    real
      p = Math::hypot(x, y),
      clam = p ? x/p : 1,    // At pole, pick lambda = 0
      slam = p ? y/p : 0,
      r = Math::hypot(z, p),
      t = r ? z/r : 0,         // At origin, pick theta = pi/2 (equator)
      u = r ? max(p/r, eps_) : 1, // Avoid the pole
      q = a/r;
    real
      q2 = Math::sq(q),
      tq = t * q,
      uq = u * q,
      uq2 = Math::sq(uq);

    // Initialize outer sum
    real vc1 = 0, vc2 = 0, vs1 = 0, vs2 = 0; // v[N + 1], v[N + 2]
    real alp, bet, w, v = 0;       // alpha, beta, temporary values of w and v
    for (int m = N; m >= 0; --m) { // m = N .. 0
      // Initialize inner sum
      real wc1 = 0, wc2 = 0, ws1 = 0, ws2 = 0; // w[N - m + 1], w[N - m + 2]
      for (int n = N; n >= m; --n) {           // n = N .. m; l = N - m .. 0
        --k;
        // alpha[l], beta[l + 1]
        alp = tq * sqrt((real(2 * n + 1) * (2 * n + 3)) /
                        (real(n - m + 1) * (n + m + 1)));
        bet = - q2 * sqrt((real(n - m + 1) * (n + m + 1) * (2 * n + 5)) /
                          (real(n - m + 2) * (n + m + 2) * (2 * n + 1)));
        w = alp * wc1 + bet * wc2 + scale_ * real(C[k]); wc2 = wc1; wc1 = w;
        w = alp * ws1 + bet * ws2 + scale_ * real(S[k]); ws2 = ws1; ws1 = w;
        //        std::cerr << "C[" << n << "," << m << "]:" << C[k]
        //                  << ",S[" << n << "," << m << "]:" << S[k] << ",\n";
      }
      // Now w1 = w[0], w2 = w[1]
      real Cv = wc1, Sv = ws1;
      if (m > 0) {
        // alpha[m], beta[m + 1]
        alp = clam * sqrt((2 * real(2 * m + 3)) / (m + 1)) * uq;
        bet = - sqrt((real(2 * m + 3) * (2 * m + 5)) /
                   (4 * real(m + 1) * (m + 2))) * uq2;
        v = alp * vc1 + bet * vc2 + Cv; vc2 = vc1; vc1 = v;
        v = alp * vs1 + bet * vs2 + Sv; vs2 = vs1; vs1 = v;
      } else {
        alp = sqrt(real(3)) * uq;          // F[1]/(q*clam) or F[1]/(q*slam)
        bet = - sqrt(real(15)/4) * uq2;    // beta[1]/q
        v = q * (Cv +  alp * (clam * vc1 + slam * vs1) + bet * vc2);
      }
    }
    if (k != 0)
      throw GeographicErr("Logic screw up");

    return v/scale_;
  }

  Math::real SphericalHarmonic::Value(int N,
                                      const std::vector<double>& C,
                                      const std::vector<double>& S,
                                      real x, real y, real z,
                                      real a,
                                      real& gradx, real& grady, real& gradz) {
    // V(x, y, z) = sum(n=0..N)[ q^(n+1) * sum(m=0..n)[
    //   (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) *
    //   Pbar[n,m](cos(theta)) ] ]
    // V(x, y, z) = sum(n=0..N)[ q^(n+1) * sum(m=0..n)[
    //   (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) *
    //   Pbar[n,m](cos(theta)) ] ]
    //
    // differentiate wrt r
    // dV/dr = sum(n=0..N)[ -(n+1)/a*q^(n+2) * sum(m=0..n)[
    //   (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) *
    //   Pbar[n,m](cos(theta)) ] ]
    //
    // differentiate wrt lambda
    // dV/dlambda = sum(n=0..N)[ q^(n+1) * sum(m=0..n)[
    //   (-C[n,m] * m*sin(m*lambda) + S[n,m] * m*cos(m*lambda)) *
    //   Pbar[n,m](cos(theta)) ] ]
    //
    // differentiate wrt theta
    // dV/dlambda = sum(n=0..N)[ q^(n+1) * sum(m=0..n)[
    //   (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) *
    //   Pbar'[n,m](cos(theta)) ] ]
    //
    // Pbar' is derivative of Pbar w.r.t. theta
    // Pbar'[n,m] = 1/u*(n*t*Pbar[n,m] - f[n,m]*Pbar[n-1,m])
    // where f[n,m] = sqrt((n-m)*(n+m)*(2*n+1)/(2*n-1))
    // Pbar'[m,m] = 1/u*(m*t*Pbar[m,m])
    //
    // Alt:
    // Pbar'[n,m] = m*t/u*Pbar[n,m] - e[n,m]*Pbar[n,m+1])
    // where e[n,m] = sqrt((n+m+1)*(n-m)/(m == 0 ? 2 : 1))
    //
    // sum(m,0,N) sum(n,m,N) q^(n-m)*C[n,m]*e[n,m]*Pbar[n,m+1]/P[m,m]
    // (m' = m+1, m=m'-1)
    // = sum(m',1,N+1) sum(n,m'-1,N)
    //   q^(n-m'+1)*C[n,m'-1]*e[n,m'-1]*Pbar[n,m'] / P[m'-1,m'-1]
    // = sum(m,1,N) * q * P[m,m]/P[m-1,m-1] *
    //      sum(n,m,N) q^(n-m) C[n,m-1]*e[n,m-1]*Pbar[n,m]/P[m,m]
    // e[n,m-1] = sqrt((n+m)*(n-m+1)/(m == 1 ? 2 : 1))
    // P[m,m]/P[m-1,m-1] = u*sqrt((2*m+1)/(2*m)) for m > 1
    // The differences in addresses C[n,m] - C[n,m-1] = N-m+1
    //
    // P'[m,m] = m/(m-1)*u*sqrt((2*m+1)/(2*m)) * P'[m-1,m-1]
    //
    // N is limited to 32766 by the requirement that the vectors C and S fit in
    // the address space of 32-bit machines.
    if (sizeof(double) * (N + 1.0) * (N + 2.0) / 2 >
        double(numeric_limits<size_t>::max()))
      throw GeographicErr("N is too large");
    size_t k = (size_t(N + 1) * size_t(N + 2)) / 2;
    if (! (C.size() == k && S.size() == k) )
      throw GeographicErr("Vectors of coefficients  are the wrong size");

    real
      p = Math::hypot(x, y),
      clam = p ? x/p : 1,    // At pole, pick lambda = 0
      slam = p ? y/p : 0,
      r = Math::hypot(z, p),
      t = r ? z/r : 0,         // At origin, pick theta = pi/2 (equator)
      u = r ? max(p/r, eps_) : 1, // Avoid the pole
      q = a/r;
    real
      q2 = Math::sq(q),
      tq = t * q,
      uq = u * q,
      uq2 = Math::sq(uq),
      tu = t / u;

    // Initialize outer sum
    real vc1 = 0, vc2 = 0, vs1 = 0, vs2 = 0;     // v[N + 1], v[N + 2]
    real vrc1 = 0, vrc2 = 0, vrs1 = 0, vrs2 = 0; // vr[N + 1], vr[N + 2]
    real vlc1 = 0, vlc2 = 0, vls1 = 0, vls2 = 0; // vl[N + 1], vl[N + 2]
    real vtc1 = 0, vtc2 = 0, vts1 = 0, vts2 = 0; // vt[N + 1], vt[N + 2]
    real wtc = 0, wts = 0;                     // previous Values of wtc1 wts1
    for (int m = N; m >= 0; --m) { // m = N .. 0
      // Initialize inner sum
      real wc1 = 0, wc2 = 0, ws1 = 0, ws2 = 0;     // w[N - m + 1], w[N - m + 2]
      real wrc1 = 0, wrc2 = 0, wrs1 = 0, wrs2 = 0; // wr[N-m+1], wr[N-m+2]
      // wt accumulates C[n,m-1]*e[n,m-1]*Pbar[n,m] (for m > 0)
      real wtc1 = 0, wtc2 = 0, wts1 = 0, wts2 = 0; // wt[N-m+1], wt[N-m+2]
      for (int n = N; n >= m; --n) {            // n = N .. m; l = N - m .. 0
        --k;
        // alpha[l], beta[l + 1]
        real w = real(2 * n + 1) / (real(n - m + 1) * (n + m + 1)),
          alp = tq * sqrt(w * (2 * n + 3)),
          bet = - q2 * sqrt(real(2 * n + 5) / (w * (n - m + 2) * (n + m + 2))),
          R = scale_ * real(C[k]);
        w = alp * wc1  + bet * wc2  +           R; wc2  = wc1 ; wc1  = w;
        w = alp * wrc1 + bet * wrc2 + (n + 1) * R; wrc2 = wrc1; wrc1 = w;
        if (m) {
          R = scale_ * real(S[k]);
          w = alp * ws1  + bet * ws2  +           R; ws2  = ws1 ; ws1  = w;
          w = alp * wrs1 + bet * wrs2 + (n + 1) * R; wrs2 = wrs1; wrs1 = w;
          // e[n,m-1]
          real e = sqrt((real(n + m) * (n - m + 1)) / real(m > 1 ? 1 : 2));
          w = alp * wtc1 + bet * wtc2 + e * scale_ * real(C[k - (N - m + 1)]);
          wtc2 = wtc1; wtc1 = w;
          w = alp * wts1 + bet * wts2 + e * scale_ * real(S[k - (N - m + 1)]);
          wts2 = wts1; wts1 = w;
        }
      }
      // Now w1 = w[0], w2 = w[1]
      real Cv = wc1, Sv = ws1;
      real Cvr = wrc1, Svr = wrs1;
      // Pbar'[n,m] = m*t/u*Pbar[n,m] - e[n,m]*Pbar[n,m+1])
      real e = uq * sqrt(real(2 * m + 3) * (m ? 1 : 2) / real(2 * m + 2)),
        Cvt = m * tu * wc1 - e * wtc,
        Svt = m * tu * ws1 - e * wts;
      wtc = wtc1; wts = wts1;   // Save values of wt[cs]1 for next time
      if (m > 0) {
        // alpha[m], beta[m + 1]
        real v,
          alp = clam * sqrt((2 * real(2 * m + 3)) / (m + 1)) * uq,
          bet = - sqrt((real(2 * m + 3) * (2 * m + 5)) /
                       (4 * real(m + 1) * (m + 2))) * uq2;
        v = alp * vc1 + bet * vc2 + Cv; vc2 = vc1; vc1 = v;
        v = alp * vs1 + bet * vs2 + Sv; vs2 = vs1; vs1 = v;
        v = alp * vrc1 + bet * vrc2 + Cvr; vrc2 = vrc1; vrc1 = v;
        v = alp * vrs1 + bet * vrs2 + Svr; vrs2 = vrs1; vrs1 = v;
        v = alp * vlc1 + bet * vlc2 + m * Sv; vlc2 = vlc1; vlc1 = v;
        v = alp * vls1 + bet * vls2 - m * Cv; vls2 = vls1; vls1 = v;
        v = alp * vtc1 + bet * vtc2 + Cvt; vtc2 = vtc1; vtc1 = v;
        v = alp * vts1 + bet * vts2 + Svt; vts2 = vts1; vts1 = v;
      } else {
        real
          alp = sqrt(real(3)) * uq,       // F[1]/(q*clam) or F[1]/(q*slam)
          bet = - sqrt(real(15)/4) * uq2, // beta[1]/q
          qs = q / scale_;
        vc1 = qs * (Cv +  alp * (clam * vc1 + slam * vs1) + bet * vc2);
        qs /= r;
        vrc1 = -qs * (Cvr +  alp * (clam * vrc1 + slam * vrs1) + bet * vrc2);
        vlc1 = qs / u * (alp * (clam * vlc1 + slam * vls1) + bet * vlc2);
        vtc1 = qs * (Cvt +  alp * (clam * vtc1 + slam * vts1) + bet * vtc2);
      }
    }
    if (k != 0)
      throw GeographicErr("Logic screw up");

    gradx = clam * (u * vrc1 + t * vtc1) - slam * vlc1;
    grady = slam * (u * vrc1 + t * vtc1) + clam * vlc1;
    gradz =         t * vrc1 - u * vtc1               ;
    return vc1;
  }

} // namespace GeographicLib

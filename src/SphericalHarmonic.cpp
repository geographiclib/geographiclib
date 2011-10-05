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

  Math::real SphericalHarmonic::Value(const std::vector<double>& C,
                                          const std::vector<double>& S,
                                          int N,
                                          work tr,    // cos(theta)
                                          work ur,    // sin(theta),
                                          work clamr, // cos(lambda)
                                          work slamr, // sin(lambda)
                                          work qr     // a/r
                                          ) {
    //
    // Evaluate sums via Clenshaw method.  See
    //    http://mathworld.wolfram.com/ClenshawRecurrenceFormula.html
    //
    // R. E. Deakin, Derivatives of the earth's potentials, Geomatics Research
    // Australasia (68), 31-60. , (June 1998)
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
    // If F[0](x) = 1 and beta(0,x) = 0, then F[1](x) = alpha(0,x) and
    //    S = alpha(0,x) * y[1] + beta(1,x) * y[2] + c[0] = y[0]
    //
    // Let F[n](theta) = P[n,m](theta) where .
    //
    // General sum
    // V(r, theta, lambda) = sum(n,0,N) sum(m,0,n)
    //   q^n * (Cnm * cos(m*lambda) + Snm * sin(m*lambda)) * P[n,m](theta)
    //
    // P[n,m] is the fully normalized associated Legendre function (usually
    // denoted Pbar)
    //
    // Switch order of sums:
    // V(r, theta, lambda)
    // = sum(m,0,N) sum(n,m,N)
    //   q^n * (Cnm * cos(m*lambda) + Snm * sin(m*lambda)) * P[n,m](theta)
    //
    // = sum(m,0,N) * [cos(m*lambda) (sum(n,m,N) q^n * Cnm * P[n,m](theta))
    //               + sin(m*lambda) (sum(n,m,N) q^n * Snm * P[n,m](theta))]
    //
    // = sum(m,0,N) * P[m,m](theta) * q^m *
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
    // alpha[l] = cos(theta) * q *
    //          sqrt((2*m+2*l+1)*(2*m+2*l+3)/((l+1)*(2*m+l+1)))
    // beta[l+1] = - q^2 * sqrt((l+1)*(2*m+l+1)*(2*m+2*l+5)/
    //                         ((l+2)*(2*m+l+2)*(2*m+2*l+1)))
    //
    // In this F[0] = 1 and beta[0] = 0 so the sum is given by y[0].
    // 
    // F[1] = q * P[m+1,m]/P[m,m] = cos(theta) * q * sqrt(2*m+3)
    // beta[1] = - q^2 * sqrt((2*m+5)/(4*(m+1)))
    //
    // Outer sum...
    //
    // V(r, theta, lambda) = sum(m,0,N) * P[m,m](theta) * q^m *
    //  [cos(m*lambda) * Sc[m] + sin(m*lambda) * Ss[m]]
    //
    // = sum(m,0,N) Sc[m] * q^m * cos(m*lambda) * P[m,m](theta)
    // + sum(m,0,N) Ss[m] * q^m * cos(m*lambda) * P[m,m](theta)
    //
    // Let F[m] = q^m * cos(m*lambda) * P[m,m](theta) [or sin(m*lambda)]
    //
    // (See Holmes + Featherstone, Eq. (13) and
    // cos((m+1)*lambda) = 2*cos(lambda)*cos(m*lambda) - cos((m-1)*lambda)
    // alpha[m] = 2*cos(lambda) * sqrt((2*m+3)/(2*(m+1))) * sin(theta) * q
    //          =   cos(lambda) * sqrt( 2*(2*m+3)/(m+1) ) * sin(theta) * q
    // beta[m+1] = - sqrt((2*m+3)*(2*m+5)/(4*(m+1)*(m+2))) * sin(theta)^2 * q^2
    //             * (if m = 0 then sqrt(2) else 1)
    //
    // F[0] = 1                                           [or 0]
    // F[1] = cos(lambda) * sqrt(3) * sin(theta) * q  [or sin(lambda)]
    // beta[1] = - sqrt(15/4) * sin(theta)^2 * q^2

    // N is limited to 32766 by the requirement that the vectors C and S fit in
    // the address space of 32-bit machines.
    if (sizeof(double) * (N + 1.0) * (N + 2.0) / 2 >
        double(std::numeric_limits<size_t>::max()))
      throw GeographicErr("N is too large");
    size_t k = (size_t(N + 1) * size_t(N + 2)) / 2;
    if (! (C.size() >= k && S.size() >= k) )
      throw GeographicErr("Vectors of coefficients too small");

    T
      t =  T(tr),
      u =  T(ur),
      clam =  T(clamr),
      slam =  T(slamr);
    T scale = std::pow(2.0, -512);
    T q = Math::hypot(t, u);
    t /= q;
    u /= q;
    q = Math::hypot(clam, slam);
    clam /= q;
    slam /= q;
    q =  T(qr);
    T
      q2 = Math::sq(q),
      tq = t * q,
      uq = u * q,
      uq2 = Math::sq(uq);

    // Initialize outer sum
    T vc1 = 0, vc2 = 0, vs1 = 0, vs2 = 0; // v[N + 1], v[N + 2]
    T a, b, y, v = 0;              // alpha, beta, temporary values of y and v
    for (int m = N; m >= 0; --m) { // m = N .. 0
      // Initialize inner sum
      T yc1 = 0, yc2 = 0, ys1 = 0, ys2 = 0; // y[N - m + 1], y[N - m + 2]
      for (int l = N - m; l >= 0; --l) {     // l = N - m .. 0
        --k;
        // alpha[l], beta[l + 1]
        a = tq * sqrt((T(2 * m + 2 * l + 1) * T(2 * m + 2 * l + 3)) /
                      (T(l + 1) * T(2 * m + l + 1)));
        b = - q2 * sqrt((T(l + 1) * T(2 * m + l + 1) * T(2 * m + 2 * l + 5)) /
                        (T(l + 2) * T(2 * m + l + 2) * T(2 * m + 2 * l + 1)));
        if (false && m == 1286 && l == 0)
          std::cerr << a << " " << b << " "
                    << yc1 << " " << yc2 << " "
                    << ys1 << " " << ys2 << "\n";
        y = a * yc1 + b * yc2 + scale * T(C[k]); yc2 = yc1; yc1 = y;
        y = a * ys1 + b * ys2 + scale * T(S[k]); ys2 = ys1; ys1 = y;
        //        if (Math::isnan(y))
        //          std::cerr << m << " " << l << "\n";
      }
      // Now y1 = y[0], y2 = y[1]
      T Cv = yc1, Sv = ys1;
      if (m > 0) {
        // alpha[m], beta[m + 1]
        a = clam * sqrt((2 * T(2 * m + 3)) / T(m + 1)) * uq;
        b = - sqrt((T(2 * m + 3) * T(2 * m + 5)) /
                   (4 * T(m + 1) * T(m + 2))) * uq2;
        v = a * vc1 + b * vc2 + Cv; vc2 = vc1; vc1 = v;
        v = a * vs1 + b * vs2 + Sv; vs2 = vs1; vs1 = v;
      } else {
        a = sqrt(T(3)) * uq;          // F[1]/clam or F[1]/slam
        b = - sqrt(T(15)/4) * uq2; // beta[1]
        v = Cv +  a * (clam * vc1 + slam * vs1) + b * vc2;
      }
    }
    if (k != 0)
      throw GeographicErr("Logic screw up");

    return work(v/scale);
  }

} // namespace GeographicLib

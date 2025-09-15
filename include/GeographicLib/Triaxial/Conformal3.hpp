/**
 * \file Conformal3.hpp
 * \brief Header for GeographicLib::Triaxial::Conformal3 class
 *
 * Copyright (c) Charles Karney (2014-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CONFORMAL3_HPP)
#define GEOGRAPHICLIB_CONFORMAL3_HPP 1

#include <GeographicLib/Angle.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

namespace GeographicLib {
  namespace Triaxial {

  /**
   * \brief Jacobi's conformal projection of a triaxial ellipsoid
   *
   * This is a conformal projection of the ellipsoid to a plane in which the
   * grid lines are straight; see Jacobi,
   *  <a href="https://books.google.com/books?id=ryEOAAAAQAAJ&pg=PA212">
   * Vorlesungen &uuml;ber Dynamik, &sect;28</a>.  The constructor takes the
   * semi-axes of the ellipsoid (which must be in order).  Member functions map
   * the ellipsoidal coordinates &omega; and &beta; separately to \e x and \e
   * y.  Jacobi's coordinates have been multiplied by \f$eb/2\f$ so that the
   * customary results are returned in the cases of a sphere or an ellipsoid of
   * revolution.  In particular the scale satisfies \f$k\ge 1\f$ with \f$k =
   * 1\f$ at \f$\cos\beta = \sin\omega = 1\f$.
   *
   * The ellipsoid is oriented so that the major principal ellipse, \f$Z=0\f$,
   * is the equator, \f$\beta=0\f$, while the minor principal ellipse,
   * \f$X=0\f$, is the meridian, \f$\omega\pm\pi/2\f$.  The four umbilic
   * points, \f$\sin\omega = \cos\beta = 0\f$, lie on median principal ellipse
   * in the plane \f$Y=0\f$.
   *
   * For more information on this projection, see \ref jacobi.
   *
   * Example of use:
   * \include example-Conformal3.cpp
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT Conformal3 {
  private:
    typedef Math::real real;
    Ellipsoid3 _t, _s, _s0;
    EllipticFunction _ex, _ey, _exs, _eys;
    real _x, _y;

    real a() const { return _t.a(); }
    real b() const { return _t.b(); }
    real c() const { return _t.c(); }
    real e2() const { return _t.e2(); }
    real k2() const { return _t.k2(); }
    real kp2() const { return _t.kp2(); }
    static real Pi(const EllipticFunction& ell, Angle phi);
    static Angle Piinv(const EllipticFunction& ell, real x);
    static real F(const EllipticFunction& ell, Angle phi);
    static Angle Finv(const EllipticFunction& ell, real x);
    static Angle omegashift(Angle omg, int dir) {
      return omg - Angle::cardinal(dir);
    }
    // Jacobi has m = 2/sqrt(lambda2 - lambda3)
    // Setting lambda2 = Cayley's k = -(b^2 \sin^2\beta + c^2 \cos^2\beta)
    //                              = -(b^2 - (b^2 - c^2) \cos^2\beta)
    // Setting lambda3 = Cayley's h = -(b^2 \cos^2\omega + a^2 \sin^2\omega)
    //                              = -(b^2 + (a^2 - b^2) \sin^2\omega)
    // lambda2 - lambda3 =
    //      - (b^2 - (b^2 - c^2) \cos^2\beta)
    //      + (b^2 + (a^2 - b^2) \sin^2\omega)
    // (a^2 - c^2) * (k'^2 * \sin^2\omega + k^2 \cos^2\beta)
    // m = 2/sqrt(a^2 - c^2) * 1/sqrt(k'^2 * \sin^2\omega + k^2 \cos^2\beta)
    real scale(Angle bet, Angle omg) const {
      using std::sqrt;
      omg = omegashift(omg, +1);
      return 1/sqrt(k2() * Math::sq(bet.c()) + kp2() * Math::sq(omg.c()));
    }
    static Ellipsoid3 EquivSphere(real x, real y,
                                  EllipticFunction& ellx,
                                  EllipticFunction& elly);
  public:
    Conformal3(const Ellipsoid3& t);
    /**
     * Constructor for a trixial ellipsoid with semi-axes.
     *
     * @param[in] a the largest semi-axis.
     * @param[in] b the middle semi-axis.
     * @param[in] c the smallest semi-axis.
     *
     * The semi-axes must satisfy \e a &ge; \e b &ge; \e c > 0 and \e a >
     * \e c.
     **********************************************************************/
    Conformal3(real a, real b, real c);
    /**
     * Alternate constructor for a triaxial ellipsoid.
     *
     * @param[in] b the middle semi-axis.
     * @param[in] e2 the smallest semi-axis.
     * @param[in] k2 the relative magnitude of \e a &minus; \e b.
     * @param[in] kp2 the relative magnitude of \e b &minus; \e c.
     *
     * This form can be used to specify a sphere.
     **********************************************************************/
    Conformal3(real b, real e2, real k2, real kp2);
    /**
     * @return the quadrant length in the \e x direction (in meters).
     **********************************************************************/
    Math::real x0() const { return _x; }
    /**
     * The \e x projection.
     *
     * @param[in] omg the longitude &omega; as an Angle.
     * @return the easting (in meters).
     **********************************************************************/
    Math::real x(Angle omg) const;
    /**
     * The \e x projection.
     *
     * @param[in] omg the longitude &omega; (in degrees).
     * @return \e x (in meters).
     **********************************************************************/
    Math::real x(real omg) const {
      return x(Angle(omg));
    }
    /**
     * The inverse \e x projection.
     *
     * @param[in] x the easting (in meters).
     * @return the longitude &omega; as an Angle.
     **********************************************************************/
    Angle omega(real x) const;
    /**
     * The inverse \e x projection.
     *
     * @param[in] x the easting (in meters).
     * @return the longitude &omega; in degrees.
     **********************************************************************/
    Math::real omegad(real x) const {
      Angle omg = omega(x);
      return real(omg);
    }
    /**
     * @return the quadrant length in the \e y direction (in meters).
     **********************************************************************/
    Math::real y0() const { return _y; }
    /**
     * The \e y projection.
     *
     * @param[in] bet the latitude &beta; as an Angle
     * @return the northing (in meters).
     **********************************************************************/
    Math::real y(Angle bet) const;
    /**
     * The \e y projection.
     *
     * @param[in] bet the latitude &beta; (in degrees).
     * @return the northing (in meters).
     **********************************************************************/
    Math::real y(real bet) const {
      return y(Angle(bet));
    }
    /**
     * The inverse \e y projection.
     *
     * @param[in] y the northing (in meters).
     * @return the latitude &beta; as an Angle.
     **********************************************************************/
    Angle beta(real y) const;
    /**
     * The inverse \e y projection.
     *
     * @param[in] y the northing (in meters).
     * @return the latitude &beta; in degrees.
     **********************************************************************/
    Math::real betad(real y) const {
      Angle bet = beta(y);
      return real(bet);
    }
    void Forward(Angle bet, Angle omg, real& x, real& y) const;
    void Forward(real bet, real omg, real& x, real& y) const {
      Forward(Angle(bet), Angle(omg), x, y);
    }
    void Forward(Angle bet, Angle omg, real& x, real& y, real& m) const;
    void Forward(real bet, real omg, real& x, real& y, real& m) const {
      Forward(Angle(bet), Angle(omg), x, y, m);
    }
    void Reverse(real x, real y, Angle& bet, Angle& omg) const;
    void Reverse(real x, real y, real& bet, real& omg) const {
      Angle beta, omga;
      Reverse(x, y, beta, omga);
      bet = real(beta); omg = real(omga);
    }
    void Reverse(real x, real y, Angle& bet, Angle& omg, real& m) const;
    void Reverse(real x, real y, real& bet, real& omg, real& m) const {
      Angle beta, omga;
      Reverse(x, y, beta, omga, m);
      bet = real(beta); omg = real(omga);
    }
    real SphereRadius() const { return _s.b(); };
    void ForwardSphere(Angle bet, Angle omg, Angle& phi, Angle& lam,
                       Angle& gamma, real& m) const;
    void ForwardSphere(real bet, real omg, real& phi, real& lam,
                       real& gamma, real& m) const {
      Angle phia, lama, gammaa;
      ForwardSphere(Angle(bet), Angle(omg), phia, lama, gammaa, m);
      phi = real(phia); lam = real(lama); gamma = real(gammaa);
    }
    void ReverseSphere(Angle phi, Angle lam, Angle& bet, Angle& omg,
                       Angle& gamma, real& m) const;
    void ReverseSphere(real phi, real lam, real& bet, real& omg,
                       real& gamma, real& m) const {
      Angle beta, omga, gammaa;
      ReverseSphere(Angle(phi), Angle(lam), beta, omga, gammaa, m);
      bet = real(beta); omg = real(omga); gamma = real(gammaa);
    }
    const Ellipsoid3& t() { return _t; }
  };

  } // namespace Triaxial
} // namespace GeographicLib
#endif  // GEOGRAPHICLIB_CONFORMAL3_HPP

#if 0
function [k2,kp2,b] = ksolve(lx, ly)
% find b, k2, kp2 s.t.
% b*K(kp2) = lx
% b*K(k2)  = ly
% K(kp2)/K(k2) = lx/ly
% lx * K(k2) - ly * K(kp2) = 0
%  diff(lx*K(k) - ly*K(kp), k2);
%  ((lx*(E(k)-K(k)*kp2)) + ly*(E(kp)-k2*K(kp))) / (2*k2*kp2)

  swap = lx < ly;
  if swap, [lx, ly] = deal(ly, lx); end
  % Now lx >= ly, k2 <= 1/2
  s = lx + ly;
  nx = lx/s; ny = ly/s;
  % Get good estimate for lx/ly large and k2 small, the K(k2) = pi/2
  % K(kp2) = ny/ny * pi/2
  n = (log(4)-log(pi))/(pi/2-log(4));
  b = exp(n*pi/2)-4^n;
  KK = nx/ny * pi/2;
  k2 = 16/(exp(n*KK) - b)^(2/n);
  k2 = min(1/2,k2);
  for i=1:100
    kp2 = 1-k2;
    [K, E] = ellipke(k2);
    [Kp, Ep] = ellipke(kp2);
    f = nx*K - ny*Kp;
    fp = ((nx*(E-kp2*K)) + ny*(Ep-k2*Kp)) / (2*k2*kp2);
    %    [k2, kp2, f]
    k2n = k2 - f/fp;
    if (k2n >= 1)
      k2 = (1 + k2)/2;
    elseif (k2n <= 0)
      k2 = k2/2;
    else
      dd = k2 - k2n;
      k2 = k2n;
      if (abs(f) < 1e-10)
        i;
        break
      end
    end
  end
  kp2 = 1 - k2;
  [K, E] = ellipke(k2);
  [Kp, Ep] = ellipke(kp2);
  b = ly/K;
  if swap, [k2, kp2] = deal(kp2, k2); end
end

Need to find k2, s.t.,
  a*K(kp) - b*K(sqrt(k)) = 0

  For estimate of root use https://arxiv.org/abs/2505.17159v4
  K(k) = 1/n*log((4/kp)^n+b);
k(K) = sqrt(1 - 16/(exp(n*K) - b)^(2/n))
  kp(K) = sqrt(16/(exp(n*K) - b)^(2/n))
        = 4/(exp(n*K) - b)^(1/n))
  where
  n=(log(4)-log(pi))/(pi/2-log(4))
  b=exp(n*pi/2)-4^n
  gradef(K(k),(E(k)-(1-k^2)*K(k))/(k*(1-k^2)));
  diff(a*K(kp) - b*K(k), k2);
  (-a*(E(kp)-k2*K(kp))-(b*(E(k)-K(k)*kp2)))/(2*k2*kp2)
#endif

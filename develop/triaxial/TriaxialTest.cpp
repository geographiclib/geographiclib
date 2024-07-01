#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include "Angle.hpp"
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"
#include "TriaxialODE.hpp"

using namespace GeographicLib;
using namespace std;

void ODEtest(Math::real a, Math::real b, Math::real c) {
    typedef Math::real real;
    typedef Triaxial::vec3 vec3;
    Triaxial t(a, b, c);
    vec3 r1{t.axes[0] * t.kp, 0, t.axes[2] * t.k},
      v1{0,1,0};
    real s12 = 4 * EllipticFunction::RG(Math::sq(t.axes[0]),
                                        Math::sq(t.axes[2]));
    vec3 r2, v2;
    int kmax = Math::digits()-4;
    TriaxialODE direct(t, r1, v1);
    if (1) {
      for (int k = 10; k <= kmax; ++k) {
        real eps = pow(real(2), -k);
        int n = direct.Position(s12, r2, v2, eps);
        cout << k << " " << n << " "
             << sqrt(Math::sq(r1[0]+r2[0]) +
                     Math::sq(r1[1]+r2[1]) +
                     Math::sq(r1[2]+r2[2])) / t.b << endl;
      }
    } else {
      int n = 0, imax = 1000;
#if GEOGRAPHICLIB_PRECISION >= 4
      real eps = pow(numeric_limits<real>::epsilon(), real(0.875));
#else
      real eps = numeric_limits<real>::epsilon() * 16;
#endif
      for (int i = 0; i < imax; ++i) {
        n += direct.Position(s12, r2, v2, eps);
        r1 = r2;
        v1 = v2;
      }
      cout << n/real(imax) << "\n";
    }
}

/*
void TriaxialTest0() {
  typedef Math::real real;
  // Triaxial t(sqrt(real(2)), 1, sqrt(real(0.5)));
  Triaxial t(1.01,1,0.8);
  // Triaxial t(sqrt(real(3)), 1, 1/sqrt(real(3)));
  real mu = real(0.0000001);
  //  cout << t.k2 << " " << t.kp2 << " " << t.e2 << "\n";
  geod_fun fa(t.k2, t.kp2, -t.e2, mu, true);
  fa.NCoeffsInv();
  cerr << fa.NCoeffs() << " "
            << fa.NCoeffsInv() << " "
            << fa.InvCounts().first << "\n";
  geod_fun fb(t.k2, t.kp2, -t.e2, mu, false);
  fb.NCoeffsInv();
  cerr << fb.NCoeffs() << " "
            << fb.NCoeffsInv() << " "
            << fb.InvCounts().first << "\n";
  for (int k = 0; k <= 360; k += 3) {
    real x = k * Math::degree(),
      u = fa.ell.F(x),
      fu = fa.fun(u),
      fx = fb.fun(x),
      uu = fa.fun.inv(fu),
      xx = fb.fun.inv(fx);
      cout << k << " " << fu << " " << fx << " " << fu - fx << " "
                << u - uu << " " << x - xx << "\n";
  }
}
*/
void TriaxialTest1(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  Triaxial t(a, b, c);
  real k2 = t.k2, kp2 = t.kp2, e2 = t.e2;
  if (0) {
    int num = 100, numk = int(round(num*k2)), numkp = num-numk;
    for (int k = -numkp; k <= numk; ++k) {
      if (k == 0) continue;
      real gam = k < 0 ?
        - kp2 * Math::sq(k / real(numkp)) : k2 * Math::sq(k / real(numk));
      if (0) {
        geod_fun fphia(k2, kp2, e2, -gam, false); fphia.NCoeffsInv();
        geod_fun fphib(k2, kp2, e2, -gam, true); fphib.NCoeffsInv();
        geod_fun fomga(kp2, k2, -e2, gam, false); fomga.NCoeffsInv();
        geod_fun fomgb(kp2, k2, -e2, gam, true); fomgb.NCoeffsInv();
         cout << k << " " << gam << " "
                   << fphia.NCoeffs() << " " << fphia.NCoeffsInv() << " "
                   << fomga.NCoeffs() << " " << fomga.NCoeffsInv() << " "
                   << fphib.NCoeffs() << " " << fphib.NCoeffsInv() << " "
                   << fomgb.NCoeffs() << " " << fomgb.NCoeffsInv() << "\n";
      } else {
        geod_fun fphia(k2, kp2, e2, -gam); fphia.NCoeffsInv();
        geod_fun fomga(kp2, k2, -e2, gam); fomga.NCoeffsInv();
         cout << k << " " << gam << " "
                   << fphia.txp() << " "
                   << fphia.NCoeffs() << " " << fphia.NCoeffsInv() << " "
                   << fomga.txp() << " "
                   << fomga.NCoeffs() << " " << fomga.NCoeffsInv() << "\n";
      }
    }
  }
  for (int k = 3; k <= 16; ++k) {
    for (int s = -1; s <= 1; s += 2) {
      real gam = pow(real(10), -k) * s;
      geod_fun fphia(k2, kp2, e2, -gam); fphia.NCoeffsInv();
      geod_fun fomga(kp2, k2, -e2, gam); fomga.NCoeffsInv();
       cout << k << " " << gam << " "
                 << fphia.HalfPeriod() << " "
                 << fphia.NCoeffs() << " " << fphia.NCoeffsInv() << " "
                 << fomga.HalfPeriod() << " "
                 << fomga.NCoeffs() << " " << fomga.NCoeffsInv() << "\n";
    }
  }
}

void DirectfunTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  typedef Angle ang;
  Triaxial t(a, b, c);
  // Circumpolar
  TriaxialLine lca(t, ang::degrees(real(1)), ang::degrees(real(0)),
                  ang::degrees(real(90)));
  if (1) {
  TriaxialLine lcb(t, ang::degrees(real(89.999)),
                   ang::degrees(real(0)),
                  ang::degrees(real(90)));
  }
  // Umbilic
  TriaxialLine lu(t, ang::degrees(real(90)), ang::degrees(real(0)),
                  ang::degrees(real(135)));
  if (1) {
  // Circumpolar
  TriaxialLine ltb(t, ang::degrees(real(90)),
                   ang::degrees(real(0.001)),
                   ang::degrees(real(180)));
  }
  TriaxialLine lta(t, ang::degrees(real(90)), ang::degrees(real(89)),
                  ang::degrees(real(180)));
}

void PositionTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  typedef Angle ang;
  Triaxial t(a, b, c);
  ang bet1, omg1, alp1, bet2, omg2, alp2;
  bet1 = ang::degrees(-0.0); omg1 = ang::degrees(105);
  alp1 = ang::degrees(90);
  //  bet1 = ang::degrees(1); omg1 = ang::degrees(2);
  //  alp1 = ang::degrees(3);
  TriaxialLine l(t, bet1, omg1, alp1);
  cout << fixed << setprecision(6);
  real ds = 1/real(10);
  //  for (int s12 = 0; s12 <= 30; ++s12) {
  {real s12 = 0.996504/ds;
    l.Position(s12*ds, bet2, omg2, alp2);
    (void) Triaxial::AngNorm(bet2, omg2, alp2);
    cout << s12 << " "
         << omg2.degrees() << "\n";
      //         << bet2.degrees() << " "
      //         << alp2.degrees() << "\n";
  }
}

Math::real HybridOLD(const Triaxial& t,
                     const Angle& bet1, const Angle& omg1,
                     const Angle& alp1,
                     const Angle& bet2, const Angle& omg2) {
  typedef Math::real real;
  typedef Angle ang;
  TriaxialLine l(t, bet1, omg1, alp1);
  real s12;
  ang bet2a, omg2a, alp2a;
  l.Hybrid(bet2, bet2a, omg2a, alp2a, s12);
  (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
  omg2a -= omg2;
  return omg2a.radians0();
}

Math::real Hybrid0(const TriaxialLineF& l,
                   const TriaxialLineF::ics ic,
                   const Angle& bet2,
                   const Angle& omg2) {
  typedef Angle ang;
  ang bet2a, omg2a, alp2a;
  (void) l.Hybrid(ic, bet2, bet2a, omg2a, alp2a);
  (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
  omg2a -= omg2;
  return omg2a.radians0();
}

Math::real HybridA(const Triaxial& t,
                   const Angle& bet1, const Angle& omg1,
                   const Angle& alp1,
                   const Angle& bet2, const Angle& omg2) {
  Triaxial::gamblk gam(t, bet1, omg1, alp1);
  TriaxialLineF l(t, gam, 0.5, 1.5);
  TriaxialLineF::ics ic(l, bet1, omg1, alp1);
  return Hybrid0(l, ic, bet2, omg2);
  /*
  cout << "AA "
       << t.a << " " << t.b << " " << t.c << " "
       << gam.gam << " "
       << b1.degrees() << " "
       << o1.degrees() << " "
       << a1.degrees() << " "
       << bet2.degrees() << " "
       << omg2.degrees() << "\n";
  cout << "BB " << alp1.degrees() << " " << d/Math::degree() << "\n";
  ang bet2a, omg2a, alp2a;
  (void) l.Hybrid(ic, bet2, bet2a, omg2a, alp2a);
  (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
  omg2a -= omg2;
  return omg2a.radians0();
  */
}

Math::real HybridB(const Triaxial& t,
                   const Angle& bet1, const Angle& omg1,
                   const Angle& alp1, const Angle& bet2,
                   Angle& bet2a, Angle& omg2a, Angle& alp2a) {
  typedef Angle ang;
  ang b1{bet1}, o1{omg1}, a1{alp1};
  Triaxial::gamblk gam(t, b1, o1, a1);
  TriaxialLineF l(t, gam, 0.5, 1.5);
  TriaxialLineF::ics ic(l, b1, o1, a1);
  TriaxialLineF::disttx d = l.Hybrid(ic, bet2, bet2a, omg2a, alp2a);
  TriaxialLineG ld(t, gam);
  TriaxialLineG::ics icd(ld, ic);
  return ld.dist(icd, d);
}

// Solve f(alp1) = 0 where alp1 is an azimuth and f(alp1) is the difference in
// lontitude on bet2 and the target longitude.
Angle findroot(const function<Math::real(const Angle&)>& f,
                  Angle xa,  Angle xb,
                  Math::real fa, Math::real fb,
                  int* countn = nullptr, int* countb = nullptr) {
  // Implement root finding method of Chandrupatla (1997)
  // https://doi.org/10.1016/s0965-9978(96)00051-8
  // Here we follow Scherer (2013), Section 6.1.7.3
  // https://doi.org/10.1007/978-3-319-00401-3

  // Here the independent variable is an ang, but the computations on this
  // variable essentially involve its conversion to radians.  There's no need
  // to worry about the angle wrapping around because (xb-xa).radians() is in
  // (0,pi).

  // require xa and xb to be normalized (the result is normalized)
  // require fa and fb to have opposite signs

  typedef Angle ang;
  ang xm;                  // The return value
  int cntn = 0, cntb = 0;
  /*
    cout << xa.degrees() << " " << fa << " "
       << xb.degrees() << " " << fb << "\n";

  int num = 360;
  //  for (int i = 0; i <= num; ++i) {
  { int i = 90;
    Math::real x = 360*i/Math::real(num),
      ff = f(ang::degrees(x));
    cout << "DAT " << x << " " << ff/Math::degree() << "\n";
  }
  return 0;
  */
  bool trip = false;
 for (Math::real t = 1/Math::real(2), ab = 0;
       cntn < 50 || GEOGRAPHICLIB_PANIC;) {
    ang xt = 2*t == 1 ?
      ang::aux(xa.y() + xb.y(),
               xa.x() + xb.x(), true) :
      (2*t < 1 ? xa - ang::radians(t * ab) :
       xb + ang::radians((1 - t) * ab)),
    /*
      ang::radians((1-t)*xa.radians() +
                                    t*xb.radians()),
    */
    /*
    ang xt((1-t) * xa.y() + t * xb.y(),
                (1-t) * xa.x() + t * xb.x(), true),
    */
      xc;
    if (fabs(hypot(xt.x(), xt.y()) - 1) > 0.001) {
      cout << "CHECK " << (2*t - 1) << " "
           << hypot(xa.x(), xa.y()) - 1 << " "
           << hypot(xb.x(), xb.y()) - 1 << " "
           << hypot(xt.x(), xt.y()) - 1 << "\n";
      ang xq(xa - ang::radians(t * ab));
      ang xr(xb + ang::radians((1 - t) * ab));
      cout << "CHECKX "
           << hypot(xq.x(), xq.y()) - 1 << " "
           << hypot(xr.x(), xr.y()) - 1 << "\n";
    }
    if (trip) {
      xm = xt;
      break;
    }
    ++cntn;
    Math::real ft = f(xt), fm, fc;
    cout << cntn << " " << xt.degrees() << " " << ft << " " << ab << " " << t << "\n";
    if (signbit(ft) == signbit(fa)) {
      xc = xa; xa = xt;
      fc = fa; fa = ft;
    } else {
      xc = xb; xb = xa; xa = xt;
      fc = fb; fb = fa; fa = ft;
    }
    if (fabs(fb) < fabs(fa)) {
      xm = xb; fm = fb;
    } else {
      xm = xa; fm = fa;
    }
    // ordering is b - a - c
    ab = (xa-xb).radians0();
    Math::real
      ca = (xc-xa).radians0(),
      cb = ca+ab,
      // Scherer has a fabs(cb).  This should be fabs(ab).
      tl = numeric_limits<Math::real>::epsilon() / fabs(ab);
    // Backward tests to deal with NaNs
    trip =  !(2 * tl < 1 && fabs(fm) > numeric_limits<Math::real>::epsilon());
    // If trip update xm one more time, then Hybrid solution is called once
    // more outside this route to update bet2, omg2, alp2, etc.
    Math::real
      xi = ab / cb,
      phi = (fa-fb) / (fc-fb);
    if ( 2 * tl < 1 && xi / (1 + sqrt(1 - xi)) < phi && phi < sqrt(xi) ) {
      t = fa/(fb-fa) * fc/(fb-fc) - ca/ab * fa/(fc-fa) * fb/(fc-fb);
      // This equation matches the pseudocode in Scherer.  His Eq (6.40) reads
      // t = fa/(fb-fa) * fc/(fb-fc) + ca/cb * fc/(fc-fa) * fb/(fb-fa); this is
      // wrong.
      t = fmin(1 - tl, fmax(tl, t));
    } else {
      t = 1/Math::real(2);
      ++cntb;
    }
  }
  if (countn) *countn += cntn;
  if (countb) *countb += cntb;
  return xm;
}

void HybridTest(Math::real a, Math::real b, Math::real c,
                Math::real bet1d, Math::real omg1d,
                Math::real bet2d, Math::real omg2d) {
  typedef Math::real real;
  typedef Angle ang;
  Triaxial t(a, b, c);
  ang bet1 = ang::degrees(bet1d),
    omg1 = ang::degrees(omg1d),
    bet2 = ang::degrees(bet2d),
    omg2 = ang::degrees(omg2d);
  cout << fixed << setprecision(6);
  /*
  for (int i = -180; i <= 180; ++i) {
  //  {  int i = 101;
    cout << i << " "
         << HybridA(t, bet1, omg1, ang::degrees(i),
         bet2,omg2) / Math::degree() << "\n";
  }
  return;
  */
  if (0) {
  cout << a << " " << b << " " << c << "\n";
  real domg[4];
  ang alp1u[4];
  alp1u[0] = ang::aux( t.kp * omg1.y(), t.k * bet1.x(), true );
  for (unsigned q = 1; q < 4; ++q) {
    alp1u[q] = alp1u[0];
    alp1u[q].setquadrant(q);
  }
  for (unsigned q = 0U; q < 4U; ++q) {
    domg[q] = HybridA(t, bet1, omg1, alp1u[q], bet2, omg2);
    cout << q << " " << alp1u[q].degrees() << " "
         << domg[q]/Math::degree() << "\n";
    if (domg[q] == 0) {
      cout << "Result " << alp1u[q].degrees() << "\n";
      return;
    }
  }
  {
    TriaxialLineF l(t, Triaxial::gamblk{}, 0.5, 1.5);
    TriaxialLineF::ics ic(l, bet1, omg1, alp1u[0]);
    for (unsigned q = 0U; q < 4U; ++q) {
      ic.setquadrant(l, q);
      real dd = Hybrid0(l, ic, bet2, omg2);
      cout << q << " " << dd/Math::degree() << "\n";
    }
  }

  ang xa, xb, xn;
  real fa, fb;
  for (unsigned q = 0U; q < 4U; ++q) {
    if (domg[q] < 0 && domg[(q+1)&3] > 0) {
      xb = alp1u[q]; xa = alp1u[(q+1)&3];
      fb = domg[q]; fa = domg[(q+1)&3];
      break;
    }
  }

  bool bisect = false;
  int countn = 0, countb = 0;
  if (bisect) {
    for (; countb < Math::digits() + 10;) {
      ++countb;
      xn = ang::aux(xb.y() + xa.y(), xb.x() + xa.x(), true);
      real ft =  HybridA(t, bet1, omg1, xn, bet2, omg2);
      cout << countb << " " << xn.degrees() << " " << ft << "\n";
      if (ft == 0 || (xa-xn).radians0() <= 0 || (xn-xb).radians0() <= 0)
        break;
      (ft > 0 ? xa : xb) = xn;
    }
  } else {
    xn = findroot(
                  [&t, &bet1, &omg1, &bet2, &omg2]
                  (const ang& x) -> real
                  { return HybridA(t, bet1, omg1, x, bet2, omg2); },
                  xa,  xb,
                  fa, fb,
                  &countn, &countb);
  }

  std::cout << xn.degrees() << " " << countn << " " << countb << "\n";
  ang alp1(xn), bet2a, alp2a, omg2a;
  real s12 = HybridB(t, bet1, omg1, alp1, bet2, bet2a, omg2a, alp2a);
  std::cout << alp1.degrees() << " " << bet2a.degrees() << " "
            << omg2a.degrees() << " " << alp2a.degrees() << " "
            << s12 << "\n";
  {
    TriaxialLine l0(t, bet1, omg1, alp1);
    ang bb,oo,aa;
    l0.Position(s12, bb,oo,aa);
    std::cout << "TRY3\n" << alp1.degrees() << " " << bb.degrees() << " "
            << oo.degrees() << " " << aa.degrees() << " "
            << s12 << "\n";
  }
  }
  TriaxialLine l = t.Inverse(bet1, omg1, bet2, omg2);
  if (1) {
    real bet1x, omg1x, alp1x, bet2x, omg2x, alp2x;
    l.pos1(bet1x, omg1x, alp1x);
    real s12 = l.Distance();
    l.Position(s12, bet2x, omg2x, alp2x);
    std::cout << bet1x << " " << omg1x << " " << alp1x << " "
              << bet2x << " " << omg2x << " " << alp2x << " "
              << s12 << "\n";
  }
}

std::pair<Math::real, Math::real>
cartdiff(const Triaxial& t,
         Angle bet1, Angle omg1, Angle alp1,
         Angle bet2, Angle omg2, Angle alp2)
{
  Triaxial::vec3 r1, r2, v1, v2;
  t.elliptocart2(bet1, omg1, alp1, r1, v1);
  t.elliptocart2(bet2, omg2, alp2, r2, v2);
  return std::pair<Math::real, Math::real>
    (Triaxial::hypot3(r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]),
     Triaxial::hypot3(v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]));
}

void InverseTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  typedef Angle ang;
  Triaxial t(a, b, c);
  cout << fixed << setprecision(6);
  real bet1d, omg1d, bet2d, omg2d, alp1d, alp2d, s12d;
  while (cin >> bet1d >> omg1d >> bet2d >> omg2d >> alp1d >> alp2d >> s12d) {
    /*
    if (bet1d == 90 && (alp1d > -90 && alp1d <= 90))
      alp1d += alp1d < 0 ? 180 : -180;
    if (bet1d ==-90 &&!(alp1d > -90 && alp1d <= 90))
      alp1d += alp1d < 0 ? 180 : -180;
    if (bet2d == 90 && (alp2d > -90 && alp2d <= 90))
      alp2d += alp2d < 0 ? 180 : -180;
    if (bet2d ==-90 &&!(alp2d > -90 && alp2d <= 90))
      alp2d += alp2d < 0 ? 180 : -180;
    */
    cout << int(bet1d) << " " << int(omg1d) << " "
         << int(bet2d) << " " << int(omg2d) << " " << flush;
    ang
      bet1(ang::degrees(bet1d)),
      omg1(ang::degrees(omg1d)),
      bet2(ang::degrees(bet2d)),
      omg2(ang::degrees(omg2d)),
      alp1(ang::degrees(alp1d)),
      alp2(ang::degrees(alp2d));
    // if (fabs(bet1d) == 90 || fabs(bet2d) == 90) {
    bool umb1 = bet1.x() == 0 && omg1.y() == 0,
      umb2 = bet2.x() == 0 && omg2.y() == 0;
    //    if (umb1 || umb2 || !( (fabs(bet1d) == 90 || fabs(bet2d) == 90) &&
    //                           !(fabs(bet1d) == 90 && fabs(bet2d) == 90) )) {
    if (0) {
    if (umb1 || umb2) {
      cout << "SKIP" << endl;
      continue;
    }
    }

    if (!(
          (fabs(bet1d) == 90 && fabs(bet2d) == 90) ||
          ((omg1d == 0 || omg1d == 180) && (omg2d == 0 || omg2d == 180)) ||
          ((fabs(bet1d) == 90 || omg1d == 0 || omg1d == 180) &&
           (fabs(bet2d) == 90 || omg2d == 0 || omg2d == 180)) ||
          (bet1d == 0 && bet2d == 0) ||
          (umb1 || umb2) ||
          (fabs(bet1d) == 90 || fabs(bet2d) == 90) ||
          ((omg1d == 0 || omg1d == 180) || (omg2d == 0 || omg2d == 180)) ||
          true
          )) {
      cout << "SKIP" << endl;
      continue;
    }
    TriaxialLine l = t.Inverse(bet1, omg1, bet2, omg2);
    ang bet1x, omg1x, alp1x, bet2x, omg2x, alp2x;
    l.pos1(bet1x, omg1x, alp1x);
    /*
    bet1x.rnd(); omg1x.rnd(); alp1x.rnd();
    Triaxial::AngNorm(bet1x, omg1x, alp1x);
    */
    real s12x = l.Distance();
    l.Position(s12x, bet2x, omg2x, alp2x);
    real ds12 = fabs(s12x - s12d);
    std::pair<real, real>
      err1 = cartdiff(t, bet1, omg1, alp1, bet1x, omg1x, alp1x),
      err2 = cartdiff(t, bet2, omg2, alp2, bet2x, omg2x, alp2x),
      err3 = cartdiff(t, bet1x, omg1x, alp1x, bet2x, omg2x, alp2x);
    if (bet1d + bet2d == 0 && err1.second >= 1.e-5) {
      alp1x.reflect(false, true);
      alp2x.reflect(false, true);
      err1 = cartdiff(t, bet1, omg1, alp1, bet1x, omg1x, alp1x);
      err2 = cartdiff(t, bet2, omg2, alp2, bet2x, omg2x, alp2x);
    }
    if (err1.first < 1e-5 &&
        err2.first < 1e-5 &&
        ds12 < 1e-5 &&
        (s12d == 0 ? err3.first < 1.e-5 && err3.second < 1.e-5 :
         err1.second < 1.e-5 && err2.second < 1.e-5)) {
      if (0)
      cout  << alp1d << " " << alp2d << " " << s12d << " OK\n"
            << setprecision(0)
           << bet1x.degrees() << " " << omg1x.degrees() << " "
           << bet2x.degrees() << " " << omg2x.degrees() << " "
            << setprecision(6)
           << alp1x.degrees() << " " << alp2x.degrees() << " "
           << s12x << " OK" << endl;
      else
        cout << "OK" << endl;
    } else {
      if (1)
      cout  << alp1d << " " << alp2d << " " << s12d << " BAD\n"
            << setprecision(1)
           << bet1x.degrees() << " " << omg1x.degrees() << " "
           << bet2x.degrees() << " " << omg2x.degrees() << " "
            << setprecision(6)
           << alp1x.degrees() << " " << alp2x.degrees() << " "
           << s12x << " BAD" << endl;
      else
        cout << "BAD" << endl;
    }
  }
}

int main() {
  try {
    Utility::set_digits();
    typedef Math::real real;
    typedef Angle ang;
   if (0) {
      Triaxial t(sqrt(2.0), 1.0, 1/sqrt(2.0));
      TriaxialLine l = t.Inverse(
                     ang::degrees(-45.0),ang::degrees(-30.0),
                     ang::degrees(-90.0),ang::degrees(0.0));
      real s12 = l.Distance();
      /*      real s12 = 0.361729;
      TriaxialLine l(t, ang::degrees(-45.0),ang::degrees(-30.0),
                      ang::degrees(135.0));
      l.SetDistance(s12);
      */
      cout << s12 << "\n";
      ang bet2, omg2, alp2;
      l.pos1(bet2, omg2, alp2);
        cout << "POS1 " << s12 << " "
             << bet2.degrees() << " "
             << omg2.degrees() << " "
             << alp2.degrees() << "\n";
      cout << fixed << setprecision(6);
      for (int i = -3; i <= 3; ++i) {
        real ss = 0 + i * 1e-6;
        l.Position(ss, bet2, omg2, alp2);
        cout << i << " " << ss << " "
             << bet2.degrees() << " "
             << omg2.degrees() << " "
             << alp2.degrees() << "\n";
      }
      for (int i = -3; i <= 3; ++i) {
        real ss = s12 + i * 1e-6;
        l.Position(ss, bet2, omg2, alp2);
        cout << i << " " << ss << " "
             << bet2.degrees() << " "
             << omg2.degrees() << " "
             << alp2.degrees() << "\n";
      }
      return 0;
    }
    if (0) {
    {
      TriaxialLine ll(Triaxial(sqrt(2.0), 1.0, 1/sqrt(2.0)),
                      ang::degrees(90), ang::degrees(180),
                      ang::degrees(158.253574));
      ll.SetDistance(0.898324);
      ang bet1x, omg1x, alp1x, bet2x, omg2x, alp2x;
      int ibet1, iomg1, ialp1, ibet2, iomg2, ialp2;
      ll.pos1(bet1x, omg1x, alp1x, &ibet1, &iomg1, &ialp1);
      real s12x = ll.Distance();
      ll.Position(s12x, bet2x, omg2x, alp2x, &ibet2, &iomg2, &ialp2);
      std::cout << "OUT\n" << bet1x.degrees() << " " << omg1x.degrees() << " " << alp1x.degrees() << " " << ibet1 << iomg1 << ialp1 << signbit(bet1x.x()) << signbit(omg1x.y()) << "\n"
                << bet2x.degrees() << " " << omg2x.degrees() << " " << alp2x.degrees() << " " << ibet2 << iomg2 << ialp2 << "\n"
                << s12x << "\n";
    }
    return 0;
        {
      TriaxialLine ll(Triaxial(sqrt(2.0), 1.0, 1/sqrt(2.0)),
                      ang::degrees(-90), ang::degrees(0),
                      ang::degrees(-87.089900));
      ll.SetDistance(2.736330);
      ang bet1x, omg1x, alp1x, bet2x, omg2x, alp2x;
      int ibet1, iomg1, ialp1, ibet2, iomg2, ialp2;
      ll.pos1(bet1x, omg1x, alp1x, &ibet1, &iomg1, &ialp1);
      real s12x = ll.Distance();
      ll.Position(s12x, bet2x, omg2x, alp2x, &ibet2, &iomg2, &ialp2);
      std::cout << "OUT\n" << bet1x.degrees() << " " << omg1x.degrees() << " " << alp1x.degrees() << " " << ibet1 << iomg1 << ialp1 << signbit(bet1x.x()) << signbit(omg1x.y()) << "\n"
                << bet2x.degrees() << " " << omg2x.degrees() << " " << alp2x.degrees() << " " << ibet2 << iomg2 << ialp2 << "\n"
                << s12x << "\n";
    }
    return 0;
    }

    real a = 6378172, b = 6378103, c = 6356753;
    if (0) {
      using std::sqrt;
      a = sqrt(real(2)); b = 1; c = 1/a;
      //      real a = 1.01, b = 1, c = 0.8;
      //      a = 6378172; b = 6378103; c = 6356753;
      //      a -= 34; b += 34;
      PositionTest(a, b, c);
      }
    if (1) {
      using std::sqrt;
      a = sqrt(real(2)); b = 1; c = 1/a;
      InverseTest(a, b, c);
    }
    if (0)
      ODEtest(a, b, c);
    /*
    if (0)
      TriaxialTest0();
    */
    if (0) {
      real a = 1.01, b = 1, c = 0.8;
      //    a = sqrt(2.0); c = 1/a;
      //    a = 1.2; c = 0.99;
      a = 6378172; b = 6378103; c = 6356753;
      a -= 34; b += 34;
      TriaxialTest1(a, b, c);
    }
    if (0) {
      real a = 1.01, b = 1, c = 0.8;
      a = 6378172; b = 6378103; c = 6356753;
      a -= 34; b += 34;
      DirectfunTest(a, b, c);
      //      a = sqrt(real(2)); b = 1; c = 1/a;
      //      DirectfunTest(a, b, c);
    }
    if (0) {
      using std::sqrt;
      real bet1 = -45, omg1 = 30, bet2 = 30, omg2 = 60;
      real a = 1.01, b = 1, c = 0.8;
      // s12 = 1.172873690527506
      // alp1 = 29.936197484667744
      //        29.936197484667770
      // alp2 = 26.341535680708272
       a = sqrt(real(2)); b = 1; c = 1/a;
      // s12 = 0.995474462124499
      // alp1 = 30.275215672294671
      //        30.275215672295103
      // alp2 = 48.296592187017524
       omg2 = 60;
      //      a = 6378172; b = 6378103; c = 6356753;
      //      a -= 34; b += 34;
      HybridTest(a, b, c, bet1, omg1, bet2, omg2);
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
}

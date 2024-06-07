#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/AuxAngle.hpp>
#include "Triaxial.hpp"
#include "TriaxialLine.hpp"

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
    if (1) {
      for (int k = 10; k <= kmax; ++k) {
        real eps = pow(real(2), -k);
        int n = t.Direct(r1, v1, s12, r2, v2, eps);
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
        n += t.Direct(r1, v1, s12, r2, v2, eps);
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
  Triaxial t(a, b, c);
  // Circumpolar
  TriaxialLine lca(t, AuxAngle::degrees(real(1)), AuxAngle::degrees(real(0)),
                  AuxAngle::degrees(real(90)));
  if (1) {
  TriaxialLine lcb(t, AuxAngle::degrees(real(89.999)),
                   AuxAngle::degrees(real(0)),
                  AuxAngle::degrees(real(90)));
  }
  // Umbilic
  TriaxialLine lu(t, AuxAngle::degrees(real(90)), AuxAngle::degrees(real(0)),
                  AuxAngle::degrees(real(135)));
  if (1) {
  // Circumpolar
  TriaxialLine ltb(t, AuxAngle::degrees(real(90)),
                   AuxAngle::degrees(real(0.001)),
                   AuxAngle::degrees(real(180)));
  }
  TriaxialLine lta(t, AuxAngle::degrees(real(90)), AuxAngle::degrees(real(89)),
                  AuxAngle::degrees(real(180)));
}

void PositionTest(Math::real a, Math::real b, Math::real c) {
  typedef Math::real real;
  Triaxial t(a, b, c);
  AuxAngle bet1, omg1, alp1, bet2, omg2, alp2;
  bet1 = AuxAngle::degrees(45); omg1 = AuxAngle::degrees(0);
  alp1 = AuxAngle::degrees(90);
  bet1 = AuxAngle::degrees(90); omg1 = AuxAngle::degrees(45);
  alp1 = AuxAngle::degrees(180);
  bet1 = AuxAngle::degrees(90); omg1 = AuxAngle::degrees(0);
  alp1 = AuxAngle::degrees(135);
  //  bet1 = AuxAngle::degrees(1); omg1 = AuxAngle::degrees(2);
  //  alp1 = AuxAngle::degrees(3);
  TriaxialLine l(t, bet1, omg1, alp1);
  cout << fixed << setprecision(14);
  for (int s12 = 0/*-5*/; s12 <= 1/*10*/; ++s12) {
    l.Position(real(s12), bet2, omg2, alp2);
    (void) Triaxial::AngNorm(bet2, omg2, alp2);
    cout << s12 << " "
         << bet2.degrees() << " "
         << omg2.degrees() << " "
         << alp2.degrees() << "\n";
  }
}

Math::real HybridOLD(const Triaxial& t,
                     const AuxAngle& bet1, const AuxAngle& omg1,
                     const AuxAngle& alp1,
                     const AuxAngle& bet2, const AuxAngle& omg2) {
  typedef Math::real real;
  TriaxialLine l(t, bet1, omg1, alp1);
  real s12;
  AuxAngle bet2a, omg2a, alp2a;
  l.Hybrid(bet2, +1, bet2a, omg2a, alp2a, s12);
  (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
  omg2a -= omg2;
  return omg2a.radians();
}

Math::real Hybrid0(const TriaxialLineF& l,
                   const TriaxialLineF::ics ic,
                   const AuxAngle& bet2,
                   const AuxAngle& omg2) {
  AuxAngle bet2a, omg2a, alp2a;
  (void) l.Hybrid(ic, bet2, bet2a, omg2a, alp2a);
  (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
  omg2a -= omg2;
  return omg2a.radians();
}

Math::real HybridA(const Triaxial& t,
                   const AuxAngle& bet1, const AuxAngle& omg1,
                   const AuxAngle& alp1,
                   const AuxAngle& bet2, const AuxAngle& omg2) {
  AuxAngle b1{bet1}, o1{omg1}, a1{alp1};
  Triaxial::gamblk gam(t, b1, o1, a1);
  TriaxialLineF l(t, gam, 0.5, 1.5);
  TriaxialLineF::ics ic(l, b1, o1, a1);
  return Hybrid0(l, ic, bet2, omg2);
  /*
  AuxAngle bet2a, omg2a, alp2a;
  (void) l.Hybrid(ic, bet2, bet2a, omg2a, alp2a);
  (void) Triaxial::AngNorm(bet2a, omg2a, alp2a);
  omg2a -= omg2;
  return omg2a.radians();
  */
}

Math::real HybridB(const Triaxial& t,
                   const AuxAngle& bet1, const AuxAngle& omg1,
                   const AuxAngle& alp1, const AuxAngle& bet2, 
                   AuxAngle& bet2a, AuxAngle& omg2a, AuxAngle& alp2a) {
  AuxAngle b1{bet1}, o1{omg1}, a1{alp1};
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
AuxAngle findroot(const function<Math::real(const AuxAngle&)>& f,
                  AuxAngle xa,  AuxAngle xb,
                  Math::real fa, Math::real fb,
                  int* countn = nullptr, int* countb = nullptr) {
  // Implement root finding method of Chandrupatla (1997)
  // https://doi.org/10.1016/s0965-9978(96)00051-8
  // Here we follow Scherer (2013), Section 6.1.7.3
  // https://doi.org/10.1007/978-3-319-00401-3

  // Here the independent variable is an AuxAngle, but the computations on this
  // variable essentially involve its conversion to radians.  There's no need
  // to worry about the angle wrapping around because (xb-xa).radians() is in
  // (0,pi).

  // require xa and xb to be normalized (the result is normalized)
  // require fa and fb to have opposite signs

  AuxAngle xm;                  // The return value
  int cntn = 0, cntb = 0;
  /*
    cout << xa.degrees() << " " << fa << " "
       << xb.degrees() << " " << fb << "\n";

  int num = 360;
  //  for (int i = 0; i <= num; ++i) {
  { int i = 90;
    Math::real x = 360*i/Math::real(num),
      ff = f(AuxAngle::degrees(x));
    cout << "DAT " << x << " " << ff/Math::degree() << "\n";
  }
  return 0;
  */
  bool trip = false;
 for (Math::real t = 1/Math::real(2), ab = 0;
       cntn < 50 || GEOGRAPHICLIB_PANIC;) {
    AuxAngle xt = 2*t == 1 ?
      AuxAngle(xa.y() + xb.y(),
               xa.x() + xb.x(), true) :
      (2*t < 1 ? xa - AuxAngle::radians(t * ab) :
       xb + AuxAngle((1 - t) * ab)),
    /*
      AuxAngle::radians((1-t)*xa.radians() +
                                    t*xb.radians()),
    */
    /*
    AuxAngle xt((1-t) * xa.y() + t * xb.y(),
                (1-t) * xa.x() + t * xb.x(), true),
    */
      xc;
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
    ab = (xa-xb).radians();
    Math::real
      ca = (xc-xa).radians(),
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
  Triaxial t(a, b, c);
  AuxAngle bet1 = AuxAngle::degrees(bet1d),
    omg1 = AuxAngle::degrees(omg1d),
    bet2 = AuxAngle::degrees(bet2d),
    omg2 = AuxAngle::degrees(omg2d);
  cout << fixed << setprecision(6);
  /*
  for (int i = -180; i <= 180; ++i) {
  //  {  int i = 101;
    cout << i << " "
         << HybridA(t, bet1, omg1, AuxAngle::degrees(i),
         bet2,omg2) / Math::degree() << "\n";
  }
  return;
  */
  cout << a << " " << b << " " << c << "\n";
  Math::real domg[4];
  AuxAngle alp1u[4];
  alp1u[0] = AuxAngle( t.kp * omg1.y(), t.k * bet1.x(), true );
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
      Math::real dd = Hybrid0(l, ic, bet2, omg2);
      cout << q << " " << dd/Math::degree() << "\n";
    }
  }

  AuxAngle xa, xb, xn;
  Math::real fa, fb;
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
      xn = AuxAngle(xb.y() + xa.y(), xb.x() + xa.x(), true);
      Math::real ft =  HybridA(t, bet1, omg1, xn, bet2, omg2);
      cout << countb << " " << xn.degrees() << " " << ft << "\n";
      if (ft == 0 || (xa-xn).radians() <= 0 || (xn-xb).radians() <= 0)
        break;
      (ft > 0 ? xa : xb) = xn;
    }
  } else {
    xn = findroot(
                  [&t, &bet1, &omg1, &bet2, &omg2]
                  (const AuxAngle& x) -> Math::real
                  { return HybridA(t, bet1, omg1, x, bet2, omg2); },
                  xa,  xb,
                  fa, fb,
                  &countn, &countb);
  }

  std::cout << xn.degrees() << " " << countn << " " << countb << "\n";
  AuxAngle alp1(xn), bet2a, alp2a, omg2a;
  Math::real s12 = HybridB(t, bet1, omg1, alp1, bet2, bet2a, omg2a, alp2a);
  std::cout << alp1.degrees() << " " << bet2a.degrees() << " "
            << omg2a.degrees() << " " << alp2a.degrees() << " "
            << s12 << "\n";
}

int main() {
  try {
    Utility::set_digits();
    typedef Math::real real;
    real a = 6378172, b = 6378103, c = 6356753;
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
      real a = 1.01, b = 1, c = 0.8;
      //      a = 6378172; b = 6378103; c = 6356753;
      //      a -= 34; b += 34;
      PositionTest(a, b, c);
    }
    if (1) {
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
       omg2 = 60-180;
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

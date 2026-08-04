#include <iostream>
#include <iomanip>
#include <functional>

#include <GeographicLib/Utility.hpp>
#include <GeographicLib/TransverseMercator.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/AuxLatitude.hpp>
#include <GeographicLib/EllipticFunction.hpp>
#include <GeographicLib/Trigfun.hpp>

using namespace std;
typedef GeographicLib::Math::real real;
typedef GeographicLib::Math::cmplx cmplx;
typedef pair<cmplx, cmplx> cmplx2;

namespace GeographicLib {
  namespace experimental {
  class TM3 {
    typedef Math::real real;
    typedef Math::cmplx cmplx;
    real _a, _f, _e2, _e12, _k0;
    cmplx _e;
    EllipticFunction _ell;
    bool _debug;
  public:
    TM3(real a, real f, real k0 = 1, bool debug = false);
    static cmplx2 TAUTOPHI(cmplx TAU) {
      return { atan(TAU), real(1) / (real(1) + Math::sq(TAU)) };
    }
    static tuple<real, real, int> sincos(real x) {
      real s = sin(x), c = cos(x);
      if (signbit(c)) { s = -s; c = -c; }
      int n = int(rint((x - atan2(s, c)) / Math::pi()));
      return { s, c, n };
    }
    static tuple<cmplx, cmplx, int> sincos(cmplx x) {
      cmplx s = sin(x), c = cos(x);
      if (signbit(c.real())) { s = -s; c = -c; }
      int n = int(rint((x.real() - atan2(s.real(), c.real())) / Math::pi()));
      return { s, c, n };
    }

    template <typename T>
    static T gd(T x) { return atan(sinh(x)); }
    template <typename T>
    static T lam(T phi) { return asinh(tan(phi)); }
    template <typename T>
    static pair<T, T> scaletan(T phi, real m) {
      auto [s, c, n] = sincos(phi);
      return {
        atan(m * s / c) + n * Math::pi(),
        m / (Math::sq(m * s) + Math::sq(c))
      };
    }
    template <typename T>
    pair<T, T>  phitopsi(T phi) const {
      return lam(phi) - _e * atanh(_e * sin(phi));
    }
    real phitopsi(real phi) const {
      return lam(phi) - (signbit(_f) ?
                         -_e.imag() * atan (_e.imag() * sin(phi)) :
                         +_e.real() * atanh(_e.real() * sin(phi)));
    }
    template <typename T>
    T phitochi(T phi) const { return gd(phitopsi(phi)); }
    real chitophi(real chi) const {
      return atan(Math::tauf(tan(chi),
                             (signbit(_f) ? -1 : 1) * abs(_e)));
    }
    real psitophi(real psi) const {
      return atan(Math::tauf(sinh(psi),
                             (signbit(_f) ? -1 : 1) * abs(_e)));
    }
    real psitotau(real psi) const {
      return Math::tauf(sinh(psi), (signbit(_f) ? -1 : 1) * abs(_e));
    }
    real tautopsi(real tau) const {
      return asinh(Math::taupf(tau, (signbit(_f) ? -1 : 1) * abs(_e)));
    }
    cmplx2 BETATOMU(cmplx BETA) const;
    // tau = tan(psi), taup = tan(chi)
    cmplx2 BETATOTAU(cmplx BETA) const;
    cmplx2 TAUTOBETA(cmplx TAU) const;
    cmplx2 TAUTOPSI(cmplx TAU) const;
    cmplx2 PSITOTAU(cmplx PSI) const;
    cmplx2 PHITOPSI(cmplx PHI, int& n) const;
    cmplx2 PSITOGEOD(cmplx PSI) const;
    cmplx2 GEODTOPSI(cmplx GEOD, int& n) const;
    cmplx2 MUTOTM(cmplx MU) const;
    cmplx2 TMTOMU(cmplx TM) const;
    cmplx2 PSITOCHI(cmplx PSI) const;
    cmplx2 CHITOPSI(cmplx CHI, int& n) const;
    cmplx2 PHITOCHI(cmplx PHI) const;
    cmplx2 CHITOPHI(cmplx CHI) const;
    cmplx2 PHITOW(cmplx PHI) const;
    cmplx2 GEODTOCHI(cmplx GEOD) const {
      int n;
      auto [PSI, DPSI] = GEODTOPSI(GEOD, n);
      auto [CHI, DCHI] = PSITOCHI(PSI);
      return { CHI + n * Math::pi(), DCHI*DPSI };
    }
    cmplx2 PHITOBETA(cmplx PHI) const {
      return scaletan(PHI, 1 - _f);
    }
    cmplx2 BETATOPHI(cmplx BETA) const {
      return scaletan(BETA, 1/(1 - _f));
    }
    cmplx2 PHITOMU(cmplx PHI) const {
      auto [BETA, DBETA] = PHITOBETA(PHI);
      auto [MU, DMU] = BETATOMU(BETA);
      return { MU, DMU * DBETA };
    }

    // ZETAP acts like chi
    // ZETA acts like mu
    // ZETA(ZETAP) gives conversion from chi to mu
    // ZETA = betatomu(phitobeta(chitophi(ZETAP)))
    // let PHI = chitophi(ZETAP), ZETAP = phitochi(PHI)
    // ZETA

    // find z = finv(w)
    // ffp(z) returns a pair f(z) and f'(z)
    // z0 = initial guess (typically f(z0) = w0 = real(w))
    // return pair of finv(w) and finv'(w)
    cmplx2
    invertRobust(const function<cmplx2(cmplx)>& ffp,
                 cmplx w, cmplx z0, int& cnt,
                 real wscale = 1, real zscale = 1, real tol = 0) const;
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y, real& gamma, real& k) const;
  };
  }
}

GeographicLib::experimental::TM3::TM3(real a, real f, real k0, bool debug)
  : _a(a)
  , _f(f)
  , _e2(_f * (2 - _f))
  , _e12(_e2 / (1 - _e2))
  , _k0(k0)
  , _e(sqrt(cmplx(_e2)))
  , _ell(-_e12)
  , _debug(debug)
{}

cmplx2 GeographicLib::experimental::TM3::BETATOMU(cmplx BETA) const {
  // _ell.E works over multiple periods
  return {
    Math::pi()/2 * _ell.E(BETA)/_ell.E(),
    Math::pi()/2 * _ell.Delta(sin(BETA), cos(BETA))/_ell.E()
  };
}
cmplx2 GeographicLib::experimental::TM3::BETATOTAU(cmplx BETA) const {
  return { tan(BETA)/(1 - _f), real(1) / (Math::sq(cos(BETA)) * (1 - _f)) };
}
cmplx2 GeographicLib::experimental::TM3::TAUTOBETA(cmplx TAU) const {
  cmplx TBETA = TAU*(1 - _f);
  return { atan(TBETA), (1 - _f) / (real(1) + Math::sq(TBETA)) };
}
cmplx2 GeographicLib::experimental::TM3::TAUTOPSI(cmplx TAU) const {
  cmplx TAU2 = Math::sq(TAU), SECPHI = sqrt(real(1) + TAU2);
  return {
    asinh(TAU) - _e * atanh(_e * TAU / SECPHI),
    (1 - _e2) * SECPHI / (real(1) + (1 - _e2) * TAU2)
  };
}
cmplx2 GeographicLib::experimental::TM3::PSITOTAU(cmplx PSI) const {
  int cnt;
  auto res = invertRobust([this]
                          (cmplx TAU) -> cmplx2
                          {return TAUTOPSI(TAU);},
                          PSI, psitotau(PSI.real()), cnt);
  if (_debug) cerr << "CNT " << cnt << "\n";
  return res;
}
// Assume PHI.real() in [-pi/2, pi/2]
cmplx2 GeographicLib::experimental::TM3::PHITOPSI(cmplx PHI, int& n) const {
  auto [s, c, n0] = sincos(PHI);
  n = n0;
  // psi/dphi = sec(phi) * (1-e^2)/(1-e^2*sin(phi)^2)
  return {
    asinh(s/c) - _e * atanh(_e * s),
    (1 - _e2) / (c * (real(1) - _e2 * Math::sq(s)))
  };
}
cmplx2 GeographicLib::experimental::TM3::PSITOGEOD(cmplx PSI) const {
  real phi = psitophi(PSI.real());
  return {
    { phi, PSI.imag() },
    _a * cos(phi) / sqrt(1 - _e2 * Math::sq(sin(phi)))
  };
}
cmplx2 GeographicLib::experimental::TM3::GEODTOPSI(cmplx GEOD, int& n) const {
  auto [s, c, n0] = sincos(GEOD.real()); n = n0;
  real phi = GEOD.real();
  return {
    { phitopsi(phi), GEOD.imag() },
    sqrt(1 - _e2 * Math::sq(sin(phi))) / (_a * cos(phi))
  };
}

cmplx2 GeographicLib::experimental::TM3::MUTOTM(cmplx MU) const {
  real A = _k0 * (1 - _f) * _a * _ell.E() / (Math::pi() / 2);
  return { MU * A, A };
}

cmplx2 GeographicLib::experimental::TM3::TMTOMU(cmplx TM) const {
  real Ainv = (Math::pi() / 2) / (_k0 * (1 - _f) * _a * _ell.E());
  return { TM * Ainv, Ainv };
}
cmplx2 GeographicLib::experimental::TM3::PSITOCHI(cmplx PSI) const {
  return { gd(PSI), real(1) / cosh(PSI) };
}
cmplx2 GeographicLib::experimental::TM3::CHITOPSI(cmplx CHI, int& n) const {
  auto [s, c, n0] = sincos(CHI); n = n0;
  return { asinh(s/c), real(1) / c };
}
cmplx2 GeographicLib::experimental::TM3::PHITOCHI(cmplx PHI) const {
  int n;
  auto [PSI, DPSI] = PHITOPSI(PHI, n);
  auto [CHI, DCHI] = PSITOCHI(PSI);
  return { CHI + n * Math::pi(), DCHI*DPSI };
}
cmplx2 GeographicLib::experimental::TM3::CHITOPHI(cmplx CHI) const {
  // _f = 0: PHI = CHI
  // _f > 0: singular points
  // A: (phi,lam) = (0, 90*(1-e))
  //    CHI = i*atanh(sin( (1-e) * pi/2 ))
  //    PHI = (0, inf)
  //    dPHI/dCHI = (0, inf) dCHI/dPHI = 0 dPHI/dCHI = inf
  // B: (phi,lam) = 0, 90, PSI = (0,pi/2)
  //    CHI = (0,inf)
  //    PSI = (0,pi/2)
  //    PHI = (pi/2, y) dPHI/dCHI = 0 dCHI/dPHI = inf
  //    where f(y) = log(coth(y/2))-e*atanh(e*cosh(y)) = 0 y in [0, acosh(1/e)]
  //      f'(y) = -1/(2*cosh(y/2)*sinh(y/2))-((e^2*sinh(y))/(1-e^2*cosh(y)^2))

  // 0 89.9999
  // CHI (0,13.95170536896) (0.08952465548519,0)
  // CNT 82
  // PHI (1.57079567315,0.9346436076028) (-1.603729148699e-12,-6.536446975358e-07)
  // CHIX (6.361151710565e-11,13.95170536896) (-3.753596088652,1529883.136456)
  // 0.0001 90
  // CHI (1.57079632674,14.39793040367) (7.668072297593e-12,-0.1398822742018)
  // CNT 58
  // PHI (1.570796326795,0.9346431892297) (4.183726301833e-07,-2.293173405618e-17)
  // CHIX (1.57079632674,14.39793040367) (2390213.718239,0.0001310117855943)
  // _f < 0: singular points
  // A: (phi,lam) = (psitophi(abs(e) * pi/2), 90, )
  //    PSI = (abs(e) * pi/2, pi/2)
  //    CHI = (pi/2, -log(tanh(abs(e) * pi/4)))
  //    PHI = (0, inf)
  // B: [phi,lam] = 0, 90, PSI = [0,pi/2]
  //    CHI = (0,inf)
  //    PHI = (0, y)
  //    where f(y) = asin(tanh(y)) + abs(e)*atanh(abs(e)*sinh(y)) = pi/2
  //          f'(y) = sech(y) + abs(e)^2 * cosh(y)/(1-abs(e)^2*sinh(y)^2)
  //                = sech(y) * (1 + abs(e)^2)/(1 - abs(e)^2*sinh(y)^2)

  static const real maxv = -log(numeric_limits<real>::epsilon());
  /*
  cout << "YYY " << CHI << " "
       << gd(Math::pi()/2*(real(1)-_e)*cmplx(real(0), real(1))) << " "
       << abs(CHI-gd(Math::pi()/2*(real(1)-_e)*cmplx(real(0), real(1))))
       << "\n";
  */

  if (1) {
    if (_f == 0) return {CHI, real(1)};
    int cnt;
    cmplx PHI0 = chitophi(fmin(Math::pi()/2-1/real(100),
                               fmax(1/real(100),CHI.real())));
    cmplx CHIX = gd(Math::pi()/2*(real(1)-_e)*cmplx(real(0), real(1)));
    if (abs(CHI - CHIX) < 10 * numeric_limits<real>::epsilon()) {
      cmplx PHI{Math::pi()/2 * 0, maxv/4};
      PHI0 = PHI;
    } else {
      if (_f > 0 &&
          abs(CHI-cmplx(real(0),atanh(sin( (1-_e.real()) * Math::pi()/2 )))) <
          10 * numeric_limits<real>::epsilon()) {
        if (_debug)
          cerr << "CHIQ " << abs(CHI-cmplx(real(0),atanh(sin( (1-_e.real()) * Math::pi()/2 )))) << "\n";
        cmplx PHI{Math::pi()/2 * 0, maxv/4};
        auto [ignore, DCHIB] = PHITOCHI(PHI);
        int n;
        auto [PSI, DPSI] = PHITOPSI(PHI, n);
        auto [CHIC, DCHIC] = PSITOCHI(PSI);
        if (_debug) {
          cerr << "PHIX " << PHI << " " << real(1)/DCHIB << " " << PSI << " " << DPSI << " " << CHIC << " " << DCHIC << " " << real(1)/(DPSI*DCHIC) << "\n";
          cerr << "PHIY " << PHI0 << " " << PHI << " "
               << chitophi(CHI.real()) << "\n";
        }
        PHI0 = PHI;
        //      return {PHI, real(1)/DCHIB};
      }
      if (CHI.imag() > maxv) {
        if (_f > 0) {
          auto ffp = [e2 = _e2]
            (real y) -> pair<real, real>
            {
              real e = sqrt(e2),
              sy = sinh(y), cy = cosh(y),
              sy2 = sinh(y/2), cy2 = cosh(y/2),
              f = log(cy2/sy2) - e*atanh(e*cy),
              fp = -1/(2*cy2*sy2) - (e2*sy)/(1 - Math::sq(e*cy));
              return {f, fp};
            };
          real ya = 0, yb = acosh(1/_e.real()), y0 = (ya + yb)/2,
            v = Trigfun::root(Trigfun::MERCATOR3,
                              ffp, real(0),
                              y0, ya, yb, 1, 1, -1);
          cmplx PHI{Math::pi()/2, v};
          auto [CHIB, DCHIB] = PHITOCHI(PHI);
          if (_debug) cerr << "VV " << PHI << " " << CHIB << " " << DCHIB << "\n";
          return {PHI, real(1)/DCHIB};
        }
        CHI = cmplx(CHI.real(), copysign(real(maxv), CHI.imag()));
      }
    }
    auto res = invertRobust([this]
                            (cmplx PHI) -> cmplx2
                            {return PHITOCHI(PHI);},
                            CHI, PHI0,
                            cnt);
    if (_debug) cerr << "CNT " << cnt << "\n";
    return res;
  } else {
    int n;
    auto [PSI, DPSI] = CHITOPSI(CHI, n);
    auto [TAU, DTAU] = PSITOTAU(PSI);
    auto [PHI, DPHI] = TAUTOPHI(TAU);
    return { PHI + n * Math::pi(), DPHI*DTAU*DPSI };
  }
}
cmplx2 GeographicLib::experimental::TM3::PHITOW(cmplx PHI) const {
  return {PHI, real(0)};
}

pair<cmplx, cmplx>
GeographicLib::experimental::TM3::invertRobust(const function<pair<cmplx, cmplx>(cmplx)>& ffp,
                  cmplx w, cmplx z0, int& cnt,
                  real wscale, real zscale, real tol) const {
  if (tol <= 0)
    tol = numeric_limits<real>::epsilon();
  real tol2 = sqrt(tol) / 100 * zscale;
  tol *= wscale;

  cmplx z = z0;
  if (_debug) cerr << "X " << w << " " << z0 << "\n";

  const int maxIter = 500;
  cnt = 0;

  auto [f, fp] = ffp(z); f -= w; ++cnt;
  if (_debug) cerr << "A " << f << " " << fp << "\n";

  for (int iter = 0; iter < maxIter; ++iter) {
    real err = abs(f);
    if (_debug) cerr << "E " << f << " " << err << "\n";
    if (!(err > tol))
      break;

    // Line search (adaptive damping)
    //    real lambda = iter < 2 ? 1/real(10) :
    //      iter < 4 ? 1/real(4) :
    //      iter < 8 ? 1/real(2) : 1;

    // Full Newton step
    cmplx dz = -f / fp;
    real lambda = fmin(real(1), zscale / (10 * abs(dz)));
    if (_debug) cerr << "B " << abs(dz) << " " << lambda << "\n";
    if (lambda == 1 && !(abs(dz) > tol + 0*tol2)) {
      z += dz;
      tie(f, fp) = ffp(z); f -= w; ++cnt;
      if (_debug) cerr << "C " << fp << "\n";
      break;
    }

    cmplx znew;

    for (int ls = 0; ls < 10; ++ls) {
      znew = z + lambda * dz;
      tie(f, fp) = ffp(znew); f -= w; ++cnt;
      if (_debug) cerr << "D " << lambda << " " << znew << " " << fp << "\n";
      real errnew = abs(f);

      // Accept if improvement
      if (!(errnew >= err))
        break;

      lambda /= 2;
    }

    z = znew;
  }

  return { z, real(1)/fp };
}

void GeographicLib::experimental::TM3::Forward(real lon0, real lat, real lon,
                  real& x, real& y, real& gamma, real& k) const {
  cmplx SCALE = real(1), DIFF;
  cmplx GEOD = cmplx( lat, lon - lon0 ) * Math::degree();
  if (0) {
    int n;
    cmplx PSI; tie(PSI, DIFF) = GEODTOPSI(GEOD, n); SCALE *= DIFF;
    cmplx TAU; tie(TAU, DIFF) = PSITOTAU(PSI); SCALE *= DIFF;
    if (_debug) cerr << "PSI/TAU " << PSI << " " << TAU << "\n";
    cmplx BETA; tie(BETA, DIFF) = TAUTOBETA(TAU); SCALE *= DIFF;
    BETA += n *Math::pi();
    cmplx MU; tie(MU, DIFF) = BETATOMU(BETA); SCALE *= DIFF;
    cmplx TM; tie(TM, DIFF) = MUTOTM(MU); SCALE *= DIFF;
    SCALE = conj(SCALE);
    x = TM.imag(); y = TM.real();
    gamma = arg(SCALE) / Math::degree(); k = abs(SCALE);
  } else {
    if (_debug) cerr << "GEOD " << GEOD << "\n";
    auto [CHI, DCHI] = GEODTOCHI(GEOD);
    if (_debug) cerr << "CHI " << CHI << " " << DCHI << "\n";
    auto [PHI, DPHI] = CHITOPHI(CHI);
    if (_debug) cerr << "PHI " << PHI << " " << DPHI << "\n";
    auto [CHIX, DCHIX] = PHITOCHI(PHI);
    if (_debug) cerr << "CHIX " << CHIX << " " << DCHIX << "\n";
    auto [MU, DMU] = PHITOMU(PHI);
    if (_debug) cerr << "MU " << MU << " " << DMU << "\n";
    auto [TM, DTM] = MUTOTM(MU);
    SCALE *= DCHI*DPHI*DMU*DTM;
    SCALE = conj(SCALE);
    x = TM.imag(); y = TM.real();
    gamma = arg(SCALE) / Math::degree(); k = abs(SCALE);
  }
}

using namespace GeographicLib;

int main() {
  //echo 45 0 | TransverseMercatorProj -e 5729577.9513082320877 1e-20 -k 1 -p 8
  // 0.00000000 4500000.00000000 0.00000000000000 1.00000000000000
  //echo 45 45 | TransverseMercatorProj -e 5729577.9513082320877 1e-20 -k 1 -p 8
  // 3147292.37309454 5473561.03172453 35.26438968275465 1.15470053837925
  // in radians  0.54930614433405519442 0.95531661812450848082
  // gd(lam(pi/4)+%i*pi/4) ->
  // 0.54930614433405484569*%i + 0.95531661812450927817
  /*
echo 0 0 90 0 | tools/GeodSolve -e 6.4e6 1/5 -E -i -p 14
0.0000000000000000000 0.0000000000000000000 9075733.72447183508203
 6.4b6/9075733.72447183508203b0*90b5;
(%o8)                       6346594.308368390178
echo 0 0 90 0 | tools/GeodSolve -e 6346594.308368390178 1/5 -E -i -p 14
echo 0 0 0 45e5 | tools/GeodSolve -e 6346594.308368390178 1/5 -E -p 14
54.3657985578814606026 0.0000000000000000000 0.0000000000000000000
ckarney@petrel:~/geographiclib/BUILD-mpfr$

  phi = 54.3657985578814606026
  beta = 48.138350966351010137
  theta = 41.758991040149728581
  mu = 45
  chi = 42.260680367805458269
  xi = 45.98909400365

starting point

lam(phi):=asinh(tan(phi));
gd(x):=atan(sinh(x));
psi(phi) := lam(phi) - e * atanh(e*sin(phi)) ;
chi(phi):=gd(psi(phi));


dpsi/dphi = sec(phi) * (1-e^2)/(1-e^2*sin(phi)^2)
chi = gd(psi)
dchi/dpsi = sech(psi)
chi = gd(lam(phi) - e * atanh(e*sin(phi)))
dchi/dphi = sec(phi) * (1-e^2)/(1-e^2*sin(phi)^2) /
(sec(phi) * cosh(e * atanh(e*sin(phi)))
 - tan(phi) * sinh(e * atanh(e*sin(phi))))
= (1-e^2)/(1-e^2*sin(phi)^2) /
(cosh(e * atanh(e*sin(phi)))
 - sin(phi) * sinh(e * atanh(e*sin(phi))))

taup = @(tau) sinh(asinh(tau) - real(e * atanh(e * tau/sqrt(1+tau^2))))
tau -> inf
sig = sinh(e*atanh(e*tau/sqrt(1+tau^2))) -> sinh(e*atanh(e))

taup -> tau*sqrt(1+sig^2) - sig*tau

psi = asinh(tau) - e * atanh(e * tau/sqrt(1+tau^2)); tau = tan(phi)
mu = pi/2 * E(beta, i*e') / E(i*e')
PSI = psi + i*lambda; scale = sec(beta) = sec(phi)*sqrt(1-e^2*sin(phi)^2)
CHI = gd(PSI)
PHI = psif^-1(PSI)
BETA = atan((1-f)*tan(PHI))
MU = pi/2 * E(BETA, i*e') / E(i*e')

pick BETA
MU = pi/2 * E(BETA, i*e') / E(i*e') -- the TM coordinates
TAU = tan(BETA)/(1-f)
PSI = asinh(TAU) - e * atanh(e * TAU/sqrt(1 + TAU^2))
    = psi + i*lambda
dPSI/dTAU = (1-e^2)*sqrt(1+TAU^2)/(1+(1-e^2)*TAU^2)
TM Eq 10
gd(Z) = atan(sinh(x)/cos(y)) + i * asinh(sin(y)/sqrt(sinh(x)^2+cos(y)^2))
[ also write im part as atanh(sin(y)/cosh(x)) ]
gd'(Z) = sech(Z) = 1/cosh(Z)
TM Eq 24
cosh(Z) = exp(i * atan(tanh(x) * tan(y))) * sqrt(sinh(x)^2 + cos(y)^2)

PSI = psi + i*lam
TM Eq 10
gd(PSI) = atan(sinh(psi)/cos(lam)) + i * asinh(sin(lam)/sqrt(sinh(psi)^2+cos(lam)^2))
gd'(PSI) = sech(PSI) = 1/cosh(PSI)
TM Eq 24
cosh(PSI) = epsip(i * atan(tanh(psi) * tan(lam))) * sqrt(sinh(psi)^2 + cos(lam)^2)

TM Eq 18
lam(Z) = =asinh(sin(x)/sqrt(sinh(y)^2+cos(x)^2))+i*atan(sinh(y)/cos(x));
[ also write im part as atanh(sin(x)/cosh(y)) ]
lam'(Z) = sec(Z) = 1/cos(Z)
TM Eq 20
cos(Z) = exp(-i * atan(tan(x) * tanh(y)) * sqrt(cos(x)^2 + sinh(y)^2)

MU = xi+i*eta
TM Eq 18
lam(MU) = =asinh(sin(xi)/sqrt(sinh(eta)^2+cos(xi)^2))+i*atan(sinh(eta)/cos(xi));
lam'(MU) = sec(MU) = 1/cos(MU)
TM Eq 20
cos(MU) = exp(-i * atan(tan(xi) * tanh(eta)) * sqrt(cos(xi)^2 + sinh(eta)^2)

Singular point PSI = 0 + i*(1-e)*pi/2
f = 1/5, e = 3/5, (1-e)*pi/2 = 36
  */
  typedef Math::real real;
  //  typedef Math::cmplx cmplx;
  Utility::set_digits();
  real a = 6.4e6, f = 1/real(50), deg = Math::degree();
  a = real(6400000); a = 1; f = 1/real(5);
  real e2 = 1/real(100);
  f = 0;
  f = 1/real(300);
  f = -1/real(4);
  f = 1/real(5);
  // f = e2 / (1+sqrt(1-e2));
  e2 = f*(2-f);
  experimental::TM3 qq(a, f, 1, true);
  cout << setprecision(13);
  cerr << setprecision(13);
  if (0) {
    real psi = Math::pi()/2*abs(sqrt(cmplx(e2))),
      phi = qq.psitophi(psi);
    cmplx psi2 = qq.phitopsi(phi);
    cout << psi << " " << phi/Math::degree() << " " << psi2 << "\n";
    cout << e2 << " " << abs(sqrt(cmplx(e2))) << "\n";
    return 0;
  }
  TransverseMercator tm(a, f, 1, true && false);
  if (0) {
    real xlim = 5, ylim = Math::pi()/2;
    int xnum = 40, ynum = 10, ndiv = 10;
    for (int ix = 0; ix <= xnum; ++ix) {
      real x = ix * xlim/xnum;
      for (int iy = 0; iy <= ynum*ndiv; ++iy) {
        real y = iy * ylim/(ynum*ndiv);
        cmplx CHI{y, x};
        auto [PHI, ignore] = qq.CHITOPHI(CHI);
        cout << CHI.imag() << " " << CHI.real() << " "
             << PHI.imag() << " " << PHI.real() << "\n";
      }
      cout << "nan nan nan nan\n";
    }
    for (int iy = 0; iy <= ynum; ++iy) {
      real y = iy * ylim/ynum;
      for (int ix = 0; ix <= xnum*ndiv; ++ix) {
        real x = ix * xlim/(xnum*ndiv);
        cmplx CHI{y, x};
        auto [PHI, ignore] = qq.CHITOPHI(CHI);
        cout << CHI.imag() << " " << CHI.real() << " "
             << PHI.imag() << " " << PHI.real() << "\n";
      }
      cout << "nan nan nan nan\n";
    }
    return 0;
  }
  if (1) {
    int num = 10; real rnum = num;
    for (int lat = 0; lat <= 80; lat += 10) {
      for (int lon = 0; lon <= 90*num; ++lon) {
        real x, y, gamma, k;
        qq.Forward(0, lat, lon/rnum, x, y, gamma, k);
        cout << lat << " " << lon/rnum << " "
             << x << " " << y << " " << gamma << " " << k << "\n";
      }
      cout << "nan nan nan nan nan nan\n";
    }
    for (int lon = 0; lon <= 90; lon += 10) {
      for (int lat = 0; lat <= 80*num; ++lat) {
        real x, y, gamma, k;
        qq.Forward(0, lat/rnum, lon, x, y, gamma, k);
        cout << lat/rnum << " " << lon << " "
             << x << " " << y << " " << gamma << " " << k << "\n";
      }
      cout << "nan nan nan nan nan nan\n";
    }
    for (int lat = 1; lat <= 9; ++lat) {
      for (int lon = 80*num; lon <= 90*num; ++lon) {
        real x, y, gamma, k;
        qq.Forward(0, lat, lon/rnum, x, y, gamma, k);
        cout << lat << " " << lon/rnum << " "
             << x << " " << y << " " << gamma << " " << k << "\n";
      }
      cout << "nan nan nan nan nan nan\n";
    }
    for (int lon = 81; lon <= 89; ++lon) {
      for (int lat = 0; lat <= 10*num; ++lat) {
        real x, y, gamma, k;
        qq.Forward(0, lat/rnum, lon, x, y, gamma, k);
        cout << lat/rnum << " " << lon << " "
             << x << " " << y << " " << gamma << " " << k << "\n";
      }
      cout << "nan nan nan nan nan nan\n";
    }
    return 0;
  }
  if (1) {
    cout << setprecision(19);
    cout << "flat " << f << " " << sqrt(cmplx(f*(2-f))) << " "
         << qq.psitophi(sqrt(abs(f*(2-f))) * Math::pi()/2)/Math::degree()
         << "\n";
    real lat, lon;
    while (cin >> lat >> lon) {
      real x, y, gamma, k;
      qq.Forward(0, lat, lon, x, y, gamma, k);
      cout << x << " " << y << " " << gamma << " " << k << "\n";
      tm.Forward(0, lat, lon, x, y, gamma, k);
      cout << x << " " << y << " " << gamma << " " << k << "\n";

    }
  } else {
    real betar, betai;
    while (cin >> betar >> betai) {
      cmplx SCALE = real(1), DIFF,
        BETA = cmplx(betar, betai) * deg;
      cout << "BETA " << BETA / deg << "\n";
      cmplx TAU; tie(TAU, DIFF) = qq.BETATOTAU(BETA); SCALE /= DIFF;
      cout << "TAU " << TAU << " " << DIFF << "\n";
      cout << "PHI " << atan(TAU) / deg << "\n";
      cmplx PSI; tie(PSI,DIFF) = qq.TAUTOPSI(TAU); SCALE /= DIFF;
      cout << "PSI " << PSI << " " << DIFF << "\n";
      cmplx GEOD; tie(GEOD, DIFF) = qq.PSITOGEOD(PSI); SCALE /= DIFF;
      cout << "GEOD " << GEOD /deg << " " << DIFF << "\n";
      cmplx MU; tie(MU, DIFF) = qq.BETATOMU(BETA); SCALE *= DIFF;
      cout << "MU " << MU /deg << " " << DIFF << "\n";
      cmplx TM; tie(TM, DIFF) = qq.MUTOTM(MU); SCALE *= DIFF;
      SCALE = conj(SCALE);
      GEOD /= deg;
      cout << TM.imag() << " " << TM.real() << " "
           << arg(SCALE) / deg << " " << abs(SCALE) << "\n";
      real x, y, gamma, k;
      tm.Forward(0, GEOD.real(), GEOD.imag(), x, y, gamma, k);
      cout << x << " " << y << " "
           << gamma << " " << k << "\n";
    }
  }
}
/*
singularity at BETA = atan((1-f)*%i) = (0, atanh(1-f))
  tan(BETA) = (1-f)*%i
  TAU = tan(PHI) = %i
  PHI = inf * %i
  SIN(PHI) = inf * %i
  atanh(e*SIN(PHI)) = pi/2 * %i
  SECPHI = 0
  PSI = asinh(TAU) - _e * atanh(_e * TAU / SECPHI)  = %pi/2*(1-e)*%i
  for e real => phi = 0, lambda = (1-e)*%pi/2
  for e imag => lambda = %pi/2, psi = %pi*abs(e)/2, phi = ...
abs(e) = 3/4, phi = 42.68333708561267059 deg

dPSI/dTAU = 0 (prop to sec(PHI) = sech(inf) = 0)
dMU/dBETA = sqrt(1 - k2 * sin(BETA)^2) = 0
k2 = -e2/(1+e2)

f = 1/5, e2 = 9/25, e12 = 9/16, e = 3/5, e' = 3/4
f = -1/4, e2 = -9/16, e12 = -9/25, e = i*3/4, e' = i*3/5

f = 1/5, e2 = 9/25, e12 = 9/16, e = 3/5, e' = 3/4
                        BETA = %i * atanh(4/5) = %i * 62.945847461890764 deg
f = -1/4, e2 = -9/16, e12 = -9/25, e = i*3/4, e' = i*3/5
                        BETA = %i * atanh(5/4) = 90deg + %i*atanh(4/5)
                        = 90deg + %i * 62.945847461890764 deg

BETA (90,62.94584746189)
TAU (2.755455298082e-17,1) (-0.45,6.888638245204e-17)
PHI (86.46303131836,1052.21291249)
SECPHI (1.305049717425e-09,2.111379559944e-08)
PSI (1.178097245096,1.570796326795) (-3.625138103958e-09,-5.864943222067e-08)
GEOD (30.99847986629,90) (5117467.24688,0)
MU (89.99999977109,45.62901999621) (1.26075147151e-08,-4.767355851386e-09)
5751622.716586 11344667.12674 107.1764438385 0.7193869519694

dPSI/dTAU = 0 (prop to sec(PHI) = sech(inf) = 0)
dMU/dBETA = sqrt(1 - k2 * sin(BETA)^2) = 0
k2 = -e2/(1+e2)

f = 1/5, e2 = 9/25, e12 = 9/16, e = 3/5, e' = 3/4
                        BETA = %i * atanh(4/5) = %i * 62.945847461890764 deg
f = -1/4, e2 = -9/16, e12 = -9/25, e = i*3/4, e' = i*3/5
                        BETA = %i * atanh(5/4) = 90deg + %i*atanh(4/5)
                        = 90deg + %i * 62.945847461890764 deg

tan(BETA)=5/4* %i
tan(PHI) = %i
PSI = pi/2*abs(_e) + %i*pi/2
*/

#include <iostream>
#include <iomanip>
#include <functional>

#include <GeographicLib/Utility.hpp>
#include <GeographicLib/TransverseMercator.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/AuxLatitude.hpp>
#include <GeographicLib/EllipticFunction.hpp>

using namespace GeographicLib;
using namespace std;
typedef Math::real real;
typedef Math::cmplx cmplx;
typedef pair<cmplx, cmplx> cmplx2;

template <typename T>
T gd(T x) { return atan(sinh(x)); }
template <typename T>
T lam(T phi) { return asinh(tan(phi)); }

class TM3 {
  typedef Math::real real;
  typedef Math::cmplx cmplx;
  real _a, _f, _e2, _e12;
  cmplx _e;
  EllipticFunction _ell;
public:
  TM3(real a, real f);
  template <typename T>
  static T gd(T x) { return atan(sinh(x)); }
  template <typename T>
  static T lam(T phi) { return asinh(tan(phi)); }
  template <typename T>
  static T scaletan(T phi, real fm1) { return atan(fm1 * tan(phi)); }
  template <typename T>

  cmplx phitopsi(T phi) {
    return lam(phi) - _e * atanh(_e * sin(phi));
  }
  real phitopsi(real phi) {
    return lam(phi) - (signbit(_f) ?
                          -_e.imag() * atan (_e.imag() * sin(phi)) :
                          +_e.real() * atanh(_e.real() * sin(phi)));
  }
  template <typename T>
  T phitochi(T phi) { return gd(phitopsi(phi)); }
  real chitophi(real chi) {
    return atan(Math::tauf(tan(chi),
                           (signbit(_f) ? -1 : 1) * abs(_e)));
  }
  real psitophi(real psi) {
    return atan(Math::tauf(sinh(psi),
                           (signbit(_f) ? -1 : 1) * abs(_e)));
  }
  // TM paper's complex chi
  cmplx CHI(real phi, real lambda)
  { return cmplx(phitopsi(phi), lambda); }
  // TM paper's complex zeta'
  cmplx ZETAP(real phi, real lambda)
  { return gd(CHI(phi, lambda)); }
  cmplx2 BETATOMU(cmplx BETA);
  // tau = tan(psi), taup = tan(chi)
  cmplx2 BETATOTAU(cmplx BETA);
  cmplx2 TAUTOPSI(cmplx TAU);
  cmplx2 PSITOGEOD(cmplx PSI);
  cmplx2 GEODTOPSI(cmplx GEOD);
  cmplx2 MUTOTM(cmplx MU);
  // ZETAP acts like chi
  // ZETA acts like mu
  // ZETA(ZETAP) gives conversion from chi to mu
  // ZETA = betatomu(phitobeta(chitophi(ZETAP)))
  // let PHI = chitophi(ZETAP), ZETAP = phitochi(PHI)
  // ZETA
  // Solve z = f(x)
  static cmplx
  root(const function<cmplx2(cmplx)>& ffp, // f(x), f'(x)
       cmplx z,
       const function<real(real)>& finvr, // f^-1(x) for real x
       real tol2 = 0);

  static bool newtonSolve(const function<cmplx(cmplx)>& f,
                          const function<cmplx(cmplx)>& fp,
                          cmplx w,
                          cmplx& z,
                          real tol,
                          int maxIter = 50);
  static cmplx invertRecursive(const function<cmplx(cmplx)>& f,
                               const function<cmplx(cmplx)>& fp,
                               const function<real(real)>& gr,
                               cmplx w,
                               cmplx z0, // carry the current best guess
                               real tol,
                               int depth = 0,
                               int maxDepth = 20);
  static cmplx invert(const function<cmplx(cmplx)>& f,
                      const function<cmplx(cmplx)>& fp,
                      const function<real(real)>& gr,
                      cmplx w,
                      real tol);
  static cmplx invertAdaptive(const function<cmplx(cmplx)>& f,
                              const function<cmplx(cmplx)>& fp,
                              const function<real(real)>& gr,
                              cmplx w_target,
                              real tol);
  static cmplx invertRobust(const function<cmplx(cmplx)>& f,
                            const function<cmplx(cmplx)>& fp,
                            const function<real(real)>& gr,
                            cmplx w,
                            real tol);
};

TM3::TM3(real a, real f)
  : _a(a)
  , _f(f)
  , _e2(_f * (2 - _f))
  , _e12(_e2 / (1 - _e2))
  , _e(sqrt(cmplx(_e2)))
  , _ell(-_e12)
{}

cmplx2 TM3::BETATOMU(cmplx BETA) {
  return cmplx2(Math::pi()/2 * _ell.E(BETA)/_ell.E(),
                Math::pi()/2 * _ell.Delta(sin(BETA), cos(BETA))/_ell.E());
}
cmplx2 TM3::BETATOTAU(cmplx BETA) {
  return cmplx2(tan(BETA)/(1 - _f),
                real(1) / (Math::sq(cos(BETA)) * (1 - _f)));
}
cmplx2 TM3::TAUTOPSI(cmplx TAU) {
  cmplx TAU2 = Math::sq(TAU), SECPHI = sqrt(real(1) + TAU2);
  return cmplx2(asinh(TAU) - _e * atanh(_e * TAU / SECPHI),
                (1 - _e2) * SECPHI / (real(1) + (1 - _e2) * TAU2));
}

cmplx2 TM3::PSITOGEOD(cmplx PSI) {
  real phi = psitophi(PSI.real());
  return cmplx2(cmplx(phi, PSI.imag()),
                cmplx(_a * cos(phi) / sqrt(1 - _e2 * Math::sq(sin(phi)))));
}
cmplx2 TM3::GEODTOPSI(cmplx GEOD) {
  real phi = GEOD.real();
  return cmplx2(cmplx(phitopsi(GEOD.real()), GEOD.imag()),
                sqrt(1 - _e2 * Math::sq(sin(phi))) / (_a * cos(phi)));
}

cmplx2 TM3::MUTOTM(cmplx MU) {
  real A = (1 - _f) * _a * _ell.E() / (Math::pi() / 2);
  return cmplx2(MU * A, cmplx(A));
}

cmplx
TM3::root(const function<cmplx2(cmplx)>& ffp,
          cmplx z,
          const function<real(real)>& finvr, // f^-1(z) for real z
          real tol2) {
  static const int maxit = 20;
  if (tol2 <= 0)
    tol2 = numeric_limits<real>::epsilon();
  cmplx x = cmplx(finvr(z.real()));
  int i = 0;
  for (; i < maxit; ++i) {
    auto [v, vp] = ffp(x);
    v -= z;
    cmplx dx = - v/vp;
    x += dx;
    if (!(norm(dx) > tol2))
      break;
  }
  return x;
}

bool TM3::newtonSolve(const function<cmplx(cmplx)>& f,
                      const function<cmplx(cmplx)>& fp,
                      cmplx w,
                      cmplx& z,
                      real tol,
                      int maxIter) {
  for (int i = 0; i < maxIter; ++i) {
    cmplx F = f(z) - w;
    if (abs(F) < tol)
      return true;

    cmplx dF = fp(z);
    if (abs(dF) < 1e-14)
      return false;

    cmplx dz = F / dF;
    z -= dz;

    if (abs(dz) < tol)
      return true;
  }
  return false;
}
cmplx TM3::invertRecursive(const function<cmplx(cmplx)>& f,
                           const function<cmplx(cmplx)>& fp,
                           const function<real(real)>& gr,
                           cmplx w,
                           cmplx z0,
                           real tol,
                           int depth,
                           int maxDepth) {
  if (depth > maxDepth)
    throw runtime_error("Max recursion depth exceeded");

  cmplx z = z0;

  // Try Newton directly
  if (newtonSolve(f, fp, w, z, tol))
    return z;

  // If imaginary part is tiny, we are stuck
  if (abs(w.imag()) < tol)
    throw runtime_error("Newton failed near real axis");

  // First midpoint
  cmplx w_mid(w.real(), w.imag() / 2);

  // Recurse toward midpoint using current guess
  cmplx z_mid = invertRecursive(f, fp, gr,
                                w_mid, z0, tol, depth + 1, maxDepth);

  // Try again from midpoint solution
  cmplx z_try = z_mid;
  if (newtonSolve(f, fp, w, z_try, tol))
    return z_try;

  // Second midpoint (closer to w), using improved guess
  cmplx w_mid2 = (w_mid + w) / real(2);

  cmplx z_mid2 = invertRecursive(f, fp, gr,
                                 w_mid2, z_mid, tol, depth + 1, maxDepth);

  z_try = z_mid2;
  if (newtonSolve(f, fp, w, z_try, tol))
    return z_try;

  throw runtime_error("Failed to converge");
}

cmplx TM3::invert(const function<cmplx(cmplx)>& f,
                  const function<cmplx(cmplx)>& fp,
                  const function<real(real)>& gr,
                  cmplx w,
                  real tol) {
  // Initial guess ONLY once
  cmplx z0(gr(w.real()), 0.0);

  return invertRecursive(f, fp, gr, w, z0, tol, 0, 20);
}

// Adaptive continuation solver

cmplx TM3::invertAdaptive(const function<cmplx(cmplx)>& f,
                     const function<cmplx(cmplx)>& fp,
                     const function<real(real)>& gr,
                     cmplx w_target,
                     real tol) {
  // Initial anchor on real axis
  cmplx w0(w_target.real(), 0.0);
  cmplx z = cmplx(gr(w0.real()), 0.0);

  cmplx direction = w_target - w0;

  real t = 0;
  real dt = 1/real(10);

  const real dt_min = real(1e-6);
  const real dt_max = 1/real(2);

  while (t < 1) {
    dt = fmin(dt, 1 - t);

    // Predictor step
    cmplx w_dot = direction;
    cmplx dzdt = w_dot / fp(z);   // dz/dt

    cmplx z_pred = z + dt * dzdt;
    cmplx w_next = w0 + (t + dt) * direction;

    cmplx z_corr = z_pred;

    // Corrector (Newton)
    bool ok = newtonSolve(f, fp, w_next, z_corr, tol);

    if (ok) {
      // Accept step
      z = z_corr;
      t += dt;

      // Adapt step size based on correction size
      real correction = abs(z_corr - z_pred);

      if (correction < tol / 10)
        dt = fmin(dt * 2, dt_max);
      else if (correction > 10 * tol)
        dt = fmax(dt / 2, dt_min);
    } else {
      // Reject step → shrink
      dt /= 2;
      if (dt < dt_min)
        throw runtime_error("Step size underflow");
    }
  }

  return z;
}

cmplx TM3::invertRobust(const function<cmplx(cmplx)>& f,
                        const function<cmplx(cmplx)>& fp,
                        const function<real(real)>& gr,
                        cmplx w,
                        real tol) {
  // Initial guess from real axis
  cmplx z(gr(w.real()));

  const int maxIter = 100;

  for (int iter = 0; iter < maxIter; ++iter) {
    cmplx F = f(z) - w;
    real err = abs(F);

    if (err < tol)
      return z;

    cmplx dF = fp(z);
    if (abs(dF) < 1e-14)
      throw runtime_error("Derivative too small");

    // Full Newton step
    cmplx step = -F / dF;

    // Line search (adaptive damping)
    real lambda = 1;
    cmplx z_trial;
    real err_trial;

    for (int ls = 0; ls < 10; ++ls) {
      z_trial = z + lambda * step;
      err_trial = abs(f(z_trial) - w);

      // Accept if improvement
      if (err_trial < err)
        break;

      lambda /= 2;
    }

    if (lambda < 1e-6)
      throw runtime_error("Line search failed");

    z = z_trial;
  }

  throw runtime_error("Max iterations exceeded");
}

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

psi = lam(phi) - e * atanh(e*sin(phi)) = psif(phi)
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
  a = real(6400000); f = 1/real(5);
  f = 1/real(5);
  TM3 qq(a, f);
  TransverseMercator tm(a, f, 1, f != 0);
  Geodesic geod(a, f, true);
  AuxLatitude aux(a, f);
  real betar, betai;
  cout << setprecision(13);
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
/*
singularity at BETA = atan((1-f)*%i) = (0, atanh(1-f))
  tan(BETA) = (1-f)*%i
  TAU = tan(PHI) = %i
  PHI = inf * %i
  SIN(PHI) = inf * %i
  atanh(e*SIN(PHI)) = pi/2 * %i
  SECPHI = 0
  PSI = asinh(TAU) - _e * atanh(_e * TAU / SECPHI)  = %pi/2*(1-e)*%i

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

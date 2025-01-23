#include <GeographicLib/Constants.hpp>

#include <array>
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

int main() {
  typedef GeographicLib::Math::real REAL;
  typedef std::array<REAL, 1> REAL1;
  //  typedef boost::numeric::odeint::bulirsch_stoer<REAL1, REAL> step;
  typedef boost::numeric::odeint::bulirsch_stoer_dense_out<REAL1, REAL> dstep;
  REAL eps = REAL(1e-15);
  dstep ds(eps);
  //  step s(eps);
  auto fun = [](const REAL1& /*y*/, REAL1& yp, REAL t) -> void {
    yp[0] = t;
  };
  REAL t0 = 0, t1 = 50;
  REAL1 y0, y1;
  y0[0] = 0;
  ds.initialize(y0, t0, 1/REAL(4));
  (void) ds.do_step(fun);
  while (ds.current_time() < t1) {
    std::cout << ds.current_time() << "\n";
    (void) ds.do_step(fun);
  }
  std::cout << ds.current_time() << "\n";
  ds.calc_state(t1, y1);
  std::cout << y1[0] << "\n";
}

#include <iostream>
#include <iomanip>
#include <GeographicLib/EllipticFunction.hpp>

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions
#  pragma warning (disable: 4127)
#endif

using namespace GeographicLib;
int main() {
  typedef GeographicLib::Math::real real;
  try {
    if (true) {
      real b = 1, a = 10,
        e2 = (a*a - b*b)/(a*a), ep2 = (a*a - b*b)/(b*b);
      EllipticFunction elle(e2, e2);
      EllipticFunction ellep(-ep2, -ep2);
      std::cout << std::fixed << std::setprecision(8);
      for (int i = 0; i <= 90; i += 5) {
        real
          beta = i * Math::degree<real>(),
          phi = atan(sqrt(1 + ep2) * tan(beta)),
          u = elle.F(phi),
          y = b * ellep.E(beta),
          M = a*elle.E();
        std::cout << (y / M) * 90 << " "
                  << i << " "
                  << phi / Math::degree<real>() << " "
                  << (u / elle.K()) * 90 << "\n";
      }
      /* Create plot with
t=load('data.txt');
plot(t(:,1),t(:,1),'-k',...
     t(:,1),t(:,2),'-kx',...
     t(:,1),t(:,4),'-ko',...
     t(:,1),t(:,3),'-k+')
axis equal;axis([0,90,0,90]);
title('meridian measures for a/b = 10');
xlabel('meridian distance (\circ)');
ylabel('measure (\circ)');
legend('distance',...
       'parametric latitude',...
       'elliptic variable',...
       'geographic latitude',...
       'Location','SouthEast');
set(gcf,'PaperSize',[4.5,4.5]);
set(gcf,'PaperPosition',[0,0,4.5,4.5]);
print meridian-measures.png -dpng

       */
      return 0;
    }

    if (false) {
      real alpha2 = 0.8, k2 = -0.4;
      EllipticFunction ellG(k2,alpha2);
      EllipticFunction ellH(k2,k2/alpha2);
      
      std::cout << std::setprecision(10);
      for (int i = -179; i <= 180; i += 10) {
        real
          phi = i * Math::degree<real>(),
          sn = sin(phi), cn = cos(phi), dn = ellG.Delta(sn, cn),
          g = ellG.G(phi),
          h = (k2/alpha2)*ellH.H(phi) + sqrt(1-k2/alpha2)/sqrt(1-alpha2)*
          atan2(sqrt(1-alpha2)*sqrt(1-k2/alpha2)*sn, dn*cn);
        
        std::cout << i << " " << g << " " << h << " " << h-g << "\n";
      }
      return 0;
    }
    if (false) {
      // For tabulated values in A+S
      real
        ASalpha = 30*Math::degree<double>(),
        k2 = Math::sq(std::sin(ASalpha)),
        alpha2 = 0.3;

      EllipticFunction ell(k2, alpha2);
      real dphi = Math::degree<real>();
      std::cout << std::fixed << std::setprecision(10);
      for (int i = 0; i <= 90; i += 15) {
        real phi = i * dphi;
        std::cout << i << " "
                  << ell.F(phi) << " "
                  << ell.E(phi) << " "
                  << ell.D(phi) << " "
                  << ell.Pi(phi) << " "
                  << ell.G(phi) << " "
                  << ell.H(phi) << "\n";
      }
      return 0;
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
  return 0;
}

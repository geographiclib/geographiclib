// Example of using the GeographicLib::EllipticFunction class

#include <iostream>
#include <iomanip>
#include <exception>
#include <cmath>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/EllipticFunction.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  using cmplx = Math::cmplx;
  try {
    EllipticFunction ell(0.1);  // parameter m = 0.1
    // See Abramowitz and Stegun, table 17.1
    cout << ell.K() << " " << ell.E() << "\n";
    double phi = 20, sn, cn;
    Math::sincosd(phi, sn ,cn);
    // See Abramowitz and Stegun, table 17.6 with
    // alpha = asin(sqrt(m)) = 18.43 deg and phi = 20 deg
    cout << ell.E(phi * Math::degree()) << " "
         << ell.E(sn, cn, ell.Delta(sn, cn))
         << "\n\n";
    // See Carlson 1995, Sec 3.
    cout << setprecision(14)
         << "RF(1,2,0)      = " << EllipticFunction::RF(1,2)         << " = "
         << EllipticFunction::RF(cmplx(1),2)                         << "\n"
         << "RF(0.5,1)      = " << EllipticFunction::RF(0.5,1)       << " = "
         << EllipticFunction::RF(cmplx(0.5),1)                       << "\n"
         << "RF(i,-i)       = "
         << EllipticFunction::RF(cmplx(0,1),cmplx(0,-1))             << "\n"
         << "RF(i-1,i)      = "
         << EllipticFunction::RF(cmplx(-1,1),cmplx(0,1))             << "\n"
         << "RF(2,3,4)      = " << EllipticFunction::RF(2,3,4)       << " = "
         << EllipticFunction::RF(cmplx(2),3,4)                       << "\n"
         << "RF(i,-i,2)     = "
         << EllipticFunction::RF(cmplx(0,1),cmplx(0,-1),2)           << "\n"
         << "RF(i-1,i,1-i)  = "
         << EllipticFunction::RF(cmplx(-1,1),cmplx(0,1),cmplx(1,-1)) << "\n\n"
         << "RC(0,1/4)      = " << EllipticFunction::RC(0,0.25)      << " = "
         << EllipticFunction::RC(cmplx(0),0.25)                      << "\n"
         << "RC(9/4,2)      = " << EllipticFunction::RC(2.25,2)      << " = "
         << EllipticFunction::RC(cmplx(2.25),2)                      << "\n"
         << "RC(0,i)        = "
         << EllipticFunction::RC(0,cmplx(0,1))                       << "\n"
         << "RC(-i,i)       = "
         << EllipticFunction::RC(cmplx(0,-1),cmplx(0,1))             << "\n"
         << "RC(1/4,-2)     = " << EllipticFunction::RC(0.25,-2)     << " = "
         << EllipticFunction::RC(cmplx(0.25),-2)                     << "\n"
         << "RC(i,-1)       = "
         << EllipticFunction::RC(cmplx(0,1),-1)                      << "\n\n"
         << "RJ(0,1,2,3)    = " << EllipticFunction::RJ(0,1,2,3)     << " = "
         << EllipticFunction::RJ(cmplx(0),1,2,3)                     << "\n"
         << "RJ(2,3,4,5)    = " << EllipticFunction::RJ(2,3,4,5)     << " = "
         << EllipticFunction::RJ(cmplx(2),3,4,5)                     << "\n"
         << "RJ(2,3,4,-1+i) = "
         << EllipticFunction::RJ(2,3,4,cmplx(-1,1))                  << "\n"
         << "RJ(i,-i,0,2)   = "
         << EllipticFunction::RJ(cmplx(0,1),cmplx(0,-1),0,2)         << "\n"
         << "RJ(-1+i,-1-i,1,2)    = "
         << EllipticFunction::RJ(cmplx(-1,1),cmplx(-1,-1),1,2)       << "\n"
         << "RJ(i,-i,0,1-i)       = "
         << EllipticFunction::RJ(cmplx(0,1),cmplx(0,-1),0,cmplx(1,-1))  << "\n"
         << "RJ(-1+i,-1-i,1,-3+i) = "
         << EllipticFunction::RJ(cmplx(-1,1),cmplx(-1,-1),1,cmplx(-3,1))<< "\n"
         << "RJ(-1+i,-2-i,-i,-1+i)= "
         << EllipticFunction::RJ(cmplx(-1,1),cmplx(-2,-1),
                                 cmplx(0,-1),cmplx(-1,1))            << "\n"
         << "RJ(2,3,4,-1/2) = " << EllipticFunction::RJ(2,3,4,-0.5)  << " = "
         << EllipticFunction::RJ(cmplx(2),3,4,-0.5)                  << "\n"
         << "RJ(3,4,2,-1/2) = " << EllipticFunction::RJ(3,4,2,-0.5)  << " = "
         << EllipticFunction::RJ(cmplx(3),4,2,-0.5)                  << "\n"
         << "RJ(4,2,3,-1/2) = " << EllipticFunction::RJ(4,2,3,-0.5)  << " = "
         << EllipticFunction::RJ(cmplx(4),2,3,-0.5)                  << "\n"
         << "RJ(2,3,4,-5)   = " << EllipticFunction::RJ(2,3,4,-5)    << " = "
         << EllipticFunction::RJ(cmplx(2),3,4,-5)                    << "\n\n"
         << "RD(0,2,1)      = " << EllipticFunction::RD(0,2,1)       << " = "
         << EllipticFunction::RD(cmplx(0),2,1)                       << "\n"
         << "RD(2,3,4)      = " << EllipticFunction::RD(2,3,4)       << " = "
         << EllipticFunction::RD(cmplx(2),3,4)                       << "\n"
         << "RD(i,-i,2)     = "
         << EllipticFunction::RD(cmplx(0,1),cmplx(0,-1),2)           << "\n"
         << "RD(0,i,-i)     = "
         << EllipticFunction::RD(0,cmplx(0,1),cmplx(0,-1))           << "\n"
         << "RD(0,i-1,i)    = "
         << EllipticFunction::RD(0,cmplx(-1,1),cmplx(0,1))           << "\n"
         << "RD(-2-i,-i,-1+i)= "
         << EllipticFunction::RD(cmplx(-2,-1),cmplx(0,-1),cmplx(-1,1)) << "\n\n"
         << "RG(0,16,16)    = " << EllipticFunction::RG(16,16)       << " = "
         << EllipticFunction::RG(cmplx(16),16)                       << "\n"
         << "RG(2,3,4)      = " << EllipticFunction::RG(2,3,4)       << " = "
         << EllipticFunction::RG(cmplx(2),3,4)                       << "\n"
         << "RG(0,i,-i)     = "
         << EllipticFunction::RG(0,cmplx(0,1),cmplx(0,-1))           << "\n"
         << "RG(i-1,i,0)    = "
         << EllipticFunction::RG(cmplx(-1,1),cmplx(0,1))             << "\n"
         << "RG(-i,i-1,i)   = "
         << EllipticFunction::RG(cmplx(0,-1),cmplx(-1,1),cmplx(0,1)) << "\n"
         << "RG(0,0.0796,4) = " << EllipticFunction::RG(0.0796,4)    << " = "
         << EllipticFunction::RG(cmplx(0.0796),4)                    << "\n";

    cout << setprecision(14) << "\n"
         << "RF(1,2,0)      = " << EllipticFunction::RF(1,2)         << " = "
         << EllipticFunction::RF(cmplx(1),2).imag()                         << "\n"
         << "RF(0.5,1)      = " << EllipticFunction::RF(0.5,1)       << " = "
         << EllipticFunction::RF(cmplx(0.5),1).imag()                       << "\n"
         << "RF(i,-i)       = "
         << EllipticFunction::RF(cmplx(0,1),cmplx(0,-1)).imag()             << "\n"
         << "RF(i-1,i)      = "
         << EllipticFunction::RF(cmplx(-1,1),cmplx(0,1)).imag()             << "\n"
         << "RF(2,3,4)      = " << EllipticFunction::RF(2,3,4)       << " = "
         << EllipticFunction::RF(cmplx(2),3,4).imag()                       << "\n"
         << "RF(i,-i,2)     = "
         << EllipticFunction::RF(cmplx(0,1),cmplx(0,-1),2).imag()           << "\n"
         << "RF(i-1,i,1-i)  = "
         << EllipticFunction::RF(cmplx(-1,1),cmplx(0,1),cmplx(1,-1)).imag() << "\n\n"
         << "RC(0,1/4)      = " << EllipticFunction::RC(0,0.25)      << " = "
         << EllipticFunction::RC(cmplx(0),0.25).imag()                      << "\n"
         << "RC(9/4,2)      = " << EllipticFunction::RC(2.25,2)      << " = "
         << EllipticFunction::RC(cmplx(2.25),2).imag()                      << "\n"
         << "RC(0,i)        = "
         << EllipticFunction::RC(0,cmplx(0,1)).imag()                       << "\n"
         << "RC(-i,i)       = "
         << EllipticFunction::RC(cmplx(0,-1),cmplx(0,1)).imag()             << "\n"
         << "RC(1/4,-2)     = " << EllipticFunction::RC(0.25,-2)     << " = "
         << EllipticFunction::RC(cmplx(0.25),-2).imag()                     << "\n"
         << "RC(i,-1)       = "
         << EllipticFunction::RC(cmplx(0,1),-1).imag()                      << "\n\n"
         << "RJ(0,1,2,3)    = " << EllipticFunction::RJ(0,1,2,3)     << " = "
         << EllipticFunction::RJ(cmplx(0),1,2,3).imag()                     << "\n"
         << "RJ(2,3,4,5)    = " << EllipticFunction::RJ(2,3,4,5)     << " = "
         << EllipticFunction::RJ(cmplx(2),3,4,5).imag()                     << "\n"
         << "RJ(2,3,4,-1+i) = "
         << EllipticFunction::RJ(2,3,4,cmplx(-1,1)).imag()                  << "\n"
         << "RJ(i,-i,0,2)   = "
         << EllipticFunction::RJ(cmplx(0,1),cmplx(0,-1),0,2).imag()         << "\n"
         << "RJ(-1+i,-1-i,1,2)    = "
         << EllipticFunction::RJ(cmplx(-1,1),cmplx(-1,-1),1,2).imag()       << "\n"
         << "RJ(i,-i,0,1-i)       = "
         << EllipticFunction::RJ(cmplx(0,1),cmplx(0,-1),0,cmplx(1,-1)).imag()  << "\n"
         << "RJ(-1+i,-1-i,1,-3+i) = "
         << EllipticFunction::RJ(cmplx(-1,1),cmplx(-1,-1),1,cmplx(-3,1)).imag()<< "\n"
         << "RJ(-1+i,-2-i,-i,-1+i)= "
         << EllipticFunction::RJ(cmplx(-1,1),cmplx(-2,-1),
                                 cmplx(0,-1),cmplx(-1,1)).imag()            << "\n"
         << "RJ(2,3,4,-1/2) = " << EllipticFunction::RJ(2,3,4,-0.5)  << " = "
         << EllipticFunction::RJ(cmplx(2),3,4,-0.5).imag()                  << "\n"
         << "RJ(3,4,2,-1/2) = " << EllipticFunction::RJ(3,4,2,-0.5)  << " = "
         << EllipticFunction::RJ(cmplx(3),4,2,-0.5).imag()                  << "\n"
         << "RJ(4,2,3,-1/2) = " << EllipticFunction::RJ(4,2,3,-0.5)  << " = "
         << EllipticFunction::RJ(cmplx(4),2,3,-0.5).imag()                  << "\n"
         << "RJ(2,3,4,-5)   = " << EllipticFunction::RJ(2,3,4,-5)    << " = "
         << EllipticFunction::RJ(cmplx(2),3,4,-5).imag()                    << "\n\n"
         << "RD(0,2,1)      = " << EllipticFunction::RD(0,2,1)       << " = "
         << EllipticFunction::RD(cmplx(0),2,1).imag()                       << "\n"
         << "RD(2,3,4)      = " << EllipticFunction::RD(2,3,4)       << " = "
         << EllipticFunction::RD(cmplx(2),3,4).imag()                       << "\n"
         << "RD(i,-i,2)     = "
         << EllipticFunction::RD(cmplx(0,1),cmplx(0,-1),2).imag()           << "\n"
         << "RD(0,i,-i)     = "
         << EllipticFunction::RD(0,cmplx(0,1),cmplx(0,-1)).imag()           << "\n"
         << "RD(0,i-1,i)    = "
         << EllipticFunction::RD(0,cmplx(-1,1),cmplx(0,1)).imag()           << "\n"
         << "RD(-2-i,-i,-1+i)= "
         << EllipticFunction::RD(cmplx(-2,-1),cmplx(0,-1),cmplx(-1,1)).imag() << "\n\n"
         << "RG(0,16,16)    = " << EllipticFunction::RG(16,16)       << " = "
         << EllipticFunction::RG(cmplx(16),16).imag()                       << "\n"
         << "RG(2,3,4)      = " << EllipticFunction::RG(2,3,4)       << " = "
         << EllipticFunction::RG(cmplx(2),3,4).imag()                       << "\n"
         << "RG(0,i,-i)     = "
         << EllipticFunction::RG(0,cmplx(0,1),cmplx(0,-1)).imag()           << "\n"
         << "RG(i-1,i,0)    = "
         << EllipticFunction::RG(cmplx(-1,1),cmplx(0,1)).imag()             << "\n"
         << "RG(-i,i-1,i)   = "
         << EllipticFunction::RG(cmplx(0,-1),cmplx(-1,1),cmplx(0,1)).imag() << "\n"
         << "RG(0,0.0796,4) = " << EllipticFunction::RG(0.0796,4)    << " = "
         << EllipticFunction::RG(cmplx(0.0796),4).imag()                    << "\n";
  }
  catch (const exception& e) {
    cout << "Caught exception: " << e.what() << "\n";
  }
}

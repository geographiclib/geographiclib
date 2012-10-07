#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <GeographicLib/MagneticModel.hpp>
#include <GeographicLib/SphericalHarmonic.hpp>
#include <GeographicLib/CircularEngine.hpp>
#include <GeographicLib/MagneticCircle.hpp>
#include <GeographicLib/NormalGravity.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Geoid.hpp>

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions and potentially
// uninitialized local variables
#  pragma warning (disable: 4127 4701)
#endif

using namespace GeographicLib;
int main() {
  typedef GeographicLib::Math::real real;
  try {
    {
      std::cout << std::fixed;
      // GravityModel egm("egm2008");
      // GravityModel egm("egm84-mod","/home/ckarney/geographiclib/gravity");
      GravityModel egm("egm2008");
      real lat, lon;
      while (std::cin >> lat >> lon) {
        real g = egm.GeoidHeight(lat, lon);
        std::cout << std::setprecision(6)
                  << std::setw(12) << lat
                  << std::setw(12) << lon
                  << std::setprecision(12) << std::setw(20)
                  << g << "\n";
        /*
          real h;
          std::cin >> h;
        real Dg01, xi, eta;
        egm.Anomaly(lat, lon, h, Dg01, xi, eta);
        Dg01 *= 1e5;
        xi *= 3600;
        eta *= 3600;
        std::cout << std::setprecision(12) << g << " " << Dg01 << " "
                  << xi << " " << eta << "\n";
        */
      }
    return 0;
    }
    {
      std::cout << std::fixed << std::setprecision(6);
      GravityModel egm("egm2008");
      real lat, lon;
      while (std::cin >> lat >> lon) {
        real
          h = 0,
          gamma = egm.ReferenceEllipsoid().SurfaceGravity(lat),
          U0 = egm.ReferenceEllipsoid().SurfacePotential();
        real deltax, deltay, deltaz;
        real T = egm.Disturbance(lat, lon, h, deltax, deltay, deltaz);
        h = T/gamma;
        real h0 = h;
        //        std::cout << "I " << U0 << " " << h << "\n";
        for (int i = 0; i < 10; ++i) {
          real gx, gy, gz;
          real W = egm.Gravity(lat, lon, h, gx, gy, gz);
          h += (W-U0)/gamma;
          if (i == 9)
          std::cout << i << " " << W << " " << h <<  " " << h - h0 << "\n";
        }
      }
      return 0;
    }
    if (false) {
    {
      std::cout << std::fixed;
      Geoid egm("egm2008-1");
      real lat0 = 89.234412, lon0 = -179.234257, dl = 180/1123, h = 0;
      for (int j = 0; j < 2000; ++j) {
        real lat = lat0 - j*dl/2;
        for (int i = 0; i < 4000; ++i) {
          real lon = lon0 + i * dl/2;
        // std::cout << std::setprecision(3) << lon << " "
        // << std::setprecision(12) << egm.GeoidHeight(lat,lon) << "\n";
          h+=egm(lat,lon);
        }
      }
      std::cout << std::setprecision(3) << 0/*lon*/ << " "
                << std::setprecision(12) << h/*c.GeoidHeight(lon)*/ << "\n";
      return 0;
    }

    {
      std::cout << std::fixed;
      GravityModel egm("egm2008","/home/ckarney/geographiclib/gravity");
      real lat = 33, lon0 = 44, dlon = 0.001;
      GravityCircle c = egm.Circle(lat, 0.0);
      for (int i = 0; i < 100000; ++i) {
        real lon = lon0 + i * dlon;
        // std::cout << std::setprecision(3) << lon << " "
        // << std::setprecision(12) << egm.GeoidHeight(lat,lon) << "\n";
        std::cout << std::setprecision(3) << lon << " "
                  << std::setprecision(12) << c.GeoidHeight(lon) << "\n";
      }
      return 0;
    }
    }
    if (true) {
      std::cout << std::setprecision(8);
      // MagneticModel mag("/scratch/WMM2010NewLinux/WMM2010ISO.COF");
      MagneticModel mag1("wmm2010");
      MagneticModel mag2("emm2010");
      MagneticModel mag3("igrf11");
      real lat, lon, h, t, bx, by, bz, bxt, byt, bzt;
      while (std::cin >> t >> lat >> lon >> h) {
        mag1(t, lat, lon, h, bx, by, bz, bxt, byt, bzt);
        std::cout << by << " " << bx << " " << -bz << " "
                  << byt << " " << bxt << " " << -bzt << "\n";
        MagneticCircle circ1(mag1.Circle(t, lat, h));
        circ1(lon, bx, by, bz, bxt, byt, bzt);
        std::cout << by << " " << bx << " " << -bz << " "
                  << byt << " " << bxt << " " << -bzt << "\n";
        /*
        mag2(t, lat, lon, h, bx, by, bz, bxt, byt, bzt);
        std::cout << by << " " << bx << " " << -bz << " "
                  << byt << " " << bxt << " " << -bzt << "\n";
        MagneticCircle circ2(mag2.Circle(t, lat, h));
        circ2(lon, bx, by, bz, bxt, byt, bzt);
        std::cout << by << " " << bx << " " << -bz << " "
                  << byt << " " << bxt << " " << -bzt << "\n";
        mag3(t, lat, lon, h, bx, by, bz, bxt, byt, bzt);
        std::cout << by << " " << bx << " " << -bz << " "
                  << byt << " " << bxt << " " << -bzt << "\n";
        MagneticCircle circ3(mag3.Circle(t, lat, h));
        circ3(lon, bx, by, bz, bxt, byt, bzt);
        std::cout << by << " " << bx << " " << -bz << " "
                  << byt << " " << bxt << " " << -bzt << "\n";
        */
      }
      return 0;
    }
    int type = 2;
    int N;
    std::string name;
    switch (type) {
    case 0:
      N = 360;
      name = "data_EGM96_360.dat";
      break;
    case 1:
      N = 2190;
      name = "data_EGM2008_2190.dat";
      break;
    case 2:
    default:
      N = 5;
      name = "harmtest.dat";
      break;
    }
    int k = ((N + 1) * (N + 2)) / 2;
    std::vector<real> C(k);
    std::vector<real> S(k);
    name = "/scratch/egm2008/harm/" + name;
    {
      std::ifstream f(name.c_str(), std::ios::binary);
      if (!f.good())
        throw GeographicErr("Cannot open coefficient file");
      f.read(reinterpret_cast<char *>(&C[0]), k * sizeof(real));
      f.read(reinterpret_cast<char *>(&S[0]), k * sizeof(real));
    }
    //    for (int i = 0; i < k; ++i)
    //      std::cout << i << " " << C[i] << " " << S[i] << "\n";
    real lat, lon;
    std::cout << std::setprecision(17);
    real a(0.9L), r(1.2L);
    SphericalHarmonic harm(C, S, N, a, SphericalHarmonic::FULL);
    std::vector<real> Z;
    while (std::cin >> lat >> lon) {
      real
        phi = Math::degree<real>() * lat,
        lam = Math::degree<real>() * lon,
        x = r * (abs(lat) == 90 ? 0 : cos(phi)) * cos(lam),
        y = r * (abs(lat) == 90 ? 0 : cos(phi)) * sin(lam),
        z = r * sin(phi);
      real
        d = 1e-7L,
        dx1 = (harm(x+d, y, z) - harm(x-d, y, z))/(2*d),
        dy1 = (harm(x, y+d, z) - harm(x, y-d, z))/(2*d),
        dz1 = (harm(x, y, z+d) - harm(x, y, z-d))/(2*d),
        dx2, dy2, dz2;
      real
        v1 = harm(x, y, z);
      real
        v2 = harm(x, y, z, dx2, dy2, dz2);
      std::cout << v1 << " " << v2 << "\n";
      std::cout << dx1 << " " << dx2 << "\n"
                << dy1 << " " << dy2 << "\n"
                << dz1 << " " << dz2 << "\n";
      CircularEngine circ1 = harm.Circle(r * cos(phi), z, false);
      CircularEngine circ2 = harm.Circle(r * cos(phi), z, true);
      real v3, dx3, dy3, dz3;
      v3 = circ2(lon, dx3, dy3, dz3);
      std::cout << v3 << " " << dx3 << " " << dy3 << " " << dz3 << "\n";
    }
    std::cout << "start timing" << std::endl;
    real phi, lam, sum, z, p, dx, dy, dz;
    lat = 33; phi = lat * Math::degree<real>();
    z = r * sin(phi);
    p = r * cos(phi);
    sum = 0;
    for (int i = 0; i < 100; ++i) {
      lam = (44 + 0.01*i) * Math::degree<real>();
      sum += harm(p * cos(lam), p * sin(lam), z, dx, dy, dz);
    }
    std::cout << "sum a " << sum << std::endl;
    CircularEngine circ(harm.Circle(p, z, true));
    sum = 0;
    for (int i = 0; i < 100000; ++i) {
      lon = (44 + 0.01*i);
      sum += circ(lon, dx, dy, dz);
    }
    std::cout << "sum b " << sum << std::endl;
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

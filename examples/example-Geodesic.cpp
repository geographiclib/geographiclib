// Example of using GeographicLib::Geodesic class
// $Id$

#include <iostream>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  Geodesic geodesic(Constants::WGS84_a(), Constants::WGS84_f());
  // Alternatively: const Geodesic& geodesic = Geodesic::WGS84;
  {
    // Sample direct calculation
    double lat1 = 43, lon1 = -75, s12 = 20e6, azi1 = -3;
    double lat2, lon2;
    geodesic.Direct(lat1, lon1, azi1, s12, lat2, lon2);
    cout << lat2 << " " << lon2 << "\n";
  }
  {
    // Sample inverse calculation
    double lat1 = 43, lon1 = -75, lat2 = 60, lon2 = 0;
    double s12;
    geodesic.Inverse(lat1, lon1, lat2, lon2, s12);
    cout << s12 << "\n";
  }
}

// Example of using GeographicLib::Geodesic class
// $Id$

#include <iostream>
#include <exception>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
    // Alternatively: const Geodesic& geod = Geodesic::WGS84;
    {
      // Sample direct calculation
      double lat1 = 43, lon1 = -75, s12 = 20e6, azi1 = -3;
      double lat2, lon2;
      geod.Direct(lat1, lon1, azi1, s12, lat2, lon2);
      cout << lat2 << " " << lon2 << "\n";
    }
    {
      // Sample inverse calculation
      double lat1 = 43, lon1 = -75, lat2 = 60, lon2 = 0;
      double s12;
      geod.Inverse(lat1, lon1, lat2, lon2, s12);
      cout << s12 << "\n";
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

// Example of using the GeographicLib::OSGB class

#include <iostream>
#include <exception>
#include <string>
#include <GeographicLib/OSGB.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    {
      // Sample forward calculation
      double lat = 55.5, lon = -1.64; // Embleton, Northumberland
      double x, y;
      OSGB::Forward(lat, lon, x, y);
      string gridref;
      OSGB::GridReference(x, y, 2, gridref);
      cout << x << " " << y << " " << gridref << "\n";
    }
    {
      // Sample reverse calculation
      string gridref = "NU2222";
      double x, y;
      int prec;
      OSGB::GridReference(gridref, x, y, prec);
      double lat, lon;
      OSGB::Reverse(x, y, lat, lon);
      cout << prec << " " << lat << " " << lon << "\n";
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

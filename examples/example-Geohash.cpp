// Example of using the GeographicLib::MGRS class
// $Id$

#include <iostream>
#include <iomanip>
#include <exception>
#include <string>
#include <GeographicLib/Geohash.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    {
      // Sample forward calculation
      double lat = 57.64911, lon = 10.40744; // Jutland (the wikipedia example)
      string geohash;
      int maxprec = Geohash::GeohashPrecision(1.0e-5);
      for (int prec = 0; prec <= maxprec; ++prec) {
        Geohash::Forward(lat, lon, prec, geohash);
        cout << prec << " " << geohash << "\n";
      }
    }
    {
      // Sample reverse calculation
      string geohash = "u4pruydqqvj";
      double lat, lon;
      cout << fixed;
      for (unsigned i = 0; i <= geohash.length(); ++i) {
        int prec;
        Geohash::Reverse(geohash.substr(0, i), lat, lon, prec);
        cout << setprecision(max(0, Geohash::DecimalPrecision(prec)))
             << prec << " " << lat << " " << lon << "\n";
      }
    }
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

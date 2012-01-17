// Example of using GeographicLib::MGRS class
// $Id$

#include <iostream>
#include <string>
#include <GeographicLib/UTMUPS.hpp>
#include <GeographicLib/MGRS.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  {
    // Sample forward calculation
    double lat = 33.3, lon = 44.4;
    int zone;
    bool northp;
    double x, y;
    UTMUPS::Forward(lat, lon, zone, northp, x, y);
    string mgrs;
    MGRS::Forward(zone, northp, x, y, 5, mgrs);
    cout << mgrs << "\n";
  }
  {
    // Sample reverse calculation
    string mgrs = "38SMB4488";
    int zone, prec;
    bool northp;
    double x, y;
    MGRS::Reverse(mgrs, zone, northp, x, y, prec);
    double lat, lon;
    UTMUPS::Reverse(zone, northp, x, y, lat, lon);
    cout << prec << " " << lat << " " << lon << "\n";
  }
}

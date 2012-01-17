// Example of using GeographicLib::UTMUPS class
// $Id$

#include <iostream>
#include <string>
#include <iomanip>
#include <GeographicLib/UTMUPS.hpp>

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
    string zonestr = UTMUPS::EncodeZone(zone, northp);
    cout << fixed << setprecision(2)
         << zonestr << " " << x << " " << y << "\n";
  }
  {
    // Sample reverse calculation
    string zonestr = "38N";
    int zone;
    bool northp;
    UTMUPS::DecodeZone(zonestr, zone, northp);
    double x = 444e3, y = 3688e3;
    double lat, lon;
    UTMUPS::Reverse(zone, northp, x, y, lat, lon);
    cout << lat << " " << lon << "\n";
  }
}

// Example of using the GeographicLib::NearestNeighbor class.  WARNING: this
// creates a file, vptree.xml or vptree.bin, in the current directory.

#include <iostream>

#include <vector>
#include <cstdlib>              // For srand, rand
#include <cmath>                // For asin
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>            // For sort
#include <GeographicLib/NearestNeighbor.hpp>
#include <GeographicLib/Geodesic.hpp>

#if !defined(GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION)
#define GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION 0
#endif

#if GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION
// If Boost serialization is available, use it.
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#endif

using namespace std;
using namespace GeographicLib;

// A structure to hold a geographic coordinate.  Also included is a field for a
// "name".  This is unused in this example.
struct pos {
  double lat, lon;
  string name;
  pos(double lat = 0, double lon = 0, const string& name = "")
    : lat(lat), lon(lon), name(name) {}
};

pos randompos() {
  double r, lat, lon;
  r = 2 * (rand() + 0.5) / (RAND_MAX + 1.0) - 1;
  lat = asin(r) / Math::degree();
  r = 2 * (rand() + 0.5) / (RAND_MAX + 1.0) - 1;
  lon = 180 * r;
  return pos(lat, lon);
}

// A class to compute the distance between 2 positions.
class DistanceCalculator {
private:
  Geodesic _geod;
public:
  explicit DistanceCalculator(const Geodesic& geod)
    : _geod(geod) {}
  double operator() (const pos& a, const pos& b) const {
    double s12;
    _geod.Inverse(a.lat, a.lon, b.lat, b.lon, s12);
    return s12;
  }
};

typedef NearestNeighbor<double, pos, DistanceCalculator> GeodesicNeighbor;

// Pick 10000 points on the ellipsoid and determine which ones are more than
// 350 km from all the others.

// In this example the NearestNeighbor object is saved to an external file and
// read back in.  This is unnecessary in this simple application, but is useful
// if many different applications need to query the same dataset.

int main() {
  try {
    // Define a distance function object
    DistanceCalculator distance(Geodesic::WGS84());
    srand(0);
    vector<pos> pts;
    int num = 10000;
    // Sample the points
    for (int i = 0; i < num; ++i) pts.push_back(randompos());
    {
      // Illustrate saving and restoring the GeodesicNeighbor
      // construct it
      GeodesicNeighbor posset(pts, distance);
      // and save it
#if GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION
      ofstream f("vptree.xml");
      boost::archive::xml_oarchive oa(f);
      oa << BOOST_SERIALIZATION_NVP(posset);
#else
      ofstream ofs("vptree.txt");
      ofs << posset << "\n";
#endif
    }
    // Construct an empty GeodesicNeighbor
    GeodesicNeighbor posset;
    // restore it from the file
    {
#if GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION
      ifstream f("vptree.xml");
      boost::archive::xml_iarchive ia(f);
      ia >> BOOST_SERIALIZATION_NVP(posset);
#else
      ifstream ifs("vptree.txt");
      ifs >> posset;
#endif
    }
    // Now use it
    vector<int> ind;
    int cnt = 0;
    double thresh = 325000;
    cout << "Points more than " << thresh/1000 << "km from their neighbors\n"
         << "latitude longitude distance\n";
    for (int i = 0; i < num; ++i) {
      // Call search with distance limits = (0, thresh].  Set exhaustive = false
      // so that the search ends as some as a neighbor is found.
      posset.Search(pts, distance, pts[i], ind, 1, thresh, 0, false);
      if (ind.size() == 0) {
        // If no neighbors in (0, thresh], search again with no upper limit and
        // with exhaustive = true (the default).
        double d = posset.Search(pts, distance, pts[i], ind, 1,
                                 numeric_limits<double>::max(), 0);
        cout << pts[i].lat << " " << pts[i].lon << " " << d << "\n";
        ++cnt;
      }
    }
    int setupcost, numsearches, searchcost, mincost, maxcost;
    double mean, sd;
    posset.Statistics(setupcost, numsearches, searchcost, mincost, maxcost,
                      mean, sd);
    int totcost = setupcost + searchcost, exhaustivecost = num * (num - 1) / 2;
    cout
      << "Number of distance calculations = " << totcost << "\n"
      << "With an exhaustive search = " << exhaustivecost << "\n"
      << "Ratio = " << double(totcost) / exhaustivecost << "\n"
      << "Efficiency improvement = "
      << 100 * (1 - double(totcost) / exhaustivecost) << "%\n";
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}

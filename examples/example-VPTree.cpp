// Example of using the GeographicLib::VPTree class

#include <iostream>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <GeographicLib/VPTree.hpp>
#include <GeographicLib/Geodesic.hpp>

using namespace std;
using namespace GeographicLib;

// A structure to hold a geographic coordinate.  This could include additional
// information such as a place name.
struct pos {
  double lat, lon;
  pos(double lat = 0, double lon = 0) : lat(lat), lon(lon) {}
};

pos randompos() {
    double r = 2 * (rand() + 0.5) / (double(RAND_MAX) + 1) - 1;
    double lat = asin(r) / Math::degree();
    r = 2 * (rand() + 0.5) / (double(RAND_MAX) + 1) - 1;
    double lon = 180 * r;
    return pos(lat, lon);
}

// A class to compute the distance between 2 positions.
class DistanceCalculator {
private:
  const Geodesic& _geod;
  // copy constructor and assignment not allowed
  DistanceCalculator(const DistanceCalculator&);
  DistanceCalculator& operator=(const DistanceCalculator&);
public:
  DistanceCalculator(const Geodesic& geod)
    : _geod(geod) {}
  double operator()(const pos& a, const pos& b) const {
    double s12;
    _geod.Inverse(a.lat, a.lon, b.lat, b.lon, s12);
    return s12;
  }
};

// Pick 10000 points on the ellipsoid and determine which ones are more than
// 350 km from all the others.

int main() {
  // Define a distance function object
  DistanceCalculator distance(Geodesic::WGS84());
  srand(0);
  vector<pos> pts;
  int num = 10000;
  // Sample the points
  for (int i = 0; i < num; ++i) pts.push_back(randompos());
  // Set up the VP tree
  VPTree<double, pos, DistanceCalculator> posset(pts, distance);
  vector<int> ind;
  int cnt = 0;
  cout << "Points more than 350km from their neighbors\n"
       << "latitude longitude distance\n";
  for (int i = 0; i < num; ++i) {
    // Call search with distance limits = (0, 350e3].  Set exhaustive = false
    // so that the search ends as some as a neighbor is found.
    posset.search(pts[i], ind, 1, 350e3, 0, false);
    if (ind.size() == 0) {
      // If no neighbors in (0, 350e3], search again with no upper limit and
      // with exhaustive = true (the default).
      posset.search(pts[i], ind, 1, numeric_limits<double>::max(), 0);
      cout << pts[i].lat << " " << pts[i].lon << " "
           << distance(pts[i], pts[ind[0]]) << "\n";
      ++cnt;
    }
  }
  posset.report(cout);
}

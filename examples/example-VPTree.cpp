// Example of using the GeographicLib::VPTree class

#include <iostream>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <GeographicLib/VPTree.hpp>
#include <GeographicLib/Geodesic.hpp>
#if GEOGRAPHICIB_HAVE_BOOST_SERIALIZATION
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#endif

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
  if (1) {
    srand(0);
  vector<pos> pts;
  int num = 10000;
  // Sample the points
  for (int i = 0; i < num; ++i) pts.push_back(randompos());
  {
    // Set up the VP tree
    VPTree<double, pos, DistanceCalculator> posset(pts, distance);
    ofstream ofs("vptree.bin", std::ios::binary);
#if GEOGRAPHICIB_HAVE_BOOST_SERIALIZATION
    {
      std::ofstream f("vptree.xml");
      boost::archive::xml_oarchive oa(f); // set up an xml archive
      oa << posset;                       // save tree to xml file vptree.xml
    }
#endif
    posset.save(ofs, true);
  }
  ifstream ifs("vptree.bin", std::ios::binary);
  VPTree<double, pos, DistanceCalculator> posset(pts, distance, ifs, true);
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
  if (0) {
    srand(1);
    vector<pos> ptsa, ptsb;
    int numa = 10000, numb = 5000;
    // Sample the points
    for (int i = 0; i < numa; ++i) ptsa.push_back(randompos());
    for (int i = 0; i < numb; ++i) ptsb.push_back(randompos());
    VPTree<double, pos, DistanceCalculator> seta(ptsa, distance);
    VPTree<double, pos, DistanceCalculator> setb(ptsb, distance);
    vector<int> ind;
    double d0 = 10e3, d1 = 100e3;
    for (int j = 0; j < numb; ++j) {
      double d = seta.search(ptsb[j], ind, 1, d0);
      if (ind.size() != 1) continue;
      int i = ind[0];
      setb.search(ptsa[i], ind, 1, d);
      if (ind[0] != j) continue;
      seta.search(ptsb[j], ind, 2, d1, false);
      if (ind.size() == 2) continue;
      setb.search(ptsa[i], ind, 2, d1, false);
      if (ind.size() == 2) continue;
      std::cout << i << " " << j << "\n";
    }
    seta.report(std::cout);
    setb.report(std::cout);
  }
}

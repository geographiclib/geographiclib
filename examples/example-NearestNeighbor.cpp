// Example of using the GeographicLib::NearestNeighbor class

#include <iostream>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <GeographicLib/NearestNeighbor.hpp>
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

typedef NearestNeighbor<double, pos, DistanceCalculator> GeodesicNeighbor;

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
    GeodesicNeighbor posset(pts, distance);
#if GEOGRAPHICIB_HAVE_BOOST_SERIALIZATION
    {
      ofstream f("vptree.xml");
      boost::archive::xml_oarchive oa(f); // set up an xml archive
      oa << BOOST_SERIALIZATION_NVP(posset);
    }
#endif
    {
      ofstream ofs("vptree.bin", ios::binary);
      posset.Save(ofs, true);
    }
  }
  GeodesicNeighbor posset;
#if GEOGRAPHICIB_HAVE_BOOST_SERIALIZATION
  {
    ifstream f("vptree.xml");
    boost::archive::xml_iarchive ia(f);
    ia >> BOOST_SERIALIZATION_NVP(posset);
  }
#else
  {
    ifstream ifs("vptree.bin", ios::binary);
    posset.Load(ifs);
  }
#endif
  vector<int> ind;
  int cnt = 0;
  cout << "Points more than 350km from their neighbors\n"
       << "latitude longitude distance\n";
  for (int i = 0; i < num; ++i) {
    // Call search with distance limits = (0, 350e3].  Set exhaustive = false
    // so that the search ends as some as a neighbor is found.
    posset.Search(pts, distance, pts[i], ind, 1, 350e3, 0, false);
    if (ind.size() == 0) {
      // If no neighbors in (0, 350e3], search again with no upper limit and
      // with exhaustive = true (the default).
      posset.Search(pts, distance, pts[i], ind, 1,
                    numeric_limits<double>::max(), 0);
      cout << pts[i].lat << " " << pts[i].lon << " "
           << distance(pts[i], pts[ind[0]]) << "\n";
      ++cnt;
    }
  }
  posset.report(cout);
  }
  if (1) {
    srand(1);
    vector<pos> ptsa, ptsb;
    int numa = 10000, numb = 5000;
    // Sample the points
    for (int i = 0; i < numa; ++i) ptsa.push_back(randompos());
    for (int i = 0; i < numb; ++i) ptsb.push_back(randompos());
    GeodesicNeighbor seta(ptsa, distance);
    GeodesicNeighbor setb(ptsb, distance);
    vector<int> ind;
    double d0 = 10e3, d1 = 100e3;
    for (int j = 0; j < numb; ++j) {
      double d = seta.Search(ptsa, distance, ptsb[j], ind, 1, d0);
      if (ind.size() != 1) continue;
      int i = ind[0];
      setb.Search(ptsb, distance, ptsa[i], ind, 1, d);
      if (ind[0] != j) continue;
      seta.Search(ptsa, distance, ptsb[j], ind, 2, d1, false);
      if (ind.size() == 2) continue;
      setb.Search(ptsb, distance, ptsa[i], ind, 2, d1, false);
      if (ind.size() == 2) continue;
      cout << i << " " << j << "\n";
    }
    seta.report(cout);
    setb.report(cout);
  }
}

// Compute a table of egm2008 geoid heights on a 1' grid.  This takes about 40
// mins on a 8-processor Intel 2.66 GHz machine using OpenMP (-DHAVE_OPENMP=1).
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#if HAVE_OPENMP
#  include <omp.h>
#endif

#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/Utility.hpp>

using namespace std;
using namespace GeographicLib;

int main(int argc, char* argv[]) {
  // Hardwired for 3 args:
  // 1 = EGM (e.g., egm2008)
  // 2 = intervals per degree
  // 3 = output file (must if pgm or gtx)
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " model ninterval output\n";
    return 1;
  }
  try {
    enum format {
      PGM = 0,
      PGM4 = 1,
      GTX = 2,
    };
    format mode = PGM4;

    std::string model(argv[1]);
    // Number of intervals per degree
    int ndeg = Utility::num<int>(std::string(argv[2]));
    std::string filename(argv[3]);
    GravityModel g(model);
    int
      nlat = 180 * ndeg + 1,
      nlon = 360 * ndeg;
    double delta = 1 / double(ndeg); // Grid spacing
    double
      offset = mode == PGM ? -108 : -128,
      scale = mode == PGM ? 0.003 : std::pow(0.5, 24);
    // Write results as floats in binary mode
    ofstream file(filename.c_str(), ios::binary);


    switch (mode) {
    case PGM:
      file << "P5\n"
           << "# Geoid file in PGM format for the GeographicLib::Geoid class\n"
           << "# Description WGS84 " << model
           << ", " << 60*delta << "-minute grid\n"
           << "# Offset " << offset << "\n"
           << "# Scale " << scale << "\n"
           << "# Origin 90N 0E\n"
           << "# AREA_OR_POINT Point\n"
           << "# Vertical_Datum WGS84\n"
           << nlon << " " << nlat << "\n" << (1<<16)-1 << "\n";
      break;
    case PGM4:
      file << std::setprecision(17)
           << "P5\n"
           << "# Geoid file in PGM format for the GeographicLib::Geoid class\n"
           << "# Description WGS84 " << model
           << ", " << 60*delta << "-minute grid\n"
           << "# Offset " << offset << "\n"
           << "# Scale " << scale << "\n"
           << "# Origin 90N 0E\n"
           << "# AREA_OR_POINT Point\n"
           << "# Vertical_Datum WGS84\n"
           << nlon << " " << nlat << "\n" << 0xffffffffu << "\n";
      break;
    case GTX:
      break;
    }

    // Compute and store results for nbatch latitudes at a time
    const int nbatch = 64;
    vector< vector<float> > Nf(mode == GTX ? nbatch : 0, vector<float>(nlon));
    vector< vector<unsigned short> >
      Ns(mode == PGM ? nbatch : 0, vector<unsigned short>(nlon));
    vector< vector<unsigned> >
      Nu(mode == PGM4 ? nbatch : 0, vector<unsigned>(nlon));

    for (int ilat0 = 0; ilat0 < nlat; ilat0 += nbatch) { // Loop over batches
      int nlat0 = min(nlat, ilat0 + nbatch);

#if HAVE_OPENMP
#  pragma omp parallel for
#endif
      for (int ilat = ilat0; ilat < nlat0; ++ilat) { // Loop over latitudes
        double lat = 90 - ilat * delta, h = 0;
        GravityCircle c(g.Circle(lat, h, GravityModel::GEOID_HEIGHT));
        for (int ilon = 0; ilon < nlon; ++ilon) { // Loop over longitudes
          double lon = ilon * delta, h = c.GeoidHeight(lon);
          switch (mode) {
          case PGM:
            Ns[ilat - ilat0][ilon]
              = std::max(0, std::min((1<<16)-1, int((h - offset)/scale + 0.5)));
            break;
          case PGM4:
            Nu[ilat - ilat0][ilon] = unsigned((h - offset)/scale + 0.5);
            break;
          case GTX:
            Nf[ilat - ilat0][ilon] = float(h);
            break;
          }
        } // longitude loop
      }   // latitude loop
      for (int ilat = ilat0; ilat < nlat0; ++ilat)
          switch (mode) {
          case PGM:
            Utility::writearray<short unsigned, short unsigned, true>
              (file, Ns[ilat - ilat0]);
            break;
          case PGM4:
            Utility::writearray<unsigned, unsigned, true>
              (file, Nu[ilat - ilat0]);
            break;
          case GTX:
            file.write(reinterpret_cast<char *>(Nf[ilat - ilat0].data()),
                       nlon * sizeof(float));
            break;
          }

    } // batch loop
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
  return 0;
}
/*
Test point
82.333333   16.250000     26.437499855168

echo 82:20 16:15 | Gravity -n egm2008 -H -p 20
26.437499855163 round to float -> 26.4375
with offset = -108; scale = 0.003
26.437499855163 -> 44812
26.4375         -> 44813

Rate of misrounding = 1/4200

for i in egm*pgm;do echo $i; pamarith -diff /usr/local/share/GeographicLib/geoids/$i $i | pgmhist;done
egm2008-1.pgm
value  count  b%      w%
-----  -----  ------  ------
    0  233245696    100%    100%
    1  55904    100%  0.024%
egm2008-2_5.pgm
value  count  b%      w%
-----  -----  ------  ------
    0  37324574    100%    100%
    1   8866    100%  0.0237%
egm2008-5.pgm
value  count  b%      w%
-----  -----  ------  ------
    0  9333323    100%    100%
    1   2197    100%  0.0235%
egm84-15.pgm
value  count  b%      w%
-----  -----  ------  ------
    0  850300   81.9%    100%
    1  187940    100%   18.1%
egm84-30.pgm
value  count  b%      w%
-----  -----  ------  ------
    0  212756   81.9%    100%
    1  47164    100%   18.1%
egm96-15.pgm
value  count  b%      w%
-----  -----  ------  ------
    0  1038064    100%    100%
    1    176    100%  0.017%
egm96-5.pgm
value  count  b%      w%
-----  -----  ------  ------
    0  9318864   99.8%    100%
    1  16656    100%  0.178%

 */

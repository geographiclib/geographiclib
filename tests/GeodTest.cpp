/**
 * \file GeodTest.cpp
 **********************************************************************/

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/GeodesicLine.hpp"
#include "GeographicLib/GeodesicExact.hpp"
#include "GeographicLib/Constants.hpp"
// Allow for mixed libraries?
namespace GeographicLibL = GeographicLib;
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>

using namespace std;
using namespace GeographicLib;

int usage(int retval) {
  ( retval ? cerr : cout ) <<
"GeodTest [ -a | -c | -t0 | -t1 | -t2 | -t3 | -h ]\n\
\n\
Check GeographicLib::Geodesic class.\n\
-a (default) accuracy test (reads test data on standard input)\n\
-E accuracy test with GeodesicExact (reads test data on standard input)\n\
-F accuracy test with GeodesicExact (reads test data on standard input\n\
   first line gives a and f)\n\
-c coverage test (reads test data on standard input)\n\
-t0 time GeodecicLine with distances using synthetic data\n\
-t1 time GeodecicLine with angles using synthetic data\n\
-t2 time Geodecic::Direct using synthetic data\n\
-t3 time Geodecic::Inverse with synthetic data\n\
\n\
-c requires an instrumented version of Geodesic.\n";
  return retval;
}

Math::extended angdiff(Math::extended a1, Math::extended a2) {
  Math::extended d = a2 - a1;
  if (d >= 180)
    d -= 360;
  else if (d < -180)
    d += 360;
  return d;
}

Math::extended azidiff(Math::extended lat,
                       Math::extended lon1, Math::extended lon2,
                       Math::extended azi1, Math::extended azi2) {
  Math::extended
    phi = lat * Math::degree<Math::extended>(),
    alpha1 = azi1 * Math::degree<Math::extended>(),
    alpha2 = azi2 * Math::degree<Math::extended>(),
    dlam = angdiff(lon1, lon2) * Math::degree<Math::extended>();
  Math::extended res = sin(alpha2-alpha1)*cos(dlam)
    -cos(alpha2-alpha1)*sin(dlam)*sin(phi)
    // -sin(alpha1)*cos(alpha2)*(1-cos(dlam))*cos(phi)*cos(phi)
    ;
  return res;
}

Math::extended dist(Math::extended a, Math::extended f,
                    Math::extended lat0, Math::extended lon0,
                    Math::extended lat1, Math::extended lon1) {
  //  typedef GeographicLibL::Math::real real;
  //  real s12;
  //  GeographicLibL::Geodesic::
  //    WGS84.Inverse(real(lat0), real(lon0), real(lat1), real(lon1), s12);
  //  return Math::extended(s12);
  a *= Math::degree<Math::extended>();
  if (abs(lat0 + lat1) > Math::extended(179.998)) {
    // Near pole, transform into polar coordinates
    Math::extended
      r0 = 90 - abs(lat0),
      r1 = 90 - abs(lat1),
      lam0 = lon0 * Math::degree<Math::extended>(),
      lam1 = lon1 * Math::degree<Math::extended>();
    return (a / (1 - f)) *
      Math::hypot
      (r0 * cos(lam0) - r1 * cos(lam1), r0 * sin(lam0) - r1 * sin(lam1));
  } else {
    // Otherwise use cylindrical formula
    Math::extended
      phi = lat0 * Math::degree<Math::extended>(),
      cphi = abs(lat0) <= 45 ? cos(phi)
      : sin((90 - abs(lat0)) * Math::degree<Math::extended>()),
      e2 = f * (2 - f),
      sinphi = sin(phi),
      n = 1/sqrt(1 - e2 * sinphi * sinphi),
      // See Wikipedia article on latitude
      degreeLon = a * cphi * n,
      degreeLat = a * (1 - e2) * n * n * n,
      dlon = angdiff(lon1, lon0),
      dlat = lat1 - lat0;
    dlat *= degreeLat;
    dlon *= degreeLon;
    return Math::hypot(dlat, dlon);
  }
}

// wreal is precision of args.
// treal is precision of test class.
// rreal is precision of reference class.
//
// err[0] error in position of point 2 for the direct problem.
// err[1] error in azimuth at point 2 for the direct problem.
// err[2] error in m12 for the direct problem & inverse (except near conjugacy)
// err[3] error in s12 for the inverse problem.
// err[4] error in the azimuths for the inverse problem scaled by m12.
// err[5] consistency of the azimuths for the inverse problem.
// err[6] area error direct & inverse (except near conjugacy)
template<class wreal, class test, class treal, class ref, class rreal>
void GeodError(const test& tgeod, const ref& rgeod,
               wreal lat1, wreal lon1, wreal azi1,
               wreal lat2, wreal lon2, wreal azi2,
               wreal s12, wreal /*a12*/, wreal m12, wreal S12,
               vector<wreal>& err) {
  treal tlat1, tlon1, tazi1, tlat2, tlon2, tazi2, ts12, tm12a, tm12b,
    tM12, tM21, tS12a, tS12b /*, ta12*/;
  rreal rlat1, rlon1, razi1, rlat2, rlon2, razi2, rm12;
  tgeod.Direct(treal(lat1), treal(lon1), treal(azi1),  treal(s12),
               tlat2, tlon2, tazi2, tm12a,
               tM12, tM21, tS12a);
  tS12a -= treal(rgeod.EllipsoidArea() * (tazi2-azi2)/720);
  tgeod.Direct(treal(lat2), treal(lon2), treal(azi2), -treal(s12),
               tlat1, tlon1, tazi1, tm12b,
               tM12, tM21, tS12b);
  tS12b -= treal(rgeod.EllipsoidArea() * (tazi1-azi1)/720);
  err[0] = max(dist(rgeod.MajorRadius(), rgeod.Flattening(),
                    lat2, lon2, tlat2, tlon2),
               dist(rgeod.MajorRadius(), rgeod.Flattening(),
                    lat1, lon1, tlat1, tlon1));
  err[1] = max(abs(azidiff(lat2, lon2, tlon2, azi2, tazi2)),
               abs(azidiff(lat1, lon1, tlon1, azi1, tazi1))) *
    rgeod.MajorRadius();
  err[2] = max(abs(tm12a - m12), abs(tm12b + m12));
  if (!Math::isnan(S12))
    err[6] = max(abs(tS12a - S12), abs(tS12b + S12)) / rgeod.MajorRadius();

  /* ta12 = */ tgeod.Inverse(treal(lat1), treal(lon1), treal(lat2), treal(lon2),
                             ts12, tazi1, tazi2, tm12a,
                             tM12, tM21, tS12a);
  tS12a -= treal(rgeod.EllipsoidArea() * ((tazi2-azi2)-(tazi1-azi1))/720);
  err[3] = abs(ts12 - s12);
  err[4] = max(abs(angdiff(azi1, tazi1)), abs(angdiff(azi2, tazi2))) *
    Math::degree<Math::extended>() * abs(m12);
  if (treal(lat1) + treal(lat2) == 0)
    err[4] = min(err[4],
                 max(abs(angdiff(azi1, tazi2)), abs(angdiff(azi2, tazi1))) *
                 Math::degree<Math::extended>() * abs(m12));
  // m12 and S12 are very sensitive with the inverse problem near conjugacy
  if (!(s12 > rgeod.MajorRadius() && m12 < 10e3)) {
    err[2] = max(err[2], wreal(abs(tm12a - m12)));
    if (!Math::isnan(S12))
      err[6] = max(err[6], wreal(abs(tS12a - S12) / rgeod.MajorRadius()));
  }
  if (s12 > rgeod.MajorRadius()) {
    rgeod.Direct(rreal(lat1), rreal(lon1), rreal(tazi1),   rreal(ts12)/2,
                 rlat2, rlon2, razi2, rm12);
    rgeod.Direct(rreal(lat2), rreal(lon2), rreal(tazi2), - rreal(ts12)/2,
                 rlat1, rlon1, razi1, rm12);
    err[5] = dist(rgeod.MajorRadius(), rgeod.Flattening(),
                  rlat1, rlon1, rlat2, rlon2);
  } else {
    rgeod.Direct(rreal(lat1), rreal(lon1), rreal(tazi1),
                 rreal(ts12) + rgeod.MajorRadius(),
                 rlat2, rlon2, razi2, rm12);
    rgeod.Direct(rreal(lat2), rreal(lon2), rreal(tazi2), rgeod.MajorRadius(),
                 rlat1, rlon1, razi1, rm12);
    err[5] = dist(rgeod.MajorRadius(), rgeod.Flattening(),
                  rlat1, rlon1, rlat2, rlon2);
    rgeod.Direct(rreal(lat1), rreal(lon1), rreal(tazi1), - rgeod.MajorRadius(),
                 rlat2, rlon2, razi2, rm12);
    rgeod.Direct(rreal(lat2), rreal(lon2), rreal(tazi2),
                 - rreal(ts12) - rgeod.MajorRadius(),
                 rlat1, rlon1, razi1, rm12);
    err[5] = max(err[5], wreal(dist(rgeod.MajorRadius(), rgeod.Flattening(),
                                    rlat1, rlon1, rlat2, rlon2)));
  }
}


int main(int argc, char* argv[]) {
  Math::real a = Constants::WGS84_a<Math::real>();
  Math::real f = Constants::WGS84_f<Math::real>();
  Math::real al =
    GeographicLibL::Constants::WGS84_a<GeographicLibL::Math::real>();
  Math::real fl =
    GeographicLibL::Constants::WGS84_f<GeographicLibL::Math::real>();
  bool timing = false;
  int timecase = 0; // 0 = line, 1 = line ang, 2 = direct, 3 = inverse
  bool accuracytest = true;
  bool coverage = false;
  bool exact = false;
  if (argc == 2) {
    string arg = argv[1];
    if (arg == "-a") {
      accuracytest = true;
      coverage = false;
      timing = false;
      exact = false;
    } else if (arg == "-E") {
      accuracytest = true;
      coverage = false;
      timing = false;
      exact = true;
    } else if (arg == "-F") {
      accuracytest = true;
      coverage = false;
      timing = false;
      exact = true;
      string s;
      getline(cin, s);
      istringstream str(s);
      str >> al >> fl;
      a = al; f = fl;
    } else if (arg == "-c") {
      accuracytest = false;
      coverage = true;
      timing = false;
      exact = false;
    } else if (arg == "-t0") {
      accuracytest = false;
      coverage = false;
      timing = true;
      timecase = 0;
      exact = false;
    } else if (arg == "-t1") {
      accuracytest = false;
      coverage = false;
      timing = true;
      timecase = 1;
      exact = false;
    } else if (arg == "-t2") {
      accuracytest = false;
      coverage = false;
      timing = true;
      timecase = 2;
      exact = false;
    } else if (arg == "-t3") {
      accuracytest = false;
      coverage = false;
      timing = true;
      timecase = 3;
      exact = false;
    } else
      return usage(arg == "-h" ? 0 : 1);
  } else if (argc > 2)
    return usage(1);

  if (timing) {
    const Geodesic& geod = Geodesic::WGS84;
    unsigned cnt = 0;
    Math::real s = 0;
    Math::real dl;
    switch (timecase) {
    case 0:
      // Time Line
      dl = 2e7/1000;
      for (int i = 0; i <= 90; ++i) {
        Math::real lat1 = i;
        for (int j = 0; j <= 180; ++j) {
          Math::real azi1 = j;
          const GeodesicLine l(geod, lat1, 0.0, azi1);
          for (int k = 0; k <= 1000; ++k) {
            Math::real s12 = dl * k;
            Math::real lat2, lon2;
            l.Position(s12, lat2, lon2);
            ++cnt;
            s += lat2;
          }
        }
      }
      cout << cnt << " " << s << "\n";
      break;
    case 1:
      // Time Line ang
      dl = 180.0/1000;
      for (int i = 0; i <= 90; ++i) {
        Math::real lat1 = i;
        for (int j = 0; j <= 180; ++j) {
          Math::real azi1 = j;
          GeodesicLine l(geod, lat1, 0.0, azi1);
          for (int k = 0; k <= 1000; ++k) {
            Math::real s12 = dl * k;
            Math::real lat2, lon2;
            l.ArcPosition(s12, lat2, lon2);
            ++cnt;
            s += lat2;
          }
        }
      }
      cout << cnt << " " << s << "\n";
      break;
    case 2:
      // Time Direct
      dl = 2e7/200;
      for (int i = 0; i <= 90; ++i) {
        Math::real lat1 = i;
        for (int j = 0; j <= 180; ++j) {
          Math::real azi1 = j;
          for (int k = 0; k <= 200; ++k) {
            Math::real s12 = dl * k;
            Math::real lat2, lon2;
            geod.Direct(lat1, 0.0, azi1, s12, lat2, lon2);
            ++cnt;
            s += lat2;
          }
        }
      }
      cout << cnt << " " << s << "\n";
      break;
    case 3:
      // Time Inverse
      for (int i = 0; i <= 179; ++i) {
        Math::real lat1 = i * 0.5;
        for (int j = -179; j <= 179; ++j) {
          Math::real lat2 = j * 0.5;
          for (int k = 0; k <= 359; ++k) {
            Math::real lon2 = k * 0.5;
            Math::real s12;
            geod.Inverse(lat1, 0.0, lat2, lon2, s12);
            ++cnt;
            s += s12;
          }
        }
      }
      cout << cnt << " " << s << "\n";
      break;
    }
  }
  else if (accuracytest || coverage) {
    const Geodesic geod(a, f);
    const GeodesicExact geode(a, f);

    const GeographicLibL::Geodesic geodl(al, fl);
    const GeographicLibL::GeodesicExact geodel(al, fl);
    const unsigned NUMERR = 7;

    cout << fixed << setprecision(2);
    vector<Math::extended> erra(NUMERR);
    vector<Math::extended> err(NUMERR, 0.0);
    vector<unsigned> errind(NUMERR);
    unsigned cnt = 0;
    string s;
    while (getline(cin, s)) {
      istringstream str(s);
      Math::extended lat1l, lon1l, azi1l, lat2l, lon2l, azi2l,
        s12l, a12l, m12l, S12l;
      if (!(str >> lat1l >> lon1l >> azi1l
                >> lat2l >> lon2l >> azi2l
                >> s12l >> a12l >> m12l))
        break;
      if (!(str >> S12l))
        S12l = Math::NaN();
      if (coverage) {
#if defined(GEOD_DIAG) && GEOD_DIAG
        Math::real
          lat1 = lat1l, lon1 = lon1l,
          lat2 = lat2l, lon2 = lon2l,
          azi1, azi2, s12, m12;
        geod.Inverse(lat1, lon1, lat2, lon2, s12, azi1, azi2, m12);
        cout << geod.coverage << " " << geod.niter << "\n";
#endif
      } else {
        exact ?
          GeodError< Math::extended,
                     GeodesicExact, Math::real,
                     GeographicLibL::GeodesicExact, GeographicLibL::Math::real >
          (geode, geodel, lat1l, lon1l, azi1l,
           lat2l, lon2l, azi2l,
           s12l, a12l, m12l, S12l,
           erra) :
          GeodError< Math::extended,
                     Geodesic, Math::real,
                     GeographicLibL::Geodesic, GeographicLibL::Math::real >
          (geod, geodl, lat1l, lon1l, azi1l,
           lat2l, lon2l, azi2l,
           s12l, a12l, m12l, S12l,
           erra);
        for (unsigned i = 0; i < NUMERR; ++i) {
          if (Math::isfinite(err[i]) && !(erra[i] <= err[i])) {
            err[i] = erra[i];
            errind[i] = cnt;
          }
        }
        ++cnt;
      }
    }
    if (accuracytest) {
      Math::extended mult = sizeof(Math::real) > sizeof(double) ? 1e12l : 1e9l;
      for (unsigned i = 0; i < NUMERR; ++i)
        cout << i << " " << mult * err[i]
             << " " << errind[i] << "\n";
    }
  }
}

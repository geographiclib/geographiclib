/**
 * \file TMTest.cpp
 * \brief Command line utility for testing transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/TransverseMercator.hpp"
#include "GeographicLib/TransverseMercatorExact.hpp"
#include "GeographicLib/Geodesic.hpp"
#include "GeographicLib/Constants.hpp"
#include <vector>
#include <algorithm>

#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#define GEOTRANSTM 1
#define WMM 1

GeographicLib::Math::real
dist(GeographicLib::Math::real a, GeographicLib::Math::real r,
     long double lat0, long double lon0,
     GeographicLib::Math::real lat1, GeographicLib::Math::real lon1) {
  using namespace GeographicLib;
  typedef Math::real real;
  real
    phi = real(lat0) * Constants::degree(),
    f = r != 0 ? 1/r : 0,
    e2 = f * (2 - f),
    sinphi = sin(phi),
    n = 1/sqrt(1 - e2 * sinphi * sinphi),
      // See Wikipedia article on latitude
    hlon = std::cos(phi) * n,
    hlat = (1 - e2) * n * n * n;
  long double dlon = (long double)(lon1) - lon0;
  if (dlon >= 180) dlon -= 360;
  else if (dlon < -180) dlon += 360;
  return a * Constants::degree() *
    Math::hypot(real((long double)(lat1) - lat0) * hlat, real(dlon) * hlon);
}

int usage(int retval) {
  ( retval ? std::cerr : std::cout ) <<
"TMTest [-s] [-d]\n\
$Id$\n\
\n\
Read in TMcoords.dat on standard input and test TransverseMercatorExact\n\
or (if -s is given) TransverseMercator.  If -d dump the error for each\n\
point; otherwise summarise errors.  If -tf, perform a timing test of the\n\
forward projection.  If -tr,  perform a timing test of the reverse\n\
projection.\n";
  return retval;
}

#if GEOTRANSTM

/***************************************************************************/
/* RSC IDENTIFIER: TRANSVERSE MERCATOR
 *
 * ABSTRACT
 *
 *    This component provides conversions between Geodetic coordinates 
 *    (latitude and longitude) and Transverse Mercator projection coordinates
 *    (easting and northing).
 *
 * ERROR HANDLING
 *
 *    This component checks parameters for valid values.  If an invalid value
 *    is found the error code is combined with the current error code using 
 *    the bitwise or.  This combining allows multiple error codes to be
 *    returned. The possible error codes are:
 *
 *       TRANMERC_NO_ERROR           : No errors occurred in function
 *       TRANMERC_LAT_ERROR          : Latitude outside of valid range
 *                                      (-90 to 90 degrees)
 *       TRANMERC_LON_ERROR          : Longitude outside of valid range
 *                                      (-180 to 360 degrees, and within
 *                                        +/-90 of Central Meridian)
 *       TRANMERC_EASTING_ERROR      : Easting outside of valid range
 *                                      (depending on ellipsoid and
 *                                       projection parameters)
 *       TRANMERC_NORTHING_ERROR     : Northing outside of valid range
 *                                      (depending on ellipsoid and
 *                                       projection parameters)
 *       TRANMERC_ORIGIN_LAT_ERROR   : Origin latitude outside of valid range
 *                                      (-90 to 90 degrees)
 *       TRANMERC_CENT_MER_ERROR     : Central meridian outside of valid range
 *                                      (-180 to 360 degrees)
 *       TRANMERC_A_ERROR            : Semi-major axis less than or equal to zero
 *       TRANMERC_INV_F_ERROR        : Inverse flattening outside of valid range
 *								  	                  (250 to 350)
 *       TRANMERC_SCALE_FACTOR_ERROR : Scale factor outside of valid
 *                                     range (0.3 to 3.0)
 *		   TRANMERC_LON_WARNING        : Distortion will result if longitude is more
 *                                     than 9 degrees from the Central Meridian
 *
 * REUSE NOTES
 *
 *    TRANSVERSE MERCATOR is intended for reuse by any application that 
 *    performs a Transverse Mercator projection or its inverse.
 *    
 * REFERENCES
 *
 *    Further information on TRANSVERSE MERCATOR can be found in the 
 *    Reuse Manual.
 *
 *    TRANSVERSE MERCATOR originated from :  
 *                      U.S. Army Topographic Engineering Center
 *                      Geospatial Information Division
 *                      7701 Telegraph Road
 *                      Alexandria, VA  22310-3864
 *
 * LICENSES
 *
 *    None apply to this component.
 *
 * RESTRICTIONS
 *
 *    TRANSVERSE MERCATOR has no restrictions.
 *
 * ENVIRONMENT
 *
 *    TRANSVERSE MERCATOR was tested and certified in the following 
 *    environments:
 *
 *    1. Solaris 2.5 with GCC, version 2.8.1
 *    2. Windows 95 with MS Visual C++, version 6
 *
 * MODIFICATIONS
 *
 *    Date              Description
 *    ----              -----------
 *    2-26-07           Original C++ Code
 *
 */

namespace GeographicLib
{


  /***************************************************************************/
  /*
   *                              DEFINES
   */

  class GeotransTM
  {
  public:

    /*
     * The constructor receives the ellipsoid
     * parameters and Tranverse Mercator projection parameters as inputs, and
     * sets the corresponding state variables. If any errors occur, an exception 
     * is thrown with a description of the error.
     *
     *    ellipsoidSemiMajorAxis     : Semi-major axis of ellipsoid, in meters    (input)
     *    ellipsoidFlattening        : Flattening of ellipsoid						        (input)
     *    centralMeridian            : Longitude in radians at the center of the  (input)
     *                                 projection
     *    latitudeOfTrueScale        : Latitude in radians at the origin of the   (input)
     *                                 projection
     *    falseEasting               : Easting/X at the center of the projection  (input)
     *    falseNorthing              : Northing/Y at the center of the projection (input)
     *    scaleFactor                : Projection scale factor                    (input) 
     */

    GeotransTM( double ellipsoidSemiMajorAxis, double invFlattening, double scaleFactor );

    /*
     * The function convertFromGeodetic converts geodetic
     * (latitude and longitude) coordinates to Transverse Mercator projection
     * (easting and northing) coordinates, according to the current ellipsoid
     * and Transverse Mercator projection coordinates.  If any errors occur, 
     * an exception is thrown with a description of the error.
     *
     *    longitude     : Longitude in radians                        (input)
     *    latitude      : Latitude in radians                         (input)
     *    easting       : Easting/X in meters                         (output)
     *    northing      : Northing/Y in meters                        (output)
     */

    void Forward(double lon0, double latitude, double longitude, double& easting, double& northing) const;



    /*
     * The function convertToGeodetic converts Transverse
     * Mercator projection (easting and northing) coordinates to geodetic
     * (latitude and longitude) coordinates, according to the current ellipsoid
     * and Transverse Mercator projection parameters.  If any errors occur, 
     * an exception is thrown with a description of the error.
     *
     *    easting       : Easting/X in meters                         (input)
     *    northing      : Northing/Y in meters                        (input)
     *    longitude     : Longitude in radians                        (output)
     *    latitude      : Latitude in radians                         (output)
     */

    void Reverse(double lon0, double easting, double northing, double& latitude, double& longitude ) const;

  private:
    
    double semiMajorAxis;
    double flattening;

    /* Ellipsoid Parameters */
    double TranMerc_es;             /* Eccentricity squared */
    double TranMerc_ebs;            /* Second Eccentricity squared */

    double TranMerc_Scale_Factor;         /* Scale factor  */

    /* Isometric to geodetic latitude parameters */
    double TranMerc_ap;
    double TranMerc_bp;
    double TranMerc_cp;
    double TranMerc_dp;
    double TranMerc_ep;

    double sphtmd( double latitude ) const;
    double sphsn( double latitude ) const;
    double sphsr( double latitude ) const;

  };

  class WMMTM {
  private:
	double K0, Eps, Epssq, R4, R4oa, Bcoeff[10], Dcoeff[6];
    int NtrmsB, NtrmsD;
    static double ATanH(double x)
    {
      return(0.5 * std::log((1 + x) / (1 - x)));
    }
  public:
    WMMTM(double scale);
    void Forward(double lon0, double latitude, double longitude,
                 double& easting, double& northing) const {}
    void Reverse(double lon0, double easting, double northing,
                 double& latitude, double& longitude,
                 double& gamma, double& k) const;
    double MajorRadius() const
    { return Constants::WGS84_a();}
    double InverseFlattening() const
    { return Constants::WGS84_r();}
  };
}
// CLASSIFICATION: UNCLASSIFIED
#endif

int main(int argc, char* argv[]) {
  using namespace GeographicLib;
  typedef Math::real real;
  bool series = false;
#if GEOTRANSTM
  bool geotrans = false;
#endif
  bool dump = false;
  bool timefor = false, timerev = false;
  for (int m = 1; m < argc; ++m) {
    std::string arg(argv[m]);
    if (arg == "-s") {
      series = true;
#if GEOTRANSTM
      geotrans = false;
    } else if (arg == "-g") {
      series = false;
      geotrans = true;
#endif
    } else if (arg == "-d") {
      dump = true;
      timefor = false;
      timerev = false;
    } else if (arg == "-tf") {
      dump = false;
      timefor = true;
      timerev = false;
    } else if (arg == "-tr") {
      dump = false;
      timefor = false;
      timerev = true;
    } else
      return usage(arg != "-h");
  }

  if (timefor || timerev) {
    real s = 0;
    int count = 0;
    real dlat = 0.015, dlon = 0.015, dx = 2e3, dy = 2e3;
    if (series) {
      const TransverseMercator& tm = TransverseMercator::UTM;
      if (timefor) {
        real x, y, gam, k;
        for (real lat = -80.0; lat <= 84.0; lat += dlat)
          for (real lon = -3.0; lon <= 3.0; lon += dlon) {
            tm.Forward(0.0, lat, lon, x, y, gam, k);
            s += k;
            ++count;
          }
      } else {
        real lat, lon, gam, k;
        for (real x = -400e3; x <= 400e3; x += dx)
          for (real y = -9000e3; y <= 9500e3; y += dy) {
            tm.Reverse(0.0, x, y, lat, lon, gam, k);
            s += k;
            ++count;
          }
      }
#if GEOTRANSTM
    } else if (geotrans) {
#if WMM
      const WMMTM tm(Constants::UTM_k0());
#else
      const GeotransTM tm(Constants::WGS84_a(),
                          Constants::WGS84_r(),
                          Constants::UTM_k0());
#endif
      if (timefor) {
        real x, y;
        for (real lat = -80.0; lat <= 84.0; lat += dlat)
          for (real lon = -3.0; lon <= 3.0; lon += dlon) {
            tm.Forward(0.0, lat, lon, x, y);
            s += x;
            ++count;
          }
      } else {
        real lat, lon, gamma, k;
        for (real x = -400e3; x <= 400e3; x += dx)
          for (real y = -9000e3; y <= 9500e3; y += dy) {
#if WMM
            tm.Reverse(0.0, x, y, lat, lon, gamma, k);
#else
            tm.Reverse(0.0, x, y, lat, lon);
#endif
            s += lat;
            ++count;
          }
      }
#endif
    } else {
      const TransverseMercatorExact tm(Constants::WGS84_a(),
                                       Constants::WGS84_r(),
                                       Constants::UTM_k0(),
                                       true);
      if (timefor) {
        real x, y, gam, k;
        for (real lat = -80.0; lat <= 84.0; lat += dlat)
          for (real lon = -3.0; lon <= 3.0; lon += dlon) {
            tm.Forward(0.0, lat, lon, x, y, gam, k);
            s += k;
            ++count;
          }
      } else {
        real lat, lon, gam, k;
        for (real x = -400e3; x <= 400e3; x += dx)
          for (real y = -9000e3; y <= 9500e3; y += dy) {
            tm.Reverse(0.0, x, y, lat, lon, gam, k);
            s += k;
            ++count;
          }
      }
    }
    std::cout << count << " " << s << "\n";
    return 0;
  }

  try {
#if GEOTRANSTM
    real minlat = (series || geotrans) ? 0 : (dump ? -100 : -15);
#else
    real minlat = series ? 0 : (dump ? -100 : -15);
#endif
    const unsigned nbins = 101;
    std::vector<real> d(nbins);
    std::vector<real> errv(nbins, 0);
    std::vector<real> errvg(nbins, 0);
    std::vector<real> errvk(nbins, 0);
    real esterr = sizeof(real) == sizeof(double) ? (series ? 3e-9 : 8e-9) :
      (series ? 4e-12 : 4e-12);
    for (unsigned i = 0; i < nbins; ++i)
      d[i] = 100e3 * i;
    d[0] = 10e3;
    d[nbins - 1] = 10001966;
    const TransverseMercator& tm = TransverseMercator::UTM;
    const TransverseMercatorExact tme(Constants::WGS84_a(),
                                      Constants::WGS84_r(),
                                      Constants::UTM_k0(),
                                      true);
#if GEOTRANSTM
#if WMM
    const WMMTM tmg(tm.CentralScale());
#else
    const GeotransTM tmg(tm.MajorRadius(),
                         tm.InverseFlattening(),
                         tm.CentralScale());
#endif
    real
      a = (series || geotrans) ? tm.MajorRadius() : tme.MajorRadius(),
      r = (series || geotrans) ? tm.InverseFlattening() : tme.InverseFlattening();
#else
    real
      a = series ? tm.MajorRadius() : tme.MajorRadius(),
      r = series ? tm.InverseFlattening() : tme.InverseFlattening();
#endif
    const Geodesic geod(a, r);
    long double lat0l, lon0l, x0l, y0l, gam0l, k0l;
    while (std::cin >> lat0l >> lon0l >> x0l >> y0l >> gam0l >> k0l) {
      real
        lat0 = lat0l,
        lon0 = lon0l,
        x0 = x0l,
        y0 = y0l,
        gam0 = gam0l,
        k0 = k0l;
      if (lat0 < minlat)
        continue;
      real azi1, azi2, s12, m12;
      real errf, errgf, errkf, errr, errgr, errkr;
      geod.Inverse(std::max(lat0,real(0)), lon0, std::max(lat0,real(0)), -lon0,
                   s12, azi1, azi2, m12);
      s12 /= 2;
      real lat, lon, x, y, gam, k;
      if (series) {
        tm.Forward(0, lat0, lon0, x, y, gam, k);
#if GEOTRANSTM
      } else if (geotrans) {
        tmg.Forward(0, lat0, lon0, x, y);
        k = 1;
        gam = 0;
#endif
      } else
        tme.Forward(0, lat0, lon0, x, y, gam, k);
      errf = real(Math::hypot((long double)(x) - x0l,
                              (long double)(y) - y0l)) / k0;
      errgf = real(std::abs((long double)(gam) - gam0));
      errkf = real(std::abs((long double)(k) - k0));
      if (series) {
        tm.Reverse(0, x0, y0, lat, lon, gam, k);
#if GEOTRANSTM
      } else if (geotrans) {
#if WMM
        tmg.Reverse(0, x0, y0, lat, lon, gam, k);
        std::cerr << std::setprecision(16) << std::fixed << lat << " " << lon << "\n";
#else
        tmg.Reverse(0, x0, y0, lat, lon);
        k = 1;
        gam = 0;
#endif
#endif
      } else
        tme.Reverse(0, x0, y0, lat, lon, gam, k);
      errr = dist(a, r, lat0l, lon0l, lat, lon);
      errgr = real(std::abs((long double)(gam) - gam0));
      errkr = real(std::abs((long double)(k) - k0));

#if WMM
      //if (geotrans)
        errf = errgf = errkf = 0;
#endif
      real
        err = std::max(errf, errr),
        errg = std::max(errgf, errgr)
        - esterr/(a * std::sin((90 - lat0) * Constants::degree())
                  * Constants::degree()),
        errk = std::max(errkf, errkr) / k0;
      if (dump)
        std::cout << std::fixed << std::setprecision(12)
                  << lat0 << " " << lon0 << " "
                  << std::scientific << std::setprecision(4)
                  << errf << " " << errr << " "
                  << errgf << " " << errgr << " "
                  << errkf << " " << errkr << "\n";
      else
        for (unsigned i = 0; i < nbins; ++i) {
          if (s12 <= d[i]) {
            errv[i] = std::max(err, errv[i]);
            errvg[i] = std::max(errg, errvg[i]);
            errvk[i] = std::max(errk, errvk[i]);
          }
        }
    }
    if (!dump)
      for (unsigned i = 0; i < nbins; ++i)
        std::cout << int(d[i]/1000) << " "
                  << errv[i] << " "
                  << errvg[i] << " "
                  << errvk[i] << "\n";
  }
  catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

#if GEOTRANSTM
// CLASSIFICATION: UNCLASSIFIED

/***************************************************************************/
/* RSC IDENTIFIER: TRANSVERSE MERCATOR
 *
 * ABSTRACT
 *
 *    This component provides conversions between Geodetic coordinates 
 *    (latitude and longitude) and Transverse Mercator projection coordinates
 *    (easting and northing).
 *
 * ERROR HANDLING
 *
 *    This component checks parameters for valid values.  If an invalid value
 *    is found the error code is combined with the current error code using 
 *    the bitwise or.  This combining allows multiple error codes to be
 *    returned. The possible error codes are:
 *
 *       TRANMERC_NO_ERROR           : No errors occurred in function
 *       TRANMERC_LAT_ERROR          : Latitude outside of valid range
 *                                      (-90 to 90 degrees)
 *       TRANMERC_LON_ERROR          : Longitude outside of valid range
 *                                      (-180 to 360 degrees, and within
 *                                        +/-90 of Central Meridian)
 *       TRANMERC_EASTING_ERROR      : Easting outside of valid range
 *                                      (depending on ellipsoid and
 *                                       projection parameters)
 *       TRANMERC_NORTHING_ERROR     : Northing outside of valid range
 *                                      (depending on ellipsoid and
 *                                       projection parameters)
 *       TRANMERC_ORIGIN_LAT_ERROR   : Origin latitude outside of valid range
 *                                      (-90 to 90 degrees)
 *       TRANMERC_CENT_MER_ERROR     : Central meridian outside of valid range
 *                                      (-180 to 360 degrees)
 *       TRANMERC_A_ERROR            : Semi-major axis less than or equal to zero
 *       TRANMERC_INV_F_ERROR        : Inverse flattening outside of valid range
 *								  	                  (250 to 350)
 *       TRANMERC_SCALE_FACTOR_ERROR : Scale factor outside of valid
 *                                     range (0.3 to 3.0)
 *		   TRANMERC_LON_WARNING        : Distortion will result if longitude is more
 *                                     than 9 degrees from the Central Meridian
 *
 * REUSE NOTES
 *
 *    TRANSVERSE MERCATOR is intended for reuse by any application that 
 *    performs a Transverse Mercator projection or its inverse.
 *    
 * REFERENCES
 *
 *    Further information on TRANSVERSE MERCATOR can be found in the 
 *    Reuse Manual.
 *
 *    TRANSVERSE MERCATOR originated from :  
 *                      U.S. Army Topographic Engineering Center
 *                      Geospatial Information Division
 *                      7701 Telegraph Road
 *                      Alexandria, VA  22310-3864
 *
 * LICENSES
 *
 *    None apply to this component.
 *
 * RESTRICTIONS
 *
 *    TRANSVERSE MERCATOR has no restrictions.
 *
 * ENVIRONMENT
 *
 *    TRANSVERSE MERCATOR was tested and certified in the following 
 *    environments:
 *
 *    1. Solaris 2.5 with GCC, version 2.8.1
 *    2. Windows 95 with MS Visual C++, version 6
 *
 * MODIFICATIONS
 *
 *    Date              Description
 *    ----              -----------
 *    2-26-07           Original C++ Code
 *
 */


/***************************************************************************/
/*
 *                               INCLUDES
 */

/*
 *    math.h      - Standard C++ math library
 *    TransverseMercator.h  - Is for prototype error checking
 *    MapProjectionCoordinates.h   - defines map projection coordinates
 *    GeodeticCoordinates.h   - defines geodetic coordinates
 *    CoordinateConversionException.h - Exception handler
 *    ErrorMessages.h  - Contains exception messages
 */


namespace GeographicLib {

  using namespace std;

  /************************************************************************/
  /*                              FUNCTIONS     
   *
   */

  GeotransTM::GeotransTM( double ellipsoidSemiMajorAxis, double invFlattening, double scaleFactor )
  {
    /*
     * The constructor receives the ellipsoid
     * parameters and Tranverse Mercator projection parameters as inputs, and
     * sets the corresponding state variables. If any errors occur, an exception 
     * is thrown with a description of the error.
     *
     *    ellipsoidSemiMajorAxis     : Semi-major axis of ellipsoid, in meters    (input)
     *    ellipsoidFlattening        : Flattening of ellipsoid						        (input)
     *    centralMeridian            : Longitude in radians at the center of the  (input)
     *                                 projection
     *    latitudeOfTrueScale        : Latitude in radians at the origin of the   (input)
     *                                 projection
     *    falseEasting               : Easting/X at the center of the projection  (input)
     *    falseNorthing              : Northing/Y at the center of the projection (input)
     *    scaleFactor                : Projection scale factor                    (input) 
     */

    double tn;         /* True Meridianal distance constant  */
    double tn2;
    double tn3;
    double tn4;
    double tn5;
    double TranMerc_b; /* Semi-minor axis of ellipsoid, in meters */

    semiMajorAxis = ellipsoidSemiMajorAxis;
    flattening = 1/invFlattening;

    TranMerc_Scale_Factor = scaleFactor;

    /* Eccentricity Squared */
    TranMerc_es = 2 * flattening - flattening * flattening;
    /* Second Eccentricity Squared */
    TranMerc_ebs = (1 / (1 - TranMerc_es)) - 1;

    TranMerc_b = semiMajorAxis * (1 - flattening);    
    /*True meridianal constants  */
    tn = (semiMajorAxis - TranMerc_b) / (semiMajorAxis + TranMerc_b);
    tn2 = tn * tn;
    tn3 = tn2 * tn;
    tn4 = tn3 * tn;
    tn5 = tn4 * tn;

    TranMerc_ap = semiMajorAxis * (1.e0 - tn + 5.e0 * (tn2 - tn3)/4.e0
                                   + 81.e0 * (tn4 - tn5)/64.e0 );
    TranMerc_bp = 3.e0 * semiMajorAxis * (tn - tn2 + 7.e0 * (tn3 - tn4)
                                          /8.e0 + 55.e0 * tn5/64.e0 )/2.e0;
    TranMerc_cp = 15.e0 * semiMajorAxis * (tn2 - tn3 + 3.e0 * (tn4 - tn5 )/4.e0) /16.0;
    TranMerc_dp = 35.e0 * semiMajorAxis * (tn3 - tn4 + 11.e0 * tn5 / 16.e0) / 48.e0;
    TranMerc_ep = 315.e0 * semiMajorAxis * (tn4 - tn5) / 512.e0;
  }


  void GeotransTM::Forward(double lon0, double latitude, double longitude, double& easting, double& northing) const
  {
    /*
     * The function convertFromGeodetic converts geodetic
     * (latitude and longitude) coordinates to Transverse Mercator projection
     * (easting and northing) coordinates, according to the current ellipsoid
     * and Transverse Mercator projection coordinates.  If any errors occur, 
     * an exception is thrown with a description of the error.
     *
     *    longitude     : Longitude in radians                        (input)
     *    latitude      : Latitude in radians                         (input)
     *    easting       : Easting/X in meters                         (output)
     *    northing      : Northing/Y in meters                        (output)
     */

    double c;       /* Cosine of latitude                          */
    double c2;
    double c3;
    double c5;
    double c7;
    double dlam;    /* Delta longitude - Difference in Longitude       */
    double eta;     /* constant - TranMerc_ebs *c *c                   */
    double eta2;
    double eta3;
    double eta4;
    double s;       /* Sine of latitude                        */
    double sn;      /* Radius of curvature in the prime vertical       */
    double t;       /* Tangent of latitude                             */
    double tan2;
    double tan3;
    double tan4;
    double tan5;
    double tan6;
    double t1;      /* Term in coordinate conversion formula - GP to Y */
    double t2;      /* Term in coordinate conversion formula - GP to Y */
    double t3;      /* Term in coordinate conversion formula - GP to Y */
    double t4;      /* Term in coordinate conversion formula - GP to Y */
    double t5;      /* Term in coordinate conversion formula - GP to Y */
    double t6;      /* Term in coordinate conversion formula - GP to Y */
    double t7;      /* Term in coordinate conversion formula - GP to Y */
    double t8;      /* Term in coordinate conversion formula - GP to Y */
    double t9;      /* Term in coordinate conversion formula - GP to Y */
    double tmd;     /* True Meridianal distance                        */

    if (longitude - lon0 > 180)
      longitude -= lon0 - 360;
    else if (longitude - lon0 <= -180)
      longitude -= lon0 + 360;
    else
      longitude -= lon0;

    longitude *= Constants::degree();
    latitude *= Constants::degree();

    dlam = longitude;
    s = sin(latitude);
    c = cos(latitude);
    c2 = c * c;
    c3 = c2 * c;
    c5 = c3 * c2;
    c7 = c5 * c2;
    t = tan (latitude);
    tan2 = t * t;
    tan3 = tan2 * t;
    tan4 = tan3 * t;
    tan5 = tan4 * t;
    tan6 = tan5 * t;
    eta = TranMerc_ebs * c2;
    eta2 = eta * eta;
    eta3 = eta2 * eta;
    eta4 = eta3 * eta;

    /* radius of curvature in prime vertical */
    sn = sphsn(latitude);

    /* True Meridianal Distances */
    tmd = sphtmd(latitude);

    /* northing */
    t1 = (tmd) * TranMerc_Scale_Factor;
    t2 = sn * s * c * TranMerc_Scale_Factor/ 2.e0;
    t3 = sn * s * c3 * TranMerc_Scale_Factor * (5.e0 - tan2 + 9.e0 * eta 
                                                + 4.e0 * eta2) /24.e0; 

    t4 = sn * s * c5 * TranMerc_Scale_Factor * (61.e0 - 58.e0 * tan2
                                                + tan4 + 270.e0 * eta - 330.e0 * tan2 * eta + 445.e0 * eta2
                                                + 324.e0 * eta3 -680.e0 * tan2 * eta2 + 88.e0 * eta4 
                                                -600.e0 * tan2 * eta3 - 192.e0 * tan2 * eta4) / 720.e0;

    t5 = sn * s * c7 * TranMerc_Scale_Factor * (1385.e0 - 3111.e0 * 
                                                tan2 + 543.e0 * tan4 - tan6) / 40320.e0;

    northing = t1 + pow(dlam, 2.e0) * t2
      + pow(dlam,4.e0) * t3 + pow(dlam,6.e0) * t4
      + pow(dlam,8.e0) * t5; 

    /* Easting */
    t6 = sn * c * TranMerc_Scale_Factor;
    t7 = sn * c3 * TranMerc_Scale_Factor * (1.e0 - tan2 + eta ) /6.e0;
    t8 = sn * c5 * TranMerc_Scale_Factor * (5.e0 - 18.e0 * tan2 + tan4
                                            + 14.e0 * eta - 58.e0 * tan2 * eta + 13.e0 * eta2 + 4.e0 * eta3 
                                            - 64.e0 * tan2 * eta2 - 24.e0 * tan2 * eta3 )/ 120.e0;
    t9 = sn * c7 * TranMerc_Scale_Factor * ( 61.e0 - 479.e0 * tan2
                                             + 179.e0 * tan4 - tan6 ) /5040.e0;

    easting = dlam * t6 + pow(dlam,3.e0) * t7 
      + pow(dlam,5.e0) * t8 + pow(dlam,7.e0) * t9;

  }


  void GeotransTM::Reverse(double lon0, double easting, double northing, double& latitude, double& longitude ) const
  {
    /*
     * The function convertToGeodetic converts Transverse
     * Mercator projection (easting and northing) coordinates to geodetic
     * (latitude and longitude) coordinates, according to the current ellipsoid
     * and Transverse Mercator projection parameters.  If any errors occur, 
     * an exception is thrown with a description of the error.
     *
     *    easting       : Easting/X in meters                         (input)
     *    northing      : Northing/Y in meters                        (input)
     *    longitude     : Longitude in radians                        (output)
     *    latitude      : Latitude in radians                         (output)
     */

    double c;       /* Cosine of latitude                          */
    double de;      /* Delta easting - Difference in Easting (easting-Fe)    */
    double dlam;    /* Delta longitude - Difference in Longitude       */
    double eta;     /* constant - TranMerc_ebs *c *c                   */
    double eta2;
    double eta3;
    double eta4;
    double ftphi;   /* Footpoint latitude                              */
    int    i;       /* Loop iterator                   */
    double s;       /* Sine of latitude                        */
    double sn;      /* Radius of curvature in the prime vertical       */
    double sr;      /* Radius of curvature in the meridian             */
    double t;       /* Tangent of latitude                             */
    double tan2;
    double tan4;
    double t10;     /* Term in coordinate conversion formula - GP to Y */
    double t11;     /* Term in coordinate conversion formula - GP to Y */
    double t12;     /* Term in coordinate conversion formula - GP to Y */
    double t13;     /* Term in coordinate conversion formula - GP to Y */
    double t14;     /* Term in coordinate conversion formula - GP to Y */
    double t15;     /* Term in coordinate conversion formula - GP to Y */
    double t16;     /* Term in coordinate conversion formula - GP to Y */
    double t17;     /* Term in coordinate conversion formula - GP to Y */
    double tmd;     /* True Meridianal distance                        */

    /*  Origin  */
    tmd = northing / TranMerc_Scale_Factor; 

    /* First Estimate */
    sr = sphsr(0.e0);
    ftphi = tmd/sr;

    for (i = 0; i < 5 ; i++)
      {
        t10 = sphtmd (ftphi);
        sr = sphsr(ftphi);
        ftphi = ftphi + (tmd - t10) / sr;
      }

    /* Radius of Curvature in the meridian */
    sr = sphsr(ftphi);

    /* Radius of Curvature in the meridian */
    sn = sphsn(ftphi);

    /* Sine Cosine terms */
    s = sin(ftphi);
    c = cos(ftphi);

    /* Tangent Value  */
    t = tan(ftphi);
    tan2 = t * t;
    tan4 = tan2 * tan2;
    eta = TranMerc_ebs * pow(c, 2.0);
    eta2 = eta * eta;
    eta3 = eta2 * eta;
    eta4 = eta3 * eta;
    de = easting;
    if (fabs(de) < 0.0001)
      de = 0.0;

    /* Latitude */
    t10 = t / (2.e0 * sr * sn * pow(TranMerc_Scale_Factor, 2.0));
    t11 = t * (5.e0  + 3.e0 * tan2 + eta - 4.e0 * pow(eta, 2.0)
               - 9.e0 * tan2 * eta) / (24.e0 * sr * pow(sn, 3.0) 
                                       * pow(TranMerc_Scale_Factor, 4.0));
    t12 = t * (61.e0 + 90.e0 * tan2 + 46.e0 * eta + 45.E0 * tan4
               - 252.e0 * tan2 * eta  - 3.e0 * eta2 + 100.e0 
               * eta3 - 66.e0 * tan2 * eta2 - 90.e0 * tan4
               * eta + 88.e0 * eta4 + 225.e0 * tan4 * eta2
               + 84.e0 * tan2* eta3 - 192.e0 * tan2 * eta4)
      / ( 720.e0 * sr * pow(sn,5.0) * pow(TranMerc_Scale_Factor, 6.0) );
    t13 = t * ( 1385.e0 + 3633.e0 * tan2 + 4095.e0 * tan4 + 1575.e0 
                * pow(t, 6.0))/ (40320.e0 * sr * pow(sn, 7.0) * pow(TranMerc_Scale_Factor, 8.0));
    latitude = ftphi - pow(de, 2.0) * t10 + pow(de, 4.0) * t11 - pow(de, 6.0) * t12 
      + pow(de, 8.0) * t13;

    t14 = 1.e0 / (sn * c * TranMerc_Scale_Factor);

    t15 = (1.e0 + 2.e0 * tan2 + eta) / (6.e0 * pow(sn, 3.0) * c * 
                                        pow(TranMerc_Scale_Factor, 3.0));

    t16 = (5.e0 + 6.e0 * eta + 28.e0 * tan2 - 3.e0 * eta2
           + 8.e0 * tan2 * eta + 24.e0 * tan4 - 4.e0 
           * eta3 + 4.e0 * tan2 * eta2 + 24.e0 
           * tan2 * eta3) / (120.e0 * pow(sn, 5.0) * c  
                             * pow(TranMerc_Scale_Factor, 5.0));

    t17 = (61.e0 +  662.e0 * tan2 + 1320.e0 * tan4 + 720.e0 
           * pow(t, 6.0)) / (5040.e0 * pow(sn, 7.0) * c 
                             * pow(TranMerc_Scale_Factor, 7.0));

    /* Difference in Longitude */
    dlam = de * t14 - pow(de, 3.0) * t15 + pow(de, 5.0) * t16 - pow(de, 7.0) * t17;

    /* Longitude */
    longitude = dlam;

    latitude /= Constants::degree();
    longitude /= Constants::degree();
    if (longitude + lon0 >= 180)
      longitude += lon0 - 360;
    else if (longitude + lon0 < -180)
      longitude += lon0 + 360;
    else
      longitude += lon0;

  }


  /************************************************************************/
  /*                              PRIVATE FUNCTIONS     
   *
   */

  double GeotransTM::sphtmd( double latitude ) const
  {
    return TranMerc_ap * latitude 
      - TranMerc_bp * sin(2.e0 * latitude) + TranMerc_cp * sin(4.e0 * latitude) 
      - TranMerc_dp * sin(6.e0 * latitude) + TranMerc_ep * sin(8.e0 * latitude);
  }


  double GeotransTM::sphsn( double  latitude ) const
  {
    return semiMajorAxis / sqrt( 1.e0 - TranMerc_es * pow(sin(latitude), 2.0));
  }


  double GeotransTM::sphsr( double latitude ) const
  {
    double denom = sqrt(1.e0 - TranMerc_es * pow(sin(latitude), 2.0));
    return semiMajorAxis * (1.e0 - TranMerc_es) / pow(denom, 3.0);
  }

  WMMTM::WMMTM( /* double ellipsoidSemiMajorAxis, double invFlattening, */
               double scaleFactor )  {
    K0=scaleFactor;
	
    //WGS84 ellipsoid
	Eps=0.081819190842621494335;
    Epssq=0.0066943799901413169961;
	R4=6367449.1458234153093;
    R4oa=0.99832429843125277950;

    Bcoeff[0]=-8.3773216405794867707E-04;
    Bcoeff[1]=-5.9058701522203651815E-08;    
    Bcoeff[2]=-1.6734826653438249274E-10;   
    Bcoeff[3]=-2.1647981104903861797E-13;   
    Bcoeff[4]=-3.7879309688396013649E-16;   
    Bcoeff[5]=-7.2367692879669040169E-19;   
    Bcoeff[6]=-1.4934544948474238627E-21;   
    Bcoeff[7]=-3.2538430893102429751E-24;   
    Bcoeff[8]=-7.3912479800652263027E-27;   
    Bcoeff[9]=-1.7344445906774504295E-29;  

    Dcoeff[0]=3.3565514691328321888E-03;
    Dcoeff[1]=6.5718731986276970786E-06;
    Dcoeff[2]=1.7646404113343270872E-08;
    Dcoeff[3]=5.3877540683992094349E-11;
    Dcoeff[4]=1.7639829049873024993E-13;
    Dcoeff[5]=6.0345363280209865351E-16;
    
    NtrmsB = NtrmsD = 6;
  }

  void WMMTM::Reverse(double lon0, double easting, double northing,
                      double& latitude, double& longitude,
                      double& gamma, double& k) const {
    double Xstar, Ystar, TwoX, TwoY;
    int K;
    double c2kx[8], s2kx[8], c2ky[8], s2ky[8];
    double U, V;
    double CU, SU, CV, SV;
    double Lam;
    double Schi, Chi, TwoChi;

    /*   These are needed for point-scale and CoM  */

    double sig3, sig4;
    double Sphi;
    double W, P, denom;


    /*  Undo offsets, scale change, and factor R4
        ---- -------  ----- ------  --- ------ --  */

    Xstar = easting /  (K0 * R4);
    Ystar = northing / (K0 * R4);


    /*  Trigonometric multiple angles
        ------------- -------- ------  */

    TwoX = 2 * Xstar;
    TwoY = 2 * Ystar;

    for (K = 0; K < NtrmsB; K++)
      {
        c2kx[K] = cosh((K+1) * TwoX);
        s2kx[K] = sinh((K+1) * TwoX);
        c2ky[K] = cos((K+1) * TwoY);
        s2ky[K] = sin((K+1) * TwoY);
      }

    /*  Second plane (x*, y*) to first plane (u, v)
        ------ ----- -------- -- ----- ----- ------  */

    U = 0;
    V = 0;

    for (K = NtrmsB - 1; K >= 0; K--)
      {
        U += Bcoeff[K] * s2kx[K] * c2ky[K];
        V += Bcoeff[K] * c2kx[K] * s2ky[K];
      }

    U += Xstar;
    V += Ystar;


    /*  First plane to sphere
        ----- ----- -- ------  */

    CU = cosh(U);
    SU = sinh(U);
    CV = cos(V);
    SV = sin(V);

    /*   Longitude from central meridian  */

    /*
    if ((abs(CV) < 10E-12) && (abs(SU) < 10E-12))
      Lam = 0;
    else
    */
      Lam = atan2(SU, CV);

    /*   Conformal latitude  */

    Schi = SV / CU;

    if (true) /* (abs(Schi) < 0.984375000L) */
      Chi = std::asin(Schi);
    else
      Chi = std::atan2(SV, sqrt((SU * SU) + (CV * CV)));
    /*  Sphere to ellipsoid
        ------ -- ---------  */

    TwoChi = 2 * Chi;
    double Phi = 0;

    for (K = NtrmsD - 1; K >= 0; K--)
      Phi += Dcoeff[K] * sin((K+1) * TwoChi);

    Phi += Chi;
    latitude = Phi;
    longitude = Lam;


    latitude /= Constants::degree();
    longitude /= Constants::degree();
    if (longitude + lon0 >= 180)
      longitude += lon0 - 360;
    else if (longitude + lon0 < -180)
      longitude += lon0 + 360;
    else
      longitude += lon0;

    /*  Point-scale and CoM
        ----------- --- ---  */

    sig3 = 0;
    sig4 = 0;

    for (K = NtrmsB - 1; K >= 0; K--)
      {
        sig3 += 2 * (K+1) * Bcoeff[K] * c2kx[K] * c2ky[K];
        sig4 += 2 * (K+1) * Bcoeff[K] * s2kx[K] * s2ky[K];
      }

    sig3++;

    Sphi = sin(Phi);
    W = sqrt(1 - (Epssq * Sphi * Sphi));
    P = exp(Eps * ATanH(Eps * Sphi));
    denom = (1 + Sphi) / P + (1 - Sphi) * P;

    k = K0 * R4oa * W * 2 / denom * CU /
      sqrt((sig3 * sig3) + (sig4 * sig4));

    gamma = atan2(Schi * SU, CV) - atan2(sig4, sig3);
    gamma /= Constants::degree();

  }

}
#endif

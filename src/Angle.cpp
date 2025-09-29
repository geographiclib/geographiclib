/**
 * \file Angle.cpp
 * \brief Implementation for the GeographicLib::Angle class.
 *
 * This file is an implementation of the methods described in
 * - C. F. F. Karney,
 *   <a href="https://doi.org/10.1080/00396265.2023.2217604">
 *   On auxiliary latitudes,</a>
 *   Survey Review 56(395), 165--180 (2024);
 *   preprint
 *   <a href="https://arxiv.org/abs/2212.05818">arXiv:2212.05818</a>.
 * .
 * Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/Angle.hpp>
#include <GeographicLib/DMS.hpp>
#include <iostream>

namespace GeographicLib {

  using namespace std;

  template<typename T>
  T AngleT<T>::rnd(T x) {
    // This value of z more-or-less matches the value z = 1/16 in
    // Math::AngRound (where the argument is in degrees).
    static const T z = 1/T(1024);
    GEOGRAPHICLIB_VOLATILE T y = fabs(x);
    GEOGRAPHICLIB_VOLATILE T w = z - y;
    // The compiler mustn't "simplify" z - (z - y) to y
    y = w > 0 ? z - w : y;
    return copysign(y, x);
  }

  template<typename T>
  void AngleT<T>::DecodeLatLon(const string& stra, const string& strb,
                           AngleT<T>& lat, AngleT<T>& lon, bool longfirst) {
    T a, b;
    DMS::flag ia, ib;
    a = T(DMS::Decode(stra, ia));
    b = T(DMS::Decode(strb, ib));
    if (ia == DMS::NONE && ib == DMS::NONE) {
      // Default to lat, long unless longfirst
      ia = longfirst ? DMS::LONGITUDE : DMS::LATITUDE;
      ib = longfirst ? DMS::LATITUDE : DMS::LONGITUDE;
    } else if (ia == DMS::NONE)
      ia = DMS::flag(DMS::LATITUDE + DMS::LONGITUDE - ib);
    else if (ib == DMS::NONE)
      ib = DMS::flag(DMS::LATITUDE + DMS::LONGITUDE - ia);
    if (ia == ib)
      throw GeographicErr("Both " + stra + " and "
                          + strb + " interpreted as "
                          + (ia == DMS::LATITUDE ? "latitudes" : "longitudes"));
    lat = AngleT<T>(ia == DMS::LATITUDE ? a : b);
    lon = AngleT<T>(ia == DMS::LATITUDE ? b : a);
  }

  template<typename T>
  AngleT<T> AngleT<T>::DecodeAzimuth(const string& azistr) {
    DMS::flag ind;
    T azi = T(DMS::Decode(azistr, ind));
    if (ind == DMS::LATITUDE)
      throw GeographicErr("Azimuth " + azistr
                          + " has a latitude hemisphere, N/S");
    return AngleT<T>(azi);
  }

  template<typename T>
  string AngleT<T>::LatLonString(AngleT<T> lat, AngleT<T> lon,
                             int prec, bool dms, char dmssep, bool longfirst) {
    string
      latstr = dms ? DMS::Encode(Math::real(T(lat)),
                                 prec, DMS::LATITUDE, dmssep) :
      DMS::Encode(Math::real(T(lat)), prec, DMS::NUMBER),
      lonstr = dms ? DMS::Encode(Math::real(T(lon)),
                                 prec, DMS::LONGITUDE, dmssep) :
      DMS::Encode(Math::real(T(lon)), prec, DMS::NUMBER);
    return
      (longfirst ? lonstr : latstr) + " " + (longfirst ? latstr : lonstr);
  }

  template<typename T>
  string AngleT<T>::AzimuthString(AngleT<T> azi, int prec, bool dms, char dmssep) {
    return dms ? DMS::Encode(Math::real(T(azi)), prec, DMS::AZIMUTH, dmssep) :
      DMS::Encode(Math::real(T(azi)), prec, DMS::NUMBER);
  }

#define GEOGRAPHICLIB_ANGLE_INSTANTIATE(T)                           \
  template T         GEOGRAPHICLIB_EXPORT AngleT<T>::rnd(T);         \
  template void      GEOGRAPHICLIB_EXPORT AngleT<T>::DecodeLatLon    \
       (const string& , const string&,AngleT<T>&, AngleT<T>&, bool); \
  template AngleT<T> GEOGRAPHICLIB_EXPORT AngleT<T>::DecodeAzimuth   \
       (const string&);                                              \
  template string    GEOGRAPHICLIB_EXPORT AngleT<T>::LatLonString    \
       (AngleT<T>, AngleT<T>, int, bool, char, bool);                \
  template string    GEOGRAPHICLIB_EXPORT AngleT<T>::AzimuthString   \
       (AngleT<T>, int, bool, char);

  // Instantiate with the standard floating type
  GEOGRAPHICLIB_ANGLE_INSTANTIATE(float)
  GEOGRAPHICLIB_ANGLE_INSTANTIATE(double)
#if GEOGRAPHICLIB_HAVE_LONG_DOUBLE
  // Instantiate if long double is distinct from double
  GEOGRAPHICLIB_ANGLE_INSTANTIATE(long double)
#endif
#if GEOGRAPHICLIB_PRECISION > 3
  // Instantiate with the high precision type
  GEOGRAPHICLIB_ANGLE_INSTANTIATE(Math::real)
#endif

} // namespace GeographicLib

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

  Math::real Angle::rnd(real x) {
    // This value of z more-or-less matches the value z = 1/16 in
    // Math::AngRound (where the argument is in degrees).
    static const real z = 1/real(1024);
    GEOGRAPHICLIB_VOLATILE real y = fabs(x);
    GEOGRAPHICLIB_VOLATILE real w = z - y;
    // The compiler mustn't "simplify" z - (z - y) to y
    y = w > 0 ? z - w : y;
    return copysign(y, x);
  }

  void Angle::DecodeLatLon(const string& stra, const string& strb,
                           Angle& lat, Angle& lon, bool longfirst) {
    real a, b;
    DMS::flag ia, ib;
    a = DMS::Decode(stra, ia);
    b = DMS::Decode(strb, ib);
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
    lat = Angle(ia == DMS::LATITUDE ? a : b);
    lon = Angle(ia == DMS::LATITUDE ? b : a);
  }

  Angle Angle::DecodeAzimuth(const string& azistr) {
    DMS::flag ind;
    real azi = DMS::Decode(azistr, ind);
    if (ind == DMS::LATITUDE)
      throw GeographicErr("Azimuth " + azistr
                          + " has a latitude hemisphere, N/S");
    return Angle(azi);
  }

  string Angle::LatLonString(Angle lat, Angle lon,
                             int prec, bool dms, char dmssep, bool longfirst) {
    string
      latstr = dms ? DMS::Encode(real(lat), prec, DMS::LATITUDE, dmssep) :
      DMS::Encode(real(lat), prec, DMS::NUMBER),
      lonstr = dms ? DMS::Encode(real(lon), prec, DMS::LONGITUDE, dmssep) :
      DMS::Encode(real(lon), prec, DMS::NUMBER);
    return
      (longfirst ? lonstr : latstr) + " " + (longfirst ? latstr : lonstr);
  }

  std::string Angle::AzimuthString(Angle azi, int prec, bool dms, char dmssep) {
    return dms ? DMS::Encode(real(azi), prec, DMS::AZIMUTH, dmssep) :
      DMS::Encode(real(azi), prec, DMS::NUMBER);
  }

} // namespace GeographicLib

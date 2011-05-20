/**
 * \file DMS.cpp
 * \brief Implementation for GeographicLib::DMS class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include "GeographicLib/DMS.hpp"

#define GEOGRAPHICLIB_DMS_CPP "$Id: DMS.cpp 6831 2010-05-31 20:51:00Z karney $"

RCSID_DECL(GEOGRAPHICLIB_DMS_CPP)
RCSID_DECL(GEOGRAPHICLIB_DMS_HPP)

namespace GeographicLib {

  using namespace std;

  const string DMS::hemispheres = "SNWE";
  const string DMS::signs = "-+";
  const string DMS::digits = "0123456789";
  const string DMS::dmsindicators = "D'\"";
  const string DMS::components[] = {"degrees", "minutes", "seconds"};

  Math::real DMS::Decode(const std::string& dms, flag& ind) {
    int sign = 1;
    unsigned
      beg = 0,
      end = unsigned(dms.size());
    while (beg < end && isspace(dms[beg]))
      ++beg;
    while (beg < end && isspace(dms[end - 1]))
      --end;
    flag ind1 = NONE;
    int k = -1;
    if (end > beg && (k = lookup(hemispheres, dms[beg])) >= 0) {
      ind1 = (k / 2) ? LONGITUDE : LATITUDE;
      sign = k % 2 ? 1 : -1;
      ++beg;
    }
    if (end > beg && (k = lookup(hemispheres, dms[end-1])) >= 0) {
      if (k >= 0) {
        if (ind1 != NONE) {
          if (toupper(dms[beg - 1]) == toupper(dms[end - 1]))
            throw GeographicErr("Repeated hemisphere indicators "
                                + str(dms[beg - 1]) + " in "
                                + dms.substr(beg - 1, end - beg + 1));
          else
            throw GeographicErr("Contradictory hemisphere indicators "
                                + str(dms[beg - 1]) + " and "
                                + str(dms[end - 1]) + " in "
                                + dms.substr(beg - 1, end - beg + 1));
        }
        ind1 = (k / 2) ? LONGITUDE : LATITUDE;
        sign = k % 2 ? 1 : -1;
        --end;
      }
    }
    if (end > beg && (k = lookup(signs, dms[beg])) >= 0) {
      if (k >= 0) {
        sign *= k ? 1 : -1;
        ++beg;
      }
    }
    if (end == beg)
      throw GeographicErr("Empty or incomplete DMS string " + dms);
    real ipieces[] = {0, 0, 0};
    real fpieces[] = {0, 0, 0};
    unsigned npiece = 0;
    real icurrent = 0;
    real fcurrent = 0;
    unsigned ncurrent = 0, p = beg;
    bool pointseen = false;
    unsigned digcount = 0;
    while (p < end) {
      char x = dms[p++];
      if ((k = lookup(digits, x)) >= 0) {
        ++ncurrent;
        if (digcount > 0)
          ++digcount;           // Count of decimal digits
        else
          icurrent = 10 * icurrent + k;
      } else if (x == '.') {
        if (pointseen)
          throw GeographicErr("Multiple decimal points in "
                              + dms.substr(beg, end - beg));
        pointseen = true;
        digcount = 1;
      } else if ((k = lookup(dmsindicators, x)) >= 0) {
        if (unsigned(k) == npiece - 1)
          throw GeographicErr("Repeated " + components[k]
                              + " component in " + dms.substr(beg, end - beg));
        else if (unsigned(k) < npiece)
          throw GeographicErr(components[k] + " component follows "
                              + components[npiece - 1] + " component in "
                              + dms.substr(beg, end - beg));
        if (ncurrent == 0)
          throw GeographicErr("Missing numbers in " + components[k]
                              + " component of " + dms.substr(beg, end - beg));
        if (digcount > 1) {
          istringstream s(dms.substr(p - digcount - 1, digcount));
          s >> fcurrent;
        }
        ipieces[k] = icurrent;
        fpieces[k] = icurrent + fcurrent;
        if (p < end) {
          npiece = k + 1;
          icurrent = fcurrent = 0;
          ncurrent = digcount = 0;
        }
      } else if (lookup(signs, x) >= 0)
        throw GeographicErr("Internal sign in DMS string "
                            + dms.substr(beg, end - beg));
      else
        throw GeographicErr("Illegal character " + str(x)
                            + " in DMS string "
                            + dms.substr(beg, end - beg));
    }
    if (lookup(dmsindicators, dms[p - 1]) < 0) {
      if (npiece >= 3)
        throw GeographicErr("Extra text following seconds in DMS string "
                            + dms.substr(beg, end - beg));
      if (ncurrent == 0)
        throw GeographicErr("Missing numbers in " + components[k]
                            + " component of " + dms.substr(beg, end - beg));
      if (digcount > 1) {
        istringstream s(dms.substr(p - digcount, digcount));
        s >> fcurrent;
      }
      ipieces[npiece] = icurrent;
      fpieces[npiece] = icurrent + fcurrent;
    }
    if (pointseen && digcount == 0)
      throw GeographicErr("Decimal point in non-terminal component of "
                          + dms.substr(beg, end - beg));
    // Note that we accept 59.999999... even though it rounds to 60.
    if (ipieces[1] >= 60)
      throw GeographicErr("Minutes " + str(fpieces[1])
                          + " not in range [0, 60)");
    if (ipieces[2] >= 60)
      throw GeographicErr("Seconds " + str(fpieces[2])
                          + " not in range [0, 60)");
    ind = ind1;
    // Assume check on range of result is made by calling routine (which might
    // be able to offer a better diagnostic).
    return real(sign) * (fpieces[0] + (fpieces[1] + fpieces[2] / 60) / 60);
  }

  Math::real DMS::Decode(const std::string& str) {
    istringstream is(str);
    real num;
    if (!(is >> num))
      throw GeographicErr("Could not read number: " + str);
    // On some platforms, is >> num gobbles final E in 1234E, so look for last
    // character which is legal as the final character in a number (digit or
    // period).
    int pos = min(int(is.tellg()), int(str.find_last_of("0123456789.")) + 1);
    if (pos != int(str.size()))
      throw GeographicErr("Extra text " + str.substr(pos) +
                          " in number " + str);
    return num;
  }

  void DMS::DecodeLatLon(const std::string& stra, const std::string& strb,
                         real& lat, real& lon) {
      real a, b;
      flag ia, ib;
      a = Decode(stra, ia);
      b = Decode(strb, ib);
      if (ia == NONE && ib == NONE) {
        // Default to lat, long
        ia = LATITUDE;
        ib = LONGITUDE;
      } else if (ia == NONE)
        ia = flag(LATITUDE + LONGITUDE - ib);
      else if (ib == NONE)
        ib = flag(LATITUDE + LONGITUDE - ia);
      if (ia == ib)
        throw GeographicErr("Both " + stra + " and "
                            + strb + " interpreted as "
                            + (ia == LATITUDE ? "latitudes" : "longitudes"));
      real
        lat1 = ia == LATITUDE ? a : b,
        lon1 = ia == LATITUDE ? b : a;
      if (! (lat1 >= -90 && lat1 <= 90))
        throw GeographicErr("Latitude " + str(lat1) + "d not in [-90d, 90d]");
      if (! (lon1 >= -180 && lon1 <= 360))
        throw GeographicErr("Latitude " + str(lon1)
                            + "d not in [-180d, 360d]");
      if (lon1 >= 180)
        lon1 -= 360;
      lat = lat1;
      lon = lon1;
  }

  Math::real DMS::DecodeAngle(const std::string& angstr) {
    DMS::flag ind;
    real ang = Decode(angstr, ind);
    if (ind != DMS::NONE)
      throw GeographicErr("Arc angle " + angstr
                          + " includes a hemisphere, N/E/W/S");
    return ang;
  }

  Math::real DMS::DecodeAzimuth(const std::string& azistr) {
    DMS::flag ind;
    real azi = Decode(azistr, ind);
    if (ind == DMS::LATITUDE)
      throw GeographicErr("Azimuth " + azistr
                          + " has a latitude hemisphere, N/S");
    if (!(azi >= -180 && azi <= 360))
      throw GeographicErr("Azimuth " + azistr + " not in range [-180,360]");
    if (azi >= 180) azi -= 360;
    return azi;
  }

  string DMS::Encode(real angle, component trailing, unsigned prec, flag ind) {
    // Assume check on range of input angle has been made by calling
    // routine (which might be able to offer a better diagnostic).
    //
    // 15 - 2 * trailing = ceiling(log10(2^53/90/60^trailing)).
    // This suffices to give full real precision for numbers in [-90,90]
    prec = min(15 - 2 * unsigned(trailing), prec);
    real scale = 1;
    for (unsigned i = 0; i < unsigned(trailing); ++i)
      scale *= 60;
    for (unsigned i = 0; i < prec; ++i)
      scale *= 10;
    if (ind == AZIMUTH)
      angle -= floor(angle/360) * 360;
    int sign = angle < 0 ? -1 : 1;
    angle *= sign;

    // Break off integer part to preserve precision in manipulation of
    // fractional part.
    real
      idegree = floor(angle),
      fdegree = floor((angle - idegree) * scale + real(0.5)) / scale;
    if (fdegree >= 1) {
      idegree += 1;
      fdegree -= 1;
    }
    real pieces[3] = {fdegree, 0, 0};
    for (unsigned i = 1; i <= unsigned(trailing); ++i) {
      real
        ip = floor(pieces[i - 1]),
        fp = pieces[i - 1] - ip;
      pieces[i] = fp * 60;
      pieces[i - 1] = ip;
    }
    pieces[0] += idegree;
    ostringstream s;
    s << fixed  << setfill('0');
    if (ind == NONE && sign < 0)
      s << '-';
    switch (trailing) {
    case DEGREE:
      if (ind != NONE)
        s << setw(1 + min(int(ind), 2) + prec + (prec ? 1 : 0));
      s << setprecision(prec) << pieces[0];
      // Don't include degree designator (d) if it is the trailing component.
      break;
    default:
      if (ind != NONE)
        s << setw(1 + min(int(ind), 2));
      s << setprecision(0) << pieces[0] << char(tolower(dmsindicators[0]));
      switch (trailing) {
      case MINUTE:
        s << setw(2 + prec + (prec ? 1 : 0)) << setprecision(prec)
          << pieces[1] <<  char(tolower(dmsindicators[1]));
        break;
      case SECOND:
        s << setw(2) << pieces[1] <<  char(tolower(dmsindicators[1]))
          << setw(2 + prec + (prec ? 1 : 0)) << setprecision(prec)
          << pieces[2] <<  char(tolower(dmsindicators[2]));
        break;
      default:
        break;
      }
    }
    if (ind != NONE && ind != AZIMUTH)
      s << hemispheres[(ind == LATITUDE ? 0 : 2) + (sign < 0 ? 0 : 1)];
    return s.str();
  }

} // namespace GeographicLib

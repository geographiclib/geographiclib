/**
 * \file DMS.cpp
 * \brief Implementation for GeographicLib::DMS class
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#include "GeographicLib/DMS.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iomanip>

namespace {
  char RCSID[] = "$Id$";
  char RCSID_H[] = DMS_HPP;
}

namespace GeographicLib {

  const std::string DMS::hemispheres = "SNWE";
  const std::string DMS::signs = "-+";
  const std::string DMS::digits = "0123456789";
  const std::string DMS::dmsindicators = "D'\"";
  const std::string DMS::components[] = {"degrees", "minutes", "seconds"};

  double DMS::Decode(const std::string& dms, flag& ind) {
    double sign = 1;
    unsigned
      beg = 0,
      end = unsigned(dms.size());
    while (beg < end && isspace(dms[beg]))
      ++beg;
    while (beg < end && isspace(dms[end - 1]))
      --end;
    ind = NONE;
    int k = -1;
    if (end > beg && (k = lookup(hemispheres, dms[beg])) >= 0) {
      ind = (k / 2) ? LONGITUDE : LATITUDE;
      sign = k % 2 ? 1 : -1;
      ++beg;
    }
    if (end > beg && (k = lookup(hemispheres, dms[end-1])) >= 0) {
      if (k >= 0) {
	if (ind != NONE) {
	  if (toupper(dms[beg - 1]) == toupper(dms[end - 1]))
	    throw std::out_of_range("Repeated hemisphere indicators "
				    + str(dms[beg - 1]) + " in "
				    + dms.substr(beg - 1, end - beg + 1));
	  else
	    throw std::out_of_range("Contradictory hemisphere indicators "
				    + str(dms[beg - 1]) + " and "
				    + str(dms[end - 1]) + " in "
				    + dms.substr(beg - 1, end - beg + 1));
	}
	ind = (k / 2) ? LONGITUDE : LATITUDE;
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
      throw std::out_of_range("Empty or incomplete DMS string " + dms);
    double ipieces[] = {0, 0, 0};
    double fpieces[] = {0, 0, 0};
    unsigned npiece = 0;
    double icurrent = 0;
    double fcurrent = 0;
    unsigned ncurrent = 0, p = beg;
    bool pointseen = false;
    double mult = 0;
    while (p < end) {
      char x = dms[p++];
      if ((k = lookup(digits, x)) >= 0) {
	++ncurrent;
	if (pointseen) {
	  mult *= 0.1;
	  fcurrent += k * mult;
	} else
	  icurrent = 10 * icurrent + k;
      } else if (x == '.') {
	if (pointseen)
	  throw std::out_of_range("Multiple decimal points in "
				  + dms.substr(beg, end - beg));
	pointseen = true;
	mult = 1;
      } else if ((k = lookup(dmsindicators, x)) >= 0) {
	if (unsigned(k) == npiece - 1)
	  throw std::out_of_range("Repeated " + components[k]
				  + " component in "
				  + dms.substr(beg, end - beg));
	else if (unsigned(k) < npiece)
	  throw std::out_of_range(components[k] + " component follows "
				  + components[npiece - 1] + " component in "
				  + dms.substr(beg, end - beg));
	if (ncurrent == 0)
	  throw std::out_of_range("Missing numbers in " + components[k]
				  + " component of "
				  + dms.substr(beg, end - beg));
	ipieces[k] = icurrent;
	fpieces[k] = icurrent + fcurrent;
	if (p < end) {
	  npiece = k + 1;
	  icurrent = 0;
	  fcurrent = 0;
	  ncurrent = 0;
	  mult = 0;
	}
      } else if (lookup(signs, x) >= 0)
	throw std::out_of_range("Internal sign in DMS string "
				+ dms.substr(beg, end - beg));
      else
	throw std::out_of_range("Illegal character " + str(x)
				+ " in DMS string "
				+ dms.substr(beg, end - beg));
    }
    if (lookup(dmsindicators, dms[p - 1]) < 0) {
      if (npiece >= 3)
	throw std::out_of_range("Extra text following seconds in DMS string "
				+ dms.substr(beg, end - beg));
      if (ncurrent == 0)
	throw std::out_of_range("Missing numbers in " + components[k]
				+ " component of "
				+ dms.substr(beg, end - beg));
      ipieces[npiece] = icurrent;
      fpieces[npiece] = icurrent + fcurrent;
    }
    if (pointseen && mult == 0)
      throw std::out_of_range("Decimal point in non-terminal component of "
			      + dms.substr(beg, end - beg));
    // Note that we accept 59.999999... even though it rounds to 60.
    if (ipieces[1] >= 60)
      throw std::out_of_range("Minutes " + str(fpieces[1])
			      + " not in range [0, 60)");
    if (ipieces[2] >= 60)
      throw std::out_of_range("Seconds " + str(fpieces[2])
			      + " not in range [0, 60)");
    // Assume check on range of result is made by calling routine (which might
    // be able to offer a better diagnostic).
    return sign * (fpieces[0] + (fpieces[1] + fpieces[2] / 60) / 60);
  }

  std::string DMS::Encode(double angle, component trailing, unsigned prec,
			     flag ind) {
    // Assume check on range of input angle has been made by calling
    // routine (which might be able to offer a better diagnostic).
    //
    // 15 - 2 * trailing = ceiling(log10(2^53/90/60^trailing)).
    // This suffices to give full double precision for numbers in [-90,90]
    prec = (std::min)(15 - 2 * unsigned(trailing), prec);
    double scale = 1;
    for (unsigned i = 0; i < unsigned(trailing); ++i)
      scale *= 60;
    for (unsigned i = 0; i < prec; ++i)
      scale *= 10;
    int sign = angle < 0 ? -1 : 1;
    angle *= sign;
    // Break off integer part to preserve precision in manipulation of
    // fractional part.
    double
      idegree = floor(angle),
      fdegree = floor((angle - idegree) * scale + 0.5) / scale;
    if (fdegree >= 1) {
      idegree += 1;
      fdegree -= 1;
    }
    double pieces[3] = {fdegree, 0, 0};
    for (unsigned i = 1; i <= unsigned(trailing); ++i) {
      double
	ip = floor(pieces[i - 1]),
	fp = pieces[i - 1] - ip;
      pieces[i] = fp * 60;
      pieces[i - 1] = ip;
    }
    pieces[0] += idegree;
    std::ostringstream s;
    s << std::fixed  << std::setfill('0');
    if (ind == NONE && sign < 0)
      s << '-';
    switch (trailing) {
    case DEGREE:
      if (ind != NONE)
	s << std::setw(1 + ind + prec + (prec ? 1 : 0));
      s << std::setprecision(prec) << pieces[0];
      break;
    default:
      if (ind != NONE)
	s << std::setw(1 + ind);
      s << std::setprecision(0) << pieces[0] << char(tolower(dmsindicators[0]));
      switch (trailing) {
      case MINUTE:
	s << std::setw(2 + prec + (prec ? 1 : 0)) << std::setprecision(prec)
	  << pieces[1] <<  char(tolower(dmsindicators[1]));
	break;
      case SECOND:
	s << std::setw(2) << pieces[1] <<  char(tolower(dmsindicators[1]))
	  << std::setw(2 + prec + (prec ? 1 : 0)) << std::setprecision(prec)
	  << pieces[2] <<  char(tolower(dmsindicators[2]));
	break;
      default:
	break;
      }
    }
    if (ind != NONE)
      s << hemispheres[(ind == LATITUDE ? 0 : 2) + (sign < 0 ? 0 : 1)];
    return s.str();
  }

} // namespace GeographicLib

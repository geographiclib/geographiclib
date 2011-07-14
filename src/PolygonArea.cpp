/**
 * \file PolygonArea.cpp
 * \brief Implementation for GeographicLib::PolygonArea class
 *
 * Copyright (c) Charles Karney (2010, 2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/PolygonArea.hpp>

#define GEOGRAPHICLIB_POLYGONAREA_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_POLYGONAREA_CPP)
RCSID_DECL(GEOGRAPHICLIB_POLYGONAREA_HPP)

namespace GeographicLib {

  using namespace std;

  void PolygonArea::AddPoint(real lat, real lon) throw() {
    if (_num == 0) {
      _lat0 = _lat1 = lat;
      _lon0 = _lon1 = lon;
    } else {
      real s12, S12, t;
      _g.GenInverse(_lat1, _lon1, lat, lon,
                    Geodesic::DISTANCE | (_polyline ? 0 : Geodesic::AREA),
                    s12, t, t, t, t, t, S12);
      _perimetersum += s12;
      if (!_polyline) {
        _areasum += S12;
        _crossings += transit(_lon1, lon);
      }
      _lat1 = lat;
      _lon1 = lon;
    }
    ++_num;
  }

  unsigned PolygonArea::Compute(bool reverse, bool sign,
                                real& perimeter, real& area) const throw() {
    real s12, S12, t;
    if (_num < 2) {
      perimeter = 0;
      if (!_polyline)
        area = 0;
      return _num;
    }
    if (_polyline) {
      perimeter = _perimetersum();
      return _num;
    }
    _g.GenInverse(_lat1, _lon1, _lat0, _lon0,
                  Geodesic::DISTANCE | Geodesic::AREA,
                  s12, t, t, t, t, t, S12);
    perimeter = _perimetersum(s12);
    Accumulator<real> tempsum(_areasum);
    tempsum += S12;
    int crossings = _crossings + transit(_lon1, _lon0);
    if (crossings & 1)
      tempsum += (tempsum < 0 ? 1 : -1) * _area0/2;
    // area is with the clockwise sense.  If !reverse convert to
    // counter-clockwise convention.
    if (!reverse)
      tempsum *= -1;
    // If sign put area in (-area0/2, area0/2], else put area in [0, area0)
    if (sign) {
      if (tempsum > _area0/2)
        tempsum -= _area0;
      else if (tempsum <= -_area0/2)
        tempsum += _area0;
    } else {
      if (tempsum >= _area0)
        tempsum -= _area0;
      else if (tempsum < 0)
        tempsum += _area0;
    }
    area = tempsum();
    return _num;
  }

  unsigned PolygonArea::Compute(bool reverse, bool sign,
                                real& perimeter, real& area,
                                real lat, real lon) const throw() {
    PolygonArea p(*this);
    p.AddPoint(lat, lon);
    return p.Compute(reverse, sign, perimeter, area);
  }

} // namespace GeographicLib

/**
 * \file Geoid.hpp
 * \brief Header for GeographicLib::Geoid class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEOID_HPP)
#define GEOGRAPHICLIB_GEOID_HPP "$Id$"

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

namespace GeographicLib {

  class Geoid {
  private:
    static inline double sq(double x) throw() { return x * x; }
    std::string _filename;
    mutable std::ifstream _file;
    double _rlonres, _rlatres;
    std::string _description;
    double _offset, _scale, _maxerror, _rmserror;
    unsigned _width, _height;
    unsigned _datastart;
    mutable std::vector< std::vector<unsigned short> > _data;
    mutable bool _cache;
    mutable unsigned _ix, _iy;
    mutable double _v00, _v01, _v10, _v11;

    double rawval(unsigned ix, unsigned iy) const {
      if (ix >= _width)
	ix -= _width;
      if (_cache)
	return double(_data[iy][ix]);
      else {
	_file.seekg(_datastart + 2 * (iy * _width + ix), std::ios::beg);
	char a, b;
	_file.get(a);
	_file.get(b);
	return double((unsigned char)(a) * 256u + (unsigned char)(b));
      }
    }
  public:
    Geoid(const std::string& filename);
    void CacheAll() const;
    double operator()(float lat, float lon) const {
      if (lon < -180 || lon > 360)
	throw std::out_of_range("Longitude not in [-180, 360]");
      if (lat < -90 || lat > 90)
	throw std::out_of_range("Latitude not in [-90, 90]");
      if (lon < 0)
	lon += 360;
      double
	fy = (90 - lat) * _rlatres,
	fx = lon * _rlatres;
      unsigned
	iy = unsigned(fy),
	ix = unsigned(fx);
      if (iy == _height - 1)
	--iy;
      fx -= ix;
      fy -= iy;
      if (!(ix == _ix && iy == _iy)) {
	_v00 = rawval(ix    , iy    );
	_v01 = rawval(ix + 1, iy    );
	_v10 = rawval(ix    , iy + 1);
	_v11 = rawval(ix + 1, iy + 1);
	_ix = ix;
	_iy = iy;
      }
      double a = (1 - fx) * _v00 + fx * _v01;
      double b = (1 - fx) * _v10 + fx * _v11;
      double c = (1 - fy) * a + fy * b;
      return _offset + _scale * c;
    }
    const std::string& Description() const { return _description; }
    const std::string& GeoidFile() const { return _filename; }
    double MaxError() const { return _maxerror; }
    double RMSError() const { return _rmserror; }
  };
} //namespace GeographicLib
#endif

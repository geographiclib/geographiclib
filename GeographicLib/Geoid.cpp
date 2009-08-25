/**
 * \file Geoid.cpp
 * \brief Implementation for GeographicLib::Geoid class
 *
 * Copyright (c) Charles Karney (2009) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/geographic/
 **********************************************************************/

#include "GeographicLib/Constants.hpp"
#include "GeographicLib/Geoid.hpp"
#include <sstream>
#include <limits>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

#define GEOGRAPHICLIB_GEOID_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GEOID_CPP)
RCSID_DECL(GEOGRAPHICLIB_GEOID_HPP)

#if !defined(GEOID_DEFAULT_PATH)
/**
 * The order of the series relating \e lambda and \e eta when expanded in
 * powers of \e fp.
 **********************************************************************/
#if defined(_MSC_VER)
#define GEOID_DEFAULT_PATH "C:/cygwin/usr/local/share/geographiclib/geoids" 
#else
#define GEOID_DEFAULT_PATH "/usr/local/share/geographiclib/geoids" 
#endif
#endif

#if defined(_MSC_VER)
// Squelch warnings about unsafe use of getenv
#pragma warning (disable: 4996)
#endif

namespace GeographicLib {

  using namespace std;

  Geoid::Geoid(const std::string& geoid, const std::string& path)
    : _a( Constants::WGS84_a() )
    , _e2( (2 - 1/Constants::WGS84_r())/Constants::WGS84_r() )
    , _degree( Constants::degree() )
    , _eps( sqrt(numeric_limits<double>::epsilon()) ) {
    string dir = path;
    if (dir.size() == 0)
      dir = GeoidPath();
    if (dir.size() == 0)
      dir = DefaultPath();
    _filename = dir + "/" + geoid + ".pgm";
    _file.open(_filename.c_str(), ios::binary);
    if (!(_file.good()))
      throw out_of_range("File not readable " + _filename);
    string s;
    if (!(getline(_file, s) && s == "P5"))
      throw out_of_range("File not in PGM format " + _filename);
    _offset = numeric_limits<double>::max();
    _scale = 0;
    _maxerror = _rmserror = -1;
    _description = "NONE";
    while (getline(_file, s)) {
      if (s.size() == 0)
	continue;
      if (s[0] == '#') {
	if (s.substr(0, 14) == "# Description ")
	  _description = s.substr(14);
	else if (s.substr(0,9) == "# Offset ") {
	  s = s.substr(9);
	  istringstream is(s);
	  if (!(is >> _offset))
	    throw out_of_range("Error reading offset " + _filename);
	} else if (s.substr(0, 8) == "# Scale ") {
	  s = s.substr(8);
	  istringstream is(s);
	  if (!(is >> _scale))
	    throw out_of_range("Error reading scale " + _filename);
	} else if (s.substr(0,19) == "# MaxBilinearError ") {
	  s = s.substr(19);
	  istringstream is(s);
	  // It's not an error if the error can't be read
	  is >> _maxerror;
	} else if (s.substr(0,19) == "# RMSBilinearError ") {
	  s = s.substr(19);
	  istringstream is(s);
	  // It's not an error if the error can't be read
	  is >> _rmserror;
	}
      } else {
	istringstream is(s);
	if (!(is >> _width >> _height))
	  throw out_of_range("Error reading raster size " + _filename);
	break;
      }
    }
    {
      unsigned maxval;
      if (!(_file >> maxval))
	throw out_of_range("Error reading maxval " + _filename);
      if (maxval != 0xffffu)
	throw out_of_range("Maxval not equal to 2^16-1 " + _filename);
      // Add 1 for whitespace after maxval
      _datastart = unsigned(_file.tellg()) + 1u;
    }
    if (_offset == numeric_limits<double>::max())
      throw out_of_range("Offset not set " + _filename);
    if (_scale == 0)
      throw out_of_range("Scale not set " + _filename);
    if (_scale < 0)
      throw out_of_range("Scale must be positive " + _filename);
    if (_height < 2 || _width < 1)
      throw out_of_range("Raster size too small " + _filename);
    _file.seekg(0, ios::end);
    if (_datastart + 2 * _width * _height != _file.tellg())
      throw out_of_range("File has the wrong length " + _filename);
    _rlonres = _width / 360.0;
    _rlatres = (_height - 1) / 180.0;
    _cache = false;
    _ix = _width;
    _iy = _height;
    _file.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
  }

  double Geoid::height(double lat, double lon, bool gradp,
		       double& gradn, double& grade) const  {
    if (lon < 0)
      lon += 360;
    double
      fy = (90 - lat) * _rlatres,
      fx = lon * _rlonres;
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
    double
      a = (1 - fx) * _v00 + fx * _v01,
      b = (1 - fx) * _v10 + fx * _v11,
      c = (1 - fy) * a + fy * b,
      h = _offset + _scale * c;
    if (gradp) {
      double
	phi = lat * _degree,
	cosphi = cos(phi),
	sinphi = sin(phi),
	n = 1/sqrt(1 - _e2 * sinphi * sinphi);
      gradn = ((1 - fx) * (_v00 - _v10) + fx * (_v01 - _v11)) *
	_rlatres / _degree / (_a * (1 - _e2) * n * n * n);
      grade = (cosphi > _eps ?
	       ((1 - fy) * (_v01 - _v00) + fy * (_v11 - _v10)) /   cosphi :
	       (sinphi > 0 ? _v11 - _v10 : _v01 - _v00) * _rlatres / _degree ) *
	_rlonres / (_degree * _a * n);
      gradn *= _scale;
      grade *= _scale;
    }
    return h;
  }

  void Geoid::CacheClear() const {
    _cache = false;
    try {
      _data.clear();
      // Use swap to release memory back to system
      std::vector< std::vector<unsigned short> > t;
      _data.swap(t);
    }
    catch (std::bad_alloc&) {
    }
  }

  void Geoid::CacheArea(double south, double west,
			double north, double east) const {
    if (south > north) {
      CacheClear();
      return;
    }
    west += west < 0 ? 360 : west >= 360 ? -360 : 0;
    east += east < 0 ? 360 : east >= 360 ? -360 : 0;
    if (east <= west)
      east += 360;
    double
      fn = (90 - north) * _rlatres,
      fs = (90 - south) * _rlatres,
      fw = west * _rlonres,
      fe = east * _rlonres;
    // Bounding indices for cached area
    unsigned
      in = unsigned(fn),
      is = unsigned(fs) + 1,
      iw = unsigned(fw),
      ie = unsigned(fe) + 1;
    if (is == _height)
      --is;
    if (in == _height - 1)
      --in;
    if (ie - iw >= _width - 1) {
      // Include entire longitude range
      iw = 0;
      ie = _width -1;
    }
    unsigned oysize = unsigned(_data.size());
    _xsize = ie - iw + 1;
    _ysize = is - in + 1;
    _xoffset = iw;
    _yoffset = in;

    try {
      _data.resize(_ysize, vector<unsigned short>(_xsize));
      for (unsigned iy = min(oysize, _ysize); iy--;)
	_data[iy].resize(_xsize);
    }
    catch (std::bad_alloc&) {
      CacheClear();
      throw out_of_range("Insufficient memory for caching " + _filename);
    }
    
    try {
      unsigned
	ie1 = min(_width - 1, ie),
	w1 = ie1 - iw + 1;
      vector<char> buf(2 * w1);
      for (unsigned iy = in; iy <= is; ++iy) {
	_file.seekg(_datastart + 2 * (iy * _width + iw), ios::beg);
	_file.read(&(buf[0]), 2 * w1);
	for (unsigned ix = 0; ix < w1; ++ix)
	  _data[iy - in][ix] =
	    (unsigned short)((unsigned char)buf[2 * ix] * 256u +
			     (unsigned char)buf[2 * ix + 1]);
      }
      if (ie1 < ie) {
	// Cached area wraps past longitude = 0
	ie1 = ie - _width;
	unsigned
	  iw1 = 0,
	  w2 = ie1 + 1;
	buf.resize(2 * w1);
	for (unsigned iy = in; iy <= is; ++iy) {
	  _file.seekg(_datastart + 2 * (iy * _width + iw1), ios::beg);
	  _file.read(&(buf[0]), 2 * w2);
	  for (unsigned ix = 0; ix < w2; ++ix)
	  _data[iy - in][ix + w1] =
	    (unsigned short)((unsigned char)buf[2 * ix] * 256u +
			     (unsigned char)buf[2 * ix + 1]);
	}
      }
      _cache = true;
    }
    catch (std::exception& e) {
      CacheClear();
      throw out_of_range(string("Error filling cache ") + e.what());
    }
  }

  std::string Geoid::DefaultPath() {
    return string(GEOID_DEFAULT_PATH);
  }

  std::string Geoid::GeoidPath() {
    string path;
    char* geoidpath = getenv("GEOID_PATH");
    if (geoidpath)
      path = string(geoidpath);
    return path;
  }
}

/*

h=(1-fy)*((1-fx)*v00+fx*v01)+
fy*((1-fx)*v10+fx*v11);

dh/dfx
(1-fy)*(v01-v00)+fy*(v11-v10);

dh/dfy
(1-fx)*(v10-v00)+fx*(v11-v01)

	fy = (90 - lat) * _rlatres,
	fx = lon * _rlonres;

dfy/dlat = -rlatres
dfx/dlon = rlonres

degree = pi/180
phi = lat * degree
lam = lon * degree

e2 = f*(2-f)

      n = 1/sqrt(1 - e2 * sin(phi)^2)

dy/dphi = a * (1-e2) * n*n*n
dx/dlon = a * cos(phi) * n

dh/dx = dh/dfx * dfx/dlon * dlon/dlam * dlam/dx

= ( (1-fy)*(v01-v00)+fy*(v11-v10) ) * rlonres / degree / (a * cos(phi) * n)


Near lat = 90, v01=v00, cos(phi) = (pi/2-phi) = fy/(rlatres/degree)
dh/dx = fy*(v11-v10) * rlonres / degree / (a * fy/(rlatres/degree) * n)
      = (v11-v10) * rlonres*rlatres / degree^2 / (a * n)

Similarly near lat = -90
dh/dx = (v01-v00) * rlonres*rlatres / degree^2 / (a * n)

dh/dy = dh/dfy * dfy/dlat * dlat/dphi * dphi/dy
= ( (1-fx)*(v10-v00)+fx*(v11-v01) ) * -rlatres / degree / (a * (1-e2) * n^3)
= ( (1-fx)*(v00-v10)+fx*(v01-v11) ) * rlatres / degree / (a * (1-e2) * n^3)

 */

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

  Geoid::Geoid(const std::string& name, bool cubic, const std::string& path)
    : _cubic(cubic)
    , _a( Constants::WGS84_a() )
    , _e2( (2 - 1/Constants::WGS84_r())/Constants::WGS84_r() )
    , _degree( Constants::degree() )
    , _eps( sqrt(numeric_limits<double>::epsilon()) ) {
    string dir = path;
    if (dir.size() == 0)
      dir = GeoidPath();
    if (dir.size() == 0)
      dir = DefaultPath();
    _filename = dir + "/" + name + ".pgm";
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
    if (_width & 1)
      throw out_of_range("Raster width is odd " + _filename);
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
    int
      iy = int(fy),
      ix = int(fx);
    if (iy == _height - 1)
      --iy;
    fx -= ix;
    fy -= iy;
    if (!(ix == _ix && iy == _iy)) {
      _ix = ix;
      _iy = iy;
      if (!_cubic) {
      _v00 = rawval(ix    , iy    );
      _v01 = rawval(ix + 1, iy    );
      _v10 = rawval(ix    , iy + 1);
      _v11 = rawval(ix + 1, iy + 1);
      } else {
	double v[16];
	for (int j = -1, k = 0; j < 3; ++j)
	  for (int i = -1; i < 3; ++i)
	    v[k++] = rawval(ix + i, iy + j);
	/*
	// This is the matrix for 124 weights
	double c[160] = {
	  -476, -709, -709, 1023, 1836, 1023, -242, -594, -594, -242,
	  1142, -381, -1778, -1749, 234, 1650, 726, 594, -396, -484,
	  -262, 1701, -950, 429, -1422, 1254, -726, 594, 396, -484,
	  -404, -611, 533, 297, -648, 429, 242, -594, 594, -242,
	  1142, -1778, -381, 1650, 234, -1749, -484, -396, 594, 726,
	  5176, -1662, -1662, -3102, -504, -3102, 1452, 396, 396, 1452,
	  1864, 3510, -1770, 1254, -288, -2706, -1452, 396, -396, 1452,
	  530, -70, -543, 198, 558, -1155, 484, -396, -594, 726,
	  -262, -950, 1701, 1254, -1422, 429, -484, 396, 594, -726,
	  1864, -1770, 3510, -2706, -288, 1254, 1452, -396, 396, -1452,
	  -1160, 2826, 2826, 1650, 1080, 1650, -1452, -396, -396, -1452,
	  -442, -106, 675, -198, 630, 1023, 484, 396, -594, -726,
	  -404, 533, -611, 429, -648, 297, -242, 594, -594, 242,
	  530, -543, -70, -1155, 558, 198, 726, -594, -396, 484,
	  -442, 675, -106, 1023, 630, -198, -726, -594, 396, 484,
	  316, -665, -665, -297, -540, -297, 242, 594, 594, 242,
	};
	_t0 = 8712;
	*/

	// This is the matrix for 012 weights;
        double c[160] = {
            0,   0,   0,    0,   0,    0,   0,   0,   0,   0, 
            9, -18, -88,    0,  96,   90,   0,   0, -60, -20, 
           -9,  18,   8,    0, -96,   30,   0,   0,  60, -20, 
            0,   0,   0,    0,   0,    0,   0,   0,   0,   0, 
            9, -88, -18,   90,  96,    0, -20, -60,   0,   0, 
          186, -42, -42, -150, -96, -150,  60,  60,  60,  60, 
           54, 162, -78,   30, -24,  -90, -60,  60, -60,  60, 
           -9, -32,  18,   30,  24,    0,  20, -60,   0,   0, 
           -9,   8,  18,   30, -96,    0, -20,  60,   0,   0, 
           54, -78, 162,  -90, -24,   30,  60, -60,  60, -60, 
          -54,  78,  78,   90, 144,   90, -60, -60, -60, -60, 
            9,  -8, -18,  -30, -24,    0,  20,  60,   0,   0, 
            0,   0,   0,    0,   0,    0,   0,   0,   0,   0, 
           -9,  18, -32,    0,  24,   30,   0,   0, -60,  20, 
            9, -18,  -8,    0, -24,  -30,   0,   0,  60,  20, 
            0,   0,   0,    0,   0,    0,   0,   0,   0,   0,
        };
        _t0 = 240;

	/*
	// This is the matrix for 011 weights
	double c[160] = {
	  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	  18, -36, -169, 0, 186, 171, 0, 0, -114, -38, 
	  -18, 36, 17, 0, -186, 57, 0, 0, 114, -38, 
	  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	  18, -169, -36, 171, 186, 0, -38, -114, 0, 0, 
	  348, -69, -69, -285, -204, -285, 114, 114, 114, 114, 
	  108, 297, -159, 57, -24, -171, -114, 114, -114, 114, 
	  -18, -59, 36, 57, 42, 0, 38, -114, 0, 0, 
	  -18, 17, 36, 57, -186, 0, -38, 114, 0, 0, 
	  108, -159, 297, -171, -24, 57, 114, -114, 114, -114, 
	  -108, 159, 159, 171, 252, 171, -114, -114, -114, -114, 
	  18, -17, -36, -57, -42, 0, 38, 114, 0, 0, 
	  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	  -18, 36, -59, 0, 42, 57, 0, 0, -114, 38, 
	  18, -36, -17, 0, -42, -57, 0, 0, 114, 38, 
	  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	};
	_t0 = 456;
	*/
	/*
	// This is the matrix for 111 weights
	double c[160] = {
	  -138, -109, -109, 195, 288, 195, -50, -90, -90, -50,
	  264, -93, -223, -345, -24, 165, 150, 90, -30, -50,
	  -24, 333, -157, 105, -156, 135, -150, 90, 30, -50,
	  -102, -131, 89, 45, -108, 105, 50, -90, 90, -50,
	  264, -223, -93, 165, -24, -345, -50, -30, 90, 150,
	  558, -171, -171, -315, -48, -315, 150, 30, 30, 150,
	  222, 351, -189, 135, -12, -285, -150, 30, -30, 150,
	  156, 43, -147, 15, 84, -255, 50, -30, -90, 150,
	  -24, -157, 333, 135, -156, 105, -50, 30, 90, -150,
	  222, -189, 351, -285, -12, 135, 150, -30, 30, -150,
	  -102, 309, 309, 165, 72, 165, -150, -30, -30, -150,
	  -96, 37, 207, -15, 96, 195, 50, 30, -90, -150,
	  -102, 89, -131, 105, -108, 45, -50, 90, -90, 50,
	  156, -147, 43, -255, 84, 15, 150, -90, -30, 50,
	  -96, 207, 37, 195, 96, -15, -150, -90, 30, 50,
	  42, -149, -149, -45, -72, -45, 50, 90, 90, 50,
	};
	_t0 = 1200;

	/*
	// This is the matrix for 123 weights
	double c[160] = {
	  -2193, -2344, -2344,  3870,  6468,  3870,  -980, -2100, -2100,  -980,
	   4818, -1926, -6266, -6810,   672,  5370,  2940,  2100, -1260, -1540,
	   -978,  6726, -3494,  2010, -4872,  4110, -2940,  2100,  1260, -1540,
	  -1647, -2456,  2024,   930, -2268,  1770,   980, -2100,  2100,  -980,
	   4818, -6266, -1926,  5370,   672, -6810, -1540, -1260,  2100,  2940,
	  16749, -5088, -5088, -9990, -1764, -9990,  4620,  1260,  1260,  4620,
	   6291, 11208, -5592,  3870,  -756, -8730, -4620,  1260, -1260,  4620,
	   2382,   146, -2514,   750,  1848, -4710,  1540, -1260, -2100,  2940,
	   -978, -3494,  6726,  4110, -4872,  2010, -1540,  1260,  2100, -2940,
	   6291, -5592, 11208, -8730,  -756,  3870,  4620, -1260,  1260, -4620,
	  -3411,  9192,  9192,  5130,  3276,  5130, -4620, -1260, -1260, -4620,
	  -1902,  -106,  3114,  -510,  2352,  4110,  1540,  1260, -2100, -2940,
	  -1647,  2024, -2456,  1770, -2268,   930,  -980,  2100, -2100,   980,
	   2382, -2514,   146, -4710,  1848,   750,  2940, -2100, -1260,  1540,
	  -1902,  3114,  -106,  4110,  2352,  -510, -2940, -2100,  1260,  1540,
 	   1167, -2624, -2624, -1170, -1932, -1170,   980,  2100,  2100,   980,
	};
	_t0 = 30240;
	*/
	for (unsigned i = 0; i < 10; ++i) {
	  _t[i] = 0;
	  for (unsigned j = 0; j < 16; ++j)
	    _t[i] += v[j] * c[10 * j + i];
	}
      }
    }
    if (!_cubic) {  
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
    } else {
      double h = _t[0] + fx * (_t[1] + fx * (_t[3] + fx * _t[6])) +
	fy * (_t[2] + fx * (_t[4] + fx * _t[7]) +
	     fy * (_t[5] + fx * _t[8] + fy * _t[9]));
      h = _offset + _scale * h / _t0;
      return h;
    }
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
    int
      in = int(fn),
      is = int(fs) + 1,
      iw = int(fw),
      ie = int(fe) + 1;
    if (is == _height)
      --is;
    if (in == _height - 1)
      --in;
    if (ie - iw >= _width - 1) {
      // Include entire longitude range
      iw = 0;
      ie = _width -1;
    }
    int oysize = int(_data.size());
    _xsize = ie - iw + 1;
    _ysize = is - in + 1;
    _xoffset = iw;
    _yoffset = in;

    try {
      _data.resize(_ysize, vector<unsigned short>(_xsize));
      for (int iy = min(oysize, _ysize); iy--;)
	_data[iy].resize(_xsize);
    }
    catch (std::bad_alloc&) {
      CacheClear();
      throw out_of_range("Insufficient memory for caching " + _filename);
    }
    
    try {
      int
	ie1 = min(_width - 1, ie),
	w1 = ie1 - iw + 1;
      vector<char> buf(2 * w1);
      for (int iy = in; iy <= is; ++iy) {
	_file.seekg(_datastart + 2 * (iy * _width + iw), ios::beg);
	_file.read(&(buf[0]), 2 * w1);
	for (int ix = 0; ix < w1; ++ix)
	  _data[iy - in][ix] =
	    (unsigned short)((unsigned char)buf[2 * ix] * 256u +
			     (unsigned char)buf[2 * ix + 1]);
      }
      if (ie1 < ie) {
	// Cached area wraps past longitude = 0
	ie1 = ie - _width;
	int
	  iw1 = 0,
	  w2 = ie1 + 1;
	buf.resize(2 * w1);
	for (int iy = in; iy <= is; ++iy) {
	  _file.seekg(_datastart + 2 * (iy * _width + iw1), ios::beg);
	  _file.read(&(buf[0]), 2 * w2);
	  for (int ix = 0; ix < w2; ++ix)
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

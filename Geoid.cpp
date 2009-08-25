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

  Geoid::Geoid(const std::string& geoid, const std::string& path) {
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

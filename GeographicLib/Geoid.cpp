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

#define GEOGRAPHICLIB_GEOID_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GEOID_CPP)
RCSID_DECL(GEOGRAPHICLIB_GEOID_HPP)

namespace GeographicLib {

  using namespace std;

  Geoid::Geoid(const std::string& filename) {
    _filename = filename;
    _file.open(filename.c_str(), ios::binary);
    if (!(_file.good()))
      throw out_of_range("File not readable " + filename);
    string s;
    if (!(getline(_file, s) && s == "P5"))
      throw out_of_range("File not in PGM format " + filename);
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
	    throw out_of_range("Error reading offset " + filename);
	} else if (s.substr(0, 8) == "# Scale ") {
	  s = s.substr(8);
	  istringstream is(s);
	  if (!(is >> _scale))
	    throw out_of_range("Error reading scale " + filename);
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
	  throw out_of_range("Error reading raster size " + filename);
	break;
      }
    }
    {
      unsigned maxval;
      if (!(_file >> maxval))
	throw out_of_range("Error reading maxval " + filename);
      if (maxval != 0xffffu)
	throw out_of_range("Maxval not equal to 2^16-1 " + filename);
      // Add 1 for whitespace after maxval
      _datastart = unsigned(_file.tellg()) + 1u;
    }
    if (_offset == numeric_limits<double>::max())
      throw out_of_range("Offset not set " + filename);
    if (_scale == 0)
      throw out_of_range("Scale not set " + filename);
    if (_scale < 0)
      throw out_of_range("Scale must be positive " + filename);
    if (_height < 2 || _width < 1)
      throw out_of_range("Raster size too small " + filename);
    _file.seekg(0, ios::end);
    if (_datastart + 2 * _width * _height != _file.tellg())
      throw out_of_range("File has the wrong length " + filename);
    _rlonres = _width / 360.0;
    _rlatres = (_height - 1) / 180.0;
    _cache = false;
    _ix = _width;
    _iy = _height;
    _file.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
  }

  void Geoid::CacheAll() const {
    try {
      _data.resize(_height, std::vector<unsigned short>(_width));
    }
    catch (std::bad_alloc&) {
      _cache = false;
      _data.clear();
      throw out_of_range("Insufficient memory for caching " + _filename);
    }
    try {
      _file.seekg(_datastart, ios::beg);
      vector<char> buf(2 * _width);
      for (unsigned iy = 0; iy < _height; ++iy) {
	_file.read(&(buf[0]), 2 * _width);
	for (unsigned ix = 0; ix < _width; ++ix)
	  _data[iy][ix] = (unsigned short)((unsigned char)buf[2 * ix] * 256u +
					   (unsigned char)buf[2 * ix + 1]);
      }
      _cache = true;
    }
    catch (std::exception& e) {
      _cache = false;
      _data.clear();
      throw out_of_range(string("Error filling cache ") + e.what());
    }
  }

  
}

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

namespace GeographicLib {

  /**
   * \brief Computing the height of the geoid
   *
   * This class evaluated the height of one of the standard geoids, EGM84,
   * EGM96, or EGM2008 by bilinear interpolation into a rectangular grid of
   * data.  These geoid models are documented in
   * - EGM84:
   *   http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html
   * - EGM96:
   *   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
   * - EGM2008:
   *   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008
   *
   * The geoids are defined in terms of spherical harmonics.  However in order
   * to provide a quick and flexible method of evaluating the geoid heights,
   * this class reads in a data file giving the geoid heights on a rectangle
   * grid and determines the geoid height by bilinear interpolation.
   *
   * See \ref geoid for details of how to install the data sets, the data
   * format, estimates of the interpolation errors, and how to use caching.
   *
   * In addition to returning the geoid height, the gradient of the geoid can
   * be calculated.  The gradient is defined as the rate of change of the geoid
   * as a function of position on the ellipsoid.  This uses the parameters for
   * the WGS84 ellipsoid.  The gradient defined in terms of the interpolated
   * heights and, because the interpolation is bilinear, the gradient has
   * discontinuities on cell boundaries.
   *
   * This class is \e not thread safe in that a single instantiation cannot be
   * safely used by multiple threads.  If multiple threads need to calculate
   * geoid heights they should all construct thread-local instantiations.
   **********************************************************************/

  class Geoid {
  private:
    std::string _filename;
    const double _a, _e2, _degree, _eps;
    mutable std::ifstream _file;
    double _rlonres, _rlatres;
    std::string _description;
    double _offset, _scale, _maxerror, _rmserror;
    unsigned _width, _height;
    unsigned _datastart;
    // Area cache
    mutable std::vector< std::vector<unsigned short> > _data;
    mutable bool _cache;
    // NE corner and extent of cache
    mutable unsigned _xoffset, _yoffset, _xsize, _ysize;
    // Cell cache
    mutable unsigned _ix, _iy;
    mutable double _v00, _v01, _v10, _v11;

    double rawval(unsigned ix, unsigned iy) const {
      if (ix >= _width)
	ix -= _width;
      if (_cache && iy >= _yoffset && iy < _yoffset + _ysize &&
	  ((ix >= _xoffset && ix < _xoffset + _xsize) ||
	   (ix + _width >= _xoffset && ix + _width < _xoffset + _xsize))) {
	return double(_data
		      [iy - _yoffset]
		      [ix >= _xoffset ?
		       ix - _xoffset :
		       ix + _width - _xoffset]);
      } else {
	_file.seekg(_datastart + 2 * (iy * _width + ix), std::ios::beg);
	char a, b;
	_file.get(a);
	_file.get(b);
	return double((unsigned char)(a) * 256u + (unsigned char)(b));
      }
    }
    double height(double lat, double lon, bool gradp,
	       double& grade, double& gradn) const;
  public:

    /**
     * Create a Geoid loading the data for geoid \e name.  The data file is
     * formed by appending ".pgm" to the name.  If \e path is specified, then
     * the file is load from that directory.  Otherwise the path is given by
     * the GEOID_PATH environment variable.  If that is undefined, a
     * compile-time default path is used (/usr/local/share/geographiclib/geoids
     * on non-Windows systems and
     * C:/cygwin/usr/local/share/geographiclib/geoids on Windows systems.
     **********************************************************************/
    Geoid(const std::string& name, const std::string& path = "");

    /**
     * Cache the data for the rectangular area defined by the four arguments \e
     * south, \e west, \e north, \e east (all in degrees).  \e east is always
     * interpreted as being east of \e west, if necessary by adding
     * 360<sup>o</sup> to its value.
     **********************************************************************/
    void CacheArea(double south, double west, double north, double east) const;

    /**
     * Cache all the data.  On most computers, this is fast for data sets with
     * grid resolution of 5' or coarser.  For a 1' grid, the required RAM is
     * 450MB; a 2.5' grid needs 72MB; and a 5' grid needs 18MB.
     **********************************************************************/
    void CacheAll() const { CacheArea(-90.0, 0.0, 90.0, 360.0); }

    /**
     * Clear the cache.
     **********************************************************************/
    void CacheClear() const;

    /**
     * Return the geoid height in meters for latitude \e lat (in [-90, 90]) and
     * longitude \e lon (in [-180,360]), both in degrees.
     **********************************************************************/
    double operator()(double lat, double lon) const {
      double gradn, grade;
      return height(lat, lon, false, gradn, grade);
    }
    /**
     * Return the geoid height in meters for latitude \e lat (in [-90, 90]) and
     * longitude \e lon (in [-180,360]), both in degrees.  In addition compute
     * the gradient of the geoid height in the northerly \e gradn and easterly
     * \e grade directions.
     **********************************************************************/
    double operator()(double lat, double lon, double& gradn, double& grade)
      const {
      return height(lat, lon, true, gradn, grade);
    }

    /**
     * Return the geoid description if availabe in the data file.  If absent,
     * return "UNKNOWN".
     **********************************************************************/
    const std::string& Description() const { return _description; }

    /**
     * Return the path name used to load the geoid data.
     **********************************************************************/
    const std::string& GeoidFile() const { return _filename; }

    /**
     * Return a estimate of the maximum interpolation and quantization error
     * (meters).  This relies on the value being stored in the data file.  If
     * the value is absent, return -1.
     **********************************************************************/
    double MaxError() const { return _maxerror; }

    /**
     * Return a estimate of the RMS interpolation and quantization error
     * (meters).  This relies on the value being stored in the data file.  If
     * the value is absent, return -1.
     **********************************************************************/
    double RMSError() const { return _rmserror; }

    /**
     * Return the compile-time default path for the geoid data files.
     **********************************************************************/
    static std::string DefaultPath();

    /**
     * Return the value of the environment variable GEOID_PATH.
     **********************************************************************/
    static std::string GeoidPath();
  };

} //namespace GeographicLib
#endif

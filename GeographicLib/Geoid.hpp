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

#include "GeographicLib/Constants.hpp"
#include <vector>
#include <string>
#include <fstream>

namespace GeographicLib {

  /**
   * \brief Computing the height of the geoid
   *
   * This class evaluated the height of one of the standard geoids, EGM84,
   * EGM96, or EGM2008 by bilinear or cubic interpolation into a rectangular
   * grid of data.  These geoid models are documented in
   * - EGM84:
   *   http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html
   * - EGM96:
   *   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
   * - EGM2008:
   *   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008
   *
   * The geoids are defined in terms of spherical harmonics.  However in order
   * to provide a quick and flexible method of evaluating the geoid heights,
   * this class evaluates the height by interpolation inot a grid of
   * precomputed values.
   *
   * See \ref geoid for details of how to install the data sets, the data
   * format, estimates of the interpolation errors, and how to use caching.
   *
   * In addition to returning the geoid height, the gradient of the geoid can
   * be calculated.  The gradient is defined as the rate of change of the geoid
   * as a function of position on the ellipsoid.  This uses the parameters for
   * the WGS84 ellipsoid.  The gradient defined in terms of the interpolated
   * heights.
   *
   * This class is \e not thread safe in that a single instantiation cannot be
   * safely used by multiple threads.  If multiple threads need to calculate
   * geoid heights they should all construct thread-local instantiations.
   **********************************************************************/

  class Geoid {
  private:
    typedef Math::real real;
    static const unsigned stencilsize = 12;
    static const unsigned nterms = ((3 + 1) * (3 + 2))/2; // for a cubic fit
    static const real c0, c0n, c0s;
    static const real c3[stencilsize * nterms];
    static const real c3n[stencilsize * nterms];
    static const real c3s[stencilsize * nterms];

    std::string _name, _dir, _filename;
    const bool _cubic;
    const real _a, _e2, _degree, _eps;
    mutable std::ifstream _file;
    real _rlonres, _rlatres;
    std::string _description, _datetime;
    real _offset, _scale, _maxerror, _rmserror;
    int _width, _height;
    unsigned long long _datastart, _swidth;
    // Area cache
    mutable std::vector< std::vector<unsigned short> > _data;
    mutable bool _cache;
    // NE corner and extent of cache
    mutable int _xoffset, _yoffset, _xsize, _ysize;
    // Cell cache
    mutable int _ix, _iy;
    mutable real _v00, _v01, _v10, _v11;
    mutable real _t[nterms];
    void filepos(int ix, int iy) const {
      _file.seekg(std::ios::streamoff(_datastart +
                                      2ULL * (unsigned(iy) * _swidth +
                                              unsigned(ix))));
    }
    real rawval(int ix, int iy) const {
      if (iy < 0) {
        iy = -iy;
        ix += _width/2;
      } else if (iy >= _height) {
        iy = 2 * (_height - 1) - iy;
        ix += _width/2;
      }
      if (ix < 0)
        ix += _width;
      else if (ix >= _width)
        ix -= _width;
      if (_cache && iy >= _yoffset && iy < _yoffset + _ysize &&
          ((ix >= _xoffset && ix < _xoffset + _xsize) ||
           (ix + _width >= _xoffset && ix + _width < _xoffset + _xsize))) {
        return real(_data[iy - _yoffset]
                    [ix >= _xoffset ? ix - _xoffset : ix + _width - _xoffset]);
      } else {
        filepos(ix, iy);
        char a, b;
        _file.get(a);
        _file.get(b);
        return real((unsigned char)(a) * 256u + (unsigned char)(b));
      }
    }
    real height(real lat, real lon, bool gradp,
                real& grade, real& gradn) const;
    Geoid(const Geoid&);        // copy constructor not allowed
    Geoid& operator=(const Geoid&); // copy assignment not allowed
  public:

    /**
     * Create a Geoid loading the data for geoid \e name.  The data file is
     * formed by appending ".pgm" to the name.  If \e path is specified (and is
     * non-empty), then the file is loaded from directory, \e path.  Otherwise
     * the path is given by the GEOID_PATH environment variable.  If that is
     * undefined, a compile-time default path is used
     * (/usr/local/share/GeographicLib/geoids on non-Windows systems and
     * C:/cygwin/usr/local/share/GeographicLib/geoids on Windows systems).  The
     * final \e cubic argument specifies whether to use bilinear (\e cubic =
     * false) or cubic (\e cubic = true, the default) interpolation.  This may
     * throw an error because the file does not exist, is unreadable, or is
     * corrupt.
     **********************************************************************/
    explicit Geoid(const std::string& name, const std::string& path = "",
                   bool cubic = true);

    /**
     * Cache the data for the rectangular area defined by the four arguments \e
     * south, \e west, \e north, \e east (all in degrees).  \e east is always
     * interpreted as being east of \e west, if necessary by adding
     * 360<sup>o</sup> to its value.  This may throw an error because of
     * insufficent memory or because of an error reading the data from the
     * file.  In this case, you can catch the error and either do nothing (you
     * will have no cache in this case) or try again with a smaller area.
     **********************************************************************/
    void CacheArea(real south, real west, real north, real east) const;

    /**
     * Cache all the data.  On most computers, this is fast for data sets with
     * grid resolution of 5' or coarser.  For a 1' grid, the required RAM is
     * 450MB; a 2.5' grid needs 72MB; and a 5' grid needs 18MB.  This may throw
     * an error because of insufficent memory or because of an error reading
     * the data from the file.  In this case, you can catch the error and
     * either do nothing (you will have no cache in this case) or try using
     * Geoid::CacheArea on a specific area.
     **********************************************************************/
    void CacheAll() const { CacheArea(real(-90), real(0),
                                      real(90), real(360)); }

    /**
     * Clear the cache.  This never throws an error.
     **********************************************************************/
    void CacheClear() const throw();

    /**
     * Return the geoid height in meters for latitude \e lat (in [-90, 90]) and
     * longitude \e lon (in [-180,360]), both in degrees.  This may throw an
     * error because of an error reading data from disk.  However, it will not
     * throw if (\e lat, \e lon) is within a successfully cached area.
     **********************************************************************/
    Math::real operator()(real lat, real lon) const {
      real gradn, grade;
      return height(lat, lon, false, gradn, grade);
    }
    /**
     * Return the geoid height in meters for latitude \e lat (in [-90, 90]) and
     * longitude \e lon (in [-180,360]), both in degrees.  In addition compute
     * the gradient of the geoid height in the northerly \e gradn and easterly
     * \e grade directions.  This may throw an error because of an error
     * reading data from disk.  However, it will not throw if (\e lat, \e lon)
     * is within a successfully cached area.
     **********************************************************************/
    Math::real operator()(real lat, real lon, real& gradn, real& grade) const {
      return height(lat, lon, true, gradn, grade);
    }

    /**
     * Return the geoid description if available in the data file.  If absent,
     * return "NONE".
     **********************************************************************/
    const std::string& Description() const throw() { return _description; }

    /**
     * Return the date of the data file.  If absent, return "UNKNOWN".
     **********************************************************************/
    const std::string& DateTime() const throw() { return _datetime; }

    /**
     * Return the full file name used to load the geoid data.
     **********************************************************************/
    const std::string& GeoidFile() const throw() { return _filename; }

    /**
     * Return the "name" used to load the geoid data (from the first argument
     * of the constructor).
     **********************************************************************/
    const std::string& GeoidName() const throw() { return _name; }

    /**
     * Return the directory used to load the geoid data.
     **********************************************************************/
    const std::string& GeoidDirectory() const throw() { return _dir; }

    /**
     * Return the interpolation method (cubic or bilinear).
     **********************************************************************/
    const std::string Interpolation() const
    { return std::string(_cubic ? "cubic" : "bilinear"); }

    /**
     * Return a estimate of the maximum interpolation and quantization error
     * (meters).  This relies on the value being stored in the data file.  If
     * the value is absent, return -1.
     **********************************************************************/
    Math::real MaxError() const throw() { return _maxerror; }

    /**
     * Return a estimate of the RMS interpolation and quantization error
     * (meters).  This relies on the value being stored in the data file.  If
     * the value is absent, return -1.
     **********************************************************************/
    Math::real RMSError() const throw() { return _rmserror; }

    /**
     * Return offset (meters) for converting pixel values to geoid heights.
     **********************************************************************/
    Math::real Offset() const throw() { return _offset; }

    /**
     * Return scale (meters) for converting pixel values to geoid
     * heights.
     **********************************************************************/
    Math::real Scale() const throw() { return _scale; }

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

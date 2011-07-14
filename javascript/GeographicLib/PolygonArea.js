/**
 * PolygonArea.js
 * Transcription of PolygonArea.[ch]pp into javascript.
 *
 * See the documentation for the C++ class.  The conversion is mostly a
 * literal conversion from C++.  However there are two javascript-ready
 * interface routines.
 *
 *   GeographicLib.PolygonArea.Area(GeographicLib.Geodesic.WGS84,
 *                                  points, polyline);
 *
 * computes the area of a polygon with vertices given by an array
 * points, each of whose elements contains lat and lon fields.  The
 * function returns an object with fields.
 *
 *   number, perimeter, area
 *
 * There is no need to "close" the polygon.  If polyline is true, that
 * the points denote a polyline and its length is returned as the
 * perimeter (and the area is not calculated).
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * $Id$
 **********************************************************************/

// Load AFTER GeographicLib/Math.js and GeographicLib/Geodesic.js
GeographicLib.PolygonArea = {};

(function() {
  var m = GeographicLib.Math;
  var a = GeographicLib.Accumulator;
  var g = GeographicLib.Geodesic;
  var p = GeographicLib.PolygonArea;

  p.transit = function(lon1, lon2) {
    // Return 1 or -1 if crossing prime meridian in east or west direction.
    // Otherwise return zero.
    lon1 = g.AngNormalize(lon1);
    lon2 = g.AngNormalize(lon2);
    // treat lon12 = -180 as an eastward geodesic, so convert to 180.
    var lon12 = -g.AngNormalize(lon1 - lon2); // In (-180, 180]
    var cross =
      lon1 < 0 && lon2 >= 0 && lon12 > 0 ? 1 :
      lon2 < 0 && lon1 >= 0 && lon12 < 0 ? -1 : 0;
    return cross;
  }

  p.PolygonArea = function(earth, polyline) {
    this._earth = earth;
    this._area0 = 4 * Math.PI * earth._c2;
    this._polyline = polyline;
    if (!this._polyline)
      this._areasum = new a.Accumulator(0);
    this._perimetersum = new a.Accumulator(0);
    this.Clear();
  }

  p.PolygonArea.prototype.Clear = function() {
    this._num = 0;
    this._crossings = 0;
    if (!this._polyline)
      this._areasum.Set(0);
    this._perimetersum.Set(0);
    this._lat0 = this._lon0 = this._lat1 = this._lon1 = 0;
  }

  p.PolygonArea.prototype.AddPoint = function(lat, lon) {
    if (this._num == 0) {
      this._lat0 = this._lat1 = lat;
      this._lon0 = this._lon1 = lon;
    } else {
      var t = this._earth.Inverse(this._lat1, this._lon1, lat, lon,
				  g.DISTANCE | (this._polyline ? 0 : g.AREA));
      this._perimetersum.Add(t.s12);
      if (!this._polyline) {
        this._areasum.Add(t.S12);
        this._crossings += p.transit(this._lon1, lon);
      }
      this._lat1 = lat;
      this._lon1 = lon;
    }
    ++this._num;
  }

  // args = perimeter, area
  p.PolygonArea.prototype.Compute = function(reverse, sign, args) {
    if (this._num < 2) {
      args.perimeter = 0;
      if (!this._polyline)
        args.area = 0;
      return this._num;
    }
    if (this._polyline) {
      args.perimeter = this._perimetersum.Sum();
      return this._num;
    }
    var t = this._earth.Inverse(this._lat1, this._lon1, this._lat0, this._lon0,
				g.DISTANCE | g.AREA);
    args.perimeter = this._perimetersum.Sum(t.s12);
    var tempsum = new a.Accumulator(this._areasum);
    tempsum.Add(t.S12);
    var crossings = this._crossings + p.transit(this._lon1, this._lon0);
    if (crossings & 1)
      tempsum.Add( (tempsum.Sum() < 0 ? 1 : -1) * this._area0/2 );
    // area is with the clockwise sense.  If !reverse convert to
    // counter-clockwise convention.
    if (!reverse)
      tempsum.Negate();
    // If sign put area in (-area0/2, area0/2], else put area in [0, area0)
    if (sign) {
      if (tempsum.Sum() > this._area0/2)
        tempsum.Add( -this._area0);
      else if (tempsum.Sum() <= -this._area0/2)
        tempsum.Add( this._area0);
    } else {
      if (tempsum.Sum() >= this._area0)
        tempsum.Add( -this._area0);
      else if (tempsum < 0)
        tempsum.Add( this._area0);
    }
    args.area = tempsum.Sum();
    return this._num;
  }

  p.Area = function(earth, points, polyline) {
    var poly = new p.PolygonArea(earth, polyline);
    for (var i = 0; i < points.length; ++i)
      poly.AddPoint(points[i].lat, points[i].lon);
    var result = {};
    result.number = poly.Compute(false, true, result);
    return result;
  }

})();

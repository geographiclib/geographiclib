/**
 * PolygonArea.js
 * Transcription of PolygonArea.[ch]pp into javascript.
 *
 * See the documentation for the C++ class.  The conversion is a literal
 * conversion from C++.
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
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
      (lon2 < 0 && lon1 >= 0 && lon12 < 0 ? -1 : 0);
    return cross;
  }

  p.PolygonArea = function(earth, polyline) {
    this._earth = earth;
    this._area0 = 4 * Math.PI * earth._c2;
    this._polyline = !polyline ? false : polyline;
    this._mask =  g.DISTANCE | (this._polyline ? 0 : g.AREA);
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
      var t = this._earth.Inverse(this._lat1, this._lon1, lat, lon, this._mask);
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

  // return number, perimeter, area
  p.PolygonArea.prototype.Compute = function(reverse, sign) {
    var vals = {number: this._num};
    if (this._num < 2) {
      vals.perimeter = 0;
      if (!this._polyline)
	vals.area = 0;
      return vals;
    }
    if (this._polyline) {
      vals.perimeter = this._perimetersum.Sum();
      return vals;
    }
    var t = this._earth.Inverse(this._lat1, this._lon1, this._lat0, this._lon0,
				this._mask);
    vals.perimeter = this._perimetersum.Sum(t.s12);
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
	tempsum.Add( -this._area0 );
      else if (tempsum.Sum() <= -this._area0/2)
	tempsum.Add( +this._area0 );
    } else {
      if (tempsum.Sum() >= this._area0)
	tempsum.Add( -this._area0 );
      else if (tempsum < 0)
	tempsum.Add( -this._area0 );
    }
    vals.area = tempsum.Sum();
    return vals;
  }

  // return number, perimeter, area
  p.TestCompute = function(lat, lon, reverse, sign) {
    var vals = {number: this._num + 1};
    if (this._num == 0) {
      vals.perimeter = 0;
      if (!this._polyline)
	vals.area = 0;
      return vals;
    }
    vals.perimeter = this._perimetersum.Sum();
    var tempsum = this._polyline ? 0 : this._areasum.Sum();
    var crossings = this._crossings;
    var t;
    for (var i = 0; i < this._polyline ? 1 : 2; ++i) {
      t = this._earth.Inverse
      (i == 0 ? this._lat1 : lat, i == 0 ? this._lon1 : lon,
       i != 0 ? this._lat0 : lat, i != 0 ? this._lon0 : lon,
       this._mask);
      vals.perimeter += t.s12;
      if (!this._polyline) {
	tempsum += t.S12;
	crossings += p.transit(i == 0 ? this._lon1 : lon,
			       i != 0 ? this._lon0 : lon);
      }
    }

    if (this._polyline)
      return vals;

    if (crossings & 1)
      tempsum += (tempsum < 0 ? 1 : -1) * this._area0/2;
    // area is with the clockwise sense.  If !reverse convert to
    // counter-clockwise convention.
    if (!reverse)
      tempsum *= -1;
    // If sign put area in (-area0/2, area0/2], else put area in [0, area0)
    if (sign) {
      if (tempsum > this._area0/2)
	tempsum -= this._area0;
      else if (tempsum <= -this._area0/2)
	tempsum += this._area0;
    } else {
      if (tempsum >= this._area0)
	tempsum -= this._area0;
      else if (tempsum < 0)
	tempsum += this._area0;
    }
    vals.area = tempsum;
    return vals;
  }

  p.Area = function(earth, points, polyline) {
    var poly = new p.PolygonArea(earth, polyline);
    for (var i = 0; i < points.length; ++i)
      poly.AddPoint(points[i].lat, points[i].lon);
    return poly.Compute(false, true);
  }

})();

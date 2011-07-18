/**
 * Interface.js
 * Javascript interface routines for GeographicLib.
 *
 **********************************************************************
 * GeographicLib.Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2, outmask);
 * GeographicLib.Geodesic.WGS84.Direct(lat1, lon1, azi1, s12, outmask);
 *
 * perform the basic geodesic calculations.  These return an object with
 * (some) of the following fields:
 *
 *   lat1 latitude of point 1
 *   lon1 longitude of point 1
 *   azi1 azimuth of line at point 1
 *   lat2 latitude of point 2
 *   lon2 longitude of point 2
 *   azi2 azimuth of line at point 2
 *   s12 distance from 1 to 2
 *   a12 arc length on auxiliary sphere from 1 to 2
 *   m12 reduced length of geodesic
 *   M12 geodesic scale 2 relative to 1
 *   M21 geodesic scale 1 relative to 2
 *   S12 area between geodesic and equator
 *
 * outmask determines which fields get included and if outmask is
 * omitted, then only the basic geodesic fields are computed.  The mask
 * is an or'ed combination of the following values
 *
 *   GeographicLib.Geodesic.LATITUDE
 *   GeographicLib.Geodesic.LONGITUDE
 *   GeographicLib.Geodesic.AZIMUTH
 *   GeographicLib.Geodesic.DISTANCE
 *   GeographicLib.Geodesic.REDUCEDLENGTH
 *   GeographicLib.Geodesic.GEODESICSCALE
 *   GeographicLib.Geodesic.AREA
 *   GeographicLib.Geodesic.ALL
 *
 **********************************************************************
 * GeographicLib.Geodesic.WGS84.Path(lat1, lon1, lat2, lon2, ds12, maxk);
 *
 * splits a geodesic line into k equal pieces which are no longer than
 * about ds12 (but k cannot exceed maxk, default 20), and returns a
 * vector of length k + 1 of objects with fields
 *
 *   lat, lon, azi
 *
 **********************************************************************
 * GeographicLib.Geodesic.WGS84.Area(points, polyline);
 *
 * computes the area of a polygon with vertices given by an array
 * points, each of whose elements contains lat and lon fields.	The
 * function returns an object with fields.
 *
 *   number, perimeter, area
 *
 * There is no need to "close" the polygon.  If polyline (default =
 * false) is true, that the points denote a polyline and its length is
 * returned as the perimeter (and the area is not calculated).
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * $Id$
 **********************************************************************/

(function() {
  var m = GeographicLib.Math;
  var g = GeographicLib.Geodesic;
  var l = GeographicLib.GeodesicLine;

  g.Geodesic.CheckPosition = function(lat, lon) {
    if (!(Math.abs(lat) <= 90))
      throw new Error("latitude " + lat + " not in [-90, 90]");
    if (!(lon >= -180 && lon <= 360))
      throw new Error("longitude " + lon + " not in [-180, 360]");
    return g.AngNormalize(lon);
  }

  g.Geodesic.CheckAzimuth = function(azi) {
    if (!(azi >= -180 && azi <= 360))
      throw new Error("longitude " + azi + " not in [-180, 360]");
  }

  g.Geodesic.CheckDistance = function(s) {
    if (!(isFinite(s)))
      throw new Error("distance " + s + " not a finite number");
  }

  g.Geodesic.prototype.Inverse = function(lat1, lon1, lat2, lon2, outmask) {
    if (!outmask) outmask = g.DISTANCE | g.AZIMUTH;
    lon1 = g.Geodesic.CheckPosition(lat1, lon1);
    lon2 = g.Geodesic.CheckPosition(lat2, lon2);

    var result = this.GenInverse(lat1, lon1, lat2, lon2, outmask);
    result.lat1 = lat1; result.lon1 = lon1;
    result.lat2 = lat2; result.lon2 = lon2;
    
    return result;
  }

  g.Geodesic.prototype.Direct = function(lat1, lon1, azi1, s12, outmask) {
    if (!outmask) outmask = g.LATITUDE | g.LONGITUDE | g.AZIMUTH;
    lon1 = g.Geodesic.CheckPosition(lat1, lon1);
    azi1 = g.Geodesic.CheckAzimuth(azi1);
    g.Geodesic.CheckDistance(s12);

    var result = this.GenDirect(lat1, lon1, azi1, false, s12, outmask);
    result.lat1 = lat1; result.lon1 = lon1;
    result.azi1 = azi1; result.s12 = s12;
    
    return result;
  }

  g.Geodesic.prototype.Path = function(lat1, lon1, lat2, lon2, ds12, maxk) {
    var t = this.Inverse(lat1, lon1, lat2, lon2);
    if (!maxk) maxk = 20;
    if (!(ds12 > 0))
      throw new Error("ds12 must be a positive number")
    var
    k = Math.max(1, Math.min(maxk, Math.ceil(t.s12/ds12))),
    points = new Array(k + 1);
    points[0] = {lat: t.lat1, lon: t.lon1, azi: t.azi1};
    points[k] = {lat: t.lat2, lon: t.lon2, azi: t.azi2};
    if (k > 1) {
      var line = new l.GeodesicLine(this, t.lat1, t.lon1, t.azi1,
				    g.LATITUDE | g.LONGITUDE | g.AZIMUTH),
      da12 = t.a12/k;
      var vals;
      for (var i = 1; i < k; ++i) {
	vals =
	line.GenPosition(true, i * da12, g.LATITUDE | g.LONGITUDE | g.AZIMUTH);
	points[i] = {lat: vals.lat2, lon: vals.lon2, azi: vals.azi2};
      }
    }
    return points;
  }

  g.Geodesic.prototype.Circle = function(lat1, lon1, azi1, s12, k) {
    if (!(Math.abs(lat1) <= 90))
      throw new Error("lat1 must be in [-90, 90]");
    if (!(lon1 >= -180 && lon1 <= 360))
      throw new Error("lon1 must be in [-180, 360]");
    if (!(azi1 >= -180 && azi1 <= 360))
      throw new Error("azi1 must be in [-180, 360]");
    if (!(isFinite(s12)))
      throw new Error("s12 must be a finite number");
    if (lon1 >= 180) lon1 -= 360;
    if (azi1 >= 180) azi1 -= 360;
    if (!k || k < 4) k = 24;
    var points = new Array(k + 1);
    var vals;
    for (var i = 0; i <= k; ++i) {
      var azi1a = azi1 + (k - i) * 360 / k; // Traverse circle counter-clocwise
      if (azi1a >= 180) azi1a -= 360;
      vals =
	this.GenDirect(lat1, lon1, azi1a, false, s12, g.LATITUDE | g.LONGITUDE);
      points[i] = {lat: vals.lat2, lon: vals.lon2};
    }
    return points;
  }

  g.Geodesic.prototype.Envelope = function(lat1, lon1, k) {
    if (!(Math.abs(lat1) <= 90))
      throw new Error("lat1 must be in [-90, 90]");
    if (!(lon1 >= -180 && lon1 <= 360))
      throw new Error("lon1 must be in [-180, 360]");
    if (lon1 >= 180) lon1 -= 360;
    if (!k || k < 4) k = 24;
    var points = new Array(k + 1);
    var vals, line, s12, j;
    for (var i = 0; i <= k; ++i) {
      var azi1 = -180 + i * 360 / k;
      line = new l.GeodesicLine(this, lat1, lon1, azi1, 
				g.LATITUDE | g.LONGITUDE | g.DISTANCE_IN |
				g.DISTANCE | g.REDUCEDLENGTH | g.GEODESICSCALE);
      vals = line.GenPosition(true, 180,
			      g.DISTANCE | g.REDUCEDLENGTH | g.GEODESICSCALE);
      j = 0;
      while (true) {
	// Solve m12(s12) = 0 by Newton's method using dm12/ds12 = M21
	s12 = vals.s12 - vals.m12/vals.M21;
	if (Math.abs(vals.m12) < 0.01 || ++j > 10)
	  break;
	vals = line.GenPosition(false, s12,
				g.DISTANCE | g.REDUCEDLENGTH | g.GEODESICSCALE);
      }
      vals = line.GenPosition(false, s12, g.LATITUDE | g.LONGITUDE);
      points[i] = {lat: vals.lat2, lon: vals.lon2};
    }
    return points;
  }

  g.Geodesic.prototype.Area = function(points, polyline) {
    return GeographicLib.PolygonArea.Area(this, points, polyline);
  }

})();


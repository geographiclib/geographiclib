Geodesic routines from GeographicLib
====================================

{@link module:GeographicLib}
{@link module:GeographicLib/Accumulator}
{@link module:GeographicLib/Accumulator.Accumulator}
{@link module:GeographicLib/Math}
{@link module:GeographicLib/Math.degree}
{@link GeographicLib}
{@link GeographicLib/Accumulator.Accumulator}


This is GeographicLib (http://geographiclib.sourceforge.net/) node.js port.

Installation:

```javascript
$ nmp install geographiclib
```

Usage:

You should simply require('geographiclib') and then use it like using js
version of GeographicLib
(http://geographiclib.sourceforge.net/scripts/geographiclib.js):

http://geographiclib.sourceforge.net/html/other.html#javascript

For library documentation refer
[here](http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geodesic.html).

Some examples below:

```javascript
var GeographlicLib = require('geographiclib'),
    geod = GeographicLib.Geodesic.WGS84; // using WGS84 ellipsoid

geo.Inverse(source_lat, source_lon, destination_lat, destination_lon);
/**
 * outputs in (numbers are just example):
 * {
 *   a12: 0.043512234,
 *   s12: 15432.779809263175,
 *   azi1: 10.450357251089693,
 *   azi2: 10.467249496997848,
 *   lat1: source_lat,
 *   lon1: source_lon,
 *   lat2: destination_lat,
 *   lon2: destination_lon
 * }
 */

geo.Direct(source_lat, source_lon, azimuth, distance);
/**
 * outputs in (numbers are just example):
 * {
 *   lon2: destination_lat,
 *   lat2: destination_lon,
 *   azi2: 10.467249496997848,
 *   a12: 0.045145534,
 *   lat1: source_lat,
 *   lon1: source_lon,
 *   azi1: azimuth,
 *   s12: distance
 * }
 */

```

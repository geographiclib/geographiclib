
````
var GeographicLib = require('./geographiclib'),
    wgs84 = GeographicLib.Geodesic.WGS84,
    international = new GeographicLib.Geodesic.Geodesic(6378388, 1/297);
geod.Inverse(1,2,3,4);
international.Inverse(1,2,3,4, GeographicLib.Geodesic.ALL);
````

See {@tutorial tutorial-geodesic}.
perform the basic geodesic calculations.  These return an object with
(some) of the following fields set:

* lat1 = &phi;<sub>1</sub>, latitude of point 1 (degrees)
* lon1 = &lambda;<sub>1</sub>, longitude of point 1 (degrees)
* azi1 = &alpha;<sub>1</sub>, azimuth of line at point 1 (degrees)
* lat2 = &phi;<sub>2</sub>, latitude of point 2 (degrees)
* lon2 = &lambda;<sub>2</sub>, longitude of point 2 (degrees)
* azi2 = &alpha;<sub>2</sub>, (forward) azimuth of line at point 2 (degrees)
* s12 = *s*<sub>12</sub>, distance from 1 to 2 (meters)
* a12 = &sigma;<sub>12</sub>, arc length on auxiliary sphere from 1 to 2
  (degrees)
* m12 = *m*<sub>12</sub>, reduced length of geodesic (meters)
* M12 = *M*<sub>12</sub>, geodesic scale at 2 relative to 1 (dimensionless)
* M21 = *M*<sub>21</sub>, geodesic scale at 1 relative to 2 (dimensionless)
* S12 = *S*<sub>12</sub>, area between geodesic and equator
  (meters<sup>2</sup>)

outmask determines which fields get included and if outmask is omitted,
then STANDARD is assumed.  The mask is an or'ed combination of the
following values

* GeographicLib.Geodesic.LATITUDE, compute latitude
* GeographicLib.Geodesic.LONGITUDE, compute longitude
* GeographicLib.Geodesic.AZIMUTH, compute azimuth
* GeographicLib.Geodesic.DISTANCE, compute distance
* GeographicLib.Geodesic.STANDARD (all of the above)
* GeographicLib.Geodesic.DISTANCE_IN, see below
* GeographicLib.Geodesic.REDUCEDLENGTH, compute reduced length
* GeographicLib.Geodesic.GEODESICSCALE, compute geodesic scales
* GeographicLib.Geodesic.AREA, compute area
* GeographicLib.Geodesic.ALL (all of the above)
* GeographicLib.Geodesic.LONG_UNROLL, see below

The input arguments are included in the returned objects.  Latitudes
outside of the allowed range [&minus;90&deg;, 90&deg;] are replaced by
NaNs.  Azimuths are reduced to the range [&minus;180&deg;, 180&deg;).

GeographicLib.Geodesic.DISTANCE_IN is a capability provided to the
GeographicLine.GeodesicLine constructor.  It allows the position on the
line to specified in terms of distance.  (Without this, the position can
only be specified in terms of the arc length.)

GeographicLib.Geodesic.LONG_UNROLL controls the treatment of longitude.
If it is not set then the lon1 and lon2 fields are both reduced to the
range [&minus;180&deg;, 180&deg;).  If it is set, then lon1 is as given
in the function call and lon2 &minus; lon1 determines how many times and
in what sense the geodesic has encircled the ellipsoid.
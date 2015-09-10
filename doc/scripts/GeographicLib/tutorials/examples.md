````
var GeographicLib = require('./geographiclib'),
    g = GeographicLib.Geodesic,
    c = GeographicLib.Constants,
    geog = new g.Geodesic(c.WGS84.a, c.WGS84.f);
geog.Inverse(-30, 20, 29.5, 199.5);
// -->
// { lat1: -30,
//   lat2: 29.5,
//   lon1: 0,
//   lon2: -5,
//   a12: 59.526700345561785,
//   s12: 6606149.877725766,
//   azi1: -5.06894183621681,9
//   azi2: -5.0437842312297185 }
geog.Inverse(-30, 0, 29.5, 179.5, g.ALL);
````

## Geodesic routines from GeographicLib

This documentation applies to version 1.46.

The documentation for other versions is available
at <tt>http://geographiclib.sourceforge.net/m.nn/js</tt> for versions
numbers <tt>m.nn</tt> &ge; 1.45.

Licensed under the MIT/X11 License; see
[LICENSE.txt](http://geographiclib.sourceforge.net/html/LICENSE.txt).

### Installation

This library is a JavaScript implementation of the geodesic routines
from [GeographicLib](http://geographiclib.sourceforge.net).  This solves the
direct and inverse geodesic problems for an ellipsoid of revolution.

The library can be used in [node](https://nodejs.org) by first
installing the
[geographiclib node package](https://www.npmjs.com/package/geographiclib)
with [npm](https://www.npmjs.com)
```bash
$ npm install geographiclib
$ node
> var GeographicLib = require("geographiclib");
```
The npm package includes a test suite.  Run this by
```bash
$ cd node_modules/geograliblib
$ npm test
```

Alternatively, you can use it in client-side JavaScript, by including in
your HTML page
```html
<script type="text/javascript"
        src="http://geographiclib.sourceforge.net/scripts/geographiclib.js">
</script>
```
Both of these prescriptions define a {@link GeographicLib} namespace.

### Examples

Now geodesic calculations can be carried out, for example,
```javascript
var geod = GeographicLib.Geodesic.WGS84, r;

// Find the distance from Wellington, NZ (41.32S, 174.81E) to
// Salamanca, Spain (40.96N, 5.50W)...
r = geod.Inverse(-41.32, 174.81, 40.96, -5.50);
console.log("The distance is " + r.s12.toFixed(3) + " m.");
// This prints "The distance is 19959679.267 m."

// Find the point 20000 km SW of Perth, Australia (32.06S, 115.74E)...
r = geod.Direct(-32.06, 115.74, 225, 20000e3);
console.log("The position is (" +
            r.lat2.toFixed(8) + ", " + r.lon2.toFixed(8) + ").");
// This prints "The position is (32.11195529, -63.95925278)."
```
Two examples of this library in use are
* [A geodesic calculator](http://geographiclib.sourceforge.net/scripts/geod-calc.html)
* [Displaying geodesics on Google
  Maps](http://geographiclib.sourceforge.net/scripts/geod-google.html)

### More information
* {@tutorial 1-geodesics}
* {@tutorial 2-interface}
* {@tutorial 3-examples}

### Implementations in various languages
* {@link http://sourceforge.net/p/geographiclib/code/ci/release/tree/
    git repository}
* C++ (complete library):
  {@link http://geographiclib.sourceforge.net/html/index.html
    documentation},
  {@link https://sourceforge.net/projects/geographiclib/files/distrib
    download};
* C (geodesic routines):
  {@link http://geographiclib.sourceforge.net/html/C/index.html
    documentation}, also included with recent versions of
  {@link https://github.com/OSGeo/proj.4/wiki
    proj.4};
* Fortran (geodesic routines):
  {@link http://geographiclib.sourceforge.net/html/Fortran/index.html
    documentation};
* Java (geodesic routines):
  {@link http://repo1.maven.org/maven2/net/sf/geographiclib/GeographicLib-Java/
    Maven Central package},
  {@link http://geographiclib.sourceforge.net/html/java/index.html
    documentation};
* JavaScript (geodesic routines):
  {@link https://www.npmjs.com/package/geographiclib
    npm package},
  {@link http://geographiclib.sourceforge.net/html/js/index.html
    documentation};
* Python (geodesic routines):
  {@link http://pypi.python.org/pypi/geographiclib
    PyPI package},
  {@link http://geographiclib.sourceforge.net/html/python/index.html
    documentation};
* Matlab/Octave (geodesic and some other routines):
  {@link http://www.mathworks.com/matlabcentral/fileexchange/50605
    Matlab Central package},
  {@link http://www.mathworks.com/matlabcentral/fileexchange/50605/content/Contents.m
    documentation};
* C# (.NET wrapper for complete C++ library):
  {@link http://geographiclib.sourceforge.net/html/NET/index.html
    documentation}.

### Change log

* Version 1.46 (released 2016-02-15)
  * Fix bugs in PolygonArea.TestEdge (problem found by threepointone).
  * Add Geodesic.DirectLine, Geodesic.ArcDirectLine,
    Geodesic.GenDirectLine, Geodesic.InverseLine,
    GeodesicLine.SetDistance, GeodesicLine.SetArc,
    GeodesicLine.GenSetDistance, GeodesicLine.s13, GeodesicLine.a13.
  * More accurate inverse solution when longitude difference is close to
    180&deg;.

### Authors

* algorithms + js code: Charles Karney (charles@karney.com)
* node.js port: Yurij Mikhalevich (0@39.yt)

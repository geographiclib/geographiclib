Test data for GeographicLib
===========================

*NOTE:* It's not necessary to download these files to use
GeographicLib.

`GeodTest.dat.gz` and `GeodTest-short.dat.gz` provide test data for
geodesics.  See

> https://geographiclib.sourceforge.io/C++/doc/geodesic.html#testgeod

`TMCoords.dat.gz` provides test data for the transverse Mercator
projection.  See

> https://geographiclib.sourceforge.io/C++/doc/transversemercator.html#testmerc

`GeoidHeights.dat.gz` provides test data for geoid heights.  See

> https://geographiclib.sourceforge.io/C++/doc/geoid.html#testgeoid

`Geod3Test-v1.txt.gz` provides test data for geodesics on a triaxial
ellipsoid with semiaxes, [a, b, c] = [√2, 1, 1/√2].  This dataset is
described in

> https://doi.org/10.5281/zenodo.12510796

`Geod3Test-*-v1.txt.gz` provide corresponding test data for geodesics on
biaxial ellipsoids.  The coordinates for the endpoints match those in
`Geod3Test-v1.txt.gz` and, note well, that these are triaxial
ellipsoidal coordinates.  The semiaxes of the ellipsoids are

* `Geod3Test-obl-v1.txt.gz`: oblate limit, [a, b, c] = [1, 1, 1/2]
* `Geod3Test-pro-v1.txt.gz`: prolate limit, [a, b, c] = [2, 1, 1]
* `Geod3Test-sph*-v1.txt.gz`: spherical limit, b = 1, a → b, and c → b,
  such that k² = (b² − c²)/(a² − c²) is given by
    * `Geod3Test-spha-v1.txt.gz`: k² = 1
    * `Geod3Test-sphb-v1.txt.gz`: k² = 2/3
    * `Geod3Test-sphc-v1.txt.gz`: k² = 1/3
    * `Geod3Test-sphd-v1.txt.gz`: k² = 0

This is a Fortran implementation of the geodesic algorithms described in

  C. F. F. Karney,
  Algorithms for geodesics,
  J. Geodesy (2012);
  http://dx.doi.org/10.1007/s00190-012-0578-z
  Addenda: http://geographiclib.sf.net/geod-addenda.html

This implementation is more-or-less a direct transcription of the C++
class Geodesic and GeodesicLine in GeographicLib

  http://geographiclib.sf.net

I have intentionally stuck with Fortran 77 (omitting features which are
now deprecated).  So the library should compile correctly with any
extant Fortran compiler.

Inventory:

  geodesic.for -- the library
  direct.for, inverse.for -- sample driver programs
  CMakeLists.txt -- the cmake file for building the driver programs

See the comments at the beginning of geodesic.for for details for the
library.

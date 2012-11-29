This is a C implementation of the geodesic algorithms described in

  C. F. F. Karney,
  Algorithms for geodesics,
  J. Geodesy (2012);
  http://dx.doi.org/10.1007/s00190-012-0578-z
  Addenda: http://geographiclib.sf.net/geod-addenda.html

This implementation is-less a direct transcription of the C++ class
Geodesic and GeodesicLine in GeographicLib

  http://geographiclib.sf.net

I have intentionally stuck with ANSI C as described in B. W. Kernigan
and D. M. Ritchie, The C Programming Language, 2nd Ed. (Prentice Hall,
1988).  So the library should compile correctly with any extant C
compiler.

Inventory:

  geodesic.h, geodesic.c -- the library
  direct.c, inverse.c -- sample driver programs
  CMakeLists.txt -- the cmake file for building the driver programs

See the comments at the beginning of geodesic.h for details for the
library.

/* A simple test of the C library for geodesics.  This solves the "inverse
 * geodesic problem" by reading in lines with lat1, lon1, lat2, lon2 and
 * printing out lines with azi1, azi2, s12 (for the WGS84 ellipsoid). */

#include <stdio.h>
#include "geodesic.h"

#if defined(_MSC_VER)
/* Squelch warnings about scanf */
#  pragma warning (disable: 4996)
#endif

int main() {
  double a = 6378137, f = 1/298.257223563; /* WGS84 */
  double lat1, lon1, azi1, lat2, lon2, azi2, s12;
  struct geod_geodesic g;

  geod_init(&g, a, f);
  while (scanf("%lf %lf %lf %lf\n", &lat1, &lon1, &lat2, &lon2) == 4) {
    geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
    printf("%.15f %.15f %.10f\n", azi1, azi2, s12);
  }
  return 0;
}

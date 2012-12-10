/* A simple test of the C library for geodesics.  This computes the area
 * of a geodesic polygon by reading in up to 100 lines with the
 * coordinates (lat, lon) for each vertex, and printing out the number
 * of points, the perimeter, and the area (for the WGS84 ellipsoid). */

#include <stdio.h>
#include "geodesic.h"

#if defined(_MSC_VER)
/* Squelch warnings about scanf */
#  pragma warning (disable: 4996)
#endif

#define MAXPTS 100
int main() {
  double a = 6378137, f = 1/298.257223563; /* WGS84 */
  double lats[MAXPTS], lons[MAXPTS], A, P;
  int n = 0;
  struct geod_geodesic g;

  while (n < MAXPTS && scanf("%lf %lf\n", &lats[n], &lons[n]) == 2)
    ++n;
  geod_init(&g, a, f);
  geod_polygonarea(&g, lats, lons, n, &A, &P);
  printf("%d %.8f %.2f\n", n, P, A);
  return 0;
}

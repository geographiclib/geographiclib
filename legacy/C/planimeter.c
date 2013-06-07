/**
 * @file planimeter.c
 * @brief A test program for geod_polygonarea()
 **********************************************************************/

#include <stdio.h>
#include "geodesic.h"

#if defined(_MSC_VER)
/* Squelch warnings about scanf */
#  pragma warning (disable: 4996)
#endif
#define MAXPTS 100000

/**
 * A simple program to compute the area of a geodesic polygon.
 *
 * This program reads in up to 100000 lines with lat, lon for each vertex
 * of a polygon.  At the end of input, the program prints the number of
 * vertices, the perimeter of the polygon and its area (for the WGS84
 * ellipsoid).
 **********************************************************************/

int main() {
  double a = 6378137, f = 1/298.257223563; /* WGS84 */
  double lats[MAXPTS], lons[MAXPTS], A, P;
  int n = 0;
  struct geod_geodesic g;

  while (n < MAXPTS && scanf("%lf %lf\n", &lats[n], &lons[n]) == 2)
    ++n;
  geod_init(&g, a, f);
  geod_polygonarea(&g, lats, lons, n, &A, &P);
  printf("%d %.8f %.3f\n", n, P, A);
  return 0;
}

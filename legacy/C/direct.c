#include <stdio.h>
#include "geodesic.h"

#if defined(_MSC_VER)
// Squelch warnings about scanf
#  pragma warning (disable: 4996)
#endif

int main() {
  double a = 6378137, f = 1/298.257223563; /* WGS84 */
  double lat1, lon1, azi1, lat2, lon2, azi2, s12;
  struct Geodesic g;

  GeodesicInit(&g, a, f);
  while (scanf("%lf %lf %lf %lf\n", &lat1, &lon1, &azi1, &s12) == 4) {
    Direct(&g, lat1, lon1, azi1, s12, &lat2, &lon2, &azi2);
    printf("%.15f %.15f %.15f\n", lat2, lon2, azi2);
    }
  return 0;
}

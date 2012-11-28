#include <stdio.h>
#include <math.h>
#include "geodesic.h"

#if defined(_MSC_VER)
// Squelch warnings about scanf
#  pragma warning (disable: 4996)
#endif

int main() {
  double a = 6378137, f = 1/298.257223563; /* WGS84 */
  double lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, S12;
  struct Geodesic g;

  GeodesicInit(&g, a, f);

  while (scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &lat1, &lon1, &azi1,
               &lat2, &lon2, &azi2, &s12, &a12, &m12, &S12 ) == 10) {
    double lat1a, lon1a, azi1a, lat2a, lon2a, azi2a, s12a, a12a, m12a, S12a;
    a12a = GenInverse(&g, lat1, lon1, lat2, lon2, &s12a, &azi1a, &azi2a,
                      &m12a, 0, 0, &S12a);
    printf("%.8f %.8f %.3f %.8f %.3f %.0f\n",
           azi1, azi2, s12, a12, m12, S12);
    printf("%.8f %.8f %.3f %.8f %.3f %.0f\n",
           azi1a, azi2a, s12a, a12a, m12a, S12a);
    a12a = GenDirect(&g, lat1, lon1, azi1, 0, s12, &lat2a, &lon2a, &azi2a,
                     0, &m12, 0, 0, &S12a);
    printf("%.8f %.8f %.8f %.8f %.3f %.0f\n",
           lat2, lon2, azi2, a12, m12, S12);
    printf("%.8f %.8f %.8f %.8f %.3f %.0f\n",
           lat2a, lon2a, azi2a, a12a, m12a, S12a);
    a12a = GenDirect(&g, lat2, lon2, azi2, 0, -s12, &lat1a, &lon1a, &azi1a,
                     0, &m12, 0, 0, &S12a);
    printf("%.8f %.8f %.8f %.8f %.3f %.0f\n",
           lat1, lon1, azi1, a12, m12, S12);
    printf("%.8f %.8f %.8f %.8f %.3f %.0f\n",
           lat1a, lon1a, azi1a, -a12a, -m12a, -S12a);
    GenDirect(&g, lat1, lon1, azi1, 1, a12, &lat2a, &lon2a, &azi2a,
              &s12a, &m12, 0, 0, &S12a);
    printf("%.8f %.8f %.8f %.3f %.3f %.0f\n",
           lat2, lon2, azi2, s12, m12, S12);
    printf("%.8f %.8f %.8f %.3f %.3f %.0f\n",
           lat2a, lon2a, azi2a, s12a, m12a, S12a);
  }
  return 0;
}

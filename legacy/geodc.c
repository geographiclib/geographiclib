#include <stdio.h>
#include "geodesicc.h"

int main() {
  double lat1, lon1, lat2, lon2;
  double a12, s12, azi1, azi2, m12, M12, M21, S12;
  double a = 6378137, f = 1/298.257223563;
  struct Geodesic g;
  GeodesicInit(&g, a, f);
  unsigned outmask = ALL;
  unsigned inverse = 0;
  if (inverse)
    while (scanf("%lf %lf %lf %lf\n", &lat1, &lon1, &lat2, &lon2) == 4) {
      a12 = GenInverse(&g, lat1, lon1, lat2, lon2, outmask,
		       &s12, &azi1, &azi2, &m12, &M12, &M21, &S12);
      printf("%.10f %.10f %.5f\n", azi1, azi2, s12);
    }
  else
    while (scanf("%lf %lf %lf %lf\n", &lat1, &lon1, &azi1, &s12) == 4) {
      double s12a,
      a12 = GenDirect(&g, lat1, lon1, azi1, 0, s12, outmask,
		      &lat2, &lon2, &azi2, &s12a, &m12, &M12, &M21, &S12);
      printf("%.10f %.10f %.10f\n", lat2, lon2, azi2);
    }
  return 0;
}

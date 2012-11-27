#include <stdio.h>
#include <math.h>
#include "geodesicc.h"

#if defined(_MSC_VER)
// Squelch warnings about scanf
#  pragma warning (disable: 4996)
#endif

int main() {
  double lat1, lon1, lat2, lon2;
  double a12, s12, azi1, azi2, m12, M12, M21, S12;
  double a = 6378137, f = 1/298.257223563;
  struct Geodesic g;
  unsigned inverse = 1, test = 0;

  GeodesicInit(&g, a, f);
  if (test) {
    double erri = 0, errd = 0;
    int i = 0, ierri = 0, ierrd = 0;
    while (scanf("%lf %lf %lf %lf %lf %lf %lf\n", &lat1, &lon1, &azi1,
                 &lat2, &lon2, &azi2, &s12) == 7) {
      double lat1a, lon1a, lat2a, lon2a, s12a, azi1a, azi2a;
      i += 1;
      Inverse(&g, lat1, lon1, lat2, lon2, &s12a, &azi1a, &azi2a);
      s12a = fabs(s12a - s12);
      if (s12a > erri) {
        erri = s12a;
        ierri = i;
      }
      Direct(&g, lat1, lon1, azi1, s12, &lat2a, &lon2a, &azi2a);
      Inverse(&g, lat2, lon2, lat2a, lon2a, &s12a, &azi1a, &azi2a);
      if (s12a > errd) {
        errd = s12a;
        ierrd = i;
      }
      Direct(&g, lat2, lon2, azi2, -s12, &lat1a, &lon1a, &azi1a);
      Inverse(&g, lat1, lon1, lat1a, lon1a, &s12a, &azi1a, &azi2a);
      if (s12a > errd) {
        errd = s12a;
        ierrd = i;
      }
    }
    printf("%d %.2f %d %.2f\n", ierrd, errd * 1e9, ierri, erri * 1e9);
  } else if (inverse)
    while (scanf("%lf %lf %lf %lf\n", &lat1, &lon1, &lat2, &lon2) == 4) {
      a12 = GenInverse(&g, lat1, lon1, lat2, lon2,
                       &s12, &azi1, &azi2, &m12, &M12, &M21, &S12);
      printf("%.15f %.15f %.10f %.15f\n", azi1, azi2, s12, a12);
    }
  else
    while (scanf("%lf %lf %lf %lf\n", &lat1, &lon1, &azi1, &s12) == 4) {
      double s12a,
        a12 = GenDirect(&g, lat1, lon1, azi1, 0, s12,
                        &lat2, &lon2, &azi2, &s12a, &m12, &M12, &M21, &S12);
      printf("%.15f %.15f %.15f %.15f\n", lat2, lon2, azi2, a12);
    }
  return 0;
}

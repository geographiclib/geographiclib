#if !defined(GEODESICC_H)
#define GEODESICC_H 1

#define GEOGRAPHICLIB_GEODESIC_ORDER 6
#define nA1_   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC1_   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC1p_  GEOGRAPHICLIB_GEODESIC_ORDER
#define nA2_   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC2_   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3_   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3x_  nA3_
#define nC3_   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC3x_  ((nC3_ * (nC3_ - 1)) / 2)
#define nC4_   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC4x_  ((nC4_ * (nC4_ + 1)) / 2)

struct Geodesic {
    double a, f, f1, e2, ep2, n, b, c2, etol2;
    double A3x[nA3x_], C3x[nC3x_], C4x[nC4x_];
};

struct GeodesicLine {
    double lat1, lon1, azi1;
    double a, f, b, c2, f1, salp0, calp0, k2,
      salp1, calp1, ssig1, csig1, dn1, stau1, ctau1, somg1, comg1,
      A1m1, A2m1, A3c, B11, B21, B31, A4, B41;
  /* index zero elements of C1a, C1pa, C2a, C3a are unused */
    double C1a[nC1_ + 1], C1pa[nC1p_ + 1], C2a[nC2_ + 1], C3a[nC3_],
      C4a[nC4_];    /* all the elements of C4a are used */
    unsigned caps;
};

    enum captype {
      CAP_NONE = 0U,
      CAP_C1   = 1U<<0,
      CAP_C1p  = 1U<<1,
      CAP_C2   = 1U<<2,
      CAP_C3   = 1U<<3,
      CAP_C4   = 1U<<4,
      CAP_ALL  = 0x1FU,
      OUT_ALL  = 0x7F80U
    };
    enum mask {
      /**
       * No capabilities, no output.
       * @hideinitializer
       **********************************************************************/
      NONE          = 0U,
      /**
       * Calculate latitude \e lat2.  (It's not necessary to include this as a
       * capability to GeodesicLine because this is included by default.)
       * @hideinitializer
       **********************************************************************/
      LATITUDE      = 1U<<7  | CAP_NONE,
      /**
       * Calculate longitude \e lon2.
       * @hideinitializer
       **********************************************************************/
      LONGITUDE     = 1U<<8  | CAP_C3,
      /**
       * Calculate azimuths \e azi1 and \e azi2.  (It's not necessary to
       * include this as a capability to GeodesicLine because this is included
       * by default.)
       * @hideinitializer
       **********************************************************************/
      AZIMUTH       = 1U<<9  | CAP_NONE,
      /**
       * Calculate distance \e s12.
       * @hideinitializer
       **********************************************************************/
      DISTANCE      = 1U<<10 | CAP_C1,
      /**
       * Allow distance \e s12 to be used as input in the direct geodesic
       * problem.
       * @hideinitializer
       **********************************************************************/
      DISTANCE_IN   = 1U<<11 | CAP_C1 | CAP_C1p,
      /**
       * Calculate reduced length \e m12.
       * @hideinitializer
       **********************************************************************/
      REDUCEDLENGTH = 1U<<12 | CAP_C1 | CAP_C2,
      /**
       * Calculate geodesic scales \e M12 and \e M21.
       * @hideinitializer
       **********************************************************************/
      GEODESICSCALE = 1U<<13 | CAP_C1 | CAP_C2,
      /**
       * Calculate area \e S12.
       * @hideinitializer
       **********************************************************************/
      AREA          = 1U<<14 | CAP_C4,
      /**
       * All capabilities.  Calculate everything.
       * @hideinitializer
       **********************************************************************/
      ALL           = OUT_ALL| CAP_ALL
    };


void GeodesicInit(struct Geodesic* g, double a, double f);
double GenInverse(const struct Geodesic* g,
		double lat1, double lon1, double lat2, double lon2,
		unsigned outmask,
		double* ps12, double* pazi1, double* pazi2,
		double* pm12, double* pM12, double* pM21, double* pS12);

void GeodesicLineInit(struct GeodesicLine* l,
                      const struct Geodesic* g,
                      double lat1, double lon1, double azi1, unsigned caps);

double GenPosition(const struct GeodesicLine* l,
                   unsigned arcmode, double s12_a12,
                   unsigned outmask,
                   double* plat2, double* plon2, double* pazi2,
                   double* ps12, double* pm12,
                   double* pM12, double* pM21,
                   double* pS12);

double GenDirect(const struct Geodesic* g,
                 double lat1, double lon1, double azi1,
                 unsigned arcmode, double s12_a12, unsigned outmask,
                 double* plat2, double* plon2, double* pazi2,
                 double* ps12, double* pm12, double* pM12, double* pM21,
                 double* pS12);

#endif

#if !defined(GEODESICC_H)
#define GEODESICC_H 1

#if defined(__cplusplus)
extern "C" {
#endif

  struct Geodesic {
    double a, f, f1, e2, ep2, n, b, c2, etol2;
    double A3x[6], C3x[15], C4x[21];
  };

  struct GeodesicLine {
    double lat1, lon1, azi1;
    double a, f, b, c2, f1, salp0, calp0, k2,
      salp1, calp1, ssig1, csig1, dn1, stau1, ctau1, somg1, comg1,
      A1m1, A2m1, A3c, B11, B21, B31, A4, B41;
    /* index zero elements of C1a, C1pa, C2a, C3a are unused */
    double C1a[6+1], C1pa[6+1], C2a[6+1], C3a[6],
      C4a[6];                   /* all the elements of C4a are used */
    unsigned caps;
  };

  void GeodesicInit(struct Geodesic* g, double a, double f);

  void Direct(const struct Geodesic* g,
              double lat1, double lon1, double azi1, double s12,
              double* plat2, double* plon2, double* pazi2);
  void Inverse(const struct Geodesic* g,
               double lat1, double lon1, double lat2, double lon2,
               double* ps12, double* pazi1, double* pazi2);

  double GenDirect(const struct Geodesic* g,
                   double lat1, double lon1, double azi1,
                   int arcmode, double s12_a12,
                   double* plat2, double* plon2, double* pazi2,
                   double* ps12, double* pm12, double* pM12, double* pM21,
                   double* pS12);
  double GenInverse(const struct Geodesic* g,
                    double lat1, double lon1, double lat2, double lon2,
                    double* ps12, double* pazi1, double* pazi2,
                    double* pm12, double* pM12, double* pM21, double* pS12);
  void GeodesicLineInit(struct GeodesicLine* l,
                        const struct Geodesic* g,
                        double lat1, double lon1, double azi1, unsigned caps);
  double GenPosition(const struct GeodesicLine* l,
                     int arcmode, double s12_a12,
                     double* plat2, double* plon2, double* pazi2,
                     double* ps12, double* pm12,
                     double* pM12, double* pM21,
                     double* pS12);

  enum mask {
    /**
     * No capabilities, no output.
     **********************************************************************/
    NONE          = 0U,
    /**
     * Calculate latitude \e lat2.  (It's not necessary to include this as a
     * capability to GeodesicLine because this is included by default.)
     **********************************************************************/
    LATITUDE      = 1U<<7  | 0U,
    /**
     * Calculate longitude \e lon2.
     **********************************************************************/
    LONGITUDE     = 1U<<8  | 1U<<3,
    /**
     * Calculate azimuths \e azi1 and \e azi2.  (It's not necessary to
     * include this as a capability to GeodesicLine because this is included
     * by default.)
     **********************************************************************/
    AZIMUTH       = 1U<<9  | 0U,
    /**
     * Calculate distance \e s12.
     **********************************************************************/
    DISTANCE      = 1U<<10 | 1U<<0,
    /**
     * Allow distance \e s12 to be used as input in the direct geodesic
     * problem.
     **********************************************************************/
    DISTANCE_IN   = 1U<<11 | 1U<<0 | 1U<<1,
    /**
     * Calculate reduced length \e m12.
     **********************************************************************/
    REDUCEDLENGTH = 1U<<12 | 1U<<0 | 1U<<2,
    /**
     * Calculate geodesic scales \e M12 and \e M21.
     **********************************************************************/
    GEODESICSCALE = 1U<<13 | 1U<<0 | 1U<<2,
    /**
     * Calculate area \e S12.
     **********************************************************************/
    AREA          = 1U<<14 | 1U<<4,
    /**
     * All capabilities.  Calculate everything.
     **********************************************************************/
    ALL           = 0x7F80U| 0x1FU
  };

#if defined(__cplusplus)
}
#endif

#endif

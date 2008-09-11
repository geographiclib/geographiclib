# $Id$

Miscellaneous code for transformed geographic information.
                                Charles Karney <charles@karney.com>
                                http://charles.karney.info/geographic

tm.mac -- maxima code for arbitrary precision transverse Mercator
    projection.

ellint.mac -- maxima code various elliptic integrals, etc. (used by
    tm.mac).

UTM-fi.txt -- extend method for transverse Mercator projection given in
    JHS 154 include higher order terms and improve the formulas for
    meridian convergence and scale.

tmseries.mac -- maxima code to generate the coefficients for the series
    in UTM-fi.txt.

tmscale.mac -- maxima code to generate series approximations for the
    scale and meridian convergence in UTM-fi.txt.

revert.mac -- maxima code for the reversion of a series (used by
    tmseries.mac).

TMcoords.dat.gz -- set of about 1/4 million test points as gzipped text
    file.  The columns are:

      1 Latitude (deg)
      2 Longitude (deg)
      3 Transverse Mercator Easting (m)
      4 Transverse Mercator Northing (m)
      5 Meridian convergence (deg)
      6 Scale

    Latitude and Longitude are randomly and uniformly sampled from an
    octant of a sphere and then rounded to the nearest 10^-12 degrees.
    Columns 3-4 are determined by the function tm1 in tm.mac and
    rounding the nearest doubles.  The parameters for the transverse
    Mercator projection are:

       a = 6378137 m
       f = 1/298.257223563
       k0 = 0.9996
       Central Meridian = 0 deg
       False Easting = 0 m
       False Northing = 0 m


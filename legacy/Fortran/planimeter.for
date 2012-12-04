* A simple test of the C library for geodesics.  This computes the area
* of a geodesic polygon by reading in up to 100 lines with the
* coordinates (lat, lon) for each vertex, and printing out the number
* of points, the perimeter, and the area (for the WGS84 ellipsoid).

      program garea
      implicit none

      interface
        subroutine area(a, f, lats, lons, n, S, P)
        integer, intent(in) :: n
        double precision, intent(in) :: a, f, lats(n), lons(n)
        double precision, intent(out) :: S, P
        end subroutine area
      end interface

      integer maxpts
      parameter (maxpts = 100)
      double precision a, f, lats(maxpts), lons(maxpts), S, P
      integer n

* WGS84 values
      a = 6378137d0
      f = 1/298.257223563d0

      n = 0
 10   continue
      if (n .ge. maxpts) go to 20
      read(*, *, end=20, err=20) lats(n + 1), lons(n + 1)
      n = n + 1
      go to 10
 20   continue
      call area(a, f, lats, lons, n, S, P)
      print 30, n, P, S
 30   format(i3, 1x, f20.8, 1x, f19.2)
      stop
      end

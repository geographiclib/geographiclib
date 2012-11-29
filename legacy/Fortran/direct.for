* A simple test of the Fortran library for geodesics.  This solves the
* "direct geodesic problem" by reading in lines with lat1, lon1, azi1,
* s12 and printing out lines with lat2, lon2, azi2 (for the WGS84
* ellipsoid).

      program geodir
      implicit none

      interface
        subroutine direct(a, f, lat1, lon1, azi1, s12a12, arcmod,
     +      lat2, lon2, azi2, omask, a12s12, m12, MM12, MM21, SS12)
        double precision, intent(in) :: a, f, lat1, lon1, azi1, s12a12
        logical, intent(in) :: arcmod
        integer, intent(in) :: omask
        double precision, intent(out) :: lat2, lon2, azi2
* optional output (depending on omask)
        double precision, intent(out) :: a12s12, m12, MM12, MM21, SS12
        end subroutine direct
      end interface

      double precision a, f, lat1, lon1, azi1, lat2, lon2, azi2, s12,
     +    dummy
      logical arcmod
      integer omask

* WGS84 values
      a = 6378137d0
      f = 1/298.257223563d0

      arcmod = .false.
      omask = 0

 10   continue
      read(*, *, end=90, err=90) lat1, lon1, azi1, s12
      call direct(a, f, lat1, lon1, azi1, s12, arcmod,
     +    lat2, lon2, azi2, omask, dummy, dummy, dummy, dummy, dummy)
      print 20, lat2, lon2, azi2
 20   format(f20.15, 1x, f20.15, 1x, f20.15)
      go to 10
 90   continue
      stop
      end

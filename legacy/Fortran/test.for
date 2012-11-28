      program geod
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

        subroutine invers(a, f, lat1, lon1, lat2, lon2,
     +      s12, azi1, azi2, omask, a12, m12, MM12, MM21, SS12)
        double precision, intent(in) :: a, f, lat1, lon1, lat2, lon2
        integer, intent(in) :: omask
        double precision, intent(out) :: s12, azi1, azi2
* optional output (depending on omask)
        double precision, intent(out) :: a12, m12, MM12, MM21, SS12
        end subroutine invers
      end interface

      double precision a, f, lat1, lon1, azi1, lat2, lon2, azi2, s12,
     +    a12, m12, SS12, lat1a, lon1a, azi1a, lat2a, lon2a,
     +    azi2a, s12a, a12a, m12a, SS12a, dummy
      logical arcmod
      integer omask

* WGS84 values
      a = 6378137d0
      f = 1/298.257223563d0
      omask = 1+2+8


 10   continue
      read(*, *, end=90, err=90) lat1, lon1, azi1, lat2, lon2, azi2,
     +    s12, a12, m12, SS12

      call invers(a, f, lat1, lon1, lat2, lon2,
     +    s12a, azi1a, azi2a, omask, a12a, m12a, dummy, dummy, SS12a)
      print 20, azi1, azi2, s12, a12, m12, SS12
      print 20, azi1a, azi2a, s12a, a12a, m12a, SS12a
 20   format(f13.8, 1x, f13.8, 1x, f12.3, 1x, f13.8, 1x, f12.3,
     +    1x, f17.0)

      arcmod = .false.
      call direct(a, f, lat1, lon1, azi1, s12, arcmod,
     +    lat2a, lon2a, azi2a, omask, a12a, m12a, dummy, dummy, SS12a)
      print 30, lat2, lon2, azi2, a12, m12, SS12
      print 30, lat2a, lon2a, azi2a, a12a, m12a, SS12a
 30   format(f13.8, 1x, f13.8, 1x, f13.8, 1x, f13.8, 1x, f12.3,
     +    1x, f17.0)

      call direct(a, f, lat2, lon2, azi2, -s12, arcmod,
     +    lat1a, lon1a, azi1a, omask, a12a, m12a, dummy, dummy, SS12a)
      print 30, lat1, lon1, azi1, a12, m12, SS12
      print 30, lat1a, lon1a, azi1a, -a12a, -m12a, -SS12a

      arcmod = .true.
      call direct(a, f, lat1, lon1, azi1, a12, arcmod,
     +    lat2a, lon2a, azi2a, omask, s12a, m12a, dummy, dummy, SS12a)
      print 40, lat2, lon2, azi2, s12, m12, SS12
      print 40, lat2a, lon2a, azi2a, s12a, m12a, SS12a
 40   format(f13.8, 1x, f13.8, 1x, f13.8, 1x, f12.3, 1x, f12.3,
     +    1x, f17.0)
      go to 10
 90   continue
      stop
      end

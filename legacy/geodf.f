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

      double precision a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2,
     +    lat1a, lon1a, lat2a, lon2a, azi1a, azi2a, s12a,
     +    a12, m12, MM12, MM21, SS12
      double precision erri, errd
      integer omask, i, ierri, ierrd
      logical invp, arcmod, test

      a = 6378137d0
      f = 1/298.257223563d0
      invp = .false.
      test = .true.
      erri = 0
      errd = 0
      i = 0
      ierri = 0
      ierrd = 0

 10   continue
      if (test) then
        read(*, *, end=90, err=90) lat1, lon1, azi1, lat2, lon2, azi2,
     +      s12
        i = i + 1
        omask = 0
        arcmod = .false.
        call invers(a, f, lat1, lon1, lat2, lon2,
     +      s12a, azi1a, azi2a,
     +      omask, a12, m12, MM12, MM21, SS12)
        s12a = abs(s12a - s12)
        if (s12a .gt. erri) then
          erri = s12a
          ierri = i
        end if
        call direct(a, f, lat1, lon1, azi1, s12, arcmod,
     +      lat2a, lon2a, azi2a,
     +      omask, a12, m12, MM12, MM21, SS12)
        call invers(a, f, lat2, lon2, lat2a, lon2a,
     +      s12a, azi1a, azi2a,
     +      omask, a12, m12, MM12, MM21, SS12)
        if (s12a .gt. errd) then
          errd = s12a
          ierrd = i
        end if
        call direct(a, f, lat2, lon2, azi2, -s12, arcmod,
     +      lat1a, lon1a, azi1a,
     +      omask, a12, m12, MM12, MM21, SS12)
        call invers(a, f, lat1, lon1, lat1a, lon1a,
     +      s12a, azi1a, azi2a,
     +      omask, a12, m12, MM12, MM21, SS12)
        if (s12a .gt. errd) then
          errd = s12a
          ierrd = i
        end if
      else if (invp) then
        read(*, *, end=90, err=90) lat1, lon1, lat2, lon2
        omask = 15
        call invers(a, f, lat1, lon1, lat2, lon2,
     +      s12, azi1, azi2,
     +      omask, a12, m12, MM12, MM21, SS12)
        print 20, azi1, azi2, s12, a12, m12, MM12, MM21, SS12
 20     format(f20.15, 1x, f20.15, 1x, f19.10, 1x, f20.15, 1x, f19.10,
     +      1x, f20.17, 1x, f20.17, 1x, f20.3)
      else
        read(*, *, end=90, err=90) lat1, lon1, azi1, s12
        omask = 15
        arcmod = .false.
        call direct(a, f, lat1, lon1, azi1, s12, arcmod,
     +      lat2, lon2, azi2,
     +      omask, a12, m12, MM12, MM21, SS12)
        print 30, lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
 30     format(f20.15, 1x, f20.15, 1x, f20.15, 1x, f20.15, 1x, f19.10,
     +      1x, f20.17, 1x, f20.17, 1x, f20.3)
      end if
      go to 10
 90   continue
      if (test) print 40, ierrd, errd * 1d9, ierri, erri * 1d9
 40   format(i7, 1x, f20.2, 1x, i7, 1x, f20.2)
      stop
      end

cb::forward
c
      program forward
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:      forward
c version:   200208.19 
c author:    stephen j. frakes
c purpose:   to compute a geodetic forward (direct problem)
c            and then display output information
c
c input parameters:
c -----------------
c
c output parameters:
c ------------------
c
c local variables and constants:
c ------------------------------
c answer           user prompt response
c arc              meridional arc distance latitude p1 to p2 (meters)      
c b                semiminor axis polar (in meters)
c baz              azimuth back (in radians)
c blimit           geodetic distance allowed on ellipsoid (in meters)
c buff             input char buffer array
c dd,dm,ds         temporary values for degrees, minutes, seconds
c dmt,d_mt         char constants for units (in meters)        
c dd_max           maximum ellipsoid distance -1 (in meters)        
c edist            ellipsoid distance (in meters)
c elips            ellipsoid choice
c esq              eccentricity squared for reference ellipsoid
c faz              azimuth forward (in radians)
c filout           output file name
c finv             reciprocal flattening
c hem              hemisphere flag for lat & lon entry  
c ierror           error condition flag with d,m,s conversion
c lgh              length of buff() array
c option           user prompt response             
c r1,r2            temporary variables    
c ss               temporary value for ellipsoid distance
c tol              tolerance for conversion of seconds    
c
c name1            name of station one
c ld1,lm1,sl1      latitude  sta one - degrees,minutes,seconds
c ald1,alm1,sl1    latitude  sta one - degrees,minutes,seconds
c lat1sn           latitude  sta one - sign (+/- 1)
c d_ns1            latitude  sta one - char ('N','S')
c lod1,lom1,slo1   longitude sta one - degrees,minutes,seconds
c alod1,alom1,slo1 longitude sta one - degrees,minutes,seconds
c lon1sn           longitude sta one - sign (+/- 1)
c d_ew1            longitude sta one - char ('E','W')
c iaz1,maz1,saz1   forward azimuth   - degrees,minutes,seconds
c isign1           forward azimuth   - flag  (+/- 1)
c azd1,azm1,saz1   forward azimuth   - degrees,minutes,seconds
c iazsn            forward azimuth   - flag  (+/- 1)
c glat1,glon1      station one       - (lat & lon in radians )
c
c name2            name of station two
c ld2,lm2,sl2      latitude  sta two - degrees,minutes,seconds
c lat2sn           latitude  sta two - sign (+/- 1)
c d_ns2            latitude  sta two - char ('N','S')
c lod2,lom2,slo2   longitude sta two - degrees,minutes,seconds
c lon2sn           longitude sta two - sign (+/- 1)
c d_ew2            longitude sta two - char ('E','W')
c iaz2,maz2,saz2   back azimuth      - degrees,minutes,seconds
c isign2           back azimuth      - flag  (+/- 1)
c glat2,glon2      station two       - (lat & lon in radians )
c
c global variables and constants:
c -------------------------------
c a                semimajor axis equatorial (in meters)
c f                flattening
c pi               constant 3.14159....
c rad              constant 180.0/pi  
c
c    this module called by:  n/a
c
c    this module calls:      elipss, getrad, dirct1, todmsp
c    gethem, trim,   bufdms, gvalr8, gvali4, fixdms, gpnarc
c    datan,  write,  read,   dabs,   open,   stop
c
c    include files used:     n/a
c
c    common blocks used:     const, elipsoid
c
c    references:             see comments within subroutines
c
c    comments:
c
c********1*********2*********3*********4*********5*********6*********7**
c::modification history
c::1990mm.dd, sjf, creation of program           
c::199412.15, bmt, creation of program on viper
c::200203.08, crs, modified by c.schwarz to correct spelling of Clarke
c::                at request of Dave Doyle
c::200207.18, rws, modified i/o & standardized program documentation
c::                added subs trim, bufdms, gethem, gvali4, gvalr8      
c::200207.23, rws, added sub gpnarc
c::200208.15, rws, fixed an error in bufdms
c::              - renamed ellips to elipss "common error" with dirct2
c::              - added FAZ & BAZ to printed output 
c::200208.19, rws, added more error flags for web interface 
c::              - added logical nowebb                      
c::200208.xx, sjf, program version number 2.0                   
c********1*********2*********3*********4*********5*********6*********7**
ce::forward
c
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
      logical  nowebb
c
      character*1  answer,option,dmt,buff(50),hem
      character*6  d_ns1, d_ew1, d_ns2, d_ew2, d_mt
      character*30 filout,name1,name2,elips
c
      integer*4    ierror
      integer*4    lgh
c
      common/const/pi,rad
      common/elipsoid/a,f
c
c     ms_unix      0 = web version
c                  1 = ms_dos or unix
c
      parameter   ( ms_unix = 0 )
c
      pi=4.d0*datan(1.d0)
      rad=180.d0/pi
      dmt='m'
      d_mt='Meters'
c
      if( ms_unix.eq.1 )then
        nowebb = .true.
      else
        nowebb = .false.
      endif
c
    1 do 2 i=1,25
        write(*,*) '  '
    2 continue
 
    5 write(*,*) '  Program Forward  -  Version 2.0 '
      write(*,*) '  '
      write(*,*) '  Ellipsoid options : '
      write(*,*) '  '
      write(*,*) '  1) GRS80 / WGS84  (NAD83) '
      write(*,*) '  2) Clarke 1866    (NAD27) '
      write(*,*) '  3) Any other ellipsoid '
      write(*,*) '  '
      write(*,*) '  Enter choice : '
      read(*,10) option
   10 format(a1)
c
      if(option.eq.'1') then
        a=6378137.d0
        f=1.d0/298.25722210088d0
        elips='GRS80 / WGS84  (NAD83)'
      elseif(option.eq.'2') then
        a=6378206.4d0
        f=1.d0/294.9786982138d0
        elips='Clarke 1866    (NAD27)'
      elseif(option.eq.'3') then
        call elipss (elips)
      else
        write(*,*) '  Enter 1, 2, or 3 !   Try again --'
        goto 5
      endif
c
      esq = f*(2.0d0-f)
c
      r1  = 0.0d0
      r2  = pi/2.0d0
      call gpnarc ( a, f, esq, pi, r1, r2, arc )
c
c     compute the geodetic limit distance (blimit), it is equal
c     to twice the distance to the pole minus one meter
c
      blimit = 2.0d0*arc-1.0d0
c
c     maximum distance allowed on ellipsoid
c
      dd_max = blimit                  
c
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
c
   15 write(*,*) '  Enter First Station '
      write(*,16)
   16 format(18x,'(Separate D,M,S by blanks or commas)')
      write(*,*) 'hDD MM SS.sssss  Latitude :        (h default = N )'
c
   11 format(50a1)
c
   22 hem='N'
      read(*,11) buff
      call trim (buff,lgh,hem)
      call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlat1 = 0
      else
        irlat1 = 1
        write(*,*) ' Invalid Latitude ... Please re-enter it '
        write(*,*) '  '
        if( nowebb )then
          goto 22
        endif
      endif
c
      ald1 = dd
      alm1 = dm
      sl1  = ds
c
      if( hem.eq.'N' ) then
        lat1sn = +1
      else
        lat1sn = -1
      endif
c
      write(*,*) 'hDDD MM SS.sssss Longitude :       (h default = W )'
c
   23 hem='W'
      read(*,11) buff
      call trim (buff,lgh,hem)
      call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlon1 = 0
      else
        irlon1 = 1 
        write(*,*) ' Invalid Longitude ... Please re-enter it '
        write(*,*) '  '
        if( nowebb )then
          goto 23
        endif
      endif
c
      alod1 = dd
      alom1 = dm
      slo1  = ds
c
      if( hem.eq.'E' ) then
        lon1sn = +1
      else
        lon1sn = -1
      endif
c
      call getrad(ald1, alm1, sl1, lat1sn, glat1)
      call getrad(alod1,alom1,slo1,lon1sn, glon1)
c
   20 write(*,*) 'DDD MM SS.sss    Forward Azimuth :     (from north)'
c
   24 hem='A'
      read(*,11) buff
      call trim (buff,lgh,hem)
      call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        iazsn  = 1
        irazi1 = 0
      else
        irazi1 = 1
        write(*,*) ' Invalid Azimuth  ... Please re-enter it '
        write(*,*) '  '
        if( nowebb )then
          goto 24
        endif
      endif
c
      azd1 = dd
      azm1 = dm
      saz1 = ds
c
      call getrad(azd1,azm1,saz1,iazsn,faz)
c
      write(*,*) 'DDDDDD.dddd      Ellipsoidal Distance : (in meters)'
c
   25 ss = 0.0d0
      read(*,*) ss
      ss = dabs(ss)
c
      if( ss.lt.dd_max )then
        edist  = ss
        irdst1 = 0
      else
        irdst1 = 1
        write(*,*) ' Invalid Distance ... Please re-enter it '
        write(*,*) '  '
        if( nowebb )then
          goto 25
        endif
        edist  = 0.001d0
      endif
c
      call dirct1 (glat1,glon1,glat2,glon2,faz,baz,edist)
c
c     set the azimuth seconds tolerance
c
      tol = 0.00005d0
c
      call todmsp(faz,iaz1,maz1,saz1,isign1)
      if(isign1.lt.0) then
        iaz1=359-iaz1
        maz1=59-maz1
        saz1=60.d0-saz1
      endif
      call fixdms ( iaz1, maz1, saz1, tol )
c
      call todmsp(baz,iaz2,maz2,saz2,isign2)
      if(isign2.lt.0) then
        iaz2=359-iaz2
        maz2=59-maz2
        saz2=60.d0-saz2
      endif
      call fixdms ( iaz2, maz2, saz2, tol )
c
      call todmsp(glat1, ld1,  lm1,  sl1,  lat1sn)
      call todmsp(glon1, lod1, lom1, slo1, lon1sn)
      call todmsp(glat2, ld2,  lm2,  sl2,  lat2sn)
      call todmsp(glon2, lod2, lom2, slo2, lon2sn)
c
      call hem_ns ( lat1sn, d_ns1 )
      call hem_ew ( lon1sn, d_ew1 )
      call hem_ns ( lat2sn, d_ns2 )
      call hem_ew ( lon2sn, d_ew2 )
c
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
      write(*,49) elips
   49 format('  Ellipsoid : ',a30)
      finv=1.d0/f
      b=a*(1.d0-f)
      write(*,50) a,b,finv
   50 format('  Equatorial axis,    a   = ',f15.4,/,
     *       '  Polar axis,         b   = ',f15.4,/,
     *       '  Inverse flattening, 1/f = ',f16.11)
c
   18 format('    LAT = ',i3,1x,i2,1x,f8.5,1x,a6)
   19 format('    LON = ',i3,1x,i2,1x,f8.5,1x,a6)
c
      write(*,*) '  '
      write(*,*) '  First  Station : '
      write(*,*) '  ---------------- '
      write(*,18) ld1, lm1, sl1, d_ns1       
      write(*,19) lod1,lom1,slo1,d_ew1
c
      write(*,*) '  '
      write(*,*) '  Second Station : '
      write(*,*) '  ---------------- '
      write(*,18) ld2, lm2, sl2, d_ns2
      write(*,19) lod2,lom2,slo2,d_ew2
c
c
   32 format('  Ellipsoidal distance     S = ',f14.4,1x,a1)
   34 format('  Forward azimuth        FAZ = ',i3,1x,i2,1x,f7.4,
     1       ' From North')
   35 format('  Back azimuth           BAZ = ',i3,1x,i2,1x,f7.4,
     1       ' From North')
c
      write(*,*) '  '
      write(*,34) iaz1,maz1,saz1
      write(*,35) iaz2,maz2,saz2
      write(*,32) edist,dmt
      write(*,*) '  '
      write(*,*) '  Do you want to save this output into a file (y/n)?'
      read(*,10) answer
c
      if( answer.eq.'Y'.or.answer.eq.'y' )then
   39   write(*,*) '  Enter the output filename : '
        read(*,40) filout
   40   format(a30)
        open(3,file=filout,status='new',err=99)
        goto 98
c
   99   write(*,*) '  File already exists, try again.'
        go to 39
c
   98   continue
        write(3,*) '  '
        write(3,49) elips
        write(3,50) a,b,finv
        write(*,*) '  Enter the First Station name : '
        read(*,40) name1
        write(*,*) '  Enter the Second Station name : '
        read(*,40) name2
c
   41   format('  First  Station : ',a30)
   42   format('  Second Station : ',a30)
   84   format('  Error:  First  Station ... Invalid Latitude  ')
   85   format('  Error:  First  Station ... Invalid Longitude ')
   86   format('  Error:  Forward Azimuth .. Invalid Entry     ')
   87   format('  Error:  Ellipsoid Distance .. Invalid Entry  ')
   88   format(1x,65(1h*))  
   91   format('          DD(0-89) MM(0-59) SS(0-59.999...)    ')
   92   format('          DDD(0-359) MM(0-59) SS(0-59.999...)  ')
   93   format('          Geodetic distance is too long        ')    
c
        write(3,*) '  '
        write(3,41) name1
        write(3,*) '  ---------------- '
c
        if( irlat1.eq.0 )then
          write(3,18) ld1, lm1, sl1, d_ns1
        else
          write(3,88)
          write(3,84)
          write(3,91)                       
          write(3,88)
        endif
c 
        if( irlon1.eq.0 )then      
          write(3,19) lod1,lom1,slo1,d_ew1
        else
          write(3,88)
          write(3,85)
          write(3,92)                         
          write(3,88)
        endif
c
        write(3,*) '  '
        write(3,42) name2
        write(3,*) '  ---------------- '
        write(3,18) ld2, lm2, sl2, d_ns2
        write(3,19) lod2,lom2,slo2,d_ew2
c
        write(3,*) '  '
        if( irazi1.eq.0 )then
          write(3,34) iaz1,maz1,saz1
        else
          write(3,88)
          write(3,86)
          write(3,92)                
          write(3,88)
        endif
c
        write(3,35) iaz2,maz2,saz2
c
        if( irdst1.eq.0 )then
          write(3,32) edist,dmt
        else
          write(3,88)
          write(3,87)
          write(3,93)            
          write(3,88)
        endif
c
        write(3,*) '  '
      endif
c
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  1) Another forward, different ellipsoid.'
      write(*,*) '  2) Same ellipsoid, two new stations.'
      write(*,*) '  3) Same ellipsoid, same First Station.'
      write(*,*) '  4) Let the Second Station be the First Station.'
      write(*,*) '  5) Done.'
      write(*,*) '  '
      write(*,*) '  Enter choice : '
      read(*,10) answer
c
      if(    answer.eq.'1')then
        goto 1
      elseif(answer.eq.'2')then
        goto 15
      elseif(answer.eq.'3')then
        goto 20
      elseif(answer.eq.'4')then
        glat1 = glat2
        glon1 = glon2
        goto 20
      else
        write(*,*) '  Thats all, folks!'
      endif
c
c     stop
      end

      subroutine bufdms (buff,lgh,hem,dd,dm,ds,ierror)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
      logical     done,flag
c
      character*1 buff(*),abuf(21)
      character*1 ch
      character*1 hem
      integer*4   ll,lgh
      integer*4   i4,id,im,is,icond,ierror
      real*8      x(5)
c
c     set the "error flag" 
c
      ierror = 0
      icond  = 0
c
c     set defaults for dd,dm,ds
c
      dd = 0.0d0
      dm = 0.0d0
      ds = 0.0d0
c
c     set default limits for "hem" flag
c
      if(     hem.eq.'N' .or. hem.eq.'S' )then
        ddmax = 90.0d0
      elseif( hem.eq.'E' .or. hem.eq.'W' )then
        ddmax = 360.0d0
      elseif( hem.eq.'A' )then
        ddmax = 360.0d0
      elseif( hem.eq.'Z' )then
        ddmax = 180.0d0
      elseif( hem.eq.'*' )then
        ddmax  = 0.0d0
        ierror = 1
      else
        ddmax = 360.0d0
      endif
c
      do 1 i=1,5
        x(i) = 0.0d0
    1 continue
c
      icolon = 0
      ipoint = 0
      icount = 0
      flag   = .true.
      jlgh   = lgh
c
      do 2 i=1,jlgh
        if( buff(i).eq.':' )then
          icolon = icolon+1
        endif
        if( buff(i).eq.'.' )then
          ipoint = ipoint+1
          flag   = .false.
        endif
        if( flag )then
          icount = icount+1
        endif
    2 continue
c
      if( ipoint.eq.1 .and. icolon.eq.0 )then
c
c       load temp buffer
c
        do 3 i=1,jlgh
          abuf(i) = buff(i)
    3   continue
        abuf(jlgh+1) = '$'
        ll = jlgh
c
        call gvalr8 (abuf,ll,r8,icond)
c
        if( icount.ge.5 )then
c
c         value is a packed decimal of ==>  DDMMSS.sssss       
c
          ss = r8/10000.0d0
          id = idint( ss )
c
          r8 = r8-10000.0d0*dble(float(id))
          ss = r8/100.0d0
          im = idint( ss )
c
          r8 = r8-100.0d0*dble(float(im))
        else
c
c         value is a decimal of ==>  .xx   X.xxx   X.  
c
          id = idint( r8 )
          r8 = (r8-id)*60.0d0
          im = idint( r8 )
          r8 = (r8-im)*60.0d0
        endif
c
c       account for rounding error
c
        is = idnint( r8*1.0d5 )
        if( is.ge.6000000 )then
           r8 = 0.0d0
           im = im+1
        endif
c
        if( im.ge.60 )then
          im = 0
          id = id+1
        endif
c
        dd = dble( float( id ) )
        dm = dble( float( im ) )
        ds = r8
      else
c
c       buff() value is a d,m,s of ==>  NN:NN:XX.xxx    
c
        k    = 0
        next = 1
        done = .false.
        ie   = jlgh
c
        do 100 j=1,5
          ib = next
          do 90 i=ib,ie
            ch   = buff(i)
            last = i
            if( i.eq.jlgh .or. ch.eq.':' )then
              if( i.eq.jlgh )then
                done = .true.
              endif
              if( ch.eq.':' )then
                last = i-1
              endif
              goto 91
            endif
   90     continue
          goto 98
c
   91     ipoint = 0
          ik     = 0
          do 92 i=next,last
            ik = ik+1
            ch = buff(i)
            if( ch.eq.'.' )then
              ipoint = ipoint+1
            endif
            abuf(ik) = buff(i) 
   92     continue
          abuf(ik+1) = '$' 
c
          ll = ik
          if( ipoint.eq.0 )then
            call gvali4 (abuf,ll,i4,icond)
            r8 = dble(float( i4 )) 
          else
            call gvalr8 (abuf,ll,r8,icond)
          endif
c
          k    = k+1
          x(k) = r8
c
   98     if( done )then
            goto 101
          endif
c
          next = last
   99     next = next+1     
          if( buff(next).eq.':' )then
            goto 99
          endif
  100   continue
c
c       load dd,dm,ds
c
  101   if( k.ge.1 )then
          dd = x(1)
        endif
c
        if( k.ge.2 )then
          dm = x(2)
        endif
c
        if( k.ge.3 )then
          ds = x(3)
        endif
      endif
c
      if( dd.gt.ddmax  .or.
     1    dm.ge.60.0d0 .or.
     1    ds.ge.60.0d0 )then
        ierror = 1
        dd = 0.0d0
        dm = 0.0d0
        ds = 0.0d0
      endif
c
      if( icond.ne.0 )then
        ierror = 1
      endif
c
      return
      end

      SUBROUTINE DIRCT1(GLAT1,GLON1,GLAT2,GLON2,FAZ,BAZ,S)
C
C *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY
C *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
C
C *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
C *** F IS THE FLATTENING OF THE REFERENCE ELLIPSOID
C *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
C *** AZIMUTHS IN RADIANS CLOCKWISE FROM NORTH
C *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A
C
C *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75
C *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CONST/PI,RAD
      COMMON/ELIPSOID/A,F
      DATA EPS/0.5D-13/
      R=1.-F
      TU=R*DSIN(GLAT1)/DCOS(GLAT1)
      SF=DSIN(FAZ)
      CF=DCOS(FAZ)
      BAZ=0.
      IF(CF.NE.0.) BAZ=DATAN2(TU,CF)*2.
      CU=1./DSQRT(TU*TU+1.)
      SU=TU*CU
      SA=CU*SF
      C2A=-SA*SA+1.
      X=DSQRT((1./R/R-1.)*C2A+1.)+1.
      X=(X-2.)/X
      C=1.-X
      C=(X*X/4.+1)/C
      D=(0.375D0*X*X-1.)*X
      TU=S/R/A/C
      Y=TU
  100 SY=DSIN(Y)
      CY=DCOS(Y)
      CZ=DCOS(BAZ+Y)
      E=CZ*CZ*2.-1.
      C=Y
      X=E*CY
      Y=E+E-1.
      Y=(((SY*SY*4.-3.)*Y*CZ*D/6.+X)*D/4.-CZ)*SY*D+TU
      IF(DABS(Y-C).GT.EPS)GO TO 100
      BAZ=CU*CY*CF-SU*SY
      C=R*DSQRT(SA*SA+BAZ*BAZ)
      D=SU*CY+CU*SY*CF
      GLAT2=DATAN2(D,C)
      C=CU*CY-SU*SY*CF
      X=DATAN2(SY*SF,C)
      C=((-3.*C2A+4.)*F+4.)*C2A*F/16.
      D=((E*CY*C+CZ)*SY*C+Y)*SA
      GLON2=GLON1+X-(1.-C)*D*F
      BAZ=DATAN2(SA,BAZ)+PI
      RETURN
      END

      subroutine elipss (elips)
      implicit double precision(a-h,o-z)
      character*1  answer
      character*30 elips
      common/elipsoid/a,f
      write(*,*) '  Other Ellipsoids.'
      write(*,*) '  -----------------'
      write(*,*) '  '
      write(*,*) '  A) Airy 1858'
      write(*,*) '  B) Airy Modified'
      write(*,*) '  C) Australian National'
      write(*,*) '  D) Bessel 1841'
      write(*,*) '  E) Clarke 1880'
      write(*,*) '  F) Everest 1830'
      write(*,*) '  G) Everest Modified'
      write(*,*) '  H) Fisher 1960'
      write(*,*) '  I) Fisher 1968'
      write(*,*) '  J) Hough 1956'
      write(*,*) '  K) International (Hayford)'
      write(*,*) '  L) Krassovsky 1938'
      write(*,*) '  M) NWL-9D (WGS 66)'
      write(*,*) '  N) South American 1969'
      write(*,*) '  O) Soviet Geod. System 1985'
      write(*,*) '  P) WGS 72'
      write(*,*) '  Q-Z) User defined.'
      write(*,*) '  '
      write(*,*) '  Enter choice : '
      read(*,10) answer
   10 format(a1)
c
      if(answer.eq.'A'.or.answer.eq.'a') then
        a=6377563.396d0
        f=1.d0/299.3249646d0
        elips='Airy 1858'
      elseif(answer.eq.'B'.or.answer.eq.'b') then
        a=6377340.189d0
        f=1.d0/299.3249646d0
        elips='Airy Modified'
      elseif(answer.eq.'C'.or.answer.eq.'c') then
        a=6378160.d0
        f=1.d0/298.25d0
        elips='Australian National'
      elseif(answer.eq.'D'.or.answer.eq.'d') then
        a=6377397.155d0
        f=1.d0/299.1528128d0
        elips='Bessel 1841'
      elseif(answer.eq.'E'.or.answer.eq.'e') then
        a=6378249.145d0
        f=1.d0/293.465d0
        elips='Clarke 1880'
      elseif(answer.eq.'F'.or.answer.eq.'f') then
        a=6377276.345d0
        f=1.d0/300.8017d0
        elips='Everest 1830'
      elseif(answer.eq.'G'.or.answer.eq.'g') then
        a=6377304.063d0
        f=1.d0/300.8017d0
        elips='Everest Modified'
      elseif(answer.eq.'H'.or.answer.eq.'h') then
        a=6378166.d0
        f=1.d0/298.3d0
        elips='Fisher 1960'
      elseif(answer.eq.'I'.or.answer.eq.'i') then
        a=6378150.d0
        f=1.d0/298.3d0
        elips='Fisher 1968'
      elseif(answer.eq.'J'.or.answer.eq.'j') then
        a=6378270.d0
        f=1.d0/297.d0
        elips='Hough 1956'
      elseif(answer.eq.'K'.or.answer.eq.'k') then
        a=6378388.d0
        f=1.d0/297.d0
        elips='International (Hayford)'
      elseif(answer.eq.'L'.or.answer.eq.'l') then
        a=6378245.d0
        f=1.d0/298.3d0
        elips='Krassovsky 1938'
      elseif(answer.eq.'M'.or.answer.eq.'m') then
        a=6378145.d0
        f=1.d0/298.25d0
        elips='NWL-9D  (WGS 66)'
      elseif(answer.eq.'N'.or.answer.eq.'n') then
        a=6378160.d0
        f=1.d0/298.25d0
        elips='South American 1969'
      elseif(answer.eq.'O'.or.answer.eq.'o') then
        a=6378136.d0
        f=1.d0/298.257d0
        elips='Soviet Geod. System 1985'
      elseif(answer.eq.'P'.or.answer.eq.'p') then
        a=6378135.d0
        f=1.d0/298.26d0
        elips='WGS 72'
      else
        elips='User defined.'
c
        write(*,*) '  Enter Equatorial axis,   a : '
        read(*,*) a
        a = dabs(a)
c
        write(*,*) '  Enter either Polar axis,   b  or '
        write(*,*) '  Reciprocal flattening,   1/f : '
        read(*,*) ss
        ss = dabs(ss)
c
        f = 0.0d0
        if( 200.0d0.le.ss .and. ss.le.310.0d0 )then
          f = 1.d0/ss   
        elseif( 6000000.0d0.lt.ss .and. ss.lt.a )then
          f = (a-ss)/a
        else
          elips = 'Error: default GRS80 used.'
          a     = 6378137.0d0
          f     = 1.0d0/298.25722210088d0
        endif
      endif
c
      return
      end

      subroutine fixdms (ideg, min, sec, tol )
c
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
c     test for seconds near 60.0-tol
c
      if( sec.ge.( 60.0d0-tol ) )then
        sec  = 0.0d0
        min  = min+1
      endif
c
c     test for minutes near 60
c
      if( min.ge.60 )then
        min  = 0
        ideg = ideg+1
      endif 
c
c     test for degrees near 360
c
      if( ideg.ge.360 )then
        ideg = 0
      endif 
c
      return
      end 

      subroutine hem_ns ( lat_sn, hem )
      implicit integer (i-n)
      character*6  hem
c
      if( lat_sn.eq.1 ) then
        hem = 'North '
      else
        hem = 'South '
      endif
c
      return
      end
      subroutine hem_ew ( lon_sn, hem )
      implicit integer (i-n)
      character*6  hem
c
      if( lon_sn.eq.1 ) then
        hem = 'East  '
      else
        hem = 'West  '
      endif
c
      return
      end
      subroutine getrad(d,m,sec,isign,val)
 
*** comvert deg, min, sec to radians
 
      implicit double precision(a-h,j-z)
      common/const/pi,rad
 
      val=(d+m/60.d0+sec/3600.d0)/rad
      val=dble(isign)*val
 
      return
      end

CB::GPNARC
C
      SUBROUTINE GPNARC (AMAX,FLAT,ESQ,PI,P1,P2,ARC)
C
C********1*********2*********3*********4*********5*********6*********7*
C
C NAME:        GPNARC
C VERSION:     200005.26
C WRITTEN BY:  ROBERT (Sid) SAFFORD
C PURPOSE:     SUBROUTINE TO COMPUTE THE LENGTH OF A MERIDIONAL ARC 
C              BETWEEN TWO LATITUDES
C
C INPUT PARAMETERS:
C -----------------
C AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
C FLAT         FLATTENING (0.0033528 ... )
C ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
C PI           3.14159...
C P1           LAT STATION 1
C P2           LAT STATION 2
C
C OUTPUT PARAMETERS:
C ------------------
C ARC          GEODETIC DISTANCE 
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C GLOBAL VARIABLES AND CONSTANTS:
C -------------------------------
C
C    MODULE CALLED BY:    GENERAL 
C
C    THIS MODULE CALLS:   
C       LLIBFORE/ OPEN,   CLOSE,  READ,   WRITE,  INQUIRE
C                 DABS,   DBLE,   FLOAT,  IABS,   CHAR,   ICHAR
C
C    INCLUDE FILES USED:
C    COMMON BLOCKS USED:  
C
C    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
C                MS-DOS Operating System
C    COMMENTS:
C********1*********2*********3*********4*********5*********6*********7*
C::MODIFICATION HISTORY
C::197507.05, RWS, VER 00 TENCOL RELEASED FOR FIELD USE
C::198311.20, RWS, VER 01 MTEN   RELEASED TO FIELD
C::198411.26, RWS, VER 07 MTEN2  RELEASED TO FIELD
C::1985xx.xx, RWS, CODE   CREATED               
C::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
C::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
C::198512.18, RWS, CODE   MODIFIED FOR MTEN3
C::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
C::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
C::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
C::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED             
C::200012.31, RWS, VER 23 MTEN5 RELEASED                                 
C********1*********2*********3*********4*********5*********6*********7*
CE::GPNARC
C ---------------------------
C     M T E N  (VERSION 3)
C     M T E N  (VERSION 5.23)
C ---------------------------
C 
      IMPLICIT REAL*8 (A-H,O-Z)
C
      LOGICAL  FLAG
C
      DATA TT/5.0D-15/
C
C     CHECK FOR A 90 DEGREE LOOKUP
C
      FLAG = .FALSE.
C
      S1 = DABS(P1)
      S2 = DABS(P2)
C
      IF( (PI/2.0D0-TT).LT.S2 .AND. S2.LT.(PI/2.0D0+TT) )THEN
        FLAG = .TRUE.
      ENDIF
C
      IF( S1.GT.TT )THEN
        FLAG = .FALSE.
      ENDIF
C
      DA = (P2-P1)
      S1 = 0.0D0
      S2 = 0.0D0
C
C     COMPUTE THE LENGTH OF A MERIDIONAL ARC BETWEEN TWO LATITUDES
C
      E2 = ESQ
      E4 = E2*E2
      E6 = E4*E2
      E8 = E6*E2
      EX = E8*E2
C
      T1 = E2*(003.0D0/4.0D0)
      T2 = E4*(015.0D0/64.0D0)
      T3 = E6*(035.0D0/512.0D0)
      T4 = E8*(315.0D0/16384.0D0)
      T5 = EX*(693.0D0/131072.0D0)
C
      A  = 1.0D0+T1+3.0D0*T2+10.0D0*T3+35.0D0*T4+126.0D0*T5
C
      IF( FLAG )THEN
        GOTO 1
      ENDIF
C
      B  = T1+4.0D0*T2+15.0D0*T3+56.0D0*T4+210.0D0*T5
      C  = T2+06.0D0*T3+28.0D0*T4+120.0D0*T5
      D  = T3+08.0D0*T4+045.0D0*T5
      E  = T4+010.0D0*T5
      F  = T5
C
      DB = DSIN(P2*2.0D0)-DSIN(P1*2.0D0)
      DC = DSIN(P2*4.0D0)-DSIN(P1*4.0D0)
      DD = DSIN(P2*6.0D0)-DSIN(P1*6.0D0)
      DE = DSIN(P2*8.0D0)-DSIN(P1*8.0D0)
      DF = DSIN(P2*10.0D0)-DSIN(P1*10.0D0)
C
C     COMPUTE THE S2 PART OF THE SERIES EXPANSION
C
      S2 = -DB*B/2.0D0+DC*C/4.0D0-DD*D/6.0D0+DE*E/8.0D0-DF*F/10.0D0
C
C     COMPUTE THE S1 PART OF THE SERIES EXPANSION
C
    1 S1 = DA*A
C
C     COMPUTE THE ARC LENGTH
C
      ARC = AMAX*(1.0D0-ESQ)*(S1+S2)
C
      RETURN
      END

      subroutine gvali4 (buff,ll,vali4,icond)
      implicit     integer (i-n)
c
      logical      plus,sign,done,error
      character*1  buff(*)
      character*1  ch
c
c     integer*2    i
c     integer*2    l1
c
      integer*4    ich,icond
      integer*4    ll    
      integer*4    vali4
c
      l1    = ll
      vali4 = 0
      icond = 0
      plus  = .true.
      sign  = .false.
      done  = .false.
      error = .false.
c
      i = 0
   10 i = i+1
      if( i.gt.l1 .or. done )then
        go to 1000
      else
        ch  = buff(i)
        ich = ichar( buff(i) )
      endif
c
      if(     ch.eq.'+' )then
c
c       enter on plus sign
c
        if( sign )then
          goto 150
        else 
          sign = .true.
          goto 10
        endif
      elseif( ch.eq.'-' )then
c
c       enter on minus sign
c
        if( sign )then
          goto 150
        else
          sign = .true.
          plus = .false.
          goto 10
        endif
      elseif( ch.ge.'0' .and. ch.le.'9' )then
        goto 100
      elseif( ch.eq.' ' )then
c
c       enter on space -- ignore leading spaces
c
        if( .not.sign )then
          goto 10
        else
          buff(i) = '0'
          ich = 48
          goto 100
        endif
      elseif( ch.eq.':' )then
c
c       enter on colon -- ignore 
c
        if( .not.sign )then
          goto 10
        else
          goto 1000
        endif
      elseif( ch.eq.'$' )then
c
c       enter on dollar "$"      
c
        done = .true.
        goto 10
      else
c
c       something wrong
c
        goto 150
      endif
c
c     enter on numeric
c
  100 vali4 = 10*vali4+(ich-48)
      sign  = .true.
      goto 10
c
c     treat illegal character
c
  150 buff(i) = '0'
      vali4 = 0
      icond = 1
c
 1000 if( .not.plus )then
        vali4 = -vali4
      endif
c
      return
      end
      subroutine gvalr8 (buff,ll,valr8,icond)
      implicit     integer (i-n)
c
      logical      plus,sign,dpoint,done
c
      character*1  buff(*)
      character*1  ch
c
c     integer*2    i, ip
c     integer*2    l1
c     integer*2    nn, num, n48
c
      integer*4    ich,icond
      integer*4    ll
c
      real*8       ten
      real*8       valr8
      real*8       zero
c
      data zero,ten/0.0d0,10.0d0/
c
      n48     =  48
      l1      =  ll
      icond   =   0
      valr8   =  zero  
      plus    = .true.
      sign    = .false.
      dpoint  = .false.
      done    = .false.
c
c     start loop thru buffer
c
      i = 0
   10 i = i+1
      if( i.gt.l1 .or. done )then
        go to 1000
      else 
        ch  = buff(i)
        nn  = ichar( ch )
        ich = nn
      endif 
c
      if(     ch.eq.'+' )then
c
c       enter on plus sign
c
        if( sign )then
          goto 150
        else
          sign = .true.
          goto 10
        endif
      elseif( ch.eq.'-' )then
c
c       enter on minus sign
c
        if( sign )then
          goto 150
        else
          sign = .true.
          plus = .false.
          goto 10
        endif
      elseif( ch.eq.'.' )then
c
c       enter on decimal point
c
        ip     = 0
        sign   = .true.
        dpoint = .true.
        goto 10
      elseif( ch.ge.'0' .and. ch.le.'9' )then
        goto 100
      elseif( ch.eq.' ' )then
c
c       enter on space
c
        if( .not.sign )then
          goto 10
        else
          buff(i) = '0'
          ich = 48
          goto 100
        endif
      elseif( ch.eq.':' .or. ch.eq.'$' )then
c
c       enter on colon or "$" sign
c
        done = .true.
        goto 10
      else
c
c       something wrong
c
        goto 150
      endif
c
c     enter on numeric
c
  100 sign = .true.
      if( dpoint )then
        ip = ip+1
      endif
c
      num   = ich
      valr8 = ten*valr8+dble(float( num-n48 ))
      goto 10
c
c     treat illegal character
c
  150 buff(i) = '0'
      valr8   =  0.0d0
      icond   =  1
c
 1000 if( dpoint )then
        valr8 =  valr8/(ten**ip)
      endif
c
      if( .not.plus )then
        valr8 = -valr8
      endif
c
      return
      end

      subroutine todmsp(val,id,im,s,isign)
 
*** convert position radians to deg,min,sec
*** range is [-pi to +pi]
 
      implicit double precision(a-h,o-z)
      common/const/pi,rad
 
    1 if(val.gt.pi) then
        val=val-pi-pi
        go to 1
      endif
 
    2 if(val.lt.-pi) then
        val=val+pi+pi
        go to 2
      endif
 
      if(val.lt.0.d0) then
        isign=-1
      else
        isign=+1
      endif
 
      s=dabs(val*rad)
      id=idint(s)
      s=(s-id)*60.d0
      im=idint(s)
      s=(s-im)*60.d0
 
*** account for rounding error
 
      is=idnint(s*1.d5)
      if(is.ge.6000000) then
        s=0.d0
        im=im+1
      endif
      if(im.ge.60) then
        im=0
        id=id+1
      endif
 
      return
      end

      subroutine trim (buff,lgh,hem)
      implicit integer (i-n)
      character*1 ch,hem
      character*1 buff(*)
      integer*4   lgh
c
      ibeg = 1
      do 10 i=1,50
        if( buff(i).ne.' ' )then
          goto 11
        endif
        ibeg = ibeg+1
   10 continue
   11 continue
      if( ibeg.ge.50 )then
        ibeg = 1
        buff(ibeg) = '0'
      endif
c
      iend = 50
      do 20 i=1,50
        j = 51-i
        if( buff(j).eq.' ' )then
          iend = iend-1
        else
          goto 21
        endif
   20 continue
   21 continue
c
      ch = buff(ibeg)
      if( hem.eq.'N' )then
        if( ch.eq.'N' .or. ch.eq.'n' .or. ch.eq.'+' )then
          hem = 'N'
          ibeg = ibeg+1
        endif
        if( ch.eq.'S' .or. ch.eq.'s' .or. ch.eq.'-' )then
          hem = 'S'
          ibeg = ibeg+1
        endif
c
c       check for wrong hemisphere entry
c
        if( ch.eq.'E' .or. ch.eq.'e' )then
          hem = '*'
          ibeg = ibeg+1
        endif
        if( ch.eq.'W' .or. ch.eq.'w' )then
          hem = '*'
          ibeg = ibeg+1
        endif
      elseif( hem.eq.'W' )then
        if( ch.eq.'E' .or. ch.eq.'e' .or. ch.eq.'+' )then
          hem = 'E'
          ibeg = ibeg+1
        endif
        if( ch.eq.'W' .or. ch.eq.'w' .or. ch.eq.'-' )then
          hem = 'W'
          ibeg = ibeg+1
        endif
c
c       check for wrong hemisphere entry
c
        if( ch.eq.'N' .or. ch.eq.'n' )then
          hem = '*'
          ibeg = ibeg+1
        endif
        if( ch.eq.'S' .or. ch.eq.'s' )then
          hem = '*'
          ibeg = ibeg+1
        endif
      elseif( hem.eq.'A' )then
        if( .not.('0'.le.ch .and. ch.le.'9') )then
          hem = '*'
          ibeg = ibeg+1
        endif
      else
c        do nothing
      endif
c
c
      do 30 i=ibeg,iend
        ch = buff(i)
c
        if(     ch.eq.':' .or. ch.eq.'.' )then
          goto 30
        elseif( ch.eq.' ' .or. ch.eq.',' )then
          buff(i) = ':'
        elseif( '0'.le.ch .and. ch.le.'9' )then
          goto 30      
        else
          buff(i) = ':'
        endif
c
   30 continue
c
c     left justify buff() array to its first character position
c     also check for a ":" char in the starting position,
c     if found!!  skip it
c
      j  = 0
      ib = ibeg
      ie = iend
c
      do 40 i=ib,ie
        if( i.eq.ibeg .and. buff(i).eq.':' )then
c
c         move the 1st position pointer to the next char &
c         do not put ":" char in buff(j) array where j=1    
c
          ibeg = ibeg+1
          goto 40
        endif
        j = j+1
        buff(j) = buff(i)
   40 continue
c
c
      lgh = iend-ibeg+1
      j   = lgh+1
      buff(j) = '$'
c
c     clean-up the rest of the buff() array
c
      do 50 i=j+1,50   
        buff(i) = ' '    
   50 continue
c
c     save a maximum of 20 characters
c
      if( lgh.gt.20 )then
        lgh = 20
        j   = lgh+1
        buff(j) = '$'
      endif
c
      return
      end

      real*8 function zvalue ( ss, tol )
c
      implicit double precision (a-z)
c
c     to check for dx,dy,dz or 
c     dn,de,du values below tol.
c
      if( dabs(ss).lt.tol )then
        ss = 0.0d0
      endif
c
      zvalue = ss
c
      return
      end


cb::inverse
c
      program inverse
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:      inverse
c version:   201105.xx
c author:    stephen j. frakes
c last mod:  dr. dennis milbert
c purpose:   to compute a geodetic inverse  
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
c b                semiminor axis polar (in meters)
c baz              azimuth back (in radians)
c buff             input char buffer array
c dd,dm,ds         temporary values for degrees, minutes, seconds
c dlon             temporary value for difference in longitude (radians)   
c dmt,d_mt         char constants for meter units         
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
c ss               temporary variable     
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
c glat1,glon1      station one       - (lat & lon in radians )
c p1,e1            standpoint one    - (lat & lon in radians )
c
c name2            name of station two
c ld2,lm2,sl2      latitude  sta two - degrees,minutes,seconds
c ald2,alm2,sl2    latitude  sta two - degrees,minutes,seconds
c lat2sn           latitude  sta two - sign (+/- 1)
c d_ns2            latitude  sta two - char ('N','S')
c lod2,lom2,slo2   longitude sta two - degrees,minutes,seconds
c alod2,alom2,slo2 longitude sta one - degrees,minutes,seconds
c lon2sn           longitude sta two - sign (+/- 1)
c d_ew2            longitude sta two - char ('E','W')
c iaz2,maz2,saz2   back azimuth      - degrees,minutes,seconds
c isign2           back azimuth      - flag  (+/- 1)
c glat2,glon2      station two       - (lat & lon in radians )
c p2,e2            forepoint two     - (lat & lon in radians )
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
c    this module calls:      elipss, getrad, inver1, todmsp
c    gethem, trim,   bufdms, gvalr8, gvali4, fixdms, gpnhri ***********
c    gethem, trim,   bufdms, gvalr8, gvali4, fixdms, invers <----------
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
c::200207.15, rws, modified i/o & standardized program documentation
c::                added subs trim, bufdms, gethem, gvali4, gvalr8      
c::200207.23, rws, replaced sub inver1 with gpnarc, gpnloa, gpnhri
c::200208.15, rws, fixed an error in bufdms
c::              - renamed ellips to elipss "common error" with dirct2
c::              - added FAZ & BAZ to printed output      
c::200208.19, rws, added more error flags for web interface code
c::              - added logical nowebb                             
c::200208.xx, sjf, program version number 2.0                   
c::201105.xx, dgm, program version number 3.0
c::              - replaced sub gpnarc, gpnloa, gpnhri with invers
c::              - tested for valid antipodal solutions (+/- 0.1 mm)
c::              - tested for polar solutions (+/- 0.1 mm)
c::              - needs improvement for long-line/antipodal boundary
c********1*********2*********3*********4*********5*********6*********7**
ce::inverse
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
c                  1 = ms_dos or unix version
c
      parameter   ( ms_unix = 0 )
c
      pi   = 4.d0*datan(1.d0)
      rad  = 180.d0/pi
      dmt  = 'm'
      d_mt = 'Meters'
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
c 
    5 write(*,*) '  Program Inverse  -  Version 3.0 '
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
        f=1.d0/298.257222100882711243162836600094d0
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
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
c
   15 write(*,*) '  Enter First Station '
      write(*,16)
   16 format(18x,'(Separate D,M,S by blanks or commas)')
      write(*,*) 'hDD MM SS.sssss  Latitude :        (h default = N )'
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
      write(*,*) '  '
      write(*,*) '  '
c
   20 write(*,*) '  Enter Second Station '
      write(*,16)
      write(*,*) 'hDD MM SS.sssss  Latitude :        (h default = N )'
c
   24 hem='N'
      read(*,11) buff
      call trim (buff,lgh,hem)
      call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlat2 = 0
      else
        irlat2 = 1
        write(*,*) ' Invalid Latitude ... Please re-enter it '
        write(*,*) '  '
        if( nowebb )then
          goto 24
        endif
      endif
c
      ald2 = dd
      alm2 = dm
      sl2  = ds
c
      if( hem.eq.'N' ) then
        lat2sn = +1
      else
        lat2sn = -1
      endif
c
      write(*,*) 'hDDD MM SS.sssss Longitude :       (h default = W )'
c
   25 hem='W'
      read(*,11) buff
      call trim (buff,lgh,hem)
      call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlon2 = 0
      else
        irlon2 = 1
        write(*,*) ' Invalid Longitude ... Please re-enter it '
        write(*,*) '  '
        if( nowebb )then
          goto 25
        endif
      endif
c
      alod2 = dd
      alom2 = dm
      slo2  = ds
c
      if( hem.eq.'E' )then
        lon2sn = +1
      else
        lon2sn = -1
      endif
c
      call getrad(ald2, alm2, sl2, lat2sn, glat2)
      call getrad(alod2,alom2,slo2,lon2sn, glon2)
c
      p1 = glat1
      e1 = glon1
      p2 = glat2
      e2 = glon2
c
      if( e1.lt.0.0d0 )then
        e1 = e1+2.0d0*pi
      endif
      if( e2.lt.0.0d0 )then
        e2 = e2+2.0d0*pi
      endif
c
c     compute the geodetic inverse
c
c ************************************************************
c *   replaced subroutine inver1 with gpnhri
c *  
c *   call inver1 (glat1,glon1,glat2,glon2,faz,baz,edist)
c *
c ************************************************************
c ************************************************************
c *   replaced subroutine inver1 with gpnhri with invers -- dgm
c *  
c *   call gpnhri (a,f,esq,pi,p1,e1,p2,e2,faz,baz,edist)
c *
c ************************************************************

      finv=1.d0/f
      call invers(a,finv,p1,e1,p2,e2,faz,baz,edist,
     *            numiter,sigmaout,dlambdaout,iflagout)

c     check for a non distance ... p1,e1 & p2,e2 equal zero ?
c
      if( edist.lt.0.00005d0 )then
        faz = 0.0d0
        baz = 0.0d0
      endif
c
c     set the tolerance (in seconds) for the azimuth conversion
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
        finv=1.d0/f
        b=a*(1.d0-f)
        write(3,50) a,b,finv
        write(*,*) '  Enter the First Station name : '
        read(*,40) name1
        write(*,*) '  Enter the Second Station name : '
        read(*,40) name2
c
   41   format('  First  Station : ',a30)
   42   format('  Second Station : ',a30)
   84   format('  Error: First  Station ... Invalid Latitude  ')
   85   format('  Error: First  Station ... Invalid Longitude ')
   86   format('  Error: Second Station ... Invalid Latitude  ')
   87   format('  Error: Second Station ... Invalid Longitude ')
   88   format(1x,65(1h*))
   91   format('         DD(0-89) MM(0-59) SS(0-59.999...)  ')
   92   format('         DDD(0-359) MM(0-59) SS(0-59.999...)  ')
c
        write(3,*) '  '
        write(3,41) name1
        write(3,*) '  ---------------- '

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
c
        if( irlat2.eq.0 )then
          write(3,18) ld2, lm2, sl2, d_ns2
        else
          write(3,88)
          write(3,86)
          write(3,91)                       
          write(3,88)
        endif
c
        if( irlon2.eq.0 )then
          write(3,19) lod2,lom2,slo2,d_ew2
        else
          write(3,88)
          write(3,87)
          write(3,92)                       
          write(3,88)
        endif
c
        write(3,*) '  '
        write(3,34) iaz1,maz1,saz1
        write(3,35) iaz2,maz2,saz2
        write(3,32) edist,dmt
        write(3,*) '  '
      endif
c
   80 write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  '
      write(*,*) '  1) Another inverse, different ellipsoid.'
      write(*,*) '  2) Same ellipsoid, two new stations.'
      write(*,*) '  3) Same ellipsoid, same First Station.'
      write(*,*) '  4) Done.'
      write(*,*) '  '
      write(*,*) '  Enter choice : '
      read(*,10) answer
c
      if(     answer.eq.'1' )then
        goto 1
      elseif( answer.eq.'2' )then
        goto 15
      elseif( answer.eq.'3' )then
        goto 20
      else
        write(*,*) '  Thats all, folks!'
      endif

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
        elips = 'User defined.'
c
        write(*,*) '  Enter Equatorial axis,   a : '
        read(*,*) a
        a  = dabs(a)
c
        write(*,*) '  Enter either Polar axis, b or '
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
c
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
*-----------------------------------------------------------------------
      subroutine invers(a,rf,b1,l1,b2,l2,faz,baz,s,it,sig,lam,kind)

*** inverse for long-line and antipodal cases.
*** latitudes may be 90 degrees exactly.
*** latitude positive north, longitude positive east, radians.
*** azimuth clockwise from north, radians.
*** original programmed by thaddeus vincenty, 1975, 1976
*** removed back side solution option, debugged, revised -- 2011may01 -- dgm
*** this version of code is interim -- antipodal boundary needs work

*** output (besides faz,baz, and s):
*** it,   iteration count
*** sig,  spherical distance on auxiliary sphere
*** lam,  longitude difference on auxiliary sphere
*** kind, solution flag:  kind=1, long-line;  kind=2, antipodal

      implicit double precision(a-h,o-z)
      double precision l1,l2,l,lam
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(tol=1.d-14, eps=1.d-15)

      boa =  1.d0-1.d0/rf
***** sinu1 = boa*dsin(b1)/dsqrt((boa*dsin(b1))**2+dcos(b1)**2)
***** cosu1 = dsqrt(-sinu1**2+1.d0)                               !*** roundoff
***** sinu2 = boa*dsin(b2)/dsqrt((boa*dsin(b2))**2+dcos(b2)**2)
***** cosu2 = dsqrt(-sinu2**2+1.d0)                               !*** roundoff

      beta1 = datan(boa*dtan(b1))             !*** better reduced latitude
      sinu1 = dsin(beta1)
      cosu1 = dcos(beta1)

      beta2 = datan(boa*dtan(b2))             !*** better reduced latitude
      sinu2 = dsin(beta2)
      cosu2 = dcos(beta2)

      l = l2-l1                               !*** longitude difference [-pi,pi]
      if(l.gt. pi) l=l-pi-pi
      if(l.lt.-pi) l=l+pi+pi
      prev = l
      test = l
      it =0
      kind=1
      lam=l                                              !*** v13   (rapp part II)

*** top of the long-line loop (kind=1)

    2 sinlam = dsin(lam)
***** if(dabs(pi-dabs(l)).lt.0.2d-10) sinlam=0.d0        !*** no--troublesome
      coslam = dcos(lam)
      temp = cosu1*sinu2-sinu1*cosu2*coslam
      sinsig = dsqrt((cosu2*sinlam)**2+temp**2)          !*** v14   (rapp part II)
      cossig = sinu1*sinu2+cosu1*cosu2*coslam            !*** v15   (rapp part II)
      sig = datan2(sinsig,cossig)
***** sinal = cosu1*cosu2*sinlam/sinsig                  !*** v17   (rapp part II)
      if(dabs(sinsig).lt.eps) then
        sinal = cosu1*cosu2*sinlam/dsign(eps,sinsig)     !*** avoid div 0
      else
        sinal = cosu1*cosu2*sinlam/sinsig
      endif
      cosal2 = -sinal**2+1.d0
***** costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig      !*** v18   (rapp part II)
      if(dabs(cosal2).lt.eps) then
        costm= -2.d0*(sinu1*sinu2/dsign(eps,cosal2))+cossig    !*** avoid div 0
      else
        costm= -2.d0*(sinu1*sinu2/cosal2)+cossig
      endif
      costm2=costm*costm
      c=((-3.d0*cosal2+4.d0)/rf+4.d0)*cosal2/rf/16.d0    !*** v10   (rapp part II)

*** entry point of the antipodal loop (kind=2)

    6 it=it+1
      d=(((2.d0*costm2-1.d0)*cossig*c+costm)*sinsig*c+sig)*(1.d0-c)/rf  !*** v11
      if(kind.eq.1) then
        lam=l+d*sinal
        if(dabs(lam-test).lt.tol) go to 100
        if(dabs(lam).gt.pi) then
          kind=2
          lam=pi
          if(l.lt.0.d0) lam=-lam
          sinal = 0.d0
          cosal2 = 1.d0
          test = 2.d0
          prev = test
          sig= pi - dabs(datan(sinu1/cosu1) + datan(sinu2/cosu2))
          sinsig = dsin(sig)
          cossig = dcos(sig)
          c=((-3.d0*cosal2+4.d0)/rf+4.d0)*cosal2/rf/16.d0
          if(dabs(sinal-prev).lt.tol) go to 100
*****     costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig
          if(dabs(cosal2).lt.eps) then
            costm= -2.d0*(sinu1*sinu2/dsign(eps,cosal2))+cossig    !*** avoid div 0
          else
            costm= -2.d0*(sinu1*sinu2/cosal2)+cossig
          endif
          costm2=costm*costm
          go to 6
        endif
        if(((lam-test)*(test-prev)).lt.0.d0.and.it.gt.5)
     *       lam=(2.d0*lam+3.d0*test+prev)/6.d0            !*** refined converge.
        prev = test
        test = lam
        go to 2
      else
        sinal= (lam-l)/d
        if(((sinal-test)*(test-prev)).lt.0.d0.and.it.gt.5)
     *       sinal=(2.d0*sinal+3.d0*test+prev)/6.d0        !*** refined converge.
        prev = test
        test = sinal
        cosal2= -sinal**2+1.d0
        sinlam= sinal*sinsig/(cosu1*cosu2)
        coslam = -dsqrt(dabs(-sinlam**2 +1.d0))
        lam= datan2(sinlam,coslam)
        temp = cosu1*sinu2-sinu1*cosu2*coslam
        sinsig = dsqrt((cosu2*sinlam)**2+temp**2)
        cossig = sinu1*sinu2+cosu1*cosu2*coslam
        sig = datan2(sinsig,cossig)
        c=((-3.d0*cosal2+4.d0)/rf+4.d0)*cosal2/rf/16.d0
        if(dabs(sinal-prev).lt.tol) go to 100
*****   costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig
        if(dabs(cosal2).lt.eps) then
          costm= -2.d0*(sinu1*sinu2/dsign(eps,cosal2))+cossig    !*** avoid div 0
        else
          costm= -2.d0*(sinu1*sinu2/cosal2)+cossig
        endif
        costm2=costm*costm
        go to 6
      endif

*** convergence

  100 if(kind.eq.2) then                        !*** antipodal
        faz= sinal/cosu1
        baz= dsqrt(-faz**2+1.d0)
        if(temp.lt.0.d0) baz=-baz
        faz=datan2(faz,baz)
        tem1=-sinal
        tem2=sinu1*sinsig-cosu1*cossig*baz
        baz=datan2(tem1,tem2)
      else                                      !*** long-line
        tem1=cosu2*sinlam
        tem2=cosu1*sinu2-sinu1*cosu2*coslam
        faz=datan2(tem1,tem2)
        tem1=-cosu1*sinlam
        tem2=sinu1*cosu2-cosu1*sinu2*coslam
        baz=datan2(tem1,tem2)
      endif
      if(faz.lt.0.d0) faz= faz+pi+pi
      if(baz.lt.0.d0) baz= baz+pi+pi

*** Helmert 1880 from Vincenty "Geodetic inverse solution between antipodal points"

      ep2 =1.d0/(boa*boa)-1.d0
      bige=dsqrt(1.d0+ep2*cosal2)                                         !*** 15
      bigf=(bige-1.d0)/(bige+1.d0)                                        !*** 16
      biga=(1.d0+bigf*bigf/4.d0)/(1.d0-bigf)                              !*** 17
      bigb=bigf*(1.d0-0.375d0*bigf*bigf)                                  !*** 18
      z   =bigb/6.d0*costm*(-3.d0+4.d0*sinsig**2)*(-3.d0+4.d0*costm2)
      dsig=bigb*sinsig*(costm+bigb/4.d0*(cossig*(-1.d0+2.d0*costm2)-z))   !*** 19
      s=(boa*a)*biga*(sig-dsig)                                           !*** 20

      return
      end

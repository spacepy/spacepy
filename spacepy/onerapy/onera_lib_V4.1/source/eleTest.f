C       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C testing the sgp4_ele routine
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      PROGRAM Test

      IMPLICIT NONE

      INTEGER*4 nmax,ii,numpoints,sysaxes,debug
      PARAMETER (nmax=100000)
      INTEGER*4 iyr, imon, iday, ihr, imin
      REAL*8 sec,epoch
      REAL*8 e1, e2, e3, e4, e5, e6
      REAL*8 startsfe, stopsfe, deltasec
      INTEGER*4 options(5), fnameLength
      PARAMETER (fnameLength = 39)
      INTEGER*4 iyrArray(nmax),idoyArray(nmax)
      INTEGER*4 iaYear(nmax),iaDoy(nmax)
      REAL*8 ut(nmax),x1(nmax),x2(nmax),x3(nmax)

      CHARACTER*500 fname,elefname
      REAL*8 trj(9,nmax) ! this is to keep the trajectory info

c     Now for the make_lstar stuff
      INTEGER*4 ntime, kext
      REAL*8 Lm(nmax),Lstar(nmax),Blocal(nmax)
      REAL*8 Bmin(nmax)
      REAL*8 XJ(nmax),MLT(nmax),maginput(25,nmax)


      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'

C epoch for ISS
      iyr = 1998
      imon = 11
      iday = 20
      ihr = 6
      imin = 50
      sec = 0.0d0
C epoch for GOES-8
      iyr = 2002
      imon = 11
      iday = 11
      ihr = 7
      imin = 10
      sec = 34.6817d0

      sysaxes = 0   ! see thoughts file for real conversion indices.  

c     compute epoch from inputs.  epoch in Julian Days.  
      CALL JDAY(iyr,imon,iday,ihr,imin,sec,epoch) 
      PRINT*,'JDAY = ',epoch

      ! set up the options array
      options(1) = 3  ! which standard element set?  
      options(2) = 1  ! for classical, x5 = omega (=1) or PI(=2)?
      options(3) = 1  ! for classical, x6 = T(=1), nu0(=2), u0(=3), l0(=4),or M0(=5)?
      options(4) = 0  ! unused
      options(5) = 0  ! unused

      startsfe = 0.0d0
      stopsfe = 5460.d0   ! oh, these are in fact in seconds.  Which version of the  
      stopsfe = 7200.d0   ! oh, these are in fact in seconds.  Which version of the 
      stopsfe = 86400.d0/4d0   ! 6 hours
      deltasec = 60.d0      ! codes am I using?    only compute for two hours (7200 sec)

      numpoints = (stopsfe - startsfe)/deltasec
C Now use GOES-8 for these two-line elements.
c 1 23051U 94022A   02315.29901252 -.00000255  00000-0  00000+0 0  4660
c 2 23051   0.3799  97.7687 0002018 118.0344 226.9733  1.00281322 38818


      if (options(1) .eq. 1) THEN
c                                ! onera elements  I'm sure this is right for ISS on 11/20/98
c         e1 = 51.5908d0        ! inclination for station
c         e2 = 323.4924d0        !  rp
c         e3 = 360.3654d0        !  ra
c         e4 = 168.3788d0         !  ascending node
c         e5 = 86.41859d0        !  arument of perigee
c         e6 = 359.7454d0        ! mean anomaly

         e1 = 0.3799d0      ! inclination for G8
         e2 = 35775.438d0            !  ap 
         e3 = 35792.454d0            !  aa 
         e4 = 97.7687d0        !  ascending node
         e5 = 118.0344d0      !  arument of perigee
         e6 = 226.9733d0        ! mean anomaly
         e6 = 180d0
         
         elefname = '/cygdrive/c/work/oneraDevel/g8oner.ele'
         elefname = elefname(1:fnameLength) 
         OPEN(UNIT=3,FILE=elefname,STATUS='UNKNOWN',ACCESS='SEQUENTIAL') 
         PRINT*,'=============output elements================'
         WRITE( *,'(6F17.6)') e1,e2,e3,e4,e5,e6
         WRITE( 3,'(6F17.6)') e1,e2,e3,e4,e5,e6
         CLOSE(3)


      ELSE IF ((options(1).eq.2)) THEN
                                ! classical elements with argument of perigee (omega,e5)
c         e1 = 1.0536d0          ! a
c         e2 = .0125362d0        !  e
c         e3 = 51.5908d0         !  i
c         e4 = 168.3788d0        !  ascending node (BigOmega)
c         e5 = 86.41859d0        !  arument of perigee  (omega)
c         e6 = epoch + (92.7127d0/(24.d0*60.d0)) ! time of perigee passage, Jdays 
c         e6 = epoch
         
         e1 = 6.6106d0      ! semimajor axis for G8
         e2 = 0.0002018d0  ! e
         e3 = 0.3799d0  ! i
         e4 = 97.7687d0  ! Omega
         e5 = 118.0344d0  ! omega  
         e6 = epoch-0.4986   ! OH, we're defining the elements at APOGEE, 1/2 a 
                               ! period away from perigee (epoch).  0.4986 is 
                                ! half what APEXRAD gives me as a P
         
         elefname = '/cygdrive/c/work/oneraDevel/g8clas.ele'
         elefname = elefname(1:fnameLength) 
         OPEN(UNIT=3,FILE=elefname,STATUS='UNKNOWN',ACCESS='SEQUENTIAL') 
         PRINT*,'=============output elements================'
         WRITE( *,'(6F17.6)') e1,e2,e3,e4,e5,e6
         WRITE( 3,'(6F17.6)') e1,e2,e3,e4,e5,e6
         CLOSE(3)


         IF (options(2).eq.1) THEN
            ! Default.  We've already set e5 = omega
         ELSE IF (options(2).eq.2) THEN
           !  Longitude of periapsis (PI), instead of omega (PI=Omega+omega)
            e5 = 254.7974d0    
         ELSE
            PRINT*,'options(2) takes values from 1-2, not',options(2)
         ENDIF

         IF (options(3).eq.1) THEN
            ! If options(3).eq.1, (default) we've already set e6 = T
         ELSE IF (options(3).eq.2) THEN
                                ! classical elements with nu0 instead of T
            e6 = 80             ! nu0  true anomaly given in degrees
            
         ELSE IF (options(3).eq.3) THEN
                                ! classical elements with u0 instead of T
            e6 = 120            ! u0  ; figure out what u0 should really be.
            
         ELSE IF (options(3).eq.4) THEN
                                ! classical elements with l0 instead of T
            e6 = 170            ! l0  ; figure out what l0 should really be.
            
         ELSE IF (options(3).eq.5) THEN
                                ! classical elements with M0 instead of T
            e6 = 359.7454d0     ! M0  ; figure out what M0 should really be.
         ELSE
            PRINT*,'options(3) takes values from 1-5, not',options(3)
         ENDIF   ! ends options(3) logic

      ELSE IF (options(1).eq.3) THEN
                                ! r,v
c         e1 = 1.0536d0*rekm     ! assume ISS is at noon GEO in km
c         e2 = 0.d0   
c         e3 = 0.d0
c         e4 = 0.0d0
c         e5 = 7.69d0           ! and going this fast in km/sc
c         e6 = 0.01d0

         e1 = 6.6106d0*rekm     ! assume G8 is at noon at apogee in km
         e2 = 0.d0   
         e3 = 0.d0
         e4 = 0.0d0
         e5 = 3.5804d0           ! and going this fast in km/sc
         e6 = 0.0d0

         elefname = '/cygdrive/c/work/oneraDevel/g8arve.ele'
         elefname = elefname(1:fnameLength) 
         OPEN(UNIT=3,FILE=elefname,STATUS='UNKNOWN',ACCESS='SEQUENTIAL') 
         PRINT*,'=============output elements================'
         WRITE( *,'(6F17.6)') e1,e2,e3,e4,e5,e6
         WRITE( 3,'(6F17.6)') e1,e2,e3,e4,e5,e6
         CLOSE(3)

      ELSE IF (options(1).eq.4) THEN
                                ! solar elements
c         e1 = 51.5908d0             ! inclination for station
c         e2 = 323.4924d0            !  A_p
c         e3 = 360.3654d0            !  A_a
c         e4 = 12.0d0            ! local time of apogee (hrs)
c         e5 = 22.0d0            ! local time of maximum inclination (hrs)
c         e6 = startsfe + (92.7127d0/(24.d0*60.d0))  ! time since perigee passage, in days
c         e6 = epoch
        ! now do it for G8
         e1 = 0.3799d0             ! inclination for G8
c         e2 = 35775.438d0            !  A_p
c         e3 = 35792.454d0            !  A_a
         e2 = 35719d0            !  A_p
         e3 = 35736d0            !  A_a
         e4 = 23.05d0            ! local time of apogee (hrs)
         e5 = 9.18d0            ! local time of maximum inclination (hrs)
         e6 = startsfe + (92.7127d0/(24.d0*60.d0))  ! time since perigee passage, in days
         e6 = epoch-0.4986   ! Oh, defining elements at apogee, 1/2 period away 
                             ! from perigee.  Get this from the classical ele. above.  

         elefname = '/cygdrive/c/work/oneraDevel/g8solr.ele'
         elefname = elefname(1:fnameLength) 
         OPEN(UNIT=3,FILE=elefname,STATUS='UNKNOWN',ACCESS='SEQUENTIAL') 
         PRINT*,'=============output elements================'
         WRITE( *,'(6F17.6)') e1,e2,e3,e4,e5,e6
         WRITE( 3,'(6F17.6)') e1,e2,e3,e4,e5,e6
         CLOSE(3)


      ELSE IF (options(1).eq.5) THEN
                                ! mean elements
c         e1 = 16.05064833d0             ! n (mean motion in revs/day
c         e2 = .0125362d0            !  e
c         e3 = 51.5908d0            !  i
c         e4 = 168.3788d0           ! Omega
c         e5 = 86.41859d0           ! omega
c         e6 = 359.7454d0 ! M0 (mean anomaly)
c         e6 = 180d0
         ! now do it for GOES-8
         e1 = 1.00281322d0   ! n (mean motion revs/day)
         e2 = 0.0002018d0  ! e
         e3 = 0.3799d0  ! i
         e4 = 97.7687d0  ! Omega
         e5 = 118.0344d0  ! omega  
         e6 = 226.9733d0  ! M0
         e6 = 180d0

         elefname = '/cygdrive/c/work/oneraDevel/g8mean.ele'
         elefname = elefname(1:fnameLength) 
         OPEN(UNIT=3,FILE=elefname,STATUS='UNKNOWN',ACCESS='SEQUENTIAL') 
         PRINT*,'=============output elements================'
         WRITE( *,'(6F17.6)') e1,e2,e3,e4,e5,e6
         WRITE( 3,'(6F17.6)') e1,e2,e3,e4,e5,e6
         CLOSE(3)


      ELSE

         PRINT*,'options(1) takes values from 1-5, not',options(1)
         
      ENDIF
      
      CALL  sgp4_ele1(sysaxes,iyr,imon,iday,ihr,imin,sec, 
     &     e1,e2,e3,e4,e5,e6,options, 
     &     startsfe,stopsfe,deltasec,
     &     iaYear,iaDoy,ut,x1,x2,x3)
      
      debug = 1
      IF (debug .eq. 1) THEN
                                ! write out test file
         fname = '/cygdrive/c/work/oneraDevel/output.trj'
         fname = fname(1:fnameLength) 
         OPEN(UNIT=1,FILE=fname,STATUS='UNKNOWN',ACCESS='SEQUENTIAL') 
         elefname = '/cygdrive/c/work/oneraDevel/output.ele'
         elefname = elefname(1:fnameLength) 
         OPEN(UNIT=2,FILE=elefname,STATUS='UNKNOWN',ACCESS='SEQUENTIAL') 
         PRINT*,'=============output elements================'
         WRITE( *,'(6F17.6)') e1,e2,e3,e4,e5,e6
         WRITE( 2,'(6F17.6)') e1,e2,e3,e4,e5,e6
         CLOSE(2)
      ENDIF
      DO ii = 1,numpoints
         
         WRITE( *,'(I4,1x,I3,1x,F17.6,1x,3F17.6)' )
     &        iaYear(ii),iaDoy(ii),ut(ii),x1(ii),x2(ii),x3(ii)
         
         IF (debug .eq. 1) THEN
                                ! also write to file for matlab inspection
            WRITE( 1,'(I4,1x,I3,1x,F17.6,1x,3F17.6)' )
     &           iaYear(ii),iaDoy(ii),ut(ii),x1(ii),x2(ii),x3(ii)
         ENDIF
      ENDDO
      IF (debug.eq.1) THEN
         CLOSE(1)
      ENDIF
      
    
c     Now use the outputs from sgp4_ele to call make_lstar

      ntime = numpoints
      ntime = 10
      kext = 2
c      sysaxes = 2   ! provided earlier
      options(1) = 1
      options(2) = 0
      options(3) = 0
      options(4) = 0
      options(5) = 0

      PRINT*, 'Calculating Lstar...'
      call make_lstar1(ntime,0,options,sysaxes,iaYear,iaDoy,ut, 
     +       x1,x2,x3, maginput, lm,lstar,blocal,bmin,xj,mlt)

      
      DO ii = 1,ntime
         WRITE( *,'(I4,1x,I3,1x,F17.6,1x,4E10.3)' )
     &        iaYear(ii),iaDoy(ii),ut(ii),lm(ii),lstar(ii),
     &        blocal(ii),bmin(ii)
         PRINT*, 'lstar = ', lstar(ii)
      ENDDO


      
      END

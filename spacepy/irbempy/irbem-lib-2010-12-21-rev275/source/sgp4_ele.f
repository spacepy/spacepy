!***************************************************************************************************
! Copyright 2001 D. Vallado, 2007 S. Bourdarie, T. Guild
!
! This file is part of IRBEM-LIB.
!
!    IRBEM-LIB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    IRBEM-LIB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
!
C       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE sgp4_ele
C      Routine to convert from many different types of 
C      orbital elements to the set used by sgp4_orb, the 
C      ONERA orbit propagator.  Then calls the orbit propagator 
C      and rotates the resultant trajectory into whichever 
C      sysaxes is specified.  Output times and trajectory returned 
C      to the calling program.  
C
C     ele_opts (element options) array description:
c      !  first element in ele_opts array corresponds to element type:
c      !  1 = ONERA; requires i (inclination; deg), A_p (altitude of perigee; km), 
c      !      A_a (altitude of apogee; km), Omega (longitude of ascending node; deg), 
c      !      omega (argument of perigee; deg), M0 (mean anomaly at epoch; deg)
c      !  2 = CLASSICAL; requires a (semimajor axis; Re), e (eccentricity), 
c      !      i (inclination; deg), Omega (longitude of ascending node; deg), 
c      !      omega (argument of perigee; deg), tsfe (time of perigee in seconds 
c      !      from epoch; sec)
c      !  3 = RV;  requires vector r (km), v (km/s)
c      !  4 = SOLAR;  requres i (inclination; deg), A_p (altitude of perigee; km), 
c      !      A_a (altitude of apogee; km), H_a (local time of apogee; hrs), 
c      !      H_i (local time of maximum inclination; hrs), 
c      !      tsfe (time of perigee in seconds from epoch; sec)
c      !  5 = MEAN;  requires n (mean motion; rev/day), e (eccentricity), 
c      !      i (inclination; deg), Omega (longitude of ascending node; deg), 
c      !      omega (argument of perigee; deg), M (mean anomaly at epoch; deg)
c      
c      ! second element in the ele_opts array switches between e5=omega and PI
c      !    ele_opts(2)=1 means e5 is omega the argument of perigee (deg)
c      !    ele_opts(2)=2 means e5 is PI the longitude of perigee (deg)
c      !    ele_opts(2) only has effect for ele_opts(1) = 1, 2, 5
c
c      ! third element in the ele_opts array switches among different e6
c      !    ele_opts(3)=1 means e6 is tsfe the time of perigee in seconds from epoch (secs)
c      !    ele_opts(3)=2 means e6 is nu0 the true anomaly at epoch (deg)
c      !    ele_opts(3)=3 means e6 is u0 the argument of latitude at epoch (deg)
c      !    ele_opts(3)=4 means e6 is l0 the true longitude at epoch (deg) 
c      !    ele_opts(3)=5 means e6 is M0 the mean anomaly at epoch
c      !    ele_opts(2) only has effect for ele_opts(1) = 1, 2, 4, 5
C
C      Calling sequence from fortran:
C           
C      CALL  sgp4_ele1(sysaxes,iyr,imon,iday,ihr,imin,sec, 
C     &     e1,e2,e3,e4,e5,e6,options, 
C     &     startsfe,stopsfe,deltasec,
C     &     iaYear,iaDoy,ut,x1,x2,x3)
C
C         with the following types:
C      INTEGER*4 sysaxes
C      INTEGER*4 iyr,imon,iday,ihr,imin
C      REAL*8 sec
C      REAL*8 e1,e2,e3,e4,e5,e6
C      INTEGER*4 ele_opts(5)
C      REAL*8 startsfe, stopsfe, deltasec
C      INTEGER*4 iaYear(ntime_max),iaDoy(ntime_max)
C      REAL*8 ut(ntime_max),x1(ntime_max),x2(ntime_max),x3(ntime_max)
C
C     As with all (most?) onera library calls, the maximum array size
C     is limited to ntime_max elements.  If startsfe/stopsfe/deltasec 
C     are given such that the trajectory exceeds ntime_max positions, 
C     only the first ntime_max will be returned.  Errors will be issued 
C     to STDOUT, but not through IDL/Matlab wrappers.   
C                              Contributed by Timothy Guild, 3.7.07
C                               timothy.guild@aero.org
C     last updated 10.9.07 by tbg
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      SUBROUTINE sgp4_ele1(sysaxes,iyr,imon,iday,ihr,imin,sec, 
     +                       e1,e2,e3,e4,e5,e6,ele_opts, 
     +                       startsfe,stopsfe,deltasec,
     +                       iaYear,iaDoy,ut,x1,x2,x3)

      IMPLICIT NONE
      INCLUDE 'ntime_max.inc'

      INTEGER*4 ii,numOutputLines,sysaxes,debug
      INTEGER*4 iyr, imon, iday, ihr, imin,idoy
      INTEGER*4 iaYear(ntime_max),iaDoy(ntime_max)
      REAL*8 ut(ntime_max),x1(ntime_max),x2(ntime_max),x3(ntime_max)
      REAL*8 sec
      REAL*8 e1, e2, e3, e4, e5, e6
      REAL*8 keepe1, keepe2, keepe3, keepe4, keepe5, keepe6
      REAL*8 o1, o2, o3, o4, o5, o6
      REAL*8 a, e, i, BigOmega, omega, tsfe, BigPi
      REAL*8 nu0,u0,l0,M0,P,E0
      REAL*8 startsfe, stopsfe, deltasec,epoch
      INTEGER*4 ele_opts(5)
      EXTERNAL GET_DOY
      INTEGER*4 GET_DOY
      real*8     xGEO(3), MLT
      REAL*8 output(10,ntime_max) ! keeps the trajectory info from SGP4_orb1
      REAL*8 xINV(3,ntime_max),xOUTV(3,ntime_max)
      INTEGER*4 sysaxesFROM,sysaxesTO

      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'

      ! TBG, 10.9.07:  save the initial input elements to pass back, as per 
      ! Sebastien's suggestion.  
      keepe1 = e1
      keepe2 = e2
      keepe3 = e3
      keepe4 = e4
      keepe5 = e5
      keepe6 = e6

      debug = 0
      IF (debug.eq.1) THEN
                                ! check the values coming in
         PRINT*,'SGP4_ELE:  INPUT ARGUMENTS:'
         PRINT*,'iyr,imon,iday,ihr,imin,sec:',iyr,imon,iday,ihr,imin,sec
         PRINT*,'e1,e2,e3,e4,e5,e6:',e1,e2,e3,e4,e5,e6
         PRINT*,'ele_opts:',ele_opts
         PRINT*,'startsfe, stopsfe, deltasec:',startsfe,stopsfe,deltasec
      ENDIF
      
c     compute epoch from inputs.  epoch in Julian Days.  
      CALL JDAY(iyr,imon,iday,ihr,imin,sec,epoch) 


      IF (ele_opts(1) .eq. 1) THEN

         PRINT*, 'Converting these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         call fifthElement(e1,e2,e3,e4,e5,e6,ele_opts)
         call sixthElement(e1,e2,e3,e4,e5,e6,ele_opts)         
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         
      ELSE IF (ele_opts(1) .eq. 2) THEN
         PRINT*, 'Converting these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         CALL classical2onera(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
         call fifthElement(e1,e2,e3,e4,e5,e6,ele_opts)
         call sixthElement(e1,e2,e3,e4,e5,e6,ele_opts) 
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6

      ELSE IF (ele_opts(1) .eq. 3) THEN
         PRINT*, 'Converting these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         CALL rv2classical(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         CALL classical2onera(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         call fifthElement(e1,e2,e3,e4,e5,e6,ele_opts)
         call sixthElement(e1,e2,e3,e4,e5,e6,ele_opts) 
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6            

      ELSE IF (ele_opts(1) .eq. 4) THEN
         PRINT*, 'Converting these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         CALL solar2classical(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         CALL classical2onera(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         call fifthElement(e1,e2,e3,e4,e5,e6,ele_opts)
         call sixthElement(e1,e2,e3,e4,e5,e6,ele_opts) 
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6

      ELSE IF (ele_opts(1) .eq. 5) THEN
         PRINT*, 'Converting these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         CALL mean2classical(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         CALL classical2onera(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6
         call fifthElement(e1,e2,e3,e4,e5,e6,ele_opts)
         call sixthElement(e1,e2,e3,e4,e5,e6,ele_opts) 
         PRINT*, 'to these ele_opts = ',ele_opts
         PRINT*, 'with these elements = ',e1,e2,e3,e4,e5,e6

      ELSE 
         PRINT*,'ele_opts(1) must be 1-5'
      ENDIF

      IF (debug .eq. 1) THEN
         PRINT*, 'elements just before SGP4_ORB1 call---------------'
         PRINT*,'iyr,imon,iday,ihr,imin,sec:',iyr,imon,iday,ihr,imin,sec
         PRINT*,'e1,e2,e3,e4,e5,e6:',e1,e2,e3,e4,e5,e6
         PRINT*,'ele_opts:',ele_opts
         PRINT*,'startsfe, stopsfe, deltasec:',startsfe,stopsfe,deltasec
      ENDIF

c     compute how many lines the output array should have here
      numOutputLines = (stopsfe - startsfe)/deltasec

      IF (numOutputLines .gt. ntime_max) THEN
         PRINT*,'SGP4_ELE:  Given > ntime_max timesteps in the orbit'
         PRINT*,'propagation.  ONERA library requires less than that,'
         PRINT*,' so SGP4_ELE has truncated your request.  '
         stopsfe = startsfe + ntime_max*deltasec ! this forces the stop time to be 
                                ! at most NTIME_MAX times, as required by ONERA lib
      ENDIF

c     Now call the orbit propagator with the epoch and proper 
C     orbital elements computed above.  Trajectory timeseries 
c     returned in the 'output' array.  
      CALL SGP4_ORB1(iyr,imon,iday,ihr,imin,sec, e1,e2,e3,e4,e5,e6, 
     +     startsfe,stopsfe,deltasec,sysaxes,output)
      
c     unpack the 'output' array into arguments make_lstar uses.  
      DO ii = 1,numOutputLines
         iaYear(ii) = INT(output(3,ii))   ! array of years
         iaDoy(ii) = GET_DOY(INT(output(3,ii)),INT(output(2,ii)),
     &        INT(output(1,ii)))   !  array of DOYs
         x1(ii) = output(8,ii)  ! input X coordinate
         x2(ii) = output(9,ii)  ! input Y coordinate
         x3(ii) = output(10,ii) ! input Z coordinate
         ut(ii) = output(4,ii)*3600d0+output(5,ii)*60d0+output(6,ii)
      ENDDO

      ! TBG, 10.9.07: now assign the original elements back to the e1-e6 variables
      e1 = keepe1
      e2 = keepe2
      e3 = keepe3
      e4 = keepe4
      e5 = keepe5
      e6 = keepe6

      END

C********************************************************************************
C Now for the subroutines which convert between element sets, and various 
C supporting subroutines
C********************************************************************************

C********************************************************************************
C     CLASSICAL2ONERA
C********************************************************************************
 
      SUBROUTINE classical2onera(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
      
      IMPLICIT NONE

      INTEGER*4 debug
      INTEGER*4 ele_opts(5)
      REAL*8 e1, e2, e3, e4, e5, e6
      REAL*8 out1, out2, out3, out4, out5, out6
      REAL*8 r_a, r_p,epoch,tsfe,P,a,e,i,A_p,A_a
      REAL*8 nu0,u0,l0,M0,bigomega,bigpi,omega

      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'

      debug = 0

C     need to convert a,   e,   i,   Omega,   omega (or PI), tsfe (or nu0,u0,l0,M0) to 
C     to              i,  A_p, A_a,  Omege,   omega (or PI), M0 (or nu0,u0,l0,tsfe)

C     First, do the first four elements.  

      i = e3  ! inclination is given
      a = e1  !  semimajor axis given
      e = e2  ! eccentricity given
      bigomega = e4  ! longitude of ascending node given

      r_a = a*(1.d0+e)  ! radius of apogee is a(1+e), in Re
      r_p = a*(1.d0-e)  ! radius of perigee is a(1-e), in Re

      A_p = (r_p-1.d0)*rekm  ! perigee altitude in km
      A_a = (r_a-1.d0)*rekm  ! apogee altitude in km  

      e1 = i
      e2 = A_p
      e3 = A_a
c      e4 = bigomega   ! pass these numbers through.  
c      e5 = omega
c      e6 = M0  

      ele_opts(1) = 1   ! save for the last two elements, this is now ONERA elements. 

      END


C********************************************************************************
C     RV2CLASSICAL   
C
C     TBG:  Wait, sgp4ext.f contains a subroutine rv2coe which does this.
C     Just manage the bookkeeping and call it in this subroutine.   
C********************************************************************************

      SUBROUTINE rv2classical(e1,e2,e3,e4,e5,e6,epoch,ele_opts)

      IMPLICIT NONE
      INTEGER*4 ii,debug,ele_opts(5)
      REAL*8 e1, e2, e3, e4, e5, e6,epoch
      REAL*8 nu0,u0,l0,M0,Period,E0,tsfe
      EXTERNAL DOT, MAG, RV2COE
      REAL*8 DOT, MAG

      ! variables for call to rv2coe
      REAL*8 R(3), V(3), P, A, Ecc, Incl, Omega, Argp, Nu, M, ArgLat,
     &         TrueLon, LonPer

      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'

      debug = 1

c     convert RV:  r(1), r(2), r(3), v(1), v(2), v(3) where r (km) and v (km/s)
c     to Classical:  a, e, i, Omega, omega (or PI), tsfe (or nu0,u0,l0,M0)
      
      R(1) = e1
      R(2) = e2
      R(3) = e3
      V(1) = e4
      V(2) = e5
      V(3) = e6

      IF (debug .eq. 1) THEN
         PRINT*, 'R, V = ',R, V
      ENDIF
C just call rv2coe from sgp4ext.f in the ONERA library
      CALL rv2coe ( R, V, P, A, Ecc, Incl, Omega, Argp, Nu,
     &                         M, ArgLat, TrueLon, LonPer )

      IF (debug .eq. 1) THEN 
         PRINT*, 'P,A,Ecc,Incl,Omega,Argp,Nu,M,ArgLat,TrueLon,LonPer',
     &        P,A,Ecc,Incl,Omega,Argp,Nu,M,ArgLat,TrueLon,LonPer
      ENDIF

      e1 = a / (rekm)  ! convert semimajor axis back to RE
      e2 = Ecc
      e3 = Incl*Rad2Deg   ! convert from Rad to degrees
      e4 = Omega*Rad2Deg  ! convert from Rad to degrees
      e5 = argp*Rad2Deg   ! convert from Rad to degrees
      e6 = M*Rad2Deg    ! mean anomaly

      ele_opts(1) = 2    !   now in classical elements
      ele_opts(2) = 1   ! with omega instead of PI
      ele_opts(3) = 5   ! and e6 = mean anomaly calculated from rv2coe

      END


C********************************************************************************
C     SOLAR2CLASSICAL
C********************************************************************************

      SUBROUTINE solar2classical(e1,e2,e3,e4,e5,e6,epoch,ele_opts)

      IMPLICIT NONE
      INTEGER*4 ii,day,epochDay,debug,ele_opts(5)
      REAL*8 e1, e2, e3, e4, e5, e6
      REAL*8 out1, out2, out3, out4, out5, out6
      REAL*8 p, a, i,BigOmega,omega,l0,u0,nu0,epoch
      REAL*8 startDate,hour,A_a, A_p,r_a,r_p,H_e,H_a,H_i
      REAL*8 Theta_a, Theta_i, e, M0, n,tsfe

      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'

      debug = 0

c     convert Solar: i, A_p, A_a, H_a, H_i, tsfe (or nu0, u0, l0, M0)
c     to Classical : a, e, i, Omega, omega, tsfe (or nu0, u0, l0, M0)

      A_a = e3
      A_p = e2
      H_a = e4
      H_i = e5
      i = e1

      CALL JDAY ( 1992,1,1,0,0,0.d0,  startDate )   ! B-3, APEXRAD
      epochDay = epoch
      day = 1 + epochDay - startDate !these are days and hours since startDate
      hour = 24.d0*DMOD(epoch,1.d0)

      r_a = rekm + A_a   ! B-4, B-5, APEXRAD
      r_p = rekm + A_p
      
      e = (r_a - r_p) / (r_a + r_p)    !  B-9, APEXRAD

      a = (r_a + r_p )/2.d0   

      IF (debug .eq. 1) THEN
         PRINT*,'SOLAR2CLASSICAL:  r_a, r_p,H_a,H_i epoch, epochDay= ',
     &        r_a,r_p,H_a,H_i,epoch,epochDay
      ENDIF 

      !  B-6, APEXRAD
c     B-6 in APEXRAD manual might have typo.  Paul thinks last term should 
c     be close to zero, not one, and the difference might be between UT and LT.
c      H_e = 6.594703 + 0.06570982463 * day + 1.00273791 * hour 
      H_e = 6.594703 + 0.06570982463 * day + 0.00273791 * hour 

      IF (debug .eq. 1) THEN
         PRINT*,'SOLAR2CLASSICAL:  day,hour,H_e = ',day,hour,H_e
      ENDIF
      Theta_a = (H_a + 12.d0) * (PI/12.d0) ! radians, B-7, B-8, APEXRAD
      Theta_i = (H_i -6.d0 ) * (PI/12.d0)
      
      IF (debug .eq. 1) THEN
         PRINT*,'S2C:  H_e, Theta_a, Theta_i = ',H_e, Theta_a, Theta_i
      ENDIF 

      IF (debug .eq. 1) THEN
         PRINT*,'S2C:  e, a, P, M0, tsfe = ',e, a, P, M0, tsfe
      ENDIF      
c     Otherwise, assume satellite at apogee at epoch (B-10, APEXRAD)
c      M0 = 180

      IF (debug .eq. 1)  THEN
         PRINT*,'SOLAR2Classical: r_a,r_p,A_a,A_p,e= ',r_a,r_p,A_a,A_p,e
      ENDIF
      n = SQRT((MU*(12.d0*3600.d0/PI)**2)/(a**3d0) ) ! B-11, APEXRAD
      BigOmega = 15.d0*(H_i - 6.d0 + H_e) ! B-12, APEXRAD
      BigOmega = DMOD(DMOD(BigOmega,360d0) + 360d0,360d0)  !  0 to 360
      IF (debug .eq. 1) THEN
         PRINT*,'S2C:  n,BigOmega = ',n, BigOmega
      ENDIF      
c     APEXRAD B-13 gives cos(omega) = cos(Theta_a)*cos(Theta_i) + sin(Theta_a)*sin(Theta_i)
c     Paul claims that the following is equivalent:
      omega = (Theta_a - Theta_i) * (180d0/PI)
      omega = DMOD(DMOD(omega,360d0) + 360d0,360d0)  !  0 to 360

      IF (debug .eq. 1) THEN
         PRINT*,'S2C: omega = ',omega
      ENDIF                     
      ! now take those mean elements and convert them to classical ones
      e1 = n
      e2 = e
      e3 = i
      e4 = BigOmega
      e5 = omega
c      e6 = M0   ! leave whatever e6 came in as.  

      ele_opts(1) = 5   ! elements are now mean elements.
      ele_opts(2) = 1   ! omega specified

      PRINT*, 'S2C: Converting mean ele 
     &     (n,e,i,Omega,omega,(tsfe,nu0,u0,l0,m0):',e1,e2,e3,e4,e5,e6
      CALL mean2classical(e1,e2,e3,e4,e5,e6,epoch,ele_opts)
      PRINT*, 'S2C: to Classical elements: ', e1,e2,e3,e4,e5,e6
      
      END

C********************************************************************************
C     MEAN2CLASSICAL
C********************************************************************************

      SUBROUTINE mean2classical(e1,e2,e3,e4,e5,e6,epoch,ele_opts)

      IMPLICIT NONE
      INTEGER*4 debug,ele_opts(5)
      REAL*8 e1, e2, e3, e4, e5, e6
      REAL*8 out1, out2, out3, out4, out5, out6
      REAL*8 p, a, i,BigOmega,omega,l0,u0,nu0,epoch
      REAL*8 n, acubed,M0,tsfe

      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'

      debug = 0

c     mean elements are n, e, i, BigOmega, omega, M0
c     classical are a, e, i, BigOmega, omega, tsfe  

      ! first convert the mean motion n into semimajor axis a
      n = e1*(TwoPi/86400.d0)  ! mean motion converted from revs/day to rads/s
      acubed = mu / (n**2.d0)  ! Eq. 6-13, SMAD
      a = acubed**(1.d0/3.d0)

      IF (debug .eq. 1) THEN
         PRINT*,'M2C:  in1, n, a = ',e1, n, a
         PRINT*,'M2C:  p, M0, tsfe = ',p, M0, epoch, tsfe
      ENDIF

      e1 = a / rekm  ! convert back to RE
c      e6 = tsfe   ! still undefined 

      ele_opts(1) = 2   ! elements now classical, leave ele_opts(2,3) as they came in.  

      END


C********************************************************************************
C     tsfe2M0:  calculate mean anomaly (M0) given time of perigee
C********************************************************************************
      SUBROUTINE tsfe2M0(tsfe,a,M0)
      
      IMPLICIT NONE
      INTEGER*4 debug
      REAL*8 a,tsfe,M0,P,negM0
c      EXTERNAL REM   ! line commented by S. Bourdarie
c      REAL*8 REM   ! line commented by S. Bourdarie

      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'      
      
      P = twopi*DSQRT(a**3/mu) ! in seconds
      M0 = -360.d0*(tsfe/P)    
      M0 = DMOD(DMOD(M0,360.d0)+360.d0,360.d0) ! DMOD = REM for positive arguments only.  

      PRINT*,'Converting tsfe (',tsfe,') to M0 (',M0,')'

      END

C********************************************************************************
C     fifthElement:  determine which e5 was passed in, convert to omega
C********************************************************************************
      SUBROUTINE fifthElement(e1,e2,e3,e4,e5,e6,ele_opts)
      
      IMPLICIT NONE
      INTEGER*4 debug,ele_opts(5)
      REAL*8 e1,e2,e3,e4,e5,e6
      REAL*8 omega, BigPi
      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'      
      
      IF ((ele_opts(1).eq.1) .or. (ele_opts(1).eq.2) .or. 
     &     (ele_opts(1).eq.5)) THEN
         
         IF (ele_opts(2).eq.1) THEN ! omega
                                ! do nothing, omega is default
         ELSE IF (ele_opts(2).eq.2) THEN ! PI, longitude of perigee
            BigPi = e5
            omega = MOD(e5 - e4 + 360D0, 360D0)
            e5 = omega
            ele_opts(2) = 1
            
            PRINT*,'Converting PI (',BigPi,') to omega (',omega,')'
            
         ELSE 
            PRINT*,'ele_opts(2) must be either 1 (omega) or 2 (PI)'
            
         ENDIF
         
      ELSE
c               ele_opts(1) is 3 or 4, both have different values for e5.           
         PRINT*,'fifthElement: no change to ele_opts(1) = 3,4'   
                                
      ENDIF     
       
      END
      
C********************************************************************************
C     sixthElement:  determine which e6 was passed in, convert to M0
C********************************************************************************
      SUBROUTINE sixthElement(e1,e2,e3,e4,e5,e6,ele_opts)
      
      IMPLICIT NONE
      INTEGER*4 debug,ele_opts(5)
      REAL*8 e1,e2,e3,e4,e5,e6
      REAL*8 M0, e, nu0, u0, l0, r_a, r_p, a,tsfe
      REAL*8 n, acubed
      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'      
      
      IF (ele_opts(1).eq.1) THEN ! ONERA elements
         
         IF (ele_opts(3).eq.1) THEN ! tsfe
            r_p = e2 + rekm ! B-4, B-5, APEXRAD
            r_a = e3 + rekm
            a = (r_a + r_p) / 2.0D0
            CALL tsfe2M0(e6,a, M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.2) THEN !  nu0
            nu0 = e6
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.3) THEN !  u0
            nu0 = MOD(e6 - e5 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.4) THEN !  l0
            nu0 = MOD(e6 - e5 - e4 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.5) THEN !  M0
            ! ONERA element default, MO.  Do nothing.  
         ELSE
            PRINT*,'ele_opts(3) must be 1-5'
         ENDIF
         
      ELSE IF (ele_opts(1).eq.2) THEN ! Classical elements

         IF (ele_opts(3).eq.1) THEN ! tsfe
            r_p = e2 + rekm ! B-4, B-5, APEXRAD
            r_a = e3 + rekm
            a = (r_a + r_p) / 2.0D0
            CALL tsfe2M0(e6, a, M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.2) THEN !  nu0
            nu0 = e6
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            PRINT*,'e,nu0,m0 = ',e,nu0,m0
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.3) THEN !  u0
            nu0 = MOD(e6 - e5 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.4) THEN !  l0
            nu0 = MOD(e6 - e5 - e4 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.5) THEN !  M0
            ! ONERA element default, MO.  Do nothing.  
         ELSE
            PRINT*,'ele_opts(3) must be 1-5'
         ENDIF

         
      ELSE IF (ele_opts(1).eq.3) THEN !  RV elements
         
         PRINT*,'sixthElement:  No effect for RV elements'


      ELSE IF (ele_opts(1).eq.4) THEN ! Solar elements
         
         IF (ele_opts(3).eq.1) THEN ! tsfe
            r_p = e2 + rekm ! B-4, B-5, APEXRAD
            r_a = e3 + rekm
            a = (r_a + r_p) / 2.0D0
            CALL tsfe2M0(e6, a, M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.2) THEN !  nu0
            nu0 = e6
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.3) THEN !  u0
            nu0 = MOD(e6 - e5 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.4) THEN !  l0
            nu0 = MOD(e6 - e5 - e4 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.5) THEN !  M0
            ! ONERA element default, MO.  Do nothing.  
         ELSE
            PRINT*,'ele_opts(3) must be 1-5'
         ENDIF


      ELSE IF (ele_opts(1).eq.5) THEN ! Mean elements
         
         IF (ele_opts(3).eq.1) THEN ! tsfe
            n = e1*(TwoPi/86400.d0) ! mean motion converted from revs/day to rads/s
            acubed = mu / (n**2.d0) ! Eq. 6-13, SMAD
            a = acubed**(1.d0/3.d0)
            CALL tsfe2M0(e6, a, M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.2) THEN !  nu0
            nu0 = e6
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.3) THEN !  u0
            nu0 = MOD(e6 - e5 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.4) THEN !  l0
            nu0 = MOD(e6 - e5 - e4 + 360D0, 360D0)
            CALL computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
            CALL meanAnom(e,nu0,M0)
            e6 = M0
            ele_opts(3) = 5
         ELSE IF (ele_opts(3).eq.5) THEN !  M0
            ! ONERA element default, MO.  Do nothing.  
         ELSE
            PRINT*,'ele_opts(3) must be 1-5'
         ENDIF


      ELSE 
         PRINT*, 'ele_opts(1) betwee 1 and 5'
      ENDIF
      
      END

C********************************************************************************
C     computeE:  calculate eccentricity given ele_opts(1)
C********************************************************************************
      SUBROUTINE computeE(e1,e2,e3,e4,e5,e6,e,ele_opts)
      
      IMPLICIT NONE
      INTEGER*4 debug,ele_opts(5)
      REAL*8 e1,e2,e3,e4,e5,e6,e,r_a,r_p
      ! variables for call to rv2coe
      REAL*8 R(3),V(3),P,A,Ecc,Incl,Omega,Argp,Nu,
     &     M,ArgLat,TrueLon,LonPer
      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'      
      
      IF ((ele_opts(1).eq.1).or.(ele_opts(1).eq.4)) THEN
         r_p = e2 + rekm
         r_a = e3 + rekm
         e = (r_a - r_p) / (r_a + r_p)
      ELSE IF ((ele_opts(1).eq.2).or.(ele_opts(1).eq.5)) THEN
         e = e2
      ELSE IF (ele_opts(1).eq.3) THEN
         R(1) = e1
         R(2) = e2
         R(3) = e3
         V(1) = e4
         V(2) = e5
         V(3) = e6
         PRINT*,'computeE: ', R, V
         CALL rv2coe( R, V, P, A, Ecc, Incl, Omega, Argp, Nu,
     &                         M, ArgLat, TrueLon, LonPer )
         e = Ecc
      ELSE 
         PRINT*,'ele_opts(1) must be 1-5'
      ENDIF   


      END

C********************************************************************************
C     meanAnom:  compute mean anomaly given e, nu0, m0
C     Derived from Kepler library routine xkep_meananom.c
c  http://pds-rings.seti.org/toolkits/kepler_102/xkep_meananom.c
C********************************************************************************
      SUBROUTINE meanAnom(e,nu0,M0)
      
      IMPLICIT NONE
      INTEGER*4 debug
      REAL*8 e,nu0,M0,eRatio,tanhalfPsi,psi
      INCLUDE 'astconst.cmn'
      INCLUDE 'astmath.cmn'      

      ! convert nu0 to radians
      nu0 = nu0 * Deg2Rad

      eRatio = DSQRT( ( 1.D0 + e ) / (1.D0 - e ) )
      tanhalfPsi = DTAN( 0.5D0 * nu0 ) / eRatio
      psi = 2.D0 * DATAN(tanhalfPsi)
      M0 = psi - e * DSIN(psi)
      M0 = M0 * Rad2Deg  ! convert to degrees for consistency 

      END


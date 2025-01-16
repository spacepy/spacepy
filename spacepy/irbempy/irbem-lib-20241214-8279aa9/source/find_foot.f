!***************************************************************************************************
! Copyright 2004 S. Bourdarie
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
C S. Bourdarie (June 2004)
c Modified S./ Bourdarie (July 2004)
c Modified from find_bm and find_mirror_point by TPO (July 2008)
c Modified to prevent overshooting into the other hemisphere ACK (Nov 2016)
c
C Routine to find foot point of field line at specified altitude and hemi
C finds foot point at specified altitude to within 1 km
c --------------------------------------------------------------------
c
        SUBROUTINE find_foot_point1(kext,options,sysaxes,iyearsat,
     &  idoy,UT,xIN1,xIN2,xIN3,stop_alt,hemi_flag,maginput,
     &  XFOOT,BFOOT,BFOOTMAG)
c
c       INPUTS have the usual meaning, except:
c      REAL*8 stop_alt - geodetic altitude of desired foot point (gdz), km
c      integer*4 hemi_flag - hemishere flag, specifies hemisphere of foot point
c      0 - same Hemisphere as start point
c      +1 - Northern Hemisphere
c      -1 - Southern Hemisphere
c      2 - opposite Hemisphere as start point
c
c      OUTPUTS
c      REAL*8 XFOOT(3) - GDZ position of foot point (alt, lat, lon)
c      REAL*8 BFOOT(3) - Magnetic field at foot point (nT, GEO)
c      REAL*8 BFOOTMAG - Magnetic field at foot point (nT)

	IMPLICIT NONE
	INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,kint,options(5)
        INTEGER*4    sysaxes
	INTEGER*4    iyearsat
	integer*4    idoy
	real*8     UT
	real*8     xIN1,xIN2,xIN3
	real*8     stop_alt
        INTEGER*4  hemi_flag
	real*8     maginput(25)
c
c Declare internal variables
	INTEGER*4    isat,iyear,Iint
        INTEGER*4    Ndays,activ,Ifail
	REAL*8     psi,mlon,tilt
	REAL*8     xGEO(3)
	real*8     alti,lati,longi
	real*8     BxGEO(3)
c
c Declare output variables	
        REAL*8     XFOOT(3),BFOOT(3),BFOOTMAG
C
	COMMON /magmod/k_ext,k_l,kint
      integer*4 int_field_select, ext_field_select
C

        if ((stop_alt.lt.0).or.(stop_alt.ge.6378.0*500.0)) then
           goto 999 ! fail, stop_alt out of range 0 to 500 Re
        endif
	kint = int_field_select ( options(5) )
	k_ext = ext_field_select ( kext )
c
        CALL INITIZE
	
	call init_fields ( kint, iyearsat, idoy, ut, options(2) )
	
	call get_coordinates ( sysaxes, xIN1, xIN2, xIN3, 
     6    alti, lati, longi, xGEO )
	    
	call set_magfield_inputs ( k_ext, maginput, ifail )
 
	if ( ifail.lt.0 ) goto 999

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then 
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
            call INIT_TS07D_TLPR
	        if ( ifail.lt.0 ) goto 999
        end if
c
c
        CALL find_foot(lati,longi,alti,stop_alt,
     &       hemi_flag,XFOOT,BFOOT,BFOOTMAG)
        return
 999    continue ! bad data
        XFOOT(1) = baddata
        XFOOT(2) = baddata
        XFOOT(3) = baddata
        BFOOT(1) = baddata
        BFOOT(2) = baddata
        BFOOT(3) = baddata
        BFOOTMAG = baddata
	END
      
       SUBROUTINE find_foot(
     &     lati,longi,alti,stop_alt,hemi_flag,
     &     XFOOT,BFOOT,BFOOTMAG)
C
c      inputs: 
c      REAL*8 lati - geodetic latitude of start point (gdz), degrees
c      REAL*8 longi - geodetic longitude of start point  (gdz), degrees
c      REAL*8 alti - geodetic altitude of start point  (gdz), km
c      REAL*8 stop_alt - geodetic altitude of desired foot point (gdz), km
c      integer*4 hemi_flag - hemishere flag, specifies hemisphere of foot point
c      0 - same Hemisphere as start point
c      +1 - Northern Hemisphere
c      -1 - Southern Hemisphere
c      2 - opposite Hemisphere as start point
c
c      outputs:
c      REAL*8 XFOOT(3) - GDZ position of foot point (alt, lat, lon)
c      REAL*8 BFOOT(3) - Magnetic field at foot point (nT, GEO)
c      REAL*8 BFOOTMAG - Magnetic field at foot point (nT)
       IMPLICIT NONE

       INTEGER*4  k_ext,k_l,kint,Ifail
       REAL*8     xx0(3)
       REAL*8     lati,longi,alti
       REAL*8     stop_alt
       INTEGER*4  hemi_flag
       REAL*8     XFOOT(3),BFOOT(3),BFOOTMAG

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
C
       call find_foot_opt ( xx0,stop_alt,hemi_flag,
     &     XFOOT,BFOOT,BFOOTMAG)

       RETURN
       END

       SUBROUTINE find_foot_opt (
     &     xx0,stop_alt,hemi_flag,
     &     XFOOT,BFOOT,BFOOTMAG)
C
c      inputs: 
c      REAL*8 xx0(3) - GEO cartesian coordinates
c      REAL*8 stop_alt - geodetic altitude of desired foot point (gdz), km
c      integer*4 hemi_flag - hemishere flag, specifies hemisphere of foot point
c      0 - same Hemisphere as start point
c      +1 - Northern Hemisphere
c      -1 - Southern Hemisphere
c      2 - opposite Hemisphere as start point
c
c      outputs:
c      REAL*8 XFOOT(3) - GDZ position of foot point (alt, lat, lon)
c      REAL*8 BFOOT(3) - Magnetic field at foot point (nT, GEO)
c      REAL*8 BFOOTMAG - Magnetic field at foot point (nT)
       IMPLICIT NONE
       INCLUDE 'variables.inc'
C
       INTEGER*4  Nreb
       PARAMETER (Nreb = 50)
C
       INTEGER*4  Ifail
       REAL*8     rr,tt
       REAL*8     xx0(3),xx(3),x1(3),x2(3)
       REAL*8     stop_alt
       INTEGER*4  hemi_flag
       REAL*8     B(3),Bl,B0,B1,B3
       REAL*8     dsreb
       
       REAL*8     XFOOT(3),BFOOT(3),BFOOTMAG

       INTEGER*4  I,J
       REAL*8     Lb
C      IFOUND is a dummy loop result variable
       integer*4  IFOUND
C
C
       CALL GEO_SM(xx0,xx)
       rr = SQRT(xx(1)*xx(1)+xx(2)*xx(2)+xx(3)*xx(3))
       tt = ACOS(xx(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
c       write(6,*)'L bete ',Lb
C
       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) THEN
          goto 999
       ENDIF

       !geo_gdz provides lat/lon/alt at x2
       call geo_gdz(xx0(1),xx0(2),xx0(3),XFOOT(2),XFOOT(3),XFOOT(1))
       if (XFOOT(1).LE.stop_alt) then
            goto 999 ! fail altitude of starting point to low
       endif


       XFOOT(1) = baddata
       XFOOT(2) = baddata
       XFOOT(3) = baddata
       BFOOT(1) = baddata
       BFOOT(2) = baddata
       BFOOT(3) = baddata
       BFOOTMAG = baddata

       dsreb = Lb/(Nreb*1.d0) ! step size

        ! introduced to prevent overshooting
        ! where one may end up in the opposite hemisphere
        ! Kellerman / Tenfjord
       if (dsreb.GT.1) THEN
        dsreb = 1
       ENDIF

C calcul du sens du depart 
C sksyst call takes southward step
       CALL sksyst(-dsreb,xx0,x1,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          goto 999
       ENDIF
       B1 = Bl
C sksyst call takes northward step
       CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          goto 999
       ENDIF
       B3 = Bl
C

c     dsreb is currently pointing north
       if (hemi_flag.eq.-1) then
          dsreb = -dsreb ! point dsreb south
       endif
       if (hemi_flag.eq.0) then
          IF (B1.GT.B3) THEN
             dsreb = -dsreb     ! point dsreb toward larger field (same hemi)
          ENDIF
       endif
       if (hemi_flag.eq.2) then
          IF (B1.LT.B3) THEN
             dsreb = -dsreb     ! point dsreb toward smaller field (toward other hemi)
          ENDIF
       endif

C
C calcul de la ligne de champ
C
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
C
       Bl = B0 ! reset to starting value
15     continue ! prepare to do loop
       IFOUND = 0
       DO J = 1,500
         B1 = Bl ! store for comparison after step
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            goto 999
         ENDIF
c
c test for completion
c need to check: alt < stop_alt, and moving to lower alt (higher B)
         ! geo_gdz provides lat/lon/alt at x2
         call geo_gdz(x2(1),x2(2),x2(3),XFOOT(2),XFOOT(3),XFOOT(1))
         if (XFOOT(1).LE.stop_alt) then
            IFOUND = 1
            goto 20 ! done with loop
         endif
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
20     CONTINUE
C

       if (IFOUND.EQ.1) then
          ! footpoint is between x1 and x2
          if (abs(XFOOT(1)-stop_alt).le.1.0) then
             !get B field at x2
             call champ(x2,BFOOT,BFOOTMAG,Ifail)
             if (Ifail.LT.0) then
                goto 999
             endif
             return
          else  ! try loop again with smaller step
             dsreb = dsreb/100.0
             goto 15
          endif
       else
          goto 999
       endif
C
100    CONTINUE
C
       RETURN
 999   CONTINUE  ! fail!
       XFOOT(1) = baddata
       XFOOT(2) = baddata
       XFOOT(3) = baddata
       BFOOT(1) = baddata
       BFOOT(2) = baddata
       BFOOT(3) = baddata
       BFOOTMAG = baddata
       RETURN
       END
C

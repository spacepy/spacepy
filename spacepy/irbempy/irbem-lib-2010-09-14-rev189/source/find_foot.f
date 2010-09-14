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
c
C Routine to find foot point of field line at specified altitude and hemi
C finds foot point at specified altitude to within 1 km
C

      REAL*4 FUNCTION find_foot_point(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  convert the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine find_foot_point1: 15 arguments
      call find_foot_point1(%VAL(argv(1)), %VAL(argv(2)),
     + %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)), %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)))

      find_foot_point = 9.9

      RETURN
      END
c

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
	INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
	INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
        REAL*8     yearsat,dec_year,a2000_ut
	REAL*8     psi,mlon,tilt
	REAL*8     xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3),xGEO(3)
	real*8     alti,lati,longi
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
	real*8     BxGEO(3)
c
c Declare output variables	
        REAL*8     XFOOT(3),BFOOT(3),BFOOTMAG
C
        COMMON/GENER/ERA,AQUAD,BQUAD
        COMMON /dip_ang/tilt
	COMMON /magmod/k_ext,k_l,kint
        COMMON /drivers/density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
     &        ,G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04,
     &         W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
	COMMON /index/activ
	COMMON /a2000_time/a2000_ut,a2000_iyear,a2000_imonth,a2000_iday
        DATA  xSUN /1.d0,0.d0,0.d0/
C

        if ((stop_alt.lt.0).or.(stop_alt.ge.6378.0*500.0)) then
           goto 999 ! fail, stop_alt out of range 0 to 500 Re
        endif
        iyear=1800
	k_ext=kext
	kint=options(5)
	IF (kint .lt. 0) THEN
	   kint=0
	   WRITE(6,*)
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)'Invalid internal field specification'
	   WRITE(6,*)'Selecting IGRF'
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)
	ENDIF
	if (kint .gt. 3) THEN
	   kint=0
	   WRITE(6,*)
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)'Invalid internal field specification'
	   WRITE(6,*)'Selecting IGRF'
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)
	ENDIF
	IF (kext .lt. 0) THEN
	   k_ext=5
	   WRITE(6,*)
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)'Invalid external field specification'
	   WRITE(6,*)'Selecting Olson-Pfitzer (quiet)'
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)
	ENDIF
	if (kext .gt. 12) THEN
	   k_ext=5
	   WRITE(6,*)
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)'Invalid external field specification'
	   WRITE(6,*)'Selecting Olson-Pfitzer (quiet)'
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)
	ENDIF
c
        CALL INITIZE
	if (kint .eq. 2) CALL JensenANDCain1960
	if (kint .eq. 3) CALL GSFC1266
        if (kint .le. 1) then
           if (options(2) .eq. 0) then	
	      if (iyearsat .ne. iyear) then
	         iyear=iyearsat
	         dec_year=iyear+0.5d0
	         CALL INIT_DTD(dec_year)
	      endif
	   else
	      if (iyearsat .ne. iyear .or.
     &        MOD(idoy*1.d0,options(2)*1.d0) .eq. 0) THEN
	         iyear=iyearsat
                 firstJanuary=JULDAY(iyear,01,01)
                 lastDecember=JULDAY(iyear,12,31)
                 currentdoy=(idoy/options(2))*options(2)
	         if (currentdoy .eq. 0) currentdoy=1
	         dec_year=iyear+currentdoy*1.d0/
     &             ((lastDecember-firstJanuary+1)*1.d0)
	         CALL INIT_DTD(dec_year)
              endif
	   endif
	endif
c
        CALL INIT_GSM(iyearsat,idoy,UT,psi)
        tilt = psi/(4.D0*ATAN(1.D0)/180.d0)
	if (sysaxes .EQ. 0) then
	    alti=xIN1
	    lati=xIN2
	    longi=xIN3
	endif
	if (sysaxes .EQ. 1) then
	    xGEO(1)=xIN1
	    xGEO(2)=xIN2
	    xGEO(3)=xIN3
	    CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
        endif
	if (sysaxes .EQ. 2) then
	    xGSM(1)=xIN1
	    xGSM(2)=xIN2
	    xGSM(3)=xIN3
	    CALL GSM_GEO(xGSM,xGEO)
	    CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	endif
	if (sysaxes .EQ. 3) then
	    xGSE(1)=xIN1
	    xGSE(2)=xIN2
	    xGSE(3)=xIN3
	    CALL GSE_GEO(xGSE,xGEO)
	    CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	endif
	if (sysaxes .EQ. 4) then
	    xSM(1)=xIN1
	    xSM(2)=xIN2
	    xSM(3)=xIN3
	    CALL SM_GEO(xSM,xGEO)
	    CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	endif
	if (sysaxes .EQ. 5) then
	    xGEI(1)=xIN1
	    xGEI(2)=xIN2
	    xGEI(3)=xIN3
	    CALL GEI_GEO(xGEI,xGEO)
	    CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	endif
	if (sysaxes .EQ. 6) then
	    xMAG(1)=xIN1
	    xMAG(2)=xIN2
	    xMAG(3)=xIN3
	    CALL MAG_GEO(xMAG,xGEO)
	    CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	endif
	if (sysaxes .EQ. 7) then
	    xMAG(1)=xIN1
	    xMAG(2)=xIN2
	    xMAG(3)=xIN3
	    CALL SPH_CAR(xMAG(1),xMAG(2),xMAG(3),xGEO)
	    CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	endif
	if (sysaxes .EQ. 8) then
	    xMAG(1)=xIN1
	    lati=xIN2
	    longi=xIN3
	    CALL RLL_GDZ(xMAG(1),lati,longi,alti)
	endif
c	   	   
c make inputs according to magn. field model chosen
c
        if (kext .eq. 1) then
c Input for MEAD
	    if (maginput(1).le.3.d0) Activ=1
	    if (maginput(1).gt.3.d0 .and. 
     &      maginput(1).lt.20.d0) Activ=2
	    if (maginput(1).ge.20.d0 .and. 
     &      maginput(1).lt.30.d0) Activ=3
	    if (maginput(1).ge.30.d0) Activ=4
c
	    if (maginput(1).lt.0.d0 .or. 
     &      maginput(1).gt.90.d0) then
               goto 999
	    endif
	endif
        if (kext .eq. 2) then
c Input for TSYG87s
	    if (maginput(1).lt.7.d0) Activ=1
	    if (maginput(1).ge.7.d0 .and. 
     &      maginput(1).lt.17.d0) Activ=2
	    if (maginput(1).ge.17.d0 .and. 
     &      maginput(1).lt.20.d0) Activ=3
	    if (maginput(1).ge.20.d0 .and. 
     &      maginput(1).lt.27.d0) Activ=4
	    if (maginput(1).ge.27.d0 .and. 
     &      maginput(1).lt.37.d0) Activ=5
	    if (maginput(1).ge.37.d0 .and. 
     &      maginput(1).lt.47.d0) Activ=6
	    if (maginput(1).ge.47.d0) Activ=7
	    if (maginput(1).ge.53.d0) Activ=8
c
	    if (maginput(1).lt.0.d0 .or. 
     &      maginput(1).gt.90.d0) then
               goto 999
	    endif
	endif
        if (kext .eq. 3) then
c Input for TSYG87l
	    if (maginput(1).lt.7.d0) Activ=1
	    if (maginput(1).ge.7.d0 .and. 
     &      maginput(1).lt.17.d0) Activ=2
	    if (maginput(1).ge.17.d0 .and. 
     &      maginput(1).lt.27.d0) Activ=3
	    if (maginput(1).ge.27.d0 .and. 
     &      maginput(1).lt.37.d0) Activ=4
	    if (maginput(1).ge.37.d0 .and. 
     &      maginput(1).lt.47.d0) Activ=5
	    if (maginput(1).ge.47.d0) Activ=6
c
	    if (maginput(1).lt.0.d0 .or. 
     &      maginput(1).gt.90.d0) then
               goto 999
	    endif
        endif
        if (kext .eq. 4) then
c Input for Tsy89
	    if (maginput(1).lt.7.d0) Activ=1
	    if (maginput(1).ge.7.d0 .and. 
     &      maginput(1).lt.17.d0) Activ=2
	    if (maginput(1).ge.17.d0 .and. 
     &      maginput(1).lt.27.d0) Activ=3
	    if (maginput(1).ge.27.d0 .and. 
     &      maginput(1).lt.37.d0) Activ=4
	    if (maginput(1).ge.37.d0 .and. 
     &      maginput(1).lt.47.d0) Activ=5
	    if (maginput(1).ge.47.d0 .and. 
     &      maginput(1).lt.57.d0) Activ=6
	    if (maginput(1).ge.57.d0) Activ=7
c
	    if (maginput(1).lt.0.d0 .or. 
     &      maginput(1).gt.90.d0) then
               goto 999
	    endif
	endif
        if (kext .eq. 6) then
c Input for OP dyn
            density=maginput(3)
	    speed=maginput(4)
	    dst_nt=maginput(2)
c
	    if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
               goto 999
	    endif
	    if (density.lt.5.d0 .or. density.gt.50.d0) then
               goto 999
	    endif
	    if (speed.lt.300.d0 .or. speed.gt.500.d0) then
               goto 999
	    endif
	endif
        if (kext .eq. 7) then
c Input for Tsy96
	    dst_nt=maginput(2)
	    Pdyn_nPa=maginput(5)
	    ByIMF_nt=maginput(6)
	    BzIMF_nt=maginput(7)
c
	    if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
               goto 999
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
               goto 999
	    endif
	    if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
               goto 999
	    endif
	    if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
               goto 999
	    endif
	endif
        if (kext .eq. 8) then
c Input for Ostapenko97
	    dst_nt=maginput(2)
	    Pdyn_nPa=maginput(5)
	    BzIMF_nt=maginput(7)
	    fkp=maginput(1)*1.d0/10.d0
	endif
        if (kext .eq. 9) then
c Input for Tsy01
	    dst_nt=maginput(2)
	    Pdyn_nPa=maginput(5)
	    ByIMF_nt=maginput(6)
	    BzIMF_nt=maginput(7)
	    G1_tsy01=maginput(8)
	    G2_tsy01=maginput(9)
c
	    if (dst_nt.lt.-50.d0 .or. dst_nt.gt.20.d0) then
               goto 999
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
               goto 999
	    endif
	    if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
               goto 999
	    endif
	    if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
               goto 999
	    endif
	    if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
               goto 999
	    endif
	    if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
               goto 999
	    endif
	endif
        if (kext .eq. 10) then
c Input for Tsy01 storm
	    dst_nt=maginput(2)
	    Pdyn_nPa=maginput(5)
	    ByIMF_nt=maginput(6)
	    BzIMF_nt=maginput(7)
	    G2_tsy01=maginput(9)
	    G3_tsy01=maginput(10)
	endif
c	   
        if (kext .eq. 11) then
c Input for Tsy04 storm
	    dst_nt=maginput(2)
	    Pdyn_nPa=maginput(5)
	    ByIMF_nt=maginput(6)
	    BzIMF_nt=maginput(7)
	    W1_tsy04=maginput(11)
	    W2_tsy04=maginput(12)
	    W3_tsy04=maginput(13)
	    W4_tsy04=maginput(14)
	    W5_tsy04=maginput(15)
	    W6_tsy04=maginput(16)
	endif
c
        if (kext .eq. 12) then
c Input for Alexeev 2000
            a2000_iyear=iyearsat
	    firstJanuary=JULDAY(a2000_iyear,01,01)
	    currentdoy=firstJanuary+idoy-1
	    CALL CALDAT(currentdoy,a2000_iyear,
     &      a2000_imonth,a2000_iday)
	    a2000_ut=UT
            density=maginput(3)
	    speed=maginput(4)
	    dst_nt=maginput(2)
	    BzIMF_nt=maginput(7)
	    Al=maginput(17)
	endif
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
       INCLUDE 'variables.inc'
C
       INTEGER*4  Nreb
       PARAMETER (Nreb = 50)
C
       INTEGER*4  k_ext,k_l,kint,Ifail
       INTEGER*4  Nrebmax
       REAL*8     rr
       REAL*8     xx0(3),xx(3),x1(3),x2(3)
       REAL*8     xmin(3)
       REAL*8     lati,longi,alti
       REAL*8     stop_alt
       INTEGER*4  hemi_flag
       REAL*8     B(3),Bl,B0,B1,B3
       REAL*8     dsreb
       
       REAL*8     XFOOT(3),BFOOT(3),BFOOTMAG

       INTEGER*4  I,J,ind
       REAL*8     Lb
       REAL*8     leI0
C
       REAL*8     pi,rad
       REAL*8     tt
       REAL*8     cste
C
       REAL*8     Bposit,Bmir
       REAL*8     sn2,sn2max,alpha
       integer*4  IFOUND ! dummy loop result variable
C
       COMMON /magmod/k_ext,k_l,kint
C
C
       pi = 4.D0*ATAN(1.D0)
       rad = pi/180.D0
C
       Nrebmax = 20*Nreb
C
       leI0 = 0.D0

C
       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
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
       XFOOT(1) = alti
       XFOOT(2) = lati
       XFOOT(3) = longi
       BFOOT(1) = B(1)
       BFOOT(2) = B(2)
       BFOOT(3) = B(3)
       BFOOTMAG = B0

       dsreb = Lb/(Nreb*1.d0) ! step size
C
C calcul du sens du depart 
C
       CALL sksyst(-dsreb,xx0,x1,Bl,Ifail) ! southward step
       IF (Ifail.LT.0) THEN
          goto 999
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xx0,x2,Bl,Ifail)! northward step
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
       ind=0
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
C
       xmin(1)=x1(1)
       xmin(2)=x1(2)
       xmin(3)=x1(3)
       Bl = B0 ! reset to starting value
15     continue ! prepare to do loop
       IFOUND = 0
       DO J = 1,100
         B1 = Bl ! store for comparison after step
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            goto 999
         ENDIF
c
c test for completion
c need to check: alt < stop_alt, and moving to lower alt (higher B)
         call geo_gdz(x2(1),x2(2),x2(3),XFOOT(2),XFOOT(3),XFOOT(1)) ! provides lat/lon/alt at x2
         if ((B1.LT.Bl).AND.(XFOOT(1).LE.stop_alt)) then
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
             call champ(x2,BFOOT,BFOOTMAG,Ifail) ! get field at x2
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

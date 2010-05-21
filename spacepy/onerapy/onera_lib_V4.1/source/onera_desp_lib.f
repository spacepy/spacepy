C----------------------------------------------------------------------------- 
C Wrappers and procedures for ONERA_DESP_LIB
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION make_lstar(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine make_Lstar: 17 arguments
      call make_lstar1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)), %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     + %VAL(argv(16)), %VAL(argv(17)))

      make_lstar = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE make_lstar1(ntime,kext,options,sysaxes,iyearsat,
     &  idoy,UT,xIN1,xIN2,xIN3,maginput,Lm,Lstar,BLOCAL,BMIN,XJ,MLT)
c
	IMPLICIT NONE
	INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    ntime_max,kext,k_ext,k_l,options(5)
	PARAMETER (ntime_max=100000)
        INTEGER*4    ntime,sysaxes
	INTEGER*4    iyearsat(ntime_max)
	integer*4    idoy(ntime_max)
	real*8     UT(ntime_max)
	real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
	real*8     maginput(25,ntime_max)
c                      1: Kp
c                      2: Dst
c                      3: dens
c                      4: velo
c                      5: Pdyn
c                      6: ByIMF
c                      7: BzIMF
c                      8: G1
c                      9: G2
c                     10: G3
c
c Declare internal variables
	INTEGER*4    isat,iyear,kint
        INTEGER*4    Ndays,activ,t_resol,r_resol,Ilflag,Ilflag_old
	INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
	INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
        REAL*8     yearsat,dec_year,a2000_ut
	REAL*8     psi,mlon,tilt
	REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
c
c Declare output variables	
        REAL*8     BLOCAL(ntime_max),BMIN(ntime_max),XJ(ntime_max)
	REAL*8     MLT(ntime_max)
        REAL*8     Lm(ntime_max),Lstar(ntime_max)
C
        COMMON/GENER/ERA,AQUAD,BQUAD
        COMMON /dip_ang/tilt
	COMMON /magmod/k_ext,k_l,kint
        COMMON /drivers/density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
     &        ,G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04,
     &         W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
	COMMON /index/activ
        COMMON /flag_L/Ilflag
	COMMON /a2000_time/a2000_ut,a2000_iyear,a2000_imonth,a2000_iday
        DATA  xSUN /1.d0,0.d0,0.d0/
C
	Ilflag=0
	Ilflag_old=Ilflag
        iyear=1800
	k_ext=kext
	if (options(3).lt.0 .or. options(3).gt.9) options(3)=0
	t_resol=options(3)+1
	r_resol=options(4)+1
	k_l=options(1)
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
c
        CALL INITIZE
	if (kint .eq. 2) CALL JensenANDCain1960
	if (kint .eq. 3) CALL GSFC1266
        DO isat = 1,ntime
c	   write(6,*)real(isat)*100./real(ntime), '% done'
c
           if (kint .le. 1) then
              if (options(2) .eq. 0) then	
	        if (iyearsat(isat) .ne. iyear) then
	           iyear=iyearsat(isat)
	           dec_year=iyear+0.5d0
	           CALL INIT_DTD(dec_year)
	        endif
	      else
	        if (iyearsat(isat) .ne. iyear .or.
     &          MOD(idoy(isat)*1.d0,options(2)*1.d0) .eq. 0) THEN
	           iyear=iyearsat(isat)
		   firstJanuary=JULDAY(iyear,01,01)
		   lastDecember=JULDAY(iyear,12,31)
		   currentdoy=(idoy(isat)/options(2))*options(2)
		   if (currentdoy .eq. 0) currentdoy=1
	           dec_year=iyear+currentdoy*1.d0/
     &             ((lastDecember-firstJanuary+1)*1.d0)
	           CALL INIT_DTD(dec_year)
                endif
	      endif
	   endif
c
           CALL INIT_GSM(iyearsat(isat),idoy(isat),UT(isat),psi)
           tilt = psi/(4.D0*ATAN(1.D0)/180.d0)
	   if (sysaxes .EQ. 0) then
	       alti=xIN1(isat)
	       lati=xIN2(isat)
	       longi=xIN3(isat)
	   endif
	   if (sysaxes .EQ. 1) then
	       xGEO(1)=xIN1(isat)
	       xGEO(2)=xIN2(isat)
	       xGEO(3)=xIN3(isat)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 2) then
	       xGSM(1)=xIN1(isat)
	       xGSM(2)=xIN2(isat)
	       xGSM(3)=xIN3(isat)
	       CALL GSM_GEO(xGSM,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 3) then
	       xGSE(1)=xIN1(isat)
	       xGSE(2)=xIN2(isat)
	       xGSE(3)=xIN3(isat)
	       CALL GSE_GEO(xGSE,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 4) then
	       xSM(1)=xIN1(isat)
	       xSM(2)=xIN2(isat)
	       xSM(3)=xIN3(isat)
	       CALL SM_GEO(xSM,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 5) then
	       xGEI(1)=xIN1(isat)
	       xGEI(2)=xIN2(isat)
	       xGEI(3)=xIN3(isat)
	       CALL GEI_GEO(xGEI,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 6) then
	       xMAG(1)=xIN1(isat)
	       xMAG(2)=xIN2(isat)
	       xMAG(3)=xIN3(isat)
	       CALL MAG_GEO(xMAG,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 7) then
	       xMAG(1)=xIN1(isat)
	       xMAG(2)=xIN2(isat)
	       xMAG(3)=xIN3(isat)
	       CALL SPH_CAR(xMAG(1),xMAG(2),xMAG(3),xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 8) then
	       xMAG(1)=xIN1(isat)
	       lati=xIN2(isat)
	       longi=xIN3(isat)
	       CALL RLL_GDZ(xMAG(1),lati,longi,alti)
	   endif
c	   
c make inputs according to magn. field model chosen
c
           if (kext .eq. 1) then
c Input for MEAD
	       if (maginput(1,isat).le.3.d0) Activ=1
	       if (maginput(1,isat).gt.3.d0 .and. 
     &         maginput(1,isat).lt.20.d0) Activ=2
	       if (maginput(1,isat).ge.20.d0 .and. 
     &         maginput(1,isat).lt.30.d0) Activ=3
	       if (maginput(1,isat).ge.30.d0) Activ=4
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 2) then
c Input for TSYG87s
	       if (maginput(1,isat).lt.7.d0) Activ=1
	       if (maginput(1,isat).ge.7.d0 .and. 
     &         maginput(1,isat).lt.17.d0) Activ=2
	       if (maginput(1,isat).ge.17.d0 .and. 
     &         maginput(1,isat).lt.20.d0) Activ=3
	       if (maginput(1,isat).ge.20.d0 .and. 
     &         maginput(1,isat).lt.27.d0) Activ=4
	       if (maginput(1,isat).ge.27.d0 .and. 
     &         maginput(1,isat).lt.37.d0) Activ=5
	       if (maginput(1,isat).ge.37.d0 .and. 
     &         maginput(1,isat).lt.47.d0) Activ=6
	       if (maginput(1,isat).ge.47.d0) Activ=7
	       if (maginput(1,isat).ge.53.d0) Activ=8
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 3) then
c Input for TSYG87l
	       if (maginput(1,isat).lt.7.d0) Activ=1
	       if (maginput(1,isat).ge.7.d0 .and. 
     &         maginput(1,isat).lt.17.d0) Activ=2
	       if (maginput(1,isat).ge.17.d0 .and. 
     &         maginput(1,isat).lt.27.d0) Activ=3
	       if (maginput(1,isat).ge.27.d0 .and. 
     &         maginput(1,isat).lt.37.d0) Activ=4
	       if (maginput(1,isat).ge.37.d0 .and. 
     &         maginput(1,isat).lt.47.d0) Activ=5
	       if (maginput(1,isat).ge.47.d0) Activ=6
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
           endif
           if (kext .eq. 4) then
c Input for Tsy89
	       if (maginput(1,isat).lt.7.d0) Activ=1
	       if (maginput(1,isat).ge.7.d0 .and. 
     &         maginput(1,isat).lt.17.d0) Activ=2
	       if (maginput(1,isat).ge.17.d0 .and. 
     &         maginput(1,isat).lt.27.d0) Activ=3
	       if (maginput(1,isat).ge.27.d0 .and. 
     &         maginput(1,isat).lt.37.d0) Activ=4
	       if (maginput(1,isat).ge.37.d0 .and. 
     &         maginput(1,isat).lt.47.d0) Activ=5
	       if (maginput(1,isat).ge.47.d0 .and. 
     &         maginput(1,isat).lt.57.d0) Activ=6
	       if (maginput(1,isat).ge.57.d0) Activ=7
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 6) then
c Input for OP dyn
               density=maginput(3,isat)
	       speed=maginput(4,isat)
	       dst_nt=maginput(2,isat)
c
	       if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (density.lt.5.d0 .or. density.gt.50.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (speed.lt.300.d0 .or. speed.gt.500.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 7) then
c Input for Tsy96
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
c
	       if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 8) then
c Input for Ostapenko97
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       BzIMF_nt=maginput(7,isat)
	       fkp=maginput(1,isat)*1.d0/10.d0
	   endif
           if (kext .eq. 9) then
c Input for Tsy01
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
	       G1_tsy01=maginput(8,isat)
	       G2_tsy01=maginput(9,isat)
c
	       if (dst_nt.lt.-50.d0 .or. dst_nt.gt.20.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
	          Lm(isat)=baddata
		  Lstar(isat)=baddata
		  XJ(isat)=baddata
		  BLOCAL(isat)=baddata
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 10) then
c Input for Tsy01 storm
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
	       G2_tsy01=maginput(9,isat)
	       G3_tsy01=maginput(10,isat)
	   endif
c
           if (kext .eq. 11) then
c Input for Tsy04 storm
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
	       W1_tsy04=maginput(11,isat)
	       W2_tsy04=maginput(12,isat)
	       W3_tsy04=maginput(13,isat)
	       W4_tsy04=maginput(14,isat)
	       W5_tsy04=maginput(15,isat)
	       W6_tsy04=maginput(16,isat)
	   endif
c
           if (kext .eq. 12) then
c Input for Alexeev 2000
               a2000_iyear=iyearsat(isat)
	       firstJanuary=JULDAY(a2000_iyear,01,01)
	       currentdoy=firstJanuary+idoy(isat)-1
	       CALL CALDAT(currentdoy,a2000_iyear,
     &         a2000_imonth,a2000_iday)
	       a2000_ut=UT(isat)
               density=maginput(3,isat)
	       speed=maginput(4,isat)
	       dst_nt=maginput(2,isat)
	       BzIMF_nt=maginput(7,isat)
	       Al=maginput(17,isat)
	   endif
c
           CALL calcul_Lstar(t_resol,r_resol,lati,longi,alti
     &     ,Lm(isat),Lstar(isat),XJ(isat),BLOCAL(isat),BMIN(isat))
           if (Ilflag_old .eq.1 .and. Lstar(isat).eq. Baddata) then
	      Ilflag=0
              CALL calcul_Lstar(t_resol,r_resol,lati,longi,alti
     &        ,Lm(isat),Lstar(isat),XJ(isat),BLOCAL(isat),BMIN(isat))
	   endif
	   Ilflag_old=Ilflag
99         continue
           CALL GDZ_GEO(lati,longi,alti
     &     ,xGEO(1),xGEO(2),xGEO(3))
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon1)
           CALL GSM_GEO(xSUN,xGEO)
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon)
           MLT(isat) = (Mlon1 - Mlon)/15.d0 + 12.d0
           IF (MLT(isat).GE.24.d0) MLT(isat) = MLT(isat) - 24.d0
           IF (MLT(isat).LT.0.d0) MLT(isat) = MLT(isat) + 24.d0
	ENDDO
	END
C----------------------------------------------------------------------------- 
C Wrapper and procedure
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION make_lstar_shell_splitting(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
      call make_lstar_shell_splitting1(%VAL(argv(1)), %VAL(argv(2)),
     + %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)), %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     + %VAL(argv(16)), %VAL(argv(17)), %VAL(argv(18)), %VAL(argv(19)))

      make_lstar_shell_splitting = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE make_lstar_shell_splitting1(ntime,Nipa,kext,options,
     & sysaxes,iyearsat,idoy,UT,xIN1,xIN2,xIN3,
     &  alpha,maginput,Lm,Lstar,BLOCAL,BMIN,XJ,MLT)
c
	IMPLICIT NONE
	INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    ntime_max,kext,k_ext,k_l,options(5),Nalp,Nipa
	PARAMETER (ntime_max=100000)
	PARAMETER (Nalp=25)
        INTEGER*4    ntime,sysaxes
	INTEGER*4    iyearsat(ntime_max)
	integer*4    idoy(ntime_max)
	real*8     UT(ntime_max)
	real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
	real*8     alpha(Nalp)
	real*8     maginput(25,ntime_max)
c                      1: Kp
c                      2: Dst
c                      3: dens
c                      4: velo
c                      5: Pdyn
c                      6: ByIMF
c                      7: BzIMF
c                      8: G1
c                      9: G2
c                     10: G3
c
c
c Declare internal variables
	INTEGER*4    isat,iyear,IPA,kint
        INTEGER*4    Ndays,activ,Ilflag,t_resol,r_resol
	INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
	INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
        REAL*8     yearsat,dec_year,a2000_ut
	REAL*8     psi,mlon,tilt,BL,BMIR,Bmin_tmp
	REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi,alti_IPA,lati_IPA,longi_IPA
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
c
c Declare output variables	
        REAL*8     BLOCAL(ntime_max,Nalp),BMIN(ntime_max)
	REAL*8     XJ(ntime_max,Nalp),MLT(ntime_max)
        REAL*8     Lm(ntime_max,Nalp),Lstar(ntime_max,Nalp)
C
        COMMON/GENER/ERA,AQUAD,BQUAD
        COMMON /dip_ang/tilt
	COMMON /magmod/k_ext,k_l,kint
        COMMON /drivers/density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
     &        ,G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04,
     &         W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
	COMMON /index/activ
        COMMON /flag_L/Ilflag
	COMMON /a2000_time/a2000_ut,a2000_iyear,a2000_imonth,a2000_iday
        DATA  xSUN /1.d0,0.d0,0.d0/
C
	Ilflag=0
        iyear=1800
	k_ext=kext
	if (options(3).lt.0 .or. options(3).gt.9) options(3)=0
	t_resol=options(3)+1
	r_resol=options(4)+1
	k_l=options(1)
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
        DO isat = 1,ntime
c	   write(6,*)real(isat)*100./real(ntime), '% done'
c	
           if (kint .le. 1) then
              if (options(2) .eq. 0) then	
	        if (iyearsat(isat) .ne. iyear) then
	           iyear=iyearsat(isat)
	           dec_year=iyear+0.5d0
	           CALL INIT_DTD(dec_year)
	        endif
	      else
	        if (iyearsat(isat) .ne. iyear .or.
     &          MOD(idoy(isat)*1.d0,options(2)*1.d0) .eq. 0) THEN
	           iyear=iyearsat(isat)
		   firstJanuary=JULDAY(iyear,01,01)
		   lastDecember=JULDAY(iyear,12,31)
		   currentdoy=(idoy(isat)/options(2))*options(2)
		   if (currentdoy .eq. 0) currentdoy=1
	           dec_year=iyear+currentdoy*1.d0/
     &             ((lastDecember-firstJanuary+1)*1.d0)
	           CALL INIT_DTD(dec_year)
                endif
	      endif
	   endif
c
           CALL INIT_GSM(iyearsat(isat),idoy(isat),UT(isat),psi)
           tilt = psi/(4.D0*ATAN(1.D0)/180.d0)
	   if (sysaxes .EQ. 0) then
	       alti=xIN1(isat)
	       lati=xIN2(isat)
	       longi=xIN3(isat)
	   endif
	   if (sysaxes .EQ. 1) then
	       xGEO(1)=xIN1(isat)
	       xGEO(2)=xIN2(isat)
	       xGEO(3)=xIN3(isat)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 2) then
	       xGSM(1)=xIN1(isat)
	       xGSM(2)=xIN2(isat)
	       xGSM(3)=xIN3(isat)
	       CALL GSM_GEO(xGSM,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 3) then
	       xGSE(1)=xIN1(isat)
	       xGSE(2)=xIN2(isat)
	       xGSE(3)=xIN3(isat)
	       CALL GSE_GEO(xGSE,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 4) then
	       xSM(1)=xIN1(isat)
	       xSM(2)=xIN2(isat)
	       xSM(3)=xIN3(isat)
	       CALL SM_GEO(xSM,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 5) then
	       xGEI(1)=xIN1(isat)
	       xGEI(2)=xIN2(isat)
	       xGEI(3)=xIN3(isat)
	       CALL GEI_GEO(xGEI,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 6) then
	       xMAG(1)=xIN1(isat)
	       xMAG(2)=xIN2(isat)
	       xMAG(3)=xIN3(isat)
	       CALL MAG_GEO(xMAG,xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 7) then
	       xMAG(1)=xIN1(isat)
	       xMAG(2)=xIN2(isat)
	       xMAG(3)=xIN3(isat)
	       CALL SPH_CAR(xMAG(1),xMAG(2),xMAG(3),xGEO)
	       CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
	   endif
	   if (sysaxes .EQ. 8) then
	       xMAG(1)=xIN1(isat)
	       lati=xIN2(isat)
	       longi=xIN3(isat)
	       CALL RLL_GDZ(xMAG(1),lati,longi,alti)
	   endif
c	   
c make inputs according to magn. field model chosen
c
           if (kext .eq. 1) then
c Input for MEAD
	       if (maginput(1,isat).le.3.d0) Activ=1
	       if (maginput(1,isat).gt.3.d0 .and. 
     &         maginput(1,isat).lt.20.d0) Activ=2
	       if (maginput(1,isat).ge.20.d0 .and. 
     &         maginput(1,isat).lt.30.d0) Activ=3
	       if (maginput(1,isat).ge.30.d0) Activ=4
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 2) then
c Input for TSYG87s
	       if (maginput(1,isat).lt.7.d0) Activ=1
	       if (maginput(1,isat).ge.7.d0 .and. 
     &         maginput(1,isat).lt.17.d0) Activ=2
	       if (maginput(1,isat).ge.17.d0 .and. 
     &         maginput(1,isat).lt.20.d0) Activ=3
	       if (maginput(1,isat).ge.20.d0 .and. 
     &         maginput(1,isat).lt.27.d0) Activ=4
	       if (maginput(1,isat).ge.27.d0 .and. 
     &         maginput(1,isat).lt.37.d0) Activ=5
	       if (maginput(1,isat).ge.37.d0 .and. 
     &         maginput(1,isat).lt.47.d0) Activ=6
	       if (maginput(1,isat).ge.47.d0) Activ=7
	       if (maginput(1,isat).ge.53.d0) Activ=8
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 3) then
c Input for TSYG87l
	       if (maginput(1,isat).lt.7.d0) Activ=1
	       if (maginput(1,isat).ge.7.d0 .and. 
     &         maginput(1,isat).lt.17.d0) Activ=2
	       if (maginput(1,isat).ge.17.d0 .and. 
     &         maginput(1,isat).lt.27.d0) Activ=3
	       if (maginput(1,isat).ge.27.d0 .and. 
     &         maginput(1,isat).lt.37.d0) Activ=4
	       if (maginput(1,isat).ge.37.d0 .and. 
     &         maginput(1,isat).lt.47.d0) Activ=5
	       if (maginput(1,isat).ge.47.d0) Activ=6
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
           endif
           if (kext .eq. 4) then
c Input for Tsy89
	       if (maginput(1,isat).lt.7.d0) Activ=1
	       if (maginput(1,isat).ge.7.d0 .and. 
     &         maginput(1,isat).lt.17.d0) Activ=2
	       if (maginput(1,isat).ge.17.d0 .and. 
     &         maginput(1,isat).lt.27.d0) Activ=3
	       if (maginput(1,isat).ge.27.d0 .and. 
     &         maginput(1,isat).lt.37.d0) Activ=4
	       if (maginput(1,isat).ge.37.d0 .and. 
     &         maginput(1,isat).lt.47.d0) Activ=5
	       if (maginput(1,isat).ge.47.d0 .and. 
     &         maginput(1,isat).lt.57.d0) Activ=6
	       if (maginput(1,isat).ge.57.d0) Activ=7
c
	       if (maginput(1,isat).lt.0.d0 .or. 
     &         maginput(1,isat).gt.90.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	    endif
           if (kext .eq. 6) then
c Input for OP dyn
               density=maginput(3,isat)
	       speed=maginput(4,isat)
	       dst_nt=maginput(2,isat)
c
	       if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (density.lt.5.d0 .or. density.gt.50.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (speed.lt.300.d0 .or. speed.gt.500.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 7) then
c Input for Tsy96
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
c
	       if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 8) then
c Input for Ostapenko97
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       BzIMF_nt=maginput(7,isat)
	       fkp=maginput(1,isat)*1.d0/10.d0
	   endif
           if (kext .eq. 9) then
c Input for Tsy01
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
	       G1_tsy01=maginput(8,isat)
	       G2_tsy01=maginput(9,isat)
c
	       if (dst_nt.lt.-50.d0 .or. dst_nt.gt.20.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	       if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
	          DO IPA=1,Nipa
		     Lm(isat,IPA)=baddata
		     Lstar(isat,IPA)=baddata
		     XJ(isat,IPA)=baddata
		     BLOCAL(isat,IPA)=baddata
		  ENDDO
		  BMIN(isat)=baddata
		  GOTO 99
	       endif
	   endif
           if (kext .eq. 10) then
c Input for Tsy01 storm
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
	       G2_tsy01=maginput(9,isat)
	       G3_tsy01=maginput(10,isat)
	   endif
c
           if (kext .eq. 11) then
c Input for Tsy04 storm
	       dst_nt=maginput(2,isat)
	       Pdyn_nPa=maginput(5,isat)
	       ByIMF_nt=maginput(6,isat)
	       BzIMF_nt=maginput(7,isat)
	       W1_tsy04=maginput(11,isat)
	       W2_tsy04=maginput(12,isat)
	       W3_tsy04=maginput(13,isat)
	       W4_tsy04=maginput(14,isat)
	       W5_tsy04=maginput(15,isat)
	       W6_tsy04=maginput(16,isat)
	   endif
c
           if (kext .eq. 12) then
c Input for Alexeev 2000
               a2000_iyear=iyearsat(isat)
	       firstJanuary=JULDAY(a2000_iyear,01,01)
	       currentdoy=firstJanuary+idoy(isat)-1
	       CALL CALDAT(currentdoy,a2000_iyear,
     &         a2000_imonth,a2000_iday)
	       a2000_ut=UT(isat)
               density=maginput(3,isat)
	       speed=maginput(4,isat)
	       dst_nt=maginput(2,isat)
	       BzIMF_nt=maginput(7,isat)
	       Al=maginput(17,isat)
	   endif
c
c Compute Bmin assuming 90° PA at S/C
	   k_l=0
           IPA=1
           CALL calcul_Lstar(t_resol,r_resol,lati
     &         ,longi,alti
     &         ,Lm(isat,IPA),Lstar(isat,IPA),XJ(isat,IPA)
     &         ,BLOCAL(isat,IPA),BMIN(isat))
	   k_l=options(1)
           DO IPA=1,Nipa
             if (alpha(IPA).ne.90.0d0) then 
                CALL find_bm(
     &           lati,longi,alti,alpha(IPA),BL,BMIR,xGEO)
	        CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati_IPA,
     &            longi_IPA,alti_IPA)
             ELSE
	         Bmir=0.d0
	         lati_IPA=lati
		 longi_IPA=longi
		 alti_IPA=alti
	     ENDIF
	     IF (Bmir.NE.baddata) THEN
	       Ilflag=0
               CALL calcul_Lstar(t_resol,r_resol,lati_IPA
     &         ,longi_IPA,alti_IPA
     &         ,Lm(isat,IPA),Lstar(isat,IPA),XJ(isat,IPA)
     &         ,BLOCAL(isat,IPA),BMIN_tmp)
             ELSE
	       Lm(isat,IPA)=baddata
	       Lstar(isat,IPA)=baddata
	       XJ(isat,IPA)=baddata
	       BLOCAL(isat,IPA)=baddata
	     ENDIF
           ENDDO
99         continue
           CALL GDZ_GEO(lati,longi,alti
     &     ,xGEO(1),xGEO(2),xGEO(3))
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon1)
           CALL GSM_GEO(xSUN,xGEO)
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon)
           MLT(isat) = (Mlon1 - Mlon)/15.d0 + 12.d0
           IF (MLT(isat).GE.24.d0) MLT(isat) = MLT(isat) - 24.d0
           IF (MLT(isat).LT.0.d0) MLT(isat) = MLT(isat) + 24.d0
	ENDDO
	END
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION drift_shell(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine make_Lstar: 17 arguments
      call drift_shell1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)), %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     + %VAL(argv(16)), %VAL(argv(17)))

      drift_shell = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE drift_shell1(kext,options,sysaxes,iyearsat,idoy,UT,
     &  xIN1,xIN2,xIN3,maginput,Lm,Lstar,BLOCAL,BMIN,XJ,posit,ind)
c
	IMPLICIT NONE
	INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,options(5)
        INTEGER*4    sysaxes
	INTEGER*4    iyearsat
	integer*4    idoy
	real*8     UT
	real*8     xIN1,xIN2,xIN3
	real*8     maginput(25)
c
c Declare internal variables
	INTEGER*4    isat,iyear,kint
        INTEGER*4    Ndays,activ
	INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
	INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
        REAL*8     yearsat,dec_year,a2000_ut
	REAL*8     psi,mlon,tilt
	REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
c
c Declare output variables	
        INTEGER*4  ind(48)
        REAL*8     BLOCAL(1000,48),BMIN,XJ
        REAL*8     Lm,Lstar
	REAL*8     posit(3,1000,48)
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
        iyear=1800
	k_ext=kext
	k_l=options(1)
	kint=options(5)
c
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
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
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
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
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
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
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
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	endif
        if (kext .eq. 6) then
c Input for OP dyn
            density=maginput(3)
	    speed=maginput(4)
	    dst_nt=maginput(2)
c
	    if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (density.lt.5.d0 .or. density.gt.50.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (speed.lt.300.d0 .or. speed.gt.500.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
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
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
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
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
	       Lm=baddata
	       Lstar=baddata
	       XJ=baddata
	       BMIN=baddata
	       RETURN
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
        CALL trace_drift_shell(lati,longi,alti
     &     ,Lm,Lstar,XJ,BLOCAL,BMIN,
     &     posit,ind)
	END

C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION trace_field_line(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine make_Lstar: 17 arguments
      call trace_field_line1(%VAL(argv(1)), %VAL(argv(2)), 
     +%VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)), %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     + %VAL(argv(16)))

      trace_field_line = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE trace_field_line1(kext,options,sysaxes,iyearsat,idoy,
     &  UT,xIN1,xIN2,xIN3,maginput,Lm,BLOCAL,BMIN,XJ,posit,ind)
c
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
	real*8     maginput(25)
c
c Declare internal variables
	INTEGER*4    isat,iyear
        INTEGER*4    Ndays,activ,i,j
	INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
	INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
        REAL*8     yearsat,dec_year,a2000_ut
	REAL*8     psi,mlon,tilt
	REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
c
c Declare output variables	
        INTEGER*4  ind
        REAL*8     BLOCAL(1000),BMIN,XJ
        REAL*8     Lm
	REAL*8     posit(3,1000)
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
        do i=1,3
	   do j=1,1000
	     posit(i,j)=baddata
	   enddo
	enddo
c
        iyear=1800
	k_ext=kext
	k_l=options(1)
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
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
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
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
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
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
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
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	endif
        if (kext .eq. 6) then
c Input for OP dyn
            density=maginput(3)
	    speed=maginput(4)
	    dst_nt=maginput(2)
c
	    if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (density.lt.5.d0 .or. density.gt.50.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (speed.lt.300.d0 .or. speed.gt.500.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
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
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
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
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
	    endif
	    if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
	       Lm=baddata
	       XJ=baddata
	       BMIN=baddata
	       ind=0
	       RETURN
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
        CALL field_line_tracing(lati,longi,alti
     &     ,Lm,XJ,BLOCAL,BMIN,posit,ind)
	END

C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION find_mirror_point(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine make_Lstar: 19 arguments
      call find_mirror_point1(%VAL(argv(1)), %VAL(argv(2)), 
     + %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)), %VAL(argv(13)), %VAL(argv(14)))

      find_mirror_point = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE find_mirror_point1(kext,options,sysaxes,iyearsat,
     &  idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,BLOCAL,BMIR,xGEO)
c
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
	real*8     alpha
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
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
	real*8     BxGEO(3)
c
c Declare output variables	
        REAL*8     BLOCAL,xGEO(3),BMIR
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
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
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
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
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
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
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
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	endif
        if (kext .eq. 6) then
c Input for OP dyn
            density=maginput(3)
	    speed=maginput(4)
	    dst_nt=maginput(2)
c
	    if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (density.lt.5.d0 .or. density.gt.50.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (speed.lt.300.d0 .or. speed.gt.500.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
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
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
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
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
	    endif
	    if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
	       xGEO(1)=baddata
	       xGEO(2)=baddata
	       xGEO(3)=baddata
	       BLOCAL=baddata
	       BMIR=baddata
	       RETURN
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
        if (alpha.eq.90.0d0) then 
          Iint=2   
          CALL CHAMP(Iint,xGEO,BxGEO,Blocal,Ifail)
	  IF (Ifail.LT.0) THEN
	    Blocal=baddata
	    Bmir=baddata
	    xGEO(1)=baddata
	    xGEO(2)=baddata
	    xGEO(3)=baddata
	  ELSE
	    BMIR=Blocal
	  ENDIF
	  RETURN
	endif	   
c
        CALL find_bm(lati,longi,alti,alpha
     &     ,BLOCAL,BMIR,xGEO)
	END



C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION find_MAGequator(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  
      call find_MAGequator1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)))

      find_MAGequator = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE FIND_MAGEQUATOR1(kext,options,sysaxes,iyearsat
     &  ,idoy,UT,xIN1,xIN2,xIN3,maginput,BMIN,posit)
c
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
	real*8     maginput(25)
c
c Declare internal variables
	INTEGER*4    isat,iyear
        INTEGER*4    Ndays,activ
	INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
	INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
        REAL*8     yearsat,dec_year,a2000_ut
	REAL*8     psi,mlon,tilt
	REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
c
c Declare output variables	
        REAL*8     BMIN
	REAL*8     posit(3)
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
	k_l=0
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
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
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
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
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
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
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
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	endif
        if (kext .eq. 6) then
c Input for OP dyn
            density=maginput(3)
	    speed=maginput(4)
	    dst_nt=maginput(2)
c
	    if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (density.lt.5.d0 .or. density.gt.50.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (speed.lt.300.d0 .or. speed.gt.500.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
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
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
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
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
	    endif
	    if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
	       posit(1)=baddata
	       posit(2)=baddata
	       posit(3)=baddata
	       BMIN=baddata
	       RETURN
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
        CALL loc_equator(
     &     lati,longi,alti,BMIN,posit)
	END
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION GET_FIELD(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  
      call GET_FIELD1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     + %VAL(argv(12)))

      GET_FIELD = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE GET_FIELD1(kext,options,sysaxes,iyearsat,idoy,UT,
     &  xIN1,xIN2,xIN3,maginput,BxGEO,Bl)
c
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
	real*8     maginput(25)
c
c Declare internal variables
	INTEGER*4    isat,iyear
        INTEGER*4    Ndays,activ,Ifail
	INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
	INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
        REAL*8     yearsat,dec_year,a2000_ut
	REAL*8     psi,mlon,tilt
	REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi
        REAL*8     ERA,AQUAD,BQUAD
	real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
	real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
        real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
c
c Declare output variables	
        REAL*8     BxGEO(3),Bl
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
	k_l=0
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
	    CALL GDZ_GEO(lati,longi,alti,xGEO(1),xGEO(2),xGEO(3))
	endif
	if (sysaxes .EQ. 1) then
	    xGEO(1)=xIN1
	    xGEO(2)=xIN2
	    xGEO(3)=xIN3
        endif
	if (sysaxes .EQ. 2) then
	    xGSM(1)=xIN1
	    xGSM(2)=xIN2
	    xGSM(3)=xIN3
	    CALL GSM_GEO(xGSM,xGEO)
	endif
	if (sysaxes .EQ. 3) then
	    xGSE(1)=xIN1
	    xGSE(2)=xIN2
	    xGSE(3)=xIN3
	    CALL GSE_GEO(xGSE,xGEO)
	endif
	if (sysaxes .EQ. 4) then
	    xSM(1)=xIN1
	    xSM(2)=xIN2
	    xSM(3)=xIN3
	    CALL SM_GEO(xSM,xGEO)
	endif
	if (sysaxes .EQ. 5) then
	    xGEI(1)=xIN1
	    xGEI(2)=xIN2
	    xGEI(3)=xIN3
	    CALL GEI_GEO(xGEI,xGEO)
	endif
	if (sysaxes .EQ. 6) then
	    xMAG(1)=xIN1
	    xMAG(2)=xIN2
	    xMAG(3)=xIN3
	    CALL MAG_GEO(xMAG,xGEO)
	endif
	if (sysaxes .EQ. 7) then
	    xMAG(1)=xIN1
	    xMAG(2)=xIN2
	    xMAG(3)=xIN3
	    CALL SPH_CAR(xMAG(1),xMAG(2),xMAG(3),xGEO)
	endif
	if (sysaxes .EQ. 8) then
	    xMAG(1)=xIN1
	    lati=xIN2
	    longi=xIN3
	    CALL RLL_GDZ(xMAG(1),lati,longi,alti)
	    CALL GDZ_GEO(lati,longi,alti,xGEO(1),xGEO(2),xGEO(3))
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
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
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
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
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
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
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
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	endif
        if (kext .eq. 6) then
c Input for OP dyn
            density=maginput(3)
	    speed=maginput(4)
	    dst_nt=maginput(2)
c
	    if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (density.lt.5.d0 .or. density.gt.50.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (speed.lt.300.d0 .or. speed.gt.500.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
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
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
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
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
	    endif
	    if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
	       Bl=baddata
	       BxGEO(1)=baddata
	       BxGEO(2)=baddata
	       BxGEO(3)=baddata
	       RETURN
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
        CALL CHAMP(xGEO,BxGEO,Bl,Ifail)
	IF (Ifail.LT.0) THEN
	   BxGEO(1)=baddata
	   BxGEO(2)=baddata
	   BxGEO(3)=baddata
	   Bl=baddata
	ENDIF
	END
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION GET_MLT(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine make_Lstar, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  
      call GET_MLT1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      GET_MLT = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE GET_MLT1(iyr,idoy,UT,xGEO,MLT)
c
	IMPLICIT NONE
	INCLUDE 'variables.inc'
C
c declare inputs
	INTEGER*4    iyr
	integer*4    idoy
	real*8     UT
	real*8     xGEO
c
c Declare internal variables
        REAL*8     dyear
	REAL*8     psi,mlon
	REAL*8     xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xTMP(3)
c
c Declare output variables	
        REAL*8     MLT
C
        DATA  xSUN /1.d0,0.d0,0.d0/
C
	dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,ut,psi)
C
        CALL geo_mag(xGEO,xMAG)
        CALL car_sph(xMAG,rM,MLAT,Mlon1)
        CALL GSM_GEO(xSUN,xTMP)
        CALL geo_mag(xTMP,xMAG)
        CALL car_sph(xMAG,rM,MLAT,Mlon)
        MLT = (Mlon1 - Mlon)/15.d0 + 12.d0
        IF (MLT.GE.24.d0) MLT = MLT - 24.d0
        IF (MLT.LT.0.d0) MLT = MLT + 24.d0

        END

c Wrapper and proceduce to access many coordinate transformation form the 
c ONERA library
c
c =======================================================================
c GEO2GSM
c
c Routine to transform Cartesian GEO to cartesian GSM coordinates
c 
c  INPUTS: iyr = integer year
c          idoy = integer day of year
c          secs = UT in seconds
c          xGEO = 3D array of cartesian position in GEO (Re)
c
c OUTPUTS: psi: angle for GSM coordinate
c         xGSM: 3D array of cartesian position in GSM  (Re)
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'geo2gsm_', $                         ;The entry point
c           iyr,idoy,secs,psi,xGEO,xGSM, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c =======================================================================
c GSM2GEO
c
c Routine to transform Cartesian GSM to cartesian GEO coordinates
c 
c  INPUTS: iyr = integer year
c          idoy = integer day of year
c          secs = UT in seconds
c          xGSM = 3D array of cartesian position in GSM (Re)
c
c OUTPUTS: psi: angle for GSM coordinate
c         xGEO: 3D array of cartesian position in GEO (Re)
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'gsm2geo_', $                         ;The entry point
c           iyr,idoy,secs,psi,xGSM,xGEO, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c
c =======================================================================
c GDZ2GEO
c
c Routine to transform GEODEZIC coordinates to cartesian GEO coordinates
c 
c  INPUTS: lati = latitude (degres)
c          longi = longitude (degres)
c          alti = altitude (km)
c
c OUTPUTS: xx = xGEO (Re)
c          yy = yGEO (Re)
c          zz = zGEO (Re)
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'gdz2geo_', $                         ;The entry point
c           lati,longi,alti,xx,yy,zz, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c
c =======================================================================
c GEO2GDZ
c
c Routine to transform cartesian GEO coordinates to GEODEZIC coordinates 
c 
c INPUTS: xx = xGEO (Re)
c          yy = yGEO (Re)
c          zz = zGEO (Re)
c
c OUTPUTS: lati = latitude (degres)
c          longi = longitude (degres)
c          alti = altitude (km)
c
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'geo2gdz_', $                         ;The entry point
c           xx,yy,zz,lati,longi,alti, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c =======================================================================
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION coord_trans(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine coord_trans1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine coord_trans1: 7 arguments
      call coord_trans1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)))

      coord_trans = 9.9

      RETURN
      END
c

C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION coord_trans_vec(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine coord_trans1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine coord_trans_vec1: 8 arguments
      call coord_trans_vec1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)))

      coord_trans_vec = 9.9

      RETURN
      END
c

C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION geo2gsm(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call geo2gsm1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)))

      geo2gsm = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE geo2gsm1(iyr,idoy,secs,psi,xGEO,xGSM)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xGSM(3)
	
	dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_GSM(xGEO,xGSM)
        end
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION gsm2geo(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine gsm2geo1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine gsm2geo: 6 arguments
      call gsm2geo1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)))

      gsm2geo = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE gsm2geo1(iyr,idoy,secs,psi,xGSM,xGEO)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xGSM(3)
	
	
	dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GSM_GEO(xGSM,xGEO)
        end

C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION geo2gse(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call geo2gse1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      geo2gse = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE geo2gse1(iyr,idoy,secs,xGEO,xGSE)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xGSE(3)
	
	dyear=iyr+0.5d0
	psi=0.d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_GSE(xGEO,xGSE)
        end
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION gse2geo(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine gsm2geo1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine gsm2geo: 6 arguments
      call gse2geo1(%VAL(argv(1)), %VAL(argv(2)),
     * %VAL(argv(3)),  %VAL(argv(4)),  %VAL(argv(5)))

      gse2geo = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE gse2geo1(iyr,idoy,secs,xGSE,xGEO)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xGSE(3)
	
	
	dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GSE_GEO(xGSE,xGEO)
        end



C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION gdz2geo(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine gdz2geo1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine gdz2geo: 6 arguments
      call GDZ_GEO(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)))

      gdz2geo = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION geo2gdz(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo_gdz, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo_gdz: 6 arguments
      call GEO_GDZ(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)))

      geo2gdz = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION geo2gei(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call geo2gei1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      geo2gei = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE geo2gei1(iyr,idoy,secs,xGEO,xGEI)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xGEI(3)
	
	dyear=iyr+0.5d0
	psi=0.d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_GEI(xGEO,xGEI)
        end
	
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION gei2geo(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call gei2geo1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      gei2geo = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE gei2geo1(iyr,idoy,secs,xGEI,xGEO)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xGEI(3)
	
	dyear=iyr+0.5d0
	psi=0.d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEI_GEO(xGEI,xGEO)
        end
	
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION geo2sm(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call geo2sm1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      geo2sm = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE geo2sm1(iyr,idoy,secs,xGEO,xSM)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xSM(3)
	
	dyear=iyr+0.5d0
	psi=0.d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_SM(xGEO,xSM)
        end
	
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION sm2geo(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call sm2geo1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      sm2geo = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE sm2geo1(iyr,idoy,secs,xSM,xGEO)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGEO(3),xSM(3)
	
	dyear=iyr+0.5d0
	psi=0.d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL SM_GEO(xSM,xGEO)
        end
	
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION gsm2sm(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call gsm2sm1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      gsm2sm = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE gsm2sm1(iyr,idoy,secs,xGSM,xSM)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGSM(3),xSM(3)
	
	dyear=iyr+0.5d0
	psi=0.d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GSM_SM(xGSM,xSM)
        end
	
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION sm2gsm(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call sm2gsm1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      sm2gsm = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE sm2gsm1(iyr,idoy,secs,xSM,xGSM)
	INTEGER*4 iyr,idoy
	REAL*8    secs,psi,dyear
	REAL*8    xGSM(3),xSM(3)
	
	dyear=iyr+0.5d0
	psi=0.d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL SM_GSM(xSM,xGSM)
        end
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION geo2mag(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call geo2mag1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)))

      geo2mag = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE geo2mag1(iyr,xGEO,xMAG)
	INTEGER*4 iyr
	REAL*8    dyear
	REAL*8    xGEO(3),xMAG(3)
	
	dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL GEO_MAG(xGEO,xMAG)
        end
	
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION mag2geo(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call mag2geo1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)))

      mag2geo = 9.9

      RETURN
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE mag2geo1(iyr,xSM,xMAG)
	INTEGER*4 iyr
	REAL*8    dyear
	REAL*8    xGEO(3),xMAG(3)
	
	dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL MAG_GEO(xMAG,xGEO)
        end
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION sph2car(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call SPH_CAR(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)), 
     &  %VAL(argv(4)))

      sph2car = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION car2sph(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call CAR_SPH(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)), 
     & %VAL(argv(4)))

      car2sph = 9.9

      RETURN
      END
c	
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION rll2gdz(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine geo2gsm1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine geo2gsm: 6 arguments
      call RLL_GDZ(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)), 
     & %VAL(argv(4)))

      rll2gdz = 9.9

      RETURN
      END
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION gse2hee(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine gse2hee1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine gse2hee1: 5 arguments
      call gse2hee1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      gse2hee = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION hee2gse(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine hee2gse1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine hee2gse1: 5 arguments
      call hee2gse1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      hee2gse = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION hae2hee(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine hae2hee1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine hae2hee1: 5 arguments
      call hae2hee1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      hae2hee = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION hee2hae(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine hee2hae1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine hee2hae1: 5 arguments
      call hee2hae1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      hee2hae = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION hae2heeq(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine hae2heeq1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine hae2heeq1: 5 arguments
      call hae2heeq1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      hae2heeq = 9.9

      RETURN
      END
c
C----------------------------------------------------------------------------- 
C Wrapper and procedure for ONERA library
C----------------------------------------------------------------------------- 

      REAL*4 FUNCTION heeq2hae(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.

c  Call subroutine heeq2hae1, converting the IDL parameters to standard FORTRAN
c  passed by reference arguments.
c
c  subroutine heeq2hae1: 5 arguments
      call heeq2hae1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)))

      heeq2hae = 9.9

      RETURN
      END
c
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - September 2005
! MODIFICATION: None
!
! DESCRIPTION: Wrapper to call fly_in_nasa_aeap1 (IN AE8_AP8.f) from IDL, converts the IDL parameters to 
!              standard FORTRAN passed by reference arguments. 
!
! INPUT: argc-> number of argument (long integer)
!        argv -> reference argument
!
! CALLING SEQUENCE: result=call_external(lib_name, 'fly_in_nasa_aeap_', ntime,sysaxes,whichm,whatf,energy,xIN1,xIN2,xIN3,flux, /f_value)
!---------------------------------------------------------------------------------------------------
      REAL*4 FUNCTION fly_in_nasa_aeap(argc, argv)   ! Called by IDL
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers

       j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.
!
      call fly_in_nasa_aeap1(%VAL(argv(1)), %VAL(argv(2)), 
     * %VAL(argv(3)),%VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),
     * %VAL(argv(7)),  %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)),
     * %VAL(argv(11)),  %VAL(argv(12)),  %VAL(argv(13)))

      fly_in_nasa_aeap = 9.9

      RETURN
      END
c
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: S. Bourdarie - March 2007 (add multi channel calculcations) - V4.1
!
! DESCRIPTION: Wrapper to call fly_in_afrl_crres1 (IN AFRL_CRRES_models.f) from IDL, converts the IDL parameters to 
!              standard FORTRAN passed by reference arguments. 
!
! INPUT: argc-> number of argument (long integer)
!        argv -> reference argument
!
! CALLING SEQUENCE: result=call_external(lib_name, 'fly_in_afrl_crres_', ntime,sysaxes,whichm,whatf,energy,xIN1,xIN2,xIN3,flux, /f_value)
!---------------------------------------------------------------------------------------------------
      REAL*4 FUNCTION fly_in_afrl_crres(argc, argv)   ! Called by IDL
!
      INTEGER   CHAR_SIZE
      PARAMETER	(CHAR_SIZE=500)

      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers
!
      j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.
c

      call fly_in_afrl_crres1(%VAL(argv(1)), %VAL(argv(2)), 
     * %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)),  %VAL(argv(9)),  %VAL(argv(10)), %VAL(argv(11)),
     * %VAL(argv(12)), %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     * %VAL(argv(16)))

      fly_in_afrl_crres = 9.9

      RETURN
      END
c
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.0
!
! CREATION: S. Bourdarie - January 2007
! MODIFICATION: None
!
! DESCRIPTION: Wrapper to call SPG4_TLE1 from IDL, converts the IDL parameters to 
!              standard FORTRAN passed by reference arguments. 
!
! INPUT: argc-> number of argument (long integer)
!        argv -> reference argument
!
! CALLING SEQUENCE: result=call_external(lib_name, 'SGP4_TLE_', runtype,startsfe,stopsfe,deltasec,InFileByte,strlenIn,OutFileByte,strlenOut, /f_value)
!---------------------------------------------------------------------------------------------------
      REAL*4 FUNCTION SGP4_TLE(argc, argv)   ! Called by IDL
!
      INTEGER   CHAR_SIZE
      PARAMETER	(CHAR_SIZE=500)

      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers
!
      j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.
c

      call SGP4_TLE1(%VAL(argv(1)), %VAL(argv(2)), 
     * %VAL(argv(3)),
     * %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  %VAL(argv(7)),
     * %VAL(argv(8)))

      SGP4_TLE = 9.9

      RETURN
      END
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.0
!
! CREATION: S. Bourdarie - January 2007
! MODIFICATION: None
!
! DESCRIPTION: Wrapper to call SPG4_ORB1 from IDL, converts the IDL parameters to 
!              standard FORTRAN passed by reference arguments. 
!
! INPUT: argc-> number of argument (long integer)
!        argv -> reference argument
!
! CALLING SEQUENCE: result=call_external(lib_name, 'SGP4_ORB_', runtype,startsfe,stopsfe,deltasec,InFileByte,strlenIn,OutFileByte,strlenOut, /f_value)
!---------------------------------------------------------------------------------------------------
      REAL*4 FUNCTION SGP4_ELE(argc, argv)   ! Called by IDL
!
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers
!
      j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.
c

      call SGP4_ELE1(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),%VAL(argv(5)), %VAL(argv(6)), %VAL(argv(7)),
     * %VAL(argv(8)),%VAL(argv(9)), %VAL(argv(10)), %VAL(argv(11)),
     * %VAL(argv(12)),%VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     * %VAL(argv(16)),%VAL(argv(17)),%VAL(argv(18)),%VAL(argv(19)),
     * %VAL(argv(20)),%VAL(argv(21)),%VAL(argv(22)),%VAL(argv(23)))

      SGP4_ELE= 9.9

      RETURN
      END
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.0
!
! CREATION: S. Bourdarie - January 2007
! MODIFICATION: None
!
! DESCRIPTION: Wrapper to call rv2coe (in sgp4ext.f) from IDL, converts the IDL parameters to 
!              standard FORTRAN passed by reference arguments. 
!
! INPUT: argc-> number of argument (long integer)
!        argv -> reference argument
!
! CALLING SEQUENCE: result=call_external(lib_name, 'RV2COE_IDL_', R, V, P, A, Ecc, Incl, Omega, Argp, Nu, M, ArgLat, TrueLon, LonPer, /f_value)
!---------------------------------------------------------------------------------------------------
      REAL*4 FUNCTION RV2COE_IDL(argc, argv)   ! Called by IDL
!
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers
!
      j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.
c

      call rv2coe(%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)),
     * %VAL(argv(4)),%VAL(argv(5)), %VAL(argv(6)), %VAL(argv(7)),
     * %VAL(argv(8)),%VAL(argv(9)), %VAL(argv(10)), %VAL(argv(11)),
     * %VAL(argv(12)),%VAL(argv(13)))

      RV2COE_IDL= 9.9

      RETURN
      END
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.1
!
! CREATION: S. Bourdarie - March 2007
! MODIFICATION: None
!
! DESCRIPTION: Wrapper to call DATE_AND_TIME2DECY from IDL, converts the IDL parameters to 
!              standard FORTRAN passed by reference arguments. 
!
! INPUT: argc-> number of argument (long integer)
!        argv -> reference argument
!
! CALLING SEQUENCE: result=call_external(lib_name, 'DATE_AND_TIME2DECY_IDL_', Year,Month,Day,hour,minute,second,decy, /f_value)
!---------------------------------------------------------------------------------------------------
      REAL*4 FUNCTION DATE_AND_TIME2DECY_IDL(argc, argv)   ! Called by IDL
!
      INCLUDE 'wrappers.inc'
c      INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers
!
      j = loc(argc)                    ! Obtains the number of arguments (argc)
                                       ! Because argc is passed by VALUE.
c

      call DATE_AND_TIME2DECY(%VAL(argv(1)), %VAL(argv(2)), 
     * %VAL(argv(3)),
     * %VAL(argv(4)),%VAL(argv(5)), %VAL(argv(6)), %VAL(argv(7)))

      DATE_AND_TIME2DECY_IDL= 9.9

      RETURN
      END
c

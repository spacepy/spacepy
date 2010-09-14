!***************************************************************************************************
! Copyright 2010 T.P. O'Brien
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
      
C-----------------------------------------------------------------------------
C     Wrappers and procedures for ONERA_DESP_LIB
C-----------------------------------------------------------------------------
      REAL*4 FUNCTION drift_bounce_orbit(argc, argv) ! Called by IDL
      INCLUDE 'wrappers.inc'
c     INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers
      
c     Call subroutine drift_bounce_orbit1, 
c     converting the IDL parameters to standard FORTRAN
c     passed by reference arguments.
c     
      call drift_bounce_orbit1(%VAL(argv(1)), %VAL(argv(2)),
     +     %VAL(argv(3)),
     *     %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  
     +     %VAL(argv(7)),  %VAL(argv(8)),  %VAL(argv(9)),  
     +     %VAL(argv(10)), %VAL(argv(11)), %VAL(argv(12)), 
     +     %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     +     %VAL(argv(16)), %VAL(argv(17)), %VAL(argv(18)),
     +     %VAL(argv(19)))
      
      drift_bounce_orbit = 9.9  ! return value (copied from make_lstar_splitting
      
      RETURN
      END

      REAL*4 FUNCTION drift_bounce_orbit2(argc, argv) ! Called by IDL
      INCLUDE 'wrappers.inc'
c     INTEGER*4 argc, argv(*)                      ! Argc and Argv are integers
      
c     Call subroutine drift_bounce_orbit1, 
c     converting the IDL parameters to standard FORTRAN
c     passed by reference arguments.
c     
      call drift_bounce_orbit2_1(%VAL(argv(1)), %VAL(argv(2)),
     +     %VAL(argv(3)),
     *     %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)),  
     +     %VAL(argv(7)),  %VAL(argv(8)),  %VAL(argv(9)),  
     +     %VAL(argv(10)), %VAL(argv(11)), %VAL(argv(12)), 
     +     %VAL(argv(13)), %VAL(argv(14)), %VAL(argv(15)),
     +     %VAL(argv(16)), %VAL(argv(17)), %VAL(argv(18)),
     +     %VAL(argv(19)), %VAL(argv(20)), %VAL(argv(21)),
     +     %VAL(argv(22)))
      
      drift_bounce_orbit2 = 9.9  ! return value (copied from make_lstar_splitting
      
      RETURN
      END
c     
c     --------------------------------------------------------------------
c     

      SUBROUTINE drift_bounce_orbit1(kext,options,sysaxes,
     &     iyearsat,idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,
     &     Lm,Lstar,BLOCAL,BMIN,BMIR,XJ,posit,ind)
      ! inputs
      INTEGER*4 kext, options(5),sysaxes,iyearsat,idoy
      REAL*8 UT,xIN1,xIN2,xIN3,alpha,maginput(25)
      ! outputs
      REAL*8 Lm,Lstar
      REAL*8 BLOCAL(1000,25),BMIN,BMIR,XJ,posit(3,1000,25)
      INTEGER*4 ind(25)
      ! locals
      REAL*8 R0,hmin,hmin_lon
      R0 = 1.D0 ! stop search at surface of Earth
      call drift_bounce_orbit2_1(kext,options,sysaxes,
     &     iyearsat,idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,
     &     R0,
     &     Lm,Lstar,BLOCAL,BMIN,BMIR,XJ,posit,ind,
     &     hmin,hmin_lon)
      end

      
      SUBROUTINE drift_bounce_orbit2_1(kext,options,sysaxes,
     &     iyearsat,idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,
     &     R0,
     &     Lm,Lstar,BLOCAL,BMIN,BMIR,XJ,posit,ind,
     &     hmin,hmin_lon)
      
c     computes posit(3,1000,25), BLOCAL(1000,25) and ind(25) in same
c     format as drift_shell1 (note 25 not 48 azimuths)
c     also provides Bmin, Bmirror and usual stuff
c     allows user to specify R0 = surface below which a particle is considered
c     lost (R0 is in RE and the usual value is 1 but for the drift loss cone, one
c     may wish to allow R0<1)
c     also provides hmin and hmin_lon, the altitude and longitude (GDZ) of the
c     minimum altitude along the drift orbit
      
      
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C     
c     declare inputs
      INTEGER*4    kext,k_ext,k_l,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     alpha
      real*8     maginput(25),R0
      REAL*8 hmin,hmin_lon
c     1: Kp
c     2: Dst
c     3: dens
c     4: velo
c     5: Pdyn
c     6: ByIMF
c     7: BzIMF
c     8: G1
c     9: G2
c     10: G3
c     
c     
c     Declare internal variables
      INTEGER*4    iyear,kint,i,j,k
      INTEGER*4    Ndays,activ,Ilflag,t_resol,r_resol
      INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
      INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
      REAL*8     yearsat,dec_year,a2000_ut
      REAL*8     psi,tilt,BL,BMIN_tmp
      REAL*8     xGEO(3),xMAG(3)
      REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
      real*8     alti,lati,longi
      REAL*8     ERA,AQUAD,BQUAD
      real*8     density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
      real*8     G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
      real*8     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
c     
c     Declare output variables
      INTEGER*4  ind(25)
      REAL*8     posit(3,1000,25)
      REAL*8     BLOCAL(1000,25)
      REAL*8     BMIN,BMIR
      REAL*8     XJ
      REAL*8     Lm,Lstar
C     
      COMMON/GENER/ERA,AQUAD,BQUAD
      COMMON /dip_ang/tilt
      COMMON /magmod/k_ext,k_l,kint
      COMMON /drivers/density,speed,dst_nt,Pdyn_nPa,ByIMF_nt,BzIMF_nt
     &     ,G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04,
     &     W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
      COMMON /index/activ
      COMMON /flag_L/Ilflag
      COMMON /a2000_time/a2000_ut,a2000_iyear,a2000_imonth,a2000_iday
C     

c     initialize outputs
      Lm=baddata
      Lstar=baddata
      XJ=baddata
      BMIN=baddata
      BMIR=baddata
      hmin = baddata
      hmin_lon = baddata
      do i=1,25
         ind(i)=0
         do j=1,1000
            BLOCAL(j,i)=baddata
            do k=1,3
               posit(k,j,i)=baddata
            enddo
         enddo
      enddo

      Ilflag=0
      iyear=1800
      k_ext=kext
      if (options(1).eq.0) options(1)=1
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
      if (kint .le. 1) then
         if (options(2) .eq. 0) then
            if (iyearsat .ne. iyear) then
               iyear=iyearsat
               dec_year=iyear+0.5d0
               CALL INIT_DTD(dec_year)
            endif
         else
            if (iyearsat .ne. iyear .or.
     &           MOD(idoy*1.d0,options(2)*1.d0) .eq. 0) THEN
               iyear=iyearsat
               firstJanuary=JULDAY(iyear,01,01)
               lastDecember=JULDAY(iyear,12,31)
               currentdoy=(idoy/options(2))*options(2)
               if (currentdoy .eq. 0) currentdoy=1
               dec_year=iyear+currentdoy*1.d0/
     &              ((lastDecember-firstJanuary+1)*1.d0)
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
c     make inputs according to magn. field model chosen
c     
      if (kext .eq. 1) then
c     Input for MEAD
         if (maginput(1).le.3.d0) Activ=1
         if (maginput(1).gt.3.d0 .and.
     &        maginput(1).lt.20.d0) Activ=2
         if (maginput(1).ge.20.d0 .and.
     &        maginput(1).lt.30.d0) Activ=3
         if (maginput(1).ge.30.d0) Activ=4
c     
         if (maginput(1).lt.0.d0 .or.
     &        maginput(1).gt.90.d0) then
            GOTO 99
         endif
      endif
      if (kext .eq. 2) then
c     Input for TSYG87s
         if (maginput(1).lt.7.d0) Activ=1
         if (maginput(1).ge.7.d0 .and.
     &        maginput(1).lt.17.d0) Activ=2
         if (maginput(1).ge.17.d0 .and.
     &        maginput(1).lt.20.d0) Activ=3
         if (maginput(1).ge.20.d0 .and.
     &        maginput(1).lt.27.d0) Activ=4
         if (maginput(1).ge.27.d0 .and.
     &        maginput(1).lt.37.d0) Activ=5
         if (maginput(1).ge.37.d0 .and.
     &        maginput(1).lt.47.d0) Activ=6
         if (maginput(1).ge.47.d0) Activ=7
         if (maginput(1).ge.53.d0) Activ=8
c     
         if (maginput(1).lt.0.d0 .or.
     &        maginput(1).gt.90.d0) then
            GOTO 99
         endif
      endif
      if (kext .eq. 3) then
c     Input for TSYG87l
         if (maginput(1).lt.7.d0) Activ=1
         if (maginput(1).ge.7.d0 .and.
     &        maginput(1).lt.17.d0) Activ=2
         if (maginput(1).ge.17.d0 .and.
     &        maginput(1).lt.27.d0) Activ=3
         if (maginput(1).ge.27.d0 .and.
     &        maginput(1).lt.37.d0) Activ=4
         if (maginput(1).ge.37.d0 .and.
     &        maginput(1).lt.47.d0) Activ=5
         if (maginput(1).ge.47.d0) Activ=6
c     
         if (maginput(1).lt.0.d0 .or.
     &        maginput(1).gt.90.d0) then
            GOTO 99
         endif
      endif
      if (kext .eq. 4) then
c     Input for Tsy89
         if (maginput(1).lt.7.d0) Activ=1
         if (maginput(1).ge.7.d0 .and.
     &        maginput(1).lt.17.d0) Activ=2
         if (maginput(1).ge.17.d0 .and.
     &        maginput(1).lt.27.d0) Activ=3
         if (maginput(1).ge.27.d0 .and.
     &        maginput(1).lt.37.d0) Activ=4
         if (maginput(1).ge.37.d0 .and.
     &        maginput(1).lt.47.d0) Activ=5
         if (maginput(1).ge.47.d0 .and.
     &        maginput(1).lt.57.d0) Activ=6
         if (maginput(1).ge.57.d0) Activ=7
c     
         if (maginput(1).lt.0.d0 .or.
     &        maginput(1).gt.90.d0) then
            GOTO 99
         endif
      endif
      if (kext .eq. 6) then
c     Input for OP dyn
         density=maginput(3)
         speed=maginput(4)
         dst_nt=maginput(2)
c     
         if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
            GOTO 99
         endif
         if (density.lt.5.d0 .or. density.gt.50.d0) then
            GOTO 99
         endif
         if (speed.lt.300.d0 .or. speed.gt.500.d0) then
            GOTO 99
         endif
      endif
      if (kext .eq. 7) then
c     Input for Tsy96
         dst_nt=maginput(2)
         Pdyn_nPa=maginput(5)
         ByIMF_nt=maginput(6)
         BzIMF_nt=maginput(7)
c     
         if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) then
            GOTO 99
         endif
         if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) then
            GOTO 99
         endif
         if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) then
            GOTO 99
         endif
         if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) then
            GOTO 99
         endif
      endif
      if (kext .eq. 8) then
c     Input for Ostapenko97
         dst_nt=maginput(2)
         Pdyn_nPa=maginput(5)
         BzIMF_nt=maginput(7)
         fkp=maginput(1)*1.d0/10.d0
      endif
      if (kext .eq. 9) then
c     Input for Tsy01
         dst_nt=maginput(2)
         Pdyn_nPa=maginput(5)
         ByIMF_nt=maginput(6)
         BzIMF_nt=maginput(7)
         G1_tsy01=maginput(8)
         G2_tsy01=maginput(9)
c     
         if (dst_nt.lt.-50.d0 .or. dst_nt.gt.20.d0) then
            GOTO 99
         endif
         if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) then
            GOTO 99
         endif
         if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) then
            GOTO 99
         endif
         if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) then
            GOTO 99
         endif
         if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) then
            GOTO 99
         endif
         if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) then
            GOTO 99
         endif
      endif
      if (kext .eq. 10) then
c     Input for Tsy01 storm
         dst_nt=maginput(2)
         Pdyn_nPa=maginput(5)
         ByIMF_nt=maginput(6)
         BzIMF_nt=maginput(7)
         G2_tsy01=maginput(9)
         G3_tsy01=maginput(10)
      endif
c     
      if (kext .eq. 11) then
c     Input for Tsy04 storm
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
c     Input for Alexeev 2000
         a2000_iyear=iyearsat
         firstJanuary=JULDAY(a2000_iyear,01,01)
         currentdoy=firstJanuary+idoy-1
         CALL CALDAT(currentdoy,a2000_iyear,
     &        a2000_imonth,a2000_iday)
         a2000_ut=UT
         density=maginput(3)
         speed=maginput(4)
         dst_nt=maginput(2)
         BzIMF_nt=maginput(7)
         Al=maginput(17)
      endif
c     
      if (alpha.ne.90.0d0) then
         CALL find_bm(
     &        lati,longi,alti,alpha,BL,BMIR,xGEO)
         CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,
     &        longi,alti)
      ELSE
         Bmir=0.0D0
      ENDIF
      IF (Bmir.NE.baddata) THEN
         Ilflag=0
         call trace_drift_bounce_orbit(t_resol,r_resol,
     &        lati,longi,alti,R0,Lm,Lstar,XJ,
     &        BLOCAL,Bmin,Bmir,posit,ind,hmin,hmin_lon)
         if (Lstar.eq.baddata) then
            hmin = baddata
            hmin_lon = baddata
         endif
         
      ENDIF
 99   continue
      
      end                       ! end subroutine drift_bounce_orbit1


c     --------------------------------------------------------------------
c     
      
      SUBROUTINE trace_drift_bounce_orbit(t_resol,r_resol,
     &     lati,longi,alti,R0,Lm,Lstar,leI0,Bposit,
     &     Bmin,Bmir,posit,Nposit,hmin,hmin_lon)
C     
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      
C     Parameters
      INTEGER*4  Nreb_def,Nder_def,Ntet_def
      PARAMETER (Nreb_def = 50, Nder_def = 25, Ntet_def = 720)
      
C     Input Variables
      INTEGER*4 t_resol,r_resol
      REAL*8     lati,longi,alti,R0
      
C     Internal Variables
      INTEGER*4  Nder,Nreb,Ntet
      INTEGER*4  k_ext,k_l,kint,n_resol
      INTEGER*4  Nrebmax
      REAL*8     rr,rr2
      REAL*8     xx0(3),xx(3),x1(3),x2(3)
      REAL*8     xmin(3)
      REAL*8     B(3),Bl,B0,B1,B3
      REAL*8     dsreb0,dsreb,smin
      
      INTEGER*4  I,J,K,Iflag,Iflag_I,Ilflag,Ifail
      INTEGER*4  Ibounce_flag
      INTEGER*4  istore         ! azimuth cursor
      REAL*8     Lb
      REAL*8     leI,leI1
      REAL*8     XY,YY
      REAL*8     aa,bb
      REAL*8     tmp_lati,tmp_longi,tmp_alti
      
C     
      REAL*8     pi,rad,tt
      REAL*8     tet(10*Nder_def),phi(10*Nder_def)
      REAL*8     tetl,tet1,dtet,lasttet
      REAL*8     somme,R02

c variables to deal with leI~0 case
      REAL*8     x1old(3)
      INTEGER*4  I0flag
C     
      REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
      
C     Output Variables       
      REAL*8     Lm,Lstar,leI0,Bmin,Bmir
      REAL*8     posit(3,20*Nreb_def,Nder_def)
      REAL*8     Bposit(20*Nreb_def,Nder_def)
      REAL*8     hmin,hmin_lon
      INTEGER*4  Nposit(Nder_def)
C     
      COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
      COMMON /calotte/tet
      COMMON /flag_L/Ilflag
      COMMON /magmod/k_ext,k_l,kint
C     
C     


      R02 = R0*R0
      Nder=Nder_def*r_resol     ! longitude steps
      Nreb=Nreb_def             ! steps along field line
      Ntet=Ntet_def*t_resol     ! latitude steps
      pi = 4.D0*ATAN(1.D0)
      rad = pi/180.D0
      dtet = pi/Ntet            ! theta (latitude) step
C     
      Nrebmax = 20*Nreb         ! maximum steps along field line
C     
C
c     initialize hmin,hmin_lon to starting point
      hmin = alti
      hmin_lon = longi

      CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
C     
      CALL GEO_SM(xx0,xx)
      rr = SQRT(xx(1)*xx(1)+xx(2)*xx(2)+xx(3)*xx(3))
      tt = ACOS(xx(3)/rr)
      Lb  = rr/SIN(tt)/SIN(tt)  ! dipole L
C     
      CALL CHAMP(xx0,B,B0,Ifail)
      Bmir = B0 ! local field at starting point
      IF (Ifail.LT.0) THEN
         Ilflag = 0
         RETURN
      ENDIF
      Bmin = B0
C     
      dsreb0 = Lb/Nreb           ! step size dipole L / Nsteps
      dsreb = dsreb0
C     
C     calcul du sens du depart
C     (compute hemisphere)
 10   
      CALL sksyst(-dsreb,xx0,x1,Bl,Ifail)
      IF (Ifail.LT.0) THEN
         Ilflag = 0
         RETURN
      ENDIF
      B1 = Bl
      
      CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
      IF (Ifail.LT.0) THEN
         Ilflag = 0
         RETURN
      ENDIF
      B3 = Bl

C     
C     attention cas equatorial
C     (equatorial special case)
      IF(B1.GT.B0 .AND. B3.GT.B0)THEN
         aa = 0.5D0*(B3+B1-2.D0*B0)
         bb = 0.5D0*(B3-B1)
         smin = -0.5D0*bb/aa
         Bmin = B0 - aa*smin*smin
         leI0 = SQRT(1.D0-Bmin/B0)*2.D0*ABS(smin*dsreb)
         Lm = (Bo/Bmin)**(1.D0/3.D0)
c     write(6,*)'L McIlwain eq ',B0,leI0,Lm
         GOTO 100
      ENDIF
      dsreb0 = dsreb ! keep smaller step size for future use
      IF (B3.GT.B1) THEN
         dsreb = -dsreb
      ENDIF
C     
C     calcul de la ligne de champ et de I
C     (compute field line and I)
      Bmin = B0
      B1 = B0
      leI = 0.D0
      DO I = 1,3
         x1(I)  = xx0(I)
      ENDDO
C     
      DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
            RETURN
         ENDIF
         IF (Bl.LT.Bmin) THEN
            xmin(1) = x2(1)
            xmin(2) = x2(2)
            xmin(3) = x2(3)
            Bmin = Bl
         ENDIF
         IF (Bl.GT.B0) GOTO 20  ! traced past southern mirror point
         x1(1) = x2(1)
         x1(2) = x2(2)
         x1(3) = x2(3)
         leI = leI + SQRT(1.D0-Bl/B0)
         
         B1 = Bl
      ENDDO
 20   CONTINUE

C     
      IF (J.GE.Nrebmax) THEN    !open field line
         Ilflag = 0
         RETURN
      ENDIF

C     
      leI = leI+0.5D0*SQRT(1.D0-B1/B0)*(B0-Bl)/(Bl-B1)
      leI = leI*ABS(dsreb)
      leI0 = leI
C     
C     calcul de L Mc Ilwain (Mc Ilwain-Hilton)
C     (compute L McIlwain (McIlwain-Hilton))
C     
      XY = leI*leI*leI*B0/Bo
      YY = 1.D0 + 1.35047D0*XY**(1.D0/3.D0)
     &     + 0.465376D0*XY**(2.D0/3.D0)
     &     + 0.0475455D0*XY
      Lm = (Bo*YY/B0)**(1.D0/3.D0)
C     
C     calcul de Bmin
c     (compute Bmin)
C     
      CALL sksyst(dsreb,xmin,x1,B3,Ifail)
      IF (Ifail.LT.0) THEN
         Ilflag = 0
         RETURN
      ENDIF
      CALL sksyst(-dsreb,xmin,x1,B1,Ifail)
      IF (Ifail.LT.0) THEN
         Ilflag = 0
         RETURN
      ENDIF
      aa = 0.5D0*(B3+B1-2.D0*Bmin)
      bb = 0.5D0*(B3-B1)
      smin = -0.5D0*bb/aa
      Bmin = Bmin - aa*smin*smin
      IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.0) THEN
         Lm = -Lm
      ENDIF
C     
 100  CONTINUE
      if (k_l .eq.0) then
         Ilflag = 0
         RETURN
      endif
      IF (ABS(Lm) .GT. 10.D0) THEN
         Ilflag = 0
         RETURN
      ENDIF
C     
C     calcul du point sur la ligne de champ a la surface de la terre du
C     cote nord
C     (compute the point nothern footpoint at Earth's surface)
C     
      DO I = 1,3
         x1(I)  = xx0(I)
      ENDDO
      dsreb = ABS(dsreb)
      DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
	    RETURN
	 ENDIF
	 rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
	 IF (rr.LT.R0) GOTO 102
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
      ENDDO
 102  CONTINUE
      smin = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
      smin = (R0-smin)/(rr-smin)
      CALL sksyst(smin*dsreb,x1,x2,Bl,Ifail)
      IF (Ifail.LT.0) THEN
         Ilflag = 0
         RETURN
      ENDIF
      rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
      tet(1) = ACOS(x2(3)/rr)
      phi(1) = ATAN2(x2(2),x2(1))
C     
C     et on tourne -> on se decale sur la surface en phi et on cherche teta
C     pour avoir leI0 et B0 constants
C     (find the thetta/phi contour on the surface of the earth that conserves
C     B0= Bmirror and leI0=I)
      call trace_bounce_orbit(x2,B0,1,dsreb0,R0,
     &     Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
      if (Ibounce_flag.ne.1) then ! try again from starting point
         call trace_bounce_orbit(xx0,B0,1,dsreb0,R0,
     &        Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
      endif
      if (Ibounce_flag.ne.1) then
         Ilflag = 0
         RETURN
      endif

      istore = 2 ! istore=1 already stored in first bounce orbit trace
      dsreb = -dsreb
      DO I = 2,Nder
         phi(I) = phi(I-1)+2.D0*pi/Nder
         Iflag_I = 0
         I0flag = 0
         IF (Ilflag.EQ.0) THEN
            tetl = tet(I-1)
            IF (I.GT.2) tetl = 2.D0*tet(I-1)-tet(I-2)
            tet1 = tetl
         ELSE
            tetl = tet(I)
            tet1 = tetl
         ENDIF
         leI1 = baddata
C     
 107     CONTINUE
         x1(1) = R0*SIN(tetl)*COS(phi(I))
         x1(2) = R0*SIN(tetl)*SIN(phi(I))
         x1(3) = R0*COS(tetl)
         Iflag = 0
         leI = baddata
         CALL CHAMP(x1,B,B1,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
            RETURN
         ENDIF
C     
         dsreb = dsreb/abs(dsreb)*dsreb0
         DO J = 1,Nrebmax
 109        CALL sksyst(dsreb,x1,x2,Bl,Ifail)
            IF (Ifail.LT.0) THEN
               Ilflag = 0
               RETURN
            ENDIF
            rr2 = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
            IF (Bl.LT.B0) THEN
               lasttet = tetl
               IF (Iflag .EQ. 0) THEN
                  CALL CHAMP(x1,B,B1,Ifail)
                  IF (Ifail.LT.0) THEN
                     Ilflag = 0
                     RETURN
                  ENDIF
                  leI = 0.5D0*SQRT(1.D0-Bl/B0)*(1.D0+(Bl-B0)/(Bl-B1))
     &                 *abs(dsreb)
                  Iflag = 1
               ELSE
                  leI = leI+SQRT(1.D0-Bl/B0)*abs(dsreb)
               ENDIF
            ENDIF
            IF (Bl.GT.B0 .AND. Iflag.EQ.1) GOTO 103
            IF (rr2.LT.R0) GOTO 103
            x1(1) = x2(1)
            x1(2) = x2(2)
            x1(3) = x2(3)
            B1 = Bl
         ENDDO
 103     CONTINUE
         if ((Iflag.eq.0).and.(J.LT.Nrebmax)) then ! never got to B < Bmirror
            leI=0 ! assume particle is equatorially mirroring
                  ! Not technically true, but effectively gives desired result
                  ! (i.e., leI=0 for field lines way inside trajectory)
            I0flag = 1
         else
            I0flag = 0
c     Pourquoi? (why? I don't know!)
            IF (rr2.LT.R0) THEN
               leI = baddata
            ENDIF
            IF (J.LT.Nrebmax .AND. rr2.GE.R0) THEN
               CALL CHAMP(x1,B,B1,Ifail)
               IF (Ifail.LT.0) THEN
                  Ilflag = 0
                  RETURN
               ENDIF
               lasttet = tetl
               leI = leI+0.5D0*SQRT(1.D0-B1/B0)*(B0-Bl)/(Bl-B1)*
     &              abs(dsreb)
            ENDIF
         endif
C     
         IF (Iflag_I .EQ.0) THEN
            IF (J.GE.Nrebmax) THEN
               tetl = tetl-dtet
            ELSE
               tetl = tetl+dtet
            ENDIF
            leI1 = leI
            tet1 = tetl
            Iflag_I = 1
            GOTO 107
         ENDIF
         IF ((leI-leI0)*(leI1-leI0) .LT. 0.D0) GOTO 108
         leI1 = leI
         tet1 = tetl
         IF (leI.LT.leI0) THEN
            tetl = tetl-dtet
         ElSE
            tetl = tetl+dtet
         ENDIF
         IF (tetl.GT.pi .OR. tetl.LT.0.D0) GOTO 108
         GOTO 107
 108     CONTINUE
         IF ((J.GE.Nrebmax .AND. leI.GT.0.D0).or.(leI.LT.0.D0)) THEN
            Ilflag = 0
            RETURN
         ENDIF

         tetl = 0.5D0*(tetl+tet1) ! both tetl and tet1 work, so average
         tet(I) = tetl
         x1(1) = R0*SIN(tet(I))*COS(phi(I))
         x1(2) = R0*SIN(tet(I))*SIN(phi(I))
         x1(3) = R0*COS(tet(I))
         CALL CHAMP(x1,B,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
            RETURN
         ENDIF
         IF (Bl.LT.B0) THEN
            Ilflag = 0
            RETURN
         ENDIF
         
! Trace bounce orbit on field line from x1
         call trace_bounce_orbit(x1,Bmir,istore,dsreb0,R0,
     &        Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
         if (Ibounce_flag.ne.1) then
            ! try again from lasttet, which worked
            tet(I) = lasttet
            x1(1) = R0*SIN(tet(I))*COS(phi(I))
            x1(2) = R0*SIN(tet(I))*SIN(phi(I))
            x1(3) = R0*COS(tet(I))
            call trace_bounce_orbit(x1,Bmir,istore,dsreb0,R0,
     &           Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
         endif
         if (Ibounce_flag.ne.1) then
            Ilflag = 0
            RETURN
         endif
         if (MOD(I-1,r_resol).eq.0) then
            istore = istore+1 ! next time, store in next array step
         endif                  ! mode(I,r_resol)
      ENDDO                     ! end of do I = 2,Nder
C     
C     calcul de somme de BdS sur la calotte nord
C     (compute the integral of BdS on the norther polar cap)
      x1(1) = 0.D0
      x1(2) = 0.D0
      x1(3) = 1.D0
      CALL CHAMP(x1,B,Bl,Ifail)
      IF (Ifail.LT.0)THEN
         Ilflag = 0
         RETURN
      ENDIF
      somme = Bl*pi*dtet*dtet/4.D0*R02
      DO I = 1,Nder
         tetl = 0.D0
         DO J = 1,Ntet
            tetl = tetl+dtet
            IF (tetl .GT. tet(I)) GOTO 111
            x1(1) = R0*SIN(tetl)*COS(phi(I))
            x1(2) = R0*SIN(tetl)*SIN(phi(I))
            x1(3) = R0*COS(tetl)
            CALL CHAMP(x1,B,Bl,Ifail)
            IF (Ifail.LT.0)THEN
               Ilflag = 0
               RETURN
            ENDIF
            somme = somme+Bl*SIN(tetl)*dtet*2.D0*pi/Nder*R02
	 ENDDO
 111     CONTINUE
      ENDDO
      if (k_l .eq.1) Lstar = 2.D0*pi*Bo/somme
      if (k_l .eq.2) Lstar = somme ! Phi and not Lstar
      IF (Lm.LT.0.D0) Lstar = -Lstar
      Ilflag = 1
C     
      END
      
      SUBROUTINE trace_bounce_orbit(xstart,Bmirror,istore,dsreb0,
     &     R0,Bposit,posit,Nposit,Iflag,hmin,hmin_lon)

      IMPLICIT NONE

c     Declare input variables
      REAL*8 xstart(3) ! GEO
      REAL*8 Bmirror  ! particle's mirror field strength
      INTEGER*4 istore ! where to store bounce orbit
      ! hmin,hmin_lon are also inputs

c     Declare output variables
      INTEGER*4  Nposit(25),Iflag ! Iflag=1 means success
      REAL*8     posit(3,1000,25)
      REAL*8     Bposit(1000,25),R0,hmin,hmin_lon

c     Declare internal variables
      INTEGER*4 i,j,k,Ifail
      REAL*8 lati,longi,alti,Bl,Bmir
      REAL*8 xmir(3)
      REAL*8 pi,rad,alpha
      REAL*8 dsreb0,dsreb,x1(3),x2(3),B1,B2,Bvec(3)
      INTEGER*4 store ! store or not?

      Iflag = 0

      if ((istore.gt.25).or.(istore.lt.1)) then
         store = 0
      else
         store = 1
      endif

      pi = 4.D0*ATAN(1.D0)
      rad = pi/180.D0

c     trace up to mirror point, leave starting pointin Bmir, xmir

      call champ(xstart,Bvec,B1,Ifail)

      dsreb = -dsreb0 ! assume starting at northern R=1 foot point

      do i = 1,3
         x1(i) = xstart(i)
      enddo

      call sksyst(dsreb,x1,x2,B2,Ifail)
      if (B2.GT.B1) then
         dsreb = -dsreb ! going wrong way
      endif

      call sksyst(dsreb,x1,x2,B2,Ifail)
      if (B2.GT.B1) then ! already at local min w/in resolution of dsreb
         ! store & return
         if (store.eq.1) then
            Bposit(1,istore) = B2
            Nposit(istore) = 1
            do k = 1,3
               posit(k,j,istore) = x1(k)
            enddo
         endif
         call check_hmin(x1(1),x1(2),x1(3),hmin,hmin_lon)
         Iflag=1
         return
      endif

      do j = 1,999
         call sksyst(dsreb,x1,x2,B2,Ifail)
         if (B2.GT.Bmirror) then ! still haven't crossed Bmirror
            do i = 1,3
               x1(i) = x2(i)
            enddo
            B1 = B2
         else
            do i = 1,3
               xmir(i) = x2(i)
            enddo
            Bmir = B2
            if (abs((Bmir/Bmirror-1.0D0)).lt.1.0D-6) then
               goto 101         ! good convergence
            endif
            dsreb = dsreb/2.D0 ! try again with smaller step (converges logarithmically)
         endif
      enddo
      return ! did not converge (Iflag=0)
 101  continue

      call check_hmin(xmir(1),xmir(2),xmir(3),hmin,hmin_lon)
      dsreb = -dsreb0

c     trace bounce trajectory
      call sksyst(dsreb,xmir,x2,B2,Ifail)
      if (B2.GT.Bmir) then
         dsreb = -dsreb ! going wrong way
      endif

c     trace to opposite mirror point
c     set up trace
      B2 = BMIR
      do k=1,3
         x1(k)=xmir(k)
      enddo
c     do trace
      do j=1,999
         if (store.eq.1) then
            Bposit(j,istore) = B2
            Nposit(istore) = j
            do k = 1,3
               posit(k,j,istore) = x1(k)
            enddo
         endif
         call sksyst(dsreb,x1,x2,B2,Ifail)
         if (Ifail.LT.0) then
            Iflag = 0
            return
         endif
         if (B2.ge.BMIR) goto 201
         do k=1,3
            x1(k)=x2(k)
         enddo
      enddo
 201  continue
      if (B2.lt.BMIR) then
         Iflag = 0 ! failed to get past BMIR
         return
      elseif (B2.gt.BMIR) then
c     finish trace to Bm: Bm between x1 and x2
         BL = B2 ! closest stored point so far (x1)
         j = j+1                ! prepare to store in next element
         do i = 1,10            ! this converges logarithmically, so reduces step size by up to 2^10
            dsreb = dsreb/2.0D0 ! restart with half step size
            call sksyst(dsreb,x1,x2,B2,Ifail)
            if (B2.LT.Bmirror) then ! still inside bounce orbit
               do k=1,3
                  x1(k) = x2(k)
               enddo
               if (store.eq.1) then ! got closer, so store x1
                  Bposit(j,istore) = B2
                  Nposit(istore) = j
                  do k=1,3
                     posit(k,j,istore) = x1(k)
                  enddo
               endif
            endif
         enddo
         call check_hmin(x1(1),x1(2),x1(3),hmin,hmin_lon)
      endif

      Iflag = 1
      END ! end sub trace_bounce_orbit

      subroutine check_hmin(x,y,z,hmin,hmin_lon)
      ! check alt/lat/lon at x,y,z and replace hmin,hmin_lon if needed
      IMPLICIT NONE
      ! inputs - x,y,z GEO
      real*8 x,y,z
      
      ! outputs
      real*8 hmin,hmin_lon

      ! locals
      real*8 lati,longi,alti

      call geo_gdz(x,y,z,lati,longi,alti)
      if (alti .lt. hmin) then
         hmin = alti
         hmin_lon = longi
      endif
      end



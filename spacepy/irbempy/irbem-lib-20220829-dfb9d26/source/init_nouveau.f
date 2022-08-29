!***************************************************************************************************
! Copyright 2000 D. Boscher, 2004 S. Bourdarie
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
       SUBROUTINE CHAMP(xGEO,BxGEO,Bl,Ifail)
C
       IMPLICIT NONE
C
       INTEGER*4   I,IOPT,Ifail,kint
       INTEGER*4   Activ,k_ext,k_l
       REAL*8      xGEO(3)
       REAL*8      xSM(3)
       REAL*8      BxSM(3),BIMF(3)
       REAL*8      Bxext(3)
       REAL*8      Bxint(3),den,vel,dst
       REAL*8      BxGEO(3),Bl,density,speed,dst_nt,fkp,fkp_osta
       REAL*8      Pdyn_nPa,ByIMF_nt,BzIMF_nt,G1_tsy01,G2_tsy01,G3_tsy01
       REAL*8      W1_tsy04,W2_tsy04,W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04
       REAL*8      Al,Al_ind,BxIMF_nt
       REAL*8      PARMOD(10),tilt,sn,Pdyn,BzIMF
C
       REAL*8      Bo,xc,yc,zc,ct,st,cp,sp
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
       COMMON /magmod/k_ext,k_l,kint
       COMMON /drivers/density,speed,dst_nt,Pdyn_nPa,BxIMF_nt,ByIMF_nt,
     &        BzIMF_nt,G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04,
     &         W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
       COMMON /index/activ
       COMMON /dip_ang/tilt
C
       Ifail=0
       IF(kint.EQ.0)THEN
         CALL IGRF(xGEO(1),xGEO(2),xGEO(3),Bxint(1),Bxint(2),Bxint(3))
       ENDIF
       IF(kint.EQ.1) THEN
         CALL DTD(xGEO(1),xGEO(2),xGEO(3),Bxint(1),Bxint(2),Bxint(3))
       ENDIF
       IF(kint.GE.2 .and. kint .le. 3) THEN
         CALL get_intfield(xGEO(1),xGEO(2),xGEO(3),
     &                    Bxint(1),Bxint(2),Bxint(3))
       ENDIF
       IF(kint .eq. 4) THEN
         CALL myOwnMagField(xGEO,Bxint)
       ENDIF
       IF(kint .eq. 5) THEN
         CALL centered_dipole(xGEO(1),xGEO(2),xGEO(3),Bxint(1),Bxint(2)
     &    ,Bxint(3))
       ENDIF
C
       Bxext(1) = 0.D0
       Bxext(2) = 0.D0
       Bxext(3) = 0.D0
c
c Case for Mead dynamic mag model (use of SM coordinates)
c inputs is magnetic activity index deduced from Kp
c This field model is only valid for rSM<15 Re
       if (k_ext .eq. 1) then
          IOPT=Activ
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.289.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL GEO_SM(xGEO,xSM)
          CALL MEAD(xSM(1),xSM(2),xSM(3),IOPT,BxSM(1),BxSM(2),BxSM(3))
          CALL SM_GEO(BxSM,Bxext)
       endif
c
c Case for Tsyganenko 87 short dynamic mag model (use of GSM coordinates)
c inputs is magnetic activity index deduced from Kp
c This field model is only valid for rGEO<30 Re
       if (k_ext .eq. 2) then
          IOPT=Activ
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.900.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL GEO_GSM(xGEO,xSM)
          CALL TSY87S(IOPT,xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),BxSM(3))
          CALL GSM_GEO(BxSM,Bxext)
       endif
c
c Case for Tsyganenko 87 long dynamic mag model (use of GSM coordinates)
c inputs is magnetic activity index deduced from Kp
c This field model is only valid for rGEO<70 Re
       if (k_ext .eq. 3) then
          IOPT=Activ
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.4900.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL GEO_GSM(xGEO,xSM)
          CALL TSY87l(IOPT,xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),BxSM(3))
          CALL GSM_GEO(BxSM,Bxext)
       endif
c
c Case for Tsyganenko 89 dynamic mag model (use of GSM coordinates)
c inputs is magnetic activity index deduced from Kp
c This field model is only valid for rGEO<70 Re
       if (k_ext .eq. 4) then
          IOPT=Activ
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.4900.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL GEO_GSM(xGEO,xSM)
          CALL T89C(IOPT,xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),BxSM(3))
          CALL GSM_GEO(BxSM,Bxext)
       endif
c Case for Olson-Pfitzer quiet mag model (use of SM coordinates)
c no parameters as input
c This field model is only valid for rGEO<15 Re
       if (k_ext .eq. 5) then
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.225.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL GEO_SM(xGEO,xSM)
          CALL BXYZMU(xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),BxSM(3))
          CALL SM_GEO(BxSM,Bxext)
       endif
c Case for Olson-Pfitzer dynamic mag model (use of SM coordinates)
c inputs are SW density, speed, and Dst index
c This field model is only valid for rGEO<60 Re
       if (k_ext .eq. 6) then
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.3600.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL GEO_SM(xGEO,xSM)
          den=density
           vel=speed
        dst=dst_nt
          CALL BDYN(den,vel,dst,xSM(1),xSM(2),xSM(3),BxSM(1),
     &    BxSM(2),BxSM(3))
          CALL SM_GEO(BxSM,Bxext)
       endif
c
c Case for rescent Tsyganenko dynamic mag model 96 (use of GSM coordinates)
c inputs is SW pressure, (nPa), DST (nT), ByIMF and BzIMF (nT)
c
       if (k_ext .eq. 7) then
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.1600.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          PARMOD(1)=Pdyn_nPa
        PARMOD(2)=dst_nt
          PARMOD(3)=ByIMF_nt
        PARMOD(4)=BzIMF_nt
          CALL GEO_GSM(xGEO,xSM)
          CALL T96_01(PARMOD,xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),
     &      BxSM(3))
          CALL GSM_GEO(BxSM,Bxext)
       endif
c
c Case for Ostapenko dynamic mag model 97 (use of SM coordinates)
c inputs is SW pressure, (nPa), DST (nT), f(Kp) and BzIMF (nT)
       if (k_ext .eq. 8) then
        dst=dst_nt
          Pdyn=Pdyn_nPa
        BzIMF=BzIMF_nt
        fkp_osta=fkp
        sn =SIN(tilt*4.D0*ATAN(1.D0)/180.d0)
          CALL GEO_SM(xGEO,xSM)
        CALL set_a(dst,Pdyn,fkp_osta,BzIMF,sn)
          CALL BOM97(xSM,BxSM)
          CALL SM_GEO(BxSM,Bxext)
       endif
c
c Case for Tsyganenko dynamic mag model 01 (use of GSM coordinates)
c inputs is SW pressure, (nPa), DST (nT), ByIMF and BzIMF (nT), G1,G2
c
       if (k_ext .eq. 9) then
          PARMOD(1)=Pdyn_nPa
        PARMOD(2)=dst_nt
          PARMOD(3)=ByIMF_nt
        PARMOD(4)=BzIMF_nt
          PARMOD(5)=G1_tsy01
        PARMOD(6)=G2_tsy01
          CALL GEO_GSM(xGEO,xSM)
        IF (xSM(1).LT.-15.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL T01_01(PARMOD,xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),
     &      BxSM(3))
          CALL GSM_GEO(BxSM,Bxext)
       endif
c
c Case for Tsyganenko dynamic mag model 01 storm (use of GSM coordinates)
c inputs is SW pressure, (nPa), DST (nT), ByIMF and BzIMF (nT), G1,G2
c
       if (k_ext .eq. 10) then
          PARMOD(1)=Pdyn_nPa
        PARMOD(2)=dst_nt
          PARMOD(3)=ByIMF_nt
        PARMOD(4)=BzIMF_nt
          PARMOD(5)=G2_tsy01
        PARMOD(6)=G3_tsy01
          CALL GEO_GSM(xGEO,xSM)
        IF (xSM(1).LT.-15.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL T01_01(PARMOD,xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),
     &      BxSM(3))
          CALL GSM_GEO(BxSM,Bxext)
       endif
c
c Case for Tsyganenko dynamic mag model 04 storm (use of GSM coordinates)
c inputs is SW pressure, (nPa), DST (nT), ByIMF and BzIMF (nT), W1-6
c
       if (k_ext .eq. 11) then
          PARMOD(1)=Pdyn_nPa
        PARMOD(2)=dst_nt
          PARMOD(3)=ByIMF_nt
        PARMOD(4)=BzIMF_nt
          PARMOD(5)=W1_tsy04
        PARMOD(6)=W2_tsy04
          PARMOD(7)=W3_tsy04
        PARMOD(8)=W4_tsy04
          PARMOD(9)=W5_tsy04
        PARMOD(10)=W6_tsy04
          CALL GEO_GSM(xGEO,xSM)
        IF (xSM(1).LT.-15.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          CALL T04_s(PARMOD,xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),
     &      BxSM(3))
          CALL GSM_GEO(BxSM,Bxext)
       endif
c
c Case for Alexeev dynamic mag model 2000 (use of GSM coordinates)
c inputs is SW pressure, (nPa), DST (nT), ByIMF and BzIMF (nT), G1,G2
c
       if (k_ext .eq. 12) then
        IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.50.D0) THEN
           Ifail=-1
           RETURN
        ENDIF
          den=density
           vel=speed
        dst=dst_nt
        BIMF(1)=BxIMF_nt
        BIMF(2)=ByIMF_nt
        BIMF(3)=BzIMF_nt
        Al_ind=Al
          CALL GEO_GSM(xGEO,xSM)
          CALL a2000(den,vel,BIMF,
     &    dst,al_ind,xSM,BxSM,ifail)
          if (Ifail .eq. 0) CALL GSM_GEO(BxSM,Bxext)
       endif
c
C
c Case for Tsyganenko and Sitnov 07d mag model (use of GSM coordinates)
c input is Pdyn
c
c      2015 version, new bess calls in external file TS07j_2015.f (previously TS07j.f)
       IF (k_ext .eq. 13) THEN 
c M. Sitnov and G. Stephens recommend limiting the valid output to 20 Re and
c less based on the current (at the time) expansion:
       IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.400.D0) THEN 
             Ifail=-1
             RETURN
          ENDIF

          CALL GEO_GSM(xGEO,xSM)
          CALL TS07D_2015 (xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),
     &      BxSM(3)) ! note that no solar wind input is necessary as it is
                    ! included through the TS07D common block
          CALL GSM_GEO(BxSM,Bxext)
       endif

c      2017 version, new bess calls implemented directly in Tsy_and_Sit_Jul2017.f code
c      Code is double precision by design
       IF (k_ext .eq. 14) THEN 
c M. Sitnov and G. Stephens recommend limiting the valid output to 20 Re and
c less based on the current (at the time) expansion:
       IF (xGEO(1)*xGEO(1)+xGEO(2)*xGEO(2)+xGEO(3)*xGEO(3)
     &    .GT.400.D0) THEN 
             Ifail=-1
             RETURN
          ENDIF

          CALL GEO_GSM(xGEO,xSM)
          CALL TS07D_JULY_2017 (xSM(1),xSM(2),xSM(3),BxSM(1),BxSM(2),
     &      BxSM(3)) ! note that no solar wind input is necessary as it is
                    ! included through the TS07D common block
          CALL GSM_GEO(BxSM,Bxext)
       endif

       Bl = 0.D0
       DO I = 1,3
         !WRITE(6,*)Bxext(I),Bxint(I)
         BxGEO(I) = Bxext(I) + Bxint(I)
         Bl = Bl + BxGEO(I)*BxGEO(I)
       ENDDO
C
       Bl = SQRT(Bl)
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE INIT_DTD(year)
C
       IMPLICIT NONE
C
       INTEGER*4 ierr
       REAL*8    year
       REAL*8    g(66),h(66)
       REAL*8    xc,yc,zc                  !Re
       REAL*8    thet,phit           !radian
       REAL*8    ct,st,cp,sp
       REAL*8    Bo
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
       COMMON /dgrf/g,h
C
       call get_igrf_coeffs(year,g,h,ierr)
c
       call get_terms(g,h,thet,phit,xc,yc,zc,Bo)
c       write(6,*)Bo
c
       CALL calcul_GH1
C
       ct = COS(thet)
       st = SIN(thet)
       cp = COS(phit)
       sp = SIN(phit)
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE INIT_GSM(iyr,iday,secs,psi)
C
       IMPLICIT NONE
C
       INTEGER*4   iyr,iday
       REAL*8      secs
       REAL*8      gst,Slong,Srasn,Sdec
       REAL*8      xS,yS,zS                     !GEI
       REAL*8      xDip,yDip,zDip               !GEO
       REAL*8      xD,yD,zD                     !GEI
       REAL*8      xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
       REAL*8      xc,yc,zc
       REAL*8      ct,st,cp,sp
       REAL*8      cgst,sgst
       REAL*8      Bo
       REAL*8      norm
       REAL*8      psi
C
       COMMON /Soleil/xS,yS,zS,cgst,sgst
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
       COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
       CALL SUN(iyr,iday,secs,gst,Slong,Srasn,Sdec)
C
       xS = COS(Srasn)*COS(Sdec)
       yS = SIN(Srasn)*COS(Sdec)
       zS = SIN(Sdec)
C
       cgst = COS(gst)
       sgst = SIN(gst)
C
       xDip = st*cp
       yDip = st*sp
       zDip = ct
C
       xD = cgst*xDip - sgst*yDip
       yD = sgst*xDip + cgst*yDip
       zD = zDip
C
       xSD = yD*zS - yS*zD
       ySD = zD*xS - zS*xD
       zSD = xD*yS - xS*yD
C
       norm =       SQRT(xSD*xSD + ySD*ySD + zSD*zSD)
C
       xSD = xSD/norm
       ySD = ySD/norm
       zSD = zSD/norm
C
       xSSD = yS*zSD - ySD*zS
       ySSD = zS*xSD - zSD*xS
       zSSD = xS*ySD - xSD*yS
C
       norm = SQRT(xSSD*xSSD + ySSD*ySSD + zSSD*zSSD)
C
       xSSD = xSSD/norm
       ySSD = ySSD/norm
       zSSD = zSSD/norm
C
       xSDD = ySD*zD - yD*zSD
       ySDD = zSD*xD - zD*xSD
       zSDD = xSD*yD - xD*ySD
C
       norm = SQRT(xSDD*xSDD + ySDD*ySDD + zSDD*zSDD)
C
       xSDD = xSDD/norm
       ySDD = ySDD/norm
       zSDD = zSDD/norm
C
       psi = ASIN(xD*xS + yD*yS + zD*zS)
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       GEO_GSM(xGEO,xGSM)
C
      IMPLICIT NONE
C
      REAL*8    xGEO(3)
      REAL*8    xGSM(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
        COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      xGEI =  cgst*xGEO(1) - sgst*xGEO(2)
      yGEI =  sgst*xGEO(1) + cgst*xGEO(2)
        zGEI =  xGEO(3)
C
      xGSM(1) = xS*xGEI + yS*yGEI + zS*zGEI
      xGSM(2) = xSD*xGEI + ySD*yGEI + zSD*zGEI
        xGSM(3) = xSSD*xGEI + ySSD*yGEI + zSSD*zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       GSM_GEO(xGSM,xGEO)
C
      IMPLICIT NONE
C
      REAL*8    xGSM(3)
      REAL*8    xGEO(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
        COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      xGEI = xS*xGSM(1) + xSD*xGSM(2) + xSSD*xGSM(3)
      yGEI = yS*xGSM(1) + ySD*xGSM(2) + ySSD*xGSM(3)
        zGEI = zS*xGSM(1) + zSD*xGSM(2) + zSSD*xGSM(3)
C
      xGEO(1) =  cgst*xGEI + sgst*yGEI
      xGEO(2) = -sgst*xGEI + cgst*yGEI
        xGEO(3) =  zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       GEO_GEI(xGEO,xGEI)
C
      IMPLICIT NONE
C
      REAL*8    xGEO(3)
      REAL*8    xGEI(3)
        REAL*8    xS,yS,zS,cgst,sgst
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
C
      xGEI(1) =  cgst*xGEO(1) - sgst*xGEO(2)
      xGEI(2) =  sgst*xGEO(1) + cgst*xGEO(2)
        xGEI(3) =  xGEO(3)
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       GEI_GEO(xGEI,xGEO)
C
      IMPLICIT NONE
C
      REAL*8    xGEI(3)
      REAL*8    xGEO(3)
        REAL*8    xS,yS,zS,cgst,sgst
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
C
      xGEO(1) =  cgst*xGEI(1) + sgst*xGEI(2)
      xGEO(2) = -sgst*xGEI(1) + cgst*xGEI(2)
        xGEO(3) =  xGEI(3)
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       GEO_SM(xGEO,xSM)
C
      IMPLICIT NONE
C
      REAL*8    xGEO(3)
      REAL*8    xSM(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
        COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      xGEI =  cgst*xGEO(1) - sgst*xGEO(2)
      yGEI =  sgst*xGEO(1) + cgst*xGEO(2)
        zGEI =  xGEO(3)
C
      xSM(1) = xSDD*xGEI + ySDD*yGEI + zSDD*zGEI
      xSM(2) = xSD*xGEI  + ySD*yGEI  + zSD*zGEI
        xSM(3) = xD*xGEI   + yD*yGEI   + zD*zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       SM_GEO(xSM,xGEO)
C
      IMPLICIT NONE
C
      REAL*8    xSM(3)
      REAL*8    xGEO(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
        COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      xGEI = xSDD*xSM(1) + xSD*xSM(2) + xD*xSM(3)
      yGEI = ySDD*xSM(1) + ySD*xSM(2) + yD*xSM(3)
        zGEI = zSDD*xSM(1) + zSD*xSM(2) + zD*xSM(3)
C
      xGEO(1) =  cgst*xGEI + sgst*yGEI
      xGEO(2) = -sgst*xGEI + cgst*yGEI
        xGEO(3) =  zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       SM_GSM(xSM,xGSM)
C
      IMPLICIT NONE
C
      REAL*8    xSM(3)
      REAL*8    xGSM(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
        COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      xGEI = xSDD*xSM(1) + xSD*xSM(2) + xD*xSM(3)
      yGEI = ySDD*xSM(1) + ySD*xSM(2) + yD*xSM(3)
        zGEI = zSDD*xSM(1) + zSD*xSM(2) + zD*xSM(3)
C
      xGSM(1) = xS*xGEI + yS*yGEI + zS*zGEI
      xGSM(2) = xSD*xGEI + ySD*yGEI + zSD*zGEI
        xGSM(3) = xSSD*xGEI + ySSD*yGEI + zSSD*zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE
     &       GSM_SM(xGSM,xSM)
C
      IMPLICIT NONE
C
      REAL*8    xGSM(3)
      REAL*8    xSM(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      COMMON /Soleil/xS,yS,zS,cgst,sgst
        COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
      xGEI = xS*xGSM(1) + xSD*xGSM(2) + xSSD*xGSM(3)
      yGEI = yS*xGSM(1) + ySD*xGSM(2) + ySSD*xGSM(3)
        zGEI = zS*xGSM(1) + zSD*xGSM(2) + zSSD*xGSM(3)
C
      xSM(1) = xSDD*xGEI + ySDD*yGEI + zSDD*zGEI
      xSM(2) = xSD*xGEI  + ySD*yGEI  + zSD*zGEI
        xSM(3) = xD*xGEI   + yD*yGEI   + zD*zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE GEO_MAG(xGEO,xMAG)
C
       IMPLICIT NONE
C
       REAL*8    xGEO(3)
       REAL*8    xMAG(3)
       REAL*8    xc,yc,zc                  !Re
       REAL*8    ct,st,cp,sp
       REAL*8    Bo
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C
       xMAG(1) =  xGEO(1)*ct*cp + xGEO(2)*ct*sp - xGEO(3)*st
       xMAG(2) = -xGEO(1)*sp    + xGEO(2)*cp
       xMAG(3) =  xGEO(1)*st*cp + xGEO(2)*st*sp + xGEO(3)*ct
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE MAG_GEO(xMAG,xGEO)
C
       IMPLICIT NONE
C
       REAL*8    xGEO(3)
       REAL*8    xMAG(3)
       REAL*8    xc,yc,zc                  !Re
       REAL*8    ct,st,cp,sp
       REAL*8    Bo
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C
       xGEO(1) =  xMAG(1)*ct*cp - xMAG(2)*sp + xMAG(3)*st*cp
       xGEO(2) =  xMAG(1)*ct*sp + xMAG(2)*cp + xMAG(3)*st*sp
       xGEO(3) = -xMAG(1)*st                 + xMAG(3)*ct
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE GEO_GSE(xGEO,xGSE)
C
       IMPLICIT NONE
C
      REAL*8    xGSE(3)
      REAL*8    xGEO(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
      REAL*8    aa,bb,y1,y2,y3
C
       COMMON /Soleil/xS,yS,zS,cgst,sgst
       COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
        aa = -0.3978D0
      bb = 0.9175D0
C
        y1 = aa*zS-bb*yS
      y2 = bb*xS
      y3 = -aa*xS
C
      xGEI =  cgst*xGEO(1) - sgst*xGEO(2)
      yGEI =  sgst*xGEO(1) + cgst*xGEO(2)
        zGEI =  xGEO(3)
C
      xGSE(1) = xS*xGEI + yS*yGEI + zS*zGEI
      xGSE(2) = y1*xGEI + y2*yGEI + y3*zGEI
        xGSE(3) = aa*yGEI + bb*zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE GSE_GEO(xGSE,xGEO)
C
       IMPLICIT NONE
C
      REAL*8    xGSE(3)
      REAL*8    xGEO(3)
      REAL*8    xGEI,yGEI,zGEI
        REAL*8    xS,yS,zS,cgst,sgst
        REAL*8    xD,yD,zD                     !GEI
        REAL*8    xSD,ySD,zSD,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
      REAL*8    aa,bb,y1,y2,y3,det
C
       COMMON /Soleil/xS,yS,zS,cgst,sgst
       COMMON /sundip/xD,yD,zD,xSD,ySD,zSD
     &      ,xSSD,ySSD,zSSD,xSDD,ySDD,zSDD
C
        aa = -0.3978D0
      bb = 0.9175D0
C
        y1 = aa*zS-bb*yS
      y2 = bb*xS
      y3 = -aa*xS
C
        det = xS*y2*bb+zS*y1*aa-xS*y3*aa-yS*y1*bb
C
      xGEI =  (y2*bb-y3*aa)*xGSE(1) - (ys*bb-zS*aa)*xGSE(2)
     &        + (yS*y3-zS*y2)*xGSE(3)
      yGEI = - y1*bb*xGSE(1) + xS*bb*xGSE(2) - (xS*y3-zS*y1)*xGSE(3)
        zGEI =   y1*aa*xGSE(1) - xS*aa*xGSE(2) + (xS*y2-yS*y1)*xGSE(3)
C
        xGEI = xGEI/det
        yGEI = yGEI/det
        zGEI = zGEI/det
C
      xGEO(1) =  cgst*xGEI + sgst*yGEI
      xGEO(2) = -sgst*xGEI + cgst*yGEI
        xGEO(3) =  zGEI
C
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE SPH_CAR(r,lati,longi,x)
C
       IMPLICIT NONE
C
       REAL*8     r,lati,longi
       REAL*8     colati
       REAL*8     x(3)
        REAL*8     pi,rad
        common /rconst/rad,pi
C
        CALL INITIZE
       colati = pi/2.d0 - lati*rad
C
       x(1) = r*SIN(colati)*COS(longi*rad)
       x(2) = r*SIN(colati)*SIN(longi*rad)
       x(3) = r*COS(colati)
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE CAR_SPH(x,r,lati,longi)
C
       IMPLICIT NONE
       REAL*8     r,lati,longi
       REAL*8     SQ
       REAL*8     x(3)
        REAL*8     pi,rad
        common /rconst/rad,pi
C
        CALL INITIZE
       SQ = x(1)*x(1) + x(2)*x(2)
       r = SQRT(SQ + x(3)*x(3))
       IF (SQ.NE.0.d0) GOTO 2
       longi = 0.d0
       IF (x(3).LT.0.d0) THEN
          lati = -90.d0
       ELSE
          lati = 90.d0
       ENDIF
       RETURN
2      SQ = SQRT(SQ)
       longi = ATAN2(x(2),x(1))/rad
       lati = 90.d0 - ATAN2(SQ,x(3))/rad
       IF (longi.LT.0.d0) longi = longi + 360.d0
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE SPH_CAR_VECT(r,lati,longi,sph,x)
C
       IMPLICIT NONE
C
       REAL*8     r,lati,longi
       REAL*8     colati
       REAL*8     st,ct,sp,cp
       REAL*8     sph(3),x(3)
        REAL*8     pi,rad
        common /rconst/rad,pi
C
        CALL INITIZE
       colati = pi/2.d0 - lati*rad
C
      st     = sin(colati)
      ct     = cos(colati)
      sp     = sin(longi*rad)
      cp     = cos(longi*rad)
      x(1) = sph(1) * st*cp + sph(2) * ct*cp - sph(3) * sp
      x(2) = sph(1) * st*sp + sph(2) * ct*sp + sph(3) * cp
      x(3) = sph(1) *  ct   - sph(2) *  st
C
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE CAR_SPH_VECT(r,lati,longi,x,sph)
C
       IMPLICIT NONE
C
       REAL*8     r,lati,longi
       REAL*8     colati
       REAL*8     st,ct,sp,cp
       REAL*8     sph(3),x(3)
        REAL*8     pi,rad
        common /rconst/rad,pi
C
        CALL INITIZE
       colati = pi/2.d0 - lati*rad
C
      st     = sin(colati)
      ct     = cos(colati)
      sp     = sin(longi*rad)
      cp     = cos(longi*rad)
      sph(1) = x(3) *  ct + x(1) * st*cp + x(2) * st*sp
      sph(2) = -st * x(3) + x(1) * ct*cp + x(2) * ct*sp
      sph(3) =                - sp * x(1)  + x(2) *  cp
C
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SUN(iyr,iday,secs,gst,Slong,Srasn,Sdec)
C
      IMPLICIT NONE
C
C program to calculate sideral, time and position of the Sun
C good for years 1901 through 2099. accuracy 0.006 degree
C input is iyr,iday (INTEGER*4s), and secs, defining universal time
C output is Greenwich mean sidereal time (gst) in degrees,
C longitude along ecliptic (Slong), and apparent right ascension
C and declination (Srasn,Sdec) of the Sun, all in degrees.
C
      INTEGER*4 iyr,iday
      INTEGER*4 iaux
      REAL*8    secs,gst,Slong,Srasn,Sdec
C
        REAL*8 aux
        REAL*8 dj,fday
      REAL*8 t,vl,g,obliq,slp,sind,cosd
        REAL*8     pi,rad
        common /rconst/rad,pi
C
      IF (iyr.LT.1901 .OR. iyr.GT.2099) RETURN
C
      fday = secs/86400.D0
      dj   = 365.D0*(iyr-1900.D0) + INT((iyr-1901)/4)
     &       + iday + fday - 0.5D0
        t    = dj/36525.D0
      aux = 279.696678D0 + 0.9856473354D0*dj
      iaux = INT(aux/360.D0)
      vl = aux-360.D0*iaux
      aux = 279.690983D0 + 0.9856473354D0*dj + 360.D0*fday + 180.D0
      iaux = INT(aux/360.D0)
      gst = aux-360.D0*iaux
      aux = 358.475845D0 + 0.985600267D0*dj
      iaux = INT(aux/360.D0)
      g = (aux-360.D0*iaux)*rad
C
        Slong = vl + (1.91946D0-0.004789D0*t)*SIN(g)
     &         + 0.020094D0*SIN(2.D0*g)
      obliq = (23.45229D0 - 0.0130125D0*t) * rad
        slp   = (Slong - 0.005686D0) * rad
      sind  = SIN(obliq) * SIN(slp)
      cosd  = SQRT(1.D0 - sind*sind)
        Sdec  = ATAN(sind/cosd) / rad
        Srasn = 180.D0
     &    - ATAN2(sind/(cosd*TAN(obliq)), -COS(slp)/cosd) / rad
C
      gst   = gst   * rad
      Slong = Slong * rad
        Sdec  = Sdec  * rad
        Srasn = Srasn * rad
C
        RETURN
        END
C
!C----------------------------------------------------------------------
!C
!      subroutine get_igrf_coeffs(year,g,h,ierr)
!c
!C      SET UP TO ACCEPT DATES BETWEEN 1965 AND 2005; COEFFICIENTS
!C      THROUGH 1985 ARE FROM DGRF MODELS COEFFICIENTS FOR 1990
!C       FROM IGRF
!C      [EOS TRANS. AGU APRIL 21, 1992, P. 182]. INTERPOLATION IS USED
!C      FOR YEARS BETWEEN DGRF MODELS, AND EXTRAPOLATION FOR YEARS
!C      BETWEEN 1990 AND 2000.
!C       Coefficients for 1995 were added by D. Boscher
!C       Modified by Paul O'Brien: ierr set, but uses nearest date limit
!c       Date limits are 1965 to 2015
!c
!c            Mauricio Peredo
!c            Hughes STX at NASA/GSFC
!c            June 29, 1995
!c               bidouille Daniel Boscher Fev 20,1997
!c               bidouille Daniel Boscher May 02,2001
!c
!C  Additional IGRF References:
!C  (J. GEOMAG. GEOELECTR.(1982), V.34, P.313-315,
!C  GEOMAGN. AND AERONOMY (1986), V.26, P.523-525).
!C
!C------INPUT PARAMETERS:
!C  iyr - YEAR NUMBER (FROM 1965 UP TO 2000)
!C----- OUTPUT PARAMETERS:
!C  g,h - coefficients for the igrf model interpolated (or extrapolated)
!c       to the epoch of iyr
!C
!C
!      IMPLICIT NONE
!C
!      REAL*8 F1,F2,DT,G(66),H(66),G65(66),H65(66),G70(66),H70(66),
!     *G75(66),H75(66),G80(66),H80(66),G85(66),H85(66),G90(66),
!     *H90(66),G95(66),H95(66),G00(66),H00(66),G05(66),H05(66),
!     *G2010(66),H2010(66),
!     *DG2010(45),DH2010(45)
!
!      INTEGER*4 N,ierr
!        REAL*8    year
!        INTEGER*4 year_error_reported ! only report year error once
!        save year_error_reported
!        DATA year_error_reported/0/ ! initially not reported
!c
!      DATA G65/0.D0,-30334.D0,-2119.D0,-1662.D0,2997.D0,1594.D0,
!     *1297.D0,-2038.D0,1292.D0,856.D0,957.D0,804.D0,479.D0,-390.D0,
!     *252.D0,-219.D0,358.D0,254.D0,-31.D0,-157.D0,-62.D0,
!     *45.D0,61.D0,8.D0,-228.D0,4.D0,1.D0,-111.D0,75.D0,-57.D0,4.D0,
!     *13.D0,-26.D0,-6.D0,13.D0,1.D0,13.D0,5.D0,-4.D0,-14.D0,0.D0,
!     *8.D0,-1.D0,11.D0,4.D0,8.D0,10.D0,2.D0,-13.D0,10.D0,-1.D0,-1.D0,
!     *5.D0,1.D0,-2.D0,
!     *-2.D0,-3.D0,2.D0,-5.D0,-2.D0,4.D0,4.D0,0.D0,2.D0,2.D0,0.D0/
!      DATA H65/0.D0,0.D0,5776.D0,0.D0,-2016.D0,114.D0,0.D0,-404.D0,
!     *240.D0,-165.D0,0.D0,148.D0,-269.D0,13.D0,-269.D0,0.D0,19.D0,
!     *128.D0,-126.D0,-97.D0,81.D0,0.D0,-11.D0,100.D0,68.D0,-32.D0,
!     *-8.D0,-7.D0,0.D0,-61.D0,-27.D0,-2.D0,6.D0,26.D0,-23.D0,-12.D0,
!     *0.D0,7.D0,-12.D0,9.D0,-16.D0,4.D0,24.D0,-3.D0,-17.D0,0.D0,-22.D0,
!     *15.D0,7.D0,-4.D0,-5.D0,10.D0,10.D0,-4.D0,1.D0,0.D0,2.D0,1.D0,
!     *2.D0,6.D0,-4.D0,0.D0,-2.D0,3.D0,0.D0,-6.D0/
!c
!      DATA G70/0.D0,-30220.D0,-2068.D0,-1781.D0,3000.D0,1611.D0,
!     *1287.D0,-2091.D0,1278.D0,838.D0,952.D0,800.D0,461.D0,-395.D0,
!     *234.D0,-216.D0,359.D0,262.D0,-42.D0,-160.D0,-56.D0,
!     *43.D0,64.D0,15.D0,-212.D0,2.D0,3.D0,-112.D0,72.D0,-57.D0,1.D0,
!     *14.D0,-22.D0,-2.D0,13.D0,-2.D0,14.D0,6.D0,-2.D0,-13.D0,-3.D0,
!     *5.D0,0.D0,11.D0,3.D0,8.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,0.D0,
!     *3.D0,1.D0,-1.D0,-3.D0,-3.D0,2.D0,-5.D0,-1.D0,6.D0,4.D0,1.D0,
!     *0.D0,3.D0,-1.D0/
!      DATA H70/0.D0,0.D0,5737.D0,0.D0,-2047.D0,25.D0,0.D0,-366.D0,
!     *251.D0,-196.D0,0.D0,167.D0,-266.D0,26.D0,-279.D0,0.D0,26.D0,
!     *139.D0,-139.D0,-91.D0,83.D0,0.D0,-12.D0,100.D0,72.D0,-37.D0,
!     *-6.D0,1.D0,0.D0,-70.D0,-27.D0,-4.D0,8.D0,23.D0,-23.D0,-11.D0,
!     *0.D0,7.D0,-15.D0,6.D0,-17.D0,6.D0,21.D0,-6.D0,-16.D0,0.D0,
!     *-21.D0,16.D0,6.D0,-4.D0,-5.D0,10.D0,11.D0,-2.D0,1.D0,0.D0,
!     *1.D0,1.D0,3.D0,4.D0,-4.D0,0.D0,-1.D0,3.D0,1.D0,-4.D0/
!c
!      DATA G75/0.D0,-30100.D0,-2013.D0,-1902.D0,3010.D0,1632.D0,
!     *1276.D0,-2144.D0,1260.D0,830.D0,946.D0,791.D0,438.D0,-405.D0,
!     *216.D0,-218.D0,356.D0,264.D0,-59.D0,-159.D0,-49.D0,
!     *45.D0,66.D0,28.D0,-198.D0,1.D0,6.D0,-111.D0,71.D0,-56.D0,1.D0,
!     *16.D0,-14.D0,0.D0,12.D0,-5.D0,14.D0,6.D0,-1.D0,-12.D0,-8.D0,
!     *4.D0,0.D0,10.D0,1.D0,7.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,-1.D0,
!     *4.D0,1.D0,-2.D0,-3.D0,-3.D0,2.D0,-5.D0,-2.D0,5.D0,4.D0,1.D0,
!     *0.D0,3.D0,-1.D0/
!      DATA H75/0.D0,0.D0,5675.D0,0.D0,-2067.D0,-68.D0,0.D0,-333.D0,
!     *262.D0,-223.D0,0.D0,191.D0,-265.D0,39.D0,-288.D0,0.D0,31.D0,
!     *148.D0,-152.D0,-83.D0,88.D0,0.D0,-13.D0,99.D0,75.D0,-41.D0,
!     *-4.D0,11.D0,0.D0,-77.D0,-26.D0,-5.D0,10.D0,22.D0,-23.D0,-12.D0,
!     *0.D0,6.D0,-16.D0,4.D0,-19.D0,6.D0,18.D0,-10.D0,-17.D0,0.D0,
!     *-21.D0,16.D0,7.D0,-4.D0,-5.D0,10.D0,11.D0,-3.D0,1.D0,0.D0,1.D0,
!     *1.D0,3.D0,4.D0,-4.D0,-1.D0,-1.D0,3.D0,1.D0,-5.D0/
!c
!      DATA G80/0.D0,-29992.D0,-1956.D0,-1997.D0,3027.D0,1663.D0,
!     *1281.D0,-2180.D0,1251.D0,833.D0,938.D0,782.D0,398.D0,-419.D0,
!     *199.D0,-218.D0,357.D0,261.D0,-74.D0,-162.D0,-48.D0,48.D0,66.D0,
!     *42.D0,-192.D0,4.D0,14.D0,-108.D0,72.D0,-59.D0,2.D0,21.D0,-12.D0,
!     *1.D0,11.D0,-2.D0,18.D0,6.D0,0.D0,-11.D0,-7.D0,4.D0,3.D0,6.D0,
!     *-1.D0,5.D0,10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,2.D0,
!     *-5.D0,-4.D0,-4.D0,2.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,3.D0,0.D0/
!      DATA H80/0.D0,0.D0,5604.D0,0.D0,-2129.D0,-200.D0,0.D0,-336.D0,
!     *271.D0,-252.D0,0.D0,212.D0,-257.D0,53.D0,-297.D0,0.D0,46.D0,
!     *150.D0,-151.D0,-78.D0,92.D0,0.D0,-15.D0,93.D0,71.D0,-43.D0,
!     *-2.D0,17.D0,0.D0,-82.D0,-27.D0,-5.D0,16.D0,18.D0,-23.D0,-10.D0,
!     *0.D0,7.D0,-18.D0,4.D0,-22.D0,9.D0,16.D0,-13.D0,-15.D0,0.D0,
!     *-21.D0,16.D0,9.D0,-5.D0,-6.D0,9.D0,10.D0,-6.D0,2.D0,0.D0,1.D0,
!     *0.D0,3.D0,6.D0,-4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0/
!c
!      DATA G85/0.D0,-29873.D0,-1905.D0,-2072.D0,3044.D0,1687.D0,
!     *1296.D0,-2208.D0,1247.D0,829.D0,936.D0,780.D0,361.D0,-424.D0,
!     *170.D0,-214.D0,355.D0,253.D0,-93.D0,-164.D0,-46.D0,
!     *53.D0,65.D0,51.D0,-185.D0,4.D0,16.D0,-102.D0,74.D0,-62.D0,3.D0,
!     *24.D0,-6.D0,4.D0,10.D0,0.D0,21.D0,6.D0,0.D0,-11.D0,-9.D0,4.D0,
!     *4.D0,4.D0,-4.D0,5.D0,10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,
!     *1.D0,-5.D0,-4.D0,-4.D0,3.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,
!     *3.D0,0.D0/
!      DATA H85/0.D0,0.D0,5500.D0,0.D0,-2197.D0,-306.D0,0.D0,-310.D0,
!     *284.D0,-297.D0,0.D0,232.D0,-249.D0,69.D0,-297.D0,0.D0,47.D0,
!     *150.D0,-154.D0,-75.D0,95.D0,0.D0,-16.D0,88.D0,69.D0,-48.D0,
!     *-1.D0,21.D0,0.D0,-83.D0,-27.D0,-2.D0,20.D0,17.D0,-23.D0,-7.D0,
!     *0.D0,8.D0,-19.D0,5.D0,-23.D0,11.D0,14.D0,-15.D0,-11.D0,0.D0,
!     *-21.D0,15.D0,9.D0,-6.D0,-6.D0,9.D0,9.D0,-7.D0,2.D0,0.D0,1.D0,
!     *0.D0,3.D0,6.D0,-4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0/
!c
!      DATA G90/0.D0,-29775.D0, -1848.D0, -2131.D0,  3059.D0,  1686.D0,
!     *      1314.D0, -2239.D0,  1248.D0,   802.D0,   939.D0,   780.D0,
!     *       325.D0,  -423.D0,   141.D0,  -214.D0,   353.D0,   245.D0,
!     *      -109.D0,  -165.D0,   -36.D0,    61.D0,    65.D0,    59.D0,
!     *      -178.D0,     3.D0,    18.D0,   -96.D0,    77.D0,   -64.D0,
!     *         2.D0,    26.D0,    -1.D0,     5.D0,     9.D0,     0.D0,
!     *        23.D0,     5.D0,    -1.D0,   -10.D0,   -12.D0,     3.D0,
!     *         4.D0,     2.D0,    -6.D0,     4.D0,     9.D0,     1.D0,
!     *       -12.D0,     9.D0,    -4.D0,    -2.D0,     7.D0,     1.D0,
!     *        -6.D0,    -3.D0,    -4.D0,     2.D0,    -5.D0,    -2.D0,
!     *         4.D0,     3.D0,     1.D0,     3.D0,     3.D0,     0.D0/
!      DATA H90/0.D0,     0.D0,  5406.D0,     0.D0, -2279.D0,  -373.D0,
!     *         0.D0,  -284.D0,   293.D0,  -352.D0,     0.D0,   247.D0,
!     *      -240.D0,    84.D0,  -299.D0,     0.D0,    46.D0,   154.D0,
!     *      -153.D0,   -69.D0,    97.D0,     0.D0,   -16.D0,    82.D0,
!     *        69.D0,   -52.D0,     1.D0,    24.D0,     0.D0,   -80.D0,
!     *       -26.D0,     0.D0,    21.D0,    17.D0,   -23.D0,    -4.D0,
!     *         0.D0,    10.D0,   -19.D0,     6.D0,   -22.D0,    12.D0,
!     *        12.D0,   -16.D0,   -10.D0,     0.D0,   -20.D0,    15.D0,
!     *        11.D0,    -7.D0,    -7.D0,     9.D0,     8.D0,    -7.D0,
!     *         2.D0,     0.D0,     2.D0,     1.D0,     3.D0,     6.D0,
!     *        -4.D0,     0.D0,    -2.D0,     3.D0,    -1.D0,    -6.D0/
!C
!      DATA G95/0.D0,-29692.D0, -1784.D0, -2200.D0,  3070.D0,  1681.D0,
!     *      1335.D0, -2267.D0,  1249.D0,   759.D0,   940.D0,   780.D0,
!     *       290.D0,  -418.D0,   122.D0,  -214.D0,   352.D0,   235.D0,
!     *      -118.D0,  -166.D0,   -17.D0,    68.D0,    67.D0,    68.D0,
!     *      -170.D0,    -1.D0,    19.D0,   -93.D0,    77.D0,   -72.D0,
!     *         1.D0,    28.D0,     5.D0,     4.D0,     8.D0,    -2.D0,
!     *        25.D0,     6.D0,    -6.D0,    -9.D0,   -14.D0,     9.D0,
!     *         6.D0,    -5.D0,    -7.D0,     4.D0,     9.D0,     3.D0,
!     *       -10.D0,     8.D0,    -8.D0,    -1.D0,    10.D0,    -2.D0,
!     *        -8.D0,    -3.D0,    -6.D0,     2.D0,    -4.D0,    -1.D0,
!     *         4.D0,     2.D0,     2.D0,     5.D0,     1.D0,     0.D0/
!      DATA H95/0.D0,     0.D0,  5306.D0,     0.D0, -2366.D0,  -413.D0,
!     *         0.D0,  -262.D0,   302.D0,  -427.D0,     0.D0,   262.D0,
!     *      -236.D0,    97.D0,  -306.D0,     0.D0,    46.D0,   165.D0,
!     *      -143.D0,   -55.D0,   107.D0,     0.D0,   -17.D0,    72.D0,
!     *        67.D0,   -58.D0,     1.D0,    36.D0,     0.D0,   -69.D0,
!     *       -25.D0,     4.D0,    24.D0,    17.D0,   -24.D0,    -6.D0,
!     *         0.D0,    11.D0,   -21.D0,     8.D0,   -23.D0,    15.D0,
!     *        11.D0,   -16.D0,    -4.D0,     0.D0,   -20.D0,    15.D0,
!     *        12.D0,    -6.D0,    -8.D0,     8.D0,     5.D0,    -8.D0,
!     *         3.D0,     0.D0,     1.D0,     0.D0,     4.D0,     5.D0,
!     *        -5.D0,    -1.D0,    -2.D0,     1.D0,    -2.D0,    -7.D0/
!C
!      DATA G00/0.0D0,-29619.4D0,-1728.2D0,-2267.7D0,3068.4D0,1670.9D0,
!     *      1339.6D0, -2288.0D0, 1252.1D0,  714.5D0, 932.3D0, 786.8D0,
!     *       250.0D0,  -403.0D0,  111.3D0, -218.8D0, 351.4D0, 222.3D0,
!     *      -130.4D0,  -168.6D0,  -12.9D0,   72.3D0,  68.2D0,  74.2D0,
!     *      -160.9D0,    -5.9D0,   16.9D0,  -90.4D0,  79.0D0, -74.0D0,
!     *         0.0D0,    33.3D0,    9.1D0,    6.9D0,   7.3D0,  -1.2D0,
!     *        24.4D0,     6.6D0,   -9.2D0,   -7.9D0, -16.6D0,   9.1D0,
!     *         7.0D0,    -7.9D0,   -7.0D0,    5.0D0,   9.4D0,   3.0D0,
!     *        -8.4D0,     6.3D0,   -8.9D0,   -1.5D0,   9.3D0,  -4.3D0,
!     *        -8.2D0,    -2.6D0,   -6.0D0,    1.7D0,  -3.1D0,  -0.5D0,
!     *         3.7D0,     1.0D0,    2.0D0,    4.2D0,   0.3D0,  -1.1D0/
!      DATA H00/0.0D0,     0.0D0, 5186.1D0,   0.0D0,-2481.6D0,-458.0D0,
!     *         0.0D0,  -227.6D0,  293.4D0, -491.1D0,   0.0D0, 272.6D0,
!     *      -231.9D0,   119.8D0, -303.8D0,    0.0D0,  43.8D0, 171.9D0,
!     *      -133.1D0,   -39.3D0,  106.3D0,    0.0D0, -17.4D0,  63.7D0,
!     *        65.1D0,   -61.2D0,    0.7D0,   43.8D0,   0.0D0, -64.6D0,
!     *       -24.2D0,     6.2D0,   24.0D0,   14.8D0, -25.4D0,  -5.8D0,
!     *         0.0D0,    11.9D0,  -21.5D0,    8.5D0, -21.5D0,  15.5D0,
!     *         8.9D0,   -14.9D0,   -2.1D0,    0.0D0, -19.7D0,  13.4D0,
!     *        12.5D0,    -6.2D0,   -8.4D0,    8.4D0,   3.8D0,  -8.2D0,
!     *         4.8D0,     0.0D0,    1.7D0,    0.0D0,   4.0D0,   4.9D0,
!     *        -5.9D0,    -1.2D0,   -2.9D0,    0.2D0,  -2.2D0,  -7.4D0/
!C
!      DATA G05/0.0D0,-29554.63D0,-1669.05D0,-2337.24D0,3047.69D0,
!     *     1657.76D0,   1336.3D0,-2305.83D0, 1246.39D0, 672.51D0,
!     *      920.55D0,   797.96D0,  210.65D0, -379.86D0, 100.0D0,
!     *       -227.D0,   354.41D0,  208.95D0, -136.54D0,-168.05D0,
!     *      -13.55D0,    73.60D0,   69.56D0,   76.74D0,-151.34D0,
!     *      -14.58D0,    14.58D0,  -86.36D0,   79.88D0, -74.46D0,
!     *       -1.65D0,    38.73D0,    12.3D0,    9.37D0,   5.42D0,
!     *        1.94D0,    24.8D0,    7.62D0,  -11.73D0,   -6.88D0,
!     *      -18.11D0,   10.17D0,    9.36D0,  -11.25D0,   -4.87D0,
!     *        5.58D0,    9.76D0,    3.58D0,   -6.94D0,    5.01D0,
!     *      -10.76D0,   -1.25D0,    8.76D0,   -6.66D0,   -9.22D0,
!     *       -2.17D0,   -6.12D0,    1.42D0,   -2.35D0,   -0.15D0,
!     *        3.06D0,    0.29D0,    2.06D0,    3.77D0,   -0.21D0,
!     *       -2.09D0/
!      DATA H05/0.0D0,    0.0D0, 5077.99D0,     0.0D0, -2594.5D0,
!     *     -515.43D0,    0.0D0, -198.86D0,  269.72D0, -524.72D0,
!     *         0.0D0,  282.07D0,-225.23D0,  145.15D0, -305.36D0,
!     *         0.0D0,   42.72D0,  180.25D0,-123.45D0,  -19.57D0,
!     *      103.85D0,    0.0D0,  -20.33D0,   54.75D0,   63.63D0,
!     *      -63.53D0,    0.24D0,   50.94D0,    0.0D0,  -61.14D0,
!     *      -22.57D0,    6.82D0,   25.35D0,   10.93D0,  -26.32D0,
!     *       -4.64D0,     0.0D0,   11.20D0,  -20.88D0,    9.83D0,
!     *      -19.71D0,   16.22D0,    7.61D0,  -12.76D0,   -0.06D0,
!     *         0.0D0,  -20.11D0,   12.67D0,    12.7D0,   -6.72D0,
!     *       -8.16D0,    8.10D0,    2.92D0,   -7.73D0,    6.01D0,
!     *         0.0D0,    2.19D0,    0.10D0,    4.46D0,    4.76D0,
!     *       -6.58D0,   -1.01D0,   -3.47D0,   -0.86D0,   -2.31D0,
!     *       -7.93D0/
!C
!      DATA G2010/
!     *      0.00d0,  !n=00 m=00-00
!     * -29496.50d0,  -1585.90d0,  !n=01 m=00-01
!     *  -2396.60d0,   3026.00d0,   1668.60d0,  !n=02 m=00-02
!     *   1339.70d0,  -2326.30d0,   1231.70d0,    634.20d0,  !n=03 m=00-03
!     *    912.60d0,    809.00d0,    166.60d0,   -357.10d0,     89.70d0,  !n=04 m=00-04
!     *   -231.10d0,    357.20d0,    200.30d0,   -141.20d0,   -163.10d0,
!     *     -7.70d0,  !n=05 m=00-05
!     *     72.80d0,     68.60d0,     76.00d0,   -141.40d0,    -22.90d0,
!     *     13.10d0,    -77.90d0,  !n=06 m=00-06
!     *     80.40d0,    -75.00d0,     -4.70d0,     45.30d0,     14.00d0,
!     *     10.40d0,      1.60d0,      4.90d0,  !n=07 m=00-07
!     *     24.30d0,      8.20d0,    -14.50d0,     -5.70d0,    -19.30d0,
!     *     11.60d0,     10.90d0,    -14.10d0,     -3.70d0,  !n=08 m=00-08
!     *      5.40d0,      9.40d0,      3.40d0,     -5.30d0,      3.10d0,
!     *    -12.40d0,     -0.80d0,      8.40d0,     -8.40d0,    -10.10d0,  !n=09 m=00-09
!     *     -2.00d0,     -6.30d0,      0.90d0,     -1.10d0,     -0.20d0,
!     *      2.50d0,     -0.30d0,      2.20d0,      3.10d0,     -1.00d0,
!     *     -2.80d0  !n=10 m=00-10
!     * /
!
!      DATA H2010/
!     *      0.00d0,  !n=00 m=00-00
!     *      0.00d0,   4945.10d0,  !n=01 m=00-01
!     *      0.00d0,  -2707.70d0,   -575.40d0,  !n=02 m=00-02
!     *      0.00d0,   -160.50d0,    251.70d0,   -536.80d0,  !n=03 m=00-03
!     *      0.00d0,    286.40d0,   -211.20d0,    164.40d0,   -309.20d0,  !n=04 m=00-04
!     *      0.00d0,     44.70d0,    188.90d0,   -118.10d0,      0.10d0,
!     *    100.90d0,  !n=05 m=00-05
!     *      0.00d0,    -20.80d0,     44.20d0,     61.50d0,    -66.30d0,
!     *      3.10d0,     54.90d0,  !n=06 m=00-06
!     *      0.00d0,    -57.80d0,    -21.20d0,      6.60d0,     24.90d0,
!     *      7.00d0,    -27.70d0,     -3.40d0,  !n=07 m=00-07
!     *      0.00d0,     10.90d0,    -20.00d0,     11.90d0,    -17.40d0,
!     *     16.70d0,      7.10d0,    -10.80d0,      1.70d0,  !n=08 m=00-08
!     *      0.00d0,    -20.50d0,     11.60d0,     12.80d0,     -7.20d0,
!     *     -7.40d0,      8.00d0,      2.20d0,     -6.10d0,      7.00d0,  !n=09 m=00-09
!     *      0.00d0,      2.80d0,     -0.10d0,      4.70d0,      4.40d0,
!     *     -7.20d0,     -1.00d0,     -4.00d0,     -2.00d0,     -2.00d0,
!     *     -8.30d0  !n=10 m=00-10
!     * /
!
!      DATA DG2010/
!     *      0.00d0,  !n=00 m=00-00
!     *     11.40d0,     16.70d0,  !n=01 m=00-01
!     *    -11.30d0,     -3.90d0,      2.70d0,  !n=02 m=00-02
!     *      1.30d0,     -3.90d0,     -2.90d0,     -8.10d0,  !n=03 m=00-03
!     *     -1.40d0,      2.00d0,     -8.90d0,      4.40d0,     -2.30d0,  !n=04 m=00-04
!     *     -0.50d0,      0.50d0,     -1.50d0,     -0.70d0,      1.30d0,
!     *      1.40d0,  !n=05 m=00-05
!     *     -0.30d0,     -0.30d0,     -0.30d0,      1.90d0,     -1.60d0,
!     *     -0.20d0,      1.80d0,  !n=06 m=00-06
!     *      0.20d0,     -0.10d0,     -0.60d0,      1.40d0,      0.30d0,
!     *      0.10d0,     -0.80d0,      0.40d0,  !n=07 m=00-07
!     *     -0.10d0,      0.10d0,     -0.50d0,      0.30d0,     -0.30d0,
!     *      0.30d0,      0.20d0,     -0.50d0,      0.20d0  !n=08 m=00-08
!     * /
!
!      DATA DH2010/
!     *      0.00d0,  !n=00 m=00-00
!     *      0.00d0,    -28.80d0,  !n=01 m=00-01
!     *      0.00d0,    -23.00d0,    -12.90d0,  !n=02 m=00-02
!     *      0.00d0,      8.60d0,     -2.90d0,     -2.10d0,  !n=03 m=00-03
!     *      0.00d0,      0.40d0,      3.20d0,      3.60d0,     -0.80d0,  !n=04 m=00-04
!     *      0.00d0,      0.50d0,      1.50d0,      0.90d0,      3.70d0,
!     *     -0.60d0,  !n=05 m=00-05
!     *      0.00d0,     -0.10d0,     -2.10d0,     -0.40d0,     -0.50d0,
!     *      0.80d0,      0.50d0,  !n=06 m=00-06
!     *      0.00d0,      0.60d0,      0.30d0,     -0.20d0,     -0.10d0,
!     *     -0.80d0,     -0.30d0,      0.20d0,  !n=07 m=00-07
!     *      0.00d0,      0.00d0,      0.20d0,      0.50d0,      0.40d0,
!     *      0.10d0,     -0.10d0,      0.40d0,      0.40d0  !n=08 m=00-08
!     * /
!c
!C
!      if ( (year.lt.1965.D0).or.(year.gt.2015.D0) ) then
!         ierr=1
!           if (year_error_reported .eq. 0) then
!              write(*,999) year
!              year_error_reported = 1
!           endif
!c         goto 300 ! continue, just throw an error
!      endif
!c
!      IF (year.LT.1970.D0) GOTO 50      !INTERPOLATE BETWEEN 1965 - 1970
!      IF (year.LT.1975.D0) GOTO 60      !INTERPOLATE BETWEEN 1970 - 1975
!      IF (year.LT.1980.D0) GOTO 70      !INTERPOLATE BETWEEN 1975 - 1980
!      IF (year.LT.1985.D0) GOTO 80          !INTERPOLATE BETWEEN 1980 - 1985
!      IF (year.LT.1990.D0) GOTO 90      !INTERPOLATE BETWEEN 1985 - 1990
!      IF (year.LT.1995.D0) GOTO 100         !INTERPOLATE BETWEEN 1990 - 1995
!      IF (year.LT.2000.D0) GOTO 110         !INTERPOLATE BETWEEN 1995 - 2000
!      IF (year.LT.2005.D0) GOTO 120         !INTERPOLATE BETWEEN 2000 - 2005
!      IF (year.LT.2010.D0) GOTO 130         !INTERPOLATE BETWEEN 2005 - 2010
!C
!C      EXTRAPOLATE BETWEEN 2010 - 2020
!C
!      DT=year-2010.D0
!      if (DT.gt.10.D0) then
!        DT = 10.D0 ! effectively 2020
!      endif
!      DO 40 N=1,66
!         G(N)=G2010(N)
!         H(N)=H2010(N)
!         IF (N.GT.45) GOTO 40
!         G(N)=G(N)+DG2010(N)*DT
!         H(N)=H(N)+DH2010(N)*DT
!40    CONTINUE
!      GOTO 300
!C
!C
!C      INTERPOLATE BETWEEEN 1965 - 1970
!C
!50    if (year.gt.1965.D0) then
!        F2=(year-1965.D0)/5.D0
!      else
!        F2=0.D0 ! effectively 1965
!      endif
!      F1=1.D0-F2
!      DO 55 N=1,66
!         G(N)=G65(N)*F1+G70(N)*F2
!55       H(N)=H65(N)*F1+H70(N)*F2
!      GOTO 300
!C
!C
!C      INTERPOLATE BETWEEN 1970 - 1975
!C
!60    F2=(year-1970.D0)/5.D0
!      F1=1.D0-F2
!      DO 65 N=1,66
!         G(N)=G70(N)*F1+G75(N)*F2
!65       H(N)=H70(N)*F1+H75(N)*F2
!      GOTO 300
!C
!C
!C      INTERPOLATE BETWEEN 1975 - 1980
!C
!70    F2=(year-1975.D0)/5.D0
!      F1=1.D0-F2
!      DO 75 N=1,66
!         G(N)=G75(N)*F1+G80(N)*F2
!75       H(N)=H75(N)*F1+H80(N)*F2
!      GOTO 300
!C
!C
!C      INTERPOLATE BETWEEN 1980 - 1985
!C
!80    F2=(year-1980.D0)/5.D0
!      F1=1.D0-F2
!      DO 85 N=1,66
!         G(N)=G80(N)*F1+G85(N)*F2
!85       H(N)=H80(N)*F1+H85(N)*F2
!      GOTO 300
!C
!C
!C      INTERPOLATE BETWEEN 1985 - 1990
!C
!90    F2=(year-1985.D0)/5.D0
!      F1=1.D0-F2
!      DO 95 N=1,66
!         G(N)=G85(N)*F1+G90(N)*F2
!95       H(N)=H85(N)*F1+H90(N)*F2
!      GOTO 300
!C
!C
!C      INTERPOLATE BETWEEN 1990 - 1995
!C
!100   F2=(year-1990.D0)/5.
!      F1=1.D0-F2
!      DO 105 N=1,66
!         G(N)=G90(N)*F1+G95(N)*F2
!105      H(N)=H90(N)*F1+H95(N)*F2
!      GOTO 300
!C
!C
!C      INTERPOLATE BETWEEN 1995 - 2000
!C
!110   F2=(year-1995.D0)/5.D0
!      F1=1.D0-F2
!      DO 115 N=1,66
!         G(N)=G95(N)*F1+G00(N)*F2
!115      H(N)=H95(N)*F1+H00(N)*F2
!      GOTO 300
!C
!C      INTERPOLATE BETWEEN 2000 - 2005
!C
!120   F2=(year-2000.D0)/5.D0
!      F1=1.D0-F2
!      DO 125 N=1,66
!         G(N)=G00(N)*F1+G05(N)*F2
!125      H(N)=H00(N)*F1+H05(N)*F2
!      GOTO 300
!C
!C      INTERPOLATE BETWEEN 2005 - 2010
!C
!130   F2=(year-2005.D0)/5.D0
!      F1=1.D0-F2
!      DO 135 N=1,66
!         G(N)=G05(N)*F1+G2010(N)*F2
!135      H(N)=H05(N)*F1+H2010(N)*F2
!      GOTO 300
!C
!C
!300      continue
!c
!C
!999      format(//1x,
!     *   '*** WARNING -- Input year = ',F7.2,/
!     *   ' is out of valid range 1965-2015. Using nearest ***'//)
!c
!      return
!      end
!C
C---------------------------------------------------------------------
C
c
      subroutine get_terms(g,h,alpha,beta,x0,y0,z0,b0)
c
      implicit none
c
      REAL*8 g(66),h(66),alpha,beta,x0,y0,z0
      REAL*8 g10,g11,h11,g20,g21,g22,h21,h22
      REAL*8 b02,b0,sq3,l0,l1,l2,e
c
      g10 = g(2)
      g11 = g(3)
      h11 = h(3)
      g20 = g(4)
      g21 = g(5)
      g22 = g(6)
      h21 = h(5)
      h22 = h(6)
c
      b02 = g10*g10 + g11*g11 + h11*h11
      b0 = sqrt(b02)
c      write(6,*)b0
      alpha = acos(-g10/b0)
      beta = atan(h11/g11)
c
      sq3 = sqrt(3.D0)
      l0 = 2.D0*g10*g20 + sq3*(g11*g21 + h11*h21)
      l1 = -g11*g20 + sq3*(g10*g21 + g11*g22 + h11*h22)
      l2 = -h11*g20 + sq3*(g10*h21 - h11*g22 + g11*h22)
      e = (l0*g10 + l1*g11 + l2*h11)/(4.*b02)
c
      z0 = (l0 - g10*e)/(3.D0*b02)
      x0 = (l1 - g11*e)/(3.D0*b02)
      y0 = (l2 - h11*e)/(3.D0*b02)
!      write(6,*)alpha*180.d0/3.14d0,beta*180.d0/3.14d0
c
      return
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE calcul_GH1
C
      IMPLICIT NONE
C
      INTEGER*4  I,J,L,M,N
      REAL*8     g(66),h(66)
      REAL*8     GH1(144),GH(144)
      REAL*8     X,F,F0
C
      COMMON /MODEL/GH1
      COMMON /dgrf/G,H
C
      L = 0
      M = 0
      DO I = 1,10
        DO J = 0,I
          M = M + 1
          L = L+1
          GH(L) = G(M+1)
          IF(J.NE.0) THEN
            L = L+1
            GH(L) = H(M+1)
          ENDIF
        ENDDO
      ENDDO
C
      GH1(1) =  0.0D0
      I=2
      F0= -1.D0
      DO 9 N=1,10
        X = N
        F0 = F0 * X * X / (4.D0 * X - 2.D0)
        F0 = F0 * (2.D0 * X - 1.D0) / X
        F = F0 * 0.5D0 * SQRT(2.D0)
        GH1(I) = GH(I-1) * F0
        I = I+1
        DO 9 M=1,N
          F = F * (X + M) / (X - M + 1.D0)
          F = F * SQRT((X - M + 1.D0) / (X + M))
          GH1(I) = GH(I-1) * F
          GH1(I+1) = GH(I) * F
          I=I+2
9     CONTINUE
C
      RETURN
      END
C
C*********************************************************************
      SUBROUTINE RLL_GDZ (rr,lati,longi,alti)
C rr adimensionne, alti en km, longi en degres, lati en degres
C
        IMPLICIT NONE
C
        REAL*8 lati,longi
        REAL*8 alti,rr
        REAL*8 CT,ST,CP,SP
        REAL*8 D
      REAL*8 bb,cc
        REAL*8 ERA,AQUAD,BQUAD
C
        COMMON/GENER/ERA,AQUAD,BQUAD
        REAL*8     pi,rad
        common /rconst/rad,pi
C
        CALL INITIZE
C
        CT = DSIN(lati*rad)
        ST = DCOS(lati*rad)
        D  = DSQRT(AQUAD-(AQUAD-BQUAD)*CT*CT)
        CP = COS(longi*rad)
        SP = SIN(longi*rad)
        bb = (BQUAD*CT*CT+AQUAD*ST*ST)/D
        cc = (BQUAD*BQUAD*CT*CT+AQUAD*AQUAD*ST*ST)/D/D-rr*rr*ERA*ERA
        alti = -bb+DSQRT(bb*bb-cc)
C
      RETURN
      END
C
C*********************************************************************
C alti en km
C lati,longi en degres
C x,y,z GEO en adimensionne
C
        SUBROUTINE GDZ_GEO(lati,longi,alti,xx,yy,zz)
C
        IMPLICIT NONE
C
        REAL*8 lati,longi
        REAL*8 alti,xx,yy,zz
        REAL*8 CT,ST,CP,SP
        REAL*8 D,RHO
        REAL*8 ERA,AQUAD,BQUAD
C
        COMMON/GENER/ERA,AQUAD,BQUAD
        REAL*8     pi,rad
        common /rconst/rad,pi
C
        CALL INITIZE
C
        CT = SIN(lati*rad)
        ST = COS(lati*rad)
        D  = SQRT(AQUAD-(AQUAD-BQUAD)*CT*CT)
        CP = COS(longi*rad)
        SP = SIN(longi*rad)
        zz = (alti+BQUAD/D)*CT/ERA
        RHO= (alti+AQUAD/D)*ST/ERA
        xx = RHO*CP
        yy = RHO*SP
C
        RETURN
        END
C
C*********************************************************************
C x,y,z en adimensionne
C lati,longi en degres
C alti en km
C
      SUBROUTINE GEO_GDZ(xx,yy,zz,lati,longi,alti)
C
      IMPLICIT NONE
C
      REAL*8 precision
      PARAMETER (precision = 1.D-5) ! increase precision to 1E-5
C
      REAL*8 lati,longi,lat0
      REAL*8 alti,xx,yy,zz,alt0
      REAL*8 D,RHO
      REAL*8 ERA,AQUAD,BQUAD
      integer*4 i
C
      COMMON/GENER/ERA,AQUAD,BQUAD
      REAL*8     pi,rad
      common /rconst/rad,pi
C
      CALL INITIZE
      longi = ATAN2(yy,xx)/rad
      RHO = SQRT(xx*xx+yy*yy)
      lati = ATAN2(zz,RHO)
      D = COS(lati)
      IF (D.LT.1.D-15) THEN  !trap for pole
        lati = lati/rad
        alti = (abs(zz)-1.D0)*SQRT(BQUAD)
      ELSE
        alti = RHO/D - 1.D0
C
        i=0
10        CONTINUE
        alt0 = alti
        lat0 = lati
        D = SQRT(AQUAD-(AQUAD-BQUAD)*sin(lat0)*sin(lat0))
        lati = ATAN2(zz*(alt0+AQUAD/D/ERA),RHO*(alt0+BQUAD/D/ERA))
        alti = RHO/cos(lati) - AQUAD/D/ERA
        i=i+1
        if (i .gt. 1000) then  !failed to converge
          alti=0.
          lati=0.
        !write(6,*)'geo2gdz ',i,alti
          return
        endif
        IF (ABS(alti - alt0).GT. precision .OR. ABS(lati -
     &     lat0) .GT. precision) GOTO 10 ! keep looping until BOTH are within precision
        alti = alti*ERA
        lati = lati/rad
        ENDIF
        !write(6,*)'geo2gdz ',i,alti
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE INITIZE
C----------------------------------------------------------------
C Initializes the parameters in COMMON/GENER/
C
C      UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
C      ERA      EARTH RADIUS FOR NORMALIZATION OF CARTESIAN
C                  COORDINATES (6371.2 KM)
C      EREQU      MAJOR HALF AXIS FOR EARTH ELLIPSOID (6378.160 KM)
C      ERPOL      MINOR HALF AXIS FOR EARTH ELLIPSOID (6356.775 KM)
C      AQUAD      SQUARE OF MAJOR HALF AXIS FOR EARTH ELLIPSOID
C      BQUAD   SQUARE OF MINOR HALF AXIS FOR EARTH ELLIPSOID
C
C ERA, EREQU and ERPOL as recommended by the INTERNATIONAL
C ASTRONOMICAL UNION .
C-----------------------------------------------------------------
        IMPLICIT NONE
      REAL*8 ERA,AQUAD,BQUAD,EREQU,ERPOL
            COMMON/GENER/ERA,AQUAD,BQUAD
      real*8 rad,pi
      common/rconst/rad,pi
      ERA=6371.2D0
c WGS84 World Geodetic System 84 (GPS)
      EREQU=6378.137D0
      ERPOL=6356.752314D0
c Older World Geodetic System
c      EREQU=6378.16D0
c      ERPOL=6356.775D0
c
      AQUAD=EREQU*EREQU
      BQUAD=ERPOL*ERPOL

       pi = 4.D0*ATAN(1.D0)
       rad = pi/180.D0

      RETURN
      END

!***************************************************************************************************
! Copyright 1990, 1991, Hedin
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
!  Contents:  MSIS FORTRAN SUBROUTINE GTD6

C-----------------------------------------------------------------------
      SUBROUTINE GTD6(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
      IMPLICIT NONE
C        Neutral Atmosphere Empirical Model from the surface to lower
C          exosphere  MSISE90 (JGR, 96, 1159-1172, 1991)
C         A.E.Hedin 4/24/90;6/3/91(add SAVE)
C         2/11/93 correct switch initialization and mks calculation
C       2/11/97 [AEH] CORRECT ERROR IN GHP6 WHEN USING METER6(.TRUE.)
C           See subroutine GHP6 to specify a pressure rather than
C           altitude.
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD or DDD (day of year from 1 to 365)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C     Note:  Ut, Local Time, and Longitude are used independently in the
C            model and are not of equal importance for every situation.  
C            For the most physically realistic calculation these three
C            variables should be consistent (STL=SEC/3600+GLONG/15).
C            F107, F107A, and AP effects are not large below 80 km 
C            and these can be set to 150., 150., and 4. respectively.
C     OUTPUT:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)                       
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C      TO GET OUTPUT IN M-3 and KG/M3:   CALL METER6(.TRUE.) 
C
C      O, H, and N set to zero below 72.5 km
C      Exospheric temperature set to average for altitudes below 120 km.
C        
C           The following is for test and special purposes:
C            TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC5(SW)
C               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C               FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C
C              To get current values of SW: CALL TRETRV5(SW)
C
      INTEGER*4 IYD,MASS,IMR,ISW,MN3,MN2,MSSL,MSS
      INTEGER*4 I,J
      REAL*8 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),D(8),T(2)
      REAL*8 TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
      REAL*8 G0,RL,DD,DB14,TR12,PTM(10),PDM(10,8)
      REAL*8 TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      REAL*8 PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),PMA(100,10)
      REAL*8 D6(8),T6(2)
      REAL*8 ZN3(5),ZN2(4),SV(25),SW(25),SWC(25)
      REAL*8 ZMIX,ALAST
      REAL*8 PAVGM(10),GSURF,RE
      REAL*8 DM04,DM16,DM28,DM32,DM40,DM01,DM14
      REAL*8 XLAT,V1,VTST,XMM,ALTT,DM28M,GLOB6S,DMC,DZ28,DMR
      REAL*8 DENSM6,TZ
c
      CHARACTER NAME(2)*4,ISDATE(3)*4,ISTIME(2)*4
      CHARACTER ISD(3)*4,IST(2)*4,NAM(2)*4
c
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO6/TN1,TN2,TN3,TGN1,TGN2,TGN3
      COMMON/LOWER6/PTM,PDM
      COMMON/PARM6/PT,PD,PS,PDL,PTL,PMA
      COMMON/DATIM6/ISD,IST,NAM
      COMMON/DATIME/ISDATE,ISTIME,NAME
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/MAVG6/PAVGM
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      SAVE
      EXTERNAL GTD6BK
      DATA MN3/5/,ZN3/32.5D0,20.D0,15.D0,10.D0,0.D0/
      DATA MN2/4/,ZN2/72.5D0,55.D0,45.D0,32.5D0/
      DATA ZMIX/62.5D0/,ALAST/99999.D0/,MSSL/-999/
      DATA SV/25*1.D0/
      IF(ISW.NE.64999) CALL TSELEC5(SV)
C      Put identification data into common/datime/
      DO 1 I=1,3
        ISDATE(I)=ISD(I)
    1 CONTINUE
      DO 2 I=1,2
        ISTIME(I)=IST(I)
        NAME(I)=NAM(I)
    2 CONTINUE
C
Ce        Test for changed input
      V1=VTST(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,1)
C       Latitude variation of gravity (none for SW(2)=0)
      XLAT=GLAT
      IF(SW(2).EQ.0.D0) XLAT=45.D0
      CALL GLATF6(XLAT,GSURF,RE)
C
      XMM=PDM(5,3)
C
C       THERMOSPHERE/UPPER MESOSPHERE [above ZN2(1)]
      ALTT=DMAX1(ALT,ZN2(1))
      MSS=MASS
Ce       Only calculate N2 in thermosphere if alt in mixed region
      IF(ALT.LT.ZMIX.AND.MASS.GT.0) MSS=28
Ce       Only calculate thermosphere if input parameters changed
Ce         or altitude above ZN2(1) in mesosphere
      IF(V1.EQ.1.D0.OR.ALT.GT.ZN2(1).OR.ALAST.GT.ZN2(1).OR.MSS.NE.MSSL)
     $ THEN
        CALL GTS6(IYD,SEC,ALTT,GLAT,GLONG,STL,F107A,F107,AP,MSS,D6,T6)
        DM28M=DM28
C         metric adjustment
        IF(IMR.EQ.1) DM28M=DM28*1.D6
        MSSL=MSS
      ENDIF
      T(1)=T6(1)
      T(2)=T6(2)
      IF(ALT.GE.ZN2(1)) THEN
        DO 5 J=1,8
          D(J)=D6(J)
    5   CONTINUE
        GOTO 10
      ENDIF
C
C       LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3(1) and ZN2(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
Ce        Only calculate nodes if input changed
       IF(V1.EQ.1.D0.OR.ALAST.GE.ZN2(1)) THEN
        TGN2(1)=TGN1(2)
        TN2(1)=TN1(5)
        TN2(2)=PMA(1,1)*PAVGM(1)/(1.D0-SW(20)*GLOB6S(PMA(1,1)))
        TN2(3)=PMA(1,2)*PAVGM(2)/(1.D0-SW(20)*GLOB6S(PMA(1,2)))
        TN2(4)=PMA(1,3)*PAVGM(3)/(1.D0-SW(20)*SW(22)*GLOB6S(PMA(1,3)))
        TGN2(2)=PAVGM(9)*PMA(1,10)*(1.D0+SW(20)*SW(22)*
     &   GLOB6S(PMA(1,10)))*TN2(4)*TN2(4)/(PMA(1,3)*PAVGM(3))**2
        TN3(1)=TN2(4)
       ENDIF
       IF(ALT.GE.ZN3(1)) GOTO 6
C
C       LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
Ce        Only calculate nodes if input changed
        IF(V1.EQ.1.D0.OR.ALAST.GE.ZN3(1)) THEN
         TGN3(1)=TGN2(2)
         TN3(2)=PMA(1,4)*PAVGM(4)/(1.D0-SW(22)*GLOB6S(PMA(1,4)))
         TN3(3)=PMA(1,5)*PAVGM(5)/(1.D0-SW(22)*GLOB6S(PMA(1,5)))
         TN3(4)=PMA(1,6)*PAVGM(6)/(1.D0-SW(22)*GLOB6S(PMA(1,6)))
         TN3(5)=PMA(1,7)*PAVGM(7)/(1.D0-SW(22)*GLOB6S(PMA(1,7)))
         TGN3(2)=PMA(1,8)*PAVGM(8)*(1.D0+SW(22)*GLOB6S(PMA(1,8)))
     $   *TN3(5)*TN3(5)/(PMA(1,7)*PAVGM(7))**2
        ENDIF
    6   CONTINUE
        IF(MASS.EQ.0) GOTO 50
Ce          Linear transition to full mixing at ZMIX from almost
Ce            full mixing at ZN2(1) to improve efficiency
        DMC=0.D0
        IF(ALT.GT.ZMIX) DMC=1.D0-(ZN2(1)-ALT)/(ZN2(1)-ZMIX)
        DZ28=D6(3)
C      ***** N2 DENSITY ****
        DMR=D6(3)/DM28M-1.D0
        D(3)=DENSM6(ALT,DM28M,XMM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
        D(3)=D(3)*(1.D0+DMR*DMC)
C      ***** HE DENSITY ****
        D(1)=0.D0
        IF(MASS.NE.4.AND.MASS.NE.48) GOTO 204
          DMR=D6(1)/(DZ28*PDM(2,1))-1.D0
          D(1)=D(3)*PDM(2,1)*(1.D0+DMR*DMC)
  204   CONTINUE
C      **** O DENSITY ****
        D(2)=0.D0
  216   CONTINUE
C      ***** O2 DENSITY ****
        D(4)=0.D0
        IF(MASS.NE.32.AND.MASS.NE.48) GOTO 232
          DMR=D6(4)/(DZ28*PDM(2,4))-1.D0
          D(4)=D(3)*PDM(2,4)*(1.D0+DMR*DMC)
  232   CONTINUE
C      ***** AR DENSITY ****
        D(5)=0.D0
        IF(MASS.NE.40.AND.MASS.NE.48) GOTO 240
          DMR=D6(5)/(DZ28*PDM(2,5))-1.D0
          D(5)=D(3)*PDM(2,5)*(1.D0+DMR*DMC)
  240   CONTINUE
C      ***** HYDROGEN DENSITY ****
        D(7)=0.D0
C      ***** ATOMIC NITROGEN DENSITY ****
        D(8)=0.D0
C
C       TOTAL MASS DENSITY
C
        IF(MASS.EQ.48) THEN
         D(6) = 1.66D-24*(4.D0*D(1)+16.D0*D(2)+28.D0*D(3)+32.D0*D(4)+
     &       40.D0*D(5)+D(7)+14.D0*D(8))  
         IF(IMR.EQ.1) D(6)=D(6)/1000.D0
        ENDIF
        T(2)=TZ
   10 CONTINUE
      GOTO 90
   50 CONTINUE
      DD=DENSM6(ALT,1.D0,0.D0,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)                
      T(2)=TZ
   90 CONTINUE
      ALAST=ALT
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GHP6(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,
     $  D,T,PRESS)
      IMPLICIT NONE
C       FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD6
C       2/11/97 [AEH] CORRECT ERROR IN GHP6 WHEN USING METER6(.TRUE.)
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD
C        SEC - UT(SEC)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        PRESS - PRESSURE LEVEL(MB)
C     OUTPUT:
C        ALT - ALTITUDE(KM) 
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      INTEGER*4 IYD,IMR,IDAY,L
      REAL*8 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),D(8),T(2),PRESS
      REAL*8 GSURF,RE,BM,RGAS,TEST,PL,ZI,CL,CL2,CD,CA,Z,XN,P,DIFF,XM
      REAL*8 G,SH
c
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      SAVE
      DATA BM/1.3806D-19/,RGAS/831.4D0/
      DATA TEST/.00043D0/
      PL=DLOG10(PRESS)
C      Initial altitude estimate
      IF(PL.GE.-5.D0) THEN
         IF(PL.GT.2.5D0) ZI=18.06D0*(3.00D0-PL)
         IF(PL.GT..75D0.AND.PL.LE.2.5D0) ZI=14.98D0*(3.08D0-PL)
         IF(PL.GT.-1.D0.AND.PL.LE..75D0) ZI=17.8D0*(2.72D0-PL)
         IF(PL.GT.-2.D0.AND.PL.LE.-1.D0) ZI=14.28D0*(3.64D0-PL)
         IF(PL.GT.-4.D0.AND.PL.LE.-2.D0) ZI=12.72D0*(4.32D0-PL)
         IF(PL.LE.-4.D0) ZI=25.3D0*(.11D0-PL)
         IDAY=MOD(IYD,1000)
         CL=GLAT/90.D0
         CL2=CL*CL
         IF(IDAY.LT.182) CD=1.D0-IDAY/91.25D0
         IF(IDAY.GE.182) CD=IDAY/91.25D0-3.D0
         CA=0.D0
         IF(PL.GT.-1.11D0.AND.PL.LE.-.23D0) CA=1.0D0
         IF(PL.GT.-.23D0) CA=(2.79D0-PL)/(2.79D0+.23D0)
         IF(PL.LE.-1.11D0.AND.PL.GT.-3.D0) CA=(-2.93D0-PL)/
     &      (-2.93D0+1.11D0)
         Z=ZI-4.87D0*CL*CD*CA-1.64D0*CL2*CA+.31D0*CA*CL
      ENDIF
      IF(PL.LT.-5.D0) Z=22.D0*(PL+4.D0)**2+110.D0
      L=0
C      ITERATION LOOP
   10 CONTINUE
        L=L+1
        CALL GTD6(IYD,SEC,Z,GLAT,GLONG,STL,F107A,F107,AP,48,D,T)
        XN=D(1)+D(2)+D(3)+D(4)+D(5)+D(7)+D(8)
        P=BM*XN*T(2)
        IF(IMR.EQ.1) P=P*1.D-6
        DIFF=PL-DLOG10(P)
        IF(ABS(DIFF).LT.TEST .OR. L.EQ.6) GOTO 20
        XM=D(6)/XN/1.66D-24
        IF(IMR.EQ.1) XM = XM*1.D3
        G=GSURF/(1.D0+Z/RE)**2
        SH=RGAS*T(2)/(XM*G)
C         New altitude estimate using scale height
        Z=Z-SH*DIFF*2.302D0
        GOTO 10
   20 CONTINUE
      IF(L.EQ.6) WRITE(6,100) PRESS,DIFF
  100 FORMAT(1X,29HGHP6 NOT CONVERGING FOR PRESS,1PE12.2,E12.2)
      ALT=Z
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLATF6(LAT,GV,REFF)
      IMPLICIT NONE
C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
C      RADIUS (REFF)
      REAL*8 LAT,GV,REFF,DGTR,C2
      SAVE
      DATA DGTR/1.74533D-2/
      C2 = COS(2.D0*DGTR*LAT)
      GV = 980.616D0*(1.D0-.0026373D0*C2)
      REFF = 2.D0*GV/(3.085462D-6 + 2.27D-9*C2)*1.D-5
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION VTST(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,IC)
      IMPLICIT NONE
C       Test if geophysical variables or switches changed and save
C       Return 0 if unchanged and 1 if changed
      INTEGER*4 IYD,IC,IYDL(2),ISW,I
      REAL*8 SEC,GLAT,GLONG,STL,F107A,F107,AP(7)
      REAL*8 SECL(2),GLATL(2),GLL(2),STLL(2)
      REAL*8 FAL(2),FL(2),APL(7,2),SWL(25,2),SWCL(25,2)
      REAL*8 SW(25),SWC(25),VTST
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      SAVE
      DATA IYDL/2*-999/,SECL/2*-999.D0/,GLATL/2*-999.D0/,GLL/2*-999.D0/
      DATA STLL/2*-999.D0/,FAL/2*-999.D0/,FL/2*-999.D0/,APL/14*-999.D0/
      DATA SWL/50*-999.D0/,SWCL/50*-999.D0/
      VTST=0.D0
      IF(IYD.NE.IYDL(IC)) GOTO 10
      IF(SEC.NE.SECL(IC)) GOTO 10
      IF(GLAT.NE.GLATL(IC)) GOTO 10
      IF(GLONG.NE.GLL(IC)) GOTO 10
      IF(STL.NE.STLL(IC)) GOTO 10
      IF(F107A.NE.FAL(IC)) GOTO 10
      IF(F107.NE.FL(IC)) GOTO 10
      DO 5 I=1,7
        IF(AP(I).NE.APL(I,IC)) GOTO 10
    5 CONTINUE
      DO 7 I=1,25
        IF(SW(I).NE.SWL(I,IC)) GOTO 10
        IF(SWC(I).NE.SWCL(I,IC)) GOTO 10
    7 CONTINUE
      GOTO 20
   10 CONTINUE
      VTST=1.D0
      IYDL(IC)=IYD
      SECL(IC)=SEC
      GLATL(IC)=GLAT
      GLL(IC)=GLONG
      STLL(IC)=STL
      FAL(IC)=F107A
      FL(IC)=F107
      DO 15 I=1,7
        APL(I,IC)=AP(I)
   15 CONTINUE
      DO 16 I=1,25
        SWL(I,IC)=SW(I)
        SWCL(I,IC)=SWC(I)
   16 CONTINUE
   20 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTS6(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
      IMPLICIT NONE
C        Neutral Thermosphere Model above 72.5 km for MSISE-90
C         A.E.Hedin 3/9/90
C         Coefficients not changed for 120km and above, but results may differ
C        by a few percent from MSIS-86 (GTS5) with introduction of a
C        latitude dependent accel. of gravity.
C         Lower thermosphere reformulated for better continuation into
C        lower atmosphere.
C        For efficiency:
C         Exospheric temperature left at average value for alt below 120km;
C         120 km gradient left at average value for alt below 72 km;
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM) (GREATER THAN 72.5 KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C     Note:  Ut, Local Time, and Longitude are used independently in the
C            model and are not of equal importance for every situation.  
C            For the most physically realistic calculation these three
C            variables should be consistent (STL=SEC/3600+GLONG/15).
C     OUTPUT:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C           The following is for test and special purposes:
C           (1) LOWER BOUND QUANTITIES IN COMMON/GTS3C/
C           (2) TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC5(SW)
C               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C               FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C
C              To get current values of SW: CALL TRETRV5(SW)
C
      INTEGER*4 MT(10),ISW,MASS,IMR,MN1,J
      INTEGER*4 IYD,I
      LOGICAL METER
      REAL*8 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),D(8),T(2)
      REAL*8 TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
      REAL*8 G0,RL,DD,DB14,TR12,PTM(10),PDM(10,8)
      REAL*8 TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      REAL*8 PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),PMA(100,10)
      REAL*8 SW(25),SWC(25),TINFG,GB,ROUT,TT(15)
      REAL*8 ALTL(8),ZN1(5)
      REAL*8 DGTR,DR,ALAST
      REAL*8 V2,VTST,TINF,DAY,ZHF,YRD,XMM,XMD,GLOB6S,ZLB,TZ
      REAL*8 GLOBE6,DENSU6,DNET6,CCOR6,DDUM
      REAL*8 G28,ZH28,ZHM28,B28,DM28
      REAL*8 G4,ZH04,ZHM04,B04,DM04,ZC04,HC04
      REAL*8 G16,ZH16,ZHM16,B16,DM16,ZC16,HC16,HCC16,ZCC16,RC16
      REAL*8 G32,ZH32,ZHM32,B32,DM32,HC32,ZC32
      REAL*8 G40,ZH40,ZHM40,B40,DM40,HC40,ZC40
      REAL*8 G1,ZH01,ZHM01,B01,DM01,HC01,ZC01,HCC01,ZCC01,RC01
      REAL*8 G14,ZH14,ZHM14,B14,DM14,HC14,ZC14,HCC14,ZCC14,RC14
c
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO6/TN1,TN2,TN3,TGN1,TGN2,TGN3
      COMMON/LOWER6/PTM,PDM
      COMMON/PARM6/PT,PD,PS,PDL,PTL,PMA
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/TTEST/TINFG,GB,ROUT,TT
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/METSEL/IMR
      SAVE
      DATA MT/48,0,4,16,28,32,40,1,49,14/
      DATA ALTL/200.D0,400.D0,160.D0,200.D0,240.D0,450.D0,
     &          320.D0,450.D0/
      DATA MN1/5/,ZN1/120.D0,110.D0,100.D0,90.D0,72.5D0/
      DATA DGTR/1.74533D-2/,DR/1.72142D-2/,ALAST/-999.D0/
Ce        Test for changed input
      V2=VTST(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,2)
C
      YRD=IYD*1.D0
      ZA=PDL(16,2)
      ZN1(1)=ZA
      DO 2 J=1,8
        D(J)=0.D0
    2 CONTINUE
Ce       TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
      IF(ALT.GT.ZN1(1)) THEN
        IF(V2.EQ.1.D0.OR.ALAST.LE.ZN1(1)) TINF=PTM(1)*PT(1)
     $  *(1.D0+SW(16)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PT))
      ELSE
        TINF=PTM(1)*PT(1)
      ENDIF
      T(1)=TINF
Ce         GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
      IF(ALT.GT.ZN1(5)) THEN
        IF(V2.EQ.1.D0.OR.ALAST.LE.ZN1(5)) G0=PTM(4)*PS(1)
     $   *(1.D0+SW(19)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PS))
      ELSE
        G0=PTM(4)*PS(1)
      ENDIF
Ce      Calculate these temperatures only if input changed
      IF(V2.EQ.1.D0)
     $  TLB=PTM(2)*(1.D0+SW(17)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,
     $  F107A,F107,AP,PD(1,4)))*PD(1,4)
       S=G0/(TINF-TLB)
Ce       Lower thermosphere temp variations not significant for
Ce        density above 300 km
       IF(ALT.LT.300.D0) THEN
        IF(V2.EQ.1.D0.OR.ALAST.GE.300.D0) THEN
         TN1(2)=PTM(7)*PTL(1,1)/(1.D0-SW(18)*GLOB6S(PTL(1,1)))
         TN1(3)=PTM(3)*PTL(1,2)/(1.D0-SW(18)*GLOB6S(PTL(1,2)))
         TN1(4)=PTM(8)*PTL(1,3)/(1.D0-SW(18)*GLOB6S(PTL(1,3)))
         TN1(5)=PTM(5)*PTL(1,4)/(1.D0-SW(18)*SW(20)*GLOB6S(PTL(1,4)))
         TGN1(2)=PTM(9)*PMA(1,9)*(1.D0+SW(18)*SW(20)*GLOB6S(PMA(1,9)))
     $   *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
        ENDIF
       ELSE
        TN1(2)=PTM(7)*PTL(1,1)
        TN1(3)=PTM(3)*PTL(1,2)
        TN1(4)=PTM(8)*PTL(1,3)
        TN1(5)=PTM(5)*PTL(1,4)
        TGN1(2)=PTM(9)*PMA(1,9)
     $  *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
       ENDIF
C
      Z0=ZN1(4)
      T0=TN1(4)
      ZLB=PTM(6)
      TR12=1.D0
C
      IF(MASS.EQ.0) GO TO 50
C       N2 variation factor at Zlb
      G28=SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107, 
     & AP,PD(1,3))
C        Variation of turbopause height
      DAY=DMOD(YRD,1000.D0)
      ZHF=PDL(25,2)
     $    *(1.D0+SW(5)*PDL(25,1)*SIN(DGTR*GLAT)*COS(DR*(DAY-PT(14))))
C
      YRD=IYD*1.D0
      T(1)=TINF
      XMM=PDM(5,3)
C
      DO 10 J = 1,10
      IF(MASS.EQ.MT(J))   GO TO 15
   10 CONTINUE
      WRITE(6,100) MASS
      GO TO 90
      
   15 IF(ALT.GT.ALTL(6).AND.MASS.NE.28.AND.MASS.NE.48) GO TO 17
C
C       **** N2 DENSITY ****
C
C      Diffusive density at Zlb
      DB28 = PDM(1,3)*EXP(G28)*PD(1,3)
C      Diffusive density at Alt
      D(3)=DENSU6(ALT,DB28,TINF,TLB, 28.D0,0.D0,T(2),ZLB,S
     &,MN1,ZN1,TN1,TGN1)
      DD=D(3)
C      Turbopause
      ZH28=PDM(3,3)*ZHF
      ZHM28=PDM(4,3)*PDL(6,2) 
      XMD=28.D0-XMM
C      Mixed density at Zlb
      B28=DENSU6(ZH28,DB28,TINF,TLB,XMD,-1.D0,TZ,ZLB,S,
     &MN1,ZN1,TN1,TGN1)
      IF(ALT.GT.ALTL(3).OR.SW(15).EQ.0.D0) GO TO 17
C      Mixed density at Alt
      DM28=DENSU6(ALT,B28,TINF,TLB,XMM,0.D0,TZ,ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(3)=DNET6(D(3),DM28,ZHM28,XMM,28.D0)
   17 CONTINUE
      GO TO (20,50,20,25,90,35,40,45,25,48),  J
   20 CONTINUE
C
C       **** HE DENSITY ****
C
C       Density variation factor at Zlb
      G4 = SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,1))
C      Diffusive density at Zlb
      DB04 = PDM(1,1)*EXP(G4)*PD(1,1)
C      Diffusive density at Alt
      D(1)=DENSU6(ALT,DB04,TINF,TLB, 4.D0,-.4D0,T(2),ZLB,S
     &,MN1,ZN1,TN1,TGN1)
      DD=D(1)
      IF(ALT.GT.ALTL(1).OR.SW(15).EQ.0.D0) GO TO 24
C      Turbopause
      ZH04=PDM(3,1)
      ZHM04=ZHM28
C      Mixed density at Zlb
      B04=DENSU6(ZH04,DB04,TINF,TLB,4.D0-XMM,-1.4D0,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM04=DENSU6(ALT,B04,TINF,TLB,XMM,0.D0,T(2),
     &ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(1)=DNET6(D(1),DM04,ZHM04,XMM,4.D0)
C      Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,1)/B04)
      ZC04=PDM(5,1)*PDL(1,2)
      HC04=PDM(6,1)*PDL(2,2)
C      Net density corrected at Alt
      D(1)=D(1)*CCOR6(ALT,RL,HC04,ZC04)
   24 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   25 CONTINUE
C
C      **** O DENSITY ****
C
C       Density variation factor at Zlb
      G16= SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,2))
C      Diffusive density at Zlb
      DB16 =  PDM(1,2)*EXP(G16)*PD(1,2)
C       Diffusive density at Alt
      D(2)=DENSU6(ALT,DB16,TINF,TLB, 16.D0,0.D0,T(2),ZLB,S
     &,MN1,ZN1,TN1,TGN1)
      DD=D(2)
      IF(ALT.GT.ALTL(2).OR.SW(15).EQ.0.D0) GO TO 34
C  Corrected from PDM(3,1) to PDM(3,2)  12/2/85
C       Turbopause
      ZH16=PDM(3,2)
      ZHM16=ZHM28
C      Mixed density at Zlb
      B16=DENSU6(ZH16,DB16,TINF,TLB,16.D0-XMM,-1.D0,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM16=DENSU6(ALT,B16,TINF,TLB,XMM,0.D0,T(2),
     &ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(2)=DNET6(D(2),DM16,ZHM16,XMM,16.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
      HC16=PDM(6,2)*PDL(4,2)
      ZC16=PDM(5,2)*PDL(3,2)
      D(2)=D(2)*CCOR6(ALT,RL,HC16,ZC16)
C       Chemistry correction
      HCC16=PDM(8,2)*PDL(14,2)
      ZCC16=PDM(7,2)*PDL(13,2)
      RC16=PDM(4,2)*PDL(15,2)
C      Net density corrected at Alt
      D(2)=D(2)*CCOR6(ALT,RC16,HCC16,ZCC16)
   34 CONTINUE
      IF(MASS.NE.48 .AND. MASS.NE.49) GO TO 90
   35 CONTINUE
C
C       **** O2 DENSITY ****
C
C       Density variation factor at Zlb
      G32= SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,5))
C      Diffusive density at Zlb
      DB32 = PDM(1,4)*EXP(G32)*PD(1,5)
C       Diffusive density at Alt
      D(4)=DENSU6(ALT,DB32,TINF,TLB, 32.D0,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
      IF(MASS.EQ.49) THEN
         DD=DD+2.*D(4)
      ELSE
         DD=D(4)
      ENDIF
      IF(ALT.GT.ALTL(4).OR.SW(15).EQ.0.D0) GO TO 39
C       Turbopause
      ZH32=PDM(3,4)
      ZHM32=ZHM28
C      Mixed density at Zlb
      B32=DENSU6(ZH32,DB32,TINF,TLB,32.D0-XMM,-1.D0,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM32=DENSU6(ALT,B32,TINF,TLB,XMM,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(4)=DNET6(D(4),DM32,ZHM32,XMM,32.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,4)/B32)
      HC32=PDM(6,4)*PDL(8,2)
      ZC32=PDM(5,4)*PDL(7,2)
C      Net density corrected at Alt
      D(4)=D(4)*CCOR6(ALT,RL,HC32,ZC32)
   39 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   40 CONTINUE
C
C       **** AR DENSITY ****
C
C       Density variation factor at Zlb
      G40= SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,6))
C      Diffusive density at Zlb
      DB40 = PDM(1,5)*EXP(G40)*PD(1,6)
C       Diffusive density at Alt
      D(5)=DENSU6(ALT,DB40,TINF,TLB, 40.D0,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
      DD=D(5)
      IF(ALT.GT.ALTL(5).OR.SW(15).EQ.0.D0) GO TO 44
C       Turbopause
      ZH40=PDM(3,5)
      ZHM40=ZHM28
C      Mixed density at Zlb
      B40=DENSU6(ZH40,DB40,TINF,TLB,40.D0-XMM,-1.D0,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM40=DENSU6(ALT,B40,TINF,TLB,XMM,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(5)=DNET6(D(5),DM40,ZHM40,XMM,40.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,5)/B40)
      HC40=PDM(6,5)*PDL(10,2)
      ZC40=PDM(5,5)*PDL(9,2)
C      Net density corrected at Alt
      D(5)=D(5)*CCOR6(ALT,RL,HC40,ZC40)
   44 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   45 CONTINUE
C
C        **** HYDROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G1 = SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,7))
C      Diffusive density at Zlb
      DB01 = PDM(1,6)*EXP(G1)*PD(1,7)
C       Diffusive density at Alt
      D(7)=DENSU6(ALT,DB01,TINF,TLB,1.D0,-.4D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
      DD=D(7)
      IF(ALT.GT.ALTL(7).OR.SW(15).EQ.0.D0) GO TO 47
C       Turbopause
      ZH01=PDM(3,6)
      ZHM01=ZHM28
C      Mixed density at Zlb
      B01=DENSU6(ZH01,DB01,TINF,TLB,1.D0-XMM,-1.4D0,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM01=DENSU6(ALT,B01,TINF,TLB,XMM,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(7)=DNET6(D(7),DM01,ZHM01,XMM,1.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,6)*ABS(PDL(18,2))/B01)
      HC01=PDM(6,6)*PDL(12,2)
      ZC01=PDM(5,6)*PDL(11,2)
      D(7)=D(7)*CCOR6(ALT,RL,HC01,ZC01)
C       Chemistry correction
      HCC01=PDM(8,6)*PDL(20,2)
      ZCC01=PDM(7,6)*PDL(19,2)
      RC01=PDM(4,6)*PDL(21,2)
C      Net density corrected at Alt
      D(7)=D(7)*CCOR6(ALT,RC01,HCC01,ZCC01)
   47 CONTINUE
   48 CONTINUE
C
C        **** ATOMIC NITROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G14 =SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,8))
C      Diffusive density at Zlb
      DB14 = PDM(1,7)*EXP(G14)*PD(1,8)
C       Diffusive density at Alt
      D(8)=DENSU6(ALT,DB14,TINF,TLB,14.D0,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
      DD=D(8)
      IF(ALT.GT.ALTL(8).OR.SW(15).EQ.0.D0) GO TO 49
C       Turbopause
      ZH14=PDM(3,7)
      ZHM14=ZHM28
C      Mixed density at Zlb
      B14=DENSU6(ZH14,DB14,TINF,TLB,14.D0-XMM,-1.D0,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM14=DENSU6(ALT,B14,TINF,TLB,XMM,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(8)=DNET6(D(8),DM14,ZHM14,XMM,14.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,7)*ABS(PDL(3,1))/B14)
      HC14=PDM(6,7)*PDL(2,1)
      ZC14=PDM(5,7)*PDL(1,1)
      D(8)=D(8)*CCOR6(ALT,RL,HC14,ZC14)
C       Chemistry correction
      HCC14=PDM(8,7)*PDL(5,1)
      ZCC14=PDM(7,7)*PDL(4,1)
      RC14=PDM(4,7)*PDL(6,1)
C      Net density corrected at Alt
      D(8)=D(8)*CCOR6(ALT,RC14,HCC14,ZCC14)
   49 CONTINUE
      IF(MASS.NE.48) GO TO 90
C
C       TOTAL MASS DENSITY
C
      D(6) = 1.66D-24*(4.*D(1)+16.D0*D(2)+28.D0*D(3)+32.D0*D(4)+
     &       40.D0*D(5)+D(7)+14.*D(8))
      DB48=1.66D-24*(4.D0*DB04+16.D0*DB16+28.D0*DB28+32.D0*DB32+
     &       40.D0*DB40+DB01+14.D0*DB14)
      GO TO 90
C       TEMPERATURE AT ALTITUDE
   50 CONTINUE
      DDUM=DENSU6(ALT,1.D0,TINF,TLB,0.D0,0.D0,T(2),ZLB,S,
     &MN1,ZN1,TN1,TGN1)
      GO TO 90
   90 CONTINUE
C       ADJUST DENSITIES FROM CGS TO KGM
      IF(IMR.EQ.1) THEN
        DO 95 I=1,8
          D(I)=D(I)*1.D6
   95   CONTINUE
        D(6)=D(6)/1000.D0
      ENDIF
      ALAST=ALT
      RETURN
  100 FORMAT(1X,'MASS', I5, '  NOT VALID')
      ENTRY METER6(METER)
      IMR=0
      IF(METER) IMR=1
      END
C-----------------------------------------------------------------------
      FUNCTION GLOBE6(YRD,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)
      IMPLICIT NONE
C       CALCULATE G(L) FUNCTION 
C       Upper Thermosphere Parameters
      INTEGER*4 IYR,ISW,NSW,J,I
c
      REAL*8 LAT, LONG,LONGL
      REAL*8 P(150),SV(25),AP(7),TINF,GB,ROUT,T(15)
      REAL*8 PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC
      REAL*8 DAY,DF,DFA,APD,APDF,APT(4),SW(25),SWC(25)
      REAL*8 YRD,SEC,TLOC,F107A,F107,XLONG,CLONG,SLONG
      REAL*8 DGTR,DR,XL,TLL,SW9,DAYL,P14,P18,P32,HR,SR,P39
      REAL*8 A,G0,SUMEX,SG0,EX
      REAL*8 C,S,C2,C4,S2
      REAL*8 CD14,CD18,CD32,CD39,F1,F2,T71,T72,T81,T82
      REAL*8 P44,P45,EXP1,EXP2,GLOBE6
c
      COMMON/TTEST/TINF,GB,ROUT,T
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/LPOLY/PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ DAY,DF,DFA,APD,APDF,APT,XLONG,CLONG,SLONG
      COMMON/LPOLYI/IYR
      SAVE
      DATA DGTR/1.74533D-2/,DR/1.72142D-2/, XL/1000.D0/,TLL/1000.D0/
      DATA SW9/1.D0/,DAYL/-1.D0/,P14/-1000.D0/,P18/-1000.D0/,
     & P32/-1000.D0/
      DATA HR/.2618D0/,SR/7.2722D-5/,SV/25*1.D0/,NSW/14/,P39/-1000.D0/
      DATA LONGL/-999.D0/
C       3hr Magnetica activity functions
      G0(A)=(A-4.D0+(P(26)-1.D0)*(A-4.D0+(EXP(-ABS(P(25))*
     &(A-4.D0))-1.D0)/ABS(P(25))))
       SUMEX(EX)=1.D0+(1.D0-EX**19)/(1.D0-EX)*EX**(.5D0)
      SG0(EX)=(G0(AP(2))+(G0(AP(3))*EX+G0(AP(4))*EX*EX+G0(AP(5))*EX**3
     $ +(G0(AP(6))*EX**4+G0(AP(7))*EX**12)*(1.D0-EX**8)/(1.D0-EX))
     $ )/SUMEX(EX)
C
      IF(ISW.NE.64999) CALL TSELEC5(SV)
      DO 10 J=1,14
       T(J)=0.D0
   10 CONTINUE
      IF(SW(9).GT.0.D0) SW9=1.D0
      IF(SW(9).LT.0.D0) SW9=-1.D0
      IYR = INT(YRD/1000.D0)
      DAY = YRD - IYR*1000.D0
      XLONG=LONG
      IF(XL.EQ.LAT)   GO TO 15
C          CALCULATE LEGENDRE POLYNOMIALS
      C = SIN(LAT*DGTR)
      S = COS(LAT*DGTR)
      C2 = C*C
      C4 = C2*C2
      S2 = S*S
      PLG(2,1) = C
      PLG(3,1) = 0.5D0*(3.D0*C2 -1.D0)
      PLG(4,1) = 0.5D0*(5.D0*C*C2-3.D0*C)
      PLG(5,1) = (35.D0*C4 - 30.D0*C2 + 3.D0)/8.D0
      PLG(6,1) = (63.D0*C2*C2*C - 70.D0*C2*C + 15.D0*C)/8.D0
      PLG(7,1) = (11.D0*C*PLG(6,1) - 5.D0*PLG(5,1))/6.D0
C     PLG(8,1) = (13.D0*C*PLG(7,1) - 6.D0*PLG(6,1))/7.D0
      PLG(2,2) = S
      PLG(3,2) = 3.D0*C*S
      PLG(4,2) = 1.5D0*(5.D0*C2-1.D0)*S
      PLG(5,2) = 2.5D0*(7.D0*C2*C-3.D0*C)*S
      PLG(6,2) = 1.875D0*(21.D0*C4 - 14.D0*C2 +1.D0)*S
      PLG(7,2) = (11.D0*C*PLG(6,2)-6.D0*PLG(5,2))/5.D0
C     PLG(8,2) = (13.D0*C*PLG(7,2)-7.D0*PLG(6,2))/6.D0
C     PLG(9,2) = (15.D0*C*PLG(8,2)-8.D0*PLG(7,2))/7.D0
      PLG(3,3) = 3.D0*S2
      PLG(4,3) = 15.D0*S2*C
      PLG(5,3) = 7.5D0*(7.D0*C2 -1.D0)*S2
      PLG(6,3) = 3.D0*C*PLG(5,3)-2.D0*PLG(4,3)
      PLG(7,3)=(11.D0*C*PLG(6,3)-7.D0*PLG(5,3))/4.D0
      PLG(8,3)=(13.D0*C*PLG(7,3)-8.D0*PLG(6,3))/5.D0
      PLG(4,4) = 15.D0*S2*S
      PLG(5,4) = 105.D0*S2*S*C 
      PLG(6,4)=(9.D0*C*PLG(5,4)-7.D0*PLG(4,4))/2.D0
      PLG(7,4)=(11.D0*C*PLG(6,4)-8.D0*PLG(5,4))/3.D0
      XL=LAT
   15 CONTINUE
      IF(TLL.EQ.TLOC)   GO TO 16
      IF(SW(7).EQ.0.D0.AND.SW(8).EQ.0.D0.AND.SW(14).EQ.0.D0) GOTO 16
      STLOC = SIN(HR*TLOC)
      CTLOC = COS(HR*TLOC)
      S2TLOC = SIN(2.D0*HR*TLOC)
      C2TLOC = COS(2.D0*HR*TLOC)
      S3TLOC = SIN(3.D0*HR*TLOC)
      C3TLOC = COS(3.D0*HR*TLOC)
      TLL = TLOC
   16 CONTINUE
      IF(LONG.NE.LONGL) THEN
        CLONG=COS(DGTR*LONG)
        SLONG=SIN(DGTR*LONG)
      ENDIF
      LONGL=LONG
      IF(DAY.NE.DAYL.OR.P(14).NE.P14) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P(18).NE.P18) CD18=COS(2.D0*DR*(DAY-P(18)))
      IF(DAY.NE.DAYL.OR.P(32).NE.P32) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P(39).NE.P39) CD39=COS(2.D0*DR*(DAY-P(39)))
      DAYL = DAY
      P14 = P(14)
      P18 = P(18)
      P32 = P(32)
      P39 = P(39)
C         F10.7 EFFECT
      DF = F107 - F107A
      DFA=F107A-150.D0
      T(1) =  P(20)*DF + P(21)*DF*DF + P(22)*DFA
     & + P(30)*DFA**2
      F1 = 1.D0 + (P(48)*DFA +P(20)*DF+P(21)*DF*DF)*SWC(1)
      F2 = 1.D0 + (P(50)*DFA+P(20)*DF+P(21)*DF*DF)*SWC(1)
C        TIME INDEPENDENT
      T(2) =
     1  (P(2)*PLG(3,1) + P(3)*PLG(5,1)+P(23)*PLG(7,1))
     & +(P(15)*PLG(3,1))*DFA*SWC(1)
     2 +P(27)*PLG(2,1)
C        SYMMETRICAL ANNUAL
      T(3) =
     1 (P(19) )*CD32
C        SYMMETRICAL SEMIANNUAL
      T(4) =
     1 (P(16)+P(17)*PLG(3,1))*CD18
C        ASYMMETRICAL ANNUAL
      T(5) =  F1*
     1  (P(10)*PLG(2,1)+P(11)*PLG(4,1))*CD14
C         ASYMMETRICAL SEMIANNUAL
      T(6) =    P(38)*PLG(2,1)*CD39
C        DIURNAL
      IF(SW(7).EQ.0.D0) GOTO 200
      T71 = (P(12)*PLG(3,2))*CD14*SWC(5)
      T72 = (P(13)*PLG(3,2))*CD14*SWC(5)
      T(7) = F2*
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2) + P(28)*PLG(6,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2) +P(29)*PLG(6,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0.D0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
      T(8) = F2*
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0.D0) GOTO 220
      T(14) = F2*
     1 ((P(40)*PLG(4,4)+(P(94)*PLG(5,4)+P(47)*PLG(7,4))*CD14*SWC(5))*
     $ S3TLOC
     2 +(P(41)*PLG(4,4)+(P(95)*PLG(5,4)+P(49)*PLG(7,4))*CD14*SWC(5))*
     $ C3TLOC)
  220 CONTINUE
C          MAGNETIC ACTIVITY BASED ON DAILY AP

      IF(SW9.EQ.-1.D0) GO TO 30
      APD=(AP(1)-4.D0)
      P44=P(44)
      P45=P(45)
      IF(P44.LT.0.D0) P44=1.D-5
      APDF = (APD+(P45-1.D0)*(APD+(EXP(-P44  *APD)-1.D0)/P44  ))
      IF(SW(9).EQ.0.D0) GOTO 40
      T(9)=APDF*(P(33)+P(46)*PLG(3,1)+P(35)*PLG(5,1)+
     $ (P(101)*PLG(2,1)+P(102)*PLG(4,1)+P(103)*PLG(6,1))*CD14*SWC(5)+
     $ (P(122)*PLG(2,2)+P(123)*PLG(4,2)+P(124)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(125))))
      GO TO 40
   30 CONTINUE
      IF(P(52).EQ.0.D0) GO TO 40
      EXP1 = EXP(-10800.D0*ABS(P(52))/(1.D0+P(139)*(45.D0-ABS(LAT))))
      IF(EXP1.GT..99999D0) EXP1=.99999D0
      EXP2 = EXP(-10800.D0*ABS(P(54)))
      IF(EXP2.GT..99999D0) EXP2=.99999D0
      IF(P(25).LT.1.D-4) P(25)=1.D-4
      APT(1)=SG0(EXP1)
      APT(3)=SG0(EXP2)
      IF(SW(9).EQ.0.D0) GOTO 40
      T(9) = APT(1)*(P(51)+P(97)*PLG(3,1)+P(55)*PLG(5,1)+
     $ (P(126)*PLG(2,1)+P(127)*PLG(4,1)+P(128)*PLG(6,1))*CD14*SWC(5)+
     $ (P(129)*PLG(2,2)+P(130)*PLG(4,2)+P(131)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(132))))
  40  CONTINUE
      IF(SW(10).EQ.0.D0.OR.LONG.LE.-1000.D0) GO TO 49
C        LONGITUDINAL
      IF(SW(11).EQ.0.D0) GOTO 230
      T(11)= (1.D0+P(81)*DFA*SWC(1))*
     $((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $ +P(104)*PLG(2,2)+P(105)*PLG(4,2)+P(106)*PLG(6,2)
     $ +SWC(5)*(P(110)*PLG(2,2)+P(111)*PLG(4,2)+P(112)*PLG(6,2))*CD14)*
     $     CLONG
     $ +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $ +P(107)*PLG(2,2)+P(108)*PLG(4,2)+P(109)*PLG(6,2)
     $ +SWC(5)*(P(113)*PLG(2,2)+P(114)*PLG(4,2)+P(115)*PLG(6,2))*CD14)*
     $  SLONG)
  230 CONTINUE
C        UT AND MIXED UT,LONGITUDE
      IF(SW(12).EQ.0.D0) GOTO 240
      T(12)=(1.D0+P(96)*PLG(2,1))*(1.D0+P(82)*DFA*SWC(1))*
     $(1.D0+P(120)*PLG(2,1)*SWC(5)*CD14)*
     $((P(69)*PLG(2,1)+P(70)*PLG(4,1)+P(71)*PLG(6,1))*
     $     COS(SR*(SEC-P(72))))
      T(12)=T(12)+SWC(11)*
     $ (P(77)*PLG(4,3)+P(78)*PLG(6,3)+P(79)*PLG(8,3))*
     $     COS(SR*(SEC-P(80))+2.D0*DGTR*LONG)*(1.D0+P(138)*DFA*SWC(1))
  240 CONTINUE
C        UT,LONGITUDE MAGNETIC ACTIVITY
      IF(SW(13).EQ.0.D0) GOTO 48
      IF(SW9.EQ.-1.D0) GO TO 45
      T(13)= APDF*SWC(11)*(1.D0+P(121)*PLG(2,1))*
     $((P( 61)*PLG(3,2)+P( 62)*PLG(5,2)+P( 63)*PLG(7,2))*
     $     COS(DGTR*(LONG-P( 64))))
     $ +APDF*SWC(11)*SWC(5)*
     $ (P(116)*PLG(2,2)+P(117)*PLG(4,2)+P(118)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(119)))
     $ + APDF*SWC(12)*
     $ (P( 84)*PLG(2,1)+P( 85)*PLG(4,1)+P( 86)*PLG(6,1))*
     $     COS(SR*(SEC-P( 76)))
      GOTO 48
   45 CONTINUE
      IF(P(52).EQ.0.D0) GOTO 48
      T(13)=APT(1)*SWC(11)*(1.D0+P(133)*PLG(2,1))*
     $((P(53)*PLG(3,2)+P(99)*PLG(5,2)+P(68)*PLG(7,2))*
     $     COS(DGTR*(LONG-P(98))))
     $ +APT(1)*SWC(11)*SWC(5)*
     $ (P(134)*PLG(2,2)+P(135)*PLG(4,2)+P(136)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(137)))
     $ +APT(1)*SWC(12)*
     $ (P(56)*PLG(2,1)+P(57)*PLG(4,1)+P(58)*PLG(6,1))*
     $     COS(SR*(SEC-P(59)))
   48 CONTINUE
   49 CONTINUE
      TINF=P(31)
      DO 50 I = 1,NSW
   50 TINF = TINF + ABS(SW(I))*T(I)
      GLOBE6 = TINF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TSELEC5(SV)
      IMPLICIT NONE
C        SET SWITCHES
C        SW FOR MAIN TERMS, SWC FOR CROSS TERMS
      INTEGER*4 I,ISW
C
      REAL*8 SW(25),SWC(25),SV(25),SVV(25),SAV(25)
c
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      SAVE
      DO 100 I = 1,25
        SAV(I)=SV(I)
        SW(I)=DMOD(SV(I),2.D0)
        IF(ABS(SV(I)).EQ.1.D0.OR.ABS(SV(I)).EQ.2.D0) THEN
          SWC(I)=1.D0
        ELSE
          SWC(I)=0.D0
        ENDIF
  100 CONTINUE
      ISW=64999
      RETURN
      ENTRY TRETRV5(SVV)
      DO 200 I=1,25
        SVV(I)=SAV(I)
  200 CONTINUE
      END
C-----------------------------------------------------------------------
      FUNCTION GLOB6S(P)
      IMPLICIT NONE
C      VERSION OF GLOBE FOR LOWER ATMOSPHERE 1/17/90
      INTEGER*4 ISW,J,I,IYR
      REAL*8 LONG,DR,DGTR,DAYL,P32,P18,P14,P39
      REAL*8 SW(25),SWC(25),P(100),T(14)
      REAL*8 PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC
      REAL*8 DAY,DF,DFA,APD,APDF,APT(4),CLONG,SLONG
      REAL*8 CD32,CD18,CD14,CD39,T71,T72,T81,T82,TT,GLOB6S
c
      COMMON/LPOLY/PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ DAY,DF,DFA,APD,APDF,APT,LONG,CLONG,SLONG
      COMMON/LPOLYI/IYR
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      SAVE
      DATA DR/1.72142D-2/,DGTR/1.74533D-2/
      DATA DAYL/-1.D0/,P32,P18,P14,P39/4*-1000.D0/
      DO 10 J=1,14
        T(J)=0.D0
   10 CONTINUE
      IF(DAY.NE.DAYL.OR.P32.NE.P(32)) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P18.NE.P(18)) CD18=COS(2.D0*DR*(DAY-P(18)))       
      IF(DAY.NE.DAYL.OR.P14.NE.P(14)) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P39.NE.P(39)) CD39=COS(2.D0*DR*(DAY-P(39)))
      DAYL=DAY
      P32=P(32)
      P18=P(18)
      P14=P(14)
      P39=P(39)
C
C       F10.7
      T(1)=P(22)*DFA
C       TIME INDEPENDENT
      T(2)=P(2)*PLG(3,1)+P(3)*PLG(5,1)+P(23)*PLG(7,1)
     $     +P(27)*PLG(2,1)+P(28)*PLG(4,1)+P(29)*PLG(6,1)
C       SYMMETRICAL ANNUAL
      T(3)=(P(19)+P(48)*PLG(3,1)+P(30)*PLG(5,1))*CD32
C       SYMMETRICAL SEMIANNUAL
      T(4)=(P(16)+P(17)*PLG(3,1)+P(31)*PLG(5,1))*CD18
C       ASYMMETRICAL ANNUAL
      T(5)=(P(10)*PLG(2,1)+P(11)*PLG(4,1)+P(36)*PLG(6,1))*CD14
C       ASYMMETRICAL SEMIANNUAL
      T(6)=(P(38)*PLG(2,1))*CD39
C        DIURNAL
      IF(SW(7).EQ.0.D0) GOTO 200
      T71 = P(12)*PLG(3,2)*CD14*SWC(5)
      T72 = P(13)*PLG(3,2)*CD14*SWC(5)
      T(7) = 
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0.D0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(47)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(49)*PLG(6,3))*CD14*SWC(5)
      T(8) = 
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0.D0) GOTO 220
      T(14) = P(40)*PLG(4,4)*S3TLOC
     $ +P(41)*PLG(4,4)*C3TLOC
  220 CONTINUE
C       MAGNETIC ACTIVITY
      IF(SW(9).EQ.0.D0) GOTO 40
      IF(SW(9).EQ.1.D0)
     $ T(9)=APDF*(P(33)+P(46)*PLG(3,1)*SWC(2))
      IF(SW(9).EQ.-1.D0)
     $ T(9)=(P(51)*APT(3)+P(97)*PLG(3,1)*APT(3)*SWC(2))
   40 CONTINUE
      IF(SW(10).EQ.0.D0.OR.SW(11).EQ.0.D0.OR.LONG.LE.-1000.D0) GO TO 49
C        LONGITUDINAL
      T(11)= (1.D0+PLG(2,1)*(P(81)*SWC(5)*COS(DR*(DAY-P(82)))
     $           +P(86)*SWC(6)*COS(2.D0*DR*(DAY-P(87))))
     $        +P(84)*SWC(3)*COS(DR*(DAY-P(85)))
     $           +P(88)*SWC(4)*COS(2.D0*DR*(DAY-P(89))))
     $ *((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $   +P(75)*PLG(2,2)+P(76)*PLG(4,2)+P(77)*PLG(6,2)
     $    )*CLONG
     $  +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $   +P(78)*PLG(2,2)+P(79)*PLG(4,2)+P(80)*PLG(6,2)
     $    )*SLONG)
   49 CONTINUE
      TT=0.D0
      DO 50 I=1,14
   50 TT=TT+ABS(SW(I))*T(I)
      GLOB6S=TT
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSU6(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,
     $  MN1,ZN1,TN1,TGN1)
      IMPLICIT NONE
C       Calculate Temperature and Density Profiles for MSIS models
C       New lower thermo polynomial 10/30/89
      INTEGER*4 MN1,MP,II,JG,LT,IERR,IFUN,N,J
      INTEGER*4 MN,K
c
      REAL*8 ZN1(MN1),TN1(MN1),TGN1(2),XS(5),YS(5),Y2OUT(5)
      REAL*8 ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,GSURF,RE,RGAS
      REAL*8 QPB(50),DV(60),ZETA,DENSU6,ZA,Z,ZG2,TT,TA,DTA
      REAL*8 Z1,Z2,T1,T2,ZG,ZGDIF,YD1,YD2,X,Y,GLB,GAMM,YI,EXPL
      REAL*8 GAMMA,DENSA,ZZ,ZL
c
      COMMON/PARMB/GSURF,RE
      COMMON/LSQV/MP,II,JG,LT,QPB,IERR,IFUN,N,J,DV
      SAVE
      DATA RGAS/831.4D0/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
CCCCCCWRITE(6,*) 'DB',ALT,DLB,TINF,TLB,XM,ALPHA,ZLB,S2,MN1,ZN1,TN1
      DENSU6=1.D0
C        Joining altitude of Bates and spline
      ZA=ZN1(1)
      Z=DMAX1(ALT,ZA)
C      Geopotential altitude difference from ZLB
      ZG2=ZETA(Z,ZLB)
C      Bates temperature
      TT=TINF-(TINF-TLB)*EXP(-S2*ZG2)
      TA=TT
      TZ=TT
      DENSU6=TZ
      IF(ALT.GE.ZA) GO TO 10
C
C       CALCULATE TEMPERATURE BELOW ZA
C      Temperature gradient at ZA from Bates profile
      DTA=(TINF-TA)*S2*((RE+ZLB)/(RE+ZA))**2
      TGN1(1)=DTA 
      TN1(1)=TA
      Z=DMAX1(ALT,ZN1(MN1))
      MN=MN1
      Z1=ZN1(1)
      Z2=ZN1(MN)
      T1=TN1(1)
      T2=TN1(MN)
C      Geopotental difference from Z1
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 20 K=1,MN
        XS(K)=ZETA(ZN1(K),Z1)/ZGDIF
        YS(K)=1.D0/TN1(K)
   20 CONTINUE
C        End node derivatives
      YD1=-TGN1(1)/(T1*T1)*ZGDIF
      YD2=-TGN1(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE6(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT6(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1.D0/Y
      DENSU6=TZ
   10 IF(XM.EQ.0.D0) GO TO 50
C
C      CALCULATE DENSITY ABOVE ZA
      GLB=GSURF/(1.D0+ZLB/RE)**2
      GAMMA=XM*GLB/(S2*RGAS*TINF)
      EXPL=EXP(-S2*GAMMA*ZG2)
      IF(EXPL.GT.50.D0.OR.TT.LE.0.D0) THEN
        EXPL=50.D0
      ENDIF
C       Density at altitude
      DENSA=DLB*(TLB/TT)**(1.D0+ALPHA+GAMMA)*EXPL
      DENSU6=DENSA
      IF(ALT.GE.ZA) GO TO 50
C
C      CALCULATE DENSITY BELOW ZA
      GLB=GSURF/(1.D0+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       integrate spline temperatures
      CALL SPLINI6(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.D0.OR.TZ.LE.0.D0) THEN
        EXPL=50.D0
      ENDIF
C       Density at altitude
      DENSU6=DENSU6*(T1/TZ)**(1.D0+ALPHA)*EXP(-EXPL)
   50 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSM6(ALT,D0,XM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
      IMPLICIT NONE
C       Calculate Temperature and Density Profiles for lower atmos.
      INTEGER*4 MN3,MN2,MP,II,JG,LT,IERR,IFUN,N,J
      INTEGER*4 MN,K
c
      REAL*8 ZN3(MN3),TN3(MN3),TGN3(2),XS(10),YS(10),Y2OUT(10)
      REAL*8 ZN2(MN2),TN2(MN2),TGN2(2)
      REAL*8 ALT,D0,XM,TZ,RGAS,GSURF,RE,TAF,QPB(50),DV(60)
      REAL*8 ZETA,DENSM6,Z,Z1,Z2,T1,Y2,ZG,ZGDIF,YD1,YD2,X,Y
      REAL*8 GLB,GAMM,YI,EXPL,T2,ZZ,ZL
c
      COMMON/PARMB/GSURF,RE
      COMMON/FIT/TAF
      COMMON/LSQV/MP,II,JG,LT,QPB,IERR,IFUN,N,J,DV
      SAVE
      DATA RGAS/831.4D0/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
      DENSM6=D0
      IF(ALT.GT.ZN2(1)) GOTO 50
C      STRATOSPHERE/MESOSPHERE TEMPERATURE
      Z=DMAX1(ALT,ZN2(MN2))
      MN=MN2
      Z1=ZN2(1)
      Z2=ZN2(MN)
      T1=TN2(1)
      T2=TN2(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 210 K=1,MN
        XS(K)=ZETA(ZN2(K),Z1)/ZGDIF
        YS(K)=1.D0/TN2(K)
  210 CONTINUE
      YD1=-TGN2(1)/(T1*T1)*ZGDIF
      YD2=-TGN2(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE6(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT6(XS,YS,Y2OUT,MN,X,Y)
C       Temperature at altitude
      TZ=1.D0/Y
      IF(XM.EQ.0.D0) GO TO 20
C
C      CALCULATE STRATOSPHERE/MESOSPHERE DENSITY 
      GLB=GSURF/(1.D0+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       Integrate temperature profile
      CALL SPLINI6(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.D0) EXPL=50.D0
C       Density at altitude
      DENSM6=DENSM6*(T1/TZ)*EXP(-EXPL)
   20 CONTINUE
      IF(ALT.GT.ZN3(1)) GOTO 50
C
C      TROPOSPHERE/STRATOSPHERE TEMPERATURE
      Z=ALT
      MN=MN3
      Z1=ZN3(1)
      Z2=ZN3(MN)
      T1=TN3(1)
      T2=TN3(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 220 K=1,MN
        XS(K)=ZETA(ZN3(K),Z1)/ZGDIF
        YS(K)=1.D0/TN3(K)
  220 CONTINUE
      YD1=-TGN3(1)/(T1*T1)*ZGDIF
      YD2=-TGN3(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE6(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT6(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1.D0/Y
      IF(XM.EQ.0.D0) GO TO 30
C
C      CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY 
C     
      GLB=GSURF/(1.D0+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C        Integrate temperature profile
      CALL SPLINI6(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.D0) EXPL=50.D0
C        Density at altitude
      DENSM6=DENSM6*(T1/TZ)*EXP(-EXPL)
   30 CONTINUE
   50 CONTINUE
      IF(XM.EQ.0.D0) DENSM6=TZ
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINE6(X,Y,N,YP1,YPN,Y2)
      IMPLICIT NONE
C        CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
C        ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
C        X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        N: SIZE OF ARRAYS X,Y
C        YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
C                 >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
C        Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
      INTEGER*4 NMAX
      PARAMETER (NMAX=100)
      INTEGER*4 N,I,K
      REAL*8 X(N),Y(N),Y2(N),U(NMAX),YP1,YPN,SIG,P,QN,UN
      SAVE
      IF(YP1.GT..99D30) THEN
        Y2(1)=0.D0
        U(1)=0.D0
      ELSE
        Y2(1)=-.5D0
        U(1)=(3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.D0
        Y2(I)=(SIG-1.D0)/P
        U(I)=(6.D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     $    /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
   11 CONTINUE
      IF(YPN.GT..99D30) THEN
        QN=0.D0
        UN=0.D0
      ELSE
        QN=.5D0
        UN=(3.D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.D0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
   12 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINT6(XA,YA,Y2A,N,X,Y)
      IMPLICIT NONE
C        CALCULATE CUBIC SPLINE INTERP VALUE
C        ADAPTED FROM NUMBERICAL RECIPES BY PRESS ET AL.
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA FOR INTERPOLATION
C        Y: OUTPUT VALUE
      INTEGER*4 N,KLO,KHI,K
c
      REAL*8 XA(N),YA(N),Y2A(N),X,Y
      REAL*8 H,A,B
      SAVE
      KLO=1
      KHI=N
    1 CONTINUE
      IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF(H.EQ.0.D0) WRITE(6,*) 'BAD XA INPUT TO SPLINT6'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     $  ((A*A*A-A)*Y2A(KLO)+(B*B*B-B)*Y2A(KHI))*H*H/6.D0
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINI6(XA,YA,Y2A,N,X,YI)
      IMPLICIT NONE
C       INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA ENDPOINT FOR INTEGRATION
C        Y: OUTPUT VALUE
      INTEGER*4 N,KLO,KHI
c
      REAL*8 XA(N),YA(N),Y2A(N),X,YI
      REAL*8 H,A,B,A2,B2,XX
      SAVE
      YI=0.D0
      KLO=1
      KHI=2
    1 CONTINUE
      IF(X.GT.XA(KLO).AND.KHI.LE.N) THEN
        XX=X
        IF(KHI.LT.N) XX=DMIN1(X,XA(KHI))
        H=XA(KHI)-XA(KLO)
        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        A2=A*A
        B2=B*B
        YI=YI+((1.D0-A2)*YA(KLO)/2.D0+B2*YA(KHI)/2.D0+
     $     ((-(1.D0+A2*A2)/4.D0+A2/2.D0)*Y2A(KLO)+
     $     (B2*B2/4.D0-B2/2.D0)*Y2A(KHI))*H*H/6.D0)*H
        KLO=KLO+1
        KHI=KHI+1
        GOTO 1
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION DNET6(DD,DM,ZHM,XMM,XM)
      IMPLICIT NONE
C       TURBOPAUSE CORRECTION FOR MSIS MODELS
C         Root mean density
C       8/20/80
C          DD - diffusive density
C          DM - full mixed density
C          ZHM - transition scale length
C          XMM - full mixed molecular weight
C          XM  - species molecular weight
C          DNET6 - combined density
      REAL*8 A,ZHM,XMM,XM,YLOG
      REAL*8 DM,DD,DNET6
c
      SAVE
      A=ZHM/(XMM-XM)
      IF(DM.GT.0.D0.AND.DD.GT.0.D0) GOTO 5
        WRITE(6,*) 'DNET6 LOG ERROR',DM,DD,XM
        IF(DD.EQ.0.D0.AND.DM.EQ.0.D0) DD=1.D0
        IF(DM.EQ.0.D0) GOTO 10
        IF(DD.EQ.0.D0) GOTO 20
    5 CONTINUE
      YLOG=A*DLOG(DM/DD)
      IF(YLOG.LT.-10.D0) GO TO 10
      IF(YLOG.GT.10.D0)  GO TO 20
        DNET6=DD*(1.D0+EXP(YLOG))**(1.D0/A)
        GO TO 50
   10 CONTINUE
        DNET6=DD
        GO TO 50
   20 CONTINUE
        DNET6=DM
        GO TO 50
   50 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR6(ALT, R,H1,ZH)
      IMPLICIT NONE
C        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C        ALT - altitude
C        R - target ratio
C        H1 - transition scale length
C        ZH - altitude of 1/2 R
      REAL*8 E,ALT,ZH,H1,EX,CCOR6,R
c
      SAVE
      E=(ALT-ZH)/H1
      IF(E.GT.70.D0) GO TO 20
      IF(E.LT.-70.D0) GO TO 10
        EX=EXP(E)
        CCOR6=R/(1.D0+EX)
        GO TO 50
   10   CCOR6=R
        GO TO 50
   20   CCOR6=0.D0
        GO TO 50
   50 CONTINUE
      CCOR6=EXP(CCOR6)
       RETURN
      END
C-----------------------------------------------------------------------
      BLOCK DATA GTD6BK
C          MSISE 90 12-MAR-90 
      INTEGER*4 IMR
c
      REAL*8  PT1(50),PT2(50),PT3(50),PA1(50),PA2(50),PA3(50)
      REAL*8  PB1(50),PB2(50),PB3(50),PC1(50),PC2(50),PC3(50)
      REAL*8  PD1(50),PD2(50),PD3(50),PE1(50),PE2(50),PE3(50)
      REAL*8  PF1(50),PF2(50),PF3(50),PG1(50),PG2(50),PG3(50)
      REAL*8  PH1(50),PH2(50),PH3(50),PI1(50),PI2(50),PI3(50)
      REAL*8  PJ1(50),PJ2(50),PJ3(50),PK1(50),PL1(50),PL2(50)
      REAL*8  PM1(50),PM2(50),PN1(50),PN2(50),PO1(50),PO2(50)
      REAL*8  PP1(50),PP2(50),PQ1(50),PQ2(50),PR1(50),PR2(50)
      REAL*8  PS1(50),PS2(50),PU1(50),PU2(50),PV1(50),PV2(50)
      REAL*8  PW1(50),PW2(50),PX1(50),PX2(50),PY1(50),PY2(50)
      REAL*8  PZ1(50),PZ2(50)
      REAL*8  PAVGM(10),PTM(10),PDM(10,8)
c
      CHARACTER NAME(2)*4,ISDATE(3)*4,ISTIME(2)*4
c
      COMMON/PARM6/PT1,PT2,PT3,PA1,PA2,PA3,
     $ PB1,PB2,PB3,PC1,PC2,PC3,
     $ PD1,PD2,PD3,PE1,PE2,PE3,
     $ PF1,PF2,PF3,PG1,PG2,PG3,
     $ PH1,PH2,PH3,PI1,PI2,PI3,
     $ PJ1,PJ2,PJ3,PK1,PL1,PL2,
     $ PM1,PM2,PN1,PN2,PO1,PO2,
     $ PP1,PP2,PQ1,PQ2,PR1,PR2,
     $ PS1,PS2,PU1,PU2,PV1,PV2,
     $ PW1,PW2,PX1,PX2,PY1,PY2,
     $ PZ1,PZ2
      COMMON/LOWER6/PTM,PDM
      COMMON/MAVG6/PAVGM
      COMMON/DATIM6/ISDATE,ISTIME,NAME
      COMMON/METSEL/IMR
      DATA IMR/0/
      DATA ISDATE/'12-M','AR-9','0   '/,ISTIME/'15:0','9:04'/
      DATA NAME/'MSIS','E 90'/
C         TEMPERATURE
      DATA PT1/
     *  9.96040D-01, 3.85528D-02, 3.03445D-03,-1.05531D-01,-6.07134D-03,
     * -5.16278D-04,-1.15622D-01, 2.02240D-03, 9.90156D-03,-1.27371D-01,
     * -3.02449D-02, 1.23512D-02,-5.26277D-03,-8.45398D+00, 0.00000D+00,
     *  1.42370D-02, 0.00000D+00, 1.25818D+02, 8.05486D-03, 1.64419D-03,
     * -6.21452D-06, 3.11701D-03, 0.00000D+00, 3.86578D-03, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,-6.41110D-06,
     *  0.00000D+00, 3.00150D+01, 5.33297D-03, 3.89146D-03, 2.04725D-03,
     *  0.00000D+00, 0.00000D+00,-1.92645D-02, 2.75905D+00, 1.47284D-03,
     *  3.41345D-04,-1.17388D-03,-3.54589D-04, 1.13139D-01, 1.69134D-01,
     *  5.08295D-03, 3.65016D-05, 4.26385D-03, 1.15102D-04, 5.11819D-03/
      DATA PT2/
     *  6.09108D-03, 4.04995D-05, 1.53049D-03, 2.41470D-05, 2.30764D-03,
     *  1.55267D-03, 1.33722D-03,-1.82318D-03,-2.63007D+02, 0.00000D+00,
     *  1.37337D-03, 9.95774D-04, 0.00000D+00,-1.08983D+02, 5.62606D-03,
     *  5.94053D-03, 1.09358D-03, 0.00000D+00,-1.33410D-02,-2.43409D-02,
     * -1.35688D-02, 3.11370D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -2.83023D+03, 8.45583D-04, 5.38706D-04, 0.00000D+00, 2.47956D+02,
     *  2.92246D-03, 0.00000D+00, 0.00000D+00, 7.47703D-05, 8.87993D-04,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -1.16540D-02,-4.49173D-03,-3.53189D-04,-1.73933D-04,-1.53218D-04,
     * -5.65411D-01, 7.77272D-03,-9.11784D+01, 6.45187D-04, 0.00000D+00/
      DATA PT3/
     * -8.37685D-04, 2.42318D-03, 4.73796D-03,-3.01801D-03,-4.23564D-03,
     * -2.48289D-03, 9.19286D-04, 2.16372D-03, 8.63968D-04, 1.89689D-03,
     *  4.15654D-03, 0.00000D+00, 1.18068D-02, 3.31190D-03, 0.00000D+00,
     *  1.20222D-03, 0.00000D+00, 0.00000D+00,-3.07246D+00, 0.00000D+00,
     *  0.00000D+00, 6.72403D-04, 1.08930D-03, 9.72278D-04, 4.68242D+00,
     * -3.15034D-04, 4.00059D-03, 5.15036D-03, 1.62989D-03, 1.08824D-03,
     *  9.95261D-04, 4.18955D+00,-3.64059D-01, 1.70182D-03, 0.00000D+00,
     *  0.00000D+00,-3.20120D+00, 0.00000D+00, 5.80206D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         HE DENSITY
      DATA PA1/
     *  1.04934D+00,-2.88362D-02,-2.07095D-01,-1.03314D-01,-7.02373D-03,
     *  1.29664D-02, 4.08853D-01,-9.19895D-03,-1.88660D-02, 1.40927D+00,
     *  1.75033D-01, 1.87351D-02, 1.10979D-01,-7.42871D+00, 0.00000D+00,
     *  2.67143D-01,-5.95979D-02, 1.05038D+02,-8.40963D-02,-6.97632D-04,
     *  2.06521D-06, 7.65306D-04, 0.00000D+00, 0.00000D+00, 1.26762D-01,
     *  1.28876D-01,-5.04479D-02,-1.30735D-02,-2.24348D-02, 0.00000D+00,
     *  0.00000D+00,-1.50832D+02,-6.29928D-03, 0.00000D+00,-4.07760D-03,
     *  0.00000D+00, 0.00000D+00, 5.25725D-02,-3.11486D+01,-3.13351D-03,
     *  2.75838D-03, 0.00000D+00, 0.00000D+00, 1.11247D-01, 1.08815D-01,
     * -4.66713D-02, 0.00000D+00,-3.29329D-03, 0.00000D+00, 1.67838D-03/
      DATA PA2/
     * -9.16691D-03, 3.45044D-05,-9.71806D-03, 0.00000D+00,-2.04672D-03,
     * -7.86899D-03,-7.98285D-03, 5.36515D-03,-5.31172D+03, 0.00000D+00,
     * -6.42781D-03,-1.71690D-03, 0.00000D+00,-6.79131D+01,-1.79912D-02,
     * -1.58305D-02,-7.12313D-03, 0.00000D+00, 2.53477D-02, 8.52960D-02,
     *  1.02163D-01, 2.95009D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -6.84625D+03,-6.19098D-03,-2.69289D-03, 0.00000D+00,-5.20231D+02,
     * -6.33463D-03, 0.00000D+00, 0.00000D+00,-6.02428D-03,-4.07077D-03,
     *  5.42264D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  4.07560D-02, 2.82288D-02, 9.08088D-03, 0.00000D+00, 0.00000D+00,
     * -4.05204D-01,-5.97931D-02,-7.31823D+01,-2.06620D-03, 0.00000D+00/
      DATA PA3/
     * -3.72723D-03,-1.88146D-02,-1.01794D-02, 8.04633D-03, 1.01090D-02,
     *  8.73253D-03, 2.38268D-02, 4.80444D-03, 1.71088D-03, 3.96369D-02,
     * -2.13809D-02, 0.00000D+00,-1.02588D-01,-5.91702D-03, 0.00000D+00,
     *  2.70923D-03, 0.00000D+00, 0.00000D+00,-1.75043D+02, 6.03489D-01,
     * -6.17589D-01, 8.38098D-03, 1.83871D-03,-7.05329D-04,-4.06644D+00,
     * -5.09347D-03,-2.84344D-02,-1.24160D-02, 1.33665D-02, 3.93410D-03,
     * -5.03723D-04,-4.57683D+00,-5.29542D-01,-4.25812D-03, 0.00000D+00,
     *  0.00000D+00, 1.91541D+01, 0.00000D+00, 3.23247D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         O DENSITY
      DATA PB1/
     *  9.31113D-01,-1.38721D-01,-1.33457D-01,-5.29542D-02,-4.44983D-03,
     *  1.35264D-02, 5.98075D-02,-3.62880D-02,-3.12798D-02, 3.72068D-01,
     *  2.95974D-02, 1.20509D-02, 5.21995D-02,-7.78888D+00, 0.00000D+00,
     *  1.18634D-01,-2.04495D-02, 1.03280D+02, 9.82432D-02, 4.77694D-04,
     *  0.00000D+00, 2.74372D-03, 0.00000D+00, 0.00000D+00, 7.57809D-02,
     *  1.71403D-01,-1.05205D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-8.73348D+00,-5.81094D-03, 0.00000D+00,-8.14944D-03,
     *  0.00000D+00, 0.00000D+00, 5.17255D-02,-1.53028D+01,-3.48932D-03,
     *  9.61771D-04, 5.57732D-03,-4.54180D-04, 9.88213D-02, 9.40456D-02,
     * -3.18797D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.32122D-03/
      DATA PB2/
     * -6.00220D-03, 2.77654D-05,-3.22019D-03, 0.00000D+00,-3.78551D-03,
     * -3.34809D-03,-1.70668D-03, 0.00000D+00, 6.36184D+03, 0.00000D+00,
     *  1.59986D-03,-3.88204D-03,-1.64825D-03,-7.47955D+01,-1.05360D-02,
     * -9.45723D-03,-1.59824D-03,-7.06730D-04,-1.68513D-02,-1.13023D-01,
     * -6.36637D-02,-1.37709D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -1.52368D+04,-5.86061D-03,-2.53108D-03, 0.00000D+00,-2.54837D+03,
     * -3.28988D-03, 0.00000D+00, 0.00000D+00,-2.76364D-03, 9.67923D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  4.34255D-02, 1.14020D-02,-6.18447D-03, 0.00000D+00, 0.00000D+00,
     * -3.02568D-01,-3.27694D-02,-6.71589D+01,-2.28340D-03, 0.00000D+00/
      DATA PB3/
     *  3.06230D-03,-4.65113D-03,-9.73421D-03, 1.28326D-02, 7.88553D-03,
     *  7.97197D-03,-1.20760D-02,-7.67547D-03,-1.20755D-03,-2.98523D-02,
     * -1.26560D-02, 0.00000D+00,-5.68350D-02,-1.53039D-02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.42911D-03,-4.01347D-03,-2.19074D-03, 3.11281D+00,
     *  3.23251D-03,-6.39523D-03,-6.63069D-03,-3.04403D-04,-4.01920D-03,
     * -1.18708D-03, 4.15211D+00,-2.01896D-01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         N2 DENSITY
      DATA PC1/
     *  1.06903D+00, 0.00000D+00, 0.00000D+00, 3.66210D-03, 0.00000D+00,
     *  1.90412D-02,-1.78929D-03, 0.00000D+00,-3.92257D-02,-1.19444D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-8.45398D+00, 0.00000D+00,
     *  2.08180D-02, 0.00000D+00, 1.39638D+02, 8.98481D-02, 0.00000D+00,
     *  0.00000D+00, 3.77113D-04, 0.00000D+00, 0.00000D+00, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-2.36325D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.43022D-03,
     * -3.99776D-06, 6.32343D-03, 5.48144D-03, 1.13139D-01, 1.69134D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PC2/
     *  0.00000D+00, 2.41470D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PC3/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         TLB
      DATA PD1/
     *  9.76619D-01, 0.00000D+00, 0.00000D+00,-2.00200D-02, 0.00000D+00,
     * -9.38391D-03,-1.95833D-03, 0.00000D+00, 1.31480D-02,-1.92414D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-8.45398D+00, 0.00000D+00,
     *  1.07674D-02, 0.00000D+00, 8.93820D+01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 5.68478D-04, 0.00000D+00, 0.00000D+00, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 4.66814D-03, 0.00000D+00, 0.00000D+00,
     *  5.11651D-05, 2.55717D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-2.60147D-03,-8.08556D-04, 1.13139D-01, 1.69134D-01,
     *  6.64196D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PD2/
     *  5.82026D-03, 2.41470D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 6.21998D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PD3/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         O2 DENSITY
      DATA PE1/
     *  9.31402D-01, 1.37976D-01, 0.00000D+00, 3.23736D-04, 0.00000D+00,
     * -9.10906D-03, 7.07506D-02, 0.00000D+00,-5.16650D-02, 6.89755D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-8.45398D+00, 0.00000D+00,
     *  2.81140D-02, 0.00000D+00, 7.36009D+01, 5.96604D-02, 0.00000D+00,
     *  0.00000D+00,-1.51792D-03, 0.00000D+00, 0.00000D+00, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 9.48758D+00, 8.84541D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 1.13139D-01, 1.69134D-01,
     *  1.45192D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PE2/
     *  1.07906D-02, 2.99942D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.48930D-02,
     * -7.87184D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -6.83420D-02,-4.41778D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.29730D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PE3/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         AR DENSITY
      DATA PF1/
     *  8.68053D-01, 2.36364D-01, 1.34306D-01, 1.03086D-02, 0.00000D+00,
     * -3.79164D-03,-1.57806D-01, 0.00000D+00,-5.87644D-02,-3.12508D-01,
     *  0.00000D+00, 4.37387D-02,-3.54091D-02,-2.23636D+01, 0.00000D+00,
     * -5.33976D-02, 0.00000D+00, 1.14091D+02, 5.17497D-02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 3.42702D+02, 1.57033D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-3.66278D-03,
     * -1.16193D-03, 0.00000D+00, 0.00000D+00, 1.13139D-01, 1.69134D-01,
     *  1.78431D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PF2/
     *  1.62864D-02, 3.16963D-05, 1.27968D-02, 0.00000D+00, 0.00000D+00,
     * -7.04599D-03, 2.07921D-03, 6.36660D-03, 2.29940D+04, 0.00000D+00,
     *  1.27833D-02,-2.08036D-03,-4.61820D-03,-6.29391D+01,-1.20745D-02,
     *  1.36675D-02, 1.36011D-02,-5.37162D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  1.92509D+04, 8.35522D-03, 4.19439D-03, 0.00000D+00, 1.20366D+04,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-1.00034D-02,-2.33267D-03,
     *  9.72374D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -2.65079D-02,-2.09125D-02,-1.09465D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.17252D-02,-7.12385D+01,-1.89428D-03, 0.00000D+00/
      DATA PF3/
     * -6.02006D-03, 1.69058D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.90646D-02,
     *  3.48971D-03, 0.00000D+00, 5.01174D-02, 5.50595D-02, 0.00000D+00,
     * -9.55897D-03, 0.00000D+00, 0.00000D+00,-1.51693D+03, 0.00000D+00,
     *  0.00000D+00, 1.29306D-02, 2.69567D-03, 0.00000D+00, 3.92243D+00,
     * -8.47690D-03, 1.16896D-02, 0.00000D+00, 1.48967D-02, 5.44521D-03,
     *  0.00000D+00, 5.64918D+00, 0.00000D+00,-7.72178D-03, 0.00000D+00,
     *  0.00000D+00,-7.34042D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          H DENSITY
      DATA PG1/
     *  1.27515D+00,-2.10472D-01,-1.77924D-01, 2.18900D-01, 2.88436D-02,
     *  1.90077D-02, 2.91001D-01, 2.17437D-02,-1.05186D-02, 4.36141D-01,
     *  1.07605D-01, 3.30755D-02, 4.00581D-02,-9.58051D+00, 0.00000D+00,
     *  1.54028D-02, 0.00000D+00, 7.34194D+01, 4.96540D-02,-5.95906D-03,
     *  3.84512D-05,-1.36000D-02, 0.00000D+00, 0.00000D+00, 1.32397D-01,
     *  2.13315D-01,-4.16610D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 1.46276D+02,-1.98408D-02, 0.00000D+00, 1.32530D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.04687D-04,
     * -1.47562D-03, 0.00000D+00, 0.00000D+00, 1.13139D-01, 1.69134D-01,
     * -1.26913D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,-6.08370D-03/
      DATA PG2/
     * -2.57587D-02, 3.19022D-05, 0.00000D+00, 0.00000D+00, 1.56644D-02,
     *  1.03640D-02, 1.05771D-03, 0.00000D+00, 3.57949D+03, 0.00000D+00,
     * -1.25672D-03, 1.52783D-03, 1.30518D-03, 7.55558D+00,-9.20341D-03,
     * -2.09142D-02,-1.34106D-02, 0.00000D+00,-4.83312D-02, 8.30900D-02,
     *  9.88009D-02,-1.41148D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -1.05513D+03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 6.73442D-03, 2.01691D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  5.98019D-02, 6.33298D-03,-1.12871D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.28604D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PG3/
     * -4.94960D-03,-1.36415D-02,-1.15039D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-5.86860D-03,-1.41732D-03, 2.13697D-03, 2.63845D+00,
     * -8.34186D-03,-1.87336D-02,-1.90870D-02,-8.03810D-03,-2.84279D-03,
     *  2.56722D-03, 1.71429D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          N DENSITY
      DATA PH1/
     *  5.73587D+01,-3.98747D-01, 0.00000D+00,-5.29554D-01,-5.82186D-03,
     *  7.14177D-02,-6.79279D-01,-1.67715D-01,-6.42434D-02,-2.11569D-01,
     * -1.59922D-01,-1.71024D-04,-1.15885D-01, 6.51603D+00, 0.00000D+00,
     * -1.76683D-01, 6.50395D-02, 1.43504D+00, 9.28208D-02, 5.11662D-03,
     *  0.00000D+00, 9.95121D-03, 0.00000D+00, 0.00000D+00, 1.32397D-01,
     *  2.13315D-01, 1.01451D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 5.67667D+01, 2.38192D-03, 0.00000D+00,-1.88240D-02,
     *  0.00000D+00, 0.00000D+00, 4.76218D-02, 2.35206D+01, 4.75901D-03,
     *  5.76162D-03, 1.51815D-02,-1.92730D-02, 1.13139D-01, 1.69134D-01,
     * -2.88771D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.18418D-03/
      DATA PH2/
     * -3.68927D-03, 3.14704D-05, 8.82198D-03, 0.00000D+00,-1.92562D-02,
     * -2.58674D-03,-2.19913D-02, 0.00000D+00, 4.38655D+03, 0.00000D+00,
     *  7.60126D-03, 2.59438D-03, 1.72310D-03, 7.79204D+01, 7.97786D-04,
     * -7.70510D-03, 1.90982D-03, 2.72707D-03, 1.01016D-02, 1.16537D-01,
     * -3.12236D-03, 1.39783D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -1.30712D+03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-3.20544D-03,-2.06970D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  1.59010D-02,-1.91427D-03,-3.42829D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-3.45379D-02, 8.94518D+01, 1.71556D-03, 0.00000D+00/
      DATA PH3/
     * -7.65278D-03,-2.08987D-04,-1.57393D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-8.60673D-03,-1.19922D-02,-6.46356D-03,-3.00107D+00,
     * -9.32511D-03,-1.50205D-02,-8.67835D-03,-7.64801D-03,-1.31495D-02,
     * -6.76720D-03,-1.82396D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C        SPARE
      DATA PI1/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-8.45398D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 1.13139D-01, 1.69134D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PI2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PI3/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          S PARAM  
      DATA PJ1/
     *  9.51363D-01,-4.67542D-02, 1.20260D-01, 0.00000D+00, 0.00000D+00,
     *  1.91357D-02, 0.00000D+00, 0.00000D+00, 1.25429D-03,-1.33240D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-8.45398D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.52317D-03, 0.00000D+00,-9.73404D-03, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-7.18482D-04, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 7.87683D-03,-2.33698D-03, 1.13139D-01, 1.69134D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PJ2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PJ3/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TURBO
      DATA PK1/
     *  9.33804D-01, 5.47446D+00, 1.53263D-01, 9.19303D-01, 1.64109D+01,
     *  4.27083D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.40925D-01,
     *  1.15897D+00, 4.71094D-01, 1.09459D+00, 5.25012D+00, 1.00000D+00,
     *  1.00000D+00, 1.03999D+00, 7.67132D-01, 1.10514D+00, 1.75636D+00,
     *  1.10845D+00, 2.33439D+00, 7.96532D-01, 4.31520D+00, 4.07300D+00,
     *  1.22807D+02, 2.39547D-01, 2.53791D-06, 8.42931D-01, 1.04192D+00,
     *  2.00202D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 9.62736D-01/
C         LOWER BOUNDARY
      DATA PTM/
     L  1.04130D+03, 3.86000D+02, 1.95000D+02, 1.66728D+01, 2.13000D+02,
     L  1.20000D+02, 2.40000D+02, 1.87000D+02,-2.00000D+00, 0.00000D+00/
      DATA PDM/
     L  2.45600D+07, 6.71072D-06, 1.00000D+02, 0.00000D+00, 1.10000D+02,
     L  1.00000D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  8.59400D+10, 5.40000D-01, 1.05000D+02,-8.00000D+00, 1.10000D+02,
     L  1.00000D+01, 9.00000D+01, 2.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  2.81000D+11, 0.00000D+00, 1.05000D+02, 2.80000D+01, 2.89500D+01,
     L  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  3.30000D+10, 2.68270D-01, 1.05000D+02, 0.00000D+00, 1.10000D+02,
     L  1.00000D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  1.33000D+09, 1.19615D-02, 1.05000D+02, 0.00000D+00, 1.10000D+02,
     L  1.00000D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  1.76100D+05, 1.00000D+00, 9.50000D+01,-8.00000D+00, 1.10000D+02,
     L  1.00000D+01, 9.00000D+01, 2.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  1.00000D+07, 1.00000D+00, 1.05000D+02,-8.00000D+00, 1.10000D+02,
     L  1.00000D+01, 9.00000D+01, 2.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  1.00000D+07, 1.00000D+00, 1.05000D+02,-8.00000D+00, 1.10000D+02,
     L  1.00000D+01, 9.00000D+01, 2.00000D+00, 0.00000D+00, 0.00000D+00/
C         TN1(2)
      DATA PL1/
     *  1.02083D+00, 4.08449D-02,-2.34582D-02, 4.38274D-04,-1.52380D-02,
     * -2.09089D-02, 4.46355D-03,-3.41250D-03,-1.12961D-02,-7.03277D-02,
     * -4.82724D-02, 0.00000D+00, 0.00000D+00,-6.20496D+00, 0.00000D+00,
     * -9.80197D-03,-1.45065D-02,-1.13226D+02, 2.28455D-02, 0.00000D+00,
     *  0.00000D+00, 4.93658D-04, 0.00000D+00, 3.79078D-03, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-8.89051D+03, 2.25900D-03, 1.76142D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.55015D-04,
     *  2.21388D-03,-5.99073D-04,-3.52331D-03, 1.13139D-01, 1.69134D-01,
     *  7.79156D-03,-1.93458D-03,-1.08596D-02,-4.39285D-04, 0.00000D+00/
      DATA PL2/
     *  3.83994D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 6.76608D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         TN1(3)
      DATA PM1/
     *  9.24880D-01, 7.41986D-02,-6.37629D-03, 6.00575D-03, 1.29382D-03,
     *  6.97550D-03,-1.70782D-03, 2.80584D-03,-8.87214D-03,-4.35703D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 4.31515D+00, 0.00000D+00,
     * -1.81474D-02,-6.06627D-02,-8.43503D+01, 8.46944D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-2.17081D-02,-2.19500D-03, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.47580D+02, 4.41585D-03, 7.80466D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 6.44155D-04,
     * -2.49166D-03, 2.90482D-03,-3.40501D-04, 1.13139D-01, 1.69134D-01,
     * -6.01460D-03,-1.63368D-03, 0.00000D+00,-4.31340D-03, 0.00000D+00/
      DATA PM2/
     *  4.53979D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-5.43660D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         TN1(4)
      DATA PN1/
     *  9.72669D-01,-4.26748D-02, 1.12876D-02,-8.44951D-03, 7.04114D-03,
     *  1.26036D-02,-3.88164D-03,-5.20509D-04,-6.09710D-04, 1.31603D-01,
     *  1.13804D-01, 0.00000D+00, 0.00000D+00,-6.15970D+00, 0.00000D+00,
     * -2.14214D-02,-6.62913D-02,-2.02884D-01, 2.35350D-02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 1.13573D-02,-1.84905D-03, 1.32397D-01,
     *  2.13315D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 1.42645D+00,-2.64405D-03,-5.57771D-04, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-2.20621D+01,-1.10313D-03,
     *  3.97063D-05, 5.47632D-05, 3.57577D-03, 1.13139D-01, 1.69134D-01,
     *  0.00000D+00, 1.18897D-03, 0.00000D+00, 7.62305D-04, 0.00000D+00/
      DATA PN2/
     * -3.52015D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-9.52550D-04,
     *  8.56253D-04, 4.33114D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.21223D-03,
     *  2.38694D-04, 9.15245D-04, 1.28385D-03, 8.67668D-04,-5.61425D-06,
     *  1.04445D+00, 3.41112D+01, 0.00000D+00,-8.40704D-01,-2.39639D+02,
     *  7.06668D-01,-2.05873D+01,-3.63696D-01, 2.39245D+01, 1.00000D+01,
     * -1.06657D-03,-7.67292D-04, 1.54534D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         TN1(5) TN2(1)
      DATA PO1/
     *  9.99368D-01, 4.33893D-02,-2.07009D-03, 1.09617D-03, 1.05440D-03,
     *  4.83408D-04, 9.77040D-04, 9.24791D-04, 4.80247D-04, 4.94737D-02,
     *  1.05985D-03, 0.00000D+00, 0.00000D+00, 2.74409D+00, 0.00000D+00,
     * -4.96656D-03,-1.51684D-02, 4.65158D+01,-7.51133D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 6.63808D-04, 1.32397D-01,
     *  2.13315D-01,-2.06652D-03,-6.32046D-03, 0.00000D+00, 0.00000D+00,
     *  5.94545D-03,-1.90958D+02, 0.00000D+00,-4.16892D-03, 0.00000D+00,
     * -1.67499D-02, 0.00000D+00, 2.58987D-03, 5.97781D+02, 0.00000D+00,
     *  0.00000D+00, 4.44890D-04, 4.66444D-04, 1.13139D-01, 1.69134D-01,
     *  0.00000D+00, 7.11360D-04, 1.32186D-02, 2.23948D-03, 0.00000D+00/
      DATA PO2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.60571D-03,
     *  6.28078D-04, 5.05469D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.57829D-03,
     * -4.00855D-04, 5.04077D-05,-1.39001D-03,-2.33406D-03,-4.81197D-04,
     *  1.46758D+00, 6.20332D+00, 0.00000D+00, 3.66476D-01,-6.19760D+01,
     *  3.09198D-01,-1.98999D+01, 0.00000D+00,-3.29933D+02, 0.00000D+00,
     * -1.10080D-03,-9.39310D-05, 1.39638D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TN2(2)
      DATA PP1/
     *  9.81637D-01,-1.41317D-03, 3.87323D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-3.58707D-02,
     * -8.63658D-03, 0.00000D+00, 0.00000D+00,-2.02226D+00, 0.00000D+00,
     * -8.69424D-03,-1.91397D-02, 8.76779D+01, 4.52188D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-7.07572D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -4.11210D-03, 3.50060D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  2.23760D-02, 0.00000D+00,-8.36657D-03, 1.61347D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.45130D-02, 0.00000D+00, 0.00000D+00/
      DATA PP2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.24152D-03,
     *  6.43365D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.33255D-03,
     *  2.42657D-03, 1.60666D-03,-1.85728D-03,-1.46874D-03,-4.79163D-06,
     *  1.22464D+00, 3.53510D+01, 0.00000D+00, 4.49223D-01,-4.77466D+01,
     *  4.70681D-01, 8.41861D+00,-2.88198D-01, 1.67854D+02, 0.00000D+00,
     *  7.11493D-04, 6.05601D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TN2(3)
      DATA PQ1/
     *  1.00422D+00,-7.11212D-03, 5.24480D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-5.28914D-02,
     * -2.41301D-02, 0.00000D+00, 0.00000D+00,-2.12219D+01, 0.00000D+00,
     * -3.28077D-03, 1.65727D-02, 1.68564D+00,-6.68154D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 8.42365D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-4.34645D-03,-1.03830D-02,-8.08279D-03, 2.16780D-02,
     *  0.00000D+00,-1.38459D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  1.45155D-02, 0.00000D+00, 7.04573D-03,-4.73204D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 1.08767D-02, 0.00000D+00, 0.00000D+00/
      DATA PQ2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 5.21769D-04,
     * -2.27387D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.26769D-03,
     *  3.16901D-03, 4.60316D-04,-1.01431D-04, 1.02131D-03, 9.96601D-04,
     *  1.25707D+00, 2.50114D+01, 0.00000D+00, 4.24472D-01,-2.77655D+01,
     *  3.44625D-01, 2.75412D+01, 0.00000D+00, 7.94251D+02, 0.00000D+00,
     *  2.45835D-03, 1.38871D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TN2(4) TN3(1)
      DATA PR1/
     *  1.01890D+00,-2.46603D-02, 1.00078D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-6.70977D-02,
     * -4.02286D-02, 0.00000D+00, 0.00000D+00,-2.29466D+01, 0.00000D+00,
     *  2.26580D-03, 2.63931D-02, 3.72625D+01,-6.39041D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.85291D-03,-7.47019D-03,-7.07265D-03, 0.00000D+00,
     *  0.00000D+00, 1.39717D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  9.58383D-03, 0.00000D+00, 9.19771D-03,-3.69121D+02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.57067D-02, 0.00000D+00, 0.00000D+00/
      DATA PR2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.92953D-03,
     * -2.77739D-03,-4.40092D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.47280D-03,
     *  2.95035D-04,-1.81246D-03, 2.81945D-03, 4.27296D-03, 9.78863D-04,
     *  1.40545D+00,-6.19173D+00, 0.00000D+00, 0.00000D+00,-7.93632D+01,
     *  4.44643D-01,-4.03085D+02, 0.00000D+00, 1.15603D+01, 0.00000D+00,
     *  2.25068D-03, 8.48557D-04,-2.98493D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TN3(2)
      DATA PS1/
     *  9.75801D-01, 3.80680D-02,-3.05198D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.85575D-02,
     *  5.04057D-02, 0.00000D+00, 0.00000D+00,-1.76046D+02, 0.00000D+00,
     * -1.48297D-03,-3.68560D-03, 3.02185D+01,-3.23338D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.15558D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 4.89620D-03, 1.44594D-02, 9.91215D-03,-1.00616D-02,
     * -8.21324D-03,-1.57757D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  1.53569D-02, 0.00000D+00, 6.63564D-03, 4.58410D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-2.51280D-02, 0.00000D+00, 0.00000D+00/
      DATA PS2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-8.73148D-04,
     * -1.29648D-03,-7.32026D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-4.68110D-03,
     * -4.66003D-03,-1.31567D-03,-7.39390D-04, 6.32499D-04,-4.65588D-04,
     * -1.29785D+00,-1.57139D+02, 0.00000D+00, 2.58350D-01,-3.69453D+01,
     *  4.10672D-01, 9.78196D+00,-1.52064D-01,-3.85084D+03, 0.00000D+00,
     * -8.52706D-04,-1.40945D-03,-7.26786D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TN3(3)
      DATA PU1/
     *  9.60722D-01, 7.03757D-02,-3.00266D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.22671D-02,
     *  4.10423D-02, 0.00000D+00, 0.00000D+00,-1.63070D+02, 0.00000D+00,
     *  5.40747D-04, 7.79481D-03, 1.44908D+02, 1.51484D-04, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.41844D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 5.77884D-03, 1.06073D-02, 5.36685D-03, 9.74319D-03,
     *  0.00000D+00,-2.88015D+03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  1.97547D-02, 0.00000D+00,-4.44902D-03,-2.92760D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 2.34419D-02, 0.00000D+00, 0.00000D+00/
      DATA PU2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-4.65325D-04,
     * -5.50628D-04, 3.31465D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.06179D-03,
     * -3.08575D-03,-7.93589D-04,-1.08629D-04, 5.95511D-04,-9.05050D-04,
     *  1.18997D+00, 4.15924D+01, 0.00000D+00,-4.72064D-01,-9.47150D+02,
     *  3.98723D-01, 1.98304D+01, 0.00000D+00, 3.73219D+03, 0.00000D+00,
     * -1.50040D-03,-1.14933D-03,-1.56769D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TN3(4)
      DATA PV1/
     *  1.03123D+00,-7.05124D-02, 8.71615D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-3.82621D-02,
     * -9.80975D-03, 0.00000D+00, 0.00000D+00, 2.89286D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 8.66153D+01, 7.91938D-04, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 4.68917D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 7.86638D-03, 9.57341D-03, 5.72268D-03, 9.90827D-03,
     *  0.00000D+00, 6.55573D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-4.00200D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 7.07457D-03, 0.00000D+00, 0.00000D+00/
      DATA PV2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.04970D-04,
     *  1.21560D-03,-8.05579D-06, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.49941D-03,
     * -4.57256D-04,-1.59311D-04, 2.96481D-04,-1.77318D-03,-6.37918D-04,
     *  1.02395D+00, 1.28172D+01, 0.00000D+00, 1.49903D-01,-2.63818D+01,
     *  0.00000D+00, 4.70628D+01,-2.22139D-01, 4.82292D-02, 0.00000D+00,
     * -8.67075D-04,-5.86479D-04, 5.32462D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TN3(5) SURFACD TDMP TSL
      DATA PW1/
     *  1.00828D+00,-9.10404D-02,-2.26549D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.32420D-02,
     * -9.08925D-03, 0.00000D+00, 0.00000D+00, 3.36105D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.24957D+01,-5.87939D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.79765D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 2.01237D+03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.75553D-02, 0.00000D+00, 0.00000D+00/
      DATA PW2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.29699D-03,
     *  1.26659D-03, 2.68402D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.17894D-03,
     *  1.48746D-03, 1.06478D-04, 1.34743D-04,-2.20939D-03,-6.23523D-04,
     *  6.36539D-01, 1.13621D+01, 0.00000D+00,-3.93777D-01, 2.38687D+03,
     *  0.00000D+00, 6.61865D+02,-1.21434D-01, 9.27608D+00, 0.00000D+00,
     *  1.68478D-04, 1.24892D-03, 1.71345D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TGN3(2) SURFACD GRAD TSLG
      DATA PX1/
     *  1.57293D+00,-6.78400D-01, 6.47500D-01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-7.62974D-02,
     * -3.60423D-01, 0.00000D+00, 0.00000D+00, 1.28358D+02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 4.68038D+01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.67898D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.90994D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 3.15706D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PX2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TGN2(1) TGN1(2)
      DATA PY1/
     *  8.66492D-01, 3.55807D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.12111D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 1.82458D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 1.01024D+02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 6.54251D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PY2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.56959D-02,
     *  1.91001D-02, 3.15971D-02, 1.00982D-02,-6.71565D-03, 2.57693D-03,
     *  1.38692D+00, 2.82132D-01, 0.00000D+00, 0.00000D+00, 3.81511D+02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          TGN3(1) TGN2(2)
      DATA PZ1/
     *  1.06029D+00,-5.25231D-02, 3.73034D-01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.31072D-02,
     * -3.88409D-01, 0.00000D+00, 0.00000D+00,-1.65295D+02, 0.00000D+00,
     * -4.38916D-02,-3.22716D-01,-8.82393D+01, 1.18458D-01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.19782D-01,-2.13801D-01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.62229D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -4.35863D-01, 0.00000D+00, 0.00000D+00,-5.37443D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-4.55788D-01, 0.00000D+00, 0.00000D+00/
      DATA PZ2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.84009D-02,
     *  3.96733D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 5.05494D-02,
     *  7.39617D-02, 1.92200D-02,-8.46151D-03,-1.34244D-02, 1.96338D-02,
     *  1.50421D+00, 1.88368D+01, 0.00000D+00, 0.00000D+00,-5.13114D+01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  5.11923D-02, 3.61225D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         MIDDLD ATMOSPHDRD AVDRAGDS
      DATA PAVGM/
     M  2.61000D+02, 2.64000D+02, 2.29000D+02, 2.17000D+02, 2.17000D+02,
     M  2.23000D+02, 2.86760D+02,-2.93940D+00, 2.50000D+00, 0.00000D+00/
      END

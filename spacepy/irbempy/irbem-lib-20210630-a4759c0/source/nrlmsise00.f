!***************************************************************************************************
! Copyright 2001 M. Picone
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
C-----------------------------------------------------------------------
      SUBROUTINE GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
      IMPLICIT NONE
C
C     NRLMSISE-00
C     -----------
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere
C
C        NEW FEATURES:
C          *Extensive satellite drag database used in model generation
C          *Revised O2 (and O) in lower thermosphere
C          *Additional nonlinear solar activity term
C          *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
C           At high altitudes (> 500 km), hot atomic oxygen or ionized
C           oxygen can become appreciable for some ranges of subroutine
C           inputs, thereby affecting drag on satellites and debris. We
C           group these species under the term "anomalous oxygen," since
C           their individual variations are not presently separable with
C           the drag data used to define this model component.
C
C        SUBROUTINES FOR SPECIAL OUTPUTS:
C
C        HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY
C        (SUBROUTINE GTD7D, OUTPUT D(6))
C           For atmospheric drag calculations at altitudes above 500 km,
C           call SUBROUTINE GTD7D to compute the "effective total mass
C           density" by including contributions from "anomalous oxygen."
C           See "NOTES ON OUTPUT VARIABLES" below on D(6).
C
C        PRESSURE GRID (SUBROUTINE GHP7)
C          See subroutine GHP7 to specify outputs at a pressure level
C          rather than at an altitude.
C
C        OUTPUT IN M-3 and KG/M3:   CALL METERS7(.TRUE.)
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
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
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES:
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C     NOTES ON OUTPUT VARIABLES:
C        TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS7(.TRUE.)
C
C        O, H, and N are set to zero below 72.5 km
C
C        T(1), Exospheric temperature, is set to global average for
C        altitudes below 120 km. The 120 km gradient is left at global
C        average value for altitudes below 72 km.
C
C        D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7
C        and GTD7D
C
C          SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
C          species labeled by indices 1-5 and 7-8 in output variable D.
C          This includes He, O, N2, O2, Ar, H, and N but does NOT include
C          anomalous oxygen (species index 9).
C
C          SUBROUTINE GTD7D -- D(6) is the "effective total mass density
C          for drag" and is the sum of the mass densities of all species
C          in this model, INCLUDING anomalous oxygen.
C
C     SWITCHES: The following is for test and special purposes:
C
C        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC7(SW),
C        WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
C        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C        FOR THE FOLLOWING VARIATIONS
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
C        To get current values of SW: CALL TRETRV7(SW)
C
      INTEGER*4 IYD,MASS,IMR,ISW,MN3,MN2,MSSL,MSS
      INTEGER*4 I,J
      REAL*8 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),D(9),T(2)
      REAL*8 TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
      REAL*8 G0,RL,DD,DB14,TR12,PTM(10),PDM(10,8)
      REAL*8 TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      REAL*8 PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),PMA(100,10)
      REAL*8 SAM(100)
      REAL*8 D6(8),T6(2)
      REAL*8 ZN3(5),ZN2(4),SV(25),SW(25),SWC(25)
      REAL*8 ZMIX,ALAST
      REAL*8 PAVGM(10),GSURF,RE
      REAL*8 DM04,DM16,DM28,DM32,DM40,DM01,DM14
      REAL*8 XLAT,V1,VTST7,XMM,ALTT,DM28M,GLOB7S,DMC,DZ28,DMR
      REAL*8 DENSM7,TZ,DS(9),TS(2)
c
      CHARACTER NAME(2)*4,ISDATE(3)*4,ISTIME(2)*4
      CHARACTER ISD(3)*4,IST(2)*4,NAM(2)*4
c
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO7/TN1,TN2,TN3,TGN1,TGN2,TGN3
      COMMON/LOWER7/PTM,PDM
      COMMON/PARM7/PT,PD,PS,PDL,PTL,PMA,SAM
      COMMON/DATIM7/ISD,IST,NAM
      COMMON/DATIME/ISDATE,ISTIME,NAME
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/MAVG7/PAVGM
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL7/IMR
      SAVE
      EXTERNAL GTD7BK
      DATA MN3/5/,ZN3/32.5D0,20.D0,15.D0,10.D0,0.D0/
      DATA MN2/4/,ZN2/72.5D0,55.D0,45.D0,32.5D0/
      DATA ZMIX/62.5D0/,ALAST/99999.D0/,MSSL/-999/
      DATA SV/25*1.D0/
      IF(ISW.NE.64999) CALL TSELEC7(SV)
C      Put identification data into common/datime/
      DO 1 I=1,3
        ISDATE(I)=ISD(I)
    1 CONTINUE
      DO 2 I=1,2
        ISTIME(I)=IST(I)
        NAME(I)=NAM(I)
    2 CONTINUE
C
C        Test for changed input
      V1=VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,1)
C       Latitude variation of gravity (none for SW(2)=0)
      XLAT=GLAT
      IF(SW(2).EQ.0.D0) XLAT=45.D0
      CALL GLATF7(XLAT,GSURF,RE)
C
      XMM=PDM(5,3)
C
C       THERMOSPHERE/MESOSPHERE (above ZN2(1))
      ALTT=DMAX1(ALT,ZN2(1))
      MSS=MASS
C       Only calculate N2 in thermosphere if alt in mixed region
      IF(ALT.LT.ZMIX.AND.MASS.GT.0) MSS=28
C       Only calculate thermosphere if input parameters changed
C         or altitude above ZN2(1) in mesosphere
      IF(V1.EQ.1.D0.OR.ALT.GT.ZN2(1).OR.ALAST.GT.ZN2(1).OR.MSS.NE.MSSL)
     $ THEN
        CALL GTS7(IYD,SEC,ALTT,GLAT,GLONG,STL,F107A,F107,AP,MSS,DS,TS)
        DM28M=DM28
C         metric adjustment
        IF(IMR.EQ.1) DM28M=DM28*1.D6
        MSSL=MSS
      ENDIF
      T(1)=TS(1)
      T(2)=TS(2)
      IF(ALT.GE.ZN2(1)) THEN
        DO 5 J=1,9
          D(J)=DS(J)
    5   CONTINUE
        GOTO 10
      ENDIF
C
C       LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3(1) and ZN2(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
C         Only calculate nodes if input changed
       IF(V1.EQ.1.D0.OR.ALAST.GE.ZN2(1)) THEN
        TGN2(1)=TGN1(2)
        TN2(1)=TN1(5)
        TN2(2)=PMA(1,1)*PAVGM(1)/(1.D0-SW(20)*GLOB7S(PMA(1,1)))
        TN2(3)=PMA(1,2)*PAVGM(2)/(1.D0-SW(20)*GLOB7S(PMA(1,2)))
        TN2(4)=PMA(1,3)*PAVGM(3)/(1.D0-SW(20)*SW(22)*GLOB7S(PMA(1,3)))
        TGN2(2)=PAVGM(9)*PMA(1,10)*(1.D0+SW(20)*SW(22)*
     &   GLOB7S(PMA(1,10)))
     $  *TN2(4)*TN2(4)/(PMA(1,3)*PAVGM(3))**2
        TN3(1)=TN2(4)
       ENDIF
       IF(ALT.GE.ZN3(1)) GOTO 6
C
C       LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
C         Only calculate nodes if input changed
        IF(V1.EQ.1.D0.OR.ALAST.GE.ZN3(1)) THEN
         TGN3(1)=TGN2(2)
         TN3(2)=PMA(1,4)*PAVGM(4)/(1.D0-SW(22)*GLOB7S(PMA(1,4)))
         TN3(3)=PMA(1,5)*PAVGM(5)/(1.D0-SW(22)*GLOB7S(PMA(1,5)))
         TN3(4)=PMA(1,6)*PAVGM(6)/(1.D0-SW(22)*GLOB7S(PMA(1,6)))
         TN3(5)=PMA(1,7)*PAVGM(7)/(1.D0-SW(22)*GLOB7S(PMA(1,7)))
         TGN3(2)=PMA(1,8)*PAVGM(8)*(1.D0+SW(22)*GLOB7S(PMA(1,8)))
     $   *TN3(5)*TN3(5)/(PMA(1,7)*PAVGM(7))**2
        ENDIF
    6   CONTINUE
        IF(MASS.EQ.0) GOTO 50
C          LINEAR TRANSITION TO FULL MIXING BELOW ZN2(1)
        DMC=0
        IF(ALT.GT.ZMIX) DMC=1.D0-(ZN2(1)-ALT)/(ZN2(1)-ZMIX)
        DZ28=DS(3)
C      ***** N2 DENSITY ****
        DMR=DS(3)/DM28M-1.D0
        D(3)=DENSM7(ALT,DM28M,XMM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
        D(3)=D(3)*(1.D0+DMR*DMC)
C      ***** HE DENSITY ****
        D(1)=0
        IF(MASS.NE.4.AND.MASS.NE.48) GOTO 204
          DMR=DS(1)/(DZ28*PDM(2,1))-1.D0
          D(1)=D(3)*PDM(2,1)*(1.D0+DMR*DMC)
  204   CONTINUE
C      **** O DENSITY ****
        D(2)=0.D0
        D(9)=0.D0
  216   CONTINUE
C      ***** O2 DENSITY ****
        D(4)=0.D0
        IF(MASS.NE.32.AND.MASS.NE.48) GOTO 232
          DMR=DS(4)/(DZ28*PDM(2,4))-1.D0
          D(4)=D(3)*PDM(2,4)*(1.D0+DMR*DMC)
  232   CONTINUE
C      ***** AR DENSITY ****
        D(5)=0.D0
        IF(MASS.NE.40.AND.MASS.NE.48) GOTO 240
          DMR=DS(5)/(DZ28*PDM(2,5))-1.D0
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
      DD=DENSM7(ALT,1.D0,0.D0,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
      T(2)=TZ
   90 CONTINUE
      ALAST=ALT
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTD7D(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,
     $ D,T)
      IMPLICIT NONE
C
C     NRLMSISE-00
C     -----------
C        This subroutine provides Effective Total Mass Density for
C        output D(6) which includes contributions from "anomalous
C        oxygen" which can affect satellite drag above 500 km.  This
C        subroutine is part of the distribution package for the
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere.  See subroutine GTD7 for more extensive comments.
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
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
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES:
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3) [includes anomalous oxygen]
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      INTEGER*4 IMR,IYD,MASS
      REAL*8 D(9),T(2),AP(7)
      REAL*8 SEC,ALT,GLAT,GLONG,STL,F107A,F107
      COMMON/METSEL7/IMR
      CALL GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C       TOTAL MASS DENSITY
C
        IF(MASS.EQ.48) THEN
         D(6) = 1.66D-24*(4.D0*D(1)+16.D0*D(2)+28.D0*D(3)+32.D0*D(4)+
     &       40.D0*D(5)+D(7)+14.D0*D(8)+16.D0*D(9))
         IF(IMR.EQ.1) D(6)=D(6)/1000.D0
         ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GHP7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,
     $  D,T,PRESS)
      IMPLICIT NONE
C       FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD7
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
C        D(9) - HOT O NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      INTEGER*4 IYD,IMR,IDAY,L,LTEST
      REAL*8 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),D(9),T(2),PRESS
      REAL*8 GSURF,RE,BM,RGAS,TEST,PL,ZI,CL,CL2,CD,CA,Z,XN,P,DIFF,XM
      REAL*8 G,SH
c
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL7/IMR
      SAVE
      DATA BM/1.3806D-19/,RGAS/831.4D0/
      DATA TEST/.00043D0/,LTEST/12/
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
     &     (-2.93D0+1.11D0)
         Z=ZI-4.87D0*CL*CD*CA-1.64D0*CL2*CA+.31D0*CA*CL
      ENDIF
      IF(PL.LT.-5.D0) Z=22.D0*(PL+4.D0)**2+110.D0
C      ITERATION LOOP
      L=0
   10 CONTINUE
        L=L+1
        CALL GTD7(IYD,SEC,Z,GLAT,GLONG,STL,F107A,F107,AP,48,D,T)
        XN=D(1)+D(2)+D(3)+D(4)+D(5)+D(7)+D(8)
        P=BM*XN*T(2)
        IF(IMR.EQ.1) P=P*1.D-6
        DIFF=PL-DLOG10(P)
        IF(ABS(DIFF).LT.TEST .OR. L.EQ.LTEST) GOTO 20
        XM=D(6)/XN/1.66D-24
        IF(IMR.EQ.1) XM = XM*1.D3
        G=GSURF/(1.D0+Z/RE)**2
        SH=RGAS*T(2)/(XM*G)
C         New altitude estimate using scale height
        IF(L.LT.6) THEN
          Z=Z-SH*DIFF*2.302D0
        ELSE
          Z=Z-SH*DIFF
        ENDIF
        GOTO 10
   20 CONTINUE
      IF(L.EQ.LTEST) WRITE(6,100) PRESS,DIFF
  100 FORMAT(1X,29HGHP7 NOT CONVERGING FOR PRESS, 1PE12.2,E12.2)
      ALT=Z
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLATF7(LAT,GV,REFF)
      IMPLICIT NONE
C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
C      RADIUS (REFF)
      REAL*8 LAT,LATL,GV,REFF,DGTR,C2
      SAVE
      DATA DGTR/1.74533D-2/
      C2 = COS(2.D0*DGTR*LAT)
      GV = 980.616D0*(1.D0-.0026373D0*C2)
      REFF = 2.D0*GV/(3.085462D-6 + 2.27D-9*C2)*1.D-5
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,IC)
      IMPLICIT NONE
C       Test if geophysical variables or switches changed and save
C       Return 0 if unchanged and 1 if changed
      INTEGER*4 IYD,IC,IYDL(2),ISW,I
      REAL*8 SEC,GLAT,GLONG,STL,F107A,F107,AP(7)
      REAL*8 SECL(2),GLATL(2),GLL(2),STLL(2)
      REAL*8 FAL(2),FL(2),APL(7,2),SWL(25,2),SWCL(25,2)
      REAL*8 SW(25),SWC(25),VTST7
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
c
      SAVE
      DATA IYDL/2*-999/,SECL/2*-999.D0/,GLATL/2*-999.D0/,GLL/2*-999.D0/
      DATA STLL/2*-999.D0/,FAL/2*-999.D0/,FL/2*-999.D0/,APL/14*-999.D0/
      DATA SWL/50*-999.D0/,SWCL/50*-999.D0/
      VTST7=0.D0
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
      VTST7=1.D0
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
      SUBROUTINE GTS7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
      IMPLICIT NONE
C
C     Thermospheric portion of NRLMSISE-00
C     See GTD7 for more extensive comments
C
C        OUTPUT IN M-3 and KG/M3:   CALL METERS7(.TRUE.)
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM) (>72.5 km)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
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
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES:
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3) [Anomalous O NOT included]
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      INTEGER*4 MT(11),ISW,MASS,IMR,MN1,J
      INTEGER*4 IYD,I
      LOGICAL METER
      REAL*8 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),D(9),T(2)
      REAL*8 TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
      REAL*8 G0,RL,DD,DB14,TR12,PTM(10),PDM(10,8)
      REAL*8 TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      REAL*8 PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),PMA(100,10)
      REAL*8 SAM(100)
      REAL*8 SW(25),SWC(25),TINFG,GB,ROUT,TT(15)
      REAL*8 ALTL(8),ZN1(5)
      REAL*8 DGTR,DR,ALAST
      REAL*8 V2,VTST7,TINF,DAY,ZHF,YRD,XMM,XMD,GLOB7S,ZLB,TZ
      REAL*8 GLOBE7,DENSU7,DNET7,CCOR7,DDUM
      REAL*8 G28,ZH28,ZHM28,B28,DM28
      REAL*8 G4,ZH04,ZHM04,B04,DM04,ZC04,HC04
      REAL*8 G16,ZH16,ZHM16,B16,DM16,ZC16,HC16,HCC16,ZCC16,RC16
      REAL*8 G32,ZH32,ZHM32,B32,DM32,HC32,ZC32
      REAL*8 G40,ZH40,ZHM40,B40,DM40,HC40,ZC40
      REAL*8 G1,ZH01,ZHM01,B01,DM01,HC01,ZC01,HCC01,ZCC01,RC01
      REAL*8 G14,ZH14,ZHM14,B14,DM14,HC14,ZC14,HCC14,ZCC14,RC14
      REAL*8 ALPHA(9),ZSHO,ZMHO,ZSHT,THO,SCALH,DB16H,G16H,RC32,ZCC32
      REAL*8 HCC232,HCC32,CCOR2,HC216,Z,T2
c
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO7/TN1,TN2,TN3,TGN1,TGN2,TGN3
      COMMON/LOWER7/PTM,PDM
      COMMON/PARM7/PT,PD,PS,PDL,PTL,PMA,SAM
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/TTEST/TINFG,GB,ROUT,TT
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/METSEL7/IMR
      SAVE
      DATA MT/48,0,4,16,28,32,40,1,49,14,17/
      DATA ALTL/200.D0,300.D0,160.D0,250.D0,240.D0,450.D0,320.D0,450.D0/
      DATA MN1/5/,ZN1/120.D0,110.D0,100.D0,90.D0,72.5D0/
      DATA DGTR/1.74533D-2/,DR/1.72142D-2/,ALAST/-999.D0/
      DATA ALPHA/-0.38D0,0.D0,0.D0,0.D0,0.17D0,0.D0,-0.38D0,0.D0,0.D0/
C        Test for changed input
      V2=VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,2)
C
      YRD=IYD*1.D0
      ZA=PDL(16,2)
      ZN1(1)=ZA
      DO 2 J=1,9
        D(J)=0.D0
    2 CONTINUE
C        TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
      IF(ALT.GT.ZN1(1)) THEN
        IF(V2.EQ.1.D0.OR.ALAST.LE.ZN1(1)) TINF=PTM(1)*PT(1)
     $  *(1.D0+SW(16)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PT))
      ELSE
        TINF=PTM(1)*PT(1)
      ENDIF
      T(1)=TINF
C          GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
      IF(ALT.GT.ZN1(5)) THEN
        IF(V2.EQ.1.D0.OR.ALAST.LE.ZN1(5)) G0=PTM(4)*PS(1)
     $   *(1.D0+SW(19)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PS))
      ELSE
        G0=PTM(4)*PS(1)
      ENDIF
C      Calculate these temperatures only if input changed
      IF(V2.EQ.1.D0 .OR. ALT.LT.300.D0)
     $  TLB=PTM(2)*(1.D0+SW(17)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,
     $  F107A,F107,AP,PD(1,4)))*PD(1,4)
       S=G0/(TINF-TLB)
C       Lower thermosphere temp variations not significant for
C        density above 300 km
       IF(ALT.LT.300.D0) THEN
        IF(V2.EQ.1.D0.OR.ALAST.GE.300.D0) THEN
         TN1(2)=PTM(7)*PTL(1,1)/(1.D0-SW(18)*GLOB7S(PTL(1,1)))
         TN1(3)=PTM(3)*PTL(1,2)/(1.D0-SW(18)*GLOB7S(PTL(1,2)))
         TN1(4)=PTM(8)*PTL(1,3)/(1.D0-SW(18)*GLOB7S(PTL(1,3)))
         TN1(5)=PTM(5)*PTL(1,4)/(1.D0-SW(18)*SW(20)*GLOB7S(PTL(1,4)))
         TGN1(2)=PTM(9)*PMA(1,9)*(1.D0+SW(18)*SW(20)*GLOB7S(PMA(1,9)))
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
      TR12=1.D0
C
      IF(MASS.EQ.0) GO TO 50
C       N2 variation factor at Zlb
      G28=SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,
     & AP,PD(1,3))
      DAY=DMOD(YRD,1000.D0)
C        VARIATION OF TURBOPAUSE HEIGHT
      ZHF=PDL(25,2)
     $    *(1.D0+SW(5)*PDL(25,1)*SIN(DGTR*GLAT)*COS(DR*(DAY-PT(14))))
      YRD=IYD*1.D0
      T(1)=TINF
      XMM=PDM(5,3)
      Z=ALT
C
      DO 10 J = 1,11
      IF(MASS.EQ.MT(J))   GO TO 15
   10 CONTINUE
      WRITE(6,100) MASS
      GO TO 90
   15 IF(Z.GT.ALTL(6).AND.MASS.NE.28.AND.MASS.NE.48) GO TO 17
C
C       **** N2 DENSITY ****
C
C      Diffusive density at Zlb
      DB28 = PDM(1,3)*EXP(G28)*PD(1,3)
C      Diffusive density at Alt
      D(3)=DENSU7(Z,DB28,TINF,TLB, 28.D0,ALPHA(3),T(2),PTM(6),S,MN1,ZN1,
     & TN1,TGN1)
      DD=D(3)
C      Turbopause
      ZH28=PDM(3,3)*ZHF
      ZHM28=PDM(4,3)*PDL(6,2)
      XMD=28.D0-XMM
C      Mixed density at Zlb
      B28=DENSU7(ZH28,DB28,TINF,TLB,XMD,ALPHA(3)-1.D0,TZ,PTM(6),S,MN1,
     & ZN1,TN1,TGN1)
      IF(Z.GT.ALTL(3).OR.SW(15).EQ.0.D0) GO TO 17
C      Mixed density at Alt
      DM28=DENSU7(Z,B28,TINF,TLB,XMM,ALPHA(3),TZ,PTM(6),S,MN1,
     & ZN1,TN1,TGN1)
C      Net density at Alt
      D(3)=DNET7(D(3),DM28,ZHM28,XMM,28.D0)
   17 CONTINUE
      GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
   20 CONTINUE
C
C       **** HE DENSITY ****
C
C       Density variation factor at Zlb
      G4 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,1))
C      Diffusive density at Zlb
      DB04 = PDM(1,1)*EXP(G4)*PD(1,1)
C      Diffusive density at Alt
      D(1)=DENSU7(Z,DB04,TINF,TLB, 4.D0,ALPHA(1),T(2),PTM(6),S,MN1,ZN1,
     & TN1,TGN1)
      DD=D(1)
      IF(Z.GT.ALTL(1).OR.SW(15).EQ.0.D0) GO TO 24
C      Turbopause
      ZH04=PDM(3,1)
C      Mixed density at Zlb
      B04=DENSU7(ZH04,DB04,TINF,TLB,4.D0-XMM,ALPHA(1)-1.D0,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM04=DENSU7(Z,B04,TINF,TLB,XMM,0.D0,T(2),PTM(6),
     &S,MN1,ZN1,TN1,TGN1)
      ZHM04=ZHM28
C      Net density at Alt
      D(1)=DNET7(D(1),DM04,ZHM04,XMM,4.D0)
C      Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,1)/B04)
      ZC04=PDM(5,1)*PDL(1,2)
      HC04=PDM(6,1)*PDL(2,2)
C      Net density corrected at Alt
      D(1)=D(1)*CCOR7(Z,RL,HC04,ZC04)
   24 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   25 CONTINUE
C
C      **** O DENSITY ****
C
C       Density variation factor at Zlb
      G16= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,2))
C      Diffusive density at Zlb
      DB16 =  PDM(1,2)*EXP(G16)*PD(1,2)
C       Diffusive density at Alt
      D(2)=DENSU7(Z,DB16,TINF,TLB, 16.D0,ALPHA(2),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(2)
      IF(Z.GT.ALTL(2).OR.SW(15).EQ.0.D0) GO TO 34
C  Corrected from PDM(3,1) to PDM(3,2)  12/2/85
C       Turbopause
      ZH16=PDM(3,2)
C      Mixed density at Zlb
      B16=DENSU7(ZH16,DB16,TINF,TLB,16.D0-XMM,ALPHA(2)-1.D0,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM16=DENSU7(Z,B16,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,
     &MN1,ZN1,TN1,TGN1)
      ZHM16=ZHM28
C      Net density at Alt
      D(2)=DNET7(D(2),DM16,ZHM16,XMM,16.D0)
C   3/16/99 Change form to match O2 departure from diff equil near 150
C   km and add dependence on F10.7
C      RL=DLOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
      RL=PDM(2,2)*PDL(17,2)*(1.D0+SW(1)*PDL(24,1)*(F107A-150.D0))
      HC16=PDM(6,2)*PDL(4,2)
      ZC16=PDM(5,2)*PDL(3,2)
      HC216=PDM(6,2)*PDL(5,2)
      D(2)=D(2)*CCOR2(Z,RL,HC16,ZC16,HC216)
C       Chemistry correction
      HCC16=PDM(8,2)*PDL(14,2)
      ZCC16=PDM(7,2)*PDL(13,2)
      RC16=PDM(4,2)*PDL(15,2)
C      Net density corrected at Alt
      D(2)=D(2)*CCOR7(Z,RC16,HCC16,ZCC16)
   34 CONTINUE
      IF(MASS.NE.48.AND.MASS.NE.49) GO TO 90
   35 CONTINUE
C
C       **** O2 DENSITY ****
C
C       Density variation factor at Zlb
      G32= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,5))
C      Diffusive density at Zlb
      DB32 = PDM(1,4)*EXP(G32)*PD(1,5)
C       Diffusive density at Alt
      D(4)=DENSU7(Z,DB32,TINF,TLB, 32.D0,ALPHA(4),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      IF(MASS.EQ.49) THEN
         DD=DD+2.D0*D(4)
      ELSE
         DD=D(4)
      ENDIF
      IF(SW(15).EQ.0.D0) GO TO 39
      IF(Z.GT.ALTL(4)) GO TO 38
C       Turbopause
      ZH32=PDM(3,4)
C      Mixed density at Zlb
      B32=DENSU7(ZH32,DB32,TINF,TLB,32.D0-XMM,ALPHA(4)-1.D0,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM32=DENSU7(Z,B32,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,
     &MN1,ZN1,TN1,TGN1)
      ZHM32=ZHM28
C      Net density at Alt
      D(4)=DNET7(D(4),DM32,ZHM32,XMM,32.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,4)/B32)
      HC32=PDM(6,4)*PDL(8,2)
      ZC32=PDM(5,4)*PDL(7,2)
      D(4)=D(4)*CCOR7(Z,RL,HC32,ZC32)
   38 CONTINUE
C      Correction for general departure from diffusive equilibrium above Zlb
      HCC32=PDM(8,4)*PDL(23,2)
      HCC232=PDM(8,4)*PDL(23,1)
      ZCC32=PDM(7,4)*PDL(22,2)
      RC32=PDM(4,4)*PDL(24,2)*(1.D0+SW(1)*PDL(24,1)*(F107A-150.D0))
C      Net density corrected at Alt
      D(4)=D(4)*CCOR2(Z,RC32,HCC32,ZCC32,HCC232)
   39 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   40 CONTINUE
C
C       **** AR DENSITY ****
C
C       Density variation factor at Zlb
      G40= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,6))
C      Diffusive density at Zlb
      DB40 = PDM(1,5)*EXP(G40)*PD(1,6)
C       Diffusive density at Alt
      D(5)=DENSU7(Z,DB40,TINF,TLB, 40.D0,ALPHA(5),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(5)
      IF(Z.GT.ALTL(5).OR.SW(15).EQ.0.D0) GO TO 44
C       Turbopause
      ZH40=PDM(3,5)
C      Mixed density at Zlb
      B40=DENSU7(ZH40,DB40,TINF,TLB,40.D0-XMM,ALPHA(5)-1.D0,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM40=DENSU7(Z,B40,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,
     &MN1,ZN1,TN1,TGN1)
      ZHM40=ZHM28
C      Net density at Alt
      D(5)=DNET7(D(5),DM40,ZHM40,XMM,40.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,5)/B40)
      HC40=PDM(6,5)*PDL(10,2)
      ZC40=PDM(5,5)*PDL(9,2)
C      Net density corrected at Alt
      D(5)=D(5)*CCOR7(Z,RL,HC40,ZC40)
   44 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   45 CONTINUE
C
C        **** HYDROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G1 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,7))
C      Diffusive density at Zlb
      DB01 = PDM(1,6)*EXP(G1)*PD(1,7)
C       Diffusive density at Alt
      D(7)=DENSU7(Z,DB01,TINF,TLB,1.D0,ALPHA(7),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(7)
      IF(Z.GT.ALTL(7).OR.SW(15).EQ.0.D0) GO TO 47
C       Turbopause
      ZH01=PDM(3,6)
C      Mixed density at Zlb
      B01=DENSU7(ZH01,DB01,TINF,TLB,1.D0-XMM,ALPHA(7)-1.D0,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM01=DENSU7(Z,B01,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,
     &MN1,ZN1,TN1,TGN1)
      ZHM01=ZHM28
C      Net density at Alt
      D(7)=DNET7(D(7),DM01,ZHM01,XMM,1.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,6)*ABS(PDL(18,2))/B01)
      HC01=PDM(6,6)*PDL(12,2)
      ZC01=PDM(5,6)*PDL(11,2)
      D(7)=D(7)*CCOR7(Z,RL,HC01,ZC01)
C       Chemistry correction
      HCC01=PDM(8,6)*PDL(20,2)
      ZCC01=PDM(7,6)*PDL(19,2)
      RC01=PDM(4,6)*PDL(21,2)
C      Net density corrected at Alt
      D(7)=D(7)*CCOR7(Z,RC01,HCC01,ZCC01)
   47 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   48 CONTINUE
C
C        **** ATOMIC NITROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G14 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,8))
C      Diffusive density at Zlb
      DB14 = PDM(1,7)*EXP(G14)*PD(1,8)
C       Diffusive density at Alt
      D(8)=DENSU7(Z,DB14,TINF,TLB,14.D0,ALPHA(8),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(8)
      IF(Z.GT.ALTL(8).OR.SW(15).EQ.0.D0) GO TO 49
C       Turbopause
      ZH14=PDM(3,7)
C      Mixed density at Zlb
      B14=DENSU7(ZH14,DB14,TINF,TLB,14.D0-XMM,ALPHA(8)-1.D0,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM14=DENSU7(Z,B14,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,
     &MN1,ZN1,TN1,TGN1)
      ZHM14=ZHM28
C      Net density at Alt
      D(8)=DNET7(D(8),DM14,ZHM14,XMM,14.D0)
C       Correction to specified mixing ratio at ground
      RL=DLOG(B28*PDM(2,7)*ABS(PDL(3,1))/B14)
      HC14=PDM(6,7)*PDL(2,1)
      ZC14=PDM(5,7)*PDL(1,1)
      D(8)=D(8)*CCOR7(Z,RL,HC14,ZC14)
C       Chemistry correction
      HCC14=PDM(8,7)*PDL(5,1)
      ZCC14=PDM(7,7)*PDL(4,1)
      RC14=PDM(4,7)*PDL(6,1)
C      Net density corrected at Alt
      D(8)=D(8)*CCOR7(Z,RC14,HCC14,ZCC14)
   49 CONTINUE
      IF(MASS.NE.48) GO TO 90
   46 CONTINUE
C
C        **** Anomalous OXYGEN DENSITY ****
C
      G16H = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,9))
      DB16H = PDM(1,8)*EXP(G16H)*PD(1,9)
      THO=PDM(10,8)*PDL(7,1)
      DD=DENSU7(Z,DB16H,THO,THO,16.D0,ALPHA(9),T2,PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      ZSHT=PDM(6,8)
      ZMHO=PDM(5,8)
      ZSHO=SCALH(ZMHO,16.D0,THO)
      D(9)=DD*EXP(-ZSHT/ZSHO*(EXP(-(Z-ZMHO)/ZSHT)-1.D0))
      IF(MASS.NE.48) GO TO 90
C
C       TOTAL MASS DENSITY
C
      D(6) = 1.66D-24*(4.D0*D(1)+16.D0*D(2)+28.D0*D(3)+32.D0*D(4)+
     &       40.D0*D(5)+D(7)+14.D0*D(8))
      DB48=1.66D-24*(4.D0*DB04+16.D0*DB16+28.D0*DB28+32.D0*DB32+
     &        40.D0*DB40+DB01+14.D0*DB14)
      GO TO 90
C       TEMPERATURE AT ALTITUDE
   50 CONTINUE
      Z=ABS(ALT)
      DDUM  = DENSU7(Z,1.D0, TINF,TLB,0.D0,0.D0,T(2),PTM(6),S,MN1
     &,ZN1,TN1,TGN1)
   90 CONTINUE
C       ADJUST DENSITIES FROM CGS TO KGM
      IF(IMR.EQ.1) THEN
        DO 95 I=1,9
          D(I)=D(I)*1.D6
   95   CONTINUE
        D(6)=D(6)/1000.D0
      ENDIF
      ALAST=ALT
      RETURN
  100 FORMAT(1X,'MASS', I5, '  NOT VALID')
      END
C-----------------------------------------------------------------------
      SUBROUTINE METERS7(METER)
      IMPLICIT NONE
C      Convert outputs to Kg & Meters if METER true
      INTEGER*4 IMR
      LOGICAL METER
      COMMON/METSEL7/IMR
      SAVE
      IMR=0
      IF(METER) IMR=1
      END
C-----------------------------------------------------------------------
      FUNCTION SCALH(ALT,XM,TEMP)
      IMPLICIT NONE
C      Calculate scale height (km)
      REAL*8 ALT,XM,TEMP,GSURF,RE,RGAS,SCALH,G
      COMMON/PARMB/GSURF,RE
      SAVE
      DATA RGAS/831.4D0/
      G=GSURF/(1.D0+ALT/RE)**2
      SCALH=RGAS*TEMP/(G*XM)
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION GLOBE7(YRD,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)
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
      REAL*8 P44,P45,EXP1,EXP2,GLOBE7
c
      COMMON/TTEST/TINF,GB,ROUT,T
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/LPOLY/PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ DAY,DF,DFA,APD,APDF,APT,XLONG,CLONG,SLONG
      COMMON/LPOLYI/IYR
      SAVE
      DATA DGTR/1.74533D-2/,DR/1.72142D-2/,XL/1000.D0/,TLL/1000.D0/
      DATA SW9/1.D0/,DAYL/-1.D0/,P14/-1000.D0/,P18/-1000.D0/,
     &P32/-1000.D0/
      DATA HR/.2618D0/,SR/7.2722D-5/,SV/25*1.D0/,NSW/14/,P39/-1000.D0/
C       3hr Magnetic activity functions
C      Eq. A24d
      G0(A)=(A-4.D0+(P(26)-1.D0)*(A-4.D0+(EXP(-ABS(P(25))*
     &(A-4.D0))-1.D0)/ABS(P(25))))
C       Eq. A24c
       SUMEX(EX)=1.D0+(1.D0-EX**19)/(1.D0-EX)*EX**(.5D0)
C       Eq. A24a
      SG0(EX)=(G0(AP(2))+(G0(AP(3))*EX+G0(AP(4))*EX*EX+G0(AP(5))*EX**3
     $ +(G0(AP(6))*EX**4+G0(AP(7))*EX**12)*(1.D0-EX**8)/(1.D0-EX))
     $ )/SUMEX(EX)
      IF(ISW.NE.64999) CALL TSELEC7(SV)
      DO 10 J=1,14
       T(J)=0.D0
   10 CONTINUE
      IF(SW(9).GT.0.D0) SW9=1.D0
      IF(SW(9).LT.0.D0) SW9=-1.D0
      IYR = INT(YRD/1000.D0)
      DAY = YRD - IYR*1000.D0
      XLONG=LONG
C      Eq. A22 (remainder of code)
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
      T(1) =  P(20)*DF*(1.D0+P(60)*DFA) + P(21)*DF*DF + P(22)*DFA
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
      APDF = APD+(P45-1.D0)*(APD+(EXP(-P44  *APD)-1.D0)/P44)
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
      IF(P(25).LT.1.D-4) P(25)=1.D-4
      APT(1)=SG0(EXP1)
C      APT(2)=SG2(EXP1)
c      APT(3)=SG0(EXP2)
C      APT(4)=SG2(EXP2)
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
     $     COS(DGTR*LONG)
     $ +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $ +P(107)*PLG(2,2)+P(108)*PLG(4,2)+P(109)*PLG(6,2)
     $ +SWC(5)*(P(113)*PLG(2,2)+P(114)*PLG(4,2)+P(115)*PLG(6,2))*CD14)*
     $  SIN(DGTR*LONG))
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
C  PARMS NOT USED: 83, 90,100,140-150
   49 CONTINUE
      TINF=P(31)
      DO 50 I = 1,NSW
   50 TINF = TINF + ABS(SW(I))*T(I)
      GLOBE7 = TINF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TSELEC7(SV)
      IMPLICIT NONE
C        SET SWITCHES
C        Output in  COMMON/CSW/SW(25),ISW,SWC(25)
C        SW FOR MAIN TERMS, SWC FOR CROSS TERMS
C
C        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC7(SV),
C        WHERE SV IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
C        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C
C        To get current values of SW: CALL TRETRV7(SW)
C
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
      ENTRY TRETRV7(SVV)
      DO 200 I=1,25
        SVV(I)=SAV(I)
  200 CONTINUE
      END
C-----------------------------------------------------------------------
      FUNCTION GLOB7S(P)
      IMPLICIT NONE
C      VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
      INTEGER*4 ISW,J,I,IYR
      REAL*8 LONG,DR,DGTR,DAYL,P32,P18,P14,P39,PSET
      REAL*8 SW(25),SWC(25),P(100),T(14)
      REAL*8 PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC
      REAL*8 DAY,DF,DFA,APD,APDF,APT(4),CLONG,SLONG
      REAL*8 CD32,CD18,CD14,CD39,T71,T72,T81,T82,TT,GLOB7S
c
      COMMON/LPOLY/PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ DAY,DF,DFA,APD,APDF,APT,LONG,CLONG,SLONG
      COMMON/LPOLYI/IYR
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      SAVE
      DATA DR/1.72142D-2/,DGTR/1.74533D-2/,PSET/2.D0/
      DATA DAYL/-1.D0/,P32,P18,P14,P39/4*-1000.D0/
C       CONFIRM PARAMETER SET
      IF(P(100).EQ.0.D0) P(100)=PSET
      IF(P(100).NE.PSET) THEN
        WRITE(6,900) PSET,P(100)
  900   FORMAT(1X,'WRONG PARAMETER SET FOR GLOB7S',3F10.1)
        STOP
      ENDIF
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
     $     +P(27)*PLG(2,1)+P(15)*PLG(4,1)+P(60)*PLG(6,1)
C       SYMMETRICAL ANNUAL
      T(3)=(P(19)+P(48)*PLG(3,1)+P(30)*PLG(5,1))*CD32
C       SYMMETRICAL SEMIANNUAL
      T(4)=(P(16)+P(17)*PLG(3,1)+P(31)*PLG(5,1))*CD18
C       ASYMMETRICAL ANNUAL
      T(5)=(P(10)*PLG(2,1)+P(11)*PLG(4,1)+P(21)*PLG(6,1))*CD14
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
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5)
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
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
     $ T(9)=(P(51)*APT(1)+P(97)*PLG(3,1)*APT(1)*SWC(2))
   40 CONTINUE
      IF(SW(10).EQ.0.D0.OR.SW(11).EQ.0.D0.OR.LONG.LE.-1000.D0) GO TO 49
C        LONGITUDINAL
      T(11)= (1.D0+PLG(2,1)*(P(81)*SWC(5)*COS(DR*(DAY-P(82)))
     $           +P(86)*SWC(6)*COS(2.D0*DR*(DAY-P(87))))
     $        +P(84)*SWC(3)*COS(DR*(DAY-P(85)))
     $           +P(88)*SWC(4)*COS(2.D0*DR*(DAY-P(89))))
     $ *((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $   +P(75)*PLG(2,2)+P(76)*PLG(4,2)+P(77)*PLG(6,2)
     $    )*COS(DGTR*LONG)
     $  +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $   +P(78)*PLG(2,2)+P(79)*PLG(4,2)+P(80)*PLG(6,2)
     $    )*SIN(DGTR*LONG))
   49 CONTINUE
      TT=0.D0
      DO 50 I=1,14
   50 TT=TT+ABS(SW(I))*T(I)
      GLOB7S=TT
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSU7(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,
     $  MN1,ZN1,TN1,TGN1)
      IMPLICIT NONE
C       Calculate Temperature and Density Profiles for MSIS models
C       New lower thermo polynomial 10/30/89
      INTEGER*4 MN1,MP,II,JG,LT,IERR,IFUN,N,J
      INTEGER*4 MN,K
c
      REAL*8 ZN1(MN1),TN1(MN1),TGN1(2),XS(5),YS(5),Y2OUT(5)
      REAL*8 ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,GSURF,RE,RGAS
      REAL*8 QPB(50),DV(60),ZETA,DENSU7,ZA,Z,ZG2,TT,TA,DTA
      REAL*8 Z1,Z2,T1,T2,ZG,ZGDIF,YD1,YD2,X,Y,GLB,GAMM,YI,EXPL
      REAL*8 GAMMA,DENSA,ZZ,ZL
c
      COMMON/PARMB/GSURF,RE
      COMMON/LSQV/MP,II,JG,LT,QPB,IERR,IFUN,N,J,DV
      SAVE
      DATA RGAS/831.4D0/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
CCCCCCWRITE(6,*) 'DB',ALT,DLB,TINF,TLB,XM,ALPHA,ZLB,S2,MN1,ZN1,TN1
      DENSU7=1.D0
C        Joining altitude of Bates and spline
      ZA=ZN1(1)
      Z=DMAX1(ALT,ZA)
C      Geopotential altitude difference from ZLB
      ZG2=ZETA(Z,ZLB)
C      Bates temperature
      TT=TINF-(TINF-TLB)*EXP(-S2*ZG2)
      TA=TT
      TZ=TT
      DENSU7=TZ
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
      CALL SPLINE7(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT7(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1.D0/Y
      DENSU7=TZ
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
      DENSU7=DENSA
      IF(ALT.GE.ZA) GO TO 50
C
C      CALCULATE DENSITY BELOW ZA
      GLB=GSURF/(1.D0+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       integrate spline temperatures
      CALL SPLINI7(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.D0.OR.TZ.LE.0.D0) THEN
        EXPL=50.D0
      ENDIF
C       Density at altitude
      DENSU7=DENSU7*(T1/TZ)**(1.D0+ALPHA)*EXP(-EXPL)
   50 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSM7(ALT,D0,XM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
      IMPLICIT NONE
C       Calculate Temperature and Density Profiles for lower atmos.
      INTEGER*4 MN3,MN2,MP,II,JG,LT,IERR,IFUN,N,J
      INTEGER*4 MN,K
c
      REAL*8 ZN3(MN3),TN3(MN3),TGN3(2),XS(10),YS(10),Y2OUT(10)
      REAL*8 ZN2(MN2),TN2(MN2),TGN2(2)
      REAL*8 ALT,D0,XM,TZ,RGAS,GSURF,RE,TAF,QPB(50),DV(60)
      REAL*8 ZETA,DENSM7,Z,Z1,Z2,T1,Y2,ZG,ZGDIF,YD1,YD2,X,Y
      REAL*8 GLB,GAMM,YI,EXPL,T2,ZZ,ZL
c
      COMMON/PARMB/GSURF,RE
      COMMON/FIT/TAF
      COMMON/LSQV/MP,II,JG,LT,QPB,IERR,IFUN,N,J,DV
      SAVE
      DATA RGAS/831.4D0/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
      DENSM7=D0
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
      CALL SPLINE7(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT7(XS,YS,Y2OUT,MN,X,Y)
C       Temperature at altitude
      TZ=1.D0/Y
      IF(XM.EQ.0.D0) GO TO 20
C
C      CALCULATE STRATOSPHERE/MESOSPHERE DENSITY
      GLB=GSURF/(1.D0+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       Integrate temperature profile
      CALL SPLINI7(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.D0) EXPL=50.D0
C       Density at altitude
      DENSM7=DENSM7*(T1/TZ)*EXP(-EXPL)
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
      CALL SPLINE7(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT7(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1.D0/Y
      IF(XM.EQ.0.D0) GO TO 30
C
C      CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY
C
      GLB=GSURF/(1.D0+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C        Integrate temperature profile
      CALL SPLINI7(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.D0) EXPL=50.D0
C        Density at altitude
      DENSM7=DENSM7*(T1/TZ)*EXP(-EXPL)
   30 CONTINUE
   50 CONTINUE
      IF(XM.EQ.0.D0) DENSM7=TZ
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINE7(X,Y,N,YP1,YPN,Y2)
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
      SUBROUTINE SPLINT7(XA,YA,Y2A,N,X,Y)
      IMPLICIT NONE
C        CALCULATE CUBIC SPLINE INTERP VALUE
C        ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
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
      IF(H.EQ.0.D0) WRITE(6,*) 'BAD XA INPUT TO SPLINT7'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     $  ((A*A*A-A)*Y2A(KLO)+(B*B*B-B)*Y2A(KHI))*H*H/6.D0
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINI7(XA,YA,Y2A,N,X,YI)
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
      FUNCTION DNET7(DD,DM,ZHM,XMM,XM)
      IMPLICIT NONE
C       TURBOPAUSE CORRECTION FOR MSIS MODELS
C         Root mean density
C       8/20/80
C          DD - diffusive density
C          DM - full mixed density
C          ZHM - transition scale length
C          XMM - full mixed molecular weight
C          XM  - species molecular weight
C          DNET7 - combined density
      REAL*8 A,ZHM,XMM,XM,YLOG
      REAL*8 DM,DD,DNET7
c
      SAVE
      A=ZHM/(XMM-XM)
      IF(DM.GT.0.D0.AND.DD.GT.0.D0) GOTO 5
        WRITE(6,*) 'DNET7 LOG ERROR',DM,DD,XM
        IF(DD.EQ.0.D0.AND.DM.EQ.0.D0) DD=1.D0
        IF(DM.EQ.0.D0) GOTO 10
        IF(DD.EQ.0.D0) GOTO 20
    5 CONTINUE
      YLOG=A*DLOG(DM/DD)
      IF(YLOG.LT.-10.D0) GO TO 10
      IF(YLOG.GT.10.D0)  GO TO 20
        DNET7=DD*(1.D0+EXP(YLOG))**(1.D0/A)
        GO TO 50
   10 CONTINUE
        DNET7=DD
        GO TO 50
   20 CONTINUE
        DNET7=DM
        GO TO 50
   50 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR7(ALT, R,H1,ZH)
      IMPLICIT NONE
C        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C        ALT - altitude
C        R - target ratio
C        H1 - transition scale length
C        ZH - altitude of 1/2 R
      REAL*8 E,ALT,ZH,H1,EX,CCOR7,R
c
      SAVE
      E=(ALT-ZH)/H1
      IF(E.GT.70.D0) GO TO 20
      IF(E.LT.-70.D0) GO TO 10
        EX=EXP(E)
        CCOR7=R/(1.D0+EX)
        GO TO 50
   10   CCOR7=R
        GO TO 50
   20   CCOR7=0.D0
        GO TO 50
   50 CONTINUE
      CCOR7=EXP(CCOR7)
       RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR2(ALT, R,H1,ZH,H2)
      IMPLICIT NONE
C       O&O2 CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
      REAL*8 E1,E2,ALT,ZH,H1,H2,EX1,EX2,CCOR2,R
c
      E1=(ALT-ZH)/H1
      E2=(ALT-ZH)/H2
      IF(E1.GT.70.D0 .OR. E2.GT.70.D0) GO TO 20
      IF(E1.LT.-70.D0 .AND. E2.LT.-70.D0) GO TO 10
        EX1=EXP(E1)
        EX2=EXP(E2)
        CCOR2=R/(1.D0+.5D0*(EX1+EX2))
        GO TO 50
   10   CCOR2=R
        GO TO 50
   20   CCOR2=0.D0
        GO TO 50
   50 CONTINUE
      CCOR2=EXP(CCOR2)
       RETURN
      END
C-----------------------------------------------------------------------
      BLOCK DATA GTD7BK
      IMPLICIT NONE
C          MSISE-00 01-FEB-02
      INTEGER*4 IMR
c
      REAL*8 PT1(50),PT2(50),PT3(50),PA1(50),PA2(50),PA3(50)
      REAL*8 PB1(50),PB2(50),PB3(50),PC1(50),PC2(50),PC3(50)
      REAL*8 PD1(50),PD2(50),PD3(50),PE1(50),PE2(50),PE3(50)
      REAL*8 PF1(50),PF2(50),PF3(50),PG1(50),PG2(50),PG3(50)
      REAL*8 PH1(50),PH2(50),PH3(50),PI1(50),PI2(50),PI3(50)
      REAL*8 PJ1(50),PJ2(50),PJ3(50),PK1(50),PL1(50),PL2(50)
      REAL*8 PM1(50),PM2(50),PN1(50),PN2(50),PO1(50),PO2(50)
      REAL*8 PP1(50),PP2(50),PQ1(50),PQ2(50),PR1(50),PR2(50)
      REAL*8 PS1(50),PS2(50),PU1(50),PU2(50),PV1(50),PV2(50)
      REAL*8 PW1(50),PW2(50),PX1(50),PX2(50),PY1(50),PY2(50)
      REAL*8 PZ1(50),PZ2(50),PAA1(50),PAA2(50)
      REAL*8  PAVGM(10),PTM(10),PDM(10,8)
c
      CHARACTER NAME(2)*4,ISDATE(3)*4,ISTIME(2)*4
c
      COMMON/PARM7/PT1,PT2,PT3,PA1,PA2,PA3,
     $ PB1,PB2,PB3,PC1,PC2,PC3,
     $ PD1,PD2,PD3,PE1,PE2,PE3,
     $ PF1,PF2,PF3,PG1,PG2,PG3,
     $ PH1,PH2,PH3,PI1,PI2,PI3,
     $ PJ1,PJ2,PJ3,PK1,PL1,PL2,
     $ PM1,PM2,PN1,PN2,PO1,PO2,
     $ PP1,PP2,PQ1,PQ2,PR1,PR2,
     $ PS1,PS2,PU1,PU2,PV1,PV2,
     $ PW1,PW2,PX1,PX2,PY1,PY2,
     $ PZ1,PZ2,PAA1,PAA2
      COMMON/LOWER7/PTM,PDM
      COMMON/MAVG7/PAVGM
      COMMON/DATIM7/ISDATE,ISTIME,NAME
      COMMON/METSEL7/IMR
      DATA IMR/0/
      DATA ISDATE/'01-F','EB-0','2   '/,ISTIME/'15:4','9:27'/
      DATA NAME/'MSIS','E-00'/
C         TEMPERATURE
      DATA PT1/
     *  9.86573D-01, 1.62228D-02, 1.55270D-02,-1.04323D-01,-3.75801D-03,
     * -1.18538D-03,-1.24043D-01, 4.56820D-03, 8.76018D-03,-1.36235D-01,
     * -3.52427D-02, 8.84181D-03,-5.92127D-03,-8.61650D+00, 0.00000D+00,
     *  1.28492D-02, 0.00000D+00, 1.30096D+02, 1.04567D-02, 1.65686D-03,
     * -5.53887D-06, 2.97810D-03, 0.00000D+00, 5.13122D-03, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,-7.27026D-06,
     *  0.00000D+00, 6.74494D+00, 4.93933D-03, 2.21656D-03, 2.50802D-03,
     *  0.00000D+00, 0.00000D+00,-2.08841D-02,-1.79873D+00, 1.45103D-03,
     *  2.81769D-04,-1.44703D-03,-5.16394D-05, 8.47001D-02, 1.70147D-01,
     *  5.72562D-03, 5.07493D-05, 4.36148D-03, 1.17863D-04, 4.74364D-03/
      DATA PT2/
     *  6.61278D-03, 4.34292D-05, 1.44373D-03, 2.41470D-05, 2.84426D-03,
     *  8.56560D-04, 2.04028D-03, 0.00000D+00,-3.15994D+03,-2.46423D-03,
     *  1.13843D-03, 4.20512D-04, 0.00000D+00,-9.77214D+01, 6.77794D-03,
     *  5.27499D-03, 1.14936D-03, 0.00000D+00,-6.61311D-03,-1.84255D-02,
     * -1.96259D-02, 2.98618D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  6.44574D+02, 8.84668D-04, 5.05066D-04, 0.00000D+00, 4.02881D+03,
     * -1.89503D-03, 0.00000D+00, 0.00000D+00, 8.21407D-04, 2.06780D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -1.20410D-02,-3.63963D-03, 9.92070D-05,-1.15284D-04,-6.33059D-05,
     * -6.05545D-01, 8.34218D-03,-9.13036D+01, 3.71042D-04, 0.00000D+00/
      DATA PT3/
     *  4.19000D-04, 2.70928D-03, 3.31507D-03,-4.44508D-03,-4.96334D-03,
     * -1.60449D-03, 3.95119D-03, 2.48924D-03, 5.09815D-04, 4.05302D-03,
     *  2.24076D-03, 0.00000D+00, 6.84256D-03, 4.66354D-04, 0.00000D+00,
     * -3.68328D-04, 0.00000D+00, 0.00000D+00,-1.46870D+02, 0.00000D+00,
     *  0.00000D+00, 1.09501D-03, 4.65156D-04, 5.62583D-04, 3.21596D+00,
     *  6.43168D-04, 3.14860D-03, 3.40738D-03, 1.78481D-03, 9.62532D-04,
     *  5.58171D-04, 3.43731D+00,-2.33195D-01, 5.10289D-04, 0.00000D+00,
     *  0.00000D+00,-9.25347D+04, 0.00000D+00,-1.99639D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         HD DENSITY
      DATA PA1/
     *  1.09979D+00,-4.88060D-02,-1.97501D-01,-9.10280D-02,-6.96558D-03,
     *  2.42136D-02, 3.91333D-01,-7.20068D-03,-3.22718D-02, 1.41508D+00,
     *  1.68194D-01, 1.85282D-02, 1.09384D-01,-7.24282D+00, 0.00000D+00,
     *  2.96377D-01,-4.97210D-02, 1.04114D+02,-8.61108D-02,-7.29177D-04,
     *  1.48998D-06, 1.08629D-03, 0.00000D+00, 0.00000D+00, 8.31090D-02,
     *  1.12818D-01,-5.75005D-02,-1.29919D-02,-1.78849D-02,-2.86343D-06,
     *  0.00000D+00,-1.51187D+02,-6.65902D-03, 0.00000D+00,-2.02069D-03,
     *  0.00000D+00, 0.00000D+00, 4.32264D-02,-2.80444D+01,-3.26789D-03,
     *  2.47461D-03, 0.00000D+00, 0.00000D+00, 9.82100D-02, 1.22714D-01,
     * -3.96450D-02, 0.00000D+00,-2.76489D-03, 0.00000D+00, 1.87723D-03/
      DATA PA2/
     * -8.09813D-03, 4.34428D-05,-7.70932D-03, 0.00000D+00,-2.28894D-03,
     * -5.69070D-03,-5.22193D-03, 6.00692D-03,-7.80434D+03,-3.48336D-03,
     * -6.38362D-03,-1.82190D-03, 0.00000D+00,-7.58976D+01,-2.17875D-02,
     * -1.72524D-02,-9.06287D-03, 0.00000D+00, 2.44725D-02, 8.66040D-02,
     *  1.05712D-01, 3.02543D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -6.01364D+03,-5.64668D-03,-2.54157D-03, 0.00000D+00, 3.15611D+02,
     * -5.69158D-03, 0.00000D+00, 0.00000D+00,-4.47216D-03,-4.49523D-03,
     *  4.64428D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  4.51236D-02, 2.46520D-02, 6.17794D-03, 0.00000D+00, 0.00000D+00,
     * -3.62944D-01,-4.80022D-02,-7.57230D+01,-1.99656D-03, 0.00000D+00/
      DATA PA3/
     * -5.18780D-03,-1.73990D-02,-9.03485D-03, 7.48465D-03, 1.53267D-02,
     *  1.06296D-02, 1.18655D-02, 2.55569D-03, 1.69020D-03, 3.51936D-02,
     * -1.81242D-02, 0.00000D+00,-1.00529D-01,-5.10574D-03, 0.00000D+00,
     *  2.10228D-03, 0.00000D+00, 0.00000D+00,-1.73255D+02, 5.07833D-01,
     * -2.41408D-01, 8.75414D-03, 2.77527D-03,-8.90353D-05,-5.25148D+00,
     * -5.83899D-03,-2.09122D-02,-9.63530D-03, 9.77164D-03, 4.07051D-03,
     *  2.53555D-04,-5.52875D+00,-3.55993D-01,-2.49231D-03, 0.00000D+00,
     *  0.00000D+00, 2.86026D+01, 0.00000D+00, 3.42722D-04, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         O DENSITY
      DATA PB1/
     *  1.02315D+00,-1.59710D-01,-1.06630D-01,-1.77074D-02,-4.42726D-03,
     *  3.44803D-02, 4.45613D-02,-3.33751D-02,-5.73598D-02, 3.50360D-01,
     *  6.33053D-02, 2.16221D-02, 5.42577D-02,-5.74193D+00, 0.00000D+00,
     *  1.90891D-01,-1.39194D-02, 1.01102D+02, 8.16363D-02, 1.33717D-04,
     *  6.54403D-06, 3.10295D-03, 0.00000D+00, 0.00000D+00, 5.38205D-02,
     *  1.23910D-01,-1.39831D-02, 0.00000D+00, 0.00000D+00,-3.95915D-06,
     *  0.00000D+00,-7.14651D-01,-5.01027D-03, 0.00000D+00,-3.24756D-03,
     *  0.00000D+00, 0.00000D+00, 4.42173D-02,-1.31598D+01,-3.15626D-03,
     *  1.24574D-03,-1.47626D-03,-1.55461D-03, 6.40682D-02, 1.34898D-01,
     * -2.42415D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 6.13666D-04/
      DATA PB2/
     * -5.40373D-03, 2.61635D-05,-3.33012D-03, 0.00000D+00,-3.08101D-03,
     * -2.42679D-03,-3.36086D-03, 0.00000D+00,-1.18979D+03,-5.04738D-02,
     * -2.61547D-03,-1.03132D-03, 1.91583D-04,-8.38132D+01,-1.40517D-02,
     * -1.14167D-02,-4.08012D-03, 1.73522D-04,-1.39644D-02,-6.64128D-02,
     * -6.85152D-02,-1.34414D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  6.07916D+02,-4.12220D-03,-2.20996D-03, 0.00000D+00, 1.70277D+03,
     * -4.63015D-03, 0.00000D+00, 0.00000D+00,-2.25360D-03,-2.96204D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  3.92786D-02, 1.31186D-02,-1.78086D-03, 0.00000D+00, 0.00000D+00,
     * -3.90083D-01,-2.84741D-02,-7.78400D+01,-1.02601D-03, 0.00000D+00/
      DATA PB3/
     * -7.26485D-04,-5.42181D-03,-5.59305D-03, 1.22825D-02, 1.23868D-02,
     *  6.68835D-03,-1.03303D-02,-9.51903D-03, 2.70021D-04,-2.57084D-02,
     * -1.32430D-02, 0.00000D+00,-3.81000D-02,-3.16810D-03, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-9.05762D-04,-2.14590D-03,-1.17824D-03, 3.66732D+00,
     * -3.79729D-04,-6.13966D-03,-5.09082D-03,-1.96332D-03,-3.08280D-03,
     * -9.75222D-04, 4.03315D+00,-2.52710D-01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C         N2 DENSITY
      DATA PC1/
     *  1.16112D+00, 0.00000D+00, 0.00000D+00, 3.33725D-02, 0.00000D+00,
     *  3.48637D-02,-5.44368D-03, 0.00000D+00,-6.73940D-02, 1.74754D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 1.74712D+02, 0.00000D+00,
     *  1.26733D-01, 0.00000D+00, 1.03154D+02, 5.52075D-02, 0.00000D+00,
     *  0.00000D+00, 8.13525D-04, 0.00000D+00, 0.00000D+00, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-2.50482D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.48894D-03,
     *  6.16053D-04,-5.79716D-04, 2.95482D-03, 8.47001D-02, 1.70147D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PC2/
     *  0.00000D+00, 2.47425D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
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
     *  9.44846D-01, 0.00000D+00, 0.00000D+00,-3.08617D-02, 0.00000D+00,
     * -2.44019D-02, 6.48607D-03, 0.00000D+00, 3.08181D-02, 4.59392D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 1.74712D+02, 0.00000D+00,
     *  2.13260D-02, 0.00000D+00,-3.56958D+02, 0.00000D+00, 1.82278D-04,
     *  0.00000D+00, 3.07472D-04, 0.00000D+00, 0.00000D+00, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 3.83054D-03, 0.00000D+00, 0.00000D+00,
     * -1.93065D-03,-1.45090D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.23493D-03, 1.36736D-03, 8.47001D-02, 1.70147D-01,
     *  3.71469D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PD2/
     *  5.10250D-03, 2.47425D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 3.68756D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00/
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
     *  1.35580D+00, 1.44816D-01, 0.00000D+00, 6.07767D-02, 0.00000D+00,
     *  2.94777D-02, 7.46900D-02, 0.00000D+00,-9.23822D-02, 8.57342D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 2.38636D+01, 0.00000D+00,
     *  7.71653D-02, 0.00000D+00, 8.18751D+01, 1.87736D-02, 0.00000D+00,
     *  0.00000D+00, 1.49667D-02, 0.00000D+00, 0.00000D+00, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-3.67874D+02, 5.48158D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 8.47001D-02, 1.70147D-01,
     *  1.22631D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PE2/
     *  8.17187D-03, 3.71617D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.10826D-03,
     * -3.13640D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -7.35742D-02,-5.00266D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 1.94965D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00/
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
     *  1.04761D+00, 2.00165D-01, 2.37697D-01, 3.68552D-02, 0.00000D+00,
     *  3.57202D-02,-2.14075D-01, 0.00000D+00,-1.08018D-01,-3.73981D-01,
     *  0.00000D+00, 3.10022D-02,-1.16305D-03,-2.07596D+01, 0.00000D+00,
     *  8.64502D-02, 0.00000D+00, 9.74908D+01, 5.16707D-02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 3.46193D+02, 1.34297D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-3.48509D-03,
     * -1.54689D-04, 0.00000D+00, 0.00000D+00, 8.47001D-02, 1.70147D-01,
     *  1.47753D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PF2/
     *  1.89320D-02, 3.68181D-05, 1.32570D-02, 0.00000D+00, 0.00000D+00,
     *  3.59719D-03, 7.44328D-03,-1.00023D-03,-6.50528D+03, 0.00000D+00,
     *  1.03485D-02,-1.00983D-03,-4.06916D-03,-6.60864D+01,-1.71533D-02,
     *  1.10605D-02, 1.20300D-02,-5.20034D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -2.62769D+03, 7.13755D-03, 4.17999D-03, 0.00000D+00, 1.25910D+04,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-2.23595D-03, 4.60217D-03,
     *  5.71794D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -3.18353D-02,-2.35526D-02,-1.36189D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.03522D-02,-6.67837D+01,-1.09724D-03, 0.00000D+00/
      DATA PF3/
     * -1.38821D-02, 1.60468D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.51574D-02,
     * -5.44470D-04, 0.00000D+00, 7.28224D-02, 6.59413D-02, 0.00000D+00,
     * -5.15692D-03, 0.00000D+00, 0.00000D+00,-3.70367D+03, 0.00000D+00,
     *  0.00000D+00, 1.36131D-02, 5.38153D-03, 0.00000D+00, 4.76285D+00,
     * -1.75677D-02, 2.26301D-02, 0.00000D+00, 1.76631D-02, 4.77162D-03,
     *  0.00000D+00, 5.39354D+00, 0.00000D+00,-7.51710D-03, 0.00000D+00,
     *  0.00000D+00,-8.82736D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          H DENSITY
      DATA PG1/
     *  1.26376D+00,-2.14304D-01,-1.49984D-01, 2.30404D-01, 2.98237D-02,
     *  2.68673D-02, 2.96228D-01, 2.21900D-02,-2.07655D-02, 4.52506D-01,
     *  1.20105D-01, 3.24420D-02, 4.24816D-02,-9.14313D+00, 0.00000D+00,
     *  2.47178D-02,-2.88229D-02, 8.12805D+01, 5.10380D-02,-5.80611D-03,
     *  2.51236D-05,-1.24083D-02, 0.00000D+00, 0.00000D+00, 8.66784D-02,
     *  1.58727D-01,-3.48190D-02, 0.00000D+00, 0.00000D+00, 2.89885D-05,
     *  0.00000D+00, 1.53595D+02,-1.68604D-02, 0.00000D+00, 1.01015D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.84552D-04,
     * -1.22181D-03, 0.00000D+00, 0.00000D+00, 8.47001D-02, 1.70147D-01,
     * -1.04927D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,-5.91313D-03/
      DATA PG2/
     * -2.30501D-02, 3.14758D-05, 0.00000D+00, 0.00000D+00, 1.26956D-02,
     *  8.35489D-03, 3.10513D-04, 0.00000D+00, 3.42119D+03,-2.45017D-03,
     * -4.27154D-04, 5.45152D-04, 1.89896D-03, 2.89121D+01,-6.49973D-03,
     * -1.93855D-02,-1.48492D-02, 0.00000D+00,-5.10576D-02, 7.87306D-02,
     *  9.51981D-02,-1.49422D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  2.65503D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 6.37110D-03, 3.24789D-04,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  6.14274D-02, 1.00376D-02,-8.41083D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.27099D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PG3/
     * -3.94077D-03,-1.28601D-02,-7.97616D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-6.71465D-03,-1.69799D-03, 1.93772D-03, 3.81140D+00,
     * -7.79290D-03,-1.82589D-02,-1.25860D-02,-1.04311D-02,-3.02465D-03,
     *  2.43063D-03, 3.63237D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C          N DENSITY
      DATA PH1/
     *  7.09557D+01,-3.26740D-01, 0.00000D+00,-5.16829D-01,-1.71664D-03,
     *  9.09310D-02,-6.71500D-01,-1.47771D-01,-9.27471D-02,-2.30862D-01,
     * -1.56410D-01, 1.34455D-02,-1.19717D-01, 2.52151D+00, 0.00000D+00,
     * -2.41582D-01, 5.92939D-02, 4.39756D+00, 9.15280D-02, 4.41292D-03,
     *  0.00000D+00, 8.66807D-03, 0.00000D+00, 0.00000D+00, 8.66784D-02,
     *  1.58727D-01, 9.74701D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 6.70217D+01,-1.31660D-03, 0.00000D+00,-1.65317D-02,
     *  0.00000D+00, 0.00000D+00, 8.50247D-02, 2.77428D+01, 4.98658D-03,
     *  6.15115D-03, 9.50156D-03,-2.12723D-02, 8.47001D-02, 1.70147D-01,
     * -2.38645D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.37380D-03/
      DATA PH2/
     * -8.41918D-03, 2.80145D-05, 7.12383D-03, 0.00000D+00,-1.66209D-02,
     *  1.03533D-04,-1.68898D-02, 0.00000D+00, 3.64526D+03, 0.00000D+00,
     *  6.54077D-03, 3.69130D-04, 9.94419D-04, 8.42803D+01,-1.16124D-02,
     * -7.74414D-03,-1.68844D-03, 1.42809D-03,-1.92955D-03, 1.17225D-01,
     * -2.41512D-02, 1.50521D+04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  1.60261D+03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-3.54403D-04,-1.87270D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  2.76439D-02, 6.43207D-03,-3.54300D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-2.80221D-02, 8.11228D+01,-6.75255D-04, 0.00000D+00/
      DATA PH3/
     * -1.05162D-02,-3.48292D-03,-6.97321D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.45546D-03,-1.31970D-02,-3.57751D-03,-1.09021D+00,
     * -1.50181D-02,-7.12841D-03,-6.64590D-03,-3.52610D-03,-1.87773D-02,
     * -2.22432D-03,-3.93895D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
C        HOT O DENSITY
      DATA PI1/
     *  6.04050D-02, 1.57034D+00, 2.99387D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.51018D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-8.61650D+00, 1.26454D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 5.50878D-03, 0.00000D+00, 0.00000D+00, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 6.23881D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 8.47001D-02, 1.70147D-01,
     * -9.45934D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
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
     *  9.56827D-01, 6.20637D-02, 3.18433D-02, 0.00000D+00, 0.00000D+00,
     *  3.94900D-02, 0.00000D+00, 0.00000D+00,-9.24882D-03,-7.94023D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 1.74712D+02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.74677D-03, 0.00000D+00, 1.54951D-02, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-6.99007D-04, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 1.24362D-02,-5.28756D-03, 8.47001D-02, 1.70147D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PJ2/
     *  0.00000D+00, 2.47425D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
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
     *  1.09930D+00, 3.90631D+00, 3.07165D+00, 9.86161D-01, 1.63536D+01,
     *  4.63830D+00, 1.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 1.28840D+00, 3.10302D-02, 1.18339D-01,
     *  1.00000D+00, 7.00000D-01, 1.15020D+00, 3.44689D+00, 1.28840D+00,
     *  1.00000D+00, 1.08738D+00, 1.22947D+00, 1.10016D+00, 7.34129D-01,
     *  1.15241D+00, 2.22784D+00, 7.95046D-01, 4.01612D+00, 4.47749D+00,
     *  1.23435D+02,-7.60535D-02, 1.68986D-06, 7.44294D-01, 1.03604D+00,
     *  1.72783D+02, 1.15020D+00, 3.44689D+00,-7.46230D-01, 9.49154D-01/
C         LOWER BOUNDARY
      DATA PTM/
     L  1.04130D+03, 3.86000D+02, 1.95000D+02, 1.66728D+01, 2.13000D+02,
     L  1.20000D+02, 2.40000D+02, 1.87000D+02,-2.00000D+00, 0.00000D+00/
      DATA PDM/
     L  2.45600D+07, 6.71072D-06, 1.00000D+02, 0.00000D+00, 1.10000D+02,
     L  1.00000D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  8.59400D+10, 1.00000D+00, 1.05000D+02,-8.00000D+00, 1.10000D+02,
     L  1.00000D+01, 9.00000D+01, 2.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  2.81000D+11, 0.00000D+00, 1.05000D+02, 2.80000D+01, 2.89500D+01,
     L  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
C
     L  3.30000D+10, 2.68270D-01, 1.05000D+02, 1.00000D+00, 1.10000D+02,
     L  1.00000D+01, 1.10000D+02,-1.00000D+01, 0.00000D+00, 0.00000D+00,
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
     L  1.00000D+06, 1.00000D+00, 1.05000D+02,-8.00000D+00, 5.50000D+02,
     L  7.60000D+01, 9.00000D+01, 2.00000D+00, 0.00000D+00, 4.00000D+03/
C         TN1(2)
      DATA PL1/
     *  1.00858D+00, 4.56011D-02,-2.22972D-02,-5.44388D-02, 5.23136D-04,
     * -1.88849D-02, 5.23707D-02,-9.43646D-03, 6.31707D-03,-7.80460D-02,
     * -4.88430D-02, 0.00000D+00, 0.00000D+00,-7.60250D+00, 0.00000D+00,
     * -1.44635D-02,-1.76843D-02,-1.21517D+02, 2.85647D-02, 0.00000D+00,
     *  0.00000D+00, 6.31792D-04, 0.00000D+00, 5.77197D-03, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-8.90272D+03, 3.30611D-03, 3.02172D-03, 0.00000D+00,
     * -2.13673D-03,-3.20910D-04, 0.00000D+00, 0.00000D+00, 2.76034D-03,
     *  2.82487D-03,-2.97592D-04,-4.21534D-03, 8.47001D-02, 1.70147D-01,
     *  8.96456D-03, 0.00000D+00,-1.08596D-02, 0.00000D+00, 0.00000D+00/
      DATA PL2/
     *  5.57917D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 9.65405D-03, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C         TN1(3)
      DATA PM1/
     *  9.39664D-01, 8.56514D-02,-6.79989D-03, 2.65929D-02,-4.74283D-03,
     *  1.21855D-02,-2.14905D-02, 6.49651D-03,-2.05477D-02,-4.24952D-02,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 1.19148D+01, 0.00000D+00,
     *  1.18777D-02,-7.28230D-02,-8.15965D+01, 1.73887D-02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.44691D-02, 2.80259D-04, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.16584D+02, 3.18713D-03, 7.37479D-03, 0.00000D+00,
     * -2.55018D-03,-3.92806D-03, 0.00000D+00, 0.00000D+00,-2.89757D-03,
     * -1.33549D-03, 1.02661D-03, 3.53775D-04, 8.47001D-02, 1.70147D-01,
     * -9.17497D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PM2/
     *  3.56082D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.00902D-02, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C         TN1(4)
      DATA PN1/
     *  9.85982D-01,-4.55435D-02, 1.21106D-02, 2.04127D-02,-2.40836D-03,
     *  1.11383D-02,-4.51926D-02, 1.35074D-02,-6.54139D-03, 1.15275D-01,
     *  1.28247D-01, 0.00000D+00, 0.00000D+00,-5.30705D+00, 0.00000D+00,
     * -3.79332D-02,-6.24741D-02, 7.71062D-01, 2.96315D-02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 6.81051D-03,-4.34767D-03, 8.66784D-02,
     *  1.58727D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 1.07003D+01,-2.76907D-03, 4.32474D-04, 0.00000D+00,
     *  1.31497D-03,-6.47517D-04, 0.00000D+00,-2.20621D+01,-1.10804D-03,
     * -8.09338D-04, 4.18184D-04, 4.29650D-03, 8.47001D-02, 1.70147D-01,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00/
      DATA PN2/
     * -4.04337D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-9.52550D-04,
     *  8.56253D-04, 4.33114D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 1.21223D-03,
     *  2.38694D-04, 9.15245D-04, 1.28385D-03, 8.67668D-04,-5.61425D-06,
     *  1.04445D+00, 3.41112D+01, 0.00000D+00,-8.40704D-01,-2.39639D+02,
     *  7.06668D-01,-2.05873D+01,-3.63696D-01, 2.39245D+01, 0.00000D+00,
     * -1.06657D-03,-7.67292D-04, 1.54534D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C         TN1(5) TN2(1)
      DATA PO1/
     *  1.00320D+00, 3.83501D-02,-2.38983D-03, 2.83950D-03, 4.20956D-03,
     *  5.86619D-04, 2.19054D-02,-1.00946D-02,-3.50259D-03, 4.17392D-02,
     * -8.44404D-03, 0.00000D+00, 0.00000D+00, 4.96949D+00, 0.00000D+00,
     * -7.06478D-03,-1.46494D-02, 3.13258D+01,-1.86493D-03, 0.00000D+00,
     * -1.67499D-02, 0.00000D+00, 0.00000D+00, 5.12686D-04, 8.66784D-02,
     *  1.58727D-01,-4.64167D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  4.37353D-03,-1.99069D+02, 0.00000D+00,-5.34884D-03, 0.00000D+00,
     *  1.62458D-03, 2.93016D-03, 2.67926D-03, 5.90449D+02, 0.00000D+00,
     *  0.00000D+00,-1.17266D-03,-3.58890D-04, 8.47001D-02, 1.70147D-01,
     *  0.00000D+00, 0.00000D+00, 1.38673D-02, 0.00000D+00, 0.00000D+00/
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
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TN2(2)
      DATA PP1/
     *  9.81637D-01,-1.41317D-03, 3.87323D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-3.58707D-02,
     * -8.63658D-03, 0.00000D+00, 0.00000D+00,-2.02226D+00, 0.00000D+00,
     * -8.69424D-03,-1.91397D-02, 8.76779D+01, 4.52188D-03, 0.00000D+00,
     *  2.23760D-02, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-7.07572D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     * -4.11210D-03, 3.50060D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-8.36657D-03, 1.61347D+01, 0.00000D+00,
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
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TN2(3)
      DATA PQ1/
     *  1.00422D+00,-7.11212D-03, 5.24480D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-5.28914D-02,
     * -2.41301D-02, 0.00000D+00, 0.00000D+00,-2.12219D+01,-1.03830D-02,
     * -3.28077D-03, 1.65727D-02, 1.68564D+00,-6.68154D-03, 0.00000D+00,
     *  1.45155D-02, 0.00000D+00, 8.42365D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-4.34645D-03, 0.00000D+00, 0.00000D+00, 2.16780D-02,
     *  0.00000D+00,-1.38459D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 7.04573D-03,-4.73204D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 1.08767D-02, 0.00000D+00, 0.00000D+00/
      DATA PQ2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-8.08279D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 5.21769D-04,
     * -2.27387D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.26769D-03,
     *  3.16901D-03, 4.60316D-04,-1.01431D-04, 1.02131D-03, 9.96601D-04,
     *  1.25707D+00, 2.50114D+01, 0.00000D+00, 4.24472D-01,-2.77655D+01,
     *  3.44625D-01, 2.75412D+01, 0.00000D+00, 7.94251D+02, 0.00000D+00,
     *  2.45835D-03, 1.38871D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TN2(4) TN3(1)
      DATA PR1/
     *  1.01890D+00,-2.46603D-02, 1.00078D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-6.70977D-02,
     * -4.02286D-02, 0.00000D+00, 0.00000D+00,-2.29466D+01,-7.47019D-03,
     *  2.26580D-03, 2.63931D-02, 3.72625D+01,-6.39041D-03, 0.00000D+00,
     *  9.58383D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.85291D-03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 1.39717D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 9.19771D-03,-3.69121D+02, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-1.57067D-02, 0.00000D+00, 0.00000D+00/
      DATA PR2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-7.07265D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.92953D-03,
     * -2.77739D-03,-4.40092D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.47280D-03,
     *  2.95035D-04,-1.81246D-03, 2.81945D-03, 4.27296D-03, 9.78863D-04,
     *  1.40545D+00,-6.19173D+00, 0.00000D+00, 0.00000D+00,-7.93632D+01,
     *  4.44643D-01,-4.03085D+02, 0.00000D+00, 1.15603D+01, 0.00000D+00,
     *  2.25068D-03, 8.48557D-04,-2.98493D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TN3(2)
      DATA PS1/
     *  9.75801D-01, 3.80680D-02,-3.05198D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.85575D-02,
     *  5.04057D-02, 0.00000D+00, 0.00000D+00,-1.76046D+02, 1.44594D-02,
     * -1.48297D-03,-3.68560D-03, 3.02185D+01,-3.23338D-03, 0.00000D+00,
     *  1.53569D-02, 0.00000D+00,-1.15558D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 4.89620D-03, 0.00000D+00, 0.00000D+00,-1.00616D-02,
     * -8.21324D-03,-1.57757D+02, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 6.63564D-03, 4.58410D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-2.51280D-02, 0.00000D+00, 0.00000D+00/
      DATA PS2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 9.91215D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-8.73148D-04,
     * -1.29648D-03,-7.32026D-05, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-4.68110D-03,
     * -4.66003D-03,-1.31567D-03,-7.39390D-04, 6.32499D-04,-4.65588D-04,
     * -1.29785D+00,-1.57139D+02, 0.00000D+00, 2.58350D-01,-3.69453D+01,
     *  4.10672D-01, 9.78196D+00,-1.52064D-01,-3.85084D+03, 0.00000D+00,
     * -8.52706D-04,-1.40945D-03,-7.26786D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TN3(3)
      DATA PU1/
     *  9.60722D-01, 7.03757D-02,-3.00266D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.22671D-02,
     *  4.10423D-02, 0.00000D+00, 0.00000D+00,-1.63070D+02, 1.06073D-02,
     *  5.40747D-04, 7.79481D-03, 1.44908D+02, 1.51484D-04, 0.00000D+00,
     *  1.97547D-02, 0.00000D+00,-1.41844D-02, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 5.77884D-03, 0.00000D+00, 0.00000D+00, 9.74319D-03,
     *  0.00000D+00,-2.88015D+03, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00,-4.44902D-03,-2.92760D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 2.34419D-02, 0.00000D+00, 0.00000D+00/
      DATA PU2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 5.36685D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-4.65325D-04,
     * -5.50628D-04, 3.31465D-04, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.06179D-03,
     * -3.08575D-03,-7.93589D-04,-1.08629D-04, 5.95511D-04,-9.05050D-04,
     *  1.18997D+00, 4.15924D+01, 0.00000D+00,-4.72064D-01,-9.47150D+02,
     *  3.98723D-01, 1.98304D+01, 0.00000D+00, 3.73219D+03, 0.00000D+00,
     * -1.50040D-03,-1.14933D-03,-1.56769D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TN3(4)
      DATA PV1/
     *  1.03123D+00,-7.05124D-02, 8.71615D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-3.82621D-02,
     * -9.80975D-03, 0.00000D+00, 0.00000D+00, 2.89286D+01, 9.57341D-03,
     *  0.00000D+00, 0.00000D+00, 8.66153D+01, 7.91938D-04, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 4.68917D-03, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 7.86638D-03, 0.00000D+00, 0.00000D+00, 9.90827D-03,
     *  0.00000D+00, 6.55573D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-4.00200D+01, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 7.07457D-03, 0.00000D+00, 0.00000D+00/
      DATA PV2/
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 5.72268D-03,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.04970D-04,
     *  1.21560D-03,-8.05579D-06, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-2.49941D-03,
     * -4.57256D-04,-1.59311D-04, 2.96481D-04,-1.77318D-03,-6.37918D-04,
     *  1.02395D+00, 1.28172D+01, 0.00000D+00, 1.49903D-01,-2.63818D+01,
     *  0.00000D+00, 4.70628D+01,-2.22139D-01, 4.82292D-02, 0.00000D+00,
     * -8.67075D-04,-5.86479D-04, 5.32462D-04, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TN3(5) SURFACE TEMP TSL
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
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TGN3(2) SURFACE GRAD TSLG
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
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TGN2(1) TGN1(2)
      DATA PY1/
     *  8.60028D-01, 3.77052D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.17570D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 7.77757D-03, 0.00000D+00,
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
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          TGN3(1) TGN2(2)
      DATA PZ1/
     *  1.06029D+00,-5.25231D-02, 3.73034D-01, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 3.31072D-02,
     * -3.88409D-01, 0.00000D+00, 0.00000D+00,-1.65295D+02,-2.13801D-01,
     * -4.38916D-02,-3.22716D-01,-8.82393D+01, 1.18458D-01, 0.00000D+00,
     * -4.35863D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00,-1.19782D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 2.62229D+01, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     *  0.00000D+00, 0.00000D+00, 0.00000D+00,-5.37443D+01, 0.00000D+00,
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
     *  0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D+00/
C          SEMIANNUAL MULT SAM
      DATA PAA1/
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00,
     *  1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00/
      DATA PAA2/
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
C         MIDDLE ATMOSPHERE AVERAGES
      DATA PAVGM/
     M  2.61000D+02, 2.64000D+02, 2.29000D+02, 2.17000D+02, 2.17000D+02,
     M  2.23000D+02, 2.86760D+02,-2.93940D+00, 2.50000D+00, 0.00000D+00/
      END

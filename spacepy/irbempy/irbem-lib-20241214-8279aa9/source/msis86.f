!***************************************************************************************************
! Copyright 1988, D. Bilitza
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
C MSIS86.FOR	D. BILITZA	10/88
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C SUBROUTINES AND FUNCTIONS:
C	GTS5, DENSS, GLOBE5, TSELEC, GLOB5L, DNET, CCOR, PRMSG5,GGM
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C
C***********************************************************************
      SUBROUTINE GTS5(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
      IMPLICIT NONE
C        MSIS-86/CIRA 1986 Neutral Thermosphere Model
C         A.E.Hedin 3/15/85;2/26/87 (Variable Names Shortened)
C         10/14/87 increase altitude limit of O mixing calculation
C             ALTL(2) from 300.0 to 400.0 km .
C     INPUT:
C        IYD - YEAR AND DAY AS YYYYDDD
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM) (GREATER THAN 85 KM)
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
C      TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
C
C          ADDITIONAL COMMENTS
C           (1) LOWER BOUND QUANTITIES IN COMMON/GTS3C/
C           (2) TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
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
C              18 - ALL T0 VAR           19 - ALL S VAR
C              20 - ALL Z0 VAR           21 - ALL NLB VAR
C              22 - ALL TR12 VAR         23 - TURBO SCALE HEIGHT VAR
C
C              To get current values of SW: CALL TRETRV(SW)
C
C !!!!!!!!!!!!!!!!!!!!!!!!!!! CHANGES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  - NAME,ISD,IST,ISDATE, and ISTIME were changed to character variables
C    in GTS5 and PRMSG5
C
C  - The variable dimension of P and AP in GLOBE5 and GLOBE5L was
C    indicted by *, rather than 1; if this does not work on your system
C    you may want to use P(150) and AP(7).
C
C  - The large data statement in PRMSG5 is now read in from file
C    MSIS86.DAT; some compilers do not allow named commons to be
C    initialized in a data statement.
C
C  - The first call to GLOBE5 should come before the common array SW(25)
C    is used in GTS5.
C
C Dieter Bilitza !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! March 87
C **********************************************************************
      INTEGER*4 IFL,IMR,MT(10),IIEE,ISW,MASS,J,I
      INTEGER*4 IYD
      LOGICAL METER
      REAL*8 DNET,ALTL(8),DD,ALT
      REAL*8 YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP(7),PT(150)
      REAL*8 SW(25),SWC(25),GGGG
      REAL*8 TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
      REAL*8 G0,RL,DB14,PTM(8),PDM(8,7)
      REAL*8 PD(150,7),PS(150),PDL(25,2),TINFG,GB,ROUT,TT(15)
      REAL*8 TINF,TR12,T(2),XMM,D(8),XMD
      REAL*8 G28,ZH28,B28,DM28,ZHM28
      REAL*8 G4,ZH04,B04,DM04,ZHM04,HC04,ZC04
      REAL*8 G16,ZH16,B16,DM16,ZHM16,HC16,ZC16,HCC16,ZCC16,RC16
      REAL*8 G32,ZH32,B32,DM32,ZHM32,HC32,ZC32
      REAL*8 G40,ZH40,B40,DM40,ZHM40,HC40,ZC40
      REAL*8 G1,ZH01,B01,DM01,ZHM01,HC01,ZC01,HCC01,ZCC01,RC01
      REAL*8 G14,ZH14,B14,DM14,ZHM14,HC14,ZC14,HCC14,ZCC14,RC14
      REAL*8 DENSS,GLOBE5,CCOR,GLOB5L,TZ,DDUM
c
      CHARACTER NAME(2)*4,ISDATE(3)*4,ISTIME(2)*4
c
      COMMON/UINR/IIEE
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     $ ,G0,RL,DD,DB14,TR12
      COMMON/LOWER5/PTM,PDM
      COMMON/PARM5/PT,PD,PS,PDL
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/TTEST/TINFG,GB,ROUT,TT
      COMMON/DATIME/ISDATE,ISTIME,NAME
      DATA MT/48,0,4,16,28,32,40,1,49,14/,IFL/0/
      DATA ALTL/200.D0,400.D0,150.D0,200.D0,240.D0,450.D0,
     &         320.D0,450.D0/
      DATA IMR/0/
c
      IF(IFL.EQ.0) THEN
	CALL PRMSG5
	IF(IIEE.GT.0) GOTO 9999
      	IFL=1
      ENDIF
      YRD=IYD*1.D0
C       Eq. A7
C!!OLD!! TINF=PTM(1)*(1.+SW(16)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,
C!!OLD!!$ AP,PT))*PT(1)
      GGGG=GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PT)
      TINF=PTM(1)*(1.D0+SW(16)*GGGG)*PT(1)
      ZA=PTM(5)*PDL(16,2)
C       Eq. A9
      T0=PTM(3)*PD(76,3)*(1.D0+SW(18)*GLOB5L(PD(76,3)))
C       Eq. A8
      TLB=PTM(2)*(1.D0+SW(17)*GLOB5L(PD(26,3)))*PD(26,3)
C       Eq. A10
      Z0=PTM(7)*(1.D0+SW(20)*GLOB5L(PD(51,3)))*PD(51,3)
C       Eq. A6
      G0=PTM(4)*PS(1)
     $ *(1.D0+SW(19)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,
     $ AP,PS))
C       Eq. A5
      S=G0/(TINF-TLB)
C       Eq. A11
      TR12=PD(101,3)*(1.D0+SW(22)*GLOB5L(PD(101,3)))
      T(1)=TINF
      IF(MASS.EQ.0) GO TO 50
C       Eq. A18  N2
      G28=SW(21)*GLOB5L(PD(1,3))
      YRD=IYD*1.D0
      T(1)=TINF
      XMM=PDM(5,3)
      DO 10 J = 1,10
      IF(MASS.EQ.MT(J))   GO TO 15
   10 CONTINUE
      WRITE(6,100) MASS
      GO TO 90
   15 IF(ALT.GT.ALTL(6).AND.MASS.NE.28.AND.MASS.NE.48) GO TO 17
C
C       **** N2 DENSITY ****
C
C       Eq. A18
      DB28 = PDM(1,3)*EXP(G28)*PD(1,3)
C       Eq. A13 - A17
      D(3)=DENSS(ALT,DB28,TINF,TLB, 28.D0,0.D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      DD=D(3)
C       Eq. A19
      ZH28=PDM(3,3)
      ZHM28=PDM(4,3)*PDL(6,2)
      XMD=28.D0-XMM
      B28=DENSS(ZH28,DB28,TINF,TLB,XMD,-1.D0,TZ,PTM(6),S,T0,ZA,Z0,TR12)
      IF(ALT.GT.ALTL(3).OR.SW(15).EQ.0.D0) GO TO 17
      DM28=DENSS(ALT,B28,TINF,TLB,XMM,0.D0,TZ,PTM(6),S,T0,ZA,Z0,TR12)
C       Eq. A12
      D(3)=DNET(D(3),DM28,ZHM28,XMM,28.D0)
   17 CONTINUE
      GO TO (20,50,20,25,90,35,40,45,25,48),  J
   20 CONTINUE
C
C       **** HE DENSITY ****
C
C       Eq. A18
      G4 = SW(21)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,1))
      DB04 = PDM(1,1)*EXP(G4)*PD(1,1)
C       Eq. A13 - A17
      D(1)=DENSS(ALT,DB04,TINF,TLB, 4.D0,-.4D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      DD=D(1)
      IF(ALT.GT.ALTL(1).OR.SW(15).EQ.0.D0) GO TO 24
C       Eq. A19
      ZH04=PDM(3,1)
      B04=DENSS(ZH04,DB04,TINF,TLB,4.D0-XMM,-1.4D0,
     $  T(2),PTM(6),S,T0,ZA,Z0,TR12)
      DM04=DENSS(ALT,B04,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,T0,ZA,Z0,TR12)
C       Eq. A12
      ZHM04=ZHM28
      D(1)=DNET(D(1),DM04,ZHM04,XMM,4.D0)
C       Eq. A20b
      RL=DLOG(B28*PDM(2,1)/B04)
C       Eq. A20a
      ZC04=PDM(5,1)*PDL(1,2)
      HC04=PDM(6,1)*PDL(2,2)
      D(1)=D(1)*CCOR(ALT,RL,HC04,ZC04)
   24 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   25 CONTINUE
C
C      **** O DENSITY ****
C
C       Eq. A18
      G16= SW(21)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,2))
      DB16 =  PDM(1,2)*EXP(G16)*PD(1,2)
C       Eq. A13 - A17
      D(2)=DENSS(ALT,DB16,TINF,TLB, 16.D0,0.D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      DD=D(2)
      IF(ALT.GT.ALTL(2).OR.SW(15).EQ.0.D0) GO TO 34
C  Corrected from PDM(3,1) to PDM(3,2)  12/2/85
C       Eq. A19
      ZH16=PDM(3,2)
      B16=DENSS(ZH16,DB16,TINF,TLB,16.D0-XMM,-1.D0,
     $  T(2),PTM(6),S,T0,ZA,Z0,TR12)
      DM16=DENSS(ALT,B16,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,T0,ZA,Z0,TR12)
C       Eq. A12
      ZHM16=ZHM28
      D(2)=DNET(D(2),DM16,ZHM16,XMM,16.D0)
C       Eq. A20b
      RL=DLOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
C       Eq. A20a
      HC16=PDM(6,2)*PDL(4,2)
      ZC16=PDM(5,2)*PDL(3,2)
      D(2)=D(2)*CCOR(ALT,RL,HC16,ZC16)
C       Eq. A21
      HCC16=PDM(8,2)*PDL(14,2)
      ZCC16=PDM(7,2)*PDL(13,2)
      RC16=PDM(4,2)*PDL(15,2)
      D(2)=D(2)*CCOR(ALT,RC16,HCC16,ZCC16)
   34 CONTINUE
      IF(MASS.NE.48 .AND. MASS.NE.49) GO TO 90
   35 CONTINUE
C
C       **** O2 DENSITY ****
C
C       Eq. A18
      G32= SW(21)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,4))
      DB32 = PDM(1,4)*EXP(G32)*PD(1,4)
C       Eq. A13 - A17
      D(4)=DENSS(ALT,DB32,TINF,TLB, 32.D0,0.D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      IF(MASS.EQ.49) THEN
         DD=DD+2.D0*D(4)
      ELSE
         DD=D(4)
      ENDIF
      IF(ALT.GT.ALTL(4).OR.SW(15).EQ.0.D0) GO TO 39
C       Eq. A19
      ZH32=PDM(3,4)
      B32=DENSS(ZH32,DB32,TINF,TLB,32.D0-XMM,-1.D0,
     $  T(2),PTM(6),S,T0,ZA,Z0,TR12)
      DM32=DENSS(ALT,B32,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,T0,ZA,Z0,TR12)
C       Eq. A12
      ZHM32=ZHM28
      D(4)=DNET(D(4),DM32,ZHM32,XMM,32.D0)
C       Eq. A20b
      RL=DLOG(B28*PDM(2,4)/B32)
C       Eq. A20a
      HC32=PDM(6,4)*PDL(8,2)
      ZC32=PDM(5,4)*PDL(7,2)
      D(4)=D(4)*CCOR(ALT,RL,HC32,ZC32)
   39 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   40 CONTINUE
C
C       **** AR DENSITY ****
C
C       Eq. A18
      G40= SW(21)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,5))
      DB40 = PDM(1,5)*EXP(G40)*PD(1,5)
C       Eq. A13 - A17
      D(5)=DENSS(ALT,DB40,TINF,TLB, 40.D0,0.D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      DD=D(5)
      IF(ALT.GT.ALTL(5).OR.SW(15).EQ.0.D0) GO TO 44
C       Eq. A19
      ZH40=PDM(3,5)
      B40=DENSS(ZH40,DB40,TINF,TLB,40.D0-XMM,-1.D0,
     $  T(2),PTM(6),S,T0,ZA,Z0,TR12)
      DM40=DENSS(ALT,B40,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,T0,ZA,Z0,TR12)
C       Eq. A12
      ZHM40=ZHM28
      D(5)=DNET(D(5),DM40,ZHM40,XMM,40.D0)
C       Eq. A20b
      RL=DLOG(B28*PDM(2,5)/B40)
C       Eq. A20a
      HC40=PDM(6,5)*PDL(10,2)
      ZC40=PDM(5,5)*PDL(9,2)
      D(5)=D(5)*CCOR(ALT,RL,HC40,ZC40)
   44 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   45 CONTINUE
C
C        **** HYDROGEN DENSITY ****
C
C       Eq. A18
      G1 = SW(21)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,6))
      DB01 = PDM(1,6)*EXP(G1)*PD(1,6)
C       Eq. A13 - A17
      D(7)=DENSS(ALT,DB01,TINF,TLB,1.D0,-.4D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      DD=D(7)
      IF(ALT.GT.ALTL(7).OR.SW(15).EQ.0.D0) GO TO 47
C       Eq. A19
      ZH01=PDM(3,6)
      B01=DENSS(ZH01,DB01,TINF,TLB,1.D0-XMM,-1.4D0,
     $  T(2),PTM(6),S,T0,ZA,Z0,TR12)
      DM01=DENSS(ALT,B01,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,T0,ZA,Z0,TR12)
C       Eq. A12
      ZHM01=ZHM28
      D(7)=DNET(D(7),DM01,ZHM01,XMM,1.D0)
C       Eq. A20b
      RL=DLOG(B28*PDM(2,6)*ABS(PDL(18,2))/B01)
C       Eq. A20a
      HC01=PDM(6,6)*PDL(12,2)
      ZC01=PDM(5,6)*PDL(11,2)
      D(7)=D(7)*CCOR(ALT,RL,HC01,ZC01)
C       Eq. A21
      HCC01=PDM(8,6)*PDL(20,2)
      ZCC01=PDM(7,6)*PDL(19,2)
      RC01=PDM(4,6)*PDL(21,2)
      D(7)=D(7)*CCOR(ALT,RC01,HCC01,ZCC01)
   47 CONTINUE
   48 CONTINUE
C
C        **** ATOMIC NITROGEN DENSITY ****
C
C       Eq. A18
      G14 = SW(21)*GLOBE5(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,7))
      DB14 = PDM(1,7)*EXP(G14)*PD(1,7)
C       Eq. A13 - A17
      D(8)=DENSS(ALT,DB14,TINF,TLB,14.D0,0.D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      DD=D(8)
      IF(ALT.GT.ALTL(8).OR.SW(15).EQ.0.D0) GO TO 49
C       Eq. A19
      ZH14=PDM(3,7)
      B14=DENSS(ZH14,DB14,TINF,TLB,14.D0-XMM,-1.D0,
     $  T(2),PTM(6),S,T0,ZA,Z0,TR12)
      DM14=DENSS(ALT,B14,TINF,TLB,XMM,0.D0,T(2),PTM(6),S,T0,ZA,Z0,TR12)
C       Eq. A12
      ZHM14=ZHM28
      D(8)=DNET(D(8),DM14,ZHM14,XMM,14.D0)
C       Eq. A20b
      RL=DLOG(B28*PDM(2,7)*ABS(PDL(3,1))/B14)
C       Eq. A20a
      HC14=PDM(6,7)*PDL(2,1)
      ZC14=PDM(5,7)*PDL(1,1)
      D(8)=D(8)*CCOR(ALT,RL,HC14,ZC14)
C       Eq. A21
      HCC14=PDM(8,7)*PDL(5,1)
      ZCC14=PDM(7,7)*PDL(4,1)
      RC14=PDM(4,7)*PDL(6,1)
      D(8)=D(8)*CCOR(ALT,RC14,HCC14,ZCC14)
   49 CONTINUE
      IF(MASS.NE.48) GO TO 90
C
C       TOTAL MASS DENSITY
C
      D(6) = 1.66D-24*(4.D0*D(1)+16.D0*D(2)+28.D0*D(3)+32.D0*D(4)+
     $       40.D0*D(5)+D(7)+14.*D(8))
      DB48=1.66D-24*(4.D0*DB04+16.D0*DB16+28.D0*DB28+32.D0*DB32+40.D0
     $        *DB40+DB01+14.D0*DB14)
      GO TO 90
   50 DDUM  = DENSS(ALT,1.D0, TINF,TLB,0.D0,0.D0,T(2),PTM(6),
     &S,T0,ZA,Z0,TR12)
      GO TO 90
   90 CONTINUE
      IF(IMR.EQ.1) THEN
        DO 95 I=1,8
          D(I)=D(I)*1.D6
   95   CONTINUE
        D(6)=D(6)/1000.D6
      ENDIF
      RETURN
  100 FORMAT(1X,'MASS', I5, '  NOT VALID')
      ENTRY METERS(METER)
      IMR=0
      IF(METER) IMR=1
9999  CONTINUE
      END
C--------------------------------------------------------------------
      FUNCTION DENSS(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,T0,ZA,Z0,TR12)
      IMPLICIT NONE
C       Calculate Temperature and Density Profiles for MSIS models
c
      INTEGER*4 MP,II,JG,LT,IERR,IFUN,N,J
      REAL*8 ZZ,ZL,ZETA,RGAS,DENSS,GSURF,RE
      REAL*8 ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,T0,ZA,Z0,TR12
      REAL*8 QPB(50),DV(60)
      REAL*8 Z,ZG2,TT,TA,ZG0,DTA,T12,ZG1,DD,CC,BB,X,X2,TAF
      REAL*8 GLB,GAMMA,DENSA,GAMM
c
      COMMON/PARMB/GSURF,RE
      COMMON/FIT/TAF
      COMMON/LSQV/MP,II,JG,LT,QPB,IERR,IFUN,N,J,DV
      DATA RGAS/831.4D0/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
c
      DENSS=1.D0
      Z=DMAX1(ALT,ZA)
C      Eq. A4a
      ZG2=ZETA(Z,ZLB)
C      Eq. A1a
      TT=TINF-(TINF-TLB)*EXP(-S2*ZG2)
      TA=TT
      TZ=TT
      DENSS=TZ
      IF(ALT.GE.ZA) GO TO 10
C      Eq. A4b
      ZG0=ZETA(Z0,ZA)
C      Eq. A2b
      DTA=(TINF-TA)*S2*((RE+ZLB)/(RE+ZA))**2
C      Eq. A3e
      T12=T0+TR12*(TA-T0)
C      Eq. A4b
      ZG1=ZETA(ALT,ZA)
C       CALCULATE TEMPERATURE BELOW ZA
C      Eq. A3a
      DD=0.666666D0*ZG0*DTA/TA**2 - 3.11111D0*(1.D0/TA-1.D0/T0)+
     $ 7.11111D0*(1.D0/T12-1.D0/T0)
C      Eq. A3b
      CC=ZG0*DTA/(2.D0*TA*TA) - (1.D0/TA-1.D0/T0) - 2.D0*DD
C      Eq. A3c
      BB=(1.D0/TA-1.D0/T0) - CC - DD
C      Eq. A3d
      X=(-(ZG1-ZG0)/ZG0)
C      Eq. A1b
      X2=X*X
      TZ=1.D0/(1.D0/T0+BB*X2+CC*X2*X2+DD*X2*X2*X2)
      DENSS=TZ
      TAF=(T12-T0)/(TA-T0)
   10 IF(XM.EQ.0.D0) GO TO 50
      IF(TA.GT.0.D0 .AND. TZ.GT.0.D0) GO TO 20
         WRITE(6,*)ALT,XM,TINF,TLB,T0,TA,II,JG,N,DV(J),IFUN,S2,ZG0,TZ
         TT=TLB
         TA=TLB
         TZ=TLB
   20 CONTINUE
C      CALCULATE DENSITY ABOVE ZA
C      Eq. A17a
      GLB=GSURF/(1.D0+ZLB/RE)**2
C      Eq. A16a
      GAMMA=XM*GLB/(S2*RGAS*TINF)
C      Eq. A13, A14a, & A15
      DENSA=DLB*(TLB/TT)**(1.D0+ALPHA+GAMMA)*EXP(-S2*GAMMA*ZG2)
      DENSS=DENSA
      IF(ALT.GE.ZA) GO TO 50
C      CALCULATE DENSITY BELOW ZA
C      Eq. A17b
      GLB=GSURF/(1.D0+ZA/RE)**2
C      Eq. A16b
      GAMM=XM*GLB*ZG0/RGAS
C      Eq. A13, A14b, & A15
      DENSS=DENSA*(TA/TZ)**(1.D0+ALPHA)*
     $ EXP(GAMM*((X-1.D0)/T0+BB*(X*X2-1.D0)/3.D0+CC*(X2*X2*X-1.D0)/5.D0+
     $ DD*(X2*X2*X2*X-1.D0)/7.D0))
   50 CONTINUE
CCCCCCWRITE(6,100)CXM,ALT,ZA,TINF,TLB,S2,T0,S1,TA,TZ,DLB,DENSA,DENSS
CC100 FORMAT(' D',1P13E10.2)
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION GLOBE5(YRD,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)
      IMPLICIT NONE
C       CALCULATE G(L) FUNCTION FOR MSIS-86/CIRA 1986
C       Upper Thermosphere Parameters
      INTEGER*4 IYR,ISW,NSW,I
C
      REAL*8 LAT, LONG
      REAL*8 P(*),SV(25),AP(*)
      REAL*8 PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC
      REAL*8 DAY,DF,DFA,APD,APDF,APT(4),SW(25),SWC(25)
      REAL*8 DGTR,DR,XL,TLL,DAYL,P14,P18,P32,HR,SR,P39
      REAL*8 T(15)
      REAL*8 YRD,SEC,TLOC,F107A,F107
      REAL*8 C,S,C2,C4,S2
      REAL*8 CD14,C2D14,CD18,CD32,CD39,F1,F2,T71,T72,T81,T82
      REAL*8 P44,P45,TINF
      REAL*8 A,EX,G0,SUMEX,SG0,GLOBE5,GB,ROUT,EXP1,EXP2
      REAL*8 XLONG,CLONG,SLONG ! not used but added for consistency with newer MSIS models
c
      COMMON/TTEST/TINF,GB,ROUT,T
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/LPOLY/PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     &             DAY,DF,DFA,APD,APDF,APT,XLONG,CLONG,SLONG
      COMMON/LPOLYI/IYR
      DATA DGTR/1.74533D-2/,DR/1.72142D-2/, XL/1000.D0/,TLL/1000.D0/
     &  DAYL/-1.D0/,P14/-1000.D0/,P18/-1000.D0/,P32/-1000.D0/
     &  HR/.2618D0/,SR/7.2722D-5/,SV/25*1.D0/,NSW/14/,P39/-1000.D0/
c
C Eq. A24d
      G0(A)=(A-4.D0+(P(26)-1.D0)*(A-4.D0+(EXP(-ABS(P(25))*
     &  (A-4.D0))-1.D0)/ABS(P(25))))
C Eq. A24c
      SUMEX(EX)=1.D0+(1.D0-EX**19)/(1.D0-EX)*EX**(.5D0)
C Eq. A24a
      SG0(EX)=(G0(AP(2))+(G0(AP(3))*EX+G0(AP(4))*EX*EX+G0(AP(5))*EX**3
     $ +(G0(AP(6))*EX**4+G0(AP(7))*EX**12)*(1.D0-EX**8)/(1.D0-EX)))
     $ /SUMEX(EX)
      IF(ISW.NE.64999) CALL TSELEC(SV)
      T(10) = 0.D0
      T(11) = 0.D0
      T(12) = 0.D0
      T(13)=0.D0
   10 CONTINUE
      IYR = INT(YRD/1000.D0)
      DAY = YRD - IYR*1000.D0
C Eq. A22 (remainder of code)
      IF(XL.EQ.LAT)   GO TO 15
C CALCULATE LEGENDRE POLYNOMIALS
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
      PLG(2,2) = S
      PLG(3,2) = 3.D0*C*S
      PLG(4,2) = 1.5D0*(5.D0*C2-1.D0)*S
      PLG(5,2) = 2.5D0*(7.D0*C2*C-3.D0*C)*S
      PLG(6,2) = 1.875D0*(21.D0*C4 - 14.D0*C2 +1.D0)*S
      PLG(7,2) = (11.D0*C*PLG(6,2)-6.D0*PLG(5,2))/5.D0
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
      STLOC = SIN(HR*TLOC)
      CTLOC = COS(HR*TLOC)
      S2TLOC = SIN(2.D0*HR*TLOC)
      C2TLOC = COS(2.D0*HR*TLOC)
      S3TLOC = SIN(3.D0*HR*TLOC)
      C3TLOC = COS(3.D0*HR*TLOC)
      TLL = TLOC
   16 CONTINUE
      IF(DAY.NE.DAYL.OR.P(14).NE.P14) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P(14).NE.P14) C2D14=COS(DR*2.D0*(DAY-P(14)))
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
     $ + P(30)*DFA**2
      F1 = 1.D0 + (P(48)*DFA +P(20)*DF+P(21)*DF*DF)*SWC(1)
      F2 = 1.D0 + (P(50)*DFA+P(20)*DF+P(21)*DF*DF)*SWC(1)
C        TIME INDEPENDENT
      T(2) =
     1  (P(2)*PLG(3,1) + P(3)*PLG(5,1)+P(23)*PLG(7,1))
     $ +(P(15)*PLG(3,1))*DFA*SWC(1)
     2 +P(27)*PLG(2,1)
C        SYMMETRICAL ANNUAL
      T(3) =
     1 (P(19) )*CD32
C        SYMMETRICAL SEMIANNUAL
      T(4) =
     1 (P(16)+P(17)*PLG(3,1))*CD18
C        ASYMMETRICAL ANNUAL
      T(5) =  F1*
     1  (P(10)*PLG(2,1) + P(11)*PLG(4,1))*CD14
C         ASYMMETRICAL SEMIANNUAL
      T(6) =    P(38)*PLG(2,1)*CD39
C        DIURNAL
      T71 = (P(12)*PLG(3,2) + P(36)*PLG(2,2))*CD14*SWC(5)
      T72 = (P(13)*PLG(3,2) + P(37)*PLG(2,2))*CD14*SWC(5)
      T(7) = F2*
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2) + P(28)*PLG(6,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2) +P(29)*PLG(6,2)
     5 + T72)*STLOC)
C        SEMIDIURNAL
      T81 = (P(24)*PLG(4,3))*CD14*SWC(5)
      T82 = (P(34)*PLG(4,3))*CD14*SWC(5)
      T(8) = F2*
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
C        TERDIURNAL
      T(14) = F2*
     1 ((P(40)*PLG(4,4)+(P(94)*PLG(5,4)+P(47)*PLG(7,4))*CD14*SWC(5))*
     $ S3TLOC
     2 +(P(41)*PLG(4,4)+(P(95)*PLG(5,4)+P(49)*PLG(7,4))*CD14*SWC(5))*
     $ C3TLOC)
C          MAGNETIC ACTIVITY BASED ON DAILY AP
      IF(SW(9).EQ.-1.D0 .AND. P(52).NE.0.D0) GO TO 30
      APD=(AP(1)-4.D0)
      P44=P(44)
      P45=P(45)
      IF(P44.LT.0.D0) P44=1.D-5
      APDF = (APD+(P45-1.D0)*(APD+(EXP(-P44  *APD)-1.D0)/P44  ))
      T(9)=APDF*(P(33)+P(46)*PLG(3,1)+P(35)*PLG(5,1)+
     $ (P(101)*PLG(2,1)+P(102)*PLG(4,1)+P(103)*PLG(6,1))*CD14*SWC(5)+
     $ (P(122)*PLG(2,2)+P(123)*PLG(4,2)+P(124)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(125))))
      GO TO 40
   30 CONTINUE
      EXP1 = EXP(-10800.D0*ABS(P(52))/(1.D0+P(139)*(45.D0-ABS(LAT))))
      IF(EXP1.GT..99999D0) EXP1=.99999D0
      EXP2 = EXP(-10800.D0*ABS(P(54)))
      IF(EXP2.GT..99999D0) EXP2=.99999D0
      IF(P(25).LT.1.D-4) P(25)=1.D-4
      APT(1)=SG0(EXP1)
      APT(3)=SG0(EXP2)
      T(9) = APT(1)*(P(51)+P(97)*PLG(3,1)+P(55)*PLG(5,1)+
     $ (P(126)*PLG(2,1)+P(127)*PLG(4,1)+P(128)*PLG(6,1))*CD14*SWC(5)+
     $ (P(129)*PLG(2,2)+P(130)*PLG(4,2)+P(131)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(132))))
  40  CONTINUE
      IF(SW(10).EQ.0.D0 .OR. LONG.LE.-1000.D0) GO TO 49
C        LONGITUDINAL
      T(11)= (1.D0+P(90)*PLG(2,1))*(1.D0+P(81)*DFA*SWC(1))*
     $((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $ +P(104)*PLG(2,2)+P(105)*PLG(4,2)+P(106)*PLG(6,2)
     $ +SWC(5)*(P(110)*PLG(2,2)+P(111)*PLG(4,2)+P(112)*PLG(6,2))*CD14)*
     $     COS(DGTR*LONG)
     $ +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $ +P(107)*PLG(2,2)+P(108)*PLG(4,2)+P(109)*PLG(6,2)
     $ +SWC(5)*(P(113)*PLG(2,2)+P(114)*PLG(4,2)+P(115)*PLG(6,2))*CD14)*
     $  SIN(DGTR*LONG))
C        UT AND MIXED UT,LONGITUDE
      T(12)=(1.D0+P(96)*PLG(2,1))*(1.D0+P(82)*DFA*SWC(1))*
     $(1.D0+P(120)*PLG(2,1)*SWC(5)*CD14)*
     $((P(69)*PLG(2,1)+P(70)*PLG(4,1)+P(71)*PLG(6,1))*
     $     COS(SR*(SEC-P(72))))
      T(12)=T(12)+SWC(11)*
     $ (P(77)*PLG(4,3)+P(78)*PLG(6,3)+P(79)*PLG(8,3))*
     $     COS(SR*(SEC-P(80))+2.D0*DGTR*LONG)*(1.D0+P(138)*DFA*SWC(1))
C        UT,LONGITUDE MAGNETIC ACTIVITY
      IF(SW(9).EQ.-1.D0 .AND. P(52).NE.0.D0) GO TO 45
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
C  PARMS NOT USED: 60,83,100,140-150
   49 TINF = 0.D0
      IF(SW(9).EQ.-1.D0) TINF=P(31)
      DO 50 I = 1,NSW
   50 TINF = TINF + ABS(SW(I))*T(I)
      GLOBE5 = TINF
      RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE TSELEC(SV)
      IMPLICIT NONE
C        SET SWITCHES
      INTEGER*4 I,ISW
C
      REAL*8 SW(25),SWC(25),SV(25),SVV(25),SAV(25)
c
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
c
      DO 100 I = 1,25
        SAV(I)=SV(I)
        SW(I)=DMOD(SV(I),2.D0)
        IF(ABS(SV(I)).GT.0.D0) THEN
          SWC(I)=1.D0
        ELSE
          SWC(I)=0.D0
        ENDIF
  100 CONTINUE
      ISW=64999
      RETURN
      ENTRY TRETRV(SVV)
      DO 200 I=1,25
        SVV(I)=SAV(I)
  200 CONTINUE
      END
C--------------------------------------------------------------------
      FUNCTION GLOB5L(P)
      IMPLICIT NONE
C      LIMITED PARAMETER VERSION OF GLOBE 9/2/82
C       CALCULATE G(L) FUNCTION FOR MSIS-86/CIRA 1986
C       Lower Thermosphere Parameters
      INTEGER*4 I,IYR,ISW
c
      REAL*8 DR,T(15),DAYL,P7,P9,P11
      REAL*8 P(*),DAY,CD7,CD9,CD11
      REAL*8 PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC
      REAL*8 DF,DFA,APD,APDF,APT(4),SW(25),SWC(25)
      REAL*8 TT,GLOB5L
      REAL*8 XLONG,CLONG,SLONG ! not used but added for consistency with newer MSIS models

c
      COMMON/LPOLY/PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ DAY,DF,DFA,APD,APDF,APT,XLONG,CLONG,SLONG
      COMMON/CSW/SW,SWC
      COMMON/CSWI/ISW
      COMMON/LPOLYI/IYR
C!!OLD!! DIMENSION P(1),T(15) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      DATA DR/1.72142D-2/,T/15*0.D0/
      DATA DAYL/-1.D0/,P7/-1000.D0/,P9/-1000.D0/,P11/-1000.D0/
      IF(DAY.NE.DAYL.OR.P7.NE.P(7)) CD7=COS(DR*(DAY-P(7)))
      IF(DAY.NE.DAYL.OR.P9.NE.P(9)) CD9=COS(2.D0*DR*(DAY-P(9)))
      IF(DAY.NE.DAYL.OR.P11.NE.P(11)) CD11=COS(DR*(DAY-P(11)))
      DAYL=DAY
      P7=P(7)
      P9=P(9)
      P11=P(11)
C
      T(1)=P(2)*DFA
      T(2)=P(4)*PLG(3,1)
      T(3)=P(6)*CD7
      T(4)=(P(8) )*CD9
      T(5)=(P(10)*PLG(2,1)+P(22)*PLG(4,1))*CD11
      T(6)=0.D0
      T(7)=(P(14)*PLG(2,2)*CTLOC+P(15)*PLG(2,2)*STLOC)
      T(8)=(P(16)*PLG(3,3)+P(18)*PLG(5,3)
     $     +(P(20)*PLG(6,3))*CD11*SWC(5)
     $     )*C2TLOC
     $     +(P(17)*PLG(3,3)+P(19)*PLG(5,3)
     $     +(P(21)*PLG(6,3))*CD11*SWC(5)
     $     )*S2TLOC
      T(14)=(P(12)*PLG(4,4)*C3TLOC
     $     +P(25)*PLG(4,4)*S3TLOC)
      IF(SW(9).EQ.1.D0)
     $ T(9)=APDF*(P(23)+P(24)*PLG(3,1)*SWC(2))
      IF(SW(9).EQ.-1.D0)
     $ T(9)=(P(3)*APT(3)+P(5)*PLG(3,1)*APT(3)*SWC(2))
C       PARMS NOT USED: 13
      TT=0.D0
      DO 50 I=1,14
   50 TT=TT+ABS(SW(I))*T(I)
      GLOB5L=TT
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DNET(DD,DM,ZHM,XMM,XM)
      IMPLICIT NONE
C       8/20/80
C       TURBOPAUSE CORRECTION FOR MSIS MODELS
C       Eq. A12b
      REAL*8 A,ZHM,XMM,XM,YLOG
      REAL*8 DM,DD,DNET
c
      A=ZHM/(XMM-XM)
C       Eq. A12a
      YLOG=A*LOG(DM/DD)
      IF(YLOG.LT.-10.D0) GO TO 10
      IF(YLOG.GT.10.D0)  GO TO 20
        DNET=DD*(1.D0+EXP(YLOG))**(1/A)
        GO TO 50
   10 CONTINUE
        DNET=DD
        GO TO 50
   20 CONTINUE
        DNET=DM
        GO TO 50
   50 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION  CCOR(ALT, R,H1,ZH)
      IMPLICIT NONE
C        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C     Eq. A20a or Eq. A21
      REAL*8 E,ALT,ZH,H1,EX,CCOR,R
c
      E=(ALT-ZH)/H1
      IF(E.GT.70.D0) GO TO 20
      IF(E.LT.-70.D0) GO TO 10
        EX=EXP(E)
        CCOR=R/(1.D0+EX)
        GO TO 50
   10   CCOR=R
        GO TO 50
   20   CCOR=0.D0
        GO TO 50
   50 CONTINUE
      CCOR=EXP(CCOR)
       RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE PRMSG5
      IMPLICIT NONE
C          CIRA     11-FEB-86
      INTEGER*4 IIEE,I
c
      REAL*8 PT1(50),PT2(50),PT3(50),PA1(50),PA2(50),PA3(50)
      REAL*8 PB1(50),PB2(50),PB3(50),PC1(50),PC2(50),PC3(50)
      REAL*8 PD1(50),PD2(50),PD3(50),PE1(50),PE2(50),PE3(50)
      REAL*8 PF1(50),PF2(50),PF3(50),PG1(50),PG2(50),PG3(50)
      REAL*8 PH1(50),PH2(50),PH3(50),PI1(50)
      REAL*8 PTM(8),PDM(8,7)
      REAL*8 GSURF,RE
c
      CHARACTER ISD(3)*4,IST(2)*4,NAME(2)*4,ISDATE(3)*4,ISTIME(2)*4
c
      COMMON/UINR/IIEE
      COMMON/PARMB/GSURF,RE
      COMMON/PARM5/PT1,PT2,PT3,PA1,PA2,PA3,
     * PB1,PB2,PB3,PC1,PC2,PC3,
     * PD1,PD2,PD3,PE1,PE2,PE3,
     * PF1,PF2,PF3,PG1,PG2,PG3,
     * PH1,PH2,PH3,PI1
      COMMON/LOWER5/PTM,PDM/DATIME/ISDATE,ISTIME,NAME
      DATA ISD/'11-F','EB-8','6   '/,IST/'18:2','3:31'/
c
      GSURF=980.665D0
      RE=6356.77D0
      IIEE=0
c
       PT1( 1)= 0.996040D+00
       PT1( 2)= 0.385528D-01
       PT1( 3)= 0.303445D-02
       PT1( 4)=-0.105531D+00
       PT1( 5)=-0.607134D-02
       PT1( 6)=-0.516278D-03
       PT1( 7)=-0.115622D+00
       PT1( 8)= 0.202240D-02
       PT1( 9)= 0.990156D-02
       PT1(10)=-0.127371D+00
       PT1(11)=-0.302449D-01
       PT1(12)= 0.123512D-01
       PT1(13)=-0.526277D-02
       PT1(14)=-0.845398D+01
       PT1(15)= 0.000000D+00
       PT1(16)= 0.142370D-01
       PT1(17)= 0.000000D+00
       PT1(18)= 0.125818D+03
       PT1(19)= 0.805486D-02
       PT1(20)= 0.164419D-02
       PT1(21)=-0.621452D-05
       PT1(22)= 0.311701D-02
       PT1(23)= 0.000000D+00
       PT1(24)= 0.386578D-02
       PT1(25)= 0.132397D+00
       PT1(26)= 0.213315D+00
       PT1(27)= 0.000000D+00
       PT1(28)= 0.000000D+00
       PT1(29)= 0.000000D+00
       PT1(30)=-0.641110D-05
       PT1(31)= 0.000000D+00
       PT1(32)= 0.300150D+02
       PT1(33)= 0.533297D-02
       PT1(34)= 0.389146D-02
       PT1(35)= 0.204725D-02
       PT1(36)= 0.000000D+00
       PT1(37)= 0.000000D+00
       PT1(38)=-0.192645D-01
       PT1(39)= 0.275905D+01
       PT1(40)= 0.147284D-02
       PT1(41)= 0.341345D-03
       PT1(42)=-0.117388D-02
       PT1(43)=-0.354589D-03
       PT1(44)= 0.113139D+00
       PT1(45)= 0.169134D+00
       PT1(46)= 0.508295D-02
       PT1(47)= 0.365016D-04
       PT1(48)= 0.426385D-02
       PT1(49)= 0.115102D-03
       PT1(50)= 0.511819D-02
       PT2( 1)= 0.609108D-02
       PT2( 2)= 0.404995D-04
       PT2( 3)= 0.153049D-02
       PT2( 4)= 0.241470D-04
       PT2( 5)= 0.230764D-02
       PT2( 6)= 0.155267D-02
       PT2( 7)= 0.133722D-02
       PT2( 8)=-0.182318D-02
       PT2( 9)=-0.263007D+03
       PT2(10)= 0.000000D+00
       PT2(11)= 0.137337D-02
       PT2(12)= 0.995774D-03
       PT2(13)= 0.000000D+00
       PT2(14)=-0.108983D+03
       PT2(15)= 0.562606D-02
       PT2(16)= 0.594053D-02
       PT2(17)= 0.109358D-02
       PT2(18)= 0.000000D+00
       PT2(19)=-0.133410D-01
       PT2(20)=-0.243409D-01
       PT2(21)=-0.135688D-01
       PT2(22)= 0.311370D+05
       PT2(23)= 0.000000D+00
       PT2(24)= 0.000000D+00
       PT2(25)= 0.000000D+00
       PT2(26)=-0.283023D+04
       PT2(27)= 0.845583D-03
       PT2(28)= 0.538706D-03
       PT2(29)= 0.000000D+00
       PT2(30)= 0.247956D+03
       PT2(31)= 0.292246D-02
       PT2(32)= 0.000000D+00
       PT2(33)= 0.000000D+00
       PT2(34)= 0.747703D-04
       PT2(35)= 0.887993D-03
       PT2(36)= 0.000000D+00
       PT2(37)= 0.000000D+00
       PT2(38)= 0.000000D+00
       PT2(39)= 0.000000D+00
       PT2(40)= 0.000000D+00
       PT2(41)=-0.116540D-01
       PT2(42)=-0.449173D-02
       PT2(43)=-0.353189D-03
       PT2(44)=-0.173933D-03
       PT2(45)=-0.153218D-03
       PT2(46)=-0.565411D+00
       PT2(47)= 0.777272D-02
       PT2(48)=-0.911784D+02
       PT2(49)= 0.645187D-03
       PT2(50)= 0.000000D+00
       PT3( 1)=-0.837685D-03
       PT3( 2)= 0.242318D-02
       PT3( 3)= 0.473796D-02
       PT3( 4)=-0.301801D-02
       PT3( 5)=-0.423564D-02
       PT3( 6)=-0.248289D-02
       PT3( 7)= 0.919286D-03
       PT3( 8)= 0.216372D-02
       PT3( 9)= 0.863968D-03
       PT3(10)= 0.189689D-02
       PT3(11)= 0.415654D-02
       PT3(12)= 0.000000D+00
       PT3(13)= 0.118068D-01
       PT3(14)= 0.331190D-02
       PT3(15)= 0.000000D+00
       PT3(16)= 0.120222D-02
       PT3(17)= 0.000000D+00
       PT3(18)= 0.000000D+00
       PT3(19)=-0.307246D+01
       PT3(20)= 0.000000D+00
       PT3(21)= 0.000000D+00
       PT3(22)= 0.672403D-03
       PT3(23)= 0.108930D-02
       PT3(24)= 0.972278D-03
       PT3(25)= 0.468242D+01
       PT3(26)=-0.315034D-03
       PT3(27)= 0.400059D-02
       PT3(28)= 0.515036D-02
       PT3(29)= 0.162989D-02
       PT3(30)= 0.108824D-02
       PT3(31)= 0.995261D-03
       PT3(32)= 0.418955D+01
       PT3(33)=-0.364059D+00
       PT3(34)= 0.170182D-02
       PT3(35)= 0.000000D+00
       PT3(36)= 0.000000D+00
       PT3(37)=-0.320120D+01
       PT3(38)= 0.000000D+00
       PT3(39)= 0.580206D-02
       PT3(40)= 0.000000D+00
       PT3(41)= 0.000000D+00
       PT3(42)= 0.000000D+00
       PT3(43)= 0.000000D+00
       PT3(44)= 0.000000D+00
       PT3(45)= 0.000000D+00
       PT3(46)= 0.000000D+00
       PT3(47)= 0.000000D+00
       PT3(48)= 0.000000D+00
       PT3(49)= 0.000000D+00
       PT3(50)= 0.000000D+00
       PA1( 1)= 0.104934D+01
       PA1( 2)=-0.288362D-01
       PA1( 3)=-0.207095D+00
       PA1( 4)=-0.103314D+00
       PA1( 5)=-0.702373D-02
       PA1( 6)= 0.129664D-01
       PA1( 7)= 0.408853D+00
       PA1( 8)=-0.919895D-02
       PA1( 9)=-0.188660D-01
       PA1(10)= 0.140927D+01
       PA1(11)= 0.175033D+00
       PA1(12)= 0.187351D-01
       PA1(13)= 0.110979D+00
       PA1(14)=-0.742871D+01
       PA1(15)= 0.000000D+00
       PA1(16)= 0.267143D+00
       PA1(17)=-0.595979D-01
       PA1(18)= 0.105038D+03
       PA1(19)=-0.840963D-01
       PA1(20)=-0.697632D-03
       PA1(21)= 0.206521D-05
       PA1(22)= 0.765306D-03
       PA1(23)= 0.000000D+00
       PA1(24)= 0.000000D+00
       PA1(25)= 0.126762D+00
       PA1(26)= 0.128876D+00
       PA1(27)=-0.504479D-01
       PA1(28)=-0.130735D-01
       PA1(29)=-0.224348D-01
       PA1(30)= 0.000000D+00
       PA1(31)= 0.000000D+00
       PA1(32)=-0.150832D+03
       PA1(33)=-0.629928D-02
       PA1(34)= 0.000000D+00
       PA1(35)=-0.407760D-02
       PA1(36)= 0.000000D+00
       PA1(37)= 0.000000D+00
       PA1(38)= 0.525725D-01
       PA1(39)=-0.311486D+02
       PA1(40)=-0.313351D-02
       PA1(41)= 0.275838D-02
       PA1(42)= 0.000000D+00
       PA1(43)= 0.000000D+00
       PA1(44)= 0.111247D+00
       PA1(45)= 0.108815D+00
       PA1(46)=-0.466713D-01
       PA1(47)= 0.000000D+00
       PA1(48)=-0.329329D-02
       PA1(49)= 0.000000D+00
       PA1(50)= 0.167838D-02
       PA2( 1)=-0.916691D-02
       PA2( 2)= 0.345044D-04
       PA2( 3)=-0.971806D-02
       PA2( 4)= 0.000000D+00
       PA2( 5)=-0.204672D-02
       PA2( 6)=-0.786899D-02
       PA2( 7)=-0.798285D-02
       PA2( 8)= 0.536515D-02
       PA2( 9)=-0.531172D+04
       PA2(10)= 0.000000D+00
       PA2(11)=-0.642781D-02
       PA2(12)=-0.171690D-02
       PA2(13)= 0.000000D+00
       PA2(14)=-0.679131D+02
       PA2(15)=-0.179912D-01
       PA2(16)=-0.158305D-01
       PA2(17)=-0.712313D-02
       PA2(18)= 0.000000D+00
       PA2(19)= 0.253477D-01
       PA2(20)= 0.852960D-01
       PA2(21)= 0.102163D+00
       PA2(22)= 0.295009D+05
       PA2(23)= 0.000000D+00
       PA2(24)= 0.000000D+00
       PA2(25)= 0.000000D+00
       PA2(26)=-0.684625D+04
       PA2(27)=-0.619098D-02
       PA2(28)=-0.269289D-02
       PA2(29)= 0.000000D+00
       PA2(30)=-0.520231D+03
       PA2(31)=-0.633463D-02
       PA2(32)= 0.000000D+00
       PA2(33)= 0.000000D+00
       PA2(34)=-0.602428D-02
       PA2(35)=-0.407077D-02
       PA2(36)= 0.542264D-02
       PA2(37)= 0.000000D+00
       PA2(38)= 0.000000D+00
       PA2(39)= 0.000000D+00
       PA2(40)= 0.000000D+00
       PA2(41)= 0.407560D-01
       PA2(42)= 0.282288D-01
       PA2(43)= 0.908088D-02
       PA2(44)= 0.000000D+00
       PA2(45)= 0.000000D+00
       PA2(46)=-0.405204D+00
       PA2(47)=-0.597931D-01
       PA2(48)=-0.731823D+02
       PA2(49)=-0.206620D-02
       PA2(50)= 0.000000D+00
       PA3( 1)=-0.372723D-02
       PA3( 2)=-0.188146D-01
       PA3( 3)=-0.101794D-01
       PA3( 4)= 0.804633D-02
       PA3( 5)= 0.101090D-01
       PA3( 6)= 0.873253D-02
       PA3( 7)= 0.238268D-01
       PA3( 8)= 0.480444D-02
       PA3( 9)= 0.171088D-02
       PA3(10)= 0.396369D-01
       PA3(11)=-0.213809D-01
       PA3(12)= 0.000000D+00
       PA3(13)=-0.102588D+00
       PA3(14)=-0.591702D-02
       PA3(15)= 0.000000D+00
       PA3(16)= 0.270923D-02
       PA3(17)= 0.000000D+00
       PA3(18)= 0.000000D+00
       PA3(19)=-0.175043D+03
       PA3(20)= 0.603489D+00
       PA3(21)=-0.617589D+00
       PA3(22)= 0.838098D-02
       PA3(23)= 0.183871D-02
       PA3(24)=-0.705329D-03
       PA3(25)=-0.406644D+01
       PA3(26)=-0.509347D-02
       PA3(27)=-0.284344D-01
       PA3(28)=-0.124160D-01
       PA3(29)= 0.133665D-01
       PA3(30)= 0.393410D-02
       PA3(31)=-0.503723D-03
       PA3(32)=-0.457683D+01
       PA3(33)=-0.529542D+00
       PA3(34)=-0.425812D-02
       PA3(35)= 0.000000D+00
       PA3(36)= 0.000000D+00
       PA3(37)= 0.191541D+02
       PA3(38)= 0.000000D+00
       PA3(39)= 0.323247D-02
       PA3(40)= 0.000000D+00
       PA3(41)= 0.000000D+00
       PA3(42)= 0.000000D+00
       PA3(43)= 0.000000D+00
       PA3(44)= 0.000000D+00
       PA3(45)= 0.000000D+00
       PA3(46)= 0.000000D+00
       PA3(47)= 0.000000D+00
       PA3(48)= 0.000000D+00
       PA3(49)= 0.000000D+00
       PA3(50)= 0.000000D+00
       PB1( 1)= 0.931113D+00
       PB1( 2)=-0.138721D+00
       PB1( 3)=-0.133457D+00
       PB1( 4)=-0.529542D-01
       PB1( 5)=-0.444983D-02
       PB1( 6)= 0.135264D-01
       PB1( 7)= 0.598075D-01
       PB1( 8)=-0.362880D-01
       PB1( 9)=-0.312798D-01
       PB1(10)= 0.372068D+00
       PB1(11)= 0.295974D-01
       PB1(12)= 0.120509D-01
       PB1(13)= 0.521995D-01
       PB1(14)=-0.778888D+01
       PB1(15)= 0.000000D+00
       PB1(16)= 0.118634D+00
       PB1(17)=-0.204495D-01
       PB1(18)= 0.103280D+03
       PB1(19)= 0.982432D-01
       PB1(20)= 0.477694D-03
       PB1(21)= 0.000000D+00
       PB1(22)= 0.274372D-02
       PB1(23)= 0.000000D+00
       PB1(24)= 0.000000D+00
       PB1(25)= 0.757809D-01
       PB1(26)= 0.171403D+00
       PB1(27)=-0.105205D-01
       PB1(28)= 0.000000D+00
       PB1(29)= 0.000000D+00
       PB1(30)= 0.000000D+00
       PB1(31)= 0.000000D+00
       PB1(32)=-0.873348D+01
       PB1(33)=-0.581094D-02
       PB1(34)= 0.000000D+00
       PB1(35)=-0.814944D-02
       PB1(36)= 0.000000D+00
       PB1(37)= 0.000000D+00
       PB1(38)= 0.517255D-01
       PB1(39)=-0.153028D+02
       PB1(40)=-0.348932D-02
       PB1(41)= 0.961771D-03
       PB1(42)= 0.557732D-02
       PB1(43)=-0.454180D-03
       PB1(44)= 0.988213D-01
       PB1(45)= 0.940456D-01
       PB1(46)=-0.318797D-01
       PB1(47)= 0.000000D+00
       PB1(48)= 0.000000D+00
       PB1(49)= 0.000000D+00
       PB1(50)= 0.232122D-02
       PB2( 1)=-0.600220D-02
       PB2( 2)= 0.277654D-04
       PB2( 3)=-0.322019D-02
       PB2( 4)= 0.000000D+00
       PB2( 5)=-0.378551D-02
       PB2( 6)=-0.334809D-02
       PB2( 7)=-0.170668D-02
       PB2( 8)= 0.000000D+00
       PB2( 9)= 0.636184D+04
       PB2(10)= 0.000000D+00
       PB2(11)= 0.159986D-02
       PB2(12)=-0.388204D-02
       PB2(13)=-0.164825D-02
       PB2(14)=-0.747955D+02
       PB2(15)=-0.105360D-01
       PB2(16)=-0.945723D-02
       PB2(17)=-0.159824D-02
       PB2(18)=-0.706730D-03
       PB2(19)=-0.168513D-01
       PB2(20)=-0.113023D+00
       PB2(21)=-0.636637D-01
       PB2(22)=-0.137709D+05
       PB2(23)= 0.000000D+00
       PB2(24)= 0.000000D+00
       PB2(25)= 0.000000D+00
       PB2(26)=-0.152368D+05
       PB2(27)=-0.586061D-02
       PB2(28)=-0.253108D-02
       PB2(29)= 0.000000D+00
       PB2(30)=-0.254837D+04
       PB2(31)=-0.328988D-02
       PB2(32)= 0.000000D+00
       PB2(33)= 0.000000D+00
       PB2(34)=-0.276364D-02
       PB2(35)= 0.967923D-02
       PB2(36)= 0.000000D+00
       PB2(37)= 0.000000D+00
       PB2(38)= 0.000000D+00
       PB2(39)= 0.000000D+00
       PB2(40)= 0.000000D+00
       PB2(41)= 0.434255D-01
       PB2(42)= 0.114020D-01
       PB2(43)=-0.618447D-02
       PB2(44)= 0.000000D+00
       PB2(45)= 0.000000D+00
       PB2(46)=-0.302568D+00
       PB2(47)=-0.327694D-01
       PB2(48)=-0.671589D+02
       PB2(49)=-0.228340D-02
       PB2(50)= 0.000000D+00
       PB3( 1)= 0.306230D-02
       PB3( 2)=-0.465113D-02
       PB3( 3)=-0.973421D-02
       PB3( 4)= 0.128326D-01
       PB3( 5)= 0.788553D-02
       PB3( 6)= 0.797197D-02
       PB3( 7)=-0.120760D-01
       PB3( 8)=-0.767547D-02
       PB3( 9)=-0.120755D-02
       PB3(10)=-0.298523D-01
       PB3(11)=-0.126560D-01
       PB3(12)= 0.000000D+00
       PB3(13)=-0.568350D-01
       PB3(14)=-0.153039D-01
       PB3(15)= 0.000000D+00
       PB3(16)= 0.000000D+00
       PB3(17)= 0.000000D+00
       PB3(18)= 0.000000D+00
       PB3(19)= 0.000000D+00
       PB3(20)= 0.000000D+00
       PB3(21)= 0.000000D+00
       PB3(22)= 0.242911D-02
       PB3(23)=-0.401347D-02
       PB3(24)=-0.219074D-02
       PB3(25)= 0.311281D+01
       PB3(26)= 0.323251D-02
       PB3(27)=-0.639523D-02
       PB3(28)=-0.663069D-02
       PB3(29)=-0.304403D-03
       PB3(30)=-0.401920D-02
       PB3(31)=-0.118708D-02
       PB3(32)= 0.415211D+01
       PB3(33)=-0.201896D+00
       PB3(34)= 0.000000D+00
       PB3(35)= 0.000000D+00
       PB3(36)= 0.000000D+00
       PB3(37)= 0.000000D+00
       PB3(38)= 0.000000D+00
       PB3(39)= 0.000000D+00
       PB3(40)= 0.000000D+00
       PB3(41)= 0.000000D+00
       PB3(42)= 0.000000D+00
       PB3(43)= 0.000000D+00
       PB3(44)= 0.000000D+00
       PB3(45)= 0.000000D+00
       PB3(46)= 0.000000D+00
       PB3(47)= 0.000000D+00
       PB3(48)= 0.000000D+00
       PB3(49)= 0.000000D+00
       PB3(50)= 0.000000D+00
       PC1( 1)= 0.106903D+01
       PC1( 2)= 0.377113D-03
       PC1( 3)= 0.000000D+00
       PC1( 4)= 0.000000D+00
       PC1( 5)= 0.000000D+00
       PC1( 6)= 0.898481D-01
       PC1( 7)=-0.236325D+02
       PC1( 8)= 0.208180D-01
       PC1( 9)= 0.139638D+03
       PC1(10)=-0.119444D+00
       PC1(11)=-0.845398D+01
       PC1(12)=-0.399776D-05
       PC1(13)= 0.000000D+00
       PC1(14)= 0.366210D-02
       PC1(15)=-0.178929D-02
       PC1(16)= 0.190412D-01
       PC1(17)=-0.392257D-01
       PC1(18)= 0.632343D-02
       PC1(19)= 0.548144D-02
       PC1(20)= 0.000000D+00
       PC1(21)= 0.000000D+00
       PC1(22)= 0.000000D+00
       PC1(23)= 0.000000D+00
       PC1(24)= 0.000000D+00
       PC1(25)=-0.243022D-02
       PC1(26)= 0.976619D+00
       PC1(27)= 0.568478D-03
       PC1(28)= 0.582026D-02
       PC1(29)= 0.000000D+00
       PC1(30)= 0.621998D-02
       PC1(31)= 0.000000D+00
       PC1(32)= 0.000000D+00
       PC1(33)= 0.107674D-01
       PC1(34)= 0.893820D+02
       PC1(35)=-0.192414D-01
       PC1(36)=-0.845398D+01
       PC1(37)= 0.000000D+00
       PC1(38)= 0.000000D+00
       PC1(39)=-0.200200D-01
       PC1(40)=-0.195833D-02
       PC1(41)=-0.938391D-02
       PC1(42)= 0.131480D-01
       PC1(43)=-0.260147D-02
       PC1(44)=-0.808556D-03
       PC1(45)= 0.511651D-04
       PC1(46)= 0.255717D-02
       PC1(47)= 0.000000D+00
       PC1(48)= 0.466814D-02
       PC1(49)= 0.664196D-02
       PC1(50)= 0.000000D+00
       PC2( 1)= 0.998594D+00
       PC2( 2)= 0.190038D-03
       PC2( 3)= 0.000000D+00
       PC2( 4)=-0.243825D-01
       PC2( 5)= 0.000000D+00
       PC2( 6)= 0.000000D+00
       PC2( 7)= 0.000000D+00
       PC2( 8)= 0.000000D+00
       PC2( 9)= 0.000000D+00
       PC2(10)= 0.522105D-01
       PC2(11)=-0.845398D+01
       PC2(12)= 0.000000D+00
       PC2(13)= 0.000000D+00
       PC2(14)= 0.000000D+00
       PC2(15)= 0.000000D+00
       PC2(16)= 0.767271D-02
       PC2(17)= 0.564539D-02
       PC2(18)=-0.270623D-02
       PC2(19)=-0.526454D-03
       PC2(20)= 0.137075D-02
       PC2(21)= 0.133060D-02
       PC2(22)= 0.000000D+00
       PC2(23)= 0.000000D+00
       PC2(24)= 0.000000D+00
       PC2(25)= 0.000000D+00
       PC2(26)= 0.949197D+00
       PC2(27)= 0.000000D+00
       PC2(28)= 0.000000D+00
       PC2(29)=-0.768008D-01
       PC2(30)= 0.000000D+00
       PC2(31)= 0.000000D+00
       PC2(32)= 0.000000D+00
       PC2(33)=-0.137993D-01
       PC2(34)=-0.140136D+01
       PC2(35)= 0.120481D+00
       PC2(36)=-0.845398D+01
       PC2(37)= 0.000000D+00
       PC2(38)= 0.000000D+00
       PC2(39)= 0.000000D+00
       PC2(40)= 0.000000D+00
       PC2(41)= 0.987746D-02
       PC2(42)= 0.175330D-02
       PC2(43)=-0.688835D-03
       PC2(44)= 0.287022D-02
       PC2(45)= 0.000000D+00
       PC2(46)= 0.000000D+00
       PC2(47)= 0.744513D-01
       PC2(48)= 0.000000D+00
       PC2(49)= 0.000000D+00
       PC2(50)= 0.000000D+00
       PC3( 1)= 0.152840D+00
       PC3( 2)= 0.000000D+00
       PC3( 3)= 0.000000D+00
       PC3( 4)= 0.116252D+01
       PC3( 5)= 0.000000D+00
       PC3( 6)= 0.000000D+00
       PC3( 7)= 0.000000D+00
       PC3( 8)= 0.000000D+00
       PC3( 9)= 0.000000D+00
       PC3(10)=-0.649190D+00
       PC3(11)=-0.845398D+01
       PC3(12)= 0.000000D+00
       PC3(13)= 0.000000D+00
       PC3(14)= 0.000000D+00
       PC3(15)= 0.000000D+00
       PC3(16)=-0.584949D-01
       PC3(17)=-0.102105D+00
       PC3(18)= 0.299153D-01
       PC3(19)=-0.486227D-01
       PC3(20)= 0.000000D+00
       PC3(21)= 0.000000D+00
       PC3(22)= 0.000000D+00
       PC3(23)= 0.000000D+00
       PC3(24)= 0.000000D+00
       PC3(25)= 0.000000D+00
       PC3(26)= 0.000000D+00
       PC3(27)= 0.000000D+00
       PC3(28)= 0.000000D+00
       PC3(29)= 0.000000D+00
       PC3(30)= 0.000000D+00
       PC3(31)= 0.000000D+00
       PC3(32)= 0.000000D+00
       PC3(33)= 0.000000D+00
       PC3(34)= 0.000000D+00
       PC3(35)= 0.000000D+00
       PC3(36)= 0.000000D+00
       PC3(37)= 0.000000D+00
       PC3(38)= 0.000000D+00
       PC3(39)= 0.000000D+00
       PC3(40)= 0.000000D+00
       PC3(41)= 0.000000D+00
       PC3(42)= 0.000000D+00
       PC3(43)= 0.000000D+00
       PC3(44)= 0.000000D+00
       PC3(45)= 0.000000D+00
       PC3(46)= 0.000000D+00
       PC3(47)= 0.000000D+00
       PC3(48)= 0.000000D+00
       PC3(49)= 0.000000D+00
       PC3(50)= 0.000000D+00
       PD1( 1)= 0.931402D+00
       PD1( 2)= 0.137976D+00
       PD1( 3)= 0.000000D+00
       PD1( 4)= 0.323736D-03
       PD1( 5)= 0.000000D+00
       PD1( 6)=-0.910906D-02
       PD1( 7)= 0.707506D-01
       PD1( 8)= 0.000000D+00
       PD1( 9)=-0.516650D-01
       PD1(10)= 0.689755D-01
       PD1(11)= 0.000000D+00
       PD1(12)= 0.000000D+00
       PD1(13)= 0.000000D+00
       PD1(14)=-0.845398D+01
       PD1(15)= 0.000000D+00
       PD1(16)= 0.281140D-01
       PD1(17)= 0.000000D+00
       PD1(18)= 0.736009D+02
       PD1(19)= 0.596604D-01
       PD1(20)= 0.000000D+00
       PD1(21)= 0.000000D+00
       PD1(22)=-0.151792D-02
       PD1(23)= 0.000000D+00
       PD1(24)= 0.000000D+00
       PD1(25)= 0.132397D+00
       PD1(26)= 0.213315D+00
       PD1(27)= 0.000000D+00
       PD1(28)= 0.000000D+00
       PD1(29)= 0.000000D+00
       PD1(30)= 0.000000D+00
       PD1(31)= 0.000000D+00
       PD1(32)= 0.948758D+01
       PD1(33)= 0.884541D-02
       PD1(34)= 0.000000D+00
       PD1(35)= 0.000000D+00
       PD1(36)= 0.000000D+00
       PD1(37)= 0.000000D+00
       PD1(38)= 0.000000D+00
       PD1(39)= 0.000000D+00
       PD1(40)= 0.000000D+00
       PD1(41)= 0.000000D+00
       PD1(42)= 0.000000D+00
       PD1(43)= 0.000000D+00
       PD1(44)= 0.113139D+00
       PD1(45)= 0.169134D+00
       PD1(46)= 0.145192D-01
       PD1(47)= 0.000000D+00
       PD1(48)= 0.000000D+00
       PD1(49)= 0.000000D+00
       PD1(50)= 0.000000D+00
       PD2( 1)= 0.107906D-01
       PD2( 2)= 0.299942D-04
       PD2( 3)= 0.000000D+00
       PD2( 4)= 0.000000D+00
       PD2( 5)= 0.000000D+00
       PD2( 6)= 0.000000D+00
       PD2( 7)= 0.000000D+00
       PD2( 8)= 0.000000D+00
       PD2( 9)= 0.000000D+00
       PD2(10)= 0.000000D+00
       PD2(11)= 0.000000D+00
       PD2(12)= 0.000000D+00
       PD2(13)= 0.000000D+00
       PD2(14)= 0.000000D+00
       PD2(15)=-0.148930D-01
       PD2(16)=-0.787184D-02
       PD2(17)= 0.000000D+00
       PD2(18)= 0.000000D+00
       PD2(19)= 0.000000D+00
       PD2(20)= 0.000000D+00
       PD2(21)= 0.000000D+00
       PD2(22)= 0.000000D+00
       PD2(23)= 0.000000D+00
       PD2(24)= 0.000000D+00
       PD2(25)= 0.000000D+00
       PD2(26)= 0.000000D+00
       PD2(27)= 0.000000D+00
       PD2(28)= 0.000000D+00
       PD2(29)= 0.000000D+00
       PD2(30)= 0.000000D+00
       PD2(31)= 0.000000D+00
       PD2(32)= 0.000000D+00
       PD2(33)= 0.000000D+00
       PD2(34)= 0.000000D+00
       PD2(35)= 0.000000D+00
       PD2(36)= 0.000000D+00
       PD2(37)= 0.000000D+00
       PD2(38)= 0.000000D+00
       PD2(39)= 0.000000D+00
       PD2(40)= 0.000000D+00
       PD2(41)=-0.683420D-01
       PD2(42)=-0.441778D-01
       PD2(43)= 0.000000D+00
       PD2(44)= 0.000000D+00
       PD2(45)= 0.000000D+00
       PD2(46)= 0.000000D+00
       PD2(47)= 0.229730D-01
       PD2(48)= 0.000000D+00
       PD2(49)= 0.000000D+00
       PD2(50)= 0.000000D+00
       PD3( 1)= 0.000000D+00
       PD3( 2)= 0.000000D+00
       PD3( 3)= 0.000000D+00
       PD3( 4)= 0.000000D+00
       PD3( 5)= 0.000000D+00
       PD3( 6)= 0.000000D+00
       PD3( 7)= 0.000000D+00
       PD3( 8)= 0.000000D+00
       PD3( 9)= 0.000000D+00
       PD3(10)= 0.000000D+00
       PD3(11)= 0.000000D+00
       PD3(12)= 0.000000D+00
       PD3(13)= 0.000000D+00
       PD3(14)= 0.000000D+00
       PD3(15)= 0.000000D+00
       PD3(16)= 0.000000D+00
       PD3(17)= 0.000000D+00
       PD3(18)= 0.000000D+00
       PD3(19)= 0.000000D+00
       PD3(20)= 0.000000D+00
       PD3(21)= 0.000000D+00
       PD3(22)= 0.000000D+00
       PD3(23)= 0.000000D+00
       PD3(24)= 0.000000D+00
       PD3(25)= 0.000000D+00
       PD3(26)= 0.000000D+00
       PD3(27)= 0.000000D+00
       PD3(28)= 0.000000D+00
       PD3(29)= 0.000000D+00
       PD3(30)= 0.000000D+00
       PD3(31)= 0.000000D+00
       PD3(32)= 0.000000D+00
       PD3(33)= 0.000000D+00
       PD3(34)= 0.000000D+00
       PD3(35)= 0.000000D+00
       PD3(36)= 0.000000D+00
       PD3(37)= 0.000000D+00
       PD3(38)= 0.000000D+00
       PD3(39)= 0.000000D+00
       PD3(40)= 0.000000D+00
       PD3(41)= 0.000000D+00
       PD3(42)= 0.000000D+00
       PD3(43)= 0.000000D+00
       PD3(44)= 0.000000D+00
       PD3(45)= 0.000000D+00
       PD3(46)= 0.000000D+00
       PD3(47)= 0.000000D+00
       PD3(48)= 0.000000D+00
       PD3(49)= 0.000000D+00
       PD3(50)= 0.000000D+00
       PE1( 1)= 0.868053D+00
       PE1( 2)= 0.236364D+00
       PE1( 3)= 0.134306D+00
       PE1( 4)= 0.103086D-01
       PE1( 5)= 0.000000D+00
       PE1( 6)=-0.379164D-02
       PE1( 7)=-0.157806D+00
       PE1( 8)= 0.000000D+00
       PE1( 9)=-0.587644D-01
       PE1(10)=-0.312508D+00
       PE1(11)= 0.000000D+00
       PE1(12)= 0.437387D-01
       PE1(13)=-0.354091D-01
       PE1(14)=-0.223636D+02
       PE1(15)= 0.000000D+00
       PE1(16)=-0.533976D-01
       PE1(17)= 0.000000D+00
       PE1(18)= 0.114091D+03
       PE1(19)= 0.517497D-01
       PE1(20)= 0.000000D+00
       PE1(21)= 0.000000D+00
       PE1(22)= 0.000000D+00
       PE1(23)= 0.000000D+00
       PE1(24)= 0.000000D+00
       PE1(25)= 0.132397D+00
       PE1(26)= 0.213315D+00
       PE1(27)= 0.000000D+00
       PE1(28)= 0.000000D+00
       PE1(29)= 0.000000D+00
       PE1(30)= 0.000000D+00
       PE1(31)= 0.000000D+00
       PE1(32)= 0.342702D+03
       PE1(33)= 0.157033D-01
       PE1(34)= 0.000000D+00
       PE1(35)= 0.000000D+00
       PE1(36)= 0.000000D+00
       PE1(37)= 0.000000D+00
       PE1(38)= 0.000000D+00
       PE1(39)= 0.000000D+00
       PE1(40)=-0.366278D-02
       PE1(41)=-0.116193D-02
       PE1(42)= 0.000000D+00
       PE1(43)= 0.000000D+00
       PE1(44)= 0.113139D+00
       PE1(45)= 0.169134D+00
       PE1(46)= 0.178431D-01
       PE1(47)= 0.000000D+00
       PE1(48)= 0.000000D+00
       PE1(49)= 0.000000D+00
       PE1(50)= 0.000000D+00
       PE2( 1)= 0.162864D-01
       PE2( 2)= 0.316963D-04
       PE2( 3)= 0.127968D-01
       PE2( 4)= 0.000000D+00
       PE2( 5)= 0.000000D+00
       PE2( 6)=-0.704599D-02
       PE2( 7)= 0.207921D-02
       PE2( 8)= 0.636660D-02
       PE2( 9)= 0.229940D+05
       PE2(10)= 0.000000D+00
       PE2(11)= 0.127833D-01
       PE2(12)=-0.208036D-02
       PE2(13)=-0.461820D-02
       PE2(14)=-0.629391D+02
       PE2(15)=-0.120745D-01
       PE2(16)= 0.136675D-01
       PE2(17)= 0.136011D-01
       PE2(18)=-0.537162D-02
       PE2(19)= 0.000000D+00
       PE2(20)= 0.000000D+00
       PE2(21)= 0.000000D+00
       PE2(22)= 0.000000D+00
       PE2(23)= 0.000000D+00
       PE2(24)= 0.000000D+00
       PE2(25)= 0.000000D+00
       PE2(26)= 0.192509D+05
       PE2(27)= 0.835522D-02
       PE2(28)= 0.419439D-02
       PE2(29)= 0.000000D+00
       PE2(30)= 0.120366D+05
       PE2(31)= 0.000000D+00
       PE2(32)= 0.000000D+00
       PE2(33)= 0.000000D+00
       PE2(34)=-0.100034D-01
       PE2(35)=-0.233267D-02
       PE2(36)= 0.972374D-02
       PE2(37)= 0.000000D+00
       PE2(38)= 0.000000D+00
       PE2(39)= 0.000000D+00
       PE2(40)= 0.000000D+00
       PE2(41)=-0.265079D-01
       PE2(42)=-0.209125D-01
       PE2(43)=-0.109465D-01
       PE2(44)= 0.000000D+00
       PE2(45)= 0.000000D+00
       PE2(46)= 0.000000D+00
       PE2(47)= 0.217252D-01
       PE2(48)=-0.712385D+02
       PE2(49)=-0.189428D-02
       PE2(50)= 0.000000D+00
       PE3( 1)=-0.602006D-02
       PE3( 2)= 0.169058D-01
       PE3( 3)= 0.000000D+00
       PE3( 4)= 0.000000D+00
       PE3( 5)= 0.000000D+00
       PE3( 6)= 0.000000D+00
       PE3( 7)= 0.000000D+00
       PE3( 8)= 0.000000D+00
       PE3( 9)= 0.000000D+00
       PE3(10)= 0.290646D-01
       PE3(11)= 0.348971D-02
       PE3(12)= 0.000000D+00
       PE3(13)= 0.501174D-01
       PE3(14)= 0.550595D-01
       PE3(15)= 0.000000D+00
       PE3(16)=-0.955897D-02
       PE3(17)= 0.000000D+00
       PE3(18)= 0.000000D+00
       PE3(19)=-0.151693D+04
       PE3(20)= 0.000000D+00
       PE3(21)= 0.000000D+00
       PE3(22)= 0.129306D-01
       PE3(23)= 0.269567D-02
       PE3(24)= 0.000000D+00
       PE3(25)= 0.392243D+01
       PE3(26)=-0.847690D-02
       PE3(27)= 0.116896D-01
       PE3(28)= 0.000000D+00
       PE3(29)= 0.148967D-01
       PE3(30)= 0.544521D-02
       PE3(31)= 0.000000D+00
       PE3(32)= 0.564918D+01
       PE3(33)= 0.000000D+00
       PE3(34)=-0.772178D-02
       PE3(35)= 0.000000D+00
       PE3(36)= 0.000000D+00
       PE3(37)=-0.734042D+02
       PE3(38)= 0.000000D+00
       PE3(39)= 0.000000D+00
       PE3(40)= 0.000000D+00
       PE3(41)= 0.000000D+00
       PE3(42)= 0.000000D+00
       PE3(43)= 0.000000D+00
       PE3(44)= 0.000000D+00
       PE3(45)= 0.000000D+00
       PE3(46)= 0.000000D+00
       PE3(47)= 0.000000D+00
       PE3(48)= 0.000000D+00
       PE3(49)= 0.000000D+00
       PE3(50)= 0.000000D+00
       PF1( 1)= 0.127515D+01
       PF1( 2)=-0.210472D+00
       PF1( 3)=-0.177924D+00
       PF1( 4)= 0.218900D+00
       PF1( 5)= 0.288436D-01
       PF1( 6)= 0.190077D-01
       PF1( 7)= 0.291001D+00
       PF1( 8)= 0.217437D-01
       PF1( 9)=-0.105186D-01
       PF1(10)= 0.436141D+00
       PF1(11)= 0.107605D+00
       PF1(12)= 0.330755D-01
       PF1(13)= 0.400581D-01
       PF1(14)=-0.958051D+01
       PF1(15)= 0.000000D+00
       PF1(16)= 0.154028D-01
       PF1(17)= 0.000000D+00
       PF1(18)= 0.734194D+02
       PF1(19)= 0.496540D-01
       PF1(20)=-0.595906D-02
       PF1(21)= 0.384512D-04
       PF1(22)=-0.136000D-01
       PF1(23)= 0.000000D+00
       PF1(24)= 0.000000D+00
       PF1(25)= 0.132397D+00
       PF1(26)= 0.213315D+00
       PF1(27)=-0.416610D-01
       PF1(28)= 0.000000D+00
       PF1(29)= 0.000000D+00
       PF1(30)= 0.000000D+00
       PF1(31)= 0.000000D+00
       PF1(32)= 0.146276D+03
       PF1(33)=-0.198408D-01
       PF1(34)= 0.000000D+00
       PF1(35)= 0.132530D-01
       PF1(36)= 0.000000D+00
       PF1(37)= 0.000000D+00
       PF1(38)= 0.000000D+00
       PF1(39)= 0.000000D+00
       PF1(40)=-0.104687D-03
       PF1(41)=-0.147562D-02
       PF1(42)= 0.000000D+00
       PF1(43)= 0.000000D+00
       PF1(44)= 0.113139D+00
       PF1(45)= 0.169134D+00
       PF1(46)=-0.126913D-01
       PF1(47)= 0.000000D+00
       PF1(48)= 0.000000D+00
       PF1(49)= 0.000000D+00
       PF1(50)=-0.608370D-02
       PF2( 1)=-0.257587D-01
       PF2( 2)= 0.319022D-04
       PF2( 3)= 0.000000D+00
       PF2( 4)= 0.000000D+00
       PF2( 5)= 0.156644D-01
       PF2( 6)= 0.103640D-01
       PF2( 7)= 0.105771D-02
       PF2( 8)= 0.000000D+00
       PF2( 9)= 0.357949D+04
       PF2(10)= 0.000000D+00
       PF2(11)=-0.125672D-02
       PF2(12)= 0.152783D-02
       PF2(13)= 0.130518D-02
       PF2(14)= 0.755558D+01
       PF2(15)=-0.920341D-02
       PF2(16)=-0.209142D-01
       PF2(17)=-0.134106D-01
       PF2(18)= 0.000000D+00
       PF2(19)=-0.483312D-01
       PF2(20)= 0.830900D-01
       PF2(21)= 0.988009D-01
       PF2(22)=-0.141148D+05
       PF2(23)= 0.000000D+00
       PF2(24)= 0.000000D+00
       PF2(25)= 0.000000D+00
       PF2(26)=-0.105513D+04
       PF2(27)= 0.000000D+00
       PF2(28)= 0.000000D+00
       PF2(29)= 0.000000D+00
       PF2(30)= 0.000000D+00
       PF2(31)= 0.000000D+00
       PF2(32)= 0.000000D+00
       PF2(33)= 0.000000D+00
       PF2(34)= 0.673442D-02
       PF2(35)= 0.201691D-02
       PF2(36)= 0.000000D+00
       PF2(37)= 0.000000D+00
       PF2(38)= 0.000000D+00
       PF2(39)= 0.000000D+00
       PF2(40)= 0.000000D+00
       PF2(41)= 0.598019D-01
       PF2(42)= 0.633298D-02
       PF2(43)=-0.112871D-02
       PF2(44)= 0.000000D+00
       PF2(45)= 0.000000D+00
       PF2(46)= 0.000000D+00
       PF2(47)=-0.128604D-01
       PF2(48)= 0.000000D+00
       PF2(49)= 0.000000D+00
       PF2(50)= 0.000000D+00
       PF3( 1)=-0.494960D-02
       PF3( 2)=-0.136415D-01
       PF3( 3)=-0.115039D-01
       PF3( 4)= 0.000000D+00
       PF3( 5)= 0.000000D+00
       PF3( 6)= 0.000000D+00
       PF3( 7)= 0.000000D+00
       PF3( 8)= 0.000000D+00
       PF3( 9)= 0.000000D+00
       PF3(10)= 0.000000D+00
       PF3(11)= 0.000000D+00
       PF3(12)= 0.000000D+00
       PF3(13)= 0.000000D+00
       PF3(14)= 0.000000D+00
       PF3(15)= 0.000000D+00
       PF3(16)= 0.000000D+00
       PF3(17)= 0.000000D+00
       PF3(18)= 0.000000D+00
       PF3(19)= 0.000000D+00
       PF3(20)= 0.000000D+00
       PF3(21)= 0.000000D+00
       PF3(22)=-0.586860D-02
       PF3(23)=-0.141732D-02
       PF3(24)= 0.213697D-02
       PF3(25)= 0.263845D+01
       PF3(26)=-0.834186D-02
       PF3(27)=-0.187336D-01
       PF3(28)=-0.190870D-01
       PF3(29)=-0.803810D-02
       PF3(30)=-0.284279D-02
       PF3(31)= 0.256722D-02
       PF3(32)= 0.171429D+01
       PF3(33)= 0.000000D+00
       PF3(34)= 0.000000D+00
       PF3(35)= 0.000000D+00
       PF3(36)= 0.000000D+00
       PF3(37)= 0.000000D+00
       PF3(38)= 0.000000D+00
       PF3(39)= 0.000000D+00
       PF3(40)= 0.000000D+00
       PF3(41)= 0.000000D+00
       PF3(42)= 0.000000D+00
       PF3(43)= 0.000000D+00
       PF3(44)= 0.000000D+00
       PF3(45)= 0.000000D+00
       PF3(46)= 0.000000D+00
       PF3(47)= 0.000000D+00
       PF3(48)= 0.000000D+00
       PF3(49)= 0.000000D+00
       PF3(50)= 0.000000D+00
       PG1( 1)= 0.573587D+02
       PG1( 2)=-0.398747D+00
       PG1( 3)= 0.000000D+00
       PG1( 4)=-0.529554D+00
       PG1( 5)=-0.582186D-02
       PG1( 6)= 0.714177D-01
       PG1( 7)=-0.679279D+00
       PG1( 8)=-0.167715D+00
       PG1( 9)=-0.642434D-01
       PG1(10)=-0.211569D+00
       PG1(11)=-0.159922D+00
       PG1(12)=-0.171024D-03
       PG1(13)=-0.115885D+00
       PG1(14)= 0.651603D+01
       PG1(15)= 0.000000D+00
       PG1(16)=-0.176683D+00
       PG1(17)= 0.650395D-01
       PG1(18)= 0.143504D+01
       PG1(19)= 0.928208D-01
       PG1(20)= 0.511662D-02
       PG1(21)= 0.000000D+00
       PG1(22)= 0.995121D-02
       PG1(23)= 0.000000D+00
       PG1(24)= 0.000000D+00
       PG1(25)= 0.132397D+00
       PG1(26)= 0.213315D+00
       PG1(27)= 0.101451D+00
       PG1(28)= 0.000000D+00
       PG1(29)= 0.000000D+00
       PG1(30)= 0.000000D+00
       PG1(31)= 0.000000D+00
       PG1(32)= 0.567667D+02
       PG1(33)= 0.238192D-02
       PG1(34)= 0.000000D+00
       PG1(35)=-0.188240D-01
       PG1(36)= 0.000000D+00
       PG1(37)= 0.000000D+00
       PG1(38)= 0.476218D-01
       PG1(39)= 0.235206D+02
       PG1(40)= 0.475901D-02
       PG1(41)= 0.576162D-02
       PG1(42)= 0.151815D-01
       PG1(43)=-0.192730D-01
       PG1(44)= 0.113139D+00
       PG1(45)= 0.169134D+00
       PG1(46)=-0.288771D-01
       PG1(47)= 0.000000D+00
       PG1(48)= 0.000000D+00
       PG1(49)= 0.000000D+00
       PG1(50)= 0.118418D-02
       PG2( 1)=-0.368927D-02
       PG2( 2)= 0.314704D-04
       PG2( 3)= 0.882198D-02
       PG2( 4)= 0.000000D+00
       PG2( 5)=-0.192562D-01
       PG2( 6)=-0.258674D-02
       PG2( 7)=-0.219913D-01
       PG2( 8)= 0.000000D+00
       PG2( 9)= 0.438655D+04
       PG2(10)= 0.000000D+00
       PG2(11)= 0.760126D-02
       PG2(12)= 0.259438D-02
       PG2(13)= 0.172310D-02
       PG2(14)= 0.779204D+02
       PG2(15)= 0.797786D-03
       PG2(16)=-0.770510D-02
       PG2(17)= 0.190982D-02
       PG2(18)= 0.272707D-02
       PG2(19)= 0.101016D-01
       PG2(20)= 0.116537D+00
       PG2(21)=-0.312236D-02
       PG2(22)= 0.139783D+05
       PG2(23)= 0.000000D+00
       PG2(24)= 0.000000D+00
       PG2(25)= 0.000000D+00
       PG2(26)=-0.130712D+04
       PG2(27)= 0.000000D+00
       PG2(28)= 0.000000D+00
       PG2(29)= 0.000000D+00
       PG2(30)= 0.000000D+00
       PG2(31)= 0.000000D+00
       PG2(32)= 0.000000D+00
       PG2(33)= 0.000000D+00
       PG2(34)=-0.320544D-02
       PG2(35)=-0.206970D-01
       PG2(36)= 0.000000D+00
       PG2(37)= 0.000000D+00
       PG2(38)= 0.000000D+00
       PG2(39)= 0.000000D+00
       PG2(40)= 0.000000D+00
       PG2(41)= 0.159010D-01
       PG2(42)=-0.191427D-02
       PG2(43)=-0.342829D-01
       PG2(44)= 0.000000D+00
       PG2(45)= 0.000000D+00
       PG2(46)= 0.000000D+00
       PG2(47)=-0.345379D-01
       PG2(48)= 0.894518D+02
       PG2(49)= 0.171556D-02
       PG2(50)= 0.000000D+00
       PG3( 1)=-0.765278D-02
       PG3( 2)=-0.208987D-03
       PG3( 3)=-0.157393D-01
       PG3( 4)= 0.000000D+00
       PG3( 5)= 0.000000D+00
       PG3( 6)= 0.000000D+00
       PG3( 7)= 0.000000D+00
       PG3( 8)= 0.000000D+00
       PG3( 9)= 0.000000D+00
       PG3(10)= 0.000000D+00
       PG3(11)= 0.000000D+00
       PG3(12)= 0.000000D+00
       PG3(13)= 0.000000D+00
       PG3(14)= 0.000000D+00
       PG3(15)= 0.000000D+00
       PG3(16)= 0.000000D+00
       PG3(17)= 0.000000D+00
       PG3(18)= 0.000000D+00
       PG3(19)= 0.000000D+00
       PG3(20)= 0.000000D+00
       PG3(21)= 0.000000D+00
       PG3(22)=-0.860673D-02
       PG3(23)=-0.119922D-01
       PG3(24)=-0.646356D-02
       PG3(25)=-0.300107D+01
       PG3(26)=-0.932511D-02
       PG3(27)=-0.150205D-01
       PG3(28)=-0.867835D-02
       PG3(29)=-0.764801D-02
       PG3(30)=-0.131495D-01
       PG3(31)=-0.676720D-02
       PG3(32)=-0.182396D+01
       PG3(33)= 0.000000D+00
       PG3(34)= 0.000000D+00
       PG3(35)= 0.000000D+00
       PG3(36)= 0.000000D+00
       PG3(37)= 0.000000D+00
       PG3(38)= 0.000000D+00
       PG3(39)= 0.000000D+00
       PG3(40)= 0.000000D+00
       PG3(41)= 0.000000D+00
       PG3(42)= 0.000000D+00
       PG3(43)= 0.000000D+00
       PG3(44)= 0.000000D+00
       PG3(45)= 0.000000D+00
       PG3(46)= 0.000000D+00
       PG3(47)= 0.000000D+00
       PG3(48)= 0.000000D+00
       PG3(49)= 0.000000D+00
       PG3(50)= 0.000000D+00
       PH1( 1)= 0.951363D+00
       PH1( 2)=-0.467542D-01
       PH1( 3)= 0.120260D+00
       PH1( 4)= 0.000000D+00
       PH1( 5)= 0.000000D+00
       PH1( 6)= 0.191357D-01
       PH1( 7)= 0.000000D+00
       PH1( 8)= 0.000000D+00
       PH1( 9)= 0.125429D-02
       PH1(10)=-0.133240D+00
       PH1(11)= 0.000000D+00
       PH1(12)= 0.000000D+00
       PH1(13)= 0.000000D+00
       PH1(14)=-0.845398D+01
       PH1(15)= 0.000000D+00
       PH1(16)= 0.000000D+00
       PH1(17)= 0.000000D+00
       PH1(18)= 0.000000D+00
       PH1(19)= 0.000000D+00
       PH1(20)= 0.000000D+00
       PH1(21)= 0.000000D+00
       PH1(22)= 0.252317D-02
       PH1(23)= 0.000000D+00
       PH1(24)=-0.973404D-02
       PH1(25)= 0.132397D+00
       PH1(26)= 0.213315D+00
       PH1(27)= 0.000000D+00
       PH1(28)= 0.000000D+00
       PH1(29)= 0.000000D+00
       PH1(30)= 0.000000D+00
       PH1(31)= 0.000000D+00
       PH1(32)= 0.000000D+00
       PH1(33)= 0.000000D+00
       PH1(34)=-0.718482D-03
       PH1(35)= 0.000000D+00
       PH1(36)= 0.000000D+00
       PH1(37)= 0.000000D+00
       PH1(38)= 0.000000D+00
       PH1(39)= 0.000000D+00
       PH1(40)= 0.000000D+00
       PH1(41)= 0.000000D+00
       PH1(42)= 0.787683D-02
       PH1(43)=-0.233698D-02
       PH1(44)= 0.113139D+00
       PH1(45)= 0.169134D+00
       PH1(46)= 0.000000D+00
       PH1(47)= 0.000000D+00
       PH1(48)= 0.000000D+00
       PH1(49)= 0.000000D+00
       PH1(50)= 0.000000D+00
       PH2( 1)= 0.000000D+00
       PH2( 2)= 0.000000D+00
       PH2( 3)= 0.000000D+00
       PH2( 4)= 0.000000D+00
       PH2( 5)= 0.000000D+00
       PH2( 6)= 0.000000D+00
       PH2( 7)= 0.000000D+00
       PH2( 8)= 0.000000D+00
       PH2( 9)= 0.000000D+00
       PH2(10)= 0.000000D+00
       PH2(11)= 0.000000D+00
       PH2(12)= 0.000000D+00
       PH2(13)= 0.000000D+00
       PH2(14)= 0.000000D+00
       PH2(15)= 0.000000D+00
       PH2(16)= 0.000000D+00
       PH2(17)= 0.000000D+00
       PH2(18)= 0.000000D+00
       PH2(19)= 0.000000D+00
       PH2(20)= 0.000000D+00
       PH2(21)= 0.000000D+00
       PH2(22)= 0.000000D+00
       PH2(23)= 0.000000D+00
       PH2(24)= 0.000000D+00
       PH2(25)= 0.000000D+00
       PH2(26)= 0.000000D+00
       PH2(27)= 0.000000D+00
       PH2(28)= 0.000000D+00
       PH2(29)= 0.000000D+00
       PH2(30)= 0.000000D+00
       PH2(31)= 0.000000D+00
       PH2(32)= 0.000000D+00
       PH2(33)= 0.000000D+00
       PH2(34)= 0.000000D+00
       PH2(35)= 0.000000D+00
       PH2(36)= 0.000000D+00
       PH2(37)= 0.000000D+00
       PH2(38)= 0.000000D+00
       PH2(39)= 0.000000D+00
       PH2(40)= 0.000000D+00
       PH2(41)= 0.000000D+00
       PH2(42)= 0.000000D+00
       PH2(43)= 0.000000D+00
       PH2(44)= 0.000000D+00
       PH2(45)= 0.000000D+00
       PH2(46)= 0.000000D+00
       PH2(47)= 0.000000D+00
       PH2(48)= 0.000000D+00
       PH2(49)= 0.000000D+00
       PH2(50)= 0.000000D+00
       PH3( 1)= 0.000000D+00
       PH3( 2)= 0.000000D+00
       PH3( 3)= 0.000000D+00
       PH3( 4)= 0.000000D+00
       PH3( 5)= 0.000000D+00
       PH3( 6)= 0.000000D+00
       PH3( 7)= 0.000000D+00
       PH3( 8)= 0.000000D+00
       PH3( 9)= 0.000000D+00
       PH3(10)= 0.000000D+00
       PH3(11)= 0.000000D+00
       PH3(12)= 0.000000D+00
       PH3(13)= 0.000000D+00
       PH3(14)= 0.000000D+00
       PH3(15)= 0.000000D+00
       PH3(16)= 0.000000D+00
       PH3(17)= 0.000000D+00
       PH3(18)= 0.000000D+00
       PH3(19)= 0.000000D+00
       PH3(20)= 0.000000D+00
       PH3(21)= 0.000000D+00
       PH3(22)= 0.000000D+00
       PH3(23)= 0.000000D+00
       PH3(24)= 0.000000D+00
       PH3(25)= 0.000000D+00
       PH3(26)= 0.000000D+00
       PH3(27)= 0.000000D+00
       PH3(28)= 0.000000D+00
       PH3(29)= 0.000000D+00
       PH3(30)= 0.000000D+00
       PH3(31)= 0.000000D+00
       PH3(32)= 0.000000D+00
       PH3(33)= 0.000000D+00
       PH3(34)= 0.000000D+00
       PH3(35)= 0.000000D+00
       PH3(36)= 0.000000D+00
       PH3(37)= 0.000000D+00
       PH3(38)= 0.000000D+00
       PH3(39)= 0.000000D+00
       PH3(40)= 0.000000D+00
       PH3(41)= 0.000000D+00
       PH3(42)= 0.000000D+00
       PH3(43)= 0.000000D+00
       PH3(44)= 0.000000D+00
       PH3(45)= 0.000000D+00
       PH3(46)= 0.000000D+00
       PH3(47)= 0.000000D+00
       PH3(48)= 0.000000D+00
       PH3(49)= 0.000000D+00
       PH3(50)= 0.000000D+00
       PI1( 1)= 0.933804D+00
       PI1( 2)= 0.547446D+01
       PI1( 3)= 0.153263D+00
       PI1( 4)= 0.919303D+00
       PI1( 5)= 0.164109D+02
       PI1( 6)= 0.427083D+01
       PI1( 7)= 0.000000D+00
       PI1( 8)= 0.000000D+00
       PI1( 9)= 0.000000D+00
       PI1(10)= 0.000000D+00
       PI1(11)= 0.000000D+00
       PI1(12)= 0.000000D+00
       PI1(13)= 0.000000D+00
       PI1(14)= 0.000000D+00
       PI1(15)= 0.000000D+00
       PI1(16)= 0.000000D+00
       PI1(17)= 0.000000D+00
       PI1(18)= 0.000000D+00
       PI1(19)= 0.000000D+00
       PI1(20)= 0.000000D+00
       PI1(21)= 0.000000D+00
       PI1(22)= 0.000000D+00
       PI1(23)= 0.000000D+00
       PI1(24)= 0.000000D+00
       PI1(25)= 0.000000D+00
       PI1(26)= 0.115897D+01
       PI1(27)= 0.471094D+00
       PI1(28)= 0.109459D+01
       PI1(29)= 0.525012D+01
       PI1(30)= 0.100000D+01
       PI1(31)= 0.100000D+01
       PI1(32)= 0.103999D+01
       PI1(33)= 0.767132D+00
       PI1(34)= 0.110514D+01
       PI1(35)= 0.175636D+01
       PI1(36)= 0.110845D+01
       PI1(37)= 0.233439D+01
       PI1(38)= 0.796532D+00
       PI1(39)= 0.431520D+01
       PI1(40)= 0.407300D+01
       PI1(41)= 0.101885D+01
       PI1(42)= 0.239547D+00
       PI1(43)= 0.253791D-05
       PI1(44)= 0.842931D+00
       PI1(45)= 0.104192D+01
       PI1(46)= 0.200202D+01
       PI1(47)= 0.100000D+01
       PI1(48)= 0.100000D+01
       PI1(49)= 0.100000D+01
       PI1(50)= 0.100000D+01
       PTM( 1)= 0.104130D+04
       PTM( 2)= 0.386000D+03
       PTM( 3)= 0.190000D+03
       PTM( 4)= 0.166728D+02
       PTM( 5)= 0.115000D+03
       PTM( 6)= 0.120000D+03
       PTM( 7)= 0.945537D+02
       PTM( 8)= 0.000000D+00
       PDM( 1, 1)= 0.245600D+08
       PDM( 1, 2)= 0.859400D+11
       PDM( 1, 3)= 0.281000D+12
       PDM( 1, 4)= 0.330000D+11
       PDM( 1, 5)= 0.133000D+10
       PDM( 1, 6)= 0.176100D+06
       PDM( 1, 7)= 0.100000D+08
       PDM( 2, 1)= 0.671072D-05
       PDM( 2, 2)= 0.540000D+00
       PDM( 2, 3)= 0.000000D+00
       PDM( 2, 4)= 0.268270D+00
       PDM( 2, 5)= 0.119615D-01
       PDM( 2, 6)= 0.100000D+01
       PDM( 2, 7)= 0.100000D+01
       PDM( 3, 1)= 0.100000D+03
       PDM( 3, 2)= 0.105000D+03
       PDM( 3, 3)= 0.105000D+03
       PDM( 3, 4)= 0.105000D+03
       PDM( 3, 5)= 0.105000D+03
       PDM( 3, 6)= 0.950000D+02
       PDM( 3, 7)= 0.105000D+03
       PDM( 4, 1)= 0.000000D+00
       PDM( 4, 2)=-0.800000D+01
       PDM( 4, 3)= 0.280000D+02
       PDM( 4, 4)= 0.000000D+00
       PDM( 4, 5)= 0.000000D+00
       PDM( 4, 6)=-0.800000D+01
       PDM( 4, 7)=-0.800000D+01
       PDM( 5, 1)= 0.110000D+03
       PDM( 5, 2)= 0.110000D+03
       PDM( 5, 3)= 0.289500D+02
       PDM( 5, 4)= 0.110000D+03
       PDM( 5, 5)= 0.110000D+03
       PDM( 5, 6)= 0.110000D+03
       PDM( 5, 7)= 0.110000D+03
       PDM( 6, 1)= 0.100000D+02
       PDM( 6, 2)= 0.100000D+02
       PDM( 6, 3)= 0.000000D+00
       PDM( 6, 4)= 0.100000D+02
       PDM( 6, 5)= 0.100000D+02
       PDM( 6, 6)= 0.100000D+02
       PDM( 6, 7)= 0.100000D+02
       PDM( 7, 1)= 0.000000D+00
       PDM( 7, 2)= 0.900000D+02
       PDM( 7, 3)= 0.000000D+00
       PDM( 7, 4)= 0.000000D+00
       PDM( 7, 5)= 0.000000D+00
       PDM( 7, 6)= 0.900000D+02
       PDM( 7, 7)= 0.900000D+02
       PDM( 8, 1)= 0.000000D+00
       PDM( 8, 2)= 0.200000D+01
       PDM( 8, 3)= 0.000000D+00
       PDM( 8, 4)= 0.000000D+00
       PDM( 8, 5)= 0.000000D+00
       PDM( 8, 6)= 0.200000D+01
       PDM( 8, 7)= 0.200000D+01
c
2395  DO 10 I=1,3
        ISDATE(I)=ISD(I)
10 	CONTINUE
      DO 20 I=1,2
        ISTIME(I)=IST(I)
20 	CONTINUE
      NAME(1)='CIRA'
      NAME(2)='-86 '
      END

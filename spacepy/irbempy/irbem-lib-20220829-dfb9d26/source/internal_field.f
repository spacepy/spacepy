!***************************************************************************************************
! Copyright 1987, E G STASSINOPOULOS AND G D MEAD, 2006 S. Bourdarie
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
      SUBROUTINE get_intfield(xGEO,yGEO,zGEO,BxGEO,ByGEO,BzGEO)
C
C!    Geomagnetic field evaluation
C                       
C     REFERENCES
C     Subroutine ALLMAG of bint.for in blxtra
C     E G STASSINOPOULOS AND G D MEAD (see below)
C
C
        REAL*8    rkm, st, ct, sph, cph, br, bt, bp
        REAL*8    alti,lati,longi
        REAL*8    xGEO,yGEO,zGEO,BxGEO,ByGEO,BzGEO
        REAL*8    BGEO(3),BSPH(3)
C
C
C     HISTORY
C *  GEOCENTRIC VERSION OF GEOMAGNETIC FIELD ROUTINE
C *  DOUBLE PRECISION DECK FOR CII-HB / DPS05
C *  LONG DECK, THROUGH NMAX=13, FIXED INDICES WITHOUT DO LOOPS
C *  EXECUTION TIME PER CALL FACTOR OF THREE LESS THAN SHORT DECK
C *  PROGRAM DESIGNED AND TESTED BY E G STASSINOPOULOS AND G D MEAD,
C *  CODE 641, NASA GODDARD SPACE FLT CTR, GREENBELT, MD 20771
C *
C **          RKM        DISTANCE IN EARTH RADII TO EARTH'S CENTRE.
C **          ST,CT      SIN & COS OF GEOCENTRIC COLATITUDE
C **          SPH,CPH    SIN & COS OF EAST LONGITUDE
C **  OUTPUT& BR,BT,BP   GEOCENTRIC FIELD COMPONENTS IN NT
C **    NOTE& FOR GREATEST EFFICIENCY, COMPLETE ALL CALCULATIONS WITH
C **            ONE MODEL AND ONE TIME BEFORE CHANGING MODELS OR TIME.
C
C
C     VARIABLES
C
        REAL*8 AR, AOR
        REAL*8 P21, P22, SP2, CP2, DP21, DP22, C2, SP3,
     &         CP3, P31, P32, P33
        REAL*8 DP31, DP32, DP33, C3, SP4, CP4, P41, P42,
     &         P43, P44, DP41
        REAL*8 DP42, DP43, DP44, C4, SP5, CP5, P51, P52,
     &         P53, P54, P55
        REAL*8 DP51, DP52, DP53, DP54, DP55, C5, SP6, CP6,
     &         P61, P62, P63
        REAL*8 P64, P65, P66, DP61, DP62, DP63, DP64, DP65,
     &         DP66, C6
        REAL*8 SP7, CP7, P71, P72, P73, P74, P75, P76, P77, 
     &         DP71, DP72
        REAL*8 DP73, DP74, DP75, DP76, DP77, C7, SP8, CP8,
     &         P81, P82, P83
        REAL*8 P84, P85, P86, P87, P88, DP81, DP82, DP83,
     &         DP84, DP85, DP86
        REAL*8 DP87, DP88, C8, SP9, CP9, P91, P92, P93, P94,
     &         P95, P96, P97
        REAL*8 P98, P99, DP91, DP92, DP93, DP94, DP95, DP96,
     &         DP97, DP98
        REAL*8 DP99, C9, SP10, CP10, P101, P102, P103, P104,
     &         P105, P106
        REAL*8 P107, P108, P109, P1010, DP101, DP102, DP103,
     &         DP104, DP105
        REAL*8 DP106, DP107, DP108, DP109, DP1010, C10, SP11,
     &         CP11, P111
        REAL*8 P112, P113, P114, P115, P116, P117, P118, P119,
     &         P1110
        REAL*8 P1111, DP111, DP112, DP113, DP114, DP115, DP116,
     &         DP117
        REAL*8 DP118, DP119, DP1110, DP1111, C11, SP12, CP12, 
     &         P121, P122
        REAL*8 P123, P124, P125, P126, P127, P128, P129, P1210, 
     &         P1211
        REAL*8 P1212, DP121, DP122, DP123, DP124, DP125, DP126,
     &         DP127
        REAL*8 DP128, DP129, DP1210, DP1211, DP1212, C12, SP13,
     &         CP13, P131
        REAL*8 P132, P133, P134, P135, P136, P137, P138, P139,
     &         P1310
        REAL*8 P1311, P1312, P1313, DP131, DP132, DP133, DP134,
     &         DP135
        REAL*8 DP136, DP137, DP138, DP139, DP1310, DP1311, DP1312,
     &         DP1313
        REAL*8 C13
        REAL*8 G(16,16)
        INTEGER*4 i, j,norder
C
      COMMON /intfield/  G,norder
        REAL*8     pi,rad
        common /rconst/rad,pi
C      
C     
      CALL geo_gdz(xGEO,yGEO,zGEO,lati,longi,alti)
      rkm=SQRT(xGEO*xGEO+yGEO*yGEO+zGEO*zGEO)
      st=SIN((90.D0-lati)*rad)
      ct=COS((90.D0-lati)*rad)
      sph=SIN(longi*rad)
      cph=COS(longi*rad)

C
      P21=CT
      P22=ST
      AR=1.0D0/RKM
      SP2=SPH
      CP2=CPH
      DP21=-P22
      DP22=P21
      AOR=AR*AR*AR
      C2=G(2,2)*CP2+G(1,2)*SP2
      BR=-(AOR+AOR)*(G(2,1)*P21+C2*P22)
      BT=AOR*(G(2,1)*DP21+C2*DP22)
      BP=AOR*(G(1,2)*CP2-G(2,2)*SP2)*P22
      IF (norder .LE. 2) GO TO 1
C                                                           N= 3
      SP3=(SP2+SP2)*CP2
      CP3=(CP2+SP2)*(CP2-SP2)
      P31=P21*P21-0.333333333D0
      P32=P21*P22
      P33=P22*P22
      DP31=-P32-P32
      DP32=P21*P21-P33
      DP33=-DP31
      AOR=AOR*AR
      C2=G(3,2)*CP2+G(1,3)*SP2
      C3=G(3,3)*CP3+G(2,3)*SP3
      BR=BR-3.0D0*AOR*(G(3,1)*P31+C2*P32+C3*P33)
      BT=BT+AOR*(G(3,1)*DP31+C2*DP32+C3*DP33)
      BP=BP-AOR*((G(3,2)*SP2-G(1,3)*CP2)*P32+2.0D0*(G(3,3)*SP3-G(2,3)*
     &   CP3)*P33)
      IF (norder .LE. 3) GO TO 1
C                                                           N= 4
      SP4=SP2*CP3+CP2*SP3
      CP4=CP2*CP3-SP2*SP3
      P41=P21*P31-0.26666666D0*P21
      DP41=P21*DP31+DP21*P31-0.26666666D0*DP21
      P42=P21*P32-0.2D0*P22
      DP42=P21*DP32+DP21*P32-0.2D0*DP22
      P43=P21*P33
      DP43=P21*DP33+DP21*P33
      P44=P22*P33
      DP44=3.0D0*P43
      AOR=AOR*AR
      C2=G(4,2)*CP2+G(1,4)*SP2
      C3=G(4,3)*CP3+G(2,4)*SP3
      C4=G(4,4)*CP4+G(3,4)*SP4
      BR=BR-4.0D0*AOR*(G(4,1)*P41+C2*P42+C3*P43+C4*P44)
      BT=BT+AOR*(G(4,1)*DP41+C2*DP42+C3*DP43+C4*DP44)
      BP=BP-AOR*((G(4,2)*SP2-G(1,4)*CP2)*P42+2.0D0*(G(4,3)*SP3-G(2,4)*
     &   CP3)*P43+3.0D0*(G(4,4)*SP4-G(3,4)*CP4)*P44)
      IF (norder .LE. 4) GO TO 1
C                                                           N= 5
      SP5=(SP3+SP3)*CP3
      CP5=(CP3+SP3)*(CP3-SP3)
      P51=P21*P41-0.25714285D0*P31
      DP51=P21*DP41+DP21*P41-0.25714285D0*DP31
      P52=P21*P42-0.22857142D0*P32
      DP52=P21*DP42+DP21*P42-0.22857142D0*DP32
      P53=P21*P43-0.14285714D0*P33
      DP53=P21*DP43+DP21*P43-0.14285714D0*DP33
      P54=P21*P44
      DP54=P21*DP44+DP21*P44
      P55=P22*P44
      DP55=4.0D0*P54
      AOR=AOR*AR
      C2=G(5,2)*CP2+G(1,5)*SP2
      C3=G(5,3)*CP3+G(2,5)*SP3
      C4=G(5,4)*CP4+G(3,5)*SP4
      C5=G(5,5)*CP5+G(4,5)*SP5
      BR=BR-5.0D0*AOR*(G(5,1)*P51+C2*P52+C3*P53+C4*P54+C5*P55)
      BT=BT+AOR*(G(5,1)*DP51+C2*DP52+C3*DP53+C4*DP54+C5*DP55)
      BP=BP-AOR*((G(5,2)*SP2-G(1,5)*CP2)*P52+2.0D0*(G(5,3)*SP3-G(2,5)*
     &   CP3)*P53+3.0D0*(G(5,4)*SP4-G(3,5)*CP4)*P54+4.0D0*(G(5,5)*SP5-
     &   G(4,5)*CP5)*P55)
      IF (norder .LE. 5) GO TO 1
C                                                           N= 6
      SP6=SP2*CP5+CP2*SP5
      CP6=CP2*CP5-SP2*SP5
      P61=P21*P51-0.25396825D0*P41
      DP61=P21*DP51+DP21*P51-0.25396825D0*DP41
      P62=P21*P52-0.23809523D0*P42
      DP62=P21*DP52+DP21*P52-0.23809523D0*DP42
      P63=P21*P53-0.19047619D0*P43
      DP63=P21*DP53+DP21*P53-0.19047619D0*DP43
      P64=P21*P54-0.11111111D0*P44
      DP64=P21*DP54+DP21*P54-0.11111111D0*DP44
      P65=P21*P55
      DP65=P21*DP55+DP21*P55
      P66=P22*P55
      DP66=5.0D0*P65
      AOR=AOR*AR
      C2=G(6,2)*CP2+G(1,6)*SP2
      C3=G(6,3)*CP3+G(2,6)*SP3
      C4=G(6,4)*CP4+G(3,6)*SP4
      C5=G(6,5)*CP5+G(4,6)*SP5
      C6=G(6,6)*CP6+G(5,6)*SP6
      BR=BR-6.0D0*AOR*(G(6,1)*P61+C2*P62+C3*P63+C4*P64+C5*P65+C6*P66)
      BT=BT+AOR*(G(6,1)*DP61+C2*DP62+C3*DP63+C4*DP64+C5*DP65+C6*DP66)
      BP=BP-AOR*((G(6,2)*SP2-G(1,6)*CP2)*P62+2.0D0*(G(6,3)*SP3-G(2,6)*
     &   CP3)*P63+3.0D0*(G(6,4)*SP4-G(3,6)*CP4)*P64+4.0D0*(G(6,5)*SP5-
     &   G(4,6)*CP5)*P65+5.0D0*(G(6,6)*SP6-G(5,6)*CP6)*P66)
      IF (norder .LE. 6) GO TO 1
C                                                           N= 7
      SP7=(SP4+SP4)*CP4
      CP7=(CP4+SP4)*(CP4-SP4)
      P71=P21*P61-0.25252525D0*P51
      DP71=P21*DP61+DP21*P61-0.25252525D0*DP51
      P72=P21*P62-0.24242424D0*P52
      DP72=P21*DP62+DP21*P62-0.24242424D0*DP52
      P73=P21*P63-0.21212121D0*P53
      DP73=P21*DP63+DP21*P63-0.21212121D0*DP53
      P74=P21*P64-0.16161616D0*P54
      DP74=P21*DP64+DP21*P64-0.16161616D0*DP54
      P75=P21*P65-0.09090909D0*P55
      DP75=P21*DP65+DP21*P65-0.09090909D0*DP55
      P76=P21*P66
      DP76=P21*DP66+DP21*P66
      P77=P22*P66
      DP77=6.0D0*P76
      AOR=AOR*AR
      C2=G(7,2)*CP2+G(1,7)*SP2
      C3=G(7,3)*CP3+G(2,7)*SP3
      C4=G(7,4)*CP4+G(3,7)*SP4
      C5=G(7,5)*CP5+G(4,7)*SP5
      C6=G(7,6)*CP6+G(5,7)*SP6
      C7=G(7,7)*CP7+G(6,7)*SP7
      BR=BR-7.0D0*AOR*(G(7,1)*P71+C2*P72+C3*P73+C4*P74+C5*P75+C6*P76+
     &   C7*P77)
      BT=BT+AOR*(G(7,1)*DP71+C2*DP72+C3*DP73+C4*DP74+C5*DP75+C6*DP76+C7*
     &   DP77)
      BP=BP-AOR*((G(7,2)*SP2-G(1,7)*CP2)*P72+2.0D0*(G(7,3)*SP3-G(2,7)*
     &   CP3)*P73+3.0D0*(G(7,4)*SP4-G(3,7)*CP4)*P74+4.0D0*(G(7,5)*SP5-
     &   G(4,7)*CP5)*P75+5.0D0*(G(7,6)*SP6-G(5,7)*CP6)*P76+6.0D0*
     &   (G(7,7)*SP7-G(6,7)*CP7)*P77)
      IF (norder .LE. 7) GO TO 1
C                                                           N= 8
      SP8=SP2*CP7+CP2*SP7
      CP8=CP2*CP7-SP2*SP7
      P81=P21*P71-0.25174825D0*P61
      DP81=P21*DP71+DP21*P71-0.25174825D0*DP61
      P82=P21*P72-0.24475524D0*P62
      DP82=P21*DP72+DP21*P72-0.24475524D0*DP62
      P83=P21*P73-0.22377622D0*P63
      DP83=P21*DP73+DP21*P73-0.22377622D0*DP63
      P84=P21*P74-0.18881118D0*P64
      DP84=P21*DP74+DP21*P74-0.18881118D0*DP64
      P85=P21*P75-0.13986013D0*P65
      DP85=P21*DP75+DP21*P75-0.13986013D0*DP65
      P86=P21*P76-0.07692307D0*P66
      DP86=P21*DP76+DP21*P76-0.07692307D0*DP66
      P87=P21*P77
      DP87=P21*DP77+DP21*P77
      P88=P22*P77
      DP88=7.0D0*P87
      AOR=AOR*AR
      C2=G(8,2)*CP2+G(1,8)*SP2
      C3=G(8,3)*CP3+G(2,8)*SP3
      C4=G(8,4)*CP4+G(3,8)*SP4
      C5=G(8,5)*CP5+G(4,8)*SP5
      C6=G(8,6)*CP6+G(5,8)*SP6
      C7=G(8,7)*CP7+G(6,8)*SP7
      C8=G(8,8)*CP8+G(7,8)*SP8
      BR=BR-8.0D0*AOR*(G(8,1)*P81+C2*P82+C3*P83+C4*P84+C5*P85+C6*P86+
     &   C7*P87+C8*P88)
      BT=BT+AOR*(G(8,1)*DP81+C2*DP82+C3*DP83+C4*DP84+C5*DP85+C6*DP86+C7*
     &   DP87+C8*DP88)
      BP=BP-AOR*((G(8,2)*SP2-G(1,8)*CP2)*P82+2.0D0*(G(8,3)*SP3-G(2,8)*
     &   CP3)*P83+3.0D0*(G(8,4)*SP4-G(3,8)*CP4)*P84+4.0D0*(G(8,5)*SP5-
     &   G(4,8)*CP5)*P85+5.0D0*(G(8,6)*SP6-G(5,8)*CP6)*P86+6.0D0*
     &   (G(8,7)*SP7-G(6,8)*CP7)*P87+7.0D0*(G(8,8)*SP8-G(7,8)*CP8)*P88)
      IF (norder .LE. 8) GO TO 1
C                                                           N= 9
      SP9=(SP5+SP5)*CP5
      CP9=(CP5+SP5)*(CP5-SP5)
      P91=P21*P81-0.25128205D0*P71
      DP91=P21*DP81+DP21*P81-0.25128205D0*DP71
      P92=P21*P82-0.24615384D0*P72
      DP92=P21*DP82+DP21*P82-0.24615384D0*DP72
      P93=P21*P83-0.23076923D0*P73
      DP93=P21*DP83+DP21*P83-0.23076923D0*DP73
      P94=P21*P84-0.20512820D0*P74
      DP94=P21*DP84+DP21*P84-0.20512820D0*DP74
      P95=P21*P85-0.16923076D0*P75
      DP95=P21*DP85+DP21*P85-0.16923076D0*DP75
      P96=P21*P86-0.12307692D0*P76
      DP96=P21*DP86+DP21*P86-0.12307692D0*DP76
      P97=P21*P87-0.06666666D0*P77
      DP97=P21*DP87+DP21*P87-0.06666666D0*DP77
      P98=P21*P88
      DP98=P21*DP88+DP21*P88
      P99=P22*P88
      DP99=8.0D0*P98
      AOR=AOR*AR
      C2=G(9,2)*CP2+G(1,9)*SP2
      C3=G(9,3)*CP3+G(2,9)*SP3
      C4=G(9,4)*CP4+G(3,9)*SP4
      C5=G(9,5)*CP5+G(4,9)*SP5
      C6=G(9,6)*CP6+G(5,9)*SP6
      C7=G(9,7)*CP7+G(6,9)*SP7
      C8=G(9,8)*CP8+G(7,9)*SP8
      C9=G(9,9)*CP9+G(8,9)*SP9
      BR=BR-9.0D0*AOR*(G(9,1)*P91+C2*P92+C3*P93+C4*P94+C5*P95+C6*P96+
     &   C7*P97+C8*P98+C9*P99)
      BT=BT+AOR*(G(9,1)*DP91+C2*DP92+C3*DP93+C4*DP94+C5*DP95+C6*DP96+C7*
     &   DP97+C8*DP98+C9*DP99)
      BP=BP-AOR*((G(9,2)*SP2-G(1,9)*CP2)*P92+2.0D0*(G(9,3)*SP3-G(2,9)*
     &   CP3)*P93+3.0D0*(G(9,4)*SP4-G(3,9)*CP4)*P94+4.0d0*(G(9,5)*SP5-
     &   G(4,9)*CP5)*P95+5.0D0*(G(9,6)*SP6-G(5,9)*CP6)*P96+6.0D0*
     &   (G(9,7)*SP7-G(6,9)*CP7)*P97+7.0D0*(G(9,8)*SP8-G(7,9)*CP8)*P98+
     &   8.0D0*(G(9,9)*SP9-G(8,9)*CP9)*P99)
      IF (norder .LE. 9) GO TO 1
C                                                           N=10
      SP10=SP2*CP9+CP2*SP9
      CP10=CP2*CP9-SP2*SP9
      P101=P21*P91-0.25098039D0*P81
      DP101=P21*DP91+DP21*P91-0.25098039D0*DP81
      P102=P21*P92-0.24705882D0*P82
      DP102=P21*DP92+DP21*P92-0.24705882D0*DP82
      P103=P21*P93-0.23529411D0*P83
      DP103=P21*DP93+DP21*P93-0.23529411D0*DP83
      P104=P21*P94-0.21568627D0*P84
      DP104=P21*DP94+DP21*P94-0.21568627D0*DP84
      P105=P21*P95-0.18823529D0*P85
      DP105=P21*DP95+DP21*P95-0.18823529D0*DP85
      P106=P21*P96-0.15294117D0*P86
      DP106=P21*DP96+DP21*P96-0.15294117D0*DP86
      P107=P21*P97-0.10980392D0*P87
      DP107=P21*DP97+DP21*P97-0.10980392D0*DP87
      P108=P21*P98-0.05882352D0*P88
      DP108=P21*DP98+DP21*P98-0.05882352D0*DP88
      P109=P21*P99
      DP109=P21*DP99+DP21*P99
      P1010=P22*P99
      DP1010=9.0D0*P109
      AOR=AOR*AR
      C2=G(10,2)*CP2+G(1,10)*SP2
      C3=G(10,3)*CP3+G(2,10)*SP3
      C4=G(10,4)*CP4+G(3,10)*SP4
      C5=G(10,5)*CP5+G(4,10)*SP5
      C6=G(10,6)*CP6+G(5,10)*SP6
      C7=G(10,7)*CP7+G(6,10)*SP7
      C8=G(10,8)*CP8+G(7,10)*SP8
      C9=G(10,9)*CP9+G(8,10)*SP9
      C10=G(10,10)*CP10+G(9,10)*SP10
      BR=BR-10.0D0*AOR*(G(10,1)*P101+C2*P102+C3*P103+C4*P104+C5*P105+
     &   C6*P106+C7*P107+C8*P108+C9*P109+C10*P1010)
      BT=BT+AOR*(G(10,1)*DP101+C2*DP102+C3*DP103+C4*DP104+C5*DP105+C6*
     &   DP106+C7*DP107+C8*DP108+C9*DP109+C10*DP1010)
      BP=BP-AOR*((G(10,2)*SP2-G(1,10)*CP2)*P102+2.0D0*(G(10,3)*SP3-
     &   G(2,10)*CP3)*P103+3.0D0*(G(10,4)*SP4-G(3,10)*CP4)*P104+4.0D0*
     &   (G(10,5)*SP5-G(4,10)*CP5)*P105+5.0D0*(G(10,6)*SP6-G(5,10)*CP6)*
     &   P106+6.0D0*(G(10,7)*SP7-G(6,10)*CP7)*P107+7.0D0*(G(10,8)*SP8-
     &   G(7,10)*CP8)*P108+8.0D0*(G(10,9)*SP9-G(8,10)*CP9)*P109+9.0D0*
     &   (G(10,10)*SP10-G(9,10)*CP10)*P1010)
      IF (norder .LE. 10) GO TO 1
C                                                           N=11
      SP11=(SP6+SP6)*CP6
      CP11=(CP6+SP6)*(CP6-SP6)
      P111=P21*P101-0.25077399D0*P91
      DP111=P21*DP101+DP21*P101-0.25077399D0*DP91
      P112=P21*P102-0.24767801D0*P92
      DP112=P21*DP102+DP21*P102-0.24767801D0*DP92
      P113=P21*P103-0.23839009D0*P93
      DP113=P21*DP103+DP21*P103-0.23839009D0*DP93
      P114=P21*P104-0.22291021D0*P94
      DP114=P21*DP104+DP21*P104-0.22291021D0*DP94
      P115=P21*P105-0.20123839D0*P95
      DP115=P21*DP105+DP21*P105-0.20123839D0*DP95
      P116=P21*P106-0.17337461D0*P96
      DP116=P21*DP106+DP21*P106-0.17337461D0*DP96
      P117=P21*P107-0.13931888D0*P97
      DP117=P21*DP107+DP21*P107-0.13931888D0*DP97
      P118=P21*P108-0.09907120D0*P98
      DP118=P21*DP108+DP21*P108-0.09907120D0*DP98
      P119=P21*P109-0.05263157D0*P99
      DP119=P21*DP109+DP21*P109-0.05263157D0*DP99
      P1110=P21*P1010
      DP1110=P21*DP1010+DP21*P1010
      P1111=P22*P1010
      DP1111=10.0D0*P1110
      AOR=AOR*AR
      C2=G(11,2)*CP2+G(1,11)*SP2
      C3=G(11,3)*CP3+G(2,11)*SP3
      C4=G(11,4)*CP4+G(3,11)*SP4
      C5=G(11,5)*CP5+G(4,11)*SP5
      C6=G(11,6)*CP6+G(5,11)*SP6
      C7=G(11,7)*CP7+G(6,11)*SP7
      C8=G(11,8)*CP8+G(7,11)*SP8
      C9=G(11,9)*CP9+G(8,11)*SP9
      C10=G(11,10)*CP10+G(9,11)*SP10
      C11=G(11,11)*CP11+G(10,11)*SP11
      BR=BR-11.0D0*AOR*(G(11,1)*P111+C2*P112+C3*P113+C4*P114+C5*P115+
     &   C6*P116+C7*P117+C8*P118+C9*P119+C10*P1110+C11*P1111)
      BT=BT+AOR*(G(11,1)*DP111+C2*DP112+C3*DP113+C4*DP114+C5*DP115+C6*
     &   DP116+C7*DP117+C8*DP118+C9*DP119+C10*DP1110+C11*DP1111)
      BP=BP-AOR*((G(11,2)*SP2-G(1,11)*CP2)*P112+2.0D0*(G(11,3)*SP3-
     &   G(2,11)*CP3)*P113+3.0D0*(G(11,4)*SP4-G(3,11)*CP4)*P114+4.0D0*
     &   (G(11,5)*SP5-G(4,11)*CP5)*P115+5.0D0*(G(11,6)*SP6-G(5,11)*CP6)*
     &   P116+6.0D0*(G(11,7)*SP7-G(6,11)*CP7)*P117+7.0D0*(G(11,8)*SP8-
     &   G(7,11)*CP8)*P118+8.0D0*(G(11,9)*SP9-G(8,11)*CP9)*P119+9.0D0*
     &   (G(11,10)*SP10-G(9,11)*CP10)*P1110+10.0D0*(G(11,11)*SP11-
     &   G(10,11)*CP11)*P1111)
      IF (norder .LE. 11) GO TO 1
C                                                           N=12
      SP12=SP2*CP11+CP2*SP11
      CP12=CP2*CP11-SP2*SP11
      P121=P21*P111-0.25062656D0*P101
      DP121=P21*DP111+DP21*P111-0.25062656D0*DP101
      P122=P21*P112-0.24812030D0*P102
      DP122=P21*DP112+DP21*P112-0.24812030D0*DP102
      P123=P21*P113-0.24060150D0*P103
      DP123=P21*DP113+DP21*P113-0.24060150D0*DP103
      P124=P21*P114-0.22807017D0*P104
      DP124=P21*DP114+DP21*P114-0.22807017D0*DP104
      P125=P21*P115-0.21052631D0*P105
      DP125=P21*DP115+DP21*P115-0.21052631D0*DP105
      P126=P21*P116-0.18796992D0*P106
      DP126=P21*DP116+DP21*P116-0.18796992D0*DP106
      P127=P21*P117-0.16040100D0*P107
      DP127=P21*DP117+DP21*P117-0.16040100D0*DP107
      P128=P21*P118-0.12781954D0*P108
      DP128=P21*DP118+DP21*P118-0.12781954D0*DP108
      P129=P21*P119-0.09022556D0*P109
      DP129=P21*DP119+DP21*P119-0.09022556D0*DP109
      P1210=P21*P1110-0.04761904D0*P1010
      DP1210=P21*DP1110+DP21*P1110-0.04761904D0*DP1010
      P1211=P21*P1111
      DP1211=P21*DP1111+DP21*P1111
      P1212=P22*P1111
      DP1212=11.0D0*P1211
      AOR=AOR*AR
      C2=G(12,2)*CP2+G(1,12)*SP2
      C3=G(12,3)*CP3+G(2,12)*SP3
      C4=G(12,4)*CP4+G(3,12)*SP4
      C5=G(12,5)*CP5+G(4,12)*SP5
      C6=G(12,6)*CP6+G(5,12)*SP6
      C7=G(12,7)*CP7+G(6,12)*SP7
      C8=G(12,8)*CP8+G(7,12)*SP8
      C9=G(12,9)*CP9+G(8,12)*SP9
      C10=G(12,10)*CP10+G(9,12)*SP10
      C11=G(12,11)*CP11+G(10,12)*SP11
      C12=G(12,12)*CP12+G(11,12)*SP12
      BR=BR-12.0D0*AOR*(G(12,1)*P121+C2*P122+C3*P123+C4*P124+C5*P125+
     &   C6*P126+C7*P127+C8*P128+C9*P129+C10*P1210+C11*P1211+C12*P1212)
      BT=BT+AOR*(G(12,1)*DP121+C2*DP122+C3*DP123+C4*DP124+C5*DP125+C6*
     &   DP126+C7*DP127+C8*DP128+C9*DP129+C10*DP1210+C11*DP1211+C12*
     &   DP1212)
      BP=BP-AOR*((G(12,2)*SP2-G(1,12)*CP2)*P122+2.0D0*(G(12,3)*SP3-
     &   G(2,12)*CP3)*P123+3.0D0*(G(12,4)*SP4-G(3,12)*CP4)*P124+4.0D0*
     &   (G(12,5)*SP5-G(4,12)*CP5)*P125+5.0D0*(G(12,6)*SP6-G(5,12)*CP6)*
     &   P126+6.0D0*(G(12,7)*SP7-G(6,12)*CP7)*P127+7.0D0*(G(12,8)*SP8-
     &   G(7,12)*CP8)*P128+8.0D0*(G(12,9)*SP9-G(8,12)*CP9)*P129+9.0D0*
     &   (G(12,10)*SP10-G(9,12)*CP10)*P1210+10.0D0*(G(12,11)*SP11-
     &   G(10,12)*CP11)*P1211+11.0D0*(G(12,12)*SP12-G(11,12)*CP12)*
     &   P1212)
      IF (norder .LE. 12) GO TO 1
C                                                           N=13
      SP13=(SP7+SP7)*CP7
      CP13=(CP7+SP7)*(CP7-SP7)
      P131=P21*P121-0.25051759D0*P111
      DP131=P21*DP121+DP21*P121-0.25051759D0*DP111
      P132=P21*P122-0.24844720D0*P112
      DP132=P21*DP122+DP21*P122-0.24844720D0*DP112
      P133=P21*P123-0.24223602D0*P113
      DP133=P21*DP123+DP21*P123-0.24223602D0*DP113
      P134=P21*P124-0.23188405D0*P114
      DP134=P21*DP124+DP21*P124-0.23188405D0*DP114
      P135=P21*P125-0.21739130D0*P115
      DP135=P21*DP125+DP21*P125-0.21739130D0*DP115
      P136=P21*P126-0.19875776D0*P116
      DP136=P21*DP126+DP21*P126-0.19875776D0*DP116
      P137=P21*P127-0.17598343D0*P117
      DP137=P21*DP127+DP21*P127-0.17598343D0*DP117
      P138=P21*P128-0.14906832D0*P118
      DP138=P21*DP128+DP21*P128-0.14906832D0*DP118
      P139=P21*P129-0.11801242D0*P119
      DP139=P21*DP129+DP21*P129-0.11801242D0*DP119
      P1310=P21*P1210-0.08281573D0*P1110
      DP1310=P21*DP1210+DP21*P1210-0.08281573D0*DP1110
      P1311=P21*P1211-0.04347826D0*P1111
      DP1311=P21*DP1211+DP21*P1211-0.04347826D0*DP1111
      P1312=P21*P1212
      DP1312=P21*DP1212+DP21*P1212
      P1313=P22*P1212
      DP1313=12.0D0*P1312
      AOR=AOR*AR
      C2=G(13,2)*CP2+G(1,13)*SP2
      C3=G(13,3)*CP3+G(2,13)*SP3
      C4=G(13,4)*CP4+G(3,13)*SP4
      C5=G(13,5)*CP5+G(4,13)*SP5
      C6=G(13,6)*CP6+G(5,13)*SP6
      C7=G(13,7)*CP7+G(6,13)*SP7
      C8=G(13,8)*CP8+G(7,13)*SP8
      C9=G(13,9)*CP9+G(8,13)*SP9
      C10=G(13,10)*CP10+G(9,13)*SP10
      C11=G(13,11)*CP11+G(10,13)*SP11
      C12=G(13,12)*CP12+G(11,13)*SP12
      C13=G(13,13)*CP13+G(12,13)*SP13
      BR=BR-13.0D0*AOR*(G(13,1)*P131+C2*P132+C3*P133+C4*P134+C5*P135+
     &   C6*P136+C7*P137+C8*P138+C9*P139+C10*P1310+C11*P1311+C12*P1312+
     &   C13*P1313)
      BT=BT+AOR*(G(13,1)*DP131+C2*DP132+C3*DP133+C4*DP134+C5*DP135+C6*
     &   DP136+C7*DP137+C8*DP138+C9*DP139+C10*DP1310+C11*DP1311+C12*
     &   DP1312+C13*DP1313)
      BP=BP-AOR*((G(13,2)*SP2-G(1,13)*CP2)*P132+2.0D0*(G(13,3)*SP3-
     &   G(2,13)*CP3)*P133+3.0D0*(G(13,4)*SP4-G(3,13)*CP4)*P134+4.0D0*
     &   (G(13,5)*SP5-G(4,13)*CP5)*P135+5.0D0*(G(13,6)*SP6-G(5,13)*CP6)*
     &   P136+6.0D0*(G(13,7)*SP7-G(6,13)*CP7)*P137+7.0D0*(G(13,8)*SP8-
     &   G(7,13)*CP8)*P138+8.0D0*(G(13,9)*SP9-G(8,13)*CP9)*P139+9.0D0*
     &   (G(13,10)*SP10-G(9,13)*CP10)*P1310+10.0D0*(G(13,11)*SP11-
     &   G(10,13)*CP11)*P1311+11.0D0*(G(13,12)*SP12-G(11,13)*CP12)*
     &   P1312+12.0D0*(G(13,13)*SP13-G(12,13)*CP13)*P1313)
      IF (norder .LE. 13) GO TO 1
C
    1 if( abs(st).lt.1.e-6 )then
        BP = -1.d31
      else       
        BP = BP / ST
      endif       
C
      BSPH(1)=BR
      BSPH(2)=BT
      BSPH(3)=BP
      call SPH_CAR_VECT(rkm,lati,longi,BSPH,BGEO)
      BxGEO=BGEO(1)
      ByGEO=BGEO(2)
      BzGEO=BGEO(3)
C
      END
c
c
c
C
C ======================================================================
C The Jensen&Cain1960 mag model has been used to generate AE8 min and max and AP8 min model.
C
C ======================================================================
c
      SUBROUTINE JensenANDCain1960
C
C!    Jensen & Cain 1960 model coefficients
C
C     VARIABLES
C
        REAL*8        gjc(27,2)
      REAL*8        G(16,16)
      REAL*8        gg(66),hh(66),thet,phit
      REAL*8        Bo,xc,yc,zc,ct,st,cp,sp
        INTEGER*4     i, j, k, n, norder
C
      COMMON /intfield/  G,norder
      COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C      
      DATA ((gjc(i,j),j=1,2),i=1,27) 
     & / 30411.2d0,    0.0d0,  2147.4d0, -5798.9d0,  2403.5d0,    0.0d0,
     &   -5125.3d0, 3312.4d0, -1338.1d0,  -157.9d0, -3151.8d0,    0.0d0,
     &    6213.0d0, 1487.0d0, -2489.8d0,  -407.5d0,  -649.6d0,   21.0d0,
     &   -4179.4d0,    0.0d0, -4529.8d0, -1182.5d0, -2179.5d0, 1000.6d0,
     &     700.8d0,   43.0d0,  -204.4d0,   138.5d0,  1625.6d0,    0.0d0,
     &   -3440.7d0,  -79.6d0, -1944.7d0,  -200.0d0,   -60.8d0,  459.7d0,
     &     277.5d0,  242.1d0,    69.7d0,  -121.8d0, -1952.3d0,    0.0d0,
     &    -485.3d0, -575.8d0,   321.2d0,  -873.5d0,  2141.3d0, -340.6d0,
     &     105.1d0,  -11.8d0,    22.7d0,  -111.6d0,   111.5d0,  -32.5d0/
C
C     CODE
C
      gg(1)=0.D0
      hh(1)=0.D0
      DO i=1,27
         gg(i+1)=-gjc(i,1)
       hh(i+1)=-gjc(i,2)
      enddo
      norder    = 7
      G(1,1) = 0.0D0
      n = 1
      DO k=2,norder
        G(k,1)     = gjc(n,1)
        DO i=2,k
          G(k,i)   = gjc(n+i-1,1)
          G(i-1,k) = gjc(n+i-1,2)
        END DO
        n = n+k
      END DO
C
      CALL get_terms(gg,hh,thet,phit,xc,yc,zc,Bo)
C
      ct = COS(thet)
      st = SIN(thet)
      cp = COS(phit)
      sp = SIN(phit)
C
      END
C
C ======================================================================
C The GSFC12/66 updated to 1970 mag model has been used to generate AP8 max model.
C
C ======================================================================
C
      SUBROUTINE GSFC1266
C
C    GSFC 12/66 model coefficients
C 
        IMPLICIT NONE
C            
        REAL*8        ggsfc(65,2), ggsfc1(65,2), ggsfc2(65,2)
        REAL*8        gg(16,16), ggt(16,16),ggtt(16,16),g(16,16)
      REAL*8        shmit(16,16),t,dum
      REAL*8        ggg(66),hhh(66),thet,phit
      REAL*8        Bo,xc,yc,zc,ct,st,cp,sp
        INTEGER*4     i, j, k, n,jj,m,norder
c
      COMMON /intfield/  G,norder
      COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C      
      DATA ((ggsfc(i,j),j=1,2),i=1,65)
     & / -30401.2d0,0.0d0,-2163.8d0,5778.2d0,-1540.1d0,0.0d0,
     & 2997.9d0,-1932.0d0,1590.3d0,202.9d0,1307.1d0,0.0d0,-1988.9d0,
     & -425.4d0,1276.8d0,227.8d0,881.2d0,-133.8d0,949.3d0,0.0d0,
     & 803.5d0,160.3d0,502.9d0,-274.3d0,-397.7d0,2.3d0,266.5d0,
     & -246.6d0,-233.5d0,0.0d0,355.7d0,5.1d0,228.4d0,117.8d0,-28.8d0,
     & -114.8d0,-157.9d0,-108.9d0,-62.2d0,82.4d0,49.2d0,0.0d0,57.5d0,
     & -12.1d0,-0.8d0,104.4d0,-238.3d0,56.6d0,-1.5d0,-23.4d0,-2.0d0,
     & -14.8d0,-108.9d0,-13.3d0,72.2d0,0.0d0,-53.7d0,-53.7d0,7.9d0,
     & -27.4d0,15.6d0,-8.1d0,-24.3d0,7.0d0,-3.6d0,24.3d0,15.5d0,
     & -22.5d0,3.6d0,-21.4d0,8.5d0,0.0d0,6.5d0,5.4d0,-9.3d0,-11.7d0,
     & -9.6d0,4.2d0,-6.1d0,-15.3d0,5.5d0,4.6d0,-8.1d0,21.9d0,13.0d0,
     & -0.7d0,7.4d0,-17.1d0,10.4d0,0.0d0,5.8d0,-22.4d0,7.5d0,13.8d0,
     & -15.1d0,6.3d0,12.1d0,-3.0d0,4.7d0,-1.9d0,0.2d0,9.0d0,1.6d0,
     & 11.5d0,0.9d0,0.1d0,0.2d0,-1.5d0,-2.9d0,0.0d0,-0.9d0,-0.1d0,
     & -2.2d0,4.5d0,0.8d0,-1.0d0,-2.8d0,2.6d0,6.4d0,-4.4d0,4.7d0,
     & -1.3d0,-0.2d0,-3.6d0,1.8d0,4.0d0,2.0d0,1.0d0,1.1d0,-2.0d0 /     
      DATA ((ggsfc1(i,j),j=1,2),i=1,65)
     & / 14.03d0,0.00d0,8.76d0,-3.71d0,-23.29d0,0.00d0,-0.09d0,
     & -14.31d0,-4.56d0,-16.62d0,-0.93d0,0.00d0,-10.62d0,5.20d0,2.31d0,
     & 2.53d0,-5.89d0,-6.98d0,1.45d0,0.00d0,0.90d0,-2.19d0,-1.75d0,
     & -0.14d0,0.66d0,1.88d0,-3.01d0,-6.52d0,1.61d0,0.00d0,0.60d0,
     & 2.24d0,3.34d0,1.59d0,-0.04d0,-2.61d0,-0.60d0,0.50d0,1.76d0,
     & -0.12d0,-0.42d0,0.00d0,0.82d0,0.05d0,0.82d0,0.09d0,2.35d0,
     & 2.55d0,0.83d0,-1.19d0,0.01d0,0.33d0,0.23d0,0.84d0,-0.57d0,
     & 0.00d0,-0.34d0,-0.96d0,-1.44d0,0.01d0,-0.90d0,0.43d0,0.03d0,
     & 0.75d0,-0.60d0,-0.33d0,-0.17d0,0.49d0,-0.64d0,0.90d0,0.35d0,
     & 0.00d0,0.50d0,-0.50d0,1.70d0,-0.21d0,-0.11d0,0.03d0,0.34d0,
     & -0.79d0,-0.07d0,0.05d0,0.43d0,0.10d0,-0.15d0,-0.36d0,-0.42d0,
     & -0.43d0,-0.10d0,0.00d0,-0.13d0,0.66d0,-1.20d0,0.54d0,0.08d0,
     & 0.03d0,-0.08d0,0.35d0,-0.39d0,-0.03d0,-0.36d0,-0.01d0,0.47d0,
     & 0.45d0,0.37d0,-0.05d0,-0.46d0,0.75d0,-0.01d0,0.00d0,-0.13d0,
     & -0.61d0,0.88d0,-0.64d0,-0.18d0,0.02d0,0.17d0,0.05d0,-0.02d0,
     & -0.63d0,0.05d0,-0.07d0,0.17d0,0.07d0,0.16d0,-0.03d0,0.31d0,
     & -0.02d0,-0.23d0,-0.45d0 /
      DATA ((ggsfc2(i,j),j=1,2),i=1,65)
     & / -0.062d0,0.000d0,0.114d0,-0.043d0,-0.154d0,0.000d0,-0.018d0,
     & 0.054d0,-0.253d0,-0.016d0,-0.123d0,0.000d0,-0.027d0,0.095d0,
     & 0.028d0,-0.007d0,-0.183d0,0.079d0,0.001d0,0.000d0,-0.044d0,
     & 0.004d0,0.017d0,0.056d0,0.007d0,-0.035d0,-0.097d0,-0.047d0,
     & 0.045d0,0.000d0,0.001d0,-0.046d0,0.075d0,0.007d0,0.008d0,
     & -0.007d0,0.015d0,0.001d0,0.056d0,-0.024d0,-0.006d0,0.000d0,
     & 0.015d0,0.020d0,0.010d0,-0.011d0,0.050d0,0.015d0,-0.011d0,
     & -0.029d0,0.026d0,0.029d0,0.023d0,-0.010d0,-0.014d0,0.000d0,
     & -0.006d0,-0.014d0,-0.034d0,0.016d0,-0.004d0,0.014d0,-0.006d0,
     & 0.005d0,-0.027d0,-0.008d0,-0.001d0,0.016d0,-0.004d0,0.011d0,
     & 0.006d0,0.000d0,0.008d0,-0.015d0,0.039d0,-0.012d0,-0.008d0,
     & 0.005d0,0.015d0,-0.011d0,-0.002d0,0.000d0,0.005d0,-0.003d0,
     & -0.008d0,-0.009d0,-0.007d0,-0.003d0,-0.005d0,0.000d0,-0.001d0,
     & 0.022d0,-0.027d0,0.007d0,0.005d0,-0.002d0,-0.007d0,0.009d0,
     & -0.006d0,0.006d0,-0.009d0,-0.001d0,0.006d0,0.009d0,0.005d0,
     & -0.004d0,-0.009d0,0.019d0,-0.003d0,0.000d0,-0.003d0,-0.012d0,
     & 0.020d0,-0.014d0,-0.008d0,0.001d0,0.007d0,0.001d0,0.001d0,
     & -0.011d0,0.001d0,-0.001d0,0.001d0,0.001d0,0.005d0,-0.001d0,
     & 0.004d0,0.001d0,-0.002d0,-0.006d0 /
C
      norder     = 11
      gg(1,1)         = 0.0D0
      ggt(1,1)        = 0.0D0
      ggtt(1,1)       = 0.0D0        
      n = 1
      DO k=2,norder
        gg(k,1)       = ggsfc(n,1)
        ggt(k,1)      = ggsfc1(n,1)
        ggtt(k,1)     = ggsfc2(n,1)
        DO j=2,k
          gg(k,j)     = ggsfc(n+j-1,1)
          gg(j-1,k)   = ggsfc(n+j-1,2)
          ggt(k,j)    = ggsfc1(n+j-1,1)
          ggt(j-1,k)  = ggsfc1(n+j-1,2)
          ggtt(k,j)   = ggsfc2(n+j-1,1)
          ggtt(j-1,k) = ggsfc2(n+j-1,2)
        ENDDO
        n = n + k
      ENDDO
C
C  The GSFC12/66 is extended from year 1960 to 1970.
C
      t = 1970.d0 - 1960.d0  
C
C Perform a Schmidt normalisation
C      
C
      shmit(1,1)       = -1.0D0
      DO n=2,16
        shmit(n,1)     = (2*n-3) * shmit(n-1,1) / (n-1)
        jj = 2
        DO m=2,n
          dum          = (n-m+1)*jj
          shmit(n,m)   = shmit(n,m-1)*SQRT(Dum/(n+m-2))
          shmit(m-1,n) = shmit(n,m)
          jj           = 1
        ENDDO
      ENDDO
C
C
      DO n=1,norder
        DO m=1,norder
          g(n,m) = (gg(n,m)+t*ggt(n,m)+t*t*ggtt(n,m)) * shmit(n,m)
        ENDDO
      ENDDO
C
      n = 1
      ggg(1)=0.D0
      hhh(1)=0.D0
      DO k=2,norder
        ggg(n+1)=-G(k,1) 
        DO i=2,k
          ggg(n+i)=-G(k,i) 
          hhh(n+i)=-G(i-1,k)
        END DO
        n = n+k
      END DO
C
      CALL get_terms(ggg,hhh,thet,phit,xc,yc,zc,Bo)
C
      ct = COS(thet)
      st = SIN(thet)
      cp = COS(phit)
      sp = SIN(phit)
C
      END

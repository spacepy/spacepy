!***************************************************************************************************
! Copyright 1986, N. Tsyganenko
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
      SUBROUTINE TSY87L(IOPT, X, Y, Z, BX, BY, BZ)
C-----------------------------------------------------------------------
C  'LONG' version of the magnetospheric magnetic field model.
C  Computes GSM components of the magnetic field of extraterrestrial
C  current systems up to geocentric radial distances about 70 Re.
C  Corresponds to magnetospheric field model (N.A.Tsyganenko, 1986)
C  based on IMP-A,C,D,E,F,G,H,I,J (years 1966-1974) and HEOS-1,2 
C  (1969-1974) satellite merged data set.
C
C  INPUT:   IOPT  specifies model version, according to ground 
C           disturbance level, as specified below
C        PS geodipole tilt angle in radians
C        X,Y,Z    GSM coordinates of the point in Earth radii
C
C  OUTPUT:  BX,BY,BZ    GSM components of the external
C              magnetic field in nanotesla
C
C  IOPT     1      2         3        4         5        6
C   KP    0,0+  1-,1,1+   2-,2,2+  3-,3,3+   4-,4,4+   >=5-  
C-----------------------------------------------------------------------

      IMPLICIT NONE

C Input
      REAL*8 ps, x, y, z,tilt
      INTEGER*4 IOPT
C Output
      REAL*8 bx, by, bz
C Local
      REAL*8 GA(32,6), PA(32), HPI, FC, RT, X1, X2
      REAL*8 C1, RRC2, DSTR, XN, RH, DY, B0, B1, B2, XN21, XN2, XNR
      REAL*8 XN22, ADLN, SPS, CPS, RPS, ZS, ZP, ZM, FY, XNX, XNX2
      REAL*8 XC1, XC2, XC22, XR2, XC12, B20, B2P, B2M, B, BP, BM, XA1
      REAL*8 XAP1, XAM1, XA2, XAP2, XAM2, XNA, XNAP, XNAM, F, FP, FM
      REAL*8 XLN1, XLNP1, XLNM1, XLN2, XLNP2, XLNM2, ALN, S0, S0P, S0M
      REAL*8 S1, S1P, S1M, S2, S2P, S2M, G1, G1P, G1M, G2, G2P, G2M
      REAL*8 EX1, EX2, Y2, Z2, YZ, XSM, ZSM, RR, RR2, ZN, BRSM, BZSM
      REAL*8 BY1, BXSM
      INTEGER*4 I, IP
      SAVE IP,PA,C1,RRC2,DSTR,XN,RH,X1,DY,B0,B1,XN21,XN2,XNR,XN22,ADLN
C
      COMMON /dip_ang/tilt
C    &-.1058,-3.221,-.00114,-.02166,-30.43,.04049,.05464,.008884,42.,
C
C errata 1.29.90 D.P. Stern
C
      DATA GA /-.09673d0,-10.63d0,1.21d0,34.57d0,-.04502d0,-.06553d0,
     :-.02952d0,.3852d0,-.03665d0,-2.084d0,.001795d0,.00638d0,-23.49d0,
     :.06082d0,.01642d0,-.02137d0,32.21d0,-.04373d0,-.02311d0,-.2832d0,
     :-.002303d0,-.000631d0,-6.397d0,-967.d0,-8650.d0,-20.55d0,5.18d0,
     :-2.796d0,2.715d0,13.58d0,8.038d0,29.21d0,-.485d0,-12.84d0,
     :1.856d0,40.06d0,-.0294d0,-.09071d0,-.02993d0,.5465d0,-.04928d0,
     :-2.453d0,.001587d0,.007402d0,-29.41d0,.08101d0,.02322d0,-.1091d0,
     :40.75d0,-.07995d0,-.03859d0,-.2755d0,-.002759d0,-.000408d0,
     :-6.189d0,-957.8d0,-7246.d0,-25.51d0,5.207d0,-4.184d0,2.641d0,
     :16.56d0,7.795d0,29.36d0,-1.132d0,-18.05d0,2.625d0,48.55d0,
     :-.004868d0,-.1087d0,-.03824d0,.8514d0,-.0522d0,-2.881d0,
     :-.000295d0,.009055d0,-29.48d0,.06394d0,.03864d0,-.2288d0,41.77d0,
     :-.05849d0,-.06443d0,-.4683d0,.001222d0,-.000519d0,-3.696d0,
     :-991.1d0,-6955.d0,-31.43d0,4.878d0,-3.151d0,3.277d0,19.19d0,
     :7.248d0,28.99d0,-1.003d0,-16.98d0,3.14d0,52.81d0,-.08625d0,
     :-.1478d0,-.03501d0,.55d0,-.07778d0,-2.97d0,.002086d0,.01275d0,
     :-26.79d0,.06328d0,.03622d0,.08345d0,39.72d0,-.06009d0,-.07825d0,
     :-.9698d0,.000178d0,-.000573d0,-.9328d0,-872.5d0,-5851.d0,
     :-39.68d0,4.902d0,-3.848d0,2.79d0,20.91d0,6.193d0,26.81d0,
     :-1.539d0,-14.29d0,3.479d0,53.36d0,-.004201d0,-.2043d0,-.03932d0,
     :.6409d0,-.1058d0,-3.221d0,-.00114d0,.02166d0,-30.43d0,.04049d0,
     :.05464d0,.008884d0,42.d0,-.01035d0,-.1053d0,-1.63d0,.003802d0,
     :-.001029d0,4.204d0,-665.6d0,-1011.d0,-43.49d0,4.514d0,-2.948d0,
     :2.99d0,21.59d0,6.005d0,22.d0,-2.581d0,-7.726d0,5.045d0,53.31d0,
     :.02262d0,-.1972d0,-.01981d0,.428d0,-.1055d0,-5.075d0,.002762d0,
     :.03277d0,-27.35d0,.04986d0,.06119d0,-.1211d0,47.48d0,-.0502d0,
     :-.1477d0,.838d0,-.01008d0,-.0057d0,9.231d0,-674.3d0,-900.d0,
     :-74.43d0,4.658d0,-3.245d0,3.39d0,21.8d0,5.62d0,25.17d0/

      DATA IP, FC, RT, X1, X2 /100, 0.3183099031D0, 30.0D0, 4.0D0,5.0D0/
      ps=tilt*4.D0*ATAN(1.D0)/180.d0

      HPI = 2.0D0 * DATAN(1.0D0)

      IF (IOPT .NE. IP) THEN
        IP = IOPT
        DO I=1,32
          PA(I) = GA(I,IP)
        END DO
        C1 = PA(29)**2
        RRC2 = PA(27)**2
        DSTR = PA(26) / RRC2 * 4.0D0
        XN = PA(28)
        RH = PA(31)
        DY = PA(30)
        B0 = PA(23)
        B1 = PA(24)
        B2 = PA(25)
        XN21 = (XN-X1)**2
        XN2 = XN - X2
        XNR = 1.0 / XN2
        XN22 = XN2**2
        ADLN = DLOG(XN22/XN21)
      END IF
      SPS = DSIN(PS)
      CPS = DCOS(PS)
      RPS = RH * SPS
C
C   COMPUTATION BEGINS HERE IF PARAMETER IOPT REMAINED UNCHANGED AFTER THE 
C   PRECEDING CALL OF THIS SUBROUTINE.
C
      ZS = Z - RPS
      ZP = Z - RT
      ZM = Z + RT
      FY = FC / (1.0D0+(Y/DY)**2.0D0)
      XNX = XN - X
      XNX2 = XNX**2.0D0
      XC1 = X - X1
      XC2 = X - X2
      XC22 = XC2**2.0D0
      XR2 = XC2 * XNR
      XC12 = XC1**2.0D0
      B20 = ZS**2 + C1
      B2P = ZP**2 + C1
      B2M = ZM**2 + C1
      B = DSQRT(B20)
      BP = DSQRT(B2P)
      BM = DSQRT(B2M)
      XA1 = XC12 + B20
      XAP1 = XC12 + B2P
      XAM1 = XC12 + B2M
      XA2 = 1.0D0 / (XC22+B20)
      XAP2 = 1.0D0 / (XC22+B2P)
      XAM2 = 1.0D0 / (XC22+B2M)
      XNA = XNX2 + B20
      XNAP = XNX2 + B2P
      XNAM = XNX2 + B2M
      F = B20 - XC22
      FP = B2P - XC22
      FM = B2M - XC22
      XLN1 = DLOG(XN21/XNA)
      XLNP1 = DLOG(XN21/XNAP)
      XLNM1 = DLOG(XN21/XNAM)
      XLN2 = XLN1 + ADLN
      XLNP2 = XLNP1 + ADLN
      XLNM2 = XLNM1 + ADLN
      ALN = 0.25D0 * (XLNP1+XLNM1-2.0D0*XLN1)
      S0 = (DATAN(XNX/B)+HPI) / B
      S0P = (DATAN(XNX/BP)+HPI) / BP
      S0M = (DATAN(XNX/BM)+HPI) / BM
      S1 = (XLN1*0.5D0+XC1*S0) / XA1
      S1P = (XLNP1*0.5D0+XC1*S0P) / XAP1
      S1M = (XLNM1*0.5D0+XC1*S0M) / XAM1
      S2 = (XC2*XA2*XLN2-XNR-F*XA2*S0) * XA2
      S2P = (XC2*XAP2*XLNP2-XNR-FP*XAP2*S0P) * XAP2
      S2M = (XC2*XAM2*XLNM2-XNR-FM*XAM2*S0M) * XAM2
      G1 = (B20*S0-0.5D0*XC1*XLN1) / XA1
      G1P = (B2P*S0P-0.5D0*XC1*XLNP1) / XAP1
      G1M = (B2M*S0M-0.5D0*XC1*XLNM1) / XAM1
      G2 = ((0.5D0*F*XLN2+2.0D0*S0*B20*XC2)*XA2+XR2) * XA2
      G2P = ((0.5D0*FP*XLNP2+2.0D0*S0P*B2P*XC2)*XAP2+XR2) * XAP2
      G2M = ((0.5D0*FM*XLNM2+2.0D0*S0M*B2M*XC2)*XAM2+XR2) * XAM2
      BX = FY * (B0*(ZS*S0-0.5D0*(ZP*S0P+ZM*S0M))+B1*
     &     (ZS*S1-0.5D0*(ZP*S1P+ZM*S1M))+B2*
     &     (ZS*S2-0.5D0*(ZP*S2P+ZM*S2M)))
      BY = 0.0D0
      BZ = FY * (B0*ALN+B1*(G1-0.5D0*(G1P+G1M))+B2*(G2-0.5D0*(G2P+G2M)))
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED.
C
      EX1 = DEXP(X/PA(32))
      EX2 = EX1**2
      Y2 = Y**2
      Z2 = Z**2
      YZ = Y * Z
      BX = BX + (EX1*PA(1)+EX2*PA(3)) * Z * CPS +
     &     (EX1*PA(2)+EX2*(PA(4)+PA(5)*Y2+PA(6)*Z2)) * SPS
      BY = (EX1*PA(7)+EX2*PA(9)) *YZ * CPS +
     &     (EX1*PA(8)+EX2*(PA(10)+PA(11)*Y2+PA(12)*Z2)) * Y * SPS
      BZ = BZ + (EX1*(PA(13)+PA(14)*Y2+PA(15)*Z2)+EX2*(PA(17)+PA(18)*Y2+
     &     PA(19)*Z2))*CPS+(EX1*PA(16)+EX2*(PA(20)+PA(21)*Y2+PA(22)*Z2))
     &     * Z * SPS
C
C   DCF AND FAC CONTRIBUTION HAS BEEN ADDED TO BX, BY, AND BZ.
C
      XSM = X * CPS - Z * SPS
      ZSM = X * SPS + Z * CPS
      Z2 = ZSM**2
      RR2 = XSM**2 + Y**2
      RR = DSQRT(RR2)
      ZN = ((RR2+Z2)/RRC2+4.0D0)**2.5D0
      BRSM = DSTR * 3.0D0 * ZSM / ZN
      BZSM = DSTR * (2.0D0*Z2-RR2+8.0D0*RRC2) / ZN
      BY1 = BRSM * Y
      BXSM = BRSM * XSM
      BX = BX + BXSM * CPS + BZSM * SPS
      BZ = BZ - BXSM * SPS + BZSM * CPS
      BY = BY + BY1
C
C   RING CURRENT FIELD HAS BEEN TAKEN INTO ACCOUNT.
C
      RETURN
      END

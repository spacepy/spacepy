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
C-----------------------------------------------------------------------
      SUBROUTINE TSY87S(IOPT,  X, Y, Z, BX, BY, BZ)
C-----------------------------------------------------------------------
C  Computes GSM components of the magnetic field of extraterrestrial 
C  current systems up to geocentric radial distances about 30 Re. 
C  Corresponds to magnetospheric field model (N.A.Tsyganenko, 1986). 
C  Based on IMP-A,C,D,E,F,G,H,I,J (years 1966-1980) and HEOS-1,2 
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
C  IOPT    1     2      3       4      5        6         7      8
C   KP    0,0+  1-,1   1+,2-   2,2+  3-,3,3+   4-,4,4+   >=5-   >=5+
C-----------------------------------------------------------------------
      
      IMPLICIT NONE

C Input
      REAL*8 x, y, z, ps,tilt
      INTEGER*4 IOPT
C Output
      REAL*8 bx, by, bz
C Local
      REAL*8 GA(24,8), PA(24), HPI, FC, RT
      REAL*8 C1, RRC2, DSTR, XN, RH, X1, DY, B0, B1, XN21, SPS, CPS, RPS
      REAL*8 ZS, ZP, ZM, FY, XNX, XNX2, XC1, XC12, B20, B2P, B2M, B, BP
      REAL*8 BM, XA1, XAP1, XAM1, XNA, XNAP, XNAM, XLN1, XLNP1, XLNM1
      REAL*8 ALN, S0, S0P, S0M, S1, S1P, S1M, G1, G1P, G1M, EX, Y2, Z2
      REAL*8 YZ, XSM, ZSM, RR, RR2, ZN, BRSM, BZSM, BY1, BXSM
      INTEGER*4 I, IP
      SAVE IP,PA,C1,RRC2,DSTR,XN,RH,X1,DY,B0,B1,XN21
C
      COMMON /dip_ang/tilt
C    &12.72,-.00867,-.001953,-.3437,-.002903,-.000999,18.41,-270.3,
C
C errata 1.29.90 D.P. Stern
C
      DATA GA /1.126d0,26.66d0,-.077d0,-.06102d0,-.06197d0,-2.048d0,
     :.00327d0,.008473d0,12.72d0,-.00867d0,-.01953d0,-.3437d0,
     :-.002903d0,-.000999d0,18.41d0,-270.3d0,-25.94d0,5.21d0,-6.2d0,
     :2.29d0,11.96d0,8.315d0,44.22d0,11.15d0,1.403d0,29.24d0,-.0693d0,
     :-.0864d0,-.07202d0,-2.068d0,.00286d0,.007438d0,16.37d0,-.02705d0,
     :-.0281d0,-.604d0,-.002256d0,.000152d0,20.2d0,-140.1d0,-29.65d0,
     :5.62d0,-5.52d0,2.02d0,14.66d0,8.06d0,27.76d0,10.94d0,1.589d0,
     :31.07d0,-.06527d0,-.07447d0,-.07632d0,-2.413d0,.002719d0,
     :.01098d0,16.2d0,-.02355d0,-.03475d0,-.4377d0,-.002169d0,
     :-.001383d0,18.7d0,-292.6d0,-35.25d0,5.29d0,-5.18d0,2.21d0,
     :14.03d0,7.66d0,17.56d0,10.9d0,1.699d0,36.28d0,-.07514d0,-.1448d0,
     :-.08049d0,-2.209d0,.000919d0,.01084d0,17.38d0,-.03516d0,
     :-.03886d0,-1.169d0,.004239d0,.000881d0,21.79d0,-162.d0,-41.87d0,
     :5.15d0,-3.62d0,2.35d0,17.26d0,7.61d0,17.99d0,10.74d0,2.141d0,
     :41.51d0,-.1518d0,-.1857d0,-.1015d0,-2.929d0,.004584d0,.01589d0,
     :18.29d0,-.02514d0,-.05927d0,-1.336d0,.00185d0,.001066d0,21.31d0,
     :-358.8d0,-47.91d0,5.13d0,-3.74d0,2.07d0,17.23d0,6.33d0,32.51d0,
     :9.73d0,2.252d0,39.35d0,-.04525d0,-.2062d0,-.1491d0,-3.059d0,
     :-.000183d0,.02614d0,15.48d0,-.02144d0,-.06608d0,-1.855d0,
     :.006199d0,-.00013d0,23.91d0,-161.d0,-51.48d0,4.61d0,-3.32d0,
     :1.68d0,15.22d0,6.68d0,.6765d0,8.007d0,2.773d0,40.95d0,.00667d0,
     :-.133d0,-.1304d0,-5.187d0,.004623d0,.03651d0,20.d0,-.03765d0,
     :-.09066d0,.5838d0,-.01462d0,-.007189d0,24.87d0,-186.1d0,-74.81d0,
     :4.57d0,-4.03d0,1.7d0,12.15d0,6.87d0,-1.746d0,8.9d0,2.919d0,
     :34.96d0,2*0.d0,-.1609d0,-5.077d0,2*0.d0,22.1d0,-.05915d0,
     :-.1051d0,.6321d0,2*0.d0,28.11d0,-330.1d0,-86.82d0,4.d0,-3.d0,
     :1.73d0,12.56d0,5.11d0,4.d0,7.866d0/

      DATA IP, FC, RT /100, 0.3183099031D0, 30.0D0/
      ps=tilt*4.D0*ATAN(1.D0)/180.d0

      HPI = 2.0D0 * DATAN(1.0D0)

      IF (IOPT .NE. IP) THEN
        IP = IOPT
        DO I=1,24
          PA(I) = GA(I,IP)         
        END DO
        C1 = PA(20)**2
        RRC2 = PA(18)**2
        DSTR = PA(17) / RRC2 * 4.0D0
        XN = PA(19)
        RH = PA(22)
        X1 = PA(23)
        DY = PA(21)
        B0 = PA(15)
        B1 = PA(16)
        XN21 = (XN-X1)**2
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
      FY = FC / (1.0D0+(Y/DY)**2)
      XNX = XN - X
      XNX2 = XNX**2
      XC1 = X - X1
      XC12 = XC1**2
      B20 = ZS**2 + C1
      B2P = ZP**2 + C1
      B2M = ZM**2 + C1
      B = DSQRT(B20)
      BP = DSQRT(B2P)
      BM = DSQRT(B2M)
      XA1 = XC12 + B20
      XAP1 = XC12 + B2P
      XAM1 = XC12 + B2M
      XNA = XNX2 + B20
      XNAP = XNX2 + B2P
      XNAM = XNX2 + B2M
      XLN1 = DLOG(XN21/XNA)
      XLNP1 = DLOG(XN21/XNAP)
      XLNM1 = DLOG(XN21/XNAM)
      ALN = 0.25D0 * (XLNP1+XLNM1-2.0D0*XLN1)
      S0 = (DATAN(XNX/B)+HPI) / B
      S0P = (DATAN(XNX/BP)+HPI) / BP
      S0M = (DATAN(XNX/BM)+HPI) / BM
      S1 = (XLN1*0.5D0+XC1*S0) / XA1
      S1P = (XLNP1*0.5D0+XC1*S0P) / XAP1
      S1M = (XLNM1*0.5D0+XC1*S0M) / XAM1
      G1 = (B20*S0-0.5D0*XC1*XLN1) / XA1
      G1P = (B2P*S0P-0.5D0*XC1*XLNP1) / XAP1
      G1M = (B2M*S0M-0.5D0*XC1*XLNM1) / XAM1
      BX = FY * (B0*(ZS*S0-0.5D0*(ZP*S0P+ZM*S0M))
     &     +B1*(ZS*S1-0.5D0*(ZP*S1P+ZM*S1M)))
      BY = 0.0D0
      BZ = FY * (B0*ALN+B1*(G1-0.5D0*(G1P+G1M)))
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      EX = DEXP(X/PA(24))
      Y2 = Y**2
      Z2 = Z**2
      YZ = Y * Z
      BX = BX + EX * (CPS*PA(1)*Z+SPS*(PA(2)+PA(3)*YZ+PA(4)*Z2))
      BY = EX * (CPS*PA(5)*YZ+SPS*Y*(PA(6)+PA(7)*Y2+PA(8)*Z2))
      BZ = BZ + EX * (CPS*(PA(9)+PA(10)*Y2+PA(11)*Z2)+SPS*Z*(PA(12)
     &     +PA(13)*Y2+PA(14)*Z2))
C
C   DCF AND FAC CONTRIBUTION HAS BEEN ADDED TO BX, BY, AND BZ
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
C   RING CURRENT FIELD HAS BEEN TAKEN INTO ACCOUNT
C
      RETURN
      END
C

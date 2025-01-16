!***************************************************************************************************
! Copyright 1989, 1992, N. Tsyganenko
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
      SUBROUTINE T89C(IOPT,X,Y,Z,BX,BY,BZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C
C   COMPUTES GSM COMPONENTS OF THE MAGNETIC FIELD PRODUCED BY EXTRA-
C   TERRESTRIAL CURRENT SYSTEMS IN THE GEOMAGNETOSPHERE. THE MODEL IS
C   VALID UP TO GEOCENTRIC DISTANCES OF 70 RE AND IS BASED ON THE MER-
C   GED IMP-A,C,D,E,F,G,H,I,J (1966-1974), HEOS-1 AND -2 (1969-1974),
C   AND ISEE-1 AND -2  SPACECRAFT DATA SET.
C
C   THIS IS A MODIFIED VERSION (T89c), WHICH REPLACED THE ORIGINAL ONE
C     IN 1992 AND DIFFERS FROM IT IN THE FOLLOWING:
C
C   (1)  ISEE-1,2 DATA WERE ADDED TO THE ORIGINAL IMP-HEOS DATASET
C   (2)  TWO TERMS WERE ADDED TO THE ORIGINAL TAIL FIELD MODES, ALLOWING
C          A MODULATION OF THE CURRENT BY THE GEODIPOLE TILT ANGLE
C
C
C  REFERENCE FOR THE ORIGINAL MODEL: N.A. TSYGANENKO, A MAGNETOSPHERIC MAGNETIC
C       FIELD MODEL WITH A WARPED TAIL CURRENT SHEET: PLANET.SPACE SCI., V.37,
C         PP.5-20, 1989.
C
C----INPUT PARAMETERS: IOPT - SPECIFIES THE GROUND DISTURBANCE LEVEL:
C
C   IOPT= 1       2        3        4        5        6      7
C                  CORRESPOND TO:
C    KP= 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  > =6-
C
C    PS - GEODIPOLE TILT ANGLE IN RADIANS
C    X, Y, Z  - GSM COORDINATES OF THE POINT IN EARTH RADII
C
C----OUTPUT PARAMETERS: BX,BY,BZ - GSM COMPONENTS OF THE MODEL MAGNETIC
C                        FIELD IN NANOTESLAS
c
c   THE PARAMETER PARMOD(10) IS A DUMMY ARRAY.  IT IS NOT USED IN THIS
C        SUBROUTINE AND IS PROVIDED JUST FOR MAKING IT COMPATIBLE WITH THE
C           NEW VERSION (4/16/96) OF THE GEOPACK SOFTWARE.
C
C   THIS RELEASE OF T89C IS DATED  FEB 12, 1996;
C--------------------------------------------------------------------------
C
C
C              AUTHOR:     NIKOLAI A. TSYGANENKO
C                          HSTX CORP./NASA GSFC
C
       DIMENSION XI(4),F(3),DER(3,30),PARAM(30,7),A(30),PARMOD(10)
c-mk-c       DOUBLE PRECISION F,DER
       real*8 F,DER,tilt
      COMMON /dip_ang/tilt
        DATA PARAM/-116.53,-10719.,42.375,59.753,-11363.,1.7844,30.268,
     * -0.35372E-01,-0.66832E-01,0.16456E-01,-1.3024,0.16529E-02,
     * 0.20293E-02,20.289,-0.25203E-01,224.91,-9234.8,22.788,7.8813,
     * 1.8362,-0.27228,8.8184,2.8714,14.468,32.177,0.01,0.0,
     * 7.0459,4.0,20.0,-55.553,-13198.,60.647,61.072,-16064.,
     * 2.2534,34.407,-0.38887E-01,-0.94571E-01,0.27154E-01,-1.3901,
     * 0.13460E-02,0.13238E-02,23.005,-0.30565E-01,55.047,-3875.7,
     * 20.178,7.9693,1.4575,0.89471,9.4039,3.5215,14.474,36.555,
     * 0.01,0.0,7.0787,4.0,20.0,-101.34,-13480.,111.35,12.386,-24699.,
     * 2.6459,38.948,-0.34080E-01,-0.12404,0.29702E-01,-1.4052,
     * 0.12103E-02,0.16381E-02,24.49,-0.37705E-01,-298.32,4400.9,18.692,
     * 7.9064,1.3047,2.4541,9.7012,7.1624,14.288,33.822,0.01,0.0,6.7442,
     * 4.0,20.0,-181.69,-12320.,173.79,-96.664,-39051.,3.2633,44.968,
     * -0.46377E-01,-0.16686,0.048298,-1.5473,0.10277E-02,0.31632E-02,
     * 27.341,-0.50655E-01,-514.10,12482.,16.257,8.5834,1.0194,3.6148,
     * 8.6042,5.5057,13.778,32.373,0.01,0.0,7.3195,4.0,20.0,-436.54,
     * -9001.0,323.66,-410.08,-50340.,3.9932,58.524,-0.38519E-01,
     * -0.26822,0.74528E-01,-1.4268,-0.10985E-02,0.96613E-02,27.557,
     * -0.56522E-01,-867.03,20652.,14.101,8.3501,0.72996,3.8149,9.2908,
     *  6.4674,13.729,28.353,0.01,0.0,7.4237,4.0,20.0,-707.77,-4471.9,
     * 432.81,-435.51,-60400.,4.6229,68.178,-0.88245E-01,-0.21002,
     * 0.11846,-2.6711,0.22305E-02,0.10910E-01,27.547,-0.54080E-01,
     * -424.23,1100.2,13.954,7.5337,0.89714,3.7813,8.2945,5.174,14.213,
     * 25.237,0.01,0.0,7.0037,4.0,20.0,-1190.4,2749.9,742.56,-1110.3,
     * -77193.,7.6727,102.05,-0.96015E-01,-0.74507,0.11214,-1.3614,
     * 0.15157E-02,0.22283E-01,23.164,-0.74146E-01,-2219.1,48253.,
     * 12.714,7.6777,0.57138,2.9633,9.3909,9.7263,11.123,21.558,0.01,
     * 0.0,4.4518,4.0,20.0/

        DATA IOP/10/
C
      ps=tilt*4.D0*ATAN(1.D0)/180.d0
c
         id=0
         IF (IOP.NE.IOPT) THEN
C
            ID=1
            IOP=IOPT
            DO 1 I=1,30
   1        A(I)=PARAM(I,IOPT)
C
         ENDIF
C
        XI(1)=X
        XI(2)=Y
        XI(3)=Z
        XI(4)=PS
         CALL T89(ID,A,XI,F,DER)
          IF (ID.EQ.1) ID=2
        BX=F(1)
        BY=F(2)
        BZ=F(3)
        RETURN
       END
C-------------------------------------------------------------------
C
          SUBROUTINE  T89 (ID, A, XI, F, DER)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C        ***  N.A. Tsyganenko ***  8-10.12.1991  ***
C
C      Calculates dependent model variables and their deriva-
C  tives for given independent variables and model parame-
C  ters.  Specifies model functions with free parameters which
C  must be determined by means of least squares fits (RMS
C  minimization procedure).
C
C      Description of parameters:
C
C  ID  - number of the data point in a set (initial assignments are performed
c        only for ID=1, saving thus CPU time)
C  A   - input vector containing model parameters;
C  XI  - input vector containing independent variables;
C  F   - output double precision vector containing
C        calculated values of dependent variables;
C  DER   - output double precision vector containing
C        calculated values for derivatives of dependent
C        variables with respect to model parameters;
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C      T89 represents external magnetospheric magnetic field
C  in Cartesian SOLAR MAGNETOSPHERIC coordinates (Tsyganenko N.A.,
C  Planet. Space Sci., 1989, v.37, p.5-20; the "T89 model" with the warped
c  tail current sheet) + A MODIFICATION ADDED IN APRIL 1992 (SEE BELOW)
C
C      Model formulas for the magnetic field components contain in total
c  30 free parameters (17 linear and 13 nonlinear parameters).
C      First 2 independent linear parameters A(1)-A(2) correspond to contribu-
c  tion from the tail current system, then follow A(3) and A(4) which are the
c  amplitudes of symmetric and antisymmetric terms in the contribution from
c  the closure currents; A(5) is the ring current amplitude. Then follow the
c coefficients A(6)-A(15) which define Chapman-Ferraro+Birkeland current field.
c    The coefficients c16-c19  (see Formula 20 in the original paper),
c   due to DivB=0 condition, are expressed through A(6)-A(15) and hence are not
c    independent ones.
c  A(16) AND A(17) CORRESPOND TO THE TERMS WHICH YIELD THE TILT ANGLE DEPEN-
C    DENCE OF THE TAIL CURRENT INTENSITY (ADDED ON APRIL 9, 1992)
C
C      Nonlinear parameters:
C
C    A(18) : DX - Characteristic scale of the Chapman-Ferraro field along the
c        X-axis
C    A(19) : ADR (aRC) - Characteristic radius of the ring current
c    A(20) : D0 - Basic half-thickness of the tail current sheet
C    A(21) : DD (GamRC)- defines rate of thickening of the ring current, as
c             we go from night- to dayside
C    A(22) : Rc - an analog of "hinging distance" entering formula (11)
C    A(23) : G - amplitude of tail current warping in the Y-direction
C    A(24) : aT - Characteristic radius of the tail current
c    A(25) : Dy - characteristic scale distance in the Y direction entering
c                 in W(x,y) in (13)
c    A(26) : Delta - defines the rate of thickening of the tail current sheet
c                 in the Y-direction (in T89 it was fixed at 0.01)
c    A(27) : Q - this parameter was fixed at 0 in the final version of T89;
c              initially it was introduced for making Dy to depend on X
c    A(28) : Sx (Xo) - enters in W(x,y) ; see (13)
c    A(29) : Gam (GamT) - enters in DT in (13) and defines rate of tail sheet
c              thickening on going from night to dayside; in T89 fixed at 4.0
c    A(30) : Dyc - the Dy parameter for closure current system; in T89 fixed
c               at 20.0
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C
         REAL*8  A(*), XI(*)
C
         DIMENSION  F(3), DER(3,30)
C
         INTEGER*4  ID, I, L
         save
         DATA A02,XLW2,YN,RPI,RT/25.D0,170.D0,30.D0,0.31830989D0,30.D0/
         DATA XD,XLD2/0.D0,40.D0/
C
C   The last four quantities define variation of tail sheet thickness along X
C
      DATA SXC,XLWC2/4.D0,50.D0/
C
C   The two quantities belong to the function WC which confines tail closure
c    current in X- and Y- direction
C
      DATA DXL/20.D0/
C
C
         IF (ID.NE.1)  GOTO  3
           DO  2  I = 1, 30
             DO  1  L = 1, 3
  1            DER(L,I) = 0.0D0
  2        CONTINUE
C
       DYC=A(30)
       DYC2=DYC**2
       DX=A(18)
       HA02=0.5D0*A02
       RDX2M=-1.D0/DX**2
       RDX2=-RDX2M
       RDYC2=1.D0/DYC2
       HLWC2M=-0.5D0*XLWC2
       DRDYC2=-2.D0*RDYC2
       DRDYC3=2.D0*RDYC2*DSQRT(RDYC2)
       HXLW2M=-0.5D0*XLW2
       ADR=A(19)
       D0=A(20)
       DD=A(21)
       RC=A(22)
       G=A(23)
       AT=A(24)
       DT=D0
       DEL=A(26)
       P=A(25)
       Q=A(27)
       SX=A(28)
       GAM=A(29)
       HXLD2M=-0.5D0*XLD2
       ADSL=0.D0
       XGHS=0.D0
       H=0.D0
       HS=0.D0
       GAMH=0.D0
       W1=-0.5D0/DX
       DBLDEL=2.D0*DEL
       W2=W1*2.D0
       W4=-1.D0/3.D0
       W3=W4/DX
       W5=-0.5D0
       W6=-3.D0
       AK1=A(1)
       AK2=A(2)
       AK3=A(3)
       AK4=A(4)
       AK5=A(5)
       AK6=A(6)
       AK7=A(7)
       AK8=A(8)
       AK9=A(9)
       AK10=A(10)
       AK11=A(11)
       AK12=A(12)
       AK13=A(13)
       AK14=A(14)
       AK15=A(15)
        AK16=A(16)
        AK17=A(17)
       SXA=0.D0
       SYA=0.D0
       SZA=0.D0
       AK610=AK6*W1+AK10*W5
       AK711=AK7*W2-AK11
       AK812=AK8*W2+AK12*W6
       AK913=AK9*W3+AK13*W4
       RDXL=1.D0/DXL
       HRDXL=0.5D0*RDXL
       A6H=AK6*0.5D0
       A9T=AK9/3.D0
       YNP=RPI/YN*0.5D0
       YND=2.D0*YN
C
  3      CONTINUE
C
           X  = XI(1)
           Y  = XI(2)
           Z  = XI(3)
           TILT=XI(4)
               TLT2=TILT**2
           SPS = DSIN(TILT)
           CPS = DSQRT (1.0D0 - SPS ** 2)
C
       X2=X*X
       Y2=Y*Y
       Z2=Z*Z
       TPS=SPS/CPS
       HTP=TPS*0.5D0
       GSP=G*SPS
       XSM=X*CPS-Z*SPS
       ZSM=X*SPS+Z*CPS
C
C   CALCULATE THE FUNCTION ZS DEFINING THE SHAPE OF THE TAIL CURRENT SHEET
C    AND ITS SPATIAL DERIVATIVES:
C
       XRC=XSM+RC
       XRC16=XRC**2+16.D0
       SXRC=DSQRT(XRC16)
       Y4=Y2*Y2
       Y410=Y4+1.D4
       SY4=SPS/Y410
       GSY4=G*SY4
       ZS1=HTP*(XRC-SXRC)
       DZSX=-ZS1/SXRC
       ZS=ZS1-GSY4*Y4
       D2ZSGY=-SY4/Y410*4.D4*Y2*Y
       DZSY=G*D2ZSGY
C
C   CALCULATE THE COMPONENTS OF THE RING CURRENT CONTRIBUTION:
C
       XSM2=XSM**2
       DSQT=DSQRT(XSM2+A02)
       FA0=0.5D0*(1.D0+XSM/DSQT)
       DDR=D0+DD*FA0
       DFA0=HA02/DSQT**3
       ZR=ZSM-ZS
       TR=DSQRT(ZR**2+DDR**2)
       RTR=1.D0/TR
       RO2=XSM2+Y2
       ADRT=ADR+TR
       ADRT2=ADRT**2
       FK=1.D0/(ADRT2+RO2)
       DSFC=DSQRT(FK)
       FC=FK**2*DSFC
       FACXY=3.0D0*ADRT*FC*RTR
       XZR=XSM*ZR
       YZR=Y*ZR
       DBXDP=FACXY*XZR
       DER(2,5)=FACXY*YZR
       XZYZ=XSM*DZSX+Y*DZSY
       FAQ=ZR*XZYZ-DDR*DD*DFA0*XSM
       DBZDP=FC*(2.D0*ADRT2-RO2)+FACXY*FAQ
       DER(1,5)=DBXDP*CPS+DBZDP*SPS
       DER(3,5)=DBZDP*CPS-DBXDP*SPS
C
C  CALCULATE THE TAIL CURRENT SHEET CONTRIBUTION:
C
       DELY2=DEL*Y2
       D=DT+DELY2
       IF (DABS(GAM).LT.1.D-6) GOTO 8
       XXD=XSM-XD
       RQD=1.D0/(XXD**2+XLD2)
       RQDS=DSQRT(RQD)
       H=0.5D0*(1.D0+XXD*RQDS)
       HS=-HXLD2M*RQD*RQDS
       GAMH=GAM*H
       D=D+GAMH
       XGHS=XSM*GAM*HS
       ADSL=-D*XGHS
   8   D2=D**2
       T=DSQRT(ZR**2+D2)
       XSMX=XSM-SX
       RDSQ2=1.D0/(XSMX**2+XLW2)
       RDSQ=DSQRT(RDSQ2)
       V=0.5D0*(1.D0-XSMX*RDSQ)
       DVX=HXLW2M*RDSQ*RDSQ2
       OM=DSQRT(DSQRT(XSM2+16.D0)-XSM)
       OMS=-OM/(OM*OM+XSM)*0.5D0
       RDY=1.D0/(P+Q*OM)
       OMSV=OMS*V
       RDY2=RDY**2
       FY=1.D0/(1.D0+Y2*RDY2)
       W=V*FY
       YFY1=2.D0*FY*Y2*RDY2
       FYPR=YFY1*RDY
       FYDY=FYPR*FY
       DWX=DVX*FY+FYDY*Q*OMSV
       YDWY=-V*YFY1*FY
       DDY=DBLDEL*Y
       ATT=AT+T
       S1=DSQRT(ATT**2+RO2)
       F5=1.D0/S1
       F7=1.D0/(S1+ATT)
       F1=F5*F7
       F3=F5**3
       F9=ATT*F3
       FS=ZR*XZYZ-D*Y*DDY+ADSL
       XDWX=XSM*DWX+YDWY
       RTT=1.D0/T
       WT=W*RTT
       BRRZ1=WT*F1
       BRRZ2=WT*F3
       DBXC1=BRRZ1*XZR
       DBXC2=BRRZ2*XZR
       DER(2,1)=BRRZ1*YZR
       DER(2,2)=BRRZ2*YZR
          DER(2,16)=DER(2,1)*TLT2
          DER(2,17)=DER(2,2)*TLT2
       WTFS=WT*FS
       DBZC1=W*F5+XDWX*F7+WTFS*F1
       DBZC2=W*F9+XDWX*F1+WTFS*F3
       DER(1,1)=DBXC1*CPS+DBZC1*SPS
       DER(1,2)=DBXC2*CPS+DBZC2*SPS
       DER(3,1)=DBZC1*CPS-DBXC1*SPS
       DER(3,2)=DBZC2*CPS-DBXC2*SPS
          DER(1,16)=DER(1,1)*TLT2
          DER(1,17)=DER(1,2)*TLT2
          DER(3,16)=DER(3,1)*TLT2
          DER(3,17)=DER(3,2)*TLT2
C
C  CALCULATE CONTRIBUTION FROM THE CLOSURE CURRENTS
C
       ZPL=Z+RT
       ZMN=Z-RT
       ROGSM2=X2+Y2
       SPL=DSQRT(ZPL**2+ROGSM2)
       SMN=DSQRT(ZMN**2+ROGSM2)
       XSXC=X-SXC
       RQC2=1.D0/(XSXC**2+XLWC2)
       RQC=DSQRT(RQC2)
       FYC=1.D0/(1.D0+Y2*RDYC2)
       WC=0.5D0*(1.D0-XSXC*RQC)*FYC
       DWCX=HLWC2M*RQC2*RQC*FYC
       DWCY=DRDYC2*WC*FYC*Y
       SZRP=1.D0/(SPL+ZPL)
       SZRM=1.D0/(SMN-ZMN)
       XYWC=X*DWCX+Y*DWCY
       WCSP=WC/SPL
       WCSM=WC/SMN
       FXYP=WCSP*SZRP
       FXYM=WCSM*SZRM
       FXPL=X*FXYP
       FXMN=-X*FXYM
       FYPL=Y*FXYP
       FYMN=-Y*FXYM
       FZPL=WCSP+XYWC*SZRP
       FZMN=WCSM+XYWC*SZRM
       DER(1,3)=FXPL+FXMN
       DER(1,4)=(FXPL-FXMN)*SPS
       DER(2,3)=FYPL+FYMN
       DER(2,4)=(FYPL-FYMN)*SPS
       DER(3,3)=FZPL+FZMN
       DER(3,4)=(FZPL-FZMN)*SPS
C
C   NOW CALCULATE CONTRIBUTION FROM CHAPMAN-FERRARO SOURCES + ALL OTHER
C
           EX=DEXP(X/DX)
           EC=EX*CPS
           ES=EX*SPS
           ECZ=EC*Z
           ESZ=ES*Z
           ESZY2=ESZ*Y2
           ESZZ2=ESZ*Z2
           ECZ2=ECZ*Z
           ESY=ES*Y
C
           DER(1,6)=ECZ
           DER(1,7)=ES
           DER(1,8)=ESY*Y
           DER(1,9)=ESZ*Z
           DER(2,10)=ECZ*Y
           DER(2,11)=ESY
           DER(2,12)=ESY*Y2
           DER(2,13)=ESY*Z2
           DER(3,14)=EC
           DER(3,15)=EC*Y2
           DER(3,6)=ECZ2*W1
           DER(3,10)=ECZ2*W5
           DER(3,7)=ESZ*W2
           DER(3,11)=-ESZ
           DER(3,8)=ESZY2*W2
           DER(3,12)=ESZY2*W6
           DER(3,9)=ESZZ2*W3
           DER(3,13)=ESZZ2*W4
C
C  FINALLY, CALCULATE NET EXTERNAL MAGNETIC FIELD COMPONENTS,
C    BUT FIRST OF ALL THOSE FOR C.-F. FIELD:
C
      SX1=AK6*DER(1,6)+AK7*DER(1,7)+AK8*DER(1,8)+AK9*DER(1,9)
      SY1=AK10*DER(2,10)+AK11*DER(2,11)+AK12*DER(2,12)+AK13*DER(2,13)
      SZ1=AK14*DER(3,14)+AK15*DER(3,15)+AK610*ECZ2+AK711*ESZ+AK812
     * *ESZY2+AK913*ESZZ2
       BXCL=AK3*DER(1,3)+AK4*DER(1,4)
       BYCL=AK3*DER(2,3)+AK4*DER(2,4)
       BZCL=AK3*DER(3,3)+AK4*DER(3,4)
       BXT=AK1*DER(1,1)+AK2*DER(1,2)+BXCL +AK16*DER(1,16)+AK17*DER(1,17)
       BYT=AK1*DER(2,1)+AK2*DER(2,2)+BYCL +AK16*DER(2,16)+AK17*DER(2,17)
       BZT=AK1*DER(3,1)+AK2*DER(3,2)+BZCL +AK16*DER(3,16)+AK17*DER(3,17)
       F(1)=BXT+AK5*DER(1,5)+SX1+SXA
       F(2)=BYT+AK5*DER(2,5)+SY1+SYA
       F(3)=BZT+AK5*DER(3,5)+SZ1+SZA
C
       RETURN
       END
c

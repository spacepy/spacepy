!***************************************************************************************************
! Copyright 1995, 1996, N. Tsyganenko
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
C----------------------------------------------------------------------
c
      SUBROUTINE T96_01 (PARMOD,X,Y,Z,BX,BY,BZ)
C
c     RELEASE DATE OF THIS VERSION:   JUNE 22, 1996.
C----------------------------------------------------------------------
C
C  WITH TWO CORRECTIONS, SUGGESTED BY T.SOTIRELIS' COMMENTS (APR.7, 1997)
C
C  (1) A "STRAY "  CLOSING PARENTHESIS WAS REMOVED IN THE S/R   R2_BIRK
C  (2) A 0/0 PROBLEM ON THE Z-AXIS WAS SIDESTEPPED (LINES 44-46 OF THE
c       DOUBLE PRECISION FUNCTION XKSI)
c--------------------------------------------------------------------
C DATA-BASED MODEL CALIBRATED BY (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
C           (2) DST (NANOTESLA),  (3) BYIMF, AND (4) BZIMF (NANOTESLA).
c THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 4 ELEMENTS
c OF THE ARRAY PARMOD(10).
C
C   THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE TILT ANGLE PS (RADIANS),
C AND   X,Y,Z -  GSM POSITION (RE)
C
c   IOPT  IS JUST A DUMMY INPUT PARAMETER, NECESSARY TO MAKE THIS SUBROUTINE
C COMPATIBLE WITH THE NEW RELEASE (APRIL 1996) OF THE TRACING SOFTWARE
C PACKAGE (GEOPACK). IOPT VALUE DOES NOT AFFECT THE OUTPUT FIELD.
c
C
c OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
C            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
C
c  (C) Copr. 1995, 1996, Nikolai A. Tsyganenko, Hughes STX, Code 695, NASA GSFC
c      Greenbelt, MD 20771, USA
c
C                            REFERENCES:
C
C               (1) N.A. TSYGANENKO AND D.P. STERN, A NEW-GENERATION GLOBAL
C           MAGNETOSPHERE FIELD MODEL  , BASED ON SPACECRAFT MAGNETOMETER DATA,
C           ISTP NEWSLETTER, V.6, NO.1, P.21, FEB.1996.
C
c              (2) N.A.TSYGANENKO,  MODELING THE EARTH'S MAGNETOSPHERIC
C           MAGNETIC FIELD CONFINED WITHIN A REALISTIC MAGNETOPAUSE,
C           J.GEOPHYS.RES., V.100, P. 5599, 1995.
C
C              (3) N.A. TSYGANENKO AND M.PEREDO, ANALYTICAL MODELS OF THE
C           MAGNETIC FIELD OF DISK-SHAPED CURRENT SHEETS, J.GEOPHYS.RES.,
C           V.99, P. 199, 1994.
C
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PDYN,DST,BYIMF,BZIMF,PS,X,Y,Z,BX,BY,BZ,QX,QY,QZ,PARMOD(10),
     *   A(9),tilt
c
      COMMON /dip_ang/tilt
      DATA PDYN0,EPS10 /2.d0,3630.7d0/
C
      DATA A/1.162d0,22.344d0,18.50d0,2.602d0,6.903d0,5.287d0,0.5790d0,
     :       0.4462d0,0.7850d0/
C
      DATA  AM0,S0,X00,DSIG/70.d0,1.08d0,5.48d0,0.005d0/
      DATA  DELIMFX,DELIMFY /20.d0,10.d0/
C
       ps=tilt*4.D0*ATAN(1.D0)/180.d0
c       
       	 idummy=iopt
       PDYN=PARMOD(1)
       DST=PARMOD(2)
       BYIMF=PARMOD(3)
       BZIMF=PARMOD(4)
C
       SPS=SIN(PS)
       PPS=PS
C
       DEPR=0.8*DST-13.*SQRT(PDYN)  
c!  DEPR is an estimate of total near-Earth
c                                         depression, based on DST and Pdyn
c                                             (usually, DEPR &lt 0 )
C
C  CALCULATE THE IMF-RELATED QUANTITIES:
C
       Bt=SQRT(BYIMF**2+BZIMF**2)

       IF (BYIMF.EQ.0.d0.AND.BZIMF.EQ.0.d0) THEN
          THETA=0.
          GOTO 1
       ENDIF
C
       THETA=ATAN2(BYIMF,BZIMF)
       IF (THETA.LE.0.D0) THETA=THETA+atan(1.0d0)*8
  1    CT=COS(THETA)
       ST=SIN(THETA)
       EPS=718.5d0*SQRT(Pdyn)*Bt*SIN(THETA/2)
C
       FACTEPS=EPS/EPS10-1
       FACTPD=SQRT(PDYN/PDYN0)-1
C
       RCAMPL=-A(1)*DEPR     
c!   RCAMPL is the amplitude of the ring current
c                  (positive and equal to abs.value of RC depression at origin)
C
       TAMPL2=A(2)+A(3)*FACTPD+A(4)*FACTEPS
       TAMPL3=A(5)+A(6)*FACTPD
       B1AMPL=A(7)+A(8)*FACTEPS
       B2AMPL=20.*B1AMPL  
c! IT IS EQUIVALENT TO ASSUMING THAT THE TOTAL CURRENT
C                           IN THE REGION 2 SYSTEM IS 40% OF THAT IN REGION 1
       RECONN=A(9)
C
       XAPPA=(PDYN/PDYN0)**0.14d0
       XAPPA3=XAPPA**3
       YS=Y*CT-Z*ST
       ZS=Z*CT+Y*ST
C
       FACTIMF=EXP(X/DELIMFX-(YS/DELIMFY)**2)
C
C  CALCULATE THE "IMF" COMPONENTS OUTSIDE THE LAYER  (HENCE BEGIN WITH "O")
C
       OIMFX=0
       OIMFY=RECONN*BYIMF*FACTIMF
       OIMFZ=RECONN*BZIMF*FACTIMF
C
       RIMFAMPL=RECONN*Bt
C
       PPS=PS
       XX=X*XAPPA
       YY=Y*XAPPA
       ZZ=Z*XAPPA
C
C  SCALE AND CALCULATE THE MAGNETOPAUSE PARAMETERS FOR THE INTERPOLATION ACROSS
C   THE BOUNDARY LAYER (THE COORDINATES XX,YY,ZZ  ARE ALREADY SCALED)
C
       X0=X00/XAPPA
       AM=AM0/XAPPA
       RHO2=Y**2+Z**2
       ASQ=AM**2
       XMXM=AM+X-X0
       IF (XMXM.LT.0) XMXM=0 
c! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
       AXX0=XMXM**2
       ARO=ASQ+RHO2
       SIGMA=SQRT((ARO+AXX0+SQRT((ARO+AXX0)**2-4*ASQ*AXX0))/(2*ASQ))
C
C   NOW, THERE ARE THREE POSSIBLE CASES:
C    (1) INSIDE THE MAGNETOSPHERE
C    (2) IN THE BOUNDARY LAYER
C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
C
      IF (SIGMA.LT.S0+DSIG) THEN  
c!  CALCULATE THE T95_06 FIELD (WITH THE
C                                POTENTIAL "PENETRATED" INTERCONNECTION FIELD):

       CALL DIPSHLD(PPS,XX,YY,ZZ,CFX,CFY,CFZ)
       CALL TAILRC96(SPS,XX,YY,ZZ,BXRC,BYRC,BZRC,BXT2,BYT2,BZT2,
     *   BXT3,BYT3,BZT3)
       CALL BIRK1TOT_02(PPS,XX,YY,ZZ,R1X,R1Y,R1Z)
       CALL BIRK2TOT_02(PPS,XX,YY,ZZ,R2X,R2Y,R2Z)
       CALL INTERCON(XX,YS*XAPPA,ZS*XAPPA,RIMFX,RIMFYS,RIMFZS)
       RIMFY=RIMFYS*CT+RIMFZS*ST
       RIMFZ=RIMFZS*CT-RIMFYS*ST
C
      FX=CFX*XAPPA3+RCAMPL*BXRC +TAMPL2*BXT2+TAMPL3*BXT3
     *  +B1AMPL*R1X +B2AMPL*R2X +RIMFAMPL*RIMFX
      FY=CFY*XAPPA3+RCAMPL*BYRC +TAMPL2*BYT2+TAMPL3*BYT3
     *  +B1AMPL*R1Y +B2AMPL*R2Y +RIMFAMPL*RIMFY
      FZ=CFZ*XAPPA3+RCAMPL*BZRC +TAMPL2*BZT2+TAMPL3*BZT3
     *  +B1AMPL*R1Z +B2AMPL*R2Z +RIMFAMPL*RIMFZ
C
C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
C
       IF (SIGMA.LT.S0-DSIG) THEN
         BX=FX
         BY=FY
         BZ=FZ
                        ELSE  
c!  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
C                                         THE INTERPOLATION REGION
       FINT=0.5d0*(1-(SIGMA-S0)/DSIG)
       FEXT=0.5d0*(1+(SIGMA-S0)/DSIG)
C
       CALL DIPOLE(PS,X,Y,Z,QX,QY,QZ)
       BX=(FX+QX)*FINT+OIMFX*FEXT -QX
       BY=(FY+QY)*FINT+OIMFY*FEXT -QY
       BZ=(FZ+QZ)*FINT+OIMFZ*FEXT -QZ
c
        ENDIF  
c!   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
C                      POSSIBILITY IS NOW THE CASE (3):
         ELSE
                CALL DIPOLE(PS,X,Y,Z,QX,QY,QZ)
                BX=OIMFX-QX
                BY=OIMFY-QY
                BZ=OIMFZ-QZ
         ENDIF
C
       RETURN
       END
C=====================================================================

      SUBROUTINE DIPSHLD(PS,X,Y,Z,BX,BY,BZ)
C
C   CALCULATES GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD DUE TO
C    SHIELDING OF THE EARTH'S DIPOLE ONLY
C
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION A1(12),A2(12)
      DATA A1 /.24777d0,-27.003d0,-.46815d0,7.0637d0,-1.5918d0,
     :   -.90317d-01,57.522d0,13.757d0,2.0100d0,10.458d0,4.5798d0,
     :   2.1695d0/
      DATA A2/-.65385d0,-18.061d0,-.40457d0,-5.0995d0,1.2846d0,
     :   .78231d-01,39.592d0,13.291d0,1.9970d0,10.062d0,4.5140d0,
     :   2.1558d0/
C
          CPS=DCOS(PS)
          SPS=DSIN(PS)
          CALL CYLHARM(A1,X,Y,Z,HX,HY,HZ)
          CALL CYLHAR1(A2,X,Y,Z,FX,FY,FZ)
C
          BX=HX*CPS+FX*SPS
          BY=HY*CPS+FY*SPS
          BZ=HZ*CPS+FZ*SPS
        RETURN
       END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PERPENDICULAR DIPOLE
C
         SUBROUTINE  CYLHARM( A, X,Y,Z, BX,BY,BZ)
C
C
C   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
C
C   An approximation for the Chapman-Ferraro field by a sum of 6 cylin-
c   drical harmonics (see pp. 97-113 in the brown GSFC notebook #1)
c
C      Description of parameters:
C
C  A   - input vector containing model parameters;
C  X,Y,Z   -  input GSM coordinates
C  BX,BY,BZ - output GSM components of the shielding field
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical harmonic
c       terms.
c  The 6 nonlinear parameters A(7)-A(12) are the corresponding scale lengths
C       for each term (see GSFC brown notebook).
c
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
	 IMPLICIT  REAL * 8  (A - H, O - Z)
C
	 DIMENSION  A(12)
C
           RHO=DSQRT(Y**2+Z**2)
            IF (RHO.LT.1.D-8) THEN
               SINFI=1.D0
               COSFI=0.D0
               RHO=1.D-8
               GOTO 1
            ENDIF
C
           SINFI=Z/RHO
           COSFI=Y/RHO
  1        SINFI2=SINFI**2
           SI2CO2=SINFI2-COSFI**2
C
             BX=0.D0
             BY=0.D0
             BZ=0.D0
C
	   DO 11 I=1,3
             DZETA=RHO/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(X/A(I+6))
             BX=BX-A(I)*XJ1*XEXP*SINFI
             BY=BY+A(I)*(2.D0*XJ1/DZETA-XJ0)*XEXP*SINFI*COSFI
             BZ=BZ+A(I)*(XJ1/DZETA*SI2CO2-XJ0*SINFI2)*XEXP
   11        CONTINUE
c
	   DO 12 I=4,6
             DZETA=RHO/A(I+6)
             XKSI=X/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(XKSI)
             BRHO=(XKSI*XJ0-(DZETA**2+XKSI-1.D0)*XJ1/DZETA)*XEXP*SINFI
             BPHI=(XJ0+XJ1/DZETA*(XKSI-1.D0))*XEXP*COSFI
             BX=BX+A(I)*(DZETA*XJ0+XKSI*XJ1)*XEXP*SINFI
             BY=BY+A(I)*(BRHO*COSFI-BPHI*SINFI)
             BZ=BZ+A(I)*(BRHO*SINFI+BPHI*COSFI)
   12        CONTINUE
C
c
         RETURN
	 END
C
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PARALLEL DIPOLE
C
         SUBROUTINE  CYLHAR1(A, X,Y,Z, BX,BY,BZ)
C
C
C   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
C
C   An approximation of the Chapman-Ferraro field by a sum of 6 cylin-
c   drical harmonics (see pages 97-113 in the brown GSFC notebook #1)
c
C      Description of parameters:
C
C  A   - input vector containing model parameters;
C  X,Y,Z - input GSM coordinates,
C  BX,BY,BZ - output GSM components of the shielding field
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C      The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical
c  harmonic terms.
c      The 6 nonlinear parameters A(7)-A(12) are the corresponding scale
c  lengths for each term (see GSFC brown notebook).
c
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
	 IMPLICIT  REAL * 8  (A - H, O - Z)
C
	 DIMENSION  A(12)
C
           RHO=DSQRT(Y**2+Z**2)
            IF (RHO.LT.1.D-10) THEN
               SINFI=1.D0
               COSFI=0.D0
               GOTO 1
            ENDIF
C
           SINFI=Z/RHO
           COSFI=Y/RHO
C
   1      BX=0.D0
          BY=0.D0
          BZ=0.D0
C
             DO 11 I=1,3
             DZETA=RHO/A(I+6)
             XKSI=X/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(XKSI)
             BRHO=XJ1*XEXP
             BX=BX-A(I)*XJ0*XEXP
             BY=BY+A(I)*BRHO*COSFI
             BZ=BZ+A(I)*BRHO*SINFI
   11        CONTINUE
c
	   DO 12 I=4,6
             DZETA=RHO/A(I+6)
             XKSI=X/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(XKSI)
             BRHO=(DZETA*XJ0+XKSI*XJ1)*XEXP
             BX=BX+A(I)*(DZETA*XJ1-XJ0*(XKSI+1.D0))*XEXP
             BY=BY+A(I)*BRHO*COSFI
             BZ=BZ+A(I)*BRHO*SINFI
   12        CONTINUE
C
         RETURN
	 END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      REAL*8 FUNCTION BES(X,K)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF (K.EQ.0) THEN
        BES=BES0(X)
        RETURN
      ENDIF
C
      IF (K.EQ.1) THEN
        BES=BES1(X)
        RETURN
      ENDIF
C
      IF (X.EQ.0.D0) THEN
        BES=0.D0
        RETURN
      ENDIF
C
      G=2.D0/X
      IF (X.LE.K) GOTO 10
C
      N=1
      XJN=BES1(X)
      XJNM1=BES0(X)
C
  1   XJNP1=G*N*XJN-XJNM1
      N=N+1
      IF (N.LT.K) GOTO 2
      BES=XJNP1
      RETURN
C
 2    XJNM1=XJN
      XJN=XJNP1
      GOTO 1
C
 10   N=24
      XJN=1.D0
      XJNP1=0.D0
      SUM=0.D0
C
  3   IF (MOD(N,2).EQ.0) SUM=SUM+XJN
      XJNM1=G*N*XJN-XJNP1
      N=N-1
C
      XJNP1=XJN
      XJN=XJNM1
      IF (N.EQ.K) BES=XJN
C
      IF (DABS(XJN).GT.1.D5) THEN
        XJNP1=XJNP1*1.D-5
        XJN=XJN*1.D-5
        SUM=SUM*1.D-5
        IF (N.LE.K) BES=BES*1.D-5
      ENDIF
C
      IF (N.EQ.0) GOTO 4
      GOTO 3
C
  4   SUM=XJN+2.D0*SUM
      BES=BES/SUM
      RETURN
      END
c-------------------------------------------------------------------
c
       REAL*8 FUNCTION BES0(X)
C
        IMPLICIT REAL*8 (A-H,O-Z)
C
        IF (DABS(X).LT.3.D0) THEN
          X32=(X/3.D0)**2
          BES0=1.D0-X32*(2.2499997D0-X32*(1.2656208D0-X32*
     *    (0.3163866D0-X32*(0.0444479D0-X32*(0.0039444D0
     *     -X32*0.00021D0)))))
        ELSE
        XD3=3.D0/X
        F0=0.79788456D0-XD3*(0.00000077D0+XD3*(0.00552740D0+XD3*
     *   (0.00009512D0-XD3*(0.00137237D0-XD3*(0.00072805D0
     *   -XD3*0.00014476D0)))))
        T0=X-0.78539816D0-XD3*(0.04166397D0+XD3*(0.00003954D0-XD3*
     *   (0.00262573D0-XD3*(0.00054125D0+XD3*(0.00029333D0
     *   -XD3*0.00013558D0)))))
        BES0=F0/DSQRT(X)*DCOS(T0)
       ENDIF
       RETURN
       END
c
c--------------------------------------------------------------------------
c
       REAL*8 FUNCTION BES1(X)
C
        IMPLICIT REAL*8 (A-H,O-Z)
C
       IF (DABS(X).LT.3.D0) THEN
        X32=(X/3.D0)**2
        BES1XM1=0.5D0-X32*(0.56249985D0-X32*(0.21093573D0-X32*
     *  (0.03954289D0-X32*(0.00443319D0-X32*(0.00031761D0
     *  -X32*0.00001109D0)))))
         BES1=BES1XM1*X
       ELSE
        XD3=3.D0/X
        F1=0.79788456D0+XD3*(0.00000156D0+XD3*(0.01659667D0+XD3*
     *   (0.00017105D0-XD3*(0.00249511D0-XD3*(0.00113653D0
     *   -XD3*0.00020033D0)))))
        T1=X-2.35619449D0+XD3*(0.12499612D0+XD3*(0.0000565D0-XD3*
     *   (0.00637879D0-XD3*(0.00074348D0+XD3*(0.00079824D0
     *   -XD3*0.00029166D0)))))
        BES1=F1/DSQRT(X)*DCOS(T1)
       ENDIF
       RETURN
       END
C------------------------------------------------------------
C
         SUBROUTINE INTERCON(X,Y,Z,BX,BY,BZ)
C
C      Calculates the potential interconnection field inside the magnetosphere,
c  corresponding to  DELTA_X = 20Re and DELTA_Y = 10Re (NB#3, p.90, 6/6/1996).
C  The position (X,Y,Z) and field components BX,BY,BZ are given in the rotated
c   coordinate system, in which the Z-axis is always directed along the BzIMF
c   (i.e. rotated by the IMF clock angle Theta)
C   It is also assumed that the IMF Bt=1, so that the components should be
c     (i) multiplied by the actual Bt, and
c     (ii) transformed to standard GSM coords by rotating back around X axis
c              by the angle -Theta.
c
C      Description of parameters:
C
C     X,Y,Z -   GSM POSITION
C      BX,BY,BZ - INTERCONNECTION FIELD COMPONENTS INSIDE THE MAGNETOSPHERE
C        OF A STANDARD SIZE (TO TAKE INTO ACCOUNT EFFECTS OF PRESSURE CHANGES,
C         APPLY THE SCALING TRANSFORMATION)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C     The 9 linear parameters are amplitudes of the "cartesian" harmonics
c     The 6 nonlinear parameters are the scales Pi and Ri entering
c    the arguments of exponents, sines, and cosines in the 9 "Cartesian"
c       harmonics (3+3)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
	 IMPLICIT  REAL * 8  (A - H, O - Z)
C
        DIMENSION A(15),RP(3),RR(3),P(3),R(3)
C
      SAVE M,RP,RR
C
      DATA A/-8.411078731d0,5932254.951d0,-9073284.93d0,-11.68794634d0,
     * 6027598.824d0,-9218378.368d0,-6.508798398d0,-11824.42793d0,
     : 18015.66212d0,
     * 7.99754043d0,13.9669886d0,90.24475036d0,16.75728834d0,
     : 1015.645781d0,
     * 1553.493216d0/
C
        DATA M/0/
C
        IF (M.NE.0) GOTO 111
        M=1
C
         P(1)=A(10)
         P(2)=A(11)
         P(3)=A(12)
         R(1)=A(13)
         R(2)=A(14)
         R(3)=A(15)
C
C
           DO 11 I=1,3
            RP(I)=1.D0/P(I)
  11        RR(I)=1.D0/R(I)
C
  111   CONTINUE
C
            L=0
C
               BX=0.d0
               BY=0.d0
               BZ=0.d0
C
c        "PERPENDICULAR" KIND OF SYMMETRY ONLY
C
               DO 2 I=1,3
                  CYPI=DCOS(Y*RP(I))
                  SYPI=DSIN(Y*RP(I))
C
                DO 2 K=1,3
                   SZRK=DSIN(Z*RR(K))
                   CZRK=DCOS(Z*RR(K))
                     SQPR=DSQRT(RP(I)**2+RR(K)**2)
                      EPR=DEXP(X*SQPR)
C
                     HX=-SQPR*EPR*CYPI*SZRK
                     HY=RP(I)*EPR*SYPI*SZRK
                     HZ=-RR(K)*EPR*CYPI*CZRK
             L=L+1
c
          BX=BX+A(L)*HX
          BY=BY+A(L)*HY
          BZ=BZ+A(L)*HZ
  2   CONTINUE
C
      RETURN
      END

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TAILRC96(SPS,X,Y,Z,BXRC,BYRC,BZRC,BXT2,BYT2,BZT2,
     *   BXT3,BYT3,BZT3)
c
c  COMPUTES THE COMPONENTS OF THE FIELD OF THE MODEL RING CURRENT AND THREE
c                   TAIL MODES WITH UNIT AMPLITUDES
C      (FOR THE RING CURRENT, IT MEANS THE DISTURBANCE OF Bz=-1nT AT ORIGIN,
C   AND FOR THE TAIL MODES IT MEANS MAXIMAL BX JUST ABOVE THE SHEET EQUAL 1 nT.
C
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION ARC(48),ATAIL2(48),ATAIL3(48)
         COMMON /WARP/ CPSS,SPSS,DPSRR,RPS,WARP,D,XS,ZS,DXSX,DXSY,DXSZ,
     *   DZSX,DZSY,DZSZ,DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW
C
         DATA ARC/-3.087699646d0,3.516259114d0,18.81380577d0,
     :  -13.95772338d0,-5.497076303d0,0.1712890838d0,2.392629189d0,
     :  -2.728020808d0,-14.79349936d0,11.08738083d0,4.388174084d0,
     :  0.2492163197d-01,0.7030375685d0,-.7966023165d0,
     :  -3.835041334d0,2.642228681d0,-0.2405352424d0,-0.7297705678d0,
     :  -0.3680255045d0,0.1333685557d0,2.795140897d0,-1.078379954d0,
     :  0.8014028630d0,0.1245825565d0,0.6149982835d0,-0.2207267314d0,
     :  -4.424578723d0,1.730471572d0,-1.716313926d0,-0.2306302941d0,
     :  -0.2450342688d0,0.8617173961d-01,1.54697858d0,
     :  -0.6569391113d0,-0.6537525353d0,0.2079417515d0,12.75434981d0,
     :  11.37659788d0,636.4346279d0,1.752483754d0,3.604231143d0,
     :  12.83078674d0,7.412066636d0,9.434625736d0,676.7557193d0,
     :  1.701162737d0,3.580307144d0,14.64298662d0/
C
        DATA ATAIL2/.8747515218d0,-.9116821411d0,2.209365387d0,
     : -2.159059518d0,-7.059828867d0,5.924671028d0,-1.916935691d0,
     : 1.996707344d0,-3.877101873d0,3.947666061d0,11.38715899d0,
     : -8.343210833d0,1.194109867d0,-1.244316975d0,3.73895491d0,
     : -4.406522465d0,-20.66884863d0,3.020952989d0,.2189908481d0,
     : -.09942543549d0,-.927225562d0,.1555224669d0,.6994137909d0,
     : -.08111721003d0,-.7565493881d0,.4686588792d0,4.266058082d0,
     : -.3717470262d0,-3.920787807d0,.02298569870d0,.7039506341d0,
     : -.5498352719d0,-6.675140817d0,.8279283559d0,-2.234773608d0,
     : -1.622656137d0,5.187666221d0,6.802472048d0,39.13543412d0,
     : 2.784722096d0,6.979576616d0,25.71716760d0,4.495005873d0,
     : 8.068408272d0,93.47887103d0,4.158030104d0,9.313492566d0,
     :57.18240483d0/
C
        DATA ATAIL3/-19091.95061d0,-3011.613928d0,20582.16203d0,
     : 4242.918430d0,-2377.091102d0,-1504.820043d0,19884.04650d0,
     : 2725.150544d0,-21389.04845d0,-3990.475093d0,2401.610097d0,
     : 1548.171792d0,-946.5493963d0,490.1528941d0,986.9156625d0,
     : -489.3265930d0,-67.99278499d0,8.711175710d0,-45.15734260d0,
     : -10.76106500d0,210.7927312d0,11.41764141d0,-178.0262808d0,
     : .7558830028d0,339.3806753d0,9.904695974d0,69.50583193d0,
     : -118.0271581d0,22.85935896d0,45.91014857d0,-425.6607164d0,
     : 15.47250738d0,118.2988915d0,65.58594397d0,-201.4478068d0,
     : -14.57062940d0,19.69877970d0,20.30095680d0,86.45407420d0,
     : 22.50403727d0,23.41617329d0,48.48140573d0,24.61031329d0,
     : 123.5395974d0,223.5367692d0,39.50824342d0,65.83385762d0,
     : 266.2948657d0/
C
       DATA RH,DR,G,D0,DELTADY/9.d0,4.d0,10.d0,2.d0,10.d0/
C
C   TO ECONOMIZE THE CODE, WE FIRST CALCULATE COMMON VARIABLES, WHICH ARE
C      THE SAME FOR ALL MODES, AND PUT THEM IN THE COMMON-BLOCK /WARP/
C
       DR2=DR*DR
       C11=DSQRT((1.D0+RH)**2+DR2)
       C12=DSQRT((1.D0-RH)**2+DR2)
       C1=C11-C12
       SPSC1=SPS/C1
       RPS=0.5d0*(C11+C12)*SPS  
c!  THIS IS THE SHIFT OF OF THE SHEET WITH RESPECT
C                            TO GSM EQ.PLANE FOR THE 3RD (ASYMPTOTIC) TAIL MODE
C
        R=DSQRT(X*X+Y*Y+Z*Z)
        SQ1=DSQRT((R+RH)**2+DR2)
        SQ2=DSQRT((R-RH)**2+DR2)
        C=SQ1-SQ2
        CS=(R+RH)/SQ1-(R-RH)/SQ2
        SPSS=SPSC1/R*C
        CPSS=DSQRT(1.D0-SPSS**2)
        DPSRR=SPS/(R*R)*(CS*R-C)/DSQRT((R*C1)**2-(C*SPS)**2)
C
        WFAC=Y/(Y**4+1.D4)   
c!   WARPING
        W=WFAC*Y**3
        WS=4.D4*Y*WFAC**2
        WARP=G*SPS*W
        XS=X*CPSS-Z*SPSS
        ZSWW=Z*CPSS+X*SPSS  
c! "WW" MEANS "WITHOUT Y-Z WARPING" (IN X-Z ONLY)
        ZS=ZSWW +WARP

        DXSX=CPSS-X*ZSWW*DPSRR
        DXSY=-Y*ZSWW*DPSRR
        DXSZ=-SPSS-Z*ZSWW*DPSRR
        DZSX=SPSS+X*XS*DPSRR
        DZSY=XS*Y*DPSRR  +G*SPS*WS  
c!  THE LAST TERM IS FOR THE Y-Z WARP
        DZSZ=CPSS+XS*Z*DPSRR        
c!      (TAIL MODES ONLY)

        D=D0+DELTADY*(Y/20.D0)**2   
c!  SHEET HALF-THICKNESS FOR THE TAIL MODES
        DDDY=DELTADY*Y*0.005D0      
c!  (THICKENS TO FLANKS, BUT NO VARIATION
C                                         ALONG X, IN CONTRAST TO RING CURRENT)
C
        DZETAS=DSQRT(ZS**2+D**2)  
c!  THIS IS THE SAME SIMPLE WAY TO SPREAD
C                                        OUT THE SHEET, AS THAT USED IN T89
        DDZETADX=ZS*DZSX/DZETAS
        DDZETADY=(ZS*DZSY+D*DDDY)/DZETAS
        DDZETADZ=ZS*DZSZ/DZETAS
C
        CALL SHLCAR3X3(ARC,X,Y,Z,SPS,WX,WY,WZ)
        CALL RINGCURR96(X,Y,Z,HX,HY,HZ)
        BXRC=WX+HX
        BYRC=WY+HY
        BZRC=WZ+HZ
C
        CALL SHLCAR3X3(ATAIL2,X,Y,Z,SPS,WX,WY,WZ)
        CALL TAILDISK(X,Y,Z,HX,HY,HZ)
        BXT2=WX+HX
        BYT2=WY+HY
        BZT2=WZ+HZ
C
        CALL SHLCAR3X3(ATAIL3,X,Y,Z,SPS,WX,WY,WZ)
        CALL TAIL87(X,Z,HX,HZ)
        BXT3=WX+HX
        BYT3=WY
        BZT3=WZ+HZ
C
        RETURN
        END
C
c********************************************************************
C
        SUBROUTINE RINGCURR96(X,Y,Z,BX,BY,BZ)
c
c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE RING CURRENT FIELD,
C        SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
C          DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN THE
C           PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996),
C            INSTEAD OF SHEARING IT IN THE SPIRIT OF THE T89 TAIL MODEL.
C
C          IN  ADDITION, INSTEAD OF 7 TERMS FOR THE RING CURRENT MODEL, WE USE
C             NOW ONLY 2 TERMS;  THIS SIMPLIFICATION ALSO GIVES RISE TO AN
C                EASTWARD RING CURRENT LOCATED EARTHWARD FROM THE MAIN ONE,
C                  IN LINE WITH WHAT IS ACTUALLY OBSERVED
C
C             FOR DETAILS, SEE NB #3, PAGES 70-73
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION F(2),BETA(2)
        COMMON /WARP/ CPSS,SPSS,DPSRR, XNEXT(3),XS,ZSWARPED,DXSX,DXSY,
     *   DXSZ,DZSX,DZSYWARPED,DZSZ,OTHER(4),ZS  
c!  ZS HERE IS WITHOUT Y-Z WARP
C

      DATA D0,DELTADX,XD,XLDX /2.d0,0.d0,0.,4.d0/  
c!  ACHTUNG !!  THE RC IS NOW
C                                            COMPLETELY SYMMETRIC (DELTADX=0)

C
        DATA F,BETA /569.895366D0,-1603.386993D0,2.722188D0,3.766875D0/
C
C  THE ORIGINAL VALUES OF F(I) WERE MULTIPLIED BY BETA(I) (TO REDUCE THE
C     NUMBER OF MULTIPLICATIONS BELOW)  AND BY THE FACTOR -0.43, NORMALIZING
C      THE DISTURBANCE AT ORIGIN  TO  B=-1nT
C
           DZSY=XS*Y*DPSRR  
c! NO WARPING IN THE Y-Z PLANE (ALONG X ONLY), AND
C                         THIS IS WHY WE DO NOT USE  DZSY FROM THE COMMON-BLOCK
           XXD=X-XD
           FDX=0.5D0*(1.D0+XXD/DSQRT(XXD**2+XLDX**2))
           DDDX=DELTADX*0.5D0*XLDX**2/DSQRT(XXD**2+XLDX**2)**3
           D=D0+DELTADX*FDX

           DZETAS=DSQRT(ZS**2+D**2)  
c!  THIS IS THE SAME SIMPLE WAY TO SPREAD
C                                        OUT THE SHEET, AS THAT USED IN T89
           RHOS=DSQRT(XS**2+Y**2)
           DDZETADX=(ZS*DZSX+D*DDDX)/DZETAS
           DDZETADY=ZS*DZSY/DZETAS
           DDZETADZ=ZS*DZSZ/DZETAS
         IF (RHOS.LT.1.D-5) THEN
            DRHOSDX=0.D0
            DRHOSDY=DSIGN(1.D0,Y)
            DRHOSDZ=0.D0
          ELSE
           DRHOSDX=XS*DXSX/RHOS
           DRHOSDY=(XS*DXSY+Y)/RHOS
           DRHOSDZ=XS*DXSZ/RHOS
         ENDIF
C
           BX=0.D0
           BY=0.D0
           BZ=0.D0
C
           DO 1 I=1,2
C
           BI=BETA(I)
C
           S1=DSQRT((DZETAS+BI)**2+(RHOS+BI)**2)
           S2=DSQRT((DZETAS+BI)**2+(RHOS-BI)**2)
           DS1DDZ=(DZETAS+BI)/S1
           DS2DDZ=(DZETAS+BI)/S2
           DS1DRHOS=(RHOS+BI)/S1
           DS2DRHOS=(RHOS-BI)/S2
C
           DS1DX=DS1DDZ*DDZETADX+DS1DRHOS*DRHOSDX
           DS1DY=DS1DDZ*DDZETADY+DS1DRHOS*DRHOSDY
           DS1DZ=DS1DDZ*DDZETADZ+DS1DRHOS*DRHOSDZ
C
           DS2DX=DS2DDZ*DDZETADX+DS2DRHOS*DRHOSDX
           DS2DY=DS2DDZ*DDZETADY+DS2DRHOS*DRHOSDY
           DS2DZ=DS2DDZ*DDZETADZ+DS2DRHOS*DRHOSDZ
C
           S1TS2=S1*S2
           S1PS2=S1+S2
           S1PS2SQ=S1PS2**2
           FAC1=DSQRT(S1PS2SQ-(2.D0*BI)**2)
           AS=FAC1/(S1TS2*S1PS2SQ)
           TERM1=1.D0/(S1TS2*S1PS2*FAC1)
           FAC2=AS/S1PS2SQ
           DASDS1=TERM1-FAC2/S1*(S2*S2+S1*(3.D0*S1+4.D0*S2))
           DASDS2=TERM1-FAC2/S2*(S1*S1+S2*(3.D0*S2+4.D0*S1))
C
           DASDX=DASDS1*DS1DX+DASDS2*DS2DX
           DASDY=DASDS1*DS1DY+DASDS2*DS2DY
           DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ
C
      BX=BX+F(I)*((2.D0*AS+Y*DASDY)*SPSS-XS*DASDZ
     *   +AS*DPSRR*(Y**2*CPSS+Z*ZS))
      BY=BY-F(I)*Y*(AS*DPSRR*XS+DASDZ*CPSS+DASDX*SPSS)
  1   BZ=BZ+F(I)*((2.D0*AS+Y*DASDY)*CPSS+XS*DASDX
     *   -AS*DPSRR*(X*ZS+Y**2*SPSS))
C
       RETURN
       END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
         SUBROUTINE TAILDISK(X,Y,Z,BX,BY,BZ)
C
c
c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT FIELD,
C        SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
C          DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN OUR
C           PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996)
C            INSTEAD OF SHEARING IT IN THE SPIRIT OF T89 TAIL MODEL.
C
C          IN  ADDITION, INSTEAD OF 8 TERMS FOR THE TAIL CURRENT MODEL, WE USE
C           NOW ONLY 4 TERMS
C
C             FOR DETAILS, SEE NB #3, PAGES 74-
C
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION F(4),BETA(4)
         COMMON /WARP/ CPSS,SPSS,DPSRR,XNEXT(3),XS,ZS,DXSX,DXSY,DXSZ,
     *    OTHER(3),DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW
C
         DATA XSHIFT /4.5d0/
C
         DATA F,BETA
     * / -745796.7338D0,1176470.141D0,-444610.529D0,-57508.01028D0,
     *   7.9250000D0,8.0850000D0,8.4712500D0,27.89500D0/
c
c  here original F(I) are multiplied by BETA(I), to economize
c    calculations
C
           RHOS=DSQRT((XS-XSHIFT)**2+Y**2)
         IF (RHOS.LT.1.D-5) THEN
            DRHOSDX=0.D0
            DRHOSDY=DSIGN(1.D0,Y)
            DRHOSDZ=0.D0
         ELSE
           DRHOSDX=(XS-XSHIFT)*DXSX/RHOS
           DRHOSDY=((XS-XSHIFT)*DXSY+Y)/RHOS
           DRHOSDZ=(XS-XSHIFT)*DXSZ/RHOS
         ENDIF
C
           BX=0.D0
           BY=0.D0
           BZ=0.D0
C
           DO 1 I=1,4
C
           BI=BETA(I)
C
           S1=DSQRT((DZETAS+BI)**2+(RHOS+BI)**2)
           S2=DSQRT((DZETAS+BI)**2+(RHOS-BI)**2)
           DS1DDZ=(DZETAS+BI)/S1
           DS2DDZ=(DZETAS+BI)/S2
           DS1DRHOS=(RHOS+BI)/S1
           DS2DRHOS=(RHOS-BI)/S2
C
           DS1DX=DS1DDZ*DDZETADX+DS1DRHOS*DRHOSDX
           DS1DY=DS1DDZ*DDZETADY+DS1DRHOS*DRHOSDY
           DS1DZ=DS1DDZ*DDZETADZ+DS1DRHOS*DRHOSDZ
C
           DS2DX=DS2DDZ*DDZETADX+DS2DRHOS*DRHOSDX
           DS2DY=DS2DDZ*DDZETADY+DS2DRHOS*DRHOSDY
           DS2DZ=DS2DDZ*DDZETADZ+DS2DRHOS*DRHOSDZ
C
           S1TS2=S1*S2
           S1PS2=S1+S2
           S1PS2SQ=S1PS2**2
           FAC1=DSQRT(S1PS2SQ-(2.D0*BI)**2)
           AS=FAC1/(S1TS2*S1PS2SQ)
           TERM1=1.D0/(S1TS2*S1PS2*FAC1)
           FAC2=AS/S1PS2SQ
           DASDS1=TERM1-FAC2/S1*(S2*S2+S1*(3.D0*S1+4.D0*S2))
           DASDS2=TERM1-FAC2/S2*(S1*S1+S2*(3.D0*S2+4.D0*S1))
C
           DASDX=DASDS1*DS1DX+DASDS2*DS2DX
           DASDY=DASDS1*DS1DY+DASDS2*DS2DY
           DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ
C
      BX=BX+F(I)*((2.D0*AS+Y*DASDY)*SPSS-(XS-XSHIFT)*DASDZ
     *   +AS*DPSRR*(Y**2*CPSS+Z*ZSWW))
C
      BY=BY-F(I)*Y*(AS*DPSRR*XS+DASDZ*CPSS+DASDX*SPSS)
  1   BZ=BZ+F(I)*((2.D0*AS+Y*DASDY)*CPSS+(XS-XSHIFT)*DASDX
     *   -AS*DPSRR*(X*ZSWW+Y**2*SPSS))

       RETURN
       END

C-------------------------------------------------------------------------
C
      SUBROUTINE TAIL87(X,Z,BX,BZ)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /WARP/ FIRST(3), RPS,WARP,D, OTHER(13)
C
C      'LONG' VERSION OF THE 1987 TAIL MAGNETIC FIELD MODEL
C              (N.A.TSYGANENKO, PLANET. SPACE SCI., V.35, P.1347, 1987)
C
C      D   IS THE Y-DEPENDENT SHEET HALF-THICKNESS (INCREASING TOWARDS FLANKS)
C      RPS  IS THE TILT-DEPENDENT SHIFT OF THE SHEET IN THE Z-DIRECTION,
C           CORRESPONDING TO THE ASYMPTOTIC HINGING DISTANCE, DEFINED IN THE
C           MAIN SUBROUTINE (TAILRC96) FROM THE PARAMETERS RH AND DR OF THE
C           T96-TYPE MODULE, AND
C      WARP  IS THE BENDING OF THE SHEET FLANKS IN THE Z-DIRECTION, DIRECTED
C           OPPOSITE TO RPS, AND INCREASING WITH DIPOLE TILT AND |Y|
C

        DATA DD/3.d0/
C
      DATA HPI,RT,XN,X1,X2,B0,B1,B2,XN21,XNR,ADLN
     * /1.5707963d0,40.d0,-10.d0,
     * -1.261d0,-0.663d0,0.391734d0,5.89715d0,24.6833d0,76.37d0,
     : -0.1071d0,0.13238005d0/
C                !!!   THESE ARE NEW VALUES OF  X1, X2, B0, B1, B2,
C                       CORRESPONDING TO TSCALE=1, INSTEAD OF TSCALE=0.6
C
C  THE ABOVE QUANTITIES WERE DEFINED AS FOLLOWS:------------------------
C       HPI=PI/2
C       RT=40.      !  Z-POSITION OF UPPER AND LOWER ADDITIONAL SHEETS
C       XN=-10.     !  INNER EDGE POSITION
C
C       TSCALE=1  !  SCALING FACTOR, DEFINING THE RATE OF INCREASE OF THE
C                       CURRENT DENSITY TAILWARDS
C
c  ATTENTION !  NOW I HAVE CHANGED TSCALE TO:  TSCALE=1.0, INSTEAD OF 0.6
c                  OF THE PREVIOUS VERSION
c
C       B0=0.391734
C       B1=5.89715 *TSCALE
C       B2=24.6833 *TSCALE**2
C
C    HERE ORIGINAL VALUES OF THE MODE AMPLITUDES (P.77, NB#3) WERE NORMALIZED
C      SO THAT ASYMPTOTIC  BX=1  AT X=-200RE
C
C      X1=(4.589  -5.85) *TSCALE -(TSCALE-1.)*XN ! NONLINEAR PARAMETERS OF THE
C                                                         CURRENT FUNCTION
C      X2=(5.187  -5.85) *TSCALE -(TSCALE-1.)*XN
c
c
C      XN21=(XN-X1)**2
C      XNR=1./(XN-X2)
C      ADLN=-DLOG(XNR**2*XN21)
C
C---------------------------------------------------------------
C
      ZS=Z -RPS +WARP
      ZP=Z-RT
      ZM=Z+RT
C
      XNX=XN-X
      XNX2=XNX**2
      XC1=X-X1
      XC2=X-X2
      XC22=XC2**2
      XR2=XC2*XNR
      XC12=XC1**2
      D2=DD**2    
c!  SQUARE OF THE TOTAL HALFTHICKNESS (DD=3Re for this mode)
      B20=ZS**2+D2
      B2P=ZP**2+D2
      B2M=ZM**2+D2
      B=DSQRT(B20)
      BP=DSQRT(B2P)
      BM=DSQRT(B2M)
      XA1=XC12+B20
      XAP1=XC12+B2P
      XAM1=XC12+B2M
      XA2=1.d0/(XC22+B20)
      XAP2=1.d0/(XC22+B2P)
      XAM2=1.d0/(XC22+B2M)
      XNA=XNX2+B20
      XNAP=XNX2+B2P
      XNAM=XNX2+B2M
      F=B20-XC22
      FP=B2P-XC22
      FM=B2M-XC22
      XLN1=DLOG(XN21/XNA)
      XLNP1=DLOG(XN21/XNAP)
      XLNM1=DLOG(XN21/XNAM)
      XLN2=XLN1+ADLN
      XLNP2=XLNP1+ADLN
      XLNM2=XLNM1+ADLN
      ALN=0.25d0*(XLNP1+XLNM1-2.d0*XLN1)
      S0=(DATAN(XNX/B)+HPI)/B
      S0P=(DATAN(XNX/BP)+HPI)/BP
      S0M=(DATAN(XNX/BM)+HPI)/BM
      S1=(XLN1*.5d0+XC1*S0)/XA1
      S1P=(XLNP1*.5d0+XC1*S0P)/XAP1
      S1M=(XLNM1*.5d0+XC1*S0M)/XAM1
      S2=(XC2*XA2*XLN2-XNR-F*XA2*S0)*XA2
      S2P=(XC2*XAP2*XLNP2-XNR-FP*XAP2*S0P)*XAP2
      S2M=(XC2*XAM2*XLNM2-XNR-FM*XAM2*S0M)*XAM2
      G1=(B20*S0-0.5d0*XC1*XLN1)/XA1
      G1P=(B2P*S0P-0.5d0*XC1*XLNP1)/XAP1
      G1M=(B2M*S0M-0.5d0*XC1*XLNM1)/XAM1
      G2=((0.5d0*F*XLN2+2.d0*S0*B20*XC2)*XA2+XR2)*XA2
      G2P=((0.5d0*FP*XLNP2+2.d0*S0P*B2P*XC2)*XAP2+XR2)*XAP2
      G2M=((0.5d0*FM*XLNM2+2.d0*S0M*B2M*XC2)*XAM2+XR2)*XAM2
      BX=B0*(ZS*S0-0.5d0*(ZP*S0P+ZM*S0M))
     : +B1*(ZS*S1-0.5d0*(ZP*S1P+ZM*S1M))
     : +B2*(ZS*S2-0.5d0*(ZP*S2P+ZM*S2M))
      BZ=B0*ALN+B1*(G1-0.5d0*(G1P+G1M))+B2*(G2-0.5d0*(G2P+G2M))
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  2x3x3=18 "CARTESIAN"
C    HARMONICS
C
         SUBROUTINE  SHLCAR3X3(A,X,Y,Z,SPS,HX,HY,HZ)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
c    harmonics (A(1)-A(36).
c  The 12 nonlinear parameters (A(37)-A(48) are the scales Pi,Ri,Qi,and Si
C   entering the arguments of exponents, sines, and cosines in each of the
C   18 "Cartesian" harmonics
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
	 IMPLICIT  REAL * 8  (A - H, O - Z)
C
         DIMENSION A(48)
C
          CPS=DSQRT(1.D0-SPS**2)
          S3PS=4.D0*CPS**2-1.D0   
c!  THIS IS SIN(3*PS)/SIN(PS)
C
           HX=0.D0
           HY=0.D0
           HZ=0.D0
           L=0
C
           DO 1 M=1,2     
c!    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(36+I)
                  Q=A(42+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
              DO 3 K=1,3
                   R=A(39+K)
                   S=A(45+K)
                   SZRK=DSIN(Z/R)
                   CZSK=DCOS(Z/S)
                   CZRK=DCOS(Z/R)
                   SZSK=DSIN(Z/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X*SQPR)
                        EQS=DEXP(X*SQQS)
C
                   DO 4 N=1,2  
c! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                  AND N=2 IS FOR THE SECOND ONE
C
                    L=L+1
                     IF (M.EQ.1) THEN
                       IF (N.EQ.1) THEN
                         DX=-SQPR*EPR*CYPI*SZRK
                         DY=EPR/P*SYPI*SZRK
                         DZ=-EPR/R*CYPI*CZRK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*CPS
                         DY=DY*CPS
                         DZ=DZ*CPS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ENDIF
                     ELSE
                       IF (N.EQ.1) THEN
                         DX=-SPS*SQQS*EQS*CYQI*CZSK
                         DY=SPS*EQS/Q*SYQI*CZSK
                         DZ=SPS*EQS/S*CYQI*SZSK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*S3PS
                         DY=DY*S3PS
                         DZ=DZ*S3PS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                       ENDIF
                 ENDIF
c
  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE
C
         RETURN
	 END

C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
       SUBROUTINE BIRK1TOT_02(PS,X,Y,Z,BX,BY,BZ)
C
C  THIS IS THE SECOND VERSION OF THE ANALYTICAL MODEL OF THE REGION 1 FIELD
C   BASED ON A SEPARATE REPRESENTATION OF THE POTENTIAL FIELD IN THE INNER AND
C   OUTER SPACE, MAPPED BY MEANS OF A SPHERO-DIPOLAR COORDINATE SYSTEM (NB #3,
C   P.91).   THE DIFFERENCE FROM THE FIRST ONE IS THAT INSTEAD OF OCTAGONAL
C   CURRENT LOOPS, CIRCULAR ONES ARE USED IN THIS VERSION FOR APPROXIMATING THE
C   FIELD IN THE OUTER REGION, WHICH IS FASTER.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION D1(3,26),D2(3,79),XI(4),C1(26),C2(79)

         COMMON /COORD11/ XX1(12),YY1(12)
         COMMON /RHDR/ RH,DR
         COMMON /LOOPDIP1/ TILT,XCENTRE(2),RADIUS(2), DIPX,DIPY
C
         COMMON /COORD21/ XX2(14),YY2(14),ZZ2(14)
         COMMON /DX1/ DX,SCALEIN,SCALEOUT
C
      DATA C1/-0.911582d-03,-0.376654d-02,-0.727423d-02,-0.270084d-02,
     * -0.123899d-02,-0.154387d-02,-0.340040d-02,-0.191858d-01,
     * -0.518979d-01,0.635061d-01,0.440680d0,-0.396570d0,0.561238d-02,
     *  0.160938d-02,-0.451229d-02,-0.251810d-02,-0.151599d-02,
     * -0.133665d-02,-0.962089d-03,-0.272085d-01,-0.524319d-01,
     *  0.717024d-01,0.523439d0,-0.405015d0,-89.5587d0,23.2806d0/

C
      DATA C2/6.04133d0,.305415d0,.606066d-02,.128379d-03,-.179406d-04,
     * 1.41714d0,-27.2586d0,-4.28833d0,-1.30675d0,35.5607d0,8.95792d0,
     : .961617d-03,
     * -.801477d-03,-.782795d-03,-1.65242d0,-16.5242d0,-5.33798d0,
     : .424878d-03,
     * .331787d-03,-.704305d-03,.844342d-03,.953682d-04,.886271d-03,
     * 25.1120d0,20.9299d0,5.14569d0,-44.1670d0,-51.0672d0,-1.87725d0,
     : 20.2998d0,
     * 48.7505d0,-2.97415d0,3.35184d0,-54.2921d0,-.838712d0,-10.5123d0,
     : 70.7594d0,
     * -4.94104d0,.106166d-03,.465791d-03,-.193719d-03,10.8439d0,
     : -29.7968d0,
     *  8.08068d0,.463507d-03,-.224475d-04,.177035d-03,-.317581d-03,
     * -.264487d-03,.102075d-03,7.71390d0,10.1915d0,-4.99797d0,
     : -23.1114d0,
     * -29.2043d0,12.2928d0,10.9542d0,33.6671d0,-9.3851d0,.174615d-03,
     : -.789777d-06,
     * .686047d-03,.460104d-04,-.345216d-02,.221871d-02,.110078d-01,
     * -.661373d-02,.249201d-02,.343978d-01,-.193145d-05,.493963d-05,
     * -.535748d-04,.191833d-04,-.100496d-03,-.210103d-03,-.232195d-02,
     * .315335d-02,-.134320d-01,-.263222d-01/
c
      DATA TILT,XCENTRE,RADIUS,DIPX,DIPY /1.00891d0,2.28397d0,
     : -5.60831d0,
     * 1.86106d0,7.83281d0,1.12541d0,0.945719d0/

      DATA DX,SCALEIN,SCALEOUT /-0.16D0,0.08D0,0.4D0/
      DATA XX1/-11.D0,2*-7.D0,2*-3.D0,3*1.D0,2*5.D0,2*9.D0/
      DATA YY1/2.D0,0.D0,4.D0,2.D0,6.D0,0.D0,4.D0,8.D0,2.D0,6.D0,0.D0,
     *  4.D0/
      DATA XX2/-10.D0,-7.D0,2*-4.D0,0.D0,2*4.D0,7.D0,10.D0,5*0.D0/
      DATA YY2/3.D0,6.D0,3.D0,9.D0,6.D0,3.D0,9.D0,6.D0,3.D0,5*0.D0/
      DATA ZZ2/2*20.D0,4.D0,20.D0,2*4.D0,3*20.D0,2.D0,3.D0,4.5D0,
     *  7.D0,10.D0/
C
      DATA RH,DR /9.D0,4.D0/   
c!  RH IS THE "HINGING DISTANCE" AND DR IS THE
C                                TRANSITION SCALE LENGTH, DEFINING THE
C                                CURVATURE  OF THE WARPING (SEE P.89, NB #2)
C
      DATA XLTDAY,XLTNGHT /78.D0,70.D0/  
c!  THESE ARE LATITUDES OF THE R-1 OVAL
C                                             AT NOON AND AT MIDNIGHT
      DATA DTET0 /0.034906d0/   
c!   THIS IS THE LATITUDINAL HALF-THICKNESS OF THE
C                                  R-1 OVAL (THE INTERPOLATION REGION BETWEEN
C                                    THE HIGH-LAT. AND THE PLASMA SHEET)
C
        TNOONN=(90.D0-XLTDAY)*0.01745329D0
        TNOONS=3.141592654D0-TNOONN     
c! HERE WE ASSUME THAT THE POSITIONS OF
C                                          THE NORTHERN AND SOUTHERN R-1 OVALS
C                                          ARE SYMMETRIC IN THE SM-COORDINATES
        DTETDN=(XLTDAY-XLTNGHT)*0.01745329D0
        DR2=DR**2
C
      SPS=DSIN(PS)
      R2=X**2+Y**2+Z**2
      R=DSQRT(R2)
      R3=R*R2
C
      RMRH=R-RH
      RPRH=R+RH
      SQM=DSQRT(RMRH**2+DR2)
      SQP=DSQRT(RPRH**2+DR2)
      C=SQP-SQM
      Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
      SPSAS=SPS/R*C/Q
      CPSAS=DSQRT(1.D0-SPSAS**2)
       XAS = X*CPSAS-Z*SPSAS
       ZAS = X*SPSAS+Z*CPSAS
        IF (XAS.NE.0.D0.OR.Y.NE.0.D0) THEN
          PAS = DATAN2(Y,XAS)
                                      ELSE
          PAS=0.D0
        ENDIF
C
      TAS=DATAN2(DSQRT(XAS**2+Y**2),ZAS)
      STAS=DSIN(TAS)
      F=STAS/(STAS**6*(1.D0-R3)+R3)**0.1666666667D0
C
      TET0=DASIN(F)
      IF (TAS.GT.1.5707963D0) TET0=3.141592654D0-TET0
      DTET=DTETDN*DSIN(PAS*0.5D0)**2
      TETR1N=TNOONN+DTET
      TETR1S=TNOONS-DTET
C
C NOW LET'S DEFINE WHICH OF THE FOUR REGIONS (HIGH-LAT., NORTHERN PSBL,
C   PLASMA SHEET, SOUTHERN PSBL) DOES THE POINT (X,Y,Z) BELONG TO:
C
       IF (TET0.LT.TETR1N-DTET0.OR.TET0.GT.TETR1S+DTET0)  LOC=1 
c! HIGH-LAT.
       IF (TET0.GT.TETR1N+DTET0.AND.TET0.LT.TETR1S-DTET0) LOC=2 
c! PL.SHEET
       IF (TET0.GE.TETR1N-DTET0.AND.TET0.LE.TETR1N+DTET0) LOC=3 
c! NORTH PSBL
       IF (TET0.GE.TETR1S-DTET0.AND.TET0.LE.TETR1S+DTET0) LOC=4 
c! SOUTH PSBL
C
       IF (LOC.EQ.1) THEN   
c! IN THE HIGH-LAT. REGION USE THE SUBROUTINE DIPOCT
C
C      print *, '  LOC=1 (HIGH-LAT)'    
c!  (test printout; disabled now)
         XI(1)=X
         XI(2)=Y
         XI(3)=Z
         XI(4)=PS
         CALL  DIPLOOP1(XI,D1)
          BX=0.D0
          BY=0.D0
          BZ=0.D0
            DO 1 I=1,26
              BX=BX+C1(I)*D1(1,I)
              BY=BY+C1(I)*D1(2,I)
  1           BZ=BZ+C1(I)*D1(3,I)
       ENDIF                                           
c!  END OF THE CASE 1
C
       IF (LOC.EQ.2) THEN
C           print *, '  LOC=2 (PLASMA SHEET)'  
c!  (test printout; disabled now)
C
         XI(1)=X
         XI(2)=Y
         XI(3)=Z
         XI(4)=PS
         CALL  CONDIP1(XI,D2)
          BX=0.D0
          BY=0.D0
          BZ=0.D0
            DO 2 I=1,79
              BX=BX+C2(I)*D2(1,I)
              BY=BY+C2(I)*D2(2,I)
  2           BZ=BZ+C2(I)*D2(3,I)
       ENDIF                                           
c!   END OF THE CASE 2
C
       IF (LOC.EQ.3) THEN
C       print *, '  LOC=3 (north PSBL)'  
c!  (test printout; disabled now)
C
         T01=TETR1N-DTET0
         T02=TETR1N+DTET0
          SQR=DSQRT(R)
          ST01AS=SQR/(R3+1.D0/DSIN(T01)**6-1.D0)**0.1666666667
          ST02AS=SQR/(R3+1.D0/DSIN(T02)**6-1.D0)**0.1666666667
          CT01AS=DSQRT(1.D0-ST01AS**2)
          CT02AS=DSQRT(1.D0-ST02AS**2)
         XAS1=R*ST01AS*DCOS(PAS)
         Y1=  R*ST01AS*DSIN(PAS)
         ZAS1=R*CT01AS
         X1=XAS1*CPSAS+ZAS1*SPSAS
         Z1=-XAS1*SPSAS+ZAS1*CPSAS 
c! X1,Y1,Z1 ARE COORDS OF THE NORTHERN
c                                                      BOUNDARY POINT
         XI(1)=X1
         XI(2)=Y1
         XI(3)=Z1
         XI(4)=PS
         CALL  DIPLOOP1(XI,D1)
          BX1=0.D0
          BY1=0.D0
          BZ1=0.D0
            DO 11 I=1,26
              BX1=BX1+C1(I)*D1(1,I) 
c!   BX1,BY1,BZ1  ARE FIELD COMPONENTS
              BY1=BY1+C1(I)*D1(2,I)  
c!  IN THE NORTHERN BOUNDARY POINT
 11           BZ1=BZ1+C1(I)*D1(3,I)  
c!
C
         XAS2=R*ST02AS*DCOS(PAS)
         Y2=  R*ST02AS*DSIN(PAS)
         ZAS2=R*CT02AS
         X2=XAS2*CPSAS+ZAS2*SPSAS
         Z2=-XAS2*SPSAS+ZAS2*CPSAS 
c! X2,Y2,Z2 ARE COORDS OF THE SOUTHERN
C                                        BOUNDARY POINT
         XI(1)=X2
         XI(2)=Y2
         XI(3)=Z2
         XI(4)=PS
         CALL  CONDIP1(XI,D2)
          BX2=0.D0
          BY2=0.D0
          BZ2=0.D0
            DO 12 I=1,79
              BX2=BX2+C2(I)*D2(1,I)
c!  BX2,BY2,BZ2  ARE FIELD COMPONENTS
              BY2=BY2+C2(I)*D2(2,I) 
c!  IN THE SOUTHERN BOUNDARY POINT
  12          BZ2=BZ2+C2(I)*D2(3,I)
C
C  NOW INTERPOLATE:
C
         SS=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
         DS=DSQRT((X-X1)**2+(Y-Y1)**2+(Z-Z1)**2)
         FRAC=DS/SS
         BX=BX1*(1.D0-FRAC)+BX2*FRAC
         BY=BY1*(1.D0-FRAC)+BY2*FRAC
         BZ=BZ1*(1.D0-FRAC)+BZ2*FRAC
C
        ENDIF                                              
c! END OF THE CASE 3
C
        IF (LOC.EQ.4) THEN
C       print *, '  LOC=4 (south PSBL)'  
c!  (test printout; disabled now)
C
         T01=TETR1S-DTET0
         T02=TETR1S+DTET0
          SQR=DSQRT(R)
          ST01AS=SQR/(R3+1.D0/DSIN(T01)**6-1.D0)**0.1666666667
          ST02AS=SQR/(R3+1.D0/DSIN(T02)**6-1.D0)**0.1666666667
          CT01AS=-DSQRT(1.D0-ST01AS**2)
          CT02AS=-DSQRT(1.D0-ST02AS**2)
         XAS1=R*ST01AS*DCOS(PAS)
         Y1=  R*ST01AS*DSIN(PAS)
         ZAS1=R*CT01AS
         X1=XAS1*CPSAS+ZAS1*SPSAS
         Z1=-XAS1*SPSAS+ZAS1*CPSAS 
c! X1,Y1,Z1 ARE COORDS OF THE NORTHERN
C                                               BOUNDARY POINT
         XI(1)=X1
         XI(2)=Y1
         XI(3)=Z1
         XI(4)=PS
         CALL  CONDIP1(XI,D2)
          BX1=0.D0
          BY1=0.D0
          BZ1=0.D0
            DO 21 I=1,79
              BX1=BX1+C2(I)*D2(1,I) 
c!  BX1,BY1,BZ1  ARE FIELD COMPONENTS
              BY1=BY1+C2(I)*D2(2,I)  
c!  IN THE NORTHERN BOUNDARY POINT
 21           BZ1=BZ1+C2(I)*D2(3,I)  
c!
C
         XAS2=R*ST02AS*DCOS(PAS)
         Y2=  R*ST02AS*DSIN(PAS)
         ZAS2=R*CT02AS
         X2=XAS2*CPSAS+ZAS2*SPSAS
         Z2=-XAS2*SPSAS+ZAS2*CPSAS 
c! X2,Y2,Z2 ARE COORDS OF THE SOUTHERN
C                                          BOUNDARY POINT
         XI(1)=X2
         XI(2)=Y2
         XI(3)=Z2
         XI(4)=PS
         CALL  DIPLOOP1(XI,D1)
          BX2=0.D0
          BY2=0.D0
          BZ2=0.D0
            DO 22 I=1,26
              BX2=BX2+C1(I)*D1(1,I) 
c!  BX2,BY2,BZ2  ARE FIELD COMPONENTS
              BY2=BY2+C1(I)*D1(2,I) 
c!     IN THE SOUTHERN BOUNDARY POINT
  22          BZ2=BZ2+C1(I)*D1(3,I)
C
C  NOW INTERPOLATE:
C
         SS=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
         DS=DSQRT((X-X1)**2+(Y-Y1)**2+(Z-Z1)**2)
         FRAC=DS/SS
         BX=BX1*(1.D0-FRAC)+BX2*FRAC
         BY=BY1*(1.D0-FRAC)+BY2*FRAC
         BZ=BZ1*(1.D0-FRAC)+BZ2*FRAC
C
        ENDIF                                        
c! END OF THE CASE 4
C
C   NOW, LET US ADD THE SHIELDING FIELD
C
        CALL  BIRK1SHLD(PS,X,Y,Z,BSX,BSY,BSZ)
        BX=BX+BSX
        BY=BY+BSY
        BZ=BZ+BSZ
         RETURN
         END
C
C------------------------------------------------------------------------------
C
C
         SUBROUTINE  DIPLOOP1(XI,D)
C
C
C      Calculates dependent model variables and their deriva-
C  tives for given independent variables and model parame-
C  ters.  Specifies model functions with free parameters which
C  must be determined by means of least squares fits (RMS
C  minimization procedure).
C
C      Description of parameters:
C
C  XI  - input vector containing independent variables;
C  D   - output double precision vector containing
C        calculated values for derivatives of dependent
C        variables with respect to LINEAR model parameters;
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
c  The  26 coefficients are moments (Z- and X-components) of 12 dipoles placed
C    inside the  R1-shell,  PLUS amplitudes of two octagonal double loops.
C     The dipoles with nonzero  Yi appear in pairs with equal moments.
c                  (see the notebook #2, pp.102-103, for details)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
         IMPLICIT  REAL * 8  (A - H, O - Z)
C
         COMMON /COORD11/ XX(12),YY(12)
         COMMON /LOOPDIP1/ TILT,XCENTRE(2),RADIUS(2),  DIPX,DIPY
         COMMON /RHDR/RH,DR
         DIMENSION XI(4),D(3,26)
C
           X = XI(1)
	   Y = XI(2)
	   Z = XI(3)
           PS= XI(4)
           SPS=DSIN(PS)
C
         DO 1 I=1,12
           R2=(XX(I)*DIPX)**2+(YY(I)*DIPY)**2
           R=DSQRT(R2)
             RMRH=R-RH
             RPRH=R+RH
             DR2=DR**2
             SQM=DSQRT(RMRH**2+DR2)
             SQP=DSQRT(RPRH**2+DR2)
             C=SQP-SQM
             Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
             SPSAS=SPS/R*C/Q
             CPSAS=DSQRT(1.D0-SPSAS**2)
         XD= (XX(I)*DIPX)*CPSAS
         YD= (YY(I)*DIPY)
         ZD=-(XX(I)*DIPX)*SPSAS
      CALL DIPXYZ(X-XD,Y-YD,Z-ZD,BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,
     *  BX1Z,BY1Z,BZ1Z)
        IF (DABS(YD).GT.1.D-10) THEN
      CALL DIPXYZ(X-XD,Y+YD,Z-ZD,BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,
     *  BX2Z,BY2Z,BZ2Z)
                                   ELSE
        BX2X=0.D0
        BY2X=0.D0
        BZ2X=0.D0
C
        BX2Z=0.D0
        BY2Z=0.D0
        BZ2Z=0.D0
                                   ENDIF
C
            D(1,I)=BX1Z+BX2Z
            D(2,I)=BY1Z+BY2Z
            D(3,I)=BZ1Z+BZ2Z
            D(1,I+12)=(BX1X+BX2X)*SPS
            D(2,I+12)=(BY1X+BY2X)*SPS
            D(3,I+12)=(BZ1X+BZ2X)*SPS
  1   CONTINUE
c
           R2=(XCENTRE(1)+RADIUS(1))**2
           R=DSQRT(R2)
             RMRH=R-RH
             RPRH=R+RH
             DR2=DR**2
             SQM=DSQRT(RMRH**2+DR2)
             SQP=DSQRT(RPRH**2+DR2)
             C=SQP-SQM
             Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
             SPSAS=SPS/R*C/Q
             CPSAS=DSQRT(1.D0-SPSAS**2)
         XOCT1= X*CPSAS-Z*SPSAS
         YOCT1= Y
         ZOCT1= X*SPSAS+Z*CPSAS
C
      CALL CROSSLP(XOCT1,YOCT1,ZOCT1,BXOCT1,BYOCT1,BZOCT1,XCENTRE(1),
     *        RADIUS(1),TILT)
            D(1,25)=BXOCT1*CPSAS+BZOCT1*SPSAS
            D(2,25)=BYOCT1
            D(3,25)=-BXOCT1*SPSAS+BZOCT1*CPSAS
C
           R2=(RADIUS(2)-XCENTRE(2))**2
           R=DSQRT(R2)
             RMRH=R-RH
             RPRH=R+RH
             DR2=DR**2
             SQM=DSQRT(RMRH**2+DR2)
             SQP=DSQRT(RPRH**2+DR2)
             C=SQP-SQM
             Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
             SPSAS=SPS/R*C/Q
             CPSAS=DSQRT(1.D0-SPSAS**2)
         XOCT2= X*CPSAS-Z*SPSAS -XCENTRE(2)
         YOCT2= Y
         ZOCT2= X*SPSAS+Z*CPSAS
            CALL CIRCLE(XOCT2,YOCT2,ZOCT2,RADIUS(2),BX,BY,BZ)
            D(1,26) =  BX*CPSAS+BZ*SPSAS
            D(2,26) =  BY
            D(3,26) = -BX*SPSAS+BZ*CPSAS
C
            RETURN
            END
c-------------------------------------------------------------------------
C
        SUBROUTINE CIRCLE(X,Y,Z,RL,BX,BY,BZ)
C
C  RETURNS COMPONENTS OF THE FIELD FROM A CIRCULAR CURRENT LOOP OF RADIUS RL
C  USES THE SECOND (MORE ACCURATE) APPROXIMATION GIVEN IN ABRAMOWITZ AND STEGUN

        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 K
cccccccc        DATA PI/3.141592654D0/
        pi=atan(1.0d0)*4
C
        RHO2=X*X+Y*Y
        RHO=DSQRT(RHO2)
        R22=Z*Z+(RHO+RL)**2
        R2=DSQRT(R22)
        R12=R22-4.D0*RHO*RL
        R32=0.5D0*(R12+R22)
        XK2=1.D0-R12/R22
        XK2S=1.D0-XK2
        DL=DLOG(1.D0/XK2S)
        K=1.38629436112d0+XK2S*(0.09666344259D0+XK2S*(0.03590092383d0+
     *     XK2S*(0.03742563713d0+XK2S*0.01451196212d0))) +DL*
     *     (0.5D0+XK2S*(0.12498593597D0+XK2S*(0.06880248576D0+
     *     XK2S*(0.03328355346D0+XK2S*0.00441787012D0))))
        E=1.D0+XK2S*(0.44325141463D0+XK2S*(0.0626060122D0+XK2S*
     *      (0.04757383546D0+XK2S*0.01736506451D0))) +DL*
     *     XK2S*(0.2499836831D0+XK2S*(0.09200180037D0+XK2S*
     *       (0.04069697526D0+XK2S*0.00526449639D0)))

        IF (RHO.GT.1.D-6) THEN
           BRHO=Z/(RHO2*R2)*(R32/R12*E-K) 
c!  THIS IS NOT EXACTLY THE B-RHO COM-
                           ELSE           
c!   PONENT - NOTE THE ADDITIONAL
           BRHO=PI*RL/R2*(RL-RHO)/R12*Z/(R32-RHO2)  
c!      DIVISION BY RHO
        ENDIF

        BX=BRHO*X
        BY=BRHO*Y
        BZ=(K-E*(R32-2.D0*RL*RL)/R12)/R2
        RETURN
        END
C-------------------------------------------------------------
C
        SUBROUTINE CROSSLP(X,Y,Z,BX,BY,BZ,XC,RL,AL)
C
c   RETURNS FIELD COMPONENTS OF A PAIR OF LOOPS WITH A COMMON CENTER AND
C    DIAMETER,  COINCIDING WITH THE X AXIS. THE LOOPS ARE INCLINED TO THE
C    EQUATORIAL PLANE BY THE ANGLE AL (RADIANS) AND SHIFTED IN THE POSITIVE
C     X-DIRECTION BY THE DISTANCE  XC.
c
        IMPLICIT REAL*8 (A-H,O-Z)
C
            CAL=DCOS(AL)
            SAL=DSIN(AL)
C
        Y1=Y*CAL-Z*SAL
        Z1=Y*SAL+Z*CAL
        Y2=Y*CAL+Z*SAL
        Z2=-Y*SAL+Z*CAL
        CALL CIRCLE(X-XC,Y1,Z1,RL,BX1,BY1,BZ1)
        CALL CIRCLE(X-XC,Y2,Z2,RL,BX2,BY2,BZ2)
        BX=BX1+BX2
        BY= (BY1+BY2)*CAL+(BZ1-BZ2)*SAL
        BZ=-(BY1-BY2)*SAL+(BZ1+BZ2)*CAL
C
        RETURN
        END

C*******************************************************************

       SUBROUTINE DIPXYZ(X,Y,Z,BXX,BYX,BZX,BXY,BYY,BZY,BXZ,BYZ,BZZ)
C
C       RETURNS THE FIELD COMPONENTS PRODUCED BY THREE DIPOLES, EACH
C        HAVING M=Me AND ORIENTED PARALLEL TO X,Y, and Z AXIS, RESP.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      X2=X**2
      Y2=Y**2
      Z2=Z**2
      R2=X2+Y2+Z2

      XMR5=30574.D0/(R2*R2*DSQRT(R2))
      XMR53=3.D0*XMR5
      BXX=XMR5*(3.D0*X2-R2)
      BYX=XMR53*X*Y
      BZX=XMR53*X*Z
C
      BXY=BYX
      BYY=XMR5*(3.D0*Y2-R2)
      BZY=XMR53*Y*Z
C
      BXZ=BZX
      BYZ=BZY
      BZZ=XMR5*(3.D0*Z2-R2)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
         SUBROUTINE  CONDIP1(XI,D)
C
C      Calculates dependent model variables and their derivatives for given
C  independent variables and model parameters.  Specifies model functions with
C  free parameters which must be determined by means of least squares fits
C  (RMS minimization procedure).
C
C      Description of parameters:
C
C  XI  - input vector containing independent variables;
C  D   - output double precision vector containing
C        calculated values for derivatives of dependent
C        variables with respect to LINEAR model parameters;
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
c  The  79 coefficients are (1) 5 amplitudes of the conical harmonics, plus
c                           (2) (9x3+5x2)x2=74 components of the dipole moments
c              (see the notebook #2, pp.113-..., for details)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
         IMPLICIT  REAL * 8  (A - H, O - Z)
C
      COMMON /DX1/ DX,SCALEIN,SCALEOUT
      COMMON /COORD21/ XX(14),YY(14),ZZ(14)
c
         DIMENSION XI(4),D(3,79),CF(5),SF(5)
C
           X = XI(1)
	   Y = XI(2)
	   Z = XI(3)
           PS= XI(4)
           SPS=DSIN(PS)
           CPS=DCOS(PS)
C
      XSM=X*CPS-Z*SPS  - DX
      ZSM=Z*CPS+X*SPS
      RO2=XSM**2+Y**2
      RO=SQRT(RO2)
C
      CF(1)=XSM/RO
      SF(1)=Y/RO
C
      CF(2)=CF(1)**2-SF(1)**2
      SF(2)=2.d0*SF(1)*CF(1)
      CF(3)=CF(2)*CF(1)-SF(2)*SF(1)
      SF(3)=SF(2)*CF(1)+CF(2)*SF(1)
      CF(4)=CF(3)*CF(1)-SF(3)*SF(1)
      SF(4)=SF(3)*CF(1)+CF(3)*SF(1)
      CF(5)=CF(4)*CF(1)-SF(4)*SF(1)
      SF(5)=SF(4)*CF(1)+CF(4)*SF(1)
C
      R2=RO2+ZSM**2
      R=DSQRT(R2)
      C=ZSM/R
      S=RO/R
      CH=DSQRT(0.5D0*(1.D0+C))
      SH=DSQRT(0.5D0*(1.D0-C))
      TNH=SH/CH
      CNH=1.D0/TNH
C
      DO 1 M=1,5
       BT=M*CF(M)/(R*S)*(TNH**M+CNH**M)
       BF=-0.5D0*M*SF(M)/R*(TNH**(M-1)/CH**2-CNH**(M-1)/SH**2)
       BXSM=BT*C*CF(1)-BF*SF(1)
       BY=BT*C*SF(1)+BF*CF(1)
       BZSM=-BT*S
C
       D(1,M)=BXSM*CPS+BZSM*SPS
       D(2,M)=BY
  1    D(3,M)=-BXSM*SPS+BZSM*CPS
C
      XSM = X*CPS-Z*SPS
      ZSM = Z*CPS+X*SPS
C
        DO 2 I=1,9
C
        IF (I.EQ.3.OR.I.EQ.5.OR.I.EQ.6) THEN
                XD =  XX(I)*SCALEIN
                YD =  YY(I)*SCALEIN
                                         ELSE
                XD =  XX(I)*SCALEOUT
                YD =  YY(I)*SCALEOUT
        ENDIF
C
         ZD =  ZZ(I)
C
      CALL DIPXYZ(XSM-XD,Y-YD,ZSM-ZD,BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,
     *  BX1Z,BY1Z,BZ1Z)
      CALL DIPXYZ(XSM-XD,Y+YD,ZSM-ZD,BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,
     *  BX2Z,BY2Z,BZ2Z)
      CALL DIPXYZ(XSM-XD,Y-YD,ZSM+ZD,BX3X,BY3X,BZ3X,BX3Y,BY3Y,BZ3Y,
     *  BX3Z,BY3Z,BZ3Z)
      CALL DIPXYZ(XSM-XD,Y+YD,ZSM+ZD,BX4X,BY4X,BZ4X,BX4Y,BY4Y,BZ4Y,
     *  BX4Z,BY4Z,BZ4Z)
C
      IX=I*3+3
      IY=IX+1
      IZ=IY+1
C
      D(1,IX)=(BX1X+BX2X-BX3X-BX4X)*CPS+(BZ1X+BZ2X-BZ3X-BZ4X)*SPS
      D(2,IX)= BY1X+BY2X-BY3X-BY4X
      D(3,IX)=(BZ1X+BZ2X-BZ3X-BZ4X)*CPS-(BX1X+BX2X-BX3X-BX4X)*SPS
C
      D(1,IY)=(BX1Y-BX2Y-BX3Y+BX4Y)*CPS+(BZ1Y-BZ2Y-BZ3Y+BZ4Y)*SPS
      D(2,IY)= BY1Y-BY2Y-BY3Y+BY4Y
      D(3,IY)=(BZ1Y-BZ2Y-BZ3Y+BZ4Y)*CPS-(BX1Y-BX2Y-BX3Y+BX4Y)*SPS
C
      D(1,IZ)=(BX1Z+BX2Z+BX3Z+BX4Z)*CPS+(BZ1Z+BZ2Z+BZ3Z+BZ4Z)*SPS
      D(2,IZ)= BY1Z+BY2Z+BY3Z+BY4Z
      D(3,IZ)=(BZ1Z+BZ2Z+BZ3Z+BZ4Z)*CPS-(BX1Z+BX2Z+BX3Z+BX4Z)*SPS
C
      IX=IX+27
      IY=IY+27
      IZ=IZ+27
C
      D(1,IX)=SPS*((BX1X+BX2X+BX3X+BX4X)*CPS+(BZ1X+BZ2X+BZ3X+BZ4X)*SPS)
      D(2,IX)=SPS*(BY1X+BY2X+BY3X+BY4X)
      D(3,IX)=SPS*((BZ1X+BZ2X+BZ3X+BZ4X)*CPS-(BX1X+BX2X+BX3X+BX4X)*SPS)
C
      D(1,IY)=SPS*((BX1Y-BX2Y+BX3Y-BX4Y)*CPS+(BZ1Y-BZ2Y+BZ3Y-BZ4Y)*SPS)
      D(2,IY)=SPS*(BY1Y-BY2Y+BY3Y-BY4Y)
      D(3,IY)=SPS*((BZ1Y-BZ2Y+BZ3Y-BZ4Y)*CPS-(BX1Y-BX2Y+BX3Y-BX4Y)*SPS)
C
      D(1,IZ)=SPS*((BX1Z+BX2Z-BX3Z-BX4Z)*CPS+(BZ1Z+BZ2Z-BZ3Z-BZ4Z)*SPS)
      D(2,IZ)=SPS*(BY1Z+BY2Z-BY3Z-BY4Z)
      D(3,IZ)=SPS*((BZ1Z+BZ2Z-BZ3Z-BZ4Z)*CPS-(BX1Z+BX2Z-BX3Z-BX4Z)*SPS)
  2   CONTINUE
C
      DO 3 I=1,5
      ZD=ZZ(I+9)
      CALL DIPXYZ(XSM,Y,ZSM-ZD,BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,BX1Z,BY1Z,
     *  BZ1Z)
      CALL DIPXYZ(XSM,Y,ZSM+ZD,BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,BX2Z,BY2Z,
     *  BZ2Z)
       IX=58+I*2
       IZ=IX+1
      D(1,IX)=(BX1X-BX2X)*CPS+(BZ1X-BZ2X)*SPS
      D(2,IX)=BY1X-BY2X
      D(3,IX)=(BZ1X-BZ2X)*CPS-(BX1X-BX2X)*SPS
C
      D(1,IZ)=(BX1Z+BX2Z)*CPS+(BZ1Z+BZ2Z)*SPS
      D(2,IZ)=BY1Z+BY2Z
      D(3,IZ)=(BZ1Z+BZ2Z)*CPS-(BX1Z+BX2Z)*SPS
C
      IX=IX+10
      IZ=IZ+10
      D(1,IX)=SPS*((BX1X+BX2X)*CPS+(BZ1X+BZ2X)*SPS)
      D(2,IX)=SPS*(BY1X+BY2X)
      D(3,IX)=SPS*((BZ1X+BZ2X)*CPS-(BX1X+BX2X)*SPS)
C
      D(1,IZ)=SPS*((BX1Z-BX2Z)*CPS+(BZ1Z-BZ2Z)*SPS)
      D(2,IZ)=SPS*(BY1Z-BY2Z)
  3   D(3,IZ)=SPS*((BZ1Z-BZ2Z)*CPS-(BX1Z-BX2Z)*SPS)
C
            RETURN
            END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
         SUBROUTINE  BIRK1SHLD(PS,X,Y,Z,BX,BY,BZ)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  The 64 linear parameters are amplitudes of the "box" harmonics.
c The 16 nonlinear parameters are the scales Pi, and Qk entering the arguments
C  of sines/cosines and exponents in each of  32 cartesian harmonics
c  N.A. Tsyganenko, Spring 1994, adjusted for the Birkeland field Aug.22, 1995
c    Revised  June 12, 1996.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(80)
      DIMENSION P1(4),R1(4),Q1(4),S1(4),RP(4),RR(4),RQ(4),RS(4)
C
      EQUIVALENCE (P1(1),A(65)),(R1(1),A(69)),(Q1(1),A(73)),
     * (S1(1),A(77))
C
      DATA A/1.174198045d0,-1.463820502d0,4.840161537d0,-3.674506864d0,
     * 82.18368896d0,-94.94071588d0,-4122.331796d0,4670.278676d0,
     : -21.54975037d0,26.72661293d0,-72.81365728d0,44.09887902d0,
     : 40.08073706d0,-51.23563510d0,1955.348537d0,-1940.971550d0,
     : 794.0496433d0,-982.2441344d0,1889.837171d0,-558.9779727d0,
     : -1260.543238d0,1260.063802d0,-293.5942373d0,344.7250789d0,
     * -773.7002492d0,957.0094135d0,-1824.143669d0,520.7994379d0,
     : 1192.484774d0,-1192.184565d0,89.15537624d0,-98.52042999d0,
     : -0.8168777675d-01,0.4255969908d-01,0.3155237661d0,
     : -0.3841755213d0,2.494553332d0,-0.6571440817d-01,-2.765661310d0,
     : 0.4331001908d0,0.1099181537d0,-0.6154126980d-01,-0.3258649260d0,
     : 0.6698439193d0,-5.542735524d0,0.1604203535d0,5.854456934d0,
     : -0.8323632049d0,3.732608869d0,-3.130002153d0,107.0972607d0,
     : -32.28483411d0,-115.2389298d0,54.45064360d0,-0.5826853320d0,
     * -3.582482231d0,-4.046544561d0,3.311978102d0,-104.0839563d0,
     : 30.26401293d0,97.29109008d0,-50.62370872d0,-296.3734955d0,
     : 127.7872523d0,5.303648988d0,10.40368955d0,69.65230348d0,
     : 466.5099509d0,1.645049286d0,3.825838190d0,11.66675599d0,
     : 558.9781177d0,1.826531343d0,2.066018073d0,25.40971369d0,
     * 990.2795225d0,2.319489258d0,4.555148484d0,9.691185703d0,
     : 591.8280358d0/
C
         BX=0.D0
         BY=0.D0
         BZ=0.D0
         CPS=DCOS(PS)
         SPS=DSIN(PS)
         S3PS=4.D0*CPS**2-1.D0
C
         DO 11 I=1,4
          RP(I)=1.D0/P1(I)
          RR(I)=1.D0/R1(I)
          RQ(I)=1.D0/Q1(I)
 11       RS(I)=1.D0/S1(I)
C
          L=0
C
           DO 1 M=1,2     
c!    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,4
                  CYPI=DCOS(Y*RP(I))
                  CYQI=DCOS(Y*RQ(I))
                  SYPI=DSIN(Y*RP(I))
                  SYQI=DSIN(Y*RQ(I))
C
                DO 3 K=1,4
                   SZRK=DSIN(Z*RR(K))
                   CZSK=DCOS(Z*RS(K))
                   CZRK=DCOS(Z*RR(K))
                   SZSK=DSIN(Z*RS(K))
                     SQPR=DSQRT(RP(I)**2+RR(K)**2)
                     SQQS=DSQRT(RQ(I)**2+RS(K)**2)
                        EPR=DEXP(X*SQPR)
                        EQS=DEXP(X*SQQS)
C
                    DO 4 N=1,2  
c! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                  AND N=2 IS FOR THE SECOND ONE
                     IF (M.EQ.1) THEN
                       IF (N.EQ.1) THEN
                         HX=-SQPR*EPR*CYPI*SZRK
                         HY=RP(I)*EPR*SYPI*SZRK
                         HZ=-RR(K)*EPR*CYPI*CZRK
                                   ELSE
                         HX=HX*CPS
                         HY=HY*CPS
                         HZ=HZ*CPS
                                   ENDIF
                     ELSE
                       IF (N.EQ.1) THEN
                         HX=-SPS*SQQS*EQS*CYQI*CZSK
                         HY=SPS*RQ(I)*EQS*SYQI*CZSK
                         HZ=SPS*RS(K)*EQS*CYQI*SZSK
                                   ELSE
                         HX=HX*S3PS
                         HY=HY*S3PS
                         HZ=HZ*S3PS
                       ENDIF
                 ENDIF
       L=L+1
c
       BX=BX+A(L)*HX
       BY=BY+A(L)*HY
  4    BZ=BZ+A(L)*HZ
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE
C
         RETURN
	 END
C
C##########################################################################
C
         SUBROUTINE BIRK2TOT_02(PS,X,Y,Z,BX,BY,BZ)
C
          IMPLICIT REAL*8 (A-H,O-Z)
C
          CALL BIRK2SHL(X,Y,Z,PS,WX,WY,WZ)
          CALL R2_BIRK(X,Y,Z,PS,HX,HY,HZ)
         BX=WX+HX
         BY=WY+HY
         BZ=WZ+HZ
         RETURN
         END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C THIS CODE IS FOR THE FIELD FROM  2x2x2=8 "CARTESIAN" HARMONICS
C
         SUBROUTINE  BIRK2SHL(X,Y,Z,PS,HX,HY,HZ)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C    The model parameters are provided to this module via common-block /A/.
C  The 16 linear parameters enter in pairs in the amplitudes of the
c       "cartesian" harmonics.
c    The 8 nonlinear parameters are the scales Pi,Ri,Qi,and Si entering the
c  arguments of exponents, sines, and cosines in each of the 8 "Cartesian"
c   harmonics
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
	 IMPLICIT  REAL * 8  (A - H, O - Z)
C
         DIMENSION P(2),R(2),Q(2),S(2)
         DIMENSION A(24)
C
         EQUIVALENCE(P(1),A(17)),(R(1),A(19)),(Q(1),A(21)),(S(1),A(23))
         DATA A/-111.6371348d0,124.5402702d0,110.3735178d0,
     : -122.0095905d0,111.9448247d0,-129.1957743d0,-110.7586562d0,
     : 126.5649012d0,-0.7865034384d0,-0.2483462721d0,0.8026023894d0,
     : 0.2531397188d0,10.72890902d0,0.8483902118d0,-10.96884315d0,
     : -0.8583297219d0,13.85650567d0,14.90554500d0,10.21914434d0,
     * 10.09021632d0,6.340382460d0,14.40432686d0,12.71023437d0,
     : 12.83966657d0/
C
            CPS=DCOS(PS)
            SPS=DSIN(PS)
            S3PS=4.D0*CPS**2-1.D0   
c!  THIS IS SIN(3*PS)/SIN(PS)
C
           HX=0.D0
           HY=0.D0
           HZ=0.D0
           L=0
C
           DO 1 M=1,2     
c!    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,2
                  CYPI=DCOS(Y/P(I))
                  CYQI=DCOS(Y/Q(I))
                  SYPI=DSIN(Y/P(I))
                  SYQI=DSIN(Y/Q(I))
C
               DO 3 K=1,2
                   SZRK=DSIN(Z/R(K))
                   CZSK=DCOS(Z/S(K))
                   CZRK=DCOS(Z/R(K))
                   SZSK=DSIN(Z/S(K))
                     SQPR=DSQRT(1.D0/P(I)**2+1.D0/R(K)**2)
                     SQQS=DSQRT(1.D0/Q(I)**2+1.D0/S(K)**2)
                        EPR=DEXP(X*SQPR)
                        EQS=DEXP(X*SQQS)
C
                   DO 4 N=1,2  
c! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                  AND N=2 IS FOR THE SECOND ONE
C
                    L=L+1
                     IF (M.EQ.1) THEN
                       IF (N.EQ.1) THEN
                         DX=-SQPR*EPR*CYPI*SZRK
                         DY=EPR/P(I)*SYPI*SZRK
                         DZ=-EPR/R(K)*CYPI*CZRK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*CPS
                         DY=DY*CPS
                         DZ=DZ*CPS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ENDIF
                     ELSE
                       IF (N.EQ.1) THEN
                         DX=-SPS*SQQS*EQS*CYQI*CZSK
                         DY=SPS*EQS/Q(I)*SYQI*CZSK
                         DZ=SPS*EQS/S(K)*CYQI*SZSK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*S3PS
                         DY=DY*S3PS
                         DZ=DZ*S3PS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                       ENDIF
                 ENDIF
c
  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE
C
         RETURN
	 END

c********************************************************************
C
       SUBROUTINE R2_BIRK(X,Y,Z,PS,BX,BY,BZ)
C
C  RETURNS THE MODEL FIELD FOR THE REGION 2 BIRKELAND CURRENT/PARTIAL RC
C    (WITHOUT SHIELDING FIELD)
C
       IMPLICIT REAL*8 (A-H,O-Z)
       SAVE PSI,CPS,SPS
       DATA DELARG/0.030D0/,DELARG1/0.015D0/,PSI/10.D0/
C
       IF (DABS(PSI-PS).GT.1.D-10) THEN
         PSI=PS
         CPS=DCOS(PS)
         SPS=DSIN(PS)
       ENDIF
C
       XSM=X*CPS-Z*SPS
       ZSM=Z*CPS+X*SPS
C
       XKS=XKSI(XSM,Y,ZSM)
      IF (XKS.LT.-(DELARG+DELARG1)) THEN
        CALL R2OUTER(XSM,Y,ZSM,BXSM,BY,BZSM)
         BXSM=-BXSM*0.02d0      
c!  ALL COMPONENTS ARE MULTIPLIED BY THE
	 BY=-BY*0.02d0          
c!  FACTOR -0.02, IN ORDER TO NORMALIZE THE
         BZSM=-BZSM*0.02d0      
c!  FIELD (SO THAT Bz=-1 nT at X=-5.3 RE, Y=Z=0)
      ENDIF
      IF (XKS.GE.-(DELARG+DELARG1).AND.XKS.LT.-DELARG+DELARG1) THEN
        CALL R2OUTER(XSM,Y,ZSM,BXSM1,BY1,BZSM1)
        CALL R2SHEET(XSM,Y,ZSM,BXSM2,BY2,BZSM2)
        F2=-0.02d0*TKSI(XKS,-DELARG,DELARG1)
        F1=-0.02d0-F2
        BXSM=BXSM1*F1+BXSM2*F2
        BY=BY1*F1+BY2*F2
        BZSM=BZSM1*F1+BZSM2*F2
      ENDIF

      IF (XKS.GE.-DELARG+DELARG1.AND.XKS.LT.DELARG-DELARG1) THEN
       CALL R2SHEET(XSM,Y,ZSM,BXSM,BY,BZSM)
         BXSM=-BXSM*0.02d0
         BY=-BY*0.02d0
         BZSM=-BZSM*0.02d0
      ENDIF
      IF (XKS.GE.DELARG-DELARG1.AND.XKS.LT.DELARG+DELARG1) THEN
        CALL R2INNER(XSM,Y,ZSM,BXSM1,BY1,BZSM1)
        CALL R2SHEET(XSM,Y,ZSM,BXSM2,BY2,BZSM2)
        F1=-0.02d0*TKSI(XKS,DELARG,DELARG1)
        F2=-0.02d0-F1
        BXSM=BXSM1*F1+BXSM2*F2
        BY=BY1*F1+BY2*F2
        BZSM=BZSM1*F1+BZSM2*F2
      ENDIF
      IF (XKS.GE.DELARG+DELARG1) THEN
         CALL R2INNER(XSM,Y,ZSM,BXSM,BY,BZSM)
         BXSM=-BXSM*0.02d0
         BY=-BY*0.02d0
         BZSM=-BZSM*0.02d0
      ENDIF
C
        BX=BXSM*CPS+BZSM*SPS
        BZ=BZSM*CPS-BXSM*SPS
C
        RETURN
        END
C
C****************************************************************

c
      SUBROUTINE R2INNER (X,Y,Z,BX,BY,BZ)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CBX(5),CBY(5),CBZ(5)
C
      DATA PL1,PL2,PL3,PL4,PL5,PL6,PL7,PL8/154.185d0,-2.12446d0,
     : .601735d-01,-.153954d-02,.355077d-04,29.9996d0,262.886d0,
     :  99.9132d0/
      DATA PN1,PN2,PN3,PN4,PN5,PN6,PN7,PN8/-8.1902d0,6.5239d0,5.504d0,
     : 7.7815d0,.8573d0,3.0986d0,.0774d0,-.038d0/
C
      CALL BCONIC(X,Y,Z,CBX,CBY,CBZ,5)
C
C   NOW INTRODUCE  ONE  4-LOOP SYSTEM:
C
       CALL LOOPS4(X,Y,Z,DBX8,DBY8,DBZ8,PN1,PN2,PN3,PN4,PN5,PN6)
C
       CALL DIPDISTR(X-PN7,Y,Z,DBX6,DBY6,DBZ6,0)
       CALL DIPDISTR(X-PN8,Y,Z,DBX7,DBY7,DBZ7,1)

C                           NOW COMPUTE THE FIELD COMPONENTS:

      BX=PL1*CBX(1)+PL2*CBX(2)+PL3*CBX(3)+PL4*CBX(4)+PL5*CBX(5)
     * +PL6*DBX6+PL7*DBX7+PL8*DBX8
      BY=PL1*CBY(1)+PL2*CBY(2)+PL3*CBY(3)+PL4*CBY(4)+PL5*CBY(5)
     * +PL6*DBY6+PL7*DBY7+PL8*DBY8
      BZ=PL1*CBZ(1)+PL2*CBZ(2)+PL3*CBZ(3)+PL4*CBZ(4)+PL5*CBZ(5)
     * +PL6*DBZ6+PL7*DBZ7+PL8*DBZ8
C
      RETURN
      END
C-----------------------------------------------------------------------

      SUBROUTINE BCONIC(X,Y,Z,CBX,CBY,CBZ,NMAX)
C
c   "CONICAL" HARMONICS
c
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION CBX(NMAX),CBY(NMAX),CBZ(NMAX)

       RO2=X**2+Y**2
       RO=SQRT(RO2)
C
       CF=X/RO
       SF=Y/RO
       CFM1=1.D0
       SFM1=0.D0
C
      R2=RO2+Z**2
      R=DSQRT(R2)
      C=Z/R
      S=RO/R
      CH=DSQRT(0.5D0*(1.D0+C))
      SH=DSQRT(0.5D0*(1.D0-C))
      TNHM1=1.D0
      CNHM1=1.D0
      TNH=SH/CH
      CNH=1.D0/TNH
C
      DO 1 M=1,NMAX
        CFM=CFM1*CF-SFM1*SF
        SFM=CFM1*SF+SFM1*CF
        CFM1=CFM
        SFM1=SFM
        TNHM=TNHM1*TNH
        CNHM=CNHM1*CNH
       BT=M*CFM/(R*S)*(TNHM+CNHM)
       BF=-0.5D0*M*SFM/R*(TNHM1/CH**2-CNHM1/SH**2)
         TNHM1=TNHM
         CNHM1=CNHM
       CBX(M)=BT*C*CF-BF*SF
       CBY(M)=BT*C*SF+BF*CF
  1    CBZ(M)=-BT*S
C
       RETURN
       END

C-------------------------------------------------------------------
C
       SUBROUTINE DIPDISTR(X,Y,Z,BX,BY,BZ,MODE)
C
C   RETURNS FIELD COMPONENTS FROM A LINEAR DISTRIBUTION OF DIPOLAR SOURCES
C     ON THE Z-AXIS.  THE PARAMETER MODE DEFINES HOW THE DIPOLE STRENGTH
C     VARIES ALONG THE Z-AXIS:  MODE=0 IS FOR A STEP-FUNCTION (Mx=const &gt 0
c         FOR Z &gt 0, AND Mx=-const &lt 0 FOR Z &lt 0)
C      WHILE MODE=1 IS FOR A LINEAR VARIATION OF THE DIPOLE MOMENT DENSITY
C       SEE NB#3, PAGE 53 FOR DETAILS.
C
C
C INPUT: X,Y,Z OF A POINT OF SPACE, AND MODE
C
        IMPLICIT REAL*8 (A-H,O-Z)
        X2=X*X
        RHO2=X2+Y*Y
        R2=RHO2+Z*Z
        R3=R2*DSQRT(R2)

        IF (MODE.EQ.0) THEN
         BX=Z/RHO2**2*(R2*(Y*Y-X2)-RHO2*X2)/R3
         BY=-X*Y*Z/RHO2**2*(2.D0*R2+RHO2)/R3
         BZ=X/R3
        ELSE
         BX=Z/RHO2**2*(Y*Y-X2)
         BY=-2.D0*X*Y*Z/RHO2**2
         BZ=X/RHO2
        ENDIF
         RETURN
         END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE R2OUTER (X,Y,Z,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DATA PL1,PL2,PL3,PL4,PL5/-34.105d0,-2.00019d0,628.639d0,
     : 73.4847d0,12.5162d0/
      DATA PN1,PN2,PN3,PN4,PN5,PN6,PN7,PN8,PN9,PN10,PN11,PN12,PN13,PN14,
     *  PN15,PN16,PN17 /.55d0,.694d0,.0031d0,1.55d0,2.8d0,.1375d0,
     : -.7d0,.2d0,.9625d0,-2.994d0,2.925d0,-1.775d0,4.3d0,-.275d0,
     : 2.7d0,.4312d0,1.55d0/
c
C    THREE PAIRS OF CROSSED LOOPS:
C
      CALL CROSSLP(X,Y,Z,DBX1,DBY1,DBZ1,PN1,PN2,PN3)
      CALL CROSSLP(X,Y,Z,DBX2,DBY2,DBZ2,PN4,PN5,PN6)
      CALL CROSSLP(X,Y,Z,DBX3,DBY3,DBZ3,PN7,PN8,PN9)
C
C    NOW AN EQUATORIAL LOOP ON THE NIGHTSIDE
C
      CALL CIRCLE(X-PN10,Y,Z,PN11,DBX4,DBY4,DBZ4)
c
c   NOW A 4-LOOP SYSTEM ON THE NIGHTSIDE
c

      CALL LOOPS4(X,Y,Z,DBX5,DBY5,DBZ5,PN12,PN13,PN14,PN15,PN16,PN17)

C---------------------------------------------------------------------

C                           NOW COMPUTE THE FIELD COMPONENTS:

      BX=PL1*DBX1+PL2*DBX2+PL3*DBX3+PL4*DBX4+PL5*DBX5
      BY=PL1*DBY1+PL2*DBY2+PL3*DBY3+PL4*DBY4+PL5*DBY5
      BZ=PL1*DBZ1+PL2*DBZ2+PL3*DBZ3+PL4*DBZ4+PL5*DBZ5

       RETURN
       END
C
C--------------------------------------------------------------------
C
       SUBROUTINE LOOPS4(X,Y,Z,BX,BY,BZ,XC,YC,ZC,R,THETA,PHI)
C
C   RETURNS FIELD COMPONENTS FROM A SYSTEM OF 4 CURRENT LOOPS, POSITIONED
C     SYMMETRICALLY WITH RESPECT TO NOON-MIDNIGHT MERIDIAN AND EQUATORIAL
C      PLANES.
C  INPUT: X,Y,Z OF A POINT OF SPACE
C        XC,YC,ZC (YC &gt 0 AND ZC &gt 0) - POSITION OF THE CENTER OF THE
C                                         1ST-QUADRANT LOOP
C        R - LOOP RADIUS (THE SAME FOR ALL FOUR)
C        THETA, PHI  -  SPECIFY THE ORIENTATION OF THE NORMAL OF THE 1ST LOOP
c      -----------------------------------------------------------

        IMPLICIT REAL*8 (A-H,O-Z)
C
          CT=DCOS(THETA)
          ST=DSIN(THETA)
          CP=DCOS(PHI)
          SP=DSIN(PHI)
C------------------------------------1ST QUADRANT:
        XS=(X-XC)*CP+(Y-YC)*SP
        YSS=(Y-YC)*CP-(X-XC)*SP
        ZS=Z-ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ1=BZSS*CT-BXSS*ST
          BX1=BXS*CP-BYS*SP
          BY1=BXS*SP+BYS*CP
C-------------------------------------2nd QUADRANT:
        XS=(X-XC)*CP-(Y+YC)*SP
        YSS=(Y+YC)*CP+(X-XC)*SP
        ZS=Z-ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ2=BZSS*CT-BXSS*ST
          BX2=BXS*CP+BYS*SP
          BY2=-BXS*SP+BYS*CP
C-------------------------------------3RD QUADRANT:
        XS=-(X-XC)*CP+(Y+YC)*SP
        YSS=-(Y+YC)*CP-(X-XC)*SP
        ZS=Z+ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ3=BZSS*CT-BXSS*ST
          BX3=-BXS*CP-BYS*SP
          BY3=BXS*SP-BYS*CP
C-------------------------------------4TH QUADRANT:
        XS=-(X-XC)*CP-(Y-YC)*SP
        YSS=-(Y-YC)*CP+(X-XC)*SP
        ZS=Z+ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ4=BZSS*CT-BXSS*ST
          BX4=-BXS*CP+BYS*SP
          BY4=-BXS*SP-BYS*CP

        BX=BX1+BX2+BX3+BX4
        BY=BY1+BY2+BY3+BY4
        BZ=BZ1+BZ2+BZ3+BZ4

         RETURN
         END
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
      SUBROUTINE R2SHEET(X,Y,Z,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DATA PNONX1,PNONX2,PNONX3,PNONX4,PNONX5,PNONX6,PNONX7,PNONX8,
     *     PNONY1,PNONY2,PNONY3,PNONY4,PNONY5,PNONY6,PNONY7,PNONY8,
     *     PNONZ1,PNONZ2,PNONZ3,PNONZ4,PNONZ5,PNONZ6,PNONZ7,PNONZ8
     */-19.0969D0,-9.28828D0,-0.129687D0,5.58594D0,22.5055D0,
     *  0.483750D-01,0.396953D-01,0.579023D-01,-13.6750D0,-6.70625D0,
     *  2.31875D0,11.4062D0,20.4562D0,0.478750D-01,0.363750D-01,
     * 0.567500D-01,-16.7125D0,-16.4625D0,-0.1625D0,5.1D0,23.7125D0,
     * 0.355625D-01,0.318750D-01,0.538750D-01/
C
C
      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,
     *  A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,
     *  A34,A35,A36,A37,A38,A39,A40,A41,A42,A43,A44,A45,A46,A47,A48,A49,
     *  A50,A51,A52,A53,A54,A55,A56,A57,A58,A59,A60,A61,A62,A63,A64,A65,
     *  A66,A67,A68,A69,A70,A71,A72,A73,A74,A75,A76,A77,A78,A79,A80
     * /8.07190D0,-7.39582D0,-7.62341D0,0.684671D0,-13.5672D0,11.6681D0,
     * 13.1154,-0.890217D0,7.78726D0,-5.38346D0,-8.08738D0,0.609385D0,
     * -2.70410D0, 3.53741D0,3.15549D0,-1.11069D0,-8.47555D0,0.278122D0,
     *  2.73514D0,4.55625D0,13.1134D0,1.15848D0,-3.52648D0,-8.24698D0,
     * -6.85710D0,-2.81369D0, 2.03795D0, 4.64383D0,2.49309D0,-1.22041D0,
     * -1.67432D0,-0.422526D0,-5.39796D0,7.10326D0,5.53730D0,-13.1918D0,
     *  4.67853D0,-7.60329D0,-2.53066D0, 7.76338D0, 5.60165D0,5.34816D0,
     * -4.56441D0,7.05976D0,-2.62723D0,-0.529078D0,1.42019D0,-2.93919D0,
     *  55.6338D0,-1.55181D0,39.8311D0,-80.6561D0,-46.9655D0,32.8925D0,
     * -6.32296D0,19.7841D0,124.731D0,10.4347D0,-30.7581D0,102.680D0,
     * -47.4037D0,-3.31278D0,9.37141D0,-50.0268D0,-533.319D0,110.426D0,
     *  1000.20D0,-1051.40D0, 1619.48D0,589.855D0,-1462.73D0,1087.10D0,
     *  -1994.73D0,-1654.12D0,1263.33D0,-260.210D0,1424.84D0,1255.71D0,
     *  -956.733D0, 219.946D0/
C
C
      DATA B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,
     *  B18,B19,B20,B21,B22,B23,B24,B25,B26,B27,B28,B29,B30,B31,B32,B33,
     *  B34,B35,B36,B37,B38,B39,B40,B41,B42,B43,B44,B45,B46,B47,B48,B49,
     *  B50,B51,B52,B53,B54,B55,B56,B57,B58,B59,B60,B61,B62,B63,B64,B65,
     *  B66,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B77,B78,B79,B80
     */-9.08427D0,10.6777D0,10.3288D0,-0.969987D0,6.45257D0,-8.42508D0,
     * -7.97464D0,1.41996D0,-1.92490D0,3.93575D0,2.83283D0,-1.48621D0,
     *0.244033D0,-0.757941D0,-0.386557D0,0.344566D0,9.56674D0,-2.5365D0,
     * -3.32916D0,-5.86712D0,-6.19625D0,1.83879D0,2.52772D0,4.34417D0,
     * 1.87268D0,-2.13213D0,-1.69134D0,-.176379D0,-.261359D0,.566419D0,
     * 0.3138D0,-0.134699D0,-3.83086D0,-8.4154D0,4.77005D0,-9.31479D0,
     * 37.5715D0,19.3992D0,-17.9582D0,36.4604D0,-14.9993D0,-3.1442D0,
     * 6.17409D0,-15.5519D0,2.28621D0,-0.891549D-2,-.462912D0,2.47314D0,
     * 41.7555D0,208.614D0,-45.7861D0,-77.8687D0,239.357D0,-67.9226D0,
     * 66.8743D0,238.534D0,-112.136D0,16.2069D0,-40.4706D0,-134.328D0,
     * 21.56D0,-0.201725D0,2.21D0,32.5855D0,-108.217D0,-1005.98D0,
     * 585.753D0,323.668D0,-817.056D0,235.750D0,-560.965D0,-576.892D0,
     * 684.193D0,85.0275D0,168.394D0,477.776D0,-289.253D0,-123.216D0,
     * 75.6501D0,-178.605D0/
C
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,
     *  C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,
     *  C34,C35,C36,C37,C38,C39,C40,C41,C42,C43,C44,C45,C46,C47,C48,C49,
     *  C50,C51,C52,C53,C54,C55,C56,C57,C58,C59,C60,C61,C62,C63,C64,C65,
     *  C66,C67,C68,C69,C70,C71,C72,C73,C74,C75,C76,C77,C78,C79,C80
     * / 1167.61D0,-917.782D0,-1253.2D0,-274.128D0,-1538.75D0,1257.62D0,
     * 1745.07D0,113.479D0,393.326D0,-426.858D0,-641.1D0,190.833D0,
     * -29.9435D0,-1.04881D0,117.125D0,-25.7663D0,-1168.16D0,910.247D0,
     * 1239.31D0,289.515D0,1540.56D0,-1248.29D0,-1727.61D0,-131.785D0,
     * -394.577D0,426.163D0,637.422D0,-187.965D0,30.0348D0,0.221898D0,
     * -116.68D0,26.0291D0,12.6804D0,4.84091D0,1.18166D0,-2.75946D0,
     * -17.9822D0,-6.80357D0,-1.47134D0,3.02266D0,4.79648D0,0.665255D0,
     * -0.256229D0,-0.857282D-1,-0.588997D0,0.634812D-1,0.164303D0,
     * -0.15285D0,22.2524D0,-22.4376D0,-3.85595D0,6.07625D0,-105.959D0,
     * -41.6698D0,0.378615D0,1.55958D0,44.3981D0,18.8521D0,3.19466D0,
     *  5.89142D0,-8.63227D0,-2.36418D0,-1.027D0,-2.31515D0,1035.38D0,
     *  2040.66D0,-131.881D0,-744.533D0,-3274.93D0,-4845.61D0,482.438D0,
     * 1567.43D0,1354.02D0,2040.47D0,-151.653D0,-845.012D0,-111.723D0,
     * -265.343D0,-26.1171D0,216.632D0/
C
c------------------------------------------------------------------
C
       XKS=XKSI(X,Y,Z)    
c!  variation across the current sheet
       T1X=XKS/DSQRT(XKS**2+PNONX6**2)
       T2X=PNONX7**3/DSQRT(XKS**2+PNONX7**2)**3
       T3X=XKS/DSQRT(XKS**2+PNONX8**2)**5 *3.493856D0*PNONX8**4
C
       T1Y=XKS/DSQRT(XKS**2+PNONY6**2)
       T2Y=PNONY7**3/DSQRT(XKS**2+PNONY7**2)**3
       T3Y=XKS/DSQRT(XKS**2+PNONY8**2)**5 *3.493856D0*PNONY8**4
C
       T1Z=XKS/DSQRT(XKS**2+PNONZ6**2)
       T2Z=PNONZ7**3/DSQRT(XKS**2+PNONZ7**2)**3
       T3Z=XKS/DSQRT(XKS**2+PNONZ8**2)**5 *3.493856D0*PNONZ8**4
C
      RHO2=X*X+Y*Y
      R=DSQRT(RHO2+Z*Z)
      RHO=DSQRT(RHO2)
C
      C1P=X/RHO
      S1P=Y/RHO
      S2P=2.D0*S1P*C1P
      C2P=C1P*C1P-S1P*S1P
      S3P=S2P*C1P+C2P*S1P
      C3P=C2P*C1P-S2P*S1P
      S4P=S3P*C1P+C3P*S1P
      CT=Z/R
      ST=RHO/R
C
      S1=FEXP(CT,PNONX1)
      S2=FEXP(CT,PNONX2)
      S3=FEXP(CT,PNONX3)
      S4=FEXP(CT,PNONX4)
      S5=FEXP(CT,PNONX5)
C
C                   NOW COMPUTE THE GSM FIELD COMPONENTS:
C
C
      BX=S1*((A1+A2*T1X+A3*T2X+A4*T3X)
     *        +C1P*(A5+A6*T1X+A7*T2X+A8*T3X)
     *        +C2P*(A9+A10*T1X+A11*T2X+A12*T3X)
     *        +C3P*(A13+A14*T1X+A15*T2X+A16*T3X))
     *    +S2*((A17+A18*T1X+A19*T2X+A20*T3X)
     *        +C1P*(A21+A22*T1X+A23*T2X+A24*T3X)
     *        +C2P*(A25+A26*T1X+A27*T2X+A28*T3X)
     *        +C3P*(A29+A30*T1X+A31*T2X+A32*T3X))
     *    +S3*((A33+A34*T1X+A35*T2X+A36*T3X)
     *        +C1P*(A37+A38*T1X+A39*T2X+A40*T3X)
     *        +C2P*(A41+A42*T1X+A43*T2X+A44*T3X)
     *        +C3P*(A45+A46*T1X+A47*T2X+A48*T3X))
     *    +S4*((A49+A50*T1X+A51*T2X+A52*T3X)
     *        +C1P*(A53+A54*T1X+A55*T2X+A56*T3X)
     *        +C2P*(A57+A58*T1X+A59*T2X+A60*T3X)
     *        +C3P*(A61+A62*T1X+A63*T2X+A64*T3X))
     *    +S5*((A65+A66*T1X+A67*T2X+A68*T3X)
     *        +C1P*(A69+A70*T1X+A71*T2X+A72*T3X)
     *        +C2P*(A73+A74*T1X+A75*T2X+A76*T3X)
     *        +C3P*(A77+A78*T1X+A79*T2X+A80*T3X))
C
C
      S1=FEXP(CT,PNONY1)
      S2=FEXP(CT,PNONY2)
      S3=FEXP(CT,PNONY3)
      S4=FEXP(CT,PNONY4)
      S5=FEXP(CT,PNONY5)
C
      BY=S1*(S1P*(B1+B2*T1Y+B3*T2Y+B4*T3Y)
     *      +S2P*(B5+B6*T1Y+B7*T2Y+B8*T3Y)
     *      +S3P*(B9+B10*T1Y+B11*T2Y+B12*T3Y)
     *      +S4P*(B13+B14*T1Y+B15*T2Y+B16*T3Y))
     *  +S2*(S1P*(B17+B18*T1Y+B19*T2Y+B20*T3Y)
     *      +S2P*(B21+B22*T1Y+B23*T2Y+B24*T3Y)
     *      +S3P*(B25+B26*T1Y+B27*T2Y+B28*T3Y)
     *      +S4P*(B29+B30*T1Y+B31*T2Y+B32*T3Y))
     *  +S3*(S1P*(B33+B34*T1Y+B35*T2Y+B36*T3Y)
     *      +S2P*(B37+B38*T1Y+B39*T2Y+B40*T3Y)
     *      +S3P*(B41+B42*T1Y+B43*T2Y+B44*T3Y)
     *      +S4P*(B45+B46*T1Y+B47*T2Y+B48*T3Y))
     *  +S4*(S1P*(B49+B50*T1Y+B51*T2Y+B52*T3Y)
     *      +S2P*(B53+B54*T1Y+B55*T2Y+B56*T3Y)
     *      +S3P*(B57+B58*T1Y+B59*T2Y+B60*T3Y)
     *      +S4P*(B61+B62*T1Y+B63*T2Y+B64*T3Y))
     *  +S5*(S1P*(B65+B66*T1Y+B67*T2Y+B68*T3Y)
     *      +S2P*(B69+B70*T1Y+B71*T2Y+B72*T3Y)
     *      +S3P*(B73+B74*T1Y+B75*T2Y+B76*T3Y)
     *      +S4P*(B77+B78*T1Y+B79*T2Y+B80*T3Y))
C
      S1=FEXP1(CT,PNONZ1)
      S2=FEXP1(CT,PNONZ2)
      S3=FEXP1(CT,PNONZ3)
      S4=FEXP1(CT,PNONZ4)
      S5=FEXP1(CT,PNONZ5)
C
      BZ=S1*((C1+C2*T1Z+C3*T2Z+C4*T3Z)
     *      +C1P*(C5+C6*T1Z+C7*T2Z+C8*T3Z)
     *      +C2P*(C9+C10*T1Z+C11*T2Z+C12*T3Z)
     *      +C3P*(C13+C14*T1Z+C15*T2Z+C16*T3Z))
     *   +S2*((C17+C18*T1Z+C19*T2Z+C20*T3Z)
     *      +C1P*(C21+C22*T1Z+C23*T2Z+C24*T3Z)
     *      +C2P*(C25+C26*T1Z+C27*T2Z+C28*T3Z)
     *      +C3P*(C29+C30*T1Z+C31*T2Z+C32*T3Z))
     *   +S3*((C33+C34*T1Z+C35*T2Z+C36*T3Z)
     *      +C1P*(C37+C38*T1Z+C39*T2Z+C40*T3Z)
     *      +C2P*(C41+C42*T1Z+C43*T2Z+C44*T3Z)
     *      +C3P*(C45+C46*T1Z+C47*T2Z+C48*T3Z))
     *   +S4*((C49+C50*T1Z+C51*T2Z+C52*T3Z)
     *      +C1P*(C53+C54*T1Z+C55*T2Z+C56*T3Z)
     *      +C2P*(C57+C58*T1Z+C59*T2Z+C60*T3Z)
     *      +C3P*(C61+C62*T1Z+C63*T2Z+C64*T3Z))
     *   +S5*((C65+C66*T1Z+C67*T2Z+C68*T3Z)
     *      +C1P*(C69+C70*T1Z+C71*T2Z+C72*T3Z)
     *      +C2P*(C73+C74*T1Z+C75*T2Z+C76*T3Z)
     *      +C3P*(C77+C78*T1Z+C79*T2Z+C80*T3Z))
C
       RETURN
       END
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      REAL*8 FUNCTION XKSI(X,Y,Z)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C   A11 - C72, R0, and DR below  ARE STRETCH PARAMETERS (P.26-27, NB# 3),
C
      DATA A11A12,A21A22,A41A42,A51A52,A61A62,B11B12,B21B22,C61C62,
     *  C71C72,R0,DR /0.305662d0,-0.383593d0,0.2677733d0,-0.097656d0,
     :  -0.636034d0,-0.359862d0,0.424706d0,-0.126366d0,0.292578d0,
     :  1.21563d0,7.50937d0/

      DATA TNOON,DTETA/0.3665191d0,0.09599309d0/ 
c! Correspond to noon and midnight
C                                         latitudes 69 and 63.5 degs, resp.
       DR2=DR*DR
C
       X2=X*X
       Y2=Y*Y
       Z2=Z*Z
       XY=X*Y
       XYZ=XY*Z
       R2=X2+Y2+Z2
       R=DSQRT(R2)
       R3=R2*R
       R4=R2*R2
       XR=X/R
       YR=Y/R
       ZR=Z/R
C
       IF (R.LT.R0) THEN
         PR=0.D0
       ELSE
         PR=DSQRT((R-R0)**2+DR2)-DR
       ENDIF
C
      F=X+PR*(A11A12+A21A22*XR+A41A42*XR*XR+A51A52*YR*YR+
     *        A61A62*ZR*ZR)
      G=Y+PR*(B11B12*YR+B21B22*XR*YR)
      H=Z+PR*(C61C62*ZR+C71C72*XR*ZR)
      G2=G*G
C
      FGH=F**2+G2+H**2
      FGH32=DSQRT(FGH)**3
      FCHSG2=F**2+G2

      IF (FCHSG2.LT.1.D-5) THEN
         XKSI=-1.D0               
c!  THIS IS JUST FOR ELIMINATING PROBLEMS
         RETURN                    
c!  ON THE Z-AXIS
      ENDIF

      SQFCHSG2=DSQRT(FCHSG2)
      ALPHA=FCHSG2/FGH32
      THETA=TNOON+0.5D0*DTETA*(1.D0-F/SQFCHSG2)
      PHI=DSIN(THETA)**2
C
      XKSI=ALPHA-PHI
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
        FUNCTION FEXP(S,A)
         IMPLICIT REAL*8 (A-H,O-Z)
          DATA E/2.718281828459D0/
          IF (A.LT.0.D0) FEXP=DSQRT(-2.D0*A*E)*S*DEXP(A*S*S)
          IF (A.GE.0.D0) FEXP=S*DEXP(A*(S*S-1.D0))
         RETURN
         END
C
C-----------------------------------------------------------------------
        FUNCTION FEXP1(S,A)
         IMPLICIT REAL*8 (A-H,O-Z)
         IF (A.LE.0.D0) FEXP1=DEXP(A*S*S)
         IF (A.GT.0.D0) FEXP1=DEXP(A*(S*S-1.D0))
         RETURN
         END
C
C************************************************************************
C
         REAL*8 FUNCTION TKSI(XKSI,XKS0,DXKSI)
         IMPLICIT REAL*8 (A-H,O-Z)
         SAVE M,TDZ3
         DATA M/0/
C
         IF (M.EQ.0) THEN
         TDZ3=2.d0*DXKSI**3
         M=1
         ENDIF
C
         IF (XKSI-XKS0.LT.-DXKSI) TKSII=0.d0
         IF (XKSI-XKS0.GE.DXKSI)  TKSII=1.d0
C
         IF (XKSI.GE.XKS0-DXKSI.AND.XKSI.LT.XKS0) THEN
           BR3=(XKSI-XKS0+DXKSI)**3
           TKSII=1.5d0*BR3/(TDZ3+BR3)
         ENDIF
C
         IF (XKSI.GE.XKS0.AND.XKSI.LT.XKS0+DXKSI) THEN
           BR3=(XKSI-XKS0-DXKSI)**3
           TKSII=1.d0+1.5d0*BR3/(TDZ3-BR3)
         ENDIF
           TKSI=TKSII
         END
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
       SUBROUTINE DIPOLE(PS,X,Y,Z,BX,BY,BZ)
C
C  CALCULATES GSM COMPONENTS OF GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH OF 1980.
C------------INPUT PARAMETERS:
C   PS - GEODIPOLE TILT ANGLE IN RADIANS, X,Y,Z - GSM COORDINATES IN RE
C------------OUTPUT PARAMETERS:
C   BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
C
C                   AUTHOR: NIKOLAI A. TSYGANENKO
C                           INSTITUTE OF PHYSICS
C                           ST.-PETERSBURG STATE UNIVERSITY
C                           STARY PETERGOF 198904
C                           ST.-PETERSBURG
C                           RUSSIA
C
      IMPLICIT NONE
C
      REAL*8 PS,X,Y,Z,BX,BY,BZ,PSI,SPS,CPS,P,U,V,T,Q
      INTEGER*4 M
      SAVE M,PSI,SPS,CPS

      DATA M,PSI/0,5.d0/
      IF(M.EQ.1.AND.ABS(PS-PSI).LT.1.d-5) GOTO 1
      SPS=SIN(PS)
      CPS=COS(PS)
      PSI=PS
      M=1
  1   P=X**2
      U=Z**2
      V=3.d0*Z*X
      T=Y**2
      Q=30574.d0/SQRT(P+T+U)**5
      BX=Q*((T+U-2.d0*P)*SPS-V*CPS)
      BY=-3.d0*Y*Q*(X*SPS+Z*CPS)
      BZ=Q*((P+T-2.d0*U)*CPS-V*SPS)
      RETURN
      END
c
c
c
c

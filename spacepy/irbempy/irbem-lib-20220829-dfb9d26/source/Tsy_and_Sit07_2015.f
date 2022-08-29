c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c Implemented by A. C. Kellerman, including improved calls to
c bessj from Jay Albert (in TS07j_2015.f renamed from TS07j.f as we now have a 2017 version of this code that implements the additional functions directly by G. K. Stephens)

      SUBROUTINE TS07D_2015 (XGSM,YGSM,ZGSM,BX,BY,BZ)
   
      IMPLICIT NONE

      REAL*8 XGSM,YGSM,ZGSM,BX,BY,BZ,PSS
      REAL*8 PDYN,TILT,AAA,SPS,CPS,BBB,PSI,CCC,A07
      INTEGER*4 M_INX,N_INX,NTOT
    
      PARAMETER(NTOT=101)

      REAL*8 BXCF,BYCF,BZCF
      REAL*8 BXR11,BYR11,BZR11
      REAL*8 BXR12,BYR12,BZR12
      REAL*8 BXR21a,BYR21a,BZR21a
      REAL*8 BXR21s,BYR21s,BZR21s

!      REAL*8, DIMENSION(5) :: BXTS,BYTS,BZTS
!      REAL*8, DIMENSION(5,4) :: BXTO,BYTO,BZTO,BXTE,BYTE,BZTE

      REAL*8 BXTS,BYTS,BZTS
      DIMENSION BXTS(5),BYTS(5),BZTS(5)
      REAL*8 BXTO,BYTO,BZTO,BXTE,BYTE,BZTE
      DIMENSION BXTO(5,4),BYTO(5,4),BZTO(5,4),BXTE(5,4),
     &BYTE(5,4),BZTE(5,4)


! for compatibility with older compilers:
!      REAL*8 BXTS(5),BYTS(5),BZTS(5)
!      REAL*8 BXTO(5,4),BYTO(5,4),BZTO(5,4),BXTE(5,4),BYTE(5,4),BZTE(5,4)

      COMMON /GEOPACK1/ AAA(10),SPS,CPS,BBB(3),PSI,CCC(18)
      COMMON /TS07D_DATA/ M_INX,N_INX,PDYN,TILT,A07(NTOT)

      PSS = PSI * 1.D0 

c      print *, '========================'
c      print *, 'M_INX = ',M_INX
c      print *, 'N_INX = ',N_INX
c      print *, 'PDYN = ',PDYN
c      print *, 'TILT = ', TILT
c      print *, 'PSI = ', PSI 
c      print *, '========================'
          
      CALL EXTERN_2015 (0,A07,NTOT,M_INX,N_INX,PSS,PDYN,XGSM,YGSM,ZGSM,
     *BXCF,BYCF,BZCF,
     *BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     *BXR11,BYR11,BZR11,BXR12,BYR12,BZR12,BXR21a,BYR21a,BZR21a,BXR21s,
     *BYR21s,BZR21s,BX,BY,BZ)

      RETURN

      END

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      SUBROUTINE EXTERN_2015 (IOPGEN,A,NTOT,
     *  M_INX,N_INX,PS,PDYN,X,Y,Z,
     *  BXCF,BYCF,BZCF,
     *  BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     *  BXR11,BYR11,BZR11,BXR12,BYR12,BZR12,
     *  BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s,
     *  BX,BY,BZ)
C
C   IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
C                                  IOPGEN=1 - DIPOLE SHIELDING ONLY
C                                  IOPGEN=2 - TAIL FIELD ONLY
C                                  IOPGEN=3 - BIRKELAND FIELD ONLY
C                                  IOPGEN=4 - RING CURRENT FIELD ONLY
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(NTOT)
C
      DIMENSION BXTS(5),BXTO(5,4),BXTE(5,4)
      DIMENSION BYTS(5),BYTO(5,4),BYTE(5,4)
      DIMENSION BZTS(5),BZTO(5,4),BZTE(5,4)
      INTEGER*4 M_INX,N_INX ! NOTE THAT THESE ARE FOR FUTURE USE
C
      COMMON /TAIL/ D  ! THE COMMON BLOCK FORWARDS TAIL SHEET THICKNESS
      COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  SCALING FACTORS FOR BIRKELAND CURRENTS
      COMMON /G/ G,TW
      COMMON /RH0/ RH0
C
      DATA A0_A,A0_S0,A0_X0 /34.586D0,1.1960D0,3.4397D0/   !   SHUE ET AL. PARAMETERS
      DATA DSIG /0.005D0/, RH2 /-5.2D0/

      XAPPA=(PDYN/2.D0)**0.155   !   0.155 is the value obtained in TS05
      XAPPA3=XAPPA**3

      D=      A(96)
      RH0=    A(97)
      G=      A(98)
      XKAPPA1=A(99)
      XKAPPA2=A(100)
      TW=     A(101)       !   THIS PARAMETER CONTROLS THE IMF-INDUCED TWISTING (ADDED 04/21/06)
c
      XX=X*XAPPA  ! pressure scaling has been reinstated here
      YY=Y*XAPPA
      ZZ=Z*XAPPA

c     print *,XAPPA,PDYN
C
      SPS=DSIN(PS)
c
      X0=A0_X0/XAPPA   ! pressure scaling has been reinstated, even though these parameters are not used in this code
      AM=A0_A/XAPPA    ! pressure scaling has been reinstated, even though these parameters are not used in this code
      S0=A0_S0
C
C   CALCULATE THE IMF CLOCK ANGLE:
C
C        IF (BYIMF.EQ.0.D0.AND.BZIMF.EQ.0.D0) THEN
C            THETA=0.D0
C         ELSE
C            THETA=DATAN2(BYIMF,BZIMF)
C            IF (THETA.LE.0.D0) THETA=THETA+6.283185307D0
C        ENDIF
C
C       CT=COS(THETA)
C       ST=SIN(THETA)
C       YS=Y*CT-Z*ST
C       ZS=Z*CT+Y*ST
C
C       STHETAH=SIN(THETA/2.)**2
C
C  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
C  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
C  OR OUTSIDE THE MAGNETOSPHERE:
C
C      FACTIMF=A(24)+A(25)*STHETAH
C
C      OIMFX=0.D0
C      OIMFY=BYIMF*FACTIMF
C      OIMFZ=BZIMF*FACTIMF
c
C =====================================================================
C  THIS FRAGMENT (BETWEEN THE ===== LINES) DISABLES THE CALCULATION OF THE MAGNETOPAUSE POSITION
C  IT SHOULD BE USED ONLY FOR THE FITTING (WE ASSUME THAT NO POINTS FROM THE SHEATH ARE PRESENT
C  IN THE DATASET, WHICH ITSELF IS STILL A QUESTION).
C
C  REMOVE IT IN THE FINAL VERSION.

C      SIGMA=0.D0
C      GOTO 1111
C======================================================================
c
C      R=SQRT(X**2+Y**2+Z**2)
C      XSS=X
C      ZSS=Z

C  1   XSOLD=XSS      !   BEGIN ITERATIVE SEARCH OF UNWARPED_2015 COORDS (TO FIND SIGMA)
C      ZSOLD=ZSS

C      RH=RH0+RH2*(ZSS/R)**2
C      SINPSAS=SPS/(1.D0+(R/RH)**3)**0.33333333D0
C      COSPSAS=DSQRT(1.D0-SINPSAS**2)
C      ZSS=X*SINPSAS+Z*COSPSAS
C      XSS=X*COSPSAS-Z*SINPSAS
C      DD=DABS(XSS-XSOLD)+DABS(ZSS-ZSOLD)
C      IF (DD.GT.1.D-6) GOTO 1
C                                END OF ITERATIVE SEARCH
C      RHO2=Y**2+ZSS**2
C      ASQ=AM**2
C      XMXM=AM+XSS-X0
C      IF (XMXM.LT.0.) XMXM=0. ! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
C      AXX0=XMXM**2
C      ARO=ASQ+RHO2
C      SIGMA=DSQRT((ARO+AXX0+SQRT((ARO+AXX0)**2-4.*ASQ*AXX0))/(2.*ASQ))

C==================================================================
C 1111 CONTINUE  !!!!!!!!!!!!  REMOVE IN THE FINAL VERSION
C==================================================================

C
C   NOW, THERE ARE THREE POSSIBLE CASES:
C    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
C    (2) IN THE BOUNDARY LAYER
C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      IF (SIGMA.LT.S0+DSIG) THEN  !  CASES (1) OR (2); CALCULATE THE MODEL FIELD
C                                   (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF (IOPGEN.LE.1) THEN

C     print *,XX,YY,ZZ
         CALL SHLCAR3X3_2015(XX,YY,ZZ,PS,CFX,CFY,CFZ)         !  DIPOLE SHIELDING FIELD
         BXCF=CFX  *XAPPA3
         BYCF=CFY  *XAPPA3
         BZCF=CFZ  *XAPPA3
      ELSE
         BXCF=0.D0
         BYCF=0.D0
         BZCF=0.D0
      ENDIF                                              !  DONE

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.2) THEN
      CALL DEFORMED_2015 (PS,XX,YY,ZZ,                !  TAIL FIELD (THREE MODES)
     *   BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE)
       ELSE
        DO 11 K=1,5

        BXTS(K)=0.D0
        BYTS(K)=0.D0
        BZTS(K)=0.D0

 11     CONTINUE

        DO 12 K=1,5

          DO 13 L=1,4

          BXTO(K,L)=0.D0
          BYTO(K,L)=0.D0
          BZTO(K,L)=0.D0

          BXTE(K,L)=0.D0
          BYTE(K,L)=0.D0
          BZTE(K,L)=0.D0

 13       CONTINUE

 12   CONTINUE

      ENDIF

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.3) THEN

         CALL BIRK_TOT_2015 (PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,BYR12,
     *   BZR12,BXR21a,BYR21a,BZR21a,BXR22a,BYR22a,BZR22a)    !   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)

         CALL BIRTOTSY_2015 (PS,XX,YY,ZZ,BXR11s,BYR11s,BZR11s,BXR12s,
     *   BYR12s,BZR12s,BXR21s,BYR21s,BZR21s,BXR22s,BYR22s,BZR22s)    !   "SYMMETRIC" BIRKELAND FIELD
C                                                                        (TWO MODES FOR R1s AND TWO MODES FOR R2s)
c                                                                        (but we actually use from here only R2s modes)
      ELSE
         BXR11=0.D0
         BYR11=0.D0
         BZR11=0.D0
         BXR12=0.D0
         BYR12=0.D0
         BZR12=0.D0
         BXR21a=0.D0
         BYR21a=0.D0
         BZR21a=0.D0
         BXR21s=0.D0
         BYR21s=0.D0
         BZR21s=0.D0
      ENDIF
C
C-----------------------------------------------------------
C
C    NOW, ADD UP ALL THE COMPONENTS:

      A_R11=A(92)
      A_R12=A(93)
      A_R21a=A(94)
      A_R21s=A(95)

      TX=0.D0
      TY=0.D0
      TZ=0.D0

C --- New tail structure -------------


        PDYN_0=2.D0   !   AVERAGE PRESSURE USED FOR NORMALIZATION

        P_FACTOR=DSQRT(PDYN/PDYN_0)-1.D0


      IND=1

        DO 911 K=1,5
         IND=IND+1
           TX=TX+(A(IND)+A(IND+45)*P_FACTOR)*BXTS(K)    !   2 - 6  &  47 - 51
           TY=TY+(A(IND)+A(IND+45)*P_FACTOR)*BYTS(K)
           TZ=TZ+(A(IND)+A(IND+45)*P_FACTOR)*BZTS(K)
 911     CONTINUE


        DO 912 K=1,5

          DO 913 L=1,4

          IND=IND+1

           TX=TX+(A(IND)+A(IND+45)*P_FACTOR)*BXTO(K,L)  !   7 -26  &  52 - 71
           TY=TY+(A(IND)+A(IND+45)*P_FACTOR)*BYTO(K,L)
           TZ=TZ+(A(IND)+A(IND+45)*P_FACTOR)*BZTO(K,L)

           TX=TX+(A(IND+20)+A(IND+65)*P_FACTOR)*BXTE(K,L) !   27 -46  &  72 - 91
           TY=TY+(A(IND+20)+A(IND+65)*P_FACTOR)*BYTE(K,L)
           TZ=TZ+(A(IND+20)+A(IND+65)*P_FACTOR)*BZTE(K,L)

 913       CONTINUE

 912   CONTINUE

      BBX=A(1)*BXCF+TX+
     * A_R11*BXR11+A_R12*BXR12+A_R21a*BXR21a+A_R21s*BXR21s

      BBY=A(1)*BYCF+TY+
     * A_R11*BYR11+A_R12*BYR12+A_R21a*BYR21a+A_R21s*BYR21s

      BBZ=A(1)*BZCF+TZ+
     * A_R11*BZR11+A_R12*BZR12+A_R21a*BZR21a+A_R21s*BZR21s
C
c   -----------------------------------------------------------
C
C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
C
C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
C
C      IF (SIGMA.LT.S0-DSIG) THEN    !  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE

C       BX=BBX
C       BY=BBY
C       BZ=BBZ
C                     ELSE           !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
C                                             THE INTERPOLATION REGION
C       FINT=0.5*(1.-(SIGMA-S0)/DSIG)
C       FEXT=0.5*(1.+(SIGMA-S0)/DSIG)
C
C       CALL DIPOLE (PS,X,Y,Z,QX,QY,QZ)
C       BX=(BBX+QX)*FINT+OIMFX*FEXT -QX
C       BY=(BBY+QY)*FINT+OIMFY*FEXT -QY
C       BZ=(BBZ+QZ)*FINT+OIMFZ*FEXT -QZ
c
C        ENDIF  !   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
C                      POSSIBILITY IS NOW THE CASE (3):
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        ELSE
C                CALL DIPOLE (PS,X,Y,Z,QX,QY,QZ)
C                BX=OIMFX-QX
C                BY=OIMFY-QY
C                BZ=OIMFZ-QZ
C        ENDIF
C

       BX=BBX
       BY=BBY
       BZ=BBZ

      RETURN
      END
C
C XXXXXXXXXXXXXXXXXXXXXXXXXXX11/15/05 16:06 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
         SUBROUTINE  SHLCAR3X3_2015(X,Y,Z,PS,BX,BY,BZ)
C
C   THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
C   REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
C   to the z=0 plane (see NB#4, p.74-74)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
c    harmonics (A(1)-A(36).
c  The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
C   entering the arguments of exponents, sines, and cosines in each of the
C   18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
C       (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(50)
      DATA A/-901.2327248,895.8011176,817.6208321,-845.5880889,
     *-83.73539535,86.58542841,336.8781402,-329.3619944,-311.2947120,
     *308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
     *-235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,
     *2.669338538,1.395023949,-.5540427503,-56.85224007,3.681827033,
     *-43.48705106,5.103131905,1.073551279,-.6673083508,12.21404266,
     *4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,
     *3.536082962,-3.222069852,9.620648151,6.082014949,27.75216226,
     *12.44199571,5.122226936,6.982039615,20.12149582,6.150973118,
     *4.663639687,15.73319647,2.303504968,5.840511214,.8385953499E-01,
     *.3477844929/
C
         P1=A(37)
         P2=A(38)
         P3=A(39)
         R1=A(40)
         R2=A(41)
         R3=A(42)
         Q1=A(43)
         Q2=A(44)
         Q3=A(45)
         S1=A(46)
         S2=A(47)
         S3=A(48)

         T1  =A(49)
         T2  =A(50)
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)
         S2PS=2.D0*CPS      !   MODIFIED HERE (INSTEAD OF SIN(3*PS) I TRY SIN(2*PS)

C
           ST1=DSIN(PS*T1)
           CT1=DCOS(PS*T1)
           ST2=DSIN(PS*T2)
           CT2=DCOS(PS*T2)

C     print *,X,Z

            X1=X*CT1-Z*ST1

C         print *,'X1=',X1

            Z1=X*ST1+Z*CT1
            X2=X*CT2-Z*ST2
            Z2=X*ST2+Z*CT2
C
C
c  MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
C
C       I=1
C
        SQPR= DSQRT(1.D0/P1**2+1.D0/R1**2)
        CYP = DCOS(Y/P1)
        SYP = DSIN(Y/P1)
        CZR = DCOS(Z1/R1)
        SZR = DSIN(Z1/R1)
c       print *,X1
        EXPR= DEXP(SQPR*X1)
        FX1 =-SQPR*EXPR*CYP*SZR
        HY1 = EXPR/P1*SYP*SZR
        FZ1 =-EXPR*CYP/R1*CZR
        HX1 = FX1*CT1+FZ1*ST1
        HZ1 =-FX1*ST1+FZ1*CT1

        SQPR= DSQRT(1.D0/P1**2+1.D0/R2**2)
        CYP = DCOS(Y/P1)
        SYP = DSIN(Y/P1)
        CZR = DCOS(Z1/R2)
        SZR = DSIN(Z1/R2)
        EXPR= DEXP(SQPR*X1)
        FX2 =-SQPR*EXPR*CYP*SZR
        HY2 = EXPR/P1*SYP*SZR
        FZ2 =-EXPR*CYP/R2*CZR
        HX2 = FX2*CT1+FZ2*ST1
        HZ2 =-FX2*ST1+FZ2*CT1

        SQPR= DSQRT(1.D0/P1**2+1.D0/R3**2)
        CYP = DCOS(Y/P1)
        SYP = DSIN(Y/P1)
        CZR = DCOS(Z1/R3)
        SZR = DSIN(Z1/R3)
        EXPR= DEXP(SQPR*X1)
        FX3 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
        HY3 = EXPR/P1*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
        FZ3 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
        HX3 = FX3*CT1+FZ3*ST1
        HZ3 =-FX3*ST1+FZ3*CT1
C
C       I=2:
C
        SQPR= DSQRT(1.D0/P2**2+1.D0/R1**2)
        CYP = DCOS(Y/P2)
        SYP = DSIN(Y/P2)
        CZR = DCOS(Z1/R1)
        SZR = DSIN(Z1/R1)
        EXPR= DEXP(SQPR*X1)
        FX4 =-SQPR*EXPR*CYP*SZR
        HY4 = EXPR/P2*SYP*SZR
        FZ4 =-EXPR*CYP/R1*CZR
        HX4 = FX4*CT1+FZ4*ST1
        HZ4 =-FX4*ST1+FZ4*CT1

        SQPR= DSQRT(1.D0/P2**2+1.D0/R2**2)
        CYP = DCOS(Y/P2)
        SYP = DSIN(Y/P2)
        CZR = DCOS(Z1/R2)
        SZR = DSIN(Z1/R2)
        EXPR= DEXP(SQPR*X1)
        FX5 =-SQPR*EXPR*CYP*SZR
        HY5 = EXPR/P2*SYP*SZR
        FZ5 =-EXPR*CYP/R2*CZR
        HX5 = FX5*CT1+FZ5*ST1
        HZ5 =-FX5*ST1+FZ5*CT1

        SQPR= DSQRT(1.D0/P2**2+1.D0/R3**2)
        CYP = DCOS(Y/P2)
        SYP = DSIN(Y/P2)
        CZR = DCOS(Z1/R3)
        SZR = DSIN(Z1/R3)
        EXPR= DEXP(SQPR*X1)
        FX6 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
        HY6 = EXPR/P2*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
        FZ6 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
        HX6 = FX6*CT1+FZ6*ST1
        HZ6 =-FX6*ST1+FZ6*CT1
C
C       I=3:
C
        SQPR= DSQRT(1.D0/P3**2+1.D0/R1**2)
        CYP = DCOS(Y/P3)
        SYP = DSIN(Y/P3)
        CZR = DCOS(Z1/R1)
        SZR = DSIN(Z1/R1)
        EXPR= DEXP(SQPR*X1)
        FX7 =-SQPR*EXPR*CYP*SZR
        HY7 = EXPR/P3*SYP*SZR
        FZ7 =-EXPR*CYP/R1*CZR
        HX7 = FX7*CT1+FZ7*ST1
        HZ7 =-FX7*ST1+FZ7*CT1

        SQPR= DSQRT(1.D0/P3**2+1.D0/R2**2)
        CYP = DCOS(Y/P3)
        SYP = DSIN(Y/P3)
        CZR = DCOS(Z1/R2)
        SZR = DSIN(Z1/R2)
        EXPR= DEXP(SQPR*X1)
        FX8 =-SQPR*EXPR*CYP*SZR
        HY8 = EXPR/P3*SYP*SZR
        FZ8 =-EXPR*CYP/R2*CZR
        HX8 = FX8*CT1+FZ8*ST1
        HZ8 =-FX8*ST1+FZ8*CT1

        SQPR= DSQRT(1.D0/P3**2+1.D0/R3**2)
        CYP = DCOS(Y/P3)
        SYP = DSIN(Y/P3)
        CZR = DCOS(Z1/R3)
        SZR = DSIN(Z1/R3)
        EXPR= DEXP(SQPR*X1)
        FX9 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
        HY9 = EXPR/P3*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
        FZ9 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
        HX9 = FX9*CT1+FZ9*ST1
        HZ9 =-FX9*ST1+FZ9*CT1


       A1=A(1)+A(2)*CPS
       A2=A(3)+A(4)*CPS
       A3=A(5)+A(6)*CPS
       A4=A(7)+A(8)*CPS
       A5=A(9)+A(10)*CPS
       A6=A(11)+A(12)*CPS
       A7=A(13)+A(14)*CPS
       A8=A(15)+A(16)*CPS
       A9=A(17)+A(18)*CPS
       BX=A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8+A9*HX9
       BY=A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8+A9*HY9
       BZ=A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8+A9*HZ9


c  MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
C
C       I=1
C
       SQQS= DSQRT(1.D0/Q1**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX1 =-SQQS*EXQS*CYQ*CZS *SPS
       HY1 = EXQS/Q1*SYQ*CZS   *SPS
       FZ1 = EXQS*CYQ/S1*SZS   *SPS
       HX1 = FX1*CT2+FZ1*ST2
       HZ1 =-FX1*ST2+FZ1*CT2

       SQQS= DSQRT(1.D0/Q1**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX2 =-SQQS*EXQS*CYQ*CZS *SPS
       HY2 = EXQS/Q1*SYQ*CZS   *SPS
       FZ2 = EXQS*CYQ/S2*SZS   *SPS
       HX2 = FX2*CT2+FZ2*ST2
       HZ2 =-FX2*ST2+FZ2*CT2

       SQQS= DSQRT(1.D0/Q1**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX3 =-SQQS*EXQS*CYQ*CZS *SPS
       HY3 = EXQS/Q1*SYQ*CZS   *SPS
       FZ3 = EXQS*CYQ/S3*SZS   *SPS
       HX3 = FX3*CT2+FZ3*ST2
       HZ3 =-FX3*ST2+FZ3*CT2
C
C       I=2
C
       SQQS= DSQRT(1.D0/Q2**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX4 =-SQQS*EXQS*CYQ*CZS *SPS
       HY4 = EXQS/Q2*SYQ*CZS   *SPS
       FZ4 = EXQS*CYQ/S1*SZS   *SPS
       HX4 = FX4*CT2+FZ4*ST2
       HZ4 =-FX4*ST2+FZ4*CT2

       SQQS= DSQRT(1.D0/Q2**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX5 =-SQQS*EXQS*CYQ*CZS *SPS
       HY5 = EXQS/Q2*SYQ*CZS   *SPS
       FZ5 = EXQS*CYQ/S2*SZS   *SPS
       HX5 = FX5*CT2+FZ5*ST2
       HZ5 =-FX5*ST2+FZ5*CT2

       SQQS= DSQRT(1.D0/Q2**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX6 =-SQQS*EXQS*CYQ*CZS *SPS
       HY6 = EXQS/Q2*SYQ*CZS   *SPS
       FZ6 = EXQS*CYQ/S3*SZS   *SPS
       HX6 = FX6*CT2+FZ6*ST2
       HZ6 =-FX6*ST2+FZ6*CT2
C
C       I=3
C
       SQQS= DSQRT(1.D0/Q3**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX7 =-SQQS*EXQS*CYQ*CZS *SPS
       HY7 = EXQS/Q3*SYQ*CZS   *SPS
       FZ7 = EXQS*CYQ/S1*SZS   *SPS
       HX7 = FX7*CT2+FZ7*ST2
       HZ7 =-FX7*ST2+FZ7*CT2

       SQQS= DSQRT(1.D0/Q3**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX8 =-SQQS*EXQS*CYQ*CZS *SPS
       HY8 = EXQS/Q3*SYQ*CZS   *SPS
       FZ8 = EXQS*CYQ/S2*SZS   *SPS
       HX8 = FX8*CT2+FZ8*ST2
       HZ8 =-FX8*ST2+FZ8*CT2

       SQQS= DSQRT(1.D0/Q3**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX9 =-SQQS*EXQS*CYQ*CZS *SPS
       HY9 = EXQS/Q3*SYQ*CZS   *SPS
       FZ9 = EXQS*CYQ/S3*SZS   *SPS
       HX9 = FX9*CT2+FZ9*ST2
       HZ9 =-FX9*ST2+FZ9*CT2

       A1=A(19)+A(20)*S2PS
       A2=A(21)+A(22)*S2PS
       A3=A(23)+A(24)*S2PS
       A4=A(25)+A(26)*S2PS
       A5=A(27)+A(28)*S2PS
       A6=A(29)+A(30)*S2PS
       A7=A(31)+A(32)*S2PS
       A8=A(33)+A(34)*S2PS
       A9=A(35)+A(36)*S2PS

       BX=BX+A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8
     *   +A9*HX9
       BY=BY+A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8
     *   +A9*HY9
       BZ=BZ+A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8
     *   +A9*HZ9
C
       RETURN
       END
c
c############################################################################
c
C
      SUBROUTINE DEFORMED_2015 (PS,X,Y,Z,
     *   BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)
C
C    CALCULATES GSM COMPONENTS OF 104 UNIT-AMPLITUDE TAIL FIELD MODES,
C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
C    WARPING IN Y-Z (DONE BY THE S/R WARPED_2015) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION BXS(5),BXO(5,4),BXE(5,4)
      DIMENSION BYS(5),BYO(5,4),BYE(5,4)
      DIMENSION BZS(5),BZO(5,4),BZE(5,4)
C
      DIMENSION BXASS(5),BXASO(5,4),BXASE(5,4)
      DIMENSION BYASS(5),BYASO(5,4),BYASE(5,4)
      DIMENSION BZASS(5),BZASO(5,4),BZASE(5,4)
C
      COMMON /RH0/ RH0
      DATA RH2,IEPS /-5.2D0,3/
C
C  RH0,RH1,RH2, AND IEPS CONTROL THE TILT-RELATED DEFORMATION OF THE TAIL FIELD
C
      SPS=DSIN(PS)
      CPS=DSQRT(1.D0-SPS**2)
      R2=X**2+Y**2+Z**2
      R=SQRT(R2)
      ZR=Z/R
      RH=RH0+RH2*ZR**2
      DRHDR=-ZR/R*2.D0*RH2*ZR
      DRHDZ= 2.D0*RH2*ZR/R
C
      RRH=R/RH
      F=1.D0/(1.D0+RRH**IEPS)**(1.D0/IEPS)
      DFDR=-RRH**(IEPS-1)*F**(IEPS+1)/RH
      DFDRH=-RRH*DFDR
c
      SPSAS=SPS*F
      CPSAS=DSQRT(1.D0-SPSAS**2)
C
      XAS=X*CPSAS-Z*SPSAS
      ZAS=X*SPSAS+Z*CPSAS
C
      FACPS=SPS/CPSAS*(DFDR+DFDRH*DRHDR)/R
      PSASX=FACPS*X
      PSASY=FACPS*Y
      PSASZ=FACPS*Z+SPS/CPSAS*DFDRH*DRHDZ
C
      DXASDX=CPSAS-ZAS*PSASX
      DXASDY=-ZAS*PSASY
      DXASDZ=-SPSAS-ZAS*PSASZ
      DZASDX=SPSAS+XAS*PSASX
      DZASDY=XAS*PSASY
      DZASDZ=CPSAS+XAS*PSASZ
      FAC1=DXASDZ*DZASDY-DXASDY*DZASDZ
      FAC2=DXASDX*DZASDZ-DXASDZ*DZASDX
      FAC3=DZASDX*DXASDY-DXASDX*DZASDY
C
C     DEFORM:
C
      CALL WARPED_2015(PS,XAS,Y,ZAS,
     *   BXASS,BYASS,BZASS,BXASO,BYASO,BZASO,BXASE,BYASE,BZASE)
C
C --- New tail structure -------------

        DO 11 K=1,5

        BXS(K)=BXASS(K)*DZASDZ-BZASS(K)*DXASDZ+BYASS(K)*FAC1
        BYS(K)=BYASS(K)*FAC2
        BZS(K)=BZASS(K)*DXASDX-BXASS(K)*DZASDX+BYASS(K)*FAC3

 11     CONTINUE

        DO 12 K=1,5

          DO 13 L=1,4

          BXO(K,L)=BXASO(K,L)*DZASDZ-BZASO(K,L)*DXASDZ
     *            +BYASO(K,L)*FAC1
          BYO(K,L)=BYASO(K,L)*FAC2
          BZO(K,L)=BZASO(K,L)*DXASDX-BXASO(K,L)*DZASDX
     *            +BYASO(K,L)*FAC3

          BXE(K,L)=BXASE(K,L)*DZASDZ-BZASE(K,L)*DXASDZ
     *            +BYASE(K,L)*FAC1
          BYE(K,L)=BYASE(K,L)*FAC2
          BZE(K,L)=BZASE(K,L)*DXASDX-BXASE(K,L)*DZASDX
     *            +BYASE(K,L)*FAC3

 13       CONTINUE

 12   CONTINUE

C ------------------------------------
C
      RETURN
      END
C
C------------------------------------------------------------------
c
C
      SUBROUTINE WARPED_2015 (PS,X,Y,Z,
     *   BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)
C
C   CALCULATES GSM COMPONENTS OF THE WARPED_2015 FIELD FOR TWO TAIL UNIT MODES.
C   THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED_2015 FIELD, COMPUTED
C   BY THE S/R "UNWARPED_2015".  THE WARPING PARAMETERS WERE TAKEN FROM THE
C   RESULTS OF GEOTAIL OBSERVATIONS (TSYGANENKO ET AL. [1998]).
C   NB # 6, P.106, OCT 12, 2000.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION BXS(5),BXO(5,4),BXE(5,4)
      DIMENSION BYS(5),BYO(5,4),BYE(5,4)
      DIMENSION BZS(5),BZO(5,4),BZE(5,4)
C
      DIMENSION BX_ASS(5),BX_ASO(5,4),BX_ASE(5,4)
      DIMENSION BY_ASS(5),BY_ASO(5,4),BY_ASE(5,4)
      DIMENSION BZ_ASS(5),BZ_ASO(5,4),BZ_ASE(5,4)
C
      COMMON /G/ G,TW
      DGDX=0.D0
      XL=20.D0
      DXLDX=0.D0

      SPS=DSIN(PS)
      RHO2=Y**2+Z**2
      RHO=DSQRT(RHO2)

      IF (Y.EQ.0.D0.AND.Z.EQ.0.D0) THEN
       PHI=0.D0
       CPHI=1.D0
       SPHI=0.D0
      ELSE
       PHI=DATAN2(Z,Y)
       CPHI=Y/RHO
       SPHI=Z/RHO
      ENDIF

      RR4L4=RHO/(RHO2**2+XL**4)

      F=PHI+G*RHO2*RR4L4*CPHI*SPS +TW*(X/10.D0)
      DFDPHI=1.D0-G*RHO2*RR4L4*SPHI*SPS
      DFDRHO=G*RR4L4**2*(3.D0*XL**4-RHO2**2)*CPHI*SPS
      DFDX=RR4L4*CPHI*SPS*(DGDX*RHO2-G*RHO*RR4L4*4.D0*XL**3*DXLDX)
     *  +TW/10.D0        !  THE LAST TERM DESCRIBES THE IMF-INDUCED TWISTING (ADDED 04/21/06)

      CF=DCOS(F)
      SF=DSIN(F)
      YAS=RHO*CF
      ZAS=RHO*SF

      CALL UNWARPED_2015 (X,YAS,ZAS,
     *   BX_ASS,BY_ASS,BZ_ASS,
     *   BX_ASO,BY_ASO,BZ_ASO,
     *   BX_ASE,BY_ASE,BZ_ASE)
C
        DO 11 K=1,5
C ------------------------------------------- Deforming symmetric modules
      BRHO_AS =  BY_ASS(K)*CF+BZ_ASS(K)*SF
      BPHI_AS = -BY_ASS(K)*SF+BZ_ASS(K)*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_ASS(K)*DFDX+BRHO_AS*DFDRHO)

        BXS(K)=BX_ASS(K)*DFDPHI
        BYS(K)=BRHO_S*CPHI-BPHI_S*SPHI
        BZS(K)=BRHO_S*SPHI+BPHI_S*CPHI

 11     CONTINUE

        DO 12 K=1,5

          DO 13 L=1,4
C -------------------------------------------- Deforming odd modules
      BRHO_AS =  BY_ASO(K,L)*CF+BZ_ASO(K,L)*SF
      BPHI_AS = -BY_ASO(K,L)*SF+BZ_ASO(K,L)*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_ASO(K,L)*DFDX+BRHO_AS*DFDRHO)

          BXO(K,L)=BX_ASO(K,L)*DFDPHI
          BYO(K,L)=BRHO_S*CPHI-BPHI_S*SPHI
          BZO(K,L)=BRHO_S*SPHI+BPHI_S*CPHI
C ------------------------------------------- Deforming even modules
      BRHO_AS =  BY_ASE(K,L)*CF+BZ_ASE(K,L)*SF
      BPHI_AS = -BY_ASE(K,L)*SF+BZ_ASE(K,L)*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_ASE(K,L)*DFDX+BRHO_AS*DFDRHO)

          BXE(K,L)=BX_ASE(K,L)*DFDPHI
          BYE(K,L)=BRHO_S*CPHI-BPHI_S*SPHI
          BZE(K,L)=BRHO_S*SPHI+BPHI_S*CPHI

 13       CONTINUE

 12   CONTINUE

      RETURN
      END
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE UNWARPED_2015 (X,Y,Z,BXS,BYS,BZS,BXO,BYO,
     * BZO,BXE,BYE,BZE)

C    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF 45 TAIL MODES WITH UNIT
C    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
C    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /TSS/ TSS(80,5)
      COMMON /TSO/ TSO(80,5,4)
      COMMON /TSE/ TSE(80,5,4)

      COMMON /TAIL/ D0

      DIMENSION BXS(5),BXO(5,4),BXE(5,4)
      DIMENSION BYS(5),BYO(5,4),BYE(5,4)
      DIMENSION BZS(5),BZO(5,4),BZE(5,4)
      
      !testing variables for Jay's correction
C      REAL*8 HXSKo,HYSKo,HZSKo,HXOKLo,HYOKLo,HZOKLo
C      REAL*8 HXEKLo,HYEKLo,HZEKLo,MR
C
C --- New tail structure -------------
C
C
C        MR=1e-7 !minimum resolution for testing Jay's optimization
        ! they are correct to ~1e-8 at worst
        DO 11 K=1,5

        CALL TAILSHT_S_2015 (K,X,Y,Z,BXSK,BYSK,BZSK)
C        CALL SHTBNORM_S (K,X,Y,Z,HXSKo,HYSKo,HZSKo) ! old routine
        CALL SHTBNORM_S_2015 (K,X,Y,Z,HXSK,HYSK,HZSK) !Jay's corrected routine

c       Below is how I went about testing the output (A.C.K.)
c       The old calls (without the j extension) remain for testing        

c        if ((HXSKo - HXSK).gt.MR.or.(HYSKo - HYSK).gt.MR.or.
c     *(HZSKo-HZSK).gt.MR) then
c            print *, 'BF(S) optimization from Jay are incorrect'
c            print *, HXSK,HXSKo,abs(HXSK-HXSKo)/HXSKo*100,'% err'
c            print *, HYSK,HYSKo,abs(HYSK-HYSKo)/HYSKo*100,'% err'
c            print *, HZSK,HZSKo,abs(HZSK-HZSKo)/HZSKo*100,'% err'
c            call sleep(1)
c        endif

        BXS(K)=BXSK+HXSK
        BYS(K)=BYSK+HYSK
        BZS(K)=BZSK+HZSK

 11     CONTINUE


        DO 12 K=1,5
          DO 13 L=1,4

          CALL TAILSHT_OE_2015 (1,K,L,X,Y,Z,BXOKL,BYOKL,BZOKL)
C          CALL SHTBNORM_O (  K,L,X,Y,Z,HXOKLo,HYOKLo,HZOKLo) !old
          CALL SHTBNORM_O_2015 (  K,L,X,Y,Z,HXOKL,HYOKL,HZOKL) !optimized

c        if ((HXOKLo - HXOKL).gt.MR .or. (HYOKLo - HYOKL).gt.MR .or.
c     * (HZOKLo-HZOKL).gt.MR) then
c            print *, 'BF(S) optimization from Jay are incorrect'
c            print *, HXOKL,HXOKLo,abs(HXOKL-HXOKLo)/HXOKLo*100,'% err'
c            print *, HYOKL,HYOKLo,abs(HYOKL-HYOKLo)/HYOKLo*100,'% err'
c            print *, HZOKL,HZOKLo,abs(HZOKL-HZOKLo)/HZOKLo*100,'% err'
c            call sleep(1)
c        endif

          BXO(K,L)=BXOKL+HXOKL
          BYO(K,L)=BYOKL+HYOKL
          BZO(K,L)=BZOKL+HZOKL

          CALL TAILSHT_OE_2015 (0,K,L,X,Y,Z,BXEKL,BYEKL,BZEKL)
c          CALL SHTBNORM_E (  K,L,X,Y,Z,HXEKLo,HYEKLo,HZEKLo) !old
          CALL SHTBNORM_E_2015 (  K,L,X,Y,Z,HXEKL,HYEKL,HZEKL) !optimized

c        if ((HXEKLo - HXEKL).gt.MR .or. (HYEKLo - HYEKL).gt.MR .or.
c     * (HZEKLo-HZEKL).gt.MR) then
c            print *, 'BF(S) optimization from Jay are incorrect'
c            print *, HXEKL,HXEKLo,abs(HXEKL-HXEKLo)/HXEKLo*100,'% err'
c            print *, HYEKL,HYEKLo,abs(HYEKL-HYEKLo)/HYEKLo*100,'% err'
c            print *, HZEKL,HZEKLo,abs(HZEKL-HZEKLo)/HZEKLo*100,'% err'
c            call sleep(1)
c        endif

          BXE(K,L)=BXEKL+HXEKL
          BYE(K,L)=BYEKL+HYEKL
          BZE(K,L)=BZEKL+HZEKL

 13       CONTINUE
 12   CONTINUE


      RETURN
      END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
c
C
      SUBROUTINE TAILSHT_S_2015 (M,X,Y,Z,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /TAIL/ D  ! THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS

C-----------------------------------------------------------------------------------
C
      RNOT=20.0        !    This can be replaced by introducing them
      DLTK=1.0         !    through the above common block


      RHO=DSQRT(X*X+Y*Y)
      CSPHI=X/RHO
      SNPHI=Y/RHO
C
      DKM=1.D0+(M-1)*DLTK
      RKM=DKM/RNOT
C
      RKMZ=RKM*Z
      RKMR=RKM*RHO

      ZD=DSQRT(Z*Z+D*D)
C
      RJ0=bessj0(RKMR)
      RJ1=bessj1(RKMR)
      REX=DEXP(RKM*ZD)
C
      BX=RKMZ*RJ1*CSPHI/ZD/REX
      BY=RKMZ*RJ1*SNPHI/ZD/REX
      BZ=RKM*RJ0/REX
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      RETURN
      END

C==========================================================================
C
       SUBROUTINE  SHTBNORM_S(K,X,Y,Z,FX,FY,FZ)
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
C
        DIMENSION AK(5)
C
        COMMON /TSS/ TSS(80,5)

c
         AK(1)=TSS(76,K)
         AK(2)=TSS(77,K)
         AK(3)=TSS(78,K)
         AK(4)=TSS(79,K)
         AK(5)=TSS(80,K)
C
c -------------------------------------------

          phi=DATAN2(Y,X)

c -------------------------------------------
c
             L=0

             FX=0.D0
             FY=0.D0
             FZ=0.D0
c
             DO 2 m1=1,15
                  m=m1-1
c
                  CMP=DCOS(m*phi)
                  SMP=DSIN(m*phi)
C
                DO 2 n=1,5
c ------------------------------------
                   RHO=dsqrt(X*X+Y*Y)
                   AKN=dabs(AK(n))
                   AKNR=AKN*RHO
c -----------------
                   CHZ=dcosh(Z*AKN)
                   SHZ=dsinh(Z*AKN)
c ------------------------------------
                   if(AKNR.lt.1.D-8) then
                   AKNRI=1.D8
                   else
                   AKNRI=1.D0/AKNR
                   end if
c -----------------
                   if(RHO.lt.1.D-8) then
                   RHOI=1.D8
                   else
                   RHOI=1.D0/RHO
                   end if
c -----------------
                   if(m.gt.2) then
                   AJM=bessj(m,AKNR)
                   AJM1=bessj(m-1,AKNR)
                   AJMD=AJM1-m*AJM*AKNRI
                   else
c                  ---------------------
                   if(m.eq.2) then
                   AJM=bessj(2,AKNR)
                   AJM1=bessj1(AKNR)
                   AJMD=AJM1-m*AJM*AKNRI
                   else
c                  ----------
                   if(m.eq.1) then
                   AJM=bessj1(AKNR)
                   AJM1=bessj0(AKNR)
                   AJMD=AJM1-AJM*AKNRI
c                  -----
                   else
                   AJM=bessj0(AKNR)
                   AJMD=-bessj1(AKNR)
c                  ----------
                   end if
c                  --------------------
                   end if
c                  --------------------
                   endif
c -----------------
                   DPDX=-Y*RHOI*RHOI
                   DPDY=X*RHOI*RHOI
c ------------------------------------
c
                   HX1=m*DPDX*SMP*SHZ*AJM
                   HX2=-AKN*X*RHOI*CMP*SHZ*AJMD
C
                   HX=HX1+HX2
c
                   HY1=m*DPDY*SMP*SHZ*AJM
                   HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD

                   HY=HY1+HY2
c
                   HZ=-AKN*CMP*CHZ*AJM
c ------------------------------------
c
              L=L+1
c
              FX=FX+HX*TSS(L,K)
              FY=FY+HY*TSS(L,K)
              FZ=FZ+HZ*TSS(L,K)

  2       CONTINUE
C
      RETURN
      END

C ===================================================================
c
C
      SUBROUTINE TAILSHT_OE_2015 (IEVO,MK,M,X,Y,Z,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /TAIL/ D0  ! THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS

C-----------------------------------------------------------------------------------
C
      RNOT=20.0     !    Rho_0 - scale parameter along the tail axis
      DLTK=1.0      !    step in Km
c
c -----------------------------------------------------------------------------------

      RHO=DSQRT(X*X+Y*Y)
C
      CSPHI=X/RHO
      SNPHI=Y/RHO
C
      phi=DATAN2(Y,X)
      CSMPHI=DCOS(m*phi)
      SNMPHI=DSIN(m*phi)
C
      DKM=1.D0+(MK-1)*DLTK
      RKM=DKM/RNOT
C
      RKMZ=RKM*Z
      RKMR=RKM*RHO
C
      ZD=DSQRT(Z*Z+D0*D0)
C
      REX=DEXP(RKM*ZD)
c
c ---- calculating Jm and its derivatives ------
c
                   if(m.gt.2) then
                   AJM=bessj(m,RKMR)
                   AJM1=bessj(m-1,RKMR)
                   AJMD=AJM1-m*AJM/RKMR
                   else
c                  --------------------
                   if(m.eq.2) then
                   AJM=bessj(2,RKMR)
                   AJM1=bessj1(RKMR)
                   AJMD=AJM1-m*AJM/RKMR
                   else
c                  --------------------
                   AJM=bessj1(RKMR)
                   AJM1=bessj0(RKMR)
                   AJMD=AJM1-AJM/RKMR
c                  --------------------
                   end if
c                  --------------------
                   endif
c -----------------------------------------
c
      if(ievo.eq.0) then
c -----------------------------------------
c calculating symmetric modes
c -----------------------------------------
c
      BRO=-M*SNMPHI*Z*AJMD/ZD/REX
      BPHI=-M*M*CSMPHI*Z*AJM/RKMR/ZD/REX
      BZ=M*SNMPHI*AJM/REX
c
c -----------------------------------------
      else
c -----------------------------------------
c calculating asymmetric modes
c -----------------------------------------
c
      BRO=M*CSMPHI*Z*AJMD/ZD/REX
      BPHI=-M*M*SNMPHI*Z*AJM/RKMR/ZD/REX
      BZ=-M*CSMPHI*AJM/REX
c
c -----------------------------------------
      end if
c
c --- transformation from cylindrical ccordinates to GSM ---
c
      BX=BRO*CSPHI-BPHI*SNPHI
      BY=BRO*SNPHI+BPHI*CSPHI
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      RETURN
      END

c ===========================================================================
C
         SUBROUTINE  SHTBNORM_O (K,L,X,Y,Z,FX,FY,FZ)
C
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
        IMPLICIT  REAL * 8  (A - H, O - Z)
C
        DIMENSION AK(5)
C
        COMMON /TSO/ TSO(80,5,4)

c
         AK(1)=TSO(76,K,L)
         AK(2)=TSO(77,K,L)
         AK(3)=TSO(78,K,L)
         AK(4)=TSO(79,K,L)
         AK(5)=TSO(80,K,L)
C
c -------------------------------------------

          phi=DATAN2(Y,X)

c -------------------------------------------
c
             L1=0
             FX=0.D0
             FY=0.D0
             FZ=0.D0
c
             DO 2 m1=1,15
                  m=m1-1
c
                  CMP=DCOS(m*phi)
                  SMP=DSIN(m*phi)
C
                DO 2 n=1,5
c ------------------------------------
                   RHO=dsqrt(X*X+Y*Y)
                   AKN=dabs(AK(n))
                   AKNR=AKN*RHO
c -----------------
                   CHZ=dcosh(Z*AKN)
                   SHZ=dsinh(Z*AKN)
c ------------------------------------
                   if(AKNR.lt.1.D-8) then
                   AKNRI=1.D8
                   else
                   AKNRI=1.D0/AKNR
                   end if
c -----------------
                   if(RHO.lt.1.D-8) then
                   RHOI=1.D8
                   else
                   RHOI=1.D0/RHO
                   end if
c -----------------
                   if(m.gt.2) then
                   AJM=bessj(m,AKNR)
                   AJM1=bessj(m-1,AKNR)
                   AJMD=AJM1-m*AJM*AKNRI
                   else
c                  ---------------------
                   if(m.eq.2) then
                   AJM=bessj(2,AKNR)
                   AJM1=bessj1(AKNR)
                   AJMD=AJM1-m*AJM*AKNRI
                   else
c                  ----------
                   if(m.eq.1) then
                   AJM=bessj1(AKNR)
                   AJM1=bessj0(AKNR)
                   AJMD=AJM1-AJM*AKNRI
c                  -----
                   else
                   AJM=bessj0(AKNR)
                   AJMD=-bessj1(AKNR)
c                  ----------
                   end if
c                  --------------------
                   end if
c                  --------------------
                   endif
c -----------------
                   DPDX=-Y*RHOI*RHOI
                   DPDY=X*RHOI*RHOI
c ------------------------------------
c
                   HX1=m*DPDX*SMP*SHZ*AJM
                   HX2=-AKN*X*RHOI*CMP*SHZ*AJMD
C
                   HX=HX1+HX2
c
                   HY1=m*DPDY*SMP*SHZ*AJM
                   HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD

                   HY=HY1+HY2
c
                   HZ=-AKN*CMP*CHZ*AJM
c ------------------------------------
c
              L1=L1+1
c
              FX=FX+HX*TSO(L1,K,L)
              FY=FY+HY*TSO(L1,K,L)
              FZ=FZ+HZ*TSO(L1,K,L)

  2       CONTINUE
C
      RETURN
      END

C
c ===========================================================================
C
         SUBROUTINE  SHTBNORM_E (K,L,X,Y,Z,FX,FY,FZ)
C
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
        IMPLICIT  REAL * 8  (A - H, O - Z)
C
        DIMENSION AK(5)
C
        COMMON /TSE/ TSE(80,5,4)

c
         AK(1)=TSE(76,K,L)
         AK(2)=TSE(77,K,L)
         AK(3)=TSE(78,K,L)
         AK(4)=TSE(79,K,L)
         AK(5)=TSE(80,K,L)
C
c -------------------------------------------

          phi=DATAN2(Y,X)

c -------------------------------------------
c
             L1=0
             FX=0.D0
             FY=0.D0
             FZ=0.D0
c
c ------------------------------------
             DO 2 m1=1,15
                  m=m1-1
c
                  CMP=DCOS(m*phi)
                  SMP=DSIN(m*phi)
C
                DO 2 n=1,5
c ------------------------------------
                   RHO=dsqrt(X*X+Y*Y)
                   AKN=dabs(AK(n))
                   AKNR=AKN*RHO
c -----------------
                   CHZ=dcosh(Z*AKN)
                   SHZ=dsinh(Z*AKN)
c ------------------------------------
                   if(AKNR.lt.1.D-8) then
                   AKNRI=1.D8
                   else
                   AKNRI=1.D0/AKNR
                   end if
c -----------------
                   if(RHO.lt.1.D-8) then
                   RHOI=1.D8
                   else
                   RHOI=1.D0/RHO
                   end if
c -----------------
                   if(m.gt.2) then
                   AJM=bessj(m,AKNR)
                   AJM1=bessj(m-1,AKNR)
                   AJMD=AJM1-m*AJM*AKNRI
                   else
c                  ---------------------
                   if(m.eq.2) then
                   AJM=bessj(2,AKNR)
                   AJM1=bessj1(AKNR)
                   AJMD=AJM1-m*AJM*AKNRI
                   else
c                  ----------
                   if(m.eq.1) then
                   AJM=bessj1(AKNR)
                   AJM1=bessj0(AKNR)
                   AJMD=AJM1-AJM*AKNRI
c                  -----
                   else
                   AJM=bessj0(AKNR)
                   AJMD=-bessj1(AKNR)
c                  ----------
                   end if
c                  --------------------
                   end if
c                  --------------------
                   endif
c -----------------
                   DPDX=-Y*RHOI*RHOI
                   DPDY=X*RHOI*RHOI
c ------------------------------------
c
                   HX1=-m*DPDX*CMP*SHZ*AJM
                   HX2=-AKN*X*RHOI*SMP*SHZ*AJMD
C
                   HX=HX1+HX2
c
                   HY1=-m*DPDY*CMP*SHZ*AJM
                   HY2=-AKN*Y*RHOI*SMP*SHZ*AJMD

                   HY=HY1+HY2
c
                   HZ=-AKN*SMP*CHZ*AJM

c
              L1=L1+1
c
              FX=FX+HX*TSE(L1,K,L)
              FY=FY+HY*TSE(L1,K,L)
              FZ=FZ+HZ*TSE(L1,K,L)

  2       CONTINUE
C
      RETURN
      END

C ===================================================================

C
      DOUBLE PRECISION FUNCTION bessj0(x)
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     *s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/
      if(Dabs(x).lt.8.D0)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=Dabs(x)
        z=8.D0/ax
        y=z**2
        xx=ax-.785398164D0
        bessj0=Dsqrt(.636619772/ax)*(Dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*Dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DOUBLE PRECISION FUNCTION bessj1(x)
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(Dabs(x).lt.8.D0)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=Dabs(x)
        z=8.D0/ax
        y=z**2
        xx=ax-2.356194491D0
        bessj1=Dsqrt(.636619772/ax)*(Dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*Dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*Dsign(1.D0,x)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.
c ===================================================================
      DOUBLE PRECISION FUNCTION bessj(n,x)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IACC=40,BIGNO=1.D10,BIGNI=1.D-10)
CU    USES bessj0,bessj1
      if(n.lt.2)pause 'bad argument n in bessj'
      ax=Dabs(x)
      if(ax.eq.0.D0)then
        bessj=0.
      else if(ax.gt.Dfloat(n))then
        tox=2.D0/ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        tox=2.D0/ax
        m=2*((n+int(Dsqrt(Dfloat(IACC*n))))/2)
        bessj=0.D0
        jsum=0
        sum=0.D0
        bjp=0.D0
        bj=1.D0
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(Dabs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.D0*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0.D0.and.mod(n,2).eq.1)bessj=-bessj
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE BIRK_TOT_2015 (PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,BZ12,
     *                          BX21,BY21,BZ21,BX22,BY22,BZ22)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
      COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.

      DATA SH11/46488.84663,-15541.95244,-23210.09824,-32625.03856,
     *-109894.4551,-71415.32808,58168.94612,55564.87578,-22890.60626,
     *-6056.763968,5091.368100,239.7001538,-13899.49253,4648.016991,
     *6971.310672,9699.351891,32633.34599,21028.48811,-17395.96190,
     *-16461.11037,7447.621471,2528.844345,-1934.094784,-588.3108359,
     *-32588.88216,10894.11453,16238.25044,22925.60557,77251.11274,
     *50375.97787,-40763.78048,-39088.60660,15546.53559,3559.617561,
     *-3187.730438,309.1487975,88.22153914,-243.0721938,-63.63543051,
     *191.1109142,69.94451996,-187.9539415,-49.89923833,104.0902848,
     *-120.2459738,253.5572433,89.25456949,-205.6516252,-44.93654156,
     *124.7026309,32.53005523,-98.85321751,-36.51904756,98.88241690,
     *24.88493459,-55.04058524,61.14493565,-128.4224895,-45.35023460,
     *105.0548704,-43.66748755,119.3284161,31.38442798,-92.87946767,
     *-33.52716686,89.98992001,25.87341323,-48.86305045,59.69362881,
     *-126.5353789,-44.39474251,101.5196856,59.41537992,41.18892281,
     *80.86101200,3.066809418,7.893523804,30.56212082,10.36861082,
     *8.222335945,19.97575641,2.050148531,4.992657093,2.300564232,
     *.2256245602,-.05841594319/

      DATA SH12/210260.4816,-1443587.401,-1468919.281,281939.2993,
     *-1131124.839,729331.7943,2573541.307,304616.7457,468887.5847,
     *181554.7517,-1300722.650,-257012.8601,645888.8041,-2048126.412,
     *-2529093.041,571093.7972,-2115508.353,1122035.951,4489168.802,
     *75234.22743,823905.6909,147926.6121,-2276322.876,-155528.5992,
     *-858076.2979,3474422.388,3986279.931,-834613.9747,3250625.781,
     *-1818680.377,-7040468.986,-414359.6073,-1295117.666,-346320.6487,
     *3565527.409,430091.9496,-.1565573462,7.377619826,.4115646037,
     *-6.146078880,3.808028815,-.5232034932,1.454841807,-12.32274869,
     *-4.466974237,-2.941184626,-.6172620658,12.64613490,1.494922012,
     *-21.35489898,-1.652256960,16.81799898,-1.404079922,-24.09369677,
     *-10.99900839,45.94237820,2.248579894,31.91234041,7.575026816,
     *-45.80833339,-1.507664976,14.60016998,1.348516288,-11.05980247,
     *-5.402866968,31.69094514,12.28261196,-37.55354174,4.155626879,
     *-33.70159657,-8.437907434,36.22672602,145.0262164,70.73187036,
     *85.51110098,21.47490989,24.34554406,31.34405345,4.655207476,
     *5.747889264,7.802304187,1.844169801,4.867254550,2.941393119,
     *.1379899178,.06607020029/

      DATA SH21/162294.6224,503885.1125,-27057.67122,-531450.1339,
     *84747.05678,-237142.1712,84133.61490,259530.0402,69196.05160,
     *-189093.5264,-19278.55134,195724.5034,-263082.6367,-818899.6923,
     *43061.10073,863506.6932,-139707.9428,389984.8850,-135167.5555,
     *-426286.9206,-109504.0387,295258.3531,30415.07087,-305502.9405,
     *100785.3400,315010.9567,-15999.50673,-332052.2548,54964.34639,
     *-152808.3750,51024.67566,166720.0603,40389.67945,-106257.7272,
     *-11126.14442,109876.2047,2.978695024,558.6019011,2.685592939,
     *-338.0004730,-81.99724090,-444.1102659,89.44617716,212.0849592,
     *-32.58562625,-982.7336105,-35.10860935,567.8931751,-1.917212423,
     *-260.2023543,-1.023821735,157.5533477,23.00200055,232.0603673,
     *-36.79100036,-111.9110936,18.05429984,447.0481000,15.10187415,
     *-258.7297813,-1.032340149,-298.6402478,-1.676201415,180.5856487,
     *64.52313024,209.0160857,-53.85574010,-98.52164290,14.35891214,
     *536.7666279,20.09318806,-309.7349530,58.54144539,67.45226850,
     *97.92374406,4.752449760,10.46824379,32.91856110,12.05124381,
     *9.962933904,15.91258637,1.804233877,6.578149088,2.515223491,
     *.1930034238,-.02261109942/

      DATA SH22/-131287.8986,-631927.6885,-318797.4173,616785.8782,
     *-50027.36189,863099.9833,47680.20240,-1053367.944,-501120.3811,
     *-174400.9476,222328.6873,333551.7374,-389338.7841,-1995527.467,
     *-982971.3024,1960434.268,297239.7137,2676525.168,-147113.4775,
     *-3358059.979,-2106979.191,-462827.1322,1017607.960,1039018.475,
     *520266.9296,2627427.473,1301981.763,-2577171.706,-238071.9956,
     *-3539781.111,94628.16420,4411304.724,2598205.733,637504.9351,
     *-1234794.298,-1372562.403,-2.646186796,-31.10055575,2.295799273,
     *19.20203279,30.01931202,-302.1028550,-14.78310655,162.1561899,
     *.4943938056,176.8089129,-.2444921680,-100.6148929,9.172262228,
     *137.4303440,-8.451613443,-84.20684224,-167.3354083,1321.830393,
     *76.89928813,-705.7586223,18.28186732,-770.1665162,-9.084224422,
     *436.3368157,-6.374255638,-107.2730177,6.080451222,65.53843753,
     *143.2872994,-1028.009017,-64.22739330,547.8536586,-20.58928632,
     *597.3893669,10.17964133,-337.7800252,159.3532209,76.34445954,
     *84.74398828,12.76722651,27.63870691,32.69873634,5.145153451,
     *6.310949163,6.996159733,1.971629939,4.436299219,2.904964304,
     *.1486276863,.06859991529/

C ====   LEAST SQUARES FITTING ONLY:
C       BX11=0.D0
C       BY11=0.D0
C       BZ11=0.D0
C       BX12=0.D0
C       BY12=0.D0
C       BZ12=0.D0
C       BX21=0.D0
C       BY21=0.D0
C       BZ21=0.D0
C       BX22=0.D0
C       BY22=0.D0
C       BZ22=0.D0
C===================================

      XKAPPA=XKAPPA1        !  FORWARDED IN BIRK_1N2_2015
      X_SC=XKAPPA1-1.1D0    !  FORWARDED IN BIRK_SHL_2015

      CALL BIRK_1N2_2015 (1,1,PS,X,Y,Z,FX11,FY11,FZ11)           !  REGION 1, MODE 1
      CALL BIRK_SHL_2015 (SH11,PS,X_SC,X,Y,Z,HX11,HY11,HZ11)
      BX11=FX11+HX11
      BY11=FY11+HY11
      BZ11=FZ11+HZ11

      CALL BIRK_1N2_2015 (1,2,PS,X,Y,Z,FX12,FY12,FZ12)           !  REGION 1, MODE 2
      CALL BIRK_SHL_2015 (SH12,PS,X_SC,X,Y,Z,HX12,HY12,HZ12)
      BX12=FX12+HX12
      BY12=FY12+HY12
      BZ12=FZ12+HZ12

      XKAPPA=XKAPPA2        !  FORWARDED IN BIRK_1N2_2015
      X_SC=XKAPPA2-1.0D0    !  FORWARDED IN BIRK_SHL_2015

      CALL BIRK_1N2_2015 (2,1,PS,X,Y,Z,FX21,FY21,FZ21)           !  REGION 2, MODE 1
      CALL BIRK_SHL_2015 (SH21,PS,X_SC,X,Y,Z,HX21,HY21,HZ21)
      BX21=FX21+HX21
      BY21=FY21+HY21
      BZ21=FZ21+HZ21

      CALL BIRK_1N2_2015 (2,2,PS,X,Y,Z,FX22,FY22,FZ22)           !  REGION 2, MODE 2
      CALL BIRK_SHL_2015 (SH22,PS,X_SC,X,Y,Z,HX22,HY22,HZ22)
      BX22=FX22+HX22
      BY22=FY22+HY22
      BZ22=FZ22+HZ22

      RETURN
      END
C
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
      SUBROUTINE BIRK_1N2_2015 (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)        !   SEE NB# 6, P.60
C
C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
C
C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A11(31),A12(31),A21(31),A22(31)
      COMMON /MODENUM/ M
      COMMON /DTHETA/ DTHETA

      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:

C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
C              TYPICAL VALUE: 0.06
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
C
      DATA BETA,RH,EPS/0.9D0,10.D0,3.D0/ ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

      DATA A11/.1618068350,-.1797957553,2.999642482,-.9322708978,
     *-.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
     *-16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
     *1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
     *-1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
     *.09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
     *2.492118385,.7113544659/
      DATA A12/.7058026940,-.2845938535,5.715471266,-2.472820880,
     *-.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
     *-212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
     *2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
     *-1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
     *.1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
     *1.212634762,.5567714182/
      DATA A21/.1278764024,-.2320034273,1.805623266,-32.37241440,
     *-.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
     *-6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
     *1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
     *-1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
     *.1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
     *1.102649543,.8867880020/
      DATA A22/.4036015198,-.3302974212,2.827730930,-45.44405830,
     *-1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
     *-233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
     *.7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
     *-1.460805289,.7719653528,-.6658988668,.2515179349E-05,
     *.02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
     *2.503482679,1.071587299,.7247997430/

      B=0.5
      RHO_0=7.0

      M=MODE
      IF (NUMB.EQ.1) THEN
          DPHI=0.055D0
          DTHETA=0.06D0
      ENDIF

      IF (NUMB.EQ.2) THEN
          DPHI=0.030D0
          DTHETA=0.09D0
      ENDIF

      Xsc=X*XKAPPA
      Ysc=Y*XKAPPA
      Zsc=Z*XKAPPA
      RHO=DSQRT(Xsc**2+Zsc**2)

      Rsc=DSQRT(Xsc**2+Ysc**2+Zsc**2)                                 !  SCALED
      RHO2=RHO_0**2

      IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
         PHI=0.D0
      ELSE
         PHI=DATAN2(-Zsc,Xsc)  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
      ENDIF

      SPHIC=DSIN(PHI)
      CPHIC=DCOS(PHI)  !  "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

      BRACK=DPHI+B*RHO2/(RHO2+1.D0)*(RHO**2-1.D0)/(RHO2+RHO**2)
      R1RH=(Rsc-1.D0)/RH
      PSIAS=BETA*PS/(1.D0+R1RH**EPS)**(1.D0/EPS)

      PHIS=PHI-BRACK*DSIN(PHI) -PSIAS
      DPHISPHI=1.D0-BRACK*DCOS(PHI)
      DPHISRHO=-2.D0*B*RHO2*RHO/(RHO2+RHO**2)**2 *DSIN(PHI)
     *   +BETA*PS*R1RH**(EPS-1.D0)*RHO/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))
      DPHISDY= BETA*PS*R1RH**(EPS-1.D0)*Ysc/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))

      SPHICS=DSIN(PHIS)
      CPHICS=DCOS(PHIS)

      XS= RHO*CPHICS
      ZS=-RHO*SPHICS

      IF (NUMB.EQ.1) THEN
        IF (MODE.EQ.1) CALL TWOCONES_2015 (A11,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES_2015 (A12,XS,Ysc,ZS,BXS,BYAS,BZS)
      ELSE
        IF (MODE.EQ.1) CALL TWOCONES_2015 (A21,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES_2015 (A22,XS,Ysc,ZS,BXS,BYAS,BZS)
      ENDIF

      BRHOAS=BXS*CPHICS-BZS*SPHICS
      BPHIAS=-BXS*SPHICS-BZS*CPHICS

      BRHO_S=BRHOAS*DPHISPHI                             *XKAPPA        ! SCALING
      BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) *XKAPPA
      BY_S=BYAS*DPHISPHI                                 *XKAPPA

      BX=BRHO_S*CPHIC-BPHI_S*SPHIC
      BY=BY_S
      BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC

      RETURN
      END
c
C=========================================================================
c
      SUBROUTINE TWOCONES_2015 (A,X,Y,Z,BX,BY,BZ)
C
C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
C

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)

      CALL ONE_CONE_2015 (A,X,Y,Z,BXN,BYN,BZN)
      CALL ONE_CONE_2015 (A,X,-Y,-Z,BXS,BYS,BZS)
      BX=BXN-BXS
      BY=BYN+BYS
      BZ=BZN+BZS

      RETURN
      END
c
C-------------------------------------------------------------------------
C
      SUBROUTINE ONE_CONE_2015(A,X,Y,Z,BX,BY,BZ)
c
c  RETURNS FIELD COMPONENTS FOR A DEFORMED_2015 CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
c    BY SIM_14.FOR.  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
c

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)

      COMMON /DTHETA/ DTHETA
      COMMON /MODENUM/ M

      DATA DR,DT/1.D-6,1.D-6/  !   JUST FOR NUMERICAL DIFFERENTIATION

      THETA0=A(31)

      RHO2=X**2+Y**2
      RHO=DSQRT(RHO2)
      R=DSQRT(RHO2+Z**2)
      THETA=DATAN2(RHO,Z)
      PHI=DATAN2(Y,X)
C
C   MAKE THE DEFORMATION OF COORDINATES:
C
       RS=R_S_2015(A,R,THETA)
       THETAS=THETA_S_2015(A,R,THETA)
       PHIS=PHI
C
C   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
C
c   MODE #M
       CALL FIALCOS_2015 (RS,THETAS,PHIS,BTAST,BFAST,M,THETA0,DTHETA) 
C
C   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
C
C      FIRST OF ALL, FIND THE DERIVATIVES:
C
       DRSDR=(R_S_2015(A,R+DR,THETA)-R_S_2015(A,R-DR,THETA))/(2.D0*DR)
       DRSDT=(R_S_2015(A,R,THETA+DT)-R_S_2015(A,R,THETA-DT))/(2.D0*DT)
       DTSDR=(THETA_S_2015(A,R+DR,THETA)-THETA_S_2015(A,R-DR,THETA))
     +  /(2.D0*DR)
   
       DTSDT=(THETA_S_2015(A,R,THETA+DT)-THETA_S_2015(A,R,THETA-DT))/
     + (2.D0*DT)
       STSST=DSIN(THETAS)/DSIN(THETA)
       RSR=RS/R

       BR     =-RSR/R*STSST*BTAST*DRSDT                 !   NB#6, P.43    BRAST DOES NOT ENTER HERE
       BTHETA = RSR*STSST*BTAST*DRSDR                  !               (SINCE IT IS ZERO IN OUR CASE)
       BPHI   = RSR*BFAST*(DRSDR*DTSDT-DRSDT*DTSDR)

       S=RHO/R
       C=Z/R
       SF=Y/RHO
       CF=X/RHO

       BE=BR*S+BTHETA*C

       BX=A(1)*(BE*CF-BPHI*SF)
       BY=A(1)*(BE*SF+BPHI*CF)
       BZ=A(1)*(BR*C-BTHETA*S)

       RETURN
       END
C
C=====================================================================================
      DOUBLE PRECISION FUNCTION R_S_2015(A,R,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
C
      R_S_2015=R+A(2)/R+A(3)*R/DSQRT(R**2+A(11)**2)+
     & A(4)*R/(R**2+A(12)**2)
     &+(A(5)+A(6)/R+A(7)*R/DSQRT(R**2+A(13)**2)+A(8)*R/(R**2+A(14)**2))*
     & DCOS(THETA)
     &+(A(9)*R/DSQRT(R**2+A(15)**2)+A(10)*R/(R**2+A(16)**2)**2)
     & *DCOS(2.D0*THETA)
C
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION THETA_S_2015(A,R,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
c
      THETA_S_2015=THETA+(A(17)+A(18)/R+A(19)/R**2
     *                +A(20)*R/DSQRT(R**2+A(27)**2))*DSIN(THETA)
     * +(A(21)+A(22)*R/DSQRT(R**2+A(28)**2)
     *                +A(23)*R/(R**2+A(29)**2))*DSIN(2.D0*THETA)
     * +(A(24)+A(25)/R+A(26)*R/(R**2+A(30)**2))*DSIN(3.D0*THETA)
C
      RETURN
      END
C
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      SUBROUTINE FIALCOS_2015(R,THETA,PHI,BTHETA,BPHI,N,THETA0,DT)
C
C  CONICAL MODEL OF BIRKELAND CURRENT FIELD; BASED ON THE OLD S/R FIALCO (OF 1990-91)
C  SEE THE OLD NOTEBOOK 1985-86-88, NOTE OF MARCH 5, BUT HERE BOTH INPUT AND OUTPUT ARE IN SPHERICAL CDS.

C  BTN, AND BPN ARE THE ARRAYS OF BTHETA AND BPHI (BTN(i), BPN(i) CORRESPOND TO i-th MODE).
C   ONLY FIRST  N  MODE AMPLITUDES ARE COMPUTED (N<=10).
C    THETA0 IS THE ANGULAR HALF-WIDTH OF THE CONE, DT IS THE ANGULAR H.-W. OF THE CURRENT LAYER

C   NOTE:  BR=0  (BECAUSE ONLY RADIAL CURRENTS ARE PRESENT IN THIS MODEL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  BTN(10),BPN(10),CCOS(10),SSIN(10)

      SINTE=DSIN(THETA)
      RO=R*SINTE
      COSTE=DCOS(THETA)
      SINFI=DSIN(PHI)
      COSFI=DCOS(PHI)
      TG=SINTE/(1.D0+COSTE)   !        TAN(THETA/2)
      CTG=SINTE/(1.D0-COSTE)  !        CTG(THETA/2)
C
C
      TETANP=THETA0+DT
      TETANM=THETA0-DT
      IF(THETA.LT.TETANM) GOTO 1
      TGP=DTAN(TETANP*0.5D0)
      TGM=DTAN(TETANM*0.5D0)
      TGM2=TGM*TGM
      TGP2=TGP*TGP
  1   CONTINUE

      COSM1=1.D0
      SINM1=0.D0
      TM=1.D0
      TGM2M=1.D0
      TGP2M=1.D0

      DO 2 M=1,N
      TM=TM*TG
      CCOS(M)=COSM1*COSFI-SINM1*SINFI
      SSIN(M)=SINM1*COSFI+COSM1*SINFI
      COSM1=CCOS(M)
      SINM1=SSIN(M)
      IF(THETA.LT.TETANM) THEN
      T=TM
      DTT=0.5D0*M*TM*(TG+CTG)
      DTT0=0.D0
      ELSE IF(THETA.LT.TETANP) THEN
      TGM2M=TGM2M*TGM2
      FC=1.D0/(TGP-TGM)
      FC1=1.D0/(2*M+1)
      TGM2M1=TGM2M*TGM
      TG21=1.D0+TG*TG
      T=FC*(TM*(TGP-TG)+FC1*(TM*TG-TGM2M1/TM))
      DTT=0.5D0*M*FC*TG21*(TM/TG*(TGP-TG)-FC1*(TM-TGM2M1/(TM*TG)))
      DTT0=0.5D0*FC*((TGP+TGM)*(TM*TG-FC1*(TM*TG-TGM2M1/TM))+
     * TM*(1.D0-TGP*TGM)-(1.D0+TGM2)*TGM2M/TM)
      ELSE
      TGP2M=TGP2M*TGP2
      TGM2M=TGM2M*TGM2
      FC=1.D0/(TGP-TGM)
      FC1=1.D0/(2*M+1)
      T=FC*FC1*(TGP2M*TGP-TGM2M*TGM)/TM
      DTT=-T*M*0.5D0*(TG+CTG)
      ENDIF

      BTN(M)=M*T*CCOS(M)/RO
  2   BPN(M)=-DTT*SSIN(M)/R

      BTHETA=BTN(N) *800.
      BPHI  =BPN(N) *800.

      RETURN
      END
C
C-------------------------------------------------------------------------
C
C
         SUBROUTINE BIRK_SHL_2015 (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
C
         IMPLICIT  REAL * 8  (A - H, O - Z)
         DIMENSION A(86)
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)

         S3PS=2.D0*CPS
C
         PST1=PS*A(85)
         PST2=PS*A(86)

         ST1=DSIN(PST1)
         CT1=DCOS(PST1)
         ST2=DSIN(PST2)
         CT2=DCOS(PST2)

         X1=X*CT1-Z*ST1
         Z1=X*ST1+Z*CT1
         X2=X*CT2-Z*ST2
         Z2=X*ST2+Z*CT2
C
         L=0
         GX=0.D0
         GY=0.D0
         GZ=0.D0
C
         DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(72+I)
                  Q=A(78+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
                DO 3 K=1,3
                   R=A(75+K)
                   S=A(81+K)
                   SZRK=DSIN(Z1/R)
                   CZSK=DCOS(Z2/S)
                   CZRK=DCOS(Z1/R)
                   SZSK=DSIN(Z2/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X1*SQPR)
                        EQS=DEXP(X2*SQQS)
C
                  DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                AND N=2 IS FOR THE SECOND ONE

                    DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                    IF (M.EQ.1) THEN
                         FX=-SQPR*EPR*CYPI*SZRK
                         FY=EPR*SYPI*SZRK/P
                         FZ=-EPR*CYPI*CZRK/R
                       IF (N.EQ.1) THEN
                         IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                         ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                         ENDIF
                       ELSE
                         IF (NN.EQ.1) THEN
                          HX=FX*CPS
                          HY=FY*CPS
                          HZ=FZ*CPS
                         ELSE
                          HX=FX*CPS*X_SC
                          HY=FY*CPS*X_SC
                          HZ=FZ*CPS*X_SC
                         ENDIF
                       ENDIF

                     ELSE                            !   M.EQ.2
                         FX=-SPS*SQQS*EQS*CYQI*CZSK
                         FY=SPS/Q*EQS*SYQI*CZSK
                         FZ=SPS/S*EQS*CYQI*SZSK
                       IF (N.EQ.1) THEN
                        IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                        ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                        ENDIF
                       ELSE
                        IF (NN.EQ.1) THEN
                         HX=FX*S3PS
                         HY=FY*S3PS
                         HZ=FZ*S3PS
                        ELSE
                         HX=FX*S3PS*X_SC
                         HY=FY*S3PS*X_SC
                         HZ=FZ*S3PS*X_SC
                        ENDIF
                       ENDIF
                  ENDIF
       L=L+1

       IF (M.EQ.1) THEN
       HXR=HX*CT1+HZ*ST1
       HZR=-HX*ST1+HZ*CT1
       ELSE
       HXR=HX*CT2+HZ*ST2
       HZR=-HX*ST2+HZ*CT2
       ENDIF

       GX=GX+HXR*A(L)
       GY=GY+HY *A(L)
  5    GZ=GZ+HZR*A(L)

  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE

      BX=GX
      BY=GY
      BZ=GZ

      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      SUBROUTINE BIRTOTSY_2015 (PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,BZ12,
     *                          BX21,BY21,BZ21,BX22,BY22,BZ22)
C
C   THIS S/R IS ALMOST IDENTICAL TO BIRK_TOT_2015, BUT IT IS FOR THE SYMMETRIC MODE, IN WHICH
C     J_parallel IS AN EVEN FUNCTION OF Ygsm.
C
C
C      IOPBS -  BIRKELAND FIELD MODE FLAG:
C         IOPBS=0 - ALL COMPONENTS
C         IOPBS=1 - REGION 1, MODES 1 & 2 (SYMMETRIC !)
C         IOPBS=2 - REGION 2, MODES 1 & 2 (SYMMETRIC !)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
      COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
C                                            (JOINT WITH  BIRK_TOT_2015  FOR THE ANTISYMMETRICAL MODE)
C
      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.
c
      DATA SH11/ 4956703.683,-26922641.21,-11383659.85,29604361.65,
     *-38919785.97,70230899.72,34993479.24,-90409215.02,30448713.69,
     *-48360257.19,-35556751.23,57136283.60,-8013815.613,30784907.86,
     * 13501620.50,-35121638.52,50297295.45,-84200377.18,-46946852.58,
     * 107526898.8,-39003263.47,59465850.17,47264335.10,-68892388.73,
     * 3375901.533,-9181255.754,-4494667.217,10812618.51,-17351920.97,
     * 27016083.00,18150032.11,-33186882.96,13340198.63,-19779685.30,
     * -17891788.15,21625767.23,16135.32442,133094.0241,-13845.61859,
     *-79159.98442,432.1215298,-85438.10368,1735.386707,41891.71284,
     * 18158.14923,-105465.8135,-11685.73823,62297.34252,-10811.08476,
     * -87631.38186,9217.499261,52079.94529,-68.29127454,56023.02269,
     * -1246.029857,-27436.42793,-11972.61726,69607.08725,7702.743803,
     * -41114.36810,12.08269108,-21.30967022,-9.100782462,18.26855933,
     * -7.000685929,26.22390883,6.392164144,-21.99351743,2.294204157,
     * -16.10023369,-1.344314750,9.342121230,148.5493329,99.79912328,
     * 70.78093196,35.23177574,47.45346891,58.44877918,139.8135237,
     * 91.96485261,6.983488815,9.055554871,19.80484284,2.860045019,
     * .8213262337E-01,-.7962186676E-05/
c
      DATA SH12/-1210748.720,-52324903.95,-14158413.33,19426123.60,
     * 6808641.947,-5138390.983,-1118600.499,-4675055.459,2059671.506,
     * -1373488.052,-114704.4353,-1435920.472,1438451.655,61199067.17,
     *  16549301.39,-22802423.47,-7814550.995,5986478.728,1299443.190,
     *  5352371.724,-2994351.520,1898553.337,203158.3658,2270182.134,
     * -618083.3112,-25950806.16,-7013783.326,9698792.575,3253693.134,
     * -2528478.464,-546323.4095,-2217735.237,1495336.589,-914647.4222,
     * -114374.1054,-1200441.634,-507068.4700,1163189.975,998411.8381,
     * -861919.3631,5252210.872,-11668550.16,-4113899.385,6972900.950,
     * -2546104.076,7704014.310,2273077.192,-5134603.198,256205.7901,
     * -589970.8086,-503821.0170,437612.8956,-2648640.128,5887640.735,
     *  2074286.234,-3519291.144,1283847.104,-3885817.147,-1145936.942,
     *  2589753.651,-408.7788403,1234.054185,739.8541716,-965.8068853,
     *  3691.383679,-8628.635819,-2855.844091,5268.500178,-1774.372703,
     *  5515.010707,1556.089289,-3665.434660,204.8672197,110.7748799,
     *  87.36036207,5.522491330,31.06364270,73.57632579,281.5331360,
     *  140.3461448,17.07537768,6.729732641,4.100970449,2.780422877,
     *  .8742978101E-01,-.1028562327E-04/
c
      DATA SH21/-67763516.61,-49565522.84,10123356.08,51805446.10,
     * -51607711.68,164360662.1,-4662006.024,-191297217.6,-7204547.103,
     * 30372354.93,-750371.9365,-36564457.17,61114395.65,45702536.50,
     * -9228894.939,-47893708.68,47290934.33,-149155112.0,4226520.638,
     * 173588334.5,7998505.443,-33150962.72,832493.2094,39892545.84,
     * -11303915.16,-8901327.398,1751557.110,9382865.820,-9054707.868,
     * 27918664.50,-788741.7146,-32481294.42,-2264443.753,9022346.503,
     * -233526.0185,-10856269.53,-244450.8850,1908295.272,185445.1967,
     * -1074202.863,41827.75224,-241553.7626,-20199.12580,123235.6084,
     * 199501.4614,-1936498.464,-178857.4074,1044724.507,121044.9917,
     * -946479.9247,-91808.28803,532742.7569,-20742.28628,120633.2193,
     * 10018.49534,-61599.11035,-98709.58977,959095.1770,88500.43489,
     * -517471.5287,-81.56122911,816.2472344,55.30711710,-454.5368824,
     * 25.74693810,-202.5007350,-7.369350794,104.9429812,58.14049362,
     * -685.5919355,-51.71345683,374.0125033,247.9296982,159.2471769,
     * 102.3151816,15.81062488,34.99767599,133.0832773,219.6475201,
     * 107.9582783,10.00264684,7.718306072,25.22866153,5.013583103,
     * .8407754233E-01,-.9613356793E-05/
c
      DATA SH22/-43404887.31,8896854.538,-8077731.036,-10247813.65,
     * 6346729.086,-9416801.212,-1921670.268,7805483.928,2299301.127,
     * 4856980.170,-1253936.462,-4695042.690,54305735.91,-11158768.10,
     * 10051771.85,12837129.47,-6380785.836,12387093.50,1687850.192,
     * -10492039.47,-5777044.862,-6916507.424,2855974.911,7027302.490,
     * -26176628.93,5387959.610,-4827069.106,-6193036.589,2511954.143,
     * -6205105.083,-553187.2984,5341386.847,3823736.361,3669209.068,
     * -1841641.700,-3842906.796,281561.7220,-5013124.630,379824.5943,
     * 2436137.901,-76337.55394,548518.2676,42134.28632,-281711.3841,
     * -365514.8666,-2583093.138,-232355.8377,1104026.712,-131536.3445,
     *  2320169.882,-174967.6603,-1127251.881,35539.82827,-256132.9284,
     * -19620.06116,131598.7965,169033.6708,1194443.500,107320.3699,
     * -510672.0036,1211.177843,-17278.19863,1140.037733,8347.612951,
     * -303.8408243,2405.771304,174.0634046,-1248.722950,-1231.229565,
     * -8666.932647,-754.0488385,3736.878824,227.2102611,115.9154291,
     * 94.34364830,3.625357304,64.03192907,109.0743468,241.4844439,
     * 107.7583478,22.36222385,6.282634037,27.79399216,2.270602235,
     * .8708605901E-01,-.1256706895E-04/
c
      XKAPPA=XKAPPA1        !  FORWARDED IN BIR1N2SY_2015
      X_SC=XKAPPA1-1.1D0    !  FORWARDED IN BIRSH_SY_2015

      CALL BIR1N2SY_2015 (1,1,PS,X,Y,Z,FX11,FY11,FZ11)           !  REGION 1, MODE 1
      CALL BIRSH_SY_2015 (SH11,PS,X_SC,X,Y,Z,HX11,HY11,HZ11)

      BX11=FX11+HX11
      BY11=FY11+HY11
      BZ11=FZ11+HZ11

      CALL BIR1N2SY_2015 (1,2,PS,X,Y,Z,FX12,FY12,FZ12)           !  REGION 1, MODE 2
      CALL BIRSH_SY_2015 (SH12,PS,X_SC,X,Y,Z,HX12,HY12,HZ12)

      BX12=FX12+HX12
      BY12=FY12+HY12
      BZ12=FZ12+HZ12

      XKAPPA=XKAPPA2        !  FORWARDED IN BIR1N2SY_2015
      X_SC=XKAPPA2-1.0D0    !  FORWARDED IN BIRSH_SY_2015

      CALL BIR1N2SY_2015 (2,1,PS,X,Y,Z,FX21,FY21,FZ21)           !  REGION 2, MODE 1
      CALL BIRSH_SY_2015 (SH21,PS,X_SC,X,Y,Z,HX21,HY21,HZ21)

      BX21=FX21+HX21
      BY21=FY21+HY21
      BZ21=FZ21+HZ21

      CALL BIR1N2SY_2015 (2,2,PS,X,Y,Z,FX22,FY22,FZ22)           !  REGION 2, MODE 2
      CALL BIRSH_SY_2015 (SH22,PS,X_SC,X,Y,Z,HX22,HY22,HZ22)

      BX22=FX22+HX22
      BY22=FY22+HY22
      BZ22=FZ22+HZ22

      RETURN
      END
C
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      SUBROUTINE BIR1N2SY_2015 (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)        !   SEE NB# 6, P.60 and NB#7, P.35-...
C
C   THIS CODE IS VERY SIMILAR TO BIRK_1N2_2015, BUT IT IS FOR THE "SYMMETRICAL" MODE, IN WHICH J_parallel
C     IS A SYMMETRIC (EVEN) FUNCTION OF Ygsm
C
C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
C
C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A11(31),A12(31),A21(31),A22(31)
      COMMON /MODENUM/ M
      COMMON /DTHETA/ DTHETA

      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:

C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
C              TYPICAL VALUE: 0.06
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
C
      DATA BETA,RH,EPS/0.9D0,10.D0,3.D0/ ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

      DATA A11/.1618068350,-.1797957553,2.999642482,-.9322708978,
     *-.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
     *-16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
     *1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
     *-1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
     *.09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
     *2.492118385,.7113544659/
      DATA A12/.7058026940,-.2845938535,5.715471266,-2.472820880,
     *-.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
     *-212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
     *2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
     *-1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
     *.1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
     *1.212634762,.5567714182/
      DATA A21/.1278764024,-.2320034273,1.805623266,-32.37241440,
     *-.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
     *-6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
     *1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
     *-1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
     *.1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
     *1.102649543,.8867880020/
      DATA A22/.4036015198,-.3302974212,2.827730930,-45.44405830,
     *-1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
     *-233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
     *.7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
     *-1.460805289,.7719653528,-.6658988668,.2515179349E-05,
     *.02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
     *2.503482679,1.071587299,.7247997430/

      B=0.5
      RHO_0=7.0

      M=MODE
      IF (NUMB.EQ.1) THEN
          DPHI=0.055D0
          DTHETA=0.06D0
      ENDIF

      IF (NUMB.EQ.2) THEN
          DPHI=0.030D0
          DTHETA=0.09D0
      ENDIF

      Xsc=X*XKAPPA
      Ysc=Y*XKAPPA
      Zsc=Z*XKAPPA
      RHO=DSQRT(Xsc**2+Zsc**2)

      Rsc=DSQRT(Xsc**2+Ysc**2+Zsc**2)                                 !  SCALED
      RHO2=RHO_0**2

      IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
         PHI=0.D0
      ELSE
         PHI=DATAN2(-Zsc,Xsc)  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
      ENDIF

      SPHIC=DSIN(PHI)
      CPHIC=DCOS(PHI)  !  "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

      BRACK=DPHI+B*RHO2/(RHO2+1.D0)*(RHO**2-1.D0)/(RHO2+RHO**2)
      R1RH=(Rsc-1.D0)/RH
      PSIAS=BETA*PS/(1.D0+R1RH**EPS)**(1.D0/EPS)

      PHIS=PHI-BRACK*DSIN(PHI) -PSIAS
      DPHISPHI=1.D0-BRACK*DCOS(PHI)
      DPHISRHO=-2.D0*B*RHO2*RHO/(RHO2+RHO**2)**2 *DSIN(PHI)
     *   +BETA*PS*R1RH**(EPS-1.D0)*RHO/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))
      DPHISDY= BETA*PS*R1RH**(EPS-1.D0)*Ysc/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))

      SPHICS=DSIN(PHIS)
      CPHICS=DCOS(PHIS)

      XS= RHO*CPHICS
      ZS=-RHO*SPHICS

      IF (NUMB.EQ.1) THEN
        IF (MODE.EQ.1) CALL TWOCONSS (A11,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONSS (A12,XS,Ysc,ZS,BXS,BYAS,BZS)
      ELSE
        IF (MODE.EQ.1) CALL TWOCONSS (A21,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONSS (A22,XS,Ysc,ZS,BXS,BYAS,BZS)
      ENDIF

      BRHOAS=BXS*CPHICS-BZS*SPHICS
      BPHIAS=-BXS*SPHICS-BZS*CPHICS

      BRHO_S=BRHOAS*DPHISPHI                             *XKAPPA        ! SCALING
      BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) *XKAPPA
      BY_S=BYAS*DPHISPHI                                 *XKAPPA

      BX=BRHO_S*CPHIC-BPHI_S*SPHIC
      BY=BY_S
      BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC

      RETURN
      END
c

C=========================================================================
c
      SUBROUTINE TWOCONSS (A,X,Y,Z,BX,BY,BZ)
C
C   DIFFERS FROM TWOCONES_2015:  THIS S/R IS FOR THE "SYMMETRIC" MODE OF BIRKELAND CURRENTS IN THAT
C                           HERE THE FIELD IS ROTATED BY 90 DEGS FOR M=1 AND BY 45 DEGS FOR M=2
C
C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
C

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
      COMMON /MODENUM/ M
      DATA HSQR2/0.707106781D0/

      IF (M.EQ.1) THEN   !   ROTATION BY 90 DEGS
       XAS = Y
       YAS =-X
      ELSE               !   ROTATION BY 45 DEGS
       XAS = (X+Y)*HSQR2
       YAS = (Y-X)*HSQR2
      ENDIF

      CALL ONE_CONE_2015 (A,XAS,YAS,Z,BXN,BYN,BZN)
      CALL ONE_CONE_2015 (A,XAS,-YAS,-Z,BXS,BYS,BZS)

      BXAS=BXN-BXS
      BYAS=BYN+BYS
      BZ=BZN+BZS

      IF (M.EQ.1) THEN   !   ROTATION BY 90 DEGS
        BX =-BYAS
        BY = BXAS
      ELSE
        BX=(BXAS-BYAS)*HSQR2
        BY=(BXAS+BYAS)*HSQR2
      ENDIF

      RETURN
      END
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
         SUBROUTINE BIRSH_SY_2015 (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
C
C   THIS S/R IS QUITE SIMILAR TO BIRK_SHL_2015, BUT IT IS FOR THE SYMMETRIC MODE OF BIRKELAND CURRENT FIELD
C     AND FOR THAT REASON THE FIELD COMPONENTS HAVE A DIFFERENT KIND OF SYMMETRY WITH RESPECT TO Y_gsm
C
         IMPLICIT  REAL * 8  (A - H, O - Z)
         DIMENSION A(86)
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)

         S3PS=2.D0*CPS
C
         PST1=PS*A(85)
         PST2=PS*A(86)

         ST1=DSIN(PST1)
         CT1=DCOS(PST1)
         ST2=DSIN(PST2)
         CT2=DCOS(PST2)

         X1=X*CT1-Z*ST1
         Z1=X*ST1+Z*CT1
         X2=X*CT2-Z*ST2
         Z2=X*ST2+Z*CT2
C
         L=0
         GX=0.D0
         GY=0.D0
         GZ=0.D0
C
         DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(72+I)
                  Q=A(78+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
                DO 3 K=1,3
                   R=A(75+K)
                   S=A(81+K)
                   SZRK=DSIN(Z1/R)
                   CZSK=DCOS(Z2/S)
                   CZRK=DCOS(Z1/R)
                   SZSK=DSIN(Z2/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X1*SQPR)
                        EQS=DEXP(X2*SQQS)
C
                  DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                AND N=2 IS FOR THE SECOND ONE

                    DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                    IF (M.EQ.1) THEN
                         FX= SQPR*EPR*SYPI*SZRK
                         FY=EPR*CYPI*SZRK/P
                         FZ= EPR*SYPI*CZRK/R
                       IF (N.EQ.1) THEN
                         IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                         ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                         ENDIF
                       ELSE
                         IF (NN.EQ.1) THEN
                          HX=FX*CPS
                          HY=FY*CPS
                          HZ=FZ*CPS
                         ELSE
                          HX=FX*CPS*X_SC
                          HY=FY*CPS*X_SC
                          HZ=FZ*CPS*X_SC
                         ENDIF
                       ENDIF

                     ELSE                            !   M.EQ.2
                         FX= SPS*SQQS*EQS*SYQI*CZSK
                         FY= SPS/Q*EQS*CYQI*CZSK
                         FZ=-SPS/S*EQS*SYQI*SZSK
                       IF (N.EQ.1) THEN
                        IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                        ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                        ENDIF
                       ELSE
                        IF (NN.EQ.1) THEN
                         HX=FX*S3PS
                         HY=FY*S3PS
                         HZ=FZ*S3PS
                        ELSE
                         HX=FX*S3PS*X_SC
                         HY=FY*S3PS*X_SC
                         HZ=FZ*S3PS*X_SC
                        ENDIF
                       ENDIF
                  ENDIF
       L=L+1

       IF (M.EQ.1) THEN
       HXR=HX*CT1+HZ*ST1
       HZR=-HX*ST1+HZ*CT1
       ELSE
       HXR=HX*CT2+HZ*ST2
       HZR=-HX*ST2+HZ*CT2
       ENDIF

       GX=GX+HXR*A(L)
       GY=GY+HY *A(L)
  5    GZ=GZ+HZR*A(L)

  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE

      BX=GX
      BY=GY
      BZ=GZ

      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
       SUBROUTINE DIPOLE_2015 (PS,X,Y,Z,BX,BY,BZ)
C
C      A DOUBLE PRECISION ROUTINE
C
C  CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH OF 2000.
C
C----INPUT PARAMETERS:
C     PS - GEODIPOLE TILT ANGLE IN RADIANS,
C     X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
C
C----OUTPUT PARAMETERS:
C     BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
C  LAST MODIFICATION: JAN. 5, 2001. THE VALUE OF THE DIPOLE MOMENT WAS UPDATED TO 2000.
C    AND A "SAVE" STATEMENT HAS BEEN ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
C    FORTRAN COMPILERS
C
C  WRITTEN BY: N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE M,PSI,SPS,CPS
      DATA M,PSI/0,5.D0/
      IF(M.EQ.1.AND.DABS(PS-PSI).LT.1.D-5) GOTO 1
      SPS=DSIN(PS)
      CPS=DCOS(PS)
      PSI=PS
      M=1
  1   P=X**2
      U=Z**2
      V=3.D0*Z*X
      T=Y**2
      Q=30115.D0/DSQRT(P+T+U)**5
      BX=Q*((T+U-2.D0*P)*SPS-V*CPS)
      BY=-3.D0*Y*Q*(X*SPS+Z*CPS)
      BZ=Q*((P+T-2.D0*U)*CPS-V*SPS)
      RETURN
      END


C(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

C
C(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
C=======================================================================================
C
      SUBROUTINE T96_MGNP_D (XN_PD,VEL,XGSM,YGSM,ZGSM,XMGNP,YMGNP,ZMGNP,
     * DIST,ID)

      IMPLICIT REAL*8 (A-H,O-Z)
C
C  DOUBLE-PRECISION VERSION !!!!!!!!   HENCE THE SUFFIX "D" IN THE NAME
C
C
C  FOR ANY POINT OF SPACE WITH GIVEN COORDINATES (XGSM,YGSM,ZGSM), THIS SUBROUTINE DEFINES
C  THE POSITION OF A POINT (XMGNP,YMGNP,ZMGNP) AT THE T96 MODEL MAGNETOPAUSE, HAVING THE
C  SAME VALUE OF THE ELLIPSOIDAL TAU-COORDINATE, AND THE DISTANCE BETWEEN THEM.  THIS IS
C  NOT THE SHORTEST DISTANCE D_MIN TO THE BOUNDARY, BUT DIST ASYMPTOTICALLY TENDS TO D_MIN,
C  AS THE OBSERVATION POINT GETS CLOSER TO THE MAGNETOPAUSE.
C
C  INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSM,YGSM,ZGSM - COORDINATES OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSM POSITION OF THE BOUNDARY POINT, HAVING THE SAME
C          VALUE OF TAU-COORDINATE AS THE OBSERVATION POINT (XGSM,YGSM,ZGSM)
C          DIST -  THE DISTANCE BETWEEN THE TWO POINTS, IN RE,
C          ID -    POSITION FLAG; ID=+1 (-1) MEANS THAT THE POINT (XGSM,YGSM,ZGSM)
C          LIES INSIDE (OUTSIDE) THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  THE PRESSURE-DEPENDENT MAGNETOPAUSE IS THAT USED IN THE T96_01 MODEL
C  (TSYGANENKO, JGR, V.100, P.5599, 1995; ESA SP-389, P.181, OCT. 1996)
C
c   AUTHOR:  N.A. TSYGANENKO
C   DATE:    AUG.1, 1995, REVISED APRIL 3, 2003.
C
C
C  DEFINE SOLAR WIND DYNAMIC PRESSURE (NANOPASCALS, ASSUMING 4% OF ALPHA-PARTICLES),
C   IF NOT EXPLICITLY SPECIFIED IN THE INPUT:

      IF (VEL.LT.0.D0) THEN
       PD=XN_PD
      ELSE
       PD=1.94D-6*XN_PD*VEL**2
C
      ENDIF
C
C  RATIO OF PD TO THE AVERAGE PRESSURE, ASSUMED EQUAL TO 2 nPa:

      RAT=PD/2.0D0
      RAT16=RAT**0.14D0

C (THE POWER INDEX 0.14 IN THE SCALING FACTOR IS THE BEST-FIT VALUE OBTAINED FROM DATA
C    AND USED IN THE T96_01 VERSION)
C
C  VALUES OF THE MAGNETOPAUSE PARAMETERS FOR  PD = 2 nPa:
C
      A0 =34.586D0
      S00=1.196D0
      X00=3.4397D0
C
C   VALUES OF THE MAGNETOPAUSE PARAMETERS, SCALED BY THE ACTUAL PRESSURE:
C
      A=A0/RAT16
      S0=S00
      X0=X00/RAT16
      XM=X0-A
C
C  (XM IS THE X-COORDINATE OF THE "SEAM" BETWEEN THE ELLIPSOID AND THE CYLINDER)
C
C     (FOR DETAILS ON THE ELLIPSOIDAL COORDINATES, SEE THE PAPER:
C      N.A.TSYGANENKO, SOLUTION OF CHAPMAN-FERRARO PROBLEM FOR AN
C      ELLIPSOIDAL MAGNETOPAUSE, PLANET.SPACE SCI., V.37, P.1037, 1989).
C
       IF (YGSM.NE.0.D0.OR.ZGSM.NE.0.D0) THEN
          PHI=DATAN2(YGSM,ZGSM)
       ELSE
          PHI=0.D0
       ENDIF
C
       RHO=DSQRT(YGSM**2+ZGSM**2)
C
       IF (XGSM.LT.XM) THEN
           XMGNP=XGSM
           RHOMGNP=A*DSQRT(S0**2-1.D0)
           YMGNP=RHOMGNP*DSIN(PHI)
           ZMGNP=RHOMGNP*DCOS(PHI)
           DIST=DSQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
           IF (RHOMGNP.GT.RHO) ID=+1
           IF (RHOMGNP.LE.RHO) ID=-1
           RETURN
       ENDIF
C
          XKSI=(XGSM-X0)/A+1.D0
          XDZT=RHO/A
          SQ1=DSQRT((1.D0+XKSI)**2+XDZT**2)
          SQ2=DSQRT((1.D0-XKSI)**2+XDZT**2)
          SIGMA=0.5D0*(SQ1+SQ2)
          TAU=0.5D0*(SQ1-SQ2)
C
C  NOW CALCULATE (X,Y,Z) FOR THE CLOSEST POINT AT THE MAGNETOPAUSE
C
          XMGNP=X0-A*(1.D0-S0*TAU)
          ARG=(S0**2-1.D0)*(1.D0-TAU**2)
          IF (ARG.LT.0.D0) ARG=0.D0
          RHOMGNP=A*DSQRT(ARG)
          YMGNP=RHOMGNP*DSIN(PHI)
          ZMGNP=RHOMGNP*DCOS(PHI)
C
C  NOW CALCULATE THE DISTANCE BETWEEN THE POINTS {XGSM,YGSM,ZGSM} AND {XMGNP,YMGNP,ZMGNP}:
C   (IN GENERAL, THIS IS NOT THE SHORTEST DISTANCE D_MIN, BUT DIST ASYMPTOTICALLY TENDS
C    TO D_MIN, AS WE ARE GETTING CLOSER TO THE MAGNETOPAUSE):
C
      DIST=DSQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
C
      IF (SIGMA.GT.S0) ID=-1   !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE
      IF (SIGMA.LE.S0) ID=+1   !  ID=+1 MEANS THAT THE POINT LIES INSIDE
C                                           THE MAGNETOSPHERE
      RETURN
      END
C

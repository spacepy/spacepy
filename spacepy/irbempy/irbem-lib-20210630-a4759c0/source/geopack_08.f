c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ c
c
c          ##########################################################################
c          #                                                                        #
c          #                          GEOPACK-2008_dp                               #
c          #                     (MAIN SET OF FORTRAN CODES)                        #
c          #               (double-precision version of 03/21/08)                   #
C          #                (IGRF coefficients updated 01/31/15)                    #
c          ##########################################################################
C
c
c  This collection of subroutines is a result of several upgrades of the original package
c  written by N. A. Tsyganenko in 1978-1979.
c
c  PREFATORY NOTE TO THE VERSION OF FEBRUARY 4, 2008:
c
c  To avoid inappropriate use of obsolete subroutines from earlier versions, a suffix 08 was
c  added to the name of each subroutine in this release.
c
c  A possibility has been added in this version to calculate vector components in the
c  "Geocentric Solar Wind" (GSW) coordinate system, which, to our knowledge, was first
c  introduced by Hones et al., Planet. Space Sci., v.34, p.889, 1986 (aka GSWM, see Appendix,
c  Tsyganenko et al., JGRA, v.103(A4), p.6827, 1998). The GSW system is analogous to the
c  standard GSM, except that its X-axis is antiparallel to the currently observed solar wind
c  flow vector, rather than aligned with the Earth-Sun line. The orientation of axes in the
c  GSW system can be uniquely defined by specifying three components (VGSEX,VGSEY,VGSEZ) of
c  the solar wind velocity, and in the case of a strictly radial anti-sunward flow (VGSEY=
c  VGSEZ=0) the GSW system becomes identical to the standard GSM, which fact was used here
c  to minimize the number of subroutines in the package. To that end, instead of the special
c  case of the GSM coordinates, this version uses a more general GSW system, and three more
c  input parameters are added in the subroutine RECALC_08, the observed values (VGSEX,VGSEY,
c  VGSEZ) of the solar wind velocity. Invoking RECALC_08 with VGSEY=VGSEZ=0 restores the
c  standard (sunward) orientation of the X axis, which allows one to easily convert vectors
c  between GSW and GSM, as well as to/from other standard and commonly used systems. For more
c  details, see the documentation file GEOPACK-2008.DOC.
c
c  Another modification allows users to have more control over the procedure of field line
c  mapping using the subroutine TRACE_08. To that end, three new input parameters were added
c  in that subroutine, allowing one to set (i) an upper limit, DSMAX, on the automatically
c  adjusted step size, (ii) a permissible step error, ERR, and (iii) maximal length, LMAX,
c  of arrays where field line point coordinates are stored. Minor changes were also made in
c  the tracing subroutine, to make it more compact and easier for understanding, and to
c  prevent the algorithm from making uncontrollable large number of multiple loops in some
c  cases with plasmoid-like field structures.
c
C  One more subroutine, named GEODGEO_08, was added to the package, allowing one to convert
c  geodetic coordinates of a point in space (altitude above the Earth's WGS84 ellipsoid and
c  geodetic latitude) to geocentric radial distance and colatitude, and vice versa.
c
C  For a complete list of modifications made earlier in previous versions, see the
c  documentation file GEOPACK-2008.DOC.
c
c----------------------------------------------------------------------------------
c
      SUBROUTINE IGRF_GSW_08 (XGSW,YGSW,ZGSW,HXGSW,HYGSW,HZGSW)
c
C  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE GEOCENTRIC SOLAR-WIND
C  (GSW) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL COEFFICIENTS
C  (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, revised 22 March, 2005)
c
C  THE GSW SYSTEM IS ESSENTIALLY SIMILAR TO THE STANDARD GSM (THE TWO SYSTEMS BECOME IDENTICAL
C  TO EACH OTHER IN THE CASE OF STRICTLY ANTI-SUNWARD SOLAR WIND FLOW). FOR A DETAILED
C  DEFINITION, SEE INTRODUCTORY COMMENTS FOR THE SUBROUTINE GSWGSE_08 .
C
C  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR, IF THE DATE/TIME (IYEAR,IDAY,IHOUR,MIN,ISEC),
C  OR THE SOLAR WIND VELOCITY COMPONENTS (VGSEX,VGSEY,VGSEZ) HAVE CHANGED, THE MODEL COEFFICIENTS
c  AND GEO-GSW ROTATION MATRIX ELEMENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC_08.
C
C-----INPUT PARAMETERS:
C
C     XGSW,YGSW,ZGSW - CARTESIAN GEOCENTRIC SOLAR-WIND COORDINATES (IN UNITS RE=6371.2 KM)
C
C-----OUTPUT PARAMETERS:
C
C     HXGSW,HYGSW,HZGSW - CARTESIAN GEOCENTRIC SOLAR-WIND COMPONENTS OF THE MAIN GEOMAGNETIC
C                           FIELD IN NANOTESLA
C
C     LAST MODIFICATION:  MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C     THIS VERSION OF THE CODE ACCEPTS DATES FROM 1965 THROUGH 2015.
c
C     AUTHOR: N. A. TSYGANENKO
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)

      DIMENSION A(14),B(14)

      CALL GEOGSW_08 (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW,-1)
      RHO2=XGEO**2+YGEO**2
      R=DSQRT(RHO2+ZGEO**2)
      C=ZGEO/R
      RHO=DSQRT(RHO2)
      S=RHO/R
      IF (S.LT.1.D-10) THEN
        CF=1.D0
        SF=0.D0
      ELSE
        CF=XGEO/RHO
        SF=YGEO/RHO
      ENDIF
C
      PP=1.D0/R
      P=PP
C
C  IN THIS VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
C    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
C      ON THE VALUE OF THE RADIAL DISTANCE R:
C
      IRP3=R+2
      NM=3+30/IRP3
      IF (NM.GT.13) NM=13

      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N

      P=1.D0
      D=0.D0
      BBR=0.D0
      BBT=0.D0
      BBF=0.D0

      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.D0
         Y=1.D0
170      Q=P
         Z=D
         BI=0.D0
         P2=0.D0
         D2=0.D0
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.D-10) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BT=BBT
      IF(S.LT.1.D-10) GOTO 210
      BF=BBF/S
      GOTO 211
210   IF(C.LT.0.) BBF=-BBF
      BF=BBF

211   HE=BR*S+BT*C
      HXGEO=HE*CF-BF*SF
      HYGEO=HE*SF+BF*CF
      HZGEO=BR*C-BT*S
C
      CALL GEOGSW_08 (HXGEO,HYGEO,HZGEO,HXGSW,HYGSW,HZGSW,1)
C
      RETURN
      END
C
c==========================================================================================
C
c
      SUBROUTINE IGRF_GEO_08 (R,THETA,PHI,BR,BTHETA,BPHI)
c
C  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE SPHERICAL GEOGRAPHIC
C  (GEOCENTRIC) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL
C  COEFFICIENTS  (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, revised: 22 March, 2005)
C
C  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR IF THE DATE (IYEAR AND IDAY) WAS CHANGED,
C  THE MODEL COEFFICIENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC_08
C
C-----INPUT PARAMETERS:
C
C   R, THETA, PHI - SPHERICAL GEOGRAPHIC (GEOCENTRIC) COORDINATES:
C   RADIAL DISTANCE R IN UNITS RE=6371.2 KM, COLATITUDE THETA AND LONGITUDE PHI IN RADIANS
C
C-----OUTPUT PARAMETERS:
C
C     BR, BTHETA, BPHI - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
C      (POSITIVE BR OUTWARD, BTHETA SOUTHWARD, BPHI EASTWARD)
C
C     LAST MODIFICATION:  MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C     THIS VERSION OF THE  CODE ACCEPTS DATES FROM 1965 THROUGH 2015.
c
C     AUTHOR: N. A. TSYGANENKO
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)

      DIMENSION A(14),B(14)

      C=DCOS(THETA)
      S=DSIN(THETA)
      CF=DCOS(PHI)
      SF=DSIN(PHI)
C
      PP=1.D0/R
      P=PP
C
C  IN THIS NEW VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
C    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
C      ON THE VALUE OF THE RADIAL DISTANCE R:
C
      IRP3=R+2
      NM=3+30/IRP3
      IF (NM.GT.13) NM=13

      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N

      P=1.D0
      D=0.D0
      BBR=0.D0
      BBT=0.D0
      BBF=0.D0

      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.D0
         Y=1.D0
170      Q=P
         Z=D
         BI=0.D0
         P2=0.D0
         D2=0.D0
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.D-5) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BTHETA=BBT
      IF(S.LT.1.D-10) GOTO 210
      BPHI=BBF/S
      RETURN
210   IF(C.LT.0.D0) BBF=-BBF
      BPHI=BBF

      RETURN
      END
C
c==========================================================================================
c
       SUBROUTINE DIP_08 (XGSW,YGSW,ZGSW,BXGSW,BYGSW,BZGSW)
C
C  CALCULATES GSW (GEOCENTRIC SOLAR-WIND) COMPONENTS OF GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH, SPECIFIED BY CALLING SUBROUTINE RECALC_08 (SHOULD BE
C  INVOKED BEFORE THE FIRST USE OF THIS ONE, OR IF THE DATE/TIME, AND/OR THE OBSERVED
C  SOLAR WIND DIRECTION, HAVE CHANGED.
C
C  THE GSW COORDINATE SYSTEM IS ESSENTIALLY SIMILAR TO THE STANDARD GSM (THE TWO SYSTEMS BECOME
C  IDENTICAL TO EACH OTHER IN THE CASE OF STRICTLY RADIAL ANTI-SUNWARD SOLAR WIND FLOW). ITS
C  DETAILED DEFINITION IS GIVEN IN INTRODUCTORY COMMENTS FOR THE SUBROUTINE GSWGSE_08 .

C--INPUT PARAMETERS: XGSW,YGSW,ZGSW - GSW COORDINATES IN RE (1 RE = 6371.2 km)
C
C--OUTPUT PARAMETERS: BXGSW,BYGSW,BZGSW - FIELD COMPONENTS IN GSW SYSTEM, IN NANOTESLA.
C
C  LAST MODIFICATION:   MARCH 21, 2008 (DOUBLE-PRECISION VERSION).
C
C  AUTHOR: N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ AA(10),SPS,CPS,BB(22)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)
C
      DIPMOM=DSQRT(G(2)**2+G(3)**2+H(3)**2)
C
      P=XGSW**2
      U=ZGSW**2
      V=3.D0*ZGSW*XGSW
      T=YGSW**2
      Q=DIPMOM/DSQRT(P+T+U)**5
      BXGSW=Q*((T+U-2.D0*P)*SPS-V*CPS)
      BYGSW=-3.D0*YGSW*Q*(XGSW*SPS+ZGSW*CPS)
      BZGSW=Q*((P+T-2.D0*U)*CPS-V*SPS)
      RETURN
      END

C*******************************************************************
c
      SUBROUTINE SUN_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C  LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA RAD/57.295779513D0/
C
      IF(IYEAR.LT.1901.OR.IYEAR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.D0
      VL=DMOD(279.696678D0+0.9856473354D0*DJ,360.D0)
      GST=DMOD(279.690983D0+.9856473354D0*DJ+360.D0*FDAY+180.D0,360.D0)/
     * RAD
      G=DMOD(358.475845D0+0.985600267D0*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946D0-0.004789D0*T)*DSIN(G)+0.020094D0
     * *DSIN(2.D0*G))/RAD
      IF(SLONG.GT.6.2831853D0) SLONG=SLONG-6.283185307D0
      IF (SLONG.LT.0.D0) SLONG=SLONG+6.283185307D0
      OBLIQ=(23.45229D0-0.0130125D0*T)/RAD
      SOB=DSIN(OBLIQ)
      SLP=SLONG-9.924D-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION DUE TO
C   EARTH'S ORBITAL MOTION
C
      SIND=SOB*DSIN(SLP)
      COSD=DSQRT(1.D0-SIND**2)
      SC=SIND/COSD
      SDEC=DATAN(SC)
      SRASN=3.141592654D0-DATAN2(DCOS(OBLIQ)/SOB*SC,-DCOS(SLP)/COSD)
      RETURN
      END
C
C================================================================================
c
      SUBROUTINE SPHCAR_08 (R,THETA,PHI,X,Y,Z,J)
C
C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICE VERSA
C    (THETA AND PHI IN RADIANS).
C
C                  J>0            J<0
C-----INPUT:   J,R,THETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,THETA,PHI
C
C  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 WHEN CONVERTING
C        FROM CARTESIAN TO SPHERICAL COORDS (I.E., FOR J<0)
C
C   LAST MOFIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C   AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)

      IF(J.GT.0) GOTO 3
      SQ=X**2+Y**2
      R=DSQRT(SQ+Z**2)
      IF (SQ.NE.0.D0) GOTO 2
      PHI=0.D0
      IF (Z.LT.0.D0) GOTO 1
      THETA=0.D0
      RETURN
  1   THETA=3.141592654D0
      RETURN
  2   SQ=DSQRT(SQ)
      PHI=DATAN2(Y,X)
      THETA=DATAN2(SQ,Z)
      IF (PHI.LT.0.D0) PHI=PHI+6.283185307D0
      RETURN
  3   SQ=R*DSIN(THETA)
      X=SQ*DCOS(PHI)
      Y=SQ*DSIN(PHI)
      Z=R*DCOS(THETA)
      RETURN
      END
C
C===========================================================================
c
      SUBROUTINE BSPCAR_08 (THETA,PHI,BR,BTHETA,BPHI,BX,BY,BZ)
C
C   CALCULATES CARTESIAN FIELD COMPONENTS FROM LOCAL SPHERICAL ONES
C
C-----INPUT:   THETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
C              BR,BTHETA,BPHI -  LOCAL SPHERICAL COMPONENTS OF THE FIELD
C-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
C
C   LAST MOFIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)

      S=DSIN(THETA)
      C=DCOS(THETA)
      SF=DSIN(PHI)
      CF=DCOS(PHI)
      BE=BR*S+BTHETA*C
      BX=BE*CF-BPHI*SF
      BY=BE*SF+BPHI*CF
      BZ=BR*C-BTHETA*S
      RETURN
      END
c
C==============================================================================
C
      SUBROUTINE BCARSP_08 (X,Y,Z,BX,BY,BZ,BR,BTHETA,BPHI)
C
CALCULATES LOCAL SPHERICAL FIELD COMPONENTS FROM THOSE IN CARTESIAN SYSTEM
C
C-----INPUT:   X,Y,Z  - CARTESIAN COMPONENTS OF THE POSITION VECTOR
C              BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD VECTOR
C-----OUTPUT:  BR,BTHETA,BPHI - LOCAL SPHERICAL COMPONENTS OF THE FIELD VECTOR
C
C  NOTE: AT THE POLES (THETA=0 OR THETA=PI) WE ASSUME PHI=0,
C        AND HENCE BTHETA=BX, BPHI=BY
C
C   WRITTEN AND ADDED TO THIS PACKAGE:  APRIL 1, 2003
c   LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C   AUTHOR:   N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)

      RHO2=X**2+Y**2
      R=DSQRT(RHO2+Z**2)
      RHO=DSQRT(RHO2)

      IF (RHO.NE.0.D0) THEN
        CPHI=X/RHO
        SPHI=Y/RHO
       ELSE
        CPHI=1.D0
        SPHI=0.D0
      ENDIF

      CT=Z/R
      ST=RHO/R

      BR=(X*BX+Y*BY+Z*BZ)/R
      BTHETA=(BX*CPHI+BY*SPHI)*CT-BZ*ST
      BPHI=BY*CPHI-BX*SPHI

      RETURN
      END
C
c=====================================================================================
C
      SUBROUTINE RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ)
C
C  1. PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
C     SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.
C
C  2. PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD
C      (IGRF MODEL)
C
C  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
C  IGRF_GEO_08, IGRF_GSW_08, DIP_08, GEOMAG_08, GEOGSW_08, MAGSW_08, SMGSW_08, GSWGSE_08,
c  GEIGEO_08, TRACE_08, STEP_08, RHAND_08.
C
C  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC_08, IF MULTIPLE CALCULATIONS ARE MADE
C    FOR THE SAME DATE/TIME AND SOLAR WIND FLOW DIRECTION.
C
C-----INPUT PARAMETERS:
C
C     IYEAR   -  YEAR NUMBER (FOUR DIGITS)
C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
C     IHOUR -  HOUR OF DAY (00 TO 23)
C     MIN   -  MINUTE OF HOUR (00 TO 59)
C     ISEC  -  SECONDS OF MINUTE (00 TO 59)
C     VGSEX,VGSEY,VGSEZ - GSE (GEOCENTRIC SOLAR-ECLIPTIC) COMPONENTS OF THE OBSERVED
C                              SOLAR WIND FLOW VELOCITY (IN KM/S)
C
C  IMPORTANT: IF ONLY QUESTIONABLE INFORMATION (OR NO INFORMATION AT ALL) IS AVAILABLE
C             ON THE SOLAR WIND SPEED, OR, IF THE STANDARD GSM AND/OR SM COORDINATES ARE
C             INTENDED TO BE USED, THEN SET VGSEX=-400.0 AND VGSEY=VGSEZ=0. IN THIS CASE,
C             THE GSW COORDINATE SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM.
C
C             IF ONLY SCALAR SPEED V OF THE SOLAR WIND IS KNOWN, THEN SETTING
C             VGSEX=-V, VGSEY=29.78, VGSEZ=0.0 WILL TAKE INTO ACCOUNT THE ~4 degs
C             ABERRATION OF THE MAGNETOSPHERE DUE TO EARTH'S ORBITAL MOTION
C
C             IF ALL THREE GSE COMPONENTS OF THE SOLAR WIND VELOCITY ARE AVAILABLE,
C             PLEASE NOTE THAT IN SOME SOLAR WIND DATABASES THE ABERRATION EFFECT
C             HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29.78 KM/S FROM VYGSE;
C             IN THAT CASE, THE UNABERRATED (OBSERVED) VYGSE VALUES SHOULD BE RESTORED
C             BY ADDING BACK THE 29.78 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST
C             BE EITHER VERIFIED WITH THE DATA ORIGINATOR OR DETERMINED BY AVERAGING
C             VGSEY OVER A SUFFICIENTLY LONG TIME INTERVAL.
C
C-----OUTPUT PARAMETERS:  NONE (ALL OUTPUT QUANTITIES ARE PLACED
C                         INTO THE COMMON BLOCKS /GEOPACK1/ AND /GEOPACK2/)
C
C    OTHER SUBROUTINES CALLED BY THIS ONE: SUN_08
C
C    AUTHOR:  N.A. TSYGANENKO
C    DATE:    DEC.1, 1991
C
C    REVISION OF JANUARY 31, 2015:
C
C     The table of IGRF coefficients was extended to include those for the epoch 2015 (igrf-12)
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
C
C    REVISION OF NOVEMBER 30, 2010:
C
C     The table of IGRF coefficients was extended to include those for the epoch 2010
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
C
C    REVISION OF NOVEMBER 15, 2007: ADDED THE POSSIBILITY TO TAKE INTO ACCOUNT THE OBSERVED
C     DEFLECTION OF THE SOLAR WIND FLOW FROM STRICTLY RADIAL DIRECTION. TO THAT END, THREE
C     GSE COMPONENTS OF THE SOLAR WIND VELOCITY WERE ADDED TO THE INPUT PARAMETERS.
C
c    CORRECTION OF MAY 9, 2006:  INTERPOLATION OF THE COEFFICIENTS (BETWEEN
C     LABELS 50 AND 105) IS NOW MADE THROUGH THE LAST ELEMENT OF THE ARRAYS
C     G(105)  AND H(105) (PREVIOUSLY MADE ONLY THROUGH N=66, WHICH IN SOME
C     CASES CAUSED RUNTIME ERRORS)
c
C    REVISION OF MAY 3, 2005:
C     The table of IGRF coefficients was extended to include those for the epoch 2005
c       the maximal order of spherical harmonics was also increased up to 13
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
c
C    REVISION OF APRIL 3, 2003:
c    The code now includes preparation of the model coefficients for the subroutines
c    IGRF_08 and GEOMAG_08. This eliminates the need for the SAVE statements, used
c    in the old versions, making the codes easier and more compiler-independent.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE ISW
C
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,
     * SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33,
     * E11,E21,E31,E12,E22,E32,E13,E23,E33
C
C  THE COMMON BLOCK /GEOPACK1/ CONTAINS ELEMENTS OF THE ROTATION MATRICES AND OTHER
C   PARAMETERS RELATED TO THE COORDINATE TRANSFORMATIONS PERFORMED BY THIS PACKAGE
C
      COMMON /GEOPACK2/ G(105),H(105),REC(105)
C
C  THE COMMON BLOCK /GEOPACK2/ CONTAINS COEFFICIENTS OF THE IGRF FIELD MODEL, CALCULATED
C    FOR A GIVEN YEAR AND DAY FROM THEIR STANDARD EPOCH VALUES. THE ARRAY REC CONTAINS
C    COEFFICIENTS USED IN THE RECURSION RELATIONS FOR LEGENDRE ASSOCIATE POLYNOMIALS.
C
      DIMENSION G65(105),H65(105),G70(105),H70(105),G75(105),H75(105),
     + G80(105),H80(105),G85(105),H85(105),G90(105),H90(105),G95(105),
     + H95(105),G00(105),H00(105),G05(105),H05(105),G10(105),H10(105),
     + G15(105),H15(105),DG15(45),DH15(45)
C
      DATA ISW /0/
c
      DATA G65/0.D0,-30334.D0,-2119.D0,-1662.D0,2997.D0,1594.D0,1297.D0,
     *-2038.D0,1292.D0,856.D0,957.D0,804.D0,479.D0,-390.D0,252.D0,
     *-219.D0,358.D0,254.D0,-31.D0,-157.D0,-62.D0,45.D0,61.D0,8.D0,
     *-228.D0,4.D0,1.D0,-111.D0,75.D0,-57.D0,4.D0,13.D0,-26.D0,-6.D0,
     *13.D0,1.D0,13.D0,5.D0,-4.D0,-14.D0,0.D0,8.D0,-1.D0,11.D0,4.D0,
     *8.D0,10.D0,2.D0,-13.D0,10.D0,-1.D0,-1.D0,5.D0,1.D0,-2.D0,-2.D0,
     *-3.D0,2.D0,-5.D0,-2.D0,4.D0,4.D0,0.D0,2.D0,2.D0,0.D0,39*0.D0/
      DATA H65/0.D0,0.D0,5776.D0,0.D0,-2016.D0,114.D0,0.D0,-404.D0,
     *240.D0,-165.D0,0.D0,148.D0,-269.D0,13.D0,-269.D0,0.D0,19.D0,
     *128.D0,-126.D0,-97.D0,81.D0,0.D0,-11.D0,100.D0,68.D0,-32.D0,-8.D0,
     *-7.D0,0.D0,-61.D0,-27.D0,-2.D0,6.D0,26.D0,-23.D0,-12.D0,0.D0,7.D0,
     *-12.D0,9.D0,-16.D0,4.D0,24.D0,-3.D0,-17.D0,0.D0,-22.D0,15.D0,7.D0,
     *-4.D0,-5.D0,10.D0,10.D0,-4.D0,1.D0,0.D0,2.D0,1.D0,2.D0,6.D0,-4.D0,
     *0.D0,-2.D0,3.D0,0.D0,-6.D0,39*0.D0/
c
      DATA G70/0.D0,-30220.D0,-2068.D0,-1781.D0,3000.D0,1611.D0,1287.D0,
     *-2091.D0,1278.D0,838.D0,952.D0,800.D0,461.D0,-395.D0,234.D0,
     *-216.D0,359.D0,262.D0,-42.D0,-160.D0,-56.D0,43.D0,64.D0,15.D0,
     *-212.D0,2.D0,3.D0,-112.D0,72.D0,-57.D0,1.D0,14.D0,-22.D0,-2.D0,
     *13.D0,-2.D0,14.D0,6.D0,-2.D0,-13.D0,-3.D0,5.D0,0.D0,11.D0,3.D0,
     *8.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,0.D0,3.D0,1.D0,-1.D0,-3.D0,
     *-3.D0,2.D0,-5.D0,-1.D0,6.D0,4.D0,1.D0,0.D0,3.D0,-1.D0,39*0.D0/
      DATA H70/0.D0,0.D0,5737.D0,0.D0,-2047.D0,25.D0,0.D0,-366.D0,
     *251.D0,-196.D0,0.D0,167.D0,-266.D0,26.D0,-279.D0,0.D0,26.D0,
     *139.D0,-139.D0,-91.D0,83.D0,0.D0,-12.D0,100.D0,72.D0,-37.D0,-6.D0,
     *1.D0,0.D0,-70.D0,-27.D0,-4.D0,8.D0,23.D0,-23.D0,-11.D0,0.D0,7.D0,
     *-15.D0,6.D0,-17.D0,6.D0,21.D0,-6.D0,-16.D0,0.D0,-21.D0,16.D0,6.D0,
     *-4.D0,-5.D0,10.D0,11.D0,-2.D0,1.D0,0.D0,1.D0,1.D0,3.D0,4.D0,-4.D0,
     *0.D0,-1.D0,3.D0,1.D0,-4.D0,39*0.D0/
c
      DATA G75/0.D0,-30100.D0,-2013.D0,-1902.D0,3010.D0,1632.D0,1276.D0,
     *-2144.D0,1260.D0,830.D0,946.D0,791.D0,438.D0,-405.D0,216.D0,
     *-218.D0,356.D0,264.D0,-59.D0,-159.D0,-49.D0,45.D0,66.D0,28.D0,
     *-198.D0,1.D0,6.D0,-111.D0,71.D0,-56.D0,1.D0,16.D0,-14.D0,0.D0,
     *12.D0,-5.D0,14.D0,6.D0,-1.D0,-12.D0,-8.D0,4.D0,0.D0,10.D0,1.D0,
     *7.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,-1.D0,4.D0,1.D0,-2.D0,-3.D0,
     *-3.D0,2.D0,-5.D0,-2.D0,5.D0,4.D0,1.D0,0.D0,3.D0,-1.D0,39*0.D0/
C
      DATA H75/0.D0,0.D0,5675.D0,0.D0,-2067.D0,-68.D0,0.D0,-333.D0,
     *262.D0,-223.D0,0.D0,191.D0,-265.D0,39.D0,-288.D0,0.D0,31.D0,
     *148.D0,-152.D0,-83.D0,88.D0,0.D0,-13.D0,99.D0,75.D0,-41.D0,-4.D0,
     *11.D0,0.D0,-77.D0,-26.D0,-5.D0,10.D0,22.D0,-23.D0,-12.D0,0.D0,
     *6.D0,-16.D0,4.D0,-19.D0,6.D0,18.D0,-10.D0,-17.D0,0.D0,-21.D0,
     *16.D0,7.D0,-4.D0,-5.D0,10.D0,11.D0,-3.D0,1.D0,0.D0,1.D0,1.D0,3.D0,
     *4.D0,-4.D0,-1.D0,-1.D0,3.D0,1.D0,-5.D0,39*0.D0/
c
      DATA G80/0.D0,-29992.D0,-1956.D0,-1997.D0,3027.D0,1663.D0,1281.D0,
     *-2180.D0,1251.D0,833.D0,938.D0,782.D0,398.D0,-419.D0,199.D0,
     *-218.D0,357.D0,261.D0,-74.D0,-162.D0,-48.D0,48.D0,66.D0,42.D0,
     *-192.D0,4.D0,14.D0,-108.D0,72.D0,-59.D0,2.D0,21.D0,-12.D0,1.D0,
     *11.D0,-2.D0,18.D0,6.D0,0.D0,-11.D0,-7.D0,4.D0,3.D0,6.D0,-1.D0,
     *5.D0,10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,2.D0,-5.D0,-4.D0,
     *-4.D0,2.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,3.D0,0.D0,39*0.D0/
C
      DATA H80/0.D0,0.D0,5604.D0,0.D0,-2129.D0,-200.D0,0.D0,-336.D0,
     *271.D0,-252.D0,0.D0,212.D0,-257.D0,53.D0,-297.D0,0.D0,46.D0,
     *150.D0,-151.D0,-78.D0,92.D0,0.D0,-15.D0,93.D0,71.D0,-43.D0,-2.D0,
     *17.D0,0.D0,-82.D0,-27.D0,-5.D0,16.D0,18.D0,-23.D0,-10.D0,0.D0,
     *7.D0,-18.D0,4.D0,-22.D0,9.D0,16.D0,-13.D0,-15.D0,0.D0,-21.D0,
     *16.D0,9.D0,-5.D0,-6.D0,9.D0,10.D0,-6.D0,2.D0,0.D0,1.D0,0.D0,3.D0,
     *6.D0,-4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0,39*0.D0/
c
      DATA G85/0.D0,-29873.D0,-1905.D0,-2072.D0,3044.D0,1687.D0,1296.D0,
     *-2208.D0,1247.D0,829.D0,936.D0,780.D0,361.D0,-424.D0,170.D0,
     *-214.D0,355.D0,253.D0,-93.D0,-164.D0,-46.D0,53.D0,65.D0,51.D0,
     *-185.D0,4.D0,16.D0,-102.D0,74.D0,-62.D0,3.D0,24.D0,-6.D0,4.D0,
     *10.D0,0.D0,21.D0,6.D0,0.D0,-11.D0,-9.D0,4.D0,4.D0,4.D0,-4.D0,5.D0,
     *10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,1.D0,-5.D0,-4.D0,-4.D0,
     *3.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,3.D0,0.D0,39*0.D0/
C
      DATA H85/0.D0,0.D0,5500.D0,0.D0,-2197.D0,-306.D0,0.D0,-310.D0,
     *284.D0,-297.D0,0.D0,232.D0,-249.D0,69.D0,-297.D0,0.D0,47.D0,
     *150.D0,-154.D0,-75.D0,95.D0,0.D0,-16.D0,88.D0,69.D0,-48.D0,-1.D0,
     *21.D0,0.D0,-83.D0,-27.D0,-2.D0,20.D0,17.D0,-23.D0,-7.D0,0.D0,8.D0,
     *-19.D0,5.D0,-23.D0,11.D0,14.D0,-15.D0,-11.D0,0.D0,-21.D0,15.D0,
     *9.D0,-6.D0,-6.D0,9.D0,9.D0,-7.D0,2.D0,0.D0,1.D0,0.D0,3.D0,6.D0,
     *-4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0,39*0.D0/
c
      DATA G90/0.D0,-29775.D0,-1848.D0,-2131.D0,3059.D0,1686.D0,1314.D0,
     *     -2239.D0,  1248.D0,  802.D0,  939.D0, 780.D0, 325.D0,-423.D0,
     *       141.D0,  -214.D0,  353.D0,  245.D0,-109.D0,-165.D0, -36.D0,
     *        61.D0,    65.D0,   59.D0, -178.D0,   3.D0,  18.D0, -96.D0,
     *        77.D0,   -64.D0,    2.D0,   26.D0,  -1.D0,   5.D0,   9.D0,
     *         0.D0,    23.D0,    5.D0,   -1.D0, -10.D0, -12.D0,   3.D0,
     *         4.D0,     2.D0,   -6.D0,    4.D0,   9.D0,   1.D0, -12.D0,
     *         9.D0,    -4.D0,   -2.D0,    7.D0,   1.D0,  -6.D0,  -3.D0,
     *        -4.D0,     2.D0,   -5.D0,   -2.D0,   4.D0,   3.D0,   1.D0,
     *         3.D0,     3.D0,    0.D0,39*0.D0/

      DATA H90/0.D0,  0.D0,5406.D0,   0.D0,-2279.D0,-373.D0,  0.D0,
     *      -284.D0,293.D0,-352.D0,   0.D0,  247.D0,-240.D0, 84.D0,
     *      -299.D0,  0.D0,  46.D0, 154.D0, -153.D0, -69.D0, 97.D0,
     *         0.D0,-16.D0,  82.D0,  69.D0,  -52.D0,   1.D0, 24.D0,
     *         0.D0,-80.D0, -26.D0,   0.D0,   21.D0,  17.D0,-23.D0,
     *        -4.D0,  0.D0,  10.D0, -19.D0,    6.D0, -22.D0, 12.D0,
     *        12.D0,-16.D0, -10.D0,   0.D0,  -20.D0,  15.D0, 11.D0,
     *        -7.D0, -7.D0,   9.D0,   8.D0,   -7.D0,   2.D0,  0.D0,
     *         2.D0,  1.D0,   3.D0,   6.D0,   -4.D0,   0.D0, -2.D0,
     *         3.D0, -1.D0,  -6.D0,39*0.D0/

      DATA G95/0.D0,-29692.D0,-1784.D0,-2200.D0,3070.D0,1681.D0,1335.D0,
     *     -2267.D0,  1249.D0,  759.D0,  940.D0, 780.D0, 290.D0,-418.D0,
     *       122.D0,  -214.D0,  352.D0,  235.D0,-118.D0,-166.D0, -17.D0,
     *        68.D0,    67.D0,   68.D0, -170.D0,  -1.D0,  19.D0, -93.D0,
     *        77.D0,   -72.D0,    1.D0,   28.D0,   5.D0,   4.D0,   8.D0,
     *        -2.D0,    25.D0,    6.D0,   -6.D0,  -9.D0, -14.D0,   9.D0,
     *         6.D0,    -5.D0,   -7.D0,    4.D0,   9.D0,   3.D0, -10.D0,
     *         8.D0,    -8.D0,   -1.D0,   10.D0,  -2.D0,  -8.D0,  -3.D0,
     *        -6.D0,     2.D0,   -4.D0,   -1.D0,   4.D0,   2.D0,   2.D0,
     *         5.D0,     1.D0,    0.D0,  39*0.D0/

      DATA H95/0.D0,  0.D0,5306.D0,  0.D0,-2366.D0,-413.D0,  0.D0,
     *      -262.D0,302.D0,-427.D0,  0.D0,  262.D0,-236.D0, 97.D0,
     *      -306.D0,  0.D0,  46.D0,165.D0, -143.D0, -55.D0,107.D0,
     *         0.D0,-17.D0,  72.D0, 67.D0,  -58.D0,   1.D0, 36.D0,
     *         0.D0,-69.D0, -25.D0,  4.D0,   24.D0,  17.D0,-24.D0,
     *        -6.D0,  0.D0,  11.D0,-21.D0,    8.D0, -23.D0, 15.D0,
     *        11.D0,-16.D0,  -4.D0,  0.D0,  -20.D0,  15.D0, 12.D0,
     *        -6.D0, -8.D0,   8.D0,  5.D0,   -8.D0,   3.D0,  0.D0,
     *         1.D0,  0.D0,   4.D0,  5.D0,   -5.D0,  -1.D0, -2.D0,
     *         1.D0, -2.D0,  -7.D0,39*0.D0/

      DATA G00/0.D0,-29619.4D0,-1728.2D0,-2267.7D0,3068.4D0,1670.9D0,
     *     1339.6D0,  -2288.D0, 1252.1D0,  714.5D0, 932.3D0, 786.8D0,
     *       250.D0,   -403.D0,  111.3D0, -218.8D0, 351.4D0, 222.3D0,
     *     -130.4D0,  -168.6D0,  -12.9D0,   72.3D0,  68.2D0,  74.2D0,
     *     -160.9D0,    -5.9D0,   16.9D0,  -90.4D0,  79.0D0, -74.0D0,
     *         0.D0,    33.3D0,    9.1D0,    6.9D0,   7.3D0,  -1.2D0,
     *       24.4D0,     6.6D0,   -9.2D0,   -7.9D0, -16.6D0,   9.1D0,
     *        7.0D0,    -7.9D0,    -7.D0,     5.D0,   9.4D0,    3.D0,
     *      - 8.4D0,     6.3D0,   -8.9D0,   -1.5D0,   9.3D0,  -4.3D0,
     *       -8.2D0,    -2.6D0,    -6.D0,    1.7D0,  -3.1D0,  -0.5D0,
     *        3.7D0,      1.D0,     2.D0,    4.2D0,   0.3D0,  -1.1D0,
     *        2.7D0,    -1.7D0,   -1.9D0,    1.5D0,  -0.1D0,   0.1D0,
     *       -0.7D0,     0.7D0,    1.7D0,    0.1D0,   1.2D0,   4.0D0,
     *       -2.2D0,    -0.3D0,    0.2D0,    0.9D0,  -0.2D0,   0.9D0,
     *       -0.5D0,     0.3D0,   -0.3D0,   -0.4D0,  -0.1D0,  -0.2D0,
     *       -0.4D0,    -0.2D0,   -0.9D0,    0.3D0,   0.1D0,  -0.4D0,
     *        1.3D0,    -0.4D0,    0.7D0,   -0.4D0,   0.3D0,  -0.1D0,
     *        0.4D0,      0.D0,    0.1D0/

      DATA H00/0.D0,   0.D0,5186.1D0,   0.D0,-2481.6D0,-458.0D0,   0.D0,
     *     -227.6D0,293.4D0,-491.1D0,   0.D0,  272.6D0,-231.9D0,119.8D0,
     *     -303.8D0,   0.D0,  43.8D0,171.9D0, -133.1D0, -39.3D0,106.3D0,
     *         0.D0,-17.4D0,  63.7D0, 65.1D0,  -61.2D0,   0.7D0, 43.8D0,
     *         0.D0,-64.6D0, -24.2D0,  6.2D0,    24.D0,  14.8D0,-25.4D0,
     *       -5.8D0,  0.0D0,  11.9D0,-21.5D0,    8.5D0, -21.5D0, 15.5D0,
     *        8.9D0,-14.9D0,  -2.1D0,  0.0D0,  -19.7D0,  13.4D0, 12.5D0,
     *       -6.2D0, -8.4D0,   8.4D0,  3.8D0,   -8.2D0,   4.8D0,  0.0D0,
     *        1.7D0,  0.0D0,   4.0D0,  4.9D0,   -5.9D0,  -1.2D0, -2.9D0,
     *        0.2D0, -2.2D0,  -7.4D0,  0.0D0,    0.1D0,   1.3D0, -0.9D0,
     *       -2.6D0,  0.9D0,  -0.7D0, -2.8D0,   -0.9D0,  -1.2D0, -1.9D0,
     *       -0.9D0,  0.0D0,  -0.4D0,  0.3D0,    2.5D0,  -2.6D0,  0.7D0,
     *        0.3D0,  0.0D0,   0.0D0,  0.3D0,   -0.9D0,  -0.4D0,  0.8D0,
     *        0.0D0, -0.9D0,   0.2D0,  1.8D0,   -0.4D0,  -1.0D0, -0.1D0,
     *        0.7D0,  0.3D0,   0.6D0,  0.3D0,   -0.2D0,  -0.5D0, -0.9D0/

      DATA G05/0.D0,  -29554.6D0,-1669.0D0,-2337.2D0,3047.7D0,1657.8D0,
     *     1336.3D0,   -2305.8D0, 1246.4D0,  672.5D0, 920.6D0, 798.0D0,
     *      210.7D0,    -379.9D0,  100.0D0, -227.0D0, 354.4D0, 208.9D0,
     *     -136.5D0,    -168.1D0,  -13.6D0,   73.6D0,  69.6D0,  76.7D0,
     *     -151.3D0,     -14.6D0,   14.6D0,  -86.4D0,  79.9D0, -74.5D0,
     *       -1.7D0,      38.7D0,   12.3D0,    9.4D0,   5.4D0,   1.9D0,
     *       24.8D0,       7.6D0,  -11.7D0,   -6.9D0, -18.1D0,  10.2D0,
     *        9.4D0,     -11.3D0,   -4.9D0,    5.6D0,   9.8D0,   3.6D0,
     *       -6.9D0,       5.0D0,  -10.8D0,   -1.3D0,   8.8D0,  -6.7D0,
     *       -9.2D0,      -2.2D0,   -6.1D0,    1.4D0,  -2.4D0,  -0.2D0,
     *        3.1D0,       0.3D0,    2.1D0,    3.8D0,  -0.2D0,  -2.1D0,
     *        2.9D0,      -1.6D0,   -1.9D0,    1.4D0,  -0.3D0,   0.3D0,
     *       -0.8D0,       0.5D0,    1.8D0,    0.2D0,   1.0D0,   4.0D0,
     *       -2.2D0,      -0.3D0,    0.2D0,    0.9D0,  -0.4D0,   1.0D0,
     *       -0.3D0,       0.5D0,   -0.4D0,   -0.4D0,   0.1D0,  -0.5D0,
     *       -0.1D0,      -0.2D0,   -0.9D0,    0.3D0,   0.3D0,  -0.4D0,
     *        1.2D0,      -0.4D0,    0.8D0,   -0.3D0,   0.4D0,  -0.1D0,
     *        0.4D0,      -0.1D0,   -0.2D0/
C
      DATA H05/0.D0,  0.0D0,5078.0D0,  0.0D0,-2594.5D0,-515.4D0,  0.0D0,
     *     -198.9D0,269.7D0,-524.7D0,  0.0D0,  282.1D0,-225.2D0,145.2D0,
     *     -305.4D0,  0.0D0,  42.7D0,180.3D0, -123.5D0, -19.6D0,103.9D0,
     *        0.0D0,-20.3D0,  54.8D0, 63.6D0,  -63.5D0,   0.2D0, 50.9D0,
     *        0.0D0,-61.1D0, -22.6D0,  6.8D0,   25.4D0,  10.9D0,-26.3D0,
     *       -4.6D0,  0.0D0,  11.2D0,-20.9D0,    9.8D0, -19.7D0, 16.2D0,
     *        7.6D0,-12.8D0,  -0.1D0,  0.0D0,  -20.1D0,  12.7D0, 12.7D0,
     *       -6.7D0, -8.2D0,   8.1D0,  2.9D0,   -7.7D0,   6.0D0,  0.0D0,
     *        2.2D0,  0.1D0,   4.5D0,  4.8D0,   -6.7D0,  -1.0D0, -3.5D0,
     *       -0.9D0, -2.3D0,  -7.9D0,  0.0D0,    0.3D0,   1.4D0, -0.8D0,
     *       -2.3D0,  0.9D0,  -0.6D0, -2.7D0,   -1.1D0,  -1.6D0, -1.9D0,
     *       -1.4D0,  0.0D0,  -0.6D0,  0.2D0,    2.4D0,  -2.6D0,  0.6D0,
     *        0.4D0,  0.0D0,   0.0D0,  0.3D0,   -0.9D0,  -0.3D0,  0.9D0,
     *        0.0D0, -0.8D0,   0.3D0,  1.7D0,   -0.5D0,  -1.1D0,  0.0D0,
     *        0.6D0,  0.2D0,   0.5D0,  0.4D0,   -0.2D0,  -0.6D0, -0.9D0/
C
      DATA G10/0.00D0,-29496.57D0,-1586.42D0,-2396.06D0,3026.34D0,
     *      1668.17D0,  1339.85D0,-2326.54D0, 1232.10D0, 633.73D0, 
     *       912.66D0,   808.97D0,  166.58D0, -356.83D0,  89.40D0,
     *      -230.87D0,   357.29D0,  200.26D0, -141.05D0,-163.17D0,
     *        -8.03D0,    72.78D0,   68.69D0,   75.92D0,-141.40D0, 
     *       -22.83D0,    13.10D0,  -78.09D0,   80.44D0, -75.00D0,
     *        -4.55D0,    45.24D0,   14.00D0,   10.46D0,   1.64D0,
     *         4.92D0,    24.41D0,    8.21D0,  -14.50D0,  -5.59D0,
     *       -19.34D0,    11.61D0,   10.85D0,  -14.05D0,  -3.54D0,
     *         5.50D0,     9.45D0,    3.45D0,   -5.27D0,   3.13D0,
     *       -12.38D0,    -0.76D0,    8.43D0,   -8.42D0, -10.08D0,
     *        -1.94D0,    -6.24D0,    0.89D0,   -1.07D0,  -0.16D0,
     *         2.45D0,    -0.33D0,    2.13D0,    3.09D0,  -1.03D0,
     *        -2.80D0,     3.05D0,   -1.48D0,   -2.03D0,   1.65D0,
     *        -0.51D0,     0.54D0,   -0.79D0,    0.37D0,   1.79D0,
     *         0.12D0,     0.75D0,    3.75D0,   -2.12D0,  -0.21D0,
     *         0.30D0,     1.04D0,   -0.63D0,    0.95D0,  -0.11D0,
     *         0.52D0,    -0.39D0,   -0.37D0,    0.21D0,  -0.77D0,
     *         0.04D0,    -0.09D0,   -0.89D0,    0.31D0,   0.42D0,
     *        -0.45D0,     1.08D0,   -0.31D0,    0.78D0,  -0.18D0,
     *         0.38D0,     0.02D0,    0.42D0,   -0.26D0,  -0.26D0/
C
      DATA H10/0.00D0,  0.00D0,4944.26D0,   0.00D0,-2708.54D0,
     *     -575.73D0,   0.00D0,-160.40D0, 251.75D0, -537.03D0, 0.00D0,
     *      286.48D0,-211.03D0, 164.46D0,-309.72D0,    0.00D0,44.58D0,
     *      189.01D0,-118.06D0,  -0.01D0, 101.04D0,    0.00D0,-20.90D0,
     *       44.18D0,  61.54D0, -66.26D0,   3.02D0,   55.40D0,  0.00D0,
     *      -57.80D0, -21.20D0,   6.54D0,  24.96D0,    7.03D0,-27.61D0,
     *       -3.28D0,   0.00D0,  10.84D0, -20.03D0,   11.83D0,-17.41D0, 
     *       16.71D0,   6.96D0, -10.74D0,   1.64D0,    0.00D0,-20.54D0,
     *       11.51D0,  12.75D0,  -7.14D0,  -7.42D0,    7.97D0,  2.14D0,
     *       -6.08D0,   7.01D0,   0.00D0,   2.73D0,   -0.10D0,  4.71D0,
     *        4.44D0,  -7.22D0,  -0.96D0,  -3.95D0,   -1.99D0, -1.97D0,
     *       -8.31D0,   0.00D0,   0.13D0,   1.67D0,   -0.66D0, -1.76D0,
     *        0.85D0,  -0.39D0,  -2.51D0,  -1.27D0,   -2.11D0, -1.94D0,
     *       -1.86D0,   0.00D0,  -0.87D0,   0.27D0,    2.13D0, -2.49D0,
     *        0.49D0,   0.59D0,   0.00D0,   0.13D0,    0.27D0, -0.86D0,
     *       -0.23D0,   0.87D0,   0.00D0,  -0.87D0,    0.30D0,  1.66D0,
     *       -0.59D0,  -1.14D0,  -0.07D0,   0.54D0,    0.10D0,  0.49D0,
     *        0.44D0,  -0.25D0,  -0.53D0,  -0.79D0/
C
      DATA G15/0.D0,-29442.0D0, -1501.0D0, -2445.1D0,  3012.9D0,
     *     1676.7D0,  1350.7D0, -2352.3D0,  1225.6D0,   582.0D0,907.6D0, 
     *      813.7D0,   120.4D0,  -334.9D0,    70.4D0,  -232.6D0,360.1D0,
     *      192.4D0,  -140.9D0,  -157.5D0,     4.1D0,    70.0D0, 67.7D0,
     *       72.7D0,  -129.9D0,   -28.9D0,    13.2D0,   -70.9D0, 81.6D0,
     *      -76.1D0,    -6.8D0,    51.8D0,    15.0D0,     9.4D0, -2.8D0,
     *        6.8D0,    24.2D0,     8.8D0,   -16.9D0,    -3.2D0,-20.6D0,
     *       13.4D0,    11.7D0,   -15.9D0,    -2.0D0,     5.4D0,  8.8D0, 
     *        3.1D0,    -3.3D0,     0.7D0,   -13.3D0,    -0.1D0,  8.7D0,
     *       -9.1D0,   -10.5D0,    -1.9D0,    -6.3D0,     0.1D0,  0.5D0,
     *       -0.5D0,     1.8D0,    -0.7D0,     2.1D0,     2.4D0, -1.8D0,
     *       -3.6D0,     3.1D0,    -1.5D0,    -2.3D0,     2.0D0, -0.8D0,
     *        0.6D0,    -0.7D0,     0.2D0,     1.7D0,    -0.2D0,  0.4D0,
     *        3.5D0,    -1.9D0,    -0.2D0,     0.4D0,     1.2D0, -0.8D0,
     *        0.9D0,     0.1D0,     0.5D0,    -0.3D0,    -0.4D0,  0.2D0,
     *       -0.9D0,     0.0D0,     0.0D0,    -0.9D0,     0.4D0,  0.5D0,
     *       -0.5D0,     1.0D0,    -0.2D0,     0.8D0,    -0.1D0,  0.3D0,
     *        0.1D0,     0.5D0,    -0.4D0,    -0.3D0/
c
      DATA H15/0.D0,     0.0D0, 4797.1D0,     0.0D0, -2845.6D0,-641.9D0,
     *        0.0D0,  -115.3D0,  244.9D0,  -538.4D0,     0.0D0, 283.3D0,
     *     -188.7D0,   180.9D0, -329.5D0,     0.0D0,    47.3D0, 197.0D0,
     *     -119.3D0,    16.0D0,  100.2D0,     0.0D0,   -20.8D0,  33.2D0,
     *       58.9D0,   -66.7D0,    7.3D0,    62.6D0,     0.0D0, -54.1D0,
     *      -19.5D0,     5.7D0,   24.4D0,     3.4D0,   -27.4D0,  -2.2D0,
     *        0.0D0,    10.1D0,  -18.3D0,    13.3D0,   -14.6D0,  16.2D0,
     *        5.7D0,    -9.1D0,    2.1D0,     0.0D0,   -21.6D0,  10.8D0,
     *       11.8D0,    -6.8D0,   -6.9D0,     7.8D0,     1.0D0,  -4.0D0,
     *        8.4D0,     0.0D0,    3.2D0,    -0.4D0,     4.6D0,   4.4D0,
     *       -7.9D0,    -0.6D0,   -4.2D0,    -2.8D0,    -1.2D0,  -8.7D0,
     *        0.0D0,    -0.1D0,    2.0D0,    -0.7D0,    -1.1D0,   0.8D0,
     *       -0.2D0,    -2.2D0,   -1.4D0,    -2.5D0,    -2.0D0,  -2.4D0,
     *        0.0D0,    -1.1D0,    0.4D0,     1.9D0,    -2.2D0,   0.3D0,
     *        0.7D0,    -0.1D0,    0.3D0,     0.2D0,    -0.9D0,  -0.1D0,
     *        0.7D0,     0.0D0,   -0.9D0,     0.4D0,     1.6D0,  -0.5D0,
     *       -1.2D0,    -0.1D0,    0.4D0,    -0.1D0,     0.4D0,   0.5D0,
     *       -0.3D0,    -0.4D0,   -0.8D0/
c
      DATA DG15/0.0D0, 10.3D0,  18.1D0,  -8.7D0,  -3.3D0,  2.1D0, 3.4D0,
     *         -5.5D0, -0.7D0, -10.1D0,  -0.7D0,   0.2D0, -9.1D0, 4.1D0,
     *         -4.3D0, -0.2D0,   0.5D0,  -1.3D0,  -0.1D0,  1.4D0, 3.9D0,
     *         -0.3D0, -0.1D0,  -0.7D0,   2.1D0,  -1.2D0,  0.3D0, 1.6D0,
     *          0.3D0, -0.2D0,  -0.5D0,   1.3D0,   0.1D0, -0.6D0,-0.8D0,
     *          0.2D0,  0.2D0,   0.0D0,  -0.6D0,   0.5D0, -0.2D0, 0.4D0,
     *          0.1D0, -0.4D0,   0.3D0/
c
      DATA DH15/0.0D0,  0.0D0, -26.6D0,   0.0D0, -27.4D0,-14.1D0, 0.0D0,
     *          8.2D0, -0.4D0,   1.8D0,   0.0D0,  -1.3D0,  5.3D0, 2.9D0,
     *         -5.2D0,  0.0D0,   0.6D0,   1.7D0,  -1.2D0,  3.4D0, 0.0D0,
     *          0.0D0,  0.0D0,  -2.1D0,  -0.7D0,   0.2D0,  0.9D0, 1.0D0,
     *          0.0D0,  0.8D0,   0.4D0,  -0.2D0,  -0.3D0, -0.6D0, 0.1D0,
     *         -0.2D0,  0.0D0,  -0.3D0,   0.3D0,   0.1D0,  0.5D0,-0.2D0,
     *         -0.3D0,  0.3D0,   0.0D0/
C
C
c      IF (VGSEY.EQ.0..AND.VGSEZ.EQ.0..AND.ISW.NE.1) THEN
c      PRINT *, ''
c      PRINT *,
c     *' RECALC_08: RADIAL SOLAR WIND --> GSW SYSTEM IDENTICAL HERE'
c      PRINT *,
c     *' TO STANDARD GSM (I.E., XGSW AXIS COINCIDES WITH EARTH-SUN LINE)'
c      PRINT *, ''
c      ISW=1
c      ENDIF

c      IF ((VGSEY.NE.0.D0.OR.VGSEZ.NE.0.D0).AND.ISW.NE.2) THEN       !  CORRECTED DEC.01, 2010
c      PRINT *, ''
c      PRINT *,
c     *' WARNING: NON-RADIAL SOLAR WIND FLOW SPECIFIED IN RECALC_08;'
c      PRINT *,
c     *' HENCE XGSW AXIS IS ASSUMED ORIENTED ANTIPARALLEL TO V_SW VECTOR'
c      PRINT *, ''
c      ISW=2
c      ENDIF
C
      IY=IYEAR
C
C  WE ARE RESTRICTED BY THE INTERVAL 1965-2020, FOR WHICH EITHER THE IGRF/DGRF COEFFICIENTS OR SECULAR VELOCITIES
c    ARE KNOWN; IF IYEAR IS OUTSIDE THIS INTERVAL, THEN THE SUBROUTINE USES THE
C      NEAREST LIMITING VALUE AND PRINTS A WARNING:
C
      IF(IY.LT.1965) THEN
       IY=1965
       WRITE (*,10) IYEAR,IY
      ENDIF

      IF(IY.GT.2020) THEN
       IY=2020
       WRITE (*,10) IYEAR,IY
      ENDIF
C
C  CALCULATE THE ARRAY REC, CONTAINING COEFFICIENTS FOR THE RECURSION RELATIONS,
C  USED IN THE IGRF SUBROUTINE FOR CALCULATING THE ASSOCIATE LEGENDRE POLYNOMIALS
C  AND THEIR DERIVATIVES:
c
      DO 20 N=1,14
         N2=2*N-1
         N2=N2*(N2-2)
         DO 20 M=1,N
            MN=N*(N-1)/2+M
20    REC(MN)=DFLOAT((N-M)*(N+M-2))/DFLOAT(N2)
C
      IF (IY.LT.1970) GOTO 50          !INTERPOLATE BETWEEN 1965 - 1970
      IF (IY.LT.1975) GOTO 60          !INTERPOLATE BETWEEN 1970 - 1975
      IF (IY.LT.1980) GOTO 70          !INTERPOLATE BETWEEN 1975 - 1980
      IF (IY.LT.1985) GOTO 80          !INTERPOLATE BETWEEN 1980 - 1985
      IF (IY.LT.1990) GOTO 90          !INTERPOLATE BETWEEN 1985 - 1990
      IF (IY.LT.1995) GOTO 100         !INTERPOLATE BETWEEN 1990 - 1995
      IF (IY.LT.2000) GOTO 110         !INTERPOLATE BETWEEN 1995 - 2000
      IF (IY.LT.2005) GOTO 120         !INTERPOLATE BETWEEN 2000 - 2005
      IF (IY.LT.2010) GOTO 130         !INTERPOLATE BETWEEN 2005 - 2010
      IF (IY.LT.2015) GOTO 140         !INTERPOLATE BETWEEN 2010 - 2015
C
C       EXTRAPOLATE BEYOND 2015:
C
      DT=DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-2015.D0
      DO 40 N=1,105
         G(N)=G15(N)
         H(N)=H15(N)
         IF (N.GT.45) GOTO 40
         G(N)=G(N)+DG15(N)*DT
         H(N)=H(N)+DH15(N)*DT
40    CONTINUE
      GOTO 300
C
C       INTERPOLATE BETWEEEN 1965 - 1970:
C
50    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1965)/5.D0
      F1=1.D0-F2
      DO 55 N=1,105
         G(N)=G65(N)*F1+G70(N)*F2
55       H(N)=H65(N)*F1+H70(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1970 - 1975:
C
60    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1970)/5.D0
      F1=1.D0-F2
      DO 65 N=1,105
         G(N)=G70(N)*F1+G75(N)*F2
65       H(N)=H70(N)*F1+H75(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1975 - 1980:
C
70    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1975)/5.D0
      F1=1.D0-F2
      DO 75 N=1,105
         G(N)=G75(N)*F1+G80(N)*F2
75       H(N)=H75(N)*F1+H80(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1980 - 1985:
C
80    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1980)/5.D0
      F1=1.D0-F2
      DO 85 N=1,105
         G(N)=G80(N)*F1+G85(N)*F2
85       H(N)=H80(N)*F1+H85(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1985 - 1990:
C
90    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1985)/5.D0
      F1=1.D0-F2
      DO 95 N=1,105
         G(N)=G85(N)*F1+G90(N)*F2
95       H(N)=H85(N)*F1+H90(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1990 - 1995:
C
100   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1990)/5.D0
      F1=1.D0-F2
      DO 105 N=1,105
         G(N)=G90(N)*F1+G95(N)*F2
105      H(N)=H90(N)*F1+H95(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1995 - 2000:
C
110   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1995)/5.D0
      F1=1.D0-F2
      DO 115 N=1,105   !  THE 2000 COEFFICIENTS (G00) GO THROUGH THE ORDER 13, NOT 10
         G(N)=G95(N)*F1+G00(N)*F2
115      H(N)=H95(N)*F1+H00(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 2000 - 2005:
C
120   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-2000)/5.D0
      F1=1.D0-F2
      DO 125 N=1,105
         G(N)=G00(N)*F1+G05(N)*F2
125      H(N)=H00(N)*F1+H05(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 2005 - 2010:
C
130   F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-2005)/5.
      F1=1.-F2
      DO 135 N=1,105
         G(N)=G05(N)*F1+G10(N)*F2
135      H(N)=H05(N)*F1+H10(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 2010 - 2015:
C
140   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25-2010)/5.
      F1=1.-F2
      DO 145 N=1,105
         G(N)=G10(N)*F1+G15(N)*F2
145      H(N)=H10(N)*F1+H15(N)*F2
      GOTO 300
C
C   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
C   THEM BY SCHMIDT NORMALIZATION FACTORS:
C
300   S=1.D0
      DO 150 N=2,14
         MN=N*(N-1)/2+1
         S=S*DFLOAT(2*N-3)/DFLOAT(N-1)
         G(MN)=G(MN)*S
         H(MN)=H(MN)*S
         P=S
         DO 150 M=2,N
            AA=1.D0
            IF (M.EQ.2) AA=2.D0
            P=P*DSQRT(AA*DFLOAT(N-M+1)/DFLOAT(N+M-2))
            MNN=MN+M-1
            G(MNN)=G(MNN)*P
150         H(MNN)=H(MNN)*P

           G_10=-G(2)
           G_11= G(3)
           H_11= H(3)
C
C  NOW CALCULATE GEO COMPONENTS OF THE UNIT VECTOR EzMAG, PARALLEL TO GEODIPOLE AXIS:
C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
C         ST0 * CL0                ST0 * SL0                CT0
C
      SQ=G_11**2+H_11**2
      SQQ=DSQRT(SQ)
      SQR=DSQRT(G_10**2+SQ)
      SL0=-H_11/SQQ
      CL0=-G_11/SQQ
      ST0=SQQ/SQR
      CT0=G_10/SQR
      STCL=ST0*CL0
      STSL=ST0*SL0
      CTSL=CT0*SL0
      CTCL=CT0*CL0
C
C  NOW CALCULATE GEI COMPONENTS (S1,S2,S3) OF THE UNIT VECTOR S = EX_GSE
C    POINTING FROM THE EARTH'S CENTER TO SUN
C
      CALL SUN_08 (IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
      S1=DCOS(SRASN)*DCOS(SDEC)
      S2=DSIN(SRASN)*DCOS(SDEC)
      S3=DSIN(SDEC)
C
C  NOW CALCULATE GEI COMPONENTS (DZ1,DZ2,DZ3) OF THE UNIT VECTOR EZGSE
C  POINTING NORTHWARD AND ORTHOGONAL TO THE ECLIPTIC PLANE, AS
C  (0,-SIN(OBLIQ),COS(OBLIQ)). FOR THE EPOCH 1978, OBLIQ = 23.44214 DEGS.
C  HERE WE USE A MORE ACCURATE TIME-DEPENDENT VALUE, DETERMINED AS:
C
      DJ=DFLOAT(365*(IY-1900)+(IY-1901)/4 +IDAY)
     * -0.5+DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      T=DJ/36525.D0
      OBLIQ=(23.45229D0-0.0130125D0*T)/57.2957795D0
      DZ1=0.D0
      DZ2=-DSIN(OBLIQ)
      DZ3=DCOS(OBLIQ)
C
C  NOW WE OBTAIN GEI COMPONENTS OF THE UNIT VECTOR EYGSE=(DY1,DY2,DY3),
C  COMPLETING THE RIGHT-HANDED SYSTEM. THEY CAN BE FOUND FROM THE VECTOR
C  PRODUCT EZGSE x EXGSE = (DZ1,DZ2,DZ3) x (S1,S2,S3):
C
      DY1=DZ2*S3-DZ3*S2
      DY2=DZ3*S1-DZ1*S3
      DY3=DZ1*S2-DZ2*S1
C
C  NOW LET'S CALCULATE GEI COMPONENTS OF THE UNIT VECTOR X = EXGSW, DIRECTED ANTIPARALLEL
C  TO THE OBSERVED SOLAR WIND FLOW. FIRST, CALCULATE ITS COMPONENTS IN GSE:
C
      V=DSQRT(VGSEX**2+VGSEY**2+VGSEZ**2)
      DX1=-VGSEX/V
      DX2=-VGSEY/V
      DX3=-VGSEZ/V
C
C  THEN IN GEI:
C
      X1=DX1*S1+DX2*DY1+DX3*DZ1
      X2=DX1*S2+DX2*DY2+DX3*DZ2
      X3=DX1*S3+DX2*DY3+DX3*DZ3
C
C  NOW CALCULATE GEI COMPONENTS (DIP1,DIP2,DIP3) OF THE UNIT VECTOR DIP = EZ_SM = EZ_MAG,
C   ALIGNED WITH THE GEODIPOLE AND POINTING NORTHWARD FROM ECLIPTIC PLANE:
C
      CGST=DCOS(GST)
      SGST=DSIN(GST)
C
      DIP1=STCL*CGST-STSL*SGST
      DIP2=STCL*SGST+STSL*CGST
      DIP3=CT0
C
C  THIS ALLOWS US TO CALCULATE GEI COMPONENTS OF THE UNIT VECTOR Y = EYGSW
C   BY TAKING THE VECTOR PRODUCT DIP x X AND NORMALIZING IT TO UNIT LENGTH:
C
      Y1=DIP2*X3-DIP3*X2
      Y2=DIP3*X1-DIP1*X3
      Y3=DIP1*X2-DIP2*X1
      Y=DSQRT(Y1*Y1+Y2*Y2+Y3*Y3)
      Y1=Y1/Y
      Y2=Y2/Y
      Y3=Y3/Y
C
C   AND GEI COMPONENTS OF THE UNIT VECTOR Z = EZGSW = EXGSW x EYGSW = X x Y:
C
      Z1=X2*Y3-X3*Y2
      Z2=X3*Y1-X1*Y3
      Z3=X1*Y2-X2*Y1
C
C   ELEMENTS OF THE MATRIX GSE TO GSW ARE THE SCALAR PRODUCTS:
C
C  E11=(EXGSE,EXGSW)  E12=(EXGSE,EYGSW)  E13=(EXGSE,EZGSW)
C  E21=(EYGSE,EXGSW)  E22=(EYGSE,EYGSW)  E23=(EYGSE,EZGSW)
C  E31=(EZGSE,EXGSW)  E32=(EZGSE,EYGSW)  E33=(EZGSE,EZGSW)
C
      E11= S1*X1 +S2*X2 +S3*X3
      E12= S1*Y1 +S2*Y2 +S3*Y3
      E13= S1*Z1 +S2*Z2 +S3*Z3
      E21=DY1*X1+DY2*X2+DY3*X3
      E22=DY1*Y1+DY2*Y2+DY3*Y3
      E23=DY1*Z1+DY2*Z2+DY3*Z3
      E31=DZ1*X1+DZ2*X2+DZ3*X3
      E32=DZ1*Y1+DZ2*Y2+DZ3*Y3
      E33=DZ1*Z1+DZ2*Z2+DZ3*Z3
C
C   GEODIPOLE TILT ANGLE IN THE GSW SYSTEM: PSI=ARCSIN(DIP,EXGSW)
C
      SPS=DIP1*X1+DIP2*X2+DIP3*X3
      CPS=DSQRT(1.D0-SPS**2)
      PSI=DASIN(SPS)
C
C   ELEMENTS OF THE MATRIX GEO TO GSW ARE THE SCALAR PRODUCTS:
C
C   A11=(EXGEO,EXGSW), A12=(EYGEO,EXGSW), A13=(EZGEO,EXGSW),
C   A21=(EXGEO,EYGSW), A22=(EYGEO,EYGSW), A23=(EZGEO,EYGSW),
C   A31=(EXGEO,EZGSW), A32=(EYGEO,EZGSW), A33=(EZGEO,EZGSW),
C
C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
C
C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
C  EXGSW=(X1,X2,X3),  EYGSW=(Y1,Y2,Y3),   EZGSW=(Z1,Z2,Z3)
C                                                           AND  THEREFORE:
C
      A11=X1*CGST+X2*SGST
      A12=-X1*SGST+X2*CGST
      A13=X3
      A21=Y1*CGST+Y2*SGST
      A22=-Y1*SGST+Y2*CGST
      A23=Y3
      A31=Z1*CGST+Z2*SGST
      A32=-Z1*SGST+Z2*CGST
      A33=Z3
C
C  NOW CALCULATE ELEMENTS OF THE MATRIX MAG TO SM (ONE ROTATION ABOUT THE GEODIPOLE AXIS);
C   THEY ARE FOUND AS THE SCALAR PRODUCTS: CFI=GM22=(EYSM,EYMAG)=(EYGSW,EYMAG),
C                                          SFI=GM23=(EYSM,EXMAG)=(EYGSW,EXMAG),
C    DERIVED AS FOLLOWS:
C
C IN GEO, THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THEIR COMPONENTS ARE:
C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
C            -ST0
C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
C            -SL0*SIN(GST)+CL0*COS(GST)
C             0
C  NOW, NOTE THAT GEI COMPONENTS OF EYSM=EYGSW WERE FOUND ABOVE AS Y1, Y2, AND Y3,
C  AND WE ONLY HAVE TO COMBINE THESE QUANTITIES INTO SCALAR PRODUCTS:
C
      EXMAGX=CT0*(CL0*CGST-SL0*SGST)
      EXMAGY=CT0*(CL0*SGST+SL0*CGST)
      EXMAGZ=-ST0
      EYMAGX=-(SL0*CGST+CL0*SGST)
      EYMAGY=-(SL0*SGST-CL0*CGST)
      CFI=Y1*EYMAGX+Y2*EYMAGY
      SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ
C
 10   FORMAT(//1X,
     *'**** RECALC_08 WARNS: YEAR IS OUT OF INTERVAL 1965-2020: IYEAR=',
     *I4,/,6X,'CALCULATIONS WILL BE DONE FOR IYEAR=',I4,/)
      RETURN
      END
c
c==================================================================================

      SUBROUTINE GSWGSE_08 (XGSW,YGSW,ZGSW,XGSE,YGSE,ZGSE,J)
C
C  THIS SUBROUTINE TRANSFORMS COMPONENTS OF ANY VECTOR BETWEEN THE STANDARD GSE
C  COORDINATE SYSTEM AND THE GEOCENTRIC SOLAR-WIND (GSW, aka GSWM), DEFINED AS FOLLOWS
C  (HONES ET AL., PLANET.SPACE SCI., V.34, P.889, 1986; TSYGANENKO ET AL., JGRA,
C  V.103(A4), P.6827, 1998):
C
C  IN THE GSW SYSTEM, X AXIS IS ANTIPARALLEL TO THE OBSERVED DIRECTION OF THE SOLAR WIND FLOW.
C  TWO OTHER AXES, Y AND Z, ARE DEFINED IN THE SAME WAY AS FOR THE STANDARD GSM, THAT IS,
C  Z AXIS ORTHOGONAL TO X AXIS, POINTS NORTHWARD, AND LIES IN THE PLANE DEFINED BY THE X-
C  AND GEODIPOLE AXIS. THE Y AXIS COMPLETES THE RIGHT-HANDED SYSTEM.
C
C  THE GSW SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM IN THE CASE OF
C   A STRICTLY RADIAL SOLAR WIND FLOW.
C
C  AUTHOR:  N. A. TSYGANENKO
C  ADDED TO 2008 VERSION OF GEOPACK: JAN 27, 2008.
C  LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C                    J>0                       J<0
C-----INPUT:   J,XGSW,YGSW,ZGSW          J,XGSE,YGSE,ZGSE
C-----OUTPUT:    XGSE,YGSE,ZGSE            XGSW,YGSW,ZGSW
C
C  IMPORTANT THINGS TO REMEMBER:
C
C   (1) BEFORE CALLING GSWGSE_08, BE SURE TO INVOKE SUBROUTINE RECALC_08, IN ORDER
C       TO DEFINE ALL NECESSARY ELEMENTS OF TRANSFORMATION MATRICES
C
C   (2) IN THE ABSENCE OF INFORMATION ON THE SOLAR WIND DIRECTION, E.G., WITH ONLY SCALAR
C       SPEED V KNOWN, THIS SUBROUTINE CAN BE USED TO CONVERT VECTORS TO ABERRATED
C       COORDINATE SYSTEM, TAKING INTO ACCOUNT EARTH'S ORBITAL SPEED OF 29 KM/S.
C       TO DO THAT, SPECIFY THE LAST 3 PARAMETERS IN RECALC_08 AS FOLLOWS:
C       VGSEX=-V, VGSEY=29.0, VGSEZ=0.0.
C
C       IT SHOULD ALSO BE KEPT IN MIND THAT IN SOME SOLAR WIND DATABASES THE ABERRATION
C       EFFECT HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29 KM/S FROM VYGSE;
C       IN THAT CASE, THE ORIGINAL VYGSE VALUES SHOULD BE RESTORED BY ADDING BACK THE
C       29 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST BE VERIFIED WITH THE DATA
C       ORIGINATOR (OR CAN BE DETERMINED BY CALCULATING THE AVERAGE VGSEY OVER
C       A SUFFICIENTLY LONG TIME INTERVAL)
C
C   (3) IF NO INFORMATION IS AVAILABLE ON THE SOLAR WIND SPEED, THEN SET VGSEX=-400.0
c       AND  VGSEY=VGSEZ=0. IN THAT CASE, THE GSW COORDINATE SYSTEM BECOMES
c       IDENTICAL TO THE STANDARD ONE.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ AAA(25),E11,E21,E31,E12,E22,E32,E13,E23,E33
C
C  DIRECT TRANSFORMATION:
C
      IF (J.GT.0) THEN
        XGSE=XGSW*E11+YGSW*E12+ZGSW*E13
        YGSE=XGSW*E21+YGSW*E22+ZGSW*E23
        ZGSE=XGSW*E31+YGSW*E32+ZGSW*E33
      ENDIF
C
C   INVERSE TRANSFORMATION: CARRIED OUT USING THE TRANSPOSED MATRIX:
C
      IF (J.LT.0) THEN
        XGSW=XGSE*E11+YGSE*E21+ZGSE*E31
        YGSW=XGSE*E12+YGSE*E22+ZGSE*E32
        ZGSW=XGSE*E13+YGSE*E23+ZGSE*E33
      ENDIF
C
      RETURN
      END
C
C========================================================================================
C
      SUBROUTINE GEOMAG_08 (XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J)
C
C    CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICE VERSA.
C
C                    J>0                       J<0
C-----INPUT:  J,XGEO,YGEO,ZGEO           J,XMAG,YMAG,ZMAG
C-----OUTPUT:    XMAG,YMAG,ZMAG           XGEO,YGEO,ZGEO
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEOMAG_08 IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR AND/OR IDAY HAVE BEEN CHANGED
C
C  NO INFORMATION IS REQUIRED HERE ON THE SOLAR WIND VELOCITY, SO ONE
C  CAN SET VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0 IN RECALC_08.
C
C   LAST MOFIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C   AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB(26)

      IF(J.GT.0) THEN
       XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
       YMAG=YGEO*CL0-XGEO*SL0
       ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
      ELSE
       XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
       YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
       ZGEO=ZMAG*CT0-XMAG*ST0
      ENDIF

      RETURN
      END
c
c=========================================================================================
c
      SUBROUTINE GEIGEO_08 (XGEI,YGEI,ZGEI,XGEO,YGEO,ZGEO,J)
C
C   CONVERTS EQUATORIAL INERTIAL (GEI) TO GEOGRAPHICAL (GEO) COORDS
C   OR VICE VERSA.
C                    J>0                J<0
C----INPUT:  J,XGEI,YGEI,ZGEI    J,XGEO,YGEO,ZGEO
C----OUTPUT:   XGEO,YGEO,ZGEO      XGEI,YGEI,ZGEI
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEIGEO_08 IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C  NO INFORMATION IS REQUIRED HERE ON THE SOLAR WIND VELOCITY, SO ONE
C  CAN SET VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0 IN RECALC_08.
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)

C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ A(13),CGST,SGST,B(19)
C
      IF(J.GT.0) THEN
       XGEO=XGEI*CGST+YGEI*SGST
       YGEO=YGEI*CGST-XGEI*SGST
       ZGEO=ZGEI
      ELSE
       XGEI=XGEO*CGST-YGEO*SGST
       YGEI=YGEO*CGST+XGEO*SGST
       ZGEI=ZGEO
      ENDIF

      RETURN
      END
C
C=======================================================================================
C
      SUBROUTINE MAGSM_08 (XMAG,YMAG,ZMAG,XSM,YSM,ZSM,J)
C
C  CONVERTS DIPOLE (MAG) TO SOLAR MAGNETIC (SM) COORDINATES OR VICE VERSA
C
C                    J>0              J<0
C----INPUT: J,XMAG,YMAG,ZMAG     J,XSM,YSM,ZSM
C----OUTPUT:    XSM,YSM,ZSM       XMAG,YMAG,ZMAG
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE MAGSM_08 IN THREE CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES, OR
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE CHANGED, AND/OR
C     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED
C
C    IMPORTANT NOTE:
C
C        A NON-STANDARD DEFINITION IS IMPLIED HERE FOR THE SOLAR MAGNETIC COORDINATE
C        SYSTEM:  IT IS ASSUMED THAT THE XSM AXIS LIES IN THE PLANE DEFINED BY THE
C        GEODIPOLE AXIS AND THE OBSERVED VECTOR OF THE SOLAR WIND FLOW (RATHER THAN
C        THE EARTH-SUN LINE).  IN ORDER TO CONVERT MAG COORDINATES TO AND FROM THE
C        STANDARD SM COORDINATES, INVOKE RECALC_08 WITH VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ A(8),SFI,CFI,B(24)
C
      IF (J.GT.0) THEN
       XSM=XMAG*CFI-YMAG*SFI
       YSM=XMAG*SFI+YMAG*CFI
       ZSM=ZMAG
      ELSE
       XMAG=XSM*CFI+YSM*SFI
       YMAG=YSM*CFI-XSM*SFI
       ZMAG=ZSM
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
       SUBROUTINE SMGSW_08 (XSM,YSM,ZSM,XGSW,YGSW,ZGSW,J)
C
C  CONVERTS SOLAR MAGNETIC (SM) TO GEOCENTRIC SOLAR-WIND (GSW) COORDINATES OR VICE VERSA.
C
C                  J>0                 J<0
C-----INPUT: J,XSM,YSM,ZSM        J,XGSW,YGSW,ZGSW
C----OUTPUT:  XGSW,YGSW,ZGSW       XSM,YSM,ZSM
C
C  ATTENTION:  SUBROUTINE RECALC_08 MUST BE INVOKED BEFORE SMGSW_08 IN THREE CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED
C
C    IMPORTANT NOTE:
C
C        A NON-STANDARD DEFINITION IS IMPLIED HERE FOR THE SOLAR MAGNETIC (SM) COORDINATE
C        SYSTEM:  IT IS ASSUMED THAT THE XSM AXIS LIES IN THE PLANE DEFINED BY THE
C        GEODIPOLE AXIS AND THE OBSERVED VECTOR OF THE SOLAR WIND FLOW (RATHER THAN
C        THE EARTH-SUN LINE).  IN ORDER TO CONVERT MAG COORDINATES TO AND FROM THE
C        STANDARD SM COORDINATES, INVOKE RECALC_08 WITH VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ A(10),SPS,CPS,B(22)

      IF (J.GT.0) THEN
       XGSW=XSM*CPS+ZSM*SPS
       YGSW=YSM
       ZGSW=ZSM*CPS-XSM*SPS
      ELSE
       XSM=XGSW*CPS-ZGSW*SPS
       YSM=YGSW
       ZSM=XGSW*SPS+ZGSW*CPS
      ENDIF

      RETURN
      END
C
C==========================================================================================
C
      SUBROUTINE GEOGSW_08 (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW,J)
C
C CONVERTS GEOGRAPHIC (GEO) TO GEOCENTRIC SOLAR-WIND (GSW) COORDINATES OR VICE VERSA.
C
C                   J>0                   J<0
C----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSW,YGSW,ZGSW
C---- OUTPUT:    XGSW,YGSW,ZGSW      XGEO,YGEO,ZGEO
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEOGSW_08 IN THREE CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES, OR
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE CHANGED, AND/OR
C     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED
C
C  NOTE: THIS SUBROUTINE CONVERTS GEO VECTORS TO AND FROM THE SOLAR-WIND GSW COORDINATE
C        SYSTEM, TAKING INTO ACCOUNT POSSIBLE DEFLECTIONS OF THE SOLAR WIND DIRECTION FROM
C        STRICTLY RADIAL.  BEFORE CONVERTING TO/FROM STANDARD GSM COORDINATES, INVOKE RECALC_08
C        WITH VGSEX=-400.0 and VGSEY=0.0, VGSEZ=0.0
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ AA(16),A11,A21,A31,A12,A22,A32,A13,A23,A33,B(9)
C
      IF (J.GT.0) THEN
       XGSW=A11*XGEO+A12*YGEO+A13*ZGEO
       YGSW=A21*XGEO+A22*YGEO+A23*ZGEO
       ZGSW=A31*XGEO+A32*YGEO+A33*ZGEO
      ELSE
       XGEO=A11*XGSW+A21*YGSW+A31*ZGSW
       YGEO=A12*XGSW+A22*YGSW+A32*ZGSW
       ZGEO=A13*XGSW+A23*YGSW+A33*ZGSW
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
      SUBROUTINE GEODGEO_08 (H,XMU,R,THETA,J)
C
C  THIS SUBROUTINE (1) CONVERTS VERTICAL LOCAL HEIGHT (ALTITUDE) H AND GEODETIC
C  LATITUDE XMU INTO GEOCENTRIC COORDINATES R AND THETA (GEOCENTRIC RADIAL
C  DISTANCE AND COLATITUDE, RESPECTIVELY; ALSO KNOWN AS ECEF COORDINATES),
C  AS WELL AS (2) PERFORMS THE INVERSE TRANSFORMATION FROM {R,THETA} TO {H,XMU}.
C
C  THE SUBROUTINE USES WORLD GEODETIC SYSTEM WGS84 PARAMETERS FOR THE EARTH'S
C  ELLIPSOID. THE ANGULAR QUANTITIES (GEO COLATITUDE THETA AND GEODETIC LATITUDE
C  XMU) ARE IN RADIANS, AND THE DISTANCES (GEOCENTRIC RADIUS R AND ALTITUDE H
C  ABOVE THE EARTH'S ELLIPSOID) ARE IN KILOMETERS.
C
C  IF J>0, THE TRANSFORMATION IS MADE FROM GEODETIC TO GEOCENTRIC COORDINATES
C   USING SIMPLE DIRECT EQUATIONS.
C  IF J<0, THE INVERSE TRANSFORMATION FROM GEOCENTRIC TO GEODETIC COORDINATES
C   IS MADE BY MEANS OF A FAST ITERATIVE ALGORITHM.
C
c-------------------------------------------------------------------------------
C                   J>0                     |            J<0
c-------------------------------------------|-----------------------------------
C--INPUT:   J        H          XMU         |    J         R          THETA
c         flag  altitude (km)  geodetic     |   flag   geocentric    spherical
c                              latitude     |         distance (km) colatitude
c                              (radians)    |                        (radians)
c-------------------------------------------|-----------------------------------
c                                           |
C----OUTPUT:         R           THETA      |          H              XMU
C                geocentric    spherical    |      altitude (km)    geodetic
C                distance (km) colatitude   |                       latitude
C                              (radians)    |                       (radians)
C-------------------------------------------------------------------------------
C
C   AUTHOR:  N. A. TSYGANENKO
c   DATE:    DEC 5, 2007
C   LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA R_EQ, BETA /6378.137D0, 6.73949674228D-3/
c
c  R_EQ is the semi-major axis of the Earth's ellipsoid, and BETA is its
c  second eccentricity squared
c
      DATA TOL /1.D-6/
c
c   Direct transformation (GEOD=>GEO):
c
      IF (J.GT.0) THEN
       COSXMU=DCOS(XMU)
       SINXMU=DSIN(XMU)
       DEN=DSQRT(COSXMU**2+(SINXMU/(1.0D0+BETA))**2)
       COSLAM=COSXMU/DEN
       SINLAM=SINXMU/(DEN*(1.0D0+BETA))
       RS=R_EQ/DSQRT(1.0D0+BETA*SINLAM**2)
       X=RS*COSLAM+H*COSXMU
       Z=RS*SINLAM+H*SINXMU
       R=DSQRT(X**2+Z**2)
       THETA=DACOS(Z/R)
      ENDIF

c
c   Inverse transformation (GEO=>GEOD):
c
      IF (J.LT.0) THEN
       N=0
       PHI=1.570796327D0-THETA
       PHI1=PHI
  1    SP=DSIN(PHI1)
       ARG=SP*(1.0D0+BETA)/DSQRT(1.0D0+BETA*(2.0D0+BETA)*SP**2)
       XMUS=DASIN(ARG)
       RS=R_EQ/DSQRT(1.0D0+BETA*DSIN(PHI1)**2)
       COSFIMS=DCOS(PHI1-XMUS)
       H=DSQRT((RS*COSFIMS)**2+R**2-RS**2)-RS*COSFIMS
       Z=RS*DSIN(PHI1)+H*DSIN(XMUS)
       X=RS*DCOS(PHI1)+H*DCOS(XMUS)
       RR=DSQRT(X**2+Z**2)
       DPHI=DASIN(Z/RR)-PHI
       PHI1=PHI1-DPHI
       N=N+1
       IF (DABS(DPHI).GT.TOL.AND.N.LT.100) GOTO 1
       XMU=XMUS
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
      SUBROUTINE RHAND_08 (X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME,INNAME)
C
C  CALCULATES THE COMPONENTS OF THE RIGHT HAND SIDE VECTOR IN THE GEOMAGNETIC FIELD
C    LINE EQUATION  (a subsidiary subroutine for the subroutine STEP_08)
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PARMOD(10)
C
C     EXNAME AND INNAME ARE NAMES OF SUBROUTINES FOR THE EXTERNAL AND INTERNAL
C     PARTS OF THE TOTAL FIELD, E.G., T96_01 AND IGRF_GSW_08
C
      COMMON /GEOPACK1/ A(12),DS3,BB(2),PSI,CC(18)

      CALL EXNAME (IOPT,PARMOD,PSI,X,Y,Z,BXGSW,BYGSW,BZGSW)
      CALL INNAME (X,Y,Z,HXGSW,HYGSW,HZGSW)

      BX=BXGSW+HXGSW
      BY=BYGSW+HYGSW
      BZ=BZGSW+HZGSW
      B=DS3/DSQRT(BX**2+BY**2+BZ**2)
      R1=BX*B
      R2=BY*B
      R3=BZ*B
      RETURN
      END
C
C===================================================================================
C
      SUBROUTINE STEP_08(X,Y,Z,DS,DSMAX,ERRIN,IOPT,PARMOD,EXNAME,INNAME)
C
C   RE-CALCULATES THE INPUT VALUES {X,Y,Z} (IN GSW COORDINATES) FOR ANY POINT ON A FIELD LINE,
C     BY MAKING A STEP ALONG THAT LINE USING RUNGE-KUTTA-MERSON ALGORITHM (G.N. Lance, Numerical
C      methods for high-speed computers, Iliffe & Sons, London 1960.)
C   DS IS A PRESCRIBED VALUE OF THE CURRENT STEP SIZE, DSMAX IS ITS UPPER LIMIT.
C   ERRIN IS A PERMISSIBLE ERROR (ITS OPTIMAL VALUE SPECIFIED IN THE S/R TRACE_08)
C     IF THE ACTUAL ERROR (ERRCUR) AT THE CURRENT STEP IS LARGER THAN ERRIN, THE STEP IS REJECTED,
C       AND THE CALCULATION IS REPEATED ANEW WITH HALVED STEPSIZE DS.
C     IF ERRCUR IS SMALLER THAN ERRIN, THE STEP IS ACCEPTED, AND THE CURRENT VALUE OF DS IS RETAINED
C       FOR THE NEXT STEP.
C     IF ERRCUR IS SMALLER THAN 0.04*ERRIN, THE STEP IS ACCEPTED, AND THE VALUE OF DS FOR THE NEXT STEP
C       IS INCREASED BY THE FACTOR 1.5, BUT NOT LARGER THAN DSMAX.
C   IOPT IS A FLAG, RESERVED FOR SPECIFYNG A VERSION OF THE EXTERNAL FIELD MODEL EXNAME.
C   ARRAY PARMOD(10) CONTAINS INPUT PARAMETERS FOR THE MODEL EXNAME.
C   EXNAME IS THE NAME OF THE SUBROUTINE FOR THE EXTERNAL FIELD MODEL.
C   INNAME IS THE NAME OF THE SUBROUTINE FOR THE INTERNAL FIELD MODEL (EITHER DIP_08 OR IGRF_GSW_08)
C
C   ALL THE ABOVE PARAMETERS ARE INPUT ONES; OUTPUT IS THE RECALCULATED VALUES OF X,Y,Z
C
C     LAST MODIFICATION:    APRIL 21, 2008 (SEE ERRATA AS OF THAT DATE)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PARMOD(10)
      COMMON /GEOPACK1/ A(12),DS3,B(21)
      EXTERNAL EXNAME,INNAME

  1   DS3=-DS/3.D0
      CALL RHAND_08 (X,Y,Z,R11,R12,R13,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND_08 (X+R11,Y+R12,Z+R13,R21,R22,R23,IOPT,PARMOD,EXNAME,
     * INNAME)
      CALL RHAND_08 (X+.5D0*(R11+R21),Y+.5D0*(R12+R22),Z+.5D0*
     *(R13+R23),R31,R32,R33,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND_08 (X+.375D0*(R11+3.D0*R31),Y+.375D0*(R12+3.D0*R32
     *),Z+.375D0*(R13+3.D0*R33),R41,R42,R43,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND_08 (X+1.5D0*(R11-3.D0*R31+4.D0*R41),Y+1.5D0*(R12-
     *3.D0*R32+4.D0*R42),Z+1.5D0*(R13-3.D0*R33+4.D0*R43),
     *R51,R52,R53,IOPT,PARMOD,EXNAME,INNAME)
      ERRCUR=DABS(R11-4.5D0*R31+4.D0*R41-.5D0*R51)+DABS(R12-4.5D0*R32
     *+4.D0*R42-.5D0*R52)+DABS(R13-4.5D0*R33+4.D0*R43-.5D0*R53)
C
C  READY FOR MAKING THE STEP, BUT CHECK THE ACCURACY; IF INSUFFICIENT,
C   REPEAT THE STEP WITH HALVED STEPSIZE:
C
      IF (ERRCUR.GT.ERRIN) THEN
      DS=DS*.5D0
      GOTO 1
      ENDIF
C
C  ACCURACY IS ACCEPTABLE, BUT CHECK IF THE STEPSIZE IS NOT TOO LARGE;
C    OTHERWISE REPEAT THE STEP WITH DS=DSMAX
C
      IF (DABS(DS).GT.DSMAX) THEN
      DS=DSIGN(DSMAX,DS)
      GOTO 1
      ENDIF
C
C  MAKING THE STEP:
C
  2   X=X+.5D0*(R11+4.D0*R41+R51)
      Y=Y+.5D0*(R12+4.D0*R42+R52)
      Z=Z+.5D0*(R13+4.D0*R43+R53)
C
C  IF THE ACTUAL ERROR IS TOO SMALL (LESS THAN 4% OF ERRIN) AND DS SMALLER
C   THAN DSMAX/1.5, THEN WE INCREASE THE STEPSIZE FOR THE NEXT STEP BY 50%
C
      IF(ERRCUR.LT.ERRIN*.04D0.AND.DS.LT.DSMAX/1.5D0) DS=DS*1.5D0
      RETURN
      END
C
C==============================================================================
C
      SUBROUTINE TRACE_08 (XI,YI,ZI,DIR,DSMAX,ERR,RLIM,R0,IOPT,PARMOD,
     * EXNAME,INNAME,XF,YF,ZF,XX,YY,ZZ,L,LMAX)
C
C  TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
C  SURFACE OR TO A MODEL LIMITING BOUNDARY.
C
C  THIS SUBROUTINE ALLOWS TWO OPTIONS:
C
C  (1) IF INNAME=IGRF_GSW_08, THEN THE IGRF MODEL WILL BE USED FOR CALCULATING
C      CONTRIBUTION FROM EARTH'S INTERNAL SOURCES. IN THIS CASE, SUBROUTINE
C      RECALC_08 MUST BE CALLED BEFORE USING TRACE_08, WITH PROPERLY SPECIFIED DATE,
C      UNIVERSAL TIME, AND SOLAR WIND VELOCITY COMPONENTS, TO CALCULATE IN ADVANCE
C      ALL QUANTITIES NEEDED FOR THE MAIN FIELD MODEL AND FOR TRANSFORMATIONS
C      BETWEEN INVOLVED COORDINATE SYSTEMS.
C
C  (2) IF INNAME=DIP_08, THEN A PURE DIPOLE FIELD WILL BE USED INSTEAD OF THE IGRF MODEL.
C      IN THIS CASE, THE SUBROUTINE RECALC_08 MUST ALSO BE CALLED BEFORE TRACE_08.
C      HERE ONE CAN CHOOSE EITHER TO
C      (a) CALCULATE DIPOLE TILT ANGLE BASED ON DATE, TIME, AND SOLAR WIND DIRECTION,
C   OR (b) EXPLICITLY SPECIFY THAT ANGLE, WITHOUT ANY REFERENCE TO DATE/UT/SOLAR WIND.
C      IN THE LAST CASE (b), THE SINE (SPS) AND COSINE (CPS) OF THE DIPOLE TILT
C      ANGLE MUST BE SPECIFIED IN ADVANCE (BUT AFTER HAVING CALLED RECALC_08) AND FORWARDED
C      IN THE COMMON BLOCK /GEOPACK1/ (IN ITS 11th AND 12th ELEMENTS, RESPECTIVELY).
C      IN THIS CASE THE ROLE OF THE SUBROUTINE RECALC_08 IS REDUCED TO ONLY CALCULATING
C      THE COMPONENTS OF THE EARTH'S DIPOLE MOMENT.
C
C------------- INPUT PARAMETERS:
C
C   XI,YI,ZI - GSW COORDS OF THE FIELD LINE STARTING POINT (IN EARTH RADII, 1 RE = 6371.2 km),
C
C   DIR - SIGN OF THE TRACING DIRECTION: IF DIR=1.0 THEN THE TRACING IS MADE ANTIPARALLEL
C     TO THE TOTAL FIELD VECTOR (E.G., FROM NORTHERN TO SOUTHERN CONJUGATE POINT);
C     IF DIR=-1.0 THEN THE TRACING PROCEEDS IN THE OPPOSITE DIRECTION, THAT IS, PARALLEL TO
C     THE TOTAL FIELD VECTOR.
C
C   DSMAX - UPPER LIMIT ON THE STEPSIZE (SETS A DESIRED MAXIMAL SPACING BETWEEN
C                 THE FIELD LINE POINTS)
C
C   ERR - PERMISSIBLE STEP ERROR. A REASONABLE ESTIMATE PROVIDING A SUFFICIENT ACCURACY FOR MOST
C         APPLICATIONS IS ERR=0.0001. SMALLER/LARGER VALUES WILL RESULT IN LARGER/SMALLER NUMBER
C         OF STEPS AND, HENCE, OF OUTPUT FIELD LINE POINTS. NOTE THAT USING MUCH SMALLER VALUES
C         OF ERR MAY REQUIRE USING A DOUBLE PRECISION VERSION OF THE ENTIRE PACKAGE.
C
C   R0 -  RADIUS OF A SPHERE (IN RE), DEFINING THE INNER BOUNDARY OF THE TRACING REGION
C         (USUALLY, EARTH'S SURFACE OR THE IONOSPHERE, WHERE R0~1.0)
C         IF THE FIELD LINE REACHES THAT SPHERE FROM OUTSIDE, ITS INBOUND TRACING IS
C         TERMINATED AND THE CROSSING POINT COORDINATES XF,YF,ZF  ARE CALCULATED.
C
C   RLIM - RADIUS OF A SPHERE (IN RE), DEFINING THE OUTER BOUNDARY OF THE TRACING REGION;
C         IF THE FIELD LINE REACHES THAT BOUNDARY FROM INSIDE, ITS OUTBOUND TRACING IS
C         TERMINATED AND THE CROSSING POINT COORDINATES XF,YF,ZF ARE CALCULATED.
C
C   IOPT - A MODEL INDEX; CAN BE USED FOR SPECIFYING A VERSION OF THE EXTERNAL FIELD
C       MODEL (E.G., A NUMBER OF THE KP-INDEX INTERVAL). ALTERNATIVELY, ONE CAN USE THE ARRAY
C       PARMOD FOR THAT PURPOSE (SEE BELOW); IN THAT CASE IOPT IS JUST A DUMMY PARAMETER.
C
C   PARMOD -  A 10-ELEMENT ARRAY CONTAINING INPUT PARAMETERS NEEDED FOR A UNIQUE
C      SPECIFICATION OF THE EXTERNAL FIELD MODEL. THE CONCRETE MEANING OF THE COMPONENTS
C      OF PARMOD DEPENDS ON A SPECIFIC VERSION OF THAT MODEL.
C
C   EXNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE EXTERNAL MAGNETIC FIELD
C    (E.G., T89, OR T96_01, ETC.).
C   INNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE INTERNAL MAGNETIC FIELD
C    (EITHER DIP_08 OR IGRF_GSW_08).
C
C   LMAX - MAXIMAL LENGTH OF THE ARRAYS XX,YY,ZZ, IN WHICH COORDINATES OF THE FIELD
C          LINE POINTS ARE STORED. LMAX SHOULD BE SET EQUAL TO THE ACTUAL LENGTH OF
C          THE ARRAYS, DEFINED IN THE MAIN PROGRAM AS ACTUAL ARGUMENTS OF THIS SUBROUTINE.
C
C-------------- OUTPUT PARAMETERS:
C
C   XF,YF,ZF - GSW COORDINATES OF THE ENDPOINT OF THE TRACED FIELD LINE.
C   XX,YY,ZZ - ARRAYS OF LENGTH LMAX, CONTAINING COORDINATES OF THE FIELD LINE POINTS.
C   L - ACTUAL NUMBER OF FIELD LINE POINTS, GENERATED BY THIS SUBROUTINE.
C
C ----------------------------------------------------------
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XX(LMAX),YY(LMAX),ZZ(LMAX), PARMOD(10)
      COMMON /GEOPACK1/ AA(12),DD,BB(21)
      EXTERNAL EXNAME,INNAME
C
      L=0
      NREV=0
      DD=DIR
C
C  INITIALIZE THE STEP SIZE AND STARTING PONT:
C
      DS=0.5D0*DIR
      X=XI
      Y=YI
      Z=ZI
c
c  here we call RHAND_08 just to find out the sign of the radial component of the field
c   vector, and to determine the initial direction of the tracing (i.e., either away
c   or towards Earth):
c
      CALL RHAND_08 (X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME,INNAME)
      AD=0.01D0
      IF (X*R1+Y*R2+Z*R3.LT.0.D0) AD=-0.01D0
C
c     |AD|=0.01 and its sign follows the rule:
c (1) if DIR=1 (tracing antiparallel to B vector) then the sign of AD is the same as of Br
c (2) if DIR=-1 (tracing parallel to B vector) then the sign of AD is opposite to that of Br
c     AD is defined in order to initialize the value of RR (radial distance at previous step):

      RR=DSQRT(X**2+Y**2+Z**2)+AD
c
  1   L=L+1
      IF(L.GT.LMAX) GOTO 7
      XX(L)=X
      YY(L)=Y
      ZZ(L)=Z
      RYZ=Y**2+Z**2
      R2=X**2+RYZ
      R=DSQRT(R2)
C
c  check if the line hit the outer tracing boundary; if yes, then terminate
c   the tracing (label 8). The outer boundary is assumed reached, when the line
c   crosses any of the 3 surfaces: (1) a sphere R=RLIM, (2) a cylinder of radius 40Re,
c   coaxial with the XGSW axis, (3) the plane X=20Re:

      IF (R.GT.RLIM.OR.RYZ.GT.1600.D0.OR.X.GT.20.D0) GOTO 8
c
c  check whether or not the inner tracing boundary was crossed from outside,
c  if yes, then calculate the footpoint position by interpolation (go to label 6):
c
      IF (R.LT.R0.AND.RR.GT.R) GOTO 6

c  check if we are moving outward, or R is still larger than 3Re; if yes, proceed further:
c
      IF (R.GE.RR.OR.R.GE.3.D0) GOTO 4
c
c  now we entered inside the sphere R=3: to avoid too large steps (and hence
c  inaccurate interpolated position of the footpoint), enforce the progressively
c  smaller stepsize values as we approach the inner boundary R=R0:
c
      FC=0.2D0
      IF(R-R0.LT.0.05D0) FC=0.05D0
      AL=FC*(R-R0+0.2D0)
      DS=DIR*AL
c
  4   XR=X
      YR=Y
      ZR=Z
c
      DRP=R-RR
      RR=R
c
      CALL STEP_08 (X,Y,Z,DS,DSMAX,ERR,IOPT,PARMOD,EXNAME,INNAME)
c
C  check the total number NREV of changes in the tracing radial direction; (NREV.GT.2) means
c   that the line started making multiple loops, in which case we stop the process:
C
      R=DSQRT(X**2+Y**2+Z**2)
      DR=R-RR
      IF (DRP*DR.LT.0.D0) NREV=NREV+1
      IF (NREV.GT.4) GOTO 8
C
      GOTO 1
c
c  find the footpoint position by interpolating between the current and previous
c   field line points:
c
  6   R1=(R0-R)/(RR-R)
      X=X-(X-XR)*R1
      Y=Y-(Y-YR)*R1
      Z=Z-(Z-ZR)*R1
      GOTO 8
  7   WRITE (*,10)
      L=LMAX
  8   XF=X
      YF=Y
      ZF=Z
C
C  replace the coordinates of the last (L-th) point in the XX,YY,ZZ arrays
C   so that they correspond to the estimated footpoint position {XF,YF,ZF},
c   satisfying:  sqrt(XF**2+YF**2+ZF**2}=R0
C
      XX(L)=XF
      YY(L)=YF
      ZZ(L)=ZF
C
      RETURN
 10   FORMAT(//,1X,'**** COMPUTATIONS IN THE SUBROUTINE TRACE_08 ARE',
     *' TERMINATED: THE NUMBER OF POINTS EXCEEDED LMAX ****'//)
      END
c
C====================================================================================
C
      SUBROUTINE SHUETAL_MGNP_08(XN_PD,VEL,BZIMF,XGSW,YGSW,ZGSW,
     *  XMGNP,YMGNP,ZMGNP,DIST,ID)
C
C  FOR ANY POINT OF SPACE WITH COORDINATES (XGSW,YGSW,ZGSW) AND SPECIFIED CONDITIONS
C  IN THE INCOMING SOLAR WIND, THIS SUBROUTINE:
C
C (1) DETERMINES IF THE POINT (XGSW,YGSW,ZGSW) LIES INSIDE OR OUTSIDE THE
C      MODEL MAGNETOPAUSE OF SHUE ET AL. (JGR-A, V.103, P. 17691, 1998).
C
C (2) CALCULATES THE GSW POSITION OF A POINT {XMGNP,YMGNP,ZMGNP}, LYING AT THE MODEL
C      MAGNETOPAUSE AND ASYMPTOTICALLY TENDING TO THE NEAREST BOUNDARY POINT WITH
C      RESPECT TO THE OBSERVATION POINT {XGSW,YGSW,ZGSW}, AS IT APPROACHES THE MAGNETO-
C      PAUSE.
C
C  INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         BZIMF - IMF BZ IN NANOTESLAS
C
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSW,YGSW,ZGSW - GSW POSITION OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSW POSITION OF THE BOUNDARY POINT
C          DIST - DISTANCE (IN RE) BETWEEN THE OBSERVATION POINT (XGSW,YGSW,ZGSW)
C                 AND THE MODEL NAGNETOPAUSE
C          ID -  POSITION FLAG:  ID=+1 (-1) MEANS THAT THE OBSERVATION POINT
C          LIES INSIDE (OUTSIDE) OF THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  OTHER SUBROUTINES USED: T96_MGNP_08
C
c          AUTHOR:  N.A. TSYGANENKO,
C          DATE:    APRIL 4, 2003.
C          LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF (VEL.LT.0.) THEN
        P=XN_PD
      ELSE
        P=1.94D-6*XN_PD*VEL**2  ! P IS THE SOLAR WIND DYNAMIC PRESSURE (IN nPa)
      ENDIF

c
c  DEFINE THE ANGLE PHI, MEASURED DUSKWARD FROM THE NOON-MIDNIGHT MERIDIAN PLANE;
C  IF THE OBSERVATION POINT LIES ON THE X AXIS, THE ANGLE PHI CANNOT BE UNIQUELY
C  DEFINED, AND WE SET IT AT ZERO:
c
      IF (YGSW.NE.0.D0.OR.ZGSW.NE.0.D0) THEN
         PHI=DATAN2(YGSW,ZGSW)
      ELSE
         PHI=0.D0
      ENDIF
C
C  FIRST, FIND OUT IF THE OBSERVATION POINT LIES INSIDE THE SHUE ET AL BDRY
C  AND SET THE VALUE OF THE ID FLAG:
C
      ID=-1
      R0=(10.22D0+1.29D0*DTANH(.184D0*(BZIMF+8.14D0)))*P**(-.15151515D0)
      ALPHA=(0.58D0-0.007D0*BZIMF)*(1.D0+0.024*DLOG(P))
      R=DSQRT(XGSW**2+YGSW**2+ZGSW**2)
      RM=R0*(2.D0/(1.D0+XGSW/R))**ALPHA
      IF (R.LE.RM) ID=+1
C
C  NOW, FIND THE CORRESPONDING T96 MAGNETOPAUSE POSITION, TO BE USED AS
C  A STARTING APPROXIMATION IN THE SEARCH OF A CORRESPONDING SHUE ET AL.
C  BOUNDARY POINT:
C
      CALL T96_MGNP_08(P,-1.D0,XGSW,YGSW,ZGSW,XMT96,YMT96,ZMT96,DIST,
     *  ID96)
C
      RHO2=YMT96**2+ZMT96**2
      R=DSQRT(RHO2+XMT96**2)
      ST=DSQRT(RHO2)/R
      CT=XMT96/R
C
C  NOW, USE NEWTON'S ITERATIVE METHOD TO FIND THE NEAREST POINT AT THE
C   SHUE ET AL.'S BOUNDARY:
C
      NIT=0

  1   T=DATAN2(ST,CT)
      RM=R0*(2.D0/(1.D0+CT))**ALPHA

      F=R-RM
      GRADF_R=1.D0
      GRADF_T=-ALPHA/R*RM*ST/(1.D0+CT)
      GRADF=DSQRT(GRADF_R**2+GRADF_T**2)

      DR=-F/GRADF**2
      DT= DR/R*GRADF_T

      R=R+DR
      T=T+DT
      ST=DSIN(T)
      CT=DCOS(T)

      DS=DSQRT(DR**2+(R*DT)**2)

      NIT=NIT+1

      IF (NIT.GT.1000) THEN
         PRINT *,
     *' BOUNDARY POINT COULD NOT BE FOUND; ITERATIONS DO NOT CONVERGE'
      ENDIF

      IF (DS.GT.1.D-4) GOTO 1

      XMGNP=R*DCOS(T)
      RHO=  R*DSIN(T)

      YMGNP=RHO*DSIN(PHI)
      ZMGNP=RHO*DCOS(PHI)

      DIST=DSQRT((XGSW-XMGNP)**2+(YGSW-YMGNP)**2+(ZGSW-ZMGNP)**2)

      RETURN
      END
C
C=======================================================================================
C
      SUBROUTINE T96_MGNP_08(XN_PD,VEL,XGSW,YGSW,ZGSW,XMGNP,YMGNP,ZMGNP,
     * DIST,ID)
C
C  FOR ANY POINT OF SPACE WITH GIVEN COORDINATES (XGSW,YGSW,ZGSW), THIS SUBROUTINE DEFINES
C  THE POSITION OF A POINT (XMGNP,YMGNP,ZMGNP) AT THE T96 MODEL MAGNETOPAUSE WITH THE
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
C         XGSW,YGSW,ZGSW - COORDINATES OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSW POSITION OF THE BOUNDARY POINT, HAVING THE SAME
C          VALUE OF TAU-COORDINATE AS THE OBSERVATION POINT (XGSW,YGSW,ZGSW)
C          DIST -  THE DISTANCE BETWEEN THE TWO POINTS, IN RE,
C          ID -    POSITION FLAG; ID=+1 (-1) MEANS THAT THE POINT (XGSW,YGSW,ZGSW)
C          LIES INSIDE (OUTSIDE) THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  THE PRESSURE-DEPENDENT MAGNETOPAUSE IS THAT USED IN THE T96_01 MODEL
C  (TSYGANENKO, JGR, V.100, P.5599, 1995; ESA SP-389, P.181, OCT. 1996)
C
c   AUTHOR:  N.A. TSYGANENKO
C   DATE:    AUG.1, 1995, REVISED APRIL 3, 2003.
C   LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C  DEFINE SOLAR WIND DYNAMIC PRESSURE (NANOPASCALS, ASSUMING 4% OF ALPHA-PARTICLES),
C   IF NOT EXPLICITLY SPECIFIED IN THE INPUT:
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
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
      A0=70.D0
      S00=1.08D0
      X00=5.48D0
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
C     (FOR DETAILS OF THE ELLIPSOIDAL COORDINATES, SEE THE PAPER:
C      N.A.TSYGANENKO, SOLUTION OF CHAPMAN-FERRARO PROBLEM FOR AN
C      ELLIPSOIDAL MAGNETOPAUSE, PLANET.SPACE SCI., V.37, P.1037, 1989).
C
       IF (YGSW.NE.0.D0.OR.ZGSW.NE.0.D0) THEN
          PHI=DATAN2(YGSW,ZGSW)
       ELSE
          PHI=0.D0
       ENDIF
C
       RHO=DSQRT(YGSW**2+ZGSW**2)
C
       IF (XGSW.LT.XM) THEN
           XMGNP=XGSW
           RHOMGNP=A*DSQRT(S0**2-1.D0)
           YMGNP=RHOMGNP*DSIN(PHI)
           ZMGNP=RHOMGNP*DCOS(PHI)
           DIST=DSQRT((XGSW-XMGNP)**2+(YGSW-YMGNP)**2+(ZGSW-ZMGNP)**2)
           IF (RHOMGNP.GT.RHO) ID=+1
           IF (RHOMGNP.LE.RHO) ID=-1
           RETURN
       ENDIF
C
          XKSI=(XGSW-X0)/A+1.D0
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
C  NOW CALCULATE THE DISTANCE BETWEEN THE POINTS {XGSW,YGSW,ZGSW} AND {XMGNP,YMGNP,ZMGNP}:
C   (IN GENERAL, THIS IS NOT THE SHORTEST DISTANCE D_MIN, BUT DIST ASYMPTOTICALLY TENDS
C    TO D_MIN, AS WE ARE GETTING CLOSER TO THE MAGNETOPAUSE):
C
      DIST=DSQRT((XGSW-XMGNP)**2+(YGSW-YMGNP)**2+(ZGSW-ZMGNP)**2)
C
      IF (SIGMA.GT.S0) ID=-1   !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE
      IF (SIGMA.LE.S0) ID=+1   !  ID=+1 MEANS THAT THE POINT LIES INSIDE
C                                           THE MAGNETOSPHERE
      RETURN
      END
C
C===================================================================================
C
c

!***************************************************************************************************
! Copyright 2007,2009 S. Bourdarie
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
!***************************************************************************************************
! CREATION: S. Bourdarie - ONERA-DESP
!
! FILE CONTENT:
!               SUBROUTINE coord_trans1: Generic coordinate transformation from one Earth or Heliospheric coordinate
!                                        system to another one
!               SUBROUTINE coord_trans_vec1: Generic coordinate transformation from one Earth or Heliospheric coordinate
!                                        system to another one (handle up to
!                                        ntime_max positions)
!
!***************************************************************************************************
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE PNM80
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009. Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C  Form the matrix of precession/nutation for a given date, IAU 1976
C  precession model, IAU 1980 nutation model.
C
C The matrix operates in the sense V(date) = RMATPN * V(J2000),
C     where the p-vector V(date) is with respect to the true
C     equatorial triad of date DATE1+DATE2 and the p-vector
C     V(J2000) is with respect to the mean equatorial triad of
C     epoch J2000.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PNM80 ( DATE1, DATE2, RMATPN )

      IMPLICIT NONE

      REAL*8 DATE1, DATE2, RMATPN(3,3)
      REAL*8 RMATP(3,3), RMATN(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Precession matrix, J2000 to date.
      CALL PMAT76 ( DATE1, DATE2, RMATP )

*  Nutation matrix.
      CALL NUTM80 ( DATE1, DATE2, RMATN )

*  Combine the matrices:  PN = N x P.
      CALL RXR ( RMATN, RMATP, RMATPN )

      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE PMAT76
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C  Precession matrix from J2000 to a specified date, IAU 1976 model.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PMAT76 ( DATE1, DATE2, RMATP )

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, RMATP(3,3)

*  Reference epoch (J2000), JD
      REAL*8 DJ00
      PARAMETER ( DJ00 = 2451545D0 )

      REAL*8 ZETA, Z, THETA, WMAT(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Precession Euler angles, J2000 to specified date.
      CALL PREC76 ( DJ00, 0D0, DATE1, DATE2, ZETA, Z, THETA )

*  Form the rotation matrix.
      CALL IR ( WMAT )
      CALL RZ ( -ZETA, WMAT )
      CALL RY ( THETA, WMAT )
      CALL RZ ( -Z, WMAT )
      CALL CR ( WMAT, RMATP )


      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE IR
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 and
C     INTEGER by INTEGER*4 by S. Bourdarie
C
C  Initialize an r-matrix to the identity matrix.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE IR ( R )
      IMPLICIT NONE

      REAL*8 R(3,3)

      INTEGER*4 I

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      CALL ZR ( R )
      DO 1 I=1,3
         R(I,I) = 1D0
 1    CONTINUE
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE ZR
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 and
C     INTEGER by INTEGER*4 by S. Bourdarie
C
C  Initialize an r-matrix to the null matrix.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ZR ( R )
      IMPLICIT NONE

      REAL*8  R(3,3)

      INTEGER*4 I, J

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 2 J=1,3
         DO 1 I=1,3
            R(I,J) = 0.D0
 1       CONTINUE
 2    CONTINUE
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE RZ
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
CRotate an r-matrix about the z-axis.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RZ ( PSI, R )
      IMPLICIT NONE

      REAL*8 PSI, R(3,3)

      REAL*8 S, C, A(3,3), W(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Matrix representing new rotation.
      S = SIN(PSI)
      C = COS(PSI)
      CALL IR ( A )
      A(1,1) = C
      A(2,1) = -S
      A(1,2) = S
      A(2,2) = C

*  Rotate.
      CALL RXR ( A, R, W )

*  Return result.
      CALL CR ( W, R )
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE RY
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C  Rotate an r-matrix about the y-axis.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RY ( THETA, R )

      IMPLICIT NONE
      REAL*8 THETA, R(3,3)

      REAL*8 S, C, A(3,3), W(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Matrix representing new rotation.
      S = SIN(THETA)
      C = COS(THETA)
      CALL IR ( A )
      A(1,1) = C
      A(3,1) = S
      A(1,3) = -S
      A(3,3) = C

*  Rotate.
      CALL RXR ( A, R, W )

*  Return result.
      CALL CR ( W, R )
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE RY
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C  Rotate an r-matrix about the x-axis.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RX ( PHI, R )
      IMPLICIT NONE

      REAL*8 PHI, R(3,3)

      REAL*8 S, C, A(3,3), W(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Matrix representing new rotation.
      S = SIN(PHI)
      C = COS(PHI)
      CALL IR ( A )
      A(2,2) = C
      A(3,2) = -S
      A(2,3) = S
      A(3,3) = C

*  Rotate.
      CALL RXR ( A, R, W )

*  Return result.
      CALL CR ( W, R )

      end


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE RXR
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C     Multiply two r-matrices.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RXR ( A, B, ATB )
      IMPLICIT NONE

      REAL*8 A(3,3), B(3,3), ATB(3,3)

      INTEGER*4 I, J, K
      REAL*8 W, WM(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 3 I=1,3
         DO 2 J=1,3
            W = 0D0
            DO 1 K=1,3
               W = W + A(I,K)*B(K,J)
 1          CONTINUE
            WM(I,J) = W
 2       CONTINUE
 3    CONTINUE
      CALL CR ( WM, ATB )
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE CR
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The changes to original routine are:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C     Subroutine CR and CP have been combined together
C
C     Copy an r-matrix.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CR ( R, C )
      IMPLICIT NONE

      REAL*8 R(3,3), C(3,3)

      INTEGER*4 I,J

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 1 I=1,3
         DO J=1,3
           C(I,J)=R(I,J)
         ENDDO
 1    CONTINUE
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE PREC76
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PREC76 ( EP01, EP02, EP11, EP12, ZETA, Z, THETA )
*+
*-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL*8 EP01, EP02, EP11, EP12, ZETA, Z, THETA

*  Arcseconds to radians
      REAL*8 DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000), JD
      REAL*8 DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      REAL*8 DJC
      PARAMETER ( DJC = 36525D0 )

      REAL*8 T0, T, TAS2R, W


*  Interval between fundamental epoch J2000.0 and beginning epoch (JC).
      T0 = ( ( EP01-DJ00 ) + EP02 ) / DJC

*  Interval over which precession required (JC).
      T = ( ( EP11-EP01 ) + ( EP12-EP02 ) ) / DJC

*  Euler angles.
      TAS2R = T * DAS2R
      W = 2306.2181D0 + (
     :       1.39656D0
     :     - 0.000139D0 * T0 ) * T0

      ZETA = ( W + ( ( 0.30188D0
     :               - 0.000344D0 * T0 )
     :               + 0.017998D0 * T ) * T ) * TAS2R

      Z = ( W + ( ( 1.09468D0
     :            + 0.000066D0 * T0 )
     :            + 0.018203D0 * T ) * T ) * TAS2R

      THETA = ( ( 2004.3109D0 + (
     :             - 0.85330D0
     :             - 0.000217D0 * T0 ) * T0 ) + ( (
     :             - 0.42665D0
     :             - 0.000217D0 * T0 )
     :             - 0.041833D0 * T ) * T ) * TAS2R
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE NUTM80
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C   Form the matrix of nutation for a given date, IAU 1980 model.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE NUTM80 ( DATE1, DATE2, RMATN )

      IMPLICIT NONE

      REAL*8 DATE1, DATE2, RMATN(3,3)

      REAL*8 DPSI, DEPS, EPSA
      REAL*8 OBL80

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Nutation components and mean obliquity.
      CALL NUT80 ( DATE1, DATE2, DPSI, DEPS )
      EPSA = OBL80 ( DATE1, DATE2 )

*  Build the rotation matrix.
      CALL NUMAT ( EPSA, DPSI, DEPS, RMATN )

      end
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE NUTM80
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C     Nutation, IAU 1980 model.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE NUT80 ( DATE1, DATE2, DPSI, DEPS )
      IMPLICIT NONE

      real*8 DATE1, DATE2, DPSI, DEPS

*  Arcseconds to radians
      real*8 DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  2Pi
      real*8 D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Units of 0.1 milliarcsecond to radians
      real*8 U2R
      PARAMETER ( U2R = DAS2R/1D4 )

*  Reference epoch (J2000), JD
      real*8 DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      real*8 DJC
      PARAMETER ( DJC = 36525D0 )

      real*8 T, EL, ELP, F, D, OM, DP, DE, ARG, S, C
      INTEGER*4 I, J

      real*8 ANPM


*  ------------------------------------------------
*  Table of multiples of arguments and coefficients
*  ------------------------------------------------
*
*  The coefficient values are in 0.1 mas units and the rates of change
*  are in mas per Julian millennium.

      REAL X(9,106)

*                Multiple of            Longitude        Obliquity
*           L    L'   F    D  Omega   coeff. of sin    coeff. of cos
*                                         1       t        1     t

      DATA ((X(I,J),I=1,9),J=1,10) /
     :      0.,  0.,  0.,  0.,  1., -171996., -1742.,  92025.,  89.,
     :      0.,  0.,  0.,  0.,  2.,    2062.,     2.,   -895.,   5.,
     :     -2.,  0.,  2.,  0.,  1.,      46.,     0.,    -24.,   0.,
     :      2.,  0., -2.,  0.,  0.,      11.,     0.,      0.,   0.,
     :     -2.,  0.,  2.,  0.,  2.,      -3.,     0.,      1.,   0.,
     :      1., -1.,  0., -1.,  0.,      -3.,     0.,      0.,   0.,
     :      0., -2.,  2., -2.,  1.,      -2.,     0.,      1.,   0.,
     :      2.,  0., -2.,  0.,  1.,       1.,     0.,      0.,   0.,
     :      0.,  0.,  2., -2.,  2.,  -13187.,   -16.,   5736., -31.,
     :      0.,  1.,  0.,  0.,  0.,    1426.,   -34.,     54.,  -1. /
      DATA ((X(I,J),I=1,9),J=11,20) /
     :      0.,  1.,  2., -2.,  2.,    -517.,    12.,    224.,  -6.,
     :      0., -1.,  2., -2.,  2.,     217.,    -5.,    -95.,   3.,
     :      0.,  0.,  2., -2.,  1.,     129.,     1.,    -70.,   0.,
     :      2.,  0.,  0., -2.,  0.,      48.,     0.,      1.,   0.,
     :      0.,  0.,  2., -2.,  0.,     -22.,     0.,      0.,   0.,
     :      0.,  2.,  0.,  0.,  0.,      17.,    -1.,      0.,   0.,
     :      0.,  1.,  0.,  0.,  1.,     -15.,     0.,      9.,   0.,
     :      0.,  2.,  2., -2.,  2.,     -16.,     1.,      7.,   0.,
     :      0., -1.,  0.,  0.,  1.,     -12.,     0.,      6.,   0.,
     :     -2.,  0.,  0.,  2.,  1.,      -6.,     0.,      3.,   0. /
      DATA ((X(I,J),I=1,9),J=21,30) /
     :      0., -1.,  2., -2.,  1.,      -5.,     0.,      3.,   0.,
     :      2.,  0.,  0., -2.,  1.,       4.,     0.,     -2.,   0.,
     :      0.,  1.,  2., -2.,  1.,       4.,     0.,     -2.,   0.,
     :      1.,  0.,  0., -1.,  0.,      -4.,     0.,      0.,   0.,
     :      2.,  1.,  0., -2.,  0.,       1.,     0.,      0.,   0.,
     :      0.,  0., -2.,  2.,  1.,       1.,     0.,      0.,   0.,
     :      0.,  1., -2.,  2.,  0.,      -1.,     0.,      0.,   0.,
     :      0.,  1.,  0.,  0.,  2.,       1.,     0.,      0.,   0.,
     :     -1.,  0.,  0.,  1.,  1.,       1.,     0.,      0.,   0.,
     :      0.,  1.,  2., -2.,  0.,      -1.,     0.,      0.,   0. /
      DATA ((X(I,J),I=1,9),J=31,40) /
     :      0.,  0.,  2.,  0.,  2.,   -2274.,    -2.,    977.,  -5.,
     :      1.,  0.,  0.,  0.,  0.,     712.,     1.,     -7.,   0.,
     :      0.,  0.,  2.,  0.,  1.,    -386.,    -4.,    200.,   0.,
     :      1.,  0.,  2.,  0.,  2.,    -301.,     0.,    129.,  -1.,
     :      1.,  0.,  0., -2.,  0.,    -158.,     0.,     -1.,   0.,
     :     -1.,  0.,  2.,  0.,  2.,     123.,     0.,    -53.,   0.,
     :      0.,  0.,  0.,  2.,  0.,      63.,     0.,     -2.,   0.,
     :      1.,  0.,  0.,  0.,  1.,      63.,     1.,    -33.,   0.,
     :     -1.,  0.,  0.,  0.,  1.,     -58.,    -1.,     32.,   0.,
     :     -1.,  0.,  2.,  2.,  2.,     -59.,     0.,     26.,   0./
      DATA ((X(I,J),I=1,9),J=41,50) /
     :      1.,  0.,  2.,  0.,  1.,     -51.,     0.,     27.,   0.,
     :      0.,  0.,  2.,  2.,  2.,     -38.,     0.,     16.,   0.,
     :      2.,  0.,  0.,  0.,  0.,      29.,     0.,     -1.,   0.,
     :      1.,  0.,  2., -2.,  2.,      29.,     0.,    -12.,   0.,
     :      2.,  0.,  2.,  0.,  2.,     -31.,     0.,     13.,   0.,
     :      0.,  0.,  2.,  0.,  0.,      26.,     0.,     -1.,   0.,
     :     -1.,  0.,  2.,  0.,  1.,      21.,     0.,    -10.,   0.,
     :     -1.,  0.,  0.,  2.,  1.,      16.,     0.,     -8.,   0.,
     :      1.,  0.,  0., -2.,  1.,     -13.,     0.,      7.,   0.,
     :     -1.,  0.,  2.,  2.,  1.,     -10.,     0.,      5.,   0. /
      DATA ((X(I,J),I=1,9),J=51,60) /
     :      1.,  1.,  0., -2.,  0.,      -7.,     0.,      0.,   0.,
     :      0.,  1.,  2.,  0.,  2.,       7.,     0.,     -3.,   0.,
     :      0., -1.,  2.,  0.,  2.,      -7.,     0.,      3.,   0.,
     :      1.,  0.,  2.,  2.,  2.,      -8.,     0.,      3.,   0.,
     :      1.,  0.,  0.,  2.,  0.,       6.,     0.,      0.,   0.,
     :      2.,  0.,  2., -2.,  2.,       6.,     0.,     -3.,   0.,
     :      0.,  0.,  0.,  2.,  1.,      -6.,     0.,      3.,   0.,
     :      0.,  0.,  2.,  2.,  1.,      -7.,     0.,      3.,   0.,
     :      1.,  0.,  2., -2.,  1.,       6.,     0.,     -3.,   0.,
     :      0.,  0.,  0., -2.,  1.,      -5.,     0.,      3.,   0. /
      DATA ((X(I,J),I=1,9),J=61,70) /
     :      1., -1.,  0.,  0.,  0.,       5.,     0.,      0.,   0.,
     :      2.,  0.,  2.,  0.,  1.,      -5.,     0.,      3.,   0.,
     :      0.,  1.,  0., -2.,  0.,      -4.,     0.,      0.,   0.,
     :      1.,  0., -2.,  0.,  0.,       4.,     0.,      0.,   0.,
     :      0.,  0.,  0.,  1.,  0.,      -4.,     0.,      0.,   0.,
     :      1.,  1.,  0.,  0.,  0.,      -3.,     0.,      0.,   0.,
     :      1.,  0.,  2.,  0.,  0.,       3.,     0.,      0.,   0.,
     :      1., -1.,  2.,  0.,  2.,      -3.,     0.,      1.,   0.,
     :     -1., -1.,  2.,  2.,  2.,      -3.,     0.,      1.,   0.,
     :     -2.,  0.,  0.,  0.,  1.,      -2.,     0.,      1.,   0. /
      DATA ((X(I,J),I=1,9),J=71,80) /
     :      3.,  0.,  2.,  0.,  2.,      -3.,     0.,      1.,   0.,
     :      0., -1.,  2.,  2.,  2.,      -3.,     0.,      1.,   0.,
     :      1.,  1.,  2.,  0.,  2.,       2.,     0.,     -1.,   0.,
     :     -1.,  0.,  2., -2.,  1.,      -2.,     0.,      1.,   0.,
     :      2.,  0.,  0.,  0.,  1.,       2.,     0.,     -1.,   0.,
     :      1.,  0.,  0.,  0.,  2.,      -2.,     0.,      1.,   0.,
     :      3.,  0.,  0.,  0.,  0.,       2.,     0.,      0.,   0.,
     :      0.,  0.,  2.,  1.,  2.,       2.,     0.,     -1.,   0.,
     :     -1.,  0.,  0.,  0.,  2.,       1.,     0.,     -1.,   0.,
     :      1.,  0.,  0., -4.,  0.,      -1.,     0.,      0.,   0. /
      DATA ((X(I,J),I=1,9),J=81,90) /
     :     -2.,  0.,  2.,  2.,  2.,       1.,     0.,     -1.,   0.,
     :     -1.,  0.,  2.,  4.,  2.,      -2.,     0.,      1.,   0.,
     :      2.,  0.,  0., -4.,  0.,      -1.,     0.,      0.,   0.,
     :      1.,  1.,  2., -2.,  2.,       1.,     0.,     -1.,   0.,
     :      1.,  0.,  2.,  2.,  1.,      -1.,     0.,      1.,   0.,
     :     -2.,  0.,  2.,  4.,  2.,      -1.,     0.,      1.,   0.,
     :     -1.,  0.,  4.,  0.,  2.,       1.,     0.,      0.,   0.,
     :      1., -1.,  0., -2.,  0.,       1.,     0.,      0.,   0.,
     :      2.,  0.,  2., -2.,  1.,       1.,     0.,     -1.,   0.,
     :      2.,  0.,  2.,  2.,  2.,      -1.,     0.,      0.,   0. /
      DATA ((X(I,J),I=1,9),J=91,100) /
     :      1.,  0.,  0.,  2.,  1.,      -1.,     0.,      0.,   0.,
     :      0.,  0.,  4., -2.,  2.,       1.,     0.,      0.,   0.,
     :      3.,  0.,  2., -2.,  2.,       1.,     0.,      0.,   0.,
     :      1.,  0.,  2., -2.,  0.,      -1.,     0.,      0.,   0.,
     :      0.,  1.,  2.,  0.,  1.,       1.,     0.,      0.,   0.,
     :     -1., -1.,  0.,  2.,  1.,       1.,     0.,      0.,   0.,
     :      0.,  0., -2.,  0.,  1.,      -1.,     0.,      0.,   0.,
     :      0.,  0.,  2., -1.,  2.,      -1.,     0.,      0.,   0.,
     :      0.,  1.,  0.,  2.,  0.,      -1.,     0.,      0.,   0.,
     :      1.,  0., -2., -2.,  0.,      -1.,     0.,      0.,   0. /
      DATA ((X(I,J),I=1,9),J=101,106) /
     :      0., -1.,  2.,  0.,  1.,      -1.,     0.,      0.,   0.,
     :      1.,  1.,  0., -2.,  1.,      -1.,     0.,      0.,   0.,
     :      1.,  0., -2.,  2.,  0.,      -1.,     0.,      0.,   0.,
     :      2.,  0.,  0.,  2.,  0.,       1.,     0.,      0.,   0.,
     :      0.,  0.,  2.,  4.,  2.,      -1.,     0.,      0.,   0.,
     :      0.,  1.,  0.,  1.,  0.,       1.,     0.,      0.,   0. /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*
*  FUNDAMENTAL ARGUMENTS in the FK5 reference system
*

*  Mean longitude of the Moon minus mean longitude of the Moon's
*  perigee.
      EL = ANPM ( ( 485866.733D0 + ( 715922.633D0 +
     :                  ( 31.310D0 + 0.064D0 * T ) * T ) * T ) * DAS2R
     :                + MOD(1325D0*T, 1D0) * D2PI )

*  Mean longitude of the Sun minus mean longitude of the Sun's perigee.
      ELP = ANPM ( ( 1287099.804D0 + ( 1292581.224D0 +
     :                   ( -0.577D0 -0.012D0 * T ) * T ) * T ) * DAS2R
     :                 + MOD(99D0*T, 1D0) * D2PI )

*  Mean longitude of the Moon minus mean longitude of the Moon's node.
      F = ANPM ( ( 335778.877D0 + ( 295263.137D0 +
     :                 ( -13.257D0 + 0.011D0 * T ) * T ) * T ) * DAS2R
     :               + MOD(1342D0*T, 1D0) * D2PI )

*  Mean elongation of the Moon from the Sun.
      D = ANPM ( ( 1072261.307D0 + ( 1105601.328D0 +
     :                 ( -6.891D0 + 0.019D0 * T ) * T ) * T ) * DAS2R
     :               + MOD(1236D0*T, 1D0) * D2PI )

*  Longitude of the mean ascending node of the lunar orbit on the
*  ecliptic, measured from the mean equinox of date.
      OM = ANPM( ( 450160.280D0 + ( -482890.539D0 +
     :                 ( 7.455D0 + 0.008D0 * T ) * T ) * T ) * DAS2R
     :               + MOD( -5D0*T, 1D0) * D2PI )

*  ---------------
*  Nutation series
*  ---------------

*  Change time argument from centuries to millennia.
      T = T / 10D0

*  Initialize nutation components.
      DP = 0D0
      DE = 0D0

*  Sum the nutation terms, ending with the biggest.
      DO 1 J=106,1,-1

*     Form argument for current term.
         ARG = DBLE(X(1,J)) * EL
     :       + DBLE(X(2,J)) * ELP
     :       + DBLE(X(3,J)) * F
     :       + DBLE(X(4,J)) * D
     :       + DBLE(X(5,J)) * OM

*     Accumulate current nutation term.
         S = DBLE(X(6,J)) + DBLE(X(7,J)) * T
         C = DBLE(X(8,J)) + DBLE(X(9,J)) * T
         IF ( S .NE. 0D0 ) DP = DP + S * SIN(ARG)
         IF ( C .NE. 0D0 ) DE = DE + C * COS(ARG)

*     Next term.
 1    CONTINUE

*  Convert results from 0.1 mas units to radians.
      DPSI = DP * U2R
      DEPS = DE * U2R

      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Function ANPM
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C     Normalize angle into the range -pi <= A < +pi.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 FUNCTION ANPM ( A )
      IMPLICIT NONE

      REAL*8 A

*  Pi
      REAL*8 DPI
      PARAMETER ( DPI = 3.141592653589793238462643D0 )

*  2Pi
      REAL*8 D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      REAL*8 W

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      W = MOD(A,D2PI)
      IF ( ABS(W) .GE. DPI ) W = W - SIGN(D2PI,A)
      ANPM = W

      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Function OBL80
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C     Mean obliquity of the ecliptic, IAU 1980 model.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 FUNCTION OBL80 ( DATE1, DATE2 )

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2

*  Arcseconds to radians
      REAL*8 DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000), JD
      REAL*8 DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      REAL*8 DJC
      PARAMETER ( DJC = 36525D0 )

      REAL*8 T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  Mean obliquity of date.
      OBL80 = DAS2R * ( 84381.448D0 +
     :                      ( -46.8150D0 +
     :                       ( -0.00059D0 +
     :                          0.001813D0 * T ) * T ) * T )

      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Subroutine NUMAT
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The only change to original routine is:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C
C     Form the matrix of nutation.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE NUMAT ( EPSA, DPSI, DEPS, RMATN )

      IMPLICIT NONE

      REAL*8 EPSA, DPSI, DEPS, RMATN(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Build the rotation matrix.
      CALL IR ( RMATN )
      CALL RX ( EPSA, RMATN )
      CALL RZ ( -DPSI, RMATN )
      CALL RX ( -(EPSA+DEPS), RMATN )

      end


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Function EQEQ94
C      This routine has been extracted from the SOFA software
C      (Standards of Fundamental Astronomy, http://www.iau-sofa.rl.ac.uk/ ).
C      by S. Bourdarie in July 2009.Therefore it does not itself
*        constitute software provided by and/or endorsed by SOFA.
C
C      The changes to original routine are:
C     DOUBLE PRECISION changed to REAL*8 by S. Bourdarie
C     Routine ANP has been included in this one by S. Bourdarie
C
C  Equation of the equinoxes, IAU 1994 model.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 FUNCTION EQEQ94 ( DATE1, DATE2 )

      IMPLICIT NONE

      REAL*8 DATE1, DATE2

*  Arcseconds to radians
      REAL*8 DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  2Pi
      REAL*8 D2PI
      PARAMETER (D2PI = 6.283185307179586476925287D0)

*  Reference epoch (J2000), JD
      REAL*8 DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      REAL*8 DJC
      PARAMETER ( DJC = 36525D0 )

      REAL*8 T, OM, DPSI, DEPS, EPS0
      REAL*8 ANPM, OBL80
      REAL*8 W

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  Longitude of the mean ascending node of the lunar orbit on the
*  ecliptic, measured from the mean equinox of date.
      OM = ANPM( ( 450160.280D0 + ( -482890.539D0 +
     :                 ( 7.455D0 + 0.008D0 * T ) * T ) * T ) * DAS2R
     :               + MOD(-5D0*T,1D0) * D2PI )

*  Nutation components and mean obliquity.
      CALL NUT80 ( DATE1, DATE2, DPSI, DEPS )
      EPS0 = OBL80 ( DATE1, DATE2 )

*  Equation of the equinoxes.
      EQEQ94 = DPSI * COS(EPS0) + DAS2R * ( 0.00264D0 * SIN(OM) +
     :                                          0.000063D0 * SIN(OM+OM))
* Normalize angle into the range 0 <= A < 2pi.
      W = MOD(EQEQ94,D2PI)
      IF ( W .LT. 0.D0 ) W = W + D2PI
      EQEQ94 = W

      end

!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4....
!
! CREATION: S. Bourdarie - July 2009
! MODIFICATION: None
!
! DESCRIPTION: Compute the dates following the J2000 method.
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!
! OUTPUT: datE1, date2 -> julday
!
! CALLING SEQUENCE: call date2dateJ2000(iyr,idoy,UT,date1,date2)
!---------------------------------------------------------------------------------------------------
      SUBROUTINE date2dateJ2000(iyr,idoy,UT,date1,date2)
      IMPLICIT NONE
c
c     input variables
      INTEGER*4 iyr,idoy
      REAL*8 UT
c
c     internal variables
      REAL*8 DJ00
      PARAMETER ( DJ00 = 2451545D0 )
      INTEGER*4 Month, Day,hour, minute, second,julday,jd
      REAL*8 eps,eps2
c
c     Output variables
      REAL*8 DATE1,DATE2

c -------------------------------------------------------------------
c     Compute the dates following J2000 method
      DATE1=DJ00   ! julday for 1st january 2000 at 12 UT (J2000 method)
      call DOY_AND_UT2DATE_AND_TIME(iyr,idoy,UT, Month, Day,
     &            hour, minute, second)

      jd=JULDAY(iyr,month,day)

      eps=2.2204460d-16
      eps2 = eps*ABS(jd)
      if (eps2 .gt. eps) eps=eps2
      DATE2 = dble(jd) + ( (Hour/24.d0 - 0.5d0) +
     &        (UT-hour*3600.d0)/86400.d0 + eps )-DATE1  !(J2000 method)
c
      end


!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4....
!
! CREATION: S. Bourdarie - July 2009
! MODIFICATION: None
!
! DESCRIPTION: Rotates coordinates from ECI J2000 to ECI true of date
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xIN -> position in ECI J2000 system (double array(3))
!
! OUTPUT: xOUT -> position in ECI True of Date (TOD) system (double array(3))
!
! CALLING SEQUENCE: call ECIJ2000_2_ECITOD(iyr,idoy,UT,xIN,xOUT)
!---------------------------------------------------------------------------------------------------
      SUBROUTINE ECIJ2000_2_ECITOD(iyr,idoy,UT,xIN,xOUT)
c
      IMPLICIT NONE
c
c     input variables
      INTEGER*4 iyr,idoy
      REAL*8 UT,xIN(3)
c
c     internal variables
      REAL*8 RMATPN(3,3)
      REAL*8 DATE1,DATE2

c     Output variables
      REAL*8 xOUT(3)

c     Compute the dates following J2000 method
      call date2dateJ2000(iyr,idoy,UT,date1,date2)
c  Form the matrix of precession/nutation for a given date, IAU 1976
c  precession model, IAU 1980 nutation model.
      call PNM80 ( DATE1, DATE2, RMATPN)


c  Rotate input coordinates to outputs
      xOUT(1)=RMATPN(1,1)*xIN(1)+RMATPN(1,2)*xIN(2)+RMATPN(1,3)*xIN(3)
      xOUT(2)=RMATPN(2,1)*xIN(1)+RMATPN(2,2)*xIN(2)+RMATPN(2,3)*xIN(3)
      xOUT(3)=RMATPN(3,1)*xIN(1)+RMATPN(3,2)*xIN(2)+RMATPN(3,3)*xIN(3)
      end

!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4....
!
! CREATION: S. Bourdarie - July 2009
! MODIFICATION: None
!
! DESCRIPTION: Rotates coordinates from ECI true of date to ECI J2000
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xIN -> position in ECI True of Date (TOD) system (double array(3))
!
! OUTPUT: xOUT -> position in ECI J2000 system (double array(3))
!
! CALLING SEQUENCE: call ECITOD_2_ECIJ2000(iyr,idoy,UT,xIN,xOUT)
!---------------------------------------------------------------------------------------------------
      SUBROUTINE ECITOD_2_ECIJ2000(iyr,idoy,UT,xIN,xOUT)
c
      IMPLICIT NONE
c
c     input variables
      INTEGER*4 iyr,idoy
      REAL*8 UT,xIN(3)
c
c     internal variables
      REAL*8 RMATPN(3,3),invRMATPN(3,3),det
      REAL*8 DATE1,DATE2
c
c     Output variables
      REAL*8 xOUT(3)

c     Compute the dates following J2000 method
      call date2dateJ2000(iyr,idoy,UT,date1,date2)
c  Form the matrix of precession/nutation for a given date, IAU 1976
c  precession model, IAU 1980 nutation model.
      call PNM80 ( DATE1, DATE2, RMATPN)

c
c      invert RMATPN matrix
      det=RMATPN(1,1)*(RMATPN(2,2)*RMATPN(3,3)-RMATPN(2,3)*RMATPN(3,2))
     &   +RMATPN(1,2)*(RMATPN(2,1)*RMATPN(3,3)-RMATPN(2,3)*RMATPN(3,1))
     &   +RMATPN(1,3)*(RMATPN(2,1)*RMATPN(3,2)-RMATPN(2,2)*RMATPN(3,1))

      invRMATPN(1,1)=(RMATPN(2,2)*RMATPN(3,3)-RMATPN(2,3)*RMATPN(3,2))
     & /det
      invRMATPN(1,2)=(RMATPN(1,3)*RMATPN(3,2)-RMATPN(1,2)*RMATPN(3,3))
     & /det
      invRMATPN(1,3)=(RMATPN(1,2)*RMATPN(2,3)-RMATPN(1,3)*RMATPN(2,2))
     & /det
      invRMATPN(2,1)=(RMATPN(2,3)*RMATPN(3,1)-RMATPN(2,1)*RMATPN(3,3))
     & /det
      invRMATPN(2,2)=(RMATPN(1,1)*RMATPN(3,3)-RMATPN(1,3)*RMATPN(3,1))
     & /det
      invRMATPN(2,3)=(RMATPN(1,3)*RMATPN(2,1)-RMATPN(1,1)*RMATPN(2,3))
     & /det
      invRMATPN(3,1)=(RMATPN(2,1)*RMATPN(3,2)-RMATPN(2,2)*RMATPN(3,1))
     & /det
      invRMATPN(3,2)=(RMATPN(1,2)*RMATPN(3,1)-RMATPN(1,1)*RMATPN(3,2))
     & /det
      invRMATPN(3,3)=(RMATPN(1,1)*RMATPN(2,2)-RMATPN(1,2)*RMATPN(2,1))
     & /det
c
c  Rotate input coordinates to outputs
      xOUT(1)=invRMATPN(1,1)*xIN(1)+invRMATPN(1,2)*xIN(2)+
     &        invRMATPN(1,3)*xIN(3)
      xOUT(2)=invRMATPN(2,1)*xIN(1)+invRMATPN(2,2)*xIN(2)+
     &        invRMATPN(2,3)*xIN(3)
      xOUT(3)=invRMATPN(3,1)*xIN(1)+invRMATPN(3,2)*xIN(2)+
     &        invRMATPN(3,3)*xIN(3)
      end

!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4....
!
! CREATION: S. Bourdarie - July 2009
! MODIFICATION: None
!
! DESCRIPTION: Rotates coordinates from TEME to ECI true of date
! TEME=True Equator Mean Equinox
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xIN -> position in TEME system (double array(3))
!
! OUTPUT: xOUT -> position in ECI TOD system (double array(3))
!
! CALLING SEQUENCE: call TEME2ECITOD(iyr,idoy,UT,xIN,xOUT)
!---------------------------------------------------------------------------------------------------
      SUBROUTINE TEME2ECITOD(iyr,idoy,UT,xIN,xOUT)
c
      IMPLICIT NONE
c
c     input variables
      INTEGER*4 iyr,idoy
      REAL*8 UT,xIN(3)

c     internal variables
      REAL*8 RMATPN(3,3),invRMATPN(3,3),det
      REAL*8 DATE1,DATE2
      REAL*8 EQEQ94,psi,S,C,A(3,3)
c
c     Output variables
      REAL*8 xOUT(3)
c
c     Compute the dates following J2000 method
      call date2dateJ2000(iyr,idoy,UT,date1,date2)
!c  Form the rotation matrix around z axis by angle -Eqe82
      psi=-EQEQ94 ( DATE1, DATE2 )
*  Matrix representing new rotation.
      S = SIN(PSI)
      C = COS(PSI)
      CALL IR ( A )
      A(1,1) = C
      A(2,1) = -S
      A(1,2) = S
      A(2,2) = C
      xOUT(1)=A(1,1)*xIN(1)+A(1,2)*xIN(2)+A(1,3)*xIN(3)
      xOUT(2)=A(2,1)*xIN(1)+A(2,2)*xIN(2)+A(2,3)*xIN(3)
      xOUT(3)=A(3,1)*xIN(1)+A(3,2)*xIN(2)+A(3,3)*xIN(3)
      end
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4....
!
! CREATION: S. Bourdarie - July 2009
! MODIFICATION: None
!
! DESCRIPTION: Rotates coordinates from ECI true of date to TEME
! TEME=True Equator Mean Equinox
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xIN -> position in ECI True of Date (TOD) system (double array(3))
!
! OUTPUT: xOUT -> position in TEME system (double array(3))
!
! CALLING SEQUENCE: call ECITOD2TEME(iyr,idoy,UT,xIN,xOUT)
!---------------------------------------------------------------------------------------------------
      SUBROUTINE ECITOD2TEME(iyr,idoy,UT,xIN,xOUT)
c
      IMPLICIT NONE
c
c     input variables
      INTEGER*4 iyr,idoy
      REAL*8 UT,xIN(3)

c     internal variables
      REAL*8 RMATPN(3,3),invRMATPN(3,3),det
      REAL*8 DATE1,DATE2
      REAL*8 EQEQ94,psi,S,C,A(3,3)
c
c     Output variables
      REAL*8 xOUT(3)
c
c     Compute the dates following J2000 method
      call date2dateJ2000(iyr,idoy,UT,date1,date2)
!c  Form the rotation matrix around z axis by angle Eqe82
      psi=EQEQ94 ( DATE1, DATE2 )
*  Matrix representing new rotation.
      S = SIN(PSI)
      C = COS(PSI)
      CALL IR ( A )
      A(1,1) = C
      A(2,1) = -S
      A(1,2) = S
      A(2,2) = C
      xOUT(1)=A(1,1)*xIN(1)+A(1,2)*xIN(2)+A(1,3)*xIN(3)
      xOUT(2)=A(2,1)*xIN(1)+A(2,2)*xIN(2)+A(2,3)*xIN(3)
      xOUT(3)=A(3,1)*xIN(1)+A(3,2)*xIN(2)+A(3,3)*xIN(3)
      end
c
c---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.0
!
! CREATION: S. Bourdarie - January 2007
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from one Earth or Heliospheric coordinate
!                                        system to another one
!
! INPUT: sysaxesIN -> designed input coordinate system (long integer)
!        sysaxesOUT> designed output coordinate system (long integer)
!        iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xIN -> position in input coordinate system (double array(3))
!
! OUTPUT: xOUT -> position in output coordinate system (double array(3))
!
! CALLING SEQUENCE: call coord_trans1(sysaxesIN,sysaxesOUT,iyr,idoy,secs,xIN,xOUT)
!---------------------------------------------------------------------------------------------------
c
      SUBROUTINE coord_trans1(sysaxesIN,sysaxesOUT,iyr,idoy,
     &   secs,xIN,xOUT)
c
      IMPLICIT NONE
c
      INTEGER*4 sysaxesIN,sysaxesOUT,iyr,idoy
      INTEGER*4 i
      REAL*8    secs,psi
      REAL*8    xIN(3),xOUT(3),xTMP(3),alti

      call initize ! sets rad, pi used by various routines

      if (sysaxesIN.EQ.sysaxesOUT) then
c         write(6,*)'sysaxesIN = sysaxesOUT!'
         do i=1,3
            xOUT(i)=xIN(i)
         enddo
         return
      endif
      if ((sysaxesIN.LT.0).or.(sysaxesIN.GT.14)) then
         write(6,*)'sysaxesIN out of range !'
         do i=1,3
            xOUT(i)=-1.D31
         enddo
         return
      endif
      if ((sysaxesOUT.LT.0).or.(sysaxesOUT.GT.14)) then
         write(6,*)'sysaxesOUT out of range !'
         do i=1,3
            xOUT(i)=-1.D31
         enddo
         return
      endif
c
c input=GDZ
      if (sysaxesIN.EQ.0) then
         if (sysaxesOUT.EQ.1) then  !GEO
            call gdz_geo(xIN(2),xIN(3),xIN(1),xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            xOUT(1)=SQRT(xTMP(1)*xTMP(1)+xTMP(2)*xTMP(2)+
     &               xTMP(3)*xTMP(3))
            xOUT(2)=xIN(2)
            xOUT(3)=xIN(3)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=GEO
      if (sysaxesIN.EQ.1) then
         if (sysaxesOUT.EQ.0) then  !GDZ
           call geo_gdz(xIN(1),xIN(2),xIN(3),xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call geo2gsm1(iyr,idoy,secs,psi,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call geo2gse1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call geo2sm1(iyr,idoy,secs,xIN,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call geo2gei1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call geo2mag1(iyr,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call CAR_SPH(xIN,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
           call geo_gdz(xIN(1),xIN(2),xIN(3),xTMP(1),xTMP(2),xTMP(3))
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call geo2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call geo2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call geo2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call geo2gei1(iyr,idoy,secs,xIN,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call geo2gei1(iyr,idoy,secs,xIN,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=GSM
      if (sysaxesIN.EQ.2) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call gsm2sm1(iyr,idoy,secs,xIN,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=GSE
      if (sysaxesIN.EQ.3) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call gse2geo1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call gse2hee1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call gse2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call gse2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call gse2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=SM
      if (sysaxesIN.EQ.4) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call sm2geo1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call sm2gsm1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call sm2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c input=GEI (which is also ECI ToD)
      if ((sysaxesIN.EQ.5).or.(sysaxesIN.EQ.12)) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call gei2geo1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call gei2geo1(iyr,idoy,secs,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.12).or.(sysaxesOUT.EQ.5)) then  !ECI ToD
            do i=1,3
               xOUT(i)=xIN(i)
            enddo
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            CALL ECITOD2TEME(iyr,idoy,secs,xIN,xOUT)
         endif
      endif
c
c
c input=MAG
      if (sysaxesIN.EQ.6) then
         if (sysaxesOUT.EQ.0) then  !GDZ
              call mag2geo1(iyr,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
              call mag2geo1(iyr,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
              call mag2geo1(iyr,xIN,xTMP)
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
              call mag2geo1(iyr,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
              call mag2geo1(iyr,xIN,xTMP)
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
              call mag2geo1(iyr,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
              call mag2geo1(iyr,xIN,xTMP)
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
              call mag2geo1(iyr,xIN,xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
              call mag2geo1(iyr,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
              call mag2geo1(iyr,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
              call mag2geo1(iyr,xIN,xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call mag2geo1(iyr,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call mag2geo1(iyr,xIN,xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=SPH
      if (sysaxesIN.EQ.7) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3)
     &             ,xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=xIN(1)
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=RLL
      if (sysaxesIN.EQ.8) then
         if (sysaxesOUT.EQ.0) then  !GDZ
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
            call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=HEE
      if (sysaxesIN.EQ.9) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            xOUT(3)=xOUT(2)
            xOUT(2)=xOUT(1)
            xOUT(1)=SQRT(xTMP(1)*xTMP(1)+xTMP(2)*xTMP(2)
     &                     +xTMP(3)*xTMP(3))
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call hee2hae1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call hee2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call hee2gse1(iyr,idoy,secs,xIN,xTMP)
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=HAE
      if (sysaxesIN.EQ.10) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            xOUT(3)=xOUT(2)
            xOUT(2)=xOUT(1)
            xOUT(1)=SQRT(xTMP(1)*xTMP(1)+xTMP(2)*xTMP(2)
     &                     +xTMP(3)*xTMP(3))
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call hae2hee1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call hae2heeq1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call hae2hee1(iyr,idoy,secs,xIN,xTMP)
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=HEEQ
      if (sysaxesIN.EQ.11) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            xOUT(3)=xOUT(2)
            xOUT(2)=xOUT(1)
            xOUT(1)=SQRT(xTMP(1)*xTMP(1)+xTMP(2)*xTMP(2)
     &                     +xTMP(3)*xTMP(3))
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call heeq2hae1(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call heeq2hae1(iyr,idoy,secs,xIN,xTMP)
            call hae2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            CALL ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=J2000
      if (sysaxesIN.EQ.13) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.14) then  !TEME
            call ECIJ2000_2_ECITOD(iyr,idoy,secs,xIN,xTMP)
            call ECITOD2TEME(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c input=TEME
      if (sysaxesIN.EQ.14) then
         if (sysaxesOUT.EQ.0) then  !GDZ
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
         endif
         if (sysaxesOUT.EQ.1) then  !GEO
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.2) then  !GSM
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.3) then  !GSE
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.4) then  !SM
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if ((sysaxesOUT.EQ.5).or.(sysaxesOUT.EQ.12)) then  !GEI (TOD)
            call TEME2ECITOD(iyr,idoy,secs,xIN,xOUT)
         endif
         if (sysaxesOUT.EQ.6) then  !MAG
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2mag1(iyr,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.7) then  !SPH
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
         endif
         if (sysaxesOUT.EQ.8) then  !RLL
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
            xOUT(2)=xTMP(1)
            xOUT(3)=xTMP(2)
         endif
         if (sysaxesOUT.EQ.9) then  !HEE
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.10) then  !HAE
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.11) then  !HEEQ
            call TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            call gei2geo1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
            do i=1,3
               xTMP(i)=xOUT(i)
            enddo
            call hae2heeq1(iyr,idoy,secs,xTMP,xOUT)
         endif
         if (sysaxesOUT.EQ.13) then  !J2000
            CALL TEME2ECITOD(iyr,idoy,secs,xIN,xTMP)
            CALL ECITOD_2_ECIJ2000(iyr,idoy,secs,xTMP,xOUT)
         endif
      endif
c
c
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE coord_trans_vec
C      Routine to vectorize coord_trans1, the single-element
C      coordinate transformation utility.
C      Calling sequence:
C
C      call coord_trans_vec(ntime,sysaxesIN,sysaxesOUT,iyear,idoy,secs,
C     +          xIN,xOUT)
C
C         with the following types:
C      INTEGER*4 sysaxesIN, sysaxesOUT
C      INTEGER*4 iyear(nmax),idoy(nmax)
C      REAL*8 secs(nmax)
C      REAL*8 xINV(3,nmax), xOUTV(3,nmax)
C      INTEGER*4 numpoints
C
C     As with all (most?) onera library calls, the maximum array size
C     is limited to ntime_max elements
C                              Contributed by Timothy Guild, 2.2.07
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

        SUBROUTINE coord_trans_vec1(ntime,sysaxesIN,sysaxesOUT,
     &   iyear,idoy,secs,xINV,xOUTV)

      IMPLICIT NONE
      INCLUDE 'ntime_max.inc'
      INCLUDE 'variables.inc'

      INTEGER*4 nmax,i,ntime, sysaxesIN, sysaxesOUT
      INTEGER*4 iyear(ntime_max),idoy(ntime_max),y,d
      REAL*8 secs(ntime_max),xINV(3,ntime_max),xOUTV(3,ntime_max)
      ! local vars
      REAL*8 xIN(3),xOUT(3),s

      ! Loop over the number of points specified, calling
      !  coord_trans1 each iteration
      DO i = 1,ntime
         y = iyear(i)
         d = idoy(i)
         s = secs(i)

        if (xINV(1,i) .eq. baddata .and. xINV(2,i) .eq. baddata
     &  .and. xINV(3,i) .eq.baddata) then
           xOUTV(1,i) = baddata
           xOUTV(2,i) = baddata
           xOUTV(3,i) = baddata
           goto 10
        endif
       xIN(1) = xINV(1,i) ! copy each array element into 3x1 array
       xIN(2) = xINV(2,i)
         xIN(3) = xINV(3,i)

         call coord_trans1( sysaxesIN,sysaxesOUT,y,d,s,xIN,xOUT)

       xOUTV(1,i) = xOUT(1)  ! copy back to 3xN array to pass back
       xOUTV(2,i) = xOUT(2)
         xOUTV(3,i) = xOUT(3)

10      continue
      ENDDO

      END


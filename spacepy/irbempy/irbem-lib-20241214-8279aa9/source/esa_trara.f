!***************************************************************************************************
! Copyright 1987, D. Bilitza, 2009 H. Evans
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
!
! CREATION: H. Evans European Space Agency, 1993-2009
!
! The method of interpolation used by ESA_TRARAP is documented in:
!   E.J. Daly et al, "Problems with Models of the Radiation Belts", 
!   IEEE Trans. Nucl. Sci, Vol 43, No.2, April 1996
!
! FILE CONTENT:
C      SUBROUTINE ESA_TRARA1 : Finds particle fluxes for given energies, mag field strength and L-shell
C      SUBROUTINE ESA_CrossP : computes cross product of two vectors
C      SUBROUTINE ESA_EVals  : returns the energy levels of the A*-8 particle model
C      SUBROUTINE ESA_LVals  : returns the L-shell values in the A*-8 particle models
C      SUBROUTINE ESA_BB0Val : returns the B/B0 and Log10(flux) values from an L-shell string
C      SUBROUTINE ESA_BVals  : returns the B and Log10(flux) values from an L-shell string
C      REAL*8 FUNCTION ESA_TRARAP     : Finds particle flux for given B and L from a model energy string
C      REAL*8 FUNCTION ESA_INTERP     : linearly interpolates between 3 points in 3d space
C      REAL*8 FUNCTION ESA_DotP       : calculates the scalar product of two vectors
C      LOGICAL FUNCTION ESA_InPoly    : Returns TRUE/FALSE if XY coords of a point are within a polygon
C                                       described by 3 points
C      REAL*8 FUNCTION ESA_XINTER     : linearly interpolates between (x1,y1) and (x2,y2)
C      INTEGER*4 FUNCTION ESA_LStrng  : returns the index in the energy string of the 
C                                       retuired L value string
C      REAL*8 FUNCTION ESA_XPHI       : computes the hybrid magnetic field coordinate
C      REAL*8 FUNCTION ESA_B0         : computes B0 for an L value using McIlwain's magnetic moment.
!
!***************************************************************************************************
      SUBROUTINE ESA_TRARA1(DESCR, MAP, FL, BoB0, E, F)
C***********************************************************************
C     B / B0 CASE     JOEL STEIN     9-15-71     X2133     KMS
C   ESA_TRARA1 DOES ENERGY VALUE SEARCH FOR FLUX CALCULATION WHEN GIVEN A
C    B AND L POINT.
*  Modified version based on Kluge and Lenhart ESOC Int. Note 78 (1971)
*  
*  Modified    APRIL 1988 E.J. DALY ESA/ESTEC/WMA
*                         can pick up geomagnetic dipole moment GMAGMO
*                         from common block and use it in place of
*                         McIlwain's value of .311653
*                         GMAGMO is passed from the geomagnetic
*                         coordinate program.
*                         Default for models should be McIlwain's value
*
*  Modified    JULY  1993 H.D.R. EVANS    ESA/ESTEC/WMA
*                         Kluge and Lenhart interpolation method in
*                         (B/B0,L) space changed to a linear polygon
*                         interpolation method in (Phi,L) space.
*
C***********************************************************************
C MAP(1) is the first word of list

      IMPLICIT NONE

      LOGICAL S0, S1, S2
C S0,S1,S2 are logical variables which indicate whether the flux for a
C particular E,B,L point has already been found in a previous call
C to ESA_TRARAP.

      REAL*8 E, F, FL, BOB0
      REAL*8 E0, E1, E2, F0, F1, F2
      REAL*8 ESA_B0, ESA_TRARAP
      REAL*4 DESCR(8)
      INTEGER*4 MAP(*), N, NB, NL, I0, I1, I2, I3, L3, IE

C      PRINT*,'ESA_TRARA1:', FL, BoB0, E

      NL = MIN(32766.0D0, ABS(FL*DBLE(DESCR(5))))
C max B/B0 allowed here is 1000 (protect against integer overflow)
C      BOB0 = MIN( BABS/ ESA_B0(FL), 1000.0D0)

      IF (BOB0 .LT. 1.0D0) BOB0 = 1.0D0
        
C force any possible representational errors to 1.0
      BoB0 = MAX(BoB0, 1.0D0)

      NB = ABS((BOB0-1.0)*DESCR(6))
C NL is the minimum of the L value or 15.999, scaled to an integer by
C THE L SCALING FACTOR
C NB IS THE DIFFERENCE BETWEEN THE INPUT B VALUE AND B EQUATORIAL,
C SCALED TO AN INTEGER BY THE B SCALING FACTOR.
      I1 = 0
      I2 = MAP(1)
      I3 = I2 + MAP(I2+1)
      L3 = MAP(I3+1)
      E1 = MAP(I1+2) / DESCR(4)
      E2 = MAP(I2+2) / DESCR(4)
      S1 = .TRUE.
      S2 = .TRUE.

C I2 IS THE NUMBER OF ELEMENTS IN THE FLUX MAP FOR THE FIRST ENERGY.
C I3 IS THE INDEX OF THE LAST ELEMENT OF THE SECOND ENERGY MAP.
C L3 IS THE LENGTH OF THE MAP FOR THE THIRD ENERGY.
C E1 IS THE ENERGY OF THE FIRST ENERGY MAP (UNSCALED)
C E2 IS THE ENERGY OF THE SECOND ENERGY MAP (UNSCALED)
C S1 AND S2 ARE TRUE TO INDICATE THAT NO FLUXES HAVE YET BEEN FOUND.
C THE DO STATEMENT LOOPS THROUGH THE ENERGIES FOR WHICH FLUXES ARE
C DESIRED AT THE GIVEN B,L POINT (BABS,FL).
 1      IF ((E .LE. E2) .OR. (L3 .EQ. 0)) GO TO 2
C THE IF STATEMENT CHECKS TO SEE IF THE INPUT ENERGY IS LESS THAN OR E
C THE ENERGY OF THE SECOND MAP, OR IF THE LENGTH OF THE THIRD MAP IS 0.
C (I.E. THERE ARE NO HIGHER ENERGIES IN THE TABLE).  IF TRUE, USE TH
C FOR THOSE TWO ENERGY MAPS TO FIND THE DESIRED FLUX AT THE DESIRED
C ENERGY.  IF FALSE, THE ZEROTH ENERGY MAP IS DEFINED TO BE TNE FIRS
C ENERGY MAP, THE FIRST BECOMES THE SECOND, AND THE SECOND BECOMES
C THE THIRD.  E0,E1,E2 ARE THE ENERGIES FOR THE ZEROTH,FIRST,AND SEC
C ENERGY MAPS.  F0,F1,F2 ARE THE FLUXES FOR THE ZEROTH, FIRST, AND
C SECOND ENERGY MAPS AT THE B,L POINT.
        I0 = I1
        I1 = I2
        I2 = I3
        I3 = I3 + L3
        L3 = MAP(I3+1)
        E0 = E1
        E1 = E2
        E2 = MAP(I2+2) / DESCR(4)
        S0 = S1
        S1 = S2
        S2 = .TRUE.
        F0 = F1
        F1 = F2
        GO TO 1
 2      IF (S1) F1 = ESA_TRARAP(DESCR, MAP(I1+1), FL, BoB0)
        IF (S2) F2 = ESA_TRARAP(DESCR, MAP(I2+1), FL, BoB0)
C THESE TWO LOGICAL IFS CALL ESA_TRARAP FOR THE FLUX FROM THE FIRST AND
C SECOND ENERGY MAPS AT THE B,L POINT IF THEY HAVE NOT ALREADY BEEN
        S1 = .FALSE.
        S2 = .FALSE.
C S1 AND S2 ARE FALSE SINCE F1 AND F2 ARE NOW FOUND.
        F = F1 + (F2-F1) * (E-E1) / (E2-E1)
C INTERPOLATE FOR THE FLUX F(IE) USING THE FLUXES AND ENERGIES FOR MAP
C ONE AND TWO.
C THE FOLLOWING COMMENTS APPLY TO THE REMAINING PROGRAM STATEMENTS.
C IF THE FLUX F2 FOR THE SECOND ENERGY MAP IS GREATER THAN ZERO, OR TH
C ZEROTH ENERGY MAP HAS NOT BEEN DEFINED, THE FINAL FLUX IS THE MAXI
C OF THE INTEROOLATED FLUX OR ZERO.  IF THE FLUX FOR THE SECOND ENER
C MAP IS EQUAL TO ZERO, AND THE ZEROTH ENERGY MAP HAS BEEN DEFINED,
C THEN INTERPOLATE FOR THE FLUX USING THE ZEROTH AND FIRST ENERGY MA
C CHOOSE THE MINIMUM OF THE TWO INTERPOLATIONS, AND THEN THE MAXIMUM
C CHOICE AND ZERO FOR THE FINAL FLUX VALUE.
        IF ((F2 .LE. 0.0D0) .AND. (I1 .NE. 0)) THEN
          IF (S0) F0 = ESA_TRARAP(DESCR, MAP(I0+1), FL, BoB0)
          S0 = .FALSE.
          F = MIN(F, F0+(F1-F0)*(E-E0)/(E1-E0))
        END IF
        F = MAX(F, 0.0D0)
C      END DO
C      PRINT*,'ESA_TRARA1:', FL, BoB0, E,F

      RETURN
      END
C
C************************************************************************ 
C 
C    REAL*8 FUNCTION ESA_TRARAP 
C 
C     PURPOSE:    This function converts the B argument to a hybrid 
C                 value, defined by Phi = ASIN( ( B - B0) / (Bmax- B0)) 
C                 and interpolates the MAP in (PHI-L) space. 
C 
C    
C
C
C     METHOD:     The conversion to the Phi-L space requires the
C                 maximum B value. This is obtained by interpolating
C                 between the maximum B values for  the L strings
C                 that subtend the L value passed to the routine.
C
C                 The Flux at the Phi point along the two L strings
C                 is determined by  locating the three model points
C                 that form a polygon (triangle) that contains the
C                 (L,Phi) point.  The interpolation is then performed
C                 by a linear interpolation using the slope of the
C                 plane in (L,Phi,lnF) space.
C
C     HISTORY:    CREATED July 1993    H.D.R. Evans ESA/ESTEC/WMA
C
C************************************************************************

      REAL*8 FUNCTION ESA_TRARAP(HEADER, MAP, L, BoB0) 

      IMPLICIT NONE 

      REAL*4 HEADER(*)
      REAL*8 BoB0, L, B, Phi 
      INTEGER*4 MAP(*)
      REAL*8 Model(3,100,2), Ls(100), Bs(100), Fluxes(100), Bm(2), bo
      REAL*8 Verts(3,3), ReqPt(3)
      INTEGER*4 ESA_LStrng, Lpos(2), Llen, Blen(2), Li, I, J
      INTEGER*4 s(2), N(2), string, other, tmp
      LOGICAL Found, EndStr
      REAL*8 ESA_XPHI, ESA_B0, ESA_XINTER, ESA_INTERP
      LOGICAL ESA_INPOLY

      IF (BoB0 .LT. 1.0) THEN
        ESA_TRARAP = 0.0
        RETURN
      END IF

      B = BoB0 * ESA_B0(L)
      CALL ESA_LVals(MAP, header, Ls, Llen)

C Find the two L strings that subtend the L value we are given...
      DO Li=1,Llen-2
        IF (L .LE. Ls(Li+1)) GO TO 6
      END DO
 6    CONTINUE

C Get the Maximum B values for the two L strings.  Allows us to
C linearly interpolate to find the maximum B value for the given 
C L value.
      DO i=1,2
        Lpos(i) = ESA_LStrng(MAP, Header, Ls(Li+i-1), Llen) + 1
        CALL ESA_BVals(MAP(Lpos(i)), Header, Bs, Fluxes, Blen(i))
        Bm(i) = Bs(Blen(i))
        Bo = ESA_B0(Ls(li+i-1))
        IF (Bo .GT. 0.0) THEN
          DO j=1,Blen(i)
            Model(1,j,i) = Ls(Li+i-1)
            Model(2,j,i) = ESA_XPHI(Bs(j)/Bo, Bm(i)/Bo) 
            Model(3,j,i) = Fluxes(j)
          END DO
        ELSE
          Blen(i) = 1 
          Model(1,1,i) = Ls(Li+i-1)
          Model(2,1,i) = -1.0
          Model(3,1,i) = 0.0
        END IF
      END DO
      Bo = ESA_B0(L)
      Phi = ESA_XPHI(B/Bo, ESA_XINTER(L,Ls(Li),Bm(1),Ls(Li+1),Bm(2))/Bo)

C Check for an invalid Phi value, e.g. if B > Bmax
      IF (Phi .LT. 0.0) THEN
        ESA_TRARAP = 0.0
        RETURN
      END IF

C If length of both strings is 1, then linear interpolation between
C points and return
      IF ((blen(1) .EQ. 1) .AND. (blen(2) .EQ. 1)) THEN
        ESA_TRARAP = (L-Model(1,1,1)) * (Model(3,1,2)-Model(3,1,1)) /
     &           (Model(1,1,2) - Model(1,1,1))
        IF (ESA_TRARAP .GT. 0.0) ESA_TRARAP = (ESA_TRARAP)
        RETURN
      END IF

C Now find the H values in both L strings that subtend the required point.
      ReqPt(1) = L
      ReqPt(2) = Phi
      S(1) = 1
      S(2) = 1
      N(1) = MIN(S(1)+1, blen(1))
      N(2) = MIN(S(2)+1, blen(2))

 10   EndStr = (s(1) .EQ. blen(1)) .AND. (s(2) .EQ. blen(2))
      IF (.NOT. EndStr) THEN
        IF (S(1) .EQ. blen(1)) THEN
C String 1 is empty, have to use string 2 now.
          string = 2
          other = 1
        ELSE IF (S(2) .EQ. blen(2)) THEN
C String 2 is empty, have to use string 1.
          string = 1
          other = 2
        ELSE
          IF (Model(2,N(1),1) .LT. Model(2,N(2),2)) THEN
            string = 1
            other = 2
          ELSE
            string = 2
            other = 1
          END IF
        END IF
        Found = ESA_INPOLY(model(1,s(string),string),
     &                 model(1,s(other),other),
     &                 model(1,N(string),string), ReqPt)
        IF (.NOT. Found) THEN
            s(string) = N(string) 
            N(string) = MIN(N(string)+1, blen(String))
        END IF
      END IF

      IF (.NOT. (Found .OR. EndStr)) GO TO 10

C Repeat this UNTIL end of both strings or poly found
C Now check for end of string condition, requires backing up and
C using a previous point.
      IF (EndStr) THEN
        IF (blen(1) .EQ. 1) THEN
          string = 2
          other = 1
        ELSE IF (blen(2) .EQ. 1) THEN
          string = 1
          other = 2
        ELSE IF (Model(2,s(string)-1,string) .LT.
     &           Model(2,s(other)-1,other)) THEN
          tmp = other
          other = string
          string = tmp
        END IF
        n(string) = s(string) - 1
      END IF
      DO j=1,3
        verts(j,1) = Model(j,s(string),string)
        verts(j,2) = Model(j,s(other),other)
        verts(j,3) = Model(j,n(string),string)
      END DO
      ESA_TRARAP = ESA_INTERP(verts, ReqPt)
      IF (ESA_TRARAP .GT. 0) ESA_TRARAP = (ESA_TRARAP)

      RETURN
      END      
C
C***********************************************************************
C
C     REAL*8 FUNCTION ESA_INTERP
C
C     PURPOSE:    Interpolates between 3 points in 3D space.
C
C     METHOD:     Constructs function of a plane containing the 3 points
C                 and calculates the Z value for the given (X,Y) point.
C
C     HISTORY:    CREATED     July 1993      H.D.R. Evans
C
C***********************************************************************

      REAL*8 FUNCTION ESA_Interp(GPnts, ReqPnt)

      IMPLICIT NONE

C      INCLUDE 'interface.h'

      REAL*8 GPnts(3,*), ReqPnt(3)
      REAL*8 v1(3), v2(3), plane(4), ESA_DotP
      INTEGER*4 I, n /3/, incx /1/, incy /1/

      ESA_Interp = 0.0
C Compute vectors in plane of three points.
      DO i=1,3
        v1(i) = GPnts(i,2) - GPnts(i,1)
        v2(i) = GPnts(i,3) - GPnts(i,1)
      END DO
C Determine normal to plane defined by v1 and v2, and plane constant.
C plane(1)x + plane(2)y + plane(3)z = plane(4)
      call ESA_CrossP(v1, v2, plane, n)
      plane(4) = ESA_DotP(n, GPnts(1,1), incx, plane, incy)
C Value we require is the Z value at the point specified
C by the solution of:
C Z = (plane(4) - Plane(1)x - Plane(2)y ) / Plane(3)
      IF (Plane(3) .NE. 0.0) THEN
        ESA_Interp = (Plane(4)-Plane(1)*ReqPnt(1)-Plane(2)*ReqPnt(2)) /
     &           Plane(3)
      ELSE
C        CALL ERR_WRITE(' Plane with 3 given points is independent of Z')
      END IF

      RETURN
      END
C
C***********************************************************************
C
C     Subroutine ESA_CrossP
C
C     PURPOSE:    Takes the cross product of the two vectors.
C
C     METHOD:     Basic vector calculations.  Vectors must be the
C                 same size.
C
C     HISTORY:    CREATED     July 1993      H.D.R. Evans ESA/ESTEC/WMA
C
C***********************************************************************

      SUBROUTINE ESA_CrossP(X, Y, Z, DIM)

      IMPLICIT NONE

C X   - 1st vector
C Y   - 2nd vector
C Z   - crossproduct = X^Y
C DIM - dimension of X, Y and Z

      REAL*8 X(*), Y(*), Z(*), magZ
      INTEGER*4 DIM, I, J, INDX

      INDX(J) = MOD(J+DIM-1, DIM) + 1
      DO I=1,DIM
        magZ = 0.0
        Z(I) = X(INDX(I+1)) * Y(INDX(I+2)) - X(INDX(I+2)) * Y(INDX(I+1))
      END DO

      RETURN
      END
C
C***********************************************************************
C
C     REAL*8 FUNCTION ESA_DotP
C
C     PURPOSE:    Returns the inner (dot) product of SX and SY.
C
C     METHOD:     Basic vector calculations.  Vectors must be the
C                 same size.
C
C     HISTORY:    CREATED     July 1993      H.D.R. Evans ESA/ESTEC/WMA
C
C***********************************************************************

      REAL*8 FUNCTION ESA_DotP(N, sx, incx, sy, incy)

      IMPLICIT NONE

      REAL*8 sx(*), sy(*)
      INTEGER*4 N, incx, incy
      INTEGER*4 I, J, k, l, POS

      pos(i, j) = (i-1) * j + 1

      k = incx
      l = incy
      ESA_DotP = 0.0
      DO i=1,n
        ESA_DotP = ESA_DotP + sx(pos(i, k)) * sy(pos(i, l))
      END DO

      return
      end
C
C***********************************************************************
C
C     LOGICAL FUNCTION ESA_InPoly
C
C     PURPOSE:    Returns .TRUE. if the (X,Y) coordinates of point Pt 
C                 is in the polygon described by the points A,B,and C
C
C     METHOD:     Pts is decomposed into the sum of the line segments 
C                 AC and BC, i.e. PTS = a * AC + b * BC.  If either a 
C                 or b is less than zero, then PTS is not subtended by 
C                 the lines AC & BC.
C
C                 This is then repeated for the AB and CB line segements 
C                 to ensure PTS isn't the other side of the AB line from C.
C
C     HISTORY:    CREATED     July 1993      H.D.R. Evans ESA/ESTEC/WMA
C
C***********************************************************************

      LOGICAL FUNCTION ESA_InPoly(A, B, C, Pt)

      IMPLICIT NONE

C      INCLUDE 'interface.h'

      REAL*8 A(2), B(2), C(2), Pt(2)
      REAL*8 M(2,2), IM(2,2), PO(2), coeff(2), DET
      INTEGER*4 I

C Initialise the matrix
      DO i=1,2
        M(1,i) = A(i) - C(i)
        M(2,i) = B(i) - C(i)
        PO(i) = Pt(i) - C(i)
      END DO
      det = M(1,1) * M(2,2) - M(1,2) * M(2,1)
      IF (det .eq. 0.0) go to 999

      IM(1,1) = M(2,2)
      IM(2,2) = M(1,1)
      IM(1,2) = -M(2,1)
      IM(2,1) = -M(1,2)
      Coeff(1) = (IM(1,1)*PO(1)+IM(1,2)*PO(2)) / DET
      Coeff(2) = (IM(2,1)*PO(1)+IM(2,2)*PO(2)) / DET
      ESA_InPoly = (coeff(1) .GE. 0.0) .AND. (Coeff(2) .GE. 0.0)
C Repeat the previous, only this time with the end point.
      DO i=1,2
        M(1,i) = B(i) - A(i)
        M(2,i) = C(i) - A(i)
        PO(i) = Pt(i) - A(i)
      END DO
      det = M(1,1) * M(2,2) - M(1,2) * M(2,1)
      IF (det .eq. 0.0) go to 999
      IM(1,1) = M(2,2)
      IM(2,2) = M(1,1)
      IM(1,2) = -M(2,1)
      IM(2,1) = -M(1,2)
      Coeff(1) = (IM(1,1)*PO(1)+IM(1,2)*PO(2)) / DET
      Coeff(2) = (IM(2,1)*PO(1)+IM(2,2)*PO(2)) / DET
      ESA_InPoly = (coeff(1) .GE. 0.0) .AND. (Coeff(2) .GE. 0.0)
     &          .AND. ESA_INPOLY
      RETURN

 999  ESA_InPoly = .FALSE.
C      CALL ERR_WRITE(' ***INPOLY*** Determinant = 0')

      END
C
C***********************************************************************
C
C     REAL*8 FUNCTION ESA_XINTER
C
C     PURPOSE:    Linearly interpolates between (x1,y1) and (x2,y2).
C
C     METHOD:     simple linear interpolation.
C
C     HISTORY:    CREATED     July 1993      H.D.R. Evans ESA/ESTEC/WMA
C
C***********************************************************************

      REAL*8 FUNCTION ESA_XINTER(X, X1, Y1, X2, Y2)

      IMPLICIT NONE

      REAL*8 X, X1, X2, Y1, Y2

      IF (X2 .NE. X1) THEN
        ESA_XINTER = y1 + (x-x1) * (y2-y1) / (x2-x1)
      ELSE
        ESA_XINTER = Y1
      END IF

      RETURN
      END

C***********************************************************************
C
C     INTEGER*4 FUNCTION ESA_LStrng
C
C     Returns the index in the Energy string (ESTR)
C     that the requested L string starts.
C     
C***********************************************************************

      INTEGER*4 FUNCTION ESA_LStrng(ESTR, Header, L, len)

      IMPLICIT NONE

      INTEGER*4 index, len
      INTEGER*4 ESTR(*)
      REAL*4 Header(*)
      REAL*8 L, MapL

C SLEN  - Position in the E string of the E string length
C LPOS  - Offset in the L string of the L value
C LSCL  - Position in the header of the L scale factor

      INTEGER*4 SLEN, LPOS, LSCL
      DATA SLEN, LPOS, LSCL /1, 1, 5/

C index = position of the first L string in the E string
      index = 3
      MapL = 0.0

      IF (L .EQ. 0.0) THEN
        ESA_LStrng = index - 1
        return
      END IF

 10   IF ((L .LE. MapL) .OR. (index .GT. ESTR(SLEN))) GO TO 20
      len = ESTR(index)
      MapL = ESTR(index+Lpos) / HEADER(LSCL)
      IF (L .GT. MapL) index = index + len
      GO TO 10
 20   ESA_LStrng = index - 1
 
      RETURN
      END      
C
C***********************************************************************
C
C     ESA_EVals
C
C     Searches through the Ax8 model for all of the energies it contains.
C
C***********************************************************************

      SUBROUTINE ESA_EVals(MAP, Header, E, Npts, iE)

      IMPLICIT NONE

      INTEGER*4 MAP(*)
      REAL*4 HEADER(*)
      REAL*8 E(*)
      INTEGER*4 NPTS, iE(*)

C ESCL   : Position in the Header of the energy scale factor
C MAPLEN : Position in the Header of the total map length

      INTEGER*4 ESCL, MAPLEN
      INTEGER*4 INDEX
      DATA ESCL, MAPLEN /4, 8/   
      
      index = 1
      Npts = 0

 10   Npts = Npts + 1
      iE(Npts) = index+1
      E(Npts) = MAP(index + 1) / HEADER(ESCL)
      index = index + MAP(index)    
      IF ((index .LE. HEADER(MAPLEN)) .AND. (MAP(index) .NE. 0)) GOTO 10

      RETURN
      END
C
C***********************************************************************
C
C     ESA_LVals
C
C     Searches through the Energy string for all of the L values 
C     it contains.
C
C***********************************************************************

      SUBROUTINE ESA_LVals(ESTR, Header, L, Npts)

      IMPLICIT NONE

      INTEGER*4 ESTR(*), NPTS
      REAL*4 HEADER(*)
      REAL*8 L(*)

C SLEN  - Position in the E string of the E string length
C LPOS  - Offset in the L string of the L value
C LSCL  - Position in the header of the L scale factor

      INTEGER*4 SLEN, LPOS, LSCL, INDEX
      DATA SLEN, LPOS, LSCL /1, 1, 5/

      index = 3
      Npts = 0

 10   IF (index .GE. ESTR(SLEN)) RETURN
      Npts = Npts + 1
      L(Npts) = ESTR(index + LPOS) / HEADER(LSCL)
      index = index + ESTR(index)
      GO TO 10
 
      END
C
C***********************************************************************
C
C     ESA_BB0Val
C
C     Returns the BoBo and Log10(flux) values contained in the L string
C     (MAP)
C
C***********************************************************************

      SUBROUTINE ESA_BB0Val(LSTR, Header, BoB0, lnFlux, Npts)

      IMPLICIT NONE

      INTEGER*4 LSTR(*), NPTS, LEN, I
      REAL*4 HEADER(*)
      REAL*8 BoB0(*), lnFlux(*)

C SLEN   - Position in the B string of the B string length
C LPOS   - Offset in the L string of the L value
C LSCL   - Position in the header of the L scale factor
C FlxSCL - Position in the Header of the lnFlux scale factor
C FlxInc - Unscaled increment in B between string values
C FlxOff - Offset in B string of the B0 lnFlux value

      INTEGER*4 SLEN, BSCL, FlxSCL, FlxInc, FlxOff
      DATA SLEN, BSCL, FlxSCL, FlxInc, FlxOff /1, 6, 7, -256, 3/

      Npts = 1
      LEN = LSTR(SLEN)
      BoB0(1) = 1.0
      lnFlux(1) = LSTR(Flxoff) / HEADER(FlxScl)

      IF (LEN .LT. 4) RETURN

      Npts = 0
      i = 4

 10   IF ((i .GT. LEN) .OR. (LSTR(i) .LE. 0)) RETURN
      Npts = Npts + 1
      lnFlux(i-2) = (LSTR(FlxOff)+FlxInc*(i-FlxOff)) / HEADER(FlxSCL)
      BoB0(i-2) = BoB0(i-FlxOff) + LSTR(i) / HEADER(BSCL)
      i = i + 1
      GO TO 10

      end
C
C***********************************************************************
C
C     ESA_BVals
C
C     Returns the B and Log10(flux) values contained in the L string
C     (MAP)
C
C***********************************************************************

      SUBROUTINE ESA_BVals(LSTR, Header, B, lnFlux, Npts)

      IMPLICIT NONE

      INTEGER*4 LSTR(*), NPTS
      REAL*4 HEADER(*)
      REAL*8 B(*), lnFlux(*)
      REAL*8 BoB0(40), L, Bo, ESA_B0 

C LPOS - Position in L string of the L value
C LSCL - Position in the Header of the L scale factor

      INTEGER*4 LPOS, LSCL, I
      DATA LPOS, LSCL /2, 5/

      L = LSTR(LPOS) / Header(LSCL)
      CALL ESA_BB0Val(LSTR, Header, BoB0, lnFlux, Npts)
      Bo = ESA_B0(L)

      DO i=1,Npts
        B(i) = BoB0(i) * Bo
        IF (B(i) .EQ. 0.0) lnFlux(i) = 0.0
      END DO
   
      return
      end
C
C****************************************************************************
C
C     REAL*8 FUNCTION ESA_XPHI
C     
C     Computes the hybrid magnetic field coordinate=
C
C                   ( B/B0 - 1     )
C     ESA_XPHI= ASIN( --------     )
C                   ( Bmax/B0 - 1  )
C
C     Where Bmax is the atmospheric cutoff value for B.
C     If Bmax = 1, then ESA_XPHI = 90. ( Arcsin( 1.0) )
C
C****************************************************************************

      REAL*8 FUNCTION ESA_XPHI(BoB0, BMax)

      IMPLICIT NONE

      REAL*8 BoB0, BMax
      REAL*8 sine, RADDEG

      RADDEG = 45.0 / ATAN(1.0)
C ARCSIN( >1.0) --- BoB0 beyond the atmospheric cut off.
      IF (Bob0 .GT. BMax) then
        ESA_XPHI = -1.0                 
        RETURN
      END IF

      IF (BMax .NE. 1.0) THEN
         SINE = (BoB0-1.0) / (BMax-1.0)
         IF ((-1.0 .LE. SINE) .AND. (SINE .LE. 1.0)) THEN
           ESA_XPHI = ASIN(SINE) * RADDEG
         ELSE
           ESA_XPHI = -1.0
         END IF
      ELSE
         ESA_XPHI = -1.0
      END IF

      return
      END
C
C****************************************************************************
C
C     REAL*8 FUNCTION ESA_B0
C
C     computes the magnetic field strength for an L shell at the magnetic
C     equator.
C
C****************************************************************************

      REAL*8 FUNCTION ESA_B0(L)
      IMPLICIT NONE
      REAL*8 L

      IF (L .GT. 0.0) THEN
        ESA_B0 = 0.311653 / (L*L*L) ! McIlwain's magnetic moment 0.311653 Gauss
      ELSE
        ESA_B0 = 0.0
      END IF

      RETURN
      END

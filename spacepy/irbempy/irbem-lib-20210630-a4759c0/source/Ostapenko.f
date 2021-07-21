!***************************************************************************************************
! Copyright 1997, Ostapenko
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
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BAZ_T(X,Y,Z,N,B)
C     BAZIS FOR TILTE-MODEL(S) MAGNETIC FIELD
C     X,Y,Z = X,Y,Z/10.D0 !  IN NORMALIZATION UNITS
      INTEGER*4 I,J,K,L,N
      REAL*8 X,Y,Z,B(3,1),SYM(3,10),ASY(3,10),JWORK(3,10)
      K=0
C 1-3:******** 1.Usym(1,3,5  ) FOR N=17 ******
C 1-4:******** 1.Usym(1,3,5,7) FOR N=29 ******
      L=3
      IF(N.EQ.29) L=4
      CALL PTNCL(X,Y,Z,7,SYM,ASY)   ! CALL FIRST TIME ONLY
      DO I=1,L
         DO J=1,3
            B(J,I+K)=SYM(J,2*I-1)
         ENDDO
      ENDDO
      K=K+L !  K=3/4
C 4- 6:******* 2.JETsym(1:3) FOR N=17 *********************
C 5-10:******* 2.JETsym(1:6) FOR N=29 *********************
      L=3
      IF(N.EQ.29) L=6
      CALL JETSYM(X,Y,Z,L,JWORK)
      DO I=1,L
         DO J=1,3
            B(J,I+K)=JWORK(J,I)
         ENDDO
      ENDDO
      K=K+L !  K=6/10
C   7,8:***** 3.Uasy(2,4)   FOR N=17 *************
C 11-13:***** 3.Uasy(2,4,6) FOR N=29 *************
      L=2
      IF(N.EQ.29) L=3
      DO I=1,L
         DO J=1,3
            B(J,I+K)=ASY(J,2*I)
         ENDDO
      ENDDO
      K=K+L !  K=8/13
C  9-12:***** 4.JETSasy(1:4) FOR N=17 *******************
C 14-22:***** 4.JETSasy(1:9) FOR N=29 *******************
      L=4
      IF(N.EQ.29) L=9
      CALL JETASY(X,Y,Z,L,JWORK)
      DO I=1,L
         DO J=1,3
            B(J,I+K)=JWORK(J,I)
         ENDDO
      ENDDO
      K=K+L !  K=12/22
C 13,14:**** 5.Usym(2,4)*sin(PSI)   FOR N=17 ************
C 23-25:**** 5.Usym(2,4,6)*sin(PSI) FOR N=29 ************
      L=2
      IF(N.EQ.29) L=3
      DO I=1,L
         DO J=1,3
            B(J,I+K)=SYM(J,2*I)
         ENDDO
      ENDDO
      K=K+L !  K=14/25
C 15-17:*** 6.Uasy(1,3,5)*sin(PSI)   FOR N=17 **********
C 26-29:*** 6.Uasy(1,3,5,7)*sin(PSI) FOR N=29 **********
      L=3
      IF(N.EQ.29) L=4
      DO I=1,L
         DO J=1,3
            B(J,I+K)=ASY(J,2*I-1)
         ENDDO
      ENDDO
      K=K+L !  K=17/29
      IF(K.EQ.N) RETURN
      WRITE(*,*) ' ERROR GENERATED IN OSTAPENKO-MALTSEV 1997'
      WRITE(*,*) ' ERROR IN MODULE BAZ_T K#MF:',K,'#',N
      STOP
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE JETSYM(X,Y,Z,L,B)
C     SYMMETRIC JET
      REAL*8 X,Y,Z,B(3,*),Z2,P2
      INTEGER*4 I,L,J
      DO J=1,L
         DO I=1,3
            B(I,J)=0.D0
         ENDDO
      ENDDO
      P2=X*X+Y*Y
      Z2=Z*Z
      B(3,1)=P2
      B(3,2)=P2*P2
      B(2,3)=-2.D0*Z*Z2    !  #3: [-2pz3,z4]
      B(3,3)=Z2*Z2
      IF(L.EQ.6) THEN
         B(3,4)=P2**3
         B(2,5)=-Z*Z2*(P2-2.D0*Z2/5.D0) !       revised 28 april 1997
         B(3,5)=Z2*Z2*(P2-2.D0*Z2/15.D0)  !     revised 28 april 1997
         B(2,6)=-3.D0*Z**5 !  #6: [-3pz5,z6]
         B(3,6)=Z2**3
      ENDIF
      DO J=1,L
         B(1,J)=X*B(2,J)
         B(2,J)=Y*B(2,J)
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE JETASY(X,Y,Z,L,B)
C     ASYMMETRIC JET ~R3 - REVISED 06/06/95
      REAL*8 X,Y,Z,B(3,*),P2,Z2,WORK
      INTEGER*4 L,I,J
      DO J=1,L
         DO I=1,3
            B(I,J)=0.D0
         ENDDO
      ENDDO
      P2=X*X+Y*Y
      Z2=Z*Z
      B(1,1)=Z
      B(3,2)=P2*X
      B(1,3)=Z2*Z
      B(1,4)=-Y*Y*Z
      B(2,4)= X*Y*Z
      B(3,4)=-X*Z2/2.D0
      IF(L.EQ.4) RETURN
C     ASYMMETRIC JET ~R5 - REVISED 16/06/95, NORMALIZATION 14/05/97
      B(3,5)=P2*P2*X       !  NORMALIZATION = 1
      WORK=3.D0*P2*Z       !  NORMALIZATION = 3
      B(1,6)=-Y*Y*WORK
      B(2,6)= X*Y*WORK
      B(3,6)=-X*Z*WORK/2.D0
      WORK=10.D0*P2*X*Z    !  NORMALIZATION = 10
      B(1,7)= X*WORK/5.D0
      B(2,7)= Y*WORK/5.D0
      B(3,7)=-Z*WORK/2.D0
      WORK=5.D0*Z2*Z       !  NORMALIZATION = 5
      B(1,8)=-Y*Y*WORK
      B(2,8)= X*Y*WORK
      B(3,8)=-X*Z*WORK/4.D0
      WORK=10.D0*X*Z*Z2
      B(1,9)= X*WORK/3.D0  !  NORMALIZATION = 10
      B(2,9)= Y*WORK/3.D0
      B(3,9)=-Z*WORK/4.D0
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PTNCL(X,Y,Z,N,SYM,ASY)
C     DERIVES grad [ P(i,j)*r^i]
      REAL*8 X,Y,Z,SYM(3,7),ASY(3,7),P(0:8,0:8),WORK,R2,C,R(-1:10),RN
      INTEGER*4 N,I,J
C
      IF(N.GT.7.OR.N.LT.1) THEN
         WRITE(*,*) ' ERROR GENERATED IN OSTAPENKO-MALTSEV 1997'
         WRITE(*,'(26H PTNCL: WRONG PARAMETER N=I2)') N
         STOP
      ENDIF
      DO J=1,N
         DO I=1,3
            SYM(I,J)=0.D0
            ASY(I,J)=0.D0
         ENDDO
      ENDDO
      SYM(3,1)=1.D0
      ASY(1,1)=1.D0
      R2=X*X+Y*Y+Z*Z
      IF(R2.EQ.0.D0) RETURN
      R(0)=1.D0
      R(1)=SQRT(R2)
      R(2)=R2
      R(-1)=1.D0/R(1)
      DO I=3,N+3
         R(I)=R(I-1)*R(1)
      ENDDO
      C=Z/R(1)
      CALL LEGNDR(C,N,P)
      DO I=2,N
C ZERO GARMONICS
         WORK=R(I-2)*(I*P(I,0)-C*P(I,1))
         SYM(1,I)=X*WORK                  ! d/dx
         SYM(2,I)=Y*WORK                  ! d/dy
         SYM(3,I)=Z*WORK+R(I-1)*P(I,1)    ! d/dz
C FIRST GARMONICS
         WORK=X*R(I-3)*((I-1)*P(I,1)-C*P(I,2))
         ASY(1,I)=X*WORK+R(I-1)*P(I,1)    ! d/dx
         ASY(2,I)=Y*WORK                  ! d/dy
         ASY(3,I)=Z*WORK+X*R(I-2)*P(I,2)  ! d/dz
      ENDDO
C SCHMIDT NORMALIZATION FOR M=1,N=2,3,...
      WORK=1.D0
      DO I=2,N
         WORK=(I-1)*WORK/(I+1)
         RN=SQRT(WORK)
         DO J=1,3
            ASY(J,I)=ASY(J,I)*RN
         ENDDO
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE LEGNDR(X,N,P)
C     P(i,0) = LEGENDRE'S POLYNOMS Pi(x)
C     AND ITS DERIVES: P(i,j)=(d/dx)^j (Pi(x))
      REAL*8 X,P(0:8,0:8)
      INTEGER*4 N,I,J
      IF(N.GT.7.OR.N.LT.1) THEN
         WRITE(*,*) ' ERROR GENERATED IN OSTAPENKO-MALTSEV 1997'
         WRITE(*,'(37H LEGENDRE POLYNOM: WRONG PARAMETER N=I2)') N
         STOP
      ENDIF
C
      P(0,0)=1.D0
      P(0,1)=0.D0
      P(1,0)=X
      P(1,1)=1.D0
      DO I=1,N-1
         P(I,I+1)=0.D0
         DO J=I+1,1,-1
            P(I+1,J)=X*P(I,J)+(I+J)*P(I,J-1)
         ENDDO
         P(I+1,0)=((2*I+1)*X*P(I,0)-I*P(I-1,0))/(I+1)
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BOM97(RE,B)
C     EXTERNAL MODEL MAGNETIC FIELD (WITH TILTE)
      INTEGER*4 MF,I,J,NA
      PARAMETER (MF=29)
      REAL*8 A(MF),RE(3),R(3),X,Y,Z,B(3),BAZIS(3,MF)
      COMMON/COEFOM97/A,NA
C
      DO J=1,3
         R(J)=RE(J)/10.D0 !  TO NORMALIZATION UNITS
         B(J)=0.D0
      ENDDO
      X=R(1)
      Y=R(2)
      Z=R(3)
      CALL BAZ_T(X,Y,Z,NA,BAZIS)
      DO I=1,NA
         DO J=1,3
            B(J)=B(J)+BAZIS(J,I)*A(I)
         ENDDO
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SET_A(Dst,Pdyn,Kp,IMFz,SN)
      INTEGER*4 MF,I,J,K,NC,L,NA
      PARAMETER (MF=29,NC=4)
      REAL*8 AA(85),C(4),SN,A(MF),Dst,Pdyn,Kp,IMFz
      COMMON/COEFOM97/A,NA
C Input parameters:
C Dst (nT), Pdyn (nPa),Kp (numeric),IMFz (nT), SN=sin(tilt of dipole)
C Kp (as key )=  0,  0+,  1-,   1,  1+,...
C Kp (numeric)=0.0, 0.3, 0.7, 1.0, 1.3,...
      DATA AA/
     ,-4.3392E+01, 2.0364E+01,-7.4532E-01,-3.2941E+00,-1.6294E+00,
     , 4.0308E+01,-1.3734E+01, 6.4463E+00, 3.4054E+00, 2.8944E+00,
     , 7.0475E+00,-3.2920E+00, 2.2886E+00, 1.0694E-01, 8.0280E-01,
     , 1.3347E+02,-6.1188E+01, 2.5064E+01, 2.5815E+00, 1.1370E+01,
     ,-3.6968E+01, 2.4557E+01,-1.2768E+01, 4.1262E+00,-4.2696E+00,
     ,-1.3010E+02, 4.6964E+01,-2.5640E+01,-1.4840E+01,-9.6417E+00,
     , 2.3632E+01, 4.9102E+00, 7.9419E+00, 4.3164E+00, 2.6925E+00,
     ,-1.7139E+00,-7.4376E-01, 2.8973E-01,-9.6157E-01,-5.0774E-01,
     , 1.2742E+00,-1.5849E+01, 7.1720E+00,-1.6217E+00,-1.0925E+01,
     ,-2.1437E+01,-1.0040E+01,-7.1910E+00,-9.9341E+00,-5.0171E+00,
     , 4.4863E+00, 3.1827E+00,-9.6490E+00,-6.9168E-01, 3.4412E+00,
     , 2.3315E+01,-6.8824E+00, 1.2323E+01, 3.5851E+00,-6.1224E+00,
     , 2.3179E+01,-2.0127E+00, 9.0600E+00,-1.8949E-01, 1.5000E-01,
     ,-1.9672E+00,-1.3847E+00, 2.9655E+00,-3.4092E+00, 1.4149E+00,
     , 1.2970E+01,-4.3018E+00, 6.4866E+00,-1.2664E-01,-2.1103E+00,
     , 5.7215E+00,-3.1930E+00, 2.2242E+00, 3.5120E-01, 9.1624E-02,
     , 4.9085E+00,-1.0754E+00, 3.5705E+00,-9.9892E-01,-6.2140E-01/
C     Input parameters: Dst (nT), Pdyn (nPa),Kp (numeric),IMFz (nT)
C     Convert from values to normalized parameters:
      C(1)=(Dst+16.9367D0)/25.2834D0
      C(2)=(Pdyn-2.278138D0)/1.882804D0
      C(3)=(Kp-2.30896D0)/1.35401D0
      C(4)=(IMFz-.0180D0)/3.7051D0
      NA=17
C     WRITE(*,'(24H NUMBER OF COEFFICIENTS=I3)') NA
      DO I=1,NA
         K=(I-1)*(NC+1)+1           !!!! REVISED 24/12/97
         A(I)=AA(K)                 !!!! REVISED 24/12/97
         DO J=1,NC
            A(I)=A(I)+AA(K+J)*C(J)  !!!! REVISED 24/12/97
         ENDDO
      ENDDO
      if( na.eq.17 )then
         L=13
      elseif( na.eq.29 )then
         L=23
      else
         WRITE(*,*) ' ERROR GENERATED IN EXT530, OSTAPENKO-MALTSEV 1997'
         WRITE(*,*) 'WRONG PARAMETER NA'
         STOP
      ENDif
      DO I=L,NA
         A(I)=A(I)*SN
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


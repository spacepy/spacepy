!***************************************************************************************************
! Copyright  2004 S. Bourdarie
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
C
C

C     Legacy wrapper w/ R0=1.0
       SUBROUTINE field_line_tracing(
     &         lati,longi,alti,Lm,leI0,Bposit,Bmin,posit,Nposit)
       IMPLICIT NONE
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 150, Ntet = 720)
       INTEGER*4  Nposit
       REAL*8     lati,longi,alti,R0
       REAL*8     Lm,leI0,Bmin,xx0(3)
       REAL*8     posit(3,20*Nreb),Bposit(20*Nreb)

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
C
       call field_line_tracing_opt2(xx0,1.D0,
     & Lm,leI0,Bposit,Bmin,posit,Nposit)

       RETURN
       END
       SUBROUTINE field_line_tracing2(
     &         lati,longi,alti,R0,Lm,leI0,Bposit,Bmin,posit,Nposit)
C
C    modified from field_line_tracing to add R0 parameter: radius (Re) of
C    reference surface
       IMPLICIT NONE
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 150, Ntet = 720)
       INTEGER*4  Nposit
       REAL*8     lati,longi,alti,R0
       REAL*8     Lm,leI0,Bmin,xx0(3)
       REAL*8     posit(3,20*Nreb),Bposit(20*Nreb)

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
C
       call field_line_tracing_opt2(xx0,R0,
     & Lm,leI0,Bposit,Bmin,posit,Nposit)

       RETURN
       END
       
       SUBROUTINE field_line_tracing_opt(
     &         xx0,Lm,leI0,Bposit,Bmin,posit,Nposit)
       IMPLICIT NONE
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 150, Ntet = 720)
       INTEGER*4  Nposit
       REAL*8     xx0(3),R0
       REAL*8     Lm,leI0,Bmin
       REAL*8     posit(3,20*Nreb),Bposit(20*Nreb)

       call field_line_tracing_opt2(xx0,1.D0,
     & Lm,leI0,Bposit,Bmin,posit,Nposit)
       RETURN
       END

       SUBROUTINE field_line_tracing_opt2(
     &         xx0,R0,Lm,leI0,Bposit,Bmin,posit,Nposit)
C
C    modified from field_line_tracing to add R0 parameter: radius (Re) of
C    reference surface
       IMPLICIT NONE
       INCLUDE 'variables.inc'
C
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 150, Ntet = 720)
C
       INTEGER*4  k_ext,k_l,kint
       INTEGER*4  Nrebmax
       REAL*8     rr,rr2
       REAL*8     xx0(3),xx(3),x1(3),x2(3)
       REAL*8     xmin(3)
       REAL*8     lati,longi,alti
       REAL*8     B(3),Bl,B0,Bmin,B1,B3
       REAL*8     dsreb,smin

       INTEGER*4  I,J,Iflag,Iflag_I,Ilflag,ind,II,Ifail
       INTEGER*4  Nposit
       REAL*8     Lm,Lstar,Lb
       REAL*8     leI,leI0,leI1
       REAL*8     XY,YY
       REAL*8     aa,bb
C
       REAL*8     tt
       REAL*8     tetl,tet1,dtet
       REAL*8     somme
C
       REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
C
       REAL*8     posit(3,20*Nreb),Bposit(20*Nreb)
       real*8     R0,R02 ! R0^2
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
       COMMON /flag_L/Ilflag
       COMMON /magmod/k_ext,k_l,kint
C
C
       R02 = R0*R0
C
       Nrebmax = 20*Nreb
c       write(*,*)'Nrebmax',Nrebmax
C
       Lm = baddata
       leI0 = 0.D0
C
       CALL GEO_SM(xx0,xx)
       rr = SQRT(xx(1)*xx(1)+xx(2)*xx(2)+xx(3)*xx(3))
       tt = ACOS(xx(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
c       write(6,*)'L bete ',Lb
C
       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
	  RETURN
       ENDIF
       Bmin = B0
C
       dsreb = Lb/Nreb
C
C calcul du sens du depart
C
       CALL sksyst(-dsreb,xx0,x1,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
	  RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
	  RETURN
       ENDIF
       B3 = Bl
C
C attention cas equatorial
C
       IF(B1.GT.B0 .AND. B3.GT.B0)THEN
         aa = 0.5D0*(B3+B1-2.D0*B0)
         bb = 0.5D0*(B3-B1)
         smin = -0.5D0*bb/aa
         Bmin = B0 - aa*smin*smin
         leI0 = SQRT(1.D0-Bmin/B0)*2.D0*ABS(smin*dsreb)
         Lm = (Bo/Bmin)**(1.D0/3.D0)
c         write(6,*)'L McIlwain eq ',B0,leI0,Lm
         GOTO 100
       ENDIF
       IF (B3.GT.B1) THEN
         dsreb = -dsreb
       ENDIF
C
C calcul de la ligne de champ et de I
C
       Bmin = B0
       leI = 0.D0
       ind=0
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
C
c       write(6,*)dsreb
       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            leI0 = baddata
            Bmin = baddata
	    RETURN
	 ENDIF
c	 write(6,*)J,x1(1),x1(2),x1(3),Bl
         IF (Bl.LT.Bmin) THEN
           xmin(1) = x2(1)
           xmin(2) = x2(2)
           xmin(3) = x2(3)
           Bmin = Bl
         ENDIF
         IF (Bl.GT.B0) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
         leI = leI + SQRT(1.D0-Bl/B0)
	 B1 = Bl
       ENDDO
20     CONTINUE
c	write(6,*)J,leI
C
       IF (J.GE.Nrebmax) THEN !open field line
        leI0 = baddata
        Bmin = baddata
	RETURN
       ENDIF
C
       leI = leI+0.5D0*SQRT(1.D0-B1/B0)*(B0-Bl)/(Bl-B1)
       leI = leI*ABS(dsreb)
       leI0 = leI
C
C calcul de L Mc Ilwain (Mc Ilwain-Hilton)
C
       XY = leI*leI*leI*B0/Bo
       YY = 1.D0 + 1.35047D0*XY**(1.D0/3.D0)
     &      + 0.465376D0*XY**(2.D0/3.D0)
     &      + 0.0475455D0*XY
       Lm = (Bo*YY/B0)**(1.D0/3.D0)
c       write(6,*)'L McIlwain ',B0,leI0,Lm
C
C calcul de Bmin
C
       CALL sksyst(dsreb,xmin,x1,B3,Ifail)
       IF (Ifail.LT.0) THEN
          Bmin = baddata
	  RETURN
       ENDIF
       CALL sksyst(-dsreb,xmin,x1,B1,Ifail)
       IF (Ifail.LT.0) THEN
          Bmin = baddata
	  RETURN
       ENDIF
       aa = 0.5D0*(B3+B1-2.D0*Bmin)
       bb = 0.5D0*(B3-B1)
       smin = -0.5D0*bb/aa
       Bmin = Bmin - aa*smin*smin
       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.D0) THEN
        Lm = -Lm
       ENDIF
C
100    CONTINUE
C
C calcul du point sur la ligne de champ a la surface de la terre du
C cote nord
C
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
       dsreb = ABS(dsreb)
       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) RETURN
	 rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
	 IF (rr.LT.R0) GOTO 102
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
102    CONTINUE
       smin = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
       smin = (R0-smin)/(rr-smin)
       CALL sksyst(smin*dsreb,x1,x2,Bl,Ifail)
       IF (Ifail.LT.0) RETURN
C
c trace la ligne de champ complete.
c
       x1(1) = x2(1)
       x1(2) = x2(2)
       x1(3) = x2(3)
       ind=1
       Nposit=ind
       posit(1,ind)=x1(1)
       posit(2,ind)=x1(2)
       posit(3,ind)=x1(3)
       Bposit(ind)=Bl
       DO J = 1,Nrebmax-1 ! this is the corrected version
                            !, gives 3000 values, A. Kellerman
         CALL sksyst(-dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) RETURN
	 ind=ind+1
	 posit(1,ind)=x2(1)
	 posit(2,ind)=x2(2)
	 posit(3,ind)=x2(3)
         Bposit(ind)=Bl
c	 write(6,*)J,x1(1),x1(2),x1(3),Bl
         rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
         IF (rr2.LT.R02) GOTO 201
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
201    CONTINUE
       Nposit=ind
       END
C

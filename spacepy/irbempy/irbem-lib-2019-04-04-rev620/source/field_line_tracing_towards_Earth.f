!***************************************************************************************************
! Copyright  2009 S. Bourdarie
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
       SUBROUTINE field_line_tracing_towards_Earth(
     &         lati,longi,alti,dsreb,posit,Nposit)
C
       IMPLICIT NONE
       
       INTEGER*4  Nreb
       PARAMETER (Nreb = 150)
       REAL*8     xx0(3)
       REAL*8     lati,longi,alti
       REAL*8     dsreb
       INTEGER*4  Nposit
       REAL*8     posit(3,20*Nreb)

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
C
       call field_line_tracing_towards_Earth_opt(
     &         xx0,dsreb,posit,Nposit)

       RETURN
       END

       SUBROUTINE field_line_tracing_towards_Earth_opt(
     &         xx0,dsreb,posit,Nposit)
C
       IMPLICIT NONE
       INCLUDE 'variables.inc'
C
       INTEGER*4  Nreb
       PARAMETER (Nreb = 150)
C
       INTEGER*4  Nrebmax
       REAL*8     rr,rr2
       REAL*8     xx0(3),xx(3),x1(3),x2(3)
       REAL*8     xmin(3)
       REAL*8     lati,longi,alti
       REAL*8     B(3),Bl,B0,Bmin,B1,B3
       REAL*8     dsreb

       INTEGER*4  J,ind,Ifail
       INTEGER*4  Nposit
C
       REAL*8     posit(3,20*Nreb)
C
       Nrebmax = 20*Nreb
C
       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) THEN
         RETURN
       ENDIF
       Bmin = B0
C
C calcul du sens du depart
C
       CALL sksyst(-dsreb,xx0,x1,Bl,Ifail)
       IF (Ifail.LT.0) THEN
         RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
       IF (Ifail.LT.0) THEN
         RETURN
       ENDIF
       B3 = Bl
C
       IF (B3.LT.B1) THEN
         dsreb = -dsreb
       ENDIF
C
c trace la ligne de champ complete.
c
       x1(1) = xx0(1)
       x1(2) = xx0(2)
       x1(3) = xx0(3)
       ind=1
       Nposit=ind
       posit(1,ind)=x1(1)
       posit(2,ind)=x1(2)
       posit(3,ind)=x1(3)
       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) RETURN
         ind=ind+1
         posit(1,ind)=x2(1)
         posit(2,ind)=x2(2)
         posit(3,ind)=x2(3)
c        write(6,*)J,x1(1),x1(2),x1(3),Bl
         rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
         IF (rr2.LT.1.) GOTO 201
         x1(1) = x2(1)
         x1(2) = x2(2)
         x1(3) = x2(3)
         if (ind .eq. Nrebmax) goto 201
       ENDDO
201    CONTINUE
       Nposit=ind
       RETURN
       END
C

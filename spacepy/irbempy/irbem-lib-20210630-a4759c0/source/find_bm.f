!***************************************************************************************************
! Copyright 2004 S. Bourdarie
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
C S. Bourdarie (June 2004)
c Modified S./ Bourdarie (July 2004)
c
C Routine to find mirror point of a trapped particle
C
       SUBROUTINE find_bm(
     &         lati,longi,alti,alpha,Bposit,Bmir,xmin)
C
       IMPLICIT NONE
       REAL*8     xx0(3),xmin(3)
       REAL*8     lati,longi,alti
       REAL*8     Bposit,Bmir
       REAL*8     alpha

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
C
       call find_bm_nalpha (xx0,1,alpha,Bposit,Bmir,xmin)

       RETURN
       END

       SUBROUTINE find_bm_nalpha(xx0,nalp,alpha,Bposit,Bmir,xmin)
C
       IMPLICIT NONE
       INCLUDE 'variables.inc'
C
       INTEGER*4  Nreb
       PARAMETER (Nreb = 50)
C
       INTEGER*4  k_ext,k_l,kint,Ifail
       INTEGER*4  Nrebmax
       REAL*8     rr,alpha(25)
       REAL*8     xx0(3),xx(3),x1(3),x2(3)
       REAL*8     xmin(3,25)
       REAL*8     B(3),Bl,B0,B1,B3
       REAL*8     dsreb

       INTEGER*4  I,J,ind,ii,nalp
       REAL*8     Lb,sn
       REAL*8     leI0
C
       REAL*8     tt
       REAL*8     cste
C
       REAL*8     Bposit,Bmir(25)
       REAL*8     sn2,sn2max
C
       COMMON /magmod/k_ext,k_l,kint
       REAL*8     pi,rad
       common /rconst/rad,pi
C
C
       Nrebmax = 20*Nreb
C
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
          Bposit=baddata
          do 50 ii=1,nalp
	    Bmir(ii)=baddata
            xmin(1,ii) = baddata
            xmin(2,ii) = baddata
            xmin(3,ii) = baddata
50        continue
	  RETURN
       ENDIF
       Bposit=B0
       
       do 100 ii=1,nalp
       
       if ( alpha(ii).eq.90.d0 ) then
           Bmir(ii) = 0.d0
	   xmin(1,ii) = xx0(1)
           xmin(2,ii) = xx0(2)
           xmin(3,ii) = xx0(3)
	   goto 100
	endif
       Bmir(ii) = B0

       sn = SIN(alpha(ii)*rad)  
       sn2=sn*sn
       Sn2max=sn2
       cste=B0/sn2
       dsreb = Lb/(Nreb*1.d0)
C
C calcul du sens du depart
C
       CALL sksyst(-dsreb,xx0,x1,Bl,Ifail)
       IF (Ifail.LT.0) THEN
	  Bmir(ii)=baddata
          xmin(1,ii) = baddata
          xmin(2,ii) = baddata
          xmin(3,ii) = baddata
	  goto 100
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
       IF (Ifail.LT.0) THEN
	  Bmir(ii)=baddata
          xmin(1,ii) = baddata
          xmin(2,ii) = baddata
          xmin(3,ii) = baddata
	  goto 100
       ENDIF
       B3 = Bl
C
       IF (B1.GT.B3) dsreb = -dsreb
C
C calcul de la ligne de champ
C
       ind=0
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
C
c       write(6,*)dsreb
       xmin(1,ii)=x1(1)
       xmin(2,ii)=x1(2)
       xmin(3,ii)=x1(3)
       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
	  Bmir(ii)=baddata
          xmin(1,ii) = baddata
          xmin(2,ii) = baddata
          xmin(3,ii) = baddata
	  goto 100
         ENDIF
c
	 sn2=Bl/cste
	 if (sn2 .GT.1.d0) GOTO 20
         IF (sn2.GT.Sn2max) THEN
           xmin(1,ii) = x2(1)
           xmin(2,ii) = x2(2)
           xmin(3,ii) = x2(3)
           Sn2max = sn2
	   Bmir(ii)=Bl
         ENDIF
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
20     CONTINUE
C
       IF (J.GE.Nrebmax) THEN !open field line
	  Bmir(ii)=baddata
          xmin(1,ii) = baddata
          xmin(2,ii) = baddata
          xmin(3,ii) = baddata
	  goto 100
       ENDIF
C
C calcul de Bmirror point
C
       DO i=1,1000
         CALL sksyst(dsreb/1000.d0,xmin(1,ii),x1,Bl,Ifail)
         IF (Ifail.LT.0) THEN
	    Bmir(ii)=baddata
            xmin(1,ii) = baddata
            xmin(2,ii) = baddata
            xmin(3,ii) = baddata
	    goto 100
         ENDIF
	 sn2=Bl/cste
	 if (sn2 .GT.1.d0) GOTO 30
         xmin(1,ii) = x1(1)
         xmin(2,ii) = x1(2)
         xmin(3,ii) = x1(3)
         Sn2max = sn2
	 Bmir(ii)=Bl
       ENDDO
30     CONTINUE
       IF (xmin(1,ii)*xmin(1,ii)+xmin(2,ii)*xmin(2,ii)+
     &      xmin(3,ii)*xmin(3,ii).LT.1.D0) THEN
	  Bmir(ii)=baddata
          xmin(1,ii) = baddata
          xmin(2,ii) = baddata
          xmin(3,ii) = baddata
	  goto 100
       ENDIF
C
100    CONTINUE
C
       RETURN
       END
C

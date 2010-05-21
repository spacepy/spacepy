C Boscher modifie pour la nieme fois le 4Feb2004
C a compiler avec IGRF.f Dipole.f Pfitzer_q.f init_nouveau.f 
C voire Pfitzer.f
C
       SUBROUTINE loc_equator(
     &         lati,longi,alti,Bmin,posit)
C
       IMPLICIT NONE
       INCLUDE 'variables.inc'
C
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 50, Ntet = 720)
C
       INTEGER*4  k_ext,k_l,kint
       INTEGER*4  Nrebmax
       REAL*8     rr
       REAL*8     xx0(3),xx(3),x1(3),x2(3)
       REAL*8     xmin(3)
       REAL*8     lati,longi,alti
       REAL*8     B(3),Bl,B0,Bmin,B1,B3
       REAL*8     dsreb,smin

       INTEGER*4  I,J,Ifail
       REAL*8     Lb
       REAL*8     aa,bb
C
       REAL*8     pi,rad,tt
       REAL*8     dtet
C
       REAL*8     posit(3)
C
       COMMON /magmod/k_ext,k_l,kint
C
C
       pi = 4.D0*ATAN(1.D0)
       rad = pi/180.D0
       dtet = pi/Ntet
C
       Nrebmax = 10*Nreb
C
       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
C           
       CALL GEO_SM(xx0,xx)
       rr = SQRT(xx(1)*xx(1)+xx(2)*xx(2)+xx(3)*xx(3))
       tt = ACOS(xx(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
C
       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) THEN
	  Bmin=baddata
          posit(1) = baddata
          posit(2) = baddata
          posit(3) = baddata
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
	  Bmin=baddata
          posit(1) = baddata
          posit(2) = baddata
          posit(3) = baddata
	  RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
       IF (Ifail.LT.0) THEN
	  Bmin=baddata
          posit(1) = baddata
          posit(2) = baddata
          posit(3) = baddata
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
         CALL sksyst(smin*dsreb,xx0,posit,Bl,Ifail)
         IF (Ifail.LT.0) THEN
	    Bmin=baddata
            posit(1) = baddata
            posit(2) = baddata
            posit(3) = baddata
	    RETURN
         ENDIF
	 RETURN
       ENDIF
       IF (B3.GT.B1) THEN
         dsreb = -dsreb
       ENDIF
C
C calcul de la ligne de champ et de I
C
       Bmin = B0
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
C
       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
	    Bmin=baddata
            posit(1) = baddata
            posit(2) = baddata
            posit(3) = baddata
	    RETURN
         ENDIF
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
	 B1 = Bl
       ENDDO
20     CONTINUE
C
       IF (J.GE.Nrebmax) THEN !open field line
          Bmin = baddata
	  posit(1)=baddata
	  posit(2)=baddata
	  posit(3)=baddata
	  RETURN
       ENDIF
c
C calcul de Bmin
C
       CALL sksyst(dsreb,xmin,x1,B3,Ifail)
       IF (Ifail.LT.0) THEN
	  Bmin=baddata
          posit(1) = baddata
          posit(2) = baddata
          posit(3) = baddata
	  RETURN
       ENDIF
       CALL sksyst(-dsreb,xmin,x1,B1,Ifail)
       IF (Ifail.LT.0) THEN
	  Bmin=baddata
          posit(1) = baddata
          posit(2) = baddata
          posit(3) = baddata
	  RETURN
       ENDIF
       aa = 0.5D0*(B3+B1-2.D0*Bmin)
       bb = 0.5D0*(B3-B1)
       smin = -0.5D0*bb/aa
       Bmin = Bmin - aa*smin*smin
       CALL sksyst(smin*dsreb,xmin,posit,Bl,Ifail)
       IF (Ifail.LT.0) THEN
	  Bmin=baddata
          posit(1) = baddata
          posit(2) = baddata
          posit(3) = baddata
	  RETURN
       ENDIF
c
       END
       

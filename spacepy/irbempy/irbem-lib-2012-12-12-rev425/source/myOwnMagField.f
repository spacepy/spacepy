      subroutine myOwnMagField(X,B)
C     User-defined B-field, must pass to/from GEO coordinates
C     This field model is a centred dipole plus a uniform,
C     southward B of 14.474nT [see Chen, Schulz, et al., JGR 1993]
C
       IMPLICIT NONE
C
       REAL*8 x(3)
       REAL*8 B(3)
       REAL*8 Bo
       REAL*8 pi,rad
       REAL*8 rr,tt,pp,Br,Bt,Bp,BxGEO,ByGEO,BzGEO
       common /rconst/rad,pi
C      COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C Uncomment the above line to use the dipole moment
C from the appropriate epoch, otherwise, use the below line
       Bo = 3.05D4
C
       rr = SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
       tt = ACOS(x(3)/rr)
       pp = ATAN2(x(2),x(1))
C
       Br = -2.D0*Bo*COS(tt)/(rr*rr*rr)
       Bt = -   Bo*SIN(tt)/(rr*rr*rr)
       Bp = 0.D0
C
       BxGEO = Br*SIN(tt)*COS(pp)
     &          + Bt*COS(tt)*COS(pp)
     &          - Bp*SIN(pp)
       ByGEO = Br*SIN(tt)*SIN(pp)
     &          + Bt*COS(tt)*SIN(pp)
     &          + Bp*COS(pp)
       BzGEO = Br*COS(tt) - Bt*SIN(tt) - 14.74
C
       B(1) = BxGEO
       B(2) = ByGEO
       B(3) = BzGEO 
       
       RETURN
       END


C       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C dipole
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE DTD(xGEO,yGEO,zGEO,BxGEO,ByGEO,BzGEO)
C
       IMPLICIT NONE
C
       REAL*8 xGEO,yGEO,zGEO
       REAL*8 xMAG,yMAG,zMAG
       REAL*8 rr,tt,pp
       REAL*8 pi
       REAL*8 Bo
       REAL*8 BxGEO,ByGEO,BzGEO
       REAL*8 BxMAG,ByMAG,BzMAG
       REAL*8 Br,Bt,Bp
       REAL*8 xc,yc,zc
       REAL*8 ct,st,cp,sp
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C
       pi = 4.D0*ATAN(1.D0)
C
       xMAG = ct*cp*(xGEO-xc) + ct*sp*(yGEO-yc) - st*(zGEO-zc)
       yMAG = -  sp*(xGEO-xc) +    cp*(yGEO-yc)
       zMAG = st*cp*(xGEO-xc) + st*sp*(yGEO-yc) + ct*(zGEO-zc)
C
       rr = SQRT(xMAG*xMAG + yMAG*yMAG + zMAG*zMAG)
       tt = ACOS(zMAG/rr)
       pp = ATAN2(yMAG,xMAG)
C
       Br = -2.D0*Bo*COS(tt)/(rr*rr*rr)
       Bt = -   Bo*SIN(tt)/(rr*rr*rr)
       Bp = 0.D0
C
       BxMAG = Br*SIN(tt)*COS(pp)
     &       + Bt*COS(tt)*COS(pp)
     &       - Bp*SIN(pp)
       ByMAG = Br*SIN(tt)*SIN(pp)
     &       + Bt*COS(tt)*SIN(pp)
     &       + Bp*COS(pp)
       BzMAG = Br*COS(tt) - Bt*SIN(tt)
C
       BxGEO = ct*cp*BxMAG - sp*ByMAG + st*cp*BzMAG
       ByGEO = ct*sp*BxMAG + cp*ByMAG + st*sp*BzMAG
       BzGEO = -st  *BxMAG            + ct   *BzMAG
C
       RETURN
       END

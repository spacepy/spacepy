!***************************************************************************************************
! Copyright , 2004 D. Boscher, S. Bourdarie
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C dipole tilte decentre
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE DTD(xGEO,yGEO,zGEO,BxGEO,ByGEO,BzGEO)
C
       IMPLICIT NONE
C
       REAL*8 xGEO,yGEO,zGEO
       REAL*8 xMAG,yMAG,zMAG
       REAL*8 rr,tt,pp
       REAL*8 Bo
       REAL*8 BxGEO,ByGEO,BzGEO
       REAL*8 BxMAG,ByMAG,BzMAG
       REAL*8 Br,Bt,Bp
       REAL*8 xc,yc,zc
       REAL*8 ct,st,cp,sp
        REAL*8     pi,rad
        common /rconst/rad,pi
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
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
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Centered dipole
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE centered_dipole(xGEO,yGEO,zGEO,BxGEO,ByGEO,BzGEO)
C
       IMPLICIT NONE
C
       REAL*8 xGEO,yGEO,zGEO
       REAL*8 Bo
       REAL*8 BxGEO,ByGEO,BzGEO
       REAL*8 rr,tt,pp
       REAL*8 Br,Bt,Bp
       REAL*8 xc,yc,zc
       REAL*8 ct,st,cp,sp
       REAL*8     pi,rad
       common /rconst/rad,pi
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C
       rr = SQRT(xGEO*xGEO+yGEO*yGEO+zGEO*zGEO)
       tt = ACOS(zGEO/rr)
       pp = ATAN2(yGEO,xGEO)
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
       BzGEO = Br*COS(tt) - Bt*SIN(tt)
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   Init centered dipole !!! BE CAREFUL : need to be called after
C                            INIT_DTD
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE INIT_CD
C
       IMPLICIT NONE
C
       REAL*8    xc,yc,zc                  !Re
       REAL*8    ct,st,cp,sp
       REAL*8    Bo
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C
C  Centered dipole !
       xc = 0.d0
       yc = 0.d0
       zc = 0.d0
       ct = 1.d0
       st = 0.d0
       cp = 1.d0
       sp = 0.d0
C
       RETURN
       END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

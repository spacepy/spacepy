c
c --------------------------------------------------------------------
c
        SUBROUTINE geo2dmag1(iyr,xGEO,xDMAG)
	    INTEGER*4 iyr
	    REAL*8    dyear
	    REAL*8    xGEO(3),xDMAG(3)

	    dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL GEO_DMAG(xGEO,xDMAG)
        end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE GEO_DMAG(xGEO,xDMAG)
C
       IMPLICIT NONE
C
       REAL*8    xGEO(3)
       REAL*8    xGEO_DEC(3)
       REAL*8    xDMAG(3)
       REAL*8    xc,yc,zc                  !Re
       REAL*8    ct,st,cp,sp
       REAL*8    Bo
C
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
C
       xGEO_DEC(1)= xGEO(1) - xc
       xGEO_DEC(2)= xGEO(2) - yc
       xGEO_DEC(3)= xGEO(3) - zc
C
       xDMAG(1) =  xGEO_DEC(1)*ct*cp + xGEO_DEC(2)*ct*sp
     & - xGEO_DEC(3)*st
       xDMAG(2) = -xGEO_DEC(1)*sp    + xGEO_DEC(2)*cp
       xDMAG(3) =  xGEO_DEC(1)*st*cp + xGEO_DEC(2)*st*sp
     & + xGEO_DEC(3)*ct
C
       RETURN
       END
C

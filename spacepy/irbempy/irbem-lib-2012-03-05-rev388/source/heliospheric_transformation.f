!***************************************************************************************************
! Copyright 2006 S. Bourdarie
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
! CREATION: S. Bourdarie - ONERA-DESP
!
! FILE CONTENT: 
!               SUBROUTINE hae2heeq1: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEEQ (Heliocentric Earth Equatorial)
!               SUBROUTINE HAE_HEEQ: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEEQ (Heliocentric Earth Equatorial)
!               SUBROUTINE heeq2hae1: Coordinate transformation from HEEQ (Heliocentric Earth Equatorial) to HAE (Heliocentric Aries Ecliptic)
!               SUBROUTINE HEEQ_HAE: Coordinate transformation from HEEQ (Heliocentric Earth Equatorial) to HAE (Heliocentric Aries Ecliptic)
!               SUBROUTINE hae2hee1: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEE (Heliocentric Earth Ecliptic)
!               SUBROUTINE HAE_HEE: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEE (Heliocentric Earth Ecliptic)
!               SUBROUTINE hee2hae1: Coordinate transformation from HEE (Heliocentric Earth Ecliptic) to HAE (Heliocentric Aries Ecliptic)
!               SUBROUTINE HEE_HAE: Coordinate transformation from HEE (Heliocentric Earth Ecliptic) to HAE (Heliocentric Aries Ecliptic)
!               SUBROUTINE gse2hee1: Coordinate transformation from GSE to HEE (Heliocentric Earth Ecliptic)
!               SUBROUTINE hee2gse1: Coordinate transformation from HEE to GSE (Heliocentric Earth Ecliptic)
!               SUBROUTINE GSE_HEE: Coordinate transformation from GSE to HEE (Heliocentric Earth Ecliptic)
!               SUBROUTINE SUN2: Calculate sideral time and position of the SUN
!
!***************************************************************************************************
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEEQ (Heliocentric Earth Equatorial)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xHAE -> position in HAE system in AU (double array(3))
!
! OUTPUT: xHEEQ -> position in HEEQ system in AU (double array(3))
!
! CALLING SEQUENCE: hae2heeq1(iyr,idoy,secs,xHAE,xHEEQ)
!---------------------------------------------------------------------------------------------------
c
        SUBROUTINE hae2heeq1(iyr,idoy,secs,xHAE,xHEEQ)
c
        IMPLICIT NONE
c
	INTEGER*4 iyr,idoy
	REAL*8    secs
	REAL*8    xHAE(3),xHEEQ(3)
	REAL*8    T0,Lamda0,st,ct,sl,cl,so,co
c
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
C
        CALL SUN2(iyr,idoy,secs,T0,Lamda0,st,ct,sl,cl,so,co)
        CALL HAE_HEEQ(xHAE,xHEEQ)
       	RETURN
	end
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEEQ (Heliocentric Earth Equatorial)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!              CALL SUN2 before calling HAE_HEEQ
!
! INPUT: HAE -> position in HAE system in AU (double array(3))
!
! OUTPUT: xHEEQ -> position in HEEQ system in AU (double array(3))
!
! CALLING SEQUENCE: HAE_HEEQ(xHAE,xHEEQ)
!---------------------------------------------------------------------------------------------------
	SUBROUTINE HAE_HEEQ(xHAE,xHEEQ)
c
        IMPLICIT NONE
c
	REAL*8 xHAE(3),xHEEQ(3)
	REAL*8 T0,Lamda0,st,ct,sl,cl,so,co
c	
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
c
	xHEEQ(1) =  (ct*co-st*cl*so)*xHAE(1) + (ct*so+st*cl*co)*xHAE(2)
     &              + st*sl*xHAE(3)
	xHEEQ(2) = (-st*co-ct*cl*so)*xHAE(1) + (ct*cl*co-st*so)*xHAE(2)
     &              + ct*sl*xHAE(3)
        xHEEQ(3) =  sl*so*xHAE(1)            - sl*co*xHAE(2)           
     &              + cl*xHAE(3)
c
       	RETURN
	END
C
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HEEQ (Heliocentric Earth Equatorial) to HAE (Heliocentric Aries Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xHEEQ -> position in HEEQ system in AU (double array(3))
!
! OUTPUT: xHAE -> position in HAE system in AU (double array(3))
!
! CALLING SEQUENCE: heeq2hae1(iyr,idoy,secs,xHEEQ,xHAE)
!---------------------------------------------------------------------------------------------------
c
        SUBROUTINE heeq2hae1(iyr,idoy,secs,xHEEQ,xHAE)
c
        IMPLICIT NONE
c
	INTEGER*4 iyr,idoy
	REAL*8    secs
	REAL*8    xHAE(3),xHEEQ(3)
	REAL*8    T0,Lamda0,st,ct,sl,cl,so,co
c
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
C
        CALL SUN2(iyr,idoy,secs,T0,Lamda0,st,ct,sl,cl,so,co)
        CALL HEEQ_HAE(xHEEQ,xHAE)
       	RETURN
	end
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HEEQ (Heliocentric Earth Equatorial) to HAE (Heliocentric Aries Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!              CALL SUN2 before calling HEEQ_HAE
!
! INPUT: HEEQ -> position in HEEQ system in AU (double array(3))
!
! OUTPUT: xHAE -> position in HAE system in AU (double array(3))
!
! CALLING SEQUENCE: HEEQ_HAE(xHEEQ,xHAE)
!---------------------------------------------------------------------------------------------------
	SUBROUTINE HEEQ_HAE(xHEEQ,xHAE)
c
        IMPLICIT NONE
c
	REAL*8 xHAE(3),xHEEQ(3)
	REAL*8 T0,Lamda0,st,ct,sl,cl,so,co
c	
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
c
	xHAE(1) =  (ct*co-st*cl*so)*xHEEQ(1) + (-st*co-ct*cl*so)*xHEEQ(2)
     &              + so*sl*xHEEQ(3)
	xHAE(2) = (ct*so+st*cl*co)*xHEEQ(1) + (ct*cl*co-st*so)*xHEEQ(2)
     &              - co*sl*xHEEQ(3)
        xHAE(3) =  sl*st*xHEEQ(1)            + sl*ct*xHEEQ(2)           
     &              + cl*xHEEQ(3)
c
       	RETURN
	END
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEE (Heliocentric Earth Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xHAE -> position in HAE system in AU (double array(3))
!
! OUTPUT: xHEE -> position in HEE system in AU (double array(3))
!
! CALLING SEQUENCE: hae2hee1(iyr,idoy,secs,xHAE,xHEE)
!---------------------------------------------------------------------------------------------------
c
        SUBROUTINE hae2hee1(iyr,idoy,secs,xHAE,xHEE)
	INTEGER*4 iyr,idoy
	REAL*8    secs
	REAL*8    xHAE(3),xHEE(3)
	REAL*8    T0,Lamda0,st,ct,sl,cl,so,co
c
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
C
        CALL SUN2(iyr,idoy,secs,T0,Lamda0,st,ct,sl,cl,so,co)
        CALL HAE_HEE(xHAE,xHEE)
       	RETURN
	end
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HAE (Heliocentric Aries Ecliptic) to HEE (Heliocentric Earth Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!              CALL SUN2 before calling HAE_HEE
!
! INPUT: HAE -> position in HAE system in AU (double array(3))
!
! OUTPUT: xHEE -> position in HEE system in AU (double array(3))
!
! CALLING SEQUENCE: HAE_HEE(xHAE,xHEE)
!---------------------------------------------------------------------------------------------------
	SUBROUTINE HAE_HEE(xHAE,xHEE)
c
        IMPLICIT NONE
c
	REAL*8 xHAE(3),xHEE(3)
	REAL*8 T0,Lamda0,rad,angle,st,ct,sl,cl,so,co
c	
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
	DATA rad /57.29577951308D0/
c
        angle=(Lamda0+180.0D0)/rad
	xHEE(1) =  COS(angle)*xHAE(1) + SIN(angle)*xHAE(2)
	xHEE(2) = -SIN(angle)*xHAE(1) + COS(angle)*xHAE(2)
        xHEE(3) =  xHAE(3)
c
       	RETURN
	END
C
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HEE (Heliocentric Earth Ecliptic) to HAE (Heliocentric Aries Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xHEE -> position in HEE system in AU (double array(3))
!
! OUTPUT: xHAE -> position in HAE system in AU (double array(3))
!
! CALLING SEQUENCE: hee2hae1(iyr,idoy,secs,xHEE,xHAE)
!---------------------------------------------------------------------------------------------------
c
        SUBROUTINE hee2hae1(iyr,idoy,secs,xHEE,xHAE)
	INTEGER*4 iyr,idoy
	REAL*8    secs
	REAL*8    xHAE(3),xHEE(3)
	REAL*8    T0,Lamda0,st,ct,sl,cl,so,co
c
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
C
        CALL SUN2(iyr,idoy,secs,T0,Lamda0,st,ct,sl,cl,so,co)
        CALL HEE_HAE(xHEE,xHAE)
       	RETURN
	end
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HEE (Heliocentric Earth Ecliptic) to HAE (Heliocentric Aries Ecliptic) 
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!              CALL SUN2 before calling HEE_HAE
!
! INPUT: xHEE -> position in HEE system in AU (double array(3))
!
! OUTPUT: xHAE -> position in HAE system in AU (double array(3))
!
! CALLING SEQUENCE: HEE_HAE(xHEE,xHAE)
!---------------------------------------------------------------------------------------------------
	SUBROUTINE HEE_HAE(xHEE,xHAE)
c
        IMPLICIT NONE
c
	REAL*8 xHAE(3),xHEE(3)
	REAL*8 T0,Lamda0,rad,angle,st,ct,sl,cl,so,co
c	
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
	DATA rad /57.29577951308D0/
c
        angle=(Lamda0+180.0D0)/rad
	xHAE(1) =  COS(angle)*xHEE(1) - SIN(angle)*xHEE(2)
	xHAE(2) =  SIN(angle)*xHEE(1) + COS(angle)*xHEE(2)
        xHAE(3) =  xHEE(3)
c
       	RETURN
	END
C
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from GSE to HEE (Heliocentric Earth Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xGSE -> position in GSE system in Re (double array(3))
!
! OUTPUT: xHEE -> position in HEE system in AU (double array(3))
!
! CALLING SEQUENCE: gse2hee1(iyr,idoy,secs,xGSE,xHEE)
!---------------------------------------------------------------------------------------------------
c
        SUBROUTINE gse2hee1(iyr,idoy,secs,xGSE,xHEE)
	INTEGER*4 iyr,idoy
	REAL*8    secs
	REAL*8    xGSE(3),xHEE(3)
        REAL*8    ERA,AQUAD,BQUAD
	REAL*8    T0,Lamda0,st,ct,sl,cl,so,co
c
        COMMON/GENER/ERA,AQUAD,BQUAD
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
C
        CALL INITIZE
	do i=1,3
	   xGSE(i)=xGSE(i)*ERA
	enddo
        CALL SUN2(iyr,idoy,secs,T0,Lamda0,st,ct,sl,cl,so,co)
        CALL GSE_HEE(xGSE,xHEE)
	do i=1,3
	   xHEE(i)=xHEE(i)/1.495985D8
	enddo
	end
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from HEE to GSE (Heliocentric Earth Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!
! INPUT: iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xHEE -> position in HEE system in AU (double array(3))
!
! OUTPUT: xGSE -> position in GSE system in Re (double array(3))
!
! CALLING SEQUENCE: hee2gse1(iyr,idoy,secs,xHEE,xGSE)
!---------------------------------------------------------------------------------------------------
c
        SUBROUTINE hee2gse1(iyr,idoy,secs,xHEE,xGSE)
	INTEGER*4 iyr,idoy
	REAL*8    secs
	REAL*8    xGSE(3),xHEE(3)
        REAL*8    ERA,AQUAD,BQUAD
	REAL*8    T0,Lamda0,st,ct,sl,cl,so,co
c
        COMMON/GENER/ERA,AQUAD,BQUAD
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
C
	do i=1,3
	   xHEE(i)=xHEE(i)*1.495985D8
	enddo
        CALL SUN2(iyr,idoy,secs,T0,Lamda0,st,ct,sl,cl,so,co)
        CALL GSE_HEE(xHEE,xGSE)
        CALL INITIZE
	do i=1,3
	   xGSE(i)=xGSE(i)/ERA
	enddo
	end
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from GSE to HEE (Heliocentric Earth Ecliptic)
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!              CALL SUN2 before calling GSE_HEE
!
! INPUT: xGSE -> position in GSE system in km (double array(3))
!
! OUTPUT: xHEE -> position in HEE system in km (double array(3))
!
! CALLING SEQUENCE: GSE_HEE(xGSE,xHEE)
!---------------------------------------------------------------------------------------------------
	SUBROUTINE GSE_HEE(xGSE,xHEE)
c
        IMPLICIT NONE
c
	REAL*8 xGSE(3),xHEE(3)
	REAL*8 T0,Lamda0,rad,st,ct,sl,cl,so,co
	REAL*8 r0,e,nu
c	
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
	DATA rad /57.29577951308D0/
C
c
        r0=1.495985D8   !km 1AU
	e=0.016709D0-0.0000418D0*T0
	nu=(Lamda0-282.94D0-1.72D0*T0)/rad
c
	xHEE(1) =  -xGSE(1) + r0*(1.D0-e*e)/(1.D0+e*COS(nu))
	xHEE(2) = - xGSE(2)
        xHEE(3) =  xGSE(3)
c
        
	RETURN
	END
C
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: None
!
! DESCRIPTION: Calculate sideral time and position of the SUN
!              From M.A. Hapgood, Space physics coordinate transformations: a user guide
!              Planet. Space Sci., Vol 40, pp 711-717, 1992
!
! INPUT: iyr -> year (long integer)
!        iday -> day of year (long integer)
!        secs -> UT in seconds (double)
!
! OUTPUT: T0 -> Time in Julian centuries (36525 days) from 12:00 UT on 1st January 2000 (double)
!         Lamda0 -> The Sun's ecliptic longitude in deg. (double)
!
! CALLING SEQUENCE: SUN2(iyr,iday,secs,T0,Lamda0)
!---------------------------------------------------------------------------------------------------
	SUBROUTINE SUN2(iyr,iday,secs,T0,Lamda0,st,ct,sl,cl,so,co)
C
	IMPLICIT NONE
C
C program to calculate sideral, time and position of the Sun
C good for years 1901 through 2099. accuracy 0.006 degree
c
c From M.A. Hapgood, Space physics coordinate transformations: a user guide
c Planet. Space Sci., Vol 40, pp 711-717, 1992
C
        INTEGER*4 JULDAY
	INTEGER*4 iyr,iday,dom,month,MJD
	REAL*8    secs,T0
	REAL*8    M,Lamda,Lamda0,rad
	REAL*8    omega,l,theta,st,ct,sl,cl,so,co
C
	DATA rad /57.29577951308D0/
c
	IF (iyr.LT.1901 .OR. iyr.GT.2099) RETURN
C
        dom=1
	month=1
        MJD=JULDAY(month,dom,iyr)+(iday-1)-2400001 !2400001=julday(11,17,1858)
	T0=(MJD-51544.5D0)/36525.0D0   !equation (3)
	M=357.528D0+35999.05D0*T0+0.04107D0*secs
	Lamda=280.46D0+36000.772D0*T0+0.04107D0*secs
	Lamda0=Lamda+(1.915D0-0.0048D0*T0)*SIN(M)+0.02D0*SIN(2.D0*M)  !equation (5)
c
	omega=(73.6667D0+0.013958D0*(MJD+3242.0D0)/365.25D0)/rad
	l=7.25D0/rad
	theta=ATAN(COS(l)*TAN(Lamda0/rad-omega))
	st=SIN(theta)
	ct=COS(theta)
	sl=SIN(l)
	cl=COS(l)
	so=SIN(omega)
	CO=COS(omega)
c
        RETURN
        END
C

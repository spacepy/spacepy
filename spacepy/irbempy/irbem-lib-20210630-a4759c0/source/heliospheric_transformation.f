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
        CALL SUN2(iyr,idoy,secs, T0,  Lamda0,st,ct,sl,cl,so,co)
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
! c	xHEEQ(1) =  (ct*co-st*cl*so)*xHAE(1) + (ct*so+st*cl*co)*xHAE(2)	!old
! 	xHEEQ(1) =  (co*ct-so*cl*st)*xHAE(1) + (-co*st-so*cl*ct)*xHAE(2)
! c	&              + st*sl*xHAE(3)	!old
!      &              + so*sl*xHAE(3)
! c	xHEEQ(2) = (-st*co-ct*cl*so)*xHAE(1) + (ct*cl*co-st*so)*xHAE(2)	!old
! 	xHEEQ(2) = (so*ct+co*cl*st)*xHAE(1) + (-so*st+co*cl*ct)*xHAE(2)
! c    &              + ct*sl*xHAE(3)	!old
!      &              + co*sl*xHAE(3)
! c      xHEEQ(3) =  sl*so*xHAE(1)            - sl*co*xHAE(2)	!old
!        xHEEQ(3) =  sl*st*xHAE(1)            + sl*ct*xHAE(2)
! c    &              + cl*xHAE(3)	!old
!      &              + cl*xHAE(3)
! c

	xHEEQ(1) =  (co*ct-so*cl*st)*xHAE(1) + (so*ct+co*cl*st)*xHAE(2)
     &              + sl*st*xHAE(3)

	xHEEQ(2) = (-co*st-so*cl*ct)*xHAE(1) + (-so*st+co*cl*ct)*xHAE(2)
     &              + sl*ct*xHAE(3)

       xHEEQ(3) =  so*sl*xHAE(1)            - co*sl*xHAE(2)
     &              + cl*xHAE(3)

 
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
        CALL SUN2(iyr,idoy,secs, T0, Lamda0,st,ct,sl,cl,so,co)
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
        CALL SUN2(iyr,idoy,secs, T0, Lamda0,st,ct,sl,cl,so,co)
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
        !angle=(Lamda0+180.0D0)/rad
        angle=(Lamda0+180)/rad		!Kellerman (this angle is simply a rotation for aries (eq point) to earth sun line (not to ecliptic long of sun as stated, checked with position and output of Lamda0
	xHEE(1) = COS(angle)*xHAE(1) + SIN(angle)*xHAE(2)
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
        SUBROUTINE hee2hae1(iyr,idoy,secs,xHEE,xHAE,TT)
	INTEGER*4 iyr,idoy
	REAL*8    secs
	REAL*8    xHAE(3),xHEE(3)
	REAL*8    T0,Lamda0,st,ct,sl,cl,so,co
c
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co
C
        CALL SUN2(iyr,idoy,secs, T0, Lamda0,st,ct,sl,cl,so,co)
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
	REAL*8 T0,Lamda0,rad,angle,st,ct,sl,cl,so,co, ang
	REAL*8 pi
	 
        COMMON /sunMAH/T0,Lamda0,st,ct,sl,cl,so,co        
	DATA rad /57.29577951308D0/
		pi = ACOS(0.0D0)*2.0D0
		angle=-(Lamda0+180)			!Kellerman
		angle=angle/rad		
	xHAE(1) =  COS(angle)*xHEE(1) + SIN(angle)*xHEE(2)
	xHAE(2) =  -SIN(angle)*xHEE(1) + COS(angle)*xHEE(2)
        xHAE(3) =  xHEE(3)
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
        CALL SUN2(iyr,idoy,secs, T0, Lamda0,st,ct,sl,cl,so,co)
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
        CALL SUN2(iyr,idoy,secs, T0, Lamda0,st,ct,sl,cl,so,co)
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
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!
! OUTPUT: T0 -> Time in Julian centuries (36525 days) from 12:00 UT on 1st January 2000 (double)
!         Lamda0 -> The Sun's ecliptic longitude in deg. (double)
!
! CALLING SEQUENCE: SUN2(iyr,idoy,secs,T0,Lamda0)
!---------------------------------------------------------------------------------------------------
	SUBROUTINE SUN2(iyr,idoy,secs,T0, Lamda0,st,ct,sl,cl,so,co)
C
	IMPLICIT NONE
C
C program to calculate sideral, time and position of the Sun
C good for years 1901 through 2099. accuracy 0.006 degree
c
c From M.A. Hapgood, Space physics coordinate transformations: a user guide
c Planet. Space Sci., Vol 40, pp 711-717, 1992
c last change A. C. Kellerman 18/02/2016
c added epoch 2000 functionality and corrected errors in the Hapgood
c formulae, i.e. hours input and errors with deg -> sin,tan functions
C
	REAL*8 QCDFTDB
	INTEGER*4 JULDAY
	INTEGER*4 iyr,idoy, dom, MJD, NST, j
	INTEGER*4    month, day, hour, mt, sc
	REAL*8 secs, T0, TT, dj, dlm, varpi, ggeo, lgeo, rgeo, lon
	REAL*8 M,Lamda,Lamda0,rad,aberr
	REAL*8 omega,l,theta,st,ct,sl,cl,so,co, londif
	REAL*8 dpi, dpi2
	 
	dpi = ACOS(0.0d0)*2.0D0
	dpi2 = 2.0D0*dpi
C
	DATA rad /57.29577951308D0/ !180 /pi
c
	IF (iyr.LT.1901 .OR. iyr.GT.2099) RETURN
C
!         dom=1
! 		month=1
!        MJD=JULDAY(month,dom,iyr)+(idoy-1)-2400001 !2400001=julday(11,17,1858)    old

	CALL DOY_AND_UT2DATE_AND_TIME(iyr,idoy,secs,month,day,hour,mt,sc)  
	TT=QCDFTDB(iyr,month,day,hour,mt,sc,0)
	!T0=(MJD-51544.5D0)/36525.0D0   !equation (3)
	
	T0=TT*1d0/(8.64d6*365.25d0)
	dj=T0*36525d0	
! 	M=357.528D0+35999.05D0*T0+0.04107D0*secs ! this is incorrect 
        !the equation is based on hours (AKellerman ref. Hapgood '92)
! 	Lamda=280.46D0+36000.772D0*T0+0.04107D0*secs	! as above
! 	Lamda0=Lamda+(1.915D0-0.0048D0*T0)*SIN(M/rad)+0.02D0*SIN(2.D0*M/rad)  !equation (5)


	dlm= ((100.4664568d0)+(mod (35999.3728565d0*T0, 360.0d0)))/rad
	varpi=((102.9373481d0)+(0.3225654d0*T0))/rad
	ggeo=dlm-varpi
	lgeo=dlm+(1.915d0*SIN(ggeo)+0.02d0*SIN(2.0d0*ggeo))/rad
	rgeo=(1.00014d0)-0.01617d0*COS(ggeo)-0.00014d0*COS(2.0d0*ggeo)
	lon=lgeo+dpi
	omega=((75.76d0)+1.397d0*T0)/rad
	l=7.25d0/rad

	aberr=20.0d0/3600.0d0/rad
	londif=mod (lon-aberr-omega, dpi2)
  
	Lamda0 = lon*rad
	
! 	omega=(73.6667D0+0.013958D0*(MJD+3242.0D0)/365.25D0)/rad
! 	l=7.25D0/rad
! 	londif=Lamda0/rad-omega
! 	londif=mod (londif, dpi2)					
	IF (londif.gt.dpi) londif=londif-dpi2		
	IF (londif.lt.-dpi) londif=londif+dpi2		
	
	theta=ATAN(COS(l)*TAN(londif))

	IF (theta.gt.dpi) theta=theta-dpi2		!Kellerman
	IF (theta.lt.-dpi) theta=theta+dpi2		!Kellerman

	IF (abs(theta-londif).lt.dpi/2.0D0) theta=theta+dpi			!Kellerman
	!WRITE(*,'((f15.5), i4, 2(f15.5))') Lamda0, idoy, theta*rad, omega*rad

	st=SIN(theta)
	ct=COS(theta)
	sl=SIN(l)
	cl=COS(l)
	so=SIN(omega)
	CO=COS(omega)	
        RETURN
        END
C

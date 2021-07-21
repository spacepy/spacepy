!***************************************************************************************************
! Copyright 1987, D. Bilitza, 2007 S. Bourdarie
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE IGRF(XXX,YYY,ZZZ,BXXX,BYYY,BZZZ)           
C-------------------------------------------------------------------
C CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
C REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61, 
C      1970.
C--------------------------------------------------------------------
C CHANGES (D. BILITZA, NOV 87):
C   - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
C   - CALCULATES DIPOL MOMENT
C--------------------------------------------------------------------
C  INPUT:  ENTRY POINT FELDG
C	 	GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
C         	GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
C         	ALT   ALTITUDE IN KM ABOVE SEA LEVEL
C
C	   ENTRY POINT FELDC
C		V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
C			X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
C			Y-AXIS POINTING TO EQUATOR AT 90 LONG.
C			Z-AXIS POINTING TO NORTH POLE
C
C	   COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED
C	     IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG.
C	
C	   COMMON /MODEL/ AND /GENER/
C		UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
C		ERA	EARTH RADIUS FOR NORMALIZATION OF CARTESIAN 
C			COORDINATES (6371.2 KM)
C		AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS FOR 
C			EARTH ELLIPSOID AS RECOMMENDED BY INTERNATIONAL 
C			ASTRONOMICAL UNION (6378.160, 6356.775 KM).
C		NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
C		TIME	YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC 
C			FIELD IS TO BE CALCULATED
C		G(M)	NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
C			M=NMAX*(NMAX+2)
C------------------------------------------------------------------------
C  OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS
C	   BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
C		  TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
C		  POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
C		  AND DOWNWARD.   
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4         NMAX,IHMAX,LAST,I,K,IH,M,IMAX,IL
      REAL*8            XXX,YYY,ZZZ
      REAL*8            BXXX,BYYY,BZZZ
      REAL*8            F,RQ,X,Y,Z,S,T
      REAL*8  		XI(3),H(144)
      REAL*8            G(144)
      COMMON/MODEL/	G  
C
C-- IS RECORDS ENTRY POINT
C
C*****ENTRY POINT  FELDG  TO BE USED WITH GEODETIC CO-ORDINATES         
      RQ=1.D0/(XXX*XXX+YYY*YYY+ZZZ*ZZZ) 
      XI(1)=XXX*RQ                                                      
      XI(2)=YYY*RQ                                                      
      XI(3)=ZZZ*RQ
C
      NMAX = 10                                                      
      IHMAX=NMAX*NMAX+1                                                 
      LAST=IHMAX+NMAX+NMAX                                              
      IMAX=NMAX+NMAX-1                                                  
      DO 8 I=IHMAX,LAST                                                 
8     H(I)=G(I)                                                         
      DO 6 K=1,3,2                                                      
      I=IMAX                                                            
      IH=IHMAX                                                          
1     IL=IH-I                                                           
      F=2.D0/((I-K+2)*1.D0)                                                 
      X=XI(1)*F                                                         
      Y=XI(2)*F                                                         
      Z=XI(3)*(F+F)                                                     
      I=I-2                                                             
      IF(I-1)5,4,2                                                      
2     DO 3 M=3,I,2                                                      
      H(IL+M+1)=G(IL+M+1)+Z*H(IH+M+1)+X*(H(IH+M+3)-H(IH+M-1))           
     A                               -Y*(H(IH+M+2)+H(IH+M-2))           
3     H(IL+M)=G(IL+M)+Z*H(IH+M)+X*(H(IH+M+2)-H(IH+M-2))                 
     A                         +Y*(H(IH+M+3)+H(IH+M-1))                 
4     H(IL+2)=G(IL+2)+Z*H(IH+2)+X*H(IH+4)-Y*(H(IH+3)+H(IH))             
      H(IL+1)=G(IL+1)+Z*H(IH+1)+Y*H(IH+4)+X*(H(IH+3)-H(IH))             
5     H(IL)=G(IL)+Z*H(IH)+2.*(X*H(IH+1)+Y*H(IH+2))                      
      IH=IL                                                             
      IF(I.GE.K)GOTO1                                                   
6     CONTINUE                                                          
      S=.5D0*H(1)+2.D0*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))                   
      T=(RQ+RQ)*SQRT(RQ)                                                
      BXXX=T*(H(3)-S*XXX)                                               
      BYYY=T*(H(4)-S*YYY)                                               
      BZZZ=T*(H(2)-S*ZZZ)                                               
      RETURN                                                            
      END                                                               
C

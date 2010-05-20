!***************************************************************************************************
! CREATION: S. Bourdarie - ONERA-DESP
!
! FILE CONTENT: 
!               SUBROUTINE coord_trans1: Generic coordinate transformation from one Earth or Heliospheric coordinate
!                                        system to another one
!               SUBROUTINE coord_trans_vec1: Generic coordinate transformation from one Earth or Heliospheric coordinate
!                                        system to another one (handle up to
!                                        100000 positions)
!
!***************************************************************************************************
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.0
!
! CREATION: S. Bourdarie - January 2007
! MODIFICATION: None
!
! DESCRIPTION: Coordinate transformation from one Earth or Heliospheric coordinate
!                                        system to another one
!
! INPUT: sysaxesIN -> designed input coordinate system (long integer)
!        sysaxesOUT> designed output coordinate system (long integer)
!        iyr -> year (long integer)
!        idoy -> day of year (long integer)
!        secs -> UT in seconds (double)
!        xIN -> position in input coordinate system (double array(3))
!
! OUTPUT: xOUT -> position in output coordinate system (double array(3))
!
! CALLING SEQUENCE: call coord_trans1(sysaxesIN,sysaxesOUT,iyr,idoy,secs,xIN,xOUT)
!---------------------------------------------------------------------------------------------------
c
        SUBROUTINE coord_trans1(sysaxesIN,sysaxesOUT,iyr,idoy,
     &   secs,xIN,xOUT)
c
        IMPLICIT NONE
c
	INTEGER*4 sysaxesIN,sysaxesOUT,iyr,idoy
	INTEGER*4 i
	REAL*8    secs,psi
	REAL*8    xIN(3),xOUT(3),xTMP(3),alti
 
        if (sysaxesIN.EQ.sysaxesOUT) then
	   write(6,*)'sysaxesIN = sysaxesOUT!'
	   do i=1,3
	      xOUT(i)=xIN(i)
	   enddo
	   return
	endif
        if (sysaxesIN.LT.0) then
	   write(6,*)'sysaxesIN out of range !'
	   do i=1,3
	      xOUT(i)=-1.D31
	   enddo
	   return
	endif
        if (sysaxesIN.GT.11) then
	   write(6,*)'sysaxesIN out of range !'
	   do i=1,3
	      xOUT(i)=-1.D31
	   enddo
	   return
	endif
        if (sysaxesOUT.LT.0) then
	   write(6,*)'sysaxesOUT out of range !'
	   do i=1,3
	      xOUT(i)=-1.D31
	   enddo
	   return
	endif
        if (sysaxesOUT.GT.11) then
	   write(6,*)'sysaxesOUT out of range !'
	   do i=1,3
	      xOUT(i)=-1.D31
	   enddo
	   return
	endif
c
c input=GDZ
	if (sysaxesIN.EQ.0) then
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call gdz_geo(xIN(2),xIN(3),xIN(1),xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=GEO
	if (sysaxesIN.EQ.1) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call geo_gdz(xIN(1),xIN(2),xIN(3),xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call geo2gsm1(iyr,idoy,secs,psi,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call geo2gse1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call geo2sm1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call geo2gei1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call geo2mag1(iyr,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call CAR_SPH(xIN,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call geo_gdz(xIN(1),xIN(2),xIN(3),xTMP(1),xTMP(2),xTMP(3))
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call geo2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call geo2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call geo2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=GSM
	if (sysaxesIN.EQ.2) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call gsm2sm1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call gsm2geo1(iyr,idoy,secs,psi,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=GSE
	if (sysaxesIN.EQ.3) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call gse2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call gse2geo1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call gse2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call gse2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call gse2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call gse2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call gse2geo1(iyr,idoy,secs,xIN,xTMP)
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call gse2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call gse2hee1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call gse2hee1(iyr,idoy,secs,xIN,xTMP)
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call gse2hee1(iyr,idoy,secs,xIN,xTMP)
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=SM
	if (sysaxesIN.EQ.4) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call sm2geo1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call sm2gsm1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call sm2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c input=GEI
	if (sysaxesIN.EQ.5) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call gei2geo1(iyr,idoy,secs,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call gei2geo1(iyr,idoy,secs,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c
c input=MAG
	if (sysaxesIN.EQ.6) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
              call mag2geo1(iyr,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
              call mag2geo1(iyr,xIN,xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
              call mag2geo1(iyr,xIN,xTMP)
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
              call mag2geo1(iyr,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
              call mag2geo1(iyr,xIN,xTMP)
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
              call mag2geo1(iyr,xIN,xTMP)
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
              call mag2geo1(iyr,xIN,xTMP)
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
              call mag2geo1(iyr,xIN,xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
              call mag2geo1(iyr,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
              call mag2geo1(iyr,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
              call mag2geo1(iyr,xIN,xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=SPH
	if (sysaxesIN.EQ.7) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3)
     &             ,xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      xOUT(1)=SQRT(xIN(1)*xIN(1)+xIN(2)*xIN(2)+xIN(3)*xIN(3))
	      xOUT(2)=xTMP(1)
	      xOUT(3)=xTMP(2)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call SPH_CAR(xIN(1),xIN(2),xIN(3),xTMP)
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=RLL
	if (sysaxesIN.EQ.8) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
              call RLL_GDZ(xIN(1),xIN(2),xIN(3),alti)
	      call gdz_geo(xIN(2),xIN(3),alti,xTMP(1),xTMP(2),xTMP(3))
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2hee1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2hae1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=HEE
	if (sysaxesIN.EQ.9) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call hee2gse1(iyr,idoy,secs,xIN,xTMP)
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      xOUT(1)=SQRT(xTMP(1)*xTMP(1)+xTMP(2)*xTMP(2)
     &                     +xTMP(3)*xTMP(3))
	      xOUT(3)=xOUT(2)
	      xOUT(2)=xOUT(1)
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call hee2hae1(iyr,idoy,secs,xIN,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call hee2hae1(iyr,idoy,secs,xIN,xTMP)
	      call hae2heeq1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	endif
c
c input=HAE
	if (sysaxesIN.EQ.10) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call hae2hee1(iyr,idoy,secs,xIN,xTMP) 
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      xOUT(1)=SQRT(xTMP(1)*xTMP(1)+xTMP(2)*xTMP(2)
     &                     +xTMP(3)*xTMP(3))
	      xOUT(3)=xOUT(2)
	      xOUT(2)=xOUT(1)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call hae2hee1(iyr,idoy,secs,xIN,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.11) then  !HEEQ
	      call hae2heeq1(iyr,idoy,secs,xIN,xOUT) 
	   endif
	endif
c
c input=HEEQ
	if (sysaxesIN.EQ.11) then
	   if (sysaxesOUT.EQ.0) then  !GDZ
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(2),xOUT(3),xOUT(1))
	   endif
	   if (sysaxesOUT.EQ.1) then  !GEO
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.2) then  !GSM
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gsm1(iyr,idoy,secs,psi,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.3) then  !GSE
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gse1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.4) then  !SM
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2sm1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.5) then  !GEI
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2gei1(iyr,idoy,secs,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.6) then  !MAG
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo2mag1(iyr,xTMP,xOUT)
	   endif
	   if (sysaxesOUT.EQ.7) then  !SPH
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call CAR_SPH(xTMP,xOUT(1),xOUT(2),xOUT(3))
	   endif
	   if (sysaxesOUT.EQ.8) then  !RLL
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call hee2gse1(iyr,idoy,secs,xTMP,xOUT) 
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call gse2geo1(iyr,idoy,secs,xTMP,xOUT)
	      do i=1,3
	         xTMP(i)=xOUT(i)
	      enddo
	      call geo_gdz(xTMP(1),xTMP(2),xTMP(3),
     &                     xOUT(1),xOUT(2),xOUT(3))
	      xOUT(1)=SQRT(xTMP(1)*xTMP(1)+xTMP(2)*xTMP(2)
     &                     +xTMP(3)*xTMP(3))
	      xOUT(3)=xOUT(2)
	      xOUT(2)=xOUT(1)
	   endif
	   if (sysaxesOUT.EQ.9) then  !HEE
	      call heeq2hae1(iyr,idoy,secs,xIN,xTMP) 
	      call hae2hee1(iyr,idoy,secs,xTMP,xOUT) 
	   endif
	   if (sysaxesOUT.EQ.10) then  !HAE
	      call heeq2hae1(iyr,idoy,secs,xIN,xOUT) 
	   endif
	endif
c
        end
C       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SUBROUTINE coord_trans_vec 
C      Routine to vectorize coord_trans1, the single-element
C      coordinate transformation utility.  
C      Calling sequence:
C           
C      call coord_trans_vec(ntime,sysaxesIN,sysaxesOUT,iyear,idoy,secs,
C     +          xIN,xOUT)
C
C         with the following types:
C      INTEGER*4 sysaxesIN, sysaxesOUT
C      INTEGER*4 iyear(nmax),idoy(nmax)
C      REAL*8 secs(nmax)
C      REAL*8 xINV(3,nmax), xOUTV(3,nmax)
C      INTEGER*4 numpoints
C
C     As with all (most?) onera library calls, the maximum array size
C     is limited to 100,000 elements
C                              Contributed by Timothy Guild, 2.2.07
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

        SUBROUTINE coord_trans_vec1(ntime,sysaxesIN,sysaxesOUT,
     &   iyear,idoy,secs,xINV,xOUTV)

      IMPLICIT NONE

      INTEGER*4 nmax,i,ntime, sysaxesIN, sysaxesOUT
      PARAMETER (nmax=100000)
      INTEGER*4 iyear(nmax),idoy(nmax),y,d
      REAL*8 secs(nmax),xINV(3,nmax),xOUTV(3,nmax)
      REAL*8 xIN(3),xOUT(3),s   ! local vars
	
      ! Loop over the number of points specified, calling 
      !  coord_trans1 each iteration
      DO i = 1,ntime
         y = iyear(i)
         d = idoy(i)
         s = secs(i)
         
	 xIN(1) = xINV(1,i) ! copy each array element into 3x1 array
	 xIN(2) = xINV(2,i)
         xIN(3) = xINV(3,i)

         call coord_trans1( sysaxesIN,sysaxesOUT,y,d,s,xIN,xOUT)

	 xOUTV(1,i) = xOUT(1)  ! copy back to 3xN array to pass back
	 xOUTV(2,i) = xOUT(2)
         xOUTV(3,i) = xOUT(3)

      ENDDO 

      END
	

!***************************************************************************************************
! Copyright 2007 T.P. O'Brien
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
      SUBROUTINE GET_HEMI1(kext,options,sysaxes,iyearsat,idoy,UT,
     &     xIN1,xIN2,xIN3,maginput,xHEMI)
c      computes magnetic hemisphere by taking a small 0.001 Re step
c      along field line and comparing B at the new point to that at the old
c      if the new B is larger, then in Northern Hemisphere (xHEMI=+1), else (xHEMI=-1)
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C     
      COMMON /magmod/k_ext,k_l,kint
c     declare inputs
      INTEGER*4    kext,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     maginput(25)

c     declare outputs
      integer*4     xHEMI

c     declare internal variables
      integer*4  i,k_ext,k_l,kint,ifail
      real*8     xGEO(3)
      real*8     alti,lati,longi
      REAL*8     BxGEO(3),B1,B2
      integer*4 int_field_select, ext_field_select ! functions to call
C     
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c     
      CALL INITIZE
      
      call init_fields ( kint, iyearsat, idoy, ut, options(2) )
      
      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3, 
     6     alti, lati, longi, xGEO )
      
      call set_magfield_inputs ( k_ext, maginput, ifail )
      
      xHEMI = 0                 ! default, indicates a failure
      
      if ( ifail.lt.0 ) RETURN
      if (k_ext .eq. 13 .or. k_ext .eq. 14) then
            call INIT_TS07D_TLPR
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
            if ( ifail.lt.0 ) RETURN
      end if

c     
      CALL CHAMP(xGEO,BxGEO,B1,Ifail)
      IF ((Ifail.LT.0).or.(B1.eq.baddata)) THEN
         return
      ENDIF
      
      do i=1,3
         xGEO(i) = xGEO(i)+BxGEO(i)/B1/1e3 ! scale so that dX is 0.001 Re
      enddo
      
      CALL CHAMP(xGEO,BxGEO,B2,Ifail)
      IF ((Ifail.LT.0).or.(B2.eq.baddata)) THEN
         return
      ENDIF
      
      if (B2.GT.B1) then
         xHEMI = +1
      else
         xHEMI = -1
      endif
      
      END
      

      SUBROUTINE GET_HEMI_MULTI(ntime,kext,options,sysaxes,iyearsat,
     & idoy,UT,xIN1,xIN2,xIN3,maginput,xHEMI)
c     calls get_hemi1 multiple times (ntime, <= ntime_max)
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      INCLUDE 'ntime_max.inc'   ! include file created by make, defines ntime_max

c     declare inputs
      INTEGER*4    ntime,kext,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat(ntime_max)
      integer*4    idoy(ntime_max)
      real*8     UT(ntime_max)
      real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
      real*8     maginput(25,ntime_max)

c     declare outputs
      integer*4     xHEMI(ntime_max)

c     declare internal variables
      integer*4  isat
c       TODO: implement the TS07D load flag here and uncomment this
c      if (kext .eq. 13) then
c        call INIT_TS07D_TLPR
c      end if
    
      do isat = 1,ntime
         call GET_HEMI1(kext,options,sysaxes,
     &    iyearsat(isat),idoy(isat),UT(isat),
     &    xIN1(isat),xIN2(isat),xIN3(isat),
     &    maginput(1,isat),xHEMI(isat))
      enddo

      end

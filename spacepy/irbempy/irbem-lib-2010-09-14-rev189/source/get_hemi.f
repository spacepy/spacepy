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
c     declare inputs
      INTEGER*4    kext,k_ext,k_l,kint,options(5)
      INTEGER*4    sysaxes,sysGEO
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
c     declare outputs
      integer*4     xHEMI

c     declare internal variables
      integer*4    i
      real*8     xGEO(3),xIN(3)
      real*8     maginput(25)
      REAL*8     BxGEO(3),Bl
      REAL*8     BxGEO2(3),B2
      
      sysGEO = 1 ! GEO
      call get_field1(kext,options,sysaxes,iyearsat,idoy,UT,
     *     xIN1,xIN2,xIN3,maginput,BxGEO,Bl)
      
      if (sysaxes.eq.sysGEO) then ! no need to convert
        xGEO(1) = xIN1
        xGEO(2) = xIN2
        xGEO(3) = xIN3
      else ! need to convert to GEO to do step
         xIN(1) = xIN1
         xIN(2) = xIN2
         xIN(3) = xIN3
         call coord_trans1(sysaxes,sysGEO,iyearsat,idoy,UT,xIN,xGEO)
      endif

!      Bl = Bl*1000.0            ! scale so that dX is 0.001 Re
      do i=1,3
         xGEO(i) = xGEO(i)+BxGEO(i)/Bl/1e3
      enddo

      call get_field1(kext,options,sysGEO,iyearsat,idoy,UT,
     *     xGEO(1),xGEO(2),xGEO(3),maginput,BxGEO2,B2)

      if (B2 .EQ. baddata) then
         xHEMI = 0
         return
      endif
      if (Bl .EQ. baddata) then
         xHEMI = 0
         return
      endif
c
      if (B2.GT.Bl) then
         xHEMI = +1
      else
         xHEMI = -1
      endif
      END
      

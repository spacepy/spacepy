!***************************************************************************************************
! Copyright 2011 T.P. O'Brien
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

      SUBROUTINE GET_Bderivs(ntime,kext,options,sysaxes,dX,
     & iyearsat,
     & idoy,UT,xIN1,xIN2,xIN3,maginput,Bgeo,Bmag,gradBmag,diffB)
C     computes derivatives of B (vector and magnitude)
C     inputs: ntime through maginput have the usual meaning, except dX
C     REAL*8 dX is the step size, in RE for the numerical derivatives (recommend 1E-3?)
C     real*8 Bgeo(3,ntime_max) - components of B in GEO, nT
C     real*8 Bmag(ntime_max) - magnitude of B in nT
C     real*8 gradBmag(3,ntime_max) - gradient of Bmag in GEO, nT/RE
C     real*8 diffB(3,3,ntime_max) - derivatives of Bgeo in GEO, nT/RE
C        diffB(i,j,t) = dB_i/dx_j for point t (t=1 to ntime)

      IMPLICIT NONE
      INCLUDE 'ntime_max.inc'   ! include file created by make, defines ntime_max
      INCLUDE 'variables.inc'
C
      COMMON /magmod/k_ext,k_l,kint

c     declare inputs
      INTEGER*4    ntime,kext,options(5)
      INTEGER*4    sysaxes
      REAL*8       dX
      INTEGER*4    iyearsat(ntime_max)
      integer*4    idoy(ntime_max)
      real*8     UT(ntime_max)
      real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
      real*8     maginput(25,ntime_max)

c     declare outputs
      real*8 Bgeo(3,ntime_max) ! components of B in GEO, nT
      real*8 Bmag(ntime_max) ! magnitude of B in nT
      real*8 gradBmag(3,ntime_max) ! gradient of Bmag in GEO, nT/RE
      real*8 diffB(3,3,ntime_max) ! derivatives of Bgeo in GEO, nT/RE

c     declare internal variables
      integer*4  isat
      integer*4  i,j,k_ext,k_l,kint,ifail
      real*8     xGEO(3),xGEOtmp(3)
      real*8     alti,lati,longi
      REAL*8     B1GEO(3),B1,BtmpGEO(3),Btmp
      integer*4 int_field_select, ext_field_select ! functions to call
C
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
      CALL INITIZE

      if (k_ext .eq. 13 .or. k_ext .eq. 14) then
        call INIT_TS07D_TLPR
      endif

      do isat = 1,ntime

         ! initialize outputs to baddata
         Bmag(isat) = baddata
         do i=1,3
            Bgeo(i,isat) = baddata
            gradBmag(i,isat) = baddata
            do j=1,3
               diffB(i,j,isat) = baddata
            enddo
         enddo

         call init_fields(kint,iyearsat(isat),idoy(isat),
     &        UT(isat),options(2))

         call get_coordinates (sysaxes,xIN1(isat),xIN2(isat),xIN3(isat),
     6        alti, lati, longi, xGEO )

         call set_magfield_inputs ( k_ext, maginput(1,isat), ifail)
         if (k_ext .eq. 13 .or. k_ext .eq. 14) then
            call INIT_TS07D_COEFFS(iyearsat(isat),idoy(isat),
     &      UT(isat),ifail)
         endif

         if (ifail.lt.0) goto 1000
c
         CALL CHAMP(xGEO,B1GEO,B1,Ifail)
         IF ((Ifail.LT.0).or.(B1.eq.baddata)) goto 1000

C copy start point to outputs
         Bmag(isat) = B1
         do i = 1,3
            Bgeo(i,isat) = B1GEO(i)
         enddo

         do j=1,3
C displace in dimension j
            do i = 1,3
               xGEOtmp(i)= xGEO(i)
            enddo
            xGEOtmp(j) = xGEOtmp(j)+dX
            CALL CHAMP(xGEOtmp,BtmpGEO,Btmp,Ifail) ! compute B at displace point
            IF ((Ifail.LT.0).or.(Btmp.eq.baddata)) THEN
               goto 1000
            ENDIF
c     compute derivatives
            gradBmag(j,isat) = (Btmp-B1)/dX
            do i = 1,3
               diffB(i,j,isat) = (BtmpGEO(i)-B1GEO(i))/dX
            enddo
         enddo ! end of j loop

 1000       continue            ! end of isat loop
         enddo

      end



      SUBROUTINE compute_grad_curv_curl(ntime,Bgeo,Bmag,gradBmag,diffB,
     & grad_par,grad_perp,grad_drift,curvature,Rcurv,curv_drift,
     & curlB,divB)
C     computes gradient factors, curvature factors, and curl of B

      IMPLICIT NONE
      INCLUDE 'ntime_max.inc'   ! include file created by make, defines ntime_max
      INCLUDE 'variables.inc'

C     all coordinates are in GEO reference frame
C     inputs: 
      integer*4 ntime           ! number of points
      real*8 Bgeo(3,ntime_max)  ! components of B in GEO, nT
      real*8 Bmag(ntime_max)    ! magnitude of B in nT
      real*8 gradBmag(3,ntime_max) ! gradient of Bmag in GEO, nT/RE
      real*8 diffB(3,3,ntime_max) ! derivatives of Bgeo in GEO, nT/RE
c     diffB(i,j,t) = dB_i/dx_j for point t (t=1 to ntime)
c     outputs:
      real*8 grad_par(ntime_max) ! gradient of Bmag along B nT/RE
      real*8 grad_perp(3,ntime_max) ! gradient of Bmag perpendicular to B nT/RE
      real*8 grad_drift(3,ntime_max) ! (bhat x grad_perp)/B, 1/RE (part of gradient drift velocity)
      real*8 curvature(3,ntime_max) ! (bhat dot grad)bhat, 1/RE (part of curvature force)
      real*8 Rcurv(ntime_max) ! 1/|curvature| RE (radius of curvature)
      real*8 curv_drift(3,ntime_max) ! (bhat x curvature), 1/RE (part of curvature drift)
      real*8 curlB(3,ntime_max) ! curl of B (nT/RE) (part of electrostatic current term)
      real*8 divB(ntime_max) ! divergence of B (nT/RE) (should be zero!)


c     internal variables
      integer*4 i,j,k,isat
      real*8 bhat(3), dbiBk ! unit vector along B and partial db_i/dB_k

      do isat = 1,ntime

         ! compute bhat and grad_par
         grad_par(isat) = 0.0
         do i = 1,3
            bhat(i) = Bgeo(i,isat)/Bmag(isat)
            grad_par(isat) = grad_par(isat) + gradBmag(i,isat)*bhat(i)
         enddo

         ! compute grad_perp
         do i = 1,3
            grad_perp(i,isat) = gradBmag(i,isat)-grad_par(isat)*bhat(i)
         enddo

         ! compute grad_drift
         grad_drift(1,isat) = (bhat(2)*grad_perp(3,isat) 
     &        - bhat(3)*grad_perp(2,isat)) / Bmag(isat)
         grad_drift(2,isat) = (bhat(3)*grad_perp(1,isat) 
     &        - bhat(1)*grad_perp(3,isat)) / Bmag(isat)
         grad_drift(3,isat) = (bhat(1)*grad_perp(2,isat) 
     &        - bhat(2)*grad_perp(1,isat)) / Bmag(isat)

         ! compute curvature (via chain rule on b = Bgeo/Bmag)
         do i = 1,3
            curvature(i,isat) = 0.0
            do k = 1,3
               ! compute db_i/dB_k
               if (i.eq.k) then
                  dbiBk = Bmag(isat)**2 - Bgeo(i,isat)**2
               else
                  dbiBk = -Bgeo(i,isat)*Bgeo(k,isat)
               endif
               dbiBk = dbiBk/Bmag(isat)**3
               do j = 1,3
                  curvature(i,isat) = curvature(i,isat)
     &                 +bhat(j)*dbiBk*diffB(k,j,isat)
               enddo
            enddo
         enddo
         
         ! compute radius of curvature = 1/|curvature|
         Rcurv(isat) = 0.0
         do i = 1,3
            Rcurv(isat) = Rcurv(isat)+curvature(i,isat)**2
         enddo
         if (Rcurv(isat).gt.0) then
            Rcurv(isat) = 1.0/sqrt(Rcurv(isat))
         else
            Rcurv(isat) = baddata
         endif
         
         ! compute curv_drift
         curv_drift(1,isat) = (bhat(2)*curvature(3,isat) 
     &        - bhat(3)*curvature(2,isat))
         curv_drift(2,isat) = (bhat(3)*curvature(1,isat) 
     &        - bhat(1)*curvature(3,isat))
         curv_drift(3,isat) = (bhat(1)*curvature(2,isat) 
     &        - bhat(2)*curvature(1,isat))

         ! compute curl of B
         curlB(1,isat) = diffB(3,2,isat)-diffB(2,3,isat)
         curlB(2,isat) = diffB(1,3,isat)-diffB(3,1,isat)
         curlB(3,isat) = diffB(2,1,isat)-diffB(1,2,isat)

         ! compute divergence of B
         divB(isat) = diffB(1,1,isat)+diffB(2,2,isat)+diffB(3,3,isat)

      enddo ! end isat loop

      end


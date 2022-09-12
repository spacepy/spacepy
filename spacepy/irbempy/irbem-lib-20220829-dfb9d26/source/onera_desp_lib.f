!***************************************************************************************************
! Copyright 2004, 2008 S. Bourdarie
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
C-----------------------------------------------------------------------------
C Wrappers and procedures for ONERA_DESP_LIB
C-----------------------------------------------------------------------------


c function returns version of fortran source code
      SUBROUTINE IRBEM_FORTRAN_VERSION1(VERSN)
        INCLUDE 'fortran_version.inc' ! include file created by make
        INTEGER*4 VERSN
        VERSN = FORTRAN_VERSION
      END


c function returns release of fortran source code
      SUBROUTINE IRBEM_FORTRAN_RELEASE1(RLS)
        INCLUDE 'fortran_release.inc'
        CHARACTER*80 RLS
        RLS = FORTRAN_RELEASE
      END


c function returns maximum size of variables
      SUBROUTINE GET_IRBEM_NTIME_MAX1(ntime_max1)
        INCLUDE 'ntime_max.inc' ! include file created by make
        INTEGER*4 ntime_max1
        ntime_max1 = ntime_max
      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE make_lstar1(ntime,kext,options,sysaxes,iyearsat,
     &  idoy,UT,xIN1,xIN2,xIN3,maginput,Lm,Lstar,BLOCAL,BMIN,XJ,MLT)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      INCLUDE 'ntime_max.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,options(5)
        INTEGER*4    ntime,sysaxes
      INTEGER*4    iyearsat(ntime_max)
      integer*4    idoy(ntime_max)
      real*8     UT(ntime_max)
      real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
      real*8     maginput(25,ntime_max)
c
c Declare internal variables
      INTEGER*4    isat,iyear,kint,ifail
        INTEGER*4    t_resol,r_resol,Ilflag,Ilflag_old
      REAL*8     mlon,mlon1
      REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT
      real*8     alti,lati,longi
c
c Declare output variables
        REAL*8     BLOCAL(ntime_max),BMIN(ntime_max),XJ(ntime_max)
      REAL*8     MLT(ntime_max)
        REAL*8     Lm(ntime_max),Lstar(ntime_max)
C
      COMMON /magmod/k_ext,k_l,kint
        COMMON /flag_L/Ilflag
        DATA  xSUN /1.d0,0.d0,0.d0/
      integer*4 int_field_select, ext_field_select
C
      Ilflag=0
      Ilflag_old=Ilflag
      if (options(3).lt.0 .or. options(3).gt.9) options(3)=0
      t_resol=options(3)+1
      r_resol=options(4)+1
      k_l=options(1)

      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
      CALL INITIZE
      if (k_ext .eq. 13 .or. k_ext .eq. 14) then !TS07D tail par init, only need it once
          call INIT_TS07D_TLPR
      end if

      DO isat = 1,ntime

          call init_fields ( kint, iyearsat(isat), idoy(isat),
     6          ut(isat), options(2) )

          call get_coordinates ( sysaxes,
     6        xIN1(isat), xIN2(isat), xIN3(isat),
     6        alti, lati, longi, xGEO )

        if (xIN1(isat) .eq. baddata .and. xIN2(isat) .eq. baddata
     &  .and. xIN3(isat) .eq.baddata) then
           Lm(isat)=baddata
           Lstar(isat)=baddata
           XJ(isat)=baddata
           BLOCAL(isat)=baddata
           BMIN(isat)=baddata
           GOTO 99
        endif

          call set_magfield_inputs ( k_ext, maginput(1,isat), ifail )

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then !TS07D coeff init
            call INIT_TS07D_COEFFS(iyearsat(isat),idoy(isat),
     &      ut(isat),ifail)
        end if

c        if (alti .le. 50.) ifail=-10 ! removed by TPO, 5/31/2011 - why would we force fail for alt<50km?
          if ( ifail.lt.0 ) then
                Lm(isat)=baddata
              Lstar(isat)=baddata
              XJ(isat)=baddata
              BLOCAL(isat)=baddata
              BMIN(isat)=baddata
              GOTO 99
             endif
c
           CALL calcul_Lstar_opt(t_resol,r_resol,XGeo
     &     ,Lm(isat),Lstar(isat),XJ(isat),BLOCAL(isat),BMIN(isat))
           if (Ilflag_old .eq.1 .and. Lstar(isat).eq. Baddata) then
            Ilflag=0
              CALL calcul_Lstar_opt(t_resol,r_resol,xGeo
     &        ,Lm(isat),Lstar(isat),XJ(isat),BLOCAL(isat),BMIN(isat))
         endif
         Ilflag_old=Ilflag

99         continue

        if (ifail .eq. -10) then
          MLT(isat) = baddata
        else
           CALL GDZ_GEO(lati,longi,alti
     &     ,xGEO(1),xGEO(2),xGEO(3))
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon1)
           CALL GSM_GEO(xSUN,xGEO)
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon)
           MLT(isat) = (Mlon1 - Mlon)/15.d0 + 12.d0
           IF (MLT(isat).GE.24.d0) MLT(isat) = MLT(isat) - 24.d0
           IF (MLT(isat).LT.0.d0) MLT(isat) = MLT(isat) + 24.d0
        endif
      ENDDO

      END
c
c --------------------------------------------------------------------
c
        SUBROUTINE make_lstar_shell_splitting1(ntime,Nipa,kext,options,
     & sysaxes,iyearsat,idoy,UT,xIN1,xIN2,xIN3,
     &  alpha,maginput,Lm,Lstar,BLOCAL,BMIN,XJ,MLT)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
        INCLUDE 'ntime_max.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,options(5),Nalp,Nipa
      PARAMETER (Nalp=25)
        INTEGER*4    ntime,sysaxes
      INTEGER*4    iyearsat(ntime_max)
      integer*4    idoy(ntime_max)
      real*8     UT(ntime_max)
      real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
      real*8     alpha(Nalp)
      real*8     maginput(25,ntime_max)
c
c
c Declare internal variables
      INTEGER*4    isat,iyear,IPA,kint,ifail
        INTEGER*4    Ilflag,t_resol,r_resol
      REAL*8     mlon,mlon1,BL,BMIR(25),Bmin_tmp
      REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT
      REAL*8     xGEOp(3,25)
      real*8     alti,lati,longi
c
c Declare output variables
        REAL*8     BLOCAL(ntime_max,Nalp),BMIN(ntime_max)
      REAL*8     XJ(ntime_max,Nalp),MLT(ntime_max)
        REAL*8     Lm(ntime_max,Nalp),Lstar(ntime_max,Nalp)
C
      COMMON /magmod/k_ext,k_l,kint
        COMMON /flag_L/Ilflag
        DATA  xSUN /1.d0,0.d0,0.d0/
      integer*4 int_field_select, ext_field_select
C
      Ilflag=0
      k_ext=kext
      if (options(3).lt.0 .or. options(3).gt.9) options(3)=0
      t_resol=options(3)+1
      r_resol=options(4)+1
      k_l=options(1)

      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
        CALL INITIZE
        if (k_ext .eq. 13 .or. k_ext .eq. 14) then !TS07D tail par init, only need it once
            call INIT_TS07D_TLPR
        end if

       DO isat = 1,ntime
        call init_fields ( kint, iyearsat(isat), idoy(isat),
     6      ut(isat), options(2) )

          call get_coordinates ( sysaxes,
     6        xIN1(isat), xIN2(isat), xIN3(isat),
     6        alti, lati, longi, xGEO )

        call set_magfield_inputs ( k_ext, maginput(1,isat), ifail )

       if (k_ext .eq. 13 .or. k_ext .eq. 14) then 
            call INIT_TS07D_COEFFS(iyearsat(isat),idoy(isat),
     &      ut(isat),ifail)
       end if


        if ( ifail.lt.0 ) then
          DO IPA=1,Nipa
                Lm(isat,IPA)=baddata
                Lstar(isat,IPA)=baddata
                XJ(isat,IPA)=baddata
                BLOCAL(isat,IPA)=baddata
            ENDDO
            BMIN(isat)=baddata
            GOTO 99
          endif
c
c Compute Bmin assuming 90deg PA at S/C
         k_l=0
           IPA=1
           CALL calcul_Lstar_opt(t_resol,r_resol,xGEO
     &         ,Lm(isat,IPA),Lstar(isat,IPA),XJ(isat,IPA)
     &         ,BLOCAL(isat,IPA),BMIN(isat))
         k_l=options(1)
           CALL find_bm_nalpha(xGEO,nipa,alpha,BL,BMIR,xGEOp)
           DO IPA=1,Nipa
           IF (Bmir(ipa).NE.baddata) THEN
             Ilflag=0
               CALL calcul_Lstar_opt(t_resol,r_resol,xGEOp(1,ipa)
     &         ,Lm(isat,IPA),Lstar(isat,IPA),XJ(isat,IPA)
     &         ,BLOCAL(isat,IPA),BMIN_tmp)
             ELSE
             Lm(isat,IPA)=baddata
             Lstar(isat,IPA)=baddata
             XJ(isat,IPA)=baddata
             BLOCAL(isat,IPA)=baddata
           ENDIF
           ENDDO

99         continue

           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon1)
           CALL GSM_GEO(xSUN,xGEO)
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon)
           MLT(isat) = (Mlon1 - Mlon)/15.d0 + 12.d0
           IF (MLT(isat).GE.24.d0) MLT(isat) = MLT(isat) - 24.d0
           IF (MLT(isat).LT.0.d0) MLT(isat) = MLT(isat) + 24.d0

      ENDDO

      END
c
C-----------------------------------------------------------------------------
c
        SUBROUTINE Lstar_Phi1(ntime,whichinv,options,iyearsat,
     &  idoy,Lstar,Phi)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
        INCLUDE 'ntime_max.inc'
C
c declare inputs
        INTEGER*4    whichinv,options(5)
        INTEGER*4    ntime
      INTEGER*4    iyearsat(ntime_max)
      integer*4    idoy(ntime_max)
c
c Declare internal variables
      INTEGER*4    isat,kint
        REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
c
c Declare output variables
        REAL*8     Phi(ntime_max),Lstar(ntime_max)
C
        COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
        REAL*8     pi,rad
        common /rconst/rad,pi
      integer*4 int_field_select, ext_field_select
C

      kint = int_field_select ( options(5) )
c
        CALL INITIZE

        DO isat = 1,ntime
          call init_fields ( kint, iyearsat(isat),idoy(isat),-1.D0,
     6      options(2) )

         if (whichinv .EQ. 1) then !Lstar to Phi
            if (Lstar(isat) .NE. baddata) then
               Phi(isat)=2.D0*pi*Bo/Lstar(isat)
            else
               Phi(isat)=baddata
            endif
         else
            if (phi(isat) .NE. baddata) then
               Lstar(isat)=2.D0*pi*Bo/Phi(isat)
            else
               Lstar(isat)=baddata
            endif
         endif
      enddo
      end
c
c --------------------------------------------------------------------
c
        SUBROUTINE drift_shell1(kext,options,sysaxes,iyearsat,idoy,UT,
     &  xIN1,xIN2,xIN3,maginput,Lm,Lstar,BLOCAL,BMIN,XJ,posit,ind)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,options(5)
        INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     maginput(25)
c
c Declare internal variables
      INTEGER*4    isat,iyear,kint,ifail
      REAL*8     xGEO(3)
      real*8     alti,lati,longi
c
c Declare output variables
        INTEGER*4  ind(48)
        REAL*8     BLOCAL(1000,48),BMIN,XJ
        REAL*8     Lm,Lstar
      REAL*8     posit(3,1000,48)
C
      COMMON /magmod/k_ext,k_l,kint
      integer*4 int_field_select, ext_field_select
C
      k_l=options(1)
c
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
        CALL INITIZE

      call init_fields ( kint, iyearsat, idoy, ut, options(2) )

      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6    alti, lati, longi, xGEO )

      call set_magfield_inputs ( k_ext, maginput, ifail )

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then
            call INIT_TS07D_TLPR
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
       end if

      if ( ifail.lt.0 ) then
             Lm=baddata
             Lstar=baddata
             XJ=baddata
             BMIN=baddata
             RETURN
          endif
c
        CALL trace_drift_shell_opt(xGeo
     &     ,Lm,Lstar,XJ,BLOCAL,BMIN,
     &     posit,ind)
      END

c --------------------------------------------------------------------
c
        SUBROUTINE trace_field_line1(kext,options,sysaxes,iyearsat,
     &  idoy,UT,xIN1,xIN2,xIN3,maginput,
     &  Lm,BLOCAL,BMIN,XJ,posit,ind)
        IMPLICIT NONE
        INTEGER*4 kext,options(5),sysaxes,iyearsat,idoy
        REAL*8 UT,xIN1,xIN2,xIN3,maginput(25),R0
        REAL*8 Lm,BLOCAL(3000),BMIN,XJ,posit(3,3000)
        INTEGER*4 ind

        R0=1.D0
        call trace_field_line2_1(kext,options,sysaxes,iyearsat,idoy,
     &  UT,xIN1,xIN2,xIN3,maginput,R0,
     &  Lm,BLOCAL,BMIN,XJ,posit,ind)
        end

        SUBROUTINE trace_field_line2_1(kext,options,sysaxes,iyearsat,
     &  idoy,UT,xIN1,xIN2,xIN3,maginput,R0,
     &  Lm,BLOCAL,BMIN,XJ,posit,ind)
C - added R0 argument - radius of reference sphere
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,kint,options(5)
        INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     maginput(25)
        real*8     R0
c
c Declare internal variables
      INTEGER*4    isat,iyear,ifail
        INTEGER*4    i,j
      REAL*8     xGEO(3)
      real*8     alti,lati,longi
c
c Declare output variables
        INTEGER*4  ind
        REAL*8     BLOCAL(3000),BMIN,XJ
        REAL*8     Lm
      REAL*8     posit(3,3000)
C
      COMMON /magmod/k_ext,k_l,kint
      integer*4 int_field_select, ext_field_select
C
        do i=1,3
         do j=1,3000
           posit(i,j)=baddata
         enddo
      enddo
c
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
        CALL INITIZE

      call init_fields ( kint, iyearsat, idoy, ut, options(2) )

      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6    alti, lati, longi, xGEO )

      call set_magfield_inputs ( k_ext, maginput, ifail )
        if (k_ext .eq. 13 .or. k_ext .eq. 14) then
            call INIT_TS07D_TLPR
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
        end if

      if ( ifail.lt.0 ) then
              Lm=baddata
             XJ=baddata
             BMIN=baddata
             ind=0
             RETURN
          endif
c
        CALL field_line_tracing_opt2(xGeo,R0
     &     ,Lm,XJ,BLOCAL,BMIN,posit,ind)
      END
c
c --------------------------------------------------------------------
c
      SUBROUTINE trace_field_line_towards_earth1(kext,options,sysaxes
     &  ,iyearsat,idoy,UT,xIN1,xIN2,xIN3,maginput,ds,posit,ind)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c declare inputs
      INTEGER*4    kext,k_ext,k_l,kint,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT,ds
      real*8     xIN1,xIN2,xIN3
      real*8     maginput(25)
c
c Declare internal variables
      INTEGER*4    isat,iyear
      INTEGER*4    Ndays,activ,i,j
      REAL*8     xGEO(3)
      real*8     alti,lati,longi
c
c Declare output variables
      INTEGER*4  ind, ifail
      REAL*8     posit(3,3000)
      integer*4 int_field_select, ext_field_select
C
C
      do i=1,3
        do j=1,3000
          posit(i,j)=baddata
        enddo
      enddo
c
      k_ext=kext
      k_l=options(1)
      kint=options(5)

      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
        CALL INITIZE

      call init_fields ( kint, iyearsat, idoy, ut, options(2) )

      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6    alti, lati, longi, xGEO )

      call set_magfield_inputs ( k_ext, maginput, ifail )

      if (k_ext .eq. 13 .or. k_ext .eq. 14) then
        call INIT_TS07D_TLPR
        call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
      end if

      if ( ifail.lt.0 ) then
             ind=0
             RETURN
          endif
c
        CALL field_line_tracing_towards_Earth_opt(xGEO,ds,posit,ind)
      END

c --------------------------------------------------------------------
c
        SUBROUTINE find_mirror_point1(kext,options,sysaxes,iyearsat,
     &  idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,BLOCAL,BMIR,xGEO)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,kint,options(5)
        INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     alpha
      real*8     maginput(25)
c
c Declare internal variables
      INTEGER*4    isat,iyear,Iint,ifail
      REAL*8     xMAG(3)
      real*8     alti,lati,longi
      real*8     BxGEO(3),xGeop(3,25),myalpha(25)
      real*8     MyBmir(25)
c
c Declare output variables
        REAL*8     BLOCAL,xGEO(3),BMIR
C
      COMMON /magmod/k_ext,k_l,kint
      integer*4 int_field_select, ext_field_select
C
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
        CALL INITIZE

      call init_fields ( kint, iyearsat, idoy, ut, options(2) )

      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6    alti, lati, longi, xGEO )

      call set_magfield_inputs ( k_ext, maginput, ifail )

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then
            call INIT_TS07D_TLPR
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
        end if

      if ( ifail.lt.0 ) then
             xGEO(1)=baddata
             xGEO(2)=baddata
             xGEO(3)=baddata
             BLOCAL=baddata
             BMIR=baddata
             RETURN
          endif
c
        if (alpha.eq.90.0d0) then
          Iint=2 ! TPO: presume this sets internal field?
c          CALL CHAMP(Iint,xGEO,BxGEO,Blocal,Ifail) ! Iint is superfluous
          CALL CHAMP(xGEO,BxGEO,Blocal,Ifail)
        IF (Ifail.LT.0) THEN
          Blocal=baddata
          Bmir=baddata
          xGEO(1)=baddata
          xGEO(2)=baddata
          xGEO(3)=baddata
        ELSE
          BMIR=Blocal
        ENDIF
        RETURN
      endif
c
       myalpha(1)=alpha
        CALL find_bm_nalpha(xGeo,1,myalpha
     &     ,BLOCAL,myBMIR,xGEOp)
        xGEO(1)=xGeop(1,1)
        xGEO(2)=xGeop(2,1)
        xGEO(3)=xGeop(3,1)
        Bmir=MyBmir(1)
      END



c --------------------------------------------------------------------
c
        SUBROUTINE FIND_MAGEQUATOR1(kext,options,sysaxes,iyearsat
     &  ,idoy,UT,xIN1,xIN2,xIN3,maginput,BMIN,posit)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,kint,options(5)
        INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     maginput(25)
c
c Declare internal variables
      INTEGER*4    isat,iyear,ifail
      REAL*8     xGEO(3)
      real*8     alti,lati,longi
c
c Declare output variables
        REAL*8     BMIN
      REAL*8     posit(3)
C
      COMMON /magmod/k_ext,k_l,kint
      integer*4 int_field_select, ext_field_select
C
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
        CALL INITIZE

      call init_fields ( kint, iyearsat, idoy, ut, options(2) )

      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6    alti, lati, longi, xGEO )

      call set_magfield_inputs ( k_ext, maginput, ifail )

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then
            call INIT_TS07D_TLPR
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
        end if

      if ( ifail.lt.0 ) then
             posit(1)=baddata
             posit(2)=baddata
             posit(3)=baddata
             BMIN=baddata
             RETURN
          endif
c
c
        CALL loc_equator_opt(xGeo,BMIN,posit)
      END
c --------------------------------------------------------------------
c
      SUBROUTINE GET_FIELD1(kext,options,sysaxes,iyearsat,idoy,UT,
     &     xIN1,xIN2,xIN3,maginput,BxGEO,Bl)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c     declare inputs
      INTEGER*4    kext,k_ext,k_l,kint,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     maginput(25)
c
c     Declare internal variables
      INTEGER*4    isat,iyear,ifail
      REAL*8     xGEO(3)
      real*8     alti,lati,longi
c
c     Declare output variables
      REAL*8     BxGEO(3),Bl
C
      COMMON /magmod/k_ext,k_l,kint
      integer*4 int_field_select, ext_field_select
C
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c

      CALL INITIZE

      call init_fields ( kint, iyearsat, idoy, ut, options(2) )

      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6     alti, lati, longi, xGEO )

      call set_magfield_inputs ( k_ext, maginput, ifail )

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then !special script to read files and
            call INIT_TS07D_TLPR
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
        end if

      if ( ifail.lt.0 ) then
         Bl=baddata
         BxGEO(1)=baddata
         BxGEO(2)=baddata
         BxGEO(3)=baddata
         RETURN
      endif
c
      CALL CHAMP(xGEO,BxGEO,Bl,Ifail)
      IF (Ifail.LT.0) THEN
         BxGEO(1)=baddata
         BxGEO(2)=baddata
         BxGEO(3)=baddata
         Bl=baddata
      ENDIF
      END
C-----------------------------------------------------------------------------

      SUBROUTINE GET_FIELD_MULTI(ntime,kext,options,sysaxes,iyearsat,
     &     idoy,UT,xIN1,xIN2,xIN3,maginput,BxGEO,Bl)
C     Call get_field1 many times (ntime, in fact, up to ntime = ntime_max)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      INCLUDE 'ntime_max.inc'   ! include file created by make, defines ntime_max
C
c     declare inputs
      INTEGER*4    ntime
      INTEGER*4    kext,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat(ntime_max)
      integer*4    idoy(ntime_max)
      real*8     UT(ntime_max)
      real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
      real*8     maginput(25,ntime_max)

c     Declare output variables
      REAL*8     BxGEO(3,ntime_max),Bl(ntime_max)
C
c     Declare internal variables
      integer*4 isat
      INTEGER*4 k_ext,k_l,kint,ifail

c    For now we call this every time
c   TODO: create a load flag so that we skip repetative loads
c      if (kext.eq.13) then
c            call INIT_TS07D_TLPR
c      endif

      do isat = 1,ntime
         call GET_FIELD1(kext,options,sysaxes,iyearsat(isat),
     &        idoy(isat),UT(isat), xIN1(isat),xIN2(isat),
     &        xIN3(isat),maginput(1,isat),BxGEO(1,isat),Bl(isat))

      enddo
      end
c
c --------------------------------------------------------------------
c
        SUBROUTINE GET_MLT1(iyr,idoy,UT,xGEO,MLT)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c declare inputs
      INTEGER*4    iyr
      integer*4    idoy
      real*8     UT
      real*8     xGEO(3)
c
c Declare internal variables
        REAL*8     dyear
      REAL*8     psi,mlon
      REAL*8     xMAG(3),xSUN(3),rM,MLAT,Mlon1
      REAL*8     xTMP(3)
c
c Declare output variables
        REAL*8     MLT
C
        DATA  xSUN /1.d0,0.d0,0.d0/
C
      dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,ut,psi)
C
        CALL geo_mag(xGEO,xMAG)
        CALL car_sph(xMAG,rM,MLAT,Mlon1)
        CALL GSM_GEO(xSUN,xTMP)
        CALL geo_mag(xTMP,xMAG)
        CALL car_sph(xMAG,rM,MLAT,Mlon)
        MLT = (Mlon1 - Mlon)/15.d0 + 12.d0
        IF (MLT.GE.24.d0) MLT = MLT - 24.d0
        IF (MLT.LT.0.d0) MLT = MLT + 24.d0

        END
C-----------------------------------------------------------------------------
c
      integer*4 function int_field_select ( kint )
      integer*4 kint

      !write(6,*)kint
      IF (kint .lt. 0) THEN
         kint=0
         WRITE(6,*)
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)'Invalid internal field specification'
         WRITE(6,*)'Selecting IGRF'
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)
      ENDIF
      if (kint .gt. 5) THEN
         kint=0
         WRITE(6,*)
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)'Invalid internal field specification'
         WRITE(6,*)'Selecting IGRF'
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)
      ENDIF

      int_field_select = kint

      return
      end

      integer*4 function ext_field_select ( kext )

        integer*4 kext

      IF (kext .lt. 0) THEN
         k_ext=5
         WRITE(6,*)
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)'Invalid external field specification'
         WRITE(6,*)'Selecting Olson-Pfitzer (quiet)'
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)
      ENDIF
      if (kext .gt. 14) THEN
         k_ext=5
         WRITE(6,*)
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)'Invalid external field specification'
         WRITE(6,*)'Selecting Olson-Pfitzer (quiet)'
         WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*)
      ENDIF

      ext_field_select = kext

      return
      end


      subroutine init_fields ( kint, iyearsat, idoy, ut, opt2 )
      INTEGER*4    firstJanuary,lastDecember,Julday,currentdoy
      INTEGER*4    kint,opt2,iyearsat,iyear
      integer*4    idoy
      real*8     UT,dec_year
      REAL*8 tilt,psi
      INTEGER*4    a2000_iyear,a2000_imonth,a2000_iday
      REAL*8      a2000_ut
      COMMON /dip_ang/tilt
      REAL*8     pi,rad

      common /rconst/rad,pi
      COMMON /a2000_time/a2000_ut,a2000_iyear,a2000_imonth,a2000_iday

      iyear = 1800

      if (kint .eq. 2) CALL JensenANDCain1960
      if (kint .eq. 3) CALL GSFC1266
c         write(6,*)real(isat)*100./real(ntime), '% done'
c
      if (kint .le. 1 .or. kint .eq. 4 .or. kint .eq. 5) then
        if (opt2 .eq. 0) then
          if (iyearsat .ne. iyear) then
            iyear=iyearsat
            dec_year=iyear+0.5d0
            if (Kint .eq. 4) then
              call myOwnMagField_init(dec_year)
            else
              CALL INIT_DTD(dec_year)
            endif
            if (kint .eq. 5) CALL INIT_CD
          endif
        else
          if (iyearsat .ne. iyear .or.
     &          MOD(idoy*1.d0,opt2*1.d0) .eq. 0) THEN
            iyear=iyearsat
            firstJanuary=JULDAY(iyear,01,01)
            lastDecember=JULDAY(iyear,12,31)
            currentdoy=(idoy/opt2)*opt2
            if (currentdoy .eq. 0) currentdoy=1
            dec_year=iyear+currentdoy*1.d0/
     &             ((lastDecember-firstJanuary+1)*1.d0)
            if (Kint .eq. 4) then
              call myOwnMagField_init(dec_year)
            else
              CALL INIT_DTD(dec_year)
            endif
          endif
        endif
      endif
c
      if ( ut.ge.0.0 ) CALL INIT_GSM(iyearsat,idoy,UT,psi)
      tilt = psi/rad

      a2000_iyear=iyearsat
      firstJanuary=JULDAY(a2000_iyear,01,01)
      currentdoy=firstJanuary+idoy-1
      CALL CALDAT(currentdoy,a2000_iyear,a2000_imonth,a2000_iday)
      a2000_ut=UT

      return
      end

      subroutine get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6    alti, lati, longi, xGEO )

      integer*4 sysaxes
      real*8 xIN1, xIN2, xIN3
      real*8 alti, lati, longi, xGEO(3)
      real*8 xGSM(3), xGSE(3), xSM(3), xGEI(3), xMAG(3)

         if (sysaxes .EQ. 0) then
             alti=xIN1
             lati=xIN2
             longi=xIN3
               CALL GDZ_GEO(lati,longi,alti,xGEO(1),xGEO(2),xGEO(3))
         endif
         if (sysaxes .EQ. 1) then
             xGEO(1)=xIN1
             xGEO(2)=xIN2
             xGEO(3)=xIN3
             CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
         endif
         if (sysaxes .EQ. 2) then
             xGSM(1)=xIN1
             xGSM(2)=xIN2
             xGSM(3)=xIN3
             CALL GSM_GEO(xGSM,xGEO)
             CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
         endif
         if (sysaxes .EQ. 3) then
             xGSE(1)=xIN1
             xGSE(2)=xIN2
             xGSE(3)=xIN3
             CALL GSE_GEO(xGSE,xGEO)
             CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
         endif
         if (sysaxes .EQ. 4) then
             xSM(1)=xIN1
             xSM(2)=xIN2
             xSM(3)=xIN3
             CALL SM_GEO(xSM,xGEO)
             CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
         endif
         if (sysaxes .EQ. 5) then
             xGEI(1)=xIN1
             xGEI(2)=xIN2
             xGEI(3)=xIN3
             CALL GEI_GEO(xGEI,xGEO)
             CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
         endif
         if (sysaxes .EQ. 6) then
             xMAG(1)=xIN1
             xMAG(2)=xIN2
             xMAG(3)=xIN3
             CALL MAG_GEO(xMAG,xGEO)
             CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
         endif
         if (sysaxes .EQ. 7) then
             xMAG(1)=xIN1
             xMAG(2)=xIN2
             xMAG(3)=xIN3
             CALL SPH_CAR(xMAG(1),xMAG(2),xMAG(3),xGEO)
             CALL GEO_GDZ(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)
         endif
         if (sysaxes .EQ. 8) then
             xMAG(1)=xIN1
             lati=xIN2
             longi=xIN3
             CALL RLL_GDZ(xMAG(1),lati,longi,alti)
             CALL GDZ_GEO(lati,longi,alti,xGEO(1),xGEO(2),xGEO(3))
         endif

         return
         end

       subroutine set_magfield_inputs ( kext, maginput, ifail )
      INCLUDE 'variables.inc'
      COMMON /index/activ
      integer*4 activ
      COMMON /drivers/density,speed,dst_nt,Pdyn_nPa,BxIMF_nt,ByIMF_nt
     &       ,BzIMF_nt,G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
     &       ,W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al
      real*8  density,speed,dst_nt,Pdyn_nPa,BxIMF_nt,ByIMF_nt,BzIMF_nt
      real*8  G1_tsy01,G2_tsy01,fkp,G3_tsy01,W1_tsy04,W2_tsy04
      real*8  W3_tsy04,W4_tsy04,W5_tsy04,W6_tsy04,Al

       integer*4 kext, ifail
       real*8 maginput(25)
c                      1: Kp
c                      2: Dst
c                      3: dens
c                      4: velo
c                      5: Pdyn
c                      6: ByIMF
c                      7: BzIMF
c                      8: G1
c                      9: G2
c                     10: G3
c
c make inputs according to magn. field model chosen
c
c     set fail flag 'on' by default
         ifail = -1
c     clear it if all tests for respective selection are ok
c
c            'none'
           if (kext .eq. 0) then
             ifail = 0
             return
         endif
c            'Olsen-Pfitzer' = default
           if (kext .eq. 5) then
             ifail = 0
             return
         endif

           if (kext .eq. 1) then
c Input for MEAD
           if (maginput(1).eq.baddata) return
             if (maginput(1).le.3.d0) Activ=1
             if (maginput(1).gt.3.d0 .and.
     &         maginput(1).lt.20.d0) Activ=2
             if (maginput(1).ge.20.d0 .and.
     &         maginput(1).lt.30.d0) Activ=3
             if (maginput(1).ge.30.d0) Activ=4
c
             if (maginput(1).lt.0.d0 .or.
     &         maginput(1).gt.90.d0) return
             ifail = 0
             return
         endif
           if (kext .eq. 2) then
c Input for TSYG87s
           if (maginput(1).eq.baddata) return
             if (maginput(1).lt.7.d0) Activ=1
             if (maginput(1).ge.7.d0 .and.
     &         maginput(1).lt.17.d0) Activ=2
             if (maginput(1).ge.17.d0 .and.
     &         maginput(1).lt.20.d0) Activ=3
             if (maginput(1).ge.20.d0 .and.
     &         maginput(1).lt.27.d0) Activ=4
             if (maginput(1).ge.27.d0 .and.
     &         maginput(1).lt.37.d0) Activ=5
             if (maginput(1).ge.37.d0 .and.
     &         maginput(1).lt.47.d0) Activ=6
             if (maginput(1).ge.47.d0) Activ=7
             if (maginput(1).ge.53.d0) Activ=8
c
             if (maginput(1).lt.0.d0 .or.
     &         maginput(1).gt.90.d0) return
             ifail = 0
             return
         endif
           if (kext .eq. 3) then
c Input for TSYG87l
           if (maginput(1).eq.baddata) return
             if (maginput(1).lt.7.d0) Activ=1
             if (maginput(1).ge.7.d0 .and.
     &         maginput(1).lt.17.d0) Activ=2
             if (maginput(1).ge.17.d0 .and.
     &         maginput(1).lt.27.d0) Activ=3
             if (maginput(1).ge.27.d0 .and.
     &         maginput(1).lt.37.d0) Activ=4
             if (maginput(1).ge.37.d0 .and.
     &         maginput(1).lt.47.d0) Activ=5
             if (maginput(1).ge.47.d0) Activ=6
c
             if (maginput(1).lt.0.d0 .or.
     &         maginput(1).gt.90.d0) return
             ifail = 0
             return
           endif
           if (kext .eq. 4) then
c Input for Tsy89
           if (maginput(1).eq.baddata) return
             if (maginput(1).lt.7.d0) Activ=1
             if (maginput(1).ge.7.d0 .and.
     &         maginput(1).lt.17.d0) Activ=2
             if (maginput(1).ge.17.d0 .and.
     &         maginput(1).lt.27.d0) Activ=3
             if (maginput(1).ge.27.d0 .and.
     &         maginput(1).lt.37.d0) Activ=4
             if (maginput(1).ge.37.d0 .and.
     &         maginput(1).lt.47.d0) Activ=5
             if (maginput(1).ge.47.d0 .and.
     &         maginput(1).lt.57.d0) Activ=6
             if (maginput(1).ge.57.d0) Activ=7
c
             if (maginput(1).lt.0.d0 .or.
     &         maginput(1).gt.90.d0) return
             ifail = 0
             return
          endif
           if (kext .eq. 6) then
c Input for OP dyn
           if (maginput(2).eq.baddata) return
           if (maginput(3).eq.baddata) return
           if (maginput(4).eq.baddata) return
               density=maginput(3)
             speed=maginput(4)
             dst_nt=maginput(2)
c
             if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) return
             if (density.lt.5.d0 .or. density.gt.50.d0) return
             if (speed.lt.300.d0 .or. speed.gt.500.d0) return
             ifail = 0
             return
         endif
           if (kext .eq. 7) then
c Input for Tsy96
           if (maginput(2).eq.baddata) return
           if (maginput(5).eq.baddata) return
           if (maginput(6).eq.baddata) return
           if (maginput(7).eq.baddata) return
             dst_nt=maginput(2)
             Pdyn_nPa=maginput(5)
             ByIMF_nt=maginput(6)
             BzIMF_nt=maginput(7)
c
             if (dst_nt.lt.-100.d0 .or. dst_nt.gt.20.d0) return
             if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.10.d0) return
             if (ByIMF_nt.lt.-10.d0 .or. ByIMF_nt.gt.10.d0) return
             if (BzIMF_nt.lt.-10.d0 .or. BzIMF_nt.gt.10.d0) return
             ifail = 0
             return
         endif
           if (kext .eq. 8) then
c Input for Ostapenko97
           if (maginput(2).eq.baddata) return
           if (maginput(5).eq.baddata) return
           if (maginput(7).eq.baddata) return
             dst_nt=maginput(2)
             Pdyn_nPa=maginput(5)
             BzIMF_nt=maginput(7)
             fkp=maginput(1)*1.d0/10.d0
             ifail = 0
             return
         endif
           if (kext .eq. 9) then
c Input for Tsy01
           if (maginput(2).eq.baddata) return
           if (maginput(5).eq.baddata) return
           if (maginput(6).eq.baddata) return
           if (maginput(7).eq.baddata) return
           if (maginput(8).eq.baddata) return
           if (maginput(9).eq.baddata) return
             dst_nt=maginput(2)
             Pdyn_nPa=maginput(5)
             ByIMF_nt=maginput(6)
             BzIMF_nt=maginput(7)
             G1_tsy01=maginput(8)
             G2_tsy01=maginput(9)
c
             if (dst_nt.lt.-50.d0 .or. dst_nt.gt.20.d0) return
             if (Pdyn_nPa.lt.0.5d0 .or. Pdyn_nPa.gt.5.d0) return
             if (ByIMF_nt.lt.-5.d0 .or. ByIMF_nt.gt.5.d0) return
             if (BzIMF_nt.lt.-5.d0 .or. BzIMF_nt.gt.5.d0) return
             if (G1_tsy01.lt.0.d0 .or. G1_tsy01.gt.10.d0) return
             if (G2_tsy01.lt.0.d0 .or. G2_tsy01.gt.10.d0) return
             ifail = 0
             return
         endif
           if (kext .eq. 10) then
c Input for Tsy01 storm
           if (maginput(2).eq.baddata) return
           if (maginput(5).eq.baddata) return
           if (maginput(6).eq.baddata) return
           if (maginput(7).eq.baddata) return
           if (maginput(9).eq.baddata) return
           if (maginput(10).eq.baddata) return
             dst_nt=maginput(2)
             Pdyn_nPa=maginput(5)
             ByIMF_nt=maginput(6)
             BzIMF_nt=maginput(7)
             G2_tsy01=maginput(9)
             G3_tsy01=maginput(10)
             ifail = 0
             return
         endif
c
           if (kext .eq. 11) then
c Input for Tsy04 storm
           if (maginput(2).eq.baddata) return
           if (maginput(5).eq.baddata) return
           if (maginput(6).eq.baddata) return
           if (maginput(7).eq.baddata) return
           if (maginput(11).eq.baddata) return
           if (maginput(12).eq.baddata) return
           if (maginput(13).eq.baddata) return
           if (maginput(14).eq.baddata) return
           if (maginput(15).eq.baddata) return
           if (maginput(16).eq.baddata) return
             dst_nt=maginput(2)
             Pdyn_nPa=maginput(5)
             ByIMF_nt=maginput(6)
             BzIMF_nt=maginput(7)
             W1_tsy04=maginput(11)
             W2_tsy04=maginput(12)
             W3_tsy04=maginput(13)
             W4_tsy04=maginput(14)
             W5_tsy04=maginput(15)
             W6_tsy04=maginput(16)
             ifail = 0
             return
         endif
c
        if (kext .eq. 12) then
c Input for Alexeev 2000
           if (maginput(2).eq.baddata) return
           if (maginput(3).eq.baddata) return
           if (maginput(4).eq.baddata) return
           if (maginput(6).eq.baddata) return
           if (maginput(7).eq.baddata) return
           if (maginput(18).eq.baddata) return
           if (maginput(17).eq.baddata) return
           density=maginput(3)
           speed=maginput(4)
           dst_nt=maginput(2)
           BxIMF_nt=maginput(18)
           ByIMF_nt=maginput(6)
           BzIMF_nt=maginput(7)
           Al=maginput(17)
           ifail = 0
           return
        endif

        if (kext .eq. 13 .or. kext .eq. 14) then
c            if (maginput(5).eq.baddata) return
c            Pdyn_nPa=maginput(5)
c           no need for solar-wind inputs anymore
c           they are read directly from the coefficient files
            ifail=0
            return
        endif

c not sure why this is here, it's not implemented in other places
c        if (kext .eq. 14) then
c Input for Mead-Tsyganenko
c           if (maginput(1).eq.baddata) return
c             fkp=maginput(1)*1.d0/10.d0
c             if (maginput(1).lt.0.d0 .or.
c     &         maginput(1).gt.90.d0) return
c           ifail = 0
c           return
c        endif
       print *, ' invalid kext'

       return
       end

c Wrapper and procedure to access many coordinate transformation form the
c ONERA library
c
c =======================================================================
c GEO2GSM
c
c Routine to transform Cartesian GEO to cartesian GSM coordinates
c
c  INPUTS: iyr = integer year
c          idoy = integer day of year
c          secs = UT in seconds
c          xGEO = 3D array of cartesian position in GEO (Re)
c
c OUTPUTS: psi: angle for GSM coordinate
c         xGSM: 3D array of cartesian position in GSM  (Re)
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'geo2gsm_', $                         ;The entry point
c           iyr,idoy,secs,psi,xGEO,xGSM, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c =======================================================================
c GSM2GEO
c
c Routine to transform Cartesian GSM to cartesian GEO coordinates
c
c  INPUTS: iyr = integer year
c          idoy = integer day of year
c          secs = UT in seconds
c          xGSM = 3D array of cartesian position in GSM (Re)
c
c OUTPUTS: psi: angle for GSM coordinate
c         xGEO: 3D array of cartesian position in GEO (Re)
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'gsm2geo_', $                         ;The entry point
c           iyr,idoy,secs,psi,xGSM,xGEO, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c
c =======================================================================
c GDZ2GEO
c
c Routine to transform GEODEZIC coordinates to cartesian GEO coordinates
c
c  INPUTS: lati = latitude (degres)
c          longi = longitude (degres)
c          alti = altitude (km)
c
c OUTPUTS: xx = xGEO (Re)
c          yy = yGEO (Re)
c          zz = zGEO (Re)
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'gdz2geo_', $                         ;The entry point
c           lati,longi,alti,xx,yy,zz, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c
c =======================================================================
c GEO2GDZ
c
c Routine to transform cartesian GEO coordinates to GEODEZIC coordinates
c
c INPUTS: xx = xGEO (Re)
c          yy = yGEO (Re)
c          zz = zGEO (Re)
c
c OUTPUTS: lati = latitude (degres)
c          longi = longitude (degres)
c          alti = altitude (km)
c
c
c CALLING SEQUENCE from IDL:
c  result = call_external(lib_name,  $            ;The sharable object file
c           'geo2gdz_', $                         ;The entry point
c           xx,yy,zz,lati,longi,alti, $               ;return values (6)
c           /f_value)                             ;function returns a float.
c
c =======================================================================
c --------------------------------------------------------------------
c
        SUBROUTINE geo2gsm1(iyr,idoy,secs,psi,xGEO,xGSM)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xGSM(3)

      dyear=iyr+0.5d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_GSM(xGEO,xGSM)
        end
c --------------------------------------------------------------------
c
        SUBROUTINE gsm2geo1(iyr,idoy,secs,psi,xGSM,xGEO)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xGSM(3)


      dyear=iyr+0.5d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GSM_GEO(xGSM,xGEO)
        end

c --------------------------------------------------------------------
c
        SUBROUTINE geo2gse1(iyr,idoy,secs,xGEO,xGSE)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xGSE(3)

        dyear=iyr+0.5d0
        psi=0.d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_GSE(xGEO,xGSE)
        end
c --------------------------------------------------------------------
c
        SUBROUTINE gse2geo1(iyr,idoy,secs,xGSE,xGEO)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xGSE(3)


      dyear=iyr+0.5d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GSE_GEO(xGSE,xGEO)
        end



c --------------------------------------------------------------------
c
        SUBROUTINE geo2gei1(iyr,idoy,secs,xGEO,xGEI)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xGEI(3)

       dyear=iyr+0.5d0
       psi=0.d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_GEI(xGEO,xGEI)
        end

c --------------------------------------------------------------------
c
        SUBROUTINE gei2geo1(iyr,idoy,secs,xGEI,xGEO)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xGEI(3)

        dyear=iyr+0.5d0
        psi=0.d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEI_GEO(xGEI,xGEO)
        end

c --------------------------------------------------------------------
c
        SUBROUTINE geo2sm1(iyr,idoy,secs,xGEO,xSM)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xSM(3)

       dyear=iyr+0.5d0
       psi=0.d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GEO_SM(xGEO,xSM)
        end

c --------------------------------------------------------------------
c
        SUBROUTINE sm2geo1(iyr,idoy,secs,xSM,xGEO)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGEO(3),xSM(3)

        dyear=iyr+0.5d0
        psi=0.d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL SM_GEO(xSM,xGEO)
        end

c --------------------------------------------------------------------
c
        SUBROUTINE gsm2sm1(iyr,idoy,secs,xGSM,xSM)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGSM(3),xSM(3)

        dyear=iyr+0.5d0
        psi=0.d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL GSM_SM(xGSM,xSM)
        end

c --------------------------------------------------------------------
c
        SUBROUTINE sm2gsm1(iyr,idoy,secs,xSM,xGSM)
      INTEGER*4 iyr,idoy
      REAL*8    secs,psi,dyear
      REAL*8    xGSM(3),xSM(3)

        dyear=iyr+0.5d0
        psi=0.d0
        call initize ! sets rad, pi used by various routines

        CALL INIT_DTD(dyear)
        CALL INIT_GSM(iyr,idoy,secs,psi)
        CALL SM_GSM(xSM,xGSM)
        end
c --------------------------------------------------------------------
c
        SUBROUTINE geo2mag1(iyr,xGEO,xMAG)
      INTEGER*4 iyr
      REAL*8    dyear
      REAL*8    xGEO(3),xMAG(3)

      dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL GEO_MAG(xGEO,xMAG)
        end

c --------------------------------------------------------------------
c
        SUBROUTINE mag2geo1(iyr,xMAG,xGEO)
      INTEGER*4 iyr
      REAL*8    dyear
      REAL*8    xGEO(3),xMAG(3)

      dyear=iyr+0.5d0
        CALL INIT_DTD(dyear)
        CALL MAG_GEO(xMAG,xGEO)
        end
c
c --------------------------------------------------------------------
c
      subroutine msis86(ntime,whichAp,DOY,UT,ALT,LAT,LONG,F107A,
     &F107,AP,Dens,Temp)
c
      IMPLICIT NONE
      INTEGER*4 ISW,ntime,whichAp,DOY(100000),I,J
      REAL*8 SV(25),STL
      REAL*8 UT(100000),ALT(100000),LAT(100000),LONG(100000)
      REAL*8 F107A(100000),F107(100000),AP(7,100000)
      REAL*8 Dens(8,100000),Temp(2,100000),APin(7),D(8),T(2)
c
      COMMON/CSWI/ISW
      DO I=1,25
         SV(I)=1.D0
      ENDDO
      if (WhichAp .EQ.2) SV(9)=-1.D0
      CALL TSELEC(SV)
      ISW=64999
      DO I=1,ntime
         STL=UT(I)/3600.D0+Long(I)/15.D0
       DO J=1,7
          APin(J)=AP(J,I)
       ENDDO
       IF (ALT(I).GE.85.D0) then
          CALL GTS5(DOY(I),UT(I),ALT(I),LAT(I),LONG(I),STL,F107A(I),
     &             F107(I),APin,48,D,T)
         ELSE
            DO J=1,8
             D(J)=-1.D31
          ENDDO
          DO J=1,2
             T(J)=-1.D31
          ENDDO
       ENDIF
         DO J=1,8
          Dens(J,I)=D(J)
       ENDDO
       DO J=1,2
          Temp(J,I)=T(J)
       ENDDO
      ENDDO
      END

c --------------------------------------------------------------------
c
      subroutine msise90(ntime,whichAp,DOY,UT,ALT,LAT,LONG,F107A,
     &F107,AP,Dens,Temp)
c
      IMPLICIT NONE
      INTEGER*4 ISW,ntime,whichAp,DOY(100000),I,J
      REAL*8 SV(25),STL
      REAL*8 UT(100000),ALT(100000),LAT(100000),LONG(100000)
      REAL*8 F107A(100000),F107(100000),AP(7,100000)
      REAL*8 Dens(8,100000),Temp(2,100000),APin(7),D(8),T(2)
c
      COMMON/CSWI/ISW
      DO I=1,25
         SV(I)=1.D0
      ENDDO
      if (WhichAp .EQ.2) SV(9)=-1.D0
      CALL TSELEC5(SV)
      ISW=64999
      DO I=1,ntime
         STL=UT(I)/3600.D0+Long(I)/15.D0
       DO J=1,7
          APin(J)=AP(J,I)
       ENDDO
       CALL GTD6(DOY(I),UT(I),ALT(I),LAT(I),LONG(I),STL,F107A(I),
     &             F107(I),APin,48,D,T)
         DO J=1,8
          Dens(J,I)=D(J)
       ENDDO
       DO J=1,2
          Temp(J,I)=T(J)
       ENDDO
      ENDDO
      END
c --------------------------------------------------------------------
c
      subroutine nrlmsise00(ntime,whichAp,DOY,UT,ALT,LAT,LONG,F107A,
     &F107,AP,Dens,Temp)
c
      IMPLICIT NONE
      INTEGER*4 ISW,ntime,whichAp,DOY(100000),I,J
      REAL*8 SV(25),STL
      REAL*8 UT(100000),ALT(100000),LAT(100000),LONG(100000)
      REAL*8 F107A(100000),F107(100000),AP(7,100000)
      REAL*8 Dens(9,100000),Temp(2,100000),APin(7),D(8),T(2)
c
      COMMON/CSWI/ISW
      DO I=1,25
         SV(I)=1.D0
      ENDDO
      if (WhichAp .EQ.2) SV(9)=-1.D0
      CALL TSELEC7(SV)
      ISW=64999
      DO I=1,ntime
         STL=UT(I)/3600.D0+Long(I)/15.D0
       DO J=1,7
          APin(J)=AP(J,I)
       ENDDO
       CALL GTD7(DOY(I),UT(I),ALT(I),LAT(I),LONG(I),STL,F107A(I),
     &             F107(I),APin,48,D,T)
         DO J=1,9
          Dens(J,I)=D(J)
       ENDDO
       DO J=1,2
          Temp(J,I)=T(J)
       ENDDO
      ENDDO
      END


c --------------------------------------------------------------------
c
c Alternate version without useless 100k time/position array
        SUBROUTINE make_lstar_shell_splitting2(Nipa,kext,options,
     & sysaxes,iyearsat,idoy,UT,xIN1,xIN2,xIN3,
     &  alpha,maginput,Lm,Lstar,BMirror,BLOCAL,BMIN,XJ,MLT)
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
C
c declare inputs
        INTEGER*4    kext,k_ext,k_l,options(5),Nalp,Nipa
      PARAMETER (Nalp=25)
        INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     alpha(Nalp)
      real*8     maginput(25)
c
c
c Declare internal variables
      INTEGER*4    iyear,IPA,kint,ifail
        INTEGER*4    Ilflag,t_resol,r_resol
      REAL*8     mlon,BL,BMIR(25),Bmin_tmp
      REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
      REAL*8     xGEOp(3,25)
      real*8     alti,lati,longi
c
c Declare output variables
        REAL*8     Bmirror(Nalp),BMIN,Blocal
      REAL*8     XJ(Nalp),MLT
        REAL*8     Lm(Nalp),Lstar(Nalp)
C
      COMMON /magmod/k_ext,k_l,kint
        COMMON /flag_L/Ilflag
        DATA  xSUN /1.d0,0.d0,0.d0/
      integer*4 int_field_select, ext_field_select
C
      Ilflag=0
      k_ext=kext
      if (options(3).lt.0 .or. options(3).gt.9) options(3)=0
      t_resol=options(3)+1
      r_resol=options(4)+1
      k_l=options(1)
      kint = int_field_select ( options(5) )
      k_ext = ext_field_select ( kext )
c
        CALL INITIZE

      call init_fields ( kint, iyearsat, idoy, ut, options(2) )

      call get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     6    alti, lati, longi, xGEO )

      call set_magfield_inputs ( k_ext, maginput, ifail )

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then !special script to read files and
            call INIT_TS07D_TLPR
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
        end if

      if ( ifail.lt.0 ) then
          DO IPA=1,Nipa
             Lm(IPA)=baddata
             Lstar(IPA)=baddata
             XJ(IPA)=baddata
             Bmirror(IPA)=baddata
          ENDDO
          BMIN=baddata
          GOTO 99
         endif
c
c      load the GDZ and GSM coordinates in unused end of maginput
       maginput(20) = alti
         maginput(21) = lati
       maginput(22) = longi
         CALL GEO_GSM(xGEO,maginput(23))

c      collect the B components for return in maginput() array
        CALL CHAMP(xGEO,maginput(17),maginput(16),Ifail)
      IF (Ifail.LT.0) THEN
         maginput(17)=baddata
         maginput(18)=baddata
         maginput(19)=baddata
         maginput(16)=baddata
      ENDIF
c
c Compute Bmin assuming 90deg PA at S/C
         k_l=0
           IPA=1
c             all returned values except Bmin are subsequently overwritten
c            this call is required to return a valid Bmin value, as it is
c            not guaranteed that input pitch angles will produce valid results
           CALL calcul_Lstar_opt(t_resol,r_resol,xGeo
     &         ,Lm(IPA),Lstar(IPA),XJ(IPA)
     &         ,Bmirror(IPA),BMIN)
         k_l=options(1)
c            for the specified location and array of pitch angles,
c              return 'Blocal' for the specified location
c              return Bmir() array for bmirror values for each pitch angle
c              return xGeop() array for location of mirror points ""  ""
           CALL find_bm_nalpha(xGEO,nipa,alpha,Blocal,BMIR,xGEOp)
           DO IPA=1,Nipa
           IF (Bmir(ipa).NE.baddata) THEN
             Ilflag=0
c    for each mirror point location in xGEOp, calc several values
               CALL calcul_Lstar_opt(t_resol,r_resol,xGEOp(1,ipa)
     &         ,Lm(IPA),Lstar(IPA),XJ(IPA)
c        bmirror is the blocal at the xGEOp mirror point,
c         (should be identical to the BMIR() array returned by find_bm_nalpha)
     &         ,Bmirror(IPA),BMIN_tmp)
             ELSE
             Lm(IPA)=baddata
             Lstar(IPA)=baddata
             XJ(IPA)=baddata
             Bmirror(IPA)=baddata
           ENDIF
           ENDDO
99         continue
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon1)
           CALL GSM_GEO(xSUN,xGEO)
           CALL geo_mag(xGEO,xMAG)
           CALL car_sph(xMAG,rM,MLAT,Mlon)
           MLT = (Mlon1 - Mlon)/15.d0 + 12.d0
           IF (MLT.GE.24.d0) MLT = MLT - 24.d0
           IF (MLT.LT.0.d0) MLT = MLT + 24.d0

      END

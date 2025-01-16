!***************************************************************************************************
! Copyright 2006, 2007 S. Bourdarie
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
!               SUBROUTINE fly_in_afrl_crres1: Allows to fly a spacecraft in AFRL CRRES ele and pro models.
!               SUBROUTINE get_crres_flux: Compute flux from AFRL CRRES models for an energy, B/B0 and L
!               SUBROUTINE Init_CRRESPRO_Quiet: Set block data for AFRL CRRESPRO QUIET model
!               SUBROUTINE Init_CRRESPRO_Active: Set block data for AFRL CRRESPRO Active model
!               SUBROUTINE Init_CRRESELE_Average: Set block data for AFRL CRRESELE Average model
!               SUBROUTINE Init_CRRESELE_WorstCase: Set block data for AFRL CRRESELE Worstcase model
!
!***************************************************************************************************
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: S. Bourdarie - March 2007 (add multi channel calculcations; add time input for coordinates transformations) - V4.1
!
! DESCRIPTION: Allows to fly a spacecraft in CRRES models (AFRL)
!
! INPUT: ntime -> maximum number of time to compute (long integer)
!        sysaxes -> define which coordinate system is provided in
!        whichm -> which model to use, 1=pro quiet 2=pro active 3=ele ave  4=ele worst case  5=ele Ap15 (long integer)
!        whatf -> which kind of flux, 1=differential 2=E range 3=integral (long integer)
!        Nene -> Number of energy channels to compute
!        energy -> energy (MeV) at which fluxes must be computed (double array [2,25])
!        iyear,idoy,UT -> times when flux are to be computed (not usefull if imput position is not in GSE, GSM, SM,GEI) (respectively long array(ntime_max), long array(ntime_max), double array(ntime_max))
!        xIN1 -> first coordinate in the chosen system (double array [ntime_max])
!        xIN2 -> second coordinate in the chosen system (double array [ntime_max])
!        xIN3 -> third coordinate in the chosen system (double array [ntime_max])
!        Ap15 -> 15 previous days average of Ap index assuming a one day delay (double array [ntime_max])
!
! OUTPUT: flux -> Computed fluxes (MeV-1 cm-2 s-1 or cm-2 s-1) (double array [ntime_max,25])
!
! CALLING SEQUENCE: CALL fly_in_afrl_crres1(ntime,sysaxes,whichm,whatf,energy,xIN1,xIN2,xIN3,flux)
!---------------------------------------------------------------------------------------------------
        SUBROUTINE fly_in_afrl_crres1(ntime,sysaxes,whichm,whatf,nene,
     &                             energy,iyear,idoy,UT,xIN1,xIN2,xIN3,
     &                             Ap15,flux,ascii_path,STRLEN)
c
	IMPLICIT NONE
        INCLUDE 'ntime_max.inc'
C
c declare inputs
        INTEGER*4  STRLEN
        BYTE       ascii_path(strlen)
c
        INTEGER*4  nene_max
	PARAMETER (nene_max=25)
        INTEGER*4  ntime,sysaxes,whichm,whatf,nene
	INTEGER*4  iyear(ntime_max),idoy(ntime_max)
	REAL*8     energy(2,nene_max)
	REAL*8     UT(ntime_max)
	real*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
c Declare internal variables
	INTEGER*4  k_ext,k_l,isat,kint
        INTEGER*4  t_resol,r_resol,Ilflag
	INTEGER*4  i,iyear_dip,idoy_dip
	REAL*8     dec_year,ut_dip
	REAL*8     xGEO(3),xMAG(3),xSUN(3),rM,MLAT,Mlon1
	REAL*8     xGSM(3),xSM(3),xGEI(3),xGSE(3)
	real*8     alti,lati,longi,psi,tilt
        REAL*8     ERA,AQUAD,BQUAD
        REAL*8     BLOCAL(ntime_max),BMIN(ntime_max),XJ(ntime_max)
        REAL*8     Lm(ntime_max),Lstar(ntime_max),BBo(ntime_max)
	REAL*8     Ap15(ntime_max)
c
c Declare output variables	
	REAL*8     flux(ntime_max,nene_max)
c
	CHARACTER*(500) afrl_crres_path
C
        COMMON/GENER/ERA,AQUAD,BQUAD
	COMMON /magmod/k_ext,k_l,kint
        COMMON /flag_L/Ilflag
        COMMON /dip_ang/tilt
        DATA  xSUN /1.d0,0.d0,0.d0/
        REAL*8     pi,rad
        common /rconst/rad,pi
C
        DO i=1,strlen
	   afrl_crres_path(i:i)=char(ascii_path(i))
	ENDDO
	afrl_crres_path=afrl_crres_path(1:strlen)
	Ilflag=0
	k_ext=5
	t_resol=3
	r_resol=0
	k_l=0
	if (whichm .lt. 1 .or. whichm .gt. 5) then
           whichm=1
	   WRITE(6,*)
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)'Invalid AFRL CRRES model specification'
	   WRITE(6,*)'Selecting crrespro quiet'
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)
	endif   
c
	if (whatf .lt. 1 .or. whatf .gt. 3) then
           whatf=1
	   WRITE(6,*)
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)'Invalid flux output specification'
	   WRITE(6,*)'Selecting differential flux'
	   WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	   WRITE(6,*)
	endif   
c
	kint=1    !set IGRF
        CALL INITIZE
	dec_year=1985.5d0
	CALL INIT_DTD(dec_year)
c
	iyear_dip=1985
        idoy_dip=182
	UT_dip=0.d0
        CALL INIT_GSM(iyear_dip,idoy_dip,UT_dip,psi)
        tilt = psi/rad
        DO isat = 1,ntime
	
	    call get_coordinates ( sysaxes, 
     6        xIN1(isat), xIN2(isat), xIN3(isat), 
     6        alti, lati, longi, xGEO )
	    
           CALL calcul_Lstar_opt(t_resol,r_resol,xGEO
     &     ,Lm(isat),Lstar(isat),XJ(isat),BLOCAL(isat),BMIN(isat))
           BBo(isat)=BLOCAL(isat)/BMIN(isat)
	   if (Lm(isat) .le. 0.D0 .and. Lm(isat) .ne. -1.D31) 
     &       Lm(isat)=-Lm(isat)
	enddo   
        call get_crres_flux(ntime,whichm,whatf,nene,
     &               energy,BBo,Lm,Ap15,flux,afrl_crres_path,strlen)
        end
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006
! MODIFICATION: S. Bourdarie - March 2007 (add multi channel calculcations) - V4.1
!
! DESCRIPTION: Compute flux from CRRES models for an energy, B/B0 and L
!
! INPUT: ntmax -> maximum number of time to compute (long integer)
!        whichm -> which model to use, 1=pro quiet 2=pro active 3=ele ave  4=ele worst case  5=ele Ap15 (long integer)
!        whatf -> which kind of flux, 1=differential 2=E range 3=integral (long integer)
!        Nene -> Number of energy channels to compute
!        energy -> energy (MeV) at which fluxes must be computed (double array [2,25])
!        BBo -> Blocal/Bequator (double array [ntime_max])
!        L -> McIlwain L (double array [ntime_max])
!        Ap15 -> 15 previous days average of Ap index assuming a one day delay (double array [ntime_max])
!
! OUTPUT: flux -> Computed fluxes (MeV-1 cm-2 s-1 or cm-2 s-1) (double array [ntime_max,25])
!
! CALLING SEQUENCE: get_crres_flux(ntmax,whichm,whatf,energy,BBo,L,flux)
!---------------------------------------------------------------------------------------------------
       SUBROUTINE get_crres_flux(ntmax,whichm,whatf,nene,
     &            energy,BBo,L,Ap15,flux,afrl_crres_path,strlen)
c       
       IMPLICIT NONE
       INCLUDE 'variables.inc'
       INCLUDE 'ntime_max.inc'
c
       INTEGER*4    STRLEN
c
       INTEGER*4   ntmax,nene,i,ie,il,ib,ie1,ie2,iap,ieny
       INTEGER*4   whichm !1=pro quiet 2=pro active 3=ele ave  4=ele worst case  5=ele Ap15
       INTEGER*4   whatf  !1=diff 2=E range 3=int
c
       INTEGER*4   Ne,Nl,Nbb0,ind
c
       REAL*8      energy(2,25)
       REAL*8      flux(ntime_max,25),BBo(ntime_max),L(ntime_max)
       REAL*8      Ap15(ntime_max)
       REAL*8      pente,cste,Flux1,Flux2
c
c crres variables
       REAL*8 c_Ec(22)
       REAL*8 c_bb0(35)
       REAL*8 c_flux(22,34,90)  !energy;bb0,Lshell
       REAL*8 c_Lshell(91)
       REAL*8 c_Ap15(7)
c
       CHARACTER*(*) afrl_crres_path
C
c
c       common /crres_model/Ne,Nbb0,Nl,c_Lshell,c_Ec,c_bb0,c_flux
       common /crres_model_int/Ne,Nbb0,Nl
       common /crres_model_dbl/c_Lshell,c_Ec,c_bb0,c_flux
c       
       data c_Ap15 /5.0D0,7.0D0,10.0D0,15.0D0,20.0D0,25.0D0,55.0D0/
c
c  init
       DO i=1,ntmax
          do ieny=1,25
             Flux(i,ieny) = baddata
          enddo
       enddo
       CALL Init_CRRES(whichm,afrl_crres_path,strlen)
c	 
       if (whatf.EQ. 1) then
	  DO i=1,ntmax
c
c search in Lbin
	       if (L(i).GE.c_Lshell(1) .AND. 
     &                   L(i).LT.c_Lshell(Nl)) then
	          do il=2,Nl
		     if (L(i).LT.c_Lshell(IL)) goto 10
	          enddo
10                continue
                  Il=Il-1
c
c search in bb0bin
                  if (BBo(i).GT.c_bb0(Nbb0)) goto 50
	          do ib=2,Nbb0
		     if (BBo(i).LE.c_bb0(IB)) goto 20
	          enddo
20                continue
                  ib=ib-1
c
c search in Ap15bin if whichm=5
                  if (whichm .EQ. 5) then 
                     if (Ap15(i).LT.c_Ap15(1) .OR. 
     &                   Ap15(i).GT.c_Ap15(7)) goto 50
	             do iap=2,7
		        if (Ap15(i).LE.c_Ap15(iap)) goto 30
	             enddo
30                   continue
                     iap=iap-1
		     ind=4+iap
		     CALL Init_CRRES(ind,afrl_crres_path,strlen)
		  endif
c
c loop on energies
                  do ieny=1,nene
		     if (energy(1,ieny).LT.c_Ec(1) .OR. 
     &                   energy(1,ieny).GT.c_Ec(Ne)) then
	                Flux(i,ieny) = baddata
                     else
	                do ie=2,Ne
	                   if (energy(1,ieny).LE.c_Ec(ie)) goto 40
	                enddo
40                      continue
c
                        if (c_flux(ie,ib,il).LE.0.D0 .OR.
     &                     c_flux(ie-1,ib,il).LE.0.D0) then
                           Flux(i,ieny) = baddata
		        else
	                   pente=(LOG(c_flux(ie,ib,il))-
     &                           LOG(c_flux(ie-1,ib,il)))/
     &                        (LOG(c_Ec(ie))-LOG(c_Ec(ie-1)))
                           cste=LOG(c_flux(ie,ib,il))-
     &                          pente*LOG(c_Ec(ie))
		           Flux(i,ieny)=exp(pente*LOG(energy(1,ieny))
     &                                  +cste)
		        endif
		     endif
		  enddo
	       else
		  do ieny=1,nene
	             Flux(i,ieny) = baddata
		  enddo
	       endif
50	       continue
	    enddo
	 endif
c	 
c 
         if (whatf.EQ. 2) then
	    DO i=1,ntmax
	       if (L(i).GE.c_Lshell(1) .AND. 
     &                 L(i).LT.c_Lshell(Nl)) then
	          do il=2,Nl
		     if (L(i).LT.c_Lshell(IL)) goto 110
	          enddo
110               continue
                  Il=Il-1
                  if (BBo(i).GT.c_bb0(Nbb0)) then
		     do ieny=1,nene
	                Flux(i,ieny) = baddata
		     enddo
		     goto 150
		  endif
	          do ib=2,Nbb0
		     if (BBo(i).LE.c_bb0(IB)) goto 120
	          enddo
120               continue
                  Ib=Ib-1
c search in Ap15bin if whichm=5
                  if (whichm .EQ. 5) then 
                     if (Ap15(i).LT.c_Ap15(1) .OR. 
     &                   Ap15(i).GT.c_Ap15(7)) then
		        do ieny=1,nene
	                   Flux(i,ieny) = baddata
		        enddo
		        goto 150
		     endif
	             do iap=2,7
		        if (Ap15(i).LE.c_Ap15(iap)) goto 125
	             enddo
125                   continue
                     iap=iap-1
		     ind=4+iap
		     CALL Init_CRRES(ind,afrl_crres_path,strlen)
		  endif
c
c loop on energies
                  do ieny=1,nene
		     if (energy(1,ieny).LT.c_Ec(1) .OR. 
     &                   energy(1,ieny).GT.c_Ec(Ne)) then
	                Flux(i,ieny) = baddata
                     else
	                do ie1=2,Ne
	                   if (energy(1,ieny).LE.c_Ec(ie1)) goto 140
	                enddo
140                     continue
	                do ie2=2,Ne
	                   if (energy(2,ieny).LE.c_Ec(ie2)) goto 145
	                enddo
145                     continue
                        if (c_flux(ie1,ib,il).LE.0.D0 .OR.
     &                      c_flux(ie1-1,ib,il).LE.0.D0) then
                           flux1=0.D0
		        else
	                   pente=(LOG(c_flux(ie1,ib,il))-
     &                         LOG(c_flux(ie1-1,ib,il)))/
     &                        (LOG(c_Ec(ie1))-LOG(c_Ec(ie1-1)))
                           cste=LOG(c_flux(ie1,ib,il))-
     &                          pente*LOG(c_Ec(ie1))
		           Flux1=exp(pente*LOG(energy(1,ieny))+cste)
		        endif
c
                        if (c_flux(ie2,ib,il).LE.0.D0 .OR.
     &                     c_flux(ie2-1,ib,il).LE.0.D0) then
                           flux2=0.D0
		        else
	                   pente=(LOG(c_flux(ie2,ib,il))-
     &                         LOG(c_flux(ie2-1,ib,il)))/
     &                        (LOG(c_Ec(ie2))-LOG(c_Ec(ie2-1)))
                           cste=LOG(c_flux(ie2,ib,il))-
     &                          pente*LOG(c_Ec(ie2))
		           Flux2=exp(pente*LOG(energy(2,ieny))+cste)
		        endif
		        if (ie1.ne.ie2) then
		           Flux(i,ieny)=(Flux1+c_flux(ie1,ib,il))*
     &                               (c_Ec(ie1)-energy(1,ieny))/2.D0
		           do ie=ie1,ie2-2
		              Flux(i,ieny)=Flux(i,ieny)+
     &                           (c_flux(ie,ib,il)+c_flux(ie+1,ib,il))*
     &                           (c_Ec(ie+1)-c_Ec(ie))/2.D0
		           enddo
		           Flux(i,ieny)=Flux(i,ieny)+
     &                        (c_flux(ie2-1,ib,il)+flux2)*
     &                        (energy(2,ieny)-c_Ec(ie2-1))/2.D0
		           Flux(i,ieny)=Flux(i,ieny)/(energy(2,ieny)-
     &                        energy(1,ieny))
                        else
		           Flux(i,ieny)=(Flux1+Flux2)/2.D0
		        endif
		        if (Flux(i,ieny).LE.0.D0) Flux(i,ieny)=baddata
		     endif
		  enddo
	       else
		  do ieny=1,nene
	             Flux(i,ieny) = baddata
		  enddo
	       endif
150	       continue
	    enddo
	 endif
c
c
         if (whatf.EQ. 3) then
	    DO i=1,ntmax
	       if (L(i).GE.c_Lshell(1) .AND. 
     &                    L(i).LT.c_Lshell(Nl)) then
	          do il=2,Nl
		     if (L(i).LT.c_Lshell(IL)) goto 210
	          enddo
210               continue
                  Il=Il-1
                  if (BBo(i).GT.c_bb0(Nbb0)) then
		     do ieny=1,nene
	                Flux(i,ieny) = baddata
		     enddo
		     goto 250
		  endif
	          do ib=2,Nbb0
		     if (BBo(i).LE.c_bb0(IB)) goto 220
	          enddo
220               continue
                  Ib=Ib-1
c search in Ap15bin if whichm=5
                  if (whichm .EQ. 5) then 
                     if (Ap15(i).LT.c_Ap15(1) .OR. 
     &                   Ap15(i).GT.c_Ap15(7)) then
		        do ieny=1,nene
	                   Flux(i,ieny) = baddata
		        enddo
		        goto 250
		     endif
	             do iap=2,7
		        if (Ap15(i).LE.c_Ap15(iap)) goto 225
	             enddo
225                  continue
                     iap=iap-1
		     ind=4+iap
		     CALL Init_CRRES(ind,afrl_crres_path,strlen)
		  endif
c
c loop on energies
                  do ieny=1,nene
		     if (energy(1,ieny).LT.c_Ec(1) .OR. 
     &                   energy(1,ieny).GT.c_Ec(Ne)) then
	                Flux(i,ieny) = baddata
                     else
	                do ie1=2,Ne
	                   if (energy(1,ieny).LE.c_Ec(ie1)) goto 240
	                enddo
240                     continue
c
                        if (c_flux(ie1,ib,il).LE.0.D0 .OR.
     &                     c_flux(ie1-1,ib,il).LE.0.D0) then
                           Flux1 = 0.D0
		        else
	                   pente=(LOG(c_flux(ie1,ib,il))-
     &                         LOG(c_flux(ie1-1,ib,il)))/
     &                        (LOG(c_Ec(ie1))-LOG(c_Ec(ie1-1)))
                           cste=LOG(c_flux(ie1,ib,il))-
     &                        pente*LOG(c_Ec(ie1))
		           Flux1=exp(pente*LOG(energy(1,ieny))+cste)
		        endif
c
		        if (ie1.ne.Ne) then
		           Flux(i,ieny)=(Flux1+c_flux(ie1,ib,il))*
     &                               (c_Ec(ie1)-energy(1,ieny))/2.D0
		           do ie=ie1,Ne-1
		              Flux(i,ieny)=Flux(i,ieny)+
     &                        (c_flux(ie,ib,il)+c_flux(ie+1,ib,il))*
     &                        (c_Ec(ie+1)-c_Ec(ie))/2.D0
		           enddo
                        else
		           Flux(i,ieny)=(Flux1+c_flux(Ne,ib,il))*
     &                               (c_Ec(Ne)-energy(1,ieny))/2.D0
		        endif
		        if (Flux(i,ieny).LE.0.D0) Flux(i,ieny)=baddata
		     endif
		  enddo
	       else
		  do ieny=1,nene
	             Flux(i,ieny) = baddata
		  enddo
	       endif
250	       continue
	    enddo
	 endif
         end 
!
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006 (from Don Brautigam)
!
! DESCRIPTION: Set block data for CRRES model
!              Valid range is: 1.5 <E < 81.3 MeV
!                             1. < bb0 < 684.6
!                             1. < L < 5.5
!
! INPUT: imod, 1=CRRESPRO Quiet, 2=CRRESPRO Active, 3=CRRESELE Average,
!              4=CRRESELE Worst case, 5=CRRESELE 5<Ap15<7, 6=CRRESELE 7<Ap15<10
!              7=CRRESELE 10<Ap15<15, 8=CRRESELE 15<Ap15<20,
!              9=CRRESELE 20<Ap15<25, 10=CRRESELE 25<Ap15<55 (long integer)
! OUTPUT: None
!
! CALLING SEQUENCE: call Init_CRRES(imod)
!---------------------------------------------------------------------------------------------------
       SUBROUTINE Init_CRRES(imod,afrl_crres_path,strlen)
c
       IMPLICIT NONE
       INCLUDE 'variables.inc'
c
       INTEGER*4    STRLEN
c
       INTEGER*4 ish,k,e,imod
       INTEGER*4 Ne,Nbb0,Nl
c
       REAL*8 Ec(22)
       REAL*8 bb0(35)
       REAL*8 flux(22,34,90)
       REAL*8 Lshell(91)
c
       CHARACTER*(*) afrl_crres_path
       CHARACTER*(1000) afrl_crres_path_tmp
C
       common /crres_model_int/Ne,Nbb0,Nl
       common /crres_model_dbl/Lshell,Ec,bb0,flux
c       
       if (imod .EQ. 1) then
          Ne=22
	  Nbb0=35
	  Nl=91
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     &'crrespro_quiet.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crrespro_quiet.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,22
	     do ish=1,90
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
          CALL Init_CRRESPRO
       endif
c
       if (imod .EQ. 2) then
          Ne=22
	  Nbb0=35
	  Nl=91
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crrespro_active.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crrespro_active.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,22
	     do ish=1,90
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
          CALL Init_CRRESPRO
       endif
c
       if (imod .EQ. 3) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Average.txt'
c         open(10,file=
c     &afrl_crres_path(1:STRLEN)//'crresele_Average.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
       if (imod .EQ. 4) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Worst_case.txt'
c          open(10,file=
c     &afrl_crres_path(1:STRLEN)//'crresele_Worst_case.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
       if (imod .EQ. 5) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Ap1.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crresele_Ap1.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
       if (imod .EQ. 6) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Ap2.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crresele_Ap2.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
       if (imod .EQ. 7) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Ap3.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crresele_Ap3.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
       if (imod .EQ. 8) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Ap4.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crresele_Ap4.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
       if (imod .EQ. 9) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Ap5.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crresele_Ap5.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
       if (imod .EQ. 10) then
          Ne=10
	  Nbb0=35
	  Nl=87
	  afrl_crres_path_tmp=afrl_crres_path(1:STRLEN)//
     & 'crresele_Ap6.txt'
c          open(10,file=afrl_crres_path(1:STRLEN)//'crresele_Ap6.txt')
          open(10,file=afrl_crres_path_tmp)
	  do e=1,10
	     do ish=1,86
	        read(10,'(34(D9.3,1x))')(flux(e,k,ish),k=1,34)
	     enddo
	  enddo
	  close(10)
c
c flux /keV in CRRESELE changed to /MeV in the lib.
	  do e=1,10
	     do ish=1,86
	        do k=1,34
	           flux(e,k,ish)=flux(e,k,ish)*1000.D0
		enddo
	     enddo
	  enddo
c
          CALL Init_CRRESELE
       endif
c
	end
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006 (from Don Brautigam)
!
! DESCRIPTION: Set block data for CRRES PRO model
!              Valid range is: 1.5 <E < 81.3 MeV
!                             1. < bb0 < 684.6
!                             1. < L < 5.5
!
! INPUT: None
! OUTPUT: None
!
! CALLING SEQUENCE: call Init_CRRESPRO
!---------------------------------------------------------------------------------------------------
       SUBROUTINE Init_CRRESPRO
c
       IMPLICIT NONE
c
       INTEGER*4 ish
       INTEGER*4 Ne,Nbb0,Nl
c
       REAL*8 Ec(22),tmp_Ec(22)
       REAL*8 bb0(35),tmp_bb0(35)
       REAL*8 flux(22,34,90)
       REAL*8 Lshell(91)
c
       common /crres_model_int/Ne,Nbb0,Nl
       common /crres_model_dbl/Lshell,Ec,bb0,flux
c       common /crres_model/Ne,Nbb0,Nl,Lshell,Ec,bb0,flux
c       
       data tmp_bb0 /1.000D0,1.004D0,1.020D0,1.046D0,1.085D0,1.140D0,
     >           1.200D0,1.300D0,1.400D0,1.520D0,1.690D0,1.880D0,
     >           2.100D0,2.400D0,2.730D0,3.130D0,3.670D0,4.350D0,
     >           5.02D0 ,6.1D0  ,7.410D0,9.088D0,11.29D0,14.22D0,
     >           18.16D0,23.56D0,31.07D0,41.47D0,57.23D0,80.25D0, 
     >           115.4D0,170.7D0,260.7D0,413.4D0,684.6D0/
c
       data tmp_Ec /1.5D0,2.1D0,2.5D0,2.9D0,3.6D0,4.3D0,5.7D0,
     >          6.8D0,8.5D0,
     >          9.7D0,10.7D0,13.2D0,16.9D0,19.4D0,26.3D0,30.9D0,36.3D0,
     >          41.1D0,47.0D0,55.0D0,65.7D0,81.3D0/
c
        do ish=1,35
	   bb0(ish)=tmp_bb0(ish)
	enddo
        do ish=1,22
	   Ec(ish)=tmp_Ec(ish)
	enddo
	do ish=1,91
           Lshell(ish)=1.0D0+(ish-1)*.05D0
	enddo
	end
c
!---------------------------------------------------------------------------------------------------
!                              Introduced in version 3.0
!
! CREATION: S. Bourdarie - May 2006 (from Don Brautigam)
!
! DESCRIPTION: Set block data for CRRES ELE Average model
!              Valid range is: 0.65 <E < 5.75 MeV
!                             1. < bb0 < 684.6
!                             2.5 < L < 6.8
!
! INPUT: None
! OUTPUT: None
!
! CALLING SEQUENCE: call Init_CRRESELE
!---------------------------------------------------------------------------------------------------
       SUBROUTINE Init_CRRESELE
c
       IMPLICIT NONE
c
       INTEGER*4 ish
       INTEGER*4 Ne,Nbb0,Nl
c
       REAL*8 Ec(22),tmp_Ec(22)
       REAL*8 bb0(35),tmp_bb0(35)
       REAL*8 flux(22,34,90)   !energy;bb0,Lshell
       REAL*8 Lshell(91)
c
c       common /crres_model/Ne,Nbb0,Nl,Lshell,Ec,bb0,flux
       common /crres_model_int/Ne,Nbb0,Nl
       common /crres_model_dbl/Lshell,Ec,bb0,flux
c       
       data tmp_bb0 /1.000D0,1.004D0,1.020D0,1.046D0,1.085D0,1.140D0,
     >           1.200D0,1.300D0,1.400D0,1.520D0,1.690D0,1.880D0,
     >           2.100D0,2.400D0,2.730D0,3.130D0,3.670D0,4.350D0,
     >           5.02D0 ,6.1D0  ,7.410D0,9.088D0,11.29D0,14.22D0,
     >           18.16D0,23.56D0,31.07D0,41.47D0,57.23D0,80.25D0, 
     >           115.4D0,170.7D0,260.7D0,413.4D0,684.6D0/
c
       data tmp_Ec /0.65D0,0.95D0,1.60D0,2.00D0,2.35D0,2.75D0,3.15D0,
     >          3.75D0,4.55D0,5.75D0,12*0.D0/
c
        do ish=1,35
	   bb0(ish)=tmp_bb0(ish)
	enddo
        do ish=1,22
	   Ec(ish)=tmp_Ec(ish)
	enddo
	do ish=1,87
           Lshell(ish)=2.5D0+(ish-1)*.05D0
	enddo
	do ish=88,91
           Lshell(ish)=0.D0
	enddo
	end

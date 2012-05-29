!***************************************************************************************************
! Copyright 2007 A. Sicard, 2008 S. Bourdarie
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
c---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.2
c fly_in_ige.f
c
c DESCRIPTION: IGE gives differential or integrated electron spectrum 
c              between 0.9 keV and 5200 keV at geostationary orbit
c              as a function of the solar cycle
c
c INPUT    :   - Year of the launch (t=0 is the minimum of the solar cycle)
c
c              - Duration of the mission (in year)
c 
c OUTPUT     : - Year in the solar cycle according to the year of launch
c
c              - Differential or integrated electron fluxes averaged on all 
c		 the duration of the mission from the year of the launch
c                in MeV-1.cm-2.s-1.sr-1 or cm-2.s-1.sr-1 respectively
c		 in 3 cases (lower, mean, upper)
c
c
c CREATION : Angelica Sicard-Piet - Decembre 2007
c            S. Bourdarie - February 2008
c---------------------------------------------------------------------------------------------------

      SUBROUTINE fly_in_ige1(launch_year,duration,whichm,whatf,nene,
     &                 energy,Lower_flux,Mean_flux,Upper_flux)

      IMPLICIT NONE
      INCLUDE 'variables.inc'

      INTEGER*4 Nenergy
      PARAMETER (Nenergy = 27)
      INTEGER*4 NsE_max
      PARAMETER (NsE_max = 100)
      INTEGER*4 launch_year,year
      INTEGER*4 year_cycle
      INTEGER*4 duration
      INTEGER*4 whichm,whatf
      INTEGER*4 I,J,IE,indE
      INTEGER*4 K,K1
      INTEGER*4 Nene,NsE
      INTEGER*4 indice,flag
      REAL*8 flux_diff_mean(Nenergy)
      REAL*8 energy(2,50),energies(Nenergy)
      REAL*8 Ecin,Emin,Emax
      REAL*8 Esmin,Esmax
      REAL*8 Es(50,NsE_max)
      REAL*8 Mean_flux(50),Upper_flux(50)
      REAL*8 Lower_flux(50)
      REAL*8 Mean_flux1(50),Upper_flux1(50)
      REAL*8 Lower_flux1(50)
      REAL*8 flux_omn,upper_flux_omn
      REAL*8 lower_flux_omn
      REAL*8 flux_omn1,upper_flux_omn1
      REAL*8 lower_flux_omn1
      REAL*8 flux_tab(Nenergy,11)
      REAL*8 p_1,c_1

      COMMON  /Param_IGE/ flux_tab,energies,Emin,Emax

      IF (launch_year .EQ. 0) then
         write(6,*)'*****************************'
         write(6,*)'Launch Year is out of range'
         write(6,*)'*****************************'
         stop
      ENDIF

      CALL calc_year_cycle(launch_year,year_cycle)

      IF (whichm .EQ. 1) 
     &    CALL search_flux_tab_V1

      IF (whichm .EQ. 2) 
     &    CALL search_flux_tab_V2

      IF (whichm .EQ. 3) 
     &    CALL search_flux_tab_V3

      IF (whichm .NE. 1 .AND. whichm .NE. 2 .AND. whichm .NE. 3) then
         write(6,*)'*****************************'
         write(6,*)'Bad Input Parameter: whichm'
         write(6,*)'Available choices: 1, 2 or 3'
         write(6,*)'*****************************'
         stop
      ENDIF 
      IF (whatf .NE. 1 .AND. whatf .NE. 2 .AND. whatf .NE. 3) then
         write(6,*)'*****************************'
         write(6,*)'Bad Input Parameter: whatf'
         write(6,*)'Available choices: 1, 2 or 3'
         write(6,*)'*****************************'
         stop
      ENDIF 

      DO i=1,Nenergy
         flux_diff_mean(i)=0.d0
      ENDDO

      DO I=1,50
         Mean_flux(i)=0.D0
         Upper_flux(i)=0.D0
         Lower_flux(i)=0.D0
         Mean_flux1(i)=0.D0
         Upper_flux1(i)=0.D0
         Lower_flux1(i)=0.D0
      ENDDO

      DO j=1,duration
           year=launch_year+(j-1)
           CALL calc_year_cycle(year,year_cycle)
           IF (year_cycle .LT. -6) year_cycle=-6
           DO i=1,Nenergy
              flux_diff_mean(i)=flux_diff_mean(i)+
     &           flux_tab(i,year_cycle+7)*1000.d0
           ENDDO
      ENDDO

      DO j=1,Nenergy
         flux_diff_mean(j) = flux_diff_mean(j)/duration
      ENDDO

      IF (nene .EQ. 0) THEN
         flag =1
         Nene=Nenergy
         IF (whatf.EQ.1) THEN
            NsE = 1
            DO I=1,Nenergy
               Energy(1,I)= Energies(I)
               Energy(2,I)= 0.d0
               Es(I,1) = Energies(I)
            ENDDO
         ENDIF
         IF (whatf.EQ.2) THEN
            NsE = NsE_max
            DO indE=1,Nenergy-1
               Energy(1,indE)= Energies(indE)
               Energy(2,indE)= Energies(indE+1)
               Esmin=Energies(indE)
               Esmax=Energies(indE+1)
               DO IE = 1,NsE
                  Es(indE,IE) = Esmin*(Esmax/Esmin)
     &                      **((IE-1.D0)/(NsE-1.D0))
              ENDDO
            ENDDO
         ENDIF
         IF (whatf.EQ.3) THEN
             NsE = NsE_max
             Esmax = Emax
             DO indE=1,Nenergy
                Energy(1,indE)=Energies(indE)
                Esmin=Energies(indE)
                DO IE = 1,NsE
                   Es(indE,IE) = Esmin*(Esmax/Esmin)
     &                      **((IE-1.D0)/(NsE-1.D0))
                ENDDO
             ENDDO
         ENDIF
      ELSE
         IF (whatf.EQ.1) THEN
            NsE = 1
            DO I=1,Nene
               Es(I,1) = Energy(1,I)
            ENDDO
         ENDIF
         IF (whatf.EQ.2) THEN
            NsE = NsE_max
            DO indE=1,Nene
               Esmin=Energy(1,indE)
               Esmax=Energy(2,indE)
               DO IE = 1,NsE
                  Es(indE,IE) = Esmin*(Esmax/Esmin)
     &                      **((IE-1.D0)/(NsE-1.D0))
               ENDDO
            ENDDO
         ENDIF
         IF (whatf.EQ.3) THEN
            NsE = NsE_max
            Esmax = Emax
            DO indE=1,Nene
               Esmin=Energy(1,indE)
	       if (Esmin .GT. Emax) then
                  DO IE = 1,NsE
                     Es(indE,IE) = Esmin
                  ENDDO
	       else   
                  DO IE = 1,NsE
                     Es(indE,IE) = Esmin*(Esmax/Esmin)
     &                      **((REAL(IE)-1.D0)/(REAL(NsE)-1.D0))
                  ENDDO
		  if (Es(indE,NsE).GT.Emax) Es(indE,NsE)=Emax
	       endif
            ENDDO
         ENDIF
      ENDIF
      DO indE=1,NeNe
         DO  IE=1,NsE
            Ecin=Es(indE,IE)
            flux_omn = 0.D0
            upper_flux_omn=0.D0
            lower_flux_omn=0.D0
            IF (Ecin .GT. Emax) THEN
               Mean_flux1(indE)=baddata
               Upper_flux1(indE)=baddata
               Lower_flux1(indE)=baddata
               GOTO 203
            ENDIF
            IF (Ecin .LT. Emin) THEN
               Mean_flux1(indE)=baddata
               Upper_flux1(indE)=baddata
               Lower_flux1(indE)=baddata
               GOTO 203
            ENDIF
 
            DO K = 2,Nenergy
               IF (energies(K).GE.Ecin.AND.energies(K).NE.Emin) GOTO 60
            ENDDO
60          CONTINUE
            K1 = K
c
            IF (flux_diff_mean(K1) .lt. 0.d0 .or.
     &         flux_diff_mean(K1-1) .lt. 0.d0 ) goto 203
   
            p_1 = LOG(flux_diff_mean(K1)/
     &         flux_diff_mean(K1-1))
            p_1 = p_1/ LOG(energies(K1)/energies(K1-1))
            c_1 = LOG(flux_diff_mean(K1))
     &         -p_1*LOG(energies(K1)) 
            flux_omn = EXP(p_1*LOG(Ecin)+c_1)
c       
            IF (Ecin .LE. energies(22)) THEN
               upper_flux_omn = flux_omn*(0.00047d0*Ecin*1.d3+1.4d0)
               lower_flux_omn = flux_omn/(0.00047d0*Ecin*1.d3+1.4d0)
            ELSE
              upper_flux_omn = flux_omn*(0.0006d0*Ecin*1.d3+1.4d0)
              lower_flux_omn = flux_omn/(0.0006d0*Ecin*1.d3+1.4d0)
            ENDIF
            IF (whatf.EQ.1)THEN
               Mean_flux1(indE) = flux_omn
               Upper_flux1(indE) = upper_flux_omn
               Lower_flux1(indE) = lower_flux_omn
            ELSE
               IF (whatf.NE.1 .AND. IE.NE.1 ) THEN
                  IF (Mean_flux1(indE) .gt. 0.d0) then
                     Mean_flux1(indE) = Mean_flux1(indE) 
     &                        +0.5D0*(flux_omn1+flux_omn)
     &                        * (Es(indE,IE)-Es(indE,IE-1))
                     Upper_flux1(indE) = Upper_flux1(indE) 
     &                        +0.5D0*(upper_flux_omn1+upper_flux_omn)
     &                        * (Es(indE,IE)-Es(indE,IE-1))
                     Lower_flux1(indE) = Lower_flux1(indE) 
     &                        +0.5D0*(lower_flux_omn1+lower_flux_omn)
     &                        * (Es(indE,IE)-Es(indE,IE-1))
                  ELSE
                     Mean_flux1(indE) = 0.5D0*(flux_omn1+flux_omn)*
     &                 (Es(indE,IE)-Es(indE,IE-1))
                     Upper_flux1(indE)=0.5D0*(upper_flux_omn1+
     &                 upper_flux_omn)*
     &                 (Es(indE,IE)-Es(indE,IE-1))
                     Lower_flux1(indE)=0.5D0*(lower_flux_omn1+
     &                 lower_flux_omn)*
     &                 (Es(indE,IE)-Es(indE,IE-1))
                  ENDIF
               ENDIF
            ENDIF
           flux_omn1 = flux_omn
           upper_flux_omn1 = upper_flux_omn
           lower_flux_omn1 = lower_flux_omn
203        CONTINUE
         ENDDO
      ENDDO

cindice=0
cDO I=1,50
c  IF (Mean_flux1(I) .NE. 0.d0) THEN
c     indice=indice+1
c     Energy(1,indice)=Energy(1,I)
c     Energy(2,indice)=Energy(2,I)
c     Mean_flux(indice)=Mean_flux1(I)
c     Upper_flux(indice)=Upper_flux1(I)
c     Lower_flux(indice)=Lower_flux1(I)
c  ENDIF
cENDDO
cNene=indice
cDO I=indice+1,50
c   Energy(1,I)=0.d0
c   Energy(2,I)=0.d0
cENDDO
c
      DO I=1,Nene
         Mean_flux(i)=Mean_flux1(I)
         Upper_flux(i)=Upper_flux1(I)
         Lower_flux(i)=Lower_flux1(I)
      ENDDO
c
c
      END
c---------------------------------------------------------------------------------------------------
c                              Introduced in version 4.2
c
c CREATION: A. Sicard-Piet - December 2007
c
c DESCRIPTION:  find the year of the solar cycle
c
c INPUT: year
c OUTPUT: year_cycle
c
c CALLING SEQUENCE: call calc_year_cycle
c---------------------------------------------------------------------------------------------------

      SUBROUTINE calc_year_cycle(year,year_cycle)

      IMPLICIT NONE

      INTEGER*4 I
      INTEGER*4 Nmin, mod_year, nsolcyc
      INTEGER*4 year,year_cycle
      PARAMETER (Nmin = 11)
      REAL*8 tab_min(Nmin)
      REAL*8 last_min,min_sol


      tab_min(1) = 1902.d0
      tab_min(2) = 1913.5d0
      tab_min(3) = 1923.583d0
      tab_min(4) = 1933.667d0
      tab_min(5) = 1944.083d0
      tab_min(6) = 1954.25d0
      tab_min(7) = 1964.75d0
      tab_min(8) = 1976.417d0
      tab_min(9) = 1986.667d0
      tab_min(10) = 1996.333d0
      tab_min(11) = 2008.92d0

      last_min=tab_min(Nmin)
      min_sol=last_min
      
      mod_year = year

C     if we're beyond the last solar minimum
C     then find the "equivalent" year in the last
C     solar cycle. We assume the solar cycle is 11
C     years long (rough estimate) and we propagate from
C     the last known solar minimum date.
      IF ( year .GE. int(last_min)) THEN
         year_cycle  = MOD( mod_year - INT(last_min), 11)
         if (year_cycle .GT. 4) year_cycle = year_cycle - 11
      ELSE
C        process a year found in the table above.
         DO i=1,Nmin 
            IF (year .GE. int(tab_min(i))
     &        .AND. year .LE. int(tab_min(i+1)))  GOTO 101
         ENDDO
101      CONTINUE
         year_cycle=year-int(tab_min(i))
         IF (year_cycle .GT. 4) THEN
            year_cycle=-(int(tab_min(i+1))-year)
         ENDIF
C        we have now found the solar cycle for the given year.
      ENDIF

      END

c---------------------------------------------------------------------------------------------------
c                              Introduced in version 4.2
c
c CREATION: A. Sicard-Piet - December 2007
c
c DESCRIPTION:  Data for IGE-2006 (energies, flux)
c
c INPUT: none
c OUTPUT: flux_tab,energies
c
c CALLING SEQUENCE: call search_flux_tab_LANL
c---------------------------------------------------------------------------------------------------

      SUBROUTINE search_flux_tab_V3

      IMPLICIT NONE

      INTEGER*4      I,J
      INTEGER*4      Nenergy
      PARAMETER (Nenergy = 27)
      REAL*8            energies(Nenergy)
      REAL*8            flux_tab(Nenergy,11)
      REAL*8            tab1(11),tab2(11),tab3(11)
      REAL*8            tab4(11),tab5(11),tab6(11)
      REAL*8            tab7(11),tab8(11),tab9(11)
      REAL*8            tab10(11),tab11(11),tab12(11)
      REAL*8            tab13(11),tab14(11),tab15(11)
      REAL*8            tab16(11),tab17(11),tab18(11)
      REAL*8            tab19(11),tab20(11),tab21(11)
      REAL*8            tab22(11),tab23(11),tab24(11)
      REAL*8            tab25(11),tab26(11),tab27(11)
      REAL*8            Emin,Emax

      COMMON  /Param_IGE/ flux_tab,energies,Emin,Emax
      
      do i=1,Nenergy
         do j=1,11
            flux_tab(i,j)=0.d0
         enddo
      enddo
       
      Emin=0.9171d-3
      Emax=5196.15d-3

      energies(1) =0.9171d-3
      energies(2) =1.204d-3
      energies(3) =1.572d-3
      energies(4) =2.051d-3
      energies(5) =2.668d-3
      energies(6) =3.468d-3
      energies(7) =4.530d-3
      energies(8) =5.899d-3
      energies(9) =7.731d-3
      energies(10) =10.15d-3
      energies(11) =13.28d-3
      energies(12) =17.35d-3
c
      energies(13)=30.d-3
      energies(14)=61.24d-3
      energies(15)=88.74d-3
      energies(16)=125.5d-3
      energies(17)=183.71d-3
      energies(18)=266.22d-3
      energies(19)=396.86d-3
      energies(20)=612.37d-3
      energies(21)=908.3d-3
      energies(22)=1284.52d-3
      energies(23)=1989.97d-3
      energies(24)=2437.21d-3
      energies(25)=3074.08d-3
      energies(26)=3968.62d-3
      energies(27)=5196.15d-3

c      indice 0 du tableau:annee -6 ... indice 10 du tableau:annee +4

c      0.92keV
        DATA tab1 /
     & 1.35D+07,1.43D+07,1.25D+07,9.80D+06,8.57D+06,6.93D+06,
     & 7.13D+06,7.87D+06,8.99D+06,9.98D+06,1.39D+07/
     
      DO I=1,11
         flux_tab(1,I)=tab1(I)
      ENDDO
      
c      1.20keV
        DATA tab2 /
     & 1.15D+07,1.21D+07,1.06D+07,8.41D+06,7.29D+06,6.01D+06,
     & 6.22D+06,6.86D+06,7.80D+06,8.58D+06,1.19D+07/
     
      DO I=1,11
         flux_tab(2,I)=tab2(I)
      ENDDO

c      1.57keV
        DATA tab3 /
     & 9.85D+06,1.04D+07,9.13D+06,7.32D+06,6.34D+06,5.28D+06,
     & 5.47D+06,6.02D+06,6.85D+06,7.47D+06,1.02D+07/
     
      DO I=1,11
         flux_tab(3,I)=tab3(I)
      ENDDO

c      2.05keV
        DATA tab4 /
     & 8.39D+06,8.84D+06,7.80D+06,6.32D+06,5.52D+06,4.63D+06,
     & 4.53D+06,4.71D+06,5.70D+06,6.27D+06,8.53D+06/
     
      DO I=1,11
         flux_tab(4,I)=tab4(I)
      ENDDO

c      2.66keV
        DATA tab5 /
     & 7.07D+06,7.42D+06,6.62D+06,5.42D+06,4.78D+06,4.01D+06,
     & 3.91D+06,4.02D+06,4.93D+06,5.36D+06,7.21D+06/
     
      DO I=1,11
         flux_tab(5,I)=tab5(I)
      ENDDO

c      3.47keV
        DATA tab6 /
     & 5.85D+06,6.06D+06,5.52D+06,4.59D+06,4.08D+06,3.40D+06,
     & 3.54D+06,3.85D+06,4.45D+06,4.70D+06,6.15D+06/
     
      DO I=1,11
         flux_tab(6,I)=tab6(I)
      ENDDO

c      4.53keV
        DATA tab7 /
     & 4.67D+06,4.76D+06,4.45D+06,3.77D+06,3.38D+06,2.79D+06,
     & 2.96D+06,3.21D+06,3.69D+06,3.85D+06,4.94D+06/
     
      DO I=1,11
         flux_tab(7,I)=tab7(I)
      ENDDO

c      5.90keV
        DATA tab8 /
     & 3.57D+06,3.58D+06,3.47D+06,2.99D+06,2.70D+06,2.19D+06,
     & 2.35D+06,2.52D+06,2.88D+06,3.00D+06,3.80D+06/
     
      DO I=1,11
         flux_tab(8,I)=tab8(I)
      ENDDO

c      7.73keV
        DATA tab9 /
     & 2.57D+06,2.54D+06,2.56D+06,2.24D+06,2.04D+06,1.61D+06,
     & 1.75D+06,1.84D+06,2.10D+06,2.20D+06,2.76D+06/
     
      DO I=1,11
         flux_tab(9,I)=tab9(I)
      ENDDO

c      10.15keV
        DATA tab10 /
     & 1.73D+06,1.68D+06,1.77D+06,1.57D+06,1.45D+06,1.10D+06,
     & 1.19D+06,1.24D+06,1.41D+06,1.51D+06,1.87D+06/
     
      DO I=1,11
         flux_tab(10,I)=tab10(I)
      ENDDO

c      13.28keV
        DATA tab11 /
     & 1.08D+06,1.05D+06,1.15D+06,1.01D+06,9.56D+05,7.09D+05,
     & 7.50D+05,7.72D+05,8.79D+05,9.68D+05,1.18D+06/
     
      DO I=1,11
         flux_tab(11,I)=tab11(I)
      ENDDO

c      17.35keV
        DATA tab12 /
     & 6.28D+05,6.10D+05,6.99D+05,6.06D+05,5.85D+05,4.28D+05,
     & 4.35D+05,4.39D+05,5.09D+05,6.17D+05,7.05D+05/
     
      DO I=1,11
         flux_tab(12,I)=tab12(I)
      ENDDO

c      30.00keV
        DATA tab13 /
     & 1.92D+05,1.89D+05,2.39D+05,2.13D+05,2.42D+05,2.14D+05,
     & 1.74D+05,1.61D+05,1.92D+05,2.14D+05,1.59D+05/
     
      DO I=1,11
         flux_tab(13,I)=tab13(I)
      ENDDO

c      61.24keV
        DATA tab14 /
     & 6.66D+04,6.42D+04,8.20D+04,7.54D+04,8.92D+04,7.76D+04,
     & 6.39D+04,5.87D+04,6.82D+04,7.30D+04,5.52D+04/
     
      DO I=1,11
         flux_tab(14,I)=tab14(I)
      ENDDO

c      88.74keV
        DATA tab15 /
     & 2.62D+04,2.48D+04,3.19D+04,3.02D+04,3.70D+04,3.18D+04,
     & 2.64D+04,2.41D+04,2.74D+04,2.83D+04,2.18D+04/
     
      DO I=1,11
         flux_tab(15,I)=tab15(I)
      ENDDO

c      125.5keV
        DATA tab16 /
     & 9.18D+03,8.52D+03,1.11D+04,1.08D+04,1.35D+04,1.17D+04,
     & 9.86D+03,9.02D+03,1.00D+04,9.90D+03,7.83D+03/
     
      DO I=1,11
         flux_tab(16,I)=tab16(I)
      ENDDO

c      183.71keV
        DATA tab17 /
     & 3.45D+03,3.15D+03,4.13D+03,4.22D+03,5.37D+03,4.71D+03,
     & 4.02D+03,3.66D+03,3.91D+03,3.65D+03,2.99D+03/
     
      DO I=1,11
         flux_tab(17,I)=tab17(I)
      ENDDO

c      266.22keV
        DATA tab18 /
     & 1.23D+03,1.10D+03,1.48D+03,1.58D+03,2.07D+03,1.79D+03,
     & 1.53D+03,1.36D+03,1.43D+03,1.27D+03,1.07D+03/
     
      DO I=1,11
         flux_tab(18,I)=tab18(I)
      ENDDO

c      396.86keV
        DATA tab19 /
     & 3.83D+02,3.23D+02,4.77D+02,5.41D+02,7.31D+02,6.22D+02,
     & 5.24D+02,4.53D+02,4.64D+02,3.87D+02,3.33D+02/
     
      DO I=1,11
         flux_tab(19,I)=tab19(I)
      ENDDO

c      612.37keV
        DATA tab20 /
     & 7.56D+01,6.33D+01,9.96D+01,1.22D+02,1.71D+02,1.45D+02,
     & 1.17D+02,9.93D+01,9.87D+01,7.70D+01,6.57D+01/
     
      DO I=1,11
         flux_tab(20,I)=tab20(I)
      ENDDO

c      908.3keV
        DATA tab21 /
     & 1.81D+01,1.57D+01,2.80D+01,3.51D+01,5.39D+01,4.28D+01,
     & 3.26D+01,2.73D+01,2.66D+01,1.97D+01,1.60D+01/
     
      DO I=1,11
         flux_tab(21,I)=tab21(I)
      ENDDO

c      1284.52keV
        DATA tab22 /
     & 4.63D+00,4.56D+00,8.68D+00,1.11D+01,1.86D+01,1.37D+01,
     & 9.84D+00,8.07D+00,7.78D+00,5.52D+00,4.25D+00/
     
      DO I=1,11
         flux_tab(22,I)=tab22(I)
      ENDDO

c      1989.97keV
        DATA tab23 /
     & 5.79D-01,6.73D-01,1.39D+00,1.81D+00,3.41D+00,2.28D+00,
     & 1.52D+00,1.22D+00,1.15D+00,7.74D-01,5.56D-01/
     
      DO I=1,11
         flux_tab(23,I)=tab23(I)
      ENDDO

c      2437.21keV
        DATA tab24 /
     & 1.96D-01,2.45D-01,5.26D-01,6.92D-01,1.37D+00,8.79D-01,
     & 5.65D-01,4.48D-01,4.21D-01,2.76D-01,1.92D-01/
     
      DO I=1,11
         flux_tab(24,I)=tab24(I)
      ENDDO

c      3074.08keV
        DATA tab25 /
     & 5.29D-02,7.23D-02,1.62D-01,2.15D-01,4.54D-01,2.76D-01,
     & 1.71D-01,1.34D-01,1.24D-01,7.93D-02,5.30D-02/
     
      DO I=1,11
         flux_tab(25,I)=tab25(I)
      ENDDO

c      3968.62keV
        DATA tab26 /
     & 9.09D-03,1.37D-02,3.21D-02,4.32D-02,9.74D-02,5.60D-02,
     & 3.31D-02,2.56D-02,2.36D-02,1.46D-02,9.36D-03/
     
      DO I=1,11
         flux_tab(26,I)=tab26(I)
      ENDDO

c      5196.15keV
        DATA tab27 /
     & 1.27D-03,2.12D-03,5.23D-03,7.12D-03,1.72D-02,9.33D-03,
     & 5.27D-03,4.01D-03,3.67D-03,2.19D-03,1.35D-03/
     
      DO I=1,11
         flux_tab(27,I)=tab27(I)
      ENDDO


      END

c---------------------------------------------------------------------------------------------------
c                              Introduced in version 4.2
c
c CREATION: A. Sicard-Piet - December 2007
c
c DESCRIPTION:  Data for POLE V2 (energies, flux)
c
c INPUT: none
c OUTPUT: flux_tab,energies
c
c CALLING SEQUENCE: call search_flux_tab_LANL
c---------------------------------------------------------------------------------------------------

      SUBROUTINE search_flux_tab_V2

      IMPLICIT NONE

      INTEGER*4      I,J
      INTEGER*4      Nenergy
      PARAMETER (Nenergy = 27)
      REAL*8            energies(Nenergy)
      REAL*8            flux_tab(Nenergy,11)
      REAL*8            tab13(11),tab14(11),tab15(11)
      REAL*8            tab16(11),tab17(11),tab18(11)
      REAL*8            tab19(11),tab20(11),tab21(11)
      REAL*8            tab22(11),tab23(11),tab24(11)
      REAL*8            tab25(11),tab26(11),tab27(11)
      REAL*8            Emin,Emax

      COMMON  /Param_IGE/ flux_tab,energies,Emin,Emax
      
      do i=1,Nenergy
         do j=1,11
            flux_tab(i,j)=0.d0
         enddo
      enddo
       
      Emin=30.d-3
      Emax=5196.15d-3


      energies(1) =0.9171d-3
      energies(2) =1.204d-3
      energies(3) =1.572d-3
      energies(4) =2.051d-3
      energies(5) =2.668d-3
      energies(6) =3.468d-3
      energies(7) =4.530d-3
      energies(8) =5.899d-3
      energies(9) =7.731d-3
      energies(10) =10.15d-3
      energies(11) =13.28d-3
      energies(12) =17.35d-3
c
      energies(13)=30.d-3
      energies(14)=61.24d-3
      energies(15)=88.74d-3
      energies(16)=125.5d-3
      energies(17)=183.71d-3
      energies(18)=266.22d-3
      energies(19)=396.86d-3
      energies(20)=612.37d-3
      energies(21)=908.3d-3
      energies(22)=1284.52d-3
      energies(23)=1989.97d-3
      energies(24)=2437.21d-3
      energies(25)=3074.08d-3
      energies(26)=3968.62d-3
      energies(27)=5196.15d-3

c      indice 0 du tableau:annee -6 ... indice 10 du tableau:annee +4

c      30.00keV
        DATA tab13 /
     & 1.92D+05,1.89D+05,2.39D+05,2.13D+05,2.42D+05,2.14D+05,
     & 1.74D+05,1.61D+05,1.92D+05,2.14D+05,1.59D+05/
     
      DO I=1,11
         flux_tab(13,I)=tab13(I)
      ENDDO
      
c      61.24keV
        DATA tab14 /
     & 6.66D+04,6.42D+04,8.20D+04,7.54D+04,8.92D+04,7.76D+04,
     & 6.39D+04,5.87D+04,6.82D+04,7.30D+04,5.52D+04/
     
      DO I=1,11
         flux_tab(14,I)=tab14(I)
      ENDDO
      
c      88.74keV
        DATA tab15 /
     & 2.62D+04,2.48D+04,3.19D+04,3.02D+04,3.70D+04,3.18D+04,
     & 2.64D+04,2.41D+04,2.74D+04,2.83D+04,2.18D+04/
     
      DO I=1,11
         flux_tab(15,I)=tab15(I)
      ENDDO
      
c      125.5keV
        DATA tab16 /
     & 9.18D+03,8.52D+03,1.11D+04,1.08D+04,1.35D+04,1.17D+04,
     & 9.86D+03,9.02D+03,1.00D+04,9.90D+03,7.83D+03/
     
      DO I=1,11
         flux_tab(16,I)=tab16(I)
      ENDDO
      
c      183.71keV
        DATA tab17 /
     & 3.45D+03,3.15D+03,4.13D+03,4.22D+03,5.37D+03,4.71D+03,
     & 4.02D+03,3.66D+03,3.91D+03,3.65D+03,2.99D+03/
     
      DO I=1,11
         flux_tab(17,I)=tab17(I)
      ENDDO
      
c      266.22keV
        DATA tab18 /
     & 1.23D+03,1.10D+03,1.48D+03,1.58D+03,2.07D+03,1.79D+03,
     & 1.53D+03,1.36D+03,1.43D+03,1.27D+03,1.07D+03/
     
      DO I=1,11
         flux_tab(18,I)=tab18(I)
      ENDDO
      
c      396.86keV
        DATA tab19 /
     & 3.83D+02,3.23D+02,4.77D+02,5.41D+02,7.31D+02,6.22D+02,
     & 5.24D+02,4.53D+02,4.64D+02,3.87D+02,3.33D+02/
     
      DO I=1,11
         flux_tab(19,I)=tab19(I)
      ENDDO
      
c      612.37keV
        DATA tab20 /
     & 7.56D+01,6.33D+01,9.96D+01,1.22D+02,1.71D+02,1.45D+02,
     & 1.17D+02,9.93D+01,9.87D+01,7.70D+01,6.57D+01/
     
      DO I=1,11
         flux_tab(20,I)=tab20(I)
      ENDDO
      
c      908.3keV
        DATA tab21 /
     & 1.81D+01,1.57D+01,2.80D+01,3.51D+01,5.39D+01,4.28D+01,
     & 3.26D+01,2.73D+01,2.66D+01,1.97D+01,1.60D+01/
     
      DO I=1,11
         flux_tab(21,I)=tab21(I)
      ENDDO
      
c      1284.52keV
        DATA tab22 /
     & 4.63D+00,4.56D+00,8.68D+00,1.11D+01,1.86D+01,1.37D+01,
     & 9.84D+00,8.07D+00,7.78D+00,5.52D+00,4.25D+00/
     
      DO I=1,11
         flux_tab(22,I)=tab22(I)
      ENDDO
      
c      1989.97keV
        DATA tab23 /
     & 5.79D-01,6.73D-01,1.39D+00,1.81D+00,3.41D+00,2.28D+00,
     & 1.52D+00,1.22D+00,1.15D+00,7.74D-01,5.56D-01/
     
      DO I=1,11
         flux_tab(23,I)=tab23(I)
      ENDDO
      
c      2437.21keV
        DATA tab24 /
     & 1.96D-01,2.45D-01,5.26D-01,6.92D-01,1.37D+00,8.79D-01,
     & 5.65D-01,4.48D-01,4.21D-01,2.76D-01,1.92D-01/
     
      DO I=1,11
         flux_tab(24,I)=tab24(I)
      ENDDO
      
c      3074.08keV
        DATA tab25 /
     & 5.29D-02,7.23D-02,1.62D-01,2.15D-01,4.54D-01,2.76D-01,
     & 1.71D-01,1.34D-01,1.24D-01,7.93D-02,5.30D-02/
     
      DO I=1,11
         flux_tab(25,I)=tab25(I)
      ENDDO
      
c      3968.62keV
        DATA tab26 /
     & 9.09D-03,1.37D-02,3.21D-02,4.32D-02,9.74D-02,5.60D-02,
     & 3.31D-02,2.56D-02,2.36D-02,1.46D-02,9.36D-03/
     
      DO I=1,11
         flux_tab(26,I)=tab26(I)
      ENDDO
      
c      5196.15keV
        DATA tab27 /
     & 1.27D-03,2.12D-03,5.23D-03,7.12D-03,1.72D-02,9.33D-03,
     & 5.27D-03,4.01D-03,3.67D-03,2.19D-03,1.35D-03/
     
      DO I=1,11
         flux_tab(27,I)=tab27(I)
      ENDDO
      
     

      END
c---------------------------------------------------------------------------------------------------
c                              Introduced in version 4.2
c
c CREATION: A. Sicard-Piet - December 2007
c
c DESCRIPTION:  Data for POLE V1 (energies, flux)
c
c INPUT: none
c OUTPUT: flux_tab,energies
c
c CALLING SEQUENCE: call search_flux_tab_LANL
c---------------------------------------------------------------------------------------------------

      SUBROUTINE search_flux_tab_V1

      IMPLICIT NONE

      INTEGER*4      I,J
      INTEGER*4      Nenergy
      PARAMETER (Nenergy = 27)
      REAL*8            energies(Nenergy)
      REAL*8            flux_tab(Nenergy,11)
      REAL*8            tab13(11),tab14(11),tab15(11)
      REAL*8            tab16(11),tab17(11),tab18(11)
      REAL*8            tab19(11),tab20(11),tab21(11)
      REAL*8            tab22(11)
      REAL*8            Emin,Emax

      COMMON  /Param_IGE/ flux_tab,energies,Emin,Emax

      do i=1,Nenergy
         do j=1,11
            flux_tab(i,j)=0.d0
         enddo
      enddo
       
      Emin=36.74d-3
      Emax=1284.52d-3


      energies(1) =0.9171d-3
      energies(2) =1.204d-3
      energies(3) =1.572d-3
      energies(4) =2.051d-3
      energies(5) =2.668d-3
      energies(6) =3.468d-3
      energies(7) =4.530d-3
      energies(8) =5.899d-3
      energies(9) =7.731d-3
      energies(10) =10.15d-3
      energies(11) =13.28d-3
      energies(12) =17.35d-3
c
      energies(13)=36.74d-3
      energies(14)=61.24d-3
      energies(15)=88.74d-3
      energies(16)=125.5d-3
      energies(17)=183.71d-3
      energies(18)=266.22d-3
      energies(19)=396.86d-3
      energies(20)=612.37d-3
      energies(21)=908.3d-3
      energies(22)=1284.52d-3
      energies(23)=1989.97d-3
      energies(24)=2437.21d-3
      energies(25)=3074.08d-3
      energies(26)=3968.62d-3
      energies(27)=5196.15d-3

c      indice 0 du tableau:annee -6 ... indice 10 du tableau:annee +4

c      36.74keV
        DATA tab13 /
     & 1.53D+05,1.49D+05,1.59D+05,1.67D+05,1.95D+05,1.72D+05,
     & 1.40D+05,1.29D+05,1.54D+05,1.69D+05,1.27D+05/
     
      DO I=1,11
         flux_tab(13,I)=tab13(I)
      ENDDO
     
c      61.24keV
        DATA tab14 /
     & 6.66D+04,6.49D+04,6.96D+04,7.52D+04,8.92D+04,7.76D+04,
     & 6.39D+04,5.87D+04,6.82D+04,7.30D+04,5.52D+04/
     
      DO I=1,11
         flux_tab(14,I)=tab14(I)
      ENDDO
     
c      88.74keV
        DATA tab15 /
     & 2.62D+04,2.55D+04,2.76D+04,3.06D+04,3.70D+04,3.18D+04,
     & 2.64D+04,2.41D+04,2.74D+04,2.83D+04,2.18D+04/
     
      DO I=1,11
         flux_tab(15,I)=tab15(I)
      ENDDO
     
c      125.5keV
        DATA tab16 /
     & 9.18D+03,8.76D+03,9.68D+04,1.10D+04,1.35D+04,1.17D+04,
     & 9.86D+03,9.02D+03,1.00D+04,9.90D+03,7.83D+03/
     
      DO I=1,11
         flux_tab(16,I)=tab16(I)
      ENDDO
     
c      183.71keV
        DATA tab17 /
     & 3.45D+03,3.19D+03,3.63D+03,4.27D+03,5.37D+03,4.71D+03,
     & 4.02D+03,3.66D+03,3.91D+03,3.65D+03,2.99D+03/
     
      DO I=1,11
         flux_tab(17,I)=tab17(I)
      ENDDO
     
c      266.22keV
        DATA tab18 /
     & 1.23D+03,1.09D+03,1.28D+03,1.56D+03,2.07D+03,1.79D+03,
     & 1.53D+03,1.36D+03,1.43D+03,1.27D+03,1.07D+03/
     
      DO I=1,11
         flux_tab(18,I)=tab18(I)
      ENDDO
     
c      396.86keV
        DATA tab19 /
     & 3.83D+02,3.01D+02,4.01D+02,5.26D+02,7.31D+02,6.22D+02,
     & 5.24D+02,4.53D+02,4.64D+02,3.87D+02,3.33D+02/
     
      DO I=1,11
         flux_tab(19,I)=tab19(I)
      ENDDO
     
c      612.37keV
        DATA tab20 /
     & 7.56D+01,5.79D+01,8.14D+01,1.18D+02,1.71D+02,1.45D+02,
     & 1.17D+02,9.93D+01,9.87D+01,7.70D+01,6.57D+01/
     
      DO I=1,11
         flux_tab(20,I)=tab20(I)
      ENDDO
     
c      908.3keV
        DATA tab21 /
     & 1.81D+01,1.46D+01,2.25D+01,3.44D+01,5.39D+01,4.28D+01,
     & 3.26D+01,2.73D+01,2.66D+01,1.97D+01,1.60D+01/
     
      DO I=1,11
         flux_tab(21,I)=tab21(I)
      ENDDO
     
c      1284.52keV
        DATA tab22 /
     & 4.63D+00,4.57D+00,6.83D+00,1.09D+01,1.86D+01,1.37D+01,
     & 9.84D+00,8.07D+00,7.78D+00,5.52D+00,4.25D+00/
     
      DO I=1,11
         flux_tab(22,I)=tab22(I)
      ENDDO

      END

!***************************************************************************************************
! Copyright 2008 S. Bourdarie
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
!***************************************************************************************************
c---------------------------------------------------------------------------------------------------
!                              Introduced in version 4.3
c fly_in_meo_gnss.f
c
c DESCRIPTION: MEO-V2 gives differential or integrated electron spectrum 
c              between 280 keV and 2240 keV at GPS navigation orbit
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
c CREATION : Sebastien Bourdarie - March 2008
c---------------------------------------------------------------------------------------------------

      SUBROUTINE fly_in_meo_gnss1(launch_year,duration,whichm,
     &                 whatf,nene,
     &                 energy,Lower_flux,Mean_flux,Upper_flux)

      IMPLICIT NONE
      INCLUDE 'variables.inc'

      INTEGER*4 Nenergy
      PARAMETER (Nenergy = 30)
      INTEGER*4 launch_year,year
      INTEGER*4 year_cycle
      INTEGER*4 duration
      INTEGER*4 whichm,whatf
      INTEGER*4 I,J,IE,indE
      INTEGER*4 K,K1
      INTEGER*4 Nene,NsE
      INTEGER*4 indice,flag
      REAL*8 flux_int_mean(Nenergy)
      REAL*8 energy(2,50),energies(Nenergy)
      REAL*8 Ecin,Emin,Emax
      REAL*8 Esmin,Esmax
      REAL*8 Es(100,20)
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
      REAL*8 eps
      PARAMETER (eps = 0.01)
      
      COMMON  /Param_meo_gnss/ flux_tab,energies,Emin,Emax

      IF (launch_year .EQ. 0) then
         write(6,*)'*****************************'
         write(6,*)'Launch Year is out of range'
         write(6,*)'*****************************'
         stop
      ENDIF

      CALL calc_year_cycle(launch_year,year_cycle)

      IF (whichm .EQ. 1) 
     &    CALL search_flux_tab_meoV1

      IF (whichm .EQ. 2) 
     &    CALL search_flux_tab_meov2

      IF (whichm .NE. 1 .AND. whichm .NE. 2) then
         write(6,*)'*****************************'
         write(6,*)'Bad Input Parameter: whichm  '
         write(6,*)'Available choices: 1, or 2   '
         write(6,*)'*****************************'
         stop
      ENDIF 
      IF (whatf .NE. 1 .AND. whatf .NE. 2 .AND. whatf .NE. 3) then
         write(6,*)'*****************************'
         write(6,*)'Bad Input Parameter: whatf   '
         write(6,*)'Available choices: 1, 2 or 3 '
         write(6,*)'*****************************'
         stop
      ENDIF 

      DO i=1,Nenergy
         flux_int_mean(i)=0.d0
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
              flux_int_mean(i)=flux_int_mean(i)+
     &           flux_tab(i,year_cycle+7)
           ENDDO
      ENDDO

      DO j=1,Nenergy
         flux_int_mean(j) = flux_int_mean(j)/duration
      ENDDO

      IF (nene .EQ. 0) THEN
         Nene=7
         IF (whatf.EQ.1) THEN
            Nse=2
            DO I=1,Nenergy
               Energy(1,I)= Energies(I)
               Es(I,1) = Energies(I)-Eps
               Es(I,2) = Energies(I)+Eps
               if (Energies(I) .EQ. Emin) Es(I,1) = Energies(I)
               if (Energies(I) .EQ. Emax) Es(I,2) = Energies(I)
            ENDDO
         ENDIF
         IF (whatf.EQ.2) THEN
            Nse=2
            DO indE=1,Nenergy-1
               Energy(1,indE)= Energies(indE)
               Energy(2,indE)= Energies(indE+1)
               Es(indE,1) = Energies(indE)
               Es(indE,2) = Energies(indE+1)
            ENDDO
         ENDIF
         IF (whatf.EQ.3) THEN
             Nse=1
             DO indE=1,Nenergy
                Energy(1,indE)=Energies(indE)
                Es(indE,1) = Energies(indE)
             ENDDO
         ENDIF
      ELSE
         IF (whatf.EQ.1) THEN
            NsE = 2
            DO I=1,Nene
               Es(I,1) = Energy(1,I)-Eps
               Es(I,2) = Energy(1,I)+Eps
               if (Energy(1,I) .EQ. Emin) Es(I,1) = Emin
               if (Energy(1,I) .EQ. Emax) Es(I,2) = Emax
            ENDDO
         ENDIF
         IF (whatf.EQ.2) THEN
            NsE = 2
            DO indE=1,Nene
               Es(indE,1) = Energy(1,indE)
               Es(indE,2) = Energy(2,indE)
            ENDDO
         ENDIF
         IF (whatf.EQ.3) THEN
            NsE = 1
            DO indE=1,Nene
               Es(indE,1) = Energy(1,indE)
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
            IF (flux_int_mean(K1) .lt. 0.d0 .or.
     &         flux_int_mean(K1-1) .lt. 0.d0 ) goto 203
   
            p_1 = LOG(flux_int_mean(K1)/
     &         flux_int_mean(K1-1))
            p_1 = p_1/ LOG(energies(K1)/energies(K1-1))
            c_1 = LOG(flux_int_mean(K1))
     &         -p_1*LOG(energies(K1)) 
            flux_omn = EXP(p_1*LOG(Ecin)+c_1)
c       
            if (whichm .EQ. 1) then
               upper_flux_omn = flux_omn*(0.002d0*Ecin*1.d3+1.76d0)
               lower_flux_omn = flux_omn/(0.002d0*Ecin*1.d3+1.76d0)
            else
               upper_flux_omn = flux_omn*(0.0016d0*Ecin*1.d3+2.4d0)
               lower_flux_omn = flux_omn/(0.0016d0*Ecin*1.d3+2.4d0)
            endif
c            
            IF (whatf .EQ. 3) THEN
               Mean_flux1(indE) = flux_omn
               Upper_flux1(indE) = upper_flux_omn
               Lower_flux1(indE) = lower_flux_omn
            ELSE
               IF (IE .EQ. 1) THEN
                  Mean_flux1(indE) = flux_omn
                  Upper_flux1(indE) = upper_flux_omn
                  Lower_flux1(indE) = lower_flux_omn
               ELSE
                  if (Mean_flux1(indE) .GT. 0.D0 .AND. 
     &               flux_omn .GT. 0.D0) THEN
                     Mean_flux1(indE)=(flux_omn1-flux_omn)/
     &                  (Es(indE,2)-Es(indE,1))
                     Upper_flux1(indE)=(upper_flux_omn1-upper_flux_omn)/
     &                 (Es(indE,2)-Es(indE,1))
                     Lower_flux1(indE)=(lower_flux_omn1-lower_flux_omn)/
     &                 (Es(indE,2)-Es(indE,1))
                  ENDIF
               ENDIF
            ENDIF
           flux_omn1 = flux_omn
           upper_flux_omn1 = upper_flux_omn
           lower_flux_omn1 = lower_flux_omn
203        CONTINUE
         ENDDO
      ENDDO
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
c                              Introduced in version 4.3
c
c CREATION: S. Bourdarie March 2008
c
c DESCRIPTION:  Data for MEO-V2 (energies, flux)
c
c INPUT: none
c OUTPUT: flux_tab,energies
c
c CALLING SEQUENCE: call search_flux_tab_meoV2
c---------------------------------------------------------------------------------------------------

      SUBROUTINE search_flux_tab_meov2

      IMPLICIT NONE

      INTEGER*4      I,J
      INTEGER*4      Nenergy
      PARAMETER (Nenergy = 30)
      REAL*8            energies(Nenergy)
      REAL*8            flux_tab(Nenergy,11)
      REAL*8            tab1(11),tab2(11),tab3(11)
      REAL*8            tab4(11),tab5(11),tab6(11)
      REAL*8            tab7(11)
      REAL*8            Emin,Emax

      COMMON  /Param_meo_gnss/ flux_tab,energies,Emin,Emax
      
      do i=1,Nenergy
         do j=1,11
            flux_tab(i,j)=0.d0
         enddo
      enddo
       
      Emin=280.d-3
      Emax=2240.d-3

      energies(1) =280.d-3
      energies(2) =400.d-3
      energies(3) =560.d-3
      energies(4) =800.d-3
      energies(5) =1120.d-3
      energies(6) =1600.d-3
      energies(7) =2240.d-3
      
c      indice 1 du tableau:annee -6 ... indice 11 du tableau:annee +4

c      280keV
      DATA tab1 /
     &  6.49D+05,5.56D+05,8.28D+05,9.51D+05,1.31D+06,1.10D+06,
     &  9.14D+05,7.89D+05,8.04D+05,6.64D+05,5.65D+05/
     
      DO I=1,11
         flux_tab(1,I)=tab1(I)
      ENDDO
      
c      400keV
      DATA tab2 /
     &  3.79D+05,3.23D+05,5.14D+05,6.18D+05,8.87D+05,7.33D+05,
     &  5.88D+05,5.00D+05,4.99D+05,3.92D+05,3.31D+05/
     
      DO I=1,11
         flux_tab(2,I)=tab2(I)
      ENDDO

c      560keV
      DATA tab3 /
     &  1.64D+05,1.43D+05,2.44D+05,3.03D+05,4.57D+05,3.67D+05,
     &  2.83D+05,2.38D+05,2.33D+05,1.76D+05,1.45D+05/
     
      DO I=1,11
         flux_tab(3,I)=tab3(I)
      ENDDO

c      800keV
      DATA tab4 /
     &  6.00D+04,5.53D+04,1.02D+05,1.29D+05,2.09D+05,1.59D+05,
     &  1.17D+05,9.70D+04,9.40D+04,6.81D+04,5.39D+04/
     
      DO I=1,11
         flux_tab(4,I)=tab4(I)
      ENDDO

c      1120keV
      DATA tab5 /
     &  2.01D+04,2.05D+04,3.98D+04,5.10D+04,8.84D+04,6.34D+04,
     &  4.47D+04,3.64D+04,3.50D+04,2.45D+04,1.86D+04/
     
      DO I=1,11
         flux_tab(5,I)=tab5(I)
      ENDDO

c      1600keV
      DATA tab6 /
     &  4.94D+03,5.65D+03,1.16D+04,1.51D+04,2.82D+04,1.90D+04,
     &  1.27D+04,1.02D+04,9.67D+03,6.53D+03,4.72D+03/
     
      DO I=1,11
         flux_tab(6,I)=tab6(I)
      ENDDO

c      2240keV
      DATA tab7 /
     &  1.12D+03,1.40D+03,3.02D+03,3.97D+03,7.91D+03,5.05D+03,
     &  3.24D+03,2.57D+03,2.41D+03,1.58D+03,1.09D+03/
     
      DO I=1,11
         flux_tab(7,I)=tab7(I)
      ENDDO

c     
      END
c---------------------------------------------------------------------------------------------------
c                              Introduced in version 4.3
c
c CREATION: S. Bourdarie March 2008
c
c DESCRIPTION:  Data for MEO-V1 (energies, flux)
c
c INPUT: none
c OUTPUT: flux_tab,energies
c
c CALLING SEQUENCE: call search_flux_tab_meoV1
c---------------------------------------------------------------------------------------------------

      SUBROUTINE search_flux_tab_meov1

      IMPLICIT NONE

      INTEGER*4      I,J
      INTEGER*4      Nenergy
      PARAMETER (Nenergy = 30)
      REAL*8            energies(Nenergy)
      REAL*8            flux_tab(Nenergy,11)
      REAL*8            tab1(11),tab2(11),tab3(11)
      REAL*8            tab4(11),tab5(11),tab6(11)
      REAL*8            tab7(11)
      REAL*8            Emin,Emax

      COMMON  /Param_meo_gnss/ flux_tab,energies,Emin,Emax
      
      do i=1,Nenergy
         do j=1,11
            flux_tab(i,j)=0.d0
         enddo
      enddo
       
      Emin=280.d-3
      Emax=1120.d-3

      energies(1) =280.d-3
      energies(2) =400.d-3
      energies(3) =560.d-3
      energies(4) =800.d-3
      energies(5) =1120.d-3
      
      
c      indice 1 du tableau:annee -6 ... indice 11 du tableau:annee +4

c      280keV
      DATA tab1 /
     &  9.17D+05,9.17D+05,9.17D+05,9.17D+05,9.17D+05,9.17D+05,
     &  9.17D+05,9.17D+05,9.17D+05,9.17D+05,9.17D+05/
     
      DO I=1,11
         flux_tab(1,I)=tab1(I)
      ENDDO
      
c      400keV
      DATA tab2 /
     &  5.54D+05,5.54D+05,5.54D+05,5.54D+05,5.54D+05,5.54D+05,
     &  5.54D+05,5.54D+05,5.54D+05,5.54D+05,5.54D+05/
     
      DO I=1,11
         flux_tab(2,I)=tab2(I)
      ENDDO

c      560keV
      DATA tab3 /
     &  2.79D+05,2.79D+05,2.79D+05,2.79D+05,2.79D+05,2.79D+05,
     &  2.79D+05,2.79D+05,2.79D+05,2.79D+05,2.79D+05/
     
      DO I=1,11
         flux_tab(3,I)=tab3(I)
      ENDDO

c      800keV
      DATA tab4 /
     &  1.17D+05,1.17D+05,1.17D+05,1.17D+05,1.17D+05,1.17D+05,
     &  1.17D+05,1.17D+05,1.17D+05,1.17D+05,1.17D+05/
     
      DO I=1,11
         flux_tab(4,I)=tab4(I)
      ENDDO

c      1120keV
      DATA tab5 /
     &  5.15D+04,5.15D+04,5.15D+04,5.15D+04,5.15D+04,5.15D+04,
     &  5.15D+04,5.15D+04,5.15D+04,5.15D+04,5.15D+04/
     
      DO I=1,11
         flux_tab(5,I)=tab5(I)
      ENDDO
c
c     
      END
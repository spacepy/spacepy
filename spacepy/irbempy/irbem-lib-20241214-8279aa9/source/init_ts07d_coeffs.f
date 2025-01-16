      ! Adapted for IRBEM by A. C. Kellerman, needs to be called on date
      ! update. TODO: check if we already loaded the current files and
      ! skip the load if so, as they currently have a cadence of 5
      ! minutes
      ! ifail returns as -1 if the files are not found, resulting in 
      ! bad data returned by each subroutine
      
      SUBROUTINE INIT_TS07D_COEFFS (iyear,idoy,dut,ifail)

      IMPLICIT NONE
     
      INTEGER*4 iyear,idoy,ifail
      REAL*8 A07,dut,COEFF_Q,COEFF_B_RMS,VXGSE,VYGSE,VZGSE
      REAL*8 PDYN,TILT
      INTEGER*4 M_INX,N_INX,NTOT,I,MODMIN
      
      INTEGER IYR,IDY,IHR,IMN,ISC
      INTEGER TS7LEN,ierr,stat,statb(13) ! file i/o
      LOGICAL OK ! I/O

      CHARACTER*18 PAR_FNAME
      CHARACTER*200 filename,FMT
      CHARACTER*80 TS7DIR
      CHARACTER*255 ts07d_env
      CHARACTER TRIM

      PARAMETER (NTOT=101)

      COMMON /TS07D_DATA/ M_INX,N_INX,PDYN,TILT,A07(NTOT)

      CALL GETENV('TS07_DATA_PATH', ts07d_env)
      if (LEN_TRIM(ts07d_env).ne.0) then 
        TS7DIR=TRIM(ts07d_env)
      else
        write(*,*) "error, TS07_DATA_PATH global variable not set"
        stop
      endif
      i=len(TS7DIR)
      do while (TS7DIR(i:i) == ' ')
        i = i - 1  
      enddo
      TS7LEN=i

      IYR=iyear
      IDY=idoy
      IHR=floor(dut/3600.D0)
      IMN=floor(((dut/3600.D0)-IHR) * 60.D0) 
      ISC=floor(dut-IHR*3600.D0-IMN*60.D0)

      MODMIN=MOD(IMN,5)
      WRITE (PAR_FNAME,'(I4,A1,I0.3,A1,I0.2,A1,I0.2,A4)')
     +iyear,'_',idoy,'_',IHR,'_',IMN-MODMIN,'.par' 

C     TS7LEN is the directory length
c      WRITE(FMT,'("(A", I0, ",A8,A18)")') TS7LEN
c      write(filename,FMT), TS7DIR(1:TS7LEN),'/Coeffs/',PAR_FNAME

      WRITE(FMT,'("(A", I0, ",A8,I4,A1,I0.3,A1,A18)")') TS7LEN
      write(filename,FMT) TS7DIR(1:TS7LEN),'/Coeffs/',iyear,'_',
     * idoy,'/',PAR_FNAME

c check that filename exists:
      INQUIRE( FILE=filename, EXIST=OK ) 
      if (OK) then
      OPEN (UNIT=1,FILE=filename,action='read') !  MODEL PARAMETER FILE FOR
        ierr=stat(filename,statb)
      if (ierr.ne.0) then
        print *, 'Could not open ',filename,' for reading.'
        stop
      endif

      READ (1,100) (A07(I),I=1,NTOT)                             !  A SPECIFIC TIME MOMENT
      READ (1,101) COEFF_Q 
      READ (1,101) COEFF_B_RMS 
      READ (1,102) M_INX
      READ (1,102) N_INX
      READ (1,101) PDYN
      READ (1,101) TILT
      CLOSE(1)                                        

c      print *, filename
c      do i=1,NTOT
c        write(*,*) 'i=',i,', A07=',A07(i)
c        call sleep(1)
c      enddo
c      write(*,*) COEFF_G, COEFF_B_RMS, M_INX, N_INX, PDYN, TILT
c      call sleep(10)

 100  FORMAT(G15.6)                                         
 101  FORMAT(7x,G15.6)                                            
 102  FORMAT(7x,I15)                                            

      
!      print *, filename
      VXGSE=-400.D0  !  GSE COMPONENTS OF SOLAR WIND VELOCITY VECTOR; THIS PARTICULAR CHOICE
      VYGSE=   0.D0  !   IS MADE IN ORDER TO MAKE THE GSW SYSTEM IDENTICAL TO THE STANDARD GSM
      VZGSE=   0.D0


      call RECALC_08 (IYR,IDY,IHR,IMN,ISC,VXGSE,VYGSE,VZGSE) ! CALCULATES TILT ANGLE AND
C                                                                  UPDATES MAIN FIELD COEFFICIENTS      

      else
        print *, 'TS07d error: No Coeff files exist for ',iyear,idoy
        print *, 'TS07d error: filename: ',filename,' does not exist'
        ifail=-1
      endif
      end

      ! Adapted for IRBEM by A. C. Kellerman, needs to be called on date
      ! update. TODO: check if we already loaded the current files and
      ! skip the load if so, as they currently have a cadence of 5
      ! minutes
      
      subroutine INIT_TS07D_COEFFS (iyear,idoy,iut)

      INCLUDE 'ts07d.inc'

      INTEGER*4 iyear,idoy
      REAL*8 A07,iut

      INTEGER IYR,IDY,IHR,IMN,ISC
      INTEGER TS7LEN,ierr,stat,statb(13) ! file i/o
      LOGICAL OK ! I/O

      CHARACTER*18 PAR_FNAME
      CHARACTER*200 filename,FMT
      CHARACTER*80 TS7DIR
      CHARACTER*255 ts07d_env

      PARAMETER (NTOT=101)

      COMMON /A07/ A07(NTOT)
      CALL GETENV('TS07_DATA_PATH', ts07d_env)
      if (LEN_TRIM(ts07d_env).ne.0) then
        TS7DIR=TRIM(ts07d_env)
      else
        TS7DIR=TRIM(TS07D_DIR)
      endif
      i=len(TS7DIR)
      do while (TS7DIR(i:i) == ' ')
        i = i - 1  
      enddo
      TS7LEN=i
      TS7DIR=TS7DIR(1:TS7LEN)//'/TS07D'
      TS7LEN = TS7LEN+6

      IYR=iyear
      IDY=idoy
      IHR=floor(iut/3600d0)
      IMN=floor((iut/3600d0)-IHR) * 60 
      ISC=iut-IHR*3600-IMN*60

      MODMIN=MOD(IMN,5)
      WRITE (PAR_FNAME,'(I4,A1,I0.3,A1,I0.2,A1,I0.2,A4)'),
     +iyear,'_',idoy,'_',IHR,'_',IMN-MODMIN,'.par' 

C     TS7LEN is the directory length
      WRITE(FMT,'("(A", I0, ",A8,A18)")') TS7LEN
      write(filename,FMT), TS7DIR(1:TS7LEN),'/Coeffs/',PAR_FNAME

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
 100  FORMAT(G15.6)                                            !  MAKE SURE TO MODIFY THE PATH
      CLOSE(1)                                        
!      print *, filename
      VXGSE=-400.  !  GSE COMPONENTS OF SOLAR WIND VELOCITY VECTOR; THIS PARTICULAR CHOICE
      VYGSE=   0.  !   IS MADE IN ORDER TO MAKE THE GSW SYSTEM IDENTICAL TO THE STANDARD GSM
      VZGSE=   0.

      call RECALC_08 (IYR,IDY,IHR,IMN,ISC,VXGSE,VYGSE,VZGSE) ! CALCULATES TILT ANGLE AND
C                                                                  UPDATES MAIN FIELD COEFFICIENTS      
      else
        print *, filename,' does not exist'
        stop
      endif
      end

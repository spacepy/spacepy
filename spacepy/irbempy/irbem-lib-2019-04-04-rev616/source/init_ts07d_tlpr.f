c     Adapted for IRBEM by A. C. Kellerman, loads tail par files
c     these only need to be loaded once
c     TODO: introduce a parameter when this is loaded, so we don't
c     call it in the GET_FIELD1 and GET_MULTI1 scripts for each call

c     revised Nov 8 2017 in preparation for 2017 TS07d model inclusion - A.C.K
      
      SUBROUTINE INIT_TS07D_TLPR

      IMPLICIT NONE

      REAL*8    ::  TSS,TSO,TSE

      INTEGER i,TS7LEN,IREAD,KREAD,KK ! dir length 

      CHARACTER*200 filename,FFMT
      CHARACTER*255 TS7DIR
      CHARACTER*255 ts07d_env

      COMMON /TSS/ TSS(80,5) ! tail pars
      COMMON /TSO/ TSO(80,5,4)
      COMMON /TSE/ TSE(80,5,4)

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

      DO 1001 IREAD=1,5
        WRITE(FFMT,'("(A", I0, ",A20,I0,A4)")') TS7LEN
        WRITE(filename,FFMT) TS7DIR(1:TS7LEN),'/TAIL_PAR/tailamebhr'
     *,IREAD,'.par'
      OPEN (UNIT=1,FILE=filename)
      READ (1,200) (TSS(KK,IREAD),KK=1,80)
 200  FORMAT(G17.10)
 1001 CLOSE(1)

      DO 1002 IREAD=1,5
        DO 1003 KREAD=1,4
        WRITE(FFMT,'("(A", I0, ",A21,I0,I0,A4)")') TS7LEN
        WRITE(filename,FFMT)TS7DIR(1:TS7LEN),'/TAIL_PAR/tailamhr_o_'
     *,IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=filename)
      READ (1,200) (TSO(KK,IREAD,KREAD),KK=1,80)
 1003   CONTINUE
 1002 CLOSE(1)

      DO 1004 IREAD=1,5
        DO 1005 KREAD=1,4
        WRITE(FFMT,'("(A", I0, ",A21,I0,I0,A4)")') TS7LEN
        WRITE(filename,FFMT)TS7DIR(1:TS7LEN),'/TAIL_PAR/tailamhr_e_'
     *,IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=filename)
      READ (1,200) (TSE(KK,IREAD,KREAD),KK=1,80)
 1005   CONTINUE
 1004 CLOSE(1)

      end

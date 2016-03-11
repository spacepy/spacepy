      ! Adapted for IRBEM by A. C. Kellerman, loads tail par files
      ! these only need to be loaded once
      ! TODO: introduce a parameter when this is loaded, so we don't
      ! call it in the GET_FIELD1 and GET_MULTI1 scripts for each call
      
      subroutine INIT_TS07D_TLPR

      INCLUDE 'ts07d.inc'

      REAL*8 TSS,TSO,TSE ! tail pars

      INTEGER i,TS7LEN ! dir length 

      CHARACTER*200 filename,FMT
      CHARACTER*80 TS7DIR
      CHARACTER*255 ts07d_env

      COMMON /TSS/ TSS(80,5) ! tail pars
      COMMON /TSO/ TSO(80,5,4)
      COMMON /TSE/ TSE(80,5,4)

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

      DO 1001 IREAD=1,5
        WRITE(FMT,'("(A", I0, ",A23,I0,A4)")') TS7LEN
        WRITE(filename,FMT) TS7DIR(1:TS7LEN),'/TSG_DYN_PAR/tailamebhr'
     *,IREAD,'.par'
      OPEN (UNIT=1,FILE=filename)
      READ (1,200) (TSS(KK,IREAD),KK=1,80)
 200  FORMAT(G17.10)
 1001 CLOSE(1)

      DO 1002 IREAD=1,5
        DO 1003 KREAD=1,4
        WRITE(FMT,'("(A", I0, ",A24,I0,I0,A4)")') TS7LEN
        WRITE(filename,FMT)TS7DIR(1:TS7LEN),'/TSG_DYN_PAR/tailamhr_o_'
     *,IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=filename)
      READ (1,200) (TSO(KK,IREAD,KREAD),KK=1,80)
 1003   CONTINUE
 1002 CLOSE(1)

      DO 1004 IREAD=1,5
        DO 1005 KREAD=1,4
        WRITE(FMT,'("(A", I0, ",A24,I0,I0,A4)")') TS7LEN
        WRITE(filename,FMT)TS7DIR(1:TS7LEN),'/TSG_DYN_PAR/tailamhr_e_'
     *,IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=filename)
      READ (1,200) (TSE(KK,IREAD,KREAD),KK=1,80)
 1005   CONTINUE
 1004 CLOSE(1)

      end

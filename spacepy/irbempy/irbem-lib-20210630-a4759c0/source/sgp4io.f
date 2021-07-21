!***************************************************************************************************
! Copyright 2004 D. Vallado
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
*     ----------------------------------------------------------------
*
*                               sgp4io.for
*
*    this file contains miscallaneous functions to read two line element
*    sets. while not formerly part of the sgp4 mathematical theory, they are
*    required for practical implemenation.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2004
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              15 aug 06  david vallado
*                           update mfe for verif time steps, consts, mfe
*    changes :
*              15 dec 05  david vallado
*                           original baseline
*       ----------------------------------------------------------------


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE TWOLINE2RVSGP4
*
*  this function converts the two line element set character string data to
*    variables and initializes the sgp4 variables. several intermediate varaibles
*    and quantities are determined. note that the result is a structure so multiple
*    satellites can be processed simaltaneously without having to reinitialize. the
*    verification mode is an important option that permits quick checks of any
*    changes to the underlying technical theory. this option works using a
*    modified tle file in which the start, stop, and delta time values are
*    included at the end of the second line of data. this only works with the
*    verification mode. the catalog mode simply propagates from -1440 to 1440 min
*    from epoch and is useful when performing entire catalog runs.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs        :
*    Numsats     - Number of satellites processed. It also becomes the record
*                  number for each satellite
*    typerun     - type of run                    verification 'V', catalog 'C', 'N'
*    whichconst  - which set of constants to use  72, 84
*
*  outputs       :
*    Code        - EOF indicator. Code = 999 when EOF reached
*    startmfe    - starttime of simulation,       min from epoch
*    stopmfe     - stoptime of simulation,        min from epoch
*    deltamin    - time step                      min
*
*  coupling      :
*    days2mdhms  - conversion of days to month, day, hour, minute, second
*    jday        - convert day month year hour minute second into julian date
*    sgp4init    - initialize the sgp4 variables
*
*  Files         :
*    Unit 10     - test.elm        input 2-line element set file
*    Unit 11     - test.bak        output file
*    Unit 15     - sgp4rec.bak     temporary file of record for 2 line element sets
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE TwoLine2RVSGP4 ( NumSats, Typerun, whichconst, 
     &                            startmfe, stopmfe, deltamin, Code )
        IMPLICIT NONE
        Character Typerun
        INTEGER*4 Code, NumSats, whichconst,i,iline,EpochYr3
        REAL*8 startmfe, stopmfe, deltamin

* ----------------------------  Locals  -------------------------------
        REAL*8 J2,TUMin,TwoPi,Rad,RadiusEarthKm,VKmPerSec, xke
        REAL*8 BC,EPDay, sec, de2ra, xpdotp, j3, j4, j3oj2 
	REAL*8 EPDay3,JDSatEpoch3
        INTEGER*4 Yr,Mon,Day,Hr,Minute,  ICrdno,nexp,bexp, error
        CHARACTER Show
        Character*130 LongStr1,LongStr2,LongStr3

        INCLUDE 'SGP4.CMN'

        COMMON /DebugHelp/ Help
        CHARACTER Help

        ! --------------------  Implementation   ----------------------
        Show = 'N'
        TwoPi         =    6.28318530717959D0
        Rad           =   57.29577951308230D0
        DE2RA         =    0.01745329251994330D0
        xpdotp        =  1440.0D0 / twopi  ! 229.1831180523293
        VKmPerSec     =  RadiusEarthKm * xke / 60.0D0

        CALL getgravconst( whichconst, tumin, radiusearthkm, xke, j2, 
     &       j3, j4, j3oj2 )

c        make sure the main program opens this file, otherwise do so here
c        ! store results in a temporary file of record
c        OPEN(15,FILE='Sgp4Rec.bak',ACCESS='DIRECT', FORM='UNFORMATTED',
c     &       RECL=1000,STATUS='UNKNOWN')

* ----------------- READ THE FIRST LINE OF ELEMENT SET ----------------
        Code = 0

        LongStr1 = ' '
   50   READ(10,'(a130)',END=999) LongStr1
        IF(LongStr1(1:1) .eq. '#') GOTO 50 ! Commented line of text, skip

        READ(LongStr1,500) ICRDNO,SatNum,SatName,EpochYr,EpDay,
     &                       NDot,NDDot,nexp,BStar,bexp,EPHTYP,ELNO
  500   FORMAT( I1,1X,I5,1X,A10,I2,D12.0,1X,D10.0,1X,
     &          F6.5,I2,1X,F6.5,I2,1X,I1,1X,I4 )

* ----------- READ THE SECOND LINE OF ELEMENT SET AND TIME ------------
        LongStr2 = ' '
   51   READ(10,'(a130)',END=999) LongStr2
        IF(LongStr2(1:1) .eq. '#') GOTO 51 ! Commented line of text, skip

        IF (Typerun.eq.'V') THEN
          READ(LongStr2,502) ICRDNO,Inclo,nodeo,Ecco,Argpo,Mo,No,REVI,
     &              startmfe, stopmfe, DeltaMin
         else
          READ(LongStr2,501) ICRDNO,Inclo,nodeo,Ecco,Argpo,Mo,No,REVI
         endif
  501   FORMAT( I1,7X,D8.0,1X,D8.0,1X,F7.7,1X,D8.0,1X,D8.0,1X,D11.0,I5)
  502   FORMAT( I1,7X,D8.0,1X,D8.0,1X,F7.7,1X,D8.0,1X,D8.0,1X,D11.0,I5,
     &          1X,F12.6,F12.6,F12.6 )

* Added Seb to process from one element to the next
*
        if (typerun .EQ. 'Y') then
           iline=0
   52      READ(10,'(a130)',END=999) LongStr3
           IF(LongStr3(1:1) .eq. '#') THEN
	      iline=iline+1
	      GOTO 52 ! Commented line of text, skip
           ENDIF
           READ(LongStr3,503) EpochYr3,EpDay3
  503      FORMAT( 18x,I2,D12.0)
           ! Temporary year fix
           IF (EpochYr.lt.50) THEN
               Yr = EpochYr + 2000
           ELSE
               Yr = EpochYr + 1900
           ENDIF
           CALL Days2MDHMS( Yr,EpDay, Mon,Day,Hr,Minute,Sec )
           CALL JDAY ( Yr,Mon,Day,Hr,Minute,Sec,  JDSatEpoch )

           IF (EpochYr3.lt.50) THEN
               Yr = EpochYr3 + 2000
           ELSE
               Yr = EpochYr3 + 1900
           ENDIF

           CALL Days2MDHMS( Yr,EpDay3, Mon,Day,Hr,Minute,Sec )
           CALL JDAY ( Yr,Mon,Day,Hr,Minute,Sec,  JDSatEpoch3 )
	   startmfe=0.D0
	   stopmfe=(JDSatEpoch3-JDSatEpoch)*1440.D0
	   if (stopmfe .EQ. 0.D0) stopmfe=1440.D0
c	   deltamin=1.D0
	   if (JDSatEpoch3.le.JDSatEpoch) then
	      READ(10,*)
	      iline=0
	      GOTO 52
	   endif
	   DO I=1,iline+1
	      backspace(10)
	   ENDDO
	ENDIF
*
* End Added Seb to process from one element to the next
*
* ---------------------- CONVERT TO INTERNAL UNITS --------------------
* ---- RADIANS, DISTANCE IN EARTH RADII, AND VELOCITY IN ER/KEMIN) ----
        NDDot  = NDDot * 10.0D0**Nexp
        NDot   = NDot / (XPDOTP*1440)
        NDDot  = NDDot / (XPDOTP*1440*1440)
        BStar  = BStar * 10.0D0**Bexp

        No     = No / XPDOTP
        a      = (No*TUMin)**(-2.0D0/3.0D0)
        Inclo  = Inclo  * DE2RA
        nodeo  = nodeo * DE2RA
        Argpo  = Argpo * DE2RA
        Mo     = Mo   * DE2RA
                                                                        
        IF (DABS(Ecco-1.0D0) .gt. 0.000001D0) THEN
            Altp= (a*(1.0D0-Ecco))-1.0D0
            Alta= (a*(1.0D0+Ecco))-1.0D0
          ELSE
            Alta= 999999.9D0
            Altp= 2.0D0* (4.0D0/(No*No)**(1.0D0/3.0D0))
          ENDIF

        ! ---- Ballistic Coefficient ----
        IF (DABS(BStar) .gt. 0.00000001D0) THEN
            BC= 1.0D0/(12.741621D0*BStar)
          ELSE
            BC= 1.111111111111111D0
          ENDIF

        ! ----------------------------------------------------------------
        ! find sgp4epoch time of element set
        ! remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
        ! and minutes from the epoch (time)
        ! ----------------------------------------------------------------

        ! ---- input start stop times manually
c        IF ((typerun.ne.'V').and.(typerun.ne.'C')) THEN
c            write(*,*) 'input start min from epoch'
c            read(*,*) startmfe
c            write(*,*) 'input stop min from epoch'
c            read(*,*) stopmfe
c            write(*,*) 'input time step in minutes'
c            read(*,*) deltamin
c          ENDIF

       ! ---- perform complete catalog evaluation
c       if (typerun .eq. 'C') THEN
c           startmfe = -1440.0D0
c           stopmfe  =  1440.0D0
c           deltamin =    20.0D0
c         ENDIF

        ! Temporary year fix
        IF (EpochYr.lt.50) THEN
            Yr = EpochYr + 2000
          ELSE
            Yr = EpochYr + 1900
          ENDIF

        CALL Days2MDHMS( Yr,EpDay, Mon,Day,Hr,Minute,Sec )
        CALL JDAY ( Yr,Mon,Day,Hr,Minute,Sec,  JDSatEpoch )

* ------------------- MAKE INITIAL PREDICTION AT EPOCH ----------------
        ! 2433281.5 - 2400000.5 = 33281.0, thus time from 1950
        CALL SGP4Init( whichconst,
     &                 SatNum,BStar, Ecco, JDSatEpoch-2433281.5D0,
     &                 Argpo,Inclo,Mo,No, nodeo, Error )

        ! ---- Write common block of data into file of record ----
c        WRITE(15,Rec=NumSats) SatName,
c     &          SatNum, ELNO  , EPHTYP, REVI  , EpochYr,
c     &          BStar , Ecco  , Inclo , nodeo, Argpo , No    , Mo    ,
c     &          NDot  , NDDot ,
c     &          alta  , altp  , a     ,
c     &          DeltaMin, JDSatEpoch, EpochDays,
c     &          Isimp , Init  , Method,
c     &          Aycof , CON41 , Cc1   , Cc4   , Cc5   , D2    , D3    ,
c     &          D4    , Delmo , Eta   , ArgpDot,Omgcof, Sinmao,
c     &          T2cof , T3cof , T4cof , T5cof , X1mth2, X7thm1, MDot  ,
c     &          nodeDot,Xlcof, Xmcof , Xnodcf,
c     &          D2201 , D2211 , D3210 , D3222 , D4410 , D4422 , D5220 ,
c     &          D5232 , D5421 , D5433 , Dedt  , Del1  , Del2  , Del3  ,
c     &          Didt  , Dmdt  , Dnodt , Domdt , E3    , Ee2   , Peo   ,
c     &          Pgho  , Pho   , Pinco , Plo   , Se2   , Se3   , Sgh2  ,
c     &          Sgh3  , Sgh4  , Sh2   , Sh3   , Si2   , Si3   , Sl2   ,
c     &          Sl3   , Sl4   , GSTo  , Xfact , Xgh2  , Xgh3  , Xgh4  ,
c     &          Xh2   , Xh3   , Xi2   , Xi3   , Xl2   , Xl3   , Xl4   ,
c     &          Xlamo , Zmol  , Zmos  , Atime , Xli   , Xni   , IRez

        IF(Error .GT. 0) THEN
            WRITE( *,*) '# *** SGP4 Model Error ***',Error
          ENDIF

c      write tle output details
c      INCLUDE 'debug8.for'

        ! ---- Fix to indicate end-of-file
        GOTO 1000
  999   Code = 999
 1000   CONTINUE
       if (typerun.EQ.'C') Code = 999
       RETURN
       END  !       SUBROUTINE TwoLine2RVSGP4





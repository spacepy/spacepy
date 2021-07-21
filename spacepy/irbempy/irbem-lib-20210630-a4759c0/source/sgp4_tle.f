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
* ---------------------------------------------------------------------
*http://www.centerforspace.com/downloads/files/pubs/AIAA%202006-6753%20Revisiting%20Spacetrack%20Report%203.pdf
*                              TESTFOR.FOR
*
*  this program tests the sgp4 propagator.
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2004
*                            by david vallado
*
*     (w) 719-573-2600, email dvallado@agi.com
*     *****************************************************************
*  current :
*            08 jan 07  Sebastien Bourdarie
*                       include into onera_desp library
*            15 aug 06  david vallado
*                        update mfe for verification time steps, constants
*  changes :
*            20 jul 05  david vallado
*                         fixes for paper, corrections from paul crawford
*             7 jul 04  david vallado
*                         fix record file and get working
*            14 may 01  david vallado
*                         2nd edition baseline
*                   97  nasa
*                         internet version
*                   80  norad
*                         original baseline
*     *****************************************************************
*
*  Files         :
*    Unit 10     - input elm file  input file for element sets
*    Unit 11     - sgp4test.out    output file
*    Unit 14     - sgp4test.dbg    debug output file
*    Unit 15     - sgp4rec.bak     temporary file of record for 2 line element sets
*
*  Uses object files:
*    Sgp4ext,
*    Sgp4io,
*    Sgp4unit

        SUBROUTINE SGP4_TLE1(runtype,startsfe,stopsfe,deltasec
     &     ,InFileByte,strlenIn,OutFileByte,strlenOut)
        IMPLICIT NONE
	INTEGER*4 runtype
        INTEGER*4 strlenIn,strlenOut
        INTEGER*4 Code, NumSats, TotalNumSats, k, error, whichconst
        INTEGER*4 i,j, Year,yr,mon,day,hr,min,doy,GET_DOY,isec
c
        BYTE      InFileByte(strlenIn),OutFileByte(strlenOut)
c
        Character typerun
        Character*500 InFileName,OutFileName

        Character*3 MonStr,Monthtitle(12)
        Real*8 ro(3),vo(3), startmfe, stopmfe, deltamin
        Real*8 startsfe, stopsfe, deltasec

        REAL*8 p, ecc, incl, node, argp, nu, m,arglat,truelon,lonper

* ----------------------------  Locals  -------------------------------
        REAL*8 J2,TwoPi,Rad,RadiusEarthKm,VKmPerSec, xke,
     &         de2ra, xpdotp, T, sec, JD, pi, j3, j4, j3oj2, tumin
	REAL*8  secs,xGEI(3),xOUT(3),alti,lati,longi,dec_y


        INCLUDE 'SGP4.CMN'

        COMMON /DebugHelp/ Help
        CHARACTER Help
        Help = 'N'

* ------------------------  Implementation   --------------------------
c   Start modif Seb  Dec 12, 2006
c
c        typerun = 'C' ! complete catlog
c        typerun = 'V' ! validation only, header satnum and r and v tsince
c        write(*,*) 'Input typerun - V, C, N, Y'
c        read(*,*) typerun
c changed input times from minute to seconds by SB. The first
c 3 lines change it back to minutes to be compliant with SGP4 code.
        startmfe=startsfe/60.0D0
        stopmfe=stopsfe/60.0D0
        deltamin=deltasec/60.0D0
        if (runtype .eq. 0) typerun='Y'
        if (runtype .eq. 1) then
	   typerun='C'
	   if (stopmfe .LE. startmfe) stopmfe=startmfe+1440.0D0
	endif
	if (deltamin .LE. 0.D0) deltamin=1.D0
c
        DO i=1,strlenIn
	   InFileName(i:i)=char(InFileByte(i))
	ENDDO
        InFileName=InFileName(1:strlenIn) ! added by TPO
        DO i=1,strlenOut
	   OutFileName(i:i)=char(OutFileByte(i))
	ENDDO
        OutFileName=OutFileName(1:strlenOut) ! added by TPO
c
c        write(*,*) 'Input whichconst - 721, 72, 84'
c        read(*,*) whichconst
        whichconst=84
c
c  End modif Seb Dec 12, 2006
c
        TwoPi         =    6.28318530717959D0
        Rad           =   57.29577951308230D0
        DE2RA         =    0.01745329251994330D0
        pi            =    3.14159265358979D0
        xpdotp        =  1440.0 / (2.0 *pi)  ! 229.1831180523293D0

        ! sgp4fix identify constants and allow alternate values
        CALL getgravconst( whichconst, tumin, radiusearthkm, xke, j2,
     &       j3, j4, j3oj2 )
        VKmPerSec     =  RadiusEarthKm * xke/60.0D0

        MonthTitle( 1)= 'Jan'
        MonthTitle( 2)= 'Feb'
        MonthTitle( 3)= 'Mar'
        MonthTitle( 4)= 'Apr'
        MonthTitle( 5)= 'May'
        MonthTitle( 6)= 'Jun'
        MonthTitle( 7)= 'Jul'
        MonthTitle( 8)= 'Aug'
        MonthTitle( 9)= 'Sep'
        MonthTitle(10)= 'Oct'
        MonthTitle(11)= 'Nov'
        MonthTitle(12)= 'Dec'

        ! ---------------- Setup files for operation ------------------
        ! 10 input 2-line element set file
c        Write(*,*) 'Input elset filename '
c        Read(*,*) InFileName
        OPEN(10,FILE = InFileName ,STATUS='OLD',
     &          ACCESS = 'SEQUENTIAL' )

        ! 11 output file
c        IF (typerun.eq.'C') THEN
c            OPEN(11,FILE = 'tforall.out' ,STATUS='UNKNOWN',
c     &              ACCESS = 'SEQUENTIAL' )
c          ELSEIF (typerun.eq.'V') THEN
c                OPEN(11,FILE = 'tforver.out' ,STATUS='UNKNOWN',
c     &              ACCESS = 'SEQUENTIAL' )
c              ELSE
c                OPEN(11,FILE = 'tfor.out' ,STATUS='UNKNOWN',
c     &              ACCESS = 'SEQUENTIAL' )
c              ENDIF
         OPEN(11,FILE = OutFileName ,STATUS='UNKNOWN',
     &           ACCESS = 'SEQUENTIAL' )

c        OPEN(14,FILE = 'sgp4test.dbg' ,STATUS='UNKNOWN',
c     &          ACCESS = 'SEQUENTIAL' )

        ! ----- 15 temporary file of record for 2 line element sets ---
c        OPEN(15,FILE = 'Sgp4Rec.bak', ACCESS = 'DIRECT',
c     &          FORM = 'UNFORMATTED', RECL = 1100, STATUS = 'UNKNOWN' )

        ! ----------------- Test simple propagation -------------------
        NumSats = 0
        DOWHILE (Code.ne.999)
            Numsats = NumSats + 1
            CALL TwoLine2RVSGP4 ( NumSats,typerun,whichconst,startmfe,
     &                            stopmfe,deltamin,Code )

c            Write(11,*) '',SatNum,' xx'
c            Write(*,*) SatNum
            ! write out epoch value
            T = 0.0D0
            CALL SGP4 ( whichconst, T, Ro, Vo, Error )
c            WRITE( 11,800 ) T, ro(1),ro(2),ro(3),vo(1),vo(2),vo(3)

            ! now initialize time variables
            T      = startmfe

            ! check so the first value isn't written twice
            IF ( DABS(T).gt.0.00000001D0 ) THEN
                T = T - DeltaMin
              ENDIF

            DOWHILE ( (T+ DeltaMin.le.stopmfe).and.(Error.eq.0) )
                T = T + DeltaMin
c                IF (T.gt.stopmfe) THEN
c                    T = stopmfe
c                  ENDIF
                CALL SGP4 ( whichconst, T, Ro, Vo, Error )

                IF (Error .gt. 0 ) THEN
                    Write(*,*) '# Error in SGP4 .. ',
     &                     Error
                  ENDIF

                IF ( error .eq. 0) THEN

c                 IF ((typerun.ne.'V').and.(typerun.ne.'C')) THEN
                 IF (typerun.ne.'V') THEN
                     JD = JDSatEpoch + T/1440.0D0
                     CALL INVJDAY( JD, Year,Mon,Day,Hr,Min, Sec )
                     IF (Year.ge.2000) THEN
                         Yr = Year - 2000
                       ELSE
                         Yr = Year - 1900
                       ENDIF
                     MonStr = MonthTitle( Mon )
c
c   Start modif Seb  Dec 12, 2006
c Force output to be Altitude, latitude, longitude rather than ECI (or GEI)
c
		     DOY = GET_DOY(Year,Mon, Day)
		     secs=3600.D0*Hr+60.D0*Min+1.D0*Sec
		     xGEI(1)=ro(1)/6371.2D0
		     xGEI(2)=ro(2)/6371.2D0
		     xGEI(3)=ro(3)/6371.2D0
             call coord_trans1(13,0,Year,DOY,secs,xGEI,xOUT)
             alti=xOUT(1)
             lati=xOUT(2)
             longi=xOUT(3)
c		     call gei2geo1(Year,DOY,secs,xGEI,xGEO)
c		     call geo_gdz(xGEO(1),xGEO(2),xGEO(3),
c     &                    lati,longi,alti)
c                     WRITE( 11,'(F17.8,I4,1x,A3,I3,I3,A1,I2,A1,F9.6,
c     &                         3F17.6,3F17.8)' )
c     &                    JD,Day,MonStr,Yr,Hr,':',Min,':',
c     &                   Sec, ro(1),ro(2),ro(3),vo(1),vo(2),vo(3)
                     isec=INT(sec)
                     CALL DATE_AND_TIME2DECY(Year,Mon, Day,
     &                    hr,min,isec,Dec_y)
                     WRITE( 11,'(I2.2,A1,I2.2,A1,I4,1x,I3,A1,I2,A1,F9.6,
     &                         1x,f13.8,1x,3F17.6)' )
     &                    Day,'/',Mon,'/',Year,Hr,':',Min,':',
     &                   Sec, Dec_y,alti,lati,longi
c
c   End modif Seb  Dec 12, 2006
c
                   ELSE
                     WRITE( 11,800 )
     &                   T,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3)
  800  FORMAT(F17.8,3F17.8,3(1X,F14.9))

c                     call rv2coe(ro, vo, p, a, ecc, incl, node, argp,
c     &                           nu, m, arglat, truelon, lonper )
c
c                     write(11,801) T ,ro(1),ro(2),ro(3),vo(1),vo(2),
c     &                             vo(3),a, ecc, incl*rad, node*rad,
c     &                             argp*rad, nu*rad, m*rad

  801  FORMAT(F17.8,3F17.8,3(1X,F13.9),f15.6,f9.6,f11.5,f11.5,f11.5,
     &          f11.5,f11.5)

                  ENDIF ! typerun

               ENDIF ! if error

             ENDDO ! propagating the case

          ENDDO ! while through file
	  close(11)

      END



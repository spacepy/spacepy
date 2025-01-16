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

      SUBROUTINE SGP4_ORB1(Yr,Mon,Day,Hr,Minute,Sec,inclination,
     &     perigee,apogee,AscendingNode,Argument_Perigee,MeanAnomaly,
     &     startsfe,stopsfe,deltasec,sysaxes,OutputArray)
        IMPLICIT NONE
        INCLUDE 'ntime_max.inc'
c   TBG, 3.27.07:  added input argument sysaxes for internal rotation.
        INTEGER*4 k, error, whichconst,strlenOut,istep,sysaxes
        INTEGER*4 i,j, Year,yr,mon,day,hr,min,doy,GET_DOY,isec
	INTEGER*4 Minute,numberOfOutputLines
c
c        BYTE      OutFileByte(strlenOut)
c
        Real*8 ro(3),vo(3), startmfe, stopmfe, deltamin
        REAL*8 p, ecc, incl, node, argp, nu, m,arglat,truelon,lonper

* ----------------------------  Locals  -------------------------------
        REAL*8 J2,TwoPi,Rad,RadiusEarthKm,VKmPerSec, xke,
     &         de2ra, xpdotp, T, sec, JD, pi, j3, j4, j3oj2, tumin
	REAL*8  secs,xGEI(3),xGEO(3),xOUT(3),alti,lati,longi,dec_y
        REAL*8 inclination,perigee,apogee,AscendingNode
	REAL*8 Argument_Perigee,MeanAnomaly,BC,mu
        Real*8 startsfe, stopsfe, deltasec
        REAL*8 outputArray(10,ntime_max)
c
c
c        Character*500 OutFileName


        INCLUDE 'SGP4.CMN'

        COMMON /DebugHelp/ Help
        CHARACTER Help
        Help = 'N'

* ------------------------  Implementation   --------------------------
c changed input times from minute to seconds by SB. The first
c 3 lines change it back to minutes to be compliant with SGP4 code.
        startmfe=startsfe/60.0D0
        stopmfe=stopsfe/60.0D0
        deltamin=deltasec/60.0D0
        numberOfOutputLines = (stopmfe - startmfe) / deltamin
c        write(*,*) 'Input whichconst - 721, 72, 84'
c        read(*,*) whichconst
        whichconst=84
	SatNum=00000
c
        TwoPi         =    6.28318530717959D0
        Rad           =   57.29577951308230D0
        DE2RA         =    0.01745329251994330D0
        pi            =    3.14159265358979D0
        xpdotp        =  1440.0 / (2.0 *pi)  ! 229.1831180523293D0
        mu     = 398600.8D0            ! in km3 / s2

        ! sgp4fix identify constants and allow alternate values
        CALL getgravconst( whichconst, tumin, radiusearthkm, xke, j2,
     &       j3, j4, j3oj2 )
        VKmPerSec     =  RadiusEarthKm * xke/60.0D0

c        DO i=1,strlenOut
c	   OutFileName(i:i)=char(OutFileByte(i))
c	ENDDO
c        OutFileName=OutFileName(1:strlenOut) ! added by TPO
c        OPEN(11,FILE = OutFileName ,STATUS='UNKNOWN',
c     &      ACCESS = 'SEQUENTIAL' )
        ! ----------------- Test simple propagation -------------------

	Inclo=inclination
	nodeo=AscendingNode
	Ecco=(apogee-perigee)/(apogee+perigee+2.D0*RadiusEarthKm)
	Argpo=Argument_Perigee
	Mo=MeanAnomaly
Cwrong	No=(24.D0/1.4D0)/(SQRT(
Cwrong     &  (0.5D0*(apogee+perigee+2.D0*RadiusEarthKm)/RadiusEarthKm)**3))
	a = 0.5D0*(apogee+perigee+2.D0*RadiusEarthKm)
	No=(43200.D0/PI) * SQRT(mu/a**3)
	a = a/RadiusEarthKm

	BSTAR=0.D0
*
* ---------------------- CONVERT TO INTERNAL UNITS --------------------
* ---- RADIANS, DISTANCE IN EARTH RADII, AND VELOCITY IN ER/KEMIN) ----
        No     = No / XPDOTP
Cxxxx        a      = (No*TUMin)**(-2.0D0/3.0D0)
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
*
        ! ----------------------------------------------------------------
        ! find sgp4epoch time of element set
        ! remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
        ! and minutes from the epoch (time)
        ! ----------------------------------------------------------------
*
        CALL JDAY ( Yr,Mon,Day,Hr,Minute,Sec,  JDSatEpoch )
*
* ------------------- MAKE INITIAL PREDICTION AT EPOCH ----------------
        ! 2433281.5 - 2400000.5 = 33281.0, thus time from 1950
*
        CALL SGP4Init( whichconst,
     &                 SatNum,BStar, Ecco, JDSatEpoch-2433281.5D0,
     &                 Argpo,Inclo,Mo,No, nodeo, Error )
        IF(Error .GT. 0) WRITE( *,*) '# *** SGP4 Model Error ***',Error

*
* ------------------- PROPAGATE ----------------
        ! 2433281.5 - 2400000.5 = 33281.0, thus time from 1950
*
        T = 0.0D0
        CALL SGP4 ( whichconst, T, Ro, Vo, Error )

        ! now initialize time variables
        T      = startmfe

        ! check so the first value isn't written twice
        IF ( DABS(T).gt.0.00000001D0 ) THEN
            T = T - DeltaMin
        ENDIF

        istep = 0
        DOWHILE ( (T.lt.stopmfe).and.(Error.eq.0) )
           istep = istep + 1
           T = T + DeltaMin
           IF (T.gt.stopmfe) THEN
               T = stopmfe
           ENDIF
           CALL SGP4 ( whichconst, T, Ro, Vo, Error )

           IF (Error .gt. 0 ) THEN
               Write(*,*) '# Error in SGP4 .. ',Error
           ENDIF

           IF ( error .eq. 0) THEN
              JD = JDSatEpoch + T/1440.0D0
              CALL INVJDAY( JD, Year,Mon,Day,Hr,Min, Sec )
              IF (Year.ge.2000) THEN
                 Yr = Year - 2000
              ELSE
                 Yr = Year - 1900
              ENDIF
c
	      DOY = GET_DOY(Year,Mon, Day)
	      secs=3600.D0*Hr+60.D0*Min+1.D0*Sec
	      xGEI(1)=ro(1)/6371.2D0
	      xGEI(2)=ro(2)/6371.2D0
	      xGEI(3)=ro(3)/6371.2D0
c     TBG 3.27.07:  Exchanged two rotations for one call to coord_trans
c     sysaxes now passed into sgp4_orb1
c     SGP4 outputs are given in TEME (True Equator Mean Equinox) coordinate
c     system (so 13 index in the lib for corrdinates conversions)
              call coord_trans1(13,sysaxes,Year,DOY,secs,xGEI,xOUT)
c	      call gei2geo1(Year,DOY,secs,xGEI,xGEO)   ! commented out by TBG
c	      call geo_gdz(xGEO(1),xGEO(2),xGEO(3),lati,longi,alti)  ! same here
              isec=INT(sec)
              CALL DATE_AND_TIME2DECY(Year,Mon, Day,
     &                    hr,min,isec,Dec_y)
c              WRITE( 11,'(I2.2,A1,I2.2,A1,I4,1x,I3,A1,I2,A1,F9.6,
c     &                         1x,f13.8,1x,3F17.6)' )
c     &                    Day,'/',Mon,'/',Year,Hr,':',Min,':',
c     &                   Sec, Dec_y,alti,lati,longi
c     Now assign trajectory info to the outputArray
              outputArray(1,istep) = Day
              outputArray(2,istep) = Mon
              outputArray(3,istep) = Year
              outputArray(4,istep) = Hr
              outputArray(5,istep) = Min
              outputArray(6,istep) = Sec
              outputArray(7,istep) = Dec_y
c              outputArray(8,istep) = alti
c              outputArray(9,istep) = lati
c              outputArray(10,istep) = longi
              outputArray(8,istep) = xOUT(1)  ! sending back the positions rotated
              outputArray(9,istep) = xOUT(2)  ! to the correct sysaxes specified.
              outputArray(10,istep) = xOUT(3)  !  TBG, 3.27.07
           ENDIF ! if error

        ENDDO ! propagating the case

c	close(11)
        AscendingNode  = nodeo / DE2RA
        Argument_Perigee  = Argpo / DE2RA
        MeanAnomaly     = Mo   / DE2RA

      END



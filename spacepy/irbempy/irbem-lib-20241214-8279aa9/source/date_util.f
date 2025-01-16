!***************************************************************************************************
! Copyright , 2005 S. Bourdarie
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
! 
!----------------------------------------------------------------------      
! NAME:
!	JULDAY
!
! PURPOSE:
!	Calculate the Julian Day Number for a given month, day, and year.
!
!
! CALLING SEQUENCE:
!	Result = JULDAY(Year,Month, Day)
!
! INPUTS:
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year.Year parameters must be valid
!               values from the civil calendar.  Years B.C.E. are represented
!               as negative integers.  Years in the common era are represented
!               as positive integers.  In particular, note that there is no
!               year 0 in the civil calendar.  1 B.C.E. (-1) is followed by
!               1 C.E. (1).
!
! OUTPUTS:
!	JULDAY returns the Julian Day Number (which begins at noon) of the
!	specified calendar date.
!-----------------------------------------------------------------------
       FUNCTION JULDAY (iyear,month,day)
!
       IMPLICIT NONE
!
       INTEGER*4 JULDAY,year,month,day,iyear
       INTEGER*4 min_calendar,max_calendar,bc,inJanFeb
       INTEGER*4 JY,JM,GREG,JA
!
!
! Gregorian Calander was adopted on Oct. 15, 1582
! skipping from Oct. 4, 1582 to Oct. 15, 1582
       GREG = 2299171  ! incorrect Julian day for Oct. 25, 1582
!
       year=iyear
       min_calendar = -4716
       max_calendar = 5000000
       IF ((year .LT. min_calendar) .OR. (year .GT. max_calendar)) THEN
	   write(6,*)'Value of Julian date is out of allowed range.'
	   stop
       ENDIF
       IF (year .eq. 0)  then
	   write(6,*)'There is no year zero in the civil calendar.'
	   write(6,*)'Value of Julian date is out of allowed range.'
	   stop
       ENDIF
!
       bc=0
       IF (year.lt.0) bc=1
       year=year+bc
!
       inJanFeb=0
       IF (month .le.2) inJanFeb=1

       JY = YEAR - inJanFeb
       JM = MONTH + (1 + 12*inJanFeb)


       JULDAY = INT(365.25d0 * JY) + INT(30.6001d0*JM) + DAY + 1720995


! Test whether to change to Gregorian Calandar.
       IF (JULDAY .GE. GREG) THEN   ! change all dates
          JA = INT(0.01d0 * JY)
	  JULDAY = JULDAY + 2 - JA +INT(0.25d0 * JA)
       ENDIF
       END 
!
!-------------------------------------------------------------------------
! NAME:
!	CALDAT
!
! PURPOSE:
!	Return the calendar date given julian day.
!	This is the inverse of the function JULDAY.
! CALLING SEQUENCE:
!	CALDAT, Julian,  Year, Month, Day
!	See also: julday, the inverse of this function.
!
! INPUTS:
!	JULIAN contains the Julian Day Number (which begins at noon) of the
!	specified calendar date.  It should be a long integer.
! OUTPUTS:
!	(Trailing parameters may be omitted if not required.)
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year.
!-------------------------------------------------------------------------
!
       SUBROUTINE CALDAT(julian, year,month, day)
!
       IMPLICIT NONE
!
       INTEGER*4  julian
       INTEGER*4  year,month, day
       INTEGER*4  min_julian,max_julian,igreg,jalpha,ja,jb,jc,jd,je
!
       min_julian = -1095
       max_julian = 1827933925
!
       IF ((julian .LT. min_julian) .OR. (julian .GT. max_julian)) THEN
          write(6,*)'Value of Julian date is out of allowed range.'
          stop
       ENDIF

       igreg = 2299161    !Beginning of Gregorian calendar

       ja=0
       IF (julian .GE. igreg) THEN    ! all are Gregorian
          jalpha = INT(((julian - 1867216) - 0.25d0) / 36524.25d0)
          ja = julian + 1 + jalpha - INT(0.25d0 * jalpha)
       ENDIF

       jb = ja + 1524
       jc = INT(6680d0 + ((jb-2439870)-122.1d0)/365.25d0)
       jd = INT(365d0 * jc + (0.25d0 * jc))
       je = INT((jb - jd) / 30.6001d0)

       day = jb - jd - INT(30.6001d0 * je)
       month = je - 1
       month = MOD((month - 1),12) + 1
       year = jc - 4715
       if (month .gt. 2) year = year - 1
       if (year .le. 0) year = year - 1

       END
! 
!----------------------------------------------------------------------      
! NAME:
!	GET_DOY
!
! PURPOSE:
!	Calculate the day of year for a given month, day, and year.
!
!
! CALLING SEQUENCE:
!	Result = GET_DOY(Year,Month, Day)
!
! INPUTS:
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year.Year parameters must be valid
!               values from the civil calendar.  Years B.C.E. are represented
!               as negative integers.  Years in the common era are represented
!               as positive integers.  In particular, note that there is no
!               year 0 in the civil calendar.  1 B.C.E. (-1) is followed by
!               1 C.E. (1).
!
! OUTPUTS:
!	GET_DOY returns the day of year of the specified calendar date.
!-----------------------------------------------------------------------
       FUNCTION GET_DOY (year,month,day)
!
       IMPLICIT NONE
!
       INTEGER*4 GET_DOY,JULDAY,year,month,day
       INTEGER*4 firstJanuary
!
       firstJanuary=JULDAY(year,01,01)
       GET_DOY=JULDAY(year,month,day)-firstJanuary+1
       END 
!----------------------------------------------------------------------      
! NAME:
!	DOY_AND_UT2DATE_AND_TIME
!
! PURPOSE:
!	Calculate month, day, and year from year and day of year.
!       Calulate time (hour, minute, second) from UT
!
! CALLING SEQUENCE:
!	CALL DOY_AND_UT2DATE_AND_TIME(Year,Doy,UT, Month, Day, hour, minute, second)
!
! INPUTS:
!	YEAR:	Number of the desired year.Year parameters must be valid
!               values from the civil calendar.  Years B.C.E. are represented
!               as negative integers.  Years in the common era are represented
!               as positive integers.  In particular, note that there is no
!               year 0 in the civil calendar.  1 B.C.E. (-1) is followed by
!               1 C.E. (1).
!
!	DOY:    The day of year of the specified calendar date.
!
!       UT:     Universal time in seconds
!
! OUTPUTS:
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	hour, minute second: time
!-----------------------------------------------------------------------
       SUBROUTINE DOY_AND_UT2DATE_AND_TIME(Year,Doy,UT, Month, Day, 
     &            hour, minute, second)
!
       IMPLICIT NONE
!
       INTEGER*4 jd,year,month,day,doy
       INTEGER*4 hour,minute,second
       INTEGER*4 firstJanuary,julday
       REAL*8    UT
!
       firstJanuary=JULDAY(year,01,01)
       jd=doy+firstJanuary-1
       CALL CALDAT(jd, year,month, day)
       hour=int(UT/3600)
       minute=int((UT-hour*3600)/60)
       second=int(UT-hour*3600.-minute*60.)
       END 
!            
!----------------------------------------------------------------------      
! NAME:
!	DECY2DATE_AND_TIME
!
! PURPOSE:
!	Calculate the date and time (yeay,month,day of month, day of year, hour, minute and second and
!         Universal Time).
!
!
! CALLING SEQUENCE:
!	CALL DECY2DATE_AND_TIME(Dec_y,Year,Month, Day, doy, hour,minute,second)
!
! INPUTS:
!       DEC_Y : Decimal year where yyyy.0d0 is January 1st at 00:00:00
!
! OUTPUTS:
!	YEAR:	Number year.Year parameters must be valid
!               values from the civil calendar.  Years B.C.E. are represented
!               as negative integers.  Years in the common era are represented
!               as positive integers.  In particular, note that there is no
!               year 0 in the civil calendar.  1 B.C.E. (-1) is followed by
!               1 C.E. (1).
!
!	MONTH:	Number  month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	DOY: Number of day of year (DOY=1 is for January 1st)
!
!       HOUR, MINUTE and SECOND: Universal time in the day
!
!      UT: Univeral time in seconds
!
!-----------------------------------------------------------------------
       SUBROUTINE DECY2DATE_AND_TIME (dec_y,year,month, day, doy, 
     &           hour,minute,second,UT)
!
       IMPLICIT NONE
!
       INTEGER*4 I
       INTEGER*4 JULDAY,year,month,day,doy,hour,minute,second
       INTEGER*4 firstJanuary,lastDecember,Ndays
       INTEGER*4 dom_n(12),dom_b(12),dom(12),tmp
!
       REAL*8  dec_y,aux,UT
!
       DATA  dom_n /31,28,31,30,31,30,31,31,30,31,30,31/
       DATA  dom_b /31,29,31,30,31,30,31,31,30,31,30,31/
!       
       year=INT(dec_y)
       firstJanuary=JULDAY(year,01,01)
       lastDecember=JULDAY(year,12,31)
       Ndays=lastDecember-firstJanuary+1
       IF (Ndays.EQ.365) THEN
          DO I=1,12
	     dom(I)=dom_n(I)
	  ENDDO
       ELSE
          DO I=1,12
	     dom(I)=dom_b(I)
	  ENDDO
       ENDIF
       aux=(dec_y-year*1.d0)*Ndays
       doy=INT(aux)+1
!
       tmp=0       
       DO I=1,12
          tmp=tmp+dom(I)
	  IF (tmp .GE. doy) GOTO 10
       ENDDO
10     CONTINUE
       month = I
       tmp=tmp-dom(I)
       day = doy-tmp
!
       aux=(aux-(doy-1)*1.d0)*24.d0
       hour=INT(aux)
       aux=(aux-hour*1.d0)*60.d0
       minute=INT(aux)
       aux=(aux-minute*1.d0)*60.d0
       second=INT(aux)
!
       UT=hour*3600.0d0+minute*60.0d0+second*1.0d0
       END 
       
!-------------------------------------------------------------------------
! NAME:
!	DATE_AND_TIME2DECY
!
! PURPOSE:
!	Return the decimal year for a given date and time.
! CALLING SEQUENCE:
!	CALL DATE_2_DECY(Year, Month, Day, hour,minute,second,decy)
!
! INPUTS:
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR,HOUR,MINUTE,SECOND.
!
! OUTPUTS:
!       decy : decimal year for a given date and time
!-------------------------------------------------------------------------
!
       SUBROUTINE DATE_AND_TIME2DECY(Year,Month,Day,hour,minute
     &,second,decy)
!
       IMPLICIT NONE
!
       INTEGER*4 Year, Month, Day, hour,minute,second
       INTEGER*4 firstJanuary,lastDecember,xut1
       INTEGER*4 julday
!
       REAL*8 decy
!
       firstJanuary=julday(Year,01,01)
       lastDecember=julday(Year,12,31)
       xut1=julday(Year,month,day)
      
       decy=Year*1.d0+(xut1-firstJanuary
     &      +(hour*3600.d0+minute*60.d0+second*1.d0)
     &       /86400.d0)/(lastDecember-firstJanuary+1.d0)

       end

!--------------------------------------------------------------------------

! This function returns TT seconds from J2000.0 for
! a given UTC date and time.
! NAME: QCDFTDB
! INPUT: year,month,day,hour,minute,sec,msec or arrays thereof
!              UTC date and time (from 1/1/1972), UT1 time for 1900-1972
! KEYWORDS: cdfepoin  calculate UTC date from CDF epoch value first
! OUTPUT: TT (~TDB) seconds from J2000.0 (1/1/2000 12:00 TDB)
! RESTRICTIONS: for array input before 2000.0 it works only if timespan
!               less than 6 months. Also does not work when array spans
!               years before and after 2000.
!
! HISTORY: created by m.fraenz@qmw.ac.uk 23/8/2001
!       03/12/2001 array input MF
!       addapted for IRBEM - A. C. Kellerman 01/02/2012
!

      FUNCTION QCDFTDB (year,month,day,hour,minute,sec,msec)
      IMPLICIT NONE
      INTEGER*4 year, month, day, hour, minute
      INTEGER*4 sec, msec, yarr(23), marr(23)
      INTEGER*4 jd, i, index0, yvar

      REAL*8 cdfepoin, deltat, theta, jsec, QCDFTDB

      ! some of this may be repetative from above
      ! calculate Julian Date from J2000 (Seidelmann 12.92):
      yvar=year+4800+(month-14)/12
      jd=(yvar*1461)/4
      jd=jd+((month-2-((month-14)/12)*12)*367)/12
      jd=jd-2451545 ! subtract jd2000.0 to keep numbers small
      yvar=yvar+100
      jd=1D0*(jd-((yvar/100)*3)/4+day-32076)+0.5d0 ! for 00:00

      ! convert to seconds and add seconds of day:
      jsec=jd*86400d0
      jsec=jsec+hour*3600d0
      jsec=jsec+minute*60d0
      jsec=jsec+sec+1D0
      jsec=jsec+msec/1000d0
      ! now add leapseconds:
      IF (year.GE.2000) then
          deltat=63.184d0
      ELSE IF (year.GE.1972) THEN
          deltat=42.184d0
      DATA yarr / 1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,
     *1982,1983,1985,1988,1990,1991,1992,1993,1994,1996,1997,1999,3000/
      DATA marr / 7,1,1,1,1,1,1,1,1,7,7,7,7,1,1,1,7,7,7,1,7,1,0 /

      i=1
100   IF (year.GE.yarr(i)) then
        deltat=deltat+1d0
        i=i+1
        goto 100
      endif
      IF ((i.LT.22).AND.(month.LT.marr(i-1))) deltat=deltat-1d0
      ELSE IF (year.GE.1900) THEN
        theta=(jsec/86400d0/36525d0)+1D0 ! Julian centuries from 1900.0
      deltat=((((((((58353.42d0*theta-232424.66d0)*theta+372919.88d0)
     **theta-303191.19d0)*theta+124906.15d0)*theta-18756.33d0)*
     *theta-2637.8d0)*theta+815.20d0)*theta+87.24)*theta-2.44d0
      ELSE
        deltat=0d0
      END IF
      jsec=jsec+deltat

      QCDFTDB = jsec
      END


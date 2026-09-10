!-----------------------------------------------------------------------------!
!           Group on Data Assimilation Development - GDAD/CPTEC/INPE          !
!-----------------------------------------------------------------------------!
!BOP
!
! !MODULE: m_time.f90
!
! !DESCRIPTION: This module contains routines and functions to manipulate time
!              periods, e.g, functions to calculate total number of hours, days, 
!              months and years between two dates, also contains routines to
!              convert julian days to gregorian day and vice and versa.
!
!
! !INTERFACE:
!
MODULE m_time
  IMPLICIT NONE
  PRIVATE
! !PUBLIC MEMBER FUNCTIONS:

  PUBLIC  :: cal2jul    ! Convert from gregorian to julian day
  PUBLIC  :: jul2cal    ! Convert fron julian to gregorian day
  PUBLIC  :: EndOfMonth ! Calculate end of month
  PUBLIC  :: noh        ! Calculate number of hours between two dates
  PUBLIC  :: nod        ! Calculate number of days between two dates
  PUBLIC  :: nom        ! Calculate number of months between two dates
  PUBLIC  :: noy        ! Calculate number of years between two dates
  PUBLIC  :: doy        ! Calculate the day of year

! !REVISION HISTORY:
! 15 Jun 2005 - J. G. de Mattos - Initial Version
! 18 Mar 2010 - J. G. de Mattos - Include Time calculation:
!                               - End day of Month [eom]
!                               - Number of hours [noh]
!                               - Number of days [nod]
!                               - Number of months [nom]
!                               - Number of year [moy]
! 23 Mar 2010 - J. G. de Mattos - Modified the call for Cal2Jul routine
!                                 was created the interface block for
!                                 use of this new Cal2Jul
! 09 May 2013 - J. G. de Mattos - Removed Bug in Number of Hours 
! 05 Feb 2014 - J. G. de Mattos - Include day of year calculation
! 12 Apr 2021 - J. G. de Mattos - standarize kinds
! 12 Apr 2021 - J. G. de Mattos - Include three ways of jul2cal return
!
! !SEE ALSO: 
! 
!
!EOP
!-----------------------------------------------------------------------------!
!

!
! KINDS
!

   INTEGER, PUBLIC, PARAMETER :: I8  = SELECTED_INT_KIND(14)
   INTEGER, PUBLIC, PARAMETER :: I4  = SELECTED_INT_KIND( 9)

   INTEGER, PUBLIC, PARAMETER :: R8  = SELECTED_REAL_KIND(15)
   INTEGER, PUBLIC, PARAMETER :: R4  = SELECTED_REAL_KIND( 6)

!
! INTERFACES
!
   INTERFACE Cal2Jul
      MODULE PROCEDURE cal2jul_, cal2jul__, cal2jul___
   END INTERFACE Cal2Jul

   INTERFACE Jul2Cal
      MODULE PROCEDURE  jul2cal_, jul2cal__, jul2cal___
   END INTERFACE Jul2Cal

CONTAINS
!
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  endOfMonth
!
! !DESCRIPTION: This function calculate the end day of month.
!
!
! !INTERFACE:
!
  FUNCTION EndOfMonth(year,month) RESULT(day)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN)  :: year
    INTEGER(kind=I4), INTENT(IN)  :: month
!
! !OUTPUT PARAMETERS:
    INTEGER(kind=I4) :: day
!
!
! !REVISION HISTORY: 
!  18 Mar 2010 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    INTEGER(kind=I4), PARAMETER, DIMENSION(12) :: dpm = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    day = dpm(month)
    IF( month.EQ.2 )THEN
      IF ( (MOD(year,4).EQ.0 .AND. MOD(year,100).NE.0).or.(MOD(year,400).EQ.0) ) day = 29
    ENDIF

  END FUNCTION
!
!EOC
!
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  noh
!
! !DESCRIPTION: This function calculate the total number of hours between two
!               dates.
!
!
! !INTERFACE:
!
  FUNCTION noh(di,df) RESULT(Nhour)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN)  :: di ! Starting Date
    INTEGER(kind=I4), INTENT(IN)  :: df ! Ending Date
!
! !OUTPUT PARAMETERS:
!
    INTEGER(I8) :: Nhour
!
!
! !REVISION HISTORY: 
!  18 Mar 2010 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

    Nhour   =  nint(mod(CAL2JUL(df),CAL2JUL(di)) * 24.0)

  END FUNCTION
!
!EOC
!
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  nod
!
! !DESCRIPTION: This function calculate the total number of days between two
!               dates.
!
!
! !INTERFACE:
!
  FUNCTION nod(di,df) RESULT(Nday)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN)  :: di ! Starting Date
    INTEGER(kind=I4), INTENT(IN)  :: df ! Ending Date
!
! !OUTPUT PARAMETERS:
!
    INTEGER(kind=I4) :: Nday
!
!
! !REVISION HISTORY: 
!  18 Mar 2010 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

    Nday   =  ABS(CAL2JUL(df) - CAL2JUL(di)) + 1

  END FUNCTION
!
!EOC
!
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  nom
!
! !DESCRIPTION: This function calculate the total number of months between two
!               dates.
!
!
! !INTERFACE:
!
  FUNCTION nom(di,df) RESULT(Nmonth)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN)  :: di ! Starting Date
    INTEGER(kind=I4), INTENT(IN)  :: df ! Ending Date
!
! !OUTPUT PARAMETERS:
!
    INTEGER(kind=I4) :: Nmonth
!
!
! !REVISION HISTORY: 
!  18 Mar 2010 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    INTEGER(kind=I4)  :: i, f

    i = INT( MOD(di,1000000)/10000 )
    f = INT( MOD(df,1000000)/10000 )

    Nmonth =  ( INT( (df - di) / 1000000 ) * 12 ) + (f-i+1)

  END FUNCTION
!
!EOC
!
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  noy
!
! !DESCRIPTION: This function calculate the total number of Years between two
!               dates.
!
!
! !INTERFACE:
!
  FUNCTION noy(di,df) RESULT(Nyear)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN)  :: di ! Starting Date
    INTEGER(kind=I4), INTENT(IN)  :: df ! Ending Date
!
! !OUTPUT PARAMETERS:
!
    INTEGER(kind=I4) :: Nyear
!
!
! !REVISION HISTORY: 
!  18 Mar 2010 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    Nyear  =   INT( ABS(df - di) / 1000000 )

  END FUNCTION
!
!EOC
!
!***************************************************************
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  doy
!
! !DESCRIPTION: This function calculate the day of the year
!
!
! !INTERFACE:
!
  FUNCTION doy(nymd, nhms) RESULT(Day)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN)  :: nymd ! year month day (yyyymmdd)  
    INTEGER(kind=I4), INTENT(IN)  :: nhms ! hour minute second (hhmnsd)
!
! !OUTPUT PARAMETERS:
!
    REAL(kind=R8) :: Day
!
!
! !REVISION HISTORY: 
!  05 Feb 2014 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

    INTEGER(kind=I4)  :: year
    INTEGER(kind=I4)  :: ymd, hms

    year   = INT ( nymd / 10000 )

    ymd = (year*10000) + 101
    hms = 0

    Day = int(cal2jul__(nymd,hms) - cal2jul__(ymd,hms)) + 1
    Day = Day + (cal2jul__(nymd,nhms) - cal2jul__(nymd,hms) )


  END FUNCTION
!
!EOC
!
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  Cal2Jul__
!
! !DESCRIPTION: This function calculate the julian day from gregorian day
!
!
! !INTERFACE:
!
  FUNCTION cal2jul_(CalDate) RESULT(julian)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN) :: CalDate

!
! !OUTPUT PARAMETERS:
!
    REAL(kind=R8)  :: julian
!
!
! !REVISION HISTORY: 
!  15 Jun 2005 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
  INTEGER(kind=I4) :: nymd
  INTEGER(kind=I4) :: nhms

  nymd = CalDate/100
  nhms = MOD(CalDate,100) * 10000

  julian = cal2jul__(nymd,nhms)

  RETURN

  END FUNCTION
!EOC
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE:  Cal2Jul__
!
! !DESCRIPTION: This function calculate the julian day from gregorian day
!
!
! !INTERFACE:
!
  FUNCTION cal2jul__(ymd,hms) RESULT(julian)
!
! !INPUT PARAMETERS:
!
    INTEGER(kind=I4), INTENT(IN) :: ymd
    INTEGER(kind=I4), INTENT(IN) :: hms
!
! !OUTPUT PARAMETERS:
!
    REAL(kind=R8)  :: julian
!
!
! !REVISION HISTORY: 
!  15 Jun 2005 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    INTEGER(kind=I4)  :: year, month, day
    INTEGER(kind=I4)  :: hour, minute, second
    

    year   = INT ( ymd / 10000 )
    month  = MOD ( ymd,  10000 ) / 100
    day    = MOD ( ymd,    100 )
    hour   = INT ( hms / 10000 )
    minute = MOD ( hms,  10000 ) / 100
    second = MOD ( hms,    100 )

    julian = cal2jul___(year,month,day,hour,minute,second)
    
    RETURN

  END FUNCTION
!
!EOC
!
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE:  Cal2Jul__
!
! !DESCRIPTION: This function calculate the julian day from gregorian day
!
!
! !INTERFACE:
!
!  FUNCTION cal2jul__(ymd,hms) RESULT(julian)
  FUNCTION cal2jul___(year, month, day, hour, minute, second) RESULT(julian)
!
! !INPUT PARAMETERS:
!
   INTEGER(kind=I4) :: year
   INTEGER(kind=I4) :: month
   INTEGER(kind=I4) :: day
   INTEGER(kind=I4) :: hour
   INTEGER(kind=I4) :: minute
   INTEGER(kind=I4) :: second
!
! !OUTPUT PARAMETERS:
!
    REAL(kind=R8)  :: julian
!
!
! !REVISION HISTORY: 
!  15 Jun 2005 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

    REAL(kind=R8)    :: A, B, C, D, E
    INTEGER(kind=I8) :: ymd
    

    IF(month < 3)THEN
       year=year-1
       month=month+12
    ENDIF

    ymd = (year*10000)+(month*100)+day
    IF(ymd>=15821015)THEN
       A = INT(year/100)
       B = INT(A/4)
       C = 2 - A + B
    ENDIF

    IF(ymd<=15821004)THEN
       C = 0
    ENDIF

    D = INT(365.25 * (year + 4716))
    E = INT(30.6001 * (month + 1))

    julian = INT(D + E + day + 0.5 + C - 1524.5) +      &
                ( hour / 24.0 ) + ( minute / (60*24) ) + &
                ( second / (60*60*24) ) 

    RETURN

  END FUNCTION
!
!EOC
!
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  Jul2Cal
!
! !DESCRIPTION: This function calculate the gregorian date from julian day.
!
!
! !INTERFACE:
!
  SUBROUTINE jul2cal_(jd,year, month, day, hour, minute, second)
!
! !INPUT PARAMETERS:
!
	 REAL (kind=R8), INTENT(IN) :: jd
!
! !OUTPUT PARAMETERS:
!
   INTEGER (kind=4)            ::  year
   INTEGER (kind=4)            ::  month
   INTEGER (kind=4)            ::  day
   INTEGER (kind=4)            ::  Hour
   INTEGER (kind=4)            ::  Minute
   INTEGER (kind=4)            ::  Second
!
!
! !REVISION HISTORY: 
!  15 Jun 2005 - J. G. de Mattos - Initial Version
!  23 Mar 2011 - J. G. de Mattos - Modified Interface
!                                  to a subroutine call
!
! !REMARKS:
!        This algorithm was adopted from Press et al.
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    INTEGER (kind=I4), PARAMETER :: Gregjd = 2299161
    INTEGER (kind=I4)            ::  j1, j2, j3, j4, j5
    INTEGER (kind=I4)            ::  Intgr
    INTEGER (kind=I4)            ::  f, tmp

    REAL (kind=R8)               ::  dayfrac, frac

    !       
    ! get the date from the Julian day number
    !       
    ! jd=2453372.25

    intgr   = INT(jd)
    frac    = real(jd - intgr,8)

    IF( intgr >= gregjd )THEN              !Gregorian calendar correction
       tmp = INT( ( (intgr - 1867216) - 0.25 ) / 36524.25 )
       j1  = intgr + 1 + tmp - INT(0.25*tmp)
    ELSE
       j1 = intgr
    ENDIF
    !       correction for half day offset

    dayfrac = frac + 0.0d0
    IF( dayfrac >= 1.0 )THEN
       dayfrac = dayfrac - 1.0d0
       j1 = j1+1
    ENDIF

    j2 = j1 + 1524
    j3 = INT( 6680.0 + ( (j2 - 2439870) - 122.1 )/365.25 )
    j4 = INT(j3*365.25)
    j5 = INT( (j2 - j4)/30.6001 )

    day   = INT(j2 - j4 - INT(j5*30.6001))
    month = INT(j5 - 1)
    IF( month > 12 ) month = month- 12
    year = INT(j3 - 4715)
    IF( month > 2 )  year = year - 1
    IF( year <= 0 ) year = year - 1

    !
    ! get time of day from day fraction
    !
    hour   = INT(dayfrac * 24.0d0)
    minute = INT((dayfrac*24.0d0 - hour)*60.0d0)
    f      = NINT( ((dayfrac*24.0d0 - hour)*60.0d0 - minute)*60.0d0)
    second = INT(f)
    f      = f-second
    IF( f > 0.5 ) second = second + 1

    IF( second == 60 )THEN
       second = 0
       minute = minute + 1
    ENDIF

    IF( minute == 60 )THEN
       minute = 0
       hour   = hour + 1
    ENDIF

    IF( hour == 24 )THEN

       hour = 0
       !
       ! this could cause a bug, but probably 
       ! will never happen in practice
       !
       day = day + 1

    ENDIF

    IF( year < 0 )THEN
       year = year * (-1)
    ENDIF

  END SUBROUTINE
!BOC
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  Jul2Cal
!
! !DESCRIPTION: This function calculate the gregorian date from julian day.
!
!
! !INTERFACE:
!
  SUBROUTINE jul2cal__(jd, gregorian)
!
! !INPUT PARAMETERS:
!
	 REAL (kind=R8), INTENT(IN) :: jd

!
! !OUTPUT PARAMETERS:
!
   INTEGER(kind=I4), INTENT(OUT) :: gregorian
!
!
! !REVISION HISTORY: 
!  23 Mar 2011 - J. G. de Mattos - created to maintain the 
!                                  original interface
!                                  of jul2cal function
! 
! !REMARKS:
!        This algorithm was adopted from Press et al.
!EOP
!-----------------------------------------------------------------------------!
!BOC

  INTEGER (kind=I4)            ::  year
  INTEGER (kind=I4)            ::  month
  INTEGER (kind=I4)            ::  day
  INTEGER (kind=I4)            ::  Hour
  INTEGER (kind=I4)            ::  Minute
  INTEGER (kind=I4)            ::  Second

  call jul2cal_(jd,year, month, day, hour, minute, second)

  gregorian=(year*1000000)+(month*10000)+(day*100)+hour

  END SUBROUTINE
!
!EOC
!-----------------------------------------------------------------------------!
!BOP
!
! !IROUTINE:  Jul2Cal
!
! !DESCRIPTION: This function calculate the gregorian date from julian day.
!
!
! !INTERFACE:
!
  SUBROUTINE jul2cal___(jd,ymd,hms)
!
! !INPUT PARAMETERS:
!
	 REAL (kind=R8), INTENT(IN) :: jd

!
! !OUTPUT PARAMETERS:
!
   INTEGER(kind=R4), INTENT(OUT) :: ymd ! year month day (yyyymmdd)
   INTEGER(kind=R4), INTENT(OUT) :: hms ! hour minute second (hhmnsd)
!
!
! !REVISION HISTORY: 
!  23 Mar 2011 - J. G. de Mattos - created to maintain the 
!                                  original interface
!                                  of jul2cal function
! 
! !REMARKS:
!        This algorithm was adopted from Press et al.
!EOP
!-----------------------------------------------------------------------------!
!BOC

   INTEGER (kind=4)            ::  year
   INTEGER (kind=4)            ::  month
   INTEGER (kind=4)            ::  day
   INTEGER (kind=4)            ::  Hour
   INTEGER (kind=4)            ::  Minute
   INTEGER (kind=4)            ::  Second

   call jul2cal_(jd, year, month, day, hour, minute, second)

   ymd = (year*10000)+(month*100)+(day)
   hms = (hour*10000)+(minute*100)+second

  END SUBROUTINE
!
!EOC
!-----------------------------------------------------------------------------!


END MODULE m_time

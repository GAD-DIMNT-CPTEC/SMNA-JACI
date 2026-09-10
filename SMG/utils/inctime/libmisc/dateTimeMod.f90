module dateTimeMod
  implicit none
  private

  public :: dateTime
  public :: today
  public :: timeDelta
  public :: isLeapYear
  public :: daysInYear
  public :: daysInMonth
  public :: deltaTime
  public :: gregorian
  public :: caldat
  !   public :: numberOfDays

  integer,  parameter :: i8  = selected_int_kind(14)
  integer,  parameter :: i4  = selected_int_kind( 9)
  integer,  parameter :: r8  = selected_real_kind(15)
  integer,  parameter :: r4  = selected_real_kind( 6)

  type timeDelta
     integer( kind = i4 ) :: years   = 0
     integer( kind = i4 ) :: months  = 0
     integer( kind = i4 ) :: days    = 0
     integer( kind = i4 ) :: hours   = 0
     integer( kind = i4 ) :: minutes = 0
     integer( kind = i4 ) :: seconds = 0
  end type timeDelta

  type dateTime
     integer( kind = i4 ) :: year     = 1900
     integer( kind = i4 ) :: month    = 1
     integer( kind = i4 ) :: day      = 1
     integer( kind = i4 ) :: hour     = 0
     integer( kind = i4 ) :: minute   = 0
     integer( kind = i4 ) :: second   = 0
     integer( kind = i4 ) :: timeZone = 0
   contains
     procedure, public  :: jdn => cal2jul
     procedure, public  :: jdn2 => julday
     procedure, public  :: jdn3 => julday2
     procedure, private :: timeAdd
     generic :: operator(+) => timeAdd
     procedure, private :: timeSub, dateDiff
     generic :: operator(-) => timeSub, dateDiff
     procedure, private :: equal
     generic :: operator(.eq.) => equal
     procedure, private :: notEqual
     generic :: operator(.ne.) => notEqual
     procedure, private :: greaterThan
     generic :: operator(.gt.) => greaterThan
     procedure, private :: greaterEqualThan
     generic :: operator(.ge.) => greaterEqualThan
     procedure, private :: lessThan
     generic :: operator(.lt.) => lessThan
     procedure, private :: lessEqualThan
     generic :: operator(.le.) => lessEqualThan
  end type dateTime



contains

  type(dateTime) function today()

    integer( kind = i4 ) :: date(9)

    call DATE_AND_TIME(Values=date)

    today%year     = date(1)
    today%month    = date(2)
    today%day      = date(3)
    today%hour     = date(5)
    today%minute   = date(6)
    today%second   = date(7)
    today%timeZone = date(4)/60.0

    return
  end function today

  pure elemental function isLeapYear(year) result(flag)
    integer, intent(in) :: year
    logical :: flag

    flag = .false.
    if (mod(year,  4) .eq. 0 .and. &
         mod(year,100) .ne. 0 .or.  &
         mod(year,400) .eq. 0) flag = .true.

  end function isLeapYear

  pure elemental function daysInYear(year) result (days)
    integer, intent(in) :: year
    integer :: days

    if (isLeapYear(year))then
       days = 366
    else
       days = 365
    endif
  end function daysInYear

  pure elemental function daysInMonth(month, year) result(days)
    integer, intent(in   ) :: year
    integer, intent(in   ) :: month
    integer                :: days

    integer, parameter :: monthDays(*) = [31, 28, 31, 30, 31, 30, &
         31, 31, 30, 31, 30, 31]

    select case (month)
    case(:0,13:)
       days = 0
    case(2)
       if (isLeapYear(year))then
          days = 29
       else
          days = monthDays(month)
       endif
    case default
       days = monthDays(month)
    end select
    return
  end function daysInMonth

  function deltaTime(days) result(delta)
    real(kind=r8), intent(in) :: days
    type(timeDelta)     :: delta

    delta%years  = int(days/365)
    delta%months = int(mod(int(days),365)/30)
    delta%days   = mod(mod(int(days),365),30)

    return
  end function deltaTime

  function numberOfDays(dt) result(days)
    class(timeDelta), intent(in) :: dt
    integer :: days

    days = dt%years*365 + dt%months*30 + dt%days

  end function numberOfDays

  function timeAdd(v1,v2) result(v3)
    class(dateTime),  intent(in) :: v1 ! date
    class(timeDelta), intent(in) :: v2 ! deltaTime
    type(dateTime) :: v3


    real(kind = r8)     :: incr
    real(kind = r8)     :: julianDay
    integer(kind = i4 ) :: value

    incr  = (float(v2%hours)/24.0) 
    incr  = incr + (float(v2%minutes)/(60.0*24.0))
    incr  = incr + (float(v2%seconds)/(60.0*60.0*24.0))

    julianDay = cal2jul(v1) + numberOfDays(v2) + incr

    call jul2cal(julianDay,v3)

    return
  end function timeAdd

  function timeSub(v1,v2) result(v3)
    class(dateTime),  intent(in) :: v1 ! date
    class(timeDelta), intent(in) :: v2 ! deltaTime
    type(dateTime) :: v3


    real(kind = r8)     :: incr
    real(kind = r8)     :: julianDay
    integer(kind = i4 ) :: value

    incr  = (float(v2%hours)/24.0) 
    incr  = incr + (float(v2%minutes)/(60.0*24.0))
    incr  = incr + (float(v2%seconds)/(60.0*60.0*24.0))

    julianDay = cal2jul(v1) - numberOfDays(v2) - incr

    call jul2cal(julianDay,v3)

    return
  end function timeSub


  function dateDiff(d1,d2) result(deltaT)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    type(timeDelta) :: deltaT

    real(kind=r8) :: days
    real(kind=r8) :: incr

    days = cal2jul(d1) - cal2jul(d2)
    
    deltaT = deltaTime(days)

    return
  end function dateDiff




  !   function timeDiff(d1,d2) result(date)
  !      class(dateTime),  intent(in) :: d1
  !      class(timeDelta), intent(in) :: d2
  !      type(dateTime) :: date
  !
  !      integer :: days
  !
  !      days = cal2jul(d1) - numberOfDays(d2)
  !
  !      date%year  = int(days/365)
  !      date%month = int(mod(int(days),365)/30)
  !      date%day   = mod(mod(int(days),365),30)
  !      
  !      return
  !   end function

  function equal(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .eq. cal2jul(d2))

    return
  end function equal

  function notequal(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .ne. cal2jul(d2))

    return
  end function notequal

  function greaterThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .gt. cal2jul(d2))

    return
  end function greaterThan

  function greaterEqualThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .ge. cal2jul(d2))

    return
  end function greaterEqualThan

  function lessThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .lt. cal2jul(d2))

    return
  end function lessThan

  function lessEqualThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .le. cal2jul(d2))

    return
  end function lessEqualThan



  function julday(d) result(jdn)
    class(dateTime), intent(in) :: d
    real(kind=r8) :: jdn

    real :: year
    real :: month

    year  = d%year
    month = d%month
    if (month < 3)then
       month = month + 12
       year  = year - 1
    endif

    jdn = d%day + int((153.0*month-457.0)/5.0) + &
         365*year + floor(year/4.0) - floor(year/100.0) + &
         floor(year/400.0) + 1721118.5

  end function julday

  function gregorian(jd)result(dt)
    real(kind=r8), intent(in) :: jd
    type(dateTime) :: dt

    real(kind=r8) :: Z, R, G, A, B, C

    Z = floor(jd - 1721118.5)
    R = jd - 1721118.5 - Z
    G = Z - .25
    A = floor(G / 36524.25)
    B = A - floor(A / 4)
    dt%year = floor((B+G) / 365.25)
    C = B + Z - floor(365.25 * dt%year)
    dt%month = int((5 * C + 456) / 153)
    dt%day   = C - int((153 * dt%month - 457) / 5) + R

    if (dt%month > 12)then
       dt%year  = dt%year + 1
       dt%month = dt%month - 12
    endif

    return
  end function gregorian



  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: lower_case - convert uppercase letters to lowercase.
  !
  ! !DESCRIPTION:
  !
  ! !INTERFACE:

  function lower_case(str) result(lstr)
    implicit none
    character(len=*), intent(in) :: str
    character(len=len(str))      :: lstr

    ! !REVISION HISTORY:
    !       13Aug96 - J. Guo        - (to do)
    !EOP
    !_______________________________________________________________________
    integer i
    integer,parameter :: iu2l=ichar('a')-ichar('A')

    lstr=str
    do i=1,len_trim(str)
       if(str(i:i).ge.'A'.and.str(i:i).le.'Z')   &
            lstr(i:i)=char(ichar(str(i:i))+iu2l)
    end do
  end function lower_case


  !-----------------------------------------------------------------------------!
  !
  !BOP
  !
  ! !IROUTINE:  cal2Jul
  !
  ! !DESCRIPTION: This function calculate the julian day from gregorian day
  !
  !
  ! !INTERFACE:
  !
  FUNCTION cal2jul(dt) RESULT(jd)
    !
    ! !INPUT PARAMETERS:
    !
    class(dateTime), intent(in) :: dt
    !
    ! !OUTPUT PARAMETERS:
    !
    REAL(kind=R8)  :: jd
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
    integer(kind=i4) :: year
    integer(kind=i4) :: month
    integer(kind=i4) :: day

    real(kind=R8)    :: A, B, C, E, F
    real(kind=r8)    :: dayFrac
    integer(kind=I8) :: ymd

    year  = dt%year
    month = dt%month
    day   = dt%day

    dayFrac = (dt%minute/60.0 + dt%hour)/24.0

    !-----------------------------------------------
    ! January and February are 13th and 14th months
    ! of the previous year
    !
    if(month .eq. 1 .or. month .eq. 2)then
       year  = year  -  1
       month = month + 12
    endif
    !-----------------------------------------------

    ymd = (year*10000)+(month*100)+day
    if(ymd>=15821015)then
       A = INT(year/100)
       B = INT(A/4)
       C = 2 - A + B
    endif

    if(ymd<=15821004)then
       c = 0
    endif

    E = INT(365.25 * (year + 4716))
    F = INT(30.6001 * (month + 1))

    jd = C + day + E + F - 1524.5 + 0.5

    jd = jd + (dt%hour/24.0)
    jd = jd + (dt%minute/(60.0*24.0))
    jd = jd + (dt%second/(60.0*60.0*24.0))

    RETURN

  END FUNCTION cal2jul
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
  SUBROUTINE jul2cal(jd, dt)
    !
    ! !INPUT PARAMETERS:
    !
    real (kind=r8), intent(in   ) :: jd
    !
    ! !OUTPUT PARAMETERS:
    !
    type(dateTime), intent(  out) :: dt
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
    integer (kind=i4), parameter :: gregjd = 2299161
    integer (kind=i4)            :: year
    integer (kind=i4)            :: month
    integer (kind=i4)            :: day
    integer (kind=i4)            :: hour
    integer (kind=i4)            :: minute
    integer (kind=i4)            :: second
    integer (kind=i4)            :: j1, j2, j3, j4, j5
    integer (kind=i4)            :: intgr
    integer (kind=i4)            :: f, tmp

    real (kind=r8)               :: dayfrac, frac

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

    dt%year   = year
    dt%month  = month
    dt%day    = day
    dt%hour   = hour
    dt%minute = minute
    dt%second = second

  END SUBROUTINE jul2cal
  !
  !EOC
  !
  !-----------------------------------------------------------------------------!

  !   function add (self, str) result(status)
  !      class(dateTime)  :: self
  !      character(len=*) :: str
  !      integer          :: status
  !      integer          :: value
  !
  !      integer :: iyear
  !      integer :: imonth
  !      integer :: iday
  !      integer :: ihour
  !      integer :: iminute
  !      integer :: isecond
  !
  !      iyear   = 0
  !      imonth  = 0
  !      iday    = 0
  !      ihour   = 0
  !      iminute = 0
  !      isecond = 0
  !
  !      j=1
  !      do i=1,len_trim(str)
  !         select case(str(i:i))
  !         case('y')
  !             read(str(j:i-1),*)value
  !             j = i + 1
  !             iyear = iyear + value
  !         case('m')
  !             read(str(j:i-1),*)value
  !             j = i + 1
  !             iyear  = iyear + int(float(value)/12.0)
  !             imonth = imonth + mod(value,12)
  !         case('d')
  !             read(str(j:i-1),*)value
  !             j = i + 1
  !             iyear = iyear + int(float(value)/365.0)
  !             iday  = iday + mod(value,365)
  !         case('h')
  !             read(str(j:i-1),*)value
  !             j = i + 1
  !             iday  = iday + int(float(value)/24.0)
  !             ihour = ihour + mod(value,24)
  !         case('n')
  !             read(str(j:i-1),*)value
  !             j = i + 1
  !             ihour   = ihour + int(float(value)/60.0)
  !             iminute = iminute+ mod(value,60)
  !         case('s')
  !             read(str(j:i-1),*)value
  !             j = i + 1
  !             iminute = iminute + int(float(value)/60.0)
  !             isecond = isecond + mod(value,60)
  !         end select   
  !      enddo
  !   
  !      value = month + imonth
  !      self%year  = self%year + iyear + int(float(value)/12.0)
  !      self%month = mod(value,12)
  !      
  !      incr  = iday 
  !      incr  = incr + (float(ihour)/24.0) 
  !      incr  = incr + (float(iminute)/(60.0*24.0))
  !      incr  = incr + (float(isecond)/(60.0*60.0*24.0))
  !
  !      julianDay = cal2jul(self%year,self%month,self%day,self%hour,self%minute,self%second)
  !      julianDay = julianDay + incr
  !
  !      call jul2cal(julianDay,self%year,self%month,self%day,self%hour,self%minute,self%second)
  !
  !      return
  !      
  !   end function

  FUNCTION julday2(dt) result(jd)
    class(dateTime), intent(in) :: dt
    REAL(R8) :: jd
    INTEGER(I4), PARAMETER :: IGREG=15+31*(10+12*1582)
    INTEGER(I4) :: ja,jm,jy
    jy=dt%year
    if (jy == 0) print*,'julday: there is no year zero'
    if (jy < 0) jy=jy+1
    if (dt%month > 2) then
       jm=dt%month+1
    else
       jy=jy-1
       jm=dt%month+13
    end if
    jd=int(365.25*jy)+int(30.6001*jm)+dt%day+1720995
    if (dt%day+31*(dt%month+12*dt%year) >= IGREG) then
       ja=int(0.01*jy)
       jd=jd+2-ja+int(0.25*ja)
    end if
  END FUNCTION julday2


  FUNCTION caldat(julian) result(dt)
    REAL(R8), INTENT(IN) :: julian
    type(dateTime) :: dt
    INTEGER(I4) :: mm,id,iyyy
    !       INTEGER(I4B), INTENT(OUT) :: mm,id,iyyy
    INTEGER(I4) :: ja,jalpha,jb,jc,jd,je
    INTEGER(I4), PARAMETER :: IGREG=2299161
    if (julian >= IGREG) then
       jalpha=int(((julian-1867216)-0.25_r4)/36524.25_r4)
       ja=julian+1+jalpha-int(0.25_r4*jalpha)
    else
       ja=julian
    end if
    jb=ja+1524
    jc=int(6680.0_r4+((jb-2439870)-122.1_r4)/365.25_r4)
    jd=365*jc+int(0.25_r4*jc)
    je=int((jb-jd)/30.6001_r4)
    id=jb-jd-int(30.6001_r4*je)
    mm=je-1
    if (mm > 12) mm=mm-12
    iyyy=jc-4715
    if (mm > 2) iyyy=iyyy-1
    if (iyyy <= 0) iyyy=iyyy-1
    dt%year = iyyy
    dt%month = mm
    dt%day = id
  END FUNCTION caldat
end module dateTimeMod

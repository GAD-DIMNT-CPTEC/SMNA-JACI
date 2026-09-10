PROGRAM inctime
! Main program
!
! Package Overview
!   Inctime is a set of routines in fortran 90 to perform 
!   calculations with dates; The main function is to calculate 
!   dates to past or to future from any date. Other operations can 
!   be performed by using the modules available in this routine. 
!   Among other duties, they are calculating the Julian day, number 
!   of years, months, days and hours between any two dates.
!
! Autor
!   João Gerd Zell de Mattos
!
! Affiliation
!   Group on Data Assimilation Development - CPTEC/INPE
!
! Date
!   March 22, 2011
! 
! Usage
!   The way to use this program is through the command line as follows:
!
!   inctime [yyyymmddhh, yyyymmdd] [<+,->nynmndnhnnns] [Form. output]
!
! where
!
!    * [yyyymmddhh, yyyymmdd] 
!                             - Initial Time			
!    * [<+,->nynmndnhnnns]    
!                             - ( -) calculate the passed date
!                             - ( +) calculate the future date (default)		
!                             - (ny) Number of year (default is 0)
!                             - (nm) Number of months (default is 0)
!                             - (nd) Number of days (default is 0)
!                             - (nh) Number of hours (default is 0)
!                             - (nn) Number of minutes (default is 0)
!                             - (ns) Number of seconds (default is 0)
!
!    * [ Form. Output ]       
!                             - Format to output date. is a template Format
!                               The format descriptors are similar to those
!                               used in the GrADS.
!
!                          * "%y4"  substitute with a 4 digit year
!                          * "%y2"  a 2 digit year
!                          * "%m1"  a 1 or 2 digit month
!                          * "%m2"  a 2 digit month
!                          * "%mc"  a 3 letter month in lower cases
!                          * "%Mc"  a 3 letter month with a leading letter in upper case
!                          * "%MC"  a 3 letter month in upper cases
!                          * "%d1"  a 1 or 2 digit day
!                          * "%d2"  a 2 digit day
!                          * "%h1"  a 1 or 2 digit hour
!                          * "%h2"  a 2 digit hour
!                          * "%h3"  a 3 digit hour (?)
!                          * "%n2"  a 2 digit minute
!                          * "%e"   a string ensemble identify
!                          * "%jd"  a julian day without hours decimals
!                          * "%jdh" a julian day with hour decimals
!                          * "%jy"  a day of current year without hours decimals
!                          * "%jyh" a day of current year with hours decimals
!
!   Can use words to compose the output format.
!
! Examples
!   * inctime 2001091000 +1d %d2/%m2/%y4
!   * inctime 2001091000 +48h30n %h2Z%d2%MC%y4
!   * inctime 2001091000 -1h30n 3B42RT.%y4%m2%d2%h2.bin
!   * inctime 2001091000 -2h45n 3B42RT.%y4%m2%d2%h2.bin
!   * inctime 2001091000 -1y3m2d1h45n ANYTHING.%y4%m2%d2%h2.ANYTHING
!   
! History
!   * 22 Mar 2011 - Joao Gerd - Initial Code
!   * 27 Jul 2011 - Joao Gerd - Bug into the end print
!   * 19 Jul 2012 - Joao Gerd - Add option to month increment/decrement
!   * 15 Apr 2013 - Joao Gerd - correct bug to incr/decr month
!   * 03 May 2013 - Joao Gerd - correct bug to incr/decr all times (ym was not 
!                               being properly initialized)
!   * 05 Feb 2014 - Joao Gerd - Add julian day
!   * 20 Jul 2014 - Joao Gerd - upgrade to new version of m_string.f90
!                             - update help banner and documentation

   USE m_time
   USE m_string
   USE m_stdio


  character(len=7)  :: myname_='inctime'

  character(len=50) :: time_incr
  character(len=14) :: time_strg
  integer           :: nymd, y, d, ym
  integer           :: nhms, h, m, s
  integer           :: incr_i, incr_y

  integer           :: iyear
  integer           :: imonth 
  integer           :: iday
  integer           :: ihour
  integer           :: iminute
  integer           :: isecond

  integer           :: year  = 1900
  integer           :: month =   01
  integer           :: day   =   15
  integer           :: hour  =   00
  integer           :: minute=   00
  integer           :: second=   00

  integer           :: ytmp
  real(kind=8)      :: incr
  real(kind=8)      :: jday
  real(kind=8)      :: jdoy
  logical           :: increment
  character(len=300):: mask,tmp,bufr
!
! Get Options 
!
   call getarg(1,bufr)
   SELECT CASE(bufr)
      CASE ('--help':'-h')
         Call usage()
         call exit(0)
      CASE DEFAULT
         IF(iargc().ne.3)THEN
            write(stderr,*)'Use: inctime [yyyymmddhh,yyyymmdd] [<+,->ndnhnmns] [output format]'
            write(stderr,*)'Try `inctime --help` for more informations.' 
            Call exit(1)
         ENDIF
   END SELECT

   !
   ! Pegando data
   !
   call getarg(1,time_strg) 
   read(time_strg,'(I4,5I2)')year,month,day,hour,minute,second

   call getarg(2,time_incr)
   call getarg(3,mask)
!
!
!
   incr      = 0.0
   iyear     = 0
   imonth    = 0
   increment = .true.
   J=1
   DO I=1,len_trim(time_incr)
      SELECT CASE(time_incr(I:I))
         case( '+' ) ! increment
            increment=.true.
            J=I+1
         case( '-' ) ! decrement
            increment=.false. 
            J=I+1
         case( 'y' ) ! Year
           read(time_incr(J:I-1),*)iyear
           J=I+1
         case( 'm' ) ! month
           read(time_incr(J:I-1),*)imonth
           J=I+1
           iyear  = iyear + int(float(imonth)/12)
           imonth = imonth - int(float(imonth)/12) * 12
         case( 'd' ) ! day
           read(time_incr(J:I-1),*)iday
           J=I+1
           incr = incr + real(float(iday),8)
         case( 'h' ) ! hour
           read(time_incr(J:I-1),*)ihour
           J=I+1
           incr = incr + real(float(ihour)/24.0,8)
         case( 'n' ) ! minute
           read(time_incr(J:I-1),*)iminute
           J=I+1
           incr = incr + real(float(iminute)/(60.0*24.0),8)
         case( 's' ) ! second
           read(time_incr(J:I-1),*)isecond
           J=I+1
           incr = incr + real(float(isecond)/(60.0*60.0*24.0),8)
      END SELECT
   ENDDO

!
! Increment/Decrement time
!

   IF (increment)THEN
      tmp  = mask

      if((imonth + month) .gt. 12)then
         iyear  = iyear + 1
         imonth = imonth - 12
      endif

      nymd =  (year + iyear)*10000 + (month + imonth)*100 + day
      nhms =  hour*10000 + minute*100 + second

      jday = cal2jul(nymd,nhms)+incr

      call jul2cal(jday,nymd,nhms)!,incr

      jdoy = doy(nymd,nhms)

      call str_template(tmp,nymd,nhms,jd=jday,doy=jdoy)
      WRITE(stdout,'(A)')TRIM(tmp)
   ELSE
      tmp  = mask

      if((month - imonth) .lt. 1)then
         iyear  = iyear + 1
         imonth = 12 - imonth
      endif

      nymd =  (year - iyear)*10000 + (month - imonth)*100 + day
      nhms =  hour*10000 + minute*100 + second

      jday = cal2jul(nymd,nhms)-incr

      call jul2cal(jday,nymd,nhms)!,incr

      jdoy = doy(nymd,nhms)

      call str_template(tmp,nymd,nhms,jd=jday,doy=jdoy)
      WRITE(stdout,'(A)')TRIM(tmp)
   ENDIF

END PROGRAM

SUBROUTINE usage()
   USE m_stdio

      write(stdout,*) "+------------------------------------------------------------------+"
      write(stdout,*) "|                          Inctime V0.01                           |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| Author: Joao Gerd - joao.gerd@cptec.inpe.br - March 2011         |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| Inctime calculate dates to past or to future from any date.      |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| To use it, type:                                                 |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| inctime [yyyymmddhh, yyyymmdd] [<+,->nynmndnhnnns] [Output Form.]|"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| Common Arguments:                                                |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| [yyyymmddhh, yyyymmdd]: Initial Time                             |"
      write(stdout,*) "| [<+,->nynmndnhnnns]   : ( -) calculate the passed date           |"
      write(stdout,*) "|                       : ( +) calculate the future date(default)  |"
      write(stdout,*) "|                       : (ny) Number of year (default is 0)       |"
      write(stdout,*) "|                       : (nd) Number of days (default is 0)       |"
      write(stdout,*) "|                       : (nm) Number of months (default is 0)     |"
      write(stdout,*) "|                       : (nh) Number of hours (default is 0)      |"
      write(stdout,*) "|                       : (nn) Number of minutes (default is 0)    |"
      write(stdout,*) "|                       : (ns) Number of seconds (default is 0)    |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| [Output Form.]        : Format to output date. is a template     |"
      write(stdout,*) "|                         Format. The format descriptors are       |"
      write(stdout,*) "|                         similar to those used in the GrADS.      |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "|                         '%y4'  substitute with a 4 digit year    |"
      write(stdout,*) "|                         '%y2'  a 2 digit year                    |"
      write(stdout,*) "|                         '%m1'  a 1 or 2 digit month              |"
      write(stdout,*) "|                         '%m2'  a 2 digit month                   |"
      write(stdout,*) "|                         '%mc'  a 3 letter month in lower cases   |"
      write(stdout,*) "|                         '%Mc'  a 3 letter month with a leading   |"
      write(stdout,*) "|                                    letter in upper case          |"
      write(stdout,*) "|                         '%MC'  a 3 letter month in upper cases   |"
      write(stdout,*) "|                         '%d1'  a 1 or 2 digit day                |"
      write(stdout,*) "|                         '%d2'  a 2 digit day                     |"
      write(stdout,*) "|                         '%h1'  a 1 or 2 digit hour               |"
      write(stdout,*) "|                         '%h2'  a 2 digit hour                    |"
      write(stdout,*) "|                         '%h3'  a 3 digit hour (?)                |"
      write(stdout,*) "|                         '%n2'  a 2 digit minute                  |"
      write(stdout,*) "|                         '%e'   a string ensemble identify        |"
      write(stdout,*) "|                         '%jd'  a julian day without hours decimal|"
      write(stdout,*) "|                         '%jdh' a julian day with hour decimals   |"
      write(stdout,*) "|                         '%jy'  a day of current year without     |"
      write(stdout,*) "|                                  hours decimals                  |"
      write(stdout,*) "|                         '%jyh' a day of current year with hours  |"
      write(stdout,*) "|                                  decimals                        |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| You can use words to compose the output format.                  |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "| Examples:                                                        |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "|       inctime 2001091000 +1d %d2/%m2/%y4                         |"
      write(stdout,*) "|       inctime 2001091000 +48h30n %h2Z%d2%MC%y4                   |"
      write(stdout,*) "|       inctime 2001091000 -1h30n 3B42RT.%y4%m2%d2%h2.bin          |"
      write(stdout,*) "|       inctime 2001091000 -2h45n 3B42RT.%y4%m2%d2%h2.bin          |"
      write(stdout,*) "|       inctime 2001091000 -1y3m2d1h45n ANYTHING.%y4%m2%d2%h2      |"
      write(stdout,*) "|                                                                  |"
      write(stdout,*) "+------------------------------------------------------------------+"

END SUBROUTINE

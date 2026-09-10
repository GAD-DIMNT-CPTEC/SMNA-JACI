!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE InputParameters

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: &
            r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
   INTEGER, PARAMETER, PUBLIC :: &
            r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers

   INTEGER, PUBLIC :: Imax, Jmax

   INTEGER, PUBLIC :: MonthBefore, MonthAfter

   REAL (KIND=r8), PUBLIC :: FactorA, FactorB

   CHARACTER (LEN=10), PUBLIC :: Date

   LOGICAL, PUBLIC :: GrADS

   REAL (KIND=r8), PUBLIC :: Undef=-999.0

   CHARACTER (LEN=5), PUBLIC :: Exts='S.unf'

   CHARACTER (LEN=9), PUBLIC :: VarName='SnowClima'

   CHARACTER (LEN=12),PUBLIC :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=13),PUBLIC :: DirModelIn='model/datain/'

   CHARACTER (LEN=16),PUBLIC :: NameLSM='ModelLandSeaMask'

   CHARACTER (LEN=6),PUBLIC :: NameAlb='Albedo'

   CHARACTER (LEN=4),PUBLIC :: VarNameS='Snow'

   CHARACTER (LEN=7),PUBLIC :: nLats='.G     '

   CHARACTER (LEN=10),PUBLIC :: mskfmt = '(      I1)'

   CHARACTER (LEN=12),PUBLIC :: TimeGrADS='  Z         '

   LOGICAL, PUBLIC :: IcePoints=.FALSE.

   CHARACTER (LEN=528), PUBLIC :: DirMain

   INTEGER, PUBLIC :: nferr=0    ! Standard Error Print Out
   INTEGER, PUBLIC :: nfinp=5    ! Standard Read In
   INTEGER, PUBLIC :: nfprt=6    ! Standard Print Out
   INTEGER, PUBLIC :: nflsm=10   ! To Read Model Land Sea Mask Data
   INTEGER, PUBLIC :: nfclm=20   ! To Read Climatological Albedo Data
   INTEGER, PUBLIC :: nfout=30   ! To Write Intepolated Climatological Snow Data
   INTEGER, PUBLIC :: nfctl=40   ! To Write Output Data Description

   PUBLIC :: InitInputParameters


CONTAINS


SUBROUTINE InitInputParameters ()

   IMPLICIT NONE

   INTEGER :: ios, mon, iadd, m

   INTEGER :: idate(4)

   INTEGER, DIMENSION (12) :: MonthLength = &
            (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

   CHARACTER (LEN=3), DIMENSION (12) :: MonthChar = (/ &
                      'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                      'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

   NAMELIST /InputDim/ Imax, Jmax, Date, GrADS, DirMain

   Imax=192
   Jmax=96
   Date='          '
   GrADS=.TRUE.
   DirMain='./ '

   OPEN (UNIT=nfinp, FILE='./'//VarName//'.nml', &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              './'//VarName//'.nml', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfinp, NML=InputDim)
   CLOSE (UNIT=nfinp)

   WRITE (UNIT=nfprt, FMT='(/,A)')  ' &InputDim'
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Imax = ', Imax
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Jmax = ', Jmax
   WRITE (UNIT=nfprt, FMT='(A)')    '     Date = '//Date
   WRITE (UNIT=nfprt, FMT='(A,L6)') '    GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')    '  DirMain = '//TRIM(DirMain)

   WRITE (nLats(3:7), '(I5.5)') Jmax

   WRITE (mskfmt(2:7), '(I6)') Imax

   ! Getting Date
   READ (date(1: 4), '(I4)') idate(4)
   READ (date(5: 6), '(I2)') idate(2)
   READ (date(7: 8), '(I2)') idate(3)
   READ (date(9:10), '(I2)') idate(1)

   WRITE (TimeGrADS(1: 2), '(I2.2)') idate(1)
   WRITE (TimeGrADS(4: 5), '(I2.2)') idate(3)
   WRITE (TimeGrADS(6: 8), '(A3)')   MonthChar(idate(2))
   WRITE (TimeGrADS(9:12), '(I4.4)') idate(4)

   ! Linear Time Interpolation Factors A and B
   mon=idate(2)
   IF (MOD(idate(4),4) == 0) MonthLength(2)=29
   MonthBefore=mon-1
   IF (idate(3) > MonthLength(mon)/2) MonthBefore=mon
   MonthAfter=MonthBefore+1
   IF (MonthBefore < 1) MonthBefore=12
   IF (MonthAfter > 12) MonthAfter=1
   iadd=MonthLength(MonthBefore)/2
   IF (MonthBefore == mon) iadd=-iadd
   FactorB=2.0_r8*REAL(idate(3)+iadd,r8)/ &
           REAL(MonthLength(MonthBefore)+MonthLength(MonthAfter),r8)
   FactorA=1.0_r8-FactorB

   WRITE (UNIT=nfprt, FMT='(/,A,3I3,I5,//,2(A,I3),//,2(A,F8.5),/)') &
         ' Climatological Snow for idate : ', idate, &
         ' MonthBefore = ', MonthBefore, ' MonthAfter = ', MonthAfter, &
         ' FactorA = ', FactorA, ' FactorB = ', FactorB

END SUBROUTINE InitInputParameters


END MODULE InputParameters

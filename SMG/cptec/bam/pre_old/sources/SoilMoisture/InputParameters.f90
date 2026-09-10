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

   REAL (KIND=r8), PUBLIC :: Undef, SeaValue, Lon0, Lat0

   INTEGER, PUBLIC :: Imax, Jmax, Idim, Jdim

   LOGICAL, PUBLIC :: PolarMean, GrADS, Linear

   LOGICAL, PUBLIC :: FlagInput(5), FlagOutput(5)

   INTEGER, DIMENSION (:,:), ALLOCATABLE, PUBLIC :: MaskInput

   CHARACTER (LEN=7), PUBLIC :: nLats='.G     '

   CHARACTER (LEN=10), PUBLIC :: mskfmt = '(      I1)'

   CHARACTER (LEN=16), PUBLIC :: NameLSMask='ModelLandSeaMask'

   CHARACTER (LEN=12), PUBLIC :: VarNameOut='SoilMoisture'

   CHARACTER (LEN=12), PUBLIC :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=13), PUBLIC :: DirModelIn='model/datain/'

   CHARACTER (LEN=32), PUBLIC :: VarName

   CHARACTER (LEN=528), PUBLIC :: DirMain

   CHARACTER (LEN=16) :: NameNML='SoilMoisture.nml'

   INTEGER, PUBLIC :: nferr=0    ! Standard Error Print Out
   INTEGER, PUBLIC :: nfinp=5    ! Standard Read In
   INTEGER, PUBLIC :: nfprt=6    ! Standard Print Out
   INTEGER, PUBLIC :: nflsm=10   ! To Read Land Sea Mask
   INTEGER, PUBLIC :: nfclm=20   ! To Read Climatological Soil Moisture Data
   INTEGER, PUBLIC :: nfout=30   ! To Write Intepolated Climatological Soil Moisture Data
   INTEGER, PUBLIC :: nfctl=40   ! To Write Output Data Description

   PUBLIC :: InitInputParameters


CONTAINS


SUBROUTINE InitInputParameters ()

   IMPLICIT NONE

   INTEGER :: ios

   NAMELIST /InputDim/ Imax, Jmax, Idim, Jdim, &
                       GrADS, Linear, VarName, DirMain

   Imax=192
   Jmax=96
   Idim=360
   Jdim=181
   GrADS=.TRUE.
   Linear=.TRUE.
   VarName='SoilMoistureClima '
   DirMain='./ '

   OPEN (UNIT=nfinp, FILE='./'//NameNML, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              './'//NameNML, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfinp, NML=InputDim)
   CLOSE (UNIT=nfinp)

   WRITE (UNIT=nfprt, FMT='(/,A)')  ' &InputDim'
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Imax = ', Imax
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Jmax = ', Jmax
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Idim = ', Idim
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Jdim = ', Jdim
   WRITE (UNIT=nfprt, FMT='(A,L6)') '    GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A,L6)') '   Linear = ', Linear
   WRITE (UNIT=nfprt, FMT='(A)')    '  VarName = '//TRIM(VarName)
   WRITE (UNIT=nfprt, FMT='(A)')    '  DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')  ' /'

   SeaValue=150.0_r8

   Undef=-999.0_r8

   ! For Linear Interpolation
   Lon0=0.0_r8  ! Start at Prime Meridian
   Lat0=90.0_r8 ! Start at North Pole

   ! For Area Weighted Interpolation
   ALLOCATE (MaskInput(Idim,Jdim))
   MaskInput=1
   PolarMean=.FALSE.
   FlagInput(1)=.TRUE.   ! Start at North Pole
   FlagInput(2)=.TRUE.   ! Start at Prime Meridian
   FlagInput(3)=.TRUE.   ! Latitudes Are at Center of Box
   FlagInput(4)=.TRUE.   ! Longitudes Are at Center of Box
   FlagInput(5)=.FALSE.  ! Regular Grid
   FlagOutput(1)=.TRUE.  ! Start at North Pole
   FlagOutput(2)=.TRUE.  ! Start at Prime Meridian
   FlagOutput(3)=.FALSE. ! Latitudes Are at North Edge of Box
   FlagOutput(4)=.TRUE.  ! Longitudes Are at Center of Box
   FlagOutput(5)=.TRUE.  ! Gaussian Grid

   WRITE (nLats(3:7), '(I5.5)') Jmax

   WRITE (mskfmt(2:7), '(I6)') Imax

END SUBROUTINE InitInputParameters


END MODULE InputParameters

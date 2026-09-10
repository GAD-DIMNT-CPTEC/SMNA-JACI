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

   INTEGER, PUBLIC :: Imax, Jmax, Idim, Jdim

   LOGICAL, PUBLIC :: PolarMean, GrADS

   LOGICAL, PUBLIC :: FlagInput(5), FlagOutput(5)

   INTEGER, DIMENSION (:,:), ALLOCATABLE, PUBLIC :: MaskInput

   CHARACTER (LEN=7), PUBLIC :: nLats='.G     '

   CHARACTER (LEN=10), PUBLIC :: mskfmt = '(      I1)'

   CHARACTER (LEN=11), PUBLIC :: NameLSM='LandSeaMask'

   CHARACTER (LEN=19), PUBLIC :: VarName='VegetationMaskClima'

   CHARACTER (LEN=16), PUBLIC :: NameLSMSSiB='ModelLandSeaMask'

   CHARACTER (LEN=14), PUBLIC :: VarNameVeg='VegetationMask'

   CHARACTER (LEN=12), PUBLIC :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=13), PUBLIC :: DirModelIn='model/datain/'

   CHARACTER (LEN=19) :: NameNML='VegetationMask.nml'

   CHARACTER (LEN=528), PUBLIC :: DirMain

   INTEGER, PUBLIC :: Undef=0

   REAL (KIND=r4), PUBLIC :: UndefG=-99.0

   INTEGER, PUBLIC :: nferr=0    ! Standard Error Print Out
   INTEGER, PUBLIC :: nfinp=5    ! Standard Read In
   INTEGER, PUBLIC :: nfprt=6    ! Standard Print Out
   INTEGER, PUBLIC :: nflsm=10   ! To Read Formatted Land Sea Mask
   INTEGER, PUBLIC :: nfcvm=20   ! To Read Unformatted Climatological Vegetation Mask
   INTEGER, PUBLIC :: nfvgm=30   ! To Write Unformatted Gaussian Grid Vegetation Mask
   INTEGER, PUBLIC :: nflsv=40   ! To Write Formatted Land Sea Mask Modified by Vegetation
   INTEGER, PUBLIC :: nfout=50   ! To Write GrADS Land Sea and Vegetation Mask
   INTEGER, PUBLIC :: nfctl=60   ! To Write GrADS Control File

   PUBLIC :: InitInputParameters

   INTEGER, PARAMETER, PUBLIC :: NumCat=13
   INTEGER, DIMENSION (NumCat), PUBLIC :: VegClass = &
            (/ 1, 1, 1, 1, 1, 1, 2, 2, 3, 2, 3, 4, 5 /)


CONTAINS


SUBROUTINE InitInputParameters ()

   IMPLICIT NONE

   INTEGER :: ios

   NAMELIST /InputDim/ Imax, Jmax, Idim, Jdim, &
                       GrADS, DirMain

   Imax=192
   Jmax=96
   Idim=360
   Jdim=180
   GrADS=.TRUE.
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
   WRITE (UNIT=nfprt, FMT='(A,I6)') '      Imax = ', Imax
   WRITE (UNIT=nfprt, FMT='(A,I6)') '      Jmax = ', Jmax
   WRITE (UNIT=nfprt, FMT='(A,I6)') '      Idim = ', Idim
   WRITE (UNIT=nfprt, FMT='(A,I6)') '      Jdim = ', Jdim
   WRITE (UNIT=nfprt, FMT='(A,L6)') '     GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')    '   DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')  ' /'

   ALLOCATE (MaskInput(Idim,Jdim))
   MaskInput=1
   PolarMean=.FALSE.
   FlagInput(1)=.TRUE.   ! Start at North Pole
   FlagInput(2)=.TRUE.   ! Start at Prime Meridian
   FlagInput(3)=.FALSE.  ! Latitudes Are at North Edge
   FlagInput(4)=.FALSE.  ! Longitudes Are at Western Edge
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

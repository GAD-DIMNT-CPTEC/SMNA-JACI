!
!  $Author: bonatti $
!  $Date: 2007/09/18 18:07:15 $
!  $Revision: 1.2 $
!
MODULE InputParameters

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: &
            r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
   INTEGER, PARAMETER, PUBLIC :: &
            r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers

   INTEGER, PUBLIC :: Idim, Jdim, Mend, Kmax, Imax, Jmax, Imx, JmaxHf, &
                      Mend1, Mend2, Mnwv0, Mnwv1, Mnwv2, Mnwv3

   INTEGER, PUBLIC :: Year, Month, Day, Hour

   REAL (KIND=r8), PUBLIC :: Undef, Lon0, Lat0, To, CO2SeaIce,SSTSeaIce, CO2OpenWater, SSTOpenWater,&
                             CO2SeaIceThreshold,SSTSeaIceThreshold, LapseRate, EMRad1, EMRad12

   LOGICAL, PUBLIC :: PolarMean, Linear, LinearGrid, GrADS

   LOGICAL, PUBLIC :: FlagInput(5), FlagOutput(5)

   INTEGER, DIMENSION (:,:), ALLOCATABLE, PUBLIC :: MaskInput

   CHARACTER (LEN=10), PUBLIC :: Trunc='T     L   '

   CHARACTER (LEN=7), PUBLIC :: nLats='.G     '

   CHARACTER (LEN=10), PUBLIC :: mskfmt = '(      I1)'

   CHARACTER (LEN=12), PUBLIC :: VarName='FLUXCO2Clima'

   CHARACTER (LEN=16), PUBLIC :: NameLSM='ModelLandSeaMask'

   CHARACTER (LEN=12), PUBLIC :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=13), PUBLIC :: DirModelIn='model/datain/'

   CHARACTER (LEN=12), PUBLIC :: DirClmFluxCO2='pre/databcs/' ! Climatological CO2 Datain Directory

   CHARACTER (LEN=12), PUBLIC :: DirClmSST='pre/databcs/' ! Climatological SST Datain Directory

   CHARACTER (LEN=11), PUBLIC :: FileClmFluxCO2='FluxCO2.bin'

   CHARACTER (LEN=11), PUBLIC :: FileClmSST='ersst.form'

   CHARACTER (LEN=16) :: NameNML='FLUXCO2Clima.nml'

   CHARACTER (LEN=7), PUBLIC :: Preffix

   CHARACTER (LEN=6), PUBLIC :: Suffix

   CHARACTER (LEN=10), PUBLIC :: DateICn

   CHARACTER (LEN=528), PUBLIC :: DirMain

   CHARACTER (LEN=3), DIMENSION (12), PUBLIC :: MonthChar = &
             (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

   INTEGER, PUBLIC :: nferr=0    ! Standard Error Print Out
   INTEGER, PUBLIC :: nfinp=5    ! Standard Read In
   INTEGER, PUBLIC :: nfprt=6    ! Standard Print Out
   INTEGER, PUBLIC :: nficn=10   ! To Read Topography from Initial Condition
   INTEGER, PUBLIC :: nflsm=20   ! To Read Formatted Land Sea Mask
   INTEGER, PUBLIC :: nfclm=30   ! To Read Formatted Climatological CO2
   INTEGER, PUBLIC :: nfclm2=31   ! To Read Formatted Climatological SST
   INTEGER, PUBLIC :: nfsto=40   ! To Write Unformatted Gaussian Grid CO2
   INTEGER, PUBLIC :: nfout=50   ! To Write GrADS Topography, Land Sea, Se Ice and Gauss CO2
   INTEGER, PUBLIC :: nfctl=60   ! To Write GrADS Control File

   PUBLIC :: InitInputParameters


CONTAINS


SUBROUTINE InitInputParameters ()

   IMPLICIT NONE

   INTEGER :: ios

   NAMELIST /InputDim/ Mend, Kmax, Idim, Jdim, &
                       SSTSeaIce, &
                       Linear, LinearGrid, GrADS, &
                       DateICn, Preffix, Suffix, DirMain

   Mend=62            ! Spectral Resolution Horizontal Truncation
   Kmax=28            ! Number of Layers of the Initial Condition for the Global Model
   Idim=360           ! Number of Longitudes For Climatological CO2 Data
   Jdim=180           ! Number of Latitudes For Climatological CO2 Data
   SSTSeaIce=-1.749_r8   ! SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
   Linear=.TRUE.         ! Flag to Bi-linear (T) or Area (F) Interpolation
   LinearGrid=.FALSE.    ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
   GrADS=.TRUE.          ! Flag for GrADS Outputs
   DateICn='yyyymmddhh'  ! Date of the Initial Condition for the Global Model
   Preffix='GANLCPT'     ! Preffix of the Initial Condition for the Global Model
   Suffix='S.unf.'       ! Suffix of the Initial Condition for the Global Model
   DirMain='./ '         ! Main Datain/Dataout Directory
   PRINT*,'passssss'
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

   WRITE (UNIT=nfprt, FMT='(/,A)')    ' &InputDim'
   WRITE (UNIT=nfprt, FMT='(A,I6)')   '       Mend = ', Mend
   WRITE (UNIT=nfprt, FMT='(A,I6)')   '       Kmax = ', Kmax
   WRITE (UNIT=nfprt, FMT='(A,I6)')   '       Idim = ', Idim
   WRITE (UNIT=nfprt, FMT='(A,I6)')   '       Jdim = ', Jdim
   WRITE (UNIT=nfprt, FMT='(A,F6.3)') '  SSTSeaIce = ', SSTSeaIce
   WRITE (UNIT=nfprt, FMT='(A,L6)')   '     Linear = ', Linear
   WRITE (UNIT=nfprt, FMT='(A,L6)')   ' LinearGrid = ', LinearGrid
   WRITE (UNIT=nfprt, FMT='(A,L6)')   '      GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')      '    DateICn = '//DateICn
   WRITE (UNIT=nfprt, FMT='(A)')      '    Preffix = '//Preffix
   WRITE (UNIT=nfprt, FMT='(A)')      '    DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')    ' /'

   CALL GetImaxJmax ()
   WRITE (UNIT=nfprt, FMT='(/,A,I6)') '   Imax = ', Imax
   WRITE (UNIT=nfprt, FMT='(A,I6,/)') '   Jmax = ', Jmax

   Mend1=Mend+1
   Mend2=Mend+2
   Mnwv2=Mend1*Mend2
   Mnwv0=Mnwv2/2
   Mnwv3=Mnwv2+2*Mend1
   Mnwv1=Mnwv3/2

   Imx=Imax+2
   JmaxHf=Jmax/2

   EMRad1=1.0_r8/6.37E6_r8
   EMRad12=EMRad1*EMRad1

   To=273.15_r8
   CO2OpenWater=0.0_r8
   CO2SeaIceThreshold=0.0_r8
   SSTOpenWater=-1.7_r8+To
   SSTSeaIceThreshold=271.2_r8

   LapseRate=0.0065_r8

   READ (DateICn, FMT='(I4,3I2)') Year, Month, Day, Hour

   Undef=-999.0_r8

   ! For Linear Interpolation
   Lon0=0.5_r8  ! Start Near Greenwhich
   Lat0=89.5_r8 ! Start Near North Pole

   ! For Area Weighted Interpolation
   ALLOCATE (MaskInput(Idim,Jdim))
   MaskInput=1
   PolarMean=.FALSE.
   FlagInput(1)=.TRUE.   ! Start at North Pole
   FlagInput(2)=.TRUE.   ! Start at Prime Meridian
   FlagInput(3)=.FALSE.  ! Latitudes Are at North Edge of Box
   FlagInput(4)=.FALSE.  ! Longitudes Are at Western Edge of Box
   FlagInput(5)=.FALSE.  ! Regular Grid
   FlagOutput(1)=.TRUE.  ! Start at North Pole
   FlagOutput(2)=.TRUE.  ! Start at Prime Meridian
   FlagOutput(3)=.FALSE. ! Latitudes Are at North Edge of Box
   FlagOutput(4)=.TRUE.  ! Longitudes Are at Center of Box
   FlagOutput(5)=.TRUE.  ! Gaussian Grid

   IF (LinearGrid) THEN
      Trunc(2:2)='L'
   ELSE
      Trunc(2:2)='Q'
   END IF
   WRITE (Trunc(3:6), FMT='(I4.4)') Mend
   WRITE (Trunc(8:10), FMT='(I3.3)') Kmax

   WRITE (nLats(3:7), '(I5.5)') Jmax

   WRITE (mskfmt(2:7), '(I6)') Imax

END SUBROUTINE InitInputParameters


SUBROUTINE GetImaxJmax ()

  IMPLICIT NONE

  INTEGER :: Nx, Nm, N2m, N3m, N5m, &
             n2, n3, n5, j, n, Check, Jfft

  INTEGER, SAVE :: Lfft=40000

  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: Ifft

  N2m=CEILING(LOG(REAL(Lfft,r8))/LOG(2.0_r8))
  N3m=CEILING(LOG(REAL(Lfft,r8))/LOG(3.0_r8))
  N5m=CEILING(LOG(REAL(Lfft,r8))/LOG(5.0_r8))
  Nx=N2m*(N3m+1)*(N5m+1)

  ALLOCATE (Ifft (Nx))
  Ifft=0

  n=0
  DO n2=1,N2m
     Jfft=(2**n2)
     IF (Jfft > Lfft) EXIT
     DO n3=0,N3m
        Jfft=(2**n2)*(3**n3)
        IF (Jfft > Lfft) EXIT
        DO n5=0,N5m
           Jfft=(2**n2)*(3**n3)*(5**n5)
           IF (Jfft > Lfft) EXIT
           n=n+1
           Ifft(n)=Jfft
        END DO
     END DO
  END DO
  Nm=n

  n=0
  DO 
     Check=0
     n=n+1
     DO j=1,Nm-1
        IF (Ifft(j) > Ifft(j+1)) THEN
           Jfft=Ifft(j)
           Ifft(j)=Ifft(j+1)
           Ifft(j+1)=Jfft
           Check=1
        END IF
     END DO
     IF (Check == 0) EXIT
  END DO

  IF (LinearGrid) THEN
     Jfft=2
  ELSE
     Jfft=3
  END IF
  Imax=Jfft*Mend+1
  DO n=1,Nm
     IF (Ifft(n) >= Imax) THEN
        Imax=Ifft(n)
        EXIT
     END IF
  END DO
  Jmax=Imax/2
  IF (MOD(Jmax, 2) /= 0) Jmax=Jmax+1

  DEALLOCATE (Ifft)

END SUBROUTINE GetImaxJmax


END MODULE InputParameters

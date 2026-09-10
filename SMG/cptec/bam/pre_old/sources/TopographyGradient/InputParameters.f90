MODULE InputParameters

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitParameters

  INTEGER, PARAMETER, PUBLIC :: r4=SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER, PUBLIC :: r8=SELECTED_REAL_KIND(15)

  INTEGER, PUBLIC :: Mend, Imax, Jmax, Kmax,Idim, Jdim

  REAL (KIND=r8), PUBLIC ::  Lon0, Lat0
  
  LOGICAL, PUBLIC :: PolarMean, GrADS, Linear
 
  LOGICAL, PUBLIC :: FlagInput(5), FlagOutput(5)

  INTEGER, DIMENSION (:,:), ALLOCATABLE, PUBLIC :: MaskInput

  INTEGER, PUBLIC :: Mend1, Mend2, Mend3, &
                     Mnwv2, Mnwv3, Mnwv0, Mnwv1, &
                     Imx,  Jmaxhf, &
                     MFactorFourier, MTrigsFourier

  INTEGER, PUBLIC :: nferr, nfprt, nfinp, nfout, nfctl,nfinpnml

  INTEGER (KIND=r4), PUBLIC :: ForecastDay

  REAL (KIND=r4), PUBLIC :: TimeOfDay

  REAL (KIND=r8), PUBLIC :: Undef, rad, dLon

  REAL (KIND=r8), PUBLIC :: EMRad   ! Earth Mean Radius (m)

  REAL (KIND=r8), PUBLIC :: EMRad1  ! 1/EMRad (1/m)

  REAL (KIND=r8), PUBLIC :: EMRad12 ! EMRad1**2 (1/m2)

  REAL (KIND=r8), PUBLIC :: EMRad2  ! EMRad**2 (m2)

  CHARACTER (LEN=12), PUBLIC :: TGrADS

  CHARACTER (LEN=33), PUBLIC :: FileInp

  CHARACTER (LEN=35), PUBLIC :: FileOut


  CHARACTER (LEN=3), DIMENSION (12), PUBLIC :: Months

  INTEGER, DIMENSION (4) :: InputDate

  CHARACTER (LEN=10), PUBLIC :: IDATE

  CHARACTER (LEN=3) :: Preffix

  CHARACTER (LEN=1) :: Grid
  CHARACTER (LEN=7),PUBLIC :: nLats='.G     '

  CHARACTER (LEN=10),PUBLIC :: VarNameT='HPRIME'
  CHARACTER (LEN=10),PUBLIC :: VarName='HPRIME'
 
  CHARACTER (LEN=528), PUBLIC :: DirMain

  CHARACTER (LEN=528),PUBLIC :: DirPreOut='pre/dataout/'

  CHARACTER (LEN=528),PUBLIC :: DirModelIn='model/datain/'

  INTEGER, PUBLIC :: nfclm=10   ! To Read Topography Data (GTOP30 or Navy)
  INTEGER, PUBLIC :: nfoub=30   ! To Write Intepolated Topography Data

CONTAINS


SUBROUTINE InitParameters ()

  IMPLICIT NONE

  NAMELIST /TopoGradNML/ Mend, Kmax,  Idim, Jdim,GrADS, Linear,IDATE, Preffix, Grid, &
                         DirMain

  Months=(/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
            'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

  ! Files Units

  nferr=0    ! Standard Error Print Out
  nfinpnml=5    ! Standard Read In
  nfprt=6    ! Standard Print Out
  nfinp=10   ! Unit To Read Input Data
  nfout=20   ! Unit To Write Output Data in GrADS Format
  nfctl=30   ! Unit To Write Control File For GrADS

  Mend=62  ! Spectral Horizontal Resolution of Output Data
  Kmax=28  ! Number of Levels of Data
  Idim=2161 ! Number of Longitudes in Navy Topography Data
  Jdim=1081 ! Number of Latitudes in Navy Topography Data
  GrADS=.TRUE. ! Flag for GrADS Outputs
  Linear=.TRUE. ! Flag for Linear (T) or Area Weighted (F) Interpolation
  IDATE='2004032600'! Date of Initial Condition
  Preffix = 'NMC' ! Preffix for Initial Condition File Name
  Grid = 'Q' ! Grid Type: Q for Quadratic and L for Linear Gaussian
  FileInp='FileInp.dat '     ! Name of File with Input Data
  FileOut='FileOut.dat '     ! Name of File with Output Data

! Imax - Number of Longitudes of Output Gaussian Grid Data
! Jmax - Number of Latitudes of Output Gaussian Grid Data

  OPEN  (UNIT=nfinpnml, FILE='TopographyGradient.nml', &
         FORM='FORMATTED', STATUS='OLD')
  READ  (UNIT=nfinpnml, NML=TopoGradNML)
  CLOSE (UNIT=nfinpnml)

  WRITE (UNIT=nfprt, FMT='(A)') ' '
  WRITE (UNIT=nfprt, NML=TopoGradNML)

  READ(IDATE(1: 4),'(I4)')InputDate(4)
  READ(IDATE(5: 6),'(I2)')InputDate(3)
  READ(IDATE(7: 8),'(I2)')InputDate(2)
  READ(IDATE(9:10),'(I2)')InputDate(1)
  
  CALL GetImaxJmax (Mend, Imax, Jmax)

  FileInp='GANL             S.unf.T     L   '
  FileInp(5:7)=Preffix
  WRITE(FileInp(8:11), FMT='(I4.4)') InputDate(4)
  WRITE(FileInp(12:13), FMT='(I2.2)') InputDate(3)
  WRITE(FileInp(14:15), FMT='(I2.2)') InputDate(2)
  WRITE(FileInp(16:17), FMT='(I2.2)') InputDate(1)
  FileInp(25:25)=Grid
  WRITE(FileInp(26:29), FMT='(I4.4)') Mend
  WRITE(FileInp(31:33), FMT='(I3.3)') Kmax

  FileOut='TopographyGradient          .G     '
  FileOut(19:28)=FileInp(8:17)
  WRITE(FileOut(31:35), FMT='(I5.5)') Jmax

  WRITE (UNIT=nfprt, FMT='(A)') ' '
  WRITE (UNIT=nfprt, FMT='(A,I5)') '  Imax = ', Imax
  WRITE (UNIT=nfprt, FMT='(A,I5)') '  Jmax = ', Jmax
  WRITE (UNIT=nfprt, FMT='(A,I6)') '  Idim = ', Idim
  WRITE (UNIT=nfprt, FMT='(A,I6)') '  Jdim = ', Jdim
  WRITE (UNIT=nfprt, FMT='(A,L6)') '  GrADS = ', GrADS
  WRITE (UNIT=nfprt, FMT='(A,L6)') '  Linear = ', Linear
  WRITE (UNIT=nfprt, FMT='(/,A)')  '  Numerical Precision (KIND): '
  WRITE (UNIT=nfprt, FMT='(A,I3)') '  r4  = ', r4
  WRITE (UNIT=nfprt, FMT='(A,I3)') '  r8  = ', r8
  WRITE (UNIT=nfprt, FMT='(A)') ' '

  WRITE (UNIT=nfprt, FMT='(/,2A)')  '  Input : ', FileInp
  WRITE (UNIT=nfprt, FMT='(2A)')    '  Output: ', FileOut
  WRITE (UNIT=nfprt, FMT='(A)') ' '

  Mend1=Mend+1
  Mend2=Mend+2
  Mend3=Mend+3
  Mnwv2=Mend1*Mend2
  Mnwv0=Mnwv2/2
  Mnwv3=Mnwv2+2*Mend1
  Mnwv1=Mnwv3/2

  Imx=Imax+2
  Jmaxhf=Jmax/2

  MFactorFourier=64
  MTrigsFourier=3*Imax/2

  rad=ATAN(1.0_r8)/45.0_r8
  dLon=360.0_r8/REAL(Imax,r8)

  EMRad=6.37E6_r8
  EMRad1=1.0_r8/EMRad
  EMRad12=EMRad1*EMRad1
  EMRad2=EMRad*EMRad

  Undef=-99999.0_r8
  
  ! For Linear Interpolation
  Lon0=0.0_r8  ! Start at Prime Meridian
  Lat0=90.0_r8-0.5_r8*(360.0_r8/REAL(Idim,r8)) ! Start at North Pole


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

   WRITE (nLats(3:7), FMT='(I5.5)') Jmax

END SUBROUTINE InitParameters


SUBROUTINE GetImaxJmax (Mend, Imax, Jmax)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Mend
  INTEGER, INTENT (OUT) :: Imax, Jmax

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

  Imax=3*Mend+1
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

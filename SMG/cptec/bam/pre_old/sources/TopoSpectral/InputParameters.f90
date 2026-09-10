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

   INTEGER, PUBLIC :: Mend, Imax, Jmax, Imx, JmaxHf, &
                      Mend1, Mend2, Mnwv0, Mnwv1, Mnwv2, Mnwv3

   REAL (KIND=r8), PUBLIC :: Undef, EMRad1, EMRad12

   CHARACTER (LEN=256), PUBLIC :: DirMain

   CHARACTER (LEN=10), PUBLIC :: VarNameT='Topography'

   CHARACTER (LEN=12), PUBLIC :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=13), PUBLIC :: DirModelIn='model/datain/'

   CHARACTER (LEN=16), PUBLIC :: NameNML='TopoSpectral.nml'

   CHARACTER (LEN=7), PUBLIC :: Trunc='.T     '

   CHARACTER (LEN=7), PUBLIC :: nLats='.G     '

   INTEGER, PUBLIC :: nferr=0    ! Standard Error Print Out
   INTEGER, PUBLIC :: nfinp=5    ! Standard Read In
   INTEGER, PUBLIC :: nfprt=6    ! Standard Print Out
   INTEGER, PUBLIC :: nftpi=10   ! To Read Reagular Grid Topography
   INTEGER, PUBLIC :: nftpo=20   ! To Write Unformatted Topography Sprectral Coefficients
   INTEGER, PUBLIC :: nfout=30   ! To Write GrADS Topography Data on a Gaussian Grid
   INTEGER, PUBLIC :: nfctl=40   ! To Write Output Data Description

   LOGICAL, PUBLIC :: LinearGrid, GrADS

   PUBLIC :: InitInputParameters


CONTAINS


SUBROUTINE InitInputParameters ()

   IMPLICIT NONE

   INTEGER :: ios

   NAMELIST /InputDim/ Mend, LinearGrid, GrADS, DirMain

   Mend=62            ! Spectral Resolution Horizontal Truncation
   LinearGrid=.FALSE. ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
   GrADS=.TRUE.       ! Flag to Get Recomposed Topography
   DirMain='./ '      ! Main Datain/Dataout Directory

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
   WRITE (UNIT=nfprt, FMT='(A,I6)')   '        Mend = ', Mend
   WRITE (UNIT=nfprt, FMT='(A,L6)')   '  LinearGrid = ', LinearGrid
   WRITE (UNIT=nfprt, FMT='(A,L6)')   '       GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')      '     DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')    ' /'

   CALL GetImaxJmax ()
   WRITE (UNIT=nfprt, FMT='(/,A,I6)') '        Imax = ', Imax
   WRITE (UNIT=nfprt, FMT='(A,I6,/)') '        Jmax = ', Jmax
   WRITE (UNIT=nfprt, FMT='(A)')      '    VarNameT = '//VarNameT

   IF (LinearGrid) THEN
      Trunc(3:3)='L'
   ELSE
      Trunc(3:3)='Q'
   END IF
   WRITE (Trunc(4:7), FMT='(I4.4)') Mend
   WRITE (nLats(3:7), '(I5.5)') Jmax

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

   Undef=-99999.0_r8

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

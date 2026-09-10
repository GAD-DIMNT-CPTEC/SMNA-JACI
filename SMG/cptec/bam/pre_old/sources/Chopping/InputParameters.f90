!
!  $Author: bonatti $
!  $Date: 2008/11/13 16:44:29 $
!  $Revision: 1.3 $
!
MODULE InputParameters

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitParameters

  ! SELECTED_INT_KIND(R):
  !   The value of the result is the kind type parameter value 
  !   of the integer type that can represent all integer values 
  !   n in the range -10**R < n < 10**R
  !   R must be a scalar of integer type(kind default)

  ! SELECTED_REAL_KIND(P,R):
  !   The value of the result is the kind type parameter value 
  !   of the real type that has a decimal precision greater than 
  !   or equal to P digits as returned by the function PRECISION 
  !   and a decimal exponent range greater than or equal to R 
  !   as returned by the function RANGE
  !   P (optional) must be a scalar of integer (kind default)
  !   R (optional) must be a scalar of integer (kind default)

  INTEGER, PARAMETER, PUBLIC :: i4=SELECTED_INT_KIND(9)
  INTEGER, PARAMETER, PUBLIC :: r4=SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER, PUBLIC :: r8=SELECTED_REAL_KIND(15)

  INTEGER, PUBLIC :: MendInp, ImaxInp, JmaxInp, KmaxInp, &
                     MendOut, ImaxOut, JmaxOut, KmaxOut, &
                     MendMin, MendCut

  INTEGER, PUBLIC :: Mnwv2Inp, Mnwv3Inp, &
                     Mend1Out, Mend2Out, Mend3Out, &
                     Mnwv2Out, Mnwv3Out, Mnwv0Out, Mnwv1Out, &
                     ImxOut,  JmaxhfOut, KmaxInpp, KmaxOutp, &
                     NTracers, Kdim, ICaseRec, ICaseDec, &
                     MFactorFourier, MTrigsFourier, Iter

  INTEGER, PUBLIC :: nferr, nfinp, nfprt, nftop, nfsig, nfnmc, &
                     nfcpt, nfozw, nftrw, nficr, nfozr, nftrr, &
                     nficw, nfozg, nftrg, nfgrd, nfctl

  INTEGER (KIND=i4), PUBLIC :: ForecastDay

  REAL (KIND=r4), PUBLIC :: TimeOfDay

  REAL (KIND=r8), PUBLIC :: cTv, SmthPerCut

  CHARACTER (LEN=128), PUBLIC :: DataCPT, DataInp, DataOut, DataTop, DataSig, DataSigInp, &
                                 GDASInp, OzonInp, TracInp, OzonOut, TracOut

  CHARACTER (LEN=1024), PUBLIC :: DirMain, DirInp, DirOut, DirTop, &
                                  DirSig, DirGrd, DGDInp, DirHome

  REAL (KIND=r8), PUBLIC :: EMRad   ! Earth Mean Radius (m)
  REAL (KIND=r8), PUBLIC :: EMRad1  ! 1/EMRad (1/m)
  REAL (KIND=r8), PUBLIC :: EMRad12 ! EMRad1**2 (1/m2)
  REAL (KIND=r8), PUBLIC :: EMRad2  ! EMRad**2 (m2)
  REAL (KIND=r8), PUBLIC :: Grav    ! Gravity Acceleration (m2/s)
  REAL (KIND=r8), PUBLIC :: Rd      ! Dry Air Gas Constant (m2/s2/K)
  REAL (KIND=r8), PUBLIC :: Rv      ! Water Vapour Air Gas Constant (m2/s2/K)
  REAL (KIND=r8), PUBLIC :: Cp      ! Dry Air Heat Capacity (m2/s2/K)
  REAL (KIND=r8), PUBLIC :: Lc      ! Latent Heat of Condensation (m2/s2)
  REAL (KIND=r8), PUBLIC :: Gamma   ! Mean Atmospheric Lapse Rate (K/m)
  REAL (KIND=r8), PUBLIC :: GEps    ! Precision For Constante Lapse Rate (No Dim)

  LOGICAL, PUBLIC :: GetOzone, GetTracers, &
                     GrADS, GrADSOnly, GDASOnly, &
                     SmoothTopo, LinearGrid

  CHARACTER (LEN=10), PUBLIC :: TruncInp

  CHARACTER (LEN=3), DIMENSION (12), PUBLIC :: MonChar

  CHARACTER (LEN=2) :: UTC

  CHARACTER (LEN=10) :: DateLabel


CONTAINS


SUBROUTINE InitParameters ()

  IMPLICIT NONE

  INTEGER :: ios

  CHARACTER (LEN=10) :: TruncOut, TrGrdOut

  CHARACTER (LEN=16) :: NCEPName

  CHARACTER (LEN=12) :: NameNML='Chopping.nml'

  LOGICAL :: ExistGANL, RmGANL,ExistGDAS

  NAMELIST /ChopNML/ MendInp,KmaxInp,MendOut, KmaxOut, MendMin, MendCut, Iter, SmthPerCut, &
                     GetOzone, GetTracers, GrADS, GrADSOnly, GDASOnly, &
                     SmoothTopo, RmGANL, LinearGrid, DateLabel, UTC, NCEPName, &
                     DirMain, DirHome

  MonChar=(/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' /)

  ! Files Units

  nferr=0    ! Standard Error Print Out
  nfinp=5    ! Standard Read In
  nfprt=6    ! Standard Print Out
  nftop=10   ! Unit To Read New Topography Field (Spectral Coefficients)
  nfsig=15   ! Unit To Read New Delta Sigma Data
  nfnmc=20   ! Unit To Read Initial Conditions From NCEP GDAS File
  nfcpt=25   ! Unit To Write Original Initial Conditions In CPTEC Format
  nfozw=30   ! Unit To Write Original Spectral Coefficients of Ozone
  nftrw=35   ! Unit To Write Original Spectral Coefficients of Other Tracers
  nficr=40   ! Unit To Read Original Initial Conditions In CPTEC Format
  nfozr=45   ! Unit To Read Original Spectral Coefficients of Ozone
  nftrr=50   ! Unit To Read Original Spectral Coefficients of Other Tracers
  nficw=55   ! Unit To Write Chopped Initial Conditions In CPTEC Format
  nfozg=60   ! Unit To Write Chopped Ozone Grid Field
  nftrg=65   ! Unit To Write Chopped Other Tracers Grid Fields
  nfgrd=70   ! Unit To Write Chopped Grid Fields For GrADS
  nfctl=75   ! Unit To Write Control File For GrADS

  MendInp =254! Spectral Horizontal Resolution of Input Data
  KmaxInp =64 ! Number of Layers of Input Data
  MendOut=213     ! Spectral Horizontal Resolution of Output Data
  KmaxOut=42      ! Number of Layers of Output Data
  MendMin=127     ! Minimum Spectral Resolution For Doing Topography Smoothing
  MendCut=0       ! Spectral Resolution Cut Off for Topography Smoothing
  Iter=10         ! Number of Iteractions in Topography Smoothing
  SmthPerCut=0.12_r8 ! Percentage for Topography Smoothing
  GetOzone=.FALSE.   ! Flag to Produce Ozone Files
  GetTracers=.FALSE. ! Flag to Produce Tracers Files
  GrADS=.TRUE.       ! Flag to Produce GrADS Files
  GrADSOnly=.FALSE.  ! Flag to Only Produce GrADS Files (Do Not Produce Inputs for Model)
  GDASOnly=.FALSE.   ! Flag to Only Produce Input CPTEC Analysis File
  SmoothTopo=.TRUE.  ! Flag to Performe Topography Smoothing
  RmGANL=.FALSE.     ! Flag to Remove GANL File if Desired
  LinearGrid=.FALSE. ! Flag to Set Linear (T) or Quadratic Gaussian Grid (F)
  DateLabel='yyyymmddhh' ! Date: yyyymmddhh or DateLabel='        hh'
                         !       If Year (yyyy), Month (mm) and Day (dd) Are Unknown
  UTC='hh'               ! UTC Hour: hh, Must Be Given if DateLabel='          ', else UTC=' '
  NCEPName='gdas1 '      ! NCEP Analysis Preffix for Input File Name
  DirMain='./ '          ! Main User Data Directory
  DirHome='./ '          ! Home User Sources Directory


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
  READ  (UNIT=nfinp, NML=ChopNML)
  CLOSE (UNIT=nfinp)

  WRITE (UNIT=nfprt, FMT='(A)')      ' '
  WRITE (UNIT=nfprt, FMT='(A)')      ' &ChopNML'
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  MendInp    = ', MendInp
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  KmaxInp    = ', KmaxInp
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  MendOut    = ', MendOut
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  KmaxOut    = ', KmaxOut
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  MendMin    = ', MendMin
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  MendCut    = ', MendCut
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  Iter       = ', Iter
  WRITE (UNIT=nfprt, FMT='(A,F7.3)') '  SmthPerCut = ', SmthPerCut
  WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GetOzone   = ', GetOzone
  WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GetTracers = ', GetTracers
  WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GrADS      = ', GrADS
  WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GrADSOnly  = ', GrADSOnly
  WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GDSAOnly   = ', GDASOnly
  WRITE (UNIT=nfprt, FMT='(A,L6)')   '  SmoothTopo = ', SmoothTopo
  WRITE (UNIT=nfprt, FMT='(A,L6)')   '  LinearGrid = ', LinearGrid
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DateLabel  = ', DateLabel
  WRITE (UNIT=nfprt, FMT='(2A)')     '  UTC        = ', UTC
  WRITE (UNIT=nfprt, FMT='(2A)')     '  NCEPName   = ', TRIM(NCEPName)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DirMain    = ', TRIM(DirMain)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DirHome    = ', TRIM(DirHome)
  WRITE (UNIT=nfprt, FMT='(A)')      ' /'

  DGDInp=TRIM(DirMain)//'pre/datain/ '
  DirInp=TRIM(DirMain)//'model/datain/ '
  DirOut=TRIM(DirMain)//'model/datain/ '
  DirTop=TRIM(DirMain)//'pre/dataTop/ '
  DirSig=TRIM(DirHome)//'sources/Chopping/ '
  DirGrd=TRIM(DirMain)//'pre/dataout/ '

  IF (DateLabel == '          ') THEN
     GDASInp=TRIM(NCEPName)//'.T'//UTC//'Z.SAnl'! Input NCEP Analysis File Name Without DateLabel
  ELSE
     UTC=DateLabel(9:10)
     GDASInp=TRIM(NCEPName)//'.T'//UTC//'Z.SAnl.'//DateLabel ! Input NCEP Analysis File Name
  END IF
  
  INQUIRE (FILE=TRIM(DGDInp)//TRIM(GDASInp), EXIST=ExistGDAS)
  
  IF(ExistGDAS)CALL GetGDASDateLabelRes ()


  WRITE (UNIT=nfprt, FMT='(/,A,I5)')   '  MendInp    = ', MendInp
  WRITE (UNIT=nfprt, FMT='(A,I5,/)')   '  KmaxInp    = ', KmaxInp

  CALL GetImaxJmax (MendInp, ImaxInp, JmaxInp)
  CALL GetImaxJmax (MendOut, ImaxOut, JmaxOut)

  IF (LinearGrid) THEN
     TruncInp='TL    L   '
     TruncOut='TL    L   '
  ELSE
     TruncInp='TQ    L   '
     TruncOut='TQ    L   '
  END IF
  WRITE (TruncInp(3:6),  FMT='(I4.4)') MendInp
  WRITE (TruncInp(8:10), FMT='(I3.3)') KmaxInp
  WRITE (TruncOut(3:6),  FMT='(I4.4)') MendOut
  WRITE (TruncOut(8:10), FMT='(I3.3)') KmaxOut

  TrGrdOut='G     L   '
  WRITE (TrGrdOut(2:6), FMT='(I5.5)') JmaxOut
  WRITE (TrGrdOut(8:10), FMT='(I3.3)') KmaxOut

  DataCPT='GANLNMC'//DateLabel//'S.unf.'//TruncInp ! Input CPTEC No Topo-Smoothed Analysis File Name
  IF (SmoothTopo) THEN
     DataInp='GANLSMT'//DateLabel//'S.unf.'//TruncInp ! Input Topo-Smoothed CPTEC Analysis File Name
     DataOut='GANLSMT'//DateLabel//'S.unf.'//TruncOut ! Output Topo-Smoothed CPTEC Analysis File Name
     OzonInp='OZONSMT'//DateLabel//'S.unf.'//TruncInp ! Input Ozone File Name
     TracInp='TRACSMT'//DateLabel//'S.unf.'//TruncInp ! Input Tracers File Name
     OzonOut='OZONSMT'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Ozone File Name
     TracOut='TRACSMT'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Tracers File Name
  ELSE
     DataInp='GANLNMC'//DateLabel//'S.unf.'//TruncInp ! Input No Topo-Smoothed CPTEC Analysis File Name
     DataOut='GANLNMC'//DateLabel//'S.unf.'//TruncOut ! Output No Topo-Smoothed CPTEC Analysis File Name
     OzonInp='OZONNMC'//DateLabel//'S.unf.'//TruncInp ! Input Ozone File Name
     TracInp='TRACNMC'//DateLabel//'S.unf.'//TruncInp ! Input Tracers File Name
     OzonOut='OZONNMC'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Ozone File Name
     TracOut='TRACNMC'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Tracers File Name
  END IF
  DataTop='Topography.'//TruncOut(1:6)         ! Input Topography File Name
  DataSigInp='DeltaSigma.'//TruncInp(7:10)     ! Delta Sigma File Input
  DataSig='DeltaSigma.'//TruncOut(7:10)        ! Delta Sigma File Output

  INQUIRE (FILE=TRIM(DirInp)//TRIM(DataInp), EXIST=ExistGANL)
  IF (ExistGANL .AND. RmGANL) THEN
     OPEN    (UNIT=nficr, FILE=TRIM(DirInp)//TRIM(DataInp))
     CLOSE   (UNIT=nficr, STATUS='DELETE')
     INQUIRE (FILE=TRIM(DirInp)//TRIM(DataInp), EXIST=ExistGANL)
     WRITE   (UNIT=nfprt, FMT='(/,A)') ' File Removed If False: '
     WRITE   (UNIT=nfprt, FMT='(A,L6)')  TRIM(DirInp)//TRIM(DataInp), ExistGANL
  END IF

  WRITE (UNIT=nfprt, FMT='(A)')      ' '
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  ImaxInp    = ', ImaxInp
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  JmaxInp    = ', JmaxInp
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  ImaxOut    = ', ImaxOut
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '  JmaxOut    = ', JmaxOut

  WRITE (UNIT=nfprt, FMT='(A)')      ' '
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DataCPT    = ', TRIM(DataCPT)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DirInp     = ', TRIM(DirInp)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DataInp    = ', TRIM(DataInp)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DirOut     = ', TRIM(DirOut)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DataOut    = ', TRIM(DataOut)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DirTop     = ', TRIM(DirTop)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DataTop    = ', TRIM(DataTop)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DirSig     = ', TRIM(DirSig)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DataSig    = ', TRIM(DataSig)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DirGrd     = ', TRIM(DirGrd)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  DGDInp     = ', TRIM(DGDInp)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  GDASInp    = ', TRIM(GDASInp)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  OzonInp    = ', TRIM(OzonInp)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  TracInp    = ', TRIM(TracInp)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  OzonOut    = ', TRIM(OzonOut)
  WRITE (UNIT=nfprt, FMT='(2A)')     '  TracOut    = ', TRIM(TracOut)

  WRITE (UNIT=nfprt, FMT='(/,A)')    '  Numerical Precision (KIND): '
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '          i4 = ', i4
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '          r4 = ', r4
  WRITE (UNIT=nfprt, FMT='(A,I5)')   '          r8 = ', r8
  WRITE (UNIT=nfprt, FMT='(A)')      ' '

  Mnwv2Inp=(MendInp+1)*(MendInp+2)
  Mnwv3Inp=Mnwv2Inp+2*(MendInp+1)

  Mend1Out=MendOut+1
  Mend2Out=MendOut+2
  Mend3Out=MendOut+3
  Mnwv2Out=Mend1Out*Mend2Out
  Mnwv0Out=Mnwv2Out/2
  Mnwv3Out=Mnwv2Out+2*Mend1Out
  Mnwv1Out=Mnwv3Out/2

  ImxOut=ImaxOut+2
  JmaxhfOut=JmaxOut/2
  KmaxInpp=KmaxInp+1
  KmaxOutp=KmaxOut+1

  NTracers=1
  Kdim=1
  ICaseRec=-1
  ICaseDec=1

  MFactorFourier=64
  MTrigsFourier=3*ImaxOut/2

  IF (MendCut <= 0 .OR. MendCut > MendOut) MendCut=MendOut

  ForecastDay=0_i4
  TimeOfDay=0.0_r4

  EMRad=6.37E6_r8
  EMRad1=1.0_r8/EMRad
  EMRad12=EMRad1*EMRad1
  EMRad2=EMRad*EMRad

  Grav=9.80665_r8
  Rd=287.05_r8
  Rv=461.50_r8
  Cp=1004.6_r8
  Lc=2.5E6_r8
  Gamma=-6.5E-3_r8
  GEps=1.E-9_r8
  cTv=Rv/Rd-1.0_r8

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


SUBROUTINE GetGDASDateLabelRes ()

  IMPLICIT NONE

  INTEGER :: ios

  INTEGER (KIND=i4), DIMENSION (4) :: Date

  REAL (KIND=r4) :: TimODay

  CHARACTER (LEN=1), DIMENSION (32) :: Descriptor

  REAL (KIND=r4) :: Header(205)

  WRITE (UNIT=nfprt, FMT='(/,A)') ' Getting Date Label from GDAS NCEP File'

  OPEN (UNIT=nfnmc, FILE=TRIM(DGDInp)//TRIM(GDASInp), FORM='UNFORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                      TRIM(TRIM(DGDInp)//TRIM(GDASInp)), &
                                      ' returned IOStat = ', ios
    STOP ' ** (Error) **'
  END IF

  ! Descriptor DateLabel (See NMC Office Note 85) (Not Used at CPTEC)
  READ (UNIT=nfnmc) Descriptor

  ! TimODay : Time of Day in Seconds
  ! Date(1) : UTC Hour
  ! Date(2) : Month
  ! Date(3) : Day
  ! Date(4) : Year
  READ  (UNIT=nfnmc) TimODay, Date, Header

  CLOSE (UNIT=nfnmc)

  ! DateLabel='yyyymmddhh'
  WRITE (DateLabel(1: 4), FMT='(I4.4)') Date(4)
  WRITE (DateLabel(5: 6), FMT='(I2.2)') Date(2)
  WRITE (DateLabel(7: 8), FMT='(I2.2)') Date(3)
  WRITE (DateLabel(9:10), FMT='(I2.2)') Date(1)

  ! UTC
  WRITE (UTC(1:2), FMT='(I2.2)') Date(1)

  ! Resolution
  MendInp=INT(Header(202))
  KmaxInp=INT(Header(203))

END SUBROUTINE GetGDASDateLabelRes


END MODULE InputParameters

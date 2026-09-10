!
!  $Author: bonatti $
!  $Date: 2008/11/13 16:44:29 $
!  $Revision: 1.3 $
!
MODULE InputParameters

  USE Parallelism, ONLY:   &
       myId

 use nemsio_module_mpi

 use nemsio_gfs

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
  INTEGER, PARAMETER, PUBLIC :: i8=SELECTED_INT_KIND(9)
  INTEGER, PARAMETER, PUBLIC :: r4=SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER, PUBLIC :: r8=SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER, PUBLIC :: r16=SELECTED_REAL_KIND(15)
  REAL (KIND=r8), PARAMETER, PUBLIC :: pai=3.14159265358979_r8
  REAL (KIND=r8), PARAMETER, PUBLIC :: twomg=1.458492e-4_r8

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

  INTEGER, PUBLIC :: ibdim_size, nproc_vert, tamBlock

  INTEGER (KIND=i4), PUBLIC :: IDVCInp
  INTEGER (KIND=i4), PUBLIC :: IDSLInp
  INTEGER (KIND=i4), PUBLIC :: ForecastDay

  REAL (KIND=r4), PUBLIC :: TimeOfDay

  REAL (KIND=r8), PUBLIC :: cTv, SmthPerCut

  CHARACTER (LEN=500), PUBLIC :: DataCPT, DataInp, DataOut, DataOup, DataTop,DataTopG, DataSig, DataSigInp, &
                                 GDASInp, OzonInp, TracInp, OzonOut, TracOut,DataGDAS

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
  REAL (KIND=r8), PUBLIC :: Gama    ! Mean Atmospheric Lapse Rate (K/m)
  REAL (KIND=r8), PUBLIC :: GEps    ! Precision For Constante Lapse Rate (No Dim)

  REAL (KIND=r8), PUBLIC :: RoCp    ! Rd over Cp
  REAL (KIND=r8), PUBLIC :: RoCp1   ! RoCp + 1
  REAL (KIND=r8), PUBLIC :: RoCpr   ! 1 / RoCp

  LOGICAL, PUBLIC :: GetOzone, GetTracers, &
                     GrADS, GrADSOnly, GDASOnly, &
                     SmoothTopo, LinearGrid, givenfouriergroups

  CHARACTER (LEN=10), PUBLIC :: TruncInp, TruncOut

  CHARACTER (LEN=3), DIMENSION (12), PUBLIC :: MonChar

  CHARACTER (LEN=2) :: UTC

  CHARACTER (LEN=10) :: DateLabel

  INCLUDE 'mpif.h'

CONTAINS


SUBROUTINE InitParameters (itertrc)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: itertrc
  INTEGER :: ios,im,jm,itrc

  CHARACTER (LEN =3) :: prefix
  
  CHARACTER (LEN=10) :: TrGrdOut
  CHARACTER (LEN=6 ) :: TrGrdOutGaus

  CHARACTER (LEN=16) :: NCEPName

  CHARACTER (LEN=255) :: NameNML='Chopping_parallel.nml'

  LOGICAL :: ExistGANL, RmGANL
 
  NAMELIST /ChopNML/ MendInp,KmaxInp,MendOut, KmaxOut, MendMin, MendCut, Iter, SmthPerCut, &
                     GetOzone, GetTracers, GrADS, GrADSOnly, GDASOnly,DataGDAS, &
                     SmoothTopo, RmGANL, LinearGrid, DateLabel, UTC, NCEPName, &
                     DirMain, DirHome, givenfouriergroups, &
                     nproc_vert, ibdim_size, tamBlock,prefix
                     

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

  MendOut=213     ! Spectral Horizontal Resolution of Output Data
  KmaxOut=42      ! Number of Layers of Output Data
  MendMin=127     ! Minimum Spectral Resolution For Doing Topography Smoothing
  MendCut=0       ! Spectral Resolution Cut Off for Topography Smoothing
  Iter=10         ! Number of Iteractions in Topography Smoothing
  nproc_vert = 1         ! Number of processors to be used in the vertical
                         ! (if givenfouriergroups set to TRUE)
  ibdim_size = 192       ! size of basic data block (ibmax)
  tamBlock = 512     ! number of fft's allocated in each block
  SmthPerCut=0.12_r8 ! Percentage for Topography Smoothing
  GetOzone=.FALSE.   ! Flag to Produce Ozone Files
  GetTracers=.FALSE. ! Flag to Produce Tracers Files
  GrADS=.TRUE.       ! Flag to Produce GrADS Files
  DataGDAS="Grid"
  GrADSOnly=.FALSE.  ! Flag to Only Produce GrADS Files (Do Not Produce Inputs for Model)
  GDASOnly=.FALSE.   ! Flag to Only Produce Input CPTEC Analysis File
  SmoothTopo=.TRUE.  ! Flag to Performe Topography Smoothing
  RmGANL=.FALSE.     ! Flag to Remove GANL File if Desired
  LinearGrid=.FALSE. ! Flag to Set Linear (T) or Quadratic Gaussian Grid (F)
  givenfouriergroups=.FALSE.! False if processor division should be automatic
  DateLabel='yyyymmddhh' ! Date: yyyymmddhh or DateLabel='        hh'
                         !       If Year (yyyy), Month (mm) and Day (dd) Are Unknown
  UTC='hh'               ! UTC Hour: hh, Must Be Given if DateLabel='          ', else UTC=' '
  NCEPName='gdas1 '      ! NCEP Analysis Preffix for Input File Name
  DirMain='./ '          ! Main User Data Directory
  DirHome='./ '          ! Home User Sources Directory
  prefix='NMC'
! MendInp : Spectral Horizontal Resolution of Input Data
! KmaxInp : Number of Layers of Input Data

  OPEN (UNIT=nfinp, FILE='./'//TRIM(NameNML), &
        FORM='FORMATTED', ACCESS='SEQUENTIAL', &
        ACTION='READ', STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
     WRITE (UNIT=nferr, FMT='(3A,I4)') &
           ' ** (Error) ** Open file ', &
             './'//TRIM(NameNML), &
           ' returned IOStat = ', ios
     STOP  ' ** (Error) **'
  END IF
  READ  (UNIT=nfinp, NML=ChopNML)
  CLOSE (UNIT=nfinp)

  IF (myid.eq.0) THEN
    WRITE (UNIT=nfprt, FMT='(A)')      ' '
    WRITE (UNIT=nfprt, FMT='(A)')      ' &ChopNML'
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  MendOut    = ', MendOut
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  KmaxOut    = ', KmaxOut
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  MendMin    = ', MendMin
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  MendCut    = ', MendCut
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  Iter       = ', Iter
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  nproc_vert = ', nproc_vert
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  ibdim_size = ', ibdim_size
    WRITE (UNIT=nfprt, FMT='(A,I5)')   '  tamBlock   = ', tamBlock
    WRITE (UNIT=nfprt, FMT='(A,F7.3)') '  SmthPerCut = ', SmthPerCut
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GetOzone   = ', GetOzone
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GetTracers = ', GetTracers
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GrADS      = ', GrADS
    WRITE (UNIT=nfprt, FMT='(A,A )')   '  DataGDAS   = ', DataGDAS
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GrADSOnly  = ', GrADSOnly
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  GDASOnly   = ', GDASOnly
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  SmoothTopo = ', SmoothTopo
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  LinearGrid = ', LinearGrid
    WRITE (UNIT=nfprt, FMT='(A,L6)')   '  givenfouriergroups = ', givenfouriergroups
    WRITE (UNIT=nfprt, FMT='(2A)')     '  DateLabel  = ', DateLabel
    WRITE (UNIT=nfprt, FMT='(2A)')     '  UTC        = ', UTC
    WRITE (UNIT=nfprt, FMT='(2A)')     '  NCEPName   = ', TRIM(NCEPName)
    WRITE (UNIT=nfprt, FMT='(2A)')     '  DirMain    = ', TRIM(DirMain)
    WRITE (UNIT=nfprt, FMT='(2A)')     '  DirHome    = ', TRIM(DirHome)
    WRITE (UNIT=nfprt, FMT='(2A)')     '  prefix     = ', TRIM(prefix)
    WRITE (UNIT=nfprt, FMT='(A)')      ' /'
  ENDIF

  DGDInp=TRIM(DirMain)//'pre/datain/ '
  DirInp=TRIM(DirMain)//'model/datain/ '
  DirOut=TRIM(DirMain)//'model/datain/ '
  DirTop=TRIM(DirMain)//'pre/dataTop/ '
!  DirTop=TRIM(DirMain)//'pre/dataout/ '
! DirSig=TRIM(DirHome)//'sources/Chopping/ '
  DirSig=TRIM(DirMain)//'pre/datain/ '
  DirGrd=TRIM(DirMain)//'pre/dataout/ '

  IF(TRIM(DataGDAS) == 'Spec')THEN
     IF (DateLabel == '          ') THEN
        GDASInp=TRIM(NCEPName)//'.T'//UTC//'Z.SAnl'! Input NCEP Analysis File Name Without DateLabel
     ELSE
        UTC=DateLabel(9:10)
        GDASInp=TRIM(NCEPName)//'.T'//UTC//'Z.SAnl.'//DateLabel ! Input NCEP Analysis File Name
     END IF
  ELSE
     IF (DateLabel == '          ') THEN
        GDASInp=TRIM(NCEPName)//'.T'//UTC//'Z.atmanl.nemsio.' ! Input NCEP Analysis File Name
     ELSE
        UTC=DateLabel(9:10)
        GDASInp=TRIM(NCEPName)//'.T'//UTC//'Z.atmanl.nemsio.'//DateLabel ! Input NCEP Analysis File Name
     END IF
  END IF

  IF(.NOT.GDASOnly)THEN
     CALL GetGDASDateLabelRes ()
  END IF

  IF (myid.eq.0) THEN
    WRITE (UNIT=nfprt, FMT='(/,A,I5)')   '  MendInp    = ', MendInp
    WRITE (UNIT=nfprt, FMT='(A,I5,/)')   '  KmaxInp    = ', KmaxInp
  END IF
  IF(TRIM(DataGDAS) == 'Spec')THEN
     CALL GetImaxJmax (MendInp, ImaxInp, JmaxInp)
     CALL GetImaxJmax (MendOut, ImaxOut, JmaxOut)

  ELSE

     IF(itertrc ==1)THEN
        DO itrc=1,40000
           CALL GetImaxJmax (itrc, im, jm)
           IF((im>ImaxInp) .and. (jm>JmaxInp))exit
        END DO
        ImaxOut=ImaxInp
        JmaxOut=JmaxInp
        MendInp=itrc
        MendOut=itrc
     ELSE
        IF(.NOT.GDASOnly)THEN
           DO itrc=1,40000
              CALL GetImaxJmax (itrc, im, jm)
              IF((im>ImaxInp) .and. (jm>JmaxInp))exit
           END DO
           MendInp=itrc
           CALL GetImaxJmax (MendOut, ImaxOut, JmaxOut)
        ELSE 
           CALL GetImaxJmax (MendOut, ImaxOut, JmaxOut)
        END IF
     END IF
  END IF  
 
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
  TrGrdOutGaus='G     '
  WRITE (TrGrdOutGaus(2:6), FMT='(I5.5)') JmaxOut
  WRITE (TrGrdOut(2:6), FMT='(I5.5)') JmaxOut
  WRITE (TrGrdOut(8:10), FMT='(I3.3)') KmaxOut

  DataCPT='GANLNMC'//DateLabel//'S.unf.'//TruncInp ! Input CPTEC No Topo-Smoothed Analysis File Name
  IF (SmoothTopo) THEN
     DataInp='GANLSMT'//DateLabel//'S.unf.'//TruncInp ! Input Topo-Smoothed CPTEC Analysis File Name
     DataOut='GANLSMT'//DateLabel//'S.unf.'//TruncOut ! Output Topo-Smoothed CPTEC Analysis File Name
     DataOup='GANLSMT'//DateLabel//'S.unf.'//TruncOut ! Output Topo-Smoothed CPTEC Analysis File Name
     OzonInp='OZONSMT'//DateLabel//'S.unf.'//TruncInp ! Input Ozone File Name
     TracInp='TRACSMT'//DateLabel//'S.unf.'//TruncInp ! Input Tracers File Name
     OzonOut='OZONSMT'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Ozone File Name
     TracOut='TRACSMT'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Tracers File Name
  ELSE
     DataInp='GANLNMC'//DateLabel//'S.unf.'//TruncInp ! Input No Topo-Smoothed CPTEC Analysis File Name
     DataOut='GANLNMC'//DateLabel//'S.unf.'//TruncOut ! Output No Topo-Smoothed CPTEC Analysis File Name
     DataOup='GANLNMC'//DateLabel//'S.unf.'//TruncOut ! Output Topo-Smoothed CPTEC Analysis File Name
     OzonInp='OZONNMC'//DateLabel//'S.unf.'//TruncInp ! Input Ozone File Name
     TracInp='TRACNMC'//DateLabel//'S.unf.'//TruncInp ! Input Tracers File Name
     OzonOut='OZONNMC'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Ozone File Name
     TracOut='TRACNMC'//DateLabel//'S.grd.'//TrGrdOut ! Grid Ouput Tracers File Name
  END IF
  DataTop='Topography2.'//TruncOut(1:6)         ! Input Topography Spec File Name
  DataTopG='Topography2.'//TRIM(TrGrdOutGaus)   ! Input Topography Grid File Name
  DataSigInp='DeltaSigma.'//TruncInp(7:10)     ! Delta Sigma File Input
  DataSig='DeltaSigma.'//TruncOut(7:10)        ! Delta Sigma File Output

  INQUIRE (FILE=TRIM(DirInp)//TRIM(DataInp), EXIST=ExistGANL)
  IF (ExistGANL .AND. RmGANL) THEN
     IF (myid.eq.0) THEN
       OPEN    (UNIT=nficr, FILE=TRIM(DirInp)//TRIM(DataInp))
       CLOSE   (UNIT=nficr, STATUS='DELETE')
       INQUIRE (FILE=TRIM(DirInp)//TRIM(DataInp), EXIST=ExistGANL)
       WRITE   (UNIT=nfprt, FMT='(/,A)') ' File Removed If False: '
       WRITE   (UNIT=nfprt, FMT='(A,L6)')  TRIM(DirInp)//TRIM(DataInp), ExistGANL
      ELSE
       ExistGANL = .FALSE.
     ENDIF
  END IF

  IF (myid.eq.0) THEN
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
    WRITE (UNIT=nfprt, FMT='(2A)')     '  DataTopG   = ', TRIM(DataTopG)
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
  END IF

  KmaxInpp=KmaxInp+1
  KmaxOutp=KmaxOut+1

  NTracers=1
  Kdim=1
  ICaseRec=-1
  ICaseDec=1

  IF (MendCut <= 0 .OR. MendCut > MendOut) MendCut=MendOut

  IDVCInp=0_i4
  IDSLInp=0_i4

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
  Gama=-6.5E-3_r8
  GEps=1.E-9_r8
  cTv=Rv/Rd-1.0_r8

  RoCp=Rd/Cp
  RoCp1=RoCp+1.0_r8
  RoCpr=1.0_r8/RoCp

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


  INTEGER (KIND=i4), DIMENSION (4) :: Date

  REAL (KIND=r4) :: TimODay

  CHARACTER (LEN=1), DIMENSION (32) :: Descriptor
  REAL (KIND=r4)   , DIMENSION (2*100+1) :: SiSl
  REAL (KIND=r4)   , DIMENSION (44)      :: Extra

  REAL (KIND=r4) :: Header(205)
!GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-
  integer             :: iret
  type(nemsio_gfile)  :: gfile
  character(255)      :: cin
  integer             :: ios
  type(nemsio_head)   :: gfshead

!GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-

  if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A)') ' Getting Date Label from GDAS NCEP File'
!GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-

  IF(TRIM(DataGDAS) == 'Spec')THEN
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
    ! Sigma Interfaces and Layers - SiSl(1:201)
    ! Extra Information           - Extra(1:44)
    !  o ID Sigma Structure       - Extra(18)
    !    = 1 Phillips
    !    = 2 Mean
    !  o ID Vertical Coordinate   - Extra(19)
    !    = 1 Sigma (0 for old files)
    !    = 2 Sigma-p
    READ  (UNIT=nfnmc) TimODay, Date, SiSl, Extra
    CLOSE (UNIT=nfnmc)

    ! DateLabel='yyyymmddhh'
    WRITE (DateLabel(1: 4), FMT='(I4.4)') Date(4)
    WRITE (DateLabel(5: 6), FMT='(I2.2)') Date(2)
    WRITE (DateLabel(7: 8), FMT='(I2.2)') Date(3)
    WRITE (DateLabel(9:10), FMT='(I2.2)') Date(1)

    ! UTC
    WRITE (UTC(1:2), FMT='(I2.2)') Date(1)

    ! Resolution
    MendInp=INT(Extra(1))
    KmaxInp=INT(Extra(2))

    IDSLInp=INT(Extra(18),i4)
    IDVCInp=INT(Extra(19),i4)


  ELSE
    !
    !Inicializa a lib nemsio
    call nemsio_init(iret=iret)
    ! Abertura do arquivo nemsio

    call nemsio_open(gfile,TRIM(DGDInp)//TRIM(GDASInp),'READ',MPI_COMM_WORLD,iret=iret)
    
    !  open (read) nemsio grid file headers
 
        call nemsio_getfilehead(gfile,                                        &
               idate=gfshead%idate, nfhour=gfshead%nfhour, nfminute=gfshead%nfminute, &
               nfsecondn=gfshead%nfsecondn, nfsecondd=gfshead%nfsecondd,      &
               version=gfshead%version, nrec=gfshead%nrec, dimx=gfshead%dimx, &
               dimy=gfshead%dimy, dimz=gfshead%dimz, jcap=gfshead%jcap,       &
               ntrac=gfshead%ntrac, ncldt=gfshead%ncldt, nsoil=gfshead%nsoil, &
               idsl=gfshead%idsl, idvc=gfshead%idvc, idvm=gfshead%idvm,       &
               idrt=gfshead%idrt, extrameta=gfshead%extrameta,                &
               nmetavari=gfshead%nmetavari, nmetavarr=gfshead%nmetavarr,      &
               nmetavarl=gfshead%nmetavarl, nmetavarr8=gfshead%nmetavarr8,    &
               nmetaaryi=gfshead%nmetaaryi, nmetaaryr=gfshead%nmetaaryr,      &
               iret=ios)
 
             if(myid.eq.0) THEN
                PRINT*,'idate=',gfshead%idate
                PRINT*,'nfhour=',gfshead%nfhour
                PRINT*,'nfminute=',gfshead%nfminute
                PRINT*,'nfsecondn=',gfshead%nfsecondn
                PRINT*,'nfsecondd=',gfshead%nfsecondd
                PRINT*,'version=',gfshead%version
                PRINT*,'nrec=',gfshead%nrec
                PRINT*,'dimx=',gfshead%dimx
                PRINT*,'dimy=',gfshead%dimy
                PRINT*,'dimz=',gfshead%dimz
                PRINT*,'jcap=',gfshead%jcap
                PRINT*,'ntrac=',gfshead%ntrac
                PRINT*,'ncldt=',gfshead%ncldt
                PRINT*,'nsoil=',gfshead%nsoil
                PRINT*,'idsl=',gfshead%idsl
                PRINT*,'idvc=',gfshead%idvc
                PRINT*,'idvm=',gfshead%idvm
                PRINT*,'idrt=',gfshead%idrt
                PRINT*,'extrameta=',gfshead%extrameta
                PRINT*,'nmetavari=',gfshead%nmetavari
                PRINT*,'nmetavarr=',gfshead%nmetavarr
                PRINT*,'nmetavarl=',gfshead%nmetavarl
                PRINT*,'nmetavarr8=',gfshead%nmetavarr8
                PRINT*,'nmetaaryi=',gfshead%nmetaaryi
                PRINT*,'nmetaaryr=',gfshead%nmetaaryr
                PRINT*,'iret=',ios

             END IF

             TimODay         =gfshead%nfsecondn
             Date   (1) =gfshead%idate(4)
             Date   (3) =gfshead%idate(3)
             Date   (2) =gfshead%idate(2)
             Date   (4) =gfshead%idate(1)

               !   SiSl
               !   Extra

               ! TimODay : Time of Day in Seconds
               ! Date(1) : UTC Hour
               ! Date(2) : Month
               ! Date(3) : Day
               ! Date(4) : Year
               ! Sigma Interfaces and Layers - SiSl(1:201)
               ! Extra Information           - Extra(1:44)
               !  o ID Sigma Structure       - Extra(18)
               !    = 1 Phillips
               !    = 2 Mean
               !  o ID Vertical Coordinate   - Extra(19)
               !    = 1 Sigma (0 for old files)
               !    = 2 Sigma-p


        !Fecha o arquivo nemsio
        
        call nemsio_close(gfile,iret=iret)
 
      !Finaliza
    
      call nemsio_finalize()


      ! DateLabel='yyyymmddhh'
      WRITE (DateLabel(1: 4), FMT='(I4.4)') Date(4)
      WRITE (DateLabel(5: 6), FMT='(I2.2)') Date(2)
      WRITE (DateLabel(7: 8), FMT='(I2.2)') Date(3)
      WRITE (DateLabel(9:10), FMT='(I2.2)') Date(1)

      ! UTC
      WRITE (UTC(1:2), FMT='(I2.2)') Date(1)

      ! Resolution
      MendInp=gfshead%jcap
      KmaxInp=gfshead%dimz
      ImaxInp=gfshead%dimx
      JmaxInp=gfshead%dimy
      IDSLInp=gfshead%idsl
      IDVCInp=gfshead%idvc

  ! Resolution

  END IF

END SUBROUTINE GetGDASDateLabelRes

END MODULE InputParameters

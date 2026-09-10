!
!  $Author: pkubota $
!  $Date: 2010/04/20 20:18:04 $
!  $Revision: 1.10 $
!
MODULE Transform


  ! LEGENDRE AND FOURIER TRANSFORMS
  !
  ! PURPOSE:
  ! Transform a set of fields from spectral reprentation to 
  ! grid point representation and backwards. Computes 
  ! derivatives on the zonal direction (longitudes).
  ! 
  ! IMPLEMENTATION GUIDELINES:
  ! Legendre Transform in matrix form.
  ! Efficiency on vector machines requires large matrices.
  ! Fields are packed together to generate large matrices
  ! prior to the transform.
  ! Version implements transform for Gaussian and Reduced Grids.
  ! OpenMP parallelism. MPI-ready data structure.
  !
  ! MODULE CREATION:
  ! InitTransform should be invoked once, prior to any other
  ! module routine. It needed to be invoked only once at each
  ! execution. It gathers model truncation and formulation info from
  ! the environment. Due to that, truncation cannot be changed during
  ! model execution.
  !
  ! GRID TO SPECTRAL USER INTERFACE:
  ! User interface contains 4 routines on the grid to
  ! spectral transform, with the following semantics:
  !
  ! CreateGridToSpec, that prepares internal data structure (arrays)
  !                   to perform Fourier and then Legendre Transforms
  !                   for a given number (input arguments) of surface and
  !                   full fields.
  !
  ! DepositGridToSpec, that gives the grid field to be transformed and
  !                    the spectral field that will receive the transformed
  !                    field. Transform is not performed by this call; it
  !                    only specifies where is input data is and where output
  !                    data will be posted.
  !
  ! DoGridToSpec, that performs the transform over all previously deposited
  !               grid fields, since last CreateGridToSpec, saving output
  !               data on the pointed positions. The number of surface grid fields
  !               and full grid fields deposited should match the information
  !               given at CreateGridToSpec. At the end of DoGridToSpec, all
  !               Spectral Fields informed by DepositGridToSpec are filled with
  !               required information.
  !
  ! DestroyGridToSpec, that removes the internal data structure created by
  !                    CreateGridToSpec.
  !
  ! GRID TO SPECTRAL USER INTERFACE ORDER OF INVOCATION:
  ! CreateGridToSpec should be invoked once, prior to any other routine.
  ! Use one DepositGridToSpec for each Grid field to be transformed.
  ! Once all Deposits are done, invoke DoGridToSpec. Finalize by invoking
  ! DestroyGridToSpec.
  !
  ! SPECTRAL TO GRID USER INTERFACE:
  ! Similar to Grid to Spec, enlarged to include zonal derivative computation.
  !
  ! CreateSpecToGrid: Similar functionality, except that the number of output grid
  !                   fields can be larger that the number of input spectral
  !                   fields, to accomodate zonal derivatives.
  ! DepositSpecToGrid: Use only to compute transform of input spectral field, 
  !                    without derivatives.
  ! DepositSpecToDelLamGrid: Use only to compute the zonal derivative
  !                    of input spectral field.
  ! DepositSpecToGridAndDelLamGrid: Use to compute both the transform and
  !                    the zonal derivative of input spectral field.
  ! DoSpecToGrid: Perform all required transforms and zonal derivatives
  !               simultaneosly.
  ! DestroySpecToGrid: Similar functionality.
  !
  ! SPECTRAL TO GRID USER INTERFACE ORDER OF INVOCATION:
  ! CreateSpecToGrid should be invoked once, prior to any other routine.
  ! Use one of DepositSpecToGrid, DepositSpecToDelLamGrid ou DepositSpecToGridAndDelLamGrid
  ! for each Grid field to be transformed, with or without zonal derivative.
  ! Once all Deposits are done, invoke DoSpecToGrid. Finalize by invoking
  ! DestroySpecToGrid.
  !
  ! OTHER FUNCTIONALITIES:
  ! NextSizeFFT returns the smallest size of FFT equal or larger than input argument.
  ! Usefull for array dimensions.


  USE InputParameters, ONLY: &
       EMRad1,   &
       r8,       &
       nfprt,    &
       tamBlock

  USE Utils, ONLY: &
       LegFuncS2F, &
       GaussWeights, &
       NoBankConflict, &
       allpolynomials

  USE Sizes, ONLY:     &
       mMax,           &
       nMax,           &
       nExtMax,        &
       mnMax,          &
       mnExtMax,       &
       mnExtMap,       &
       nodeHasM,       &
       nodeHasJ_f,     &
       lm2m,           &
       myMMax,         &
       myMNMax,        &
       myMNExtMax,     &
       myMMap,         &
       myNMap,         &
       myMNMap,        &
       myMExtMap,      &
       myNExtMap,      &
       myMNExtMap,     &
       havesurf,       &
       myfirstlev,     &
       mylastlev,      &
       myfirstlat_f,   &
       mylastlat_f,    &
       MsPerProc,      &
       firstlatinproc_f, &
       lastlatinproc_f,  &
       nlatsinproc_f,    &
       messages_f,     &
       messages_g,     &
       messproc_f,     &
       messproc_g,     &
       nrecs_f,        &
       nrecs_g,        &
       kfirst_four,    &
       klast_four,     &
       myfirstlat_f,   &
       myfirstlat,     &
       mylastlat_f,    &
       mylastlat,      &
       myfirstlon,     &
       mylastlon,      &
       iMax,           &
       jMax,           &
       jMaxHalf,       &
       kMaxloc,        &
       ibMax,          &
       ibMaxPerJB,     &
       jbMax,          &
       ibPerIJ,        &
       jbPerIJ,        &
       jMinPerM,       &
       jMaxPerM,       &
       iMaxPerJ,       &
       myJMax_f,         &
       JMaxlocal_f,    &
       MMaxlocal,      &
       ThreadDecomp

  USE Parallelism, ONLY: &
       myid,             &
       myid_four,        &
       mygroup_four,     &
       COMM_FOUR,        &
       maxnodes,         &
       maxNodes_four

  USE Communications, ONLY: &
       bufsend,          &
       bufrec,           &
       dimrecbuf,        &
       dimsendbuf

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  PRIVATE
  PUBLIC :: InitTransform
  PUBLIC :: CreateSpecToGrid
  PUBLIC :: DepositSpecToGrid
  PUBLIC :: DepositSpecToDelLamGrid
  PUBLIC :: DepositSpecToGridAndDelLamGrid
  PUBLIC :: DoSpecToGrid
  PUBLIC :: DestroySpecToGrid
  PUBLIC :: CreateGridToSpec
  PUBLIC :: DepositGridToSpec
  PUBLIC :: DepositGridToSpec_PK
  PUBLIC :: DoGridToSpec
  PUBLIC :: DestroyGridToSpec
  PUBLIC :: NextSizeFFT
  PUBLIC :: Clear_Transform

  ! Legendre Transform
  ! 
  ! Internal representation of Spectral Fields
  !
  ! Array Spe(dlmn, dv), using Spe(myMNExtMax, 2*nVertSpec)
  ! dlmn and dv avoid memory bank conflict
  !
  ! second dimension of Spe maps real and imaginary verticals 
  ! of all fields to be transformed to an 1D structure 
  ! in array element order of:(Vertical, Field, RealImaginary).
  ! 
  ! first dimension of Spe maps (lm,n) into an 1D structure st
  ! n goes faster than lm; for each lm, all n are represented.
  ! lm are stored sequentially.
  ! For each lm, first store all n st n+m is even, followed by
  ! all n st n+m is odd.
  !
  !  lm2m(1)
  ! S(n,v)
  !  even
  !
  !  lm2m(1)
  ! S(n,v)
  !  odd
  !
  !  lm2m(2)
  ! S(n,v)
  !  even
  !
  !  lm2m(2)
  ! S(n,v)
  !  odd
  !
  ! .......
  !
  !  lm2m(myMMax)
  ! S(n,v)
  !  even
  !
  !  lm2m(myMMax)
  ! S(n,v)
  !  odd
  !
  ! PAD
  !
  ! nEven(lm) has the number of spectral coeficients (real or imag) for m+n even
  ! nOdd(lm) has the number of spectral coeficients (real or imag) for m+n odd
  !
  ! dnEven(lm) is the dimensioning of nEven(lm) to avoid bank conflicts
  ! dnOdd(lm) is the dimensioning of nEven(lm) to avoid bank conflicts
  !
  ! firstNEven(lm) points to the first row of even n's for lm
  ! firstNOdd(lm) points to the first row of odd n's for lm
  ! 
  ! lmnExtMap(lmn), lmn = 1,...,myMNExtMax maps the external representation
  ! of Spectral Extended Fields to the transform internal representation,
  ! for a single vertical
  !
  ! lmnMap(lmn), lmn = 1,...,myMNMax maps the external representation
  ! of Spectral Regular Fields to the transform internal representation,
  ! for a single vertical
  !
  ! lastv stores the last used vertical of the real part of S; 
  ! first free vertical (of the real part) is lastv+1
  !
  ! nVertSpec stores the distance (in verticals) from the real to the imaginary part
  !
  ! first vertical (of the real part) of each stored field is thisv


  ! Fourier Output of Matrix Multiplications
  !
  ! Arrays FoEv(dv, djh) and FoOd(dv, djh)
  !
  ! first dimension of FoEv, FoOd maps real and imaginary verticals 
  ! of all fields to an 1D structure 
  ! in array element order of:(Vertical, Field, RealImaginary),
  ! identical to the second dimension of Spe.
  ! 
  ! second dimension of FEO contains latitudes, from 1 to
  ! jMaxHf on even fields, and from jMax to jMaxHf+1 on odd
  ! fields. Not all latitudes are used; for each m, latitudes
  ! jMinPerM(m) to jMaxHf on even fields and jMax-jMinPerM(m)+1 to
  ! jMaxHf on odd fields. Remaining latitudes are null.


  ! Fourier Fields
  !
  ! Array Four(dvdlj, dip1), using (nVertGrid*myJMax_f, mMax+1)
  ! dvdlj and di avoid memory bank conflict. nVertGrid is the
  ! total number of verticals stored (sum over all fields,
  ! including lambda derivatives in Spectral to Grids transforms)
  !
  ! first dimension of Fou maps verticals of all fields 
  ! and latitudes stored by this process into an 1D structure 
  ! in array element order of:(Vertical, Field, latitude).
  ! 
  ! Second dimension of Fou contains real and imaginary
  ! fourier coefficients, in this order. A surplus value is
  ! required for the FFT (null values).

  ! Internal representation of Grid Fields
  !
  ! Array Grid(dvdlj, di), using Grid(nVertGrid*myJMax_f, iMax)
  ! dvdlj and di avoid memory bank conflict. nVertGrid is the
  ! total number of verticals stored (sum over all fields,
  ! including lambda derivatives in Spectral to Grids transforms)
  !
  ! first dimension of Grid maps verticals of all fields 
  ! and latitudes stored by this process into an 1D structure 
  ! in array element order of:(Vertical, Field, latitude).
  ! 
  ! Second dimension of Grid is longitude.

  INTEGER              :: mGlob ! global variables for OpenMp
  INTEGER              :: nGlob  
  INTEGER              :: iGlob 
  INTEGER              :: jGlob 
  INTEGER              :: ibGlob 
  INTEGER              :: ksg
  INTEGER              :: kountg
  INTEGER              :: ipar2g
  INTEGER              :: ipar3g

  INTERFACE DepositSpecToGrid
     MODULE PROCEDURE Deposit1D, Deposit2D
  END INTERFACE

  INTERFACE DepositSpecToDelLamGrid
     MODULE PROCEDURE DepDL1D, DepDL2D
  END INTERFACE

  INTERFACE DepositSpecToGridAndDelLamGrid
     MODULE PROCEDURE DepDLG1D, DepDLG2D
  END INTERFACE

  INTERFACE DestroySpecToGrid
     MODULE PROCEDURE Destroy
  END INTERFACE

  INTERFACE DepositGridToSpec
     MODULE PROCEDURE Deposit1D, Deposit2D
  END INTERFACE

  INTERFACE DepositGridToSpec_PK
     MODULE PROCEDURE Deposit1D_pk, Deposit2D_PK
  END INTERFACE

  INTERFACE DestroyGridToSpec
     MODULE PROCEDURE Destroy
  END INTERFACE

  REAL(KIND=r8), ALLOCATABLE :: Spec(:,:)
  REAL(KIND=r8), ALLOCATABLE :: Four(:,:)
  INTEGER,       ALLOCATABLE :: mnodes(:)
  INTEGER,       ALLOCATABLE :: requests(:)
  INTEGER,       ALLOCATABLE :: requestr(:)
  INTEGER,       ALLOCATABLE :: status(:)
  INTEGER,       ALLOCATABLE :: stat(:,:)

  TYPE p1d
     REAL(KIND=r8), POINTER :: p(:)
  END TYPE p1d

  TYPE p2d
     REAL(KIND=r8), POINTER :: p(:,:)
  END TYPE p2d

  TYPE p3d
     REAL(KIND=r8), POINTER :: p(:,:,:)
  END TYPE p3d

  INTEGER                :: nSpecFields        ! total # Spectral Fields in transform
  INTEGER                :: usedSpecFields     ! # Spectral Fields currently deposited
  INTEGER                :: lastUsedSpecVert   ! at internal array Spec
  LOGICAL,   ALLOCATABLE :: surfSpec(:)        ! TRUE iff Surface Spectral Field
  INTEGER,   ALLOCATABLE :: prevSpec(:)        ! prior to first real vertical of this field at all internal arrays
  TYPE(p1d), ALLOCATABLE :: Spec1d(:)          ! points to Surface Spectral Field
  TYPE(p2d), ALLOCATABLE :: Spec2d(:)          ! points to Full Spectral Field

  INTEGER                :: nGridFields        ! total # Grid Fields in transform
  INTEGER                :: usedGridFields     ! # Grid Fields currently deposited
  INTEGER                :: lastUsedGridVert   ! at internal array Grid
  LOGICAL,   ALLOCATABLE :: surfGrid(:)        ! TRUE iff Surface Grid Field
  INTEGER,   ALLOCATABLE :: prevGrid(:)        ! prior to first real vertical of this field at all internal arrays
  TYPE(p2d), ALLOCATABLE :: Grid2d(:)          ! points to Surface Grid Field
  TYPE(p3d), ALLOCATABLE :: Grid3d(:)          ! points to Full Grid Field
  LOGICAL,   ALLOCATABLE :: fieldForDelLam(:)  ! TRUE iff this field position stores Lambda Derivative

  LOGICAL                :: willDelLam
  INTEGER                :: usedDelLamFields   ! last position at Grid array used for Lambda Der. fields
  INTEGER                :: lastUsedDelLamVert ! at Grid array
  INTEGER,   ALLOCATABLE :: prevVertDelLamSource(:)  ! source of Lambda Derivative


  INTEGER              :: dlmn
  INTEGER, ALLOCATABLE :: nEven(:)
  INTEGER, ALLOCATABLE :: dnEven(:)
  INTEGER, ALLOCATABLE :: firstNEven(:)
  INTEGER, ALLOCATABLE :: nOdd(:)
  INTEGER, ALLOCATABLE :: dnOdd(:)
  INTEGER, ALLOCATABLE :: firstNOdd(:)
  INTEGER, ALLOCATABLE :: lmnExtMap(:)
  INTEGER, ALLOCATABLE :: lmnMap(:)
  INTEGER, ALLOCATABLE :: lmnZero(:)


  INTEGER              :: nVertSpec
  INTEGER              :: nVertGrid
  INTEGER              :: nFull_g  
  INTEGER              :: nSurf_g  
  INTEGER              :: nFull_s  
  INTEGER              :: nSurf_s  
  INTEGER              :: dv
  INTEGER              :: dvdlj
  INTEGER              :: di
  INTEGER              :: dip1
  INTEGER, ALLOCATABLE :: previousJ(:)

  INTEGER              :: djh
  INTEGER              :: dvjh

  !  HIDDEN DATA, FFT SIZE INDEPENDENT:
  !  nBase=SIZE(Base)
  !  Base are the bases for factorization of n; base 4 should come
  !       before base 2 to improve FFT efficiency
  !  Permutation defines order of bases when computing FFTs
  !  sinXX, cosYY are constants required for computing FFTs

  INTEGER, PARAMETER   :: nBase=4
  INTEGER, PARAMETER   :: Base(nBase) = (/ 4, 2, 3, 5 /)
  INTEGER, PARAMETER   :: Permutation(nBase) = (/ 2, 3, 1, 4 /)
  REAL(KIND=r8)                 :: sin60
  REAL(KIND=r8)                 :: sin36
  REAL(KIND=r8)                 :: sin72
  REAL(KIND=r8)                 :: cos36
  REAL(KIND=r8)                 :: cos72

  !  HIDDEN DATA, FFT SIZE DEPENDENT:
  !  For each FFT Size:
  !     nLong is FFT size;
  !     nFactors, Factors(:), nTrig, Trig(:) are size dependent constants;
  !     firstLat is the first FFT latitude for this block
  !     lastLat is the last FFT latitude for this block
  !  BlockFFT has size dependent data for all required FFTs;
  !  nBlockFFT is BlockFFT size

  TYPE MultiFFT
     INTEGER           :: nLong
     INTEGER           :: nFactors
     INTEGER, POINTER  :: Factors(:)
     INTEGER           :: nTrigs
     REAL(KIND=r8),    POINTER  :: Trigs(:)
     INTEGER           :: firstLat
     INTEGER           :: lastLat
  END TYPE MultiFFT

  TYPE(MultiFFT), ALLOCATABLE, TARGET :: BlockFFT(:)
  INTEGER :: nBlockFFT

  REAL(KIND=r8),    ALLOCATABLE :: LS2F(:,:)
  REAL(KIND=r8),    ALLOCATABLE :: LF2S(:,:)
  REAL(KIND=r8),    ALLOCATABLE :: consIm(:)
  REAL(KIND=r8),    ALLOCATABLE :: consRe(:)

  LOGICAL, PARAMETER :: dumpLocal=.FALSE.


  ! CALLING TREE (internal)
  !
  ! InitTransform ! MpiMappings
  !               ! CreateFFT    ! Factorize
  !                              ! TrigFactors
  ! CreateSpecToGrid 
  !
  ! DepositSpecToGrid ! Deposit1D or Deposit2D
  !
  ! DepositSpecToDelLamGrid ! Deposit1D or Deposit2D
  !
  ! DepositSpecToGridAndDelLamGrid ! Deposit1D or Deposit2D
  !
  ! DoSpecToGrid ! DepositSpec
  !              ! SpecToFour
  !              ! DelLam
  !              ! InvFFTTrans
  !              ! WithdrawGrid
  !
  ! DestroySpecToGrid ! Destroy
  !
  ! CreateGridToSpec
  !
  ! DepositGridToSpec ! Deposit1D or Deposit2D
  !
  ! DoGridToSpec ! DepositGrid
  !              ! DirFFTTrans
  !              ! FourToSpec
  !              ! WithdrawSpectral
  !
  ! DestroyGridToSpec ! Destroy
  !
  !






CONTAINS



  SUBROUTINE InitTransform()
    INTEGER :: m
    INTEGER :: n
    INTEGER :: lm
    INTEGER :: lmn
    INTEGER :: i
    INTEGER :: j
    REAL(KIND=r8)    :: radi
    INTEGER :: iBlockFFT
    INTEGER :: lFirst, lLast
    TYPE(MultiFFT), POINTER :: p
    CHARACTER(LEN=10) :: c0, c1, c2, c3
    CHARACTER(LEN=*), PARAMETER :: h="**(InitTransform)**"

    ! Mapping (lm,n) into array Spec, for transform;
    ! for each lm:
    !    map all even components of n, 
    !    followed by all odd components of n.
    ! Limits of even and odd components at array Spec:
    ! nEven(myMMax): How many n with (n+m) even for this lm
    ! nOdd (myMMax): How many n with (n+m) odd  for this lm

    ALLOCATE(nEven(myMMax))
    ALLOCATE(nOdd(myMMax))
    DO lm = 1, myMMax
       m = lm2m(lm)
       nEven(lm) = (mMax - m + 3)/2
       nOdd(lm)  = (mMax - m + 2)/2
    END DO

    ! Spec Bounds:
    ! dnEven(myMMax): nEven without memory bank conflicts
    ! dnOdd (myMMax): nOdd  without memory bank conflicts

    ALLOCATE(dnEven(myMMax))
    ALLOCATE(dnOdd(myMMax))
    dnEven = NoBankConflict(nEven)
    dnOdd  = NoBankConflict(nOdd)

    ! for each lm, find:
    ! first n even, first n odd;
    ! all n even will be mapped into sequential positions,
    ! from firstNEven(lm)
    ! all n odd will be mapped into sequential positions,
    ! from firstNOdd(lm)

    ALLOCATE(firstNEven(myMMax))
    ALLOCATE(firstNOdd(myMMax))
    firstNEven(1) = 1
    firstNOdd(1)  = nEven(1) + 1
    DO lm = 2, myMMax
       firstNEven(lm) = firstNOdd(lm-1) + nOdd(lm-1)
       firstNOdd (lm) = firstNEven(lm)  + nEven(lm)
    END DO

    ! size of first dimension of array Spec:
    ! space for all even and odd values, avoiding bank conflict;
    ! garantees space for each lm matrix;

    dlmn = SUM(nEven(1:myMMax)) + SUM(nOdd(1:myMMax-1)) + dnOdd(myMMax)
    dlmn = NoBankConflict(dlmn)

    ! lmnExtMap: maps external representation of Extended Spectral fields
    !            into internal representation for transforms (array Spec)

    ALLOCATE(lmnExtMap(myMNExtMax))
    DO lmn = 1, myMNExtMax
       lm = myMExtMap(lmn)
       n  = myNExtMap(lmn)
       m  = lm2m(lm)
       IF (MOD(m+n,2) == 0) THEN
          lmnExtMap(lmn) = firstNEven(lm) + (n - m)/2
       ELSE
          lmnExtMap(lmn) = firstNOdd(lm) + (n - m - 1)/2
       END IF
    END DO

    ! lmnExtMap: maps external representation of Regular Spectral fields
    !            into internal representation for transforms (array Spec)

    ALLOCATE(lmnMap(myMNMax))
    DO lmn = 1, myMNMax
       lm = myMMap(lmn)
       n  = myNMap(lmn)
       m  = lm2m(lm)
       IF (MOD(m+n,2) == 0) THEN
          lmnMap(lmn) = firstNEven(lm) + (n - m)/2
       ELSE
          lmnMap(lmn) = firstNOdd(lm) + (n - m - 1)/2
       END IF
    END DO

    ! lmnZero: null positions for Regular Spectral fields into
    !          internal representations for transforms

    ALLOCATE(lmnZero(myMMax))
    DO lm = 1, myMMax
       lmnZero(lm) = lmnExtMap(myMNExtMap(lm,mMax+1))
    END DO

    ! dimensions

    di = NoBankConflict(iMax)
    dip1 = NoBankConflict(iMax+1)
    djh = NoBankConflict(jMaxHalf)


    ! previousJ: how many latitudes are stored at array Four
    !            before latitude j

    ALLOCATE(previousJ(jMax))
    DO j = myfirstlat_f,mylastlat_f
       previousJ(j       ) = j-myfirstlat_f
    END DO

    ! consIm, consRe : values to be used in dellam

    ALLOCATE(consIm(mMax))
    ALLOCATE(consRe(mMax))
    DO m = 1,mMax
       consIm(m) = EMRad1*REAL(m-1,r8)
       consRe(m) = - consIm(m)
    END DO
    ! LS2F: Associated Legendre Functions (latitude, mn)
    !       for spectral to fourier computation

    ALLOCATE(LS2F(djh, dlmn))
    IF (allpolynomials) THEN
       DO lmn = 1, myMNExtMax
          lm = myMExtMap(lmn)
          n  = myNExtMap(lmn)
          m  = lm2m(lm)
          DO j = 1, jMaxHalf
             LS2F(j,lmnExtMap(lmn)) = LegFuncS2F(j,mnExtMap(m,n))
          END DO
       END DO
     ELSE
       DO lmn = 1, myMNExtMax
          m = myMExtMap(lmn)
          n  = myNExtMap(lmn)
          DO j = 1, jMaxHalf
             LS2F(j,lmnExtMap(lmn)) = LegFuncS2F(j,mymnExtMap(m,n))
          END DO
       END DO
    END IF
    LS2F(jMaxHalf+1:djh,                 :) = 0.0_r8
    LS2F(:,              myMNExtMax+1:dlmn) = 0.0_r8

    ! LF2S: Associated Legendre Functions (latitude, mn)
    !       for fourier to spectral computation

    ALLOCATE(LF2S(dlmn, djh))
    DO j = 1, jMaxHalf
       DO lmn = 1, myMNExtMax
          LF2S(lmnExtMap(lmn),j) =  LS2F(j,lmnExtMap(lmn)) * GaussWeights(j)
       END DO
    END DO
    LF2S(:,                 jMaxHalf+1:djh) = 0.0_r8
    LF2S(myMNExtMax+1:dlmn,           :   ) = 0.0_r8

    ! FFT size independent data

    radi=ATAN(1.0_r8)/45.0_r8
    sin60=SIN(60.0_r8*radi)
    sin36=SIN(36.0_r8*radi)
    sin72=SIN(72.0_r8*radi)
    cos36=COS(36.0_r8*radi)
    cos72=COS(72.0_r8*radi)

    ! count FFT Blocks (one size FFT per block)

    nBlockFFT = 1
    DO i = myfirstlat_f+1,mylastlat_f
       IF (iMaxPerJ(i).ne.iMaxPerJ(i-1)) THEN
          nBlockFFT = nBlockFFT+1
       END IF
    END DO

    ! Block FFT size 
    ! Block FFT first Lat position on array Four
    ! Block FFT last Lat position on array Four

    ALLOCATE (BlockFFT(nBlockFFT))
    iBlockFFT = 1
    p => BlockFFT(iBlockFFT)
    p%nLong = iMaxPerJ(myfirstlat_f)
    p%firstLat = 1
    DO j = myfirstlat_f+1,mylastlat_f
       IF (iMaxPerJ(j) /= p%nLong) THEN
          p%lastLat = j - myfirstlat_f
          iBlockFFT = iBlockFFT + 1
          p => BlockFFT(iBlockFFT)
          p%nLong = iMaxPerJ(j)
          p%firstLat = j+1-myfirstlat_f
       END IF
    END DO
    IF (iBlockFFT /= nBlockFFT) THEN
       WRITE(nfprt,"(a,' iBlockFFT (',i5,') /= nBlockFFT (',i5,')')") &
            h, iBlockFFT, nBLockFFT
       STOP h
    ELSE
       p%lastLat = mylastlat_f+1-myfirstlat_f
    END IF

    ! fill internal FFT data for each FFT Block

    DO iBlockFFT = 1, nBlockFFT
       p => BlockFFT(iBlockFFT)
       CALL CreateFFT(p%nLong, p%Factors, p%nFactors, p%Trigs, p%nTrigs)
    END DO

    ! debug dumping

    IF (dumpLocal) THEN
       WRITE(nfprt,"(a,' starts dumping internal data')") h
       WRITE(nfprt,"(a,' dlmn      =',i10)") h, dlmn
       WRITE(nfprt,"(a,' di        =',i10)") h, di
       WRITE(nfprt,"(a,' dip1      =',i10)") h, dip1
       WRITE(nfprt,"(a,' djh       =',i10)") h, djh
       WRITE(c0,"(i10)") nBlockFFT
       WRITE(nfprt,"(a,' There are ',a,' FFT Blocks')") h, TRIM(ADJUSTL(c0))
       WRITE(nfprt,"(a,' Block #; FFT Size(longitudes); latitudes')") h
       DO iBlockFFT = 1, nBlockFFT
          p => BlockFFT(iBlockFFT)
          lFirst = p%firstLat
          lLast = p%lastLat
          IF (lLast == lFirst + 1) THEN
             WRITE(c0,"(i10)") (lFirst+1)/2
             WRITE(c1,"(i10)") jMax - (lLast+1)/2 + 1
             WRITE(nfprt,"(a,i8,i20,5x,a,' and ',a)") h, iBlockFFT, p%nLong, &
                  TRIM(ADJUSTL(c0)), TRIM(ADJUSTL(c1))
          ELSE
             WRITE(c0,"(i10)") (lFirst+1)/2
             WRITE(c1,"(i10)") (lLast+1)/2
             WRITE(c2,"(i10)") jMax - (lFirst+1)/2 + 1
             WRITE(c3,"(i10)") jMax - (lLast+1)/2 + 1
             WRITE(nfprt,"(a,i8,i20,5x,a,':',a,' and ',a,':',a)") h, iBlockFFT, p%nLong, &
                  TRIM(ADJUSTL(c0)), TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c3)), TRIM(ADJUSTL(c2))
          END IF
       END DO
    END IF
  END SUBROUTINE InitTransform





  SUBROUTINE CreateSpecToGrid(nFullSpec, nSurfSpec, nFullGrid, nSurfGrid)
    INTEGER, INTENT(IN) :: nFullSpec
    INTEGER, INTENT(IN) :: nSurfSpec
    INTEGER, INTENT(IN) :: nFullGrid
    INTEGER, INTENT(IN) :: nSurfGrid
    CHARACTER(LEN=*), PARAMETER :: h="**(CreateSpecToGrid)**"
    INTEGER :: nsusp, nsugr

    IF (havesurf) THEN
       nsusp = nSurfSpec
       nsugr = nSurfGrid
      ELSE
       nsusp = 0
       nsugr = 0
    ENDIF
    nfull_g = nFullGrid
    nsurf_g = nSurfGrid
    nfull_s = nFullSpec
    nsurf_s = nSurfSpec
    nSpecFields = nFullSpec + nsusp
    usedSpecFields = 0
    lastUsedSpecVert = 0
    ALLOCATE(surfSpec(nSpecFields))
    ALLOCATE(prevSpec(nSpecFields))
    ALLOCATE(Spec1d(nSpecFields))
    ALLOCATE(Spec2d(nSpecFields))
    nVertSpec = nFullSpec*kMaxloc + nsusp

    nGridFields = nFullGrid + nSurfGrid
    usedGridFields = 0
    lastUsedGridVert = 0
    ALLOCATE(surfGrid(nGridFields))
    ALLOCATE(prevGrid(nGridFields))
    ALLOCATE(Grid2d(nGridFields))
    ALLOCATE(Grid3d(nGridFields))
    ALLOCATE(fieldForDelLam(nGridFields))
    nVertGrid = nFullGrid*kMaxloc + nsugr

    willDelLam = .FALSE.
    usedDelLamFields = nFullSpec + nSurfSpec
    lastUsedDelLamVert = nVertSpec
    ALLOCATE(prevVertDelLamSource(nGridFields))

    dv = NoBankConflict(2*nVertSpec)
    dvjh = NoBankConflict(nVertSpec*jMaxHalf)
    dvdlj = NoBankConflict(nVertGrid*myJMax_f)

    ALLOCATE (Spec(dlmn, dv))
    ALLOCATE (Four(dvdlj, dip1))
    IF (.NOT.ALLOCATED(requests)) THEN
       ALLOCATE(requests(0:maxnodes))
       ALLOCATE(requestr(0:maxnodes))
       ALLOCATE(status(MPI_STATUS_SIZE))
       ALLOCATE(stat(MPI_STATUS_SIZE,maxnodes))
    END IF

    IF (dumpLocal) THEN
       WRITE(nfprt,"(a,' Spec: n, used, lastVert, nVert=',4i5)") &
            h, nSpecFields, usedSpecFields, lastUsedSpecVert, nVertSpec
       WRITE(nfprt,"(a,' Grid: n, used, lastVert, nVert=',4i5)") &
            h, nGridFields, usedGridFields, lastUsedGridVert, nVertGrid
       WRITE(nfprt,"(a,' DelLam: used, lastVert=',2i5)") &
            h, usedDelLamFields, lastUsedDelLamVert
    END IF
  END SUBROUTINE CreateSpecToGrid


  SUBROUTINE Deposit1D_PK (ArgSpec, ArgGrid, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgGrid(:,:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgDelLam(:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(Deposit1D)**"
    LOGICAL :: getGrid, getDelLam
    INTEGER :: lusf
    INTEGER :: pusv
    INTEGER :: lusv

    INTEGER :: lugf
    INTEGER :: pugv
    INTEGER :: lugv

    INTEGER :: ludlf
    INTEGER :: pudlv
    INTEGER :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid   = PRESENT(ArgGrid)
    getDelLam = PRESENT(ArgDelLam)

    ! critical section to update data structure counters 

    IF (.not.havesurf) THEN

    !
    !  this processor won't transform surface fields
    !  ( it only needs to know which are the correponding grid fields )
    !
       IF (getGrid) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          surfGrid(lugf) = .TRUE.
          prevGrid(lugf) = -1
          NULLIFY(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgGrid
          fieldForDelLam(lugf) = .FALSE.
          prevVertDelLamSource(lugf) = -1
          IF (getDelLam) THEN
             ludlf = usedDelLamFields + 1
             usedDelLamFields = ludlf
             surfGrid(ludlf) = .TRUE.
             prevGrid(ludlf) = -1
             NULLIFY(Grid3d(ludlf)%p)
             Grid2d(ludlf)%p => ArgDelLam
             fieldForDelLam(lugf) = .FALSE.
             prevVertDelLamSource(lugf) = -1
          END IF
       ELSE IF (getDelLam) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          surfGrid(lugf) = .TRUE.
          prevGrid(lugf) = -1
          NULLIFY(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgDelLam
          fieldForDelLam(lugf) = .FALSE. 
          prevVertDelLamSource(lugf) = - 1
       ELSE
          WRITE(nfprt,"(a, ' no Grid or DelLam field required')") h
          STOP h
       END IF

    ELSE 

       lusf = usedSpecFields + 1
       usedSpecFields = lusf
       pusv = lastUsedSpecVert 
       lusv = pusv + 1
       lastUsedSpecVert = lusv
       IF (getGrid) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          pugv = lastUsedGridVert
          lugv = pugv + 1
          lastUsedGridVert = lugv
          IF (getDelLam) THEN
             ludlf = usedDelLamFields + 1
             usedDelLamFields = ludlf
             pudlv = lastUsedDelLamVert 
             ludlv = pudlv + 1
             lastUsedDelLamVert = ludlv
          END IF
       ELSE IF (getDelLam) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          pugv = lastUsedGridVert
          lugv = pugv + 1
          lastUsedGridVert = lugv
       ELSE
          WRITE(nfprt,"(a, ' no Grid or DelLam field required')") h
          STOP h
       END IF
   
       ! deposit spectral field

       IF (lusf <= nSpecFields) THEN
          surfSpec(lusf) = .TRUE.
          prevSpec(lusf) = pusv
          NULLIFY(Spec2d(lusf)%p)
          Spec1d(lusf)%p => ArgSpec
       ELSE
          WRITE(nfprt,"(a, ' too many spectral fields')") h
          WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
          STOP h
       END IF

       IF (dumpLocal) THEN
          WRITE(nfprt,"(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
               h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
       END IF

       IF (getGrid) THEN
   
          ! ArgGrid Field is required; 
          ! Store ArgGrid Field info at first available position of the Grid array

          IF (lugf <= nSpecFields) THEN
             surfGrid(lugf) = .TRUE.
             prevGrid(lugf) = pugv
             NULLIFY(Grid3d(lugf)%p)
             Grid2d(lugf)%p => ArgGrid
             fieldForDelLam(lugf) = .FALSE.
             prevVertDelLamSource(lugf) = -1
          ELSE
             WRITE(nfprt,"(a, ' too many grid fields')") h
             WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
             STOP h
          END IF

          IF (dumpLocal) THEN
             WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                  h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
             WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                  h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
          END IF

          IF (getDelLam) THEN
   
             ! ArgGrid and ArgDelLam are required; 
             ! Store ArgDelLam info at first available position of the Grid array
             ! Mark field for DelLam

             willDelLam = .TRUE.
             IF (ludlf <= nGridFields) THEN
                surfGrid(ludlf) = .TRUE.
                prevGrid(ludlf) = pudlv
                NULLIFY(Grid3d(ludlf)%p)
                Grid2d(ludlf)%p => ArgDelLam
                fieldForDelLam(ludlf) = .TRUE.
                prevVertDelLamSource(ludlf) = prevGrid(lugf)
             ELSE
                WRITE(nfprt,"(a, ' too many grid fields including DelLam')") h
                WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
                STOP h
             END IF

             IF (dumpLocal) THEN
                WRITE(nfprt,"(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                     h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
                WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                     h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
             END IF
          END IF

       ELSE IF (getDelLam) THEN

          ! ArgDelLam Field is required; Grid is not
          ! Store ArgDelLam info at first available position of the Grid array
          ! Mark field for DelLam

          willDelLam = .TRUE.
          IF (lugf <= nSpecFields) THEN
             surfGrid(lugf) = .TRUE.
             prevGrid(lugf) = pugv
             NULLIFY(Grid3d(lugf)%p)
             Grid2d(lugf)%p => ArgDelLam
             fieldForDelLam(lugf) = .TRUE.
             prevVertDelLamSource(lugf) = prevGrid(lugf)
          ELSE
             WRITE(nfprt,"(a, ' too many grid fields')") h
             WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
             STOP h
          END IF

          IF (dumpLocal) THEN
             WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                  h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
             WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                  h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
          END IF
       END IF

    END IF 


  END SUBROUTINE Deposit1D_PK

  SUBROUTINE Deposit1D (ArgSpec, ArgGrid, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgGrid(:,:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgDelLam(:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(Deposit1D)**"
    LOGICAL :: getGrid, getDelLam
    INTEGER :: lusf
    INTEGER :: pusv
    INTEGER :: lusv

    INTEGER :: lugf
    INTEGER :: pugv
    INTEGER :: lugv

    INTEGER :: ludlf
    INTEGER :: pudlv
    INTEGER :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid   = PRESENT(ArgGrid)
    getDelLam = PRESENT(ArgDelLam)

    ! critical section to update data structure counters 

    IF (.not.havesurf) THEN

    !
    !  this processor won't transform surface fields
    !  ( it only needs to know which are the correponding grid fields )
    !
       IF (getGrid) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          surfGrid(lugf) = .TRUE.
          prevGrid(lugf) = -1
          NULLIFY(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgGrid
          fieldForDelLam(lugf) = .FALSE.
          prevVertDelLamSource(lugf) = -1
          IF (getDelLam) THEN
             ludlf = usedDelLamFields + 1
             usedDelLamFields = ludlf
             surfGrid(ludlf) = .TRUE.
             prevGrid(ludlf) = -1
             NULLIFY(Grid3d(ludlf)%p)
             Grid2d(ludlf)%p => ArgDelLam
             fieldForDelLam(lugf) = .FALSE.
             prevVertDelLamSource(lugf) = -1
          END IF
       ELSE IF (getDelLam) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          surfGrid(lugf) = .TRUE.
          prevGrid(lugf) = -1
          NULLIFY(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgDelLam
          fieldForDelLam(lugf) = .FALSE. 
          prevVertDelLamSource(lugf) = - 1
       ELSE
          WRITE(nfprt,"(a, ' no Grid or DelLam field required')") h
          STOP h
       END IF

    ELSE 

       lusf = usedSpecFields + 1
       usedSpecFields = lusf
       pusv = lastUsedSpecVert 
       lusv = pusv + 1
       lastUsedSpecVert = lusv
       IF (getGrid) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          pugv = lastUsedGridVert
          lugv = pugv + 1
          lastUsedGridVert = lugv
          IF (getDelLam) THEN
             ludlf = usedDelLamFields + 1
             usedDelLamFields = ludlf
             pudlv = lastUsedDelLamVert 
             ludlv = pudlv + 1
             lastUsedDelLamVert = ludlv
          END IF
       ELSE IF (getDelLam) THEN
          lugf = usedGridFields+1
          usedGridFields = lugf
          pugv = lastUsedGridVert
          lugv = pugv + 1
          lastUsedGridVert = lugv
       ELSE
          WRITE(nfprt,"(a, ' no Grid or DelLam field required')") h
          STOP h
       END IF
   
       ! deposit spectral field

       IF (lusf <= nSpecFields) THEN
          surfSpec(lusf) = .TRUE.
          prevSpec(lusf) = pusv
          NULLIFY(Spec2d(lusf)%p)
          Spec1d(lusf)%p => ArgSpec
       ELSE
          WRITE(nfprt,"(a, ' too many spectral fields')") h
          WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
          STOP h
       END IF

       IF (dumpLocal) THEN
          WRITE(nfprt,"(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
               h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
       END IF

       IF (getGrid) THEN
   
          ! ArgGrid Field is required; 
          ! Store ArgGrid Field info at first available position of the Grid array

          IF (lugf <= nSpecFields) THEN
             surfGrid(lugf) = .TRUE.
             prevGrid(lugf) = pugv
             NULLIFY(Grid3d(lugf)%p)
             Grid2d(lugf)%p => ArgGrid
             fieldForDelLam(lugf) = .FALSE.
             prevVertDelLamSource(lugf) = -1
          ELSE
             WRITE(nfprt,"(a, ' too many grid fields')") h
             WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
             STOP h
          END IF

          IF (dumpLocal) THEN
             WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                  h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
             WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                  h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
          END IF

          IF (getDelLam) THEN
   
             ! ArgGrid and ArgDelLam are required; 
             ! Store ArgDelLam info at first available position of the Grid array
             ! Mark field for DelLam

             willDelLam = .TRUE.
             IF (ludlf <= nGridFields) THEN
                surfGrid(ludlf) = .TRUE.
                prevGrid(ludlf) = pudlv
                NULLIFY(Grid3d(ludlf)%p)
                Grid2d(ludlf)%p => ArgDelLam
                fieldForDelLam(ludlf) = .TRUE.
                prevVertDelLamSource(ludlf) = prevGrid(lugf)
             ELSE
                WRITE(nfprt,"(a, ' too many grid fields including DelLam')") h
                WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
                STOP h
             END IF

             IF (dumpLocal) THEN
                WRITE(nfprt,"(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                     h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
                WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                     h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
             END IF
          END IF

       ELSE IF (getDelLam) THEN

          ! ArgDelLam Field is required; Grid is not
          ! Store ArgDelLam info at first available position of the Grid array
          ! Mark field for DelLam

          willDelLam = .TRUE.
          IF (lugf <= nSpecFields) THEN
             surfGrid(lugf) = .TRUE.
             prevGrid(lugf) = pugv
             NULLIFY(Grid3d(lugf)%p)
             Grid2d(lugf)%p => ArgDelLam
             fieldForDelLam(lugf) = .TRUE.
             prevVertDelLamSource(lugf) = prevGrid(lugf)
          ELSE
             WRITE(nfprt,"(a, ' too many grid fields')") h
             WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
             STOP h
          END IF

          IF (dumpLocal) THEN
             WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                  h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
             WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                  h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
          END IF
       END IF

    END IF 


  END SUBROUTINE Deposit1D



  SUBROUTINE Deposit2D (ArgSpec, ArgGrid, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:,:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgGrid(:,:,:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgDelLam(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(Deposit2D)**"
    LOGICAL :: getGrid, getDelLam
    INTEGER :: lusf
    INTEGER :: pusv
    INTEGER :: lusv

    INTEGER :: lugf
    INTEGER :: pugv
    INTEGER :: lugv

    INTEGER :: ludlf
    INTEGER :: pudlv
    INTEGER :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid   = PRESENT(ArgGrid)
    getDelLam = PRESENT(ArgDelLam)

    ! critical section to update data structure counters 

    lusf = usedSpecFields + 1
    usedSpecFields = lusf
    pusv = lastUsedSpecVert 
    lusv = pusv + kMaxloc
    lastUsedSpecVert = lusv
    IF (getGrid) THEN
       lugf = usedGridFields+1
       usedGridFields = lugf
       pugv = lastUsedGridVert
       lugv = pugv + kMaxloc
       lastUsedGridVert = lugv
       IF (getDelLam) THEN
          ludlf = usedDelLamFields + 1
          usedDelLamFields = ludlf
          pudlv = lastUsedDelLamVert 
          ludlv = pudlv + kMaxloc
          lastUsedDelLamVert = ludlv
       END IF
    ELSE IF (getDelLam) THEN
       lugf = usedGridFields+1
       usedGridFields = lugf
       pugv = lastUsedGridVert
       lugv = pugv + kMaxloc
       lastUsedGridVert = lugv
    ELSE
       WRITE(nfprt,"(a, ' no Grid or DelLam field required')") h
       STOP h
    END IF

    ! deposit spectral field

    IF (lusf <= nSpecFields) THEN
       surfSpec(lusf) = .FALSE.
       prevSpec(lusf) = pusv
       NULLIFY(Spec1d(lusf)%p)
       Spec2d(lusf)%p => ArgSpec
    ELSE
       WRITE(nfprt,"(a, ' too many spectral fields')") h
       WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
       STOP h
    END IF

    IF (dumpLocal) THEN
       WRITE(nfprt,"(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
            h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
    END IF

    IF (getGrid) THEN

       ! ArgGrid Field is required; 
       ! Store ArgGrid Field info at first available position of the Grid array

       IF (lugf <= nfull_s+nsurf_s) THEN
          surfGrid(lugf) = .FALSE.
          prevGrid(lugf) = pugv
          NULLIFY(Grid2d(lugf)%p)
          Grid3d(lugf)%p => ArgGrid
          fieldForDelLam(lugf) = .FALSE.
          prevVertDelLamSource(lugf) = -1
       ELSE
          WRITE(nfprt,"(a, ' too many grid fields')") h
          WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          STOP h
       END IF

       IF (dumpLocal) THEN
          WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
               h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
               h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
       END IF

       IF (getDelLam) THEN

          ! ArgGrid and ArgDelLam are required; 
          ! Store ArgDelLam info at first available position of the Grid array
          ! Mark field for DelLam

          willDelLam = .TRUE.
          IF (ludlf <= nGridFields) THEN
             surfGrid(ludlf) = .FALSE.
             prevGrid(ludlf) = pudlv
             NULLIFY(Grid2d(ludlf)%p)
             Grid3d(ludlf)%p => ArgDelLam
             fieldForDelLam(ludlf) = .TRUE.
             prevVertDelLamSource(ludlf) = prevGrid(lugf)
          ELSE
             WRITE(nfprt,"(a, ' too many grid fields including DelLam')") h
             WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
             STOP h
          END IF

          IF (dumpLocal) THEN
             WRITE(nfprt,"(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                  h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
             WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                  h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
          END IF
       END IF

    ELSE IF (getDelLam) THEN

       ! ArgDelLam Field is required; Grid is not
       ! Store ArgDelLam info at first available position of the Grid array
       ! Mark field for DelLam

       willDelLam = .TRUE.
       IF (lugf <= nfull_s+nsurf_s) THEN
          surfGrid(lugf) = .FALSE.
          prevGrid(lugf) = pugv
          NULLIFY(Grid2d(lugf)%p)
          Grid3d(lugf)%p => ArgDelLam
          fieldForDelLam(lugf) = .TRUE.
          prevVertDelLamSource(lugf) = prevGrid(lugf)
       ELSE
          WRITE(nfprt,"(a, ' too many grid fields')") h
          WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          STOP h
       END IF

       IF (dumpLocal) THEN
          WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
               h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
               h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
       END IF
    END IF
  END SUBROUTINE Deposit2D


  SUBROUTINE Deposit2D_PK (ArgSpec, ArgGrid, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:,:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgGrid(:,:,:)
    REAL(KIND=r8), TARGET, OPTIONAL, INTENT(IN) :: ArgDelLam(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(Deposit2D)**"
    LOGICAL :: getGrid, getDelLam
    INTEGER :: lusf
    INTEGER :: pusv
    INTEGER :: lusv

    INTEGER :: lugf
    INTEGER :: pugv
    INTEGER :: lugv

    INTEGER :: ludlf
    INTEGER :: pudlv
    INTEGER :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid   = PRESENT(ArgGrid)
    getDelLam = PRESENT(ArgDelLam)

    ! critical section to update data structure counters 

    lusf = usedSpecFields + 1
    usedSpecFields = lusf
    pusv = lastUsedSpecVert 
    lusv = pusv + kMaxloc
    lastUsedSpecVert = lusv
    IF (getGrid) THEN
       lugf = usedGridFields+1
       usedGridFields = lugf
       pugv = lastUsedGridVert
       lugv = pugv + kMaxloc
       lastUsedGridVert = lugv
       IF (getDelLam) THEN
          ludlf = usedDelLamFields + 1
          usedDelLamFields = ludlf
          pudlv = lastUsedDelLamVert 
          ludlv = pudlv + kMaxloc
          lastUsedDelLamVert = ludlv
       END IF
    ELSE IF (getDelLam) THEN
       lugf = usedGridFields+1
       usedGridFields = lugf
       pugv = lastUsedGridVert
       lugv = pugv + kMaxloc
       lastUsedGridVert = lugv
    ELSE
       WRITE(nfprt,"(a, ' no Grid or DelLam field required')") h
       STOP h
    END IF

    ! deposit spectral field

    IF (lusf <= nSpecFields) THEN
       surfSpec(lusf) = .FALSE.
       prevSpec(lusf) = pusv
       NULLIFY(Spec1d(lusf)%p)
       Spec2d(lusf)%p => ArgSpec
    ELSE
       WRITE(nfprt,"(a, ' too many spectral fields')") h
       WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
       STOP h
    END IF

    IF (dumpLocal) THEN
       WRITE(nfprt,"(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
            h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
    END IF

    IF (getGrid) THEN

       ! ArgGrid Field is required; 
       ! Store ArgGrid Field info at first available position of the Grid array

       IF (lugf <= nfull_g+nsurf_g) THEN
          surfGrid(lugf) = .FALSE.
          prevGrid(lugf) = pugv
          NULLIFY(Grid2d(lugf)%p)
          Grid3d(lugf)%p => ArgGrid
          fieldForDelLam(lugf) = .FALSE.
          prevVertDelLamSource(lugf) = -1
       ELSE
          WRITE(nfprt,"(a, ' too many grid fields')") h
          WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          STOP h
       END IF

       IF (dumpLocal) THEN
          WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
               h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
               h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
       END IF

       IF (getDelLam) THEN

          ! ArgGrid and ArgDelLam are required; 
          ! Store ArgDelLam info at first available position of the Grid array
          ! Mark field for DelLam

          willDelLam = .TRUE.
          IF (ludlf <= nGridFields) THEN
             surfGrid(ludlf) = .FALSE.
             prevGrid(ludlf) = pudlv
             NULLIFY(Grid2d(ludlf)%p)
             Grid3d(ludlf)%p => ArgDelLam
             fieldForDelLam(ludlf) = .TRUE.
             prevVertDelLamSource(ludlf) = prevGrid(lugf)
          ELSE
             WRITE(nfprt,"(a, ' too many grid fields including DelLam')") h
             WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
             STOP h
          END IF

          IF (dumpLocal) THEN
             WRITE(nfprt,"(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
                  h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
             WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
                  h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
          END IF
       END IF

    ELSE IF (getDelLam) THEN

       ! ArgDelLam Field is required; Grid is not
       ! Store ArgDelLam info at first available position of the Grid array
       ! Mark field for DelLam

       willDelLam = .TRUE.
       IF (lugf <= nfull_s+nsurf_s) THEN
          surfGrid(lugf) = .FALSE.
          prevGrid(lugf) = pugv
          NULLIFY(Grid2d(lugf)%p)
          Grid3d(lugf)%p => ArgDelLam
          fieldForDelLam(lugf) = .TRUE.
          prevVertDelLamSource(lugf) = prevGrid(lugf)
       ELSE
          WRITE(nfprt,"(a, ' too many grid fields')") h
          WRITE(nfprt,"(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          STOP h
       END IF

       IF (dumpLocal) THEN
          WRITE(nfprt,"(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
               h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          WRITE(nfprt,"(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
               h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
       END IF
    END IF
  END SUBROUTINE Deposit2D_PK


  SUBROUTINE DepDLG1D (ArgSpec, ArgGrid, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgGrid(:,:)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgDelLam(:,:)
    CALL Deposit1D (ArgSpec, ArgGrid, ArgDelLam)
  END SUBROUTINE DepDLG1D





  SUBROUTINE DepDLG2D (ArgSpec, ArgGrid, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:,:)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgGrid(:,:,:)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgDelLam(:,:,:)
    CALL Deposit2D (ArgSpec, ArgGrid, ArgDelLam)
  END SUBROUTINE DepDLG2D





  SUBROUTINE DepDL1D (ArgSpec, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgDelLam(:,:)
    CALL Deposit1D (ArgSpec, ArgDelLam=ArgDelLam)
  END SUBROUTINE DepDL1D





  SUBROUTINE DepDL2D (ArgSpec, ArgDelLam)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgSpec(:,:)
    REAL(KIND=r8), TARGET, INTENT(IN) :: ArgDelLam(:,:,:)
    CALL Deposit2D (ArgSpec, ArgDelLam=ArgDelLam)
  END SUBROUTINE DepDL2D





  SUBROUTINE DoSpecToGrid()
    CHARACTER(LEN=*), PARAMETER :: h="**(DoSpecToGrid)**"
    INTEGER :: mnFirst
    INTEGER :: mnLast
    INTEGER :: mnExtFirst
    INTEGER :: mnExtLast
    INTEGER :: mFirst
    INTEGER :: mLast
    INTEGER :: jbFirst
    INTEGER :: jbLast
    INTEGER :: iBlockFFT
    INTEGER :: first
    INTEGER :: last
    INTEGER :: FFTfirst 
    INTEGER :: FFTLast
    INTEGER :: FFTSize
    INTEGER :: firstBlock
    INTEGER :: lastBlock
    INTEGER :: sizeBlock
    INTEGER :: i
    INTEGER :: j
    INTEGER :: iFirst
    INTEGER :: iLast
    TYPE(MultiFFT), POINTER :: p


    ! start OMP parallelism


    CALL ThreadDecomp(1, mymnMax, mnFirst, mnLast, "DoSpecToGrid")
    CALL ThreadDecomp(1, mymnExtMax, mnExtFirst, mnExtLast, "DoSpecToGrid")
    CALL ThreadDecomp(1, mymMax, mFirst, mLast, "DoSpecToGrid")
    CALL ThreadDecomp(1, jbMax, jbFirst, jbLast, "DoSpecToGrid")
    CALL ThreadDecomp(1, dip1, iFirst, iLast, "DoSpecToGrid")

    ! nullify Four
    IF (nfull_g.gt.0.or.havesurf) THEN

    DO i = iFirst, iLast
       DO j = 1, dvdlj
          Four(j,i) = 0.0_r8
       END DO
    END DO

    ! ingest spectral fields

    CALL DepositSpec(mnFirst, mnLast, mnExtFirst, mnExtLast, mFirst, mLast)

    ! Fourier from Spectral
    !$OMP BARRIER

    CALL SpecToFour()

    ! DelLam where required

    IF (willDelLam) THEN
       CALL ThreadDecomp(1, mMax, mFirst, mLast, "DoSpecToGrid")
       CALL DelLam(mFirst, mLast)
    END IF

    ! FFT Fourier to Grid
    !$OMP BARRIER

    DO iBlockFFT = 1, nBlockFFT
       p => BlockFFT(iBlockFFT)
       first = (p%firstLat-1) * nVertGrid + 1
       last  = (p%lastLat   ) * nVertGrid 
       CALL ThreadDecomp(first, last, FFTFirst, FFTLast, "DoSpecToGrid")
       FFTSize  = FFTLast - FFTFirst + 1
       IF (FFTSize.le.0) CYCLE
       IF(tamBlock == 0) THEN
          CALL InvFFTTrans (Four(FFTfirst ,1), dvdlj, dip1, p%nLong, FFTSize, &
               p%Trigs, p%nTrigs, p%Factors, p%nFactors)       
       ELSE
          firstBlock = FFTFirst
          DO
             lastBlock = MIN(firstBlock+tamBlock-1, FFTLast)
             sizeBlock = lastBlock-firstBlock+1
             CALL InvFFTTrans (Four(firstBlock,1), dvdlj, dip1, p%nLong, sizeBlock, &
                  p%Trigs, p%nTrigs, p%Factors, p%nFactors)
             firstBlock = firstBlock + tamBlock
             IF (firstBlock > FFTLast) THEN
                EXIT
             END IF
          END DO
       END IF

    END DO

    END IF
    !$OMP BARRIER

    ! Withdraw Grid Fields

    CALL WithdrawGrid

  END SUBROUTINE DoSpecToGrid






  SUBROUTINE DepositSpec(mnFirst, mnLast, mnExtFirst, mnExtLast, mFirst, mLast)
    INTEGER, INTENT(IN) :: mnFirst
    INTEGER, INTENT(IN) :: mnLast
    INTEGER, INTENT(IN) :: mnExtFirst
    INTEGER, INTENT(IN) :: mnExtLast
    INTEGER, INTENT(IN) :: mFirst
    INTEGER, INTENT(IN) :: mLast
    INTEGER :: is
    INTEGER :: imn
    INTEGER :: iv
    INTEGER :: lastv
    REAL(KIND=r8), POINTER :: s1(:)
    REAL(KIND=r8), POINTER :: s2(:,:)

    DO is = 1, nSpecFields

       IF (surfSpec(is)) THEN
          s1 => Spec1d(is)%p
          lastv = prevSpec(is) + 1
          IF (SIZE(s1,1) == 2*myMNExtMax) THEN
             !CDIR NODEP
             DO imn = mnExtFirst, mnExtLast
                Spec(lmnExtMap(imn), lastv)           = s1(2*imn-1)
                Spec(lmnExtMap(imn), lastv+nVertSpec) = s1(2*imn  )
             END DO
          ELSE 
             !CDIR NODEP
             DO imn = mnFirst, mnLast
                Spec(lmnMap(imn), lastv)           = s1(2*imn-1)
                Spec(lmnMap(imn), lastv+nVertSpec) = s1(2*imn  )
             END DO
             !CDIR NODEP
             DO imn = mFirst, mLast
                Spec(lmnZero(imn), lastv)           = 0.0_r8
                Spec(lmnZero(imn), lastv+nVertSpec) = 0.0_r8
             END DO
          END IF

       ELSE
          s2 => Spec2d(is)%p
          lastv = prevSpec(is)
          IF (SIZE(s2,1) == 2*myMNExtMax) THEN
             DO iv = 1, kMaxloc
                !CDIR NODEP
                DO imn = mnExtFirst, mnExtlast
                   Spec(lmnExtMap(imn), iv+lastv)           = s2(2*imn-1, iv)
                   Spec(lmnExtMap(imn), iv+lastv+nVertSpec) = s2(2*imn  , iv)
                END DO
             END DO
          ELSE 
             DO iv = 1, kMaxloc
                !CDIR NODEP
                DO imn = mnFirst, mnLast
                   Spec(lmnMap(imn), iv+lastv)           = s2(2*imn-1, iv)
                   Spec(lmnMap(imn), iv+lastv+nVertSpec) = s2(2*imn  , iv)
                END DO
             END DO
             DO iv = 1, kMaxloc
                !CDIR NODEP
                DO imn = mFirst, mLast
                   Spec(lmnZero(imn), iv+lastv)           = 0.0_r8
                   Spec(lmnZero(imn), iv+lastv+nVertSpec) = 0.0_r8
                END DO
             END DO
          END IF
       END IF
    END DO

    !$OMP SINGLE
    DO imn = myMNExtMax+1, dlmn
       Spec(imn,:) = 0.0_r8
    END DO
    DO iv = 2*nVertSpec+1, dv
       Spec(:,iv) = 0.0_r8
    END DO
    !$OMP END SINGLE
  END SUBROUTINE DepositSpec





  SUBROUTINE SpecToFour()

    INTEGER :: lm
    INTEGER :: m
    INTEGER :: j
    INTEGER :: v
    INTEGER :: js
    INTEGER :: jn
    INTEGER :: kn
    INTEGER :: kdim
    INTEGER :: ldim
    INTEGER :: ks
    INTEGER :: k
    INTEGER :: inin
    INTEGER :: inis
    INTEGER :: inin1
    INTEGER :: inis1
    INTEGER :: mm
    INTEGER :: comm
    INTEGER :: ierr
    INTEGER :: index
    INTEGER :: jfirst
    INTEGER :: jlast
    REAL(KIND=r8) :: FoEv(dv, djh)
    REAL(KIND=r8) :: FoOd(dv, djh)


    CALL ThreadDecomp(1,myjmax_f,jfirst,jlast,"SpecToFour")
    kdim = nVertSpec*myjMax_f*2*MMaxlocal
    ldim = nVertSpec*jMaxlocal_f*2*myMMax
    !$OMP SINGLE
    mGlob = 0
    IF (dimsendbuf.lt.ldim*maxnodes_four) THEN
       dimsendbuf = ldim*maxnodes_four
       DEALLOCATE (bufsend)
       ALLOCATE (bufsend(dimsendbuf))
    ENDIF
    IF (dimrecbuf.lt.kdim*maxnodes_four) THEN
       dimrecbuf = kdim*maxnodes_four
       DEALLOCATE (bufrec)
       ALLOCATE (bufrec(dimrecbuf))
    ENDIF
    bufsend = 0.0_r8
    bufrec  = 0.0_r8
    !$OMP END SINGLE
   
    DO
       !$OMP CRITICAL(lmcrit)
       mGlob = mGlob + 1
       lm = mGlob
       !$OMP END CRITICAL(lmcrit)
       IF (lm > myMMax) EXIT

       m = lm2m(lm)

       ! Spectral to Fourier Even/Odd

       FoEv = 0.0_r8 
       CALL mmt (LS2F(jMinPerM(m),FirstNEven(lm)), djh, &
            Spec(FirstNEven(lm),1), dlmn, &
            FoEv(1,jMinPerM(m)), dv, djh, &
            jMaxHalf-jMinPerM(m)+1, 2*nVertSpec, nEven(lm))
       FoOd = 0.0_r8
       CALL mmt (LS2F(jMinPerM(m),FirstNOdd(lm)), djh, &
            Spec(FirstNOdd(lm),1), dlmn, &
            FoOd(1,jMinPerM(m)), dv, djh, &
            jMaxHalf-jMinPerM(m)+1, 2*nVertSpec, nOdd(lm))

       ! Fourier Even/Odd to Full Fourier

       DO jn = jMinPerM(m), jMaxHalf
          js = jmax-jn+1
          kn = nodeHasJ_f(jn)
          ks = nodeHasJ_f(js)
          inin = (jn-firstlatinproc_f(kn)+nlatsinproc_f(kn)*2*(lm-1))*nVertSpec
          inin1 = inin + nlatsinproc_f(kn)*nVertSpec
          inis = (js-firstlatinproc_f(ks)+nlatsinproc_f(ks)*2*(lm-1))*nVertSpec
          inis1 = inis + nlatsinproc_f(ks)*nVertSpec
          !CDIR NODEP
          DO v = 1, nVertSpec
             bufsend(inin +v+kn*ldim) = FoEv(v,          jn) + FoOd(v,          jn)
             bufsend(inin1+v+kn*ldim) = FoEv(v+nVertSpec,jn) + FoOd(v+nVertSpec,jn)
             bufsend(inis +v+ks*ldim) = FoEv(v,          jn) - FoOd(v,          jn)
             bufsend(inis1+v+ks*ldim) = FoEv(v+nVertSpec,jn) - FoOd(v+nVertSpec,jn)
          END DO
       END DO
    END DO
    !$OMP BARRIER
    !$OMP SINGLE
    comm = COMM_FOUR 
    requestr(myid_four) = MPI_REQUEST_NULL
    requests(myid_four) = MPI_REQUEST_NULL
    DO k=0,MaxNodes_four-1
       IF(k.ne.myid_four) THEN
          CALL MPI_IRECV(bufrec(1+k*kdim),2*myjmax_f*nvertspec*MsPerProc(k), &
               MPI_DOUBLE_PRECISION,k,97,comm,requestr(k),ierr)
       ENDIF
    ENDDO
    DO k=0,MaxNodes_four-1
       IF(k.ne.myid_four) THEN
          CALL MPI_ISEND(bufsend(1+k*ldim),2*nlatsinproc_f(k)*nvertspec*mymmax, &
               MPI_DOUBLE_PRECISION,k,97,comm,requests(k),ierr)
       ENDIF
    ENDDO
    !$OMP END SINGLE
    mm = 0
    DO m=1,Mmax
       IF (myid_four.eq.nodeHasM(m,mygroup_four)) THEN
          mm = mm+1
          inin = (myjmax_f*2*(mm-1))*nVertSpec
          inin1 = inin + myjmax_f*nVertSpec
          DO j=jfirst,jlast
              DO v = 1, nVertSpec
                 Four((j-1)*nVertGrid + v, 2*m-1) = &
                                bufsend((j-1)*nVertSpec+v+inin+myid_four*ldim)
                 Four((j-1)*nVertGrid + v, 2*m  ) = &
                                bufsend((j-1)*nVertSpec+v+inin1+myid_four*ldim)
             END DO
          END DO
       ENDIF
    END DO
    DO k=1,MaxNodes_four-1
       !$OMP BARRIER
       !$OMP SINGLE
       CALL MPI_WAITANY(MaxNodes_four,requestr(0),index,status,ierr)
       !$OMP END SINGLE
       ks = status(MPI_SOURCE)
       mm = 0
       DO m=1,Mmax
          IF (ks.eq.nodeHasM(m,mygroup_four)) THEN
             mm = mm+1
             inin = (myjmax_f*2*(mm-1))*nVertSpec
             inin1 = inin + myjmax_f*nVertSpec
             DO j=jfirst,jlast
                 DO v = 1, nVertSpec
                    Four((j-1)*nVertGrid + v, 2*m-1) = &
                                   bufrec((j-1)*nVertSpec+v+inin+ks*kdim)
                    Four((j-1)*nVertGrid + v, 2*m  ) = &
                                   bufrec((j-1)*nVertSpec+v+inin1+ks*kdim)
                END DO
             END DO
          ENDIF
       END DO
    ENDDO
    !$OMP SINGLE
    CALL MPI_WAITALL(maxnodes_four,requests(0),stat,ierr)
    !$OMP END SINGLE
  END SUBROUTINE SpecToFour





  SUBROUTINE DelLam(mFirst, mLast)
    INTEGER, INTENT(IN) :: mFirst
    INTEGER, INTENT(IN) :: mLast

    INTEGER :: m
    INTEGER :: j
    INTEGER :: v
    REAL(KIND=r8)    :: auxRe
    INTEGER :: ig
    INTEGER :: vBaseFrom
    INTEGER :: vBaseTo
    INTEGER :: jAux
    INTEGER :: vMax

    DO m = mFirst, mLast


       DO ig = 1, nGridFields
          IF (fieldForDelLam(ig)) THEN
             vBaseFrom = prevVertDelLamSource(ig)
             vBaseTo   = prevGrid(ig)

             IF (surfGrid(ig)) THEN
                IF (havesurf) THEN
                   vMax = 1
                 ELSE
                   vMax = 0
                END IF
             ELSE
                vMax = kMaxloc
             END IF
             !CDIR NODEP
             DO v = 1, vMax
                DO j = 1, myjMax_f
                   jAux = (j-1)*nVertGrid + v
                   auxRe = consIm(m) * Four(jAux + vBaseFrom, 2*m-1)
                   Four(jAux + vBaseTo, 2*m-1) = consRe(m) * Four(jAux + vBaseFrom, 2*m  )
                   Four(jAux + vBaseTo, 2*m  ) = auxRe
                END DO
             END DO
          END IF
       END DO
    END DO
  END SUBROUTINE DelLam





  SUBROUTINE WithdrawGrid

    INTEGER :: ig
    INTEGER :: n
    INTEGER :: m  
    INTEGER :: l
    INTEGER :: k
    INTEGER :: kl
    INTEGER :: ns
    INTEGER :: j
    INTEGER :: i
    INTEGER :: ipar
    INTEGER :: ib
    INTEGER :: v0
    INTEGER :: proc
    INTEGER :: comm
    INTEGER :: index
    INTEGER :: ierr
    INTEGER :: ibr(nrecs_f+1)
    INTEGER :: ibn(nrecs_f+1)
    INTEGER :: ibs(nrecs_g+1)
    REAL(KIND=r8), POINTER :: g2(:,:)
    REAL(KIND=r8), POINTER :: g3(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(WithdrawGrid)**"

    comm = MPI_COMM_WORLD
    m = 0
    ibr(1) = 1
    DO n = 1, nrecs_f
       ibn(n) = 0
       DO ipar = m+1,messproc_f(2,n)
          ibn(n) = ibn(n) + 1 + messages_f(2,ipar)-messages_f(1,ipar)
       END DO
       ibr(n+1) = ibr(n) + ibn(n) * nvertgrid
       m = messproc_f(2,n)
    END DO
    m = 0
    ibs(1) = 1
    DO n = 1, nrecs_g
    ib = 0
       DO ipar = m+1,messproc_g(2,n)
          ib = ib + 1 + messages_g(2,ipar)-messages_g(1,ipar)
       END DO
       k = klast_four(messproc_g(1,n))-kfirst_four(messproc_g(1,n))+1
       ns = 0
       IF (kfirst_four(messproc_g(1,n)).eq.1) ns = nsurf_g
       ibs(n+1) = ibs(n) + ib * (k*nfull_g+ns)
       m = messproc_g(2,n)
    END DO
    !$OMP SINGLE
    IF (dimsendbuf.lt.ibs(nrecs_g+1)-1) THEN
       dimsendbuf = ibs(nrecs_g+1)
       DEALLOCATE (bufsend)
       ALLOCATE (bufsend(dimsendbuf))
    ENDIF
    IF (dimrecbuf.lt.ibr(nrecs_f+1)-1) THEN
       dimrecbuf = ibr(nrecs_f+1) 
       DEALLOCATE (bufrec)
       ALLOCATE (bufrec(dimrecbuf))
    ENDIF
    DO n = 1, nrecs_g
       proc = messproc_g(1,n)
       CALL MPI_IRECV(bufsend(ibs(n)),ibs(n+1)-ibs(n),MPI_DOUBLE_PRECISION, &
                     proc,78,comm,requestr(n),ierr)
    END DO
    nglob = 1
    iglob = 0
    ibglob = 0
    !$OMP END SINGLE
    DO
       !$OMP CRITICAL(nigcrit)
       iglob = iglob + 1
       IF (iglob.gt.nGridFields) THEN
          iglob = 1
          nglob = nglob + 1
       ENDIF
       n = nglob
       IF (nglob.le.nrecs_f) THEN
          ig = iglob
          ib = ibglob
          kl = kmaxloc
          IF (surfGrid(ig)) THEN
             IF (havesurf) THEN 
                kl = 1
              ELSE
                kl = 0
              ENDIF
          ENDIF
          ibglob = ibglob + kl * ibn(n)
       ENDIF
       !$OMP END CRITICAL(nigcrit)
       IF (n.gt.nrecs_f) EXIT
       m = messproc_f(2,n-1)
       IF (surfGrid(ig)) THEN
          IF (havesurf) THEN
             DO ipar = m+1,messproc_f(2,n)
                j = messages_f(3,ipar)
                v0 = previousJ(j)*nVertGrid + prevGrid(ig)
                DO i = messages_f(1,ipar),messages_f(2,ipar)
                   ib = ib + 1
                   bufrec(ib) = Four(v0+1, i)
                END DO
             END DO
          END IF
       ELSE
          DO ipar = m+1,messproc_f(2,n)
             j = messages_f(3,ipar)
             v0 = previousJ(j)*nVertGrid + prevGrid(ig)
             DO k=1,kmaxloc
                DO i = messages_f(1,ipar),messages_f(2,ipar)
                   ib = ib + 1
                   bufrec(ib) = Four(k+v0, i) 
                END DO
             END DO
          END DO
       END IF
    END DO
    !$OMP BARRIER
    !$OMP SINGLE
    DO n = 1, nrecs_f
       proc = messproc_f(1,n)
       CALL MPI_ISEND(bufrec(ibr(n)),ibr(n+1)-ibr(n),MPI_DOUBLE_PRECISION,proc,78,&
                     comm,requests(n),ierr)
    END DO
    !$OMP END SINGLE
    iglob = 1
    jglob = max(myfirstlat,myfirstlat_f)
    !$OMP BARRIER
    !
    ! local values
    ! ------------
    IF (max(myfirstlat,myfirstlat_f).le.min(mylastlat,mylastlat_f)) THEN
       DO 
          !$OMP CRITICAL(loccrit)
          ig = iglob
          DO 
             IF (iglob.gt.nGridFields) EXIT
             IF (.not.surfGrid(iglob).or.myfirstlev.eq.1) EXIT
             iglob = iglob + 1
          ENDDO
          ig = iglob
          IF (iglob.le.nGridFields) THEN
             j = jglob
             IF (jglob.ge.min(mylastlat,mylastlat_f)) THEN
                iglob = iglob + 1
                jglob = max(myfirstlat,myfirstlat_f)
              ELSE
                jglob = jglob + 1
             ENDIF
          ENDIF
          !$OMP END CRITICAL(loccrit)
          IF (ig.gt.nGridFields) EXIT
          IF (surfGrid(ig)) THEN
             g2 => Grid2d(ig)%p
             v0 = previousJ(j)*nVertGrid + prevGrid(ig)
             DO i = myfirstlon(j),mylastlon(j)
                g2(ibPerIJ(i,j),jbPerIJ(i,j)) = Four(v0+1,i)
             END DO
          ELSE
             g3 => Grid3d(ig)%p
             v0 = previousJ(j)*nVertGrid + prevGrid(ig)
             DO k = myfirstlev,mylastlev
                v0 = v0 + 1
                DO i = myfirstlon(j),mylastlon(j)
                   g3(ibPerIJ(i,j),k,jbPerIJ(i,j)) = Four(v0, i)
                END DO
             END DO
          END IF
       END DO
    END IF

    !
    !$OMP SINGLE
    kountg = 1
    !$OMP END SINGLE
    IF (nrecs_g.gt.0) THEN
       DO 
          !$OMP SINGLE
          CALL MPI_WAITANY(nrecs_g,requestr(1),index,status,ierr)
          ksg = status(MPI_SOURCE)
          DO l=1,nrecs_g
             IF (ksg.eq.messproc_g(1,l)) THEN
                nglob = l
                ibglob = ibs(nglob) - 1
                mglob = messproc_g(2,nglob-1)
                EXIT
             ENDIF
          ENDDO
          iglob = 1
          ipar2g = mglob + 1
          ipar3g = mglob + 1
          !$OMP END SINGLE
          n = nglob
          m = mglob
          kl =  klast_four(ksg) - kfirst_four(ksg) + 1
          DO 
             !$OMP CRITICAL(reccrit)
             ig = iglob
             DO 
                IF (iglob.gt.nGridFields) EXIT
                IF (.not.surfGrid(iglob).or.kfirst_four(ksg).eq.1) EXIT 
                iglob = iglob + 1
             ENDDO
             ig = iglob
             IF (iglob.le.nGridFields) THEN
                ib = ibglob
                IF (surfGrid(ig)) THEN
                   ipar = ipar2g
                   ibglob = ibglob + messages_g(2,ipar) - messages_g(1,ipar) + 1
                   ipar2g = ipar2g + 1
                   IF (ipar2g.gt.messproc_g(2,n)) THEN
                      ipar2g = m + 1
                      iglob = iglob + 1
                   ENDIF
                 ELSE
                   ipar = ipar3g
                   ibglob = ibglob + kl * (messages_g(2,ipar) - messages_g(1,ipar) + 1)
                   ipar3g = ipar3g + 1
                   IF (ipar3g.gt.messproc_g(2,n)) THEN
                      ipar3g = m + 1
                      iglob = iglob + 1
                   ENDIF
                ENDIF
             ENDIF
             !$OMP END CRITICAL(reccrit)
             IF (ig.gt.nGridFields) EXIT
             IF (surfGrid(ig)) THEN
                g2 => Grid2d(ig)%p
                j = messages_g(3,ipar)
                v0 = previousJ(j)*nVertGrid + prevGrid(ig)
                DO i = messages_g(1,ipar),messages_g(2,ipar)
                   ib = ib + 1
                   g2(ibperij(i,j),jbperij(i,j)) = bufsend(ib)
                END DO
             ELSE
                g3 => Grid3d(ig)%p
                j = messages_g(3,ipar)
                v0 = previousJ(j)*nVertGrid + prevGrid(ig)
                DO k=kfirst_four(ksg),klast_four(ksg)
                   DO i = messages_g(1,ipar),messages_g(2,ipar)
                      ib = ib + 1
                      g3(ibperij(i,j),k,jbperij(i,j)) = bufsend(ib)
                   END DO
                END DO
             END IF
          END DO
          !$OMP SINGLE
          kountg = kountg + 1
          !$OMP END SINGLE
          IF (kountg.GT.nrecs_g) EXIT
       END DO
    END IF

    !$OMP SINGLE
    CALL MPI_WAITALL(nrecs_f,requests(1),stat,ierr)
    !$OMP END SINGLE

  END SUBROUTINE WithdrawGrid





  SUBROUTINE CreateGridToSpec(nFull, nSurf)
    INTEGER, INTENT(IN) :: nFull
    INTEGER, INTENT(IN) :: nSurf
    CHARACTER(LEN=*), PARAMETER :: h="**(CreateGridToSpec)**"
    INTEGER :: nsu

    IF (havesurf) THEN
       nsu = nSurf
      ELSE
       nsu = 0
    ENDIF
    nFull_g = nFull  
    nSurf_g = nSurf  

    nSpecFields = nFull + nsu
    usedSpecFields = 0
    lastUsedSpecVert = 0
    IF (.NOT.ALLOCATED(surfSpec))ALLOCATE(surfSpec(nSpecFields))
    IF (.NOT.ALLOCATED(prevSpec))ALLOCATE(prevSpec(nSpecFields))
    IF (.NOT.ALLOCATED(Spec1d))ALLOCATE(Spec1d(nSpecFields))
    IF (.NOT.ALLOCATED(Spec2d))ALLOCATE(Spec2d(nSpecFields))
    nVertSpec = nFull*kMaxloc + nsu

    nGridFields = nFull + nSurf
    usedGridFields = 0
    lastUsedGridVert = 0
    IF (.NOT.ALLOCATED(surfGrid)) ALLOCATE(surfGrid(nGridFields))
    IF (.NOT.ALLOCATED(prevGrid)) ALLOCATE(prevGrid(nGridFields))
    IF (.NOT.ALLOCATED(Grid2d)) ALLOCATE(Grid2d(nGridFields))
    IF (.NOT.ALLOCATED(Grid3d)) ALLOCATE(Grid3d(nGridFields))
    IF (.NOT.ALLOCATED(fieldForDelLam)) ALLOCATE(fieldForDelLam(nGridFields))
    nVertGrid = nFull*kMaxloc + nsu

    willDelLam = .FALSE.
    usedDelLamFields = nSpecFields
    lastUsedDelLamVert = nVertSpec
    IF (.NOT.ALLOCATED(prevVertDelLamSource)) ALLOCATE(prevVertDelLamSource(nGridFields))

    dv = NoBankConflict(2*nVertSpec)
    dvjh = NoBankConflict(nVertSpec*jMaxHalf)
    dvdlj = NoBankConflict(nVertGrid*myJMax_f)

    IF (.NOT.ALLOCATED(Spec)) ALLOCATE (Spec(dlmn, dv))
    IF (.NOT.ALLOCATED(Four)) ALLOCATE (Four(dvdlj, dip1))

    ALLOCATE(mnodes(0:maxnodes_four-1))
    IF (.NOT.ALLOCATED(requests)) THEN
       ALLOCATE(requests(0:maxnodes))
       ALLOCATE(requestr(0:maxnodes))
    END IF
    IF (.NOT.ALLOCATED(status)) THEN
       ALLOCATE(status(MPI_STATUS_SIZE))
    END IF
    IF (.NOT.ALLOCATED(stat)) THEN
       ALLOCATE(stat(MPI_STATUS_SIZE,maxnodes))
    END IF


    IF (dumpLocal) THEN
       WRITE(nfprt,"(a,' Spec: n, used, lastVert, nVert=',4i5)") &
            h, nSpecFields, usedSpecFields, lastUsedSpecVert, nVertSpec
       WRITE(nfprt,"(a,' Grid: n, used, lastVert, nVert=',4i5)") &
            h, nGridFields, usedGridFields, lastUsedGridVert, nVertGrid
    END IF

  END SUBROUTINE CreateGridToSpec




  SUBROUTINE DoGridToSpec()
    INTEGER :: mnFirst
    INTEGER :: mnLast
    INTEGER :: mnExtFirst
    INTEGER :: mnExtLast
    INTEGER :: jbFirst
    INTEGER :: jbLast
    INTEGER :: iBlockFFT
    INTEGER :: first 
    INTEGER :: last
    INTEGER :: FFTfirst 
    INTEGER :: FFTLast
    INTEGER :: FFTSize
    INTEGER :: firstBlock
    INTEGER :: lastBlock
    INTEGER :: sizeBlock
    TYPE(MultiFFT), POINTER :: p
    CHARACTER(LEN=*), PARAMETER :: h="**(DoGridToSpec)**"

    ! all fields were deposited?

    IF (usedSpecFields /= nSpecFields) THEN
       WRITE(nfprt,"(a, ' not all spectral fields were deposited')") h
       STOP h
    ELSE IF (usedGridFields /= nGridFields) THEN
       WRITE(nfprt,"(a, ' not all gauss fields were deposited')") h
       STOP h
    END IF

    ! start OMP parallelism

    CALL ThreadDecomp(1, mymnMax, mnFirst, mnLast, "DoGridToSpec")
    CALL ThreadDecomp(1, mymnExtMax, mnExtFirst, mnExtLast, "DoGridToSpec")
    CALL ThreadDecomp(1, jbMax, jbFirst, jbLast, "DoGridToSpec")


    ! deposit all grid fields

    CALL DepositGrid
    !$OMP BARRIER

    ! FFT Grid to Fourier

    DO iBlockFFT = 1, nBlockFFT
       p => BlockFFT(iBlockFFT)
       first = (p%firstLat-1) * nVertGrid + 1
       last  = (p%lastLat   ) * nVertGrid 
       CALL ThreadDecomp(first, last, FFTFirst, FFTLast, "DoGridToSpec")
       FFTSize  = FFTLast - FFTFirst + 1
       IF (FFTSize.le.0) CYCLE
       IF(tamBlock == 0) THEN
          CALL DirFFTTrans (Four(FFTfirst ,1), dvdlj, dip1, p%nLong, FFTSize, &
               p%Trigs, p%nTrigs, p%Factors, p%nFactors)       
       ELSE
          firstBlock = FFTFirst
          DO
             lastBlock = MIN(firstBlock+tamBlock-1, FFTLast)
             sizeBlock = lastBlock-firstBlock+1
             CALL DirFFTTrans (Four(firstBlock ,1), dvdlj, dip1, p%nLong, sizeBlock, &
                               p%Trigs, p%nTrigs, p%Factors, p%nFactors)
             firstBlock = firstBlock + tamBlock
             IF (firstBlock > FFTLast) THEN
                EXIT
             END IF
          END DO
       END IF
       
       !CALL DirFFTTrans (Four(FFTfirst ,1), dvdlj, dip1, p%nLong, FFTSize, &
       !    p%Trigs, p%nTrigs, p%Factors, p%nFactors)

    END DO

    !$OMP BARRIER

    ! Fourier to Spectral

    CALL FourToSpec()
    !$OMP BARRIER

    ! retrieve Spectral fields

    CALL WithdrawSpectral(mnFirst, mnLast, mnExtFirst, mnExtLast)
  END SUBROUTINE DoGridToSpec





  SUBROUTINE DepositGrid

    INTEGER :: ig
    INTEGER :: n
    INTEGER :: m  
    INTEGER :: l
    INTEGER :: k
    INTEGER :: kl
    INTEGER :: ns
    INTEGER :: j
    INTEGER :: i
    INTEGER :: ipar
    INTEGER :: ib
    INTEGER :: v0
    INTEGER :: proc
    INTEGER :: comm
    INTEGER :: index
    INTEGER :: ierr
    INTEGER :: ibr(nrecs_f+1)
    INTEGER :: ibs(nrecs_g+1)
    REAL(KIND=r8), POINTER :: g2(:,:)
    REAL(KIND=r8), POINTER :: g3(:,:,:)

    comm = MPI_COMM_WORLD
    m = 0
    ibr(1) = 1
    DO n = 1, nrecs_f
       ib = 0
       DO ipar = m+1,messproc_f(2,n)
          ib = ib + 1 + messages_f(2,ipar)-messages_f(1,ipar)
       END DO
       ibr(n+1) = ibr(n) + ib * nvertgrid
       m = messproc_f(2,n)
    END DO
    m = 0
    ibs(1) = 1
    DO n = 1, nrecs_g
       ib = 0
       DO ipar = m+1,messproc_g(2,n)
          ib = ib + 1 + messages_g(2,ipar)-messages_g(1,ipar)
       END DO
       k = klast_four(messproc_g(1,n))-kfirst_four(messproc_g(1,n))+1
       ns = 0
       IF (kfirst_four(messproc_g(1,n)).eq.1) ns = nsurf_g
       ibs(n+1) = ibs(n) + ib * (k*nfull_g+ns)
       m = messproc_g(2,n)
    END DO
    !$OMP SINGLE
    IF (dimsendbuf.lt.ibs(nrecs_g+1)-1) THEN
       dimsendbuf = ibs(nrecs_g+1)
       DEALLOCATE (bufsend)
       ALLOCATE (bufsend(dimsendbuf))
    ENDIF
    IF (dimrecbuf.lt.ibr(nrecs_f+1)-1) THEN
       dimrecbuf = ibr(nrecs_f+1) 
       DEALLOCATE (bufrec)
       ALLOCATE (bufrec(dimrecbuf))
    ENDIF
    DO n = 1, nrecs_f
       proc = messproc_f(1,n)
       CALL MPI_IRECV(bufrec(ibr(n)),ibr(n+1)-ibr(n),MPI_DOUBLE_PRECISION, &
                     proc,77,comm,requestr(n),ierr)
    END DO
    mglob = 0
    ibglob = 0
    ipar2g = mglob + 1
    ipar3g = mglob + 1
    kountg = 0
    !$OMP END SINGLE
    DO
       !$OMP SINGLE
       iglob = 1
       kountg = kountg + 1
       !$OMP END SINGLE
       IF (kountg.GT.nrecs_g) EXIT
       n = kountg
       proc = messproc_g(1,n)
       kl = klast_four(proc) - kfirst_four(proc) + 1
       DO 
          !$OMP CRITICAL(sndcrit)
          ig = iglob
          DO 
            IF (iglob.gt.nGridFields) EXIT
            IF (.not.surfGrid(iglob).or.kfirst_four(proc).eq.1) EXIT
            iglob = iglob + 1
          ENDDO
          ig = iglob 
          IF (iglob.le.nGridFields) THEN
             ib = ibglob
             IF (surfGrid(ig)) THEN
                ipar = ipar2g
                ipar2g = ipar2g + 1
                IF (ipar2g.gt.messproc_g(2,n)) THEN
                   ipar2g = mglob + 1
                   iglob = iglob + 1
                ENDIF
                ibglob = ibglob + messages_g(2,ipar) - messages_g(1,ipar) + 1
              ELSE
                ipar = ipar3g
                ipar3g = ipar3g + 1
                IF (ipar3g.gt.messproc_g(2,n)) THEN
                   ipar3g = mglob + 1
                   iglob = iglob + 1
                ENDIF
                ibglob = ibglob + (messages_g(2,ipar) - messages_g(1,ipar) + 1) * kl
             ENDIF
          ENDIF
          !$OMP END CRITICAL(sndcrit)
          IF (ig.gt.nGridFields) EXIT
          IF (surfGrid(ig)) THEN
             g2 => Grid2d(ig)%p
             j = messages_g(3,ipar)
             DO i = messages_g(1,ipar),messages_g(2,ipar)
                 ib = ib + 1
                 bufsend(ib) = g2(ibperij(i,j),jbperij(i,j))
             END DO
          ELSE
             g3 => Grid3d(ig)%p
             j = messages_g(3,ipar)
             DO k=kfirst_four(proc),klast_four(proc)
                DO i = messages_g(1,ipar),messages_g(2,ipar)
                   ib = ib + 1
                   bufsend(ib) = g3(ibperij(i,j),k,jbperij(i,j))
                END DO
             END DO
          END IF
       END DO
       !$OMP BARRIER
       !$OMP SINGLE
       CALL MPI_ISEND(bufsend(ibs(n)),ibs(n+1)-ibs(n),MPI_DOUBLE_PRECISION,proc,77,&
                     comm,requests(n),ierr)
       mglob = messproc_g(2,n)
       ipar2g = mglob + 1
       ipar3g = mglob + 1
       !$OMP END SINGLE
    END DO
    !$OMP SINGLE
    iglob = 1
    jglob = max(myfirstlat,myfirstlat_f)
    !$OMP END SINGLE
    !
    ! local values
    ! ------------
    IF (max(myfirstlat,myfirstlat_f).le.min(mylastlat,mylastlat_f)) THEN
       DO
          !$OMP CRITICAL(loccrit)
          ig = iglob
          DO
             IF (iglob.gt.nGridFields) EXIT
             IF (.not.surfGrid(iglob).or.myfirstlev.eq.1) EXIT
             iglob = iglob + 1
          ENDDO
          ig = iglob
          IF (iglob.le.nGridFields) THEN
             j = jglob
             IF (jglob.ge.min(mylastlat,mylastlat_f)) THEN
                iglob = iglob + 1
                jglob = max(myfirstlat,myfirstlat_f)
              ELSE
                jglob = jglob + 1
             ENDIF
          ENDIF
          !$OMP END CRITICAL(loccrit)
          IF (ig.gt.nGridFields) EXIT
    
          IF (surfGrid(ig)) THEN
             g2 => Grid2d(ig)%p
             v0 = previousJ(j)*nVertGrid + prevGrid(ig)
             DO i = myfirstlon(j),mylastlon(j)
                Four(v0+1, i) = g2(ibPerIJ(i,j),jbPerIJ(i,j))
             END DO
          ELSE
             g3 => Grid3d(ig)%p
             v0 = previousJ(j)*nVertGrid + prevGrid(ig)
             DO k = myfirstlev,mylastlev
                v0 = v0 + 1
                DO i = myfirstlon(j),mylastlon(j)
                   Four(v0, i) = g3(ibPerIJ(i,j),k,jbPerIJ(i,j))
                END DO
             END DO
          END IF
       END DO
    END IF

    ! Nullify remaining Fourier

    !$OMP SINGLE
    DO i = 2*iMax+1, dip1
       Four(:,i) = 0.0_r8
    END DO

    DO j = nVertGrid*myJMax_f+1, dvdlj
       Four(j,:) = 0.0_r8
    END DO
    kountg = 0
    !$OMP END SINGLE

    DO 
       !$OMP SINGLE
       kountg = kountg + 1
       !$OMP END SINGLE
       IF (kountg.GT.nrecs_f) EXIT
       !$OMP SINGLE
       CALL MPI_WAITANY(nrecs_f,requestr(1),index,status,ierr)
       ksg = status(MPI_SOURCE)
       DO l=1,nrecs_f
          IF (ksg.eq.messproc_f(1,l)) THEN
             nglob = l
             ibglob = ibr(nglob) - 1
             mglob = messproc_f(2,nglob-1)
             EXIT
          ENDIF
       ENDDO
       iglob = 1
       ipar2g = mglob + 1
       ipar3g = mglob + 1
       !$OMP END SINGLE
       n = nglob
       m = mglob
!      kl =  klast_four(ksg) - kfirst_four(ksg) + 1
       kl =  kmaxloc
       DO
          !$OMP CRITICAL(reccrit)
          ig = iglob
          DO 
             IF (iglob.gt.nGridFields) EXIT
             IF (.not.surfGrid(iglob).or.havesurf) EXIT
             iglob = iglob + 1
          ENDDO
          ig = iglob
          IF (iglob.le.nGridFields) THEN
             ib = ibglob
             IF (surfGrid(ig)) THEN
                ipar = ipar2g
                ibglob = ibglob + messages_f(2,ipar) - messages_f(1,ipar) + 1
                ipar2g = ipar2g + 1
                IF (ipar2g.gt.messproc_f(2,n)) THEN
                   ipar2g = m + 1
                   iglob = iglob + 1
                ENDIF
              ELSE
                ipar = ipar3g
                ibglob = ibglob + kl * (messages_f(2,ipar) - messages_f(1,ipar) + 1)
                ipar3g = ipar3g + 1
                IF (ipar3g.gt.messproc_f(2,n)) THEN
                   ipar3g = m + 1
                   iglob = iglob + 1
                ENDIF
             ENDIF
          ENDIF
          !$OMP END CRITICAL(reccrit)
          IF (ig.gt.nGridFields) EXIT

          IF (surfGrid(ig)) THEN
             j = messages_f(3,ipar)
             v0 = previousJ(j)*nVertGrid + prevGrid(ig)
             DO i = messages_f(1,ipar),messages_f(2,ipar)
                ib = ib + 1
                Four(v0+1, i) = bufrec(ib)
             END DO
          ELSE
             j = messages_f(3,ipar)
             v0 = previousJ(j)*nVertGrid + prevGrid(ig)
             DO k=1,kmaxloc
                DO i = messages_f(1,ipar),messages_f(2,ipar)
                   ib = ib + 1
                   Four(k+v0, i) = bufrec(ib)
                END DO
             END DO
          END IF
       END DO
    END DO

    !$OMP SINGLE
    CALL MPI_WAITALL(nrecs_g,requests(1),stat,ierr)
    !$OMP END SINGLE

  END SUBROUTINE DepositGrid




  SUBROUTINE FourToSpec()
    INTEGER :: lm
    INTEGER :: m
    INTEGER :: mn
    INTEGER :: j
    INTEGER :: jl
    INTEGER :: js
    INTEGER :: v
    INTEGER :: k
    INTEGER :: kdim
    INTEGER :: ldim
    INTEGER :: kn
    INTEGER :: ks
    INTEGER :: inin
    INTEGER :: inin1
    INTEGER :: inis
    INTEGER :: inis1
    INTEGER :: comm
    INTEGER :: ierr
    INTEGER :: index
    REAL(KIND=r8) :: FoEv(djh, dv)
    REAL(KIND=r8) :: FoOd(djh, dv)

    kdim = nVertSpec*myjMax_f*2*MMaxlocal
    ldim = nVertSpec*jMaxlocal_f*2*myMMax
    !$OMP SINGLE
    mnodes = 0
    mGlob = 0
    IF (dimsendbuf.lt.ldim*maxnodes_four) THEN
       dimsendbuf = ldim*maxnodes_four
       DEALLOCATE (bufsend)
       ALLOCATE (bufsend(dimsendbuf))
    ENDIF
    IF (dimrecbuf.lt.kdim*maxnodes_four) THEN
       dimrecbuf = kdim*maxnodes_four
       DEALLOCATE (bufrec)
       ALLOCATE (bufrec(dimrecbuf))
    ENDIF
    !$OMP END SINGLE
    DO 
       !$OMP CRITICAL(lmcrit1)
       mGlob = mGlob + 1
       m = mGlob
       IF (m.le.MMax) THEN
          kn = NodeHasM(m,mygroup_four)
          mnodes(kn) = mnodes(kn) + 1
          lm = mnodes(kn)
       ENDIF
       !$OMP END CRITICAL(lmcrit1)
       IF (m > MMax) EXIT
       IF (kn.ne.myid_four) THEN
          !CDIR NODEP
          DO j = MAX(jMinPerM(m),myfirstlat_f),MIN(mylastlat_f,jMaxPerM(m))
             DO v = 1, nVertSpec
                jl = j-myfirstlat_f+1
                bufrec((myjmax_f*2*(lm-1)+(jl-1))*nVertGrid + v+kn*kdim) = &
                                               Four((jl-1)*nVertGrid + v, 2*m-1)
                bufrec((myjmax_f*(2*lm-1)+(jl-1))*nVertGrid + v+kn*kdim) = &
                                               Four((jl-1)*nVertGrid + v, 2*m)
             END DO
          END DO
         ELSE
          DO j = MAX(jMinPerM(m),myfirstlat_f),MIN(mylastlat_f,jMaxPerM(m))
             DO v = 1, nVertSpec
                jl = j-myfirstlat_f+1
                bufsend((myjmax_f*2*(lm-1)+(jl-1))*nVertGrid + v+kn*ldim) = &
                                               Four((jl-1)*nVertGrid + v, 2*m-1)
                bufsend((myjmax_f*(2*lm-1)+(jl-1))*nVertGrid + v+kn*ldim) = &
                                               Four((jl-1)*nVertGrid + v, 2*m)
             END DO
          END DO
       END IF
    END DO
    !$OMP BARRIER
    !$OMP SINGLE
    comm = COMM_FOUR 
    requestr(myid_four) = MPI_REQUEST_NULL
    requests(myid_four) = MPI_REQUEST_NULL
    DO k=0,MaxNodes_four-1
       IF(k.ne.myid_four) THEN
          CALL MPI_IRECV(bufsend(1+k*ldim),2*nlatsinproc_f(k)*nvertspec*mymmax, &
               MPI_DOUBLE_PRECISION,k,98,comm,requestr(k),ierr)
       ENDIF
    ENDDO
    DO k=0,MaxNodes_four-1
       IF(k.ne.myid_four) THEN
          CALL MPI_ISEND(bufrec(1+k*kdim),2*myjmax_f*nvertspec*MsPerProc(k), &
               MPI_DOUBLE_PRECISION,k,98,comm,requests(k),ierr)
       ENDIF
    ENDDO
    DO k=1,MaxNodes_four-1
       CALL MPI_WAITANY(MaxNodes_four,requestr(0),index,status,ierr)
    ENDDO
    DO k=0,MaxNodes_four-1
       IF(k.ne.myid_four) THEN
          CALL MPI_WAIT(requests(k),status,ierr)
       ENDIF
    ENDDO
    !$OMP END SINGLE
    FoEv = 0.0_r8
    FoOd = 0.0_r8

    !$OMP SINGLE
    mGlob = 0
    !$OMP END SINGLE
    DO
       !$OMP CRITICAL(lmcrit2)
       mGlob = mGlob + 1
       lm = mGlob
       !$OMP END CRITICAL(lmcrit2)
       IF (lm > myMMax) EXIT
       m = lm2m(lm)

       !CDIR NODEP
       DO j = jMinPerM(m), jMaxHalf
          js = jMax + 1 - j
          kn = nodeHasJ_f(j)
          ks = nodeHasJ_f(js)
          inin = (nlatsinproc_f(kn)*2*(lm-1)+j-firstlatinproc_f(kn))*nVertSpec
          inin1 = inin + nlatsinproc_f(kn)*nVertSpec
          inis = (nlatsinproc_f(ks)*2*(lm-1)+js-firstlatinproc_f(ks))*nVertSpec
          inis1 = inis + nlatsinproc_f(ks)*nVertSpec
          DO v = 1, nVertSpec
             FoEv(j,v          ) = bufsend(inin+v+kn*ldim) +bufsend(inis+v+ks*ldim)
             FoEv(j,v+nVertSpec) = bufsend(inin1+v+kn*ldim)+bufsend(inis1+v+ks*ldim)
             FoOd(j,v          ) = bufsend(inin+v+kn*ldim) -bufsend(inis+v+ks*ldim)
             FoOd(j,v+nVertSpec) = bufsend(inin1+v+kn*ldim)-bufsend(inis1+v+ks*ldim)
          END DO
       END DO
       DO j = 1, dv
          DO mn = FirstNEven(lm), FirstNEven(lm) + nEven(lm) - 1
             Spec(mn,j) = 0.0_r8
          END DO
       END DO
       CALL mmd (&
            LF2S(FirstNEven(lm),jMinPerM(m)), dlmn, &
            FoEv(jMinPerM(m),1), djh, &
            Spec(FirstNEven(lm),1), dlmn, dv, &
            nEven(lm), 2*nVertSpec, jMaxHalf-jMinPerM(m)+1)
       DO j = 1, dv
          DO mn = FirstNOdd(lm), FirstNOdd(lm) + nOdd(lm) - 1
             Spec(mn,j) = 0.0_r8
          END DO
       END DO
       CALL mmd (&
            LF2S(FirstNOdd(lm),jMinPerM(m)), dlmn, &
            FoOd(jMinPerM(m),1), djh, &
            Spec(FirstNOdd(lm),1), dlmn, dv, &
            nOdd(lm), 2*nVertSpec, jMaxHalf-jMinPerM(m)+1)
    END DO
  END SUBROUTINE FourToSpec





  SUBROUTINE WithdrawSpectral(mnFirst, mnLast, mnExtFirst, mnExtLast)
    INTEGER, INTENT(IN) :: mnFirst
    INTEGER, INTENT(IN) :: mnLast
    INTEGER, INTENT(IN) :: mnExtFirst
    INTEGER, INTENT(IN) :: mnExtLast

    INTEGER :: is
    INTEGER :: imn
    INTEGER :: iv
    INTEGER :: lastv
    REAL(KIND=r8), POINTER :: s1(:)
    REAL(KIND=r8), POINTER :: s2(:,:)


    DO is = 1, nSpecFields

       IF (surfSpec(is)) THEN
          s1 => Spec1d(is)%p
          lastv = prevSpec(is) + 1
          IF (SIZE(s1,1) == 2*myMNExtMax) THEN
             !CDIR NODEP
             DO imn = mnExtFirst, mnExtLast
                s1(2*imn-1) = Spec(lmnExtMap(imn), lastv)
                s1(2*imn  ) = Spec(lmnExtMap(imn), lastv+nVertSpec)
             END DO
          ELSE 
             !CDIR NODEP
             DO imn = mnFirst, mnLast
                s1(2*imn-1) = Spec(lmnMap(imn), lastv)
                s1(2*imn  ) = Spec(lmnMap(imn), lastv+nVertSpec)
             END DO
          END IF
       ELSE
          s2 => Spec2d(is)%p
          lastv = prevSpec(is)
          IF (SIZE(s2,1) == 2*myMNExtMax) THEN
             DO iv = 1, kMaxloc
                !CDIR NODEP
                DO imn = mnExtFirst, mnExtlast
                   s2(2*imn-1, iv) = Spec(lmnExtMap(imn), iv+lastv)
                   s2(2*imn  , iv) = Spec(lmnExtMap(imn), iv+lastv+nVertSpec) 
                END DO
             END DO
          ELSE 
             DO iv = 1, kMaxloc
                !CDIR NODEP
                DO imn = mnFirst, mnLast
                   s2(2*imn-1, iv) = Spec(lmnMap(imn), iv+lastv)       
                   s2(2*imn  , iv) = Spec(lmnMap(imn), iv+lastv+nVertSpec) 
                END DO
             END DO
          END IF
       END IF
    END DO

  END SUBROUTINE WithdrawSpectral


  SUBROUTINE Destroy()
    DEALLOCATE(Spec)
    DEALLOCATE(Four)
    IF (ALLOCATED(mnodes)) THEN
       DEALLOCATE(mnodes)
    END IF
    IF (ALLOCATED(surfSpec)) THEN
       DEALLOCATE(surfSpec)
    END IF
    IF (ALLOCATED(prevSpec)) THEN
       DEALLOCATE(prevSpec)
    END IF
    IF (ALLOCATED(Spec1d)) THEN
       DEALLOCATE(Spec1d)
    END IF
    IF (ALLOCATED(Spec2d)) THEN
       DEALLOCATE(Spec2d)
    END IF

    IF (ALLOCATED(surfGrid)) THEN
       DEALLOCATE(surfGrid)
    END IF
    IF (ALLOCATED(prevGrid)) THEN
       DEALLOCATE(prevGrid)
    END IF
    IF (ALLOCATED(Grid2d)) THEN
       DEALLOCATE(Grid2d)
    END IF
    IF (ALLOCATED(Grid3d)) THEN
       DEALLOCATE(Grid3d)
    END IF
    IF (ALLOCATED(fieldForDelLam)) THEN
       DEALLOCATE(fieldForDelLam)
    END IF
    IF (ALLOCATED(prevVertDelLamSource)) THEN
       DEALLOCATE(prevVertDelLamSource)
    END IF
  END SUBROUTINE Destroy






  SUBROUTINE mmd (A, lda, B, ldb, C, ldc, tdc, ni, nj, nk)
    INTEGER, INTENT(IN ) :: ni , nj,  nk
    INTEGER, INTENT(IN ) :: lda, ldb, ldc, tdc
    REAL(KIND=r8),    INTENT(IN ) :: A(lda, nk)
    REAL(KIND=r8),    INTENT(IN ) :: B(ldb, nj)
    REAL(KIND=r8),    INTENT(INOUT) :: C(ldc, tdc)
    INTEGER i, j, k
    DO i = 1, ni
       DO j = 1, nj
          DO k = 1, nk
             C(i,j)=C(i,j)+A(i,k)*B(k,j)
          END DO
       END DO
    END DO
  END SUBROUTINE mmd




  SUBROUTINE mmt (A, lda, B, ldb, C, ldc, tdc, ni, nj, nk)
    INTEGER, INTENT(IN ) :: ni , nj,  nk
    INTEGER, INTENT(IN ) :: lda, ldb, ldc, tdc
    REAL(KIND=r8),    INTENT(IN ) :: A(lda, nk)
    REAL(KIND=r8),    INTENT(IN ) :: B(ldb, nj)
    REAL(KIND=r8),    INTENT(INOUT) :: C(ldc, tdc)
    INTEGER i, j, k
    CHARACTER(LEN=*), PARAMETER :: h="**(mmt)**"
    DO i = 1, ni
       DO j = 1, nj
          DO k = 1, nk
             C(j,i)=C(j,i)+A(i,k)*B(k,j)
          END DO
       END DO
    END DO
  END SUBROUTINE mmt




  !InvFFT: Computes inverse FFT of 'lot' sequences of 'n+1' input real data
  !        as rows of 'fin', dimensioned 'fin(ldin,tdin)'. Input data is
  !        kept unchanged. Input values 'fin(n+2:ldin,:)' and 
  !        'fin(:,lot+1:tdin)' are not visited. 
  !        Outputs 'lot' sequences of 'n' real data
  !        as rows of 'fout', dimensioned 'fout(ldout,tdout)'. Output
  !        values 'fout(n+1:ldout,:)' and 'fout(:,lot+1:tdout)' are set to 0.



  SUBROUTINE InvFFTTrans (fInOut, ldInOut, tdInOut, n, lot, Trigs, nTrigs, Factors, nFactors)
    INTEGER, INTENT(IN ) :: ldInOut, tdInOut
    REAL(KIND=r8), INTENT(INOUT ) :: fInOut (ldInOut ,tdInOut)
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(IN ) :: lot
    INTEGER, INTENT(IN ) :: nTrigs
    REAL(KIND=r8),    INTENT(IN ) :: Trigs(nTrigs)
    INTEGER, INTENT(IN ) :: nFactors
    INTEGER, INTENT(IN ) :: Factors(nFactors)

    INTEGER :: nh
    INTEGER :: nfax
    INTEGER :: la
    INTEGER :: k
    LOGICAL :: ab2cd
    CHARACTER(LEN=*), PARAMETER :: h="**(InvFFTTrans)**"
    REAL(KIND=r8) :: a(lot,n/2)
    REAL(KIND=r8) :: b(lot,n/2)
    REAL(KIND=r8) :: c(lot,n/2)
    REAL(KIND=r8) :: d(lot,n/2)


    nfax=Factors(1)
    nh=n/2

    CALL SplitFourTrans (fInOut, a, b, ldInOut, tdInOut, n, nh, lot, Trigs, nTrigs)

    la=1
    ab2cd=.TRUE.
    DO k=1,nfax
       IF (ab2cd) THEN
          CALL OnePass (a, b, c, d, lot, nh, Factors(k+1), la, Trigs, nTrigs)
          ab2cd=.FALSE.
       ELSE
          CALL OnePass (c, d, a, b, lot, nh, Factors(k+1), la, Trigs, nTrigs)
          ab2cd=.TRUE.
       END IF
       la=la*Factors(k+1)
    END DO

    IF (ab2cd) THEN
       CALL JoinGridTrans (a, b, fInOut, ldInOut, tdInOut, nh, lot)
    ELSE
       CALL JoinGridTrans (c, d, fInOut, ldInOut, tdInOut, nh, lot)
    END IF
  END SUBROUTINE InvFFTTrans


  ! Split Fourier Fields


  SUBROUTINE SplitFourTrans (fin, a, b, ldin, tdin, n, nh, lot, Trigs, nTrigs)
    INTEGER, INTENT(IN ) :: ldin
    INTEGER, INTENT(IN ) :: tdin
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(IN ) :: nh
    INTEGER, INTENT(IN ) :: lot
    REAL(KIND=r8),    INTENT(IN ) :: fin(ldin, tdin)
    REAL(KIND=r8),    INTENT(OUT) :: a  (lot , nh)
    REAL(KIND=r8),    INTENT(OUT) :: b  (lot , nh)
    INTEGER, INTENT(IN ) :: nTrigs
    REAL(KIND=r8),    INTENT(IN ) :: Trigs(nTrigs)

    INTEGER :: i, j
    REAL(KIND=r8)    :: c, s

    !CDIR NODEP
    DO i = 1, lot
       a(i,1)=fin(i,1)+fin(i,n+1)
       b(i,1)=fin(i,1)-fin(i,n+1)
    END DO

    DO j = 2, (nh+1)/2
       c=Trigs(n+2*j-1)
       s=Trigs(n+2*j  )
       !CDIR NODEP
       DO i = 1, lot
          a(i,j     )=   (fin(i,2*j-1)+fin(i,n+3-2*j)) &
               -      (s*(fin(i,2*j-1)-fin(i,n+3-2*j)) &
               +       c*(fin(i,2*j  )+fin(i,n+4-2*j)))
          a(i,nh+2-j)=   (fin(i,2*j-1)+fin(i,n+3-2*j)) &
               +      (s*(fin(i,2*j-1)-fin(i,n+3-2*j)) &
               +       c*(fin(i,2*j  )+fin(i,n+4-2*j)))
          b(i,j     )=(c*(fin(i,2*j-1)-fin(i,n+3-2*j)) &
               -       s*(fin(i,2*j  )+fin(i,n+4-2*j)))&
               +         (fin(i,2*j  )-fin(i,n+4-2*j))
          b(i,nh+2-j)=(c*(fin(i,2*j-1)-fin(i,n+3-2*j)) &
               -       s*(fin(i,2*j  )+fin(i,n+4-2*j)))&
               -         (fin(i,2*j  )-fin(i,n+4-2*j))
       END DO
    END DO
    IF ( (nh>=2) .AND. (MOD(nh,2)==0) ) THEN
       !CDIR NODEP
       DO i = 1, lot
          a(i,nh/2+1)= 2.0_r8*fin(i,nh+1)
          b(i,nh/2+1)=-2.0_r8*fin(i,nh+2)
       END DO
    END IF
  END SUBROUTINE SplitFourTrans




  !JoinGrid: Merge fundamental algorithm complex output into 
  !          sequences of real numbers



  SUBROUTINE JoinGridTrans (a, b, fout, ldout, tdout, nh, lot)
    INTEGER, INTENT(IN ) :: ldout, tdout
    INTEGER, INTENT(IN ) :: nh
    INTEGER, INTENT(IN ) :: lot
    REAL(KIND=r8),    INTENT(OUT) :: fout(ldout,tdout)
    REAL(KIND=r8),    INTENT(IN ) :: a   (lot  ,nh)
    REAL(KIND=r8),    INTENT(IN ) :: b   (lot  ,nh)

    INTEGER :: i, j

    DO j = 1, nh
       DO i = 1, lot
          fout(i,2*j-1)=a(i,j)
          fout(i,2*j  )=b(i,j)
       END DO
    END DO

    !    fout(:,2*nh+1:tdout)=0.0_r8
    !    fout(lot+1:ldout, :)=0.0_r8
  END SUBROUTINE JoinGridTrans


  !OnePass: single pass of fundamental algorithm



  SUBROUTINE OnePass (a, b, c, d, lot, nh, ifac, la, Trigs, nTrigs)
    INTEGER, INTENT(IN ) :: lot
    INTEGER, INTENT(IN ) :: nh         ! = PROD(factor(1:K))
    REAL(KIND=r8),    INTENT(IN ) :: a(lot,nh)
    REAL(KIND=r8),    INTENT(IN ) :: b(lot,nh)
    REAL(KIND=r8),    INTENT(OUT) :: c(lot,nh)
    REAL(KIND=r8),    INTENT(OUT) :: d(lot,nh)
    INTEGER, INTENT(IN ) :: ifac       ! = factor(k)
    INTEGER, INTENT(IN ) :: la         ! = PROD(factor(1:k-1))
    INTEGER, INTENT(IN ) :: nTrigs
    REAL(KIND=r8),    INTENT(IN ) :: Trigs(nTrigs)

    INTEGER :: m
    INTEGER :: jump
    INTEGER :: i, j, k
    INTEGER :: ia, ja
    INTEGER :: ib, jb, kb
    INTEGER :: ic, jc, kc
    INTEGER :: id, jd, kd
    INTEGER :: ie, je, ke
    REAL(KIND=r8)    :: c1, s1
    REAL(KIND=r8)    :: c2, s2
    REAL(KIND=r8)    :: c3, s3
    REAL(KIND=r8)    :: c4, s4
    REAL(KIND=r8)    :: wka, wkb
    REAL(KIND=r8)    :: wksina, wksinb
    REAL(KIND=r8)    :: wkaacp, wkbacp
    REAL(KIND=r8)    :: wkaacm, wkbacm

    m=nh/ifac
    jump=(ifac-1)*la

    ia=  0
    ib=  m
    ic=2*m
    id=3*m
    ie=4*m

    ja=  0
    jb=  la
    jc=2*la
    jd=3*la
    je=4*la

    IF (ifac == 2) THEN
       DO j = 1, la
          !CDIR NODEP
          DO i = 1, lot
             c(i,j+ja)=a(i,j+ia)+a(i,j+ib)
             c(i,j+jb)=a(i,j+ia)-a(i,j+ib)
             d(i,j+ja)=b(i,j+ia)+b(i,j+ib)
             d(i,j+jb)=b(i,j+ia)-b(i,j+ib)
          END DO
       END DO
       DO k = la, m-1, la
          kb=k+k
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          ja=ja+jump
          jb=jb+jump
          DO j = k+1, k+la
             !CDIR NODEP
             DO i = 1, lot
                wka      =a(i,j+ia)-a(i,j+ib)
                c(i,j+ja)=a(i,j+ia)+a(i,j+ib)
                wkb      =b(i,j+ia)-b(i,j+ib)
                d(i,j+ja)=b(i,j+ia)+b(i,j+ib)
                c(i,j+jb)=c1*wka-s1*wkb
                d(i,j+jb)=s1*wka+c1*wkb
             END DO
          END DO
       END DO
    ELSEIF (ifac == 3) THEN
       DO j = 1, la
          !CDIR  NODEP
          DO i = 1, lot
             wka      =       a(i,j+ib)+a(i,j+ic)
             wksina   =sin60*(a(i,j+ib)-a(i,j+ic))
             wkb      =       b(i,j+ib)+b(i,j+ic)
             wksinb   =sin60*(b(i,j+ib)-b(i,j+ic))
             c(i,j+ja)=       a(i,j+ia)+wka
             c(i,j+jb)=      (a(i,j+ia)-0.5_r8*wka)-wksinb
             c(i,j+jc)=      (a(i,j+ia)-0.5_r8*wka)+wksinb
             d(i,j+ja)=       b(i,j+ia)+wkb
             d(i,j+jb)=      (b(i,j+ia)-0.5_r8*wkb)+wksina
             d(i,j+jc)=      (b(i,j+ia)-0.5_r8*wkb)-wksina
          END DO
       END DO
       DO k = la, m-1, la
          kb=k+k
          kc=kb+kb
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          c2=Trigs(kc+1)
          s2=Trigs(kc+2)
          ja=ja+jump
          jb=jb+jump
          jc=jc+jump
          DO j = k+1, k+la
             !CDIR NODEP
             DO i = 1, lot
                wka      =       a(i,j+ib)+a(i,j+ic)
                wksina   =sin60*(a(i,j+ib)-a(i,j+ic))
                wkb      =       b(i,j+ib)+b(i,j+ic)
                wksinb   =sin60*(b(i,j+ib)-b(i,j+ic))
                c(i,j+ja)=       a(i,j+ia)+wka
                d(i,j+ja)=       b(i,j+ia)+wkb
                c(i,j+jb)=c1*  ((a(i,j+ia)-0.5_r8*wka)-wksinb) &
                     -    s1*  ((b(i,j+ia)-0.5_r8*wkb)+wksina)
                d(i,j+jb)=s1*  ((a(i,j+ia)-0.5_r8*wka)-wksinb) &
                     +    c1*  ((b(i,j+ia)-0.5_r8*wkb)+wksina)
                c(i,j+jc)=c2*  ((a(i,j+ia)-0.5_r8*wka)+wksinb) &
                     -    s2*  ((b(i,j+ia)-0.5_r8*wkb)-wksina)
                d(i,j+jc)=s2*  ((a(i,j+ia)-0.5_r8*wka)+wksinb) &
                     +    c2*  ((b(i,j+ia)-0.5_r8*wkb)-wksina)
             END DO
          END DO
       END DO
    ELSEIF (ifac == 4) THEN
       DO j = 1, la
          !CDIR NODEP
          DO i = 1, lot
             wkaacp   =        a(i,j+ia)+a(i,j+ic)
             wkaacm   =        a(i,j+ia)-a(i,j+ic)
             wkbacp   =        b(i,j+ia)+b(i,j+ic)
             wkbacm   =        b(i,j+ia)-b(i,j+ic)
             c(i,j+ja)=wkaacp+(a(i,j+ib)+a(i,j+id))
             c(i,j+jc)=wkaacp-(a(i,j+ib)+a(i,j+id))
             d(i,j+jb)=wkbacm+(a(i,j+ib)-a(i,j+id))
             d(i,j+jd)=wkbacm-(a(i,j+ib)-a(i,j+id))
             d(i,j+ja)=wkbacp+(b(i,j+ib)+b(i,j+id))
             d(i,j+jc)=wkbacp-(b(i,j+ib)+b(i,j+id))
             c(i,j+jb)=wkaacm-(b(i,j+ib)-b(i,j+id))
             c(i,j+jd)=wkaacm+(b(i,j+ib)-b(i,j+id))
          END DO
       END DO
       DO k = la, m-1, la
          kb=k+k
          kc=kb+kb
          kd=kc+kb
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          c2=Trigs(kc+1)
          s2=Trigs(kc+2)
          c3=Trigs(kd+1)
          s3=Trigs(kd+2)
          ja=ja+jump
          jb=jb+jump
          jc=jc+jump
          jd=jd+jump
          DO j = k+1, k+la
             !CDIR NODEP
             DO i = 1, lot
                wkaacp   =            a(i,j+ia)+a(i,j+ic)
                wkbacp   =            b(i,j+ia)+b(i,j+ic)
                wkaacm   =            a(i,j+ia)-a(i,j+ic)
                wkbacm   =            b(i,j+ia)-b(i,j+ic)
                c(i,j+ja)=    wkaacp+(a(i,j+ib)+a(i,j+id))
                d(i,j+ja)=    wkbacp+(b(i,j+ib)+b(i,j+id))
                c(i,j+jc)=c2*(wkaacp-(a(i,j+ib)+a(i,j+id))) &
                     -    s2*(wkbacp-(b(i,j+ib)+b(i,j+id))) 
                d(i,j+jc)=s2*(wkaacp-(a(i,j+ib)+a(i,j+id))) &
                     +    c2*(wkbacp-(b(i,j+ib)+b(i,j+id)))
                c(i,j+jb)=c1*(wkaacm-(b(i,j+ib)-b(i,j+id))) &
                     -    s1*(wkbacm+(a(i,j+ib)-a(i,j+id)))
                d(i,j+jb)=s1*(wkaacm-(b(i,j+ib)-b(i,j+id))) &
                     +    c1*(wkbacm+(a(i,j+ib)-a(i,j+id)))
                c(i,j+jd)=c3*(wkaacm+(b(i,j+ib)-b(i,j+id))) &
                     -    s3*(wkbacm-(a(i,j+ib)-a(i,j+id)))
                d(i,j+jd)=s3*(wkaacm+(b(i,j+ib)-b(i,j+id))) &
                     +    c3*(wkbacm-(a(i,j+ib)-a(i,j+id)))
             END DO
          END DO
       END DO
    ELSEIF (ifac == 5) THEN
       DO j = 1, la
          !CDIR NODEP
          DO i = 1, lot
             c(i,j+ja)=a(i,j+ia)+(a(i,j+ib)+a(i,j+ie))+(a(i,j+ic)+a(i,j+id))
             d(i,j+ja)=    b(i,j+ia)&
                  +       (b(i,j+ib)+b(i,j+ie))&
                  +       (b(i,j+ic)+b(i,j+id))
             c(i,j+jb)=   (a(i,j+ia)&
                  + cos72*(a(i,j+ib)+a(i,j+ie))&
                  - cos36*(a(i,j+ic)+a(i,j+id)))&
                  -(sin72*(b(i,j+ib)-b(i,j+ie))&
                  + sin36*(b(i,j+ic)-b(i,j+id)))
             c(i,j+je)=   (a(i,j+ia)&
                  + cos72*(a(i,j+ib)+a(i,j+ie))&
                  - cos36*(a(i,j+ic)+a(i,j+id)))&
                  +(sin72*(b(i,j+ib)-b(i,j+ie))&
                  + sin36*(b(i,j+ic)-b(i,j+id)))
             d(i,j+jb)=   (b(i,j+ia)&
                  + cos72*(b(i,j+ib)+b(i,j+ie))&
                  - cos36*(b(i,j+ic)+b(i,j+id)))&
                  +(sin72*(a(i,j+ib)-a(i,j+ie))&
                  + sin36*(a(i,j+ic)-a(i,j+id)))
             d(i,j+je)=   (b(i,j+ia)&
                  + cos72*(b(i,j+ib)+b(i,j+ie))&
                  - cos36*(b(i,j+ic)+b(i,j+id)))&
                  -(sin72*(a(i,j+ib)-a(i,j+ie))&
                  + sin36*(a(i,j+ic)-a(i,j+id)))
             c(i,j+jc)=   (a(i,j+ia)&
                  - cos36*(a(i,j+ib)+a(i,j+ie))&
                  + cos72*(a(i,j+ic)+a(i,j+id)))&
                  -(sin36*(b(i,j+ib)-b(i,j+ie))&
                  - sin72*(b(i,j+ic)-b(i,j+id)))
             c(i,j+jd)=   (a(i,j+ia)&
                  - cos36*(a(i,j+ib)+a(i,j+ie))&
                  + cos72*(a(i,j+ic)+a(i,j+id)))&
                  +(sin36*(b(i,j+ib)-b(i,j+ie))&
                  - sin72*(b(i,j+ic)-b(i,j+id)))
             d(i,j+jc)=   (b(i,j+ia)&
                  - cos36*(b(i,j+ib)+b(i,j+ie))&
                  + cos72*(b(i,j+ic)+b(i,j+id)))&
                  +(sin36*(a(i,j+ib)-a(i,j+ie))&
                  - sin72*(a(i,j+ic)-a(i,j+id)))
             d(i,j+jd)=   (b(i,j+ia)&
                  - cos36*(b(i,j+ib)+b(i,j+ie))&
                  + cos72*(b(i,j+ic)+b(i,j+id)))&
                  -(sin36*(a(i,j+ib)-a(i,j+ie))&
                  - sin72*(a(i,j+ic)-a(i,j+id)))
          END DO
       END DO
       DO k = la, m-1, la
          kb=k+k
          kc=kb+kb
          kd=kc+kb
          ke=kd+kb
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          c2=Trigs(kc+1)
          s2=Trigs(kc+2)
          c3=Trigs(kd+1)
          s3=Trigs(kd+2)
          c4=Trigs(ke+1)
          s4=Trigs(ke+2)
          ja=ja+jump
          jb=jb+jump
          jc=jc+jump
          jd=jd+jump
          je=je+jump
          DO j = k+1, k+la
             !CDIR NODEP
             DO i = 1, lot
                c(i,j+ja)=     a(i,j+ia)&
                     +        (a(i,j+ib)+a(i,j+ie))&
                     +        (a(i,j+ic)+a(i,j+id))
                d(i,j+ja)=     b(i,j+ia)&
                     +        (b(i,j+ib)+b(i,j+ie))&
                     +        (b(i,j+ic)+b(i,j+id))
                c(i,j+jb)=c1*((a(i,j+ia)&
                     +  cos72*(a(i,j+ib)+a(i,j+ie))&
                     -  cos36*(a(i,j+ic)+a(i,j+id)))&
                     - (sin72*(b(i,j+ib)-b(i,j+ie))&
                     +  sin36*(b(i,j+ic)-b(i,j+id))))&
                     -    s1*((b(i,j+ia)&
                     +  cos72*(b(i,j+ib)+b(i,j+ie))&
                     -  cos36*(b(i,j+ic)+b(i,j+id)))&
                     + (sin72*(a(i,j+ib)-a(i,j+ie))&
                     +  sin36*(a(i,j+ic)-a(i,j+id))))
                d(i,j+jb)=s1*((a(i,j+ia)&
                     +  cos72*(a(i,j+ib)+a(i,j+ie))&
                     -  cos36*(a(i,j+ic)+a(i,j+id)))&
                     - (sin72*(b(i,j+ib)-b(i,j+ie))&
                     +  sin36*(b(i,j+ic)-b(i,j+id))))&
                     +    c1*((b(i,j+ia)&
                     +  cos72*(b(i,j+ib)+b(i,j+ie))&
                     -  cos36*(b(i,j+ic)+b(i,j+id)))&
                     + (sin72*(a(i,j+ib)-a(i,j+ie)) &
                     +  sin36*(a(i,j+ic)-a(i,j+id))))
                c(i,j+je)=c4*((a(i,j+ia)&
                     +  cos72*(a(i,j+ib)+a(i,j+ie))&
                     -  cos36*(a(i,j+ic)+a(i,j+id)))&
                     + (sin72*(b(i,j+ib)-b(i,j+ie)) &
                     +  sin36*(b(i,j+ic)-b(i,j+id)))) &
                     -    s4*((b(i,j+ia)&
                     +  cos72*(b(i,j+ib)+b(i,j+ie))&
                     -  cos36*(b(i,j+ic)+b(i,j+id)))&
                     - (sin72*(a(i,j+ib)-a(i,j+ie))&
                     +  sin36*(a(i,j+ic)-a(i,j+id))))
                d(i,j+je)=s4*((a(i,j+ia)&
                     +  cos72*(a(i,j+ib)+a(i,j+ie))&
                     -  cos36*(a(i,j+ic)+a(i,j+id)))&
                     + (sin72*(b(i,j+ib)-b(i,j+ie))&
                     +  sin36*(b(i,j+ic)-b(i,j+id))))&
                     +    c4*((b(i,j+ia)&
                     +  cos72*(b(i,j+ib)+b(i,j+ie))&
                     -  cos36*(b(i,j+ic)+b(i,j+id))) &
                     - (sin72*(a(i,j+ib)-a(i,j+ie))&
                     +  sin36*(a(i,j+ic)-a(i,j+id))))
                c(i,j+jc)=c2*((a(i,j+ia)&
                     -  cos36*(a(i,j+ib)+a(i,j+ie))&
                     +  cos72*(a(i,j+ic)+a(i,j+id)))&
                     - (sin36*(b(i,j+ib)-b(i,j+ie))&
                     -  sin72*(b(i,j+ic)-b(i,j+id))))&
                     -    s2*((b(i,j+ia)&
                     -  cos36*(b(i,j+ib)+b(i,j+ie))&
                     +  cos72*(b(i,j+ic)+b(i,j+id)))&
                     + (sin36*(a(i,j+ib)-a(i,j+ie))&
                     -  sin72*(a(i,j+ic)-a(i,j+id)))) 
                d(i,j+jc)=s2*((a(i,j+ia)&
                     -  cos36*(a(i,j+ib)+a(i,j+ie))&
                     +  cos72*(a(i,j+ic)+a(i,j+id)))&
                     - (sin36*(b(i,j+ib)-b(i,j+ie))&
                     -  sin72*(b(i,j+ic)-b(i,j+id))))&
                     +    c2*((b(i,j+ia)&
                     -  cos36*(b(i,j+ib)+b(i,j+ie))&
                     +  cos72*(b(i,j+ic)+b(i,j+id)))&
                     + (sin36*(a(i,j+ib)-a(i,j+ie))&
                     -  sin72*(a(i,j+ic)-a(i,j+id))))
                c(i,j+jd)=c3*((a(i,j+ia)&
                     -  cos36*(a(i,j+ib)+a(i,j+ie))&
                     +  cos72*(a(i,j+ic)+a(i,j+id)))&
                     + (sin36*(b(i,j+ib)-b(i,j+ie))&
                     -  sin72*(b(i,j+ic)-b(i,j+id))))&
                     -     s3*((b(i,j+ia)&
                     -  cos36*(b(i,j+ib)+b(i,j+ie))&
                     +  cos72*(b(i,j+ic)+b(i,j+id)))&
                     - (sin36*(a(i,j+ib)-a(i,j+ie))&
                     -  sin72*(a(i,j+ic)-a(i,j+id))))
                d(i,j+jd)=s3*((a(i,j+ia)&
                     -  cos36*(a(i,j+ib)+a(i,j+ie))&
                     +  cos72*(a(i,j+ic)+a(i,j+id)))&
                     + (sin36*(b(i,j+ib)-b(i,j+ie))&
                     -  sin72*(b(i,j+ic)-b(i,j+id))))&
                     +    c3*((b(i,j+ia)&
                     -  cos36*(b(i,j+ib)+b(i,j+ie))&
                     +  cos72*(b(i,j+ic)+b(i,j+id)))&
                     - (sin36*(a(i,j+ib)-a(i,j+ie))&
                     -  sin72*(a(i,j+ic)-a(i,j+id))))
             END DO
          END DO
       END DO
    ENDIF
  END SUBROUTINE OnePass






  !CreateFFT: Allocates and computes intermediate values used by
  !           the FFT for sequences of size nIn. If size
  !           is not in the form 2 * 2**i * 3**j * 5**k, with at
  !           least one of i,j,k/=0, stops and prints (stdout) the
  !           next possible size.



  SUBROUTINE CreateFFT(nIn, Factors, nFactors, Trigs, nTrigs)
    INTEGER, INTENT(IN)  :: nIn
    INTEGER, POINTER     :: Factors(:)
    INTEGER, INTENT(OUT)  :: nFactors
    REAL(KIND=r8),    POINTER     :: Trigs(:)
    INTEGER, INTENT(OUT)  :: nTrigs
    CHARACTER(LEN=15), PARAMETER :: h="**(CreateFFT)**" ! header
    CALL Factorize  (nIn, Factors,nBase,Base)
    nFactors = SIZE(Factors)
    CALL TrigFactors(nIn, Trigs)
    nTrigs = SIZE(Trigs)
  END SUBROUTINE CreateFFT



  !DestroyFFT: Dealocates input area, returning NULL pointers



  SUBROUTINE DestroyFFT(Factors, Trigs)
    INTEGER, POINTER :: Factors(:)
    REAL(KIND=r8),    POINTER :: Trigs(:)
    CHARACTER(LEN=16), PARAMETER :: h="**(DestroyFFT)**" ! header
    DEALLOCATE(Factors); NULLIFY(Factors)
    DEALLOCATE(Trigs  ); NULLIFY(Trigs  )
  END SUBROUTINE DestroyFFT



  !Factorize: Factorizes nIn/2 in powers of 4, 3, 2, 5, if possible.
  !           Otherwise, stops with error message



  SUBROUTINE Factorize (nIn, Factors,nBase,Base)
    INTEGER, INTENT(IN ) :: nIn
    ! As variaveis Base e nBase sao globais ao modulo, portanto seria desnecessaria
    ! a passagem por parametros. Entretanto para rodar com a chave de compilacao 
    ! "-openmp" ligada no TX7, essas variaveis sao parametros de entrada. (???)
    INTEGER, INTENT(IN ) :: nBase
    INTEGER, INTENT(IN ) :: Base(nBase)
    INTEGER, POINTER     :: Factors(:)
    CHARACTER(LEN=15), PARAMETER :: h="**(Factorize)**" ! header
    CHARACTER(LEN=15) :: c ! Character representation of integer
    INTEGER :: Powers(nBase)
    INTEGER :: nOut
    INTEGER :: sumPowers
    INTEGER :: ifac
    INTEGER :: i
    INTEGER :: j
    INTEGER :: left ! portion of nOut/2 yet to be factorized

    nOut= NextSizeFFT(nIn)

    IF (nIn /= nOut) THEN
       WRITE(c,"(i15)") nIn
       WRITE(nfprt,"(a,' FFT size = ',a,' not factorizable ')")&
            h, TRIM(ADJUSTL(c))
       WRITE(c,"(i15)") nOut
       WRITE(nfprt,"(a,' Next factorizable FFT size is ',a)")&
            h, TRIM(ADJUSTL(c))
       STOP
    END IF

    ! Loop over evens, starting from nOut, getting factors of nOut/2

    left = nOut/2
    Powers = 0

    ! factorize nOut/2

    DO i = 1, nBase 
       DO
          IF (MOD(left, Base(i)) == 0) THEN
             Powers(i) = Powers(i) + 1
             left = left / Base(i)
          ELSE
             EXIT
          END IF
       END DO
    END DO

    sumPowers=SUM(Powers)
    ALLOCATE (Factors(sumPowers+1))
    Factors(1)=sumPowers
    ifac = 1
    DO i = 1, nBase
       j = Permutation(i)
       Factors(ifac+1:ifac+Powers(j)) = Base(j)
       ifac = ifac + Powers(j)
    END DO
  END SUBROUTINE Factorize



  !TrigFactors: Sin and Cos required to compute FFT of size nIn



  SUBROUTINE TrigFactors (nIn, Trigs)
    INTEGER, INTENT(IN) :: nIn
    REAL(KIND=r8), POINTER       :: Trigs(:)
    INTEGER :: nn, nh, i
    REAL(KIND=r8)    :: pi, del, angle

    nn =  nIn  / 2
    nh = (nn+1) / 2
    ALLOCATE (Trigs(2*(nn+nh)))

    pi = 2.0_r8 * ASIN(1.0_r8)
    del = (2 * pi) / REAL(nn,r8)

    DO i = 1, 2*nn, 2
       angle = 0.5_r8 * REAL(i-1,r8) * del
       Trigs(i  ) = COS(angle)
       Trigs(i+1) = SIN(angle)
    END DO

    del = 0.5_r8 * del
    DO i = 1, 2*nh, 2
       angle = 0.5_r8 * REAL(i-1,r8) * del
       Trigs(2*nn+i  ) = COS(angle)
       Trigs(2*nn+i+1) = SIN(angle)
    END DO
  END SUBROUTINE TrigFactors



  !NextSizeFFT: Smallest integer >= input in the form 2 * 2**i * 3**j * 5**k



  FUNCTION NextSizeFFT(nIn) RESULT(nOUT)
    INTEGER, INTENT(IN ) :: nIn
    INTEGER              :: nOut
    CHARACTER(LEN=22), PARAMETER :: h="**(NextSizeFFT)**" ! header
    REAL(KIND=r8), PARAMETER :: limit = HUGE(nIn)-1   ! Maximum representable integer
    CHARACTER(LEN=15) :: charNIn ! Character representation of nIn
    INTEGER :: i 
    INTEGER :: left ! portion of nOut/2 yet to be factorized

    ! nOut = positive even 

    IF (nIn <= 0) THEN
       WRITE (charNIn,"(i15)") nIn
       WRITE(nfprt,"(a,' Meaningless FFT size='a)")h, TRIM(ADJUSTL(charNIn))
       STOP
    ELSE IF (MOD(nIn,2) == 0) THEN
       nOut = nIn
    ELSE
       nOut = nIn+1
    END IF

    ! Loop over evens, starting from nOut, looking for
    ! next factorizable even/2

    DO
       left = nOut/2

       ! factorize nOut/2

       DO i = 1, nBase 
          DO
             IF (MOD(left, Base(i)) == 0) THEN
                left = left / Base(i)
             ELSE
                EXIT
             END IF
          END DO
       END DO

       IF (left == 1) THEN
          EXIT
       ELSE IF (nOut < limit) THEN
          nOut = nOut + 2
       ELSE
          WRITE (charNIn,"(i15)") nIn
          WRITE(nfprt,"(a,' Next factorizable FFT size > ',a,&
               &' is not representable in this machine')")&
               h, TRIM(ADJUSTL(charNIn))
          STOP
       END IF
    END DO
  END FUNCTION NextSizeFFT



  !DirFFT: Computes direct FFT of 'lot' sequences of 'n' input real data
  !        as rows of 'fin', dimensioned 'fin(ldin,tdin)'. Input data is
  !        kept unchanged. Input values 'fin(n+1:ldin,:)' and 
  !        'fin(:,lot+1:tdin)' are not visited. 
  !        Outputs 'lot' sequences of 'n+1' real data
  !        as rows of 'fout', dimensioned 'fout(ldout,tdout)'. Output
  !        values 'fout(:,lot+1:tdout)' and 'fout(n+2:ldout,:)' are set to 0.



  SUBROUTINE DirFFTTrans (fInOut, ldInOut, tdInOut, n, lot, Trigs, nTrigs, Factors, nFactors)
    INTEGER, INTENT(IN ) :: ldInOut, tdInOut
    REAL(KIND=r8),    INTENT(INOUT) :: fInOut (ldInOut ,tdInOut)
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(IN ) :: lot
    INTEGER, INTENT(IN ) :: nTrigs
    REAL(KIND=r8),    INTENT(IN ) :: Trigs(nTrigs)
    INTEGER, INTENT(IN ) :: nFactors
    INTEGER, INTENT(IN ) :: Factors(nFactors)

    INTEGER :: nh
    INTEGER :: nfax
    INTEGER :: la
    INTEGER :: k
    LOGICAL :: ab2cd
    CHARACTER(LEN=12), PARAMETER :: h="**(Dir)**"
    REAL(KIND=r8) :: a(lot,n/2)
    REAL(KIND=r8) :: b(lot,n/2)
    REAL(KIND=r8) :: c(lot,n/2)
    REAL(KIND=r8) :: d(lot,n/2)

    nfax=Factors(1)
    nh=n/2

    CALL SplitGridTrans (fInOut, a, b, ldInOut, tdInOut, nh, lot)

    la=1
    ab2cd=.TRUE.
    DO k=1,nfax
       IF (ab2cd) THEN
          CALL OnePass (a, b, c, d, lot, nh, Factors(k+1), la, Trigs, nTrigs)
          ab2cd=.FALSE.
       ELSE
          CALL OnePass (c, d, a, b, lot, nh, Factors(k+1), la, Trigs, nTrigs)
          ab2cd=.TRUE.
       END IF
       la=la*Factors(k+1)
    END DO

    IF (ab2cd) THEN
       CALL JoinFourTrans (a, b, fInOut, ldInOut, tdInOut, n, nh, lot, Trigs, nTrigs)
    ELSE
       CALL JoinFourTrans (c, d, fInOut, ldInOut, tdInOut, n, nh, lot, Trigs, nTrigs)
    END IF
  END SUBROUTINE DirFFTTrans



  !SplitGrid: Split space domain real input into complex pairs to
  !           feed fundamental algorithm



  SUBROUTINE SplitGridTrans (fin, a, b, ldin, tdin, nh, lot)
    INTEGER, INTENT(IN ) :: ldin
    INTEGER, INTENT(IN ) :: tdin
    INTEGER, INTENT(IN ) :: nh
    INTEGER, INTENT(IN ) :: lot
    REAL(KIND=r8),    INTENT(IN ) :: fin(ldin, tdin)
    REAL(KIND=r8),    INTENT(OUT) :: a  (lot , nh)
    REAL(KIND=r8),    INTENT(OUT) :: b  (lot , nh)

    INTEGER :: i, j

    DO j = 1, nh
       DO i = 1, lot
          a(i,j)=fin(i,2*j-1)
          b(i,j)=fin(i,2*j  )
       END DO
    END DO

  END SUBROUTINE SplitGridTrans



  !JoinFour: Unscramble frequency domain complex output from fundamental
  !          algorithm into real sequences



  SUBROUTINE JoinFourTrans (a, b, fout, ldout, tdout, n, nh, lot, Trigs, nTrigs)
    INTEGER, INTENT(IN ) :: ldout, tdout
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(IN ) :: nh
    INTEGER, INTENT(IN ) :: lot
    REAL(KIND=r8),    INTENT(OUT) :: fout(ldout,tdout)
    REAL(KIND=r8),    INTENT(IN ) :: a   (lot  ,nh)
    REAL(KIND=r8),    INTENT(IN ) :: b   (lot  ,nh)
    INTEGER, INTENT(IN ) :: nTrigs
    REAL(KIND=r8),    INTENT(IN ) :: Trigs(nTrigs)

    INTEGER :: i, j
    REAL(KIND=r8)    :: scale, scalh
    REAL(KIND=r8)    :: c, s

    scale=1.0_r8/REAL(n,r8)
    scalh=0.5_r8*scale

    !cdir nodep
    DO i = 1, lot
       fout(i,  1)=scale*(a(i,1)+b(i,1))
       fout(i,n+1)=scale*(a(i,1)-b(i,1))
       fout(i,  2)=0.0_r8
    END DO

    DO j = 2, (nh+1)/2
       c=Trigs(n+2*j-1)
       s=Trigs(n+2*j  )
       !cdir nodep
       DO i = 1, lot
          fout(i,  2*j-1)=scalh*(   (a(i,j     )+a(i,nh+2-j)) &
               +                 (c*(b(i,j     )+b(i,nh+2-j))&
               +                  s*(a(i,j     )-a(i,nh+2-j))))
          fout(i,n+3-2*j)=scalh*(   (a(i,j     )+a(i,nh+2-j)) &
               -                 (c*(b(i,j     )+b(i,nh+2-j))&
               +                  s*(a(i,j     )-a(i,nh+2-j))))
          fout(i,    2*j)=scalh*((c*(a(i,j     )-a(i,nh+2-j))&
               -                  s*(b(i,j     )+b(i,nh+2-j)))&
               +                    (b(i,nh+2-j)-b(i,j     )))
          fout(i,n+4-2*j)=scalh*((c*(a(i,j     )-a(i,nh+2-j))&
               -                  s*(b(i,j     )+b(i,nh+2-j)))&
               -                    (b(i,nh+2-j)-b(i,j     )))
       END DO
    END DO

    IF ((nh>=2) .AND. (MOD(nh,2)==0)) THEN
       !cdir nodep
       DO i = 1, lot
          fout(i,nh+1)= scale*a(i,nh/2+1)
          fout(i,nh+2)=-scale*b(i,nh/2+1)
       END DO
    END IF

  END SUBROUTINE JoinFourTrans


SUBROUTINE Clear_Transform()

 IF (ALLOCATED(Spec))  DEALLOCATE ( Spec)
 IF (ALLOCATED(Four))  DEALLOCATE ( Four)
 IF (ALLOCATED(mnodes)) DEALLOCATE ( mnodes)
 IF (ALLOCATED(requests)) DEALLOCATE ( requests)
 IF (ALLOCATED(requestr)) DEALLOCATE ( requestr)
 IF (ALLOCATED(status)) DEALLOCATE ( status)
 IF (ALLOCATED(stat)) DEALLOCATE ( stat)

 IF (ALLOCATED(surfSpec)) DEALLOCATE ( surfSpec )! TRUE iff Surface Spectral Field
 IF (ALLOCATED(prevSpec)) DEALLOCATE ( prevSpec )! prior to first real vertical of this field at all internal arrays
 IF (ALLOCATED(Spec1d)) DEALLOCATE ( Spec1d )    ! points to Surface Spectral Field
 IF (ALLOCATED(Spec2d)) DEALLOCATE ( Spec2d )    ! points to Full Spectral Field

 IF (ALLOCATED(surfGrid)) DEALLOCATE ( surfGrid )! TRUE iff Surface Grid Field
 IF (ALLOCATED(prevGrid)) DEALLOCATE ( prevGrid )! prior to first real vertical of this field at all internal arrays
 IF (ALLOCATED(Grid2d)) DEALLOCATE ( Grid2d )    ! points to Surface Grid Field
 IF (ALLOCATED(Grid3d)) DEALLOCATE ( Grid3d )    ! points to Full Grid Field
 IF (ALLOCATED(fieldForDelLam)) DEALLOCATE ( fieldForDelLam )  ! TRUE iff this field position stores Lambda Derivative

 IF (ALLOCATED(prevVertDelLamSource)) DEALLOCATE ( prevVertDelLamSource )  ! source of Lambda Derivative


 IF (ALLOCATED(nEven)) DEALLOCATE ( nEven )
 IF (ALLOCATED(dnEven)) DEALLOCATE ( dnEven )
 IF (ALLOCATED(firstNEven)) DEALLOCATE ( firstNEven )
 IF (ALLOCATED(nOdd)) DEALLOCATE ( nOdd )
 IF (ALLOCATED(dnOdd)) DEALLOCATE ( dnOdd )
 IF (ALLOCATED(firstNOdd)) DEALLOCATE ( firstNOdd )
 IF (ALLOCATED(lmnExtMap)) DEALLOCATE ( lmnExtMap )
 IF (ALLOCATED(lmnMap)) DEALLOCATE ( lmnMap )
 IF (ALLOCATED(lmnZero)) DEALLOCATE ( lmnZero )


 IF (ALLOCATED(previousJ)) DEALLOCATE ( previousJ )


 IF (ALLOCATED(BlockFFT)) DEALLOCATE ( BlockFFT )
 IF (ALLOCATED(BlockFFT)) DEALLOCATE ( BlockFFT )
 IF (ALLOCATED(LS2F)) DEALLOCATE ( LS2F )
 IF (ALLOCATED(LF2S)) DEALLOCATE ( LF2S )
 IF (ALLOCATED(consIm)) DEALLOCATE ( consIm )
 IF (ALLOCATED(consRe)) DEALLOCATE ( consRe )
END SUBROUTINE Clear_Transform



END MODULE Transform






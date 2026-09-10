!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Transform </br></br>
!#
!# **Brief**: LEGENDRE AND FOURIER TRANSFORMS </br>
!# 
!# PURPOSE:
!# <ul type="disc">
!#  <li>Transforms a set of fields from spectral reprentation to grid point
!# representation and backwards. </li>
!#  <li>Computes derivatives on the zonal direction (longitudes). </li>
!# </ul>
!#  
!# IMPLEMENTATION GUIDELINES:
!# <ul type="disc">
!#  <li>Legendre Transform in matrix form. </li>
!#  <li>Efficiency on vector machines requires large matrices. </li>
!#  <li>Fields are packed together to generate large matrices prior to the transform. </li>
!#  <li>Version implements transform for Gaussian and Reduced Grids. </li>
!#  <li>OpenMP parallelism. </li>
!#  <li>MPI-ready data structure. </li>
!# </ul>
!# 
!# MODULE CREATION:
!# <ul type="disc">
!#  <li>InitTransform should be invoked once, prior to any other module routine. </li>
!#  <li>It needed to be invoked only once at each execution. </li>
!#  <li>It gathers model truncation and formulation info from the environment. </li>
!#  <li>Due to that, truncation cannot be changed during model execution. </li>
!# </ul>
!# 
!# GRID TO SPECTRAL USER INTERFACE: </br>
!# User interface contains 4 routines on the grid to spectral transform, with
!# the following semantics:
!# 
!# <ul type="disc">
!#  <li>CreateGridToSpec, <ul type="disc"> that prepares internal data structure (arrays) to perform
!# Fourier and then Legendre Transforms for a given number (input arguments) of
!# surface and full fields. </ul></li>
!#
!#  <li>DepositGridToSpec, <ul type="disc"> that gives the grid field to be transformed and the
!# spectral field that will receive the transformed field. Transform is not
!# performed by this call; it only specifies where is input data is and where
!# output data will be posted. </ul></li>
!#
!#  <li>DoGridToSpec, <ul type="disc"> that performs the transform over all previously deposited 
!#               grid fields, since last CreateGridToSpec, saving output
!#               data on the pointed positions. The number of surface grid fields
!#               and full grid fields deposited should match the information 
!#               given at CreateGridToSpec. At the end of DoGridToSpec, all 
!#               Spectral Fields informed by DepositGridToSpec are filled with 
!#               required information. </ul></li>
!#
!#  <li>DestroyGridToSpec, <ul type="disc"> that removes the internal data structure created by 
!#                    CreateGridToSpec.  </ul></li>
!# </ul>
!# 
!# GRID TO SPECTRAL USER INTERFACE ORDER OF INVOCATION: 
!# <ul type="disc">
!#  <li>CreateGridToSpec should be invoked once, prior to any other routine. </li>
!#  <li>Use one DepositGridToSpec for each Grid field to be transformed. </li>
!#  <li>Once all Deposits are done, invoke DoGridToSpec. </li>
!#  <li>Finalize by invoking DestroyGridToSpec. </li>
!# </ul>
!# 
!# SPECTRAL TO GRID USER INTERFACE: 
!# <ul type="disc">
!#  <li>Similar to Grid to Spec: <ul type="disc"> enlarged to include zonal derivative computation. </ul></li>
!#  <li>CreateSpecToGrid: <ul type="disc"> Similar functionality, except that the number of output grid 
!#                   fields can be larger that the number of input spectral 
!#                   fields, to accomodate zonal derivatives. </ul></li>
!#  <li>DepositSpecToGrid: <ul type="disc"> Use only to compute transform of input spectral field,
!#                    without derivatives. </ul></li>
!#  <li>DepositSpecToDelLamGrid: <ul type="disc"> Use only to compute the zonal derivative 
!#                    of input spectral field.  </ul></li>
!#  <li>DepositSpecToGridAndDelLamGrid: <ul type="disc"> Use to compute both the transform and 
!#                    the zonal derivative of input spectral field.  </ul></li>
!#  <li>DoSpecToGrid: <ul type="disc"> Perform all required transforms and zonal derivatives 
!#               simultaneosly.  </ul></li>
!#  <li>DestroySpecToGrid: <ul type="disc"> Similar functionality. </ul></li>
!# </ul>
!# 
!# SPECTRAL TO GRID USER INTERFACE ORDER OF INVOCATION: 
!# <ul type="disc">
!#  <li>CreateSpecToGrid should be invoked once, prior to any other routine. </li>
!#  <li>Use one of DepositSpecToGrid, DepositSpecToDelLamGrid ou DepositSpecToGridAndDelLamGrid 
!# for each Grid field to be transformed, with or without zonal derivative. </li>
!#  <li>Once all Deposits are done, invoke DoSpecToGrid. </li>
!#  <li>Finalize by invoking DestroySpecToGrid. </li>
!# </ul>
!# 
!# OTHER FUNCTIONALITIES: </br>
!# <ul type="disc">
!#  <li>NextSizeFFT returns the smallest size of FFT equal or larger than input argument. </li>
!#  <li>Usefull for array dimensions. </li>
!# </ul>  </br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!# 
!# **Author**: Paulo Kubota </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>20-04-2010 - Paulo Kubota - version: 1.10.0 </li>
!#  <li>26-04-2019 - Denis Eiras  - version: 2.0.0 - some adaptations for modularizing Chopping </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endtodo
!#
!# @documentation
!#
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2> </br>
!# **&#9993;**<mailto:atende.cptec@inpe.br> </br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---


module Mod_Transform

  use Mod_InputParameters, only : &
    EMRad1

  use Mod_Utils, only : &
    LegFuncS2F, &
    GaussWeights, &
    NoBankConflict, &
    allpolynomials

  use Mod_Sizes, only : &
    mMax, &
    nMax, &
    nExtMax, &
    mnMax, &
    mnExtMax, &
    mnExtMap, &
    nodeHasM, &
    nodeHasJ_f, &
    lm2m, &
    myMMax, &
    myMNMax, &
    myMNExtMax, &
    myMMap, &
    myNMap, &
    myMNMap, &
    myMExtMap, &
    myNExtMap, &
    myMNExtMap, &
    havesurf, &
    myfirstlev, &
    mylastlev, &
    myfirstlat_f, &
    mylastlat_f, &
    MsPerProc, &
    firstlatinproc_f, &
    lastlatinproc_f, &
    nlatsinproc_f, &
    messages_f, &
    messages_g, &
    messproc_f, &
    messproc_g, &
    nrecs_f, &
    nrecs_g, &
    kfirst_four, &
    klast_four, &
    myfirstlat_f, &
    myfirstlat, &
    mylastlat_f, &
    mylastlat, &
    myfirstlon, &
    mylastlon, &
    iMax, &
    jMax, &
    jMaxHalf, &
    kMaxloc, &
    ibMax, &
    ibMaxPerJB, &
    jbMax, &
    ibPerIJ, &
    jbPerIJ, &
    jMinPerM, &
    jMaxPerM, &
    iMaxPerJ, &
    myJMax_f, &
    JMaxlocal_f, &
    MMaxlocal, &
    ThreadDecomp

  use Mod_Parallelism_Group_Chopping, only : &
    myid &
    , maxnodes &
    , mpiCommGroup

  use Mod_Parallelism_Fourier, only : &
    myid_four &
    , mygroup_four &
    , COMM_FOUR &
    , maxNodes_four

  use Mod_Communications, only : &
    bufsend, &
    bufrec, &
    dimrecbuf, &
    dimsendbuf

  implicit none

  include 'mpif.h'
  include 'precision.h'
  include 'messages.h'

  private
  public :: InitTransform
  public :: CreateSpecToGrid
  public :: DepositSpecToGrid
  public :: DepositSpecToDelLamGrid
  public :: DepositSpecToGridAndDelLamGrid
  public :: DoSpecToGrid
  public :: DestroySpecToGrid
  public :: CreateGridToSpec
  public :: DepositGridToSpec
  public :: DepositGridToSpec_PK
  public :: DoGridToSpec
  public :: DestroyGridToSpec
  public :: NextSizeFFT
  public :: Clear_Transform

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
  ! lmnExtMap(lmn), lmn = 1,..., myMNExtMax maps the external representation
  ! of Spectral Extended Fields to the transform internal representation,
  ! for a single vertical
  ! 
  ! lmnMap(lmn), lmn = 1,..., myMNMax maps the external representation
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
  ! Array Grid(dvdlj, di), using Grid(nVertGrid*myJMax_f, iMax)dirInp = './'                        ! input data directory
  ! dvdlj and di avoid memory bank conflict. nVertGrid is the
  ! total number of verticals stored (sum over all fields,
  ! including lambda derivatives in Spectral to Grids transforms)
  ! 
  ! first dimension of Grid maps verticals of all fields
  ! and latitudes stored by this process into an 1D structure
  ! in array element order of:(Vertical, Field, latitude).
  ! 
  ! Second dimension of Grid is longitude.

  integer :: mGlob ! global variables for OpenMp
  integer :: nGlob
  integer :: iGlob
  integer :: jGlob
  integer :: ibGlob
  integer :: ksg
  integer :: kountg
  integer :: ipar2g
  integer :: ipar3g

  interface DepositSpecToGrid
    module procedure Deposit1D, Deposit2D
  end interface

  interface DepositSpecToDelLamGrid
    module procedure DepDL1D, DepDL2D
  end interface

  interface DepositSpecToGridAndDelLamGrid
    module procedure DepDLG1D, DepDLG2D
  end interface

  interface DestroySpecToGrid
    module procedure Destroy
  end interface

  interface DepositGridToSpec
    module procedure Deposit1D, Deposit2D
  end interface

  interface DepositGridToSpec_PK
    module procedure Deposit1D_pk, Deposit2D_PK
  end interface

  interface DestroyGridToSpec
    module procedure Destroy
  end interface

  real(kind = p_r8), allocatable :: Spec(:, :)
  real(kind = p_r8), allocatable :: Four(:, :)
  integer, allocatable :: mnodes(:)
  integer, allocatable :: requests(:)
  integer, allocatable :: requestr(:)
  integer, allocatable :: status(:)
  integer, allocatable :: stat(:, :)

  TYPE p1d
    real(kind = p_r8), pointer :: p(:)
  end TYPE p1d

  TYPE p2d
    real(kind = p_r8), pointer :: p(:, :)
  end TYPE p2d

  TYPE p3d
    real(kind = p_r8), pointer :: p(:, :, :)
  end TYPE p3d

  integer :: nSpecFields        
  !# total Spectral Fields in transform
  integer :: usedSpecFields     
  !# Spectral Fields currently deposited
  integer :: lastUsedSpecVert   
  !# at internal array Spec
  logical, allocatable :: surfSpec(:)        
  !# true iff Surface Spectral Field
  integer, allocatable :: prevSpec(:)        
  !# prior to first real vertical of this field at all internal arrays
  TYPE(p1d), allocatable :: Spec1d(:)          
  !# points to Surface Spectral Field
  TYPE(p2d), allocatable :: Spec2d(:)          
  !# points to Full Spectral Field

  integer :: nGridFields        
  !# total Grid Fields in transform
  integer :: usedGridFields     
  !# Grid Fields currently deposited
  integer :: lastUsedGridVert   
  !# at internal array Grid
  logical, allocatable :: surfGrid(:)        
  !# true iff Surface Grid Field
  integer, allocatable :: prevGrid(:)        
  !# prior to first real vertical of this field at all internal arrays
  TYPE(p2d), allocatable :: Grid2d(:)          
  !# points to Surface Grid Field
  TYPE(p3d), allocatable :: Grid3d(:)          
  !# points to Full Grid Field
  logical, allocatable :: fieldForDelLam(:)  
  !# true iff this field position stores Lambda Derivative

  logical :: willDelLam
  integer :: usedDelLamFields   
  !# last position at Grid array used for Lambda Der. fields
  integer :: lastUsedDelLamVert 
  !# at Grid array
  integer, allocatable :: prevVertDelLamSource(:)  
  !# source of Lambda Derivative

  integer :: dlmn
  integer, allocatable :: nEven(:)
  integer, allocatable :: dnEven(:)
  integer, allocatable :: firstNEven(:)
  integer, allocatable :: nOdd(:)
  integer, allocatable :: dnOdd(:)
  integer, allocatable :: firstNOdd(:)
  integer, allocatable :: lmnExtMap(:)
  integer, allocatable :: lmnMap(:)
  integer, allocatable :: lmnZero(:)

  integer :: nVertSpec
  integer :: nVertGrid
  integer :: nFull_g
  integer :: nSurf_g
  integer :: nFull_s
  integer :: nSurf_s
  integer :: dv
  integer :: dvdlj
  integer :: di
  integer :: dip1
  integer, allocatable :: previousJ(:)

  integer :: djh
  integer :: dvjh

  !  HIDDEN DATA, FFT SIZE INDEPENDENT:
  !  nBase=SIZE(Base)
  !  Base are the bases for factorization of n; base 4 should come
  !       before base 2 to improve FFT efficiency
  !  Permutation defines order of bases when computing FFTs
  !  sinXX, cosYY are constants required for computing FFTs

  integer, parameter :: nBase = 4
  integer, parameter :: Base(nBase) = (/ 4, 2, 3, 5 /)
  integer, parameter :: Permutation(nBase) = (/ 2, 3, 1, 4 /)
  real(kind = p_r8) :: sin60
  real(kind = p_r8) :: sin36
  real(kind = p_r8) :: sin72
  real(kind = p_r8) :: cos36
  real(kind = p_r8) :: cos72

  !  HIDDEN DATA, FFT SIZE DEPENDENT:
  !  For each FFT Size:
  !     nLong is FFT size;
  !     nFactors, Factors(:), nTrig, Trig(:) are size dependent constants;
  !     firstLat is the first FFT latitude for this block
  !     lastLat is the last FFT latitude for this block
  !  BlockFFT has size dependent data for all required FFTs;
  !  nBlockFFT is BlockFFT size

  TYPE MultiFFT
    integer :: nLong
    integer :: nFactors
    integer, pointer :: Factors(:)
    integer :: nTrigs
    real(kind = p_r8), pointer :: Trigs(:)
    integer :: firstLat
    integer :: lastLat
  end TYPE MultiFFT

  TYPE(MultiFFT), allocatable, target :: BlockFFT(:)
  integer :: nBlockFFT

  real(kind = p_r8), allocatable :: LS2F(:, :)
  real(kind = p_r8), allocatable :: LF2S(:, :)
  real(kind = p_r8), allocatable :: consIm(:)
  real(kind = p_r8), allocatable :: consRe(:)

  logical, parameter :: dumpLocal = .false.


  ! CALLING TREE (internal)
  !
  ! InitTransform ! MpiMappings
  !               ! CreateFFT    ! Factorize
  !                              ! TrigFactors
  ! CreateSpecToGrid
  !
  ! DepositSpecToGrid ! Deposit1D or Deposit2D
  !
  ! DepositSpecToDelLamGrid ! DepDL1D or DepDL2D
  !
  ! DepositSpecToGridAndDelLamGrid ! DepDLG1D or DepDLG2D
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
contains


  subroutine InitTransform()
    !# Initializes Transform
    !# ---
    !# @info
    !# **Brief:** 
    !# <ul type="disc">
    !#  <li>InitTransform should be invoked once, prior to any other module routine. </li>
    !#  <li>It needed to be invoked only once at each execution. </li>
    !#  <li>It gathers model truncation and formulation info from the environment. </li>
    !#  <li>Due to that, truncation cannot be changed during model execution. </li>
    !# </ul></br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer :: m
    integer :: n
    integer :: lm
    integer :: lmn
    integer :: i
    integer :: j
    real(kind = p_r8) :: radi
    integer :: iBlockFFT
    integer :: lFirst, lLast
    TYPE(MultiFFT), pointer :: p
    character(len = 10) :: c0, c1, c2, c3
    character(len = *), parameter :: h = "**(InitTransform)**"

    ! Mapping (lm,n) into array Spec, for transform;
    ! for each lm:
    !    map all even components of n,
    !    followed by all odd components of n.
    ! Limits of even and odd components at array Spec:
    ! nEven(myMMax): How many n with (n+m) even for this lm
    ! nOdd (myMMax): How many n with (n+m) odd  for this lm

    allocate(nEven(myMMax))
    allocate(nOdd(myMMax))
    do lm = 1, myMMax
      m = lm2m(lm)
      nEven(lm) = (mMax - m + 3) / 2
      nOdd(lm) = (mMax - m + 2) / 2
    end do

    ! Spec Bounds:
    ! dnEven(myMMax): nEven without memory bank conflicts
    ! dnOdd (myMMax): nOdd  without memory bank conflicts

    allocate(dnEven(myMMax))
    allocate(dnOdd(myMMax))
    dnEven = NoBankConflict(nEven)
    dnOdd = NoBankConflict(nOdd)

    ! for each lm, find:
    ! first n even, first n odd;
    ! all n even will be mapped into sequential positions,
    ! from firstNEven(lm)
    ! all n odd will be mapped into sequential positions,
    ! from firstNOdd(lm)

    allocate(firstNEven(myMMax))
    allocate(firstNOdd(myMMax))
    firstNEven(1) = 1
    firstNOdd(1) = nEven(1) + 1
    do lm = 2, myMMax
      firstNEven(lm) = firstNOdd(lm - 1) + nOdd(lm - 1)
      firstNOdd (lm) = firstNEven(lm) + nEven(lm)
    end do

    ! size of first dimension of array Spec:
    ! space for all even and odd values, avoiding bank conflict;
    ! garantees space for each lm matrix;

    dlmn = sum(nEven(1:myMMax)) + sum(nOdd(1:myMMax - 1)) + dnOdd(myMMax)
    dlmn = NoBankConflict(dlmn)

    ! lmnExtMap: maps external representation of Extended Spectral fields
    !            into internal representation for transforms (array Spec)

    allocate(lmnExtMap(myMNExtMax))
    do lmn = 1, myMNExtMax
      lm = myMExtMap(lmn)
      n = myNExtMap(lmn)
      m = lm2m(lm)
      if (mod(m + n, 2) == 0) then
        lmnExtMap(lmn) = firstNEven(lm) + (n - m) / 2
      else
        lmnExtMap(lmn) = firstNOdd(lm) + (n - m - 1) / 2
      end if
    end do

    ! lmnExtMap: maps external representation of Regular Spectral fields
    !            into internal representation for transforms (array Spec)

    allocate(lmnMap(myMNMax))
    do lmn = 1, myMNMax
      lm = myMMap(lmn)
      n = myNMap(lmn)
      m = lm2m(lm)
      if (mod(m + n, 2) == 0) then
        lmnMap(lmn) = firstNEven(lm) + (n - m) / 2
      else
        lmnMap(lmn) = firstNOdd(lm) + (n - m - 1) / 2
      end if
    end do

    ! lmnZero: null positions for Regular Spectral fields into
    !          internal representations for transforms

    allocate(lmnZero(myMMax))
    do lm = 1, myMMax
      lmnZero(lm) = lmnExtMap(myMNExtMap(lm, mMax + 1))
    end do

    ! dimensions

    di = NoBankConflict(iMax)
    dip1 = NoBankConflict(iMax + 1)
    djh = NoBankConflict(jMaxHalf)


    ! previousJ: how many latitudes are stored at array Four
    !            before latitude j

    allocate(previousJ(jMax))
    do j = myfirstlat_f, mylastlat_f
      previousJ(j) = j - myfirstlat_f
    end do

    ! consIm, consRe : values to be used in dellam

    allocate(consIm(mMax))
    allocate(consRe(mMax))
    do m = 1, mMax
      consIm(m) = EMRad1 * real(m - 1, p_r8)
      consRe(m) = - consIm(m)
    end do
    ! LS2F: Associated Legendre Functions (latitude, mn)
    !       for spectral to fourier computation

    allocate(LS2F(djh, dlmn))
    if (allpolynomials) then
      do lmn = 1, myMNExtMax
        lm = myMExtMap(lmn)
        n = myNExtMap(lmn)
        m = lm2m(lm)
        do j = 1, jMaxHalf
          LS2F(j, lmnExtMap(lmn)) = LegFuncS2F(j, mnExtMap(m, n))
        end do
      end do
    else
      do lmn = 1, myMNExtMax
        m = myMExtMap(lmn)
        n = myNExtMap(lmn)
        do j = 1, jMaxHalf
          LS2F(j, lmnExtMap(lmn)) = LegFuncS2F(j, mymnExtMap(m, n))
        end do
      end do
    end if
    LS2F(jMaxHalf + 1:djh, :) = 0.0_p_r8
    LS2F(:, myMNExtMax + 1:dlmn) = 0.0_p_r8

    ! LF2S: Associated Legendre Functions (latitude, mn)
    !       for fourier to spectral computation

    allocate(LF2S(dlmn, djh))
    do j = 1, jMaxHalf
      do lmn = 1, myMNExtMax
        LF2S(lmnExtMap(lmn), j) = LS2F(j, lmnExtMap(lmn)) * GaussWeights(j)
      end do
    end do
    LF2S(:, jMaxHalf + 1:djh) = 0.0_p_r8
    LF2S(myMNExtMax + 1:dlmn, :) = 0.0_p_r8

    ! FFT size independent data

    radi = atan(1.0_p_r8) / 45.0_p_r8
    sin60 = sin(60.0_p_r8 * radi)
    sin36 = sin(36.0_p_r8 * radi)
    sin72 = sin(72.0_p_r8 * radi)
    cos36 = cos(36.0_p_r8 * radi)
    cos72 = cos(72.0_p_r8 * radi)

    ! count FFT Blocks (one size FFT per block)

    nBlockFFT = 1
    do i = myfirstlat_f + 1, mylastlat_f
      if (iMaxPerJ(i).ne.iMaxPerJ(i - 1)) then
        nBlockFFT = nBlockFFT + 1
      end if
    end do

    ! Block FFT size
    ! Block FFT first Lat position on array Four
    ! Block FFT last Lat position on array Four

    allocate (BlockFFT(nBlockFFT))
    iBlockFFT = 1
    p => BlockFFT(iBlockFFT)
    p%nLong = iMaxPerJ(myfirstlat_f)
    p%firstLat = 1
    do j = myfirstlat_f + 1, mylastlat_f
      if (iMaxPerJ(j) /= p%nLong) then
        p%lastLat = j - myfirstlat_f
        iBlockFFT = iBlockFFT + 1
        p => BlockFFT(iBlockFFT)
        p%nLong = iMaxPerJ(j)
        p%firstLat = j + 1 - myfirstlat_f
      end if
    end do
    if (iBlockFFT /= nBlockFFT) then
      write(p_nfprt, "(a,' iBlockFFT (',i5,') /= nBlockFFT (',i5,')')") &
        h, iBlockFFT, nBLockFFT
      stop h
    else
      p%lastLat = mylastlat_f + 1 - myfirstlat_f
    end if

    ! fill internal FFT data for each FFT Block

    do iBlockFFT = 1, nBlockFFT
      p => BlockFFT(iBlockFFT)
      call CreateFFT(p%nLong, p%Factors, p%nFactors, p%Trigs, p%nTrigs)
    end do

    ! debug dumping

    if (dumpLocal) then
      write(p_nfprt, "(a,' starts dumping internal data')") h
      write(p_nfprt, "(a,' dlmn      =',i10)") h, dlmn
      write(p_nfprt, "(a,' di        =',i10)") h, di
      write(p_nfprt, "(a,' dip1      =',i10)") h, dip1
      write(p_nfprt, "(a,' djh       =',i10)") h, djh
      write(c0, "(i10)") nBlockFFT
      write(p_nfprt, "(a,' There are ',a,' FFT Blocks')") h, trim(adjustl(c0))
      write(p_nfprt, "(a,' Block #; FFT size(longitudes); latitudes')") h
      do iBlockFFT = 1, nBlockFFT
        p => BlockFFT(iBlockFFT)
        lFirst = p%firstLat
        lLast = p%lastLat
        if (lLast == lFirst + 1) then
          write(c0, "(i10)") (lFirst + 1) / 2
          write(c1, "(i10)") jMax - (lLast + 1) / 2 + 1
          write(p_nfprt, "(a,i8,i20,5x,a,' and ',a)") h, iBlockFFT, p%nLong, &
            trim(adjustl(c0)), trim(adjustl(c1))
        else
          write(c0, "(i10)") (lFirst + 1) / 2
          write(c1, "(i10)") (lLast + 1) / 2
          write(c2, "(i10)") jMax - (lFirst + 1) / 2 + 1
          write(c3, "(i10)") jMax - (lLast + 1) / 2 + 1
          write(p_nfprt, "(a,i8,i20,5x,a,':',a,' and ',a,':',a)") h, iBlockFFT, p%nLong, &
            trim(adjustl(c0)), trim(adjustl(c1)), trim(adjustl(c3)), trim(adjustl(c2))
        end if
      end do
    end if
  end subroutine InitTransform


  subroutine CreateSpecToGrid(nFullSpec, nSurfSpec, nFullGrid, nSurfGrid)
    !# Creates Spectral To Grid
    !# ---
    !# @info
    !# **Brief:** Similar functionality to CreateGridToSpec, except that the number
    !# of output grid fields can be larger that the number of input spectral
    !# fields, to accomodate zonal derivatives.  </br>    
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: nFullSpec
    integer, intent(in) :: nSurfSpec
    integer, intent(in) :: nFullGrid
    integer, intent(in) :: nSurfGrid
    character(len = *), parameter :: h = "**(CreateSpecToGrid)**"
    integer :: nsusp, nsugr

    if (havesurf) then
      nsusp = nSurfSpec
      nsugr = nSurfGrid
    else
      nsusp = 0
      nsugr = 0
    endif
    nfull_g = nFullGrid
    nsurf_g = nSurfGrid
    nfull_s = nFullSpec
    nsurf_s = nSurfSpec
    nSpecFields = nFullSpec + nsusp
    usedSpecFields = 0
    lastUsedSpecVert = 0
    allocate(surfSpec(nSpecFields))
    allocate(prevSpec(nSpecFields))
    allocate(Spec1d(nSpecFields))
    allocate(Spec2d(nSpecFields))
    nVertSpec = nFullSpec * kMaxloc + nsusp

    nGridFields = nFullGrid + nSurfGrid
    usedGridFields = 0
    lastUsedGridVert = 0
    allocate(surfGrid(nGridFields))
    allocate(prevGrid(nGridFields))
    allocate(Grid2d(nGridFields))
    allocate(Grid3d(nGridFields))
    allocate(fieldForDelLam(nGridFields))
    nVertGrid = nFullGrid * kMaxloc + nsugr

    willDelLam = .false.
    usedDelLamFields = nFullSpec + nSurfSpec
    lastUsedDelLamVert = nVertSpec
    allocate(prevVertDelLamSource(nGridFields))

    dv = NoBankConflict(2 * nVertSpec)
    dvjh = NoBankConflict(nVertSpec * jMaxHalf)
    dvdlj = NoBankConflict(nVertGrid * myJMax_f)

    allocate (Spec(dlmn, dv))
    allocate (Four(dvdlj, dip1))
    if (.not.allocated(requests)) then
      allocate(requests(0:maxnodes))
      allocate(requestr(0:maxnodes))
      allocate(status(MPI_STATUS_SIZE))
      allocate(stat(MPI_STATUS_SIZE, maxnodes))
    end if

    if (dumpLocal) then
      write(p_nfprt, "(a,' Spec: n, used, lastVert, nVert=',4i5)") &
        h, nSpecFields, usedSpecFields, lastUsedSpecVert, nVertSpec
      write(p_nfprt, "(a,' Grid: n, used, lastVert, nVert=',4i5)") &
        h, nGridFields, usedGridFields, lastUsedGridVert, nVertGrid
      write(p_nfprt, "(a,' DelLam: used, lastVert=',2i5)") &
        h, usedDelLamFields, lastUsedDelLamVert
    end if
  end subroutine CreateSpecToGrid


  subroutine Deposit1D_PK (ArgSpec, ArgGrid, ArgDelLam)
    !# Deposit 1D_PK
    !# ---
    !# @info
    !# **Brief:** Gives the grid field to be transformed and the spectral field
    !# that will receive the transformed field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    real(kind = p_r8), target, intent(in) :: ArgSpec(:)
    real(kind = p_r8), target, optional, intent(in) :: ArgGrid(:, :)
    real(kind = p_r8), target, optional, intent(in) :: ArgDelLam(:, :)
    character(len = *), parameter :: h = "**(Deposit1D)**"
    logical :: getGrid, getDelLam
    integer :: lusf
    integer :: pusv
    integer :: lusv

    integer :: lugf
    integer :: pugv
    integer :: lugv

    integer :: ludlf
    integer :: pudlv
    integer :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid = present(ArgGrid)
    getDelLam = present(ArgDelLam)

    ! critical section to update data structure counters

    if (.not.havesurf) then

      !
      !  this processor won't transform surface fields
      !  ( it only needs to know which are the correponding grid fields )
      !
      if (getGrid) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        surfGrid(lugf) = .true.
        prevGrid(lugf) = -1
        nullify(Grid3d(lugf)%p)
        Grid2d(lugf)%p => ArgGrid
        fieldForDelLam(lugf) = .false.
        prevVertDelLamSource(lugf) = -1
        if (getDelLam) then
          ludlf = usedDelLamFields + 1
          usedDelLamFields = ludlf
          surfGrid(ludlf) = .true.
          prevGrid(ludlf) = -1
          nullify(Grid3d(ludlf)%p)
          Grid2d(ludlf)%p => ArgDelLam
          fieldForDelLam(lugf) = .false.
          prevVertDelLamSource(lugf) = -1
        end if
      else if (getDelLam) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        surfGrid(lugf) = .true.
        prevGrid(lugf) = -1
        nullify(Grid3d(lugf)%p)
        Grid2d(lugf)%p => ArgDelLam
        fieldForDelLam(lugf) = .false.
        prevVertDelLamSource(lugf) = - 1
      else
        write(p_nfprt, "(a, ' no Grid or DelLam field required')") h
        stop h
      end if

    else

      lusf = usedSpecFields + 1
      usedSpecFields = lusf
      pusv = lastUsedSpecVert
      lusv = pusv + 1
      lastUsedSpecVert = lusv
      if (getGrid) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        pugv = lastUsedGridVert
        lugv = pugv + 1
        lastUsedGridVert = lugv
        if (getDelLam) then
          ludlf = usedDelLamFields + 1
          usedDelLamFields = ludlf
          pudlv = lastUsedDelLamVert
          ludlv = pudlv + 1
          lastUsedDelLamVert = ludlv
        end if
      else if (getDelLam) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        pugv = lastUsedGridVert
        lugv = pugv + 1
        lastUsedGridVert = lugv
      else
        write(p_nfprt, "(a, ' no Grid or DelLam field required')") h
        stop h
      end if

      ! deposit spectral field

      if (lusf <= nSpecFields) then
        surfSpec(lusf) = .true.
        prevSpec(lusf) = pusv
        nullify(Spec2d(lusf)%p)
        Spec1d(lusf)%p => ArgSpec
      else
        write(p_nfprt, "(a, ' too many spectral fields')") h
        write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
        stop h
      end if

      if (dumpLocal) then
        write(p_nfprt, "(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
          h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
      end if

      if (getGrid) then

        ! ArgGrid Field is required;
        ! Store ArgGrid Field info at first available position of the Grid array

        if (lugf <= nSpecFields) then
          surfGrid(lugf) = .true.
          prevGrid(lugf) = pugv
          nullify(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgGrid
          fieldForDelLam(lugf) = .false.
          prevVertDelLamSource(lugf) = -1
        else
          write(p_nfprt, "(a, ' too many grid fields')") h
          write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          stop h
        end if

        if (dumpLocal) then
          write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
            h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
            h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
        end if

        if (getDelLam) then

          ! ArgGrid and ArgDelLam are required;
          ! Store ArgDelLam info at first available position of the Grid array
          ! Mark field for DelLam

          ! TODO check pudlv not initialized
          willDelLam = .true.
          if (ludlf <= nGridFields) then
            surfGrid(ludlf) = .true.
            prevGrid(ludlf) = pudlv
            nullify(Grid3d(ludlf)%p)
            Grid2d(ludlf)%p => ArgDelLam
            fieldForDelLam(ludlf) = .true.
            prevVertDelLamSource(ludlf) = prevGrid(lugf)
          else
            write(p_nfprt, "(a, ' too many grid fields including DelLam')") h
            write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
            stop h
          end if

          if (dumpLocal) then
            write(p_nfprt, "(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
              h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
            write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
              h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
          end if
        end if

      else if (getDelLam) then

        ! ArgDelLam Field is required; Grid is not
        ! Store ArgDelLam info at first available position of the Grid array
        ! Mark field for DelLam

        willDelLam = .true.
        if (lugf <= nSpecFields) then
          surfGrid(lugf) = .true.
          prevGrid(lugf) = pugv
          nullify(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgDelLam
          fieldForDelLam(lugf) = .true.
          prevVertDelLamSource(lugf) = prevGrid(lugf)
        else
          write(p_nfprt, "(a, ' too many grid fields')") h
          write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          stop h
        end if

        if (dumpLocal) then
          write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
            h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
            h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
        end if
      end if

    end if

  end subroutine Deposit1D_PK

  subroutine Deposit1D (ArgSpec, ArgGrid, ArgDelLam)
    !# Deposit 1D
    !# ---
    !# @info
    !# **Brief:** Gives the input field to be transformed and the output field
    !# that will receive the transformed field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    real(kind = p_r8), target, intent(in) :: ArgSpec(:)
    real(kind = p_r8), target, optional, intent(in) :: ArgGrid(:, :)
    real(kind = p_r8), target, optional, intent(in) :: ArgDelLam(:, :)
    character(len = *), parameter :: h = "**(Deposit1D)**"
    logical :: getGrid, getDelLam
    integer :: lusf
    integer :: pusv
    integer :: lusv

    integer :: lugf
    integer :: pugv
    integer :: lugv

    integer :: ludlf
    integer :: pudlv
    integer :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid = present(ArgGrid)
    getDelLam = present(ArgDelLam)

    ! critical section to update data structure counters

    if (.not.havesurf) then

      !
      !  this processor won't transform surface fields
      !  ( it only needs to know which are the correponding grid fields )
      !
      if (getGrid) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        surfGrid(lugf) = .true.
        prevGrid(lugf) = -1
        nullify(Grid3d(lugf)%p)
        Grid2d(lugf)%p => ArgGrid
        fieldForDelLam(lugf) = .false.
        prevVertDelLamSource(lugf) = -1
        if (getDelLam) then
          ludlf = usedDelLamFields + 1
          usedDelLamFields = ludlf
          surfGrid(ludlf) = .true.
          prevGrid(ludlf) = -1
          nullify(Grid3d(ludlf)%p)
          Grid2d(ludlf)%p => ArgDelLam
          fieldForDelLam(lugf) = .false.
          prevVertDelLamSource(lugf) = -1
        end if
      else if (getDelLam) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        surfGrid(lugf) = .true.
        prevGrid(lugf) = -1
        nullify(Grid3d(lugf)%p)
        Grid2d(lugf)%p => ArgDelLam
        fieldForDelLam(lugf) = .false.
        prevVertDelLamSource(lugf) = - 1
      else
        write(p_nfprt, "(a, ' no Grid or DelLam field required')") h
        stop h
      end if

    else

      lusf = usedSpecFields + 1
      usedSpecFields = lusf
      pusv = lastUsedSpecVert
      lusv = pusv + 1
      lastUsedSpecVert = lusv
      if (getGrid) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        pugv = lastUsedGridVert
        lugv = pugv + 1
        lastUsedGridVert = lugv
        if (getDelLam) then
          ludlf = usedDelLamFields + 1
          usedDelLamFields = ludlf
          pudlv = lastUsedDelLamVert
          ludlv = pudlv + 1
          lastUsedDelLamVert = ludlv
        end if
      else if (getDelLam) then
        lugf = usedGridFields + 1
        usedGridFields = lugf
        pugv = lastUsedGridVert
        lugv = pugv + 1
        lastUsedGridVert = lugv
      else
        write(p_nfprt, "(a, ' no Grid or DelLam field required')") h
        stop h
      end if

      ! deposit spectral field

      if (lusf <= nSpecFields) then
        surfSpec(lusf) = .true.
        prevSpec(lusf) = pusv
        nullify(Spec2d(lusf)%p)
        Spec1d(lusf)%p => ArgSpec
      else
        write(p_nfprt, "(a, ' too many spectral fields')") h
        write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
        stop h
      end if

      if (dumpLocal) then
        write(p_nfprt, "(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
          h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
      end if

      if (getGrid) then

        ! ArgGrid Field is required;
        ! Store ArgGrid Field info at first available position of the Grid array

        if (lugf <= nSpecFields) then
          surfGrid(lugf) = .true.
          prevGrid(lugf) = pugv
          nullify(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgGrid
          fieldForDelLam(lugf) = .false.
          prevVertDelLamSource(lugf) = -1
        else
          write(p_nfprt, "(a, ' too many grid fields')") h
          write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          stop h
        end if

        if (dumpLocal) then
          write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
            h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
            h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
        end if

        if (getDelLam) then

          ! ArgGrid and ArgDelLam are required;
          ! Store ArgDelLam info at first available position of the Grid array
          ! Mark field for DelLam

          ! TODO check pudlv not initialized
          willDelLam = .true.
          if (ludlf <= nGridFields) then
            surfGrid(ludlf) = .true.
            prevGrid(ludlf) = pudlv
            nullify(Grid3d(ludlf)%p)
            Grid2d(ludlf)%p => ArgDelLam
            fieldForDelLam(ludlf) = .true.
            prevVertDelLamSource(ludlf) = prevGrid(lugf)
          else
            write(p_nfprt, "(a, ' too many grid fields including DelLam')") h
            write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
            stop h
          end if

          if (dumpLocal) then
            write(p_nfprt, "(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
              h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
            write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
              h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
          end if
        end if

      else if (getDelLam) then

        ! ArgDelLam Field is required; Grid is not
        ! Store ArgDelLam info at first available position of the Grid array
        ! Mark field for DelLam

        willDelLam = .true.
        if (lugf <= nSpecFields) then
          surfGrid(lugf) = .true.
          prevGrid(lugf) = pugv
          nullify(Grid3d(lugf)%p)
          Grid2d(lugf)%p => ArgDelLam
          fieldForDelLam(lugf) = .true.
          prevVertDelLamSource(lugf) = prevGrid(lugf)
        else
          write(p_nfprt, "(a, ' too many grid fields')") h
          write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
          stop h
        end if

        if (dumpLocal) then
          write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
            h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
          write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
            h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
        end if
      end if

    end if

  end subroutine Deposit1D


  subroutine Deposit2D (ArgSpec, ArgGrid, ArgDelLam)
    !# Deposit 2D
    !# ---
    !# @info
    !# **Brief:** Gives the input field to be transformed and the output field
    !# that will receive the transformed field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    real(kind = p_r8), target, intent(in) :: ArgSpec(:, :)
    real(kind = p_r8), target, optional, intent(in) :: ArgGrid(:, :, :)
    real(kind = p_r8), target, optional, intent(in) :: ArgDelLam(:, :, :)
    character(len = *), parameter :: h = "**(Deposit2D)**"
    logical :: getGrid, getDelLam
    integer :: lusf
    integer :: pusv
    integer :: lusv

    integer :: lugf
    integer :: pugv
    integer :: lugv

    integer :: ludlf
    integer :: pudlv
    integer :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid = present(ArgGrid)
    getDelLam = present(ArgDelLam)

    ! critical section to update data structure counters

    lusf = usedSpecFields + 1
    usedSpecFields = lusf
    pusv = lastUsedSpecVert
    lusv = pusv + kMaxloc
    lastUsedSpecVert = lusv
    if (getGrid) then
      lugf = usedGridFields + 1
      usedGridFields = lugf
      pugv = lastUsedGridVert
      lugv = pugv + kMaxloc
      lastUsedGridVert = lugv
      if (getDelLam) then
        ludlf = usedDelLamFields + 1
        usedDelLamFields = ludlf
        pudlv = lastUsedDelLamVert
        ludlv = pudlv + kMaxloc
        lastUsedDelLamVert = ludlv
      end if
    else if (getDelLam) then
      lugf = usedGridFields + 1
      usedGridFields = lugf
      pugv = lastUsedGridVert
      lugv = pugv + kMaxloc
      lastUsedGridVert = lugv
    else
      write(p_nfprt, "(a, ' no Grid or DelLam field required')") h
      stop h
    end if

    ! deposit spectral field

    if (lusf <= nSpecFields) then
      surfSpec(lusf) = .false.
      prevSpec(lusf) = pusv
      nullify(Spec1d(lusf)%p)
      Spec2d(lusf)%p => ArgSpec
    else
      write(p_nfprt, "(a, ' too many spectral fields')") h
      write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
      stop h
    end if

    if (dumpLocal) then
      write(p_nfprt, "(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
        h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
    end if

    if (getGrid) then

      ! ArgGrid Field is required;
      ! Store ArgGrid Field info at first available position of the Grid array

      if (lugf <= nfull_s + nsurf_s) then
        surfGrid(lugf) = .false.
        prevGrid(lugf) = pugv
        nullify(Grid2d(lugf)%p)
        Grid3d(lugf)%p => ArgGrid
        fieldForDelLam(lugf) = .false.
        prevVertDelLamSource(lugf) = -1
      else
        write(p_nfprt, "(a, ' too many grid fields')") h
        write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
        stop h
      end if

      if (dumpLocal) then
        write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
          h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
        write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
          h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
      end if

      if (getDelLam) then

        ! ArgGrid and ArgDelLam are required;
        ! Store ArgDelLam info at first available position of the Grid array
        ! Mark field for DelLam

        ! TODO check pudlv not initialized
        willDelLam = .true.
        if (ludlf <= nGridFields) then
          surfGrid(ludlf) = .false.
          prevGrid(ludlf) = pudlv
          nullify(Grid2d(ludlf)%p)
          Grid3d(ludlf)%p => ArgDelLam
          fieldForDelLam(ludlf) = .true.
          prevVertDelLamSource(ludlf) = prevGrid(lugf)
        else
          write(p_nfprt, "(a, ' too many grid fields including DelLam')") h
          write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
          stop h
        end if

        if (dumpLocal) then
          write(p_nfprt, "(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
            h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
          write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
            h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
        end if
      end if

    else if (getDelLam) then

      ! ArgDelLam Field is required; Grid is not
      ! Store ArgDelLam info at first available position of the Grid array
      ! Mark field for DelLam

      willDelLam = .true.
      if (lugf <= nfull_s + nsurf_s) then
        surfGrid(lugf) = .false.
        prevGrid(lugf) = pugv
        nullify(Grid2d(lugf)%p)
        Grid3d(lugf)%p => ArgDelLam
        fieldForDelLam(lugf) = .true.
        prevVertDelLamSource(lugf) = prevGrid(lugf)
      else
        write(p_nfprt, "(a, ' too many grid fields')") h
        write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
        stop h
      end if

      if (dumpLocal) then
        write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
          h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
        write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
          h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
      end if
    end if
  end subroutine Deposit2D


  subroutine Deposit2D_PK (ArgSpec, ArgGrid, ArgDelLam)
    !# Deposit 2D_PK
    !# ---
    !# @info
    !# **Brief:** Gives the grid field to be transformed and the spectral field
    !# that will receive the transformed field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    real(kind = p_r8), target, intent(in) :: ArgSpec(:, :)
    real(kind = p_r8), target, optional, intent(in) :: ArgGrid(:, :, :)
    real(kind = p_r8), target, optional, intent(in) :: ArgDelLam(:, :, :)
    character(len = *), parameter :: h = "**(Deposit2D)**"
    logical :: getGrid, getDelLam
    integer :: lusf
    integer :: pusv
    integer :: lusv

    integer :: lugf
    integer :: pugv
    integer :: lugv

    integer :: ludlf
    integer :: pudlv
    integer :: ludlv

    ! Grid field required?
    ! DelLam field required?

    getGrid = present(ArgGrid)
    getDelLam = present(ArgDelLam)

    ! critical section to update data structure counters

    lusf = usedSpecFields + 1
    usedSpecFields = lusf
    pusv = lastUsedSpecVert
    lusv = pusv + kMaxloc
    lastUsedSpecVert = lusv
    if (getGrid) then
      lugf = usedGridFields + 1
      usedGridFields = lugf
      pugv = lastUsedGridVert
      lugv = pugv + kMaxloc
      lastUsedGridVert = lugv
      if (getDelLam) then
        ludlf = usedDelLamFields + 1
        usedDelLamFields = ludlf
        pudlv = lastUsedDelLamVert
        ludlv = pudlv + kMaxloc
        lastUsedDelLamVert = ludlv
      end if
    else if (getDelLam) then
      lugf = usedGridFields + 1
      usedGridFields = lugf
      pugv = lastUsedGridVert
      lugv = pugv + kMaxloc
      lastUsedGridVert = lugv
    else
      write(p_nfprt, "(a, ' no Grid or DelLam field required')") h
      stop h
    end if

    ! deposit spectral field

    if (lusf <= nSpecFields) then
      surfSpec(lusf) = .false.
      prevSpec(lusf) = pusv
      nullify(Spec1d(lusf)%p)
      Spec2d(lusf)%p => ArgSpec
    else
      write(p_nfprt, "(a, ' too many spectral fields')") h
      write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lusf, nSpecFields
      stop h
    end if

    if (dumpLocal) then
      write(p_nfprt, "(a,' usedSpecFields, lastUsedSpecVert, surfSpec, prevSpec=', 2i5,l2,i5)") &
        h, lusf, lusv, surfSpec(lusf), prevSpec(lusf)
    end if

    if (getGrid) then

      ! ArgGrid Field is required;
      ! Store ArgGrid Field info at first available position of the Grid array

      if (lugf <= nfull_g + nsurf_g) then
        surfGrid(lugf) = .false.
        prevGrid(lugf) = pugv
        nullify(Grid2d(lugf)%p)
        Grid3d(lugf)%p => ArgGrid
        fieldForDelLam(lugf) = .false.
        prevVertDelLamSource(lugf) = -1
      else
        write(p_nfprt, "(a, ' too many grid fields')") h
        write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
        stop h
      end if

      if (dumpLocal) then
        write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
          h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
        write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
          h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
      end if

      if (getDelLam) then

        ! ArgGrid and ArgDelLam are required;
        ! Store ArgDelLam info at first available position of the Grid array
        ! Mark field for DelLam

        willDelLam = .true.
        if (ludlf <= nGridFields) then
          surfGrid(ludlf) = .false.
          prevGrid(ludlf) = pudlv
          nullify(Grid2d(ludlf)%p)
          Grid3d(ludlf)%p => ArgDelLam
          fieldForDelLam(ludlf) = .true.
          prevVertDelLamSource(ludlf) = prevGrid(lugf)
        else
          write(p_nfprt, "(a, ' too many grid fields including DelLam')") h
          write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, ludlf, nGridFields
          stop h
        end if

        if (dumpLocal) then
          write(p_nfprt, "(a,' usedDelLamFields, lastUsedDelLamVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
            h, ludlf, ludlv, surfGrid(ludlf), prevGrid(ludlf)
          write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
            h, fieldForDelLam(ludlf), prevVertDelLamSource(ludlf)
        end if
      end if

    else if (getDelLam) then

      ! ArgDelLam Field is required; Grid is not
      ! Store ArgDelLam info at first available position of the Grid array
      ! Mark field for DelLam

      willDelLam = .true.
      if (lugf <= nfull_s + nsurf_s) then
        surfGrid(lugf) = .false.
        prevGrid(lugf) = pugv
        nullify(Grid2d(lugf)%p)
        Grid3d(lugf)%p => ArgDelLam
        fieldForDelLam(lugf) = .true.
        prevVertDelLamSource(lugf) = prevGrid(lugf)
      else
        write(p_nfprt, "(a, ' too many grid fields')") h
        write(p_nfprt, "(a, ' used=',i5,'; declared=',i5)") h, lugf, nSpecFields
        stop h
      end if

      if (dumpLocal) then
        write(p_nfprt, "(a,' usedGridFields, lastUsedGridVert, surfGrid, prevGrid=', 2i5,l2,i5)") &
          h, lugf, lugv, surfGrid(lugf), prevGrid(lugf)
        write(p_nfprt, "(a,' fieldForDelLam, prevVertDelLamSource=', l2,i5)") &
          h, fieldForDelLam(lugf), prevVertDelLamSource(lugf)
      end if
    end if
  end subroutine Deposit2D_PK


  subroutine DepDLG1D (ArgSpec, ArgGrid, ArgDelLam)
    !# Deposit DLG1D
    !# ---
    !# @info
    !# **Brief:** Use to compute both the transform and the zonal derivative of
    !# input spectral field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin 
    real(kind = p_r8), target, intent(in) :: ArgSpec(:)
    real(kind = p_r8), target, intent(in) :: ArgGrid(:, :)
    real(kind = p_r8), target, intent(in) :: ArgDelLam(:, :)
    call Deposit1D (ArgSpec, ArgGrid, ArgDelLam)
  end subroutine DepDLG1D


  subroutine DepDLG2D (ArgSpec, ArgGrid, ArgDelLam)
    !# Deposit DLG2D
    !# ---
    !# @info
    !# **Brief:** Use to compute both the transform and the zonal derivative of
    !# input spectral field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin 
    real(kind = p_r8), target, intent(in) :: ArgSpec(:, :)
    real(kind = p_r8), target, intent(in) :: ArgGrid(:, :, :)
    real(kind = p_r8), target, intent(in) :: ArgDelLam(:, :, :)
    call Deposit2D (ArgSpec, ArgGrid, ArgDelLam)
  end subroutine DepDLG2D


  subroutine DepDL1D (ArgSpec, ArgDelLam)
    !# Deposit DL1D
    !# ---
    !# @info
    !# **Brief:** Use only to compute the zonal derivative of input spectral field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    real(kind = p_r8), target, intent(in) :: ArgSpec(:)
    real(kind = p_r8), target, intent(in) :: ArgDelLam(:, :)
    call Deposit1D (ArgSpec, ArgDelLam = ArgDelLam)
  end subroutine DepDL1D


  subroutine DepDL2D (ArgSpec, ArgDelLam)
    !# Deposit DL2D
    !# ---
    !# @info
    !# **Brief:** Use only to compute the zonal derivative of input spectral field. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    real(kind = p_r8), target, intent(in) :: ArgSpec(:, :)
    real(kind = p_r8), target, intent(in) :: ArgDelLam(:, :, :)
    call Deposit2D (ArgSpec, ArgDelLam = ArgDelLam)
  end subroutine DepDL2D


  subroutine DoSpecToGrid(tamBlock)
    !# Do Spectral to Grid
    !# ---
    !# @info
    !# **Brief:** Perform all required transforms and zonal derivatives simultaneosly. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: tamBlock
    character(len = *), parameter :: h = "**(DoSpecToGrid)**"
    integer :: mnFirst
    integer :: mnLast
    integer :: mnExtFirst
    integer :: mnExtLast
    integer :: mFirst
    integer :: mLast
    integer :: jbFirst
    integer :: jbLast
    integer :: iBlockFFT
    integer :: first
    integer :: last
    integer :: FFTfirst
    integer :: FFTLast
    integer :: FFTSize
    integer :: firstBlock
    integer :: lastBlock
    integer :: sizeBlock
    integer :: i
    integer :: j
    integer :: iFirst
    integer :: iLast
    TYPE(MultiFFT), pointer :: p


    ! start OMP parallelism

    call ThreadDecomp(1, mymnMax, mnFirst, mnLast, "DoSpecToGrid")
    call ThreadDecomp(1, mymnExtMax, mnExtFirst, mnExtLast, "DoSpecToGrid")
    call ThreadDecomp(1, mymMax, mFirst, mLast, "DoSpecToGrid")
    call ThreadDecomp(1, jbMax, jbFirst, jbLast, "DoSpecToGrid")
    call ThreadDecomp(1, dip1, iFirst, iLast, "DoSpecToGrid")

    ! nullify Four
    if (nfull_g.gt.0.or.havesurf) then

      do i = iFirst, iLast
        do j = 1, dvdlj
          Four(j, i) = 0.0_p_r8
        end do
      end do

      ! ingest spectral fields

      call DepositSpec(mnFirst, mnLast, mnExtFirst, mnExtLast, mFirst, mLast)

      ! Fourier from Spectral
      !$OMP BARRIER

      call SpecToFour()

      ! DelLam where required

      if (willDelLam) then
        call ThreadDecomp(1, mMax, mFirst, mLast, "DoSpecToGrid")
        call DelLam(mFirst, mLast)
      end if

      ! FFT Fourier to Grid
      !$OMP BARRIER

      do iBlockFFT = 1, nBlockFFT
        p => BlockFFT(iBlockFFT)
        first = (p%firstLat - 1) * nVertGrid + 1
        last = (p%lastLat) * nVertGrid
        call ThreadDecomp(first, last, FFTFirst, FFTLast, "DoSpecToGrid")
        FFTSize = FFTLast - FFTFirst + 1
        if (FFTSize.le.0) cycle
        if(tamBlock == 0) then
          call InvFFTTrans (Four(FFTfirst, 1), dvdlj, dip1, p%nLong, FFTSize, &
            p%Trigs, p%nTrigs, p%Factors, p%nFactors)
        else
          firstBlock = FFTFirst
          do
            lastBlock = min(firstBlock + tamBlock - 1, FFTLast)
            sizeBlock = lastBlock - firstBlock + 1
            call InvFFTTrans (Four(firstBlock, 1), dvdlj, dip1, p%nLong, sizeBlock, &
              p%Trigs, p%nTrigs, p%Factors, p%nFactors)
            firstBlock = firstBlock + tamBlock
            if (firstBlock > FFTLast) then
              EXIT
            end if
          end do
        end if

      end do

    end if
    !$OMP BARRIER

    ! Withdraw Grid Fields

    call WithdrawGrid

  end subroutine DoSpecToGrid


  subroutine DepositSpec(mnFirst, mnLast, mnExtFirst, mnExtLast, mFirst, mLast)
    !# Deposit Spectral
    !# ---
    !# @info
    !# **Brief:** Ingest spectral fields.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: mnFirst
    integer, intent(in) :: mnLast
    integer, intent(in) :: mnExtFirst
    integer, intent(in) :: mnExtLast
    integer, intent(in) :: mFirst
    integer, intent(in) :: mLast
    integer :: is
    integer :: imn
    integer :: iv
    integer :: lastv
    real(kind = p_r8), pointer :: s1(:)
    real(kind = p_r8), pointer :: s2(:, :)

    do is = 1, nSpecFields

      if (surfSpec(is)) then
        s1 => Spec1d(is)%p
        lastv = prevSpec(is) + 1
        if (size(s1, 1) == 2 * myMNExtMax) then
          !CDIR NODEP
          do imn = mnExtFirst, mnExtLast
            Spec(lmnExtMap(imn), lastv) = s1(2 * imn - 1)
            Spec(lmnExtMap(imn), lastv + nVertSpec) = s1(2 * imn)
          end do
        else
          !CDIR NODEP
          do imn = mnFirst, mnLast
            Spec(lmnMap(imn), lastv) = s1(2 * imn - 1)
            Spec(lmnMap(imn), lastv + nVertSpec) = s1(2 * imn)
          end do
          !CDIR NODEP
          do imn = mFirst, mLast
            Spec(lmnZero(imn), lastv) = 0.0_p_r8
            Spec(lmnZero(imn), lastv + nVertSpec) = 0.0_p_r8
          end do
        end if

      else
        s2 => Spec2d(is)%p
        lastv = prevSpec(is)
        if (size(s2, 1) == 2 * myMNExtMax) then
          do iv = 1, kMaxloc
            !CDIR NODEP
            do imn = mnExtFirst, mnExtlast
              Spec(lmnExtMap(imn), iv + lastv) = s2(2 * imn - 1, iv)
              Spec(lmnExtMap(imn), iv + lastv + nVertSpec) = s2(2 * imn, iv)
            end do
          end do
        else
          do iv = 1, kMaxloc
            !CDIR NODEP
            do imn = mnFirst, mnLast
              Spec(lmnMap(imn), iv + lastv) = s2(2 * imn - 1, iv)
              Spec(lmnMap(imn), iv + lastv + nVertSpec) = s2(2 * imn, iv)
            end do
          end do
          do iv = 1, kMaxloc
            !CDIR NODEP
            do imn = mFirst, mLast
              Spec(lmnZero(imn), iv + lastv) = 0.0_p_r8
              Spec(lmnZero(imn), iv + lastv + nVertSpec) = 0.0_p_r8
            end do
          end do
        end if
      end if
    end do

    !$OMP SINGLE
    do imn = myMNExtMax + 1, dlmn
      Spec(imn, :) = 0.0_p_r8
    end do
    do iv = 2 * nVertSpec + 1, dv
      Spec(:, iv) = 0.0_p_r8
    end do
    !$OMP END SINGLE
  end subroutine DepositSpec


  subroutine SpecToFour()
    !# Spectral to Fourier
    !# ---
    !# @info
    !# **Brief:** Spectral to Fourier.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer :: lm
    integer :: m
    integer :: j
    integer :: v
    integer :: js
    integer :: jn
    integer :: kn
    integer :: kdim
    integer :: ldim
    integer :: ks
    integer :: k
    integer :: inin
    integer :: inis
    integer :: inin1
    integer :: inis1
    integer :: mm
    integer :: comm
    integer :: ierr
    integer :: index
    integer :: jfirst
    integer :: jlast
    real(kind = p_r8) :: FoEv(dv, djh)
    real(kind = p_r8) :: FoOd(dv, djh)

    call ThreadDecomp(1, myjmax_f, jfirst, jlast, "SpecToFour")
    kdim = nVertSpec * myjMax_f * 2 * MMaxlocal
    ldim = nVertSpec * jMaxlocal_f * 2 * myMMax
    !$OMP SINGLE
    mGlob = 0
    if (dimsendbuf.lt.ldim * maxnodes_four) then
      dimsendbuf = ldim * maxnodes_four
      deallocate (bufsend)
      allocate (bufsend(dimsendbuf))
    endif
    if (dimrecbuf.lt.kdim * maxnodes_four) then
      dimrecbuf = kdim * maxnodes_four
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    bufsend = 0.0_p_r8
    bufrec = 0.0_p_r8
    !$OMP END SINGLE

    do
      !$OMP CRITICAL(lmcrit)
      mGlob = mGlob + 1
      lm = mGlob
      !$OMP END CRITICAL(lmcrit)
      if (lm > myMMax) EXIT

      m = lm2m(lm)

      ! Spectral to Fourier Even/Odd

      FoEv = 0.0_p_r8
      call mmt (LS2F(jMinPerM(m), FirstNEven(lm)), djh, &
        Spec(FirstNEven(lm), 1), dlmn, &
        FoEv(1, jMinPerM(m)), dv, djh, &
        jMaxHalf - jMinPerM(m) + 1, 2 * nVertSpec, nEven(lm))
      FoOd = 0.0_p_r8
      call mmt (LS2F(jMinPerM(m), FirstNOdd(lm)), djh, &
        Spec(FirstNOdd(lm), 1), dlmn, &
        FoOd(1, jMinPerM(m)), dv, djh, &
        jMaxHalf - jMinPerM(m) + 1, 2 * nVertSpec, nOdd(lm))

      ! Fourier Even/Odd to Full Fourier

      do jn = jMinPerM(m), jMaxHalf
        js = jmax - jn + 1
        kn = nodeHasJ_f(jn)
        ks = nodeHasJ_f(js)
        inin = (jn - firstlatinproc_f(kn) + nlatsinproc_f(kn) * 2 * (lm - 1)) * nVertSpec
        inin1 = inin + nlatsinproc_f(kn) * nVertSpec
        inis = (js - firstlatinproc_f(ks) + nlatsinproc_f(ks) * 2 * (lm - 1)) * nVertSpec
        inis1 = inis + nlatsinproc_f(ks) * nVertSpec
        !CDIR NODEP
        do v = 1, nVertSpec
          bufsend(inin + v + kn * ldim) = FoEv(v, jn) + FoOd(v, jn)
          bufsend(inin1 + v + kn * ldim) = FoEv(v + nVertSpec, jn) + FoOd(v + nVertSpec, jn)
          bufsend(inis + v + ks * ldim) = FoEv(v, jn) - FoOd(v, jn)
          bufsend(inis1 + v + ks * ldim) = FoEv(v + nVertSpec, jn) - FoOd(v + nVertSpec, jn)
        end do
      end do
    end do
    !$OMP BARRIER
    !$OMP SINGLE
    comm = COMM_FOUR
    requestr(myid_four) = MPI_REQUEST_NULL
    requests(myid_four) = MPI_REQUEST_NULL
    do k = 0, MaxNodes_four - 1
      if(k.ne.myid_four) then
        call MPI_IRECV(bufrec(1 + k * kdim), 2 * myjmax_f * nvertspec * MsPerProc(k), &
          MPI_DOUBLE_PRECISION, k, 97, comm, requestr(k), ierr)
      endif
    enddo
    do k = 0, MaxNodes_four - 1
      if(k.ne.myid_four) then
        call MPI_ISEND(bufsend(1 + k * ldim), 2 * nlatsinproc_f(k) * nvertspec * mymmax, &
          MPI_DOUBLE_PRECISION, k, 97, comm, requests(k), ierr)
      endif
    enddo
    !$OMP END SINGLE
    mm = 0
    do m = 1, Mmax
      if (myid_four.eq.nodeHasM(m, mygroup_four)) then
        mm = mm + 1
        inin = (myjmax_f * 2 * (mm - 1)) * nVertSpec
        inin1 = inin + myjmax_f * nVertSpec
        do j = jfirst, jlast
          do v = 1, nVertSpec
            Four((j - 1) * nVertGrid + v, 2 * m - 1) = &
              bufsend((j - 1) * nVertSpec + v + inin + myid_four * ldim)
            Four((j - 1) * nVertGrid + v, 2 * m) = &
              bufsend((j - 1) * nVertSpec + v + inin1 + myid_four * ldim)
          end do
        end do
      endif
    end do
    do k = 1, MaxNodes_four - 1
      !$OMP BARRIER
      !$OMP SINGLE
      call MPI_WAITANY(MaxNodes_four, requestr(0), index, status, ierr)
      !$OMP END SINGLE
      ks = status(MPI_SOURCE)
      mm = 0
      do m = 1, Mmax
        if (ks.eq.nodeHasM(m, mygroup_four)) then
          mm = mm + 1
          inin = (myjmax_f * 2 * (mm - 1)) * nVertSpec
          inin1 = inin + myjmax_f * nVertSpec
          do j = jfirst, jlast
            do v = 1, nVertSpec
              Four((j - 1) * nVertGrid + v, 2 * m - 1) = &
                bufrec((j - 1) * nVertSpec + v + inin + ks * kdim)
              Four((j - 1) * nVertGrid + v, 2 * m) = &
                bufrec((j - 1) * nVertSpec + v + inin1 + ks * kdim)
            end do
          end do
        endif
      end do
    enddo
    !$OMP SINGLE
    call MPI_WAITALL(maxnodes_four, requests(0), stat, ierr)
    !$OMP END SINGLE
  end subroutine SpecToFour


  subroutine DelLam(mFirst, mLast)
    !# DelLam...?
    !# ---
    !# @info
    !# **Brief:** DelLam...?  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: mFirst
    integer, intent(in) :: mLast

    integer :: m
    integer :: j
    integer :: v
    real(kind = p_r8) :: auxRe
    integer :: ig
    integer :: vBaseFrom
    integer :: vBaseTo
    integer :: jAux
    integer :: vMax

    do m = mFirst, mLast

      do ig = 1, nGridFields
        if (fieldForDelLam(ig)) then
          vBaseFrom = prevVertDelLamSource(ig)
          vBaseTo = prevGrid(ig)

          if (surfGrid(ig)) then
            if (havesurf) then
              vMax = 1
            else
              vMax = 0
            end if
          else
            vMax = kMaxloc
          end if
          !CDIR NODEP
          do v = 1, vMax
            do j = 1, myjMax_f
              jAux = (j - 1) * nVertGrid + v
              auxRe = consIm(m) * Four(jAux + vBaseFrom, 2 * m - 1)
              Four(jAux + vBaseTo, 2 * m - 1) = consRe(m) * Four(jAux + vBaseFrom, 2 * m)
              Four(jAux + vBaseTo, 2 * m) = auxRe
            end do
          end do
        end if
      end do
    end do
  end subroutine DelLam


  subroutine WithdrawGrid
    !# Withdraw Grid
    !# ---
    !# @info
    !# **Brief:** Withdraw Grid Fields.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer :: ig
    integer :: n
    integer :: m
    integer :: l
    integer :: k
    integer :: kl
    integer :: ns
    integer :: j
    integer :: i
    integer :: ipar
    integer :: ib
    integer :: v0
    integer :: proc
    integer :: comm
    integer :: index
    integer :: ierr
    integer :: ibr(nrecs_f + 1)
    integer :: ibn(nrecs_f + 1)
    integer :: ibs(nrecs_g + 1)
    real(kind = p_r8), pointer :: g2(:, :)
    real(kind = p_r8), pointer :: g3(:, :, :)
    character(len = *), parameter :: h = "**(WithdrawGrid)**"

    comm = mpiCommGroup
    m = 0
    ibr(1) = 1
    do n = 1, nrecs_f
      ibn(n) = 0
      do ipar = m + 1, messproc_f(2, n)
        ibn(n) = ibn(n) + 1 + messages_f(2, ipar) - messages_f(1, ipar)
      end do
      ibr(n + 1) = ibr(n) + ibn(n) * nvertgrid
      m = messproc_f(2, n)
    end do
    m = 0
    ibs(1) = 1
    do n = 1, nrecs_g
      ib = 0
      do ipar = m + 1, messproc_g(2, n)
        ib = ib + 1 + messages_g(2, ipar) - messages_g(1, ipar)
      end do
      k = klast_four(messproc_g(1, n)) - kfirst_four(messproc_g(1, n)) + 1
      ns = 0
      if (kfirst_four(messproc_g(1, n)).eq.1) ns = nsurf_g
      ibs(n + 1) = ibs(n) + ib * (k * nfull_g + ns)
      m = messproc_g(2, n)
    end do
    !$OMP SINGLE
    if (dimsendbuf.lt.ibs(nrecs_g + 1) - 1) then
      dimsendbuf = ibs(nrecs_g + 1)
      deallocate (bufsend)
      allocate (bufsend(dimsendbuf))
    endif
    if (dimrecbuf.lt.ibr(nrecs_f + 1) - 1) then
      dimrecbuf = ibr(nrecs_f + 1)
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    do n = 1, nrecs_g
      proc = messproc_g(1, n)
      call MPI_IRECV(bufsend(ibs(n)), ibs(n + 1) - ibs(n), MPI_DOUBLE_PRECISION, &
        proc, 78, comm, requestr(n), ierr)
    end do
    nglob = 1
    iglob = 0
    ibglob = 0
    !$OMP END SINGLE
    do
      !$OMP CRITICAL(nigcrit)
      iglob = iglob + 1
      if (iglob.gt.nGridFields) then
        iglob = 1
        nglob = nglob + 1
      endif
      n = nglob
      if (nglob.le.nrecs_f) then
        ig = iglob
        ib = ibglob
        kl = kmaxloc
        if (surfGrid(ig)) then
          if (havesurf) then
            kl = 1
          else
            kl = 0
          endif
        endif
        ibglob = ibglob + kl * ibn(n)
      endif
      !$OMP END CRITICAL(nigcrit)
      if (n.gt.nrecs_f) EXIT
      m = messproc_f(2, n - 1)
      if (surfGrid(ig)) then
        if (havesurf) then
          do ipar = m + 1, messproc_f(2, n)
            j = messages_f(3, ipar)
            v0 = previousJ(j) * nVertGrid + prevGrid(ig)
            do i = messages_f(1, ipar), messages_f(2, ipar)
              ib = ib + 1
              bufrec(ib) = Four(v0 + 1, i)
            end do
          end do
        end if
      else
        do ipar = m + 1, messproc_f(2, n)
          j = messages_f(3, ipar)
          v0 = previousJ(j) * nVertGrid + prevGrid(ig)
          do k = 1, kmaxloc
            do i = messages_f(1, ipar), messages_f(2, ipar)
              ib = ib + 1
              bufrec(ib) = Four(k + v0, i)
            end do
          end do
        end do
      end if
    end do
    !$OMP BARRIER
    !$OMP SINGLE
    do n = 1, nrecs_f
      proc = messproc_f(1, n)
      call MPI_ISEND(bufrec(ibr(n)), ibr(n + 1) - ibr(n), MPI_DOUBLE_PRECISION, proc, 78, &
        comm, requests(n), ierr)
    end do
    !$OMP END SINGLE
    iglob = 1
    jglob = max(myfirstlat, myfirstlat_f)
    !$OMP BARRIER
    !
    ! local values
    ! ------------
    if (max(myfirstlat, myfirstlat_f).le.min(mylastlat, mylastlat_f)) then
      do
        !$OMP CRITICAL(loccrit)
        ig = iglob
        do
          if (iglob.gt.nGridFields) EXIT
          if (.not.surfGrid(iglob).or.myfirstlev.eq.1) EXIT
          iglob = iglob + 1
        enddo
        ig = iglob
        if (iglob.le.nGridFields) then
          j = jglob
          if (jglob.ge.min(mylastlat, mylastlat_f)) then
            iglob = iglob + 1
            jglob = max(myfirstlat, myfirstlat_f)
          else
            jglob = jglob + 1
          endif
        endif
        !$OMP END CRITICAL(loccrit)
        if (ig.gt.nGridFields) EXIT
        if (surfGrid(ig)) then
          g2 => Grid2d(ig)%p
          v0 = previousJ(j) * nVertGrid + prevGrid(ig)
          do i = myfirstlon(j), mylastlon(j)
            g2(ibPerIJ(i, j), jbPerIJ(i, j)) = Four(v0 + 1, i)
          end do
        else
          g3 => Grid3d(ig)%p
          v0 = previousJ(j) * nVertGrid + prevGrid(ig)
          do k = myfirstlev, mylastlev
            v0 = v0 + 1
            do i = myfirstlon(j), mylastlon(j)
              g3(ibPerIJ(i, j), k, jbPerIJ(i, j)) = Four(v0, i)
            end do
          end do
        end if
      end do
    end if

    !
    !$OMP SINGLE
    kountg = 1
    !$OMP END SINGLE
    if (nrecs_g.gt.0) then
      do
        !$OMP SINGLE
        call MPI_WAITANY(nrecs_g, requestr(1), index, status, ierr)
        ksg = status(MPI_SOURCE)
        do l = 1, nrecs_g
          if (ksg.eq.messproc_g(1, l)) then
            nglob = l
            ibglob = ibs(nglob) - 1
            mglob = messproc_g(2, nglob - 1)
            EXIT
          endif
        enddo
        iglob = 1
        ipar2g = mglob + 1
        ipar3g = mglob + 1
        !$OMP END SINGLE
        n = nglob
        m = mglob
        kl = klast_four(ksg) - kfirst_four(ksg) + 1
        do
          !$OMP CRITICAL(reccrit)
          ig = iglob
          do
            if (iglob.gt.nGridFields) EXIT
            if (.not.surfGrid(iglob).or.kfirst_four(ksg).eq.1) EXIT
            iglob = iglob + 1
          enddo
          ig = iglob
          if (iglob.le.nGridFields) then
            ib = ibglob
            if (surfGrid(ig)) then
              ipar = ipar2g
              ibglob = ibglob + messages_g(2, ipar) - messages_g(1, ipar) + 1
              ipar2g = ipar2g + 1
              if (ipar2g.gt.messproc_g(2, n)) then
                ipar2g = m + 1
                iglob = iglob + 1
              endif
            else
              ipar = ipar3g
              ibglob = ibglob + kl * (messages_g(2, ipar) - messages_g(1, ipar) + 1)
              ipar3g = ipar3g + 1
              if (ipar3g.gt.messproc_g(2, n)) then
                ipar3g = m + 1
                iglob = iglob + 1
              endif
            endif
          endif
          !$OMP END CRITICAL(reccrit)
          if (ig.gt.nGridFields) EXIT
          if (surfGrid(ig)) then
            g2 => Grid2d(ig)%p
            j = messages_g(3, ipar)
            v0 = previousJ(j) * nVertGrid + prevGrid(ig)
            do i = messages_g(1, ipar), messages_g(2, ipar)
              ib = ib + 1
              g2(ibperij(i, j), jbperij(i, j)) = bufsend(ib)
            end do
          else
            g3 => Grid3d(ig)%p
            j = messages_g(3, ipar)
            v0 = previousJ(j) * nVertGrid + prevGrid(ig)
            do k = kfirst_four(ksg), klast_four(ksg)
              do i = messages_g(1, ipar), messages_g(2, ipar)
                ib = ib + 1
                g3(ibperij(i, j), k, jbperij(i, j)) = bufsend(ib)
              end do
            end do
          end if
        end do
        !$OMP SINGLE
        kountg = kountg + 1
        !$OMP END SINGLE
        if (kountg.GT.nrecs_g) EXIT
      end do
    end if

    !$OMP SINGLE
    call MPI_WAITALL(nrecs_f, requests(1), stat, ierr)
    !$OMP END SINGLE

  end subroutine WithdrawGrid


  subroutine CreateGridToSpec(nFull, nSurf)
    !# Creates Grid To Spectral
    !# ---
    !# @info
    !# **Brief:** Prepares internal data structure (arrays) to perform Fourier
    !# and then Legendre Transforms for a given number (input arguments) of
    !# surface and full fields. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: nFull
    integer, intent(in) :: nSurf
    character(len = *), parameter :: h = "**(CreateGridToSpec)**"
    integer :: nsu

    if (havesurf) then
      nsu = nSurf
    else
      nsu = 0
    endif
    nFull_g = nFull
    nSurf_g = nSurf

    nSpecFields = nFull + nsu
    usedSpecFields = 0
    lastUsedSpecVert = 0
    if (.not.allocated(surfSpec))allocate(surfSpec(nSpecFields))
    if (.not.allocated(prevSpec))allocate(prevSpec(nSpecFields))
    if (.not.allocated(Spec1d))allocate(Spec1d(nSpecFields))
    if (.not.allocated(Spec2d))allocate(Spec2d(nSpecFields))
    nVertSpec = nFull * kMaxloc + nsu

    nGridFields = nFull + nSurf
    usedGridFields = 0
    lastUsedGridVert = 0
    if (.not.allocated(surfGrid)) allocate(surfGrid(nGridFields))
    if (.not.allocated(prevGrid)) allocate(prevGrid(nGridFields))
    if (.not.allocated(Grid2d)) allocate(Grid2d(nGridFields))
    if (.not.allocated(Grid3d)) allocate(Grid3d(nGridFields))
    if (.not.allocated(fieldForDelLam)) allocate(fieldForDelLam(nGridFields))
    nVertGrid = nFull * kMaxloc + nsu

    willDelLam = .false.
    usedDelLamFields = nSpecFields
    lastUsedDelLamVert = nVertSpec
    if (.not.allocated(prevVertDelLamSource)) allocate(prevVertDelLamSource(nGridFields))

    dv = NoBankConflict(2 * nVertSpec)
    dvjh = NoBankConflict(nVertSpec * jMaxHalf)
    dvdlj = NoBankConflict(nVertGrid * myJMax_f)

    if (.not.allocated(Spec)) allocate (Spec(dlmn, dv))
    if (.not.allocated(Four)) allocate (Four(dvdlj, dip1))

    allocate(mnodes(0:maxnodes_four - 1))
    if (.not.allocated(requests)) then
      allocate(requests(0:maxnodes))
      allocate(requestr(0:maxnodes))
    end if
    if (.not.allocated(status)) then
      allocate(status(MPI_STATUS_SIZE))
    end if
    if (.not.allocated(stat)) then
      allocate(stat(MPI_STATUS_SIZE, maxnodes))
    end if

    if (dumpLocal) then
      write(p_nfprt, "(a,' Spec: n, used, lastVert, nVert=',4i5)") &
        h, nSpecFields, usedSpecFields, lastUsedSpecVert, nVertSpec
      write(p_nfprt, "(a,' Grid: n, used, lastVert, nVert=',4i5)") &
        h, nGridFields, usedGridFields, lastUsedGridVert, nVertGrid
    end if

  end subroutine CreateGridToSpec


  subroutine DoGridToSpec(tamBlock)
    !# Do Grid To Spectral
    !# ---
    !# @info
    !# **Brief:** 
    !# <ul type="disc">
    !#  <li>Performs the transform over all previously deposited grid fields,
    !# since last CreateGridToSpec, saving output data on the pointed positions. </li>
    !#  <li>The number of surface grid fields and full grid fields deposited should
    !# match the information given at CreateGridToSpec. </li>
    !#  <li>At the end of DoGridToSpec, all Spectral Fields informed by
    !# DepositGridToSpec are filled with required information. </li>
    !# </ul> </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: tamBlock

    integer :: mnFirst
    integer :: mnLast
    integer :: mnExtFirst
    integer :: mnExtLast
    integer :: jbFirst
    integer :: jbLast
    integer :: iBlockFFT
    integer :: first
    integer :: last
    integer :: FFTfirst
    integer :: FFTLast
    integer :: FFTSize
    integer :: firstBlock
    integer :: lastBlock
    integer :: sizeBlock
    TYPE(MultiFFT), pointer :: p
    character(len = *), parameter :: h = "**(DoGridToSpec)**"

    ! all fields were deposited?

    if (usedSpecFields /= nSpecFields) then
      write(p_nfprt, "(a, ' not all spectral fields were deposited')") h
      stop h
    else if (usedGridFields /= nGridFields) then
      write(p_nfprt, "(a, ' not all gauss fields were deposited')") h
      stop h
    end if

    ! start OMP parallelism

    call ThreadDecomp(1, mymnMax, mnFirst, mnLast, "DoGridToSpec")
    call ThreadDecomp(1, mymnExtMax, mnExtFirst, mnExtLast, "DoGridToSpec")
    call ThreadDecomp(1, jbMax, jbFirst, jbLast, "DoGridToSpec")


    ! deposit all grid fields

    call DepositGrid
    !$OMP BARRIER

    ! FFT Grid to Fourier

    do iBlockFFT = 1, nBlockFFT
      p => BlockFFT(iBlockFFT)
      first = (p%firstLat - 1) * nVertGrid + 1
      last = (p%lastLat) * nVertGrid
      call ThreadDecomp(first, last, FFTFirst, FFTLast, "DoGridToSpec")
      FFTSize = FFTLast - FFTFirst + 1
      if (FFTSize.le.0) cycle
      if(tamBlock == 0) then
        call DirFFTTrans (Four(FFTfirst, 1), dvdlj, dip1, p%nLong, FFTSize, &
          p%Trigs, p%nTrigs, p%Factors, p%nFactors)
      else
        firstBlock = FFTFirst
        do
          lastBlock = min(firstBlock + tamBlock - 1, FFTLast)
          sizeBlock = lastBlock - firstBlock + 1
          call DirFFTTrans (Four(firstBlock, 1), dvdlj, dip1, p%nLong, sizeBlock, &
            p%Trigs, p%nTrigs, p%Factors, p%nFactors)
          firstBlock = firstBlock + tamBlock
          if (firstBlock > FFTLast) then
            EXIT
          end if
        end do
      end if

      !CALL DirFFTTrans (Four(FFTfirst ,1), dvdlj, dip1, p%nLong, FFTSize, &
      !    p%Trigs, p%nTrigs, p%Factors, p%nFactors)

    end do

    !$OMP BARRIER

    ! Fourier to Spectral

    call FourToSpec()
    !$OMP BARRIER

    ! retrieve Spectral fields

    call WithdrawSpectral(mnFirst, mnLast, mnExtFirst, mnExtLast)
  end subroutine DoGridToSpec


  subroutine DepositGrid
    !# Deposit Grid
    !# ---
    !# @info
    !# **Brief:** Deposit all grid fields.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer :: ig
    integer :: n
    integer :: m
    integer :: l
    integer :: k
    integer :: kl
    integer :: ns
    integer :: j
    integer :: i
    integer :: ipar
    integer :: ib
    integer :: v0
    integer :: proc
    integer :: comm
    integer :: index
    integer :: ierr
    integer :: ibr(nrecs_f + 1)
    integer :: ibs(nrecs_g + 1)
    real(kind = p_r8), pointer :: g2(:, :)
    real(kind = p_r8), pointer :: g3(:, :, :)

    comm = mpiCommGroup
    m = 0
    ibr(1) = 1
    do n = 1, nrecs_f
      ib = 0
      do ipar = m + 1, messproc_f(2, n)
        ib = ib + 1 + messages_f(2, ipar) - messages_f(1, ipar)
      end do
      ibr(n + 1) = ibr(n) + ib * nvertgrid
      m = messproc_f(2, n)
    end do
    m = 0
    ibs(1) = 1
    do n = 1, nrecs_g
      ib = 0
      do ipar = m + 1, messproc_g(2, n)
        ib = ib + 1 + messages_g(2, ipar) - messages_g(1, ipar)
      end do
      k = klast_four(messproc_g(1, n)) - kfirst_four(messproc_g(1, n)) + 1
      ns = 0
      if (kfirst_four(messproc_g(1, n)).eq.1) ns = nsurf_g
      ibs(n + 1) = ibs(n) + ib * (k * nfull_g + ns)
      m = messproc_g(2, n)
    end do
    !$OMP SINGLE
    if (dimsendbuf.lt.ibs(nrecs_g + 1) - 1) then
      dimsendbuf = ibs(nrecs_g + 1)
      deallocate (bufsend)
      allocate (bufsend(dimsendbuf))
    endif
    if (dimrecbuf.lt.ibr(nrecs_f + 1) - 1) then
      dimrecbuf = ibr(nrecs_f + 1)
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    do n = 1, nrecs_f
      proc = messproc_f(1, n)
      call MPI_IRECV(bufrec(ibr(n)), ibr(n + 1) - ibr(n), MPI_DOUBLE_PRECISION, &
        proc, 77, comm, requestr(n), ierr)
    end do
    mglob = 0
    ibglob = 0
    ipar2g = mglob + 1
    ipar3g = mglob + 1
    kountg = 0
    !$OMP END SINGLE
    do
      !$OMP SINGLE
      iglob = 1
      kountg = kountg + 1
      !$OMP END SINGLE
      if (kountg.GT.nrecs_g) EXIT
      n = kountg
      proc = messproc_g(1, n)
      kl = klast_four(proc) - kfirst_four(proc) + 1
      do
        !$OMP CRITICAL(sndcrit)
        ig = iglob
        do
          if (iglob.gt.nGridFields) EXIT
          if (.not.surfGrid(iglob).or.kfirst_four(proc).eq.1) EXIT
          iglob = iglob + 1
        enddo
        ig = iglob
        ! TODO check warning ipar not initialized
        ! it was initialized above, but could be = 0 here ?
        if (iglob.le.nGridFields) then
          ib = ibglob
          if (surfGrid(ig)) then
            ipar = ipar2g
            ipar2g = ipar2g + 1
            if (ipar2g.gt.messproc_g(2, n)) then
              ipar2g = mglob + 1
              iglob = iglob + 1
            endif
            ibglob = ibglob + messages_g(2, ipar) - messages_g(1, ipar) + 1
          else
            ipar = ipar3g
            ipar3g = ipar3g + 1
            if (ipar3g.gt.messproc_g(2, n)) then
              ipar3g = mglob + 1
              iglob = iglob + 1
            endif
            ibglob = ibglob + (messages_g(2, ipar) - messages_g(1, ipar) + 1) * kl
          endif
        endif
        !$OMP END CRITICAL(sndcrit)
        if (ig.gt.nGridFields) EXIT
        ! TODO check warning ipar not initialized
        if (surfGrid(ig)) then
          g2 => Grid2d(ig)%p
          j = messages_g(3, ipar)
          do i = messages_g(1, ipar), messages_g(2, ipar)
            ib = ib + 1
            bufsend(ib) = g2(ibperij(i, j), jbperij(i, j))
          end do
        else
          g3 => Grid3d(ig)%p
          j = messages_g(3, ipar)
          do k = kfirst_four(proc), klast_four(proc)
            do i = messages_g(1, ipar), messages_g(2, ipar)
              ib = ib + 1
              bufsend(ib) = g3(ibperij(i, j), k, jbperij(i, j))
            end do
          end do
        end if
      end do
      !$OMP BARRIER
      !$OMP SINGLE
      call MPI_ISEND(bufsend(ibs(n)), ibs(n + 1) - ibs(n), MPI_DOUBLE_PRECISION, proc, 77, &
        comm, requests(n), ierr)
      mglob = messproc_g(2, n)
      ipar2g = mglob + 1
      ipar3g = mglob + 1
      !$OMP END SINGLE
    end do
    !$OMP SINGLE
    iglob = 1
    jglob = max(myfirstlat, myfirstlat_f)
    !$OMP END SINGLE
    !
    ! local values
    ! ------------
    if (max(myfirstlat, myfirstlat_f).le.min(mylastlat, mylastlat_f)) then
      do
        !$OMP CRITICAL(loccrit)
        ig = iglob
        do
          if (iglob.gt.nGridFields) EXIT
          if (.not.surfGrid(iglob).or.myfirstlev.eq.1) EXIT
          iglob = iglob + 1
        enddo
        ig = iglob
        if (iglob.le.nGridFields) then
          j = jglob
          if (jglob.ge.min(mylastlat, mylastlat_f)) then
            iglob = iglob + 1
            jglob = max(myfirstlat, myfirstlat_f)
          else
            jglob = jglob + 1
          endif
        endif
        !$OMP END CRITICAL(loccrit)
        if (ig.gt.nGridFields) EXIT

        if (surfGrid(ig)) then
          g2 => Grid2d(ig)%p
          v0 = previousJ(j) * nVertGrid + prevGrid(ig)
          do i = myfirstlon(j), mylastlon(j)
            Four(v0 + 1, i) = g2(ibPerIJ(i, j), jbPerIJ(i, j))
          end do
        else
          g3 => Grid3d(ig)%p
          v0 = previousJ(j) * nVertGrid + prevGrid(ig)
          do k = myfirstlev, mylastlev
            v0 = v0 + 1
            do i = myfirstlon(j), mylastlon(j)
              Four(v0, i) = g3(ibPerIJ(i, j), k, jbPerIJ(i, j))
            end do
          end do
        end if
      end do
    end if

    ! Nullify remaining Fourier

    !$OMP SINGLE
    do i = 2 * iMax + 1, dip1
      Four(:, i) = 0.0_p_r8
    end do

    do j = nVertGrid * myJMax_f + 1, dvdlj
      Four(j, :) = 0.0_p_r8
    end do
    kountg = 0
    !$OMP END SINGLE

    do
      !$OMP SINGLE
      kountg = kountg + 1
      !$OMP END SINGLE
      if (kountg.GT.nrecs_f) EXIT
      !$OMP SINGLE
      call MPI_WAITANY(nrecs_f, requestr(1), index, status, ierr)
      ksg = status(MPI_SOURCE)
      do l = 1, nrecs_f
        if (ksg.eq.messproc_f(1, l)) then
          nglob = l
          ibglob = ibr(nglob) - 1
          mglob = messproc_f(2, nglob - 1)
          EXIT
        endif
      enddo
      iglob = 1
      ipar2g = mglob + 1
      ipar3g = mglob + 1
      !$OMP END SINGLE
      n = nglob
      m = mglob
      !      kl =  klast_four(ksg) - kfirst_four(ksg) + 1
      kl = kmaxloc
      do
        !$OMP CRITICAL(reccrit)
        ig = iglob
        do
          if (iglob.gt.nGridFields) EXIT
          if (.not.surfGrid(iglob).or.havesurf) EXIT
          iglob = iglob + 1
        enddo
        ig = iglob
        ! TODO check warning ipar not initialized
        ! it was initialized above, but could be = 0 here ?
        if (iglob.le.nGridFields) then
          ib = ibglob
          if (surfGrid(ig)) then
            ipar = ipar2g
            ibglob = ibglob + messages_f(2, ipar) - messages_f(1, ipar) + 1
            ipar2g = ipar2g + 1
            if (ipar2g.gt.messproc_f(2, n)) then
              ipar2g = m + 1
              iglob = iglob + 1
            endif
          else
            ipar = ipar3g
            ibglob = ibglob + kl * (messages_f(2, ipar) - messages_f(1, ipar) + 1)
            ipar3g = ipar3g + 1
            if (ipar3g.gt.messproc_f(2, n)) then
              ipar3g = m + 1
              iglob = iglob + 1
            endif
          endif
        endif
        !$OMP END CRITICAL(reccrit)
        if (ig.gt.nGridFields) EXIT
        ! TODO check warning ipar not initialized
        if (surfGrid(ig)) then
          j = messages_f(3, ipar)
          v0 = previousJ(j) * nVertGrid + prevGrid(ig)
          do i = messages_f(1, ipar), messages_f(2, ipar)
            ib = ib + 1
            Four(v0 + 1, i) = bufrec(ib)
          end do
        else
          j = messages_f(3, ipar)
          v0 = previousJ(j) * nVertGrid + prevGrid(ig)
          do k = 1, kmaxloc
            do i = messages_f(1, ipar), messages_f(2, ipar)
              ib = ib + 1
              Four(k + v0, i) = bufrec(ib)
            end do
          end do
        end if
      end do
    end do

    !$OMP SINGLE
    call MPI_WAITALL(nrecs_g, requests(1), stat, ierr)
    !$OMP END SINGLE

  end subroutine DepositGrid


  subroutine FourToSpec()
    !# Fourier To Spectral
    !# ---
    !# @info
    !# **Brief:** Fourier To Spectral.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer :: lm
    integer :: m
    integer :: mn
    integer :: j
    integer :: jl
    integer :: js
    integer :: v
    integer :: k
    integer :: kdim
    integer :: ldim
    integer :: kn
    integer :: ks
    integer :: inin
    integer :: inin1
    integer :: inis
    integer :: inis1
    integer :: comm
    integer :: ierr
    integer :: index
    real(kind = p_r8) :: FoEv(djh, dv)
    real(kind = p_r8) :: FoOd(djh, dv)

    kdim = nVertSpec * myjMax_f * 2 * MMaxlocal
    ldim = nVertSpec * jMaxlocal_f * 2 * myMMax
    !$OMP SINGLE
    mnodes = 0
    lm = 0 ! not initialized, initializing with 0, because mnodes = 0 above
    ! TODO check how to initialize kn
    mGlob = 0
    if (dimsendbuf.lt.ldim * maxnodes_four) then
      dimsendbuf = ldim * maxnodes_four
      deallocate (bufsend)
      allocate (bufsend(dimsendbuf))
    endif
    if (dimrecbuf.lt.kdim * maxnodes_four) then
      dimrecbuf = kdim * maxnodes_four
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    !$OMP END SINGLE
    do
      !$OMP CRITICAL(lmcrit1)
      mGlob = mGlob + 1
      m = mGlob
      if (m.le.MMax) then
        kn = NodeHasM(m, mygroup_four)
        mnodes(kn) = mnodes(kn) + 1
        lm = mnodes(kn)
      endif
      !$OMP END CRITICAL(lmcrit1)
      if (m > MMax) EXIT
      if (kn.ne.myid_four) then
        !CDIR NODEP
        do j = max(jMinPerM(m), myfirstlat_f), min(mylastlat_f, jMaxPerM(m))
          do v = 1, nVertSpec
            jl = j - myfirstlat_f + 1
            bufrec((myjmax_f * 2 * (lm - 1) + (jl - 1)) * nVertGrid + v + kn * kdim) = &
              Four((jl - 1) * nVertGrid + v, 2 * m - 1)
            bufrec((myjmax_f * (2 * lm - 1) + (jl - 1)) * nVertGrid + v + kn * kdim) = &
              Four((jl - 1) * nVertGrid + v, 2 * m)
          end do
        end do
      else
        do j = max(jMinPerM(m), myfirstlat_f), min(mylastlat_f, jMaxPerM(m))
          do v = 1, nVertSpec
            jl = j - myfirstlat_f + 1
            bufsend((myjmax_f * 2 * (lm - 1) + (jl - 1)) * nVertGrid + v + kn * ldim) = &
              Four((jl - 1) * nVertGrid + v, 2 * m - 1)
            bufsend((myjmax_f * (2 * lm - 1) + (jl - 1)) * nVertGrid + v + kn * ldim) = &
              Four((jl - 1) * nVertGrid + v, 2 * m)
          end do
        end do
      end if
    end do
    !$OMP BARRIER
    !$OMP SINGLE
    comm = COMM_FOUR
    requestr(myid_four) = MPI_REQUEST_NULL
    requests(myid_four) = MPI_REQUEST_NULL
    do k = 0, MaxNodes_four - 1
      if(k.ne.myid_four) then
        call MPI_IRECV(bufsend(1 + k * ldim), 2 * nlatsinproc_f(k) * nvertspec * mymmax, &
          MPI_DOUBLE_PRECISION, k, 98, comm, requestr(k), ierr)
      endif
    enddo
    do k = 0, MaxNodes_four - 1
      if(k.ne.myid_four) then
        call MPI_ISEND(bufrec(1 + k * kdim), 2 * myjmax_f * nvertspec * MsPerProc(k), &
          MPI_DOUBLE_PRECISION, k, 98, comm, requests(k), ierr)
      endif
    enddo
    do k = 1, MaxNodes_four - 1
      call MPI_WAITANY(MaxNodes_four, requestr(0), index, status, ierr)
    enddo
    do k = 0, MaxNodes_four - 1
      if(k.ne.myid_four) then
        call MPI_WAIT(requests(k), status, ierr)
      endif
    enddo
    !$OMP END SINGLE
    FoEv = 0.0_p_r8
    FoOd = 0.0_p_r8

    !$OMP SINGLE
    mGlob = 0
    !$OMP END SINGLE
    do
      !$OMP CRITICAL(lmcrit2)
      mGlob = mGlob + 1
      lm = mGlob
      !$OMP END CRITICAL(lmcrit2)
      if (lm > myMMax) EXIT
      m = lm2m(lm)

      !CDIR NODEP
      do j = jMinPerM(m), jMaxHalf
        js = jMax + 1 - j
        kn = nodeHasJ_f(j)
        ks = nodeHasJ_f(js)
        inin = (nlatsinproc_f(kn) * 2 * (lm - 1) + j - firstlatinproc_f(kn)) * nVertSpec
        inin1 = inin + nlatsinproc_f(kn) * nVertSpec
        inis = (nlatsinproc_f(ks) * 2 * (lm - 1) + js - firstlatinproc_f(ks)) * nVertSpec
        inis1 = inis + nlatsinproc_f(ks) * nVertSpec
        do v = 1, nVertSpec
          FoEv(j, v) = bufsend(inin + v + kn * ldim) + bufsend(inis + v + ks * ldim)
          FoEv(j, v + nVertSpec) = bufsend(inin1 + v + kn * ldim) + bufsend(inis1 + v + ks * ldim)
          FoOd(j, v) = bufsend(inin + v + kn * ldim) - bufsend(inis + v + ks * ldim)
          FoOd(j, v + nVertSpec) = bufsend(inin1 + v + kn * ldim) - bufsend(inis1 + v + ks * ldim)
        end do
      end do
      do j = 1, dv
        do mn = FirstNEven(lm), FirstNEven(lm) + nEven(lm) - 1
          Spec(mn, j) = 0.0_p_r8
        end do
      end do
      call mmd (&
        LF2S(FirstNEven(lm), jMinPerM(m)), dlmn, &
        FoEv(jMinPerM(m), 1), djh, &
        Spec(FirstNEven(lm), 1), dlmn, dv, &
        nEven(lm), 2 * nVertSpec, jMaxHalf - jMinPerM(m) + 1)
      do j = 1, dv
        do mn = FirstNOdd(lm), FirstNOdd(lm) + nOdd(lm) - 1
          Spec(mn, j) = 0.0_p_r8
        end do
      end do
      call mmd (&
        LF2S(FirstNOdd(lm), jMinPerM(m)), dlmn, &
        FoOd(jMinPerM(m), 1), djh, &
        Spec(FirstNOdd(lm), 1), dlmn, dv, &
        nOdd(lm), 2 * nVertSpec, jMaxHalf - jMinPerM(m) + 1)
    end do
  end subroutine FourToSpec


  subroutine WithdrawSpectral(mnFirst, mnLast, mnExtFirst, mnExtLast)
    !# Withdraw Spectral
    !# ---
    !# @info
    !# **Brief:** Retrieve Spectral fields.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: mnFirst
    integer, intent(in) :: mnLast
    integer, intent(in) :: mnExtFirst
    integer, intent(in) :: mnExtLast

    integer :: is
    integer :: imn
    integer :: iv
    integer :: lastv
    real(kind = p_r8), pointer :: s1(:)
    real(kind = p_r8), pointer :: s2(:, :)

    do is = 1, nSpecFields

      if (surfSpec(is)) then
        s1 => Spec1d(is)%p
        lastv = prevSpec(is) + 1
        if (size(s1, 1) == 2 * myMNExtMax) then
          !CDIR NODEP
          do imn = mnExtFirst, mnExtLast
            s1(2 * imn - 1) = Spec(lmnExtMap(imn), lastv)
            s1(2 * imn) = Spec(lmnExtMap(imn), lastv + nVertSpec)
          end do
        else
          !CDIR NODEP
          do imn = mnFirst, mnLast
            s1(2 * imn - 1) = Spec(lmnMap(imn), lastv)
            s1(2 * imn) = Spec(lmnMap(imn), lastv + nVertSpec)
          end do
        end if
      else
        s2 => Spec2d(is)%p
        lastv = prevSpec(is)
        if (size(s2, 1) == 2 * myMNExtMax) then
          do iv = 1, kMaxloc
            !CDIR NODEP
            do imn = mnExtFirst, mnExtlast
              s2(2 * imn - 1, iv) = Spec(lmnExtMap(imn), iv + lastv)
              s2(2 * imn, iv) = Spec(lmnExtMap(imn), iv + lastv + nVertSpec)
            end do
          end do
        else
          do iv = 1, kMaxloc
            !CDIR NODEP
            do imn = mnFirst, mnLast
              s2(2 * imn - 1, iv) = Spec(lmnMap(imn), iv + lastv)
              s2(2 * imn, iv) = Spec(lmnMap(imn), iv + lastv + nVertSpec)
            end do
          end do
        end if
      end if
    end do

  end subroutine WithdrawSpectral


  subroutine Destroy()
    !# Destroy
    !# ---
    !# @info
    !# **Brief:** Removes the internal data structure created by CreateGridToSpec
    !# and created by CreateSpecToGrid.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    deallocate(Spec)
    deallocate(Four)
    if (allocated(mnodes)) then
      deallocate(mnodes)
    end if
    if (allocated(surfSpec)) then
      deallocate(surfSpec)
    end if
    if (allocated(prevSpec)) then
      deallocate(prevSpec)
    end if
    if (allocated(Spec1d)) then
      deallocate(Spec1d)
    end if
    if (allocated(Spec2d)) then
      deallocate(Spec2d)
    end if

    if (allocated(surfGrid)) then
      deallocate(surfGrid)
    end if
    if (allocated(prevGrid)) then
      deallocate(prevGrid)
    end if
    if (allocated(Grid2d)) then
      deallocate(Grid2d)
    end if
    if (allocated(Grid3d)) then
      deallocate(Grid3d)
    end if
    if (allocated(fieldForDelLam)) then
      deallocate(fieldForDelLam)
    end if
    if (allocated(prevVertDelLamSource)) then
      deallocate(prevVertDelLamSource)
    end if
  end subroutine Destroy


  subroutine mmd (A, lda, B, ldb, C, ldc, tdc, ni, nj, nk)
    !# mmd...?
    !# ---
    !# @info
    !# **Brief:** mmd...?.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: lda, ldb, ldc, tdc
    real(kind = p_r8), intent(in) :: A(lda, nk)
    real(kind = p_r8), intent(in) :: B(ldb, nj)
    real(kind = p_r8), intent(inout) :: C(ldc, tdc)
    integer i, j, k
    do i = 1, ni
      do j = 1, nj
        do k = 1, nk
          C(i, j) = C(i, j) + A(i, k) * B(k, j)
        end do
      end do
    end do
  end subroutine mmd


  subroutine mmt (A, lda, B, ldb, C, ldc, tdc, ni, nj, nk)
    !# mmt...?
    !# ---
    !# @info
    !# **Brief:** mmt...?.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: lda, ldb, ldc, tdc
    real(kind = p_r8), intent(in) :: A(lda, nk)
    real(kind = p_r8), intent(in) :: B(ldb, nj)
    real(kind = p_r8), intent(inout) :: C(ldc, tdc)
    integer i, j, k
    character(len = *), parameter :: h = "**(mmt)**"
    do i = 1, ni
      do j = 1, nj
        do k = 1, nk
          C(j, i) = C(j, i) + A(i, k) * B(k, j)
        end do
      end do
    end do
  end subroutine mmt


  subroutine InvFFTTrans (fInOut, ldInOut, tdInOut, n, lot, Trigs, nTrigs, Factors, nFactors)
    !# Computes inverse FFT of 'lot' sequences of 'n+1' input real data as rows
    !# of 'fin', dimensioned 'fin(ldin,tdin)'
    !# ---
    !# @info
    !# **Brief:** Input data is kept unchanged. Input values 'fin(n+2:ldin,:)'
    !# and 'fin(:,lot+1:tdin)' are not visited. Outputs 'lot' sequences of 'n'
    !# real data as rows of 'fout', dimensioned 'fout(ldout,tdout)'. Output
    !# values 'fout(n+1:ldout,:)' and 'fout(:,lot+1:tdout)' are set to 0.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ldInOut, tdInOut
    real(kind = p_r8), intent(inout) :: fInOut (ldInOut, tdInOut)
    integer, intent(in) :: n
    integer, intent(in) :: lot
    integer, intent(in) :: nTrigs
    real(kind = p_r8), intent(in) :: Trigs(nTrigs)
    integer, intent(in) :: nFactors
    integer, intent(in) :: Factors(nFactors)

    integer :: nh
    integer :: nfax
    integer :: la
    integer :: k
    logical :: ab2cd
    character(len = *), parameter :: h = "**(InvFFTTrans)**"
    real(kind = p_r8) :: a(lot, n / 2)
    real(kind = p_r8) :: b(lot, n / 2)
    real(kind = p_r8) :: c(lot, n / 2)
    real(kind = p_r8) :: d(lot, n / 2)

    nfax = Factors(1)
    nh = n / 2

    call SplitFourTrans (fInOut, a, b, ldInOut, tdInOut, n, nh, lot, Trigs, nTrigs)

    la = 1
    ab2cd = .true.
    do k = 1, nfax
      if (ab2cd) then
        call OnePass (a, b, c, d, lot, nh, Factors(k + 1), la, Trigs, nTrigs)
        ab2cd = .false.
      else
        call OnePass (c, d, a, b, lot, nh, Factors(k + 1), la, Trigs, nTrigs)
        ab2cd = .true.
      end if
      la = la * Factors(k + 1)
    end do

    if (ab2cd) then
      call JoinGridTrans (a, b, fInOut, ldInOut, tdInOut, nh, lot)
    else
      call JoinGridTrans (c, d, fInOut, ldInOut, tdInOut, nh, lot)
    end if
  end subroutine InvFFTTrans


  subroutine SplitFourTrans (fin, a, b, ldin, tdin, n, nh, lot, Trigs, nTrigs)
    !# Split Fourier Fields
    !# ---
    !# @info
    !# **Brief:** Split Fourier Fields.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ldin
    integer, intent(in) :: tdin
    integer, intent(in) :: n
    integer, intent(in) :: nh
    integer, intent(in) :: lot
    real(kind = p_r8), intent(in) :: fin(ldin, tdin)
    real(kind = p_r8), intent(OUT) :: a  (lot, nh)
    real(kind = p_r8), intent(OUT) :: b  (lot, nh)
    integer, intent(in) :: nTrigs
    real(kind = p_r8), intent(in) :: Trigs(nTrigs)

    integer :: i, j
    real(kind = p_r8) :: c, s

    !CDIR NODEP
    do i = 1, lot
      a(i, 1) = fin(i, 1) + fin(i, n + 1)
      b(i, 1) = fin(i, 1) - fin(i, n + 1)
    end do

    do j = 2, (nh + 1) / 2
      c = Trigs(n + 2 * j - 1)
      s = Trigs(n + 2 * j)
      !CDIR NODEP
      do i = 1, lot
        a(i, j) = (fin(i, 2 * j - 1) + fin(i, n + 3 - 2 * j)) &
          - (s * (fin(i, 2 * j - 1) - fin(i, n + 3 - 2 * j)) &
            + c * (fin(i, 2 * j) + fin(i, n + 4 - 2 * j)))
        a(i, nh + 2 - j) = (fin(i, 2 * j - 1) + fin(i, n + 3 - 2 * j)) &
          + (s * (fin(i, 2 * j - 1) - fin(i, n + 3 - 2 * j)) &
            + c * (fin(i, 2 * j) + fin(i, n + 4 - 2 * j)))
        b(i, j) = (c * (fin(i, 2 * j - 1) - fin(i, n + 3 - 2 * j)) &
          - s * (fin(i, 2 * j) + fin(i, n + 4 - 2 * j)))&
          + (fin(i, 2 * j) - fin(i, n + 4 - 2 * j))
        b(i, nh + 2 - j) = (c * (fin(i, 2 * j - 1) - fin(i, n + 3 - 2 * j)) &
          - s * (fin(i, 2 * j) + fin(i, n + 4 - 2 * j)))&
          - (fin(i, 2 * j) - fin(i, n + 4 - 2 * j))
      end do
    end do
    if ((nh>=2) .and. (mod(nh, 2)==0)) then
      !CDIR NODEP
      do i = 1, lot
        a(i, nh / 2 + 1) = 2.0_p_r8 * fin(i, nh + 1)
        b(i, nh / 2 + 1) = -2.0_p_r8 * fin(i, nh + 2)
      end do
    end if
  end subroutine SplitFourTrans


  subroutine JoinGridTrans (a, b, fout, ldout, tdout, nh, lot)
    !# Join Grid Fields
    !# ---
    !# @info
    !# **Brief:** Merge fundamental algorithm complex output into sequences of
    !# real numbers.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ldout, tdout
    integer, intent(in) :: nh
    integer, intent(in) :: lot
    real(kind = p_r8), intent(OUT) :: fout(ldout, tdout)
    real(kind = p_r8), intent(in) :: a   (lot, nh)
    real(kind = p_r8), intent(in) :: b   (lot, nh)

    integer :: i, j

    do j = 1, nh
      do i = 1, lot
        fout(i, 2 * j - 1) = a(i, j)
        fout(i, 2 * j) = b(i, j)
      end do
    end do

    !    fout(:,2*nh+1:tdout)=0.0_p_r8
    !    fout(lot+1:ldout, :)=0.0_p_r8
  end subroutine JoinGridTrans


  subroutine OnePass (a, b, c, d, lot, nh, ifac, la, Trigs, nTrigs)
    !# One Pass
    !# ---
    !# @info
    !# **Brief:** Single pass of fundamental algorithm.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: lot
    integer, intent(in) :: nh         ! = PROD(factor(1:K))
    real(kind = p_r8), intent(in) :: a(lot, nh)
    real(kind = p_r8), intent(in) :: b(lot, nh)
    real(kind = p_r8), intent(OUT) :: c(lot, nh)
    real(kind = p_r8), intent(OUT) :: d(lot, nh)
    integer, intent(in) :: ifac       ! = factor(k)
    integer, intent(in) :: la         ! = PROD(factor(1:k-1))
    integer, intent(in) :: nTrigs
    real(kind = p_r8), intent(in) :: Trigs(nTrigs)

    integer :: m
    integer :: jump
    integer :: i, j, k
    integer :: ia, ja
    integer :: ib, jb, kb
    integer :: ic, jc, kc
    integer :: id, jd, kd
    integer :: ie, je, ke
    real(kind = p_r8) :: c1, s1
    real(kind = p_r8) :: c2, s2
    real(kind = p_r8) :: c3, s3
    real(kind = p_r8) :: c4, s4
    real(kind = p_r8) :: wka, wkb
    real(kind = p_r8) :: wksina, wksinb
    real(kind = p_r8) :: wkaacp, wkbacp
    real(kind = p_r8) :: wkaacm, wkbacm

    m = nh / ifac
    jump = (ifac - 1) * la

    ia = 0
    ib = m
    ic = 2 * m
    id = 3 * m
    ie = 4 * m

    ja = 0
    jb = la
    jc = 2 * la
    jd = 3 * la
    je = 4 * la

    if (ifac == 2) then
      do j = 1, la
        !CDIR NODEP
        do i = 1, lot
          c(i, j + ja) = a(i, j + ia) + a(i, j + ib)
          c(i, j + jb) = a(i, j + ia) - a(i, j + ib)
          d(i, j + ja) = b(i, j + ia) + b(i, j + ib)
          d(i, j + jb) = b(i, j + ia) - b(i, j + ib)
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        c1 = Trigs(kb + 1)
        s1 = Trigs(kb + 2)
        ja = ja + jump
        jb = jb + jump
        do j = k + 1, k + la
          !CDIR NODEP
          do i = 1, lot
            wka = a(i, j + ia) - a(i, j + ib)
            c(i, j + ja) = a(i, j + ia) + a(i, j + ib)
            wkb = b(i, j + ia) - b(i, j + ib)
            d(i, j + ja) = b(i, j + ia) + b(i, j + ib)
            c(i, j + jb) = c1 * wka - s1 * wkb
            d(i, j + jb) = s1 * wka + c1 * wkb
          end do
        end do
      end do
    ELSEIF (ifac == 3) then
      do j = 1, la
        !CDIR  NODEP
        do i = 1, lot
          wka = a(i, j + ib) + a(i, j + ic)
          wksina = sin60 * (a(i, j + ib) - a(i, j + ic))
          wkb = b(i, j + ib) + b(i, j + ic)
          wksinb = sin60 * (b(i, j + ib) - b(i, j + ic))
          c(i, j + ja) = a(i, j + ia) + wka
          c(i, j + jb) = (a(i, j + ia) - 0.5_p_r8 * wka) - wksinb
          c(i, j + jc) = (a(i, j + ia) - 0.5_p_r8 * wka) + wksinb
          d(i, j + ja) = b(i, j + ia) + wkb
          d(i, j + jb) = (b(i, j + ia) - 0.5_p_r8 * wkb) + wksina
          d(i, j + jc) = (b(i, j + ia) - 0.5_p_r8 * wkb) - wksina
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        kc = kb + kb
        c1 = Trigs(kb + 1)
        s1 = Trigs(kb + 2)
        c2 = Trigs(kc + 1)
        s2 = Trigs(kc + 2)
        ja = ja + jump
        jb = jb + jump
        jc = jc + jump
        do j = k + 1, k + la
          !CDIR NODEP
          do i = 1, lot
            wka = a(i, j + ib) + a(i, j + ic)
            wksina = sin60 * (a(i, j + ib) - a(i, j + ic))
            wkb = b(i, j + ib) + b(i, j + ic)
            wksinb = sin60 * (b(i, j + ib) - b(i, j + ic))
            c(i, j + ja) = a(i, j + ia) + wka
            d(i, j + ja) = b(i, j + ia) + wkb
            c(i, j + jb) = c1 * ((a(i, j + ia) - 0.5_p_r8 * wka) - wksinb) &
              - s1 * ((b(i, j + ia) - 0.5_p_r8 * wkb) + wksina)
            d(i, j + jb) = s1 * ((a(i, j + ia) - 0.5_p_r8 * wka) - wksinb) &
              + c1 * ((b(i, j + ia) - 0.5_p_r8 * wkb) + wksina)
            c(i, j + jc) = c2 * ((a(i, j + ia) - 0.5_p_r8 * wka) + wksinb) &
              - s2 * ((b(i, j + ia) - 0.5_p_r8 * wkb) - wksina)
            d(i, j + jc) = s2 * ((a(i, j + ia) - 0.5_p_r8 * wka) + wksinb) &
              + c2 * ((b(i, j + ia) - 0.5_p_r8 * wkb) - wksina)
          end do
        end do
      end do
    ELSEIF (ifac == 4) then
      do j = 1, la
        !CDIR NODEP
        do i = 1, lot
          wkaacp = a(i, j + ia) + a(i, j + ic)
          wkaacm = a(i, j + ia) - a(i, j + ic)
          wkbacp = b(i, j + ia) + b(i, j + ic)
          wkbacm = b(i, j + ia) - b(i, j + ic)
          c(i, j + ja) = wkaacp + (a(i, j + ib) + a(i, j + id))
          c(i, j + jc) = wkaacp - (a(i, j + ib) + a(i, j + id))
          d(i, j + jb) = wkbacm + (a(i, j + ib) - a(i, j + id))
          d(i, j + jd) = wkbacm - (a(i, j + ib) - a(i, j + id))
          d(i, j + ja) = wkbacp + (b(i, j + ib) + b(i, j + id))
          d(i, j + jc) = wkbacp - (b(i, j + ib) + b(i, j + id))
          c(i, j + jb) = wkaacm - (b(i, j + ib) - b(i, j + id))
          c(i, j + jd) = wkaacm + (b(i, j + ib) - b(i, j + id))
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        kc = kb + kb
        kd = kc + kb
        c1 = Trigs(kb + 1)
        s1 = Trigs(kb + 2)
        c2 = Trigs(kc + 1)
        s2 = Trigs(kc + 2)
        c3 = Trigs(kd + 1)
        s3 = Trigs(kd + 2)
        ja = ja + jump
        jb = jb + jump
        jc = jc + jump
        jd = jd + jump
        do j = k + 1, k + la
          !CDIR NODEP
          do i = 1, lot
            wkaacp = a(i, j + ia) + a(i, j + ic)
            wkbacp = b(i, j + ia) + b(i, j + ic)
            wkaacm = a(i, j + ia) - a(i, j + ic)
            wkbacm = b(i, j + ia) - b(i, j + ic)
            c(i, j + ja) = wkaacp + (a(i, j + ib) + a(i, j + id))
            d(i, j + ja) = wkbacp + (b(i, j + ib) + b(i, j + id))
            c(i, j + jc) = c2 * (wkaacp - (a(i, j + ib) + a(i, j + id))) &
              - s2 * (wkbacp - (b(i, j + ib) + b(i, j + id)))
            d(i, j + jc) = s2 * (wkaacp - (a(i, j + ib) + a(i, j + id))) &
              + c2 * (wkbacp - (b(i, j + ib) + b(i, j + id)))
            c(i, j + jb) = c1 * (wkaacm - (b(i, j + ib) - b(i, j + id))) &
              - s1 * (wkbacm + (a(i, j + ib) - a(i, j + id)))
            d(i, j + jb) = s1 * (wkaacm - (b(i, j + ib) - b(i, j + id))) &
              + c1 * (wkbacm + (a(i, j + ib) - a(i, j + id)))
            c(i, j + jd) = c3 * (wkaacm + (b(i, j + ib) - b(i, j + id))) &
              - s3 * (wkbacm - (a(i, j + ib) - a(i, j + id)))
            d(i, j + jd) = s3 * (wkaacm + (b(i, j + ib) - b(i, j + id))) &
              + c3 * (wkbacm - (a(i, j + ib) - a(i, j + id)))
          end do
        end do
      end do
    ELSEIF (ifac == 5) then
      do j = 1, la
        !CDIR NODEP
        do i = 1, lot
          c(i, j + ja) = a(i, j + ia) + (a(i, j + ib) + a(i, j + ie)) + (a(i, j + ic) + a(i, j + id))
          d(i, j + ja) = b(i, j + ia)&
            + (b(i, j + ib) + b(i, j + ie))&
            + (b(i, j + ic) + b(i, j + id))
          c(i, j + jb) = (a(i, j + ia)&
            + cos72 * (a(i, j + ib) + a(i, j + ie))&
            - cos36 * (a(i, j + ic) + a(i, j + id)))&
            - (sin72 * (b(i, j + ib) - b(i, j + ie))&
              + sin36 * (b(i, j + ic) - b(i, j + id)))
          c(i, j + je) = (a(i, j + ia)&
            + cos72 * (a(i, j + ib) + a(i, j + ie))&
            - cos36 * (a(i, j + ic) + a(i, j + id)))&
            + (sin72 * (b(i, j + ib) - b(i, j + ie))&
              + sin36 * (b(i, j + ic) - b(i, j + id)))
          d(i, j + jb) = (b(i, j + ia)&
            + cos72 * (b(i, j + ib) + b(i, j + ie))&
            - cos36 * (b(i, j + ic) + b(i, j + id)))&
            + (sin72 * (a(i, j + ib) - a(i, j + ie))&
              + sin36 * (a(i, j + ic) - a(i, j + id)))
          d(i, j + je) = (b(i, j + ia)&
            + cos72 * (b(i, j + ib) + b(i, j + ie))&
            - cos36 * (b(i, j + ic) + b(i, j + id)))&
            - (sin72 * (a(i, j + ib) - a(i, j + ie))&
              + sin36 * (a(i, j + ic) - a(i, j + id)))
          c(i, j + jc) = (a(i, j + ia)&
            - cos36 * (a(i, j + ib) + a(i, j + ie))&
            + cos72 * (a(i, j + ic) + a(i, j + id)))&
            - (sin36 * (b(i, j + ib) - b(i, j + ie))&
              - sin72 * (b(i, j + ic) - b(i, j + id)))
          c(i, j + jd) = (a(i, j + ia)&
            - cos36 * (a(i, j + ib) + a(i, j + ie))&
            + cos72 * (a(i, j + ic) + a(i, j + id)))&
            + (sin36 * (b(i, j + ib) - b(i, j + ie))&
              - sin72 * (b(i, j + ic) - b(i, j + id)))
          d(i, j + jc) = (b(i, j + ia)&
            - cos36 * (b(i, j + ib) + b(i, j + ie))&
            + cos72 * (b(i, j + ic) + b(i, j + id)))&
            + (sin36 * (a(i, j + ib) - a(i, j + ie))&
              - sin72 * (a(i, j + ic) - a(i, j + id)))
          d(i, j + jd) = (b(i, j + ia)&
            - cos36 * (b(i, j + ib) + b(i, j + ie))&
            + cos72 * (b(i, j + ic) + b(i, j + id)))&
            - (sin36 * (a(i, j + ib) - a(i, j + ie))&
              - sin72 * (a(i, j + ic) - a(i, j + id)))
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        kc = kb + kb
        kd = kc + kb
        ke = kd + kb
        c1 = Trigs(kb + 1)
        s1 = Trigs(kb + 2)
        c2 = Trigs(kc + 1)
        s2 = Trigs(kc + 2)
        c3 = Trigs(kd + 1)
        s3 = Trigs(kd + 2)
        c4 = Trigs(ke + 1)
        s4 = Trigs(ke + 2)
        ja = ja + jump
        jb = jb + jump
        jc = jc + jump
        jd = jd + jump
        je = je + jump
        do j = k + 1, k + la
          !CDIR NODEP
          do i = 1, lot
            c(i, j + ja) = a(i, j + ia)&
              + (a(i, j + ib) + a(i, j + ie))&
              + (a(i, j + ic) + a(i, j + id))
            d(i, j + ja) = b(i, j + ia)&
              + (b(i, j + ib) + b(i, j + ie))&
              + (b(i, j + ic) + b(i, j + id))
            c(i, j + jb) = c1 * ((a(i, j + ia)&
              + cos72 * (a(i, j + ib) + a(i, j + ie))&
              - cos36 * (a(i, j + ic) + a(i, j + id)))&
              - (sin72 * (b(i, j + ib) - b(i, j + ie))&
                + sin36 * (b(i, j + ic) - b(i, j + id))))&
              - s1 * ((b(i, j + ia)&
                + cos72 * (b(i, j + ib) + b(i, j + ie))&
                - cos36 * (b(i, j + ic) + b(i, j + id)))&
                + (sin72 * (a(i, j + ib) - a(i, j + ie))&
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + jb) = s1 * ((a(i, j + ia)&
              + cos72 * (a(i, j + ib) + a(i, j + ie))&
              - cos36 * (a(i, j + ic) + a(i, j + id)))&
              - (sin72 * (b(i, j + ib) - b(i, j + ie))&
                + sin36 * (b(i, j + ic) - b(i, j + id))))&
              + c1 * ((b(i, j + ia)&
                + cos72 * (b(i, j + ib) + b(i, j + ie))&
                - cos36 * (b(i, j + ic) + b(i, j + id)))&
                + (sin72 * (a(i, j + ib) - a(i, j + ie)) &
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            c(i, j + je) = c4 * ((a(i, j + ia)&
              + cos72 * (a(i, j + ib) + a(i, j + ie))&
              - cos36 * (a(i, j + ic) + a(i, j + id)))&
              + (sin72 * (b(i, j + ib) - b(i, j + ie)) &
                + sin36 * (b(i, j + ic) - b(i, j + id)))) &
              - s4 * ((b(i, j + ia)&
                + cos72 * (b(i, j + ib) + b(i, j + ie))&
                - cos36 * (b(i, j + ic) + b(i, j + id)))&
                - (sin72 * (a(i, j + ib) - a(i, j + ie))&
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + je) = s4 * ((a(i, j + ia)&
              + cos72 * (a(i, j + ib) + a(i, j + ie))&
              - cos36 * (a(i, j + ic) + a(i, j + id)))&
              + (sin72 * (b(i, j + ib) - b(i, j + ie))&
                + sin36 * (b(i, j + ic) - b(i, j + id))))&
              + c4 * ((b(i, j + ia)&
                + cos72 * (b(i, j + ib) + b(i, j + ie))&
                - cos36 * (b(i, j + ic) + b(i, j + id))) &
                - (sin72 * (a(i, j + ib) - a(i, j + ie))&
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            c(i, j + jc) = c2 * ((a(i, j + ia)&
              - cos36 * (a(i, j + ib) + a(i, j + ie))&
              + cos72 * (a(i, j + ic) + a(i, j + id)))&
              - (sin36 * (b(i, j + ib) - b(i, j + ie))&
                - sin72 * (b(i, j + ic) - b(i, j + id))))&
              - s2 * ((b(i, j + ia)&
                - cos36 * (b(i, j + ib) + b(i, j + ie))&
                + cos72 * (b(i, j + ic) + b(i, j + id)))&
                + (sin36 * (a(i, j + ib) - a(i, j + ie))&
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + jc) = s2 * ((a(i, j + ia)&
              - cos36 * (a(i, j + ib) + a(i, j + ie))&
              + cos72 * (a(i, j + ic) + a(i, j + id)))&
              - (sin36 * (b(i, j + ib) - b(i, j + ie))&
                - sin72 * (b(i, j + ic) - b(i, j + id))))&
              + c2 * ((b(i, j + ia)&
                - cos36 * (b(i, j + ib) + b(i, j + ie))&
                + cos72 * (b(i, j + ic) + b(i, j + id)))&
                + (sin36 * (a(i, j + ib) - a(i, j + ie))&
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
            c(i, j + jd) = c3 * ((a(i, j + ia)&
              - cos36 * (a(i, j + ib) + a(i, j + ie))&
              + cos72 * (a(i, j + ic) + a(i, j + id)))&
              + (sin36 * (b(i, j + ib) - b(i, j + ie))&
                - sin72 * (b(i, j + ic) - b(i, j + id))))&
              - s3 * ((b(i, j + ia)&
                - cos36 * (b(i, j + ib) + b(i, j + ie))&
                + cos72 * (b(i, j + ic) + b(i, j + id)))&
                - (sin36 * (a(i, j + ib) - a(i, j + ie))&
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + jd) = s3 * ((a(i, j + ia)&
              - cos36 * (a(i, j + ib) + a(i, j + ie))&
              + cos72 * (a(i, j + ic) + a(i, j + id)))&
              + (sin36 * (b(i, j + ib) - b(i, j + ie))&
                - sin72 * (b(i, j + ic) - b(i, j + id))))&
              + c3 * ((b(i, j + ia)&
                - cos36 * (b(i, j + ib) + b(i, j + ie))&
                + cos72 * (b(i, j + ic) + b(i, j + id)))&
                - (sin36 * (a(i, j + ib) - a(i, j + ie))&
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
          end do
        end do
      end do
    endif
  end subroutine OnePass


  subroutine CreateFFT(nIn, Factors, nFactors, Trigs, nTrigs)
    !# Creates FFT
    !# ---
    !# @info
    !# **Brief:** Allocates and computes intermediate values used by the FFT for
    !# sequences of size nIn. If size is not in the form 2 * 2**i * 3**j * 5**k,
    !# with at least one of i,j,k/=0, stops and prints (stdout) the next possible
    !# size.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: nIn
    integer, pointer :: Factors(:)
    integer, intent(OUT) :: nFactors
    real(kind = p_r8), pointer :: Trigs(:)
    integer, intent(OUT) :: nTrigs
    character(len = 15), parameter :: h = "**(CreateFFT)**" ! header
    call Factorize  (nIn, Factors, nBase, Base)
    nFactors = size(Factors)
    call TrigFactors(nIn, Trigs)
    nTrigs = size(Trigs)
  end subroutine CreateFFT


  subroutine DestroyFFT(Factors, Trigs)
    !# Destroy FFT
    !# ---
    !# @info
    !# **Brief:** Dealocates input area, returning NULL pointers.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, pointer :: Factors(:)
    real(kind = p_r8), pointer :: Trigs(:)
    character(len = 16), parameter :: h = "**(DestroyFFT)**" ! header
    deallocate(Factors); nullify(Factors)
    deallocate(Trigs); nullify(Trigs)
  end subroutine DestroyFFT


  subroutine Factorize (nIn, Factors, nBase, Base)
    !# Factorizes
    !# ---
    !# @info
    !# **Brief:** Factorizes nIn/2 in powers of 4, 3, 2, 5, if possible. Otherwise,
    !# stops with error message  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: nIn
    ! As variaveis Base e nBase sao globais ao modulo, portanto seria desnecessaria
    ! a passagem por parametros. Entretanto para rodar com a chave de compilacao
    ! "-openmp" ligada no TX7, essas variaveis sao parametros de entrada. (???)
    integer, intent(in) :: nBase
    integer, intent(in) :: Base(nBase)
    integer, pointer :: Factors(:)
    character(len = 15), parameter :: h = "**(Factorize)**" ! header
    character(len = 15) :: c ! character representation of integer
    integer :: Powers(nBase)
    integer :: nOut
    integer :: sumPowers
    integer :: ifac
    integer :: i
    integer :: j
    integer :: left ! portion of nOut/2 yet to be factorized

    nOut = NextSizeFFT(nIn)

    if (nIn /= nOut) then
      write(c, "(i15)") nIn
      write(p_nfprt, "(a,' FFT size = ',a,' not factorizable ')")&
        h, trim(adjustl(c))
      write(c, "(i15)") nOut
      write(p_nfprt, "(a,' Next factorizable FFT size is ',a)")&
        h, trim(adjustl(c))
      stop
    end if

    ! Loop over evens, starting from nOut, getting factors of nOut/2

    left = nOut / 2
    Powers = 0

    ! factorize nOut/2

    do i = 1, nBase
      do
        if (mod(left, Base(i)) == 0) then
          Powers(i) = Powers(i) + 1
          left = left / Base(i)
        else
          EXIT
        end if
      end do
    end do

    sumPowers = sum(Powers)
    allocate (Factors(sumPowers + 1))
    Factors(1) = sumPowers
    ifac = 1
    do i = 1, nBase
      j = Permutation(i)
      Factors(ifac + 1:ifac + Powers(j)) = Base(j)
      ifac = ifac + Powers(j)
    end do
  end subroutine Factorize


  subroutine TrigFactors (nIn, Trigs)
    !# Trig Factors
    !# ---
    !# @info
    !# **Brief:** Sin and Cos required to compute FFT of size nIn. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: nIn
    real(kind = p_r8), pointer :: Trigs(:)
    integer :: nn, nh, i
    real(kind = p_r8) :: pi, del, angle

    nn = nIn / 2
    nh = (nn + 1) / 2
    allocate (Trigs(2 * (nn + nh)))

    pi = 2.0_p_r8 * asin(1.0_p_r8)
    del = (2 * pi) / real(nn, p_r8)

    do i = 1, 2 * nn, 2
      angle = 0.5_p_r8 * real(i - 1, p_r8) * del
      Trigs(i) = cos(angle)
      Trigs(i + 1) = sin(angle)
    end do

    del = 0.5_p_r8 * del
    do i = 1, 2 * nh, 2
      angle = 0.5_p_r8 * real(i - 1, p_r8) * del
      Trigs(2 * nn + i) = cos(angle)
      Trigs(2 * nn + i + 1) = sin(angle)
    end do
  end subroutine TrigFactors


  function NextSizeFFT(nIn) result(nOUT)
    !# Next Size FFT
    !# ---
    !# @info
    !# **Brief:** Returns the smallest size of FFT equal or larger than input
    !# argument. Smallest integer >= input in the form 2 * 2**i * 3**j * 5**k  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: nIn
    integer :: nOut
    character(len = 22), parameter :: h = "**(NextSizeFFT)**" ! header
    real(kind = p_r8), parameter :: limit = huge(nIn) - 1   ! Maximum representable integer
    character(len = 15) :: charNIn ! character representation of nIn
    integer :: i
    integer :: left ! portion of nOut/2 yet to be factorized

    ! nOut = positive even

    if (nIn <= 0) then
      write (charNIn, "(i15)") nIn
      write(p_nfprt, "(a,' Meaningless FFT size='a)")h, trim(adjustl(charNIn))
      stop
    else if (mod(nIn, 2) == 0) then
      nOut = nIn
    else
      nOut = nIn + 1
    end if

    ! Loop over evens, starting from nOut, looking for
    ! next factorizable even/2

    do
      left = nOut / 2

      ! factorize nOut/2

      do i = 1, nBase
        do
          if (mod(left, Base(i)) == 0) then
            left = left / Base(i)
          else
            EXIT
          end if
        end do
      end do

      if (left == 1) then
        EXIT
      else if (nOut < limit) then
        nOut = nOut + 2
      else
        write (charNIn, "(i15)") nIn
        write(p_nfprt, "(a,' Next factorizable FFT size > ',a,&
          &' is not representable in this machine')")&
          h, trim(adjustl(charNIn))
        stop
      end if
    end do
  end function NextSizeFFT


  subroutine DirFFTTrans (fInOut, ldInOut, tdInOut, n, lot, Trigs, nTrigs, Factors, nFactors)
    !# Computes direct FFT of 'lot' sequences of 'n' input real data as rows of
    !# 'fin', dimensioned 'fin(ldin,tdin)'
    !# ---
    !# @info
    !# **Brief:** Input data is kept unchanged. Input values 'fin(n+1:ldin,:)'
    !# and 'fin(:,lot+1:tdin)' are not visited. Outputs 'lot' sequences of 'n+1'
    !# real data as rows of 'fout', dimensioned 'fout(ldout,tdout)'. Output values
    !# 'fout(:,lot+1:tdout)' and 'fout(n+2:ldout,:)' are set to 0. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ldInOut, tdInOut
    real(kind = p_r8), intent(inout) :: fInOut (ldInOut, tdInOut)
    integer, intent(in) :: n
    integer, intent(in) :: lot
    integer, intent(in) :: nTrigs
    real(kind = p_r8), intent(in) :: Trigs(nTrigs)
    integer, intent(in) :: nFactors
    integer, intent(in) :: Factors(nFactors)

    integer :: nh
    integer :: nfax
    integer :: la
    integer :: k
    logical :: ab2cd
    character(len = 12), parameter :: h = "**(Dir)**"
    real(kind = p_r8) :: a(lot, n / 2)
    real(kind = p_r8) :: b(lot, n / 2)
    real(kind = p_r8) :: c(lot, n / 2)
    real(kind = p_r8) :: d(lot, n / 2)

    nfax = Factors(1)
    nh = n / 2

    call SplitGridTrans (fInOut, a, b, ldInOut, tdInOut, nh, lot)

    la = 1
    ab2cd = .true.
    do k = 1, nfax
      if (ab2cd) then
        call OnePass (a, b, c, d, lot, nh, Factors(k + 1), la, Trigs, nTrigs)
        ab2cd = .false.
      else
        call OnePass (c, d, a, b, lot, nh, Factors(k + 1), la, Trigs, nTrigs)
        ab2cd = .true.
      end if
      la = la * Factors(k + 1)
    end do

    if (ab2cd) then
      call JoinFourTrans (a, b, fInOut, ldInOut, tdInOut, n, nh, lot, Trigs, nTrigs)
    else
      call JoinFourTrans (c, d, fInOut, ldInOut, tdInOut, n, nh, lot, Trigs, nTrigs)
    end if
  end subroutine DirFFTTrans


  subroutine SplitGridTrans (fin, a, b, ldin, tdin, nh, lot)
    !# Split Grid Fields
    !# ---
    !# @info
    !# **Brief:** Split space domain real input into complex pairs to feed
    !# fundamental algorithm. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ldin
    integer, intent(in) :: tdin
    integer, intent(in) :: nh
    integer, intent(in) :: lot
    real(kind = p_r8), intent(in) :: fin(ldin, tdin)
    real(kind = p_r8), intent(OUT) :: a  (lot, nh)
    real(kind = p_r8), intent(OUT) :: b  (lot, nh)

    integer :: i, j

    do j = 1, nh
      do i = 1, lot
        a(i, j) = fin(i, 2 * j - 1)
        b(i, j) = fin(i, 2 * j)
      end do
    end do

  end subroutine SplitGridTrans


  subroutine JoinFourTrans (a, b, fout, ldout, tdout, n, nh, lot, Trigs, nTrigs)
    !# Join Fourier Fields
    !# ---
    !# @info
    !# **Brief:** Unscramble frequency domain complex output from fundamental
    !# algorithm into real sequences. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ldout, tdout
    integer, intent(in) :: n
    integer, intent(in) :: nh
    integer, intent(in) :: lot
    real(kind = p_r8), intent(OUT) :: fout(ldout, tdout)
    real(kind = p_r8), intent(in) :: a   (lot, nh)
    real(kind = p_r8), intent(in) :: b   (lot, nh)
    integer, intent(in) :: nTrigs
    real(kind = p_r8), intent(in) :: Trigs(nTrigs)

    integer :: i, j
    real(kind = p_r8) :: scale, scalh
    real(kind = p_r8) :: c, s

    scale = 1.0_p_r8 / real(n, p_r8)
    scalh = 0.5_p_r8 * scale

    !cdir nodep
    do i = 1, lot
      fout(i, 1) = scale * (a(i, 1) + b(i, 1))
      fout(i, n + 1) = scale * (a(i, 1) - b(i, 1))
      fout(i, 2) = 0.0_p_r8
    end do

    do j = 2, (nh + 1) / 2
      c = Trigs(n + 2 * j - 1)
      s = Trigs(n + 2 * j)
      !cdir nodep
      do i = 1, lot
        fout(i, 2 * j - 1) = scalh * ((a(i, j) + a(i, nh + 2 - j)) &
          + (c * (b(i, j) + b(i, nh + 2 - j))&
            + s * (a(i, j) - a(i, nh + 2 - j))))
        fout(i, n + 3 - 2 * j) = scalh * ((a(i, j) + a(i, nh + 2 - j)) &
          - (c * (b(i, j) + b(i, nh + 2 - j))&
            + s * (a(i, j) - a(i, nh + 2 - j))))
        fout(i, 2 * j) = scalh * ((c * (a(i, j) - a(i, nh + 2 - j))&
          - s * (b(i, j) + b(i, nh + 2 - j)))&
          + (b(i, nh + 2 - j) - b(i, j)))
        fout(i, n + 4 - 2 * j) = scalh * ((c * (a(i, j) - a(i, nh + 2 - j))&
          - s * (b(i, j) + b(i, nh + 2 - j)))&
          - (b(i, nh + 2 - j) - b(i, j)))
      end do
    end do

    if ((nh>=2) .and. (mod(nh, 2)==0)) then
      !cdir nodep
      do i = 1, lot
        fout(i, nh + 1) = scale * a(i, nh / 2 + 1)
        fout(i, nh + 2) = -scale * b(i, nh / 2 + 1)
      end do
    end if

  end subroutine JoinFourTrans


  subroutine Clear_Transform()
    !# Cleans Transform
    !# ---
    !# @info
    !# **Brief:** Cleans Transform. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin

    if (allocated(Spec))  deallocate (Spec)
    if (allocated(Four))  deallocate (Four)
    if (allocated(mnodes)) deallocate (mnodes)
    if (allocated(requests)) deallocate (requests)
    if (allocated(requestr)) deallocate (requestr)
    if (allocated(status)) deallocate (status)
    if (allocated(stat)) deallocate (stat)

    if (allocated(surfSpec)) deallocate (surfSpec)! true iff Surface Spectral Field
    if (allocated(prevSpec)) deallocate (prevSpec)! prior to first real vertical of this field at all internal arrays
    if (allocated(Spec1d)) deallocate (Spec1d)    ! points to Surface Spectral Field
    if (allocated(Spec2d)) deallocate (Spec2d)    ! points to Full Spectral Field

    if (allocated(surfGrid)) deallocate (surfGrid)! true iff Surface Grid Field
    if (allocated(prevGrid)) deallocate (prevGrid)! prior to first real vertical of this field at all internal arrays
    if (allocated(Grid2d)) deallocate (Grid2d)    ! points to Surface Grid Field
    if (allocated(Grid3d)) deallocate (Grid3d)    ! points to Full Grid Field
    if (allocated(fieldForDelLam)) deallocate (fieldForDelLam)  ! true iff this field position stores Lambda Derivative

    if (allocated(prevVertDelLamSource)) deallocate (prevVertDelLamSource)  ! source of Lambda Derivative

    if (allocated(nEven)) deallocate (nEven)
    if (allocated(dnEven)) deallocate (dnEven)
    if (allocated(firstNEven)) deallocate (firstNEven)
    if (allocated(nOdd)) deallocate (nOdd)
    if (allocated(dnOdd)) deallocate (dnOdd)
    if (allocated(firstNOdd)) deallocate (firstNOdd)
    if (allocated(lmnExtMap)) deallocate (lmnExtMap)
    if (allocated(lmnMap)) deallocate (lmnMap)
    if (allocated(lmnZero)) deallocate (lmnZero)

    if (allocated(previousJ)) deallocate (previousJ)

    if (allocated(BlockFFT)) deallocate (BlockFFT)
    if (allocated(BlockFFT)) deallocate (BlockFFT)
    if (allocated(LS2F)) deallocate (LS2F)
    if (allocated(LF2S)) deallocate (LF2S)
    if (allocated(consIm)) deallocate (consIm)
    if (allocated(consRe)) deallocate (consRe)
  end subroutine Clear_Transform

end module Mod_Transform

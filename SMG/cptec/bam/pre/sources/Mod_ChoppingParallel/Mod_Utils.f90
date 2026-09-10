!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Utils </br></br>
!#
!# **Brief**: Module computes and stores Associated Legendre Functions and Epslon. </br>
!# 
!# Module exports two routines:
!# <ul type="disc">
!#  <li>CreateAssocLegFunc initializes module and compute functions; </li>
!#  <li>DestroyAssocLegFunc destroys module; </li>
!# </ul>
!# 
!# Module usage: </br>
!# CreateAssocLegFunc should be invoked once, prior to any other </br>
!# module routine. It computes and hides the function values. </br>
!# DestroyAssocLegFunc destroys all internal info. </br>
!# 
!# Module use values exported by Sizes and procedures from Auxiliary </br>
!#
!# CreateAssocLegFunc </br>
!# DestroyAssocLegFunc </br>
!# 
!# CreateGaussQuad    ------------------| CreateLegPol </br>
!# 
!# CreateGridValues </br>
!# 
!# DestroyGaussQuad   ------------------| DestroyLegPol </br>
!# 
!# iminv </br>
!# 
!# Rg                 ------------------| Balanc </br>
!#                                  | 
!#                                  | Orthes </br>
!#                                  | 
!#                                  | Ortran </br>
!#                                  | 
!#                                  | Hqr2 ------ Hqr3 </br>
!#                                  | 
!#                                  | Balbak </br>
!#                                  | 
!#                                  | Znorma </br>
!# 
!# Tql2 </br>
!# Tred2 </br>
!# tmstmp2 </br>
!# InitTimeStamp </br>
!# 
!# TimeStamp  --------------------------| caldat </br>
!# 
!# IBJBtoIJ_R (Interface) </br>
!# IJtoIBJB_R (Interface) </br>
!# IBJBtoIJ_I (Interface) </br>
!# IJtoIBJB_I (Interface) </br>
!# 
!# SplineIJtoIBJB_R2D (Interface) ------| CyclicCubicSpline </br>
!# 
!# SplineIBJBtoIJ_R2D (Interface) ------| CyclicCubicSpline </br>
!# 
!# LinearIJtoIBJB_R2D (Interface) ------| CyclicLinear </br>
!# 
!# LinearIBJBtoIJ_R2D (Interface) ------| CyclicLinear </br>
!# 
!# NearestIJtoIBJB_I2D (Interface)------| CyclicNearest_i </br>
!# 
!# NearestIJtoIBJB_R2D (Interface)------| CyclicNearest_r </br>
!# 
!# NearestIBJBtoIJ_I2D (Interface)------| CyclicNearest_i </br>
!# 
!# NearestIBJBtoIJ_R2D (Interface)------| CyclicNearest_r </br>
!# 
!# FreqBoxIJtoIBJB_I2D (Interface)------| CyclicFreqBox_i </br>
!# 
!# FreqBoxIJtoIBJB_R2D (Interface)------| CyclicFreqBox_r </br>
!# 
!# SeaMaskIJtoIBJB_R2D (Interface)------| CyclicSeaMask_r </br>
!# 
!# SeaMaskIBJBtoIJ_R2D (Interface)------| CyclicSeaMask_r </br>
!# 
!# AveBoxIJtoIBJB_R2D (Interface) ------| CyclicAveBox_r </br>
!# 
!# AveBoxIBJBtoIJ_R2D (Interface) ------| CyclicAveBox_r </br>
!# 
!# vfirec     --------------------------| vfinit </br>
!# 
!#------------------------------------------------------------------- </br>
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
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>07-04-2011 - Paulo Kubota   - version: 1.18.0 </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 2.0.0 - some adaptations for modularizing Chopping </li>
!#  <li>09-10-2019 - Eduardo Khamis - version: 2.1.0 - changing for operational Chopping </li>
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


module Mod_Utils

  use Mod_Sizes, only : &
    mnMax, &
    mnMap, &
    mMax, &
    nMax, &
    nExtMax, &
    mnExtMax, &
    mymnExtMax, &
    iMax, &
    jMax, &
    kMax, &
    jMaxHalf, &
    jbMax, &
    ibMax, &
    ibPerIJ, &
    jbPerIJ, &
    iPerIJB, &
    jPerIJB, &
    ibMaxPerJB, &
    iMaxPerJ, &
    lm2m, &
    mExtMap, &
    nExtMap, &
    mymMax, &
    mymnExtMap, &
    mymExtMap, &
    mynExtMap, &
    myfirstlat, &
    myfirstlon, &
    mylastlon, &
    mylastlat, &
    mnExtMap

  use Mod_Parallelism_Group_Chopping, only : &
    myId, &
    msgOutMaster, &
    fatalError

  use Mod_InputParameters, only : pai, twomg, emrad, &
    ImaxOut, JmaxOut, KmaxInpp, KmaxOutp, &
    Rd, Cp, GEps, Gama, Grav, Rv

  use Mod_InputArrays, only : DelSigmaOut, SigLayerOut, SigInterOut, &
    DelSInp, SigLInp, SigIInp, &
    gTopoInp, gPsfcInp, gTvirInp, gPresInp, &
    gTopoOut, gPsfcOut, gLnPsOut

  implicit none

  private
  include 'precision.h'
  include 'messages.h'

  !  LEGANDRE POLINOMIAL AND ITS ROOTS
  !
  !  Module exports four routines:
  !     CreateLegPol  initializes module;
  !     DestroyLegPol destroys module;
  !     LegPol        computes polinomial
  !     LegPolRoots   computes roots of even degree Legandre Pol
  !
  !  Module does not export (or require) any data value.

  public :: NewSigma
  public :: SigmaInp
  public :: NewPs
  public :: CreateLegPol
  public :: DestroyLegPol
  public :: LegPol
  public :: LegPolRootsandWeights
  public :: CreateGridValues
  public :: CreateGaussQuad
  public :: DestroyGaussQuad
  public :: GaussColat
  public :: SinGaussColat
  public :: CosGaussColat
  public :: AuxGaussColat
  public :: GaussPoints
  public :: GaussWeights
  public :: colrad
  public :: colrad2D
  public :: cos2lat
  public :: ercossin
  public :: fcor
  public :: cosiv
  public :: allpolynomials
  public :: coslatj
  public :: sinlatj
  public :: coslat
  public :: sinlat
  public :: coslon
  public :: sinlon
  public :: longit
  public :: rcl
  public :: rcs2
  public :: lonrad
  public :: lati
  public :: long
  public :: cosz
  public :: cos2d
  public :: NoBankConflict
  public :: CreateAssocLegFunc
  public :: DestroyAssocLegFunc
  public :: Reset_Epslon_To_Local
  public :: Epslon
  public :: LegFuncS2F
  public :: IBJBtoIJ
  public :: IJtoIBJB
  public :: SplineIJtoIBJB
  public :: SplineIBJBtoIJ
  public :: LinearIJtoIBJB
  public :: LinearIBJBtoIJ
  public :: NearestIBJBtoIJ
  public :: NearestIJtoIBJB
  public :: FreqBoxIJtoIBJB
  public :: SeaMaskIBJBtoIJ
  public :: SeaMaskIJtoIBJB
  public :: AveBoxIBJBtoIJ
  public :: AveBoxIJtoIBJB
  public :: CyclicNearest_r
  public :: CyclicLinear
  public :: CyclicLinear_ABS
  public :: CyclicLinear_inter
  public :: Clear_Utils
  
  interface CyclicLinear_inter
    module procedure CyclicLinear_r4, CyclicLinear_r8
  end interface

  interface IBJBtoIJ
    module procedure IBJBtoIJ_R, IBJBtoIJ_I
  end interface

  interface IJtoIBJB
    module procedure &
      IJtoIBJB_R, IJtoIBJB_I, &
      IJtoIBJB3_R, IJtoIBJB3_I
  end interface

  interface SplineIBJBtoIJ
    module procedure SplineIBJBtoIJ_R2D
  end interface

  interface SplineIJtoIBJB
    module procedure SplineIJtoIBJB_R2D
  end interface

  interface LinearIBJBtoIJ
    module procedure LinearIBJBtoIJ_R2D
  end interface

  interface LinearIJtoIBJB
    module procedure LinearIJtoIBJB_R2D
  end interface

  interface NearestIBJBtoIJ
    module procedure NearestIBJBtoIJ_I2D, NearestIBJBtoIJ_R2D
  end interface

  interface NearestIJtoIBJB
    module procedure &
      NearestIJtoIBJB_I2D, NearestIJtoIBJB_R2D, &
      NearestIJtoIBJB_I3D, NearestIJtoIBJB_R3D
  end interface

  interface SeaMaskIBJBtoIJ
    module procedure  SeaMaskIBJBtoIJ_R2D
  end interface

  interface SeaMaskIJtoIBJB
    module procedure SeaMaskIJtoIBJB_R2D
  end interface

  interface FreqBoxIJtoIBJB
    module procedure FreqBoxIJtoIBJB_I2D, FreqBoxIJtoIBJB_R2D
  end interface

  interface AveBoxIBJBtoIJ
    module procedure AveBoxIBJBtoIJ_R2D
  end interface

  interface AveBoxIJtoIBJB
    module procedure AveBoxIJtoIBJB_R2D
  end interface

  interface NoBankConflict
    module procedure NoBankConflictS, NoBankConflictV
  end interface

  !  Module usage:
  !     CreateGaussQuad  should be invoked once, before any other routine
  !                      of this module, to set up maximum degree
  !                      of base functions (say, n);
  !                      it computes and hides n Gaussian Points and Weights
  !                      over interval [-1:1];
  !                      it also creates and uses Mod LegPol
  !     DestroyGaussQuad destrois hidden data structure and leaves module ready
  !                      for re-start, if desired, with another maximum degree.

  real(kind = p_r8), allocatable :: GaussColat(:)
  real(kind = p_r8), allocatable :: SinGaussColat(:)
  real(kind = p_r8), allocatable :: CosGaussColat(:)
  real(kind = p_r8), allocatable :: AuxGaussColat(:)
  real(kind = p_r8), allocatable :: GaussPoints(:)
  real(kind = p_r8), allocatable :: GaussWeights(:)
  real(kind = p_r8), allocatable :: auxpol(:, :)
  real(kind = p_r8), allocatable :: colrad(:)
  real(kind = p_r8), allocatable :: rcs2(:)
  real(kind = p_r8), allocatable :: colrad2D(:, :)
  real(kind = p_r8), allocatable :: cos2lat(:, :)
  real(kind = p_r8), allocatable :: ercossin(:, :)
  real(kind = p_r8), allocatable :: fcor(:, :)
  real(kind = p_r8), allocatable :: cosiv(:, :)
  real(kind = p_r8), allocatable :: coslatj(:)
  real(kind = p_r8), allocatable :: sinlatj(:)
  real(kind = p_r8), allocatable :: coslat(:, :)
  real(kind = p_r8), allocatable :: sinlat(:, :)
  real(kind = p_r8), allocatable :: coslon(:, :)
  real(kind = p_r8), allocatable :: sinlon(:, :)
  real(kind = p_r8), allocatable :: longit(:, :)
  real(kind = p_r8), allocatable :: rcl(:, :)
  real(kind = p_r8), allocatable :: lonrad(:, :)
  real(kind = p_r8), allocatable :: lati(:)
  real(kind = p_r8), allocatable :: long(:)
  real(kind = p_r8), allocatable :: cosz(:)
  real(kind = p_r8), allocatable :: cos2d(:, :)

  !  Module Hided Data:
  !     maxDegree is the degree of the base functions (n)
  !     created specifies if module was created or not

  logical :: created = .false.
  logical :: reducedgrid = .false.
  logical :: allpolynomials
  integer :: maxDegree = -1

  real(kind = p_r8), allocatable :: Epslon(:)
  real(kind = p_r8), allocatable :: LegFuncS2F(:, :)
  real(kind = p_r8), allocatable :: Square(:)
  real(kind = p_r8), allocatable :: Den(:)


  !  Module Hidden data

  integer :: nAuxPoly
  real(kind = p_r16), allocatable, dimension(:) :: AuxPoly1, AuxPoly2
  logical, parameter :: dumpLocal = .false.

  ! Index mappings to/from diagonal from/to column for
  ! 'extended' spectral representations

  integer, allocatable :: ExtDiagPerCol(:)   ! diag=DiagPerCol(col )
  integer, allocatable :: ExtColPerDiag(:)   ! col =ColPerDiag(diag)

  ! date are always in the form yyyymmddhh ( year, month, day, hour). hour in
  ! is in 0-24 form

  integer, private :: JulianDayInitIntegration

contains


  subroutine NewSigma(kMaxOut)
    !# Computes Delta Sigma Values 
    !# ---
    !# @info
    !# **Brief:** Sigma Interface and Layer Values.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    implicit none
    integer, intent(in) :: kMaxOut

    integer :: k

    real (kind = p_r8) :: RdoCp, CpoRd, RdoCp1

    ! Compute New Sigma Interface Values
    SigInterOut(1) = 1.0_p_r8
    do k = 1, KmaxOut
      SigInterOut(k + 1) = SigInterOut(k) - DelSigmaOut(k)
    end do
    SigInterOut(KmaxOutp) = 0.0_p_r8

    ! Compute New Sigma Layer Values
    RdoCp = Rd / Cp
    CpoRd = Cp / Rd
    RdoCp1 = RdoCp + 1.0_p_r8
    do k = 1, KmaxOut
      SigLayerOut(k) = ((SigInterOut(k)**RdoCp1 - SigInterOut(k + 1)**RdoCp1) / &
        (RdoCp1 * (SigInterOut(k) - SigInterOut(k + 1))))**CpoRd
    end do

  end subroutine NewSigma


  subroutine SigmaInp(kMaxInp)
    !# Computes Delta Sigma Values 
    !# ---
    !# @info
    !# **Brief:** Sigma Interface and Layer Values.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    implicit none
    integer, intent(in) :: kMaxInp

    ! Given Delta Sigma Values Computes
    ! Sigma Interface and Layer Values

    real (kind = p_r8) :: RdoCp, CpoRd, RdoCp1

    integer :: k

    ! Compute New Sigma Interface Values
    SigIInp(1) = 1.0_p_r8
    do k = 1, KmaxInp
      SigIInp(k + 1) = SigIInp(k) - DelSInp(k)
    end do
    SigIInp(KmaxInpp) = 0.0_p_r8

    ! Compute New Sigma Layer Values
    RdoCp = Rd / Cp
    CpoRd = Cp / Rd
    RdoCp1 = RdoCp + 1.0_p_r8
    do k = 1, KmaxInp
      SigLInp(k) = ((SigIInp(k)**RdoCp1 - SigIInp(k + 1)**RdoCp1) / &
        (RdoCp1 * (SigIInp(k) - SigIInp(k + 1))))**CpoRd
    end do

  end subroutine SigmaInp


  subroutine NewPs(kMaxInp)
    !# Computes a New Surface Pressure given a New Orography
    !# ---
    !# @info
    !# **Brief:** The New Pressure is computed assuming a Hydrostatic Balance
    !# and a Constant Temperature Lapse Rate. Below Ground, the Lapse Rate is
    !# assumed to be -6.5 K/km.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    implicit none
    integer, intent(in) :: kMaxInp

    integer :: ls, k, i, j

    real (kind = p_r8) :: goRd, GammaLR, m1, m2

    real (kind = p_r8), dimension (Ibmax, Jbmax) :: Zu

    goRd = Grav / Rd

    !write(p_nfprt,*) myid, 'MINVAL(gTopoInp) ', MINVAL(gTopoInp), 'MAXVAL(gTopoInp)= ',MAXVAL(gTopoInp)
    !write(p_nfprt,*) myid, 'MINVAL(gTvirInp) ', MINVAL(gTvirInp), 'MAXVAL(gTvirInp)= ',MAXVAL(gTvirInp)
    !write(p_nfprt,*) myid, 'MINVAL(gPsfcInp) ', MINVAL(gPsfcInp), 'MAXVAL(gPsfcInp)= ',MAXVAL(gPsfcInp)
    !write(p_nfprt,*) myid, 'MINVAL(gPresInp) ', MINVAL(gPresInp), 'MAXVAL(gPresInp)= ',MAXVAL(gPresInp)
    !write(p_nfprt,*) myid, 'MINVAL(gTopoOut) ', MINVAL(gTopoOut), 'MAXVAL(gTopoOut)= ',MAXVAL(gTopoOut)

    ! Compute Surface Pressure Below the Original Ground
    ls = 0
    k = 1
    GammaLR = Gama
    do j = 1, Jbmax
      do i = 1, Ibmaxperjb(j)
        Zu(i, j) = gTopoInp(i, j) - gTvirInp(i, k, j) / GammaLR * &
          ((gPsfcInp(i, j) / gPresInp(i, k, j))**(-GammaLR / goRd) - 1.0_p_r8)
        if (gTopoOut(i, j) <= Zu(i, j)) then
          if (abs(GammaLR) > GEps) then
            gPsfcOut(i, j) = gPresInp(i, k, j) * (1.0_p_r8 + GammaLR / gTvirInp(i, k, j) * &
              (gTopoOut(i, j) - Zu(i, j)))**(-goRd / GammaLR)
          else
            gPsfcOut(i, j) = gPresInp(i, k, j) * exp(-goRd / gTvirInp(i, k, j) * &
              (gTopoOut(i, j) - Zu(i, j)))
          end if
        else
          gPsfcOut(i, j) = 0.0_p_r8
          ls = ls + 1
        end if
      end do
    end do

    ! Compute Surface Pressure Above the Original Ground
    do k = 2, KmaxInp
      if (ls > 0) then
        do j = 1, Jbmax
          do i = 1, Ibmaxperjb(j)
            if (gPsfcOut(i, j) == 0.0_p_r8) then
              GammaLR = -goRd * log(gTvirInp(i, k - 1, j) / gTvirInp(i, k, j)) / &
                log(gPresInp(i, k - 1, j) / gPresInp(i, k, j))
              if (abs(GammaLR) > GEps) then
                Zu(i, j) = Zu(i, j) - gTvirInp(i, k, j) / GammaLR * &
                  ((gPresInp(i, k - 1, j) / gPresInp(i, k, j))** &
                    (-GammaLR / goRd) - 1.0_p_r8)
              else
                Zu(i, j) = Zu(i, j) + gTvirInp(i, k, j) / &
                  goRd * log(gPresInp(i, k - 1, j) / gPresInp(i, k, j))
              end if
              if (gTopoOut(i, j) <= Zu(i, j)) then
                if (abs(GammaLR) > GEps) then
                  gPsfcOut(i, j) = gPresInp(i, k, j) * (1.0_p_r8 + GammaLR / gTvirInp(i, k, j) * &
                    (gTopoOut(i, j) - Zu(i, j)))**(-goRd / GammaLR)
                else
                  gPsfcOut(i, j) = gPresInp(i, k, j) * &
                    exp(-goRd / gTvirInp(i, k, j) * (gTopoOut(i, j) - Zu(i, j)))
                end if
                ls = ls - 1
              end if
            end if
          end do
        end do
      end if
    end do

    ! Compute Surface Pressure Over the Top
    if (ls > 0) then
      k = KmaxInp
      GammaLR = 0.0_p_r8
      do j = 1, Jbmax
        do i = 1, Ibmaxperjb(j)
          if (gPsfcOut(i, j) == 0.0_p_r8) then
            gPsfcOut(i, j) = gPresInp(i, k, j) * &
              exp(-goRd / gTvirInp(i, k, j) * (gTopoOut(i, j) - Zu(i, j)))
          end if
        end do
      end do
    end if

    ! ln(ps) in cBar (from mBar)
    m1 = 1.e7_p_r8
    m2 = 0.
    do j = 1, Jbmax
      do i = 1, Ibmaxperjb(j)
        m1 = min(m1, gPsfcOut(i, j))
        m2 = max(m2, gPsfcOut(i, j))
        gLnPsOut(i, j) = log(gPsfcOut(i, j) / 10.0_p_r8)
      end do
    end do
    !write(p_nfprt,*) myid, 'minnewps ', m1, 'maxnewps ',m2
  end subroutine NewPs


  subroutine CreateAssocLegFunc(allpolynomials)
    !# Initializes module and compute functions
    !# ---
    !# @info
    !# **Brief:** Should be invoked once, before any other routine, to compute
    !# and store Associated Legendre Functions at Gaussian Points for all
    !# Legendre Orders and Degrees (defined at Sizes). </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin  
    logical, intent(in) :: allpolynomials
    integer :: m, n, mn, j, mp, np, jp, mglobalp
    character(len = *), parameter :: h = "**(CreateAssocLegFunc)**"

    if (allpolynomials) then
      allocate (Epslon(mnExtMax))
      allocate (LegFuncS2F(jMaxHalf, mnExtMax))
    else
      allocate (Epslon(mymnExtMax))
      allocate (LegFuncS2F(jMaxHalf, mymnExtMax))
      allocate (auxpol(jMaxHalf, mMax))
    endif

    allocate (Square(nExtMax))
    allocate (Den(nExtMax))
    do n = 1, nExtMax
      Square(n) = real((n - 1) * (n - 1), p_r8)
      Den(n) = 1.0_p_r8 / (4.0_p_r8 * Square(n) - 1.0_p_r8)
    end do
    if (allpolynomials) then
      do mn = 1, mnExtMax
        m = mExtMap(mn)
        n = nExtMap(mn)
        Epslon(mn) = sqrt((Square(n) - Square(m)) * Den(n))
      end do
    else
      do mn = 1, mymnExtMax
        m = lm2m(mymExtMap(mn))
        n = mynExtMap(mn)
        Epslon(mn) = sqrt((Square(n) - Square(m)) * Den(n))
      end do
    end if

    if (allpolynomials) then
      LegFuncS2F(1:jMaxHalf, mnExtMap(1, 1)) = sqrt(0.5_p_r8)
      do m = 2, mMax
        do j = 1, jMaxHalf
          LegFuncS2F(j, mnExtMap(m, m)) = &
            sqrt(1.0_p_r8 + 0.5_p_r8 / real(m - 1, p_r8)) * &
              SinGaussColat(j) * &
              LegFuncS2F(j, mnExtMap(m - 1, m - 1))
        end do
      end do
      !$OMP PARALLEL DO PRIVATE(mp,jp)
      do mp = 1, mMax
        do jp = 1, jMaxHalf
          LegFuncS2F(jp, mnExtMap(mp, mp + 1)) = &
            sqrt(1.0_p_r8 + 2.0_p_r8 * real(mp, p_r8)) * &
              GaussPoints(jp) * &
              LegFuncS2F(jp, mnExtMap(mp, mp))
        end do
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO PRIVATE(mp,jp,np)
      do mp = 1, mMax
        do np = mp + 2, nExtMax
          do jp = 1, jMaxHalf
            LegFuncS2F(jp, mnExtMap(mp, np)) = &
              (GaussPoints(jp) * &
                LegFuncS2F(jp, mnExtMap(mp, np - 1)) - &
                Epslon(mnExtMap(mp, np - 1)) * &
                  LegFuncS2F(jp, mnExtMap(mp, np - 2))   &
                ) / Epslon(mnExtMap(mp, np))
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      call Reset_Epslon_To_Local ()
    else
      auxpol(1:jMaxHalf, 1) = sqrt(0.5_p_r8)
      do m = 2, mMax
        do j = 1, jMaxHalf
          auxpol(j, m) = &
            sqrt(1.0_p_r8 + 0.5_p_r8 / real(m - 1, p_r8)) * &
              SinGaussColat(j) * &
              auxpol(j, m - 1)
        end do
      end do
      !$OMP PARALLEL DO PRIVATE(mp,jp,mglobalp)
      do mp = 1, mymMax
        mglobalp = lm2m(mp)
        do jp = 1, jMaxHalf
          LegFuncS2F(jp, mymnExtMap(mp, mglobalp)) = &
            auxpol(jp, mglobalp)
        end do
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO PRIVATE(mp,jp,mglobalp)
      do mp = 1, mymMax
        mglobalp = lm2m(mp)
        do jp = 1, jMaxHalf
          LegFuncS2F(jp, mymnExtMap(mp, mglobalp + 1)) = &
            sqrt(1.0_p_r8 + 2.0_p_r8 * real(mglobalp, p_r8)) * &
              GaussPoints(jp) * &
              LegFuncS2F(jp, mymnExtMap(mp, mglobalp))
        end do
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO PRIVATE(mp,jp,mglobalp,np)
      do mp = 1, mymMax
        mglobalp = lm2m(mp)
        do np = mglobalp + 2, nExtMax
          do jp = 1, jMaxHalf
            LegFuncS2F(jp, mymnExtMap(mp, np)) = &
              (GaussPoints(jp) * &
                LegFuncS2F(jp, mymnExtMap(mp, np - 1)) - &
                Epslon(mymnExtMap(mp, np - 1)) * &
                  LegFuncS2F(jp, mymnExtMap(mp, np - 2))   &
                ) / Epslon(mymnExtMap(mp, np))
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine CreateAssocLegFunc


  subroutine Reset_Epslon_To_Local ()
    !# Resets Epslon To Local
    !# ---
    !# @info
    !# **Brief:** Converts Epslon to local mpi structure. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin  
    integer :: m, n, mn

    do mn = 1, mymnExtMax
      m = lm2m(mymExtMap(mn))
      n = mynExtMap(mn)
      Epslon(mn) = sqrt((Square(n) - Square(m)) * Den(n))
    end do

  end subroutine Reset_Epslon_To_Local


  subroutine DestroyAssocLegFunc()
    !# Destroys all internal info
    !# ---
    !# @info
    !# **Brief:** Deallocates all stored values. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin 
    character(len = *), parameter :: h = "**(DestroyAssocLegFunc)**"
    deallocate (Epslon, LegFuncS2F)
    deallocate (Square, Den)
  end subroutine DestroyAssocLegFunc


  !  AUXILIARY PROCEDURES
  !  Module exports two routines:
  !     NoBankConflict  

  function NoBankConflictS(s) result(p)
    !# Avoids memory bank conflicts
    !# ---
    !# @info
    !# **Brief:** Given an input integer (size of an array in any dimension)
    !# returns the next integer that should dimension the array to avoid memory
    !# bank conflicts. The vector version returns a vector of integers, given a
    !# vector of input integers. Module does not require any other module and
    !# does not export any value.</br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin 
    integer, intent(in) :: s
    integer :: p
    if ((mod(s, 2)==0) .and. (s/=0)) then
      p = s + 1
    else
      p = s
    end if
  end function NoBankConflictS
 
 
  function NoBankConflictV(s) result(p)
    !# Avoids memory bank conflicts
    !# ---
    !# @info
    !# **Brief:** Given an input integer (size of an array in any dimension)
    !# returns the next integer that should dimension the array to avoid memory
    !# bank conflicts. The vector version returns a vector of integers, given a
    !# vector of input integers. Module does not require any other module and
    !# does not export any value.</br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin 
    integer, intent(in) :: s(:)
    integer :: p(size(s))
    where ((mod(s, 2)==0) .and. (s/=0))
      p = s + 1
    elsewhere
      p = s
    end where
  end function NoBankConflictV

  !  GAUSSIAN POINTS AND WEIGHTS FOR QUADRATURE
  !  OVER LEGENDRE POLINOMIALS BASE FUNCTIONS

  !  Module exports two routines:
  !     CreateGaussQuad   initializes module
  !     DestroyGaussQuad  destroys module

  !  Module export two arrays:
  !     GaussPoints and GaussWeights
  !
  !  Module uses Module LegPol (Legandre Polinomials);
  !  transparently to the user, that does not have
  !  to create and/or destroy LegPol

  subroutine CreateGaussQuad (degreeGiven)
    !# Initializes module
    !# ---
    !# @info
    !# **Brief:** Computes and hides 'degreeGiven' Gaussian Points and Weights
    !# over interval [-1:1]; Creates and uses Mod LegPol (Legandre Polinomials).</br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin 
    integer, intent(in) :: degreeGiven
    real(kind = p_r16), allocatable :: FVals(:)
    character(len = *), parameter :: h = "**(CreateGaussQuad)**"
    character(len = 256) :: line
    character(len = 10) :: c1
    integer :: j

    !check invocation sequence and input data

    if (degreeGiven <=0) then
      write(c1, "(i10)") degreeGiven
      write(p_nfprt, "(a,' invoked with degree ',a)") h, trim(adjustl(c1))
      stop
    else
      maxDegree = degreeGiven
    end if

    !allocate areas
    if (allocated(GaussColat))deallocate(GaussColat)
    if (allocated(SinGaussColat))deallocate(SinGaussColat)
    if (allocated(CosGaussColat))deallocate(CosGaussColat)
    if (allocated(AuxGaussColat))deallocate(AuxGaussColat)
    if (allocated(FVals))deallocate(FVals)
    if (allocated(GaussPoints))deallocate(GaussPoints)
    if (allocated(GaussWeights))deallocate(GaussWeights)
    if (allocated(colrad))deallocate(colrad)
    if (allocated(rcs2))deallocate(rcs2)

    allocate (GaussColat(maxDegree / 2))
    allocate (SinGaussColat(maxDegree / 2))
    allocate (CosGaussColat(maxDegree / 2))
    allocate (AuxGaussColat(maxDegree / 2))
    allocate (FVals(maxDegree / 2))
    allocate (GaussPoints(maxDegree))
    allocate (GaussWeights(maxDegree))
    allocate (colrad(maxDegree))
    allocate (rcs2(maxDegree / 2))

    !create ModLegPol

    call CreateLegPol (maxDegree)

    !Gaussian Points are the roots of legandre polinomial of degree maxDegree

    call LegPolRootsandWeights(maxDegree)
    GaussWeights(maxDegree / 2 + 1:maxDegree) = GaussWeights(maxDegree / 2:1:-1)
    GaussColat = acos(CosGaussColat)
    SinGaussColat = sin(GaussColat)
    AuxGaussColat = 1.0_p_r16 / (SinGaussColat * SinGaussColat)
    GaussPoints(1:maxDegree / 2) = CosGaussColat
    GaussPoints(maxDegree / 2 + 1:maxDegree) = -GaussPoints(maxDegree / 2:1:-1)

    do j = 1, maxDegree / 2
      colrad(j) = GaussColat(j)
      colrad(maxDegree + 1 - j) = pai - GaussColat(j)
      rcs2(j) = AuxGaussColat(j)
    end do

  end subroutine CreateGaussQuad


  subroutine CreateGridValues
    !# Creates Grid Values
    !# ---
    !# @info
    !# **Brief:** Creates Grid Values.</br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin 
    integer :: i, j, ib, jb, jhalf
    real (kind = p_r8) :: rd
    allocate (colrad2D(ibMax, jbMax))
    allocate (rcl     (ibMax, jbMax))
    rcl = 0.0_p_r8
    allocate (coslatj (jMax))
    allocate (sinlatj (jMax))
    allocate (coslat  (ibMax, jbMax))
    allocate (sinlat  (ibMax, jbMax))
    allocate (coslon  (ibMax, jbMax))
    allocate (sinlon  (ibMax, jbMax))
    allocate (cos2lat (ibMax, jbMax))
    allocate (ercossin(ibMax, jbMax))
    allocate (fcor    (ibMax, jbMax))
    allocate (cosiv   (ibMax, jbMax))
    allocate (longit  (ibMax, jbMax))
    allocate (lonrad  (ibMax, jbMax))
    allocate (lati    (jMax))
    allocate (long    (iMax))
    allocate (cosz    (jMax))
    allocate (cos2d   (ibMax, jbMax))

    rd = 45.0_p_r8 / atan(1.0_p_r8)
    do j = 1, jMax
      lati(j) = 90.0_p_r8 - colrad(j) * rd
    end do
    do i = 1, imax
      long(i) = (i - 1) * 360.0_p_r8 / real(iMax, p_r8)
    enddo

    !$OMP PARALLEL DO PRIVATE(jhalf)
    do j = 1, jMax / 2
      jhalf = jMax - j + 1
      sinlatj (j) = cos(colrad(j))
      coslatj (j) = sin(colrad(j))
      sinlatj (jhalf) = - sinlatj(j)
      coslatj (jhalf) = coslatj(j)
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ib,i,j,jhalf)
    do jb = 1, jbMax
      do ib = 1, ibMaxPerJB(jb)
        j = jPerIJB(ib, jb)
        i = iPerIJB(ib, jb)
        jhalf = min(j, jMax - j + 1)
        colrad2D(ib, jb) = colrad(j)
        sinlat  (ib, jb) = sinlatj(j)
        coslat  (ib, jb) = coslatj(j)
        rcl     (ib, jb) = rcs2(jhalf)
        longit  (ib, jb) = (i - 1) * pai * 2.0_p_r8 / iMaxPerJ(j)
        lonrad  (ib, jb) = (i - 1) * 360.0_p_r8 / real(iMaxPerJ(j), p_r8)
        sinlon  (ib, jb) = sin(longit(ib, jb))
        coslon  (ib, jb) = cos(longit(ib, jb))
        cos2lat (ib, jb) = 1.0_p_r8 / rcl(ib, jb)
        !**(JP)** mudei
        !          ercossin(ib,jb)=sinlat(ib,jb)*rcl(ib,jb)/er
        !          fcor    (ib,jb)=twomg*sinlat(ib,jb)
        !          cosiv   (ib,jb)=1./coslat(ib,jb)
        ! para a forma antiga, so para bater binario no euleriano;
        ! impacto no SemiLagrangeano eh desconhecido
        ercossin(ib, jb) = cos(colrad(j)) * rcl(ib, jb) / emrad
        fcor    (ib, jb) = twomg * cos(colrad(j))
        cosiv   (ib, jb) = sqrt(rcl(ib, jb))
        !**(JP)** fim de alteracao
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine CreateGridValues


  subroutine DestroyGaussQuad
    !# Destroys module
    !# ---
    !# @info
    !# **Brief:** 
    !# <ul type="disc">
    !#  <li>Destroy module internal data; </li>
    !#  <li>Destroys module LegPol; </li>
    !#  <li>Get ready for new module usage. </li>
    !# </ul>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin 
    character(len = *), parameter :: h = "**(DestroyGaussQuad)**"

    !check invocation sequence

    maxDegree = -1

    !destroy ModLegPol

    call DestroyLegPol()

    !deallocate module areas
    !
    deallocate (GaussColat)
    deallocate (SinGaussColat)
    deallocate (AuxGaussColat)
    deallocate (GaussPoints)
    deallocate (GaussWeights)
  end subroutine DestroyGaussQuad


  !  Module usage:
  !     CreateLegPol should be invoked once, before any other routine;
  !     LegPol can be invoked after CreateLegPol, as much as required;
  !     LegPolRoots can be invoked after CreateLegPol, as much as required;
  !     DestroyLegPol should be invoked at the end of the computation.
  !     If this maximum degree has to be changed, module should be destroied and
  !     created again.
  !
  !  Module Hided Data:
  !     AuxPoly1, AuxPoly2 are constants used to evaluate LegPol
  !     nAuxPoly is the maximum degree of the Polinomial to be computed
  !     created specifies if module was created or not

  subroutine CreateLegPol (maxDegree)
    !# Does module initialization
    !# ---
    !# @info
    !# **Brief:**  Should be invoked once, before any other routine. Estabilishes
    !# the maximum degree for which LegPol and LegPolRoots should be invoked. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: maxDegree
    integer :: i
    !
    !  Check Correction
    !
    if (maxDegree <= 0) then
      write (0, "('**(CreateLegPol)** maxDegree <= 0; maxDegree =')") maxDegree
      stop
    else if (created) then
      write (0, "('**(CreateLegPol)** invoked twice without destruction')")
      stop
    end if
    !
    !  Allocate and Compute Constants
    !
    nAuxPoly = maxDegree
    if (allocated(AuxPoly1))deallocate(AuxPoly1)
    if (allocated(AuxPoly2))deallocate(AuxPoly2)

    allocate (AuxPoly1(nAuxPoly))
    allocate (AuxPoly2(nAuxPoly))
    do i = 1, nAuxPoly
      AuxPoly1(i) = real(2 * i - 1, p_r16) / real(i, p_r16)
      AuxPoly2(i) = real(1 - i, p_r16) / real(i, p_r16)
    end do
    created = .true.
  end subroutine CreateLegPol


  subroutine DestroyLegPol
    !# Destroys module initialization
    !# ---
    !# @info
    !# **Brief:** Should be invoked at the end of the computation. Required if
    !# maximum degree of LegPol (Legandre Polinomials) has to be changed. In this case, invoke
    !# DestroyLegPol and CreateLegPol with new degree. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    if (.not. created) then
      write (0, "('**(DestroyLegPol)** invoked without initialization')")
      stop
    end if
    deallocate (AuxPoly1, AuxPoly2)
    created = .false.
  end subroutine DestroyLegPol


  function LegPol (degree, Col)
    !# Computes LegPol (Legandre Polinomial)
    !# ---
    !# @info
    !# **Brief:**  
    !# <ul type="disc">
    !#  <li>Can be invoked after CreateLegPol, as much as required; </li>
    !#  <li>Computes Legandre Polinomial of degree 'degree' at a vector of
    !#      colatitudes 'Col'; </li>
    !#  <li>Degree should be <= maximum degree estabilished by CreateLegPol; </li>
    !#  <li>Deals with colatitudes; </li>
    !#  <li>Takes colatitudes as abcissas; </li>
    !#  <li>Colatitudes are expressed in radians, with 0 at the North Pole, pi/2 at the Equator and pi at the South Pole. </li>
    !# </ul>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: degree
    real(kind = p_r16), dimension(:), intent(in) :: Col
    real(kind = p_r16), dimension(size(Col)) :: LegPol
    !
    !  Auxiliary Variables
    !
    real(kind = p_r16), dimension(size(Col)) :: P0     ! Polinomial of degree i - 2
    real(kind = p_r16), dimension(size(Col)) :: P1     ! Polinomial of degree i - 1
    real(kind = p_r16), dimension(size(Col)) :: X      ! Cosine of colatitude
    integer :: iDegree                   ! loop index
    integer :: left                      ! loop iterations before unrolling
    !
    !  Check Correctness
    !
    if (.not. created) then
      write (0, "('**(LegPol)** invoked without initialization')")
      stop
    end if
    !
    !  Case degree >=2 and degree <= maximum degree
    !
    if ((degree >= 2) .and. (degree <= nAuxPoly)) then
      !
      !  Initialization
      !
      left = mod(degree - 1, 6)
      X = cos(Col)
      P0 = 1.0_p_r16
      P1 = X
      !
      !  Apply recurrence relation
      !     Loop upper bound should be 'degree';
      !     it is not due to unrolling
      !
      do iDegree = 2, left + 1
        LegPol = X * P1 * AuxPoly1(iDegree) + P0 * AuxPoly2(iDegree)
        P0 = P1; P1 = LegPol
      end do
      !
      !  Unroll recurrence relation, to speed up computation
      !
      do iDegree = left + 2, degree, 6
        LegPol = X * P1 * AuxPoly1(iDegree) + P0 * AuxPoly2(iDegree)
        P0 = X * LegPol * AuxPoly1(iDegree + 1) + P1 * AuxPoly2(iDegree + 1)
        P1 = X * P0 * AuxPoly1(iDegree + 2) + LegPol * AuxPoly2(iDegree + 2)
        LegPol = X * P1 * AuxPoly1(iDegree + 3) + P0 * AuxPoly2(iDegree + 3)
        P0 = X * LegPol * AuxPoly1(iDegree + 4) + P1 * AuxPoly2(iDegree + 4)
        P1 = X * P0 * AuxPoly1(iDegree + 5) + LegPol * AuxPoly2(iDegree + 5)
      end do
      LegPol = P1
      !
      !  Case degree == 0
      !
    else if (degree == 0) then
      LegPol = 1.0_p_r16
      !
      !  Case degree == 1
      !
    else if (degree == 1) then
      LegPol = cos(Col)
      !
      !  Case degree <= 0 or degree > maximum degree
      !
    else
      write(p_nfprt, "('**(LegPol)** invoked with degree ',i6,&
        &' out of bounds')") degree
      stop
    end if
  end function LegPol


  subroutine LegPolRootsandWeights (degree)
    !# Computes the roots of the Legandre Polinomial
    !# ---
    !# @info
    !# **Brief:** 
    !# <ul type="disc">
    !#  <li>Can be invoked after CreateLegPol, as much as required; </li>
    !#  <li>Computes the roots of the Legandre Polinomial of even degree 'degree'
    !#      and respective Gaussian Weights; </li>
    !#  <li>Deal with colatitudes; </li>
    !#  <li>Produces colatitudes as roots; </li>
    !#  <li>Colatitudes are expressed in radians, with 0 at the North Pole, pi/2
    !#      at the Equator and pi at the South Pole. </li>
    !# </ul>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: degree
    real(kind = p_r16), dimension(degree / 2) :: Col
    real(kind = p_r16), allocatable, dimension(:) :: XSearch, FSearch
    real(kind = p_r16), dimension(degree / 2) :: Pol    ! Polinomial of degree i
    real(kind = p_r16), dimension(degree / 2) :: P0     ! Polinomial of degree i-2
    real(kind = p_r16), dimension(degree / 2) :: P1     ! Polinomial of degree i-1
    real(kind = p_r16), dimension(degree / 2) :: X
    real(kind = p_r16), dimension(degree / 2) :: XC
    logical, allocatable, dimension(:) :: Mask
    integer :: i                               ! loop index
    integer :: nPoints                         ! to start bissection
    integer, parameter :: multSearchStart = 4 ! * factor to start bissection
    real(kind = p_r16) :: step
    real(kind = p_r16) :: pi, scale
    integer, parameter :: itmax = 10     ! maximum number of newton iterations
    integer, parameter :: nDigitsOut = 2 ! precision digits of gaussian points:
    ! machine epsilon - nDigitsOut
    real(kind = p_r16) :: rootPrecision   ! relative error in gaussian points

    integer :: halfDegree
    integer :: iDegree, it               ! loop index
    integer :: left                      ! loop iterations before unrolling
    !
    !  Check Correctness
    !
    if (.not. created) then
      write (0, "('**(LegPolRootsandWeights)** invoked without initialization')")
      stop
    else if ((degree <= 0) .or. (degree > nAuxPoly)) then
      write (0, "('**(LegPolRootsandWeights)** invoked with degree ',i6,&
        &' out of bounds')") degree
      stop
    else if (mod(degree, 2) .ne. 0) then
      write (0, "('**(LegPolRootsandWeights)** invoked with odd degree ',i6)") degree
      stop
    end if
    !
    !  Initialize Constants
    !
    pi = 4.0_p_r16 * atan(1.0_p_r16)
    rootPrecision = epsilon(1.0_p_r16) * 10.0_p_r16**(nDigitsOut)
    !
    !  LegPolRoots uses root simmetry with respect to pi/2.
    !  It finds all roots in the interval [0,pi/2]
    !  Remaining roots are simmetric
    !
    halfDegree = degree / 2
    !
    !  bissection method to find roots:
    !  get equally spaced points in interval [0,pi/2]
    !  to find intervals containing roots
    !
    nPoints = multSearchStart * halfDegree
    step = pi / (2.0_p_r16 * real(nPoints, p_r16))
    allocate (XSearch(nPoints))
    allocate (FSearch(nPoints))
    allocate (Mask   (nPoints - 1))
    do i = 1, nPoints
      XSearch(i) = step * real(i - 1, p_r16)
    end do
    FSearch = LegPol (degree, XSearch)
    !
    !  select intervals containing roots
    !
    Mask = FSearch(1:nPoints - 1) * FSearch(2:nPoints) < 0.0
    !
    !  are there enough intervals?
    !
    if (count(Mask) .ne. halfDegree) then
      write(p_nfprt, "('**(LegPolRoots)** ',i6,' bracketing intervals to find '&
        &,i6,' roots')") count(Mask), halfDegree
      stop
    end if
    !
    !  extract intervals containing roots
    !
    Col = 0.5_p_r16 * (pack(XSearch(1:nPoints - 1), Mask) + &
      pack(XSearch(2:nPoints), Mask))
    deallocate (XSearch)
    deallocate (FSearch)
    deallocate (Mask)
    !
    !    loop while there is a root to be found
    !
    it = 1
    X = cos(Col)
    scale = 2.0_p_r16 / real(degree * degree, p_r16)
    left = mod(degree - 1, 6)
    !   WRITE (UNIT=*, FMT='(/,I6,2I5,1P2E16.8)') degree, p_r16, p_r16, rootPrecision,ATAN(1.0_p_r16)/REAL(degree,p_r16)
    do
      if (it.gt.itmax) then
        write (0, "('**(LegPolRootsandWeights)** failed to converge  ',i6,&
          &' itmax')") itmax
        EXIT
      end if
      !
      !   initialization
      !
      it = it + 1
      P0 = 1.0_p_r16
      P1 = X
      !
      !   Apply recurrence relation
      !     Loop upper bound should be 'degree';
      !     it is not due to unrolling
      !
      do iDegree = 2, left + 1
        Pol = X * P1 * AuxPoly1(iDegree) + P0 * AuxPoly2(iDegree)
        P0 = P1; P1 = Pol
      end do
      !
      !   Unroll recurrence relation, to speed up computation
      !
      do iDegree = left + 2, degree, 6
        Pol = X * P1 * AuxPoly1(iDegree) + P0 * AuxPoly2(iDegree)
        P0 = X * Pol * AuxPoly1(iDegree + 1) + P1 * AuxPoly2(iDegree + 1)
        P1 = X * P0 * AuxPoly1(iDegree + 2) + Pol * AuxPoly2(iDegree + 2)
        Pol = X * P1 * AuxPoly1(iDegree + 3) + P0 * AuxPoly2(iDegree + 3)
        P0 = X * Pol * AuxPoly1(iDegree + 4) + P1 * AuxPoly2(iDegree + 4)
        P1 = X * P0 * AuxPoly1(iDegree + 5) + Pol * AuxPoly2(iDegree + 5)
      end do
      XC = P1 * (1.0_p_r16 - X * X) / (real(degree, p_r16) * (P0 - X * P1))
      X = X - XC
      if (maxval(abs(XC / X)).lt.rootPrecision) then
        GaussWeights(1:halfDegree) = real(scale * (1.0_p_r16 - X * X) / (P0 * P0), kind = p_r8)
        GaussWeights(halfDegree + 1:Degree) = GaussWeights(halfDegree:1:-1)
        EXIT
      end if
    end do
    !
    CosGaussColat = real(X, kind = p_r8)
    !
  end subroutine LegPolRootsandWeights


  subroutine IBJBtoIJ_R(var_in, var_out)
    !# Maps (ib,jb) into (i,j)
    !# ---
    !# @info
    !# **Brief:** Maps (ib,jb) into (i,j) </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: var_in (:, :)
    real(kind = p_r8), intent(OUT) :: var_out(:, :)
    integer :: i
    integer :: j
    integer :: ib
    integer :: jb
    !$OMP PARALLEL DO PRIVATE(i,j,ib)
    do jb = 1, jbmax
      do ib = 1, ibmaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb) - myfirstlat + 1
        var_out(i, j) = var_in(ib, jb)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine IBJBtoIJ_R


  subroutine IJtoIBJB_R(var_in, var_out)
    !# Maps (i,j) into (ib,jb)
    !# ---
    !# @info
    !# **Brief:** Maps (i,j) into (ib,jb) </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: var_in (iMax, jMax)
    real(kind = p_r8), intent(OUT) :: var_out(ibMax, jbMax)
    integer :: i
    integer :: j
    integer :: ib
    integer :: jb
    !$OMP PARALLEL DO PRIVATE(i,j,ib)
    do jb = 1, jbMax
      do ib = 1, ibMaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb)
        var_out(ib, jb) = var_in(i, j)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine IJtoIBJB_R


  subroutine IJtoIBJB3_R(var_in, var_out)
    !# Maps IJ into IBJB3_R
    !# ---
    !# @info
    !# **Brief:** 3D version by hmjb </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: var_in (iMax, kMax, jMax)
    real(kind = p_r8), intent(OUT) :: var_out(ibMax, kMax, jbMax)
    integer :: i
    integer :: j
    integer :: ib
    integer :: jb
    !$OMP PARALLEL DO PRIVATE(i,j,ib)
    do jb = 1, jbMax
      do ib = 1, ibMaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb)
        var_out(ib, :, jb) = var_in(i, :, j)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine IJtoIBJB3_R


  subroutine IBJBtoIJ_I(var_in, var_out)
    !# Maps (ib,jb) into (i,j)
    !# ---
    !# @info
    !# **Brief:** Maps (ib,jb) into (i,j) </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: var_in (:, :)
    integer, intent(OUT) :: var_out(:, :)
    integer :: i
    integer :: j
    integer :: ib
    integer :: jb
    !$OMP PARALLEL DO PRIVATE(i,j,ib)
    do jb = 1, jbmax
      do ib = 1, ibmaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb) - myfirstlat + 1
        var_out(i, j) = var_in(ib, jb)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine IBJBtoIJ_I


  subroutine IJtoIBJB_I(var_in, var_out)
    !# Maps (i,j) into (ib,jb)
    !# ---
    !# @info
    !# **Brief:** Maps (i,j) into (ib,jb) </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: var_in (:, :)
    integer, intent(OUT) :: var_out(:, :)
    integer :: i
    integer :: j
    integer :: ib
    integer :: jb
    !$OMP PARALLEL DO PRIVATE(i,j,ib)
    do jb = 1, jbMax
      do ib = 1, ibMaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb)
        var_out(ib, jb) = var_in(i, j)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine IJtoIBJB_I
  
  
  subroutine IJtoIBJB3_I(var_in, var_out)
    !# 3D version by hmjb
    !# ---
    !# @info
    !# **Brief:** 3D version by hmjb </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: var_in (iMax, kMax, jMax)
    integer, intent(OUT) :: var_out(ibMax, kMax, jbMax)
    integer :: i
    integer :: j
    integer :: ib
    integer :: jb
    !$OMP PARALLEL DO PRIVATE(i,j,ib)
    do jb = 1, jbMax
      do ib = 1, ibMaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb)
        var_out(ib, :, jb) = var_in(i, :, j)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine IJtoIBJB3_I


  subroutine SplineIJtoIBJB_R2D(FieldIn, FieldOut)
    !# Spline IJ into IBJB_R2D
    !# ---
    !# @info
    !# **Brief:** Spline IJ into IBJB_R2D. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (iMax, jMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(ibMax, jbMax)

    real(kind = p_r8) :: FOut(iMax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**SplineIJtoIBJB**"

    print *, h

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,ib,jb,i,fout)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicCubicSpline(iMax, iMaxPerJ(j), &
          FieldIn(1, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine SplineIJtoIBJB_R2D


  subroutine SplineIBJBtoIJ_R2D(FieldIn, FieldOut)
    !# Spline IBJB into IJ_R2D
    !# ---
    !# @info
    !# **Brief:** Spline IBJB into IJ_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (ibMax, jbMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(iMax, jMax)

    integer :: i
    integer :: iFirst
    integer :: ib
    integer :: j
    integer :: jl
    integer :: jb

    character(len = *), parameter :: h = "**SplineIBJBtoIJ**"

    print *, h
    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
      do j = myfirstlat, mylastlat
        jl = j - myfirstlat + 1
        iFirst = ibPerIJ(1, j)
        jb = jbPerIJ(1, j)
        !         CALL CyclicCubicSpline(iMaxPerJ(j), iMax, &
        !              FieldIn(iFirst,jb), FieldOut(1,jl))
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb) - myfirstlat + 1
          FieldOut(i, j) = FieldIn(ib, jb)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine SplineIBJBtoIJ_R2D


  subroutine CyclicCubicSpline (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Cubic Spline interpolation on cyclic
    !# ---
    !# @info
    !# **Brief:** Cubic Spline interpolation on cyclic, equally spaced data over [0:2pi].
    !# <ul type="disc">
    !#  <li>Input data: <ul type="disc"><li>DimIn: How many input points: abcissae
    !#  are supposed to be x(i), i=1,...,DimIn, with x(i) = 2*pi*(i-1)/DimIn. </li>
    !#                                  <li>DimOut: How many output points: abcissae
    !#  are supposed to be x(i), i=1,...,DimOut, with x(i) = 2*pi*(i-1)/DimOut. </li>
    !#                                  <li>FieldIn: function values at the DimIn abcissae. </li></ul></li>
    !#  <li>Output data: <ul type="disc"> FieldOut: function values at the DimOut abcissae </ul></li>
    !#  <li>Requirements: <ul type="disc"> DimIn >= 4. Program Stops if DimIn < 4. </ul></li>
    !# </ul>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real(kind = p_r8), intent(in) :: FieldIn(DimIn)
    real(kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: alfa(DimIn + 2)
    real(kind = p_r8) :: beta(DimIn + 2)
    real(kind = p_r8) :: gama(DimIn + 2)
    real(kind = p_r8) :: delta(DimIn + 2)
    real(kind = p_r8) :: der(DimIn + 2)
    real(kind = p_r8) :: c1(DimIn + 2)
    real(kind = p_r8) :: c2(DimIn + 2)
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: dxm
    real(kind = p_r8) :: dxm2
    real(kind = p_r8) :: dx
    real(kind = p_r8) :: dx2
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hIn2
    real(kind = p_r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 4

    if (DimIn < 4) then
      stop "**(CyclicCubicSpline)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 2) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)
    hIn2 = hIn * hIn

    ! tridiagonal system initialization:
    ! indices of alfa, beta, gama and delta are equation numbers,
    ! which are derivatives on a diagonal system

    alfa(2) = 0.0_p_r8
    beta(2) = 6.0_p_r8
    gama(2) = 0.0_p_r8
    delta(2) = (6.0_p_r8 / hIn2) * &
      (FieldIn(1) - 2.0_p_r8 * FieldIn(2) + FieldIn(3))
    do iIn = 3, DimIn - 1
      alfa(iIn) = 1.0_p_r8
      beta(iIn) = 4.0_p_r8
      gama(iIn) = 1.0_p_r8
      delta(iIn) = (6.0_p_r8 / hIn2) * &
        (FieldIn(iIn - 1) - 2.0_p_r8 * FieldIn(iIn) + FieldIn(iIn + 1))
    end do
    alfa(DimIn) = 1.0_p_r8
    beta(DimIn) = 4.0_p_r8
    gama(DimIn) = 1.0_p_r8
    delta(DimIn) = (6.0_p_r8 / hIn2) * &
      (FieldIn(DimIn - 1) - 2.0_p_r8 * FieldIn(DimIn) + FieldIn(1))
    alfa(DimIn + 1) = 0.0_p_r8
    beta(DimIn + 1) = 6.0_p_r8
    gama(DimIn + 1) = 0.0_p_r8
    delta(DimIn + 1) = (6.0_p_r8 / hIn2) * &
      (FieldIn(DimIn) - 2.0_p_r8 * FieldIn(1) + FieldIn(2))

    ! backward elimination

    do iIn = 3, DimIn
      beta(iIn) = beta(iIn) - (gama(iIn - 1) * alfa(iIn)) / beta(iIn - 1)
      delta(iIn) = delta(iIn) - (delta(iIn - 1) * alfa(iIn)) / beta(iIn - 1)
    end do

    ! forward substitution

    der(DimIn + 1) = delta(DimIn + 1) / beta(DimIn + 1)
    do iIn = DimIn, 2, -1
      der(iIn) = (delta(iIn) - gama(iIn) * der(iIn + 1)) / beta(iIn)
    end do

    der(1) = 2.0_p_r8 * der(2) - der(3)
    der(DimIn + 2) = 2.0_p_r8 * der(DimIn + 1) - der(DimIn)

    ! interpolation

    do iIn = 1, DimIn + 2
      c1(iIn) = der(iIn) / (6.0_p_r8 * hIn)
    end do
    do iIn = 1, DimIn
      c2(iIn) = FieldIn(iIn) / hIn - der(iIn) * (hIn / 6.0_p_r8)
    end do
    c2(DimIn + 1) = FieldIn(1) / hIn - der(DimIn + 1) * (hIn / 6.0_p_r8)
    c2(DimIn + 2) = FieldIn(2) / hIn - der(DimIn + 2) * (hIn / 6.0_p_r8)
    do iOut = ifirst, ilast
      iIn = int(real(iOut - 1, p_r8) * ratio) + 2
      dxm = real(iOut - 1, p_r8) * hOut - real(iIn - 2, p_r8) * hIn
      dxm2 = dxm * dxm
      dx = real(iIn - 1, p_r8) * hIn - real(iOut - 1, p_r8) * hOut
      dx2 = dx * dx
      FieldOut(iOut) = &
        (dxm2 * c1(iIn) + c2(iIn)) * dxm + &
          (dx2 * c1(iIn - 1) + c2(iIn - 1)) * dx
    end do
  end subroutine CyclicCubicSpline


  subroutine LinearIJtoIBJB_R2D(FieldIn, FieldOut)
    !# Linear IJ into IBJB_R2D
    !# ---
    !# @info
    !# **Brief:** Linear IJ into IBJB_R2D. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (iMax, jMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(ibMax, jbMax)

    real(kind = p_r8) :: Fout(imax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**LinearIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,i,ib,fout,jb)
      do j = myfirstlat, mylastlat
        ifirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicLinear(iMax, iMaxPerJ(j), &
          FieldIn(1, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine LinearIJtoIBJB_R2D


  subroutine LinearIBJBtoIJ_R2D(FieldIn, FieldOut)  
    !# Linear IBJB into IJ_R2D
    !# ---
    !# @info
    !# **Brief:** Linear IBJB into IJ_R2D_R2D. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (ibMax, jbMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(iMax, jMax)

    integer :: i
    integer :: iFirst
    integer :: ib
    integer :: j
    integer :: jl
    integer :: jb

    character(len = *), parameter :: h = "**LinearIBJBtoIJ**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
      do j = myfirstlat, mylastlat
        jl = j - myfirstlat + 1
        iFirst = ibPerIJ(1, j)
        jb = jbPerIJ(1, j)
        !         CALL CyclicLinear(iMaxPerJ(j), iMax, &
        !              FieldIn(iFirst,jb), FieldOut(1,jl))
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb) - myfirstlat + 1
          FieldOut(i, j) = FieldIn(ib, jb)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine LinearIBJBtoIJ_R2D
 

  subroutine CyclicLinear (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Linear interpolation on cyclic
    !# ---
    !# @info
    !# **Brief:** Linear interpolation on cyclic, equally spaced data over [0:2pi].
    !# <ul type="disc">
    !#  <li>Input data: <ul type="disc"><li>DimIn: How many input points: abcissae
    !# are supposed to be x(i), i=1,...,DimIn, with x(i) = 2*pi*(i-1)/DimIn. </li>
    !#                                  <li>DimOut: How many output points: abcissae
    !# are supposed to be x(i), i=1,...,DimOut, with x(i) = 2*pi*(i-1)/DimOut. </li>
    !#                                  <li>FieldIn: function values at the DimIn abcissae. </li></ul></li>
    !#  <li>Output data: <ul type="disc"> FieldOut: function values at the DimOut abcissae. </ul></li>
    !#  <li>Requirements: <ul type="disc">DimIn >= 1. Program Stops if DimIn < 1. </ul></li>
    !# </ul>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real(kind = p_r8), intent(in) :: FieldIn(DimIn)
    real(kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: c(DimIn + 1)
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: dxm
    real(kind = p_r8) :: dx
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hInInv
    real(kind = p_r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicLinear)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hInInv = 1.0_p_r8 / hIn
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    ! interpolation

    do iIn = 1, DimIn
      c(iIn) = FieldIn(iIn) * hInInv
    end do
    c(DimIn + 1) = FieldIn(1) * hInInv

    do iOut = ifirst, ilast
      iIn = int(real(iOut - 1, p_r8) * ratio) + 2
      dxm = real(iOut - 1, p_r8) * hOut - real(iIn - 2, p_r8) * hIn
      dx = real(iIn - 1, p_r8) * hIn - real(iOut - 1, p_r8) * hOut
      FieldOut(iOut) = dxm * c(iIn) + dx * c(iIn - 1)
    end do
  end subroutine CyclicLinear


  subroutine CyclicLinear_r4 (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Linear interpolation on cyclic r4
    !# ---
    !# @info
    !# **Brief:** Linear interpolation on cyclic r4 </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real(kind = p_r8), intent(in) :: FieldIn(DimIn)
    real(kind = p_r4), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: c(DimIn + 1)
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: dxm
    real(kind = p_r8) :: dx
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hInInv
    real(kind = p_r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicLinear)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hInInv = 1.0_p_r8 / hIn
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    ! interpolation

    do iIn = 1, DimIn
      c(iIn) = FieldIn(iIn) * hInInv
    end do
    c(DimIn + 1) = FieldIn(1) * hInInv

    do iOut = ifirst, ilast
      iIn = int(real(iOut - 1, p_r8) * ratio) + 2
      dxm = real(iOut - 1, p_r8) * hOut - real(iIn - 2, p_r8) * hIn
      dx = real(iIn - 1, p_r8) * hIn - real(iOut - 1, p_r8) * hOut
      FieldOut(iOut) = dxm * c(iIn) + dx * c(iIn - 1)
    end do
  end subroutine CyclicLinear_r4


  subroutine CyclicLinear_r8 (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Linear interpolation on cyclic r8
    !# ---
    !# @info
    !# **Brief:** Linear interpolation on cyclic r8 </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real(kind = p_r8), intent(in) :: FieldIn(DimIn)
    real(kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: c(DimIn + 1)
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: dxm
    real(kind = p_r8) :: dx
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hInInv
    real(kind = p_r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicLinear)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hInInv = 1.0_p_r8 / hIn
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    ! interpolation

    do iIn = 1, DimIn
      c(iIn) = FieldIn(iIn) * hInInv
    end do
    c(DimIn + 1) = FieldIn(1) * hInInv

    do iOut = ifirst, ilast
      iIn = int(real(iOut - 1, p_r8) * ratio) + 2
      dxm = real(iOut - 1, p_r8) * hOut - real(iIn - 2, p_r8) * hIn
      dx = real(iIn - 1, p_r8) * hIn - real(iOut - 1, p_r8) * hOut
      FieldOut(iOut) = dxm * c(iIn) + dx * c(iIn - 1)
    end do
  end subroutine CyclicLinear_r8
  
  
  subroutine CyclicLinear_ABS (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Linear interpolation on cyclic ABS
    !# ---
    !# @info
    !# **Brief:** Linear interpolation on cyclic ABS </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real(kind = p_r8), intent(in) :: FieldIn(DimIn)
    real(kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: c(DimIn + 1)
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: dxm
    real(kind = p_r8) :: dx
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hInInv
    real(kind = p_r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicLinear_ABS)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hInInv = 1.0_p_r8 / hIn
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    ! interpolation

    do iIn = 1, DimIn
      c(iIn) = abs(FieldIn(iIn)) * hInInv
    end do
    c(DimIn + 1) = abs(FieldIn(1)) * hInInv
    do iOut = ifirst, ilast
      iIn = int(real(iOut - 1, p_r8) * ratio) + 2
      dxm = real(iOut - 1, p_r8) * hOut - real(iIn - 2, p_r8) * hIn
      dx = real(iIn - 1, p_r8) * hIn - real(iOut - 1, p_r8) * hOut
      FieldOut(iOut) = dxm * c(iIn) + dx * c(iIn - 1)
    end do
  end subroutine CyclicLinear_ABS


  subroutine NearestIJtoIBJB_I2D(FieldIn, FieldOut)
    !# Nearest IJ into IBJB_I2D
    !# ---
    !# @info
    !# **Brief:** Nearest IJ into IBJB_I2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: FieldIn (iMax, jMax)
    integer, intent(OUT) :: FieldOut(ibMax, jbMax)

    integer :: Fout(iMax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**NearestIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicNearest_i(iMax, iMaxPerJ(j), &
          FieldIn(1, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine NearestIJtoIBJB_I2D


  subroutine NearestIJtoIBJB_R2D(FieldIn, FieldOut)
    !# Nearest IJ into IBJB_R2D
    !# ---
    !# @info
    !# **Brief:** Nearest IJ into IBJB_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (iMax, jMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(ibMax, jbMax)

    real(kind = p_r8) :: Fout(iMax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**NearestIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicNearest_r(iMax, iMaxPerJ(j), &
          FieldIn(1, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine NearestIJtoIBJB_R2D

  
  subroutine NearestIJtoIBJB_I3D(FieldIn, FieldOut)
    !# Nearest IJ into IBJB_I3D
    !# ---
    !# @info
    !# **Brief:** 3D version by hmjb </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: FieldIn (iMax, kMax, jMax)
    integer, intent(OUT) :: FieldOut(ibMax, kMax, jbMax)

    integer :: FOut(iMax)
    integer :: i, k
    integer :: ifirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**NearestIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb,k)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        do k = 1, kMax
          call CyclicNearest_i(iMax, iMaxPerJ(j), &
            FieldIn(1, k, j), FOut, ifirst, ilast)
          do i = ifirst, ilast
            ib = ibperij(i, j)
            jb = jbperij(i, j)
            FieldOut(ib, k, jb) = Fout(i)
          end do
        enddo
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, :, jb) = FieldIn(i, :, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine NearestIJtoIBJB_I3D


  subroutine NearestIJtoIBJB_R3D(FieldIn, FieldOut)
    !# Nearest IJ into IBJB_R3D
    !# ---
    !# @info
    !# **Brief:** 3D version by hmjb </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (iMax, kMax, jMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(ibMax, kMax, jbMax)

    real(kind = p_r8) :: FOut(iMax)
    integer :: i, k
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**NearestIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb,k)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        do k = 1, kMax
          call CyclicNearest_r(iMax, iMaxPerJ(j), &
            FieldIn(1, k, j), FOut, ifirst, ilast)
          do i = ifirst, ilast
            ib = ibperij(i, j)
            jb = jbperij(i, j)
            FieldOut(ib, k, jb) = Fout(i)
          end do
        enddo
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, :, jb) = FieldIn(i, :, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine NearestIJtoIBJB_R3D


  subroutine NearestIBJBtoIJ_I2D(FieldIn, FieldOut)
    !# Nearest IBJB into IJ_I2D
    !# ---
    !# @info
    !# **Brief:** Nearest IBJB into IJ_I2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: FieldIn (ibMax, jbMax)
    integer, intent(OUT) :: FieldOut(iMax, jMax)

    integer :: i
    integer :: iFirst
    integer :: ib
    integer :: j
    integer :: jl
    integer :: jb

    character(len = *), parameter :: h = "**NearestIBJBtoIJ**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
      do j = myfirstlat, mylastlat
        jl = j - myfirstlat + 1
        iFirst = ibPerIJ(1, j)
        jb = jbPerIJ(1, j)
        !         CALL CyclicNearest_i(iMaxPerJ(j), iMax, &
        !              FieldIn(iFirst,jb), FieldOut(1,jl))
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb) - myfirstlat + 1
          FieldOut(i, j) = FieldIn(ib, jb)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine NearestIBJBtoIJ_I2D


  subroutine NearestIBJBtoIJ_R2D(FieldIn, FieldOut)
    !# Nearest IBJB into IJ_R2D
    !# ---
    !# @info
    !# **Brief:** Nearest IBJB into IJ_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (ibMax, jbMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(iMax, jMax)

    integer :: i
    integer :: iFirst
    integer :: ib
    integer :: j
    integer :: jl
    integer :: jb

    character(len = *), parameter :: h = "**NearestIBJBtoIJ**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
      do j = myfirstlat, mylastlat
        jl = j - myfirstlat + 1
        iFirst = ibPerIJ(1, j)
        jb = jbPerIJ(1, j)
        !         CALL CyclicNearest_r(iMaxPerJ(j), iMax, &
        !              FieldIn(iFirst,jb), FieldOut(1,jl))
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb) - myfirstlat + 1
          FieldOut(i, j) = FieldIn(ib, jb)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine NearestIBJBtoIJ_R2D


  subroutine CyclicNearest_i(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Nearest interpolation on cyclic i
    !# ---
    !# @info
    !# **Brief:** Nearest interpolation on cyclic i </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    integer, intent(in) :: FieldIn (DimIn)
    integer, intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hInInv
    real(kind = p_r8) :: hOut
    real(kind = p_r8) :: difalfa
    real(kind = p_r8) :: alfaIn
    real(kind = p_r8) :: alfaOut
    integer :: mplon(DimOut)

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicNearest)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)
    hInInv = 1.0_p_r8 / hIn

    ! interpolation

    do iOut = ifirst, ilast
      alfaOut = (iOut - 1) * hOut
      difalfa = 1000E+12
      alfaIn = 0.0
      do iIn = 1, DimIn
        difalfa = min(abs(alfaIn - alfaOut), difalfa)
        alfaIn = alfaIn + hIn
      end do
      alfaIn = 0.0
      do iIn = 1, DimIn
        if (abs(alfaIn - alfaOut) == difalfa) then
          mplon(iOut) = iIn
          FieldOut(iOut) = FieldIn(mplon(iOut))
        end if
        alfaIn = alfaIn + hIn
      end do
    end do
    return
  end subroutine CyclicNearest_i


  subroutine CyclicNearest_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Nearest interpolation on cyclic r
    !# ---
    !# @info
    !# **Brief:** Nearest interpolation on cyclic r </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real   (kind = p_r8), intent(in) :: FieldIn (DimIn)
    real   (kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real   (kind = p_r8) :: ratio
    real   (kind = p_r8) :: pi
    real   (kind = p_r8) :: hIn
    real   (kind = p_r8) :: hOut
    real   (kind = p_r8) :: difalfa
    real   (kind = p_r8) :: difinout
    real   (kind = p_r8) :: alfaIn
    real   (kind = p_r8) :: alfaOut
    integer :: lonout

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicNearest)** ERROR: Few input data points"
    end if

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    ! interpolation

    do iOut = ifirst, ilast
      alfaOut = (iOut - 1) * hOut
      alfaIn = 0.0
      difalfa = min (alfaOut, 2.0_p_r8 * pi - alfaOut)
      lonout = 1
      do iIn = 2, DimIn
        alfaIN = (iIn - 1) * hIn
        difinout = abs(alfaIn - alfaOut)
        if (difinout.lt.difalfa) then
          difalfa = difinout
          lonout = iIn
        endif
      end do
      FieldOut(iOut) = FieldIn(lonout)
    end do
    return
  end subroutine CyclicNearest_r


  subroutine CyclicNearest_r2(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Nearest interpolation on cyclic r2
    !# ---
    !# @info
    !# **Brief:** Nearest interpolation on cyclic r2 </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real   (kind = p_r8), intent(in) :: FieldIn (DimIn)
    real   (kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real   (kind = p_r8) :: ratio
    real   (kind = p_r8) :: pi
    real   (kind = p_r8) :: hIn
    real   (kind = p_r8) :: hInInv
    real   (kind = p_r8) :: hOut
    real   (kind = p_r8) :: difalfa
    real   (kind = p_r8) :: alfaIn
    real   (kind = p_r8) :: alfaOut
    integer :: mplon(DimOut)

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicNearest)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)
    hInInv = 1.0_p_r8 / hIn

    ! interpolation

    do iOut = ifirst, ilast
      alfaOut = (iOut - 1) * hOut
      difalfa = 1000E+12
      alfaIn = 0.0
      do iIn = 1, DimIn
        difalfa = min(abs(alfaIn - alfaOut), difalfa)
        alfaIn = alfaIn + hIn
      end do
      alfaIn = 0.0
      do iIn = 1, DimIn
        if (abs(alfaIn - alfaOut) == difalfa) then
          mplon(iOut) = iIn
          FieldOut(iOut) = FieldIn(mplon(iOut))
        end if
        alfaIn = alfaIn + hIn
      end do
    end do
    return
  end subroutine CyclicNearest_r2


  subroutine FreqBoxIJtoIBJB_I2D(FieldIn, FieldOut)
    !# FreqBox IJ into IBJB_I2D
    !# ---
    !# @info
    !# **Brief:** FreqBox IJ into IBJB_I2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: FieldIn (iMax, jMax)
    integer, intent(OUT) :: FieldOut(ibMax, jbMax)

    integer :: FOut(iMax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**FreqBoxIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicFreqBox_i(iMax, iMaxPerJ(j), &
          FieldIn(1:, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine FreqBoxIJtoIBJB_I2D


  subroutine FreqBoxIJtoIBJB_R2D(FieldIn, FieldOut)
    !# FreqBox IJ into IBJB_R2D
    !# ---
    !# @info
    !# **Brief:** FreqBox IJ into IBJB_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (iMax, jMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(ibMax, jbMax)

    real(kind = p_r8) :: FOut(iMax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**FreqBoxIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicFreqBox_r(iMax, iMaxPerJ(j), &
          FieldIn(1:, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine FreqBoxIJtoIBJB_R2D


  subroutine CyclicFreqBox_i(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# FreqBox on cyclic i
    !# ---
    !# @info
    !# **Brief:** FreqBox on cyclic i </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    implicit none
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    integer, intent(in) :: FieldIn (DimIn)
    integer, intent(OUT) :: FieldOut(DimOut)
    integer, parameter :: ncat = 13 !number of catagories found

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hOut
    integer :: ici
    integer :: ico
    integer :: ioi
    integer :: ioo
    integer :: i
    integer :: k
    logical :: flgin  (5)
    logical :: flgout (5)
    real   (kind = p_r8) :: dwork  (2 * (DimIn + DimOut + 2))
    real   (kind = p_r8) :: wtlon  (DimIn + DimOut + 2)
    integer :: mplon  (DimIn + DimOut + 2, 2)
    real   (kind = p_r8) :: work   (ncat, DimOut)
    real   (kind = p_r8) :: wrk2   (DimOut)
    integer :: undef = 0
    real   (kind = p_r8) :: dof
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: lond
    real   (kind = p_r8) :: wln
    integer :: lni
    integer :: lno
    integer :: nc
    integer :: mm
    integer :: n
    integer :: nn
    integer :: nd
    integer :: mdist  (7)
    integer :: ndist  (7)
    real   (kind = p_r8) :: fm
    integer :: nx
    integer :: kl
    real   (kind = p_r8) :: b      (5)
    real   (kind = p_r8) :: fr
    real   (kind = p_r8) :: cmx
    integer :: kmx
    real   (kind = p_r8) :: fq
    real   (kind = p_r8) :: fmk
    real   (kind = p_r8) :: frk
    integer :: iq
    integer :: nq
    integer :: ns
    integer :: nxk
    integer :: klass  (ncat)
    data KLASS/6*1, 2, 2, 3, 2, 3, 4, 5/

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicFreqBox)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicFreqBox)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1) = .true.
    flgin (2) = .true.
    flgin (3) = .false.
    flgin (4) = .true.
    flgin (5) = .true.
    flgout(1) = .true.
    flgout(2) = .true.
    flgout(3) = .false.
    flgout(4) = .true.
    flgout(5) = .true.

    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi = DimIn + DimOut + 2
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    if (flgin(5) .or. flgin(4)) then
      ici = 0
      dof = 0.5_p_r8
    else
      ici = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimIn
      dwork(i + ioi) = (dof + dble(i - 1)) * hIn
    end do
    !
    !     output grid longitudes
    !
    ioo = 2 * DimIn + DimOut + 3
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    if (flgout(5) .or. flgout(4)) then
      ico = 0
      dof = 0.5_p_r8
    else
      ico = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimOut
      dwork(i + ioo) = (dof + dble(i - 1)) * hOut
    end do
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1 = 1
    i2 = 1
    i3 = 1
    do
      if (dwork(i1 + ioi) == dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1) then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if(i1 - ici > DimIn / 2)mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i2 = i2 + 1
        i3 = i3 + 1
      else if (dwork(i1 + ioi) < dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i3 = i3 + 1
      else
        dwork(i3) = dwork(i2 + ioo)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i2 = i2 + 1
        i3 = i3 + 1
      end if
      if ((i1 > DimIn) .or. (i2 > DimOut)) EXIT
    end do

    if (i1 > DimIn) i1 = 1
    if (i2 > DimOut) i2 = 1
    do
      if (i2 /= 1) then
        dwork(i3) = dwork(i2 + ioo)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = 1
        if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + 1
          if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
        end if
        mplon(i3 - 1, 2) = i2 - ico
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
          if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
        end if
        i2 = i2 + 1
        if (i2 > DimOut)i2 = 1
        i3 = i3 + 1
      end if
      if (i1 /= 1)then
        dwork(i3) = dwork(i1 + ioi)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = i1 - ici
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
          if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
        end if
        mplon(i3 - 1, 2) = 1
        if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + 1
          if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
        end if
        i1 = i1 + 1
        if (i1 > DimIn)i1 = 1
        i3 = i3 + 1
      end if
      if (i1 == 1 .and. i2 == 1) EXIT
    end do
    wtlon(i3 - 1) = 2.0_p_r8 * pi + dwork(1) - dwork(i3 - 1)
    mplon(i3 - 1, 1) = 1
    if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
    if (.not.flgin(2)) then
      mplon(i3 - 1, 1) = DimIn / 2 + 1
      if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
    end if
    mplon(i3 - 1, 2) = 1
    if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
    if (.not.flgout(2)) then
      mplon(i3 - 1, 2) = DimOut / 2 + 1
      if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
    end if
    lond = i3 - 1

    ! interpolation

    mdist(1:7) = 0
    ndist(1:7) = 0
    FieldOut(ifirst:ilast) = 0.0_p_r8
    wrk2  (1:DimOut) = 0.0_p_r8
    work(1:ncat, 1:DimOut) = 0.0_p_r8

    do i = 1, lond

      wln = wtlon(i)
      lni = mplon(i, 1)
      lno = mplon(i, 2)

      if (FieldIn(lni) == undef)    cycle

      nc = FieldIn(lni)

      if (nc > ncat.or.lno > DimOut) then
        write(p_nfprt, *)nc, lno, i, lni
        stop 'ERROR in nc,lno,i and lni point'
      end if

      if (nc.lt.1.or.lno.lt.1) then
        write(p_nfprt, *)nc, lno, i, lni
        stop 'ERROR in nc,lno,i and lni point'
      end if
      work(nc, lno) = work(nc, lno) + wln
      wrk2(lno) = wrk2(lno) + wln
    end do

    fq = 1.0_p_r8
    nd = 0
    ns = 0

    do i = ifirst, ilast

      FieldOut(i) = undef

      if (wrk2(i) == 0.0_p_r8) cycle

      fm = 0.0_p_r8
      nx = undef
      mm = 0
      nn = 1
      b(1) = 0.0_p_r8
      b(2) = 0.0_p_r8
      b(3) = 0.0_p_r8
      b(4) = 0.0_p_r8
      b(5) = 0.0_p_r8

      do n = 1, ncat
        fr = work(n, i) / wrk2(i)

        if (fm < fr) then
          fm = fr
          nx = n
        end if

        kl = klass(n)
        b(kl) = b(kl) + fr

        if (fr > 0.5_p_r8) nn = 0
        if (work(n, i).ne.0.0_p_r8) mm = mm + 1
      end do

      cmx = 0.0_p_r8
      kmx = 0

      do k = 1, 5
        if (b(k) > cmx) then
          cmx = b(k)
          kmx = k
        end if
      end do

      if (klass(nx) == kmx) then
        FieldOut(i) = nx
        nd = nd + 1
        if (fm.ne.0.0_p_r8.and.fm < fq) then
          fq = fm
          iq = i
          nq = nx
        end if

      else

        fmk = 0.0_p_r8

        nxk = 1 ! avoid not initializated
        do n = 1, ncat
          if (klass(n).ne.kmx) cycle
          frk = work(n, i) / wrk2(i)
          if (fmk.lt.frk) then
            fmk = frk
            nxk = n
          end if
        end do

        FieldOut(i) = nxk
        ns = ns + 1
        if (fmk.ne.0.0_p_r8.and.fm.lt.fq) then
          fq = fmk
          iq = i
          nq = nxk
        end if
      end if

      if (mm.gt.7.and.mm.gt.0)mm = 7

      mdist(mm) = mdist(mm) + 1
      ndist(mm) = ndist(mm) + nn
    end do
    return
  end subroutine CyclicFreqBox_i


  subroutine CyclicFreqBox_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# FreqBox on cyclic r
    !# ---
    !# @info
    !# **Brief:** FreqBox on cyclic r </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    implicit none
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real   (kind = p_r8), intent(in) :: FieldIn (DimIn)
    real   (kind = p_r8), intent(OUT) :: FieldOut(DimOut)
    integer, parameter :: ncat = 13 !number of catagories found

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hOut
    integer :: ici
    integer :: ico
    integer :: ioi
    integer :: ioo
    integer :: i
    integer :: k
    logical :: flgin  (5)
    logical :: flgout (5)
    real   (kind = p_r8) :: dwork  (2 * (DimIn + DimOut + 2))
    real   (kind = p_r8) :: wtlon  (DimIn + DimOut + 2)
    integer :: mplon  (DimIn + DimOut + 2, 2)
    real   (kind = p_r8) :: work   (ncat, DimOut)
    real   (kind = p_r8) :: wrk2   (DimOut)
    real   (kind = p_r8) :: undef = 0.0
    real   (kind = p_r8) :: dof
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: lond
    real   (kind = p_r8) :: wln
    integer :: lni
    integer :: lno
    integer :: nc
    integer :: mm
    integer :: n
    integer :: nn
    integer :: nd
    integer :: mdist  (7)
    integer :: ndist  (7)
    real   (kind = p_r8) :: fm
    integer :: nx
    integer :: kl
    real   (kind = p_r8) :: b      (5)
    real   (kind = p_r8) :: fr
    real   (kind = p_r8) :: cmx
    integer :: kmx
    real   (kind = p_r8) :: fq
    real   (kind = p_r8) :: fmk
    real   (kind = p_r8) :: frk
    integer :: iq
    integer :: nq
    integer :: ns
    integer :: nxk
    integer :: klass  (ncat)
    data KLASS/6*1, 2, 2, 3, 2, 3, 4, 5/

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicLinear)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1) = .true.
    flgin (2) = .true.
    flgin (3) = .false.
    flgin (4) = .true.
    flgin (5) = .true.
    flgout(1) = .true.
    flgout(2) = .true.
    flgout(3) = .false.
    flgout(4) = .true.
    flgout(5) = .true.
    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi = DimIn + DimOut + 2
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    if (flgin(5) .or. flgin(4)) then
      ici = 0
      dof = 0.5_p_r8
    else
      ici = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimIn
      dwork(i + ioi) = (dof + dble(i - 1)) * hIn
    end do
    !
    !     output grid longitudes
    !
    ioo = 2 * DimIn + DimOut + 3
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    if (flgout(5) .or. flgout(4)) then
      ico = 0
      dof = 0.5_p_r8
    else
      ico = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimOut
      dwork(i + ioo) = (dof + dble(i - 1)) * hOut
    end do
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1 = 1
    i2 = 1
    i3 = 1
    do
      if (dwork(i1 + ioi) == dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1) then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if(i1 - ici > DimIn / 2)mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i2 = i2 + 1
        i3 = i3 + 1
      else if (dwork(i1 + ioi) < dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i3 = i3 + 1
      else
        dwork(i3) = dwork(i2 + ioo)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i2 = i2 + 1
        i3 = i3 + 1
      end if
      if ((i1 > DimIn) .or. (i2 > DimOut)) EXIT
    end do

    if (i1 > DimIn) i1 = 1
    if (i2 > DimOut) i2 = 1
    do
      if (i2 /= 1) then
        dwork(i3) = dwork(i2 + ioo)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = 1
        if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + 1
          if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
        end if
        mplon(i3 - 1, 2) = i2 - ico
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
          if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
        end if
        i2 = i2 + 1
        if (i2 > DimOut)i2 = 1
        i3 = i3 + 1
      end if
      if (i1 /= 1)then
        dwork(i3) = dwork(i1 + ioi)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = i1 - ici
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
          if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
        end if
        mplon(i3 - 1, 2) = 1
        if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + 1
          if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
        end if
        i1 = i1 + 1
        if (i1 > DimIn)i1 = 1
        i3 = i3 + 1
      end if
      if (i1 == 1 .and. i2 == 1) EXIT
    end do
    wtlon(i3 - 1) = 2.0_p_r8 * pi + dwork(1) - dwork(i3 - 1)
    mplon(i3 - 1, 1) = 1
    if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
    if (.not.flgin(2)) then
      mplon(i3 - 1, 1) = DimIn / 2 + 1
      if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
    end if
    mplon(i3 - 1, 2) = 1
    if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
    if (.not.flgout(2)) then
      mplon(i3 - 1, 2) = DimOut / 2 + 1
      if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
    end if
    lond = i3 - 1

    ! interpolation

    mdist(1:7) = 0
    ndist(1:7) = 0
    FieldOut(ifirst:ilast) = 0.0_p_r8
    wrk2  (1:DimOut) = 0.0_p_r8
    work(1:ncat, 1:DimOut) = 0.0_p_r8

    do i = 1, lond

      wln = wtlon(i)
      lni = mplon(i, 1)
      lno = mplon(i, 2)

      if (FieldIn(lni) == undef)    cycle

      nc = int(FieldIn(lni))

      if (nc > ncat.or.lno > DimOut) then
        write(p_nfprt, *)nc, lno, i, lni
        stop 'ERROR in nc,lno,i and lni point'
      end if

      if (nc.lt.1.or.lno.lt.1) then
        write(p_nfprt, *)nc, lno, i, lni
        stop 'ERROR in nc,lno,i and lni point'
      end if
      work(nc, lno) = work(nc, lno) + wln
      wrk2(lno) = wrk2(lno) + wln
    end do

    fq = 1.0_p_r8
    nd = 0
    ns = 0

    do i = ifirst, ilast

      FieldOut(i) = undef

      if (wrk2(i) == 0.0_p_r8) cycle

      fm = 0.0_p_r8
      nx = undef
      mm = 0
      nn = 1
      b(1) = 0.0_p_r8
      b(2) = 0.0_p_r8
      b(3) = 0.0_p_r8
      b(4) = 0.0_p_r8
      b(5) = 0.0_p_r8

      do n = 1, ncat
        fr = work(n, i) / wrk2(i)

        if (fm < fr) then
          fm = fr
          nx = n
        end if

        kl = klass(n)
        b(kl) = b(kl) + fr

        if (fr > 0.5_p_r8) nn = 0
        if (work(n, i).ne.0.0_p_r8) mm = mm + 1
      end do

      cmx = 0.0_p_r8
      kmx = 0

      do k = 1, 5
        if (b(k) > cmx) then
          cmx = b(k)
          kmx = k
        end if
      end do

      if (klass(nx) == kmx) then
        FieldOut(i) = nx
        nd = nd + 1
        if (fm.ne.0.0_p_r8.and.fm < fq) then
          fq = fm
          iq = i
          nq = nx
        end if

      else

        fmk = 0.0_p_r8

        ! TODO - check this initialization (before was no initialization)
        nxk = 1
        do n = 1, ncat
          if (klass(n).ne.kmx) cycle
          frk = work(n, i) / wrk2(i)
          if (fmk.lt.frk) then
            fmk = frk
            nxk = n
          end if
        end do

        FieldOut(i) = nxk
        ns = ns + 1
        if (fmk.ne.0.0_p_r8.and.fm.lt.fq) then
          fq = fmk
          iq = i
          nq = nxk
        end if
      end if

      if (mm.gt.7.and.mm.gt.0)mm = 7

      mdist(mm) = mdist(mm) + 1
      ndist(mm) = ndist(mm) + nn
    end do
    return
  end subroutine CyclicFreqBox_r


  subroutine SeaMaskIJtoIBJB_R2D(FieldIn, FieldOut)
    !# Sea Mask IJ into IBJB_R2D
    !# ---
    !# @info
    !# **Brief:** Sea Mask IJ into IBJB_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (iMax, jMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(ibMax, jbMax)

    real(kind = p_r8) :: FOut(iMax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**SeaMaskIJtoIBJB**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicSeaMask_r(iMax, iMaxPerJ(j), &
          FieldIn(1, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine SeaMaskIJtoIBJB_R2D


  subroutine SeaMaskIBJBtoIJ_R2D(FieldIn, FieldOut)
    !# Sea Mask IBJB into IJ_R2D
    !# ---
    !# @info
    !# **Brief:** Sea Mask IBJB into IJ_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (ibMax, jbMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(iMax, jMax)

    integer :: i
    integer :: iFirst
    integer :: ib
    integer :: j
    integer :: jl
    integer :: jb

    character(len = *), parameter :: h = "**SeaMaskIBJBtoIJ**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
      do j = myfirstlat, mylastlat
        jl = j - myfirstlat + 1
        iFirst = ibPerIJ(1, j)
        jb = jbPerIJ(1, j)
        !         CALL CyclicSeaMask_r(iMaxPerJ(j), iMax, &
        !              FieldIn(iFirst,jb), FieldOut(1,jl))
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb) - myfirstlat + 1
          FieldOut(i, j) = FieldIn(ib, jb)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine SeaMaskIBJBtoIJ_R2D


  subroutine CyclicSeaMask_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# Sea Mask on cyclic r
    !# ---
    !# @info
    !# **Brief:** Sea Mask on cyclic r </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real   (kind = p_r8), intent(in) :: FieldIn (DimIn)
    real   (kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hOut
    integer :: ici
    integer :: ico
    integer :: ioi
    integer :: ioo
    integer :: i
    logical :: flgin  (5)
    logical :: flgout (5)
    real   (kind = p_r8) :: dwork  (2 * (DimIn + DimOut + 2))
    real   (kind = p_r8) :: wtlon  (DimIn + DimOut + 2)
    integer :: mplon  (DimIn + DimOut + 2, 2)
    real   (kind = p_r8) :: work   (DimOut)
    real   (kind = p_r8) :: undef = 290.0_p_r8
    real   (kind = p_r8) :: dof
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: lond
    real   (kind = p_r8) :: wln
    integer :: lni
    integer :: lno

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicLinear)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1) = .true.
    flgin (2) = .true.
    flgin (3) = .false.
    flgin (4) = .true.
    flgin (5) = .true.
    flgout(1) = .true.
    flgout(2) = .true.
    flgout(3) = .false.
    flgout(4) = .true.
    flgout(5) = .true.

    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi = DimIn + DimOut + 2
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    if (flgin(5) .or. flgin(4)) then
      ici = 0
      dof = 0.5_p_r8
    else
      ici = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimIn
      dwork(i + ioi) = (dof + dble(i - 1)) * hIn
    end do
    !
    !     output grid longitudes
    !
    ioo = 2 * DimIn + DimOut + 3
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    if (flgout(5) .or. flgout(4)) then
      ico = 0
      dof = 0.5_p_r8
    else
      ico = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimOut
      dwork(i + ioo) = (dof + (i - 1)) * hOut
    end do
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1 = 1
    i2 = 1
    i3 = 1
    do
      if (dwork(i1 + ioi) == dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1) then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if(i1 - ici > DimIn / 2)mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i2 = i2 + 1
        i3 = i3 + 1
      else if (dwork(i1 + ioi) < dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i3 = i3 + 1
      else
        dwork(i3) = dwork(i2 + ioo)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i2 = i2 + 1
        i3 = i3 + 1
      end if
      if ((i1 > DimIn) .or. (i2 > DimOut)) EXIT
    end do

    if (i1 > DimIn) i1 = 1
    if (i2 > DimOut) i2 = 1
    do
      if (i2 /= 1) then
        dwork(i3) = dwork(i2 + ioo)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = 1
        if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + 1
          if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
        end if
        mplon(i3 - 1, 2) = i2 - ico
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
          if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
        end if
        i2 = i2 + 1
        if (i2 > DimOut)i2 = 1
        i3 = i3 + 1
      end if
      if (i1 /= 1)then
        dwork(i3) = dwork(i1 + ioi)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = i1 - ici
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
          if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
        end if
        mplon(i3 - 1, 2) = 1
        if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + 1
          if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
        end if
        i1 = i1 + 1
        if (i1 > DimIn)i1 = 1
        i3 = i3 + 1
      end if
      if (i1 == 1 .and. i2 == 1) EXIT
    end do
    wtlon(i3 - 1) = 2.0_p_r8 * pi + dwork(1) - dwork(i3 - 1)
    mplon(i3 - 1, 1) = 1
    if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
    if (.not.flgin(2)) then
      mplon(i3 - 1, 1) = DimIn / 2 + 1
      if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
    end if
    mplon(i3 - 1, 2) = 1
    if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
    if (.not.flgout(2)) then
      mplon(i3 - 1, 2) = DimOut / 2 + 1
      if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
    end if
    lond = i3 - 1

    ! interpolation

    FieldOut(1:DimOut) = 0.0_p_r8
    work(1:DimOut) = 0.0_p_r8

    do i = 1, lond
      lni = mplon(i, 1)
      if (FieldIn(lni) < 0.0_p_r8) then ! OBS 0.0 Kelvin valor minino de
        ! temperatura usada na interpolacao
        wln = wtlon(i)
        lno = mplon(i, 2)
        FieldOut(lno) = FieldOut(lno) + FieldIn(lni) * wln
        work(lno) = work(lno) + wln
      end if
    end do

    do i = ifirst, ilast
      if (work(i) == 0.0_p_r8)then
        FieldOut(i) = undef
      else
        FieldOut(i) = FieldOut(i) / work(i)
      end if
    end do

  end subroutine CyclicSeaMask_r


  subroutine AveBoxIJtoIBJB_R2D(FieldIn, FieldOut)
    !# AveBox IJ into IBJB_R2D
    !# ---
    !# @info
    !# **Brief:** AveBox IJ into IBJB_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (iMax, jMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(ibMax, jbMax)

    real(kind = p_r8) :: FOut(iMax)
    integer :: i
    integer :: iFirst
    integer :: ilast
    integer :: ib
    integer :: j
    integer :: jb

    character(len = *), parameter :: h = "**AveBoxIJtoIBJB_R2D**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
      do j = myfirstlat, mylastlat
        iFirst = myfirstlon(j)
        ilast = mylastlon(j)
        call CyclicAveBox_r(iMax, iMaxPerJ(j), &
          FieldIn(1, j), FOut, ifirst, ilast)
        do i = ifirst, ilast
          ib = ibperij(i, j)
          jb = jbperij(i, j)
          FieldOut(ib, jb) = Fout(i)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          FieldOut(ib, jb) = FieldIn(i, j)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine AveBoxIJtoIBJB_R2D


  subroutine AveBoxIBJBtoIJ_R2D(FieldIn, FieldOut)
    !# AveBox IBJB into IJ_R2D
    !# ---
    !# @info
    !# **Brief:** AveBox IBJB into IJ_R2D </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    real(kind = p_r8), intent(in) :: FieldIn (ibMax, jbMax)
    real(kind = p_r8), intent(OUT) :: FieldOut(iMax, jMax)

    integer :: i
    integer :: iFirst
    integer :: ib
    integer :: j
    integer :: jl
    integer :: jb

    character(len = *), parameter :: h = "**AveBoxIBJBtoIJ_R2D**"

    if (reducedGrid) then
      !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
      do j = myfirstlat, mylastlat
        jl = j - myfirstlat + 1
        iFirst = ibPerIJ(1, j)
        jb = jbPerIJ(1, j)
        !         CALL CyclicAveBox_r(iMaxPerJ(j), iMax, &
        !              FieldIn(iFirst,jb), FieldOut(1,jl))
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ib,i,j)
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb) - myfirstlat + 1
          FieldOut(i, j) = FieldIn(ib, jb)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine AveBoxIBJBtoIJ_R2D


  subroutine CyclicAveBox_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    !# AveBox on cyclic r
    !# ---
    !# @info
    !# **Brief:** AveBox on cyclic r </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    integer, intent(in) :: DimIn
    integer, intent(in) :: DimOut
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real   (kind = p_r8), intent(in) :: FieldIn (DimIn)
    real   (kind = p_r8), intent(OUT) :: FieldOut(DimOut)

    integer :: iIn
    integer :: iOut
    integer :: iRatio
    real(kind = p_r8) :: ratio
    real(kind = p_r8) :: pi
    real(kind = p_r8) :: hIn
    real(kind = p_r8) :: hOut
    integer :: ici
    integer :: ico
    integer :: ioi
    integer :: ioo
    integer :: i
    logical :: flgin  (5)
    logical :: flgout (5)
    real   (kind = p_r8) :: dwork  (2 * (DimIn + DimOut + 2))
    real   (kind = p_r8) :: wtlon  (DimIn + DimOut + 2)
    integer :: mplon  (DimIn + DimOut + 2, 2)
    real   (kind = p_r8) :: work   (DimOut)
    real   (kind = p_r8) :: undef = -999.0_p_r8
    real   (kind = p_r8) :: dof
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: lond
    real   (kind = p_r8) :: wln
    integer :: lni
    integer :: lno

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn / DimOut
    if (iRatio * DimOut == DimIn) then
      do iOut = ifirst, ilast
        iIn = iRatio * (iOut - 1) + 1
        FieldOut(iOut) = FieldIn(iIn)
      end do
      return
    end if

    ! protection: input data size should be at least 1

    if (DimIn < 1) then
      stop "**(CyclicLinear)** ERROR: Few input data points"
    end if

    ! protection: output data should fit into input data + 2 intervals

    ratio = real(DimIn, p_r8) / real(DimOut, p_r8)
    iIn = int(real(DimOut - 1, p_r8) * ratio) + 2
    if (iIn > DimIn + 1) then
      stop "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    end if

    ! initialization

    pi = 4.0_p_r8 * atan(1.0_p_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1) = .true.
    flgin (2) = .true.
    flgin (3) = .false.
    flgin (4) = .true.
    flgin (5) = .true.
    flgout(1) = .true.
    flgout(2) = .true.
    flgout(3) = .false.
    flgout(4) = .true.
    flgout(5) = .true.

    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi = DimIn + DimOut + 2
    hIn = (2.0_p_r8 * pi) / real(DimIn, p_r8)
    if (flgin(5) .or. flgin(4)) then
      ici = 0
      dof = 0.5_p_r8
    else
      ici = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimIn
      dwork(i + ioi) = (dof + dble(i - 1)) * hIn
    end do
    !
    !     output grid longitudes
    !
    ioo = 2 * DimIn + DimOut + 3
    hOut = (2.0_p_r8 * pi) / real(DimOut, p_r8)

    if (flgout(5) .or. flgout(4)) then
      ico = 0
      dof = 0.5_p_r8
    else
      ico = 1
      dof = 0.0_p_r8
    end if
    do i = 1, DimOut
      dwork(i + ioo) = (dof + (i - 1)) * hOut
    end do
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1 = 1
    i2 = 1
    i3 = 1
    do
      if (dwork(i1 + ioi) == dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1) then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if(i1 - ici > DimIn / 2)mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i2 = i2 + 1
        i3 = i3 + 1
      else if (dwork(i1 + ioi) < dwork(i2 + ioo)) then
        dwork(i3) = dwork(i1 + ioi)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i1 = i1 + 1
        i3 = i3 + 1
      else
        dwork(i3) = dwork(i2 + ioo)
        if (i3 /= 1)then
          wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if (.not.flgin(2)) then
            mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
            if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flgout(2)) then
            mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
            if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
          end if
        end if
        i2 = i2 + 1
        i3 = i3 + 1
      end if
      if ((i1 > DimIn) .or. (i2 > DimOut)) EXIT
    end do

    if (i1 > DimIn) i1 = 1
    if (i2 > DimOut) i2 = 1
    do
      if (i2 /= 1) then
        dwork(i3) = dwork(i2 + ioo)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = 1
        if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + 1
          if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
        end if
        mplon(i3 - 1, 2) = i2 - ico
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + i2 - ico
          if (i2 - ico > DimOut / 2) mplon(i3 - 1, 2) = i2 - ico - DimOut / 2
        end if
        i2 = i2 + 1
        if (i2 > DimOut)i2 = 1
        i3 = i3 + 1
      end if
      if (i1 /= 1)then
        dwork(i3) = dwork(i1 + ioi)
        wtlon(i3 - 1) = dwork(i3) - dwork(i3 - 1)
        mplon(i3 - 1, 1) = i1 - ici
        if (.not.flgin(2)) then
          mplon(i3 - 1, 1) = DimIn / 2 + i1 - ici
          if (i1 - ici > DimIn / 2) mplon(i3 - 1, 1) = i1 - ici - DimIn / 2
        end if
        mplon(i3 - 1, 2) = 1
        if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
        if (.not.flgout(2)) then
          mplon(i3 - 1, 2) = DimOut / 2 + 1
          if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
        end if
        i1 = i1 + 1
        if (i1 > DimIn)i1 = 1
        i3 = i3 + 1
      end if
      if (i1 == 1 .and. i2 == 1) EXIT
    end do
    wtlon(i3 - 1) = 2.0_p_r8 * pi + dwork(1) - dwork(i3 - 1)
    mplon(i3 - 1, 1) = 1
    if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn
    if (.not.flgin(2)) then
      mplon(i3 - 1, 1) = DimIn / 2 + 1
      if (.not.(flgin(4) .or. flgin(5))) mplon(i3 - 1, 1) = DimIn / 2
    end if
    mplon(i3 - 1, 2) = 1
    if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut
    if (.not.flgout(2)) then
      mplon(i3 - 1, 2) = DimOut / 2 + 1
      if (.not.(flgout(4) .or. flgout(5))) mplon(i3 - 1, 2) = DimOut / 2
    end if
    lond = i3 - 1

    ! interpolation

    FieldOut(1:DimOut) = 0.0_p_r8
    work(1:DimOut) = 0.0_p_r8

    do i = 1, lond
      lni = mplon(i, 1)
      if (FieldIn(lni) /= undef) then
        wln = wtlon(i)
        lno = mplon(i, 2)
        FieldOut(lno) = FieldOut(lno) + FieldIn(lni) * wln
        work(lno) = work(lno) + wln
      end if
    end do

    do i = ifirst, ilast
      if (work(i) == 0.0_p_r8)then
        FieldOut(i) = undef
      else
        FieldOut(i) = FieldOut(i) / work(i)
      end if
    end do

  end subroutine CyclicAveBox_r


  subroutine Clear_Utils()
    !# Cleans utils
    !# ---
    !# @info
    !# **Brief:** Cleans utils. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    deallocate(GaussColat)
    deallocate(SinGaussColat)
    deallocate(CosGaussColat)
    deallocate(AuxGaussColat)
    deallocate(GaussPoints)
    deallocate(GaussWeights)
    deallocate(auxpol)
    deallocate(colrad)
    deallocate(rcs2)
    deallocate(colrad2D)
    deallocate(cos2lat)
    deallocate(ercossin)
    deallocate(fcor)
    deallocate(cosiv)
    deallocate(coslatj)
    deallocate(sinlatj)
    deallocate(coslat)
    deallocate(sinlat)
    deallocate(coslon)
    deallocate(sinlon)
    deallocate(longit)
    deallocate(rcl)
    deallocate(lonrad)
    deallocate(lati)
    deallocate(long)
    deallocate(cosz)
    deallocate(cos2d)

    !  DEALLOCATE(Epslon  )
    !  DEALLOCATE(LegFuncS2F  )
    !  DEALLOCATE(Square  )
    !  DEALLOCATE(Den  )

    deallocate (AuxPoly1, AuxPoly2)

    ! CALL DestroyGaussQuad()
    !  DEALLOCATE (GaussColat)
    !  DEALLOCATE (SinGaussColat)
    !  DEALLOCATE (AuxGaussColat)
    !  DEALLOCATE (GaussPoints)
    !  DEALLOCATE (GaussWeights)

    ! DEALLOCATE (AuxPoly1, AuxPoly2)
    created = .false.

    !DEALLOCATE (ExtDiagPerCol  )   ! diag=DiagPerCol(col )
    !  DEALLOCATE (ExtColPerDiag  )   ! col =ColPerDiag(diag)
    created = .false.
    reducedgrid = .false.
    maxDegree = -1

  end subroutine Clear_Utils

end module Mod_Utils

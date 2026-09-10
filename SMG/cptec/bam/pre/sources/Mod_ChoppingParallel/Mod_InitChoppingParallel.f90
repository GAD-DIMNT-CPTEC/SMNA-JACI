!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_InitChoppingParallel </br></br>
!#
!# **Brief**: Module used for Initialize MPI. </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!# 
!# **Author**: Paulo Bonatti </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2008 - Paulo Bonatti  - version: 1.3.0 </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 2.0.0 - some adaptations for modularizing Chopping </li>
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

module Mod_InitChoppingParallel

  use Mod_Parallelism_Group_Chopping, only : &
    myId &
    , maxNodes

  use Mod_Parallelism_Fourier, only : &
    myId_four &
    , mygroup_four &
    , CreateFourierGroup &
    , maxnodes_four

  use Mod_InputParameters, only : &
    ChoppingNameListData

  use Mod_Sizes, only : ibMax, jbMax, kMax, kmaxloc_out, Kmaxloc_In, &
    iMax, jMax, jMaxHalf, jMinPerM, iMaxPerJ, &
    mMax, nExtMax, mnMax, mnExtMax, mnExtMap, &
    mymnMax, mymnExtMax, kmaxloc, jbMax_ext, &
    RegisterBasicSizes, RegisterOtherSizes, &
    ngroups_four, npperg_four, first_proc_four, map_four, nprocmax_four, &
    nrecs_gr, nsends_gr, nrecs_f, nrecs_g, nrecs_diag, nsends_diag, &
    SpectralDomainDecomp, GridDecomposition, VerticalGroups, &
    ReshapeVerticalGroups

  use Mod_Utils, only : CreateGaussQuad, &
    GaussPoints, GaussWeights, CreateAssocLegFunc, &
    LegFuncS2F, CreateGridValues, DestroyAssocLegFunc, allpolynomials

  use Mod_SpecDynamics, only : InitDZtoUV, InitUvtodz

  use Mod_Transform, only : InitTransform

  use Mod_Communications, only : Set_Communic_buffer, exchange_ftog, Exchange_diag

  implicit none

  private

  public :: InitAll


contains


  !  subroutine InitAll (nlon, trunc_in, trunc_out, vert, vert_out, givenfouriergroups, nproc_vert, ibdim_size)
  !  (ImaxOut, dtoChp%mEndInp, dtoChp%mEndOut, dtoChp%kMaxInp, dtoChp%kMaxOut, dtoChp%givenfouriergroups, dtoChp%nproc_vert, dtoChp%ibdim_size)
  subroutine InitAll (nlon, dtoChp)
    !# Initializes All
    !# ---
    !# @info
    !# **Brief:** Initializes All. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endin
    integer, intent(in) :: nlon
    type(ChoppingNameListData), intent(inout) :: dtoChp

    character(len = *), parameter :: h = "**(InitAll)**"
    integer :: wavesFactor
    integer :: nLat
    integer :: j
    integer, allocatable :: mPerLat(:)
    integer, allocatable :: iMaxPerLat(:)

    nLat = nLon / 2

    ! verify consistency and register basic sizes
    call RegisterBasicSizes(dtoChp%mEndInp, dtoChp%mEndOut, nLat, nLon, dtoChp%kMaxInp)
    ! Initialize GaussQuad (Gaussian points and weights)
    call CreateGaussQuad(jMax)
    ! Spectral Domain Decomposition
    call VerticalGroups(dtoChp%givenfouriergroups, dtoChp%nproc_vert, dtoChp%kMaxOut, kMaxloc_out)

    call ReshapeVerticalGroups(kmax, kMaxloc_in)

    kmaxloc = kMaxloc_in
    call SpectralDomainDecomp()

    ! Initialize AssocLegFunc (computes associated legendre functions)
    allpolynomials = .false.
    call CreateAssocLegFunc(allpolynomials)

    ! computes mPerLat
    if (allocated(mPerLat))then
      deallocate (mPerLat)
      allocate (mPerLat(jMax))
    else
      allocate (mPerLat(jMax))
    end if

    mPerLat = mMax

    ! compute iMaxPerLat

    if (allocated(iMaxPerLat))then
      deallocate (iMaxPerLat)
      allocate(iMaxPerLat(jMax))
    else
      allocate(iMaxPerLat(jMax))
    end if

    iMaxPerLat = iMax

    ! finish building sizes

    call RegisterOtherSizes(iMaxPerLat, mPerLat)

    call CreateFourierGroup(mygroup_four, myid_four)
    ibmax = dtoChp%ibdim_size
    call GridDecomposition(ibmax, jbmax, maxnodes)
    call Exchange_ftog(nrecs_f, nrecs_g)
    call Exchange_diag(nrecs_diag, nsends_diag)

    ! for now on, all global constants defined at sizes can be used

    call CreateGridValues

    ! initialize dztouv and uvtodz

    call InitDZtoUV()
    call InitUvtodz()

    call Set_Communic_buffer()

    call InitTransform()

    ! deallocate legendre functions (already stored in transform,
    ! epslon already used in initializations)

    call DestroyAssocLegFunc()

  end subroutine InitAll

end module Mod_InitChoppingParallel

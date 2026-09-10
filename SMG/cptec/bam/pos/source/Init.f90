!
!  $Author: bonatti $
!  $Date: 2008/11/13 16:44:29 $
!  $Revision: 1.3 $
!
MODULE Init

  USE Parallelism, ONLY:   &
       myId,               &
       myId_four,          &
       mygroup_four,       &
       CreateFourierGroup, &
       maxnodes_four,      &
       maxNodes

  USE Sizes, ONLY : ibMax, jbMax,kMax, kmaxloc_out, Kmaxloc_In, &
       iMax, jMax, jMaxHalf, jMinPerM, iMaxPerJ, lmaxloc, &
       mMax, nExtMax, mnMax, mnExtMax, mnExtMap, &
       mymnMax, mymnExtMax, kmaxloc, jbMax_ext, &
       RegisterBasicSizes, RegisterOtherSizes, &
       ngroups_four,npperg_four,first_proc_four,map_four,nprocmax_four, &
       nrecs_gr, nsends_gr, nrecs_f, nrecs_g, nrecs_diag, nsends_diag, & 
       SpectralDomainDecomp, GridDecomposition, VerticalGroups, &
       ReshapeVerticalGroups

  USE Utils, ONLY: CreateGaussQuad, &
       GaussPoints, GaussWeights, CreateAssocLegFunc, &
       LegFuncS2F, CreateGridValues, DestroyAssocLegFunc, allpolynomials

  USE SpecDynamics, ONLY: InitDZtoUV, InitUvtodz, InitGozrim

  USE Transform, ONLY: InitTransform

  USE Communications, ONLY: Set_Communic_buffer,exchange_ftog,Exchange_diag

  USE Constants, ONLY: givenfouriergroups, nproc_vert, ibdim_size,nfprt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitAll


CONTAINS


  SUBROUTINE InitAll (trunc, nlon, nlat, vert_in, vert_out)
    INTEGER, INTENT(IN)           :: trunc
    INTEGER, INTENT(IN)           :: nlon
    INTEGER, INTENT(IN)           :: nlat
    INTEGER, INTENT(IN)           :: vert_in
    INTEGER, INTENT(IN)           :: vert_out

    CHARACTER(LEN=*), PARAMETER :: h="**(InitAll)**"
    INTEGER :: wavesFactor
    INTEGER :: j
    INTEGER, ALLOCATABLE :: mPerLat(:)
    INTEGER, ALLOCATABLE :: iMaxPerLat(:)

    ! verify consistency and register basic sizes

    CALL RegisterBasicSizes(trunc, nLat, nLon, vert_in, vert_out)

    ! Initialize GaussQuad (Gaussian points and weights)

    CALL CreateGaussQuad(jMax)

    ! Spectral Domain Decomposition

    CALL VerticalGroups(givenfouriergroups,nproc_vert,vert_out,lMaxloc)
    CALL ReshapeVerticalGroups(kmax,kMaxloc)
    CALL SpectralDomainDecomp()

    ! Initialize AssocLegFunc (computes associated legendre functions)

    allpolynomials = .false.
    CALL CreateAssocLegFunc(allpolynomials)

    ! computes mPerLat

    ALLOCATE (mPerLat(jMax))
    mPerLat = mMax

    ! compute iMaxPerLat

    ALLOCATE(iMaxPerLat(jMax))
    iMaxPerLat=iMax

    ! finish building sizes

    CALL RegisterOtherSizes(iMaxPerLat, mPerLat, myid, maxnodes)
    CALL CreateFourierGroup(mygroup_four,myid_four)
    ibmax = ibdim_size
    CALL GridDecomposition(ibmax,jbmax,maxnodes,myid)
    CALL Exchange_ftog(nrecs_f,nrecs_g)
    CALL Exchange_diag(nrecs_diag,nsends_diag)

    ! for now on, all global constants defined at sizes can be used

    CALL CreateGridValues

    ! initialize dztouv and uvtodz

    CALL InitDZtoUV()
    CALL InitUvtodz()
    CALL InitGozrim()

    CALL Set_Communic_buffer()

    CALL InitTransform()

    ! deallocate legendre functions (already stored in transform,
    ! epslon already used in initializations)

    CALL DestroyAssocLegFunc()

  END SUBROUTINE InitAll

END MODULE Init

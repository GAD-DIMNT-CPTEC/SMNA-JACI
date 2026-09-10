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
       iMax, jMax, jMaxHalf, jMinPerM, iMaxPerJ, &
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

  USE SpecDynamics, ONLY: InitDZtoUV, InitUvtodz

  USE Transform, ONLY: InitTransform

  USE Communications, ONLY: Set_Communic_buffer,exchange_ftog,Exchange_diag

  USE InputParameters, ONLY: givenfouriergroups, nproc_vert, ibdim_size,nfprt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitAll


CONTAINS


  SUBROUTINE InitAll (trunc_in, trunc_out, vert, nlon, vert_out)
    INTEGER, INTENT(IN)           :: trunc_in
    INTEGER, INTENT(IN)           :: trunc_out
    INTEGER, INTENT(IN)           :: vert
    INTEGER, INTENT(IN)           :: vert_out
    INTEGER, INTENT(IN)           :: nlon

    CHARACTER(LEN=*), PARAMETER :: h="**(InitAll)**"
    INTEGER :: wavesFactor
    INTEGER :: nLat
    INTEGER :: j
    INTEGER, ALLOCATABLE :: mPerLat(:)
    INTEGER, ALLOCATABLE :: iMaxPerLat(:)

    nLat = nLon/2

    ! verify consistency and register basic sizes

    CALL RegisterBasicSizes(trunc_in, trunc_out, nLat, nLon, vert)
    ! Initialize GaussQuad (Gaussian points and weights)

    CALL CreateGaussQuad(jMax)

    ! Spectral Domain Decomposition

    CALL VerticalGroups(givenfouriergroups,nproc_vert,vert_out,kMaxloc_out)

    CALL ReshapeVerticalGroups(kmax,kMaxloc_in)

    kmaxloc = kMaxloc_in
    CALL SpectralDomainDecomp()

    ! Initialize AssocLegFunc (computes associated legendre functions)

    allpolynomials = .false.
    CALL CreateAssocLegFunc(allpolynomials)

    ! computes mPerLat
    IF (ALLOCATED(mPerLat))THEN
       DEALLOCATE (mPerLat)
       ALLOCATE (mPerLat(jMax))
    ELSE
       ALLOCATE (mPerLat(jMax))
    END IF

    mPerLat = mMax

    ! compute iMaxPerLat

    IF (ALLOCATED(iMaxPerLat))THEN
       DEALLOCATE (iMaxPerLat)
       ALLOCATE(iMaxPerLat(jMax))
    ELSE
       ALLOCATE(iMaxPerLat(jMax))
    END IF

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

    CALL Set_Communic_buffer()

    CALL InitTransform()

    ! deallocate legendre functions (already stored in transform,
    ! epslon already used in initializations)

    CALL DestroyAssocLegFunc()

  END SUBROUTINE InitAll

END MODULE Init

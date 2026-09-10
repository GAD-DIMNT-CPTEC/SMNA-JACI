!
!  $Author: pkubota $
!  $Date: 2010/04/20 20:18:04 $
!  $Revision: 1.16 $
!
MODULE Sizes   ! Version 0 of Nov 25th, 2001



  ! SINGLE REPOSITORY OF SIZES AND MAPPINGS OF GRID, FOURIER AND SPECTRAL
  ! FIELD REPRESENTATIONS. SHOULD BE USED BY ALL MODEL MODULES THAT REQUIRE
  ! SUCH INFORMATION.
  !
  ! Module exports three procedures:
  ! RegisterBasicSizes: set all repository values for spectral representation.
  !                     It also sets latitudes and maximum number of longitudes
  !                     for grid representation.
  ! RegisterOtherSizes: set remaining repository values.
  ! DumpSizes:          detailed dumping of repository values
  !
  ! Module exports a set of values and arrays, described bellow.
  !
  ! Notation used throughout the code:
  ! m  Legendre Order (also Fourier wave number). Goes from 1 to mMax.
  ! n  Legendre Degree. Goes from m to nMax on regular spectral fields
  !    and from m to nExtMax on extended spectral fields.
  ! mn index to store the spectral triangle (m,n) into a single dimension.
  ! k  vertical index, from 1 to kMax
  ! j  latitude index, from 1 to jMax (full field)
  !    or from 1 to jMaxHalf (hemisphere)
  ! i  longitude index, from 1 to iMax (gaussian grid) or
  !    from 1 to iMaxPerJ(j) (at latitude j of the reduced grid)
  ! ib index of a longitude on a block of longitudes packed in one
  !    dimension of grids;
  ! jb index for which block of longitudes to take, packed in another
  !    dimension of grids

  ! SPECTRAL REPRESENTATION:
  ! A regular spectral field should be dimensioned (2*mnMax,kMax), where
  ! the first dimension accomodates pairs of (real,imaginary) spectral
  ! coefficients for a fixed vertical and the second dimension varies
  ! verticals.
  ! Extended spectral fields should use (2*mnExtMax,kMax).
  ! This module provides mapping functions to map (m,n) into mn:
  ! For a regular spectral field with complex entries, Legendre Order m
  ! and Degree n has real coefficient at position 2*mnMap(m,n)-1 and
  ! imaginary coefficient at position 2*mnMap(m,n). Inverse mappings 
  ! (from mn to m and n) are also provided. Maps for the extended spectral 
  ! field are also provided.

  ! GRID REPRESENTATION:
  ! Since the number of longitudes per latitude may vary with latitude
  ! (on reduced grids), sets of latitudes (with all longitudes) 
  ! are packed together in the first dimension of Grids. 
  ! Near the pole, many latitudes are packed; near the equator, just a few. 
  ! Second dimension is vertical.
  ! Third dimension is the number of packed latitudes required to represent
  ! a full field.
  ! A grid field should be dimensioned (ibMax, kMax, jbMax).
  ! This module provides mapping functions to map latitude j and longitude
  ! i into (ib,jb).
  ! Map ibPerIJ(i,j) gives index of first dimension that stores longitude
  ! i of latitude j. Map jbPerIJ(i,j) gives index of third dimension that
  ! stores longitude i of latitude j. Consequently, the point of longitude
  ! i, latitude j and vertical k is stored at (ibPerIJ(i,j), k, jbPerIJ(i,j)).
  ! Inverse mappings (from (ib,jb) to (i,j)) are also provided.

  ! FOURIER REPRESENTATION:
  ! For the moment, Fourier fields are represented externally to the
  ! transform. That should not happen in the future.
  ! First dimension contains pairs of (real,imaginary) fourier 
  ! coefficients. Second dimension is latitude. Third dimension is
  ! vertical. A full fourier field should be dimensioned 
  ! (iMax+1, jMax, kMax)

  USE InputParameters, ONLY: &
       r8,nfprt

  USE Parallelism, ONLY:&
       MsgDump,         &
       MsgOne,          &
       FatalError,      &
       mygroup_four,    &
       maxnodes,        &
       myId,            &
       maxnodes_four,   &
       myId_four,       &
       unitDump

  IMPLICIT NONE


  ! INTERNAL DATA STRUCTURE:
  ! Provision for MPI computation:
  !
  ! Domain Decomposition:
  ! There are maxNodes MPI processes, numbered from 0 to maxNodes - 1
  ! The number of this MPI process is myId, 0 <= myId <= maxNodes - 1
  ! SPECTRAL REPRESENTATION:
  ! mMax is the maximum Legendre Order (also maximum Fourier wave number).
  !      it is set to model truncation + 1 (Legendre Order 0 is m=1)
  ! nMax is the maximum Legendre Degree for regular spectral fields.
  !      it is set to model truncation + 1 (Legendre Degree 1 is n=1)
  ! mnMax is the amount of spectral coeficients per vertical, for 
  !       regular spectral fields. It is the number of points in the
  !       regular triangular spectral plane.
  ! mMap indexed by mn in [1:mnMax]; returns m at this mn for regular
  !      spectral fields;
  ! nMap indexed by mn in [1:mnMax]; returns n at this mn for regular
  !      spectral fields;
  ! mnMap indexed by (m,n); returns mn that stores spectral coefficient
  !       (m,n), when spectral coefficients are real (e.g. Associated 
  !       Legendre Functions). For complex spectral coefficients, 
  !       mn=2*mnMap(m,n) for the imaginary part and mn=2*mnMap(m,n)-1
  !       for the real part.
  ! Remaining variables with Ext extension have the same meaning, but
  ! for the extended spectral space (where nExtMax=trunc+2), that is,
  ! where each Legendre Order has one more Legendre Degree than usual.

  ! Spectral coefficients are spread across processes.
  ! The set of values of m is partitioned among processes.
  ! Consequently, all values of n for some values of m are stored
  ! at each process. Each value of m is stored at a single process.
  !
  ! Array nodeHasM(m) has which MPI process stores m, 1 <= m <= mMax
  !
  ! Variable myMMax (0 <= myMMax <= mMax) has how many m's are
  ! stored at this node.
  !
  ! Array lm2m(l), l=1,...,myMMax has the m's stored at this
  ! node. That is,
  ! it maps local m (lm) into global m (m).


  INTEGER              :: mMax=-1
  INTEGER              :: mMax_in=-1
  INTEGER              :: nMax=-1
  INTEGER              :: mnMax=-1
  INTEGER              :: mnMax_out=-1
  INTEGER, ALLOCATABLE :: mMap(:)
  INTEGER, ALLOCATABLE :: nMap(:)
  INTEGER, ALLOCATABLE :: mnMap(:,:)
  INTEGER, ALLOCATABLE :: mnMap_out(:,:)
  INTEGER              :: nExtMax=-1
  INTEGER              :: mnExtMax=-1
  INTEGER, ALLOCATABLE :: mExtMap(:)
  INTEGER, ALLOCATABLE :: nExtMap(:)
  INTEGER, ALLOCATABLE :: mnExtMap(:,:)

  INTEGER, ALLOCATABLE :: nodeHasM(:,:)
  INTEGER, ALLOCATABLE :: lm2m(:)
  INTEGER, ALLOCATABLE :: msInProc(:,:)
  INTEGER, ALLOCATABLE :: msPerProc(:)

  INTEGER              :: myMMax
  INTEGER              :: mMaxLocal
  INTEGER              :: mnMaxLocal
  INTEGER              :: mnExtMaxLocal
  INTEGER              :: myMNMax
  INTEGER              :: myMNExtMax
  LOGICAL              :: HaveM1

  INTEGER, ALLOCATABLE :: myMMap(:)
  INTEGER, ALLOCATABLE :: myNMap(:)
  INTEGER, ALLOCATABLE :: mnsPerProc(:)
  INTEGER, ALLOCATABLE :: mnsExtPerProc(:)
  INTEGER, ALLOCATABLE :: myMNMap(:,:)
  INTEGER, ALLOCATABLE :: myMExtMap(:)
  INTEGER, ALLOCATABLE :: myNExtMap(:)
  INTEGER, ALLOCATABLE :: myMNExtMap(:,:)

  !
  !  Spectral representation and division to be used in the semi-implicit part
  !  -------------------------------------------------------------------------

  !
  !  Spectral communicators to spread surface field
  !  ----------------------------------------------

  INTEGER              :: ncomm_spread 
  INTEGER, ALLOCATABLE :: comm_spread(:,:)
  INTEGER, ALLOCATABLE :: ms_spread(:,:)

  ! LATITUDES REPRESENTATION:
  ! jMax is the number of latitudes for full Fourier and Grid representations.
  ! jMaxHalf is the number of latitudes at each hemisphere.
  ! mMaxPerJ is an array indexed by latitudes that stores the maximum value
  !          of m at each latitude. Consequently, the contribution of 
  !          latitude j for the Legendre Transform should be taken only
  !          up to Legendre Order mMaxPerJ(j). By the same token, 
  !          FFTs at latitude j should be computed only up to 
  !          Fourier wave number mMaxPerJ(j).
  !          For regular grids, mMaxPerJ(j) = mMax for all j.
  ! jMinPerM is the inverse mapping - an array indexed by Legendre Orders (m)
  !          containing the smallest latitude (value of j) that has that
  !          order. Latitudes that contain Legendre Order (and Fourier
  !          wave number) m are jMinPerM(m) <= j <= jMax-jMinPerM(m)+1
  !          For regular grids, jMinPerM(m) = 1 for all m.

  


  INTEGER              :: jMax=-1
  INTEGER              :: jMaxHalf=-1
  INTEGER, ALLOCATABLE :: mMaxPerJ(:)
  INTEGER, ALLOCATABLE :: jMinPerM(:)
  INTEGER, ALLOCATABLE :: jMaxPerM(:)



  ! LONGITUDES REPRESENTATION:
  ! iMax is the maximum number of longitudes per latitude;
  !      it is the number of longitudes per latitude for all latitudes
  !      on regular grids;
  !      it is only the maximum number of longitudes per latitude on 
  !      reduced grids; it is the actual number of longitudes per latitude
  !      close to the equator on reduced grids;
  ! ijMax is the number of horizontal grid points at regular or reduced
  !       grids.
  ! iMaxPerJ is the actual number of longitudes per latitude on regular
  !          and reduced grids; latitude j has iMaxPerJ(j) longitudes.



  INTEGER              :: iMax=-1
  INTEGER              :: ijMax=-1
  INTEGER              :: ijMaxGauQua=-1
  INTEGER, ALLOCATABLE :: iMaxPerJ(:)



  ! GRID REPRESENTATION:
  ! All longitudes of a set of latitudes are packed together in the
  ! first dimension of grids. That decreases the waste of memory when
  ! comparing to store one latitude only, for the case of reduced
  ! grid. 
  ! ibMax is the maximum number of longitudes packed into the first dimension
  !       of grids. The actual number of longitudes vary with the third
  !       dimension of the grid representation.
  ! jbMax is the number of sets of longitudes required to store an
  !       entire field.
  ! ibPerIJ maps longitude i and latitude j into the first dimension
  !         of grid representation. It is indexed ibPerIJ(i,j).
  ! jbPerIJ maps longitude i and latitude j into the third dimension
  !         of grid representation. It is indexed jbPerIJ(i,j).
  ! iPerIJB gives which longitude is stored at first dimension index
  !         i and third dimension index j of Grid representations.
  !         It is indexed iPerIJB(i,j)
  ! jPerIJB gives which latitude is stored at first dimension index
  !         i and third dimension index j of Grid representations.
  !         It is indexed jPerIJB(i,j)
  ! ibMaxPerJB gives how many latitudes are actually stored at third
  !            dimension jb. Since the number of longitudes vary with
  !            latitudes, the amount of space actually used in the first
  !            dimension of grid representations vary with the third
  !            dimension. Array ibMaxPerJB, indexed by jb, accounts for
  !            such variation.

  ! Grid Point Decomposition
  !
  ! Blocks of surface points are spreaded across processes.
  ! Each process has all longitudes of a set of latitudes (block)
  !

  !
  !  decomposition in Fourier space
  !
  
  INTEGER, ALLOCATABLE :: ibMaxPerJB(:)
  INTEGER, ALLOCATABLE :: firstlatinproc_f(:)
  INTEGER, ALLOCATABLE :: lastlatinproc_f(:)
  INTEGER, ALLOCATABLE :: nlatsinproc_f(:)
  INTEGER, ALLOCATABLE :: nodeHasJ_f(:)
  INTEGER, ALLOCATABLE :: kfirst_four(:)
  INTEGER, ALLOCATABLE :: klast_four(:)
  INTEGER, ALLOCATABLE :: nlevperg_four(:)
  INTEGER, ALLOCATABLE :: npperg_four(:)
  INTEGER, ALLOCATABLE :: map_four(:,:)
  INTEGER, ALLOCATABLE :: first_proc_four(:)
  INTEGER, ALLOCATABLE :: nlatsinproc_d(:)
  INTEGER, ALLOCATABLE :: messages_f(:,:)
  INTEGER, ALLOCATABLE :: messages_g(:,:)
  INTEGER, ALLOCATABLE :: messproc_f(:,:)
  INTEGER, ALLOCATABLE :: messproc_g(:,:)
  INTEGER, ALLOCATABLE :: nodeHasJ(:)
  INTEGER, ALLOCATABLE :: nset(:)

  INTEGER              :: myfirstlat_f
  INTEGER              :: mylastlat_f
  INTEGER              :: myfirstlev
  INTEGER              :: mylastlev
  INTEGER              :: myJMax_f
  INTEGER              :: ngroups_four
  INTEGER              :: nprocmax_four
  INTEGER              :: JMaxlocal_f
  INTEGER              :: kMaxloc 
  INTEGER              :: kMaxloc_in 
  INTEGER              :: kMaxloc_out
  INTEGER              :: nrecs_f 
  INTEGER              :: nrecs_g 
  LOGICAL              :: havesurf

  !
  ! grid decomposition 
  !
  
  INTEGER              :: ibMax=-1
  INTEGER              :: jbMax=-1
  INTEGER              :: jbMax_ext
  INTEGER              :: myfirstlat
  INTEGER              :: mylastlat
  INTEGER              :: nrecs_gr 
  INTEGER              :: nsends_gr
  INTEGER              :: nrecs_diag
  INTEGER              :: nsends_diag
  INTEGER              :: myfirstlat_diag
  INTEGER              :: mylastlat_diag
  INTEGER              :: myJMax_d
  INTEGER              :: jovlap
  INTEGER, ALLOCATABLE :: firstandlastlat(:,:)
  INTEGER, ALLOCATABLE :: myfirstlon(:)
  INTEGER, ALLOCATABLE :: mylastlon(:)
  INTEGER, ALLOCATABLE :: firstlon(:,:)
  INTEGER, ALLOCATABLE :: lastlon(:,:)
  INTEGER, ALLOCATABLE :: firstlat(:)
  INTEGER, ALLOCATABLE :: lastlat(:)
  INTEGER, ALLOCATABLE :: ibPerIJ(:,:)
  INTEGER, ALLOCATABLE :: jbPerIJ(:,:)
  INTEGER, ALLOCATABLE :: iPerIJB(:,:)
  INTEGER, ALLOCATABLE :: jPerIJB(:,:)
  INTEGER, ALLOCATABLE :: pointsinproc(:)
  INTEGER, ALLOCATABLE :: myrecsgr(:,:)
  INTEGER, ALLOCATABLE :: myrecspr(:,:)
  INTEGER, ALLOCATABLE :: mysendsgr(:,:)
  INTEGER, ALLOCATABLE :: mysendspr(:,:)
  INTEGER, ALLOCATABLE :: myrecs_diag(:,:)
  INTEGER, ALLOCATABLE :: myrecspr_diag(:,:)
  INTEGER, ALLOCATABLE :: mysends_diag(:,:)
  INTEGER, ALLOCATABLE :: mysendspr_diag(:,:)
  INTEGER, ALLOCATABLE :: gridmap(:,:)

  INTEGER              :: kMax=-1

  REAL(KIND=r8), ALLOCATABLE :: ci(:)      ! 1 - sigma each level (level 1 at surface)
  REAL(KIND=r8), ALLOCATABLE :: si(:)      ! sigma
  REAL(KIND=r8), ALLOCATABLE :: del(:)     ! layer thickness (in sigma)
  REAL(KIND=r8), ALLOCATABLE :: delcl(:)     ! layer thickness (in cl)
  REAL(KIND=r8), ALLOCATABLE :: rdel2(:)
  REAL(KIND=r8), ALLOCATABLE :: sl(:)      ! sigma at layer midpoint
  REAL(KIND=r8), ALLOCATABLE :: cl(:)      ! 1.0 - sl
  REAL(KIND=r8), ALLOCATABLE :: rpi(:)     ! 'pi' ratios at adjacent layers
  LOGICAL, PARAMETER, PRIVATE :: dumpLocal=.FALSE.

CONTAINS



  SUBROUTINE RegisterBasicSizes(trunc_in, trunc_out, nLat, nLon, vert)
    INTEGER, INTENT(IN) :: trunc_in
    INTEGER, INTENT(IN) :: trunc_out
    INTEGER, INTENT(IN) :: nLat
    INTEGER, INTENT(IN) :: nLon
    INTEGER, INTENT(IN) :: vert
    INTEGER :: m, n, mn
    CHARACTER(LEN=*), PARAMETER :: h="**(RegisterBasicSizes)**"

    jMax     = nLat
    jMaxHalf = nLat/2

    iMax     = nLon

    kMax     = vert

    mMax_in  = trunc_in + 1
    mMax     = trunc_out + 1
    nMax     = mMax
    nExtMax  = mMax + 1
    mnExtMax = (nExtMax+2)*(nExtMax-1)/2
    ALLOCATE (mExtMap(mnExtMax))
    ALLOCATE (nExtMap(mnExtMax))
    ALLOCATE (mnExtMap(mMax,nExtMax))
    nExtMap = -1  ! flag mapping error
    mExtMap = -1  ! flag mapping error
    mnExtMap = -1 ! flag mapping error
    mn = 0
    DO m = 1, mMax
       DO n = m, mMax+1
          mn = mn + 1
          mnExtMap(m,n) = mn
          mExtMap(mn)   = m
          nExtMap(mn)   = n
       END DO
    END DO
    mnMax = (mmax_in * (mMax_in+1))/2
    mnMax_out = (mmax * (mMax+1))/2
    if (myid.eq.0) write(*,*) ' mnMax ',mnMax,' mnMax_out ',mnMax_out
    ALLOCATE (mnMap(mMax_in,mMax_in))
    ALLOCATE (mnMap_out(mMax,mMax))
    ALLOCATE (mMap(mnMax))
    ALLOCATE (nMap(mnMax))
    mnMap = -1  ! flag mapping error
    mnMap_out = -1  ! flag mapping error
    mMap = -1   ! flag mapping error
    nMap = -1   ! flag mapping error
    mn = 0
    DO m = 1, mMax_in
       DO n = m, mMax_in
          mn = mn + 1
          mnMap(m,n) = mn
          mMap(mn)   = m
          nMap(mn)   = n
       END DO
    END DO
    mn = 0
    DO m = 1, mMax
       DO n = m, mMax
          mn = mn + 1
          mnMap_out(m,n) = mn
       END DO
    END DO
    ijMaxGauQua = iMax*jMax

  END SUBROUTINE RegisterBasicSizes


  SUBROUTINE VerticalGroups(givenfouriergroups,nproc_vert,kmax_in, &
                            kmaxloc_io)
    LOGICAL, INTENT(IN) :: givenfouriergroups
    INTEGER, INTENT(IN) :: nproc_vert
    INTEGER, INTENT(IN) :: kmax_in
    INTEGER, INTENT(OUT):: kmaxloc_io
    CHARACTER(LEN=*), PARAMETER :: h="**(RegisterOtherSizes)**"
    INTEGER :: i, n, nlp, nrest, ipar, ng, npsq, nl, nn, np, nlg
    INTEGER :: next, k, m, nprest, l
    !
    !  Compute vertical decomposition of fourier groups
    !
    IF (givenfouriergroups) THEN
       ng = nproc_vert
     ELSE
       npsq = SQRT(REAL(maxnodes))
       DO nl=1,kmax_in
          nn = nl*maxnodes/kmax_in
          IF(nn.ge.npsq) EXIT
       ENDDO
       ng = kmax_in / nl
       IF (nl*ng.lt.kmax_in) ng = ng + 1
    ENDIF
    ngroups_four = ng
    nl = kmax_in / ng
    IF (nl*ng.lt.kmax_in) nl = nl + 1
    nrest = nl * ng - kmax_in
    ALLOCATE(kfirst_four(0:maxnodes))
    ALLOCATE(klast_four(0:maxnodes))
    ALLOCATE(npperg_four(ng))
    ALLOCATE(nlevperg_four(ng))
    ALLOCATE(first_proc_four(ng))
    np = 0
    DO i=1,ng
       IF (i.le.nrest) THEN 
          nlg = nl - 1
        ELSE
          nlg = nl
       ENDIF
       nlevperg_four(i) = nlg
       npperg_four(i) = nlg * maxnodes / kmax_in
       np = np + npperg_four(i)
    ENDDO
    nprest = maxnodes - np
    next = nprest / ng
    nprest = nprest - ng * next
    DO i=1,ng
       IF (i.le.nprest) THEN
          npperg_four(i) = npperg_four(i) + next + 1
        ELSE
          npperg_four(i) = npperg_four(i) + next
       ENDIF
    ENDDO
    ALLOCATE(nset(ng))
    nprocmax_four = MAXVAL(npperg_four)
    ALLOCATE(map_four(ng,0:nprocmax_four-1))
    n = 0
    nset = 0
    DO
       DO i=1,ngroups_four
          IF (nset(i).lt.npperg_four(i)) THEN
             IF (myid.eq.n) THEN
                mygroup_four = i
                maxnodes_four = npperg_four(i)
                myid_four = nset(i)
             ENDIF
             map_four(i,nset(i)) = n
             nset(i) = nset(i) + 1
             n = n + 1
          ENDIF
       END DO
       IF (n.ge.maxnodes) EXIT
    END DO
    m = 0
    DO i=1,ngroups_four
       DO k=0,npperg_four(i)-1
          l = map_four(i,k)
          kfirst_four(l) = m + 1
          klast_four(l) = m + nlevperg_four(i)
       ENDDO
       m = m + nlevperg_four(i)
    ENDDO
    DO i=1,ngroups_four
       first_proc_four(i) = i - 1
    ENDDO
    kmaxloc_io = klast_four(myid) - kfirst_four(myid) + 1
    myfirstlev = kfirst_four(myid)
    mylastlev = klast_four(myid)

  END SUBROUTINE VerticalGroups

  SUBROUTINE ReshapeVerticalGroups(kmax_in, kmaxloc_io)
    INTEGER, INTENT(IN) :: kmax_in
    INTEGER, INTENT(OUT):: kmaxloc_io
    CHARACTER(LEN=*), PARAMETER :: h="**(RegisterOtherSizes)**"
    INTEGER :: i, n, nlp, nrest, ipar, ng, npsq, nl, nn, np, nlg
    INTEGER :: next, k, m, nprest, l
    !
    !  Compute vertical decomposition of fourier groups
    !
    ng = ngroups_four
    nl = kmax_in / maxnodes
    np = 0
    DO i=1,ng
       np = np + npperg_four(i) * nl
       nlevperg_four(i) = npperg_four(i) * nl
    ENDDO
    nrest = kmax_in - np
    DO
       DO i=ng,1,-1
          np = min (1,nrest)
          nlevperg_four(i) = nlevperg_four(i) + np
          nrest = nrest - np
       ENDDO
       IF (nrest.eq.0) EXIT
    ENDDO
    IF (nlevperg_four(1).eq.0) THEN
       nlevperg_four(1) = 1
       nlevperg_four(ng) = nlevperg_four(ng) - 1
    ENDIF
    m = 0
    DO i=1,ngroups_four
       DO k=0,npperg_four(i)-1
          l = map_four(i,k)
          kfirst_four(l) = m + 1
          klast_four(l) = m + nlevperg_four(i)
       ENDDO
       m = m + nlevperg_four(i)
    ENDDO
    kmaxloc_io = klast_four(myid) - kfirst_four(myid) + 1
    myfirstlev = kfirst_four(myid)
    mylastlev = klast_four(myid)

  END SUBROUTINE ReshapeVerticalGroups



  SUBROUTINE RegisterOtherSizes(iMaxPerLat, mPerLat, myid, maxnodes)
    INTEGER, INTENT(IN) :: iMaxPerLat(jMax)
    INTEGER, INTENT(IN) :: mPerLat(jMax)
    INTEGER, INTENT(IN) :: myid
    INTEGER, INTENT(IN) :: maxnodes
    CHARACTER(LEN=*), PARAMETER :: h="**(RegisterOtherSizes)**"
    INTEGER :: MinLatPerBlk
    !$ INTEGER, EXTERNAL ::  OMP_GET_MAX_THREADS
    INTEGER :: i, j, k, l, jk, m, ib, jb, jh, cnt, imp2
    INTEGER :: meanl, meanp, jlast, jused, jfirst
    INTEGER :: mfirst(maxnodes_four), mlast(maxnodes_four)
    INTEGER :: npoints(maxnodes_four)
    INTEGER :: mlast1(maxnodes_four), npoints1(maxnodes_four)
    INTEGER :: maxpoints, maxpointsold, nproc, nlp, nrest
    LOGICAL :: improved, done

    ALLOCATE (mMaxPerJ(jMax))
    mMaxPerJ = mPerLat
    ALLOCATE (iMaxPerJ(jMax))
    iMaxPerJ = iMaxPerLat
    IF (iMax < MAXVAL(iMaxPerJ)) THEN
       STOP ' imax and imaxperj disagree'
    END IF
    ijMax = SUM(iMaxPerJ)

    ALLOCATE (firstlatinproc_f(0:maxNodes_four-1))
    ALLOCATE (lastlatinproc_f(0:maxNodes_four-1))
    ALLOCATE (nlatsinproc_f(0:maxNodes_four-1))
    ALLOCATE (nodeHasJ_f(jMax))

    nproc = maxnodes_four
       
       nlp = jmax / nproc
       nrest = jmax - nproc * nlp
       DO i=1,nproc-nrest
          mfirst(i) = (i-1)*nlp + 1
          mlast(i) = mfirst(i) + nlp - 1
          npoints(i) = imaxperj(1)*nlp
       ENDDO
       DO i=nproc-nrest+1,nproc
          mfirst(i) = mlast(i-1) + 1
          mlast(i) = mfirst(i) + nlp
          npoints(i) = imaxperj(1)*(nlp+1)
       ENDDO

    IF (ANY(npoints(:) <= 0)) THEN
       CALL FatalError(h//" Too many MPI processes; "//&
            &"there are processes with 0 latitudes")
       STOP
    END IF

    myfirstlat_f = mfirst(myid_four+1)
    mylastlat_f = mlast(myid_four+1)
    jMaxlocal_f = 0
    nodeHasJ_f=-1
    DO k=0,maxNodes_four-1
       firstlatinproc_f(k) = mfirst(k+1)
       lastlatinproc_f(k) = mlast(k+1)
       nlatsinproc_f(k) = mlast(k+1)-mfirst(k+1)+1
       DO j=firstlatinproc_f(k),lastlatinproc_f(k)
          nodeHasJ_f(j) = k
       ENDDO
    ENDDO
    jMaxlocal_f = MAXVAL(nlatsinproc_f)
    myJMax_f = mylastlat_f - myfirstlat_f + 1

    ALLOCATE (jMinPerM(mMax))
    ALLOCATE (jMaxPerM(mMax))
    jMinPerM = jMaxHalf
    DO j = 1, jMaxHalf
       m = mMaxPerJ(j)
       jMinPerM(1:m) = MIN(j, jMinPerM(1:m))
    END DO
    jMaxPerM = jMax - jMinPerM + 1

  END SUBROUTINE RegisterOtherSizes

  SUBROUTINE SpectralDomainDecomp()

    INTEGER :: m
    INTEGER :: n
    INTEGER :: mn
    INTEGER :: mBase
    INTEGER :: mMid
    INTEGER :: lm
    INTEGER :: mm
    INTEGER :: mng
    INTEGER :: k 
    INTEGER :: rest
    INTEGER :: is
    INTEGER :: i
    INTEGER :: MaxN, ns, ip, inc, mngiv, np, npl
    INTEGER :: ierr
    INTEGER, ALLOCATABLE :: ini(:), sends(:)
    CHARACTER(LEN=8) :: c0
    CHARACTER(LEN=*), PARAMETER :: h="**(SpectralDomainDecomp)**"



    ! DOMAIN DECOMPOSITION OF FOURIER WAVE NUMBERS (m's)



    ! any MPI process has at most mMaxLocal m's (may have less)
    ! mMaxLocal is used to dimension arrays over all MPI processes

    mMaxLocal = mMax / maxNodes_four
    IF (mMax-mMaxLocal*maxNodes_four.ne.0) mMaxLocal = mMaxLocal + 1


    ! msPerProc: how many m's at each MPI process

    ALLOCATE(msPerProc(0:maxNodes_four-1))

    ! msInProc: which m's are at each MPI process
    ! note that indexing is restricted to
    ! msInProc(1:msPerProc(pId), pId), with pId=0,maxNodes-1

    ALLOCATE(msInProc(mMaxLocal,0:maxNodes_four-1))

    ! nodeHasM: which process has a particular value of m.
    ! Values of nodeHasM for p processes are:
    ! (0, 1, 2, ..., p-1, p-1, ..., 2, 1, 0, 0, 1, 2, ...)
    ! This distribution tries to spread evenly load across
    ! processes, leaving uneven load to the smaller m's.

    ALLOCATE(nodeHasM(mMax,ngroups_four))

    ! domain decomposition of m's

    DO i=1,ngroups_four
       maxN = npperg_four(i)
       mm = 1
       DO mBase = 1, mMax, 2*maxN
          mMid = mBase+maxN-1
          DO m = mBase, MIN(mMid, mMax)
             nodeHasM(m,i) = m - mBase
          END DO
          mm = mm + 1
          DO m = mMid+1, MIN(mMid+maxN, mMax)
             nodeHasM(m,i) = maxN + mMid - m
          END DO
          mm = mm + 1
       END DO
    END DO
    IF (nodeHasM(1,mygroup_four)==myid_four) THEN
       haveM1 = .TRUE.
     ELSE
       haveM1 = .FALSE.
    ENDIF
    havesurf = .FALSE.
    IF (mygroup_four.eq.1) havesurf = .TRUE.
    mm = 1
    msPerProc = 0
    i = mygroup_four
    DO mBase = 1, mMax, 2*maxNodes_four
       mMid = mBase+maxNodes_four-1
       DO m = mBase, MIN(mMid, mMax)
          msPerProc(nodeHasM(m,i)) = msPerProc(nodeHasM(m,i)) + 1
          msInProc(mm,nodeHasM(m,i)) = m
       END DO
       mm = mm + 1
       DO m = mMid+1, MIN(mMid+maxNodes_four, mMax)
          msPerProc(nodeHasM(m,i)) = msPerProc(nodeHasM(m,i)) + 1
          msInProc(mm,nodeHasM(m,i)) = m
       END DO
       mm = mm + 1
    END DO

    ! current parallelism restricts the number of MPI processes 
    ! to truncation + 1

    IF (ANY(msPerProc <= 0)) THEN
       CALL FatalError(h//" Too many MPI processes; "//&
            &"there are processes with 0 Fourier waves")
       STOP
    END IF

    ! myMMax: scalar containing how many m's at this MPI process

    myMMax = msPerProc(myId_four)

    ! lm2m: maps local indexing of m to global indexing of m,
    ! that is, maps (1:myMMax) to (1:mMax)

    ALLOCATE(lm2m(myMMax))
    lm2m(1:mymmax) = msInProc(1:mymmax,myid_four)




    ! DOMAIN DECOMPOSITION OF SPECTRAL COEFFICIENTS (mn's)
    ! all mn's of a single m belongs to  a unique process




    ! mnsPerProc: how many mn's at each MPI process

    ALLOCATE(mnsPerProc(0:maxNodes_four-1))

    ! mnsExtPerProc: how many mnExt's at each MPI process

    ALLOCATE(mnsExtPerProc(0:maxNodes_four-1))

    ! domain decomposition of mn's and mnExt's

    mm = 1
    mnsPerProc = 0
    mnsExtPerProc = 0
    i = mygroup_four
    DO mBase = 1, mMax, 2*maxNodes_four
       mMid = mBase+maxNodes_four-1
       DO m = mBase, MIN(mMid, mMax)
          mnsPerProc(nodeHasM(m,i)) =  mnsPerProc(nodeHasM(m,i)) + mmax - m + 1
          mnsExtPerProc(nodeHasM(m,i)) =  mnsExtPerProc(nodeHasM(m,i)) + mmax - m + 2
       END DO
       mm = mm + 1
       DO m = mMid+1, MIN(mMid+maxNodes_four, mMax)
          mnsPerProc(nodeHasM(m,i)) =  mnsPerProc(nodeHasM(m,i)) + mmax - m + 1
          mnsExtPerProc(nodeHasM(m,i)) =  mnsExtPerProc(nodeHasM(m,i)) + mmax - m + 2
       END DO
       mm = mm + 1
    END DO

    ! any MPI process has at most mnMaxLocal mn's (may have less) and
    ! at most mnExtMaxLocal mnExt's (may have less).
    ! mnMaxLocal and mnExtMaxLocal are used to dimension arrays over all MPI processes

    mnMaxLocal = MAXVAL(mnsPerProc)
    mnExtMaxLocal = MAXVAL(mnsExtPerProc)

    ! myMNMax: scalar containing how many mn's at this MPI process
    ! myMNExtMax: scalar containing how many mnExt's at this MPI process

    myMNMax = mnsPerProc(myId_four)
    myMNExtMax = mnsExtPerProc(myId_four)




    ! MAPPINGS OF LOCAL INDICES mn TO (localm,n) 




    ! Mapping Local pairs (lm,n) to 1D for Regular Spectral:
    ! (1) myMMap(mn): which lm is stored at this position
    ! (2) myNMap(mn): which  n is stored at this position
    ! (3) myMNMap(lm,n): position storing pair (lm,n)

    ALLOCATE (myMNMap(myMMax,nMax))
    ALLOCATE (myMMap(myMNMax))
    ALLOCATE (myNMap(myMNMax))
    myMNMap = -1  ! flag mapping error
    myMMap = -1   ! flag mapping error
    myNMap = -1   ! flag mapping error
    mn = 0
    DO lm = 1, myMMax
       DO n = lm2m(lm), mMax
          mn = mn + 1
          myMNMap(lm,n) = mn
          myMMap(mn)    = lm
          myNMap(mn)    = n
       END DO
    END DO


    ! Mapping Local pairs (lm,n) to 1D for Extended Spectral:
    ! (1) myMExtMap(mn): which lm is stored at this position
    ! (2) myNExtMap(mn): which  n is stored at this position
    ! (3) myMNExtMap(lm,n): position storing pair (lm,n)

    ALLOCATE (myMExtMap(myMNExtMax))
    ALLOCATE (myNExtMap(myMNExtMax))
    ALLOCATE (myMNExtMap(myMMax,nExtMax))
    myMExtMap = -1  ! flag mapping error
    myNExtMap = -1  ! flag mapping error
    myMNExtMap = -1 ! flag mapping error
    mn = 0
    DO lm = 1, myMMax
       DO n = lm2m(lm), mMax+1
          mn = mn + 1
          myMNExtMap(lm,n) = mn
          myMExtMap(mn)    = lm
          myNExtMap(mn)    = n
       END DO
    END DO

  !   spectral communicators for surface field replication 
  !

    IF (ngroups_four.gt.1) THEN
       IF (mygroup_four.EQ.1) THEN
          ALLOCATE (comm_spread(maxnodes,2))
          ALLOCATE (ms_spread(mymmax,ngroups_four))
          ALLOCATE (ini(0:nprocmax_four))
          comm_spread = 0
          mng = 0
          DO n=2,ngroups_four
             ini = 0
             DO m=1,mymmax
                npl = nodeHasM(lm2m(m),n)
                np = map_four(n,npl)
                IF (ini(npl).eq.0) THEN
                   mng = mng + 1
                   ini(npl) = mng
                ENDIF
                mn = ini(npl)
                comm_spread(mn,1) = np
                comm_spread(mn,2) = comm_spread(mn,2) + MMax + 1 - lm2m(m)
                ms_spread(m,n) = mn
             ENDDO
          ENDDO
          comm_spread(:,2) = 2 * comm_spread(:,2) 
          ncomm_spread = mng
        ELSE
          ALLOCATE (comm_spread(npperg_four(1),2))
          ALLOCATE (ms_spread(mymmax,1))
          ALLOCATE (ini(0:nprocmax_four))
          comm_spread = 0
          mng = 0
          ini = 0
          DO m=1,mymmax
             npl = nodeHasM(lm2m(m),1)
             np = map_four(1,npl)
             IF (ini(npl).eq.0) THEN
                mng = mng + 1
                ini(npl) = mng
             ENDIF
             mn = ini(npl)
             comm_spread(mn,1) = np
             comm_spread(mn,2) = comm_spread(mn,2) + MMax + 1 - lm2m(m)
             ms_spread(m,1) = mn
          ENDDO
          comm_spread(:,2) = 2 * comm_spread(:,2) 
          ncomm_spread = mng
       ENDIF
    ENDIF
          
  END SUBROUTINE SpectralDomainDecomp


  SUBROUTINE ThreadDecomp(firstInd, lastInd, minInd, maxInd, msg)
    INTEGER, INTENT(IN ) :: firstInd
    INTEGER, INTENT(IN ) :: lastInd
    INTEGER, INTENT(OUT) :: minInd
    INTEGER, INTENT(OUT) :: maxInd
    CHARACTER(LEN=*), INTENT(IN) :: msg
    INTEGER :: chunk
    INTEGER :: left
    INTEGER :: nTrd
    INTEGER :: iTrd
    INTEGER :: length
    LOGICAL :: inParallel
    LOGICAL :: op
    CHARACTER(LEN=*), PARAMETER :: h="**(ThreadDecomp)**"
    !$ INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
    !$ INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
    !$ LOGICAL, EXTERNAL :: OMP_IN_PARALLEL

    inParallel = .FALSE.
    nTrd=1
    iTrd=0
    !$ inParallel = OMP_IN_PARALLEL()
    IF (inParallel) THEN
       !$ nTrd = OMP_GET_NUM_THREADS()
       !$ iTrd = OMP_GET_THREAD_NUM()
       length = lastInd - firstInd + 1
       chunk = length/nTrd
       left  = length - chunk*nTrd
       IF (iTrd < left) THEN
          minInd =  iTrd   *(chunk+1) + firstInd
          maxInd = (iTrd+1)*(chunk+1) + firstInd - 1
       ELSE
          minInd =  iTrd   *(chunk) + left + firstInd
          maxInd = (iTrd+1)*(chunk) + left + firstInd - 1
       END IF
    ELSE
       minInd = firstInd
       maxInd = lastInd
    END IF

    IF (dumpLocal) THEN

       ! since using unitDump directly (instead of using MsgDump),
       ! check if unitDump in open
       
       INQUIRE(unitDump, opened=op)
       IF (.NOT. op) THEN
          CALL FatalError(h//" unitDump not opened; CreateParallelism not invoked")
       END IF

       IF (inParallel) THEN
          WRITE(unitDump,"(a,' thread ',i2,' got [',i8,':',i8,&
               &'] from [',i8,':',i8,'] in parallel region at ',a)") &
               h, iTrd, minInd, maxInd, firstInd, lastInd, msg
       ELSE
          WRITE(unitDump,"(a,' kept domain [',i8,':',i8,&
               &'] since not in parallel region at ',a)") &
               h, minInd, maxInd, msg
       END IF
    END IF
  END SUBROUTINE ThreadDecomp


  SUBROUTINE ThreadDecompms(m, myms, nms)
    INTEGER, INTENT(IN ) :: m 
    INTEGER, INTENT(INOUT) :: myms(m)
    INTEGER, INTENT(OUT) :: nms   
    INTEGER :: i, j, k
    INTEGER :: nTrd
    INTEGER :: iTrd
    LOGICAL :: inParallel
    !$ INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
    !$ INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
    !$ LOGICAL, EXTERNAL :: OMP_IN_PARALLEL

    inParallel = .FALSE.
    nTrd=1
    iTrd=0
    !$ inParallel = OMP_IN_PARALLEL()
    IF (inParallel) THEN
       !$ nTrd = OMP_GET_NUM_THREADS()
       !$ iTrd = OMP_GET_THREAD_NUM() + 1
       nms = 0
       i = 1
       j = 1
       DO k=1,m
          IF (i.eq.iTrd) THEN
            nms = nms + 1
            myms(nms) = k
          ENDIF
          i = i + j
          IF (i.eq.nTrd+1) THEN
            j = -1
            i = nTrd
          ENDIF
          IF (i.eq.0) THEN
            j = 1
            i = 1
          ENDIF
       ENDDO
    ELSE
       DO k=1,m
          myms(k) = k
       ENDDO
       nms = m
    END IF

  END SUBROUTINE ThreadDecompms

  SUBROUTINE GridDecomposition(ibmax,jbmax,nproc,myid)
    INTEGER, INTENT(IN) :: ibmax
    INTEGER, INTENT(OUT) :: jbmax
    INTEGER, INTENT(IN) :: nproc
    INTEGER, INTENT(IN) :: myid
    INTEGER :: np2, np1, np, ngroups, npperg, nrest, npr
    INTEGER :: ngptotal, lat, usedinlat, next, ndim, jbdim
    INTEGER :: ngpperproc, n, imp2, ndim_f, iovmax, iadd
    INTEGER :: ipar, i1, i2, il1, il2, j1, j2, ij, ij1, ind, aux(4)
    INTEGER :: firstextlat, lastextlat, iex, mygridpoints
    INTEGER :: i, k, iproc, j, jlo(1), ib, jb, ijb, jh
    INTEGER, ALLOCATABLE :: lon(:), iovlap(:)
    INTEGER, ALLOCATABLE :: nprocsingroup(:)
    INTEGER, ALLOCATABLE :: npointsingroup(:)
    INTEGER, ALLOCATABLE :: firstlatingroup(:)
    INTEGER, ALLOCATABLE :: lastlatingroup(:)
    INTEGER, ALLOCATABLE :: firstloningroup(:)
    INTEGER, ALLOCATABLE :: lastloningroup(:)
    INTEGER, ALLOCATABLE :: jbmaxingroup(:)
    INTEGER, ALLOCATABLE :: procingroup(:)
    INTEGER, ALLOCATABLE :: iaux(:,:)
    REAL (KIND=r8) :: hj, ifirst, ilast
    CHARACTER(LEN=*), PARAMETER :: h="**(GridDecomposition)**"

    IF (nproc.le.3) THEN
       ngroups = nproc
       npperg = 1
       nrest = 0
       iex = 0
     ELSE 
       IF (nproc.le.5) THEN 
          ngroups = 1
          npperg = nproc-1
          nrest = 0
          iex = 1
        ELSE
          np2 = nproc / 2
          np = SQRT(REAL(np2))
          IF (np*(np+1).lt.np2) THEN
             ngroups = np+1
           ELSE
             ngroups = np
          ENDIF
          npperg = nproc / ngroups
          nrest = nproc - npperg * ngroups
          IF (npperg.le.8) THEN
             iex = 1
           ELSE
             iex = 2
          ENDIF
       ENDIF
    ENDIF
    hj = ACOS(-1._r8) / jmax
    ALLOCATE (gridmap(1:iMax,1:jMax))
    ALLOCATE (procingroup(0:nproc-1))
    ALLOCATE (nprocsingroup(ngroups+2*iex))
    nprocsingroup(1)= 1
    nprocsingroup(ngroups+2*iex)= 1
    IF (iex.eq.2) THEN
       npr = SQRT(REAL(npperg))
       nprocsingroup(iex) = npr + 1
       nprocsingroup(iex+1) = npperg - npr - 2
       nprocsingroup(ngroups+3) = npr + 1
       nprocsingroup(ngroups+iex) = npperg - npr - 2
     ELSE
       IF (iex.eq.1) THEN
          nprocsingroup(2) = npperg - 1
          nprocsingroup(ngroups+1) = npperg - 1
       ENDIF
    ENDIF
    nprocsingroup(2+iex:ngroups+iex-1) = npperg
    DO k=1,nrest
       nprocsingroup(1+iex+k) = nprocsingroup(1+iex+k) + 1
    ENDDO
    ngroups = ngroups +2*iex
    i = -1
    DO k=1,ngroups
       procingroup(i+1:i+nprocsingroup(k))=k
       i = i+nprocsingroup(k)
    ENDDO
    ALLOCATE (npointsingroup(ngroups))
    ALLOCATE (pointsinproc(0:nproc-1))
    ALLOCATE (firstlatingroup(ngroups))
    ALLOCATE (firstloningroup(ngroups))
    ALLOCATE (lastlatingroup(ngroups))
    ALLOCATE (lastloningroup(ngroups))
    ALLOCATE (jbmaxingroup(ngroups))
    ngptotal = SUM(imaxperj(1:jmax))
    ngpperproc = ngptotal / nproc
    nrest = ngptotal - ngpperproc * nproc
    pointsinproc(0:nrest-1) = ngpperproc + 1
    pointsinproc(nrest:nproc-1) = ngpperproc
    lat = 1
    usedinlat = 0
    iproc = 0
    DO k=1,ngroups
       next = MIN(nrest,nprocsingroup(k))
       npointsingroup(k) = ngpperproc * nprocsingroup(k) + next
       nrest = nrest - next
       firstlatingroup(k) = lat
       firstloningroup(k) = usedinlat + 1
       np = imaxperj(lat) - usedinlat
       DO 
          IF (np.ge.npointsingroup(k)) EXIT
          lat = lat + 1
          np = np + imaxperj(lat)
       ENDDO
       lastlatingroup(k) = lat
       usedinlat = imaxperj(lat) - np + npointsingroup(k)
       lastloningroup(k) = usedinlat
       IF (usedinlat.eq.imaxperj(lat)) THEN
          lat = lat +1
          usedinlat = 0
       ENDIF 
    ENDDO
    jbmaxingroup = lastlatingroup - firstlatingroup + 1
    jbdim = MAXVAL(jbmaxingroup)
    jovlap = 0
    ndim = 4*jbdim+2*jovlap*npperg
    ALLOCATE (mysendsgr(4,ndim))
    ALLOCATE (mysendspr(2,nproc))
    ALLOCATE (myrecsgr(4,ndim))
    ALLOCATE (myrecspr(2,nproc))
    ALLOCATE (firstlat(0:nproc-1))
    ALLOCATE (firstlon(jbdim,0:nproc-1))
    ALLOCATE (lastlat(0:nproc-1))
    ALLOCATE (lastlon(jbdim,0:nproc-1))
    ALLOCATE (lon(jbdim))

    iproc = 0
    DO k=1,ngroups
       firstlon(1,iproc) = firstloningroup(k)
       firstlon(2:jbmaxingroup(k),iproc) = 1
       lastlon(1:jbmaxingroup(k),iproc) = 0
       lon(1:jbmaxingroup(k)) = firstlon(1:jbmaxingroup(k),iproc)
       DO n=1,nprocsingroup(k)
          DO np=1,pointsinproc(iproc)
             jlo = MINLOC(REAL(lon(1:jbmaxingroup(k))) &
             / REAL(imaxperj(firstlatingroup(k):lastlatingroup(k))) )
             j = jlo(1)
             lastlon(j,iproc) = lon(j)
             gridmap(lon(j),j+firstlatingroup(k)-1) = iproc
             lon(j) = lon(j)+1
             IF (j.eq.jbmaxingroup(k).and.lon(j).gt.lastloningroup(k)) &
                lon(j) = imaxperj(lastlatingroup(k))+1
          ENDDO
          IF (lastlon(1,iproc).eq.0) THEN 
             firstlat(iproc) = firstlatingroup(k)+1
           ELSE
             firstlat(iproc) = firstlatingroup(k)
          ENDIF
          IF (lastlon(jbmaxingroup(k),iproc).eq.0) THEN 
             lastlat(iproc) = lastlatingroup(k)-1
           ELSE 
             lastlat(iproc) = lastlatingroup(k)
          ENDIF
          iproc = iproc + 1
          IF (iproc.ne.nproc) THEN
             firstlon(1:jbmaxingroup(k),iproc) = lon(1:jbmaxingroup(k))
             lastlon(1:jbmaxingroup(k),iproc) = 0
          ENDIF
       ENDDO
    ENDDO
    mygridpoints = pointsinproc(myid)

    jbmax = mygridpoints / ibmax 
    if(jbmax*ibmax.lt.mygridpoints) jbmax = jbmax + 1
    ALLOCATE (iPerIJB(ibmax,jbmax))
    ALLOCATE (jPerIJB(ibmax,jbmax))
    ALLOCATE (ibMaxPerJB(jbmax))

    myfirstlat = firstlat(myid)
    mylastlat = lastlat(myid)
    firstextlat = MAX(myfirstlat-jovlap,1)
    lastextlat = MIN(mylastlat+jovlap,jmax)
    ALLOCATE (myfirstlon(firstextlat:lastextlat))
    ALLOCATE (mylastlon(firstextlat:lastextlat))
    ALLOCATE (iovlap(firstextlat:lastextlat))
    DO j=firstextlat,lastextlat
       i = nint(jovlap*imaxperj(j)/(imax*SIN((j-.5_r8)*hj)))
       iovlap(j) = i + 1
    enddo
       iovmax = 0

    ALLOCATE (ibPerIJ(1-iovmax:iMax+iovmax,-1:jMax+2 ))
    ALLOCATE (jbPerIJ(1-iovmax:iMax+iovmax,-1:jMax+2 ))
    ibPerIJ = 0
    jbPerIJ = 0

    myfirstlon(myfirstlat:mylastlat) = &
      firstlon(1+myfirstlat-firstlatingroup(procingroup(myid)): &
      jbmaxingroup(procingroup(myid))-lastlatingroup(procingroup(myid)) &
      + mylastlat,myid)
    mylastlon(myfirstlat:mylastlat) = &
      lastlon(1+myfirstlat-firstlatingroup(procingroup(myid)): &
      jbmaxingroup(procingroup(myid))-lastlatingroup(procingroup(myid)) &
      + mylastlat,myid)
    ifirst = MINVAL(REAL(myfirstlon(myfirstlat:mylastlat)) &
                     / REAL(imaxperj(myfirstlat:mylastlat)) )
    ilast  = MAXVAL(REAL(mylastlon(myfirstlat:mylastlat)) &
                     / REAL(imaxperj(myfirstlat:mylastlat)) )

    !  interior domain
    !  ---------------
    ijb = 0
    DO j=myfirstlat,mylastlat
       DO i=myfirstlon(j),mylastlon(j)
          ib = MOD(ijb,ibmax)+1
          jb = ijb / ibmax + 1
          ijb = ijb + 1
          iPerIJB(ib,jb) = i
          jPerIJB(ib,jb) = j
          ibPerIJ(i,j) = ib
          jbPerIJ(i,j) = jb
       ENDDO
    ENDDO
    ibMaxPerJB(1:jbmax-1) = ibmax
    ibMaxPerJB(jbmax) = ib


    ! define messages to be exchanged between fourier and grid computations
    ! ---------------------------------------------------------------------
    ndim_f = jMaxlocal_f*(MAXVAL(nprocsingroup)+3)
    ALLOCATE (messages_f(4,ndim_f))
    ALLOCATE (messproc_f(2,0:nproc))
    ALLOCATE (messages_g(4,ndim_f))
    ALLOCATE (messproc_g(2,0:nproc))
    ipar = 1
    DO j=myfirstlat_f,mylastlat_f
       messages_f(1,ipar) = 1
       messages_f(3,ipar) = j
       messages_f(4,ipar) = gridmap(1,j)
       DO i=2,imaxperj(j)
          IF (gridmap(i,j).ne.messages_f(4,ipar)) THEN
             messages_f(2,ipar) = i-1
             IF(messages_f(4,ipar).ne.myid) ipar = ipar + 1
             messages_f(1,ipar) = i
             messages_f(3,ipar) = j
             messages_f(4,ipar) = gridmap(i,j)
          ENDIF
          IF (i.eq.imaxperj(j)) THEN
             messages_f(2,ipar) = i
             IF(messages_f(4,ipar).ne.myid) ipar = ipar + 1
          ENDIF
       ENDDO
    ENDDO
    ipar = ipar - 1
    IF (ipar.gt.ndim_f) THEN
       WRITE(nfprt,*) ' ndim_f, ipar  ',ndim_f,ipar
       WRITE(nfprt,"(a, ' dimensioning of segment messages insufficient')") h
       STOP h
    ENDIF

    !sort messages by processors
    !---------------------------
    DO i=2,ipar
       DO j=i,2,-1
          IF(messages_f(4,j).lt.messages_f(4,j-1)) THEN
             aux = messages_f(:,j-1)
             messages_f(:,j-1) = messages_f(:,j)
             messages_f(:,j) = aux
            ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO
    messproc_f(2,0) = 0
    IF (ipar.gt.0) THEN
       messproc_f(1,1) = messages_f(4,1)
       nrecs_f = 1
       DO i=2,ipar
          IF(messages_f(4,i).ne.messages_f(4,i-1)) THEN
             nrecs_f = nrecs_f + 1
             messproc_f(2,nrecs_f-1) = i - 1
             messproc_f(1,nrecs_f) = messages_f(4,i)
          ENDIF
       ENDDO
       messproc_f(2,nrecs_f) = ipar
     ELSE
       nrecs_f = 0
    ENDIF

    !  set communication structure for grid diagnostics
    !  ------------------------------------------------
    ALLOCATE(firstandlastlat(2,0:maxnodes-1))
    ALLOCATE (nlatsinproc_d(0:maxNodes-1))
    lat = jmax / maxnodes
    nrest = jmax - lat * maxnodes
    np = nrest / 2
    np1 = nrest - np
    n = 0
    DO i=0,maxnodes-1
       firstandlastlat(1,i) = n + 1
       IF (i.lt.np.or.i.ge.maxnodes-np1) THEN
          n = n + lat + 1
         ELSE
          n = n + lat
       ENDIF
       firstandlastlat(2,i) = n
    ENDDO
    nlatsinproc_d = firstandlastlat(2,:) - firstandlastlat(1,:) + 1
    myfirstlat_diag = firstandlastlat(1,myid)
    mylastlat_diag = firstandlastlat(2,myid)
    myJMax_d = mylastlat_diag - myfirstlat_diag + 1
         
    ndim = (5+MAXVAL(nprocsingroup))*(lat+1)
    ALLOCATE (mysends_diag(4,ndim))
    ALLOCATE (mysendspr_diag(2,0:nproc))
    ALLOCATE (myrecs_diag(4,ndim))
    ALLOCATE (myrecspr_diag(2,0:nproc))

    ij = 0
    DO j=myfirstlat_diag,mylastlat_diag
       ij = ij + 1
       myrecs_diag(1,ij) = 1
       myrecs_diag(3,ij) = j
       myrecs_diag(4,ij) = gridmap(1,j)
       DO i=2,imaxperj(j)
          IF (gridmap(i,j).ne.myrecs_diag(4,ij)) THEN
             myrecs_diag(2,ij) = i-1
             ij = ij + 1
             myrecs_diag(1,ij) = i
             myrecs_diag(3,ij) = j
             myrecs_diag(4,ij) = gridmap(i,j)
          ENDIF
       ENDDO
       myrecs_diag(2,ij) = imaxperj(j)
    ENDDO

    !sort messages by processors
    !---------------------------
    DO i=2,ij   
       DO j=i,2,-1
          IF(myrecs_diag(4,j).lt.myrecs_diag(4,j-1)) THEN
             aux = myrecs_diag(:,j-1)
             myrecs_diag(:,j-1) = myrecs_diag(:,j)
             myrecs_diag(:,j) = aux
            ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO
    n = 0
    i1 = 0
    DO i=1,ij 
       IF(myrecs_diag(4,i).eq.myid) THEN
         i1 = i
         n = n + 1
       ENDIF
    ENDDO
    ij = ij - n
    DO i=i1-n+1,ij
       myrecs_diag(:,i) = myrecs_diag(:,i+n)
    ENDDO
    myrecspr_diag(2,0) = 0
    IF (ij.gt.0) THEN
       myrecspr_diag(1,1) = myrecs_diag(4,1)
       nrecs_diag = 1
       DO i=2,ij
          IF(myrecs_diag(4,i).ne.myrecs_diag(4,i-1)) THEN
             nrecs_diag = nrecs_diag + 1
             myrecspr_diag(2,nrecs_diag-1) = i - 1
             myrecspr_diag(1,nrecs_diag) = myrecs_diag(4,i)
          ENDIF
       ENDDO
       myrecspr_diag(2,nrecs_diag) = ij
     ELSE
       nrecs_diag = 0
    ENDIF

  END SUBROUTINE GridDecomposition

END MODULE Sizes

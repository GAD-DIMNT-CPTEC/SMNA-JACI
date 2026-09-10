!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Sizes </br></br>
!#
!# **Brief**: Deals with parallelism for Modules requires more than one processor,
!# but not all. Overrides Mod_Parallelism maxNodes and mpiMasterProc. </br>
!#
!# Single repository of sizes and mappings of grid, fourier and spectral field
!# representations. Should be used by all model modules that require such information. </br>
!#
!# Module exports three procedures:
!# <ul type="disc">
!#  <li>RegisterBasicSizes: <ul type="disc"><li> set all repository values for spectral
!# representation. It also sets latitudes and maximum number of longitudes for
!# grid representation. </li></ul></li>
!#
!#  <li>RegisterOtherSizes: <ul type="disc"><li> set remaining repository values. </li></ul></li>
!#
!#  <li>DumpSizes: <ul type="disc"><li> detailed dumping of repository values. </li></ul></li>
!# </ul>
!# 
!# Module exports a set of values and arrays, described bellow. </br>
!# 
!# Notation used throughout the code: 
!# <ul type="disc">
!#  <li>m  Legendre Order (also Fourier wave number). Goes from 1 to mMax. </li>
!#  <li>n  Legendre Degree. Goes from m to nMax on regular spectral fields and from
!# m to nExtMax on extended spectral fields. </li>
!#  <li>mn index to store the spectral triangle (m,n) into a single dimension. </li>
!#  <li>k  vertical index, from 1 to kMax. </li>
!#  <li>j  latitude index, from 1 to jMax (full field) or from 1 to jMaxHalf (hemisphere). </li>
!#  <li>i  longitude index, from 1 to iMax (gaussian grid) or from 1 to iMaxPerJ(j)
!# (at latitude j of the reduced grid). </li>
!#  <li>ib index of a longitude on a block of longitudes packed in one dimension of
!# grids. </li>
!#  <li>jb index for which block of longitudes to take, packed in another dimension
!# of grids. </li>
!# </ul>
!# 
!# SPECTRAL REPRESENTATION: 
!# <ul type="disc">
!#  <li>A regular spectral field should be dimensioned (2*mnMax,kMax), where the
!# first dimension accomodates pairs of (real,imaginary) spectral coefficients
!# for a fixed vertical and the second dimension varies verticals. </li>
!#  <li>Extended spectral fields should use (2*mnExtMax,kMax). </li>
!#  <li> This module provides mapping functions to map (m,n) into mn:
!# <ul type="disc">
!#  <li>For a regular spectral field with complex entries, Legendre Order m and
!# Degree n has real coefficient at position 2*mnMap(m,n)-1 and </li>
!#  <li>imaginary coefficient at position 2*mnMap(m,n). Inverse mappings (from mn to
!# m and n) are also provided. Maps for the extended spectral field are also provided. </li>
!# </ul> </li>
!# </ul>
!#
!# GRID REPRESENTATION: </br>
!# <ul type="disc">
!#  <li>Since the number of longitudes per latitude may vary with latitude (on
!# reduced grids), sets of latitudes (with all longitudes) are packed together
!# in the first dimension of Grids. </li>
!#  <li>Near the pole, many latitudes are packed; near the equator, just a few.
!# Second dimension is vertical. </li>
!#  <li>Third dimension is the number of packed latitudes required to represent
!# a full field. </li>
!#  <li>A grid field should be dimensioned (ibMax, kMax, jbMax). </li>
!#  <li>This module provides mapping functions to map latitude j and longitude i
!# into (ib,jb). </li>
!#  <li>Map ibPerIJ(i,j) gives index of first dimension that stores longitude i
!# of latitude j. Map jbPerIJ(i,j) gives index of third dimension that stores
!# longitude i of latitude j. </li>
!#  <li>Map jbPerIJ(i,j) gives index of third dimension that stores longitude i
!# of latitude j. </li>
!#  <li>Consequently, the point of longitude i, latitude j and vertical k is
!# stored at (ibPerIJ(i,j), k, jbPerIJ(i,j)). </li>
!#  <li>Inverse mappings (from (ib,jb) to (i,j)) are also provided. </li>
!# </ul>
!# 
!# FOURIER REPRESENTATION: </br>
!# For the moment, Fourier fields are represented externally to the transform.
!# That should not happen in the future. </br>
!# First dimension contains pairs of (real,imaginary) fourier coefficients.
!# Second dimension is latitude. Third dimension is vertical. A full fourier
!# field should be dimensioned (iMax+1, jMax, kMax). </br></br>
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
!#  <li>20-04-2010 - Paulo Kubota - version: 1.16.0 </li>
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


module Mod_Sizes   ! Version 0 of Nov 25th, 2001

  use Mod_Parallelism_Group_Chopping, only : &
    msgDump &
    , msgOutMaster &
    , fatalError &
    , maxnodes &
    , myId &
    , getUnitDump

  use Mod_Parallelism_Fourier, only : &
    mygroup_four &
    , maxnodes_four &
    , myId_four

  implicit none

  include 'mpif.h'
  include 'precision.h'
  include 'messages.h'

  ! INTERNAL DATA STRUCTURE:
  ! Provision for MPI computation:
  
  ! Domain Decomposition:
  ! There are maxNodes MPI processes, numbered from 0 to maxNodes - 1
  ! The number of this MPI process is myId, 0 <= myId <= maxNodes - 1
  ! SPECTRAL REPRESENTATION:
  integer :: mMax = -1
  !# is the maximum Legendre Order (also maximum Fourier wave number).
  !      it is set to model truncation + 1 (Legendre Order 0 is m=1)
  integer :: mMax_in = -1
  !# nMax is the maximum Legendre Degree for regular spectral fields.
  !      it is set to model truncation + 1 (Legendre Degree 1 is n=1)
  integer :: nMax = -1
  integer :: mnMax = -1
  !# is the amount of spectral coeficients per vertical, for regular spectral
  !# fields. It is the number of points in the regular triangular spectral plane.
  integer :: mnMax_out = -1
  integer, allocatable :: mMap(:)
  !# mMap indexed by mn in [1:mnMax]; returns m at this mn for regular spectral
  !# fields;
  integer, allocatable :: nMap(:)
  !# nMap indexed by mn in [1:mnMax]; returns n at this mn for regular spectral
  !# fields;
  integer, allocatable :: mnMap(:, :)
  !# mnMap indexed by (m,n); returns mn that stores spectral coefficient (m,n),
  !# when spectral coefficients are real (e.g. Associated Legendre Functions).
  
  ! For complex spectral coefficients, mn=2*mnMap(m,n) for the imaginary part
  ! and mn=2*mnMap(m,n)-1 for the real part.
  
  integer, allocatable :: mnMap_out(:, :)
  
  ! Remaining variables with Ext extension have the same meaning, but
  ! for the extended spectral space (where nExtMax=trunc+2), that is,
  ! where each Legendre Order has one more Legendre Degree than usual.

  ! Spectral coefficients are spread across processes.
  ! The set of values of m is partitioned among processes.
  ! Consequently, all values of n for some values of m are stored
  ! at each process. Each value of m is stored at a single process.
  integer :: nExtMax = -1
  integer :: mnExtMax = -1
  integer, allocatable :: mExtMap(:)
  integer, allocatable :: nExtMap(:)
  integer, allocatable :: mnExtMap(:, :)
  integer, allocatable :: nodeHasM_out(:, :)
  integer, allocatable :: nodeHasM(:, :)
  !# Array nodeHasM(m) has which MPI process stores m, 1 <= m <= mMax
  integer, allocatable :: lm2m(:)
  integer, allocatable :: msInProc(:, :)
  integer, allocatable :: msInProc_out(:, :)
  integer, allocatable :: msPerProc(:)
  integer, allocatable :: msPerProc_out(:)
  integer :: myMMax
  !# Variable myMMax (0 <= myMMax <= mMax) has how many m's are stored at this
  !# node.
  
  ! Array lm2m(l), l=1,...,myMMax has the m's stored at this node. That is, it
  ! maps local m (lm) into global m (m).
  
  integer :: myMMax_out
  integer :: mMaxLocal
  integer :: mnMaxLocal
  integer :: mnExtMaxLocal
  integer :: myMNMax
  integer :: myMNExtMax
  logical :: HaveM1

  integer, allocatable :: myMMap(:)
  integer, allocatable :: myNMap(:)
  integer, allocatable :: mnsPerProc(:)
  integer, allocatable :: mnsExtPerProc(:)
  integer, allocatable :: myMNMap(:, :)
  integer, allocatable :: myMExtMap(:)
  integer, allocatable :: myNExtMap(:)
  integer, allocatable :: myMNExtMap(:, :)

  !  Spectral representation and division to be used in the semi-implicit part
  !  -------------------------------------------------------------------------

  !  Spectral communicators to spread surface field
  !  ----------------------------------------------

  integer :: ncomm_spread
  integer, allocatable :: comm_spread(:, :)
  integer, allocatable :: ms_spread(:, :)

  ! LATITUDES REPRESENTATION:
  integer :: jMax = -1
  !# jMax is the number of latitudes for full Fourier and Grid representations.
  integer :: jMaxHalf = -1
  !# jMaxHalf is the number of latitudes at each hemisphere.
  integer, allocatable :: mMaxPerJ(:)
  !# mMaxPerJ is an array indexed by latitudes that stores the maximum value of
  !# m at each latitude. 
  
  ! Consequently, the contribution of latitude j for the Legendre Transform
  ! should be taken only up to Legendre Order mMaxPerJ(j). By the same token,
  ! FFTs at latitude j should be computed only up to Fourier wave number
  ! mMaxPerJ(j). For regular grids, mMaxPerJ(j) = mMax for all j.
  
  integer, allocatable :: jMinPerM(:)
  !# jMinPerM is the inverse mapping - an array indexed by Legendre Orders (m)
  !# containing the smallest latitude (value of j) that has that order. 
  
  ! Latitudes that contain Legendre Order (and Fourier wave number) m are
  ! jMinPerM(m) <= j <= jMax-jMinPerM(m)+1. For regular grids, jMinPerM(m) = 1
  ! for all m.
  
  integer, allocatable :: jMaxPerM(:)

  ! LONGITUDES REPRESENTATION:
  integer :: iMax = -1
  !# iMax is the maximum number of longitudes per latitude; 
  
  ! it is the number of longitudes per latitude for all latitudes on regular grids;
  ! it is only the maximum number of longitudes per latitude on reduced grids; it
  ! is the actual number of longitudes per latitude close to the equator on reduced grids;
  
  integer :: ijMax = -1
  !# ijMax is the number of horizontal grid points at regular or reduced grids.
  integer :: ijMaxGauQua = -1
  integer, allocatable :: iMaxPerJ(:)
  !# iMaxPerJ is the actual number of longitudes per latitude on regular and reduced grids;
  !# latitude j has iMaxPerJ(j) longitudes.

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
  
  ! Blocks of surface points are spreaded across processes.
  ! Each process has all longitudes of a set of latitudes (block)

  !  decomposition in Fourier space

  integer, allocatable :: ibMaxPerJB(:)
  integer, allocatable :: firstlatinproc_f(:)
  integer, allocatable :: lastlatinproc_f(:)
  integer, allocatable :: nlatsinproc_f(:)
  integer, allocatable :: nodeHasJ_f(:)
  integer, allocatable :: kfirst_four(:)
  integer, allocatable :: klast_four(:)
  integer, allocatable :: nlevperg_four(:)
  integer, allocatable :: npperg_four(:)
  integer, allocatable :: map_four(:, :)
  integer, allocatable :: first_proc_four(:)
  integer, allocatable :: nlatsinproc_d(:)
  integer, allocatable :: messages_f(:, :)
  integer, allocatable :: messages_g(:, :)
  integer, allocatable :: messproc_f(:, :)
  integer, allocatable :: messproc_g(:, :)
  integer, allocatable :: nodeHasJ(:)
  integer, allocatable :: nset(:)

  integer :: myfirstlat_f
  integer :: mylastlat_f
  integer :: myfirstlev
  integer :: mylastlev
  integer :: myJMax_f
  integer :: ngroups_four
  integer :: nprocmax_four
  integer :: JMaxlocal_f
  integer :: kMaxloc
  integer :: kMaxloc_in
  integer :: kMaxloc_out
  integer :: nrecs_f
  integer :: nrecs_g
  logical :: havesurf

  ! grid decomposition

  integer :: ibMax = -1
  integer :: jbMax = -1
  integer :: jbMax_ext
  integer :: myfirstlat
  integer :: mylastlat
  integer :: nrecs_gr
  integer :: nsends_gr
  integer :: nrecs_diag
  integer :: nsends_diag
  integer :: myfirstlat_diag
  integer :: mylastlat_diag
  integer :: myJMax_d
  integer :: jovlap
  integer, allocatable :: firstandlastlat(:, :)
  integer, allocatable :: myfirstlon(:)
  integer, allocatable :: mylastlon(:)
  integer, allocatable :: firstlon(:, :)
  integer, allocatable :: lastlon(:, :)
  integer, allocatable :: firstlat(:)
  integer, allocatable :: lastlat(:)
  integer, allocatable :: ibPerIJ(:, :)
  integer, allocatable :: jbPerIJ(:, :)
  integer, allocatable :: iPerIJB(:, :)
  integer, allocatable :: jPerIJB(:, :)
  integer, allocatable :: pointsinproc(:)
  integer, allocatable :: myrecsgr(:, :)
  integer, allocatable :: myrecspr(:, :)
  integer, allocatable :: mysendsgr(:, :)
  integer, allocatable :: mysendspr(:, :)
  integer, allocatable :: myrecs_diag(:, :)
  integer, allocatable :: myrecspr_diag(:, :)
  integer, allocatable :: mysends_diag(:, :)
  integer, allocatable :: mysendspr_diag(:, :)
  integer, allocatable :: gridmap(:, :)

  integer :: kMax = -1

  real(kind = p_r8), allocatable :: ci(:)      
  !# 1 - sigma each level (level 1 at surface)
  real(kind = p_r8), allocatable :: si(:)      
  !# sigma
  real(kind = p_r8), allocatable :: del(:)     
  !# layer thickness (in sigma)
  real(kind = p_r8), allocatable :: delcl(:)     
  !# layer thickness (in cl)
  real(kind = p_r8), allocatable :: rdel2(:)
  real(kind = p_r8), allocatable :: sl(:)      
  !# sigma at layer midpoint
  real(kind = p_r8), allocatable :: cl(:)      
  !# 1.0 - sl
  real(kind = p_r8), allocatable :: rpi(:)     
  !# 'pi' ratios at adjacent layers
  logical, parameter, private :: dumpLocal = .false.

contains


  subroutine RegisterBasicSizes(trunc_in, trunc_out, nLat, nLon, vert)
    !# Records Basic Sizes
    !# ---
    !# @info
    !# **Brief:** Records Basic Sizes.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin  
    integer, intent(in) :: trunc_in
    integer, intent(in) :: trunc_out
    integer, intent(in) :: nLat
    integer, intent(in) :: nLon
    integer, intent(in) :: vert
    integer :: m, n, mn
    character(len = *), parameter :: h = "**(RegisterBasicSizes)**"

    jMax = nLat
    jMaxHalf = nLat / 2

    iMax = nLon

    kMax = vert

    mMax_in = trunc_in + 1
    mMax = trunc_out + 1
    nMax = mMax
    nExtMax = mMax + 1
    mnExtMax = (nExtMax + 2) * (nExtMax - 1) / 2
    if (allocated(mExtMap))then
      deallocate (mExtMap)
      allocate (mExtMap(mnExtMax))
    else
      allocate (mExtMap(mnExtMax))
    end if
    if (allocated(nExtMap))then
      deallocate (nExtMap)
      allocate (nExtMap(mnExtMax))
    else
      allocate (nExtMap(mnExtMax))
    end if
    if (allocated(mnExtMap))then
      deallocate (mnExtMap)
      allocate (mnExtMap(mMax, nExtMax))
    else
      allocate (mnExtMap(mMax, nExtMax))
    end if

    nExtMap = -1  
    ! flag mapping error
    mExtMap = -1  
    ! flag mapping error
    mnExtMap = -1 
    ! flag mapping error
    mn = 0
    do m = 1, mMax
      do n = m, mMax + 1
        mn = mn + 1
        mnExtMap(m, n) = mn
        mExtMap(mn) = m
        nExtMap(mn) = n
      end do
    end do
    mnMax = (mmax_in * (mMax_in + 1)) / 2
    mnMax_out = (mmax * (mMax + 1)) / 2
    if (myid.eq.0) write(*, *) ' mnMax ', mnMax, ' mnMax_out ', mnMax_out
    if (allocated(mnMap))deallocate  (mnMap)
    if (allocated(mnMap_out))deallocate  (mnMap_out)
    if (allocated(mMap))deallocate  (mMap)
    if (allocated(nMap))deallocate  (nMap)
    allocate (mnMap(mMax_in, mMax_in))
    allocate (mnMap_out(mMax, mMax))
    allocate (mMap(mnMax))
    allocate (nMap(mnMax))
    mnMap = -1  
    ! flag mapping error
    mnMap_out = -1  
    ! flag mapping error
    mMap = -1   
    ! flag mapping error
    nMap = -1   
    ! flag mapping error
    mn = 0
    do m = 1, mMax_in
      do n = m, mMax_in
        mn = mn + 1
        mnMap(m, n) = mn
        mMap(mn) = m
        nMap(mn) = n
      end do
    end do
    mn = 0
    do m = 1, mMax
      do n = m, mMax
        mn = mn + 1
        mnMap_out(m, n) = mn
      end do
    end do
    ijMaxGauQua = iMax * jMax

  end subroutine RegisterBasicSizes


  subroutine VerticalGroups(givenfouriergroups, nproc_vert, kmax_in, &
    !# Computes vertical decomposition of fourier groups
    !# ---
    !# @info
    !# **Brief:** Computes vertical decomposition of fourier groups.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin  
    kmaxloc_io)
    logical, intent(in) :: givenfouriergroups
    integer, intent(in) :: nproc_vert
    integer, intent(in) :: kmax_in
    integer, intent(OUT) :: kmaxloc_io
    character(len = *), parameter :: h = "**(RegisterOtherSizes)**"
    integer :: i, n, nlp, nrest, ipar, ng, npsq, nl, nn, np, nlg
    integer :: next, k, m, nprest, l

    if (givenfouriergroups) then
      ng = nproc_vert
    else
      npsq = sqrt(real(maxNodes))
      do nl = 1, kmax_in
        nn = nl * maxNodes / kmax_in
        if(nn.ge.npsq) EXIT
      enddo
      ng = kmax_in / nl
      if (nl * ng.lt.kmax_in) ng = ng + 1
    endif
    ngroups_four = ng
    nl = kmax_in / ng
    if (nl * ng.lt.kmax_in) nl = nl + 1
    nrest = nl * ng - kmax_in
    if (allocated(kfirst_four))deallocate  (kfirst_four)
    if (allocated(klast_four))deallocate  (klast_four)
    if (allocated(npperg_four))deallocate  (npperg_four)
    if (allocated(nlevperg_four))deallocate  (nlevperg_four)
    if (allocated(first_proc_four))deallocate  (first_proc_four)
    allocate(kfirst_four    (0:maxNodes))
    allocate(klast_four     (0:maxNodes))
    allocate(npperg_four    (ng))
    allocate(nlevperg_four  (ng))
    allocate(first_proc_four(ng))
    np = 0
    do i = 1, ng
      if (i.le.nrest) then
        nlg = nl - 1
      else
        nlg = nl
      endif
      nlevperg_four(i) = nlg
      npperg_four(i) = nlg * maxNodes / kmax_in
      np = np + npperg_four(i)
    enddo
    nprest = maxNodes - np
    next = nprest / ng
    nprest = nprest - ng * next
    do i = 1, ng
      if (i.le.nprest) then
        npperg_four(i) = npperg_four(i) + next + 1
      else
        npperg_four(i) = npperg_four(i) + next
      endif
    enddo
    if (allocated(nset))deallocate  (nset)
    allocate(nset(ng))
    nprocmax_four = maxval(npperg_four)
    if (allocated(map_four))deallocate  (map_four)
    allocate(map_four(ng, 0:nprocmax_four - 1))
    n = 0
    nset = 0
    do
      do i = 1, ngroups_four
        if (nset(i).lt.npperg_four(i)) then
          if (myid.eq.n) then
            mygroup_four = i
            maxnodes_four = npperg_four(i)
            myid_four = nset(i)
          endif
          map_four(i, nset(i)) = n
          nset(i) = nset(i) + 1
          n = n + 1
        endif
      end do
      if (n.ge.maxNodes) EXIT
    end do
    m = 0
    do i = 1, ngroups_four
      do k = 0, npperg_four(i) - 1
        l = map_four(i, k)
        kfirst_four(l) = m + 1
        klast_four(l) = m + nlevperg_four(i)
      enddo
      m = m + nlevperg_four(i)
    enddo
    do i = 1, ngroups_four
      first_proc_four(i) = i - 1
    enddo
    kmaxloc_io = klast_four(myid) - kfirst_four(myid) + 1
    myfirstlev = kfirst_four(myid)
    mylastlev = klast_four(myid)

  end subroutine VerticalGroups

  subroutine ReshapeVerticalGroups(kmax_in, kmaxloc_io)
    !# Reshapes vertical decomposition of fourier groups
    !# ---
    !# @info
    !# **Brief:** Reshapes vertical decomposition of fourier groups.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin  
    integer, intent(in) :: kmax_in
    integer, intent(OUT) :: kmaxloc_io
    character(len = *), parameter :: h = "**(RegisterOtherSizes)**"
    integer :: i, n, nlp, nrest, ipar, ng, npsq, nl, nn, np, nlg
    integer :: next, k, m, nprest, l

    ng = ngroups_four
    nl = kmax_in / maxNodes
    np = 0
    do i = 1, ng
      np = np + npperg_four(i) * nl
      nlevperg_four(i) = npperg_four(i) * nl
    enddo
    nrest = kmax_in - np
    do
      do i = ng, 1, -1
        np = min (1, nrest)
        nlevperg_four(i) = nlevperg_four(i) + np
        nrest = nrest - np
      enddo
      if (nrest.eq.0) EXIT
    enddo
    if (nlevperg_four(1).eq.0) then
      nlevperg_four(1) = 1
      nlevperg_four(ng) = nlevperg_four(ng) - 1
    endif
    m = 0
    do i = 1, ngroups_four
      do k = 0, npperg_four(i) - 1
        l = map_four(i, k)
        kfirst_four(l) = m + 1
        klast_four(l) = m + nlevperg_four(i)
      enddo
      m = m + nlevperg_four(i)
    enddo
    kmaxloc_io = klast_four(myid) - kfirst_four(myid) + 1
    myfirstlev = kfirst_four(myid)
    mylastlev = klast_four(myid)

  end subroutine ReshapeVerticalGroups


  subroutine RegisterOtherSizes(iMaxPerLat, mPerLat)
    !# Records Other Sizes
    !# @info
    !# **Brief:** Records Other Sizes.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin 
    integer, intent(in) :: iMaxPerLat(jMax)
    integer, intent(in) :: mPerLat(jMax)
    character(len = *), parameter :: h = "**(RegisterOtherSizes)**"
    integer :: MinLatPerBlk
    !$ INTEGER, EXTERNAL ::  OMP_GET_MAX_THREADS
    integer :: i, j, k, l, jk, m, ib, jb, jh, cnt, imp2
    integer :: meanl, meanp, jlast, jused, jfirst
    integer :: mfirst(maxnodes_four), mlast(maxnodes_four)
    integer :: npoints(maxnodes_four)
    integer :: mlast1(maxnodes_four), npoints1(maxnodes_four)
    integer :: maxpoints, maxpointsold, nproc, nlp, nrest
    logical :: improved, done

    allocate (mMaxPerJ(jMax))
    mMaxPerJ = mPerLat
    allocate (iMaxPerJ(jMax))
    iMaxPerJ = iMaxPerLat
    if (iMax < maxval(iMaxPerJ)) then
      stop ' imax and imaxperj disagree'
    end if
    ijMax = sum(iMaxPerJ)

    allocate (firstlatinproc_f(0:maxNodes_four - 1))
    allocate (lastlatinproc_f(0:maxNodes_four - 1))
    allocate (nlatsinproc_f(0:maxNodes_four - 1))
    allocate (nodeHasJ_f(jMax))

    nproc = maxnodes_four

    nlp = jmax / nproc
    nrest = jmax - nproc * nlp
    do i = 1, nproc - nrest
      mfirst(i) = (i - 1) * nlp + 1
      mlast(i) = mfirst(i) + nlp - 1
      npoints(i) = imaxperj(1) * nlp
    enddo
    do i = nproc - nrest + 1, nproc
      mfirst(i) = mlast(i - 1) + 1
      mlast(i) = mfirst(i) + nlp
      npoints(i) = imaxperj(1) * (nlp + 1)
    enddo

    if (any(npoints(:) <= 0)) then
      call fatalError(h, " Too many MPI processes; " // &
        &"there are processes with 0 latitudes")
    end if

    myfirstlat_f = mfirst(myid_four + 1)
    mylastlat_f = mlast(myid_four + 1)
    jMaxlocal_f = 0
    nodeHasJ_f = -1
    do k = 0, maxNodes_four - 1
      firstlatinproc_f(k) = mfirst(k + 1)
      lastlatinproc_f(k) = mlast(k + 1)
      nlatsinproc_f(k) = mlast(k + 1) - mfirst(k + 1) + 1
      do j = firstlatinproc_f(k), lastlatinproc_f(k)
        nodeHasJ_f(j) = k
      enddo
    enddo
    jMaxlocal_f = maxval(nlatsinproc_f)
    myJMax_f = mylastlat_f - myfirstlat_f + 1

    allocate (jMinPerM(mMax))
    allocate (jMaxPerM(mMax))
    jMinPerM = jMaxHalf
    do j = 1, jMaxHalf
      m = mMaxPerJ(j)
      jMinPerM(1:m) = min(j, jMinPerM(1:m))
    end do
    jMaxPerM = jMax - jMinPerM + 1

    ! OpenMP parallelism

    if (dumpLocal) then
      call msgDump(h, ' Dump at the end ')
      !CALL DumpSizes()
    end if
  end subroutine RegisterOtherSizes
  
  subroutine SpectralDomainDecomp()
    !# Domain Decomposition of Fourier Wave Numbers (m's)
    !# @info
    !# **Brief:** Domain Decomposition of Fourier Wave Numbers (m's). </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer :: m
    integer :: n
    integer :: mn
    integer :: mBase
    integer :: mMid
    integer :: lm
    integer :: mm
    integer :: mng
    integer :: k
    integer :: rest
    integer :: is
    integer :: i
    integer :: MaxN, ns, ip, inc, mngiv, np, npl
    integer :: ierr
    integer, allocatable :: ini(:), sends(:)
    character(len = 8) :: c0
    character(len = *), parameter :: h = "**(SpectralDomainDecomp)**"

    ! any MPI process has at most mMaxLocal m's (may have less)
    ! mMaxLocal is used to dimension arrays over all MPI processes

    mMaxLocal = mMax / maxNodes_four
    if (mMax - mMaxLocal * maxNodes_four.ne.0) mMaxLocal = mMaxLocal + 1

    ! msPerProc: how many m's at each MPI process

    allocate(msPerProc    (0:maxNodes_four - 1))

    allocate(msPerProc_out(0:maxNodes_four - 1))

    ! msInProc: which m's are at each MPI process
    ! note that indexing is restricted to
    ! msInProc(1:msPerProc(pId), pId), with pId=0,maxNodes-1

    allocate(msInProc(mMaxLocal, 0:maxNodes_four - 1))
    allocate(msInProc_out(mMaxLocal, 0:maxNodes_four - 1))

    ! nodeHasM: which process has a particular value of m.
    ! Values of nodeHasM for p processes are:
    ! (0, 1, 2, ..., p-1, p-1, ..., 2, 1, 0, 0, 1, 2, ...)
    ! This distribution tries to spread evenly load across
    ! processes, leaving uneven load to the smaller m's.

    allocate(nodeHasM(mMax, ngroups_four))
    allocate(nodeHasM_out(mMax, ngroups_four))

    ! domain decomposition of m's

    do i = 1, ngroups_four
      maxN = npperg_four(i)
      mm = 1
      do mBase = 1, mMax, 2 * maxN
        mMid = mBase + maxN - 1
        do m = mBase, min(mMid, mMax)
          nodeHasM(m, i) = m - mBase
        end do
        mm = mm + 1
        do m = mMid + 1, min(mMid + maxN, mMax)
          nodeHasM(m, i) = maxN + mMid - m
        end do
        mm = mm + 1
      end do
    end do
    if (nodeHasM(1, mygroup_four)==myid_four) then
      haveM1 = .true.
    else
      haveM1 = .false.
    endif
    havesurf = .false.
    if (mygroup_four.eq.1) havesurf = .true.
    !---------
    mm = 1
    msPerProc = 0
    i = mygroup_four
    do mBase = 1, mMax, 2 * maxNodes_four
      mMid = mBase + maxNodes_four - 1
      do m = mBase, min(mMid, mMax)
        msPerProc(nodeHasM(m, i)) = msPerProc(nodeHasM(m, i)) + 1
        msInProc(mm, nodeHasM(m, i)) = m
      end do
      mm = mm + 1
      do m = mMid + 1, min(mMid + maxNodes_four, mMax)
        msPerProc(nodeHasM(m, i)) = msPerProc(nodeHasM(m, i)) + 1
        msInProc(mm, nodeHasM(m, i)) = m
      end do
      mm = mm + 1
    end do
    !---------
    mm = 1
    msPerProc_out = 0
    i = mygroup_four
    do mBase = 1, mMax, 2 * maxNodes_four
      mMid = mBase + maxNodes_four - 1
      do m = mBase, min(mMid, mMax)
        msPerProc_out(nodeHasM(m, i)) = msPerProc_out(nodeHasM(m, i)) + 1
        msInProc_out(mm, nodeHasM(m, i)) = m
      end do
      mm = mm + 1
      do m = mMid + 1, min(mMid + maxNodes_four, mMax)
        msPerProc_out(nodeHasM(m, i)) = msPerProc_out(nodeHasM(m, i)) + 1
        msInProc_out(mm, nodeHasM(m, i)) = m
      end do
      mm = mm + 1
    end do

    !---------

    ! current parallelism restricts the number of MPI processes
    ! to truncation + 1

    if (any(msPerProc <= 0)) then
      call FatalError(h, " Too many MPI processes; " // &
        &"there are processes with 0 Fourier waves")
    end if

    ! myMMax: scalar containing how many m's at this MPI process

    myMMax = msPerProc(myId_four)
    myMMax_out = msPerProc_out(myId_four)

    ! lm2m: maps local indexing of m to global indexing of m,
    ! that is, maps (1:myMMax) to (1:mMax)

    allocate(lm2m(myMMax))
    lm2m(1:mymmax) = msInProc(1:mymmax, myid_four)

    ! DOMAIN DECOMPOSITION OF SPECTRAL COEFFICIENTS (mn's)
    ! all mn's of a single m belongs to a unique process

    ! mnsPerProc: how many mn's at each MPI process

    allocate(mnsPerProc(0:maxNodes_four - 1))

    ! mnsExtPerProc: how many mnExt's at each MPI process

    allocate(mnsExtPerProc(0:maxNodes_four - 1))

    ! domain decomposition of mn's and mnExt's

    mm = 1
    mnsPerProc = 0
    mnsExtPerProc = 0
    i = mygroup_four
    do mBase = 1, mMax, 2 * maxNodes_four
      mMid = mBase + maxNodes_four - 1
      do m = mBase, min(mMid, mMax)
        mnsPerProc(nodeHasM(m, i)) = mnsPerProc(nodeHasM(m, i)) + mmax - m + 1
        mnsExtPerProc(nodeHasM(m, i)) = mnsExtPerProc(nodeHasM(m, i)) + mmax - m + 2
      end do
      mm = mm + 1
      do m = mMid + 1, min(mMid + maxNodes_four, mMax)
        mnsPerProc(nodeHasM(m, i)) = mnsPerProc(nodeHasM(m, i)) + mmax - m + 1
        mnsExtPerProc(nodeHasM(m, i)) = mnsExtPerProc(nodeHasM(m, i)) + mmax - m + 2
      end do
      mm = mm + 1
    end do

    ! any MPI process has at most mnMaxLocal mn's (may have less) and at most
    ! mnExtMaxLocal mnExt's (may have less).
    ! mnMaxLocal and mnExtMaxLocal are used to dimension arrays over all MPI processes

    mnMaxLocal = maxval(mnsPerProc)
    mnExtMaxLocal = maxval(mnsExtPerProc)

    ! myMNMax: scalar containing how many mn's at this MPI process
    ! myMNExtMax: scalar containing how many mnExt's at this MPI process

    myMNMax = mnsPerProc(myId_four)
    myMNExtMax = mnsExtPerProc(myId_four)

    ! MAPPINGS OF LOCAL INDICES mn TO (localm,n)

    ! Mapping Local pairs (lm,n) to 1D for Regular Spectral:
    ! (1) myMMap(mn): which lm is stored at this position
    ! (2) myNMap(mn): which  n is stored at this position
    ! (3) myMNMap(lm,n): position storing pair (lm,n)

    allocate (myMNMap(myMMax, nMax))
    allocate (myMMap(myMNMax))
    allocate (myNMap(myMNMax))
    myMNMap = -1  ! flag mapping error
    myMMap = -1   ! flag mapping error
    myNMap = -1   ! flag mapping error
    mn = 0
    do lm = 1, myMMax
      do n = lm2m(lm), mMax
        mn = mn + 1
        myMNMap(lm, n) = mn
        myMMap(mn) = lm
        myNMap(mn) = n
      end do
    end do

    ! Mapping Local pairs (lm,n) to 1D for Extended Spectral:
    ! (1) myMExtMap(mn): which lm is stored at this position
    ! (2) myNExtMap(mn): which  n is stored at this position
    ! (3) myMNExtMap(lm,n): position storing pair (lm,n)

    allocate (myMExtMap(myMNExtMax))
    allocate (myNExtMap(myMNExtMax))
    allocate (myMNExtMap(myMMax, nExtMax))
    myMExtMap = -1  ! flag mapping error
    myNExtMap = -1  ! flag mapping error
    myMNExtMap = -1 ! flag mapping error
    mn = 0
    do lm = 1, myMMax
      do n = lm2m(lm), mMax + 1
        mn = mn + 1
        myMNExtMap(lm, n) = mn
        myMExtMap(mn) = lm
        myNExtMap(mn) = n
      end do
    end do

    !   spectral communicators for surface field replication

    if (ngroups_four.gt.1) then
      if (mygroup_four.eq.1) then
        allocate (comm_spread(maxNodes, 2))
        allocate (ms_spread(mymmax, ngroups_four))
        allocate (ini(0:nprocmax_four))
        comm_spread = 0
        mng = 0
        do n = 2, ngroups_four
          ini = 0
          do m = 1, mymmax
            npl = nodeHasM(lm2m(m), n)
            np = map_four(n, npl)
            if (ini(npl).eq.0) then
              mng = mng + 1
              ini(npl) = mng
            endif
            mn = ini(npl)
            comm_spread(mn, 1) = np
            comm_spread(mn, 2) = comm_spread(mn, 2) + MMax + 1 - lm2m(m)
            ms_spread(m, n) = mn
          enddo
        enddo
        comm_spread(:, 2) = 2 * comm_spread(:, 2)
        ncomm_spread = mng
      else
        allocate (comm_spread(npperg_four(1), 2))
        allocate (ms_spread(mymmax, 1))
        allocate (ini(0:nprocmax_four))
        comm_spread = 0
        mng = 0
        ini = 0
        do m = 1, mymmax
          npl = nodeHasM(lm2m(m), 1)
          np = map_four(1, npl)
          if (ini(npl).eq.0) then
            mng = mng + 1
            ini(npl) = mng
          endif
          mn = ini(npl)
          comm_spread(mn, 1) = np
          comm_spread(mn, 2) = comm_spread(mn, 2) + MMax + 1 - lm2m(m)
          ms_spread(m, 1) = mn
        enddo
        comm_spread(:, 2) = 2 * comm_spread(:, 2)
        ncomm_spread = mng
      endif
    endif

  end subroutine SpectralDomainDecomp


  subroutine ThreadDecomp(firstInd, lastInd, minInd, maxInd, msg)
    !# Thread Decomposition
    !# @info
    !# **Brief:** Thread Decomposition. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: firstInd
    integer, intent(in) :: lastInd
    integer, intent(OUT) :: minInd
    integer, intent(OUT) :: maxInd
    character(len = *), intent(in) :: msg
    integer :: chunk
    integer :: left
    integer :: nTrd
    integer :: iTrd
    integer :: length
    logical :: inParallel
    logical :: op
    character(len = *), parameter :: h = "**(ThreadDecomp)**"
    !$ INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
    !$ INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
    !$ LOGICAL, EXTERNAL :: OMP_IN_PARALLEL

    inParallel = .false.
    nTrd = 1
    iTrd = 0
    !$ inParallel = OMP_IN_PARALLEL()
    if (inParallel) then
      !$ nTrd = OMP_GET_NUM_THREADS()
      !$ iTrd = OMP_GET_THREAD_NUM()
      length = lastInd - firstInd + 1
      chunk = length / nTrd
      left = length - chunk * nTrd
      if (iTrd < left) then
        minInd = iTrd * (chunk + 1) + firstInd
        maxInd = (iTrd + 1) * (chunk + 1) + firstInd - 1
      else
        minInd = iTrd * (chunk) + left + firstInd
        maxInd = (iTrd + 1) * (chunk) + left + firstInd - 1
      end if
    else
      minInd = firstInd
      maxInd = lastInd
    end if

    if (dumpLocal) then

      ! since using unitDump directly (instead of using msgDump),
      ! check if unitDump in open

      inquire(getUnitDump(), opened = op)
      if (.not. op) then
        call FatalError(h, " unitDump not opened; CreateParallelism not invoked")
      end if

      if (inParallel) then
        write(getUnitDump(), "(a,' thread ',i2,' got [',i8,':',i8,&
          &'] from [',i8,':',i8,'] in parallel region at ',a)") &
          h, iTrd, minInd, maxInd, firstInd, lastInd, msg
      else
        write(getUnitDump(), "(a,' kept domain [',i8,':',i8,&
          &'] since not in parallel region at ',a)") &
          h, minInd, maxInd, msg
      end if
    end if
  end subroutine ThreadDecomp


  subroutine ThreadDecompms(m, myms, nms)
    !# Thread Decomposition ms
    !# @info
    !# **Brief:** Thread Decomposition ms. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: m
    integer, intent(inout) :: myms(m)
    integer, intent(OUT) :: nms
    integer :: i, j, k
    integer :: nTrd
    integer :: iTrd
    logical :: inParallel
    !$ INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
    !$ INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
    !$ LOGICAL, EXTERNAL :: OMP_IN_PARALLEL

    inParallel = .false.
    nTrd = 1
    iTrd = 0
    !$ inParallel = OMP_IN_PARALLEL()
    if (inParallel) then
      !$ nTrd = OMP_GET_NUM_THREADS()
      !$ iTrd = OMP_GET_THREAD_NUM() + 1
      nms = 0
      i = 1
      j = 1
      do k = 1, m
        if (i.eq.iTrd) then
          nms = nms + 1
          myms(nms) = k
        endif
        i = i + j
        if (i.eq.nTrd + 1) then
          j = -1
          i = nTrd
        endif
        if (i.eq.0) then
          j = 1
          i = 1
        endif
      enddo
    else
      do k = 1, m
        myms(k) = k
      enddo
      nms = m
    end if

  end subroutine ThreadDecompms

  subroutine GridDecomposition(ibmax, jbmax, nproc)
    !# Grid Decomposition
    !# @info
    !# **Brief:** Grid Decomposition. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    integer, intent(in) :: ibmax
    integer, intent(OUT) :: jbmax
    integer, intent(in) :: nproc
    integer :: np2, np1, np, ngroups, npperg, nrest, npr
    integer :: ngptotal, lat, usedinlat, next, ndim, jbdim
    integer :: ngpperproc, n, imp2, ndim_f, iovmax, iadd
    integer :: ipar, i1, i2, il1, il2, j1, j2, ij, ij1, ind, aux(4)
    integer :: firstextlat, lastextlat, iex, mygridpoints
    integer :: i, k, iproc, j, jlo(1), ib, jb, ijb, jh
    integer, allocatable :: lon(:), iovlap(:)
    integer, allocatable :: nprocsingroup(:)
    integer, allocatable :: npointsingroup(:)
    integer, allocatable :: firstlatingroup(:)
    integer, allocatable :: lastlatingroup(:)
    integer, allocatable :: firstloningroup(:)
    integer, allocatable :: lastloningroup(:)
    integer, allocatable :: jbmaxingroup(:)
    integer, allocatable :: procingroup(:)
    integer, allocatable :: iaux(:, :)
    real (kind = p_r8) :: hj, ifirst, ilast
    character(len = *), parameter :: h = "**(GridDecomposition)**"

    if (nproc.le.3) then
      ngroups = nproc
      npperg = 1
      nrest = 0
      iex = 0
    else
      if (nproc.le.5) then
        ngroups = 1
        npperg = nproc - 1
        nrest = 0
        iex = 1
      else
        np2 = nproc / 2
        np = sqrt(real(np2))
        if (np * (np + 1).lt.np2) then
          ngroups = np + 1
        else
          ngroups = np
        endif
        npperg = nproc / ngroups
        nrest = nproc - npperg * ngroups
        if (npperg.le.8) then
          iex = 1
        else
          iex = 2
        endif
      endif
    endif
    hj = acos(-1._p_r8) / jmax
    allocate (gridmap(1:iMax, 1:jMax))
    allocate (procingroup(0:nproc - 1))
    allocate (nprocsingroup(ngroups + 2 * iex))
    nprocsingroup(1) = 1
    nprocsingroup(ngroups + 2 * iex) = 1
    if (iex.eq.2) then
      npr = sqrt(real(npperg))
      nprocsingroup(iex) = npr + 1
      nprocsingroup(iex + 1) = npperg - npr - 2
      nprocsingroup(ngroups + 3) = npr + 1
      nprocsingroup(ngroups + iex) = npperg - npr - 2
    else
      if (iex.eq.1) then
        nprocsingroup(2) = npperg - 1
        nprocsingroup(ngroups + 1) = npperg - 1
      endif
    endif
    nprocsingroup(2 + iex:ngroups + iex - 1) = npperg
    do k = 1, nrest
      nprocsingroup(1 + iex + k) = nprocsingroup(1 + iex + k) + 1
    enddo
    ngroups = ngroups + 2 * iex
    i = -1
    do k = 1, ngroups
      procingroup(i + 1:i + nprocsingroup(k)) = k
      i = i + nprocsingroup(k)
    enddo
    allocate (npointsingroup(ngroups))
    allocate (pointsinproc(0:nproc - 1))
    allocate (firstlatingroup(ngroups))
    allocate (firstloningroup(ngroups))
    allocate (lastlatingroup(ngroups))
    allocate (lastloningroup(ngroups))
    allocate (jbmaxingroup(ngroups))

    ngptotal = sum(imaxperj(1:jmax))
    ngpperproc = ngptotal / nproc
    nrest = ngptotal - ngpperproc * nproc
    pointsinproc(0:nrest - 1) = ngpperproc + 1
    pointsinproc(nrest:nproc - 1) = ngpperproc
    lat = 1
    usedinlat = 0
    iproc = 0
    do k = 1, ngroups
      next = min(nrest, nprocsingroup(k))
      npointsingroup(k) = ngpperproc * nprocsingroup(k) + next
      nrest = nrest - next
      firstlatingroup(k) = lat
      firstloningroup(k) = usedinlat + 1
      np = imaxperj(lat) - usedinlat
      do
        if (np.ge.npointsingroup(k)) EXIT
        lat = lat + 1
        np = np + imaxperj(lat)
      enddo
      lastlatingroup(k) = lat
      usedinlat = imaxperj(lat) - np + npointsingroup(k)
      lastloningroup(k) = usedinlat
      if (usedinlat.eq.imaxperj(lat)) then
        lat = lat + 1
        usedinlat = 0
      endif
    enddo
    jbmaxingroup = lastlatingroup - firstlatingroup + 1
    jbdim = maxval(jbmaxingroup)
    jovlap = 0

    ndim = 4 * jbdim + 2 * jovlap * npperg
    allocate (mysendsgr(4, ndim))
    allocate (mysendspr(2, nproc))
    allocate (myrecsgr(4, ndim))
    allocate (myrecspr(2, nproc))
    allocate (firstlat(0:nproc - 1))
    allocate (firstlon(jbdim, 0:nproc - 1))
    allocate (lastlat(0:nproc - 1))
    allocate (lastlon(jbdim, 0:nproc - 1))
    allocate (lon(jbdim))

    iproc = 0
    do k = 1, ngroups
      firstlon(1, iproc) = firstloningroup(k)
      firstlon(2:jbmaxingroup(k), iproc) = 1
      lastlon(1:jbmaxingroup(k), iproc) = 0
      lon(1:jbmaxingroup(k)) = firstlon(1:jbmaxingroup(k), iproc)
      do n = 1, nprocsingroup(k)
        do np = 1, pointsinproc(iproc)
          jlo = minloc(real(lon(1:jbmaxingroup(k))) &
            / real(imaxperj(firstlatingroup(k):lastlatingroup(k))))
          j = jlo(1)
          lastlon(j, iproc) = lon(j)
          gridmap(lon(j), j + firstlatingroup(k) - 1) = iproc
          lon(j) = lon(j) + 1
          if (j.eq.jbmaxingroup(k).and.lon(j).gt.lastloningroup(k)) &
            lon(j) = imaxperj(lastlatingroup(k)) + 1
        enddo
        if (lastlon(1, iproc).eq.0) then
          firstlat(iproc) = firstlatingroup(k) + 1
        else
          firstlat(iproc) = firstlatingroup(k)
        endif
        if (lastlon(jbmaxingroup(k), iproc).eq.0) then
          lastlat(iproc) = lastlatingroup(k) - 1
        else
          lastlat(iproc) = lastlatingroup(k)
        endif
        iproc = iproc + 1
        if (iproc.ne.nproc) then
          firstlon(1:jbmaxingroup(k), iproc) = lon(1:jbmaxingroup(k))
          lastlon(1:jbmaxingroup(k), iproc) = 0
        endif
      enddo
    enddo
    mygridpoints = pointsinproc(myid)

    jbmax = mygridpoints / ibmax
    if(jbmax * ibmax.lt.mygridpoints) jbmax = jbmax + 1
    allocate (iPerIJB(ibmax, jbmax));iPerIJB = 0
    allocate (jPerIJB(ibmax, jbmax));jPerIJB = 0
    allocate (ibMaxPerJB(jbmax))   ;ibMaxPerJB = 0

    myfirstlat = firstlat(myid)
    mylastlat = lastlat(myid)
    firstextlat = max(myfirstlat - jovlap, 1)
    lastextlat = min(mylastlat + jovlap, jmax)

    allocate (myfirstlon(firstextlat:lastextlat))
    allocate (mylastlon(firstextlat:lastextlat))
    allocate (iovlap(firstextlat:lastextlat))

    do j = firstextlat, lastextlat
      i = nint(jovlap * imaxperj(j) / (imax * sin((j - .5_p_r8) * hj)))
      iovlap(j) = i + 1
    enddo
    iovmax = 0

    allocate (ibPerIJ(1 - iovmax:iMax + iovmax, -1:jMax + 2))
    allocate (jbPerIJ(1 - iovmax:iMax + iovmax, -1:jMax + 2))
    ibPerIJ = 0
    jbPerIJ = 0

    myfirstlon(myfirstlat:mylastlat) = &
      firstlon(1 + myfirstlat - firstlatingroup(procingroup(myid)):&
        jbmaxingroup(procingroup(myid)) - lastlatingroup(procingroup(myid)) &
          + mylastlat, myid)
    mylastlon(myfirstlat:mylastlat) = &
      lastlon(1 + myfirstlat - firstlatingroup(procingroup(myid)):&
        jbmaxingroup(procingroup(myid)) - lastlatingroup(procingroup(myid)) &
          + mylastlat, myid)
    ifirst = minval(real(myfirstlon(myfirstlat:mylastlat)) &
      / real(imaxperj(myfirstlat:mylastlat)))
    ilast = maxval(real(mylastlon(myfirstlat:mylastlat)) &
      / real(imaxperj(myfirstlat:mylastlat)))

    !  interior domain
    !  ---------------
    ijb = 0
    do j = myfirstlat, mylastlat
      do i = myfirstlon(j), mylastlon(j)
        ib = mod(ijb, ibmax) + 1
        jb = ijb / ibmax + 1
        ijb = ijb + 1
        iPerIJB(ib, jb) = i
        jPerIJB(ib, jb) = j
        ibPerIJ(i, j) = ib
        jbPerIJ(i, j) = jb
      enddo
    enddo
    ibMaxPerJB(1:jbmax - 1) = ibmax
    ibMaxPerJB(jbmax) = ib

    ! define messages to be exchanged between fourier and grid computations
    ! ---------------------------------------------------------------------
    ndim_f = jMaxlocal_f * (maxval(nprocsingroup) + 3)
    allocate (messages_f(4, ndim_f))
    allocate (messproc_f(2, 0:nproc))
    allocate (messages_g(4, ndim_f))
    allocate (messproc_g(2, 0:nproc))
    ipar = 1
    do j = myfirstlat_f, mylastlat_f
      messages_f(1, ipar) = 1
      messages_f(3, ipar) = j
      messages_f(4, ipar) = gridmap(1, j)
      do i = 2, imaxperj(j)
        if (gridmap(i, j).ne.messages_f(4, ipar)) then
          messages_f(2, ipar) = i - 1
          if(messages_f(4, ipar).ne.myid) ipar = ipar + 1
          messages_f(1, ipar) = i
          messages_f(3, ipar) = j
          messages_f(4, ipar) = gridmap(i, j)
        endif
        if (i.eq.imaxperj(j)) then
          messages_f(2, ipar) = i
          if(messages_f(4, ipar).ne.myid) ipar = ipar + 1
        endif
      enddo
    enddo
    ipar = ipar - 1
    if (ipar.gt.ndim_f) then
      write(p_nfprt, *) ' ndim_f, ipar  ', ndim_f, ipar
      write(p_nfprt, "(a, ' dimensioning of segment messages insufficient')") h
      stop h
    endif

    !sort messages by processors
    !---------------------------
    do i = 2, ipar
      do j = i, 2, -1
        if(messages_f(4, j).lt.messages_f(4, j - 1)) then
          aux = messages_f(:, j - 1)
          messages_f(:, j - 1) = messages_f(:, j)
          messages_f(:, j) = aux
        else
          EXIT
        endif
      enddo
    enddo
    messproc_f(2, 0) = 0
    if (ipar.gt.0) then
      messproc_f(1, 1) = messages_f(4, 1)
      nrecs_f = 1
      do i = 2, ipar
        if(messages_f(4, i).ne.messages_f(4, i - 1)) then
          nrecs_f = nrecs_f + 1
          messproc_f(2, nrecs_f - 1) = i - 1
          messproc_f(1, nrecs_f) = messages_f(4, i)
        endif
      enddo
      messproc_f(2, nrecs_f) = ipar
    else
      nrecs_f = 0
    endif

    !  set communication structure for grid diagnostics
    !  ------------------------------------------------
    allocate(firstandlastlat(2, 0:maxNodes - 1))
    allocate (nlatsinproc_d(0:maxNodes - 1))
    lat = jmax / maxNodes
    nrest = jmax - lat * maxNodes
    np = nrest / 2
    np1 = nrest - np
    n = 0
    do i = 0, maxNodes - 1
      firstandlastlat(1, i) = n + 1
      if (i.lt.np.or.i.ge.maxNodes - np1) then
        n = n + lat + 1
      else
        n = n + lat
      endif
      firstandlastlat(2, i) = n
    enddo
    nlatsinproc_d = firstandlastlat(2, :) - firstandlastlat(1, :) + 1
    myfirstlat_diag = firstandlastlat(1, myid)
    mylastlat_diag = firstandlastlat(2, myid)
    myJMax_d = mylastlat_diag - myfirstlat_diag + 1

    ndim = (5 + maxval(nprocsingroup)) * (lat + 1)
    allocate (mysends_diag(4, ndim))
    allocate (mysendspr_diag(2, 0:nproc))
    allocate (myrecs_diag(4, ndim))
    allocate (myrecspr_diag(2, 0:nproc))

    ij = 0
    do j = myfirstlat_diag, mylastlat_diag
      ij = ij + 1
      myrecs_diag(1, ij) = 1
      myrecs_diag(3, ij) = j
      myrecs_diag(4, ij) = gridmap(1, j)
      do i = 2, imaxperj(j)
        if (gridmap(i, j).ne.myrecs_diag(4, ij)) then
          myrecs_diag(2, ij) = i - 1
          ij = ij + 1
          myrecs_diag(1, ij) = i
          myrecs_diag(3, ij) = j
          myrecs_diag(4, ij) = gridmap(i, j)
        endif
      enddo
      myrecs_diag(2, ij) = imaxperj(j)
    enddo

    !sort messages by processors
    !---------------------------
    do i = 2, ij
      do j = i, 2, -1
        if(myrecs_diag(4, j).lt.myrecs_diag(4, j - 1)) then
          aux = myrecs_diag(:, j - 1)
          myrecs_diag(:, j - 1) = myrecs_diag(:, j)
          myrecs_diag(:, j) = aux
        else
          EXIT
        endif
      enddo
    enddo
    n = 0
    i1 = 0
    do i = 1, ij
      if(myrecs_diag(4, i).eq.myid) then
        i1 = i
        n = n + 1
      endif
    enddo
    ij = ij - n
    do i = i1 - n + 1, ij
      myrecs_diag(:, i) = myrecs_diag(:, i + n)
    enddo
    myrecspr_diag(2, 0) = 0
    if (ij.gt.0) then
      myrecspr_diag(1, 1) = myrecs_diag(4, 1)
      nrecs_diag = 1
      do i = 2, ij
        if(myrecs_diag(4, i).ne.myrecs_diag(4, i - 1)) then
          nrecs_diag = nrecs_diag + 1
          myrecspr_diag(2, nrecs_diag - 1) = i - 1
          myrecspr_diag(1, nrecs_diag) = myrecs_diag(4, i)
        endif
      enddo
      myrecspr_diag(2, nrecs_diag) = ij
    else
      nrecs_diag = 0
    endif

  end subroutine GridDecomposition

  subroutine Clear_Sizes()
    !# Cleans Sizes
    !# @info
    !# **Brief:** Cleans Sizes. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2010 <br>
    !# @endin
    deallocate (mMap)
    deallocate (nMap)
    deallocate (mnMap)
    deallocate (mnMap_out)
    deallocate (mExtMap)
    deallocate (nExtMap)
    deallocate (mnExtMap)
    deallocate (nodeHasM_out)
    deallocate (nodeHasM)
    deallocate (lm2m)
    deallocate (msInProc)
    deallocate (msInProc_out)
    deallocate (msPerProc)
    deallocate (msPerProc_out)

    deallocate (myMMap)
    deallocate (myNMap)
    deallocate (mnsPerProc)
    deallocate (mnsExtPerProc)
    deallocate (myMNMap)
    deallocate (myMExtMap)
    deallocate (myNExtMap)
    deallocate (myMNExtMap)

    deallocate (comm_spread)
    deallocate (ms_spread)

    deallocate (mMaxPerJ)
    deallocate (jMinPerM)
    deallocate (jMaxPerM)

    deallocate (iMaxPerJ)

    deallocate (ibMaxPerJB)
    deallocate (firstlatinproc_f)
    deallocate (lastlatinproc_f)
    deallocate (nlatsinproc_f)
    deallocate (nodeHasJ_f)
    deallocate (kfirst_four)
    deallocate (klast_four)
    deallocate (nlevperg_four)
    deallocate (npperg_four)
    deallocate (map_four)
    deallocate (first_proc_four)
    deallocate (nlatsinproc_d)
    deallocate (messages_f)
    deallocate (messages_g)
    deallocate (messproc_f)
    deallocate (messproc_g)
    deallocate (nset)
    !
    ! grid decomposition
    !
    deallocate (firstandlastlat)
    deallocate (myfirstlon)
    deallocate (mylastlon)
    deallocate (firstlon)
    deallocate (lastlon)
    deallocate (firstlat)
    deallocate (lastlat)
    deallocate (ibPerIJ)
    deallocate (jbPerIJ)
    deallocate (iPerIJB)
    deallocate (jPerIJB)
    deallocate (pointsinproc)
    deallocate (myrecsgr)
    deallocate (myrecspr)
    deallocate (mysendsgr)
    deallocate (mysendspr)
    deallocate (myrecs_diag)
    deallocate (myrecspr_diag)
    deallocate (mysends_diag)
    deallocate (mysendspr_diag)
    deallocate (gridmap)

    ! DEALLOCATE ( si  )      ! sigma
    ! DEALLOCATE ( del  )     ! layer thickness (in sigma)
    ! DEALLOCATE ( delcl  )     ! layer thickness (in cl)
    ! DEALLOCATE ( rdel2  )
    ! DEALLOCATE ( sl  )      ! sigma at layer midpoint
    ! DEALLOCATE ( cl  )      ! 1.0 - sl
    ! DEALLOCATE ( rpi  )     ! 'pi' ratios at adjacent layers

    mMax = -1
    mMax_in = -1
    nMax = -1
    mnMax = -1
    mnMax_out = -1
    nExtMax = -1
    mnExtMax = -1
    myMMax = -1
    myMMax_out = -1
    mMaxLocal = -1
    mnMaxLocal = -1
    mnExtMaxLocal = -1
    myMNMax = -1
    myMNExtMax = -1
    HaveM1 = .true.

    ncomm_spread = -1

    jMax = -1
    jMaxHalf = -1

    iMax = -1
    ijMax = -1
    ijMaxGauQua = -1

    myfirstlat_f = -0
    mylastlat_f = -0
    myfirstlev = -0
    mylastlev = -0
    myJMax_f = -0
    ngroups_four = -0
    nprocmax_four = -0
    JMaxlocal_f = -0
    kMaxloc = -0
    kMaxloc_in = -0
    kMaxloc_out = -0
    nrecs_f = -0
    nrecs_g = -0
    havesurf = .false.
    !
    ! grid decomposition
    !
    ibMax = -1
    jbMax = -1
    jbMax_ext = -0
    myfirstlat = -0
    mylastlat = -0
    nrecs_gr = -0
    nsends_gr = -0
    nrecs_diag = -0
    nsends_diag = -0
    myfirstlat_diag = -0
    mylastlat_diag = -0
    myJMax_d = -0
    jovlap = -0

    kMax = -1

  end subroutine Clear_Sizes

end module Mod_Sizes

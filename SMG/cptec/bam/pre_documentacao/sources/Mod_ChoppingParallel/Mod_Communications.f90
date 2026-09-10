!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Communications </br></br>
!#
!# **Brief**: Module used for MPI communication </br></br>
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
!#  <li>07-04-2011 - Paulo Kubota   - version: 1.0.0 </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 2.0.0 - some adaptations for modularizing Chopping </li>
!#  <li>14-02-2020 - Eduardo Khamis - version: 2.0.0 </li>
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


module Mod_Communications

  use Mod_Parallelism_Group_Chopping, only : &
    myId &
    , maxNodes &
    , mpiCommGroup

  use Mod_Parallelism_Fourier, only : &
    maxNodes_four &
    , myId_four &
    , COMM_FOUR &
    , mygroup_four
    
  !  USE Dumpgraph, ONLY:   &
  !       dumpgra, Writectl

  use Mod_Sizes, only : &
    mymmax, &
    mymnmax, &
    mymnextmax, &
    kmax, &
    kmaxloc, &
    mmax, &
    mmap, &
    mnmax_out, &
    mnextmax, &
    mnmaxlocal, &
    mnextmaxlocal, &
    ibmax, &
    jbmax, &
    jbMax_ext, &
    ijmax, &
    imax, &
    jmax, &
    ibmaxperjb, &
    imaxperj, &
    jbperij, &
    ibperij, &
    Msperproc, &
    Msinproc, &
    mnsPerProc, &
    mnsExtPerProc, &
    NodehasM, &
    lm2m, &
    myfirstlat, &
    mylastlat, &
    mysendsgr, &
    mysendspr, &
    myrecsgr, &
    myrecspr, &
    firstlat, &
    lastlat, &
    firstlon, &
    lastlon, &
    myfirstlon, &
    mylastlon, &
    messages_f, &
    messproc_f, &
    messages_g, &
    messproc_g, &
    nrecs_diag, &
    nsends_diag, &
    myfirstlat_diag, &
    mylastlat_diag, &
    myjmax_d, &
    firstandlastlat, &
    myrecs_diag, &
    myrecspr_diag, &
    mysends_diag, &
    mysendspr_diag, &
    havesurf, &
    myfirstlev, &
    map_four, &
    ngroups_four, &
    kfirst_four, &
    klast_four, &
    nlevperg_four, &
    first_proc_four, &
    ncomm_spread, &
    comm_spread, &
    ms_spread, &
    nlatsinproc_d, &
    gridmap, &
    pointsinproc

  use Mod_Utils, only : &
    CyclicNearest_r, &
    CyclicLinear_ABS, &
    CyclicLinear

  implicit none
  include 'mpif.h'
  include 'precision.h'
  include 'messages.h'

  private

  public :: Collect_Grid_Red
  public :: Collect_Grid_Sur
  public :: Collect_Grid_His
  public :: Collect_Grid_Full
  public :: Collect_Grid_d
  public :: Collect_Gauss
  public :: Collect_Spec
  public :: Collect_Spec_Ext
  public :: Exchange_ftog
  public :: Exchange_diag
  public :: Set_Communic_buffer
  public :: Spread_surf_Spec
  public :: p2d
  public :: Clear_Communications
  real(kind = p_r8), public, allocatable :: bufrec(:)
  real(kind = p_r8), public, allocatable :: bufsend(:)
  integer, allocatable :: isbrec(:)
  integer, allocatable :: isbsend(:)
  integer, allocatable :: ilrecbuf(:)
  integer, allocatable :: ilsendbuf(:)
  integer, public :: dimrecbuf
  integer, public :: dimsendbuf
  TYPE p2d
    real(kind = p_r8), pointer :: p(:, :)
  end TYPE p2d

contains

  subroutine Set_Communic_buffer
    !# Sets Communication with buffer
    !# ---
    !# @info
    !# **Brief:** Sets Communication with buffer. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endinfo

    integer :: ndim

    ndim =  kmax * jmax * mmax * 2_p_r8 * 8_p_r8 / maxnodes

    dimrecbuf = ndim
    dimsendbuf = ndim
    allocate (bufrec(dimrecbuf))
    allocate (bufsend(dimsendbuf))

  end subroutine Set_Communic_buffer

  subroutine Collect_Grid_Red(field, fieldglob)
    !# Collects Grid Red
    !# ---
    !# @info
    !# **Brief:** Collects Grid Red. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    ! 
    !   Processor 0 has output in fieldglob
    ! 
    real(kind = p_r8), intent(in) :: field(ibMax * jbMax)
    real(kind = p_r8), intent(OUT) :: fieldglob(ijmax)

    integer :: ij, j, i, ii
    integer :: comm
    integer :: ierr
    integer :: index
    integer :: request
    integer :: requestr(0:MaxNodes - 1)
    integer :: ini(0:MaxNodes - 1)
    integer :: status(MPI_STATUS_SIZE)

    comm = mpiCommGroup
    if (myid.ne.0) then
      call MPI_ISEND(field, pointsinproc(myid), MPI_DOUBLE_PRECISION, 0, &
        91, comm, request, ierr)
      call MPI_WAIT(request, status, ierr)
    else
      requestr(0) = MPI_REQUEST_NULL
      ini(0) = 1
      ij = 1 + pointsinproc(0)
      do ii = 1, MaxNodes - 1
        ini(ii) = ij
        call MPI_IRECV(bufrec(ij), pointsinproc(ii), MPI_DOUBLE_PRECISION, ii, 91, &
          comm, requestr(ii), ierr)
        ij = ij + pointsinproc(ii)
      enddo
      bufrec(1:pointsinproc(0)) = field(1:pointsinproc(0))
      do ii = 1, MaxNodes - 1
        call MPI_WAITANY(MaxNodes - 1, requestr(1), index, status, ierr)
      end do
      ii = 1
      do j = 1, jmax
        do i = 1, imaxperj(j)
          ij = gridmap(i, j)
          fieldglob(ii) = bufrec(ini(ij))
          ii = ii + 1
          ini(ij) = ini(ij) + 1
        end do
      end do
    end if

  end subroutine Collect_Grid_Red


  subroutine Collect_Grid_His(field, fieldglob, ngpts, ngptslocal, nproc, nf, &
    ngptsperproc, mapglobal)
    !# Collects Grid His
    !# ---
    !# @info
    !# **Brief:** Collects Grid His </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nproc 
    !# destination processor
    integer, intent(in) :: nf
    integer, intent(in) :: ngpts
    integer, intent(in) :: ngptslocal
    integer, intent(in) :: mapglobal(ngpts)
    real(kind = p_r8), intent(in) :: field(ngptslocal, nf)
    integer, intent(in) :: ngptsperproc(0:maxnodes - 1)
    real(kind = p_r8), intent(OUT) :: fieldglob(ngpts, nf)

    integer :: ij, i, n, i1, i2
    integer :: comm
    integer :: ierr
    integer :: index
    integer :: request
    integer :: requestr(0:MaxNodes - 1)
    integer :: ini(0:MaxNodes)
    integer :: status(MPI_STATUS_SIZE)

    if (dimrecbuf.lt.ngpts * nf) then
      dimrecbuf = ngpts * nf
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    comm = mpiCommGroup
    if (myid.ne.nproc) then
      call MPI_ISEND(field, ngptslocal * nf, MPI_DOUBLE_PRECISION, nproc, 92, comm, request, ierr)
      call MPI_WAIT(request, status, ierr)
    else
      requestr(nproc) = MPI_REQUEST_NULL
      ini(0) = 0
      ij = 1
      do i = 0, MaxNodes - 1
        if (i.ne.nproc) then
          call MPI_IRECV(bufrec(ij), ngptsperproc(i) * nf, MPI_DOUBLE_PRECISION, i, 92, &
            comm, requestr(i), ierr)
        endif
        ini(i + 1) = ini(i) + ngptsperproc(i)
        ij = ij + ngptsperproc(i) * nf
      enddo
      i1 = ini(nproc) + 1
      i2 = ini(nproc + 1)
      fieldglob(mapglobal(i1:i2), :) = field(1:ngptslocal, :)
      do i = 1, MaxNodes - 1
        call MPI_WAITANY(MaxNodes, requestr(0), index, status, ierr)
        ij = status(MPI_SOURCE)
        i1 = ini(ij) * nf
        !CDIR NODEP
        do n = 1, nf
          fieldglob(mapglobal(ini(ij) + 1:ini(ij + 1)), n) = &
            bufrec(i1 + 1:i1 + ngptsperproc(ij))
          i1 = i1 + ngptsperproc(ij)
        enddo
      enddo
    end if

  end subroutine Collect_Grid_His


  subroutine Collect_Grid_Sur(field, fieldglob, nproc)
    !# Collects Grid Surface Fields
    !# ---
    !# @info
    !# **Brief:** Collects grid surface fields </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nproc 
    ! destination processor
    real(kind = p_r8), intent(in) :: field(imax, myjMax_d)
    real(kind = p_r8), intent(OUT) :: fieldglob(imax, jmax)

    integer :: ij, i
    integer :: comm
    integer :: ierr
    integer :: request
    integer :: requestr(0:MaxNodes - 1)
    integer :: status(MPI_STATUS_SIZE, maxnodes)

    comm = mpiCommGroup
    if (myid.ne.nproc) then
      if (myjmax_d.gt.0) then
        call MPI_ISEND(field, imax * myjmax_d, MPI_DOUBLE_PRECISION, nproc, 93, comm, request, ierr)
        call MPI_WAIT(request, status, ierr)
      endif
    else
      requestr = MPI_REQUEST_NULL
      if(myjmax_d.gt.0) fieldglob(:, myfirstlat_diag:mylastlat_diag) = field(:, :)
      ij = 1
      do i = 0, MaxNodes - 1
        if (i.ne.nproc.and.nlatsinproc_d(i).gt.0) then
          call MPI_IRECV(fieldglob(1, ij), nlatsinproc_d(i) * imax, MPI_DOUBLE_PRECISION, i, 93, &
            comm, requestr(i), ierr)
        endif
        ij = ij + nlatsinproc_d(i)
      enddo
      call MPI_WAITALL(MaxNodes, requestr(0), status, ierr)
    end if

  end subroutine Collect_Grid_Sur


  subroutine Collect_Grid_d(field, fieldglob, levs, nproc)
    !# Collects Grid d
    !# ---
    !# @info
    !# **Brief:** Collects Grid d </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    
    integer, intent(in) :: nproc 
    !# destination processor
    integer, intent(in) :: levs
    real(kind = p_r8), intent(in) :: field(imax, myjMax_d, levs)
    real(kind = p_r8), intent(OUT) :: fieldglob(imax, jmax, levs)

    integer :: ij, i, nr, k, j, ip
    integer :: comm
    integer :: ierr
    integer :: request
    integer :: requestr(0:MaxNodes - 1)
    integer :: ini(0:MaxNodes)
    integer :: index
    integer :: status(MPI_STATUS_SIZE)

    comm = mpiCommGroup
    if (myid.ne.nproc) then
      if (myjmax_d.gt.0) then
        call MPI_ISEND(field, imax * myjmax_d * levs, MPI_DOUBLE_PRECISION, nproc, 93, comm, request, ierr)
        call MPI_WAIT(request, status, ierr)
      endif
    else
      requestr = MPI_REQUEST_NULL
      if(myjmax_d.gt.0) fieldglob(:, myfirstlat_diag:mylastlat_diag, :) = field(:, :, :)
      ij = 1
      nr = 0
      do i = 0, MaxNodes - 1
        ini(i) = ij
        if (i.ne.nproc.and.nlatsinproc_d(i).gt.0) then
          call MPI_IRECV(bufrec(ij), nlatsinproc_d(i) * imax * levs, &
            MPI_DOUBLE_PRECISION, i, 93, comm, requestr(i), ierr)
          nr = nr + 1
          ij = ij + nlatsinproc_d(i) * imax * levs
        endif
      enddo
      do i = 1, nr
        call MPI_WAITANY(MaxNodes, requestr(0), index, status, ierr)
        ij = status(MPI_SOURCE)
        ip = ini(ij) - 1
        do k = 1, levs
          !CDIR NODEP
          do j = firstandlastlat(1, ij), firstandlastlat(2, ij)
            fieldglob(1:imax, j, k) = bufrec(ip + 1:ip + imax)
            ip = ip + imax
          enddo
        enddo
      enddo
    end if

  end subroutine Collect_Grid_d


  subroutine Collect_Grid_Full(field, fieldglob, levs, nproc)
    !# Collects Full Grid
    !# ---
    !# @info
    !# **Brief:** Collects Full Grid </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    
    integer, intent(in) :: nproc 
    !# destination processor
    integer, intent(in) :: levs
    real(kind = p_r8), intent(in) :: field(ibmax, levs, jbMax)
    real(kind = p_r8), intent(OUT) :: fieldglob(imax, jmax * levs)

    integer :: j, i, k, m, n, l, ic, iold, ks, j1, jc, js
    integer :: comm
    integer :: ierr
    integer :: index
    integer :: requests(nsends_diag + 1)
    integer :: requestr(0:maxnodes)
    integer :: status(MPI_STATUS_SIZE)
    integer :: stat(MPI_STATUS_SIZE, nsends_diag)
    integer :: ib(0:maxnodes)

    comm = mpiCommGroup
    if (dimrecbuf.lt.imax * jmax * levs) then
      dimrecbuf = imax * jmax * levs
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    if (myid.ne.nproc) then
      js = myfirstlat_diag - 1
      jc = (mylastlat_diag - js)
    else
      js = 0
      jc = jmax
    endif
    ib(0) = 1
    m = 0
    do k = 1, nrecs_diag
      ic = 0
      do i = m + 1, myrecspr_diag(2, k)
        ic = ic + myrecs_diag(2, i) - myrecs_diag(1, i) + 1
      enddo
      m = myrecspr_diag(2, k)
      ib(k) = ib(k - 1) + ic * levs
      call MPI_IRECV(bufrec(ib(k - 1)), ib(k) - ib(k - 1), MPI_DOUBLE_PRECISION, &
        myrecspr_diag(1, k), 88, comm, requestr(k), ierr)
    enddo
    m = 0
    ic = 0
    iold = 0
    do k = 1, nsends_diag
      do l = m + 1, mysendspr_diag(2, k)
        j = mysends_diag(3, l)
        do i = mysends_diag(1, l), mysends_diag(2, l)
          bufsend(ic + 1:ic + levs) = field(ibperij(i, j), :, jbperij(i, j))
          ic = ic + levs
        enddo
      enddo
      call MPI_ISEND(bufsend(iold + 1), ic - iold, MPI_DOUBLE_PRECISION, &
        mysendspr_diag(1, k), 88, comm, requests(k), ierr)
      m = mysendspr_diag(2, k)
      iold = ic
    enddo
    do j = max(myfirstlat, myfirstlat_diag), min(mylastlat, mylastlat_diag)
      j1 = j - js
      do k = 1, levs
        do i = myfirstlon(j), mylastlon(j)
          fieldglob(i, j1) = field(ibperij(i, j), k, jbperij(i, j))
        enddo
        j1 = j1 + jc
      enddo
    enddo
    do k = 1, nrecs_diag
      call MPI_WAITANY(nrecs_diag, requestr(1), index, status, ierr)
      ks = status(MPI_SOURCE)
      n = 1 
      !# avoiding uninitializated
      do l = 1, nrecs_diag
        if (ks.eq.myrecspr_diag(1, l)) then
          n = l
          ic = ib(n - 1) - 1
          m = myrecspr_diag(2, n - 1)
          EXIT
        endif
      enddo
      do l = m + 1, myrecspr_diag(2, n)
        j = myrecs_diag(3, l) - js
        do i = myrecs_diag(1, l), myrecs_diag(2, l)
          j1 = j
          do ks = 1, levs
            fieldglob(i, j1) = bufrec(ic + ks)
            j1 = j1 + jc
          end do
          ic = ic + levs
        end do
      end do
    end do
    if(nsends_diag.gt.0) call MPI_WAITALL(nsends_diag, requests(1), stat, ierr)
    if (myid.ne.nproc) then
      ic = (mylastlat_diag - myfirstlat_diag + 1) * imax * levs

      if (ic.gt.0) call MPI_ISEND(fieldglob, ic, MPI_DOUBLE_PRECISION, &
        nproc, 89, comm, requests(1), ierr)
      call MPI_WAIT(requests(1), status, ierr)
    else
      ib(0) = 1
      requestr = MPI_REQUEST_NULL
      n = 0
      do k = 0, maxnodes - 1
        if (k.ne.myid) then
          ic = (firstandlastlat(2, k) - firstandlastlat(1, k) + 1) * imax * levs
          if (ic.gt.0) then
            call MPI_IRECV(bufrec(ib(k)), ic, MPI_DOUBLE_PRECISION, &
              k, 89, comm, requestr(k), ierr)
            n = n + 1
          endif
        else
          ic = 0
        endif
        ib(k + 1) = ib(k) + ic
      enddo
      do k = 1, n
        call MPI_WAITANY(MaxNodes, requestr(0), index, status, ierr)
        ks = status(MPI_SOURCE)
        ic = ib(ks) - 1
        do l = 1, levs
          j1 = (l - 1) * jmax
          do j = firstandlastlat(1, ks), firstandlastlat(2, ks)
            Fieldglob(:, j1 + j) = bufrec(ic + 1:ic + imax)
            ic = ic + imax
          enddo
        enddo
      enddo
    endif

  end subroutine Collect_Grid_Full

  subroutine Collect_Gauss(gauss, gauss_out, nf)
    !# Collects Gaussian
    !# ---
    !# @info
    !# **Brief:** Collects Gaussian </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    
    integer, intent(in) :: nf
    real(kind = p_r8), intent(in) :: gauss(ibmax, nf, jbmax)
    real(kind = p_r8), intent(OUT) :: gauss_out(imax, myjmax_d, nf)
    real(kind = p_r8) :: saux(imax)

    integer :: j, i, k, m, n, l, ic, iold, ks, j1
    integer :: comm
    integer :: ierr
    integer :: index
    integer :: requests(nsends_diag + 1)
    integer :: requestr(0:maxnodes)
    integer :: status(MPI_STATUS_SIZE)
    integer :: stat(MPI_STATUS_SIZE, nsends_diag)
    integer :: ib(0:maxnodes)

    comm = mpiCommGroup
    if (dimrecbuf.lt.imax * jmax * nf) then
      dimrecbuf = imax * jmax * nf
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif

    if (dimsendbuf.lt.ibmax * jbmax * nf) then
      dimsendbuf = ibmax * jbmax * nf
      deallocate (bufsend)
      allocate (bufsend(dimsendbuf))
    endif

    ib(0) = 1
    m = 0
    do k = 1, nrecs_diag
      ic = 0
      do i = m + 1, myrecspr_diag(2, k)
        ic = ic + myrecs_diag(2, i) - myrecs_diag(1, i) + 1
      enddo
      m = myrecspr_diag(2, k)
      ib(k) = ib(k - 1) + ic * nf
      call MPI_IRECV(bufrec(ib(k - 1)), ib(k) - ib(k - 1), MPI_DOUBLE_PRECISION, &
        myrecspr_diag(1, k), 88, comm, requestr(k), ierr)
    enddo
    m = 0
    ic = 1
    iold = 1
    do k = 1, nsends_diag
      do n = 1, nf
        do l = m + 1, mysendspr_diag(2, k)
          j = mysends_diag(3, l)
          do i = mysends_diag(1, l), mysends_diag(2, l)
            bufsend(ic) = gauss(ibperij(i, j), n, jbperij(i, j))
            ic = ic + 1
          enddo
        enddo
      enddo
      call MPI_ISEND(bufsend(iold), ic - iold, MPI_DOUBLE_PRECISION, &
        mysendspr_diag(1, k), 88, comm, requests(k), ierr)
      m = mysendspr_diag(2, k)
      iold = ic
    enddo
    do k = 1, nf
      do j = max(myfirstlat, myfirstlat_diag), min(mylastlat, mylastlat_diag)
        j1 = j - myfirstlat_diag + 1
        do i = myfirstlon(j), mylastlon(j)
          gauss_out(i, j1, k) = gauss(ibperij(i, j), k, jbperij(i, j))
        enddo
      enddo
    enddo
    do k = 1, nrecs_diag
      call MPI_WAITANY(nrecs_diag, requestr(1), index, status, ierr)
      ks = status(MPI_SOURCE)
      do l = 1, nrecs_diag
        if (ks.eq.myrecspr_diag(1, l)) then
          n = l
          ic = ib(n - 1)
          m = myrecspr_diag(2, n - 1)
          EXIT
        endif
      enddo
      do ks = 1, nf
        do l = m + 1, myrecspr_diag(2, n)
          j = myrecs_diag(3, l) - myfirstlat_diag + 1
          do i = myrecs_diag(1, l), myrecs_diag(2, l)
            gauss_out(i, j, ks) = bufrec(ic)
            ic = ic + 1
          end do
        end do
      end do
    end do
    if(nsends_diag.gt.0) call MPI_WAITALL(nsends_diag, requests(1), stat, ierr)
    !   IF (reducedgrid) THEN
    !      DO k=1,nf
    !         DO j=myfirstlat_diag,mylastlat_diag
    !            j1 = j-myfirstlat_diag+1
    !            saux(1:imaxperj(j)) = gauss_out(1:imaxperj(j),j1,k)
    !            CALL CyclicLinear(iMaxPerJ(j), iMax, &
    !                               saux,gauss_out(1,j1,k),1,imax)
    !         ENDDO
    !      ENDDO
    !   ENDIF

  end subroutine Collect_Gauss


  subroutine Collect_Spec(field, fieldglob, levs, levsg, nproc)
    !# Collects Spectral Coefficients
    !# ---
    !# @info
    !# **Brief:** Collects Spectral Coefficients </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nproc 
    !# destination processor (has to be the first processor of one fourier group)
    integer, intent(in) :: levs
    integer, intent(in) :: levsg
    real(kind = p_r8), intent(in) :: field(2 * mymnmax, levs)
    real(kind = p_r8), intent(OUT) :: fieldglob(2 * mnmax_out, levsg)

    character(len = *), parameter :: h = "**(Collect_Spec)**"
    integer :: j, i, m, mn, mnloc, ns, l, lev, kdim, kp, kl, ll
    integer :: comm
    integer :: ierr
    integer :: index
    integer :: request
    integer :: requestr(0:MaxNodes)
    integer :: statu(MPI_STATUS_SIZE)
    integer :: status(MPI_STATUS_SIZE, maxnodes)
    !
    !   Collect inside fourier groups (to first processor in each group)
    !
    if (.not.any(first_proc_four.eq.nproc)) then
      write(p_nfprt, *) ' nproc ', nproc
      write(p_nfprt, "(a, ' Spectral fields should be collected to a first processor in a fourier group')") h
      stop h
    else if (levsg.eq.1.and..not.havesurf) then
      write(p_nfprt, *) ' myid  ', myid
      write(p_nfprt, "(a, ' should not be calling collect_spec of surface field')") h
      stop h
      !   ELSE IF (levsg.ne.1.and.levsg.ne.kmax) THEN
      !      WRITE(p_nfprt,*) ' levsg ',levsg
      !      WRITE(p_nfprt,"(a, ' collect_spec should be used for a global or a surface spectral field')") h
      !      STOP h
    end if
    comm = COMM_FOUR
    kdim = 2 * mnmaxlocal * levs
    if (dimrecbuf.lt.kdim * maxnodes_four) then
      dimrecbuf = kdim * maxnodes_four
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    if (myid_four.ne.0) then
      call MPI_ISEND(field, 2 * levs * mymnmax, MPI_DOUBLE_PRECISION, 0, 95, comm, request, ierr)
      call MPI_WAIT(request, status, ierr)
    else
      requestr(0) = MPI_REQUEST_NULL
      do i = 1, MaxNodes_four - 1
        call MPI_IRECV(bufrec(1 + i * kdim), 2 * mnsPerProc(i) * levs, &
          MPI_DOUBLE_PRECISION, i, 95, comm, requestr(i), ierr)
      enddo
      mnloc = 0
      mn = 0
      kl = myfirstlev - 1
      do m = 1, Mmax
        ns = 2 * (Mmax - m + 1)
        if(NodeHasM(m, mygroup_four).eq.0) then
          do l = 1, ns
            fieldglob(mn + l, kl + 1:kl + levs) = field(mnloc + l, 1:levs)
          enddo
          mnloc = mnloc + ns
        endif
        mn = mn + ns
      enddo
      do i = 1, MaxNodes_four - 1
        call MPI_WAITANY(MaxNodes_four, requestr(0), index, statu, ierr)
        j = statu(MPI_SOURCE)
        do lev = 1, levs
          mnloc = 2 * mnsPerProc(j) * (lev - 1)
          mn = 0
          do m = 1, Mmax
            ns = 2 * (Mmax - m + 1)
            if(NodeHasM(m, mygroup_four).eq.j) then
              do l = 1, ns
                fieldglob(mn + l, kl + lev) = bufrec(mnloc + l + j * kdim)
              enddo
              mnloc = mnloc + ns
            endif
            mn = mn + ns
          enddo
        enddo
      enddo

      if (levsg.eq.1.or.Ngroups_four.eq.1) return
      !   Collect Global Field
      !
      comm = mpiCommGroup
      if (myid.ne.nproc) then
        call MPI_ISEND(fieldglob(1, kl + 1), 2 * levs * mnmax_out, MPI_DOUBLE_PRECISION, &
          nproc, 96, comm, request, ierr)
        call MPI_WAIT(request, status, ierr)
      else
        requestr(1:Ngroups_four) = MPI_REQUEST_NULL
        do i = 1, Ngroups_four
          kp = first_proc_four(i)
          kl = kfirst_four(kp)
          ll = nlevperg_four(i)
          if (kp.ne.nproc) then
            call MPI_IRECV(fieldglob(1, kl), 2 * ll * mnmax_out, &
              MPI_DOUBLE_PRECISION, kp, 96, comm, requestr(i), ierr)
          endif
        enddo
        call MPI_WAITALL(Ngroups_four, requestr(1), status, ierr)
      endif
    endif

  end subroutine Collect_Spec


  subroutine Spread_surf_Spec(field)
    !# Spreads surf Spec
    !# ---
    !# @info
    !# **Brief:** Spreads surf Spec </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    real(kind = p_r8), intent(inout) :: field(2 * mymnmax)

    character(len = *), parameter :: h = "**(Spread_surf_Spec)**"
    integer :: i, m, mng, len, np, n
    integer :: ini(maxnodes)
    integer :: comm
    integer :: ierr
    integer :: requestr(0:MaxNodes)
    integer :: status(MPI_STATUS_SIZE, maxnodes)

    !
    if (ngroups_four.eq.1) return
    comm = mpiCommGroup
    requestr = MPI_REQUEST_NULL
    ini(1) = 0
    do n = 1, ncomm_spread
      ini(n + 1) = ini(n) + comm_spread(n, 2)
    enddo
    if (mygroup_four.eq.1) then
      do n = 2, ngroups_four
        mng = 0
        do m = 1, mymmax
          len = 2 * (MMax + 1 - lm2m(m))
          np = ms_spread(m, n)
          bufsend(ini(np) + 1:ini(np) + len) = field(mng + 1:mng + len)
          ini(np) = ini(np) + len
          mng = mng + len
        enddo
      enddo
      ini(1) = 1
      do n = 1, ncomm_spread
        ini(n + 1) = ini(n) + comm_spread(n, 2)
        call MPI_ISEND(bufsend(ini(n)), comm_spread(n, 2), MPI_DOUBLE_PRECISION, comm_spread(n, 1), 75, comm, requestr(n), ierr)
      enddo
      call MPI_WAITALL(ncomm_spread, requestr(1), status, ierr)
    else
      do i = 1, ncomm_spread
        call MPI_IRECV(bufrec(ini(i) + 1), comm_spread(i, 2), &
          MPI_DOUBLE_PRECISION, comm_spread(i, 1), 75, comm, requestr(i), ierr)
      enddo
      call MPI_WAITALL(ncomm_spread, requestr(1), status, ierr)
      mng = 0
      do m = 1, mymmax
        len = 2 * (MMax + 1 - lm2m(m))
        np = ms_spread(m, 1)
        field(mng + 1:mng + len) = bufrec(ini(np) + 1:ini(np) + len)
        ini(np) = ini(np) + len
        mng = mng + len
      enddo
    endif

  end subroutine Spread_surf_Spec


  subroutine Collect_Spec_Ext(field, fieldglob, levs, levsg, nproc)
    !# Collects Spectral Coefficients Extended
    !# ---
    !# @info
    !# **Brief:** Collects Spectral Coefficients Extended </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nproc 
    ! destination processor
    integer, intent(in) :: levs
    integer, intent(in) :: levsg
    real(kind = p_r8), intent(in) :: field(2 * mymnextmax, levs)
    real(kind = p_r8), intent(OUT) :: fieldglob(2 * mnextmax, levsg)

    character(len = *), parameter :: h = "**(Collect_Spec_Ext)**"
    integer :: j, i, m, mn, mnloc, ns, l, lev, kdim, kp, kl, ll
    integer :: comm
    integer :: ierr
    integer :: index
    integer :: request
    integer :: requestr(0:MaxNodes)
    integer :: statu(MPI_STATUS_SIZE)
    integer :: status(MPI_STATUS_SIZE, maxnodes)
    ! 
    !   Collect inside fourier groups (to first processor in each group)
    ! 
    if (.not.any(first_proc_four.eq.nproc)) then
      write(p_nfprt, *) ' nproc ', nproc
      write(p_nfprt, "(a, ' Spectral fields should be collected to a first processor in a fourier group')") h
      stop h
    else if (levs.eq.1.and..not.havesurf) then
      write(p_nfprt, *) ' myid  ', myid
      write(p_nfprt, "(a, ' should not be calling collect_spec_ext of surface field')") h
      stop h
    else if (levsg.ne.1.and.levsg.ne.kmax) then
      write(p_nfprt, *) ' levsg ', levsg
      write(p_nfprt, "(a, ' collect_spec_ext should be used for a global or a surface spectral field')") h
      stop h
    end if
    comm = COMM_FOUR
    kdim = 2 * mnextmaxlocal * levs
    if (dimrecbuf.lt.kdim * maxnodes_four) then
      dimrecbuf = kdim * maxnodes_four
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    if (myid_four.ne.0) then
      call MPI_ISEND(field, 2 * levs * mymnextmax, MPI_DOUBLE_PRECISION, &
        0, 95, comm, request, ierr)
      call MPI_WAIT(request, status, ierr)
    else
      requestr(0) = MPI_REQUEST_NULL
      do i = 1, MaxNodes_four - 1
        call MPI_IRECV(bufrec(1 + i * kdim), 2 * mnsExtPerProc(i) * levs, &
          MPI_DOUBLE_PRECISION, i, 95, comm, requestr(i), ierr)
      enddo
      mnloc = 0
      mn = 0
      kl = myfirstlev - 1
      do m = 1, Mmax
        ns = 2 * (Mmax - m + 2)
        if(NodeHasM(m, mygroup_four).eq.0) then
          do l = 1, ns
            fieldglob(mn + l, kl + 1:kl + levs) = field(mnloc + l, :)
          enddo
          mnloc = mnloc + ns
        endif
        mn = mn + ns
      enddo
      do i = 1, MaxNodes_four - 1
        call MPI_WAITANY(MaxNodes_four, requestr(0), index, statu, ierr)
        j = statu(MPI_SOURCE)
        do lev = 1, levs
          mnloc = 2 * mnsExtPerProc(j) * (lev - 1)
          mn = 0
          do m = 1, Mmax
            ns = 2 * (Mmax - m + 2)
            if(NodeHasM(m, mygroup_four).eq.j) then
              do l = 1, ns
                fieldglob(mn + l, kl + lev) = bufrec(mnloc + l + j * kdim)
              enddo
              mnloc = mnloc + ns
            endif
            mn = mn + ns
          enddo
        enddo
      enddo

      if (levsg.eq.1.or.Ngroups_four.eq.1) return
      !   Collect Global Field
      ! 
      comm = mpiCommGroup
      if (myid.ne.nproc) then
        call MPI_ISEND(fieldglob(1, kl + 1), 2 * levs * mnextmax, MPI_DOUBLE_PRECISION, &
          nproc, 96, comm, request, ierr)
        call MPI_WAIT(request, status, ierr)
      else
        requestr(1:Ngroups_four) = MPI_REQUEST_NULL
        do i = 1, Ngroups_four
          kp = first_proc_four(i)
          kl = kfirst_four(kp)
          ll = nlevperg_four(i)
          if (kp.ne.nproc) then
            call MPI_IRECV(fieldglob(1, kl), 2 * ll * mnextmax, &
              MPI_DOUBLE_PRECISION, kp, 96, comm, requestr(i), ierr)
          endif
        enddo
        call MPI_WAITALL(Ngroups_four, requestr(1), status, ierr)
      endif
    endif

  end subroutine Collect_Spec_Ext

  subroutine Exchange_ftog(nrecs_f, nrecs_g)
    !# Exchanges ftog
    !# ---
    !# @info
    !# **Brief:** Exchanges ftog </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nrecs_f
    integer, intent(OUT) :: nrecs_g

    integer :: i, m, k
    integer :: comm
    integer :: ierr
    integer :: ns(0:MaxNodes - 1)
    integer :: requestr(0:MaxNodes - 1)
    integer :: requests(0:MaxNodes - 1)
    integer :: status(MPI_STATUS_SIZE, maxnodes)

    comm = mpiCommGroup
    requestr(myid) = MPI_REQUEST_NULL
    requests(myid) = MPI_REQUEST_NULL
    messproc_g(2, 0) = 0
    do i = 0, MaxNodes - 1
      if (i.ne.myid) then
        call MPI_IRECV(messproc_g(2, i + 1), 1, MPI_INTEGER, i, 18, &
          comm, requestr(i), ierr)
      else
        messproc_g(2, i + 1) = 0
      endif
    enddo
    m = 0
    k = 1
    do i = 0, MaxNodes - 1
      if (i.ne.myid) then
        if (k.le.nrecs_f.and.i.eq.messproc_f(1, k)) then
          ns(i) = messproc_f(2, k) - m
          m = messproc_f(2, k)
          k = k + 1
        else
          ns(i) = 0
        endif
        call MPI_ISEND(ns(i), 1, MPI_INTEGER, i, 18, comm, requests(i), ierr)
      endif
    enddo
    call MPI_WAITALL(MaxNodes, requestr(0), status, ierr)
    call MPI_WAITALL(MaxNodes, requests(0), status, ierr)
    k = 0
    m = 0
    do i = 0, MaxNodes - 1
      if (messproc_g(2, i + 1).ne.0) then
        k = k + 1
        m = m + messproc_g(2, i + 1)
        messproc_g(2, k) = m
        messproc_g(1, k) = i
      endif
    enddo
    nrecs_g = k
    m = 0
    do i = 1, nrecs_g
      ns(i) = messproc_g(2, i) - m
      call MPI_IRECV(messages_g(1, m + 1), 4 * ns(i), MPI_INTEGER, messproc_g(1, i), 19, &
        comm, requestr(i), ierr)
      m = messproc_g(2, i)
    enddo
    m = 0
    do i = 1, nrecs_f
      ns(i) = messproc_f(2, i) - m
      call MPI_ISEND(messages_f(1, m + 1), ns(i) * 4, MPI_INTEGER, messproc_f(1, i), 19, &
        comm, requests(i), ierr)
      m = messproc_f(2, i)
    enddo
    if (nrecs_g.gt.0) call MPI_WAITALL(nrecs_g, requestr(1), status, ierr)
    if (nrecs_f.gt.0) call MPI_WAITALL(nrecs_f, requests(1), status, ierr)

  end subroutine Exchange_ftog

  subroutine Exchange_diag(nrecs_diag, nsends_diag)
    !# Exchanges diag
    !# ---
    !# @info
    !# **Brief:** Exchanges diag </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nrecs_diag
    integer, intent(OUT) :: nsends_diag

    integer :: i, m, k
    integer :: comm
    integer :: ierr
    integer :: ns(0:MaxNodes)
    integer :: requestr(0:MaxNodes - 1)
    integer :: requests(0:MaxNodes - 1)
    integer :: status(MPI_STATUS_SIZE, maxnodes)

    comm = mpiCommGroup
    requestr(myid) = MPI_REQUEST_NULL
    requests(myid) = MPI_REQUEST_NULL
    mysendspr_diag(2, 0) = 0
    do i = 0, MaxNodes - 1
      if (i.ne.myid) then
        call MPI_IRECV(mysendspr_diag(2, i + 1), 1, MPI_INTEGER, i, 18, &
          comm, requestr(i), ierr)
      else
        mysendspr_diag(2, i + 1) = 0
      endif
    enddo
    m = 0
    k = 1
    do i = 0, MaxNodes - 1
      if (i.ne.myid) then
        if (k.le.nrecs_diag.and.i.eq.myrecspr_diag(1, k)) then
          ns(i) = myrecspr_diag(2, k) - m
          m = myrecspr_diag(2, k)
          k = k + 1
        else
          ns(i) = 0
        endif
        call MPI_ISEND(ns(i), 1, MPI_INTEGER, i, 18, comm, requests(i), ierr)
      endif
    enddo
    call MPI_WAITALL(MaxNodes, requestr(0), status, ierr)
    call MPI_WAITALL(MaxNodes, requests(0), status, ierr)
    k = 0
    m = 0
    do i = 0, MaxNodes - 1
      if (mysendspr_diag(2, i + 1).ne.0) then
        k = k + 1
        m = m + mysendspr_diag(2, i + 1)
        mysendspr_diag(2, k) = m
        mysendspr_diag(1, k) = i
      endif
    enddo
    nsends_diag = k
    m = 0
    do i = 1, nsends_diag
      ns(i) = mysendspr_diag(2, i) - m
      call MPI_IRECV(mysends_diag(1, m + 1), 4 * ns(i), MPI_INTEGER, mysendspr_diag(1, i), 19, &
        comm, requestr(i), ierr)
      m = mysendspr_diag(2, i)
    enddo
    m = 0
    do i = 1, nrecs_diag
      ns(i) = myrecspr_diag(2, i) - m
      call MPI_ISEND(myrecs_diag(1, m + 1), ns(i) * 4, MPI_INTEGER, myrecspr_diag(1, i), 19, &
        comm, requests(i), ierr)
      m = myrecspr_diag(2, i)
    enddo
    if (nsends_diag.gt.0) call MPI_WAITALL(nsends_diag, requestr(1), status, ierr)
    if (nrecs_diag.gt.0) call MPI_WAITALL(nrecs_diag, requests(1), status, ierr)

  end subroutine Exchange_diag

  subroutine Exchange_Hallos(nrec, nsend, nscalars)
    !# Exchanges Hallos
    !# ---
    !# @info
    !# **Brief:** Exchanges Hallos </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nrec
    integer, intent(in) :: nscalars
    integer, intent(OUT) :: nsend

    integer :: i, m, k, irec, nlen, isnd
    integer :: comm
    integer :: ierr
    integer :: ns(0:MaxNodes - 1)
    integer :: requestr(0:MaxNodes - 1)
    integer :: requests(0:MaxNodes - 1)
    integer :: status(MPI_STATUS_SIZE, maxnodes)

    comm = mpiCommGroup
    requestr(myid) = MPI_REQUEST_NULL
    requests(myid) = MPI_REQUEST_NULL
    do i = 0, MaxNodes - 1
      if (i.ne.myid) then
        call MPI_IRECV(mysendspr(2, i + 1), 1, MPI_INTEGER, i, 15, &
          comm, requestr(i), ierr)
      else
        mysendspr(2, i + 1) = 0
      endif
    enddo
    m = 0
    k = 1
    do i = 0, MaxNodes - 1
      if (i.ne.myid) then
        if (k.le.nrec.and.i.eq.myrecspr(1, k)) then
          ns(i) = myrecspr(2, k) - m
          m = myrecspr(2, k)
          k = k + 1
        else
          ns(i) = 0
        endif
        call MPI_ISEND(ns(i), 1, MPI_INTEGER, i, 15, comm, requests(i), ierr)
      endif
    enddo
    call MPI_WAITALL(MaxNodes, requestr(0), status, ierr)
    call MPI_WAITALL(MaxNodes, requests(0), status, ierr)
    k = 0
    m = 0
    do i = 0, MaxNodes - 1
      if (mysendspr(2, i + 1).ne.0) then
        k = k + 1
        m = m + mysendspr(2, i + 1)
        mysendspr(2, k) = m
        mysendspr(1, k) = i
      endif
    enddo
    nsend = k
    m = 0
    do i = 1, nsend
      ns(i) = mysendspr(2, i) - m
      call MPI_IRECV(mysendsgr(1, m + 1), 4 * ns(i), MPI_INTEGER, mysendspr(1, i), 16, &
        comm, requestr(i), ierr)
      m = mysendspr(2, i)
    enddo
    allocate (isbrec (nrec))
    allocate (ilrecbuf(nrec))
    allocate (isbsend(nsend))
    allocate (ilsendbuf(nsend))
    m = 0
    irec = 1
    do i = 1, nrec
      isbrec(i) = irec
      ns(i) = myrecspr(2, i) - m
      call MPI_ISEND(myrecsgr(1, m + 1), ns(i) * 4, MPI_INTEGER, myrecspr(1, i), 16, &
        comm, requests(i), ierr)
      nlen = 0
      do k = m + 1, myrecspr(2, i)
        nlen = nlen + myrecsgr(2, k) - myrecsgr(1, k) + 1
      enddo
      ilrecbuf(i) = nlen * (1 + kmax * (nscalars + 4))
      irec = ilrecbuf(i) + irec
      m = myrecspr(2, i)
    enddo
    call MPI_WAITALL(nsend, requestr(1), status, ierr)
    call MPI_WAITALL(nrec, requests(1), status, ierr)

    m = 0
    isnd = 1
    do i = 1, nsend
      isbsend(i) = isnd
      nlen = 0
      do k = m + 1, mysendspr(2, i)
        nlen = nlen + mysendsgr(2, k) - mysendsgr(1, k) + 1
      enddo
      ilsendbuf(i) = nlen * (1 + kmax * (nscalars + 4))
      isnd = ilsendbuf(i) + isnd
      m = mysendspr(2, i)
    enddo

  end subroutine Exchange_Hallos

  subroutine Exchange_Fields (u, v, t, q, lps, fgpass_scalar, adr, nscalars, &
    nrec, nsend)
    !# Exchanges Fields
    !# ---
    !# @info
    !# **Brief:** Exchanges Fields </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nsend
    integer, intent(in) :: nrec
    integer, intent(in) :: nscalars
    integer, intent(in) :: adr
    real(kind = p_r8), intent(inout) :: u(ibMax, kmax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: v(ibMax, kmax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: t(ibMax, kmax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: q(ibMax, kmax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: &
      fgpass_scalar(ibMax, kmax, jbMax_ext, nscalars, 2)
    real(kind = p_r8), intent(inout) :: lps(ibMax, jbMax_ext)
    integer :: index
    integer :: statu(MPI_STATUS_SIZE)
    integer :: status(MPI_STATUS_SIZE, nsend)
    integer :: requests(nsend)
    integer :: requestr(nrec)

    integer :: j, ns, ibr, ibs, jbr, k
    integer :: m, i, kr, ks, n
    integer :: comm, ierr

    if (dimrecbuf.lt.isbrec(nrec) + ilrecbuf(nrec)) then
      dimrecbuf = isbrec(nrec) + ilrecbuf(nrec)
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    if (dimsendbuf.lt.isbsend(nsend) + ilsendbuf(nsend)) then
      dimsendbuf = isbsend(nsend) + ilsendbuf(nsend)
      deallocate (bufsend)
      allocate (bufsend(dimsendbuf))
    endif
    comm = mpiCommGroup
    do k = 1, nrec
      call MPI_IRECV(bufrec(isbrec(k)), ilrecbuf(k), MPI_DOUBLE_PRECISION, &
        myrecspr(1, k), 75, comm, requestr(k), ierr)
    enddo
    m = 1
    do k = 1, nsend
      ibs = isbsend(k) - 1
      do ns = m, mysendspr(2, k)
        j = mysendsgr(3, ns)
        do i = mysendsgr(1, ns), mysendsgr(2, ns)
          jbr = jbperij(i, j)
          ibr = ibperij(i, j)
          bufsend(ibs + 1:ibs + kmax) = u(ibr, :, jbr)
          ibs = ibs + kmax
          bufsend(ibs + 1:ibs + kmax) = v(ibr, :, jbr)
          ibs = ibs + kmax
          bufsend(ibs + 1:ibs + kmax) = t(ibr, :, jbr)
          ibs = ibs + kmax
          bufsend(ibs + 1:ibs + kmax) = q(ibr, :, jbr)
          ibs = ibs + kmax + 1
          bufsend(ibs) = lps(ibr, jbr)
          do n = 1, nscalars
            bufsend(ibs + 1:ibs + kmax) = fgpass_scalar(ibr, :, jbr, n, adr)
            ibs = ibs + kmax
          enddo
        enddo
      enddo
      call MPI_ISEND(bufsend(isbsend(k)), ilsendbuf(k), MPI_DOUBLE_PRECISION, &
        mysendspr(1, k), 75, comm, requests(k), ierr)
      m = mysendspr(2, k) + 1
    enddo
    do k = 1, nrec
      call MPI_WAITANY(nrec, requestr, index, statu, ierr)
      kr = statu(MPI_SOURCE)
      do j = 1, nrec
        if(myrecspr(1, j).eq.kr) then
          ks = j
          exit
        endif
      enddo
      ibs = isbrec(ks) - 1
      if(ks.eq.1) then
        m = 1
      else
        m = myrecspr(2, ks - 1) + 1
      endif
      do ns = m, myrecspr(2, ks)
        j = myrecsgr(3, ns)
        do i = myrecsgr(1, ns), myrecsgr(2, ns)
          jbr = jbperij(i, j)
          ibr = ibperij(i, j)
          u(ibr, :, jbr) = bufrec(ibs + 1:ibs + kmax)
          ibs = ibs + kmax
          v(ibr, :, jbr) = bufrec(ibs + 1:ibs + kmax)
          ibs = ibs + kmax
          t(ibr, :, jbr) = bufrec(ibs + 1:ibs + kmax)
          ibs = ibs + kmax
          q(ibr, :, jbr) = bufrec(ibs + 1:ibs + kmax)
          ibs = ibs + kmax + 1
          lps(ibr, jbr) = bufrec(ibs)
          do n = 1, nscalars
            fgpass_scalar(ibr, :, jbr, n, adr) = bufrec(ibs + 1:ibs + kmax)
            ibs = ibs + kmax
          enddo
        enddo
      enddo
    enddo
    call MPI_WAITALL(nsend, requests, status, ierr)

  end subroutine Exchange_Fields


  subroutine Exchange_Winds (u, v, w, um, vm, nrec, nsend)
    !# Exchanges Winds
    !# ---
    !# @info
    !# **Brief:** Exchanges Winds </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    integer, intent(in) :: nsend
    integer, intent(in) :: nrec
    real(kind = p_r8), intent(inout) :: u(ibMax, kmax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: v(ibMax, kmax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: w(ibMax, kmax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: um(ibMax, jbMax_ext)
    real(kind = p_r8), intent(inout) :: vm(ibMax, jbMax_ext)
    integer :: index
    integer :: statu(MPI_STATUS_SIZE)
    integer :: status(MPI_STATUS_SIZE, nsend)
    integer :: requests(nsend)
    integer :: requestr(nrec)

    integer :: j, ns, ibr, ibs, jbr, k
    integer :: m, i, kr, ks
    integer :: comm, ierr

    comm = mpiCommGroup
    if (dimrecbuf.lt.isbrec(nrec) + ilrecbuf(nrec)) then
      dimrecbuf = isbrec(nrec) + ilrecbuf(nrec)
      deallocate (bufrec)
      allocate (bufrec(dimrecbuf))
    endif
    if (dimsendbuf.lt.isbsend(nsend) + ilsendbuf(nsend)) then
      dimsendbuf = isbsend(nsend) + ilsendbuf(nsend)
      deallocate (bufsend)
      allocate (bufsend(dimsendbuf))
    endif

    do k = 1, nrec
      call MPI_IRECV(bufrec(isbrec(k)), ilrecbuf(k), MPI_DOUBLE_PRECISION, &
        myrecspr(1, k), 76, comm, requestr(k), ierr)
    enddo
    m = 1
    do k = 1, nsend
      ibs = isbsend(k) - 1
      do ns = m, mysendspr(2, k)
        j = mysendsgr(3, ns)
        do i = mysendsgr(1, ns), mysendsgr(2, ns)
          jbr = jbperij(i, j)
          ibr = ibperij(i, j)
          bufsend(ibs + 1:ibs + kmax) = u(ibr, :, jbr)
          ibs = ibs + kmax
          bufsend(ibs + 1:ibs + kmax) = v(ibr, :, jbr)
          ibs = ibs + kmax
          bufsend(ibs + 1:ibs + kmax) = w(ibr, :, jbr)
          ibs = ibs + kmax + 1
          bufsend(ibs) = um(ibr, jbr)
          ibs = ibs + 1
          bufsend(ibs) = vm(ibr, jbr)
        enddo
      enddo

      call MPI_ISEND(bufsend(isbsend(k)), ibs - isbsend(k) + 1, MPI_DOUBLE_PRECISION, &
        mysendspr(1, k), 76, comm, requests(k), ierr)
      m = mysendspr(2, k) + 1
    enddo
    do k = 1, nrec
      call MPI_WAITANY(nrec, requestr, index, statu, ierr)
      kr = statu(MPI_SOURCE)

      do j = 1, nrec
        if(myrecspr(1, j).eq.kr) then
          ks = j
          exit
        endif
      enddo
      ibs = isbrec(ks) - 1
      if(ks.eq.1) then
        m = 1
      else
        m = myrecspr(2, ks - 1) + 1
      endif
      do ns = m, myrecspr(2, ks)
        j = myrecsgr(3, ns)
        do i = myrecsgr(1, ns), myrecsgr(2, ns)
          jbr = jbperij(i, j)
          ibr = ibperij(i, j)
          u(ibr, :, jbr) = bufrec(ibs + 1:ibs + kmax)
          ibs = ibs + kmax
          v(ibr, :, jbr) = bufrec(ibs + 1:ibs + kmax)
          ibs = ibs + kmax
          w(ibr, :, jbr) = bufrec(ibs + 1:ibs + kmax)
          ibs = ibs + kmax + 1
          um(ibr, jbr) = bufrec(ibs)
          ibs = ibs + 1
          vm(ibr, jbr) = bufrec(ibs)
        enddo
      enddo
    enddo
    call MPI_WAITALL(nsend, requests, status, ierr)

  end subroutine Exchange_Winds

  subroutine Clear_Communications()
    !# Cleans Communications
    !# ---
    !# @info
    !# **Brief:** Cleans Communications </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 </br>
    !# @endinfo
    deallocate (bufrec)
    deallocate (bufsend)
    !  DEALLOCATE (isbrec)
    !  DEALLOCATE (isbsend)
    !  DEALLOCATE (ilrecbuf)
    !  DEALLOCATE (ilsendbuf)
  end subroutine Clear_Communications

end module Mod_Communications

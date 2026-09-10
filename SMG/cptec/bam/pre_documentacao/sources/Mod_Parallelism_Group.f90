!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Parallelism_Group </br></br>
!#
!# **Brief**: Deals with parallelism for Modules requires more than one processor,
!# but not all. </br></br>
!# 
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!#
!# **Author**: Denis Eiras </br>
!#
!# **Version**: 1.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>03-09-2019 - Denis Eiras - version 1.0.0 - Create </li>
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

module Mod_Parallelism_Group

  use Mod_Parallelism, only : &
    myId &         ! MPI process rank
    , maxNodesWorld => maxNodes &
    , getUnitDump

  use Mod_String_Functions

  use Mod_Messages, only : &
  msgOut &
  , msgDumpModMessages => msgDump

  implicit none

  private

  type ParallelGroup
    character(len = 20) :: groupName

    integer, public :: mpiCommGroup
    integer, public :: maxNodes      
    !# MPI processes in the group
    integer, public :: mpiMasterProc 
    !# Master processor
    integer, allocatable :: groupRanks(:)
  end type

  type(ParallelGroup), allocatable :: groups(:)

  character(len = *), parameter :: headerMsg = 'Parallel Group        | '

  ! public procedures
  public :: createMpiCommGroup
  public :: getMaxNodes
  public :: getMpiCommGroup
  public :: getMpiMasterProc
  public :: getIsMasterProc           
  !#  MPI process rank of master in this group

  public :: fatalErrorGroup
  public :: msgOutMasterGroup
  public :: msgDump

  include 'mpif.h'
  include 'precision.h'
  include 'messages.h'


contains

  function getMaxNodes(name_group) result(maxRanks)
    !# Gets Max Nodes
    !# ---
    !# @info
    !# **Brief:** Gets Max Nodes. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: name_group
    integer :: maxRanks

    maxRanks = groups(getGroupIndex(trim(name_group)))%maxNodes

  end function getMaxNodes


  function getMpiMasterProc(name_group) result(mpiMasterProc)
    !# Gets Mpi Master Processor
    !# ---
    !# @info
    !# **Brief:** Gets Mpi Master Processor. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: name_group
    integer :: mpiMasterProc

    mpiMasterProc = groups(getGroupIndex(trim(name_group)))%mpiMasterProc

  end function getMpiMasterProc


  function getMpiCommGroup(name_group) result(mpiCommGroup)
    !# Gets Mpi Comm Group
    !# ---
    !# @info
    !# **Brief:** Gets Mpi Comm Group. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: name_group
    integer :: mpiCommGroup

    mpiCommGroup = groups(getGroupIndex(trim(name_group)))%mpiCommGroup

  end function getMpiCommGroup


  function getIsMasterProc(nameGroup) result(isMaster)
    !# Returns if the Master Processor is running
    !# ---
    !# @info
    !# **Brief:** Returns if the Master Processor is running. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    character(len = *) :: nameGroup

    logical :: isMaster
    isMaster = (myId .eq. getMpiMasterProc(nameGroup))
  end function getIsMasterProc


  function getGroupIndex(nameGroup) result(idxGroup)
    !# Gets Group Index
    !# ---
    !# @info
    !# **Brief:** Gets Group Index. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    character(len = *) :: nameGroup
    integer :: idxGroup

    do idxGroup = 1, size(groups)
      if(trim(nameGroup) == trim(groups(idxGroup)%groupName)) then
        return
      endif
    enddo
    idxGroup = -1
  end function getGroupIndex


  subroutine createMpiCommGroup(name_group, mpiMasterProc_group, maxNodes_group)
    !# Creates Mpi Communicator Group
    !# ---
    !# @info
    !# **Brief:** Creates Mpi Communicator Group. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    character(len = *) :: name_group
    integer, intent(in) :: mpiMasterProc_group
    integer, intent(in) :: maxNodes_group

    integer :: worldGroup, thisGroup, ierr, groupSize, groupRank, i, newGrpIdx
    type(ParallelGroup), allocatable :: groups_last(:)

    if (allocated(groups)) then
      newGrpIdx = size(groups) + 1
      allocate(groups_last(newGrpIdx - 1))
      groups_last = groups
      deallocate(groups)
      allocate(groups(newGrpIdx))
      groups(1:newGrpIdx - 1) = groups_last
    else
      newGrpIdx = 1
      allocate(groups(newGrpIdx))
    endif

    groups(newGrpIdx)%groupName = trim(name_group)
    groups(newGrpIdx)%mpiMasterProc = mpiMasterProc_group
    groups(newGrpIdx)%maxNodes = maxNodes_group
    allocate(groups(newGrpIdx)%groupRanks(maxNodes_group))

    do i = 1, groups(newGrpIdx)%maxNodes
      groups(newGrpIdx)%groupRanks(i) = i - 1 + groups(newGrpIdx)%mpiMasterProc
    enddo

    call MPI_Comm_group(MPI_COMM_WORLD, worldGroup, ierr)
    call MPI_Group_incl(worldGroup, groups(newGrpIdx)%maxNodes, groups(newGrpIdx)%groupRanks, thisGroup, ierr)
    !    call MPI_Comm_create_group(MPI_COMM_WORLD, thisGroup, 99999, groups(newGrpIdx)%mpiCommGroup, ierr)

    call MPI_Comm_create(MPI_COMM_WORLD, thisGroup, groups(newGrpIdx)%mpiCommGroup, ierr)

    !    thisGroup = myId / groups(newGrpIdx)%maxNodes
    !    call MPI_COMM_SPLIT(MPI_COMM_WORLD, thisGroup, myId, groups(newGrpIdx)%mpiCommGroup, ierr)

    groupRank = -1
    groupSize = -1
    ! If this rank isn't in the new communicator, it will be MPI_COMM_NULL.
    ! Using MPI_COMM_NULL for MPI_Comm_rank or MPI_Comm_size is erroneous
    if (MPI_COMM_NULL .ne. groups(newGrpIdx)%mpiCommGroup) then
      call MPI_Comm_rank(groups(newGrpIdx)%mpiCommGroup, groupRank, ierr)
      call MPI_Comm_size(groups(newGrpIdx)%mpiCommGroup, groupSize, ierr)
      call msgOutMasterGroup(headerMsg, "Group " // trim(groups(newGrpIdx)%groupName) // &
        " created with " // intToStr(groups(newGrpIdx)%maxNodes) // " elements", trim(name_group))
      call MPI_BARRIER(groups(newGrpIdx)%mpiCommGroup, ierr)
    else
      call msgOut(headerMsg, "Error creating Group " // trim(groups(newGrpIdx)%groupName))
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  end subroutine createMpiCommGroup


  subroutine msgOutMasterGroup(header, message, nameGroup)
    !# Sends messages only in master processor of Group
    !# ---
    !# @info
    !# **Brief:** Sends messages only in master processor of Group. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: header
    character(len = *), intent(in) :: message
    character(len = *) :: nameGroup

    if (getIsMasterProc(nameGroup)) then
      call msgOut(header, message)
    end if
  end subroutine msgOutMasterGroup


  subroutine msgDump(header, message)
    !# Sends messages only in master processor of Group
    !# ---
    !# @info
    !# **Brief:** Sends messages only in master processor of Group. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: header
    character(len = *) :: message

    call msgDumpModMessages(getUnitDump(), header, message)

  end subroutine msgDump


  subroutine fatalErrorGroup(header, message, mpiCommGroupName, mpiCommGroup)
    !# Dumps error message everywhere in Group and destroy group parallelism
    !# ---
    !# @info
    !# **Brief:** Dumps error message everywhere in Group and destroy group parallelism. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: header
    character(len = *), intent(in) :: message
    character(len = *), intent(in) :: mpiCommGroupName
    integer, intent(in) :: mpiCommGroup
    integer :: ierror = -1
    integer :: ierr
    character(len = *), parameter :: headerFatal = "***(FATAL ERROR IN MPI GROUP)***"

    call msgOut(header, headerFatal // " MPI GROUP = " // mpiCommGroupName // ' ' // message)
    call msgDump(header, headerFatal // " MPI GROUP = " // mpiCommGroupName // ' ' // message)
    ! only abort this mpi communicator !
    ! call MPI_ABORT(mpiCommGroup, ierror, ierr)
    ! call MPI_Comm_free(mpiCommGroup)
    
  end subroutine fatalErrorGroup

end module Mod_Parallelism_Group


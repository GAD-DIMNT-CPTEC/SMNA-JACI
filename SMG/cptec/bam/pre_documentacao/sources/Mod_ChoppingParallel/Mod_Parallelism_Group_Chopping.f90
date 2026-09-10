!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Parallelism_Group_Chopping </br></br>
!#
!# **Brief**: Deals with parallelism for Modules requires more than one processor,
!# but not all. Overrides Mod_Parallelism maxNodes and mpiMasterProc. </br></br>
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
!#  <li>03-06-2019 - Denis Eiras = version 1.0.0 - Create </li>
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

module Mod_Parallelism_Group_Chopping

  use Mod_Parallelism, only : &
    myId &         
    ! MPI process rank
    , getMyIdString &
    , parf_bcast_int_scalar &
    , parf_bcast_int_1d &
    , parf_bcast_char &
    , parf_bcast_char_vec &
    , parf_send_char &
    , parf_get_char &
    , parf_barrier &
    , parf_get_noblock_char &
    , parf_send_noblock_char &
    , parf_wait_any_nostatus &
    , getUnitDump

  use Mod_Parallelism_Group, only : &
    msgOutMasterGroup &
    , msgDumpGroup => msgDump &
    , getMaxNodes &
    , getMpiCommGroup &
    , getMpiMasterProc &
    , getIsMasterProc &   
    !  MPI process rank of master in this group
    , fatalErrorGroup

  implicit none

  private

  character(len = *), parameter :: nameGroup = "Chopping"

  ! public variables
  integer, public :: mpiCommGroup
  integer, public :: maxNodes
  integer, public :: mpiMasterProc
  public :: myId

  ! public procedures
  public :: initializeParallelismVariables
  public :: msgOutMaster
  public :: isMasterProc           
  !  MPI process rank of master in this group
  public :: getMyIdString
  public :: fatalError
  public :: msgDump
  public :: parf_bcast_int_scalar
  public :: parf_bcast_int_1d
  public :: parf_bcast_char
  public :: parf_bcast_char_vec
  public :: parf_send_char
  public :: parf_get_char
  public :: parf_barrier
  public :: parf_get_noblock_char
  public :: parf_send_noblock_char
  public :: parf_wait_any_nostatus
  public :: getUnitDump

  include 'mpif.h'
  include 'precision.h'
  include 'messages.h'


contains

  subroutine initializeParallelismVariables()
    implicit none

    maxNodes = getMaxNodes(nameGroup)
    mpiCommGroup = getMpiCommGroup(nameGroup)
    mpiMasterProc = getMpiMasterProc(nameGroup)

  end subroutine initializeParallelismVariables


  subroutine msgOutMaster(header, message)
    !# Sends messages only in master processor of Chopping
    !# ---
    !# @info
    !# **Brief:** Sends messages only in master processor of Chopping. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jun/2019 <br>
    !# @endin
    character(len = *), intent(in) :: header, message

    call msgOutMasterGroup(header, message, nameGroup)
  end subroutine msgOutMaster


  subroutine msgDump(header, message)
    !# Dumps messages of Chopping
    !# ---
    !# @info
    !# **Brief:** Dumps messages of Chopping. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: sep/2019 <br>
    !# @endin  
    character(len = *), intent(in) :: header, message

    call msgDumpGroup(header, message)
  end subroutine msgDump


  function isMasterProc() result(isMaster)
    !# Returns if the Master Processor is running
    !# ---
    !# @info
    !# **Brief:** Returns if the Master Processor is running. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 <br>
    !# @endin    
    implicit none
    logical :: isMaster

    isMaster = getIsMasterProc(nameGroup)
  end function isMasterProc


    subroutine fatalError(header, message)
    !# Dumps error message
    !# ---
    !# @info
    !# **Brief:** Dumps error message everywhere in Chopping Group and destroy
    !# group parallelism. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: sep/2019 <br>
    !# @endin  
      character(len = *), intent(in) :: header
      character(len = *), intent(in) :: message
      
      call fatalErrorGroup(header, message, nameGroup, mpiCommGroup)
    end subroutine fatalError

end module Mod_Parallelism_Group_Chopping


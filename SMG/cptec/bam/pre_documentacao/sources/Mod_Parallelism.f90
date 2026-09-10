!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Parallelism </br></br>
!#
!# **Brief**: Deals with parallelism </br></br>
!# 
!# Exports: </br>
!# <ul type="disc">
!#  <li>MPI procedures for:
!#     <ul type="disc">
!#       <li>initialization. </li>
!#       <li>normal finalization. </li>
!#       <li>error finalization. </li>
!#     </ul>
!#  <li>MPI process data:
!#     <ul type="disc">
!#       <li>how many processes in MPI computation. </li>
!#       <li>which processes are these. </li>
!#     </ul>
!#  <li>Creates one file per process for printouts, unscrambling outputs. Provides
!#      procedures for:
!#     <ul type="disc">
!#       <li>writing messages to output and to the file. </li>
!#       <li>writing messages only to the file. </li>
!#     </ul>
!# </ul>
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
!# **Version**: 1.0.1
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2008 - Paulo Kubota   - version: 1.0.0 - from Chopping Parallel </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 1.0.1 - adds simplified routines
!#      that incapsulates MPI routines</li>
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

! TODO list
! - create getter (like getUnitDump) instead of public variables maxNodes, myId, mpiMasterProc


module Mod_Parallelism

  use Mod_String_Functions, only : &
    Replace_Text

  use Mod_FileManager, only : &
    getFileUnit, &
    openfile 

  use Mod_Messages, only : &
    msgOut &
    , msgDump

  implicit none

  private
  integer :: unitDump

  ! public data
  integer, public :: maxNodes      
  !# MPI processes in the computation
  integer, public :: myId          
  !# MPI process rank
  integer, public :: mpiMasterProc 
  !# Master processor - must set it


  ! public procedures
  public :: CreateParallelism
  public :: DestroyParallelism
  public :: parf_bcast_int_scalar
  public :: parf_bcast_int_1d
  public :: parf_bcast_char
  public :: parf_bcast_char_vec
  public :: parf_send_char
  public :: parf_get_char
  public :: parf_barrier
  public :: isMasterProc
  public :: getMyIdString
  public :: parf_get_noblock_char
  public :: parf_send_noblock_char
  public :: parf_wait_any_nostatus
  public :: getUnitDump
  public :: fatalError



  ! private data
  character(len = 4) :: cNProc      
  !# maxNodes in characters, for printing
  character(len = 4) :: cThisProc   
  !# myId in characters, for printing
  character(len = *), parameter :: header = "Paralellism           | "

  include 'mpif.h'
  include 'precision.h'

contains


  subroutine CreateParallelism()
    !# Creates Parallelism
    !# ---
    !# @info
    !# **Brief:** Initiates MPI Parallelism: number of processors and threads. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer :: ierror
    integer :: CompName, numThreads
    character(len = MPI_MAX_PROCESSOR_NAME) :: Pname
    character(len = 4) :: caux
    character(len = 14) :: nameDump

    !$    INTEGER,   EXTERNAL   ::  OMP_GET_MAX_THREADS

    ! inicialize MPI
    call MPI_INIT(ierror)
    if (ierror /= MPI_SUCCESS) then
      write(caux, "(i4)") ierror
      call fatalError(header, " MPI_INIT returns " // caux)
    end if

    ! number of MPI processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, maxNodes, ierror)
    if (ierror /= MPI_SUCCESS) then
      write(caux, "(i4)") ierror
      call fatalError(header, " MPI_COMM_SIZE returns " // caux)
    end if

    ! number of this process on MPI computation
    call MPI_COMM_RANK(MPI_COMM_WORLD, myId, ierror)
    if (ierror /= MPI_SUCCESS) then
      write(caux, "(i4)") ierror
      call fatalError(header, " MPI_COMM_RANK returns " // caux)
    end if

    ! process name
    call MPI_GET_PROCESSOR_NAME(Pname, CompName, ierror)
    if (ierror /= MPI_SUCCESS) then
      write(caux, "(i4)") ierror
      call fatalError(header, " MPI_GET_PROCESSOR_NAME retorna " // caux)
    end if

    ! OMP parallelism
    numThreads = 1
    !$   numThreads = OMP_GET_MAX_THREADS()

    ! Number of processes for file name
    write(cNProc, "(i4.4)") maxNodes * numThreads
    write(cThisProc, "(i4.4)") myId

    ! select unit for dumping
    ! TODO - use File manager to open
    !    unitDump = 30
    ! unitDump = getFileUnit()

    ! generate dump file name and open file
    nameDump = "Dump." // cThisProc // "." // cNProc
    unitDump = openFile(nameDump, 'formatted', 'sequential', -1, 'write', 'replace')

    ! Number of processes for file name
    write(cNProc, "(i4)") maxNodes
    write(cThisProc, "(i4)") myId
    write (caux, "(i4)") numThreads

    ! tell the world I'm alive

    if (isMasterProc()) then
      call msgOut(header, "Process " // trim(adjustl(cThisProc)) // " (" // &
        &Pname(1:CompName) // ") among " // trim(adjustl(cNProc)) // &
        &" processes with " // trim(adjustl(caux)) // &
        &" threads is alive")
      call msgDump(getUnitDump(), header, "Process " // trim(adjustl(cThisProc)) // " (" // &
        &Pname(1:CompName) // ") among " // trim(adjustl(cNProc)) // &
        &" processes with " // trim(adjustl(caux)) // &
        &" threads is alive")
    endif
  end subroutine CreateParallelism


  function getUnitDump() result(unitForDump)
    !# Gets Unit Dump
    !# ---
    !# @info
    !# **Brief:** Gets Unit Dump. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    integer :: unitForDump

    unitForDump = unitDump
  end function getUnitDump


  subroutine DestroyParallelism(message)
    !# Destroys Parallelism
    !# ---
    !# @info
    !# **Brief:** Destroys Parallelism. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: message
    integer ierror
    call msgOut(header, "Process " // trim(adjustl(cThisProc)) // " " // message)
    call msgDump(getUnitDump(), header, "Process " // trim(adjustl(cThisProc)) // " " // message)
    close(unitDump)
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    call MPI_FINALIZE(ierror)
  end subroutine DestroyParallelism
  

  subroutine parf_bcast_int_scalar(buff, source_host)
    !# MPI routine for broadcast integer
    !# ---
    !# @info
    !# **Brief:** MPI routine for broadcast integer. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: buff
    integer, intent(in) :: source_host
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_BCAST(buff, 1, MPI_INTEGER, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8)') rank, ierr
      call fatalError(header, "Error in parf_bcast_int_scalar: rank, ierr=" // &
        trim(string))
    endif
  end subroutine parf_bcast_int_scalar


  subroutine parf_bcast_int_1d(buff, buff_len, source_host)
    !# MPI routine for broadcast integer array
    !# ---
    !# @info
    !# **Brief:** MPI routine for broadcast integer array. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer(kind = p_i8), intent(in) :: buff_len
    integer, intent(inout) :: buff(buff_len)
    integer, intent(in) :: source_host
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_BCAST(buff, buff_len, MPI_INTEGER, source_host, &
      MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8)') rank, ierr
      call fatalError(header, "Error in parf_bcast_int_1d: rank, ierr=" // trim(string))
    endif
  end subroutine parf_bcast_int_1d


  subroutine parf_bcast_char(buff, buff_len, source_host)
    !# MPI routine for broadcast character
    !# ---
    !# @info
    !# **Brief:** MPI routine for broadcast character. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer(kind = p_i8), intent(in) :: buff_len
    integer, intent(in) :: source_host
    character(len = buff_len), intent(inout) :: buff
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_BCAST(buff, buff_len, MPI_CHARACTER, source_host, &
      MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8)') rank, ierr
      call fatalError(header, "Error in parf_bcast_char: rank, ierr=" // trim(string))
    endif
  end subroutine parf_bcast_char


  subroutine parf_bcast_char_vec(buff, buff_len, buff_size, source_host)
    !# MPI routine for broadcast character array
    !# ---
    !# @info
    !# **Brief:** MPI routine for broadcast character array. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer(kind = p_i8), intent(in) :: buff_len
    integer(kind = p_i8), intent(in) :: buff_size
    integer, intent(in) :: source_host
    character(len = buff_len), intent(inout) :: buff(buff_size)
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_BCAST(buff, buff_len * buff_size, MPI_CHARACTER, source_host, &
      MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8)') rank, ierr
      call fatalError(header, "Error in parf_bcast_char_vec: rank, ierr=" // &
        trim(string))
    endif
  end subroutine parf_bcast_char_vec


  subroutine parf_barrier(ibarrier)
    !# MPI routine for barrier
    !# ---
    !# @info
    !# **Brief:** MPI routine for barrier. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: ibarrier
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (ierr /= MPI_SUCCESS) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I6,X,I8)') rank, ibarrier, ierr
      call fatalError(header, "Error in parf_barrier: " // trim(string))
    endif
  end subroutine parf_barrier


  function isMasterProc() result(isMaster)
    !# Returns if the Master Processor is running
    !# ---
    !# @info
    !# **Brief:** Returns if the Master Processor is running. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    logical :: isMaster
    isMaster = (myId .eq. mpiMasterProc)
  end function isMasterProc


  subroutine parf_send_char(buff, buff_len, dest_host, tag)
    !# MPI routine for sending character to another processor
    !# ---
    !# @info
    !# **Brief:** MPI routine for sending character to another processor. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer(kind = p_i8), intent(in) :: buff_len
    character(len = buff_len), intent(inout) :: buff
    integer, intent(in) :: dest_host
    integer, intent(in) :: tag
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_send(buff, buff_len, MPI_CHARACTER, dest_host, tag, &
      MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8)') rank, ierr
      call fatalError(header, "Error in parf_send_char: rank, ierr=" // trim(string))
    endif
  end subroutine parf_send_char


  subroutine parf_get_char(buff, buff_len, source_host, tag)
    !# MPI routine for receiving character from another processor
    !# ---
    !# @info
    !# **Brief:** MPI routine for receiving character from another processor. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer(kind = p_i8), intent(in) :: buff_len
    character(len = buff_len), intent(out) :: buff
    integer, intent(in) :: source_host
    integer, intent(in) :: tag
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string
    integer :: status(MPI_STATUS_SIZE)

    call MPI_recv(buff, buff_len, MPI_CHARACTER, source_host, tag, &
      MPI_COMM_WORLD, status, ierr)

    if(ierr /= MPI_SUCCESS) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8)') rank, ierr
      call fatalError(header, "Error in parf_get_block: rank, ierr=" // trim(string))
    endif
  end subroutine parf_get_char


  function getMyIdString() result(myIdString)
    !# Returns current rank as String
    !# ---
    !# @info
    !# **Brief:** Returns current rank as String. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    character(len = 4) :: myIdString         ! myId in characters, for printing
    write(myIdString, "(i4)") myId
  end function getMyIdString


  subroutine parf_wait_any_nostatus(total, request, number)
    !# Calls MPI_Waitany
    !# ---
    !# @info
    !# **Brief:** Waits for any non blockind send (isend). </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    
    ! Arguments:
    integer, intent(inout) :: request(*)
    integer, intent(in) :: total
    integer, intent(out) :: number
    
    ! Local Variables:
    integer :: status(MPI_STATUS_SIZE)
    integer :: ierr

    call MPI_Waitany(total, request, number, status, ierr)
    if(ierr /= MPI_Success) call fatalError(header, "Error in parf_wait_any")

  end subroutine parf_wait_any_nostatus


  subroutine parf_get_noblock_char(buff, buff_len, source_host, tag, request)
    !# Calls MPI_Irecv, MPI_Comm_rank
    !# ---
    !# @info
    !# **Brief:** Waits for get non blockind char. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo

    ! Arguments:
    integer, intent(in) :: buff_len, source_host, tag
    character, intent(out) :: buff(buff_len)
    integer, intent(out) :: request
    
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_Irecv(buff, buff_len, MPI_CHARACTER, source_host, tag, &
      MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_Success) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8,X,I8)') rank, request, ierr
      call fatalError(header, "Error in parf_get_noblock: rank, request, ierr=" // trim(string))
    endif

  end subroutine parf_get_noblock_char


  subroutine parf_send_noblock_char(buff, buff_len, dest_host, tag, request)
    !# Calls MPI_Isend, MPI_Comm_rank
    !# ---
    !# @info
    !# **Brief:** Waits for send non blockind char. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    
    ! Arguments:
    integer, intent(in) :: buff_len, dest_host, tag
    character, intent(in) :: buff(buff_len)
    integer, intent(out) :: request
    
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len = 20) :: string

    call MPI_Isend(buff, buff_len, MPI_CHARACTER, dest_host, tag, &
      MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_Success) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
      write(string, FMT = '(I6.6,X,I8,X,I8)') rank, request, ierr
      call fatalError(header, "Error in parf_send_noblock: rank, request, ierr=" // &
        trim(string))
    endif

  end subroutine parf_send_noblock_char


  subroutine fatalError(header, message)
    !# Dump error message everywhere and destroy parallelism
    !# ---
    !# @info
    !# **Brief:** WDump error message everywhere and destroy parallelism. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: header 
    character(len = *), intent(in) :: message
    integer :: ierror = -1
    integer :: ierr
    character(len = *), parameter :: headerFatal = " ***(FATAL ERROR)*** "

    call msgOut(header,  headerFatal // message)
    call msgDump(getUnitDump(), header, headerFatal // message)
    call MPI_ABORT(MPI_COMM_WORLD, ierror, ierr)
    call MPI_Comm_free(MPI_COMM_WORLD)
    call MPI_FINALIZE(ierr)
    stop
  end subroutine fatalError

end module Mod_Parallelism


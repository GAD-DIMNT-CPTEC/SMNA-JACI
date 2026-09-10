!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Messages_Parallel </br></br>
!#
!# **Brief**: Module responsible for sending a warning message to Master 
!# Processor, messages only in master processor and dump messages only in 
!# master processor. </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.1 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti   - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita            - version: 1.1.1 </li>
!#  <li>12-09-2019 - Denis Eiras       - version: 2.0.0 </li>
!#  <li>27-11-2019 - Eduardo Khamis    - version: 2.0.1 </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time</li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time</li>
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

module Mod_Messages_Parallel
  
  use Mod_Messages, only : &
    msgOut &
    , msgWarningOut &
    , msgDump

  use Mod_Parallelism, only : &
    isMasterProc &
    , getUnitDump
    

  implicit none
  
  ! public procedures
  public :: msgWarningOutMaster
  public :: msgOutMaster
  public :: msgDumpMaster

 
  include 'messages.h'
  
  contains

  
  subroutine  msgWarningOutMaster(h, message, isDebug)
    !# Sends a warning message to Master Processor
    !# ---
    !# @info
    !# **Brief:** Sends a warning message to Master Processor. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: h, message
    logical, optional :: isDebug

    if (isMasterProc()) then
      call msgOut(h // "***** WARNING ***** ", message, isDebug)
    endif

  end subroutine msgWarningOutMaster


  subroutine msgOutMaster(h, message, isDebug)
    !# Sends messages only in master processor
    !# ---
    !# @info
    !# **Brief:** Sends messages only in master processor. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: h, message
    logical, optional :: isDebug

    if (isMasterProc()) then
      call msgOut(h, message, isDebug)
    endif
  end subroutine msgOutMaster


  subroutine msgDumpMaster(h, message)
    !# Dumps messages only in master processor
    !# ---
    !# @info
    !# **Brief:** Dumps messages only in master processor. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: h, message

    if (isMasterProc()) then
      call msgDump(getUnitDump(), h, message)
    endif
  end subroutine msgDumpMaster


end module Mod_Messages_Parallel

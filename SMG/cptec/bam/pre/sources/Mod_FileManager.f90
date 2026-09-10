!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_FileManager </br></br>
!#
!# **Brief**: Module responsible for file units generation and opening files </br></br>
!#
!# Task 6061 - using subroutines to return erros to Mod_Pre </br></br>
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
!# **Version**: 1.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>12-11-2019 - Denis Eiras     - version: 1.1.0 </li>
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

module Mod_FileManager

  use Mod_Messages, only : &
    msgWarningOut &
    , msgInLineFormatOut &
    , msgNewLine &
    , msgOut

  include "files.h"
  include "messages.h"
  private

  logical :: op                       
  !# io unit number opened or not
  public :: getFileUnit
  public :: openFile
  public :: fileExists
  public :: openFileErrorMsg

  character (len=*), parameter :: header = 'File Manager          |'
contains

  
  function getFileUnit(initialUnit) result(file_unit)
    !# Gets a fresh new file Unit
    !# ---
    !# @info
    !# **Brief:** Gets a fresh new file Unit for using in openFile method. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    implicit none
    integer, intent(in), optional :: initialUnit
    integer :: file_unit, firstUnitThisProc, lastUnitThisProc

    if(.not. present(initialUnit)) then
      firstUnitThisProc = firstUnit
    else
      firstUnitThisProc = firstUnit + initialUnit
    endif
    lastUnitThisProc = firstUnitThisProc + maxUnitsPerRank

    do file_unit = firstUnitThisProc, lastUnitThisProc
      inquire(file_unit, opened = op)
      if (.not. op) exit
    end do
    if (file_unit >= lastUnitThisProc) then
      call msgWarningOut(header, " all i/o units in use !")
      
    end if
  end function getFileUnit


  function openFile(fileName, openForm, openAcess, recSize, openAction, openStatus, initialUnit) result(fileUnit)
    !# Gets a fresh new file Unit for opening files for reading or writing
    !# ---
    !# @info
    !# **Brief:** Gets a fresh new file Unit for opening files for reading or
    !# writing. It sends a error message in case of error when opening </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    implicit none
    character (len = *), intent(in) :: fileName
    !# fileName is the path +name of the file
    character (len = *), intent(in) :: openForm
    !# openForm form  formatted or unformatted
    character (len = *), intent(in) :: openAcess
    !# openAcess acess sequential or direct
    integer, intent(in)             :: recSize
    !# recSize the size of the record, in case of direct acess. When using sequential, use -1
    character (len = *), intent(in) :: openAction
    !# openAction Action read or write
    character (len = *), intent(in) :: openStatus
    !# openStatus status old for reading, replace for writing and so on
    integer, intent(in), optional   :: initialUnit
    integer :: fileUnit, ios

    ! call msgOut(header, ' Opening file ' // trim(fileName) // ' for ' // openAction)
    if(openAction == 'read' .and. .not. fileExists(fileName)) then
      call msgWarningOut(header, 'File not found: ' // fileName)
      fileUnit = -1
      return
    endif

    fileUnit = getFileUnit(initialUnit)
    if(recSize .eq. -1) then
      open (file = fileName, unit = fileUnit, form = openForm, access = openAcess, &
        action = openAction, status = openStatus, iostat = ios)
    else
      open (file = fileName, unit = fileUnit, form = openForm, access = openAcess, recl = recSize, &
        action = openAction, status = openStatus, iostat = ios)
    end if

    call openFileErrorMsg(ios, fileName)
    if (ios /= 0) then
      fileUnit = -1
    endif

  end function openFile


  subroutine openFileErrorMsg(ios, fileName)
    !# Sends Open Files Error Message
    !# ---
    !# @info
    !# **Brief:** Sends a message error in case of error when opening. </br>
    !# &bull; Jose P. Bonatti </br>
    !# **Date**: nov/2008 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: ios
    character (len = *), intent(in) :: fileName

    if (ios /= 0) then
      call msgInLineFormatOut(header // ' Error opening file ' // trim(fileName) // ' Error code: ', '(A)')
      call msgInLineFormatOut(ios, '(I6)')
      call msgNewLine()
    end if

  end subroutine openFileErrorMsg


  function fileExists(filename) result(res)
    !# Checks if file exists
    !# ---
    !# @info
    !# **Brief:** Checks for file existence. </br>
    !# &bull; Jose P. Bonatti </br>
    !# **Date**: nov/2008 </br>
    !# @endinfo
    implicit none
    character(len = *), intent(in) :: filename
    logical :: res

    ! Check if the file exists
    inquire(file = trim(filename), exist = res)
  end function fileExists


end module Mod_FileManager

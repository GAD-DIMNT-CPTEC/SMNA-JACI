!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Messages </br></br>
!#
!# **Brief**: Module responsible for sending newline to unit in all processors, 
!# messages to unit and to output in all processors using character formatter, 
!# without advancing line, sets on/off debug messages, dump messsage only at 
!# dump file and warning message. </br></br>
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
!#  <li>01-04-2018 - Daniel M. Lamosa  - version: 2.0.0 </li>
!#  <li>12-09-2019 - Denis Eiras       - version: 2.0.1 </li>
!#  <li>25-09-2019 - Eduardo Khamis    - version: 2.0.1 </li>
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

module Mod_Messages

  use Mod_String_Functions, only : &
    Replace_Text


  implicit none

  ! public procedures

  public :: msgWarningOut
  public :: msgOut
  public :: msgDump
  public :: setDebugMode
  public :: msgInLineFormatOut
  public :: msgInLineFormatToUnit
  public :: msgNewLine
  public :: msgNewLineToUnit
    
  private
  logical :: isDebugMode
  
  include 'messages.h'
  include 'precision.h'

  interface msgInLineFormatToUnit
    module procedure msgInLineFormatToUnit_Real_8
    module procedure msgInLineFormatToUnit_Integer
    module procedure msgInLineFormatToUnit_IntegerArray
    module procedure msgInLineFormatToUnit_Char
  end interface

  interface msgInLineFormatOut
    module procedure msgInLineFormatOut_Real_8
    module procedure msgInLineFormatOut_Integer
    module procedure msgInLineFormatOut_IntegerArray
    module procedure msgInLineFormatOut_Char
  end interface

  contains


  subroutine msgNewLineToUnit(writeUnit, isDebug)
    !# Sends newline to unit in all processors
    !# ---
    !# @info
    !# **Brief:** Sends newline to unit in all processors. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: writeUnit
    logical, optional :: isDebug

    ! show only debug messages, if isDebugMode is on and isDebug parameter is on
    if(present(isDebug)) then
      if(.not. isDebugMode .or. .not. isDebug) return
    endif

    write(writeUnit, fmt='(A)', advance='yes') ''

  end subroutine msgNewLineToUnit


  subroutine msgNewLine(isDebug)
    !# Sends newline to out in all processors
    !# ---
    !# @info
    !# **Brief:** Sends newline to out in all processors. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    logical, optional :: isDebug

    call msgNewLineToUnit(p_nfprt, isDebug)

  end subroutine msgNewLine


  ! ================================== interface msgInLineFormatToUnit ==================================

  subroutine msgInLineFormatToUnit_Real_8(writeUnit, messageReal_8, formatString,  isDebug)
    !# Sends messages to unit in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to unit in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: writeUnit
    real (kind = p_r8), intent(in) :: messageReal_8
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    ! show only debug messages, if isDebugMode is on and isDebug parameter is on
    if(present(isDebug)) then
      if(.not. isDebugMode .or. .not. isDebug) return
    endif

    write(writeUnit, fmt=trim(formatString), advance='no') messageReal_8

  end subroutine msgInLineFormatToUnit_Real_8

  
  subroutine msgInLineFormatToUnit_Integer(writeUnit, messageInteger, formatString,  isDebug)
    !# Sends messages to unit in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to unit in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: writeUnit
    integer, intent(in) :: messageInteger
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    ! show only debug messages, if isDebugMode is on and isDebug parameter is on
    if(present(isDebug)) then
      if(.not. isDebugMode .or. .not. isDebug) return
    endif

    write(writeUnit, fmt=trim(formatString), advance='no') messageInteger

  end subroutine msgInLineFormatToUnit_Integer

  
  subroutine msgInLineFormatToUnit_IntegerArray(writeUnit, messageIntegerArray, formatString,  isDebug)
    !# Sends messages to unit in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to unit in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: writeUnit
    integer, intent(in) :: messageIntegerArray(:)
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    ! show only debug messages, if isDebugMode is on and isDebug parameter is on
    if(present(isDebug)) then
      if(.not. isDebugMode .or. .not. isDebug) return
    endif

    write(writeUnit, fmt=trim(formatString), advance='no') messageIntegerArray

  end subroutine msgInLineFormatToUnit_IntegerArray


  subroutine msgInLineFormatToUnit_Char(writeUnit, messageChar, formatString,  isDebug)
    !# Sends messages to unit in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to unit in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: writeUnit
    character(len = *), intent(in) :: messageChar
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    ! show only debug messages, if isDebugMode is on and isDebug parameter is on
    if(present(isDebug)) then
      if(.not. isDebugMode .or. .not. isDebug) return
    endif

    write(writeUnit, fmt=trim(formatString), advance='no') trim(messageChar)

  end subroutine msgInLineFormatToUnit_Char

  
  ! ================================== interface msgInLineFormatOut ==================================

  subroutine msgInLineFormatOut_Real_8(messageReal_8, formatString,  isDebug)
    !# Sends messages to output in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to output in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    real (kind = p_r8), intent(in) :: messageReal_8
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    call msgInLineFormatToUnit(p_nfprt, messageReal_8, formatString,  isDebug)

  end subroutine msgInLineFormatOut_Real_8


  subroutine msgInLineFormatOut_Integer(messageInteger, formatString,  isDebug)
    !# Sends messages to output in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to output in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: messageInteger
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    call msgInLineFormatToUnit(p_nfprt, messageInteger, formatString,  isDebug)

  end subroutine msgInLineFormatOut_Integer


  subroutine msgInLineFormatOut_IntegerArray(messageIntegerArray, formatString,  isDebug)
    !# Sends messages to output in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to output in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: messageIntegerArray(:)
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    call msgInLineFormatToUnit(p_nfprt, messageIntegerArray, formatString,  isDebug)

  end subroutine msgInLineFormatOut_IntegerArray

  
  subroutine msgInLineFormatOut_Char(messageChar, formatString,  isDebug)
    !# Sends messages to output in all processors using character formatter, without
    !# advancing line
    !# ---
    !# @info
    !# **Brief:** Sends messages to output in all processors using character formatter,
    !# without advancing line. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: messageChar
    character(len = *), intent(in) :: formatString
    logical, optional :: isDebug

    call msgInLineFormatToUnit(p_nfprt, messageChar, formatString,  isDebug)

  end subroutine msgInLineFormatOut_Char


  subroutine setDebugMode(isDebug)
    !# Sets on/off debug messages
    !# ---
    !# @info
    !# **Brief:** Sets on/off debug messages. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    logical, intent(in) :: isDebug

    isDebugMode = isDebug
  end subroutine setDebugMode


  subroutine msg(writeUnit, h, message, isDebug)
    !# Internal module routine for send messages
    !# ---
    !# @info
    !# **Brief:** Internal module routine for send messages. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: writeUnit
    character(len = *), intent(in) :: h, message
    logical, optional :: isDebug

    integer :: thisThread
    integer, parameter :: maxLineLen = 160 ! stdout record size
    character(len = len(h) + len(message) + 2) :: fullMsg !full message
    integer :: fullMsgLen, lineCount, lines, first, last, idxNewLine
    !$ INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

    ! show only debug messages, if isDebugMode is on and isDebug parameter is on
    if(present(isDebug)) then
      if(.not. isDebugMode .or. .not. isDebug) return
    endif

    thisThread = 0
    !$ thisThread = OMP_GET_THREAD_NUM()
    if (thisThread == 0) then
      fullMsg = trim(h) // ' ' // trim(message)
      fullMsgLen = len_trim(fullMsg)
      lineCount = fullMsgLen / maxLineLen
      if(lineCount * maxLineLen < fullMsgLen)then
        lineCount = lineCount + 1
      end if
      do lines = 1, lineCount
        first = (lines - 1) * maxLineLen + 1
        last = min(lines * maxLineLen, fullMsgLen)

        idxNewLine = index(fullMsg(first:last), "$")
        if (idxNewLine .ne. 0) then
          write(writeUnit, *) fullMsg(first:idxNewLine - 2)
          write(writeUnit, *)
          write(writeUnit, *) fullMsg(idxNewLine + 1:last)
        else
          write(writeUnit, *) fullMsg(first:last)
        end if
        ! CALL flush(writeUnit)
      end do
    end if
  end subroutine msg


  subroutine msgDump(unitDump, h, message)
    !# Dumps messsage only at dump file
    !# ---
    !# @info
    !# **Brief:** Dumps messsage only at dump file. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: unitDump
    character(len = *), intent(in) :: h, message

    call msg(unitDump, h, message)
  end subroutine msgDump


  subroutine msgOut(h, message, isDebug)
    !# Sends messages in all processors
    !# ---
    !# @info
    !# **Brief:** Sends messages in all processors. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: h, message
    logical, optional :: isDebug

    call msg(p_nfprt, h, message, isDebug)
  end subroutine msgOut

  
  subroutine  msgWarningOut(h, message, isDebug)
    !# Sends a warning message
    !# ---
    !# @info
    !# **Brief:** Sends a warning message. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    character(len = *), intent(in) :: h, message
    logical, optional :: isDebug

    call msg(p_nfprt, h, "***** WARNING ***** " // message, isDebug)
  end subroutine msgWarningOut


end module Mod_Messages

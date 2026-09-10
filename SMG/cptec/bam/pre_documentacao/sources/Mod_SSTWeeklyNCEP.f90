!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SSTWeeklyNCEP </br></br>
!#
!# **Brief**: Module responsible for correcting the weekly sea surface temperature (SST)
!! by topography and generate the sea ice mask (1: Sea Ice 0: No Ice) </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/datain/gdas1.T00Z.sstgrd.YYYYMMDDHHS (Ex.: pre/datain/gdas1.T00Z.sstgrd.2015043000)
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; pre/dataout/SSTWeekly.YYYYMMDD.ctl  (Ex.: pre/dataout/SSTWeekly.20150430.ctl) </br>
!# &bull; pre/dataout/SSTWeekly.YYYYMMDD
!# </br></br>
!# 
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.1 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti     - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita              - version: 1.1.1 </li>
!#  <li>08-05-2019 - Eduardo Khamis      - version: 2.0.1 </li>
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

module Mod_SSTWeeklyNCEP

  use Mod_FileManager
  use Mod_Namelist,     only : varCommonNameListData

  ! First Point of Initial Data is near South Pole and near Greenwhich
  ! First Point of Output  Data is near North Pole and near Greenwhich

  implicit none
  private
  include 'pre.h'
  include 'files.h'
  include 'precision.h'
  include 'messages.h'

  public :: getNameSSTWeeklyNCEP, initSSTWeeklyNCEP, generateSSTWeeklyNCEP, shouldRunSSTWeeklyNCEP

  integer :: month
  character (len = 12) :: timeGrADS = '  Z         '
  character (len = 3), dimension(12) :: monthChar = &
    (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
      'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC' /)

  real (kind = p_r4), dimension (:, :), allocatable :: sst

  type SSTWeeklyNCEPNameListData
    logical :: autoDims = .true. 
    !# Use auto dimensions accordding to fileSST filename
    integer :: xDim = 4320  
    !# Number of Longitudes
    integer :: yDim = 2160  
    !# Number of Latitudes
    character (len = maxPathLength) :: varName = 'SSTWeekly'         
    !# sst output file name prefix
    character (len = maxPathLength) :: fileSST = 'gdas1.T  Z.sstgrd' 
    !# sst input file name prefix
  end type SSTWeeklyNCEPNameListData

! TODO which one of above should we use?
!  namelist /SSTWeeklyNCEPNameList/ xDim, yDim, grads, date, dirPreIn, dirPreOut ! have to do all these in the namelist
!  namelist /SSTWeeklyNCEPNameList/ date, dirPreIn, dirPreOut, xDim, yDim, grads ! looks like this works too
!  namelist /SSTWeeklyNCEPNameList/ xDim, yDim, grads

  type(varCommonNameListData)     :: varCommon
  type(SSTWeeklyNCEPNameListData) :: var
  namelist /SSTWeeklyNCEPNameList/   var


contains


  function getNameSSTWeeklyNCEP() result(returnModuleName)
    !# Returns SSTWeeklyNCEP Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SSTWeeklyNCEP Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    returnModuleName = "SSTWeeklyNCEP"
  end function getNameSSTWeeklyNCEP


  subroutine initSSTWeeklyNCEP(nameListFileUnit, varCommon_)
    !# Initialization of SSTWeeklyNCEP module
    !# ---
    !# @info
    !# **Brief:** Initialization of SSTWeeklyNCEP module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
    
    read(unit = nameListFileUnit, nml = SSTWeeklyNCEPNameList)
    varCommon = varCommon_

  end subroutine initSSTWeeklyNCEP


  function shouldRunSSTWeeklyNCEP() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunSSTWeeklyNCEP


  function getOutFileName() result(sstWeeklyNCEPOutFilename)
    !# Gets SSTWeeklyNCEP Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SSTWeeklyNCEP Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: sstWeeklyNCEPOutFilename

    sstWeeklyNCEPOutFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.' // varCommon%date(1:8)
  end function getOutFileName


  function generateSSTWeeklyNCEP() result(isExecOk)
    !# Generates SSTWeeklyNCEP
    !# ---
    !# @info
    !# **Brief:** Generates SSTWeeklyNCEP. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
!    integer, intent(in) :: xMax_    !max value of longitude
!    integer, intent(in) :: yMax_    !max value of latitude
    integer :: inFileUnit, outFileUnit, outRecSize
    character (len = maxPathLength) :: fileSSTInFilename, sstOutFilename
    logical :: isExecOk
    integer :: sstYear, sstMon, sstDay

    isExecOk = .false.

    if (var%autoDims) then

      read(varCommon%date(1:4), *) sstYear
      read(varCommon%date(5:6), *) sstMon
      read(varCommon%date(7:8), *) sstDay

      if ( (sstYear < 2015) .or. ( sstYear == 2015 .and. sstMon == 1 .and. sstDay <= 14 )) then
        var%xDim = 360
        var%yDim = 180
      else
        var%xDim = 4320
        var%yDim = 2160
      endif
    endif

    
    allocate (sst(var%xDim, var%yDim))

    var%fileSST(8:9) = varCommon%date(9:10)

    fileSSTInFilename = trim(varCommon%dirPreIn) // trim(var%fileSST) // '.' // varCommon%date
    inFileUnit = openFile(trim(fileSSTInFilename), 'unformatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return
    read  (unit = inFileUnit) sst
    close (unit = inFileUnit)

    inquire (iolength = outRecSize) sst
    sstOutFilename = getOutFileName()
    outFileUnit = openFile(trim(sstOutFilename), 'unformatted', 'direct', outRecSize, 'write', 'replace')
    if(outFileUnit < 0) return
    write (unit = outFileUnit, REC = 1) sst
    close (unit = outFileUnit)

    if (varCommon%grads) then
      call generateGrads()
    end if

    deallocate (sst)

    isExecOk = .true.
  end function generateSSTWeeklyNCEP


  subroutine generateGrads()
    !# Generates Grads
    !# ---
    !# @info
    !# **Brief:** Generates Grads. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxPathLength) :: ctlPathFileName
    character (len = maxPathLength) :: gradsBaseName
    character (len = maxPathLength) :: ctlFileName
    integer :: ctlFileUnit      ! Temp File Unit

    timeGrADS(1:2) = varCommon%date(9:10)
    timeGrADS(4:5) = varCommon%date(7:8)
    timeGrADS(9:12) = varCommon%date(1:4)
    read (varCommon%date(5:6), fmt = '(I2)') month
    timeGrADS(6:8) = monthChar(month)

    gradsBaseName = trim(var%varName) // '.' // varCommon%date(1:8)
    ctlFileName = trim(gradsBaseName) // '.ctl'
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(ctlFileName)
    ctlFileUnit = openFile(ctlPathFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if(ctlFileUnit < 0) return
    write (unit = ctlFileUnit, fmt = '(A)') 'DSET ^' // &
      trim(varCommon%dirPreOut) // trim(var%varName) // '.' // varCommon%date(1:8)
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'UNDEF -999.0'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'TITLE NCEP Weekly SST'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
      'XDEF ', var%xDim, ' LINEAR ', 0.5_p_r4, 360.0_p_r4 / real(var%xDim, p_r4)
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
      'YDEF ', var%yDim, ' LINEAR ', -89.5_p_r4, 180.0_p_r4 / real(var%yDim - 1, p_r4)
    write (unit = ctlFileUnit, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
    write (unit = ctlFileUnit, fmt = '(A)') 'TDEF 1 LINEAR ' // timeGrADS // ' 6HR'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'VARS 1'
    write (unit = ctlFileUnit, fmt = '(A)') 'SSTW 0 99 NCEP Weekly SST [K]'
    write (unit = ctlFileUnit, fmt = '(A)') 'ENDVARS'
    close (unit = ctlFileUnit)
  end subroutine generateGrads


  subroutine printNameList()    
    !# Prints NameList
    !# ---
    !# @info
    !# **Brief:** Prints NameList. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo

    implicit none
    write (unit = p_nfprt, fmt = '(/,A)')  ' &SSTWeeklyNCEPNameList'
    write (unit = p_nfprt, fmt = '(A,I6)') '    Idim = ', var%xDim
    write (unit = p_nfprt, fmt = '(A,I6)') '    Jdim = ', var%yDim
    write (unit = p_nfprt, fmt = '(A,L6)') '   GrADS = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)')    '    Date = ' // varCommon%date
    write (unit = p_nfprt, fmt = '(A,/)')  ' /'
  end subroutine printNamelist

end module Mod_SSTWeeklyNCEP

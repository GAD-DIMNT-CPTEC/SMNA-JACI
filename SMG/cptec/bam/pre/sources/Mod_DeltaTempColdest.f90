!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_DeltaTempColdest </br></br>
!#
!# **Brief**: Linear (.true.) or area weighted (.false.) Interpolation interpolate
!# regular to gaussian regular input data is assumed to be oriented with the north
!# pole and greenwich as the first point gaussian output data is interpolated to
!# be oriented with the north pole and greenwich as the first point input for the
!# AGCM is assumed ieee 32 bits big endian. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/DeltaTempColdestClima.dat </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/DeltaTempColdest.G00450 </br>
!# &bull; pre/dataout/DeltaTempColdest.G00450 </br>
!# &bull; pre/dataout/DeltaTempColdest.G00450.ctl </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 <br><br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti  - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita           - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel M. Lamosa - version: 2.0.0 </li>
!#  <li>31-01-2020 - Denis Eiras      - version: 2.1.0 </li>
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
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2> <br>
!# **&#9993;**<mailto:atende.cptec@inpe.br> <br><br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---

module Mod_DeltaTempColdest

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_LinearInterpolation, only : gLatsL => LatOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, initAreaInterpolation, doAreaInterpolation

  implicit none

  public generateDeltaTempColdest
  public initDeltaTempColdest
  public getNameDeltaTempColdest
  public shouldRunDeltaTempColdest

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'

  !input variables
  type DeltaTempColdestNameListData
    integer :: xDim  
    !# Number of Longitudes
    integer :: yDim  
    !# Number of Latitudes
    character(len = maxPathLength) :: varName = 'DeltaTempColdestClima'  
    !# input prefix file name
    character(len = maxPathLength) :: varNameOut = 'DeltaTempColdes'     
    !# output prefix file name
    logical :: linear ! flag for interpolation type (true = linear, false = area)
  end type DeltaTempColdestNameListData

  type(varCommonNameListData) :: varCommon
  type(DeltaTempColdestNameListData) :: var
  namelist /DeltaTempColdestNameList/ var

  !internal variables
  real(kind = p_r8), parameter :: p_Lat0 = 90.0_p_r8 
  !# Start at North Pole
  real(kind = p_r8), parameter :: p_Lon0 = 0.0_p_r8  
  !# Start at Prime Meridian
  character(len = 7) :: nLats = '.G     '                  
  !# posfix land sea mask file
  real(kind = p_r4), dimension(:, :, :), allocatable :: deltaTempColdestIn     !
  real(kind = p_r8), dimension(:, :, :), allocatable :: deltaTempColdestInput  !
  real(kind = p_r4), dimension(:, :, :), allocatable :: deltaTempColdestOut    !
  real(kind = p_r8), dimension(:, :, :), allocatable :: deltaTempColdestOutput !

contains


  function getNameDeltaTempColdest() result(returnModuleName)
    !# Returns DeltaTempColdest Module Name
    !# ---
    !# @info
    !# **Brief:** Returns DeltaTempColdest Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName ! variable for store module name

    returnModuleName = "DeltaTempColdest"
  end function getNameDeltaTempColdest


  subroutine initDeltaTempColdest(nameListFileUnit, varCommon_)
    !# Initializes DeltaTempColdest module
    !# ---
    !# @info
    !# **Brief:** Initializes DeltaTempColdest module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit ! file unit of namelist PRE_run.nml
    type(varCommonNameListData), intent(in) :: varCommon_ ! variable of type varCommonNameListData for managing common variables at PRE_run.nml

    read(unit = nameListFileUnit, nml = DeltaTempColdestNameList)
    varCommon = varCommon_

    write(nLats(3:7), '(i5.5)') varCommon%yMax
  end subroutine initDeltaTempColdest


  function shouldRunDeltaTempColdest() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun ! result variable for store if method should run

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunDeltaTempColdest


  function getOutFileName() result(outFileName)
    !# Gets DeltaTempColdest output file name
    !# ---
    !# @info
    !# **Brief:** Gets DeltaTempColdest output file name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: outFileName

    outFileName = trim(varCommon%dirModelIn) // trim(var%varNameOut) // nLats
  end function getOutFileName


  subroutine allocateData()
    !# Allocates data
    !# ---
    !# @info
    !# **Brief:** Allocates matrixes based in namelist. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    allocate(deltaTempColdestIn(var%xDim, var%yDim, 1))
    allocate(deltaTempColdestInput(var%xDim, var%yDim, 1))
    allocate(deltaTempColdestOut(varCommon%xMax, varCommon%yMax, 1))
    allocate(deltaTempColdestOutput(varCommon%xMax, varCommon%yMax, 1))
  end subroutine allocateData


  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Deallocates matrixes. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    deallocate(deltaTempColdestIn)
    deallocate(deltaTempColdestInput)
    deallocate(deltaTempColdestOut)
    deallocate(deltaTempColdestOutput)
  end subroutine deallocateData


  function readDeltaTempColdest() result(isReadOk)
    !# Reads Delta Temp Coldest
    !# ---
    !# @info
    !# **Brief:** Reads Delta Temp Coldest, input file. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none

    logical :: isReadOk 
    !# return value (if function was executed sucessfully)
    integer :: inRecSize 
    !# size of file record
    character (len = maxPathLength) :: inFilename ! input filename
    integer :: inFileUnit

    isReadOk = .false.
    inquire (iolength = inRecSize) deltaTempColdestIn

    inFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return

    read(unit = inFileUnit, rec = 1) deltaTempColdestIn
    deltaTempColdestInput = real(deltaTempColdestIn, p_r8)
    close(unit = inFileUnit)
    isReadOk = .true.
  end function readDeltaTempColdest


  function writeDeltaTempColdest(outFilename) result(isWriteOk)
    !# Writes Delta Temp Coldest
    !# ---
    !# @info
    !# **Brief:** Generates ZoRL Input File for the BAM Model. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isWriteOk  
    !# return value (if function was executed sucessfully)
    character(len=*), intent(in) :: outFilename
    integer :: lRecOut  
    !# size of file record
    integer :: datFileUnit

    isWriteOk = .false.
    inquire(iolength = lRecOut) deltaTempColdestOut
    datFileUnit = openFile(trim(outFilename), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(datFileUnit < 0) return

    write(unit = datFileUnit, rec = 1) deltaTempColdestOut
    close(unit = datFileUnit)
    isWriteOk = .true.

  end function writeDeltaTempColdest


  function generateGrads() result(isGradsOk)
    !# Generates Grads files
    !# ---
    !# @info
    !# **Brief:** Generates .ctl file to check output. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isGradsOk ! return value (if function was executed sucessfully)
    character (len = maxPathLength) :: ctlPathFileName ! ctl Grads File name with path
    character (len = maxPathLength) :: gradsBaseName ! dat GrADS base name without extension
    integer :: gradsFileUnit ! File Unit for all grads files

    ! dat file is the same of getOutputFile; not generated here
    isGradsOk = .false.
    gradsBaseName = trim(var%varNameOut) // nLats
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.ctl'

    if(.not. writeDeltaTempColdest(trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.dat')) return

    ! Write .ctl file -------------------------------------------------------------
    gradsFileUnit = openFile(trim(ctlPathFileName), 'formatted', 'sequential', -1, 'write', 'replace')
    if(gradsFileUnit < 0) return

    write(unit = gradsFileUnit, fmt = '(a)') 'DSET ^' // trim(gradsBaseName) // '.dat'
    write(unit = gradsFileUnit, fmt = '(a)') '*'
    write(unit = gradsFileUnit, fmt = '(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit = gradsFileUnit, fmt = '(a)') '*'
    write(unit = gradsFileUnit, fmt = '(a)') 'UNDEF -999.0'
    write(unit = gradsFileUnit, fmt = '(a)') '*'
    write(unit = gradsFileUnit, fmt = '(a)') 'TITLE DeltaTempColdest on a Gaussian Grid'
    write(unit = gradsFileUnit, fmt = '(a)') '*'
    write(unit = gradsFileUnit, fmt = '(a,i5,a,f8.3,f15.10)') 'XDEF ', varCommon%xMax, ' LINEAR ', &
      0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write(unit = gradsFileUnit, fmt = '(a,i5,a)') 'YDEF ', varCommon%yMax, ' LEVELS '
    if(var%linear) then
      write(unit = gradsFileUnit, fmt = '(8f10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write(unit = gradsFileUnit, fmt = '(8f10.5)') gLatsA(varCommon%yMax:1:-1)
    end if
    write(unit = gradsFileUnit, fmt = '(a)') 'ZDEF 1 LEVELS 1000'
    write(unit = gradsFileUnit, fmt = '(a)') 'TDEF 12 LINEAR JAN2005 1MO'
    write(unit = gradsFileUnit, fmt = '(a)') 'VARS 1'
    write(unit = gradsFileUnit, fmt = '(a)') 'deltat 0 99 DeltaTempColdest [C]'
    write(unit = gradsFileUnit, fmt = '(a)') 'ENDVARS'
    close(unit = gradsFileUnit)
    isGradsOk = .true.

  end function generateGrads


  function generateDeltaTempColdest() result(isExecOk)
    !# Generates Delta Temp Coldest
    !# ---
    !# @info
    !# **Brief:** Generates Delta Temp Coldest output. This subroutine is the main
    !3 method for use this module.  </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk ! return value (if function was executed sucessfully)
    integer :: k
    logical :: flagInput(5), flagOutput(5)

    isExecOk = .false.
    flagInput(1) = .true.   ! Start at North Pole
    flagInput(2) = .true.   ! Start at Prime Meridian
    flagInput(3) = .true.   ! Latitudes Are at Center of Box
    flagInput(4) = .true.   ! Longitudes Are at Center of Box
    flagInput(5) = .false.  ! Regular Grid
    flagOutput(1) = .true.  ! Start at North Pole
    flagOutput(2) = .true.  ! Start at Prime Meridian
    flagOutput(3) = .false. ! Latitudes Are at North Edge of Box
    flagOutput(4) = .true.  ! Longitudes Are at Center of Box
    flagOutput(5) = .true.  ! Gaussian Grid

    if(var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_Lat0, p_Lon0)
    else
      call initAreaInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if
    call allocateData()
    if(.not.readDeltaTempColdest()) return

    ! Interpolate Input Regular Grid Deep Soil DeltaTempColdest To Gaussian Grid on Output
    do k = 1, 1
      if(var%linear) then
        call doLinearInterpolation(deltaTempColdestInput(:, :, k), deltaTempColdestOutput(:, :, k))
      else
        call doAreaInterpolation(deltaTempColdestInput(:, :, k), deltaTempColdestOutput(:, :, k))
      end if
      deltaTempColdestOut(:, :, k) = real(deltaTempColdestOutput(:, :, k), p_r4)
    end do

    if(.not. writeDeltaTempColdest(getOutFileName())) return

    if(varCommon%grads) then
      if(.not. generateGrads()) return
    end if

    call deallocateData()
    isExecOk = .true.

  end function generateDeltaTempColdest


end module Mod_DeltaTempColdest

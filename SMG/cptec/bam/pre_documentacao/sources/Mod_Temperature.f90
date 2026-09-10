!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Temperature </br></br>
!#
!# **Brief**: Module responsible for reading the data from the
!# TemperatureClima.dat file and interpolating this data in the grid and
!# resolution used by the BAM model </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/dataout/TemperatureClima.dat
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; model/datain/Temperature.GZZZZZ (Ex.: model/datain/Temperature.G00450) </br>
!# &bull; pre/dataout/Temperature.GZZZZZ </br>
!# &bull; pre/dataout/Temperature.GZZZZZ.ctl
!# </br></br>
!# 
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti   - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita            - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel M. Lamosa  - version: 2.0.0 </li>
!#  <li>05-02-2020 - Denis Eiras       - version: 2.1.0 </li>
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

module Mod_Temperature

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_LinearInterpolation, only : gLatsL => LatOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, initAreaInterpolation, doAreaInterpolation

  implicit none

  public generateTemperature
  public initTemperature
  public getNameTemperature
  public shouldRunTemperature

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'

  !input variables
  type TemperatureNameListData
    integer :: xDim  
    !# Number of Longitudes
    integer :: yDim  
    !# Number of Latitudes
    character(len = maxPathLength) :: varName = 'TemperatureClima'  
    !# input prefix file name
    character(len = maxPathLength) :: varNameOut = 'Temperature'    
    !# output prefix file name
    logical :: linear = .true.  
    !# flag for interpolation type (true = linear, false = area)
  end type TemperatureNameListData

  type(varCommonNameListData) :: varCommon
  type(TemperatureNameListData) :: var
  namelist /TemperatureNameList/ var

  !internal variables
  real(kind = p_r8), parameter :: p_Lat0 = 90.0_p_r8 
  !# Start at North Pole
  real(kind = p_r8), parameter :: p_Lon0 = 0.0_p_r8  
  !# Start at Prime Meridian
  character(len = 7) :: nLats = '.G     '                  
  !# posfix land sea mask file
  logical :: polarMean
  logical :: flaginput(5), flagoutput(5)
  integer, dimension (:, :), allocatable, public :: maskinput

  real(kind = p_r4), dimension(:, :, :), allocatable :: temperatureIn     
  real(kind = p_r8), dimension(:, :, :), allocatable :: temperatureInput  
  real(kind = p_r4), dimension(:, :, :), allocatable :: temperatureOut    
  real(kind = p_r8), dimension(:, :, :), allocatable :: temperatureOutput 

contains


  function getNameTemperature() result(returnModuleName)
    !# Returns Temperature Module Name
    !# ---
    !# @info
    !# **Brief:** Returns Temperature Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName 
    !# variable for store module name

    returnModuleName = "Temperature"
  end function getNameTemperature


  subroutine initTemperature(nameListFileUnit, varCommon_)
    !# Initialization of Temperature module
    !# ---
    !# @info
    !# **Brief:** Initialization of Temperature module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    
    implicit none
    integer, intent(in) :: nameListFileUnit 
    !# file unit of namelist PRE_run.nml
    type(varCommonNameListData), intent(in) :: varCommon_ 
    !# variable of type varCommonNameListData for managing common variables at PRE_run.nml

    read(unit = nameListFileUnit, nml = TemperatureNameList)
    varCommon = varCommon_

    write(nLats(3:7), '(i5.5)') varCommon%yMax
    ! For Area Weighted Interpolation
    allocate (maskInput(var%xDim, var%yDim))
    maskinput = 1
    polarmean = .false.
    flaginput(1) = .true.   ! start at north pole
    flaginput(2) = .true.   ! start at prime meridian
    flaginput(3) = .true.   ! latitudes are at center of box
    flaginput(4) = .true.   ! longitudes are at center of box
    flaginput(5) = .false.  ! regular grid
    flagoutput(1) = .true.  ! start at north pole
    flagoutput(2) = .true.  ! start at prime meridian
    flagoutput(3) = .false. ! latitudes are at north edge of box
    flagoutput(4) = .true.  ! longitudes are at center of box
    flagoutput(5) = .true.  ! gaussian grid
  end subroutine initTemperature


  function shouldRunTemperature() result(shouldRun)
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
    logical :: shouldRun 
    !# result variable for store if method should run

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunTemperature


  function getOutFileName() result(outFileName)
    !# Gets temperature output file name
    !# ---
    !# @info
    !# **Brief:** Gets temperature output file name. </br>
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
    allocate(temperatureIn(var%xDim, var%yDim, 12))
    allocate(temperatureInput(var%xDim, var%yDim, 12))
    allocate(temperatureOut(varCommon%xMax, varCommon%yMax, 12))
    allocate(temperatureOutput(varCommon%xMax, varCommon%yMax, 12))
  end subroutine allocateData


  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Dellocates matrixes. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    deallocate(temperatureIn)
    deallocate(temperatureInput)
    deallocate(temperatureOut)
    deallocate(temperatureOutput)
  end subroutine deallocateData


  function readTemperature() result(isReadOk)
    !# Reads Delta Temp Coldest
    !# ---
    !# @info
    !# **Brief:** Read Delta Temp Coldest, input file. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none

    logical :: isReadOk 
    !# return value (if function was executed sucessfully)
    integer :: inRecSize 
    !# size of file record
    character (len = maxPathLength) :: inFilename 
    !# input filename
    integer :: inFileUnit

    isReadOk = .false.
    inquire (iolength = inRecSize) temperatureIn

    inFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return

    read(unit = inFileUnit, rec = 1) temperatureIn
    temperatureInput = real(temperatureIn, p_r8)
    close(unit = inFileUnit)
    isReadOk = .true.
  end function readTemperature


  function writeTemperature(outFilename) result(isWriteOk)
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
    character(len = *), intent(in) :: outFilename
    integer :: lRecOut  
    !# size of file record
    integer :: datFileUnit

    isWriteOk = .false.
    inquire(iolength = lRecOut) temperatureOut
    datFileUnit = openFile(trim(outFilename), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(datFileUnit < 0) return

    write(unit = datFileUnit, rec = 1) temperatureOut
    close(unit = datFileUnit)
    isWriteOk = .true.

  end function writeTemperature


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
    logical :: isGradsOk 
    !# return value (if function was executed sucessfully)
    character (len = maxPathLength) :: ctlPathFileName 
    !# ctl Grads File name with path
    character (len = maxPathLength) :: gradsBaseName 
    !# dat GrADS base name without extension
    integer :: gradsFileUnit 
    !# File Unit for all grads files

    ! dat file is the same of getOutputFile; not generated here
    isGradsOk = .false.
    gradsBaseName = trim(var%varNameOut) // nLats
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.ctl'

    if(.not. writeTemperature(trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.dat')) return

    ! Write .ctl file -------------------------------------------------------------
    gradsFileUnit = openFile(trim(ctlPathFileName), 'formatted', 'sequential', -1, 'write', 'replace')
    if(gradsFileUnit < 0) return

    write(unit = gradsFileUnit, fmt = '(a)') 'DSET ^' // trim(gradsBaseName) // '.dat'
    write(unit = gradsFileUnit, fmt = '(a)') '*'
    write(unit = gradsFileUnit, fmt = '(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit = gradsFileUnit, fmt = '(a)') '*'
    write(unit = gradsFileUnit, fmt = '(a)') 'UNDEF -999.0'
    write(unit = gradsFileUnit, fmt = '(a)') '*'
    write(unit = gradsFileUnit, fmt = '(a)') 'TITLE Temperature on a Gaussian Grid'
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
    write(unit = gradsFileUnit, fmt = '(a)') 'deltat 0 99 Temperature [C]'
    write(unit = gradsFileUnit, fmt = '(a)') 'ENDVARS'
    close(unit = gradsFileUnit)
    isGradsOk = .true.

  end function generateGrads


  function generateTemperature() result(isExecOk)
    !# Generates Delta Temp Coldest
    !# ---
    !# @info
    !# **Brief:** Generates Delta Temp Coldest output. This subroutine is the main
    !# method for use this module. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk 
    !# return value (if function was executed sucessfully)
    integer :: k
    logical :: flagInput(5), flagOutput(5)

    isExecOk = .false.

    if(var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_Lat0, p_Lon0)
    else
      call initAreaInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if
    call allocateData()
    if(.not.readTemperature()) return

    ! Interpolate Input Regular Grid Deep Soil Temperature To Gaussian Grid on Output
    do k = 1, 12
      if(var%linear) then
        call doLinearInterpolation(temperatureInput(:, :, k), temperatureOutput(:, :, k))
      else
        call doAreaInterpolation(temperatureInput(:, :, k), temperatureOutput(:, :, k))
      end if
      temperatureOut(:, :, k) = real(temperatureOutput(:, :, k), p_r4)
    end do

    if(.not. writeTemperature(getOutFileName())) return

    if(varCommon%grads) then
      if(.not. generateGrads()) return
    end if

    call deallocateData()
    isExecOk = .true.

  end function generateTemperature


end module Mod_Temperature

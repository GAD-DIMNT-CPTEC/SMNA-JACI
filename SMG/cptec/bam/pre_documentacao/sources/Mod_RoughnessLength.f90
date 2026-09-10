!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_RoughnessLength </br></br>
!#
!# **Brief**: Module responsible for reading the map with the global roughness
!# distribution and interpolating for the resolution of the model</br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/dataout/RoughnessLengthClima.dat
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; model/datain/RoughnessLength.GZZZZZ (Ex.: pre/dataout/LandSeaMask.G00450)</br>
!# &bull; pre/dataout/RoughnessLength.GZZZZZ.ctl
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
!#  <li>27-01-2020 - Eduardo Khamis    - version: 2.1.0 </li>
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

module Mod_RoughnessLength

  use Mod_LinearInterpolation, only: gLatsL=>latOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation,   only: gLatsA=>gLats,  initAreaInterpolation,   doAreaInterpolation
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut
  
  implicit none
  
  public :: generateRoughnessLength
  public :: getNameRoughnessLength
  public :: initRoughnessLength
  public :: shouldRunRoughnessLength

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'
  
  !parameters
  real(kind=p_r8), parameter :: p_Lat0 = 90.0_p_r8    
  !# Start at North Pole
  real(kind=p_r8), parameter :: p_Lon0 = 0.0_p_r8     
  !# Start at Prime Meridian
  
  !input variables
  type RoughnessLengthNameListData
    integer :: xDim   
    !# Number of Longitudes 
    integer :: yDim   
    !# Number of Latitudes 
    logical :: linear 
    !# flag for interpolation type (true = linear, false = area)
    character(len=maxPathLength) :: varName='RoughnessLengthClima'  
    !# Input name file
    character(len=maxPathLength) :: varNameOut='RoughnessLength'    
    !# Output name file
  end type RoughnessLengthNameListData

  type(varCommonNameListData)       :: varCommon
  type(RoughnessLengthNameListData) :: var
  namelist /RoughnessLengthNameList/   var  

  !set values are default
  character(len=7)   :: nLats='.G     '                 
  !# posfix of roughness Length file
  
  !internal variables
  real(kind=p_r4), dimension(:,:), allocatable :: roughLengthIn
  real(kind=p_r8), dimension(:,:), allocatable :: roughLengthInput
  real(kind=p_r4), dimension(:,:), allocatable :: roughLengthOut
  real(kind=p_r8), dimension(:,:), allocatable :: roughLengthOutput
  logical :: flagInput(5), flagOutput(5)
  
  character(len=*), parameter :: header = 'Roughness Length          | '


  contains

  
  function getNameRoughnessLength() result(returnModuleName)
    !# Returns RoughnessLength Module Name
    !# ---
    !# @info
    !# **Brief:** Returns RoughnessLength Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "RoughnessLength"
  end function getNameRoughnessLength


  subroutine initRoughnessLength(nameListFileUnit, varCommon_)
    !# Initialization of RoughnessLength module
    !# ---
    !# @info
    !# **Brief:** Initialization of RoughnessLength module. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = RoughnessLengthNameList)
    varCommon = varCommon_
    write (nLats(3:7), '(I5.5)') varCommon%yMax

    flagInput(1)=.true.   ! Start at North Pole
    flagInput(2)=.true.   ! Start at Prime Meridian
    flagInput(3)=.true.   ! Latitudes Are at Center of Box
    flagInput(4)=.true.   ! Longitudes Are at Center of Box
    flagInput(5)=.false.  ! Regular Grid
    flagOutput(1)=.true.  ! Start at North Pole
    flagOutput(2)=.true.  ! Start at Prime Meridian
    flagOutput(3)=.false. ! Latitudes Are at North Edge of Box
    flagOutput(4)=.true.  ! Longitudes Are at Center of Box
    flagOutput(5)=.true.  ! Gaussian Grid

  end subroutine initRoughnessLength


  function shouldRunRoughnessLength() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunRoughnessLength


  function getOutFileName() result(roughnessLengthOutFilename)
    !# Gets RoughnessLength Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets RoughnessLength Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: roughnessLengthOutFilename

    roughnessLengthOutFilename = trim(varCommon%dirModelIn) // trim(var%varNameOut) // nLats
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
    allocate(roughLengthIn(var%xDim,var%yDim))
    allocate(roughLengthInput(var%xDim,var%yDim))
    allocate(roughLengthOut(varCommon%xMax,varCommon%yMax))
    allocate(roughLengthOutput(varCommon%xMax,varCommon%yMax))    
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
    deallocate(roughLengthIn)
    deallocate(roughLengthInput)
    deallocate(roughLengthOut)
    deallocate(roughLengthOutput)    
  end subroutine deallocateData
  
  
  function readRoughLength() result(isReadOk)
    !# Reads In Input Roughness Length
    !# ---
    !# @info
    !# **Brief:** Reads Roughness Length, input file. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReadOk 
    character(len = maxPathLength) :: inFilename 
    integer :: lRecIn 
    integer :: inFileUnit 
 
    isReadOk = .false.   
    inquire(iolength=lRecIn) roughLengthIn
    inFilename = trim(varCommon%dirPreOut)//trim(var%varName)//'.dat'
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', lRecIn, 'read', 'old')
    if (inFileUnit < 0) return
    read(unit=inFileUnit, rec=1) roughLengthIn
    roughLengthInput = real(roughLengthIn, p_r8)
    close(unit=inFileUnit)
    isReadOk = .true.

  end function readRoughLength
 
 
  function writeRoughLength(outputFile) result(isWriteOk)
    !# Writes Rough Length
    !# ---
    !# @info
    !# **Brief:** Generates ZoRL Input File for the BAM Model. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    
    implicit none
    logical :: isWriteOk
    character(len=*), intent(in) :: outputFile 
    !# complete path to write file
    integer :: lRecOut
    integer :: outFileUnit

    isWriteOk = .false.    
    inquire(iolength=lRecOut) roughLengthOut
    outFileUnit = openFile(outputFile, 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if (outFileUnit < 0) return
    write(unit=outFileUnit, rec=1) roughLengthOut
    close(unit=outFileUnit)
    isWriteOk = .true.

  end function writeRoughLength

  
  function generateGrads() result(isGradsOk)
    !# Generates Grads files 
    !# ---
    !# @info
    !# **Brief:** Generates .ctl and .dat files to check output. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isGradsOk
    integer :: gradsFileUnit
 
    isGradsOk = .false.
    ! Save the same file in pre/dataout directory to view in grads  
    !if (.not. writeRoughLength(trim(varCommon%dirPreOut)//trim(var%varNameOut)//nLats)) return
    
    ! Write GrADS Control File --------------------------------------------------------------
    gradsFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varNameOut)//nLats//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if (gradsFileUnit < 0) return    

    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'//trim(var%varNameOut)//nLats
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'UNDEF -999.0'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE Roughness Length on a Gaussian Grid'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',varCommon%xMax,' LINEAR ', &
                            0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write(unit=gradsFileUnit, fmt='(a,i5,a)') 'YDEF ',varCommon%yMax,' LEVELS '
    if(var%linear) then
      write(unit=gradsFileUnit, fmt='(8f10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write(unit=gradsFileUnit, fmt='(8f10.5)') gLatsA(varCommon%yMax:1:-1)
    end if
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF 1 LEVELS 1000'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 1 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS 1'
    write(unit=gradsFileUnit, fmt='(a)') 'RGHL 0 99 Roughness Length [cm]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)
    ! ---------------------------------------------------------------------------------------
    isGradsOk = .true.

  end function generateGrads
  

  function generateRoughnessLength() result(isExecOk)
    !# Generates Rough Length 
    !# ---
    !# @info
    !# **Brief:** Generates Rough Length output. This subroutine is the main method
    !# for use this module. Only file name of namelist is needed to use it. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# &bull; Eduardo Khamis - changing subroutine to function </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk

    isExecOk = .false.    
    if(var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_Lat0, p_Lon0)
    else
      call initAreaInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if
    call allocateData()
    if (.not. readRoughLength()) then
      call msgWarningOut(header, "Error reading RoughnessLengthClima.dat file")
      return
    end if
    
    ! Interpolate Input Regular Rough Length To Gaussian Grid on Output
    if(var%linear) then
      call doLinearInterpolation(roughLengthInput, roughLengthOutput)
    else
      call doAreaInterpolation(roughLengthInput, roughLengthOutput)
    end if
    roughLengthOut = real(roughLengthOutput, p_r4)
    
    if (.not. writeRoughLength(trim(varCommon%dirModelIn)//trim(var%varNameOut)//nLats)) then
      call msgWarningOut(header, "Error writing RoughnessLength file")
      return
    end if
    
    if(varCommon%grads) then
      if (.not. generateGrads()) then
        call msgWarningOut(header, "Error while generating grads files")
        return
      end if
    end if
    
    call deallocateData()
    isExecOk = .true.

  end function generateRoughnessLength
  
  
end module Mod_RoughnessLength

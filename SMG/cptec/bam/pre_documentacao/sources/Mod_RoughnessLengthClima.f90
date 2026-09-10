!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_RoughnessLengthClima </br></br>
!#
!# **Brief**: Module responsible for reading the data set described in the 
!# zorlng.form file and generating the map with the global roughness distribution </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/databcs/zorlng.form
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; pre/dataout/RoughnessLengthClima.dat </br>
!# &bull; pre/dataout/RoughnessLengthClima.ctl
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

module Mod_RoughnessLengthClima

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut

  implicit none
  
  public :: generateRoughnessLengthClima
  public :: getNameRoughnessLengthClima
  public :: initRoughnessLengthClima
  public :: shouldRunRoughnessLengthClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !input variables ---------------------------------------------------------------------------
  type RoughnessLengthClimaNameListData
    integer :: xDim  
    !# Number of Longitudes in SSiB Vegetation Mask Data
    integer :: yDim  
    !# Number of Latitudes in SSiB Vegetation Mask Data
    character(len=maxPathLength) :: varName='RoughnessLengthClima'    
    !# output prefix file name
    character(len=maxPathLength) :: fileBCs='zorlng.form'             
    !# input file name
  end type RoughnessLengthClimaNameListData

  type(varCommonNameListData)            :: varCommon
  type(RoughnessLengthClimaNameListData) :: var
  namelist /RoughnessLengthClimaNameList/   var  
  
  !internal variables
  real(kind=p_r4), dimension (:,:), allocatable :: roughLength
  
  ! aux variables
  integer :: xDim1

  character(len=*), parameter :: header = 'Roughness Length Clima          | '
  
  
  contains

  
  function getNameRoughnessLengthClima() result(returnModuleName)
    !# Returns RoughnessLengthClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns RoughnessLengthClima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "RoughnessLengthClima"
  end function getNameRoughnessLengthClima


  subroutine initRoughnessLengthClima(nameListFileUnit, varCommon_)
    !# Initialization of RoughnessLengthClima module
    !# ---
    !# @info
    !# **Brief:** Initialization of RoughnessLengthClima module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = RoughnessLengthClimaNameList)
    varCommon = varCommon_

    xDim1 = var%xDim - 1

  end subroutine initRoughnessLengthClima


  function shouldRunRoughnessLengthClima() result(shouldRun)
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
  end function shouldRunRoughnessLengthClima


  function getOutFileName() result(roughnessLengthClimaOutFilename)
    !# Gets RoughnessLengthClima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets RoughnessLengthClima Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: roughnessLengthClimaOutFilename

    roughnessLengthClimaOutFilename = trim(varCommon%dirModelIn) // trim(var%varName) // '.dat'
  end function getOutFileName


  subroutine allocateMatrixes()
    !# Allocates matrix 
    !# ---
    !# @info
    !# **Brief:** Allocates matrix roughLength. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinf
    allocate(roughLength(var%xDim, var%yDim))
  end subroutine allocateMatrixes


  subroutine deallocateMatrixes()
    !# Deallocates matrix
    !# ---
    !# @info
    !# **Brief:** Deallocates matrix roughLength. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinf
    deallocate(roughLength)
  end subroutine deallocateMatrixes


  function readRoughLength() result(isReadOk)
    !# Reads sib mask file
    !# ---
    !# @info
    !# **Brief:** Reads tgdeep file and load in roughLength. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinf
    implicit none
    logical :: isReadOk
    character(len = maxPathLength) :: inFilename 
    integer :: inFileUnit 

    isReadOk = .false.
    inFilename = trim(varCommon%dirBCs)//trim(var%fileBCs)
    inFileUnit = openFile(trim(inFilename), 'formatted', 'sequential', -1, 'read', 'old')
    if (inFileUnit < 0) return
    read(unit=inFileUnit, fmt='(5e15.8)') roughLength
    close(unit=inFileUnit)
    isReadOk = .true.

  end function readRoughLength

  
  function writeRoughLength(outputFile) result(isWriteOk)
    !# Writes Deep Soil Temp
    !# ---
    !# @info
    !# **Brief:** Writes a file with deep soil temp. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinf
    
    
    implicit none
    logical :: isWriteOk
    character(len=*), intent(in) :: outputFile 
    !# complete path to write file
    integer :: lRec
    integer :: outFileUnit

    isWriteOk = .false.
    inquire(iolength=lRec) roughLength(1:xDim1,:)
    outFileUnit = openFile(outputFile, 'unformatted', 'direct', lRec, 'write', 'replace')
    if (outFileUnit < 0) return
    write(unit=outFileUnit, rec=1) roughLength(1:xDim1,:)
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
    !# @endinf
    implicit none
    logical :: isGradsOk
    integer :: gradsFileUnit
 
    isGradsOk = .false.  
    ! Save the same file in pre/dataout directory to view in grads
    !if (.not. writeRoughLength(trim(varCommon%dirPreOut)//trim(var%varName)//'.dat')) return
    
    ! create .ctl file ------------------------------------------------------------------
    gradsFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varName)//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if (gradsFileUnit < 0) return

    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'//trim(var%varName)//'.dat'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'UNDEF -999.0'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE CLimatological Roughness Length'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',xDim1,' LINEAR ', &
                            0.0_p_r4, 360.0_p_r4 / real(xDim1, p_r4)
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'YDEF ',var%yDim,' LINEAR ',  &
                            -90.0_p_r4, 180.0_p_r4 / real(var%yDim - 1, p_r4)
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF 1 LEVELS 1000'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 1 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS 1'
    write(unit=gradsFileUnit, fmt='(a)') 'RGHL 0 99 Roughness Length [cm]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)
    isGradsOk = .true.

  end function generateGrads

  
  function generateRoughnessLengthClima() result(isExecOk)
    !# Generates CLimatological Roughness Length 
    !# ---
    !# @info
    !# **Brief:** Generates Rough Length Clima output. This subroutine is the main
    !# method for use this module. Only file name of namelist is needed to use it. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# &bull; Eduardo Khamis   - changing subroutine to function </br>
    !# **Date**: jan/2020 </br>
    !# @endinf
    implicit none
    logical :: isExecOk

    isExecOk = .false.    
    call allocateMatrixes()
    if (.not. readRoughLength()) then
      call msgWarningOut(header, "Error reading zorlng.form file")
      return
    end if
    if (.not. writeRoughLength(trim(varCommon%dirPreOut)//trim(var%varName)//'.dat')) then
      call msgWarningOut(header, "Error writing RoughnessLengthClima file")
      return
    end if
    
    ! Generate grads output for debug or test
    if(varCommon%grads) then
      if (.not. generateGrads()) then
        call msgWarningOut(header, "Error while generating grads files")
        return
      end if
    end if
    
    call deallocateMatrixes()
    isExecOk = .true.    
  end function generateRoughnessLengthClima
  
 
end module Mod_RoughnessLengthClima

!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_DeepSoilTemperatureClima </br></br>
!#
!# **Brief**: Module responsible for generating a global map with deep soil temperature </br>
!#
!# First Point of Input and Output Data is at North Pole and Greenwhich. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/databcs/tgdeep.form </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/DeepSoilTemperatureClima.dat </br>
!# &bull; pre/dataout/DeepSoilTemperatureClima.ctl
!# </br></br>
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
!#  <li>01-04-2018 - Daniel Lamosa    - version: 2.0.0 </li>
!#  <li>27-01-2020 - Eduardo Khamis   - version: 2.1.0 </li>
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


module Mod_DeepSoilTemperatureClima

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut

  implicit none
  
  public :: generateDeepSoilTemperatureClima
  public :: getNameDeepSoilTemperatureClima
  public :: initDeepSoilTemperatureClima
  public :: shouldRunDeepSoilTemperatureClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !input variables

  type DeepSoilTemperatureClimaNameListData
    integer :: xDim  
    !# Number of Longitudes in SSiB Vegetation Mask Data
    integer :: yDim  
    !# Number of Latitudes in SSiB Vegetation Mask Data
    character(len=maxPathLength) :: varName='DeepSoilTemperatureClima'    
    !# output prefix file name
    character(len=maxPathLength) :: fileBCs='tgdeep.form'
    !# ?                
  end type DeepSoilTemperatureClimaNameListData

  type(varCommonNameListData)                :: varCommon
  type(DeepSoilTemperatureClimaNameListData) :: var
  namelist /DeepSoilTemperatureClimaNameList/   var  
  
  ! internal variables
  real(kind=p_r4), dimension (:,:), allocatable :: deepSoilTemp
    
  character(len=*), parameter :: header = 'Deep Soil Temperature Clima          | '


  contains

  
  function getNameDeepSoilTemperatureClima() result(returnModuleName)
    !# Returns DeepSoilTemperatureClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns DeepSoilTemperatureClima Module Name </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "DeepSoilTemperatureClima"
  end function getNameDeepSoilTemperatureClima

  subroutine initDeepSoilTemperatureClima(nameListFileUnit, varCommon_)
    !# Initializes the DeepSoilTemperatureClima module
    !# ---
    !# @info
    !# **Brief:** Initializes the DeepSoilTemperatureClima module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo  
    implicit none 
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = DeepSoilTemperatureClimaNameList)
    varCommon = varCommon_
  end subroutine initDeepSoilTemperatureClima

  function shouldRunDeepSoilTemperatureClima() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it
    !# does not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunDeepSoilTemperatureClima


  function getOutFileName() result(deepSoilTemperatureClimaOutFilename)
    !# Gets DeepSoilTemperatureClima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets DeepSoilTemperatureClima Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: deepSoilTemperatureClimaOutFilename

    deepSoilTemperatureClimaOutFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
  end function getOutFileName


   subroutine allocateMatrixes()
    !# Allocates matrix
    !# ---
    !# @info
    !# **Brief:** Allocates matrix DeepSoilTemp. </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo
    allocate(DeepSoilTemp(var%xDim, var%yDim))
  end subroutine allocateMatrixes


  subroutine deallocateMatrixes()
    !# Deallocates matrix
    !# ---
    !# @info
    !# **Brief:** Deallocates matrix DeepSoilTemp. </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo
    deallocate(DeepSoilTemp)
  end subroutine deallocateMatrixes
  
  
  function readDeepSoilTemp() result(isReadOk)
    !# Reads sib mask file
    !# ---
    !# @info
    !# **Brief:** Reads tgdeep file and load in DeepSoilTemp. </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo
    implicit none
    logical :: isReadOk 
    character(len = maxPathLength) :: inFilename 
    integer :: inFileUnit 
 
    isReadOk = .false.   
    inFilename = trim(varCommon%dirBCs)//trim(var%fileBCs)
    inFileUnit = openFile(trim(inFilename), 'formatted', 'sequential', -1, 'read', 'old')
    if (inFileUnit < 0) return
    read(unit=inFileUnit, fmt='(5e15.8)') deepSoilTemp
    close(unit=inFileUnit)
    isReadOk = .true.

  end function readDeepSoilTemp


  function writeDeepSoilTemp(outputFile) result(isWriteOk)    
    !# Writes Deep Soil Temp
    !# ---
    !# @info
    !# **Brief:** Writes a file with deep soil temp. </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo
    implicit none
    logical :: isWriteOk 
    character(len=*), intent(in) :: outputFile 
    !# complete path to write file
    integer :: lRec
    integer :: outFileUnit
 
    isWriteOk = .false.   
    inquire(iolength=lRec) deepSoilTemp
    outFileUnit = openFile(outputFile, 'unformatted', 'direct', lRec, 'write', 'replace')
    if (outFileUnit < 0) return
    write(unit=outFileUnit, rec=1) deepSoilTemp
    close(unit=outFileUnit)
    isWriteOk = .true.

  end function writeDeepSoilTemp
  

  function generateGrads() result(isGradsOk)
    !# Generates Grads files
    !# ---
    !# @info
    !# **Brief:** Generates .ctl and .dat files to check output </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo
    implicit none
    logical :: isGradsOk
    integer :: gradsFileUnit
 
    isGradsOk = .false.    
    gradsFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varName)//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    ! create .ctl file
    if (gradsFileUnit < 0) return    

    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'//trim(var%varName)//'.dat'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'UNDEF -999.0'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE CLimatological Deep Soil Temperature'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',var%xDim,' LINEAR ', &
                            0.0_p_r4, 360.0_p_r4 / real(var%xDim, p_r4)
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'YDEF ',var%yDim,' LINEAR ', &
                            -90.0_p_r4, 180.0_p_r4 / real(var%yDim - 1, p_r4)
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF 1 LEVELS 1000'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 1 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS 1'
    write(unit=gradsFileUnit, fmt='(a)') 'DSTP 0 99 Deep Soil Temperature [K]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)
    ! -----------------------------------------------------------------------------------
    isGradsOk = .true.

  end function generateGrads


  function generateDeepSoilTemperatureClima() result(isExecOk)
    !# Generates Deep Soil Temperature Clima
    !# ---
    !# @info
    !# **Brief:** Generates Deep Soil Temperature Clima output. This subroutine
    !# is the main method for use this module. Only file name of namelist is
    !# needed to use it. </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# &bull; Eduardo Khamis   - changing subroutine to function </br>
    !# **Date**: jan/2020 <br>
    !# @endinfo
    implicit none
    logical :: isExecOk

    isExecOk = .false.    
    call allocateMatrixes()
    if (.not. readDeepSoilTemp()) then
      call msgWarningOut(header, "Error reading tgdeep.form file")
      return
    end if
    if (.not. writeDeepSoilTemp(trim(varCommon%dirPreOut)//trim(var%varName)//'.dat')) then
      call msgWarningOut(header, "Error writing DeepSoilTemperatureClima file")
      return
    end if
    ! Generate grads output for debug or test
    if (varCommon%grads .and. .not. generateGrads()) then
      call msgWarningOut(header, "Error while generating grads files")
      return
    end if
    call deallocateMatrixes()
    isExecOk = .true.

  end function generateDeepSoilTemperatureClima
  
  
end module Mod_DeepSoilTemperatureClima

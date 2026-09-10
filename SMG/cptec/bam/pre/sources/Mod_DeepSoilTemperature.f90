!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_DeepSoilTemperature </br></br>
!#
!# **Brief**: Module responsible for interpolating the deep soil temperature 
!# field for the model grid </br>
!#
!# Linear (.TRUE.) or area weighted (.FALSE.) interpolation. Interpolate regular
!# to gaussian. Regular input data is assumed to be oriented with the north pole
!# and greenwich as the first point. Gaussian output data is interpolated to be
!# oriented with the north pole and greenwich as the first point input for the 
!# AGCM is assumed IEEE 32 bits big endian. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/DeepSoilTemperatureClima.dat </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/DeepSoilTemperature.G00450 </br>
!# &bull; pre/dataout/DeepSoilTemperature.G00450.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita          - version: 1.1.1 </li>
!#  <li>21-03-2019 - Denis Eiras     - version: 2.0.0 </li>
!#  <li>20-05-2019 - Eduardo Khamis  - version: 2.0.1 </li>
!#  <li>27-01-2020 - Eduardo Khamis  - version: 2.1.0 </li>
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


module Mod_DeepSoilTemperature

  use Mod_LinearInterpolation, only: gLatsL=>LatOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation,   only: gLatsA=>gLats,  initAreaInterpolation,   doAreaInterpolation
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut
  
  implicit none
  
  public :: generateDeepSoilTemperature
  public :: getNameDeepSoilTemperature
  public :: initDeepSoilTemperature
  public :: shouldRunDeepSoilTemperature
  
  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  ! set values are default
  character(len=7)   :: nLats='.G     '                   
  !# posfix land sea mask file
  real(kind=p_r8), parameter :: p_Lat0 = 90.0_p_r8 
  !# Start at North Pole        check: in some modules is = 89.5
  real(kind=p_r8), parameter :: p_Lon0 = 0.0_p_r8  
  !# Start at Prime Meridian    check: in some modules is = 0.5
  
  ! input variables
  type DeepSoilTemperatureNameListData
    integer :: xDim   
    !# Number of Longitudes 
    integer :: yDim   
    !# Number of Latitudes 
    logical :: linear 
    !# flag for interpolation type (true = linear, false = area) 
    character(len=maxPathLength) :: varName='DeepSoilTempClima'       
    character(len=maxPathLength) :: varNameOut='DeepSoilTemperature'  
  end type DeepSoilTemperatureNameListData

  type(varCommonNameListData)           :: varCommon
  type(DeepSoilTemperatureNameListData) :: var
  namelist /DeepSoilTemperatureNameList/   var  

  ! internal variables
  real(kind=p_r4), dimension(:,:), allocatable :: deepSoilTempIn
  real(kind=p_r8), dimension(:,:), allocatable :: deepSoilTempInput
  real(kind=p_r4), dimension(:,:), allocatable :: deepSoilTempOut
  real(kind=p_r8), dimension(:,:), allocatable :: deepSoilTempOutput
  logical :: flagInput(5), flagOutput(5)

  ! temporary variables
  integer :: ios 
  !# io status
  
  character(len=*), parameter :: header = 'Deep Soil Temperature          | '


  contains
  
  
  function getNameDeepSoilTemperature() result(returnModuleName)
    !# Returns DeepSoilTemperature Module Name
    !# ---
    !# @info
    !# **Brief:** Returns DeepSoilTemperature Module Name. </br>
    !# **Author**: Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "DeepSoilTemperature"
  end function getNameDeepSoilTemperature


  subroutine initDeepSoilTemperature(nameListFileUnit, varCommon_)
    !# Initializes the DeepSoilTemperature module
    !# ---
    !# @info
    !# **Brief:** Initializes the DeepSoilTemperature module, defined in PRE_run.nml. </br>
    !# **Author**: Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = DeepSoilTemperatureNameList)
    varCommon = varCommon_
    write (nLats(3:7), '(I5.5)') varCommon%yMax

    flagInput(1)=.true.   
    ! Start at North Pole
    flagInput(2)=.true.   
    ! Start at Prime Meridian
    flagInput(3)=.true.   
    ! Latitudes Are at Center of Box
    flagInput(4)=.true.   
    ! Longitudes Are at Center of Box
    flagInput(5)=.false.  
    ! Regular Grid
    flagOutput(1)=.true.  
    ! Start at North Pole
    flagOutput(2)=.true.  
    ! Start at Prime Meridian
    flagOutput(3)=.false. 
    ! Latitudes Are at North Edge of Box
    flagOutput(4)=.true.  
    ! Longitudes Are at Center of Box
    flagOutput(5)=.true.  
    ! Gaussian Grid

  end subroutine initDeepSoilTemperature


  function shouldRunDeepSoilTemperature() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Author**: Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunDeepSoilTemperature


  function getOutFileName() result(deepSoilTemperatureOutFilename)
    !# Gets DeepSoilTemperature Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets DeepSoilTemperature Out Filename. </br>
    !# **Author**: Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: deepSoilTemperatureOutFilename

    deepSoilTemperatureOutFilename = trim(varCommon%dirModelIn) // trim(var%varNameOut) // nLats
  end function getOutFileName


  subroutine allocateData()
    !# Allocates data
    !# ---
    !# @info
    !# **Brief:** Allocates matrixes based in namelist. </br>
    !# **Author**: Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    allocate(deepSoilTempIn(var%xDim,var%yDim))
    allocate(deepSoilTempInput(var%xDim,var%yDim))
    allocate(deepSoilTempOut(varCommon%xMax,varCommon%yMax))
    allocate(deepSoilTempOutput(varCommon%xMax,varCommon%yMax))    
  end subroutine allocateData

  
  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Dellocates matrixes based in namelist. </br>
    !# **Author**: Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    deallocate(deepSoilTempIn)
    deallocate(deepSoilTempInput)
    deallocate(deepSoilTempOut)
    deallocate(deepSoilTempOutput)    
  end subroutine deallocateData

  
  function readDeepSoilTempClima() result(isReadOk)
    !# Reads Deep Soil Temperature Clima
    !# ---
    !# @info
    !# **Brief:** Reads Deep Soil Temperature Clima, input file. </br>
    !# **Author**: Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReadOk 
    character(len = maxPathLength) :: inFilename 
    integer :: lRecIn 
    integer :: inFileUnit 
 
    isReadOk = .false.   
    inquire(iolength=lRecIn) deepSoilTempIn
    inFilename = trim(varCommon%dirPreOut)//trim(var%varName)//'.dat'
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', lRecIn, 'read', 'old')
    if (inFileUnit < 0) return
    read(unit=inFileUnit, rec=1) deepSoilTempIn
    deepSoilTempInput = real(deepSoilTempIn, p_r8)
    close(unit=inFileUnit)
    isReadOk = .true.

  end function readDeepSoilTempClima
  
  
  function writeDeepSoilTemperature(outputFile) result(isWriteOk)
    !# Writes Deep Soil Temperature
    !# ---
    !# @info
    !# **Brief:** Generates Tg1, Tg2, Tg3 input file for the AGCM Model at this
    !# climatological input it is assumed that Tg1=Tg2=Tg3. </br>
    !# **Author**: Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isWriteOk
    character(len=*), intent(in) :: outputFile 
    ! complete path to write file
    integer :: lRecOut
    integer :: outFileUnit
    
    isWriteOk = .false.
    inquire(iolength=lRecOut) deepSoilTempOut
    outFileUnit = openFile(outputFile, 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if (outFileUnit < 0) return
    write(unit=outFileUnit, rec=1) deepSoilTempOut
    write(unit=outFileUnit, rec=2) deepSoilTempOut
    write(unit=outFileUnit, rec=3) deepSoilTempOut
    close(unit=outFileUnit)
    isWriteOk = .true.

  end function writeDeepSoilTemperature  


  function generateGrads() result(isGradsOk)
    !# Generates Grads files
    !# ---
    !# @infos
    !# **Brief:** Generates .ctl and .dat files to check output. </br>
    !# **Author**: Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isGradsOk
    integer :: gradsFileUnit

    isGradsOk = .false.
    
    ! Save the same file in pre/dataout directory to view in grads
    !if (.not. writeDeepSoilTemperature(trim(varCommon%dirPreOut)//trim(var%varNameOut)//nLats)) then
    !  call msgWarningOut(header, "Error writing DeepSoilTemperature file")
    !  return
    !end if
    
    ! Write GrADS control file ---------------------------------------------
    gradsFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varNameOut)//nLats//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if (gradsFileUnit < 0) return    
    
    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'//trim(var%varNameOut)//nLats
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,1pg12.5)') 'UNDEF ', p_undef
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE Soil Temperature on a Gaussian Grid'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') &
          'XDEF ',varCommon%xMax,' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write(unit=gradsFileUnit, fmt='(a,i5,a)') 'YDEF ',varCommon%yMax,' LEVELS '
    if(var%linear) then
      write(unit=gradsFileUnit, fmt='(8f10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write(unit=gradsFileUnit, fmt='(8f10.5)') gLatsA(varCommon%yMax:1:-1)
    end if
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF 1 LEVELS 1000'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 1 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS 3'
    write(unit=gradsFileUnit, fmt='(a)') 'DST1 0 99 Deep Soil Temperature [K]'
    write(unit=gradsFileUnit, fmt='(a)') 'DST2 0 99 Deep Soil Temperature [K]'
    write(unit=gradsFileUnit, fmt='(a)') 'DST3 0 99 Deep Soil Temperature [K]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)

    isGradsOk = .true.

  end function generateGrads
  

  function generateDeepSoilTemperature() result(isExecOk)
    !# Generates Deep Soil Temperature
    !# ---
    !# @info
    !# **Brief:** Generates Deep Soil Temperature output. This subroutine is the
    !# main method for use this module. Only file name of namelist is needed to use it. </br>
    !# **Author**: Eduardo Khamis - changing subroutine to function </br>
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
    if (.not. readDeepSoilTempClima()) then
      call msgWarningOut(header, "Error reading DeepSoilTemperatureClima.dat file")
      return
    end if
    
    ! Interpolate input regular grid deep soil temperature to gaussian grid on output
    if(var%linear) then
      call doLinearInterpolation(deepSoilTempInput, deepSoilTempOutput)
    else
      call doAreaInterpolation(deepSoilTempInput, deepSoilTempOutput)
    end if
    deepSoilTempOut = real(deepSoilTempOutput, p_r4)
    
    if (.not. writeDeepSoilTemperature(trim(varCommon%dirModelIn)//trim(var%varNameOut)//nLats)) then
      call msgWarningOut(header, "Error writing DeepSoilTemperature file")
      return
    end if
    
    if(varCommon%grads .and. .not. generateGrads()) then
      call msgWarningOut(header, "Error while generating grads files")
      return
    end if
    
    call deallocateData()
    isExecOk = .true.
  end function generateDeepSoilTemperature
  

end module Mod_DeepSoilTemperature

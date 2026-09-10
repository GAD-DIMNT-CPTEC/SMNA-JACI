!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SoilMoisture </br></br>
!#
!# **Brief**: Module responsible for reading the SoilMoistureClima.dat file and
!# the ModelLandSeaMask.GZZZZZ land / sea mask and performs the intersection and
!# interpolation of the fields. In processing, data is generated in the grid and
!# model resolution.. </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/dataout/SoilMoisture Clima.dat </br>
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ (Ex.: pre/dataout/ModelLandSeaMask.G00450)
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/SoilMoisture.GZZZZZ (Ex.: pre/dataout/SoilMoisture.G00450)
!# &bull; pre/dataout/SoilMoisture.GZZZZZ.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <table>
!# <ul type="disc">
!# <tr><td><li>13-11-2004 </li></td> <td> - Jose P. Bonatti </td>                               <td> - version: 1.0.0 <td></tr>
!# <tr><td><li>01-08-2007 </li></td> <td> - Tomita </td>                                        <td> - version: 1.1.1 <td></tr>
!# <tr><td><li>01-04-2018 </li></td> <td> - Daniel M. Lamosa/Barbara A. G. P. Yamada  </td><td> - version: 2.0.0 <td></tr>
!# <tr><td><li>12-03-2020 </li></td> <td> - Eduardo Khamis </td>                                <td> - version: 2.1.0 <td></tr>
!# </ul>
!# </table>
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

module Mod_SoilMoisture

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut
  use Mod_LinearInterpolation, only: gLatsL=>latOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation,   only: gLatsA=>gLats,  initAreaInterpolation,   doAreaInterpolation
  
  implicit none
  
  public :: initSoilMoisture
  public :: generateSoilMoisture
  public :: getnameSoilMoisture
  public :: shouldRunSoilMoisture

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !parameters
  real(kind=p_r8), parameter :: p_seaValue = 150.0_p_r8  
  real(kind=p_r8), parameter :: p_lat0     = 90.0_p_r8   
  !# Start at North Pole
  real(kind=p_r8), parameter :: p_lon0     = 0.0_p_r8    
  !# Start at Prime Meridian
  
  !input variables
  type SoilMoistureNameListData
    integer :: xDim                                                  
    !# Number of Longitudes in Vegetation Mask Data  
    integer :: yDim                                                  
    !# Number of Latitudes  in Vegetation Mask Data
    logical :: linear = .true.                                       
    !# flag for interpolation type (true = linear, false = area) 
    character(len=maxPathLength) :: varName    = 'SoilMoistureClima' 
    character(len=maxPathLength) :: varNameOut = 'SoilMoisture'      
    character(len=maxPathLength) :: nameLSMask = 'ModelLandSeaMask'  
  end type SoilMoistureNameListData

  type(varCommonNameListData)    :: varCommon
  type(SoilMoistureNameListData) :: var
  namelist /SoilMoistureNameList/   var  

  !set values are default
  character(len=7)   :: nLats      = '.G     '           
  !# posfix land sea mask file
  character(len=10)  :: mskfmt     = '(      I1)'        
  !# posfix 
  
  !internal variables
  integer,         dimension(:,:), allocatable :: landSeaMask
  real(kind=p_r4), dimension(:,:), allocatable :: soilMoistureIn
  real(kind=p_r8), dimension(:,:), allocatable :: soilMoistureInput
  real(kind=p_r4), dimension(:,:), allocatable :: soilMoistureOut
  real(kind=p_r8), dimension(:,:), allocatable :: soilMoistureOutput
  !integer,         dimension(:,:), allocatable :: maskInput
  logical                                      :: flagInput(5), flagOutput(5)

  character(len=*), parameter :: header = 'Soil Moisture          | '


  contains


  function getNameSoilMoisture() result(returnModuleName)
    !# Returns SoilMoisture Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SoilMoisture Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo  
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "SoilMoisture"
  end function getNameSoilMoisture


  subroutine initSoilMoisture(nameListFileUnit, varCommon_)
    !# Initialization of SoilMoisture module
    !# ---
    !# @info
    !# **Brief:** Initialization of SoilMoisture module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo 
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
    read(unit = nameListFileUnit, nml = SoilMoistureNameList)
    varCommon = varCommon_
    
    write(mskfmt(2:7), '(i6)') varCommon%xMax
    write(nLats(3:7), '(i5.5)') varCommon%yMax

    ! For Area Weighted Interpolation
    ! allocate(maskInput(var%xDim,var%yDim))
    ! maskInput=1

    flagInput(1)=.true.   ! Start at North Pole
    flagInput(2)=.true.   ! Start at Prime Meridian
    flagInput(3)=.true.   ! Latitudes Are at North Edge
    flagInput(4)=.true.   ! Longitudes Are at Western Edge
    flagInput(5)=.false.  ! Regular Grid
    flagOutput(1)=.true.  ! Start at North Pole
    flagOutput(2)=.true.  ! Start at Prime Meridian
    flagOutput(3)=.false. ! Latitudes Are at North Edge of Box
    flagOutput(4)=.true.  ! Longitudes Are at Center of Box
    flagOutput(5)=.true.  ! Gaussian Grid

  end subroutine initSoilMoisture


  function shouldRunSoilMoisture() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo 
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunSoilMoisture


  function getOutFileName() result(soilMoistureOutFilename)
    !# Gets SoilMoisture Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SoilMoisture Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo 
    implicit none
    character(len = maxPathLength) :: soilMoistureOutFilename

    soilMoistureOutFilename = trim(varCommon%dirModelIn) // trim(var%varNameOut) // nLats
  end function getOutFileName

  
  subroutine allocateData()
    !# Allocates data
    !# ---
    !# @info
    !# **Brief:** Allocates matrixes based in namelist. </br>
    !# **Authors**: </br>
    !# &bull; Jose P. Bonatti </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo 
    implicit none

    allocate(soilMoistureIn(var%xDim, var%yDim))
    allocate(soilMoistureInput(var%xDim, var%yDim))
    allocate(soilMoistureOut(varCommon%xMax, varCommon%yMax))
    allocate(soilMoistureOutput(varCommon%xMax, varCommon%yMax))
    allocate(landSeaMask(varCommon%xMax, varCommon%yMax))
  end subroutine allocateData

  
  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Deallocates matrixes. </br>
    !# **Authors**: </br>
    !# &bull; Jose P. Bonatti </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo 
    implicit none

    deallocate(soilMoistureIn)
    deallocate(soilMoistureInput)
    deallocate(soilMoistureOut)
    deallocate(soilMoistureOutput)
    deallocate(landSeaMask)
  end subroutine deallocateData

  
  function readLandSeaMask() result(isReadOk)
    !# Reads Land Sea Mask
    !# ---
    !# @info
    !# **Brief:** Reads Land Sea Mask, input file. </br>
    !# **Authors**: </br>
    !# &bull; Jose P. Bonatti </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReadOk
    character(len = maxPathLength) :: inFileName 
    integer :: inFileUnit 

    isReadOk = .false.

    inFileName = trim(varCommon%dirPreOut)//trim(var%nameLSMask)//nLats
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'sequential', -1, 'read', 'old')
    if (inFileUnit < 0) return
    read(unit=inFileUnit) landSeaMask
    close(unit=inFileUnit)

    isReadOk = .true.
  end function readLandSeaMask

  
  function generateGrads() result(isGradsOk)
    !# Generates Grads files
    !# ---
    !# @info
    !# **Brief:** Generates .ctl file to check output. This routine differ to 
    !# another, because the generation output file differ too. Ctl file put in
    !# output directory. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isGradsOk
    integer :: ctlFileUnit

    isGradsOk = .false.
    
    ! Write GrADS Control File
    ctlFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varNameOut)//nLats//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if (ctlFileUnit < 0) return

    write(unit=ctlFileUnit, fmt='(a)') 'DSET ^'//trim(var%varNameOut)//nLats
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a,1pg12.5)') 'UNDEF ', p_undef
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a)') 'TITLE Soil Moisture on a Gaussian Grid'
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',varCommon%xMax,' LINEAR ', &
                            0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write(unit=ctlFileUnit, fmt='(a,i5,a)') 'YDEF ',varCommon%yMax,' LEVELS '
    if(var%linear) then
      write(unit=ctlFileUnit, fmt='(8f10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write(unit=ctlFileUnit, fmt='(8f10.5)') gLatsA(varCommon%yMax:1:-1)
    end if
    write(unit=ctlFileUnit, fmt='(a)') 'ZDEF  1 LEVELS 1000'
    write(unit=ctlFileUnit, fmt='(a)') 'TDEF 12 LINEAR JAN2005 1MO'
    write(unit=ctlFileUnit, fmt='(a)') 'VARS  1'
    write(unit=ctlFileUnit, fmt='(a)') 'SOMO  0 99 SoilMoisture [cm]'
    write(unit=ctlFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=ctlFileUnit)

    isGradsOk = .true.
  end function generateGrads
  
  
  function generateSoilMoisture() result(isExecOk)
    !# Generates Soil Moisture
    !# ---
    !# @info
    !# **Brief:** Generates Soil Moisture output. This subroutine is the main 
    !# method for use this module. Only file name of namelist is needed to use it. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk
    
    ! internal variables
    integer :: inFileUnit, outFileUnit
    integer :: lRecIn  
    integer :: lRecOut 
    integer :: month   
    !# counter month
    character (len = maxPathLength) :: soilMoistureInFilename, soilMoistureOutFilename

    isExecOk = .false.
    
    if(var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_Lat0, p_Lon0)
    else
      call initAreaInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if
    call allocateData()
    if (.not. readlandSeaMask()) then
      call msgWarningOut(header, "Error reading LandSeaMask file")
      return
    end if
    
    ! Open soilMoistureIn file to read ------------------------------------------
    inquire(iolength=lRecIn) soilMoistureIn
    soilMoistureInFilename = trim(varCommon%dirPreOut)//trim(var%varName)//'.dat'
    inFileUnit = openFile(trim(soilMoistureInFilename), 'unformatted', 'direct', lRecIn, 'read', 'old')
    if (inFileUnit < 0) return
    ! --------------------------------------------------------------------------
    
    ! Open soilMoistureOut to write --------------------------------------------
    inquire(iolength=lRecOut) soilMoistureOut
    soilMoistureOutFilename = trim(varCommon%dirModelIn)//trim(var%varNameOut)//nLats
    outFileUnit = openFile(trim(soilMoistureOutFilename), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if (outFileUnit < 0) return
    ! --------------------------------------------------------------------------
    
    ! Read In Input SoilMoisture
    do month=1, 12

      read(unit=inFileUnit, rec=month) soilMoistureIn
      soilMoistureInput = real(soilMoistureIn, p_r8)
  
      ! Interpolate Input Regular SoilMoisture To Gaussian Grid Output
      if (var%linear) then
        call doLinearInterpolation(soilMoistureInput, soilMoistureOutput)
      else
        call doAreaInterpolation(soilMoistureInput, soilMoistureOutput)
      end if

      ! Adjust Value Over Sea (150)
      where(landSeaMask == 0) soilMoistureOutput = p_seaValue
      soilMoistureOut = real(soilMoistureOutput, p_r4)
  
      ! Write Out Adjusted Interpolated SoilMoisture
      write (unit=outFileUnit, rec=month) soilMoistureOut
    end do
    
    ! Close files  
    close (unit=inFileUnit)
    close (unit=outFileUnit)
    
    if(varCommon%grads) then
      if (.not. generateGrads()) return
    end if
    
    call deallocateData()
    isExecOk = .true.
  end function generateSoilMoisture
  

end module Mod_SoilMoisture

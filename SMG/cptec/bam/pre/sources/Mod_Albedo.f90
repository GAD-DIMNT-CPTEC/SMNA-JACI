!# @info
!# ---
!# INPE/CPTEC, DIMNT, Division of Numerical Modeling of the Earth System
!# ---
!# </br>
!#
!# **Module**: Mod_Albedo </br></br>
!#
!# **Brief**: Module responsible for generating Albedo files </br></br>
!#
!# This subroutine reads the file AlbedoClima.dat, interpolates to the model
!# resolution and generates two files, Albedo.GZZZZZ and Albedo.GZZZZZ.ctl, which
!# will be used for visualization in GrADS and for the stage of processing of
!# accumulated snow. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/AlbedoClima.dat </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/Albedo.GZZZZZ (Ex.: pre/dataout/Albedo.G00450) </br>
!# &bull; pre/dataout/Albedo.GZZZZZ.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.1 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita          - version: 1.1.1.1 </li>
!#  <li>10-12-2018 - Daniel Lamosa   - version: 2.0.0 - Module creation </li>
!#  <li>25-03-2019 - Denis Eiras     - version: 2.0.1 - Parallel Pre Program </li>
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

module Mod_Albedo

  use Mod_LinearInterpolation, only : gLatsL => latOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, initAreaInterpolation, doAreaInterpolation
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  public :: generateAlbedo
  public :: getNameAlbedo
  public :: initAlbedo
  public :: shouldRunAlbedo

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'


  type AlbedoNameListData
  !# Type for namelist
    integer :: xDim
    !# Number of Longitudes
    integer :: yDim
    !# Number of Latitudes
    integer :: zDim
    logical :: linear
    ! set values are default
    character(len = maxPathLength) :: varName = 'AlbedoClima'
    !# output prefix file name
    character(len = maxPathLength) :: varNameOut = 'Albedo'
    !# Output name file
  end type AlbedoNameListData

  type(varCommonNameListData) :: varCommon
  !# common namelist variable
  type(AlbedoNameListData)    :: var
  !# variable for use in namelist
  namelist /AlbedoNameList/      var

  ! internal variables
  real (kind = p_r4), dimension (:, :), allocatable :: AlbedoIn
  real (kind = p_r4), dimension (:, :), allocatable :: AlbedoOut
  real (kind = p_r8), dimension (:, :), allocatable :: AlbedoInput
  real (kind = p_r8), dimension (:, :), allocatable :: AlbedoOutput
  logical :: flagInput(5)
  !# Input  grid flags
  logical :: flagOutput(5)
  !# Output grid flags

  ! other variables
  character (len = 7) :: nLats = '.G     '
  character(len=*), parameter :: header = 'Albedo                | '


contains

  function getNameAlbedo() result(returnModuleName)
    !# Returns Albedo Module Name
    !# ---
    !# @info
    !# **Brief:** Returns Albedo Module Name. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    !# name of the module

    returnModuleName = "Albedo"
  end function getNameAlbedo


  subroutine initAlbedo(nameListFileUnit, varCommon_)
    !# Initializes Albedo module
    !# ---
    !# @info
    !# **Brief:** Initializes Albedo module. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    !# namelist file Unit comming from Pre_Program
    type(varCommonNameListData), intent(in) :: varCommon_
    !# common namelist variable
    read(unit = nameListFileUnit, nml = AlbedoNameList)
    varCommon = varCommon_
    write (nLats(3:7), '(I5.5)') varCommon%yMax
  end subroutine initAlbedo


  function shouldRunAlbedo() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it
    !# does not generated its out files and was not marked to run. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun
    !# true if module should run (there's output files)

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunAlbedo


  function getOutFileName() result(albedoOutFilename)
    !# Gets Albedo Out file name
    !# ---
    !# @info
    !# **Brief:** Get Albedo Out file name. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: albedoOutFilename
    !# output file name

    write (nLats(3:7), '(I5.5)') varCommon%yMax
    albedoOutFilename = trim(varCommon%dirPreOut) // trim(var%varNameOut) // nLats
  end function getOutFileName


  function generateAlbedo() result(isExecOk)
    !# Generates Albedo
    !# ---
    !# @info
    !# **Brief:** Generates Albedo output </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    implicit none

    logical :: isExecOk
    !# true if execution was sucessfull
    real(kind = p_r8), parameter :: p_Lat0 = 90.0_p_r8
    !# Start at North Pole
    real(kind = p_r8), parameter :: p_Lon0 = 0.0_p_r8
    !# Start at Prime Meridian
    integer :: inFileUnit, outFileUnit, inRecSize, outRecSize, month
    character (len = maxPathLength) :: albedoOutFilename, albedoInFilename

    ! Linear (.TRUE.) or Area Weighted (.FALSE.) Interpolation
    ! Interpolatp_Lat0e Regular To Gaussian
    ! Regular Input Data is Assumed to be Oriented with
    ! the North Pole and Greenwich as the First Point
    ! Gaussian Output Data is Interpolated to be Oriented with
    ! the North Pole and Greenwich as the First Point
    ! Input for the AGCM is Assumed IEEE 32 Bits Big Endian
    ! p_Lat0 and p_Lon0 could be diferent for diferent modules

    isExecOk = .false.
    if (var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_Lat0, p_Lon0)
    else
      flagInput(1) = .true.
      ! Start at North Pole
      flagInput(2) = .true.
      ! Start at Prime Meridian
      flagInput(3) = .true.
      ! Latitudes Are at Center of Box
      flagInput(4) = .true.
      ! Longitudes Are at Center of Box
      flagInput(5) = .false.
      ! Regular Grid
      flagOutput(1) = .true.
      ! Start at North Pole
      flagOutput(2) = .true.
      ! Start at Prime Meridian
      flagOutput(3) = .false.
      ! Latitudes Are at North Edge of Box
      flagOutput(4) = .true.
      ! Longitudes Are at Center of Box
      flagOutput(5) = .true.
      ! Gaussian Grid
      call initAreaInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    endif
    call allocateData()

    inquire (iolength = inRecSize) AlbedoIn
    albedoInFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
    inFileUnit = openFile(trim(albedoInFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return

    inquire (iolength = outRecSize) AlbedoOut
    albedoOutFilename = getOutFileName()
    outFileUnit = openFile(trim(albedoOutFilename), 'unformatted', 'direct', outRecSize, 'write', 'replace')
    if(outFileUnit < 0) return

    ! Read In Input Albedo
    do month = 1, 12
      read (unit = inFileUnit, rec = month) AlbedoIn
      AlbedoInput = real(Albedoin, p_r8)
      ! Interpolate Input Regular Grid Albedo To Gaussian Grid on Output
      if (var%linear) then
        call doLinearInterpolation (AlbedoInput, AlbedoOutput)
      else
        call doAreaInterpolation (AlbedoInput, AlbedoOutput)
      endif
      AlbedoOut = real(AlbedoOutput, p_r4)
      ! Write Out Adjusted Interpolated Albedo
      write (unit = outFileUnit, rec = month) AlbedoOut
    enddo

    close (unit = inFileUnit)
    close (unit = outFileUnit)

    if(varCommon%grads) then
      if(.not. generateGrads()) return
    endif
    call deallocateData()

    isExecOk = .true.
  end function generateAlbedo


  function generateGrads() result(isExecOk)
    !# Generates Grads files
    !# ---
    !# @info
    !# **Brief:** Generates .ctl and .dat files to check output Grads files. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo

    implicit none
    logical :: isExecOk
    !# true if execution was sucessfull
    character (len = maxPathLength) :: ctlPathFileName
    character (len = maxPathLength) :: gradsBaseName
    character (len = maxPathLength) :: ctlFileName
    integer :: ctlFileUnit
    !# Temp File Unit

    isExecOk = .false.
    gradsBaseName = trim(var%varNameOut) // nLats
    ctlFileName = trim(gradsBaseName) // '.ctl'
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(ctlFileName)

    ctlFileUnit = openFile(ctlPathFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if(ctlFileUnit < 0) return

    write (unit = ctlFileUnit, fmt = '(A)') 'DSET ^' // trim(gradsBaseName)
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,1PG12.5)') 'UNDEF ', -999.0_p_r8
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'TITLE Albedo on a Gaussian Grid'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
      'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write (unit = ctlFileUnit, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
    if (var%linear) then
      write (unit = ctlFileUnit, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write (unit = ctlFileUnit, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
    endif
    write (unit = ctlFileUnit, fmt = '(A)') 'ZDEF  1 LEVELS 1000'
    write (unit = ctlFileUnit, fmt = '(A)') 'TDEF 12 LINEAR JAN2005 1MO'
    write (unit = ctlFileUnit, fmt = '(A)') 'VARS  1'
    write (unit = ctlFileUnit, fmt = '(A)') 'ALBE  0 99 Albedo [No Dim]'
    write (unit = ctlFileUnit, fmt = '(A)') 'ENDVARS'
    close (unit = ctlFileUnit)

    close(unit = ctlFileUnit)
    isExecOk = .true.
  end function generateGrads


  subroutine printNameList()
    !# Prints Namelist of Albedo module
    !# ---
    !# @info
    !# **Brief:** Prints Namelist of Albedo module, defined in PRE_run.nml. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    write (unit = p_nfprt, FMT = '(/,A)')  ' &AlbedoNameList'
    write (unit = p_nfprt, FMT = '(A,I6)') '      xDim = ', var%xDim
    write (unit = p_nfprt, FMT = '(A,I6)') '      yDim = ', var%yDim
    write (unit = p_nfprt, FMT = '(A,I6)') '      zDim = ', var%zDim
    write (unit = p_nfprt, FMT = '(A,I6)') '      xMax = ', varCommon%xMax
    write (unit = p_nfprt, FMT = '(A,I6)') '      yMax = ', varCommon%yMax
    write (unit = p_nfprt, FMT = '(A,L6)') '     grads = ', varCommon%grads
    write (unit = p_nfprt, FMT = '(A,L6)') '    linear = ', var%linear

    write (unit = p_nfprt, FMT = '(A)') '  dirPreIn = ' // trim(varCommon%dirPreIn)
    write (unit = p_nfprt, FMT = '(A)') ' dirPreOut = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, FMT = '(A)') '   varName = ' // trim(var%varName)
    write (unit = p_nfprt, FMT = '(A)') '   varNameOut = ' // trim(var%varNameOut)
  end subroutine printNameList


  subroutine allocateData()
    !# Allocates data
    !# ---
    !# @info
    !# **Brief:** Allocates matrixes based in namelist. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    implicit none
    allocate(AlbedoIn(var%xDim, var%yDim))
    allocate(AlbedoInput(var%xDim, var%yDim))
    allocate(AlbedoOut(varCommon%xMax, varCommon%yMax))
    allocate(AlbedoOutput(varCommon%xMax, varCommon%yMax))
  end subroutine allocateData


  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Deallocates matrixes. </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    implicit none
    deallocate(AlbedoIn)
    deallocate(AlbedoInput)
    deallocate(AlbedoOut)
    deallocate(AlbedoOutput)
  end subroutine deallocateData


end module Mod_Albedo

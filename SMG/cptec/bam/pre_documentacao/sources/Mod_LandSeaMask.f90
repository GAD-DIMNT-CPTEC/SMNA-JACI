!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_LandSeaMask </br></br>
!#
!# **Brief**: Module responsible for generating the mask land/sea interpolated 
!# to the Gaussian grid of the model. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/WaterNavy.dat
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/LandSeaMask.GZZZZZ (ex.: pre/dataout/LandSeaMask.G00450) </br>
!# &bull; pre/dataout/LandSeaMaskNavy.GZZZZZ.dat </br>
!# &bull; pre/dataout/LandSeaMaskNavy.GZZZZZ.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.1 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-08-2011 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita          - version: 1.1.1 </li>
!#  <li>21-03-2019 - Denis Eiras     - version: 2.0.0 </li>
!#  <li>29-05-2019 - Eduardo Khamis  - version: 2.0.1 </li>
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


module Mod_LandSeaMask

  use Mod_AreaInterpolation, only : gLats, initAreaInterpolation, &
    doAreaInterpolation
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  public :: getNameLandSeaMask, initLandSeaMask, generateLandSeaMask, shouldRunLandSeaMask

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !  input parameters

  type LandSeaMaskNameListData
    integer :: xDim
    integer :: yDim
    character (len = 32) :: varNameIn = 'WaterNavy '
    character (len = 32) :: varNameOut = 'LandSeaMask '
  end type LandSeaMaskNameListData

  type(varCommonNameListData)   :: varCommon
  type(LandSeaMaskNameListData) :: var
  namelist /LandSeaMaskNameList/   var


  ! internal variables

  logical :: flagInput(5), flagOutput(5)
  character (len = 7) :: nLats = '.G     '
  character (len = 10) :: mskfmt = '(      I1)'

  ! Horizontal Area Interpolator
  ! Interpolate Regular To Gaussian
  ! Regular Input Data is Assumed to be Oriented with
  ! the North Pole and Greenwich as the First Point
  ! Gaussian Output Data is Interpolated to be Oriented with
  ! the North Pole and Greenwich as the First Point
  ! Input for the AGCM is Assumed IEEE 32 Bits Big Endian
  integer :: j, lRecIn, lRecOut

  integer, dimension (:, :), allocatable :: landSeaMaskOut
  real (kind = p_r4), dimension (:, :), allocatable :: waterIn
  real (kind = p_r4), dimension (:, :), allocatable :: waterOut
  real (kind = p_r8), dimension (:, :), allocatable :: waterInput
  real (kind = p_r8), dimension (:, :), allocatable :: waterOutput


contains


  function getNameLandSeaMask() result(returnModuleName)
    !# Returns LandSeaMask Module Name
    !# ---
    !# @info
    !# **Brief:** Returns LandSeaMask Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    
    returnModuleName = "LandSeaMask"
  end function getNameLandSeaMask


  subroutine initLandSeaMask(nameListFileUnit, varCommon_)
    !# Initializates LandSeaMask module
    !# ---
    !# @info
    !# **Brief:** Initializates LandSeaMask module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = LandSeaMaskNameList)
    varCommon = varCommon_

    write (nLats(3:7), '(I5.5)') varCommon%yMax
    write (mskfmt(2:7), '(I6)') varCommon%xMax
  end subroutine initLandSeaMask


  function shouldRunLandSeaMask() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunLandSeaMask


  function getOutFileName() result(landSeaMaskOutFilename)
    !# Gets LandSeaMask Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets LandSeaMask Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: landSeaMaskOutFilename

    landSeaMaskOutFilename = trim(varCommon%dirPreOut) // trim(var%varNameOut) // nLats
  end function getOutFileName


  subroutine printNameList()
    !# Prints NameList
    !# ---
    !# @info
    !# **Brief:** Prints NameList. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    write (unit = *, fmt = '(/,A)')  ' &LandSeaMaskNameList'
    write (unit = *, fmt = '(A,I6)') '        Imax = ', varCommon%xMax
    write (unit = *, fmt = '(A,I6)') '        Jmax = ', varCommon%yMax
    write (unit = *, fmt = '(A,I6)') '        Idim = ', var%xDim
    write (unit = *, fmt = '(A,I6)') '        Jdim = ', var%yDim
    write (unit = *, fmt = '(A,L6)') '       GrADS = ', varCommon%grads
    write (unit = *, fmt = '(A)')    '   VarNameIn = ' // trim(var%varNameIn)
    write (unit = *, fmt = '(A)')    '  VarNameOut = ' // trim(var%varNameOut)
    write (unit = *, fmt = '(A)')    '   DirPreOut = ' // trim(varCommon%dirPreOut)
    write (unit = *, fmt = '(A,/)')  ' /'
  end subroutine printNameList


  function generateLandSeaMask() result(isExecOk)
    !# Generates LandSeaMask
    !# ---
    !# @info
    !# **Brief:** Generates LandSeaMask. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer :: nfclm ! To read Water Percentage data (WPD)
    integer :: nfoua ! To write Intepolated Land Sea Mask (LSM)
    integer :: nfoub ! To write Intepolated WPD and LSM
    integer :: nfctl ! To write Output data Description
    character(len = maxPathLength) :: inputWaterFileName 
    character(len = maxPathLength) :: gradsWaterOutFileName 
    character(len = maxPathLength) :: ctlWaterOutFileName 
    logical :: isExecOk

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

    call initAreaInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)

    allocate (waterInput(var%xDim, var%yDim))
    allocate (waterOutput(varCommon%xMax, varCommon%yMax))
    allocate (waterIn(var%xDim, var%yDim))
    allocate (waterOut(varCommon%xMax, varCommon%yMax))
    allocate (landSeaMaskOut(varCommon%xMax, varCommon%yMax))

    ! Read In Input Water

    inquire (iolength = lRecIn) waterIn(:, 1)
    inputWaterFileName = trim(varCommon%dirPreOut) // trim(var%varNameIn) // '.dat'
    nfclm = openFile(trim(inputWaterFileName), 'unformatted', 'direct', lRecIn, 'read', 'old')
    if(nfclm < 0) return

    do j = 1, var%yDim
      read (unit = nfclm, rec = j) waterIn(:, j)
    end do
    close (unit = nfclm)
    waterInput = real(waterIn, p_r8)

    ! Interpolate Input Water To Output Water

    call doAreaInterpolation (waterInput, waterOutput)
    waterOut = real(waterOutput, p_r4)

    ! Generate Land Sea Mask (= 1 over Land and = 0 over Sea)
    where (waterOutput < 50.0_p_r8)
      landSeaMaskOut = 1
    elsewhere
      landSeaMaskOut = 0
    endwhere

    nfoua = openFile(trim(getOutFileName()), 'unformatted', 'sequential', -1, 'write', 'replace')
    if(nfoua < 0) return

    write (unit = nfoua) landSeaMaskOut
    close (unit = nfoua)

    if (varCommon%grads) then
      ! Write Out Adjusted Interpolated Water
      inquire (iolength = lRecOut) waterOut
      gradsWaterOutFileName = trim(varCommon%dirPreOut) // trim(var%varNameOut) // nLats // '.dat'
      nfoub = openFile(trim(gradsWaterOutFileName), 'unformatted', 'direct', lRecOut, 'write', 'replace')
      if(nfoub < 0) return

      write (unit = nfoub, rec = 1) waterOut
      write (unit = nfoub, rec = 2) real(landSeaMaskOut, p_r4)
      close (unit = nfoub)

      ! Write Grads Control File
      ctlWaterOutFileName = trim(varCommon%dirPreOut) // trim(var%varNameOut) // nLats // '.ctl'

      nfctl = openFile(trim(ctlWaterOutFileName), 'formatted', 'sequential', -1, 'write', 'replace')
      if(nfctl < 0) return
      write (unit = nfctl, fmt = '(A)') 'DSET ^' // &
        trim(var%varNameOut) // nLats // '.dat'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'TITLE Land Sea Mask on a Gaussian Grid'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
        'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0 / real(varCommon%xMax, p_r8)
      write (unit = nfctl, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
      write (unit = nfctl, fmt = '(8F10.5)') gLats(varCommon%yMax:1:-1)
      write (unit = nfctl, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
      write (unit = nfctl, fmt = '(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      write (unit = nfctl, fmt = '(A)') 'VARS 2'
      write (unit = nfctl, fmt = '(A)') 'WPER  0 99 Percentage of Water [%]'
      write (unit = nfctl, fmt = '(A)') 'LSMK  0 99 Land Sea Mask [1-Land 0-Sea]'
      write (unit = nfctl, fmt = '(A)') 'ENDVARS'
      close (unit = nfctl)
    end if

    isExecOk = .true.
  end function generateLandSeaMask


end module Mod_LandSeaMask

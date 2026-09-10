!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_VarTopo </br></br>
!#
!# **Brief**: Module responsible for generating the average topography and the 
!# variance of topography to the Gaussian grid of the model </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/TopoNavy.dat </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/TopoVariance.GZZZZZ (Ex.: model/datain/TopoVariance.G00450) </br>
!# &bull; pre/dataout/Topography.GZZZZZ </br>
!# &bull; pre/dataout/VarTopoNavy.GZZZZZ.ctl </br>
!# &bull; pre/dataout/VarTopoNavy.GZZZZZ.dat or (pre/dataout/VarTopoNavy.GZZZZZ)
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 1.2.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>18-09-2007 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita          - version: 1.1.1 </li>
!#  <li>03-05-2019 - Eduardo Khamis  - version: 1.2.0 </li>
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

module Mod_VarTopo

  use Mod_LinearInterpolation, only : gLatsL => LatOut, &
    initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, &
    initAreaInterpolation, doAreaInterpolation
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : &
    msgOut

  implicit none

  public :: initVarTopo, getNameVarTopo, generateVarTopo, shouldRunVarTopo

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !  input parameters

  integer, parameter :: maxNumberOfRegions = 100

  type VarTopoNameListData
    integer             :: xDim
    integer             :: yDim
    logical             :: linear = .false.
    character(len = maxPathLength) :: varName = 'TopoNavy'     
    !# file name prefix for writing Navy Topography data
    character(len = maxPathLength) :: varNameG = 'VarTopoNavy' 
    !# file name prefix for writing Navy Topography data in grADS
    character(len = 50) :: regionName(maxNumberOfRegions)
    real                :: latIni(maxNumberOfRegions)
    real                :: latEnd(maxNumberOfRegions)
    real                :: lonIni(maxNumberOfRegions)
    real                :: lonEnd(maxNumberOfRegions)
    logical             :: useTopoCutOff(maxNumberOfRegions) 
    real                :: topoHeightCutOff(maxNumberOfRegions)
    logical             :: useShapiroFilter(maxNumberOfRegions) 
    real                :: shapiroCoef(maxNumberOfRegions) 
  end type VarTopoNameListData

  type(varCommonNameListData) :: varCommon
  type(VarTopoNameListData)   :: var
  namelist /VarTopoNameList/     var

  ! internal variables

  real (kind = p_r8) :: lon0, lat0
  logical :: flagInput(5), flagOutput(5)

  character (len = 7) :: nLats = '.G     '
  character (len = 10) :: varNameT = 'Topography'
  character (len = 12) :: varNameV = 'TopoVariance'

  real (kind = p_r8) :: varMin = 0.0_p_r8   ! Minimum of Height variance in m2
  !real (kind = p_r8) :: varMax = 1.6E5_p_r8 ! Maximum of Height variance in m2

  !   character (len = 528) :: dirMain
  !   character (len = 11) :: nameNML='VarTopo.nml'

  integer :: nfclm  ! To read Topography data (GTOP30 or Navy)
  integer :: nfoua  ! To write Intepoloated Topography Variance data
  integer :: nfoub  ! To write Intepolated Topography data
  integer :: nfouc  ! To write GrADS Topography and Variance data
  integer :: nfctl  ! To write Output data Description

  ! Horizontal Area Interpolator
  ! Interpolate Regular To Gaussian
  ! Regular Input Data is Assumed to be Oriented with
  ! the North Pole and Greenwich as the First Point

  ! Set Undefined Value for Input Data at Locations which
  ! are not to be Included in Interpolation

  integer :: j, lRecIn, lRecOut

  real(kind = p_r4), dimension (:, :), allocatable :: topoIn
  real(kind = p_r4), dimension (:, :), allocatable :: fieldOut
  real(kind = p_r8), dimension (:, :), allocatable :: topoInput
  real(kind = p_r8), dimension (:, :), allocatable :: topoOutput
  real(kind = p_r8), dimension (:, :), allocatable :: varTopInput
  real(kind = p_r8), dimension (:, :), allocatable :: varTopOutput

  character(len = *), parameter :: headerMsg = 'Var Topo              |'

contains


  function getNameVarTopo() result(returnModuleName)
    !# Returns VarTopo Module Name
    !# ---
    !# @info
    !# **Brief:** Returns VarTopo Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "VarTopo"
  end function getNameVarTopo


  function shouldRunVarTopo() result(shouldRun)
    !# Returns Returns true if Module Should Run as a dependency
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

    shouldRun = .not. fileExists(getOutTopographyFileName()) .or. .not. fileExists(getOutTopographyVarianceFileName())
  end function shouldRunVarTopo


  function getOutTopographyFileName() result(topoOutFilename)
    !# Gets Topography Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets Topography Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: topoOutFilename

    topoOutFilename = trim(varCommon%dirPreOut) // varNameT // nLats
  end function getOutTopographyFileName


  function getOutTopographyVarianceFileName() result(varTopoOutFilename)
    !# Gets Topography Variance Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets Topography Variance Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: varTopoOutFilename

    varTopoOutFilename = trim(varCommon%dirModelIn) // varNameV // nLats
  end function getOutTopographyVarianceFileName


  subroutine initVarTopo(nameListFileUnit, varCommon_)
    !# Initializes VarTopo
    !# ---
    !# @info
    !# **Brief:** Initializes VarTopo. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = VarTopoNameList)
    varCommon = varCommon_

    ! For Linear Interpolation
    lon0 = 0.0_p_r8  ! Start at Prime Meridian
    lat0 = 90.0_p_r8 - 0.5_p_r8 * (360.0_p_r8 / real(var%xDim, p_r8)) ! Start at North Pole

    ! For Area Weighted Interpolation
    flagInput(1) = .true.   ! Start at North Pole
    flagInput(2) = .true.   ! Start at Prime Meridian
    flagInput(3) = .false.  ! Latitudes Are at North Edge of Box
    flagInput(4) = .false.  ! Longitudes Are at Western Edge of Box
    flagInput(5) = .false.  ! Regular Grid
    flagOutput(1) = .true.  ! Start at North Pole
    flagOutput(2) = .true.  ! Start at Prime Meridian
    flagOutput(3) = .false. ! Latitudes Are at North Edge of Box
    flagOutput(4) = .true.  ! Longitudes Are at Center of Box
    flagOutput(5) = .true.  ! Gaussian Grid

    write (nLats(3:7), fmt = '(I5.5)') varCommon%yMax

  end subroutine initVarTopo


  subroutine printNameList()
    !# Prints NameList
    !# ---
    !# @info
    !# **Brief:** Prints NameList. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: ago/2019 </br>
    !# @endinfo
    implicit none

    write (unit = p_nfprt, fmt = '(/,A)')  ' &VarTopoNameList'
    write (unit = p_nfprt, fmt = '(A,I6)') '      xMax = ', varCommon%xMax
    write (unit = p_nfprt, fmt = '(A,I6)') '      yMax = ', varCommon%yMax
    write (unit = p_nfprt, fmt = '(A,I6)') '      xDim = ', var%xDim
    write (unit = p_nfprt, fmt = '(A,I6)') '      yDim = ', var%yDim
    write (unit = p_nfprt, fmt = '(A,L6)') '     grads = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A,L6)') '    linear = ', var%linear
    write (unit = p_nfprt, fmt = '(A)')    '   varName = ' // trim(var%varName)
    write (unit = p_nfprt, fmt = '(A)')    '  varNameG = ' // trim(var%varNameG)
    write (unit = p_nfprt, fmt = '(A)')    '   dirPreOut  = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, fmt = '(A)')    '   dirModelIn = ' // trim(varCommon%dirModelIn)
    write (unit = p_nfprt, fmt = '(A,/)')  ' /'

  end subroutine printNameList


  function generateVarTopo() result(isExecOk)
    !# Generates VarTopo
    !# ---
    !# @info
    !# **Brief:** Generates VarTopo. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk

    real (kind = p_r8) :: dlat, dlon
    real (kind = p_r8) :: aa, bb, cc, aa1, aa2, bb1, bb2, ss1, ss2
    integer :: i
    real (kind = p_r8), dimension (:), allocatable :: lat2
    real (kind = p_r8), dimension (:), allocatable :: lon2
    real (kind = p_r8), dimension (:, :), allocatable :: orog  !filtrado
    real (kind = p_r8), dimension (:, :), allocatable :: orog2 !antes de filtrar
    integer :: idxRegion

    isExecOk = .false.

    if (var%linear) then
      call initLinearInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, lat0, lon0)
    else
      call initAreaInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if

    allocate (topoIn(var%xDim, var%yDim))
    allocate (fieldOut(varCommon%xMax, varCommon%yMax))
    allocate (topoInput(var%xDim, var%yDim))
    allocate (topoOutput(varCommon%xMax, varCommon%yMax))
    allocate (varTopInput(var%xDim, var%yDim))
    allocate (varTopOutput(varCommon%xMax, varCommon%yMax))

    allocate (lat2(var%yDim))
    allocate (lon2(var%xDim))
    allocate (orog(var%xDim, var%yDim))
    allocate (orog2(var%xDim, var%yDim))

    ! Read In Input Topo
    inquire (iolength = lRecIn) topoIn(:, 1)
    nfclm = openFile(trim(varCommon%dirPreOut) // trim(var%varName) // '.dat', &
      'unformatted', 'direct', lRecIn, 'read', 'old')
    if(nfclm < 0) return
    do j = 1, var%yDim
      read (unit = nfclm, rec = j) topoIn(:, j)
    end do
    close (unit = nfclm)
    topoInput = real(topoIn, p_r8)

    orog = 0.0_p_r8
    orog2 = 0.0_p_r8
    dlat = 0.0
    dlon = 0.0
    lat2 = 0.0
    lon2 = 0.0
    ss1 = 0.0_p_r8
    ss2 = 0.0_p_r8

    !  INPUT FULL TOPO TopogIn (in r4) e Topog (em p_r8)
    orog = TopoInput
    dlat = 180.0_p_r8 / float(var%yDim)
    dlon = 360.0_p_r8 / float(var%xDim)
    lat2(1) = 90.0_p_r8
    lat2(var%yDim) = -90.0_p_r8
    lon2(1) = 0.0_p_r8

    do j = 2, var%yDim / 2
      lat2(j) = 90.0_p_r8 - dlat * J
    enddo

    i = 0
    do j = var%yDim / 2 + 1, var%yDim
      i = i + 1
      lat2(j) = -dlat * i
    enddo

    do i = 2, var%xDim
      lon2(i) = 0.0 + dlon * (i - 1)
    enddo

    ! loop over all regions configured at VarTopo namelist (var%region(#region)% ...)
    do idxRegion = 1, maxNumberOfRegions
      if(trim(var%regionName(idxRegion)) == 'UNDEF') cycle

      ! Apply Topo Cut Off
      if(var%useTopoCutOff(idxRegion)) then
        call msgOut(headerMsg, 'Applying Topo Cut off over region ' //  trim(var%regionName(idxRegion)))
        do j = 1, var%yDim
          aa = lat2(j)
          aa1 = var%latIni(idxRegion)
          aa2 = var%latEnd(idxRegion)
          do i = 1, var%xDim
            bb = lon2(i)
            bb1 = 360.0_p_r8 + var%lonIni(idxRegion)
            bb2 = 360_p_r8 + var%lonEnd(idxRegion)

            if(aa.ge.aa1.and.aa.le.aa2)then
              if(bb.ge.bb1.and.bb.le.bb2)then
                !if(orog(i, j).ge.20.0_p_r8)orog(i, j) = 20.0_p_r8
                if(orog(i, j) .ge. var%topoHeightCutOff(idxRegion)) then
                  orog(i, j) = var%topoHeightCutOff(idxRegion)
                end if
              endif
            endif

          enddo
        end do
       endif

      ! Apply ShapiroFilter
      if(var%useShapiroFilter(idxRegion)) then
        call msgOut(headerMsg, 'Applying Shapiro filter over region ' //  trim(var%regionName(idxRegion)))

        orog2 = orog
        do j = 2, var%yDim - 1
          aa = lat2(j)
          aa1 = var%latIni(idxRegion)
          aa2 = var%latEnd(idxRegion)
          do i = 2, var%xDim -1
            bb = lon2(i)
            bb1 = 360.0_p_r8 + var%lonIni(idxRegion)
            bb2 = 360_p_r8 + var%lonEnd(idxRegion)

            if(aa.ge.aa1.and.aa.le.aa2)then
              if(bb.ge.bb1.and.bb.le.bb2)then
                  ss1 = orog2(i, j)
                  ss2 = orog2(i - 1, j) + orog2(i + 1, j) + orog2(i, j - 1) + orog2(i, j + 1)
                  orog(i, j) = ss1 + var%shapiroCoef(idxRegion) * (ss2 - 4.0_p_r8 * ss1)
              endif
            endif

          enddo
        enddo
      endif

    enddo

    !BACK TO Topog
    topoInput = orog

    varTopInput = topoInput * topoInput

    ! Interpolate Input Regular Grid Albedo To Gaussian Grid on Output
    if (var%linear) then
      call doLinearInterpolation (topoInput, topoOutput)
      call doLinearInterpolation (varTopInput, varTopOutput)
    else
      call doAreaInterpolation (topoInput, topoOutput)
      call doAreaInterpolation (varTopInput, varTopOutput)
    end if
    varTopOutput = max(varTopOutput - topoOutput * topoOutput, varMin)
    !varTopOutput=min(varTopOutput,varMax)

    inquire (iolength = lRecOut) fieldOut

    ! Write Out Adjusted Interpolated Topography Variance (m2)
    nfoua = openFile(getOutTopographyVarianceFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(nfoua < 0) return
    fieldOut = real(varTopOutput, p_r4)
    write (unit = nfoua, rec = 1) fieldOut
    close (unit = nfoua)

    ! Write Out Adjusted Interpolated Topography (m)
    nfoub = openFile(getOutTopographyFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(nfoub < 0) return
    fieldOut = real(topoOutput, p_r4)
    write (unit = nfoub, rec = 1) fieldOut
    close (unit = nfoub)

    if (varCommon%grads) then
      call generateGrads()
    endif

    isExecOk = .true.
  end function generateVarTopo


  subroutine generateGrads()
    !# Generates Grads
    !# ---
    !# @info
    !# **Brief:** Generates Grads. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none

    ! Write Out Adjusted Interpolated Topography and Variance
    nfouc = openFile(trim(varCommon%dirPreOut) // trim(var%varNameG) // nLats, &
      'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(nfouc < 0) return
    fieldOut = real(topoOutput, p_r4)
    write (unit = nfouc, rec = 1) fieldOut
    fieldOut = real(varTopOutput, p_r4)
    write (unit = nfouc, rec = 2) fieldOut
    close (unit = nfouc)

    ! Write GrADS Control File
    nfctl = openFile(trim(varCommon%dirPreOut) // trim(var%varNameG) // nLats // '.ctl', &
      'formatted', 'sequential', -1, 'write', 'replace')
    if(nfctl < 0) return
    write (unit = nfctl, fmt = '(A)') 'DSET ^' // &
      trim(var%varNameG) // nLats
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A)') 'TITLE Topography and Variance on a Gaussian Grid'
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
      'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write (unit = nfctl, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
    if (var%linear) then
      write (unit = nfctl, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write (unit = nfctl, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
    end if
    write (unit = nfctl, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
    write (unit = nfctl, fmt = '(A)') 'TDEF 1 LINEAR JAN2005 1MO'
    write (unit = nfctl, fmt = '(A)') 'VARS 2'
    write (unit = nfctl, fmt = '(A)') 'TOPO 0 99 Topography [m]'
    write (unit = nfctl, fmt = '(A)') 'VART 0 99 Variance of Topography [m2]'
    write (unit = nfctl, fmt = '(A)') 'ENDVARS'
    close (unit = nfctl)
  end subroutine generateGrads


end module Mod_VarTopo

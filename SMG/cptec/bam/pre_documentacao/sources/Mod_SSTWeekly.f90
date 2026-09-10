!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SSTWeekly </br></br>
!#
!# **Brief**: Module responsible for correcting the weekly sea surface temperature 
!! (SST) by topography and generate the sea ice mask (1: Sea Ice 0: No Ice) </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/databcs/sstaoi.form (if climWindow=.true.) </br>
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ  (Ex.: pre/dataout/ModelLandSeaMask.G00450) </br>
!# &bull; model/datain/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ  (Ex.: model/datain/GANLNMC2015043000S.unf.TQ0299L064) </br>
!# &bull; pre/dataout/SSTWeekly.YYYYMMDD
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; model/datain/SSTWeeklyYYYYMMDD.GZZZZZ </br>
!# &bull; pre/dataout/SSTWeeklyYYYYMMDD.GZZZZZ </br>
!# &bull; pre/dataout/SSTWeeklyYYYYMMDD.GZZZZZ.ctl (if grads=.true.)
!# </br></br>
!# 
!# SSTWeekly20150430.G00450 from pre/dataout and model/datain are different </br></br>
!# 
!# The SSTWeeklyYYYYMMDD.GZZZZZ file that is in the model/datain directory is the 
!# file that the model will use in its execution, has the specific format (unformatted, 
!# direct access, 64-byte integer) and the SSTWeeklyYYYYMMDD.GZZZZZ file that is in 
!# the pre/dataout directory is specific to viewing in GrADS software. </br></br>
!# 
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.1 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti         - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita                  - version: 1.1.1 </li>
!#  <li>08-05-2019 - Eduardo Khamis          - version: 2.0.1 </li>
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

module Mod_SSTWeekly

  use Mod_FastFourierTransform, only : createFFT, destroyFFT
  use Mod_LegendreTransform, only : createGaussRep, createSpectralRep, createLegTrans, destroyLegendreObjects
  use Mod_SpectralGrid, only : transp, specCoef2Grid
  use Mod_LinearInterpolation, only : gLatsL => latOut, &
    initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, &
    initAreaInterpolation, doAreaInterpolation
  use Mod_Get_xMax_yMax, only : getxMaxyMax
  use Mod_FileManager
  use Mod_Namelist,             only : varCommonNameListData

  implicit none
  private
  include 'pre.h'
  include 'files.h'
  include 'precision.h'

  public :: getNameSSTWeekly, initSSTWeekly, generateSSTWeekly, shouldRunSSTWeekly

  ! input parameters

  integer :: xmx, yMaxHf, &
    mEnd1, mEnd2, mnwv0, mnwv1, mnwv2, mnwv3
  integer :: year, month, day, hour

  real (kind = p_r8) ::  p_lon0, p_lat0, p_to, sstOpenWater, &
    sstSeaIceThreshold, lapseRate

  logical :: flagInput(5), flagOutput(5)

  character (len = 12) :: gradsTime = '  Z         '
  character (len = 10) :: trunc = 'T     L   '
  character (len = 7) :: nLats = '.G     '
  character (len = 10) :: mskfmt = '(      I1)'
  character (len = 9) :: varName = 'SSTWeekly'
  character (len = 16) :: nameLSM = 'ModelLandSeaMask'
!  character (len = 528) :: dirPreOut != 'pre/dataout/'
!  character (len = 528) :: dirModelIn != 'model/datain/'
!  character (len = 528) :: dirClmSST != 'pre/databcs/' ! Climatological SST Datain Directory
  !  character (len = 12) :: dirPreOut != 'pre/dataout/'
  !  character (len = 13) :: dirModelIn != 'model/datain/'
  !  character (len = 12) :: dirClmSST != 'pre/databcs/' ! Climatological SST Datain Directory
  character (len = 11) :: fileClmSST = 'sstaoi.form'
!  character (len = 7) :: preffix
!  character (len = 6) :: suffix
!  character (len = 10) :: dateICn
  character (len = 3), dimension (12) :: monthChar = &
    (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
      'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

  integer, dimension (12) :: monthLength = &
    (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  !  integer :: p_nferr = 0    ! Standard Error print Out
  !  integer :: p_nfinp = 5    ! Standard read in
  integer :: p_nfprt = 6    ! Standard print Out
  !  integer :: p_nficn = 10   ! To read Topography from Initial Condition
  !  integer :: p_nflsm = 20   ! To read Formatted Land Sea Mask
  !  integer :: p_nfsti = 30   ! To read Unformatted 1x1 SST
  !  integer :: p_nfclm = 40   ! To read Formatted Climatological SST
  !  integer :: p_nfsto = 50   ! To write Unformatted Gaussian Grid SST
  !  integer :: p_nfout = 60   ! To write GrADS Topography, Land Sea, Se Ice and Gauss SST
  !  integer :: p_nfctl = 70   ! To write GrADS Control file

  ! local

  ! Reads the Mean Weekly 1x1 SST Global From NCEP,
  ! Interpolates it Using Area Weigth Into a Gaussian Grid

  integer :: j, i, js, jn, js1, jn1, ja, jb
!  integer :: ios
  integer :: forecastDay

  real (kind = p_r4) :: timeOfDay
  real (kind = p_r8) :: rgSSTMax, rgSSTMin, ggSSTMax, ggSSTMin, mgSSTMax, mgSSTMin

  integer :: iCnDate(4), currentDate(4)
  integer, dimension (:, :), allocatable :: landSeaMask, seaIceMask

  real (kind = p_r4), dimension (:), allocatable :: coefTopIn
  real (kind = p_r4), dimension (:, :), allocatable :: sstWklIn, WrOut
  real (kind = p_r8), dimension (:), allocatable :: coefTop
  real (kind = p_r8), dimension (:, :), allocatable :: topog, sstClim, &
    sstIn, sstGaus, seaIceFlagIn, seaIceFlagOut
  character (len = 6), dimension (:), allocatable :: sstLabel

  type SSTWeeklyNameListData
    integer            :: zmax                 
    !# Number of Layers 
    integer            :: xDim                 
    !# Number of Longitudes
    integer            :: yDim                 
    !# Number of Latitudes
    real(kind = p_r8)  :: sstSeaIce = -1.749   
    !# SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
    real(kind = p_r8)  :: latClimSouth =-50.0  
    !# Southern Latitude For Climatological SST Data
    real(kind = p_r8)  :: latClimNorth = 60.0  
    !# Northern Latitude For Climatological SST Data
    logical            :: climWindow = .false. 
    !# Flag to Climatological SST Data Window
    logical            :: linear = .true.      
    !# Flag for Linear (T) or Area Weighted (F) Interpolation
    logical            :: linearGrid = .false. 
    !# Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
    character(len = 7) :: preffix = 'GANLCPT'  
    !# Preffix of the Initial Condition for the Global Model
    character(len = 6) :: suffix = 'S.unf.'    
    !# Suffix of the Initial Condition for the Global Model
  end type SSTWeeklyNameListData

  type(varCommonNameListData) :: varCommon
  type(SSTWeeklyNameListData) :: var
  namelist /SSTWeeklyNameList/   var

!  namelist /SSTWeeklyNameList/ mEnd, zMax, xDim, yDim, &
!      sstSeaIce, latClimSouth, latClimNorth, &
!      climWindow, linear, linearGrid, grads, &
!      dateICn, preffix, suffix, dirPreOut, dirModelIn, dirClmSST


contains


  function getNameSSTWeekly() result(returnModuleName)
    !# Returns SSTWeekly Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SSTWeekly Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    returnModuleName = "SSTWeekly"
  end function getNameSSTWeekly


  subroutine initSSTWeekly(nameListFileUnit, varCommon_)
    !# Initialization of SSTWeekly module
    !# ---
    !# @info
    !# **Brief:** Initialization of SSTWeekly module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    
    !    integer :: ios

    !    mEnd = 62            ! Spectral Resolution Horizontal Truncation
    !    zMax = 28            ! Number of Layers of the Initial Condition for the Global Model
    !    xDim = 360           ! Number of Longitudes For Climatological SST data
    !    yDim = 180           ! Number of Latitudes For Climatological SST data
    !    sstSeaIce = -1.749_p_r8   ! SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
    !    latClimSouth = -50.0_p_r8 ! Southern Latitude For Climatological SST data
    !    latClimNorth = 60.0_p_r8  ! Northern Latitude For Climatological SST data
    !    climWindow = .false.    ! Flag to Climatological SST data Window
    !    linear = .true.         ! Flag to Bi-linear (T) or Area (F) Interpolation
    !    linearGrid = .false.    ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
    !    grads = .true.          ! Flag for GrADS Outputs
    !    dateICn = 'yyyymmddhh'  ! Date of the Initial Condition for the Global Model
    !    preffix = 'GANLCPT'     ! Preffix of the Initial Condition for the Global Model
    !    suffix = 'S.unf.'       ! Suffix of the Initial Condition for the Global Model
    !    dirMain = './ '         ! Main Datain/Dataout Directory

    !    open (unit = p_nfinp, file = './' // nameNML, &
    !      form = 'formatted', access = 'sequential', &
    !      action = 'read', status = 'old', iostat = ios)
    !    if (ios /= 0) then
    !      write (unit = p_nferr, fmt = '(3A,I4)') &
    !        ' ** (Error) ** open file ', &
    !        './' // nameNML, &
    !        ' returned iostat = ', ios
    !      stop  ' ** (Error) **'
    !    end if
    !    read  (unit = p_nfinp, nml = inputDim)
    !    close (unit = p_nfinp)

    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
    
    read(unit = nameListFileUnit, nml = SSTWeeklyNameList)
    varCommon = varCommon_

    call getxMaxyMax(varCommon%mEnd, varCommon%xMax, varCommon%yMax)

    mEnd1 = varCommon%mEnd + 1
    mEnd2 = varCommon%mEnd + 2
    mnwv2 = mEnd1 * mEnd2
    mnwv0 = mnwv2 / 2
    mnwv3 = mnwv2 + 2 * mEnd1
    mnwv1 = mnwv3 / 2

    xmx = varCommon%xMax + 2
    yMaxHf = varCommon%yMax / 2

    p_to = 273.15_p_r8
    sstOpenWater = -1.7_p_r8 + p_to
    sstSeaIceThreshold = 271.2_p_r8

    lapseRate = 0.0065_p_r8

    read (varCommon%date, fmt = '(I4,3I2)') year, month, day, hour
    if (mod(year, 4) == 0) monthLength(2) = 29

    ! For Linear Interpolation
    p_lon0 = 0.5_p_r8  ! Start Near Greenwhich
    p_lat0 = 89.5_p_r8 ! Start Near North Pole

    ! For Area Weighted Interpolation
    if (var%linearGrid) then
      trunc(2:2) = 'L'
    else
      trunc(2:2) = 'Q'
    end if
    write (trunc(3:6), fmt = '(I4.4)') varCommon%mEnd
    write (trunc(8:10), fmt = '(I3.3)') var%zMax

    write (nLats(3:7), '(I5.5)') varCommon%yMax

    write (mskfmt(2:7), '(I6)') varCommon%xMax

    write (gradsTime(1:2), fmt = '(I2.2)') hour
    write (gradsTime(4:5), fmt = '(I2.2)') day
    write (gradsTime(6:8), fmt = '(A3)') monthChar(month)
    write (gradsTime(9:12), fmt = '(I4.4)') year
  end subroutine initSSTWeekly


  function shouldRunSSTWeekly() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunSSTWeekly


  function getOutFileName() result(sstWeeklyOutFilename)
    !# Gets SSTWeekly Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SSTWeekly Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: sstWeeklyOutFilename

    sstWeeklyOutFilename = trim(varCommon%dirModelIn) // trim(varName) // varCommon%date(1:8) // nLats
  end function getOutFileName


  function generateSSTWeekly() result(isExecOk)
    !# Generates SSTWeekly
    !# ---
    !# @info
    !# **Brief:** Generates SSTWeekly. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    
    implicit none
!    integer,             intent(in) :: xMax_  ! max value of longitude
!    integer,             intent(in) :: yMax_  ! max value of latitude
    integer :: inFileUnit, outFileUnit, outRecSize, lRecIn
    character (len = maxPathLength) :: sstWklInFilename, landSeaMaskInFilename, coefTopInFilename, wrOutFilename
    logical :: isExecOk

    isExecOk = .false.
!    xMax = xMax_
!    yMax = yMax_
!   TODO adjust common variables like varCommon%xMax
    write (unit = p_nfprt, fmt = '(/,A,I6)') '   Imax = ', varCommon%xMax
    write (unit = p_nfprt, fmt = '(A,I6,/)') '   Jmax = ', varCommon%yMax

    call createSpectralRep(mEnd1, mEnd2, mnwv1)
    call createGaussRep(varCommon%yMax, yMaxHf)
    call createFFT(varCommon%xMax)
    call createLegTrans(mnwv0, mnwv1, mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf)
    if (var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_lat0, p_lon0)
    else
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
      call initAreaInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if

    allocate (landSeaMask(varCommon%xMax, varCommon%yMax), seaIceMask(varCommon%xMax, varCommon%yMax))
    allocate (coefTopIn(mnwv2), sstWklIn(var%xDim, var%yDim))
    allocate (coefTop(mnwv2), topog(varCommon%xMax, varCommon%yMax), sstClim(var%xDim, var%yDim))
    allocate (sstIn(var%xDim, var%yDim), sstGaus(varCommon%xMax, varCommon%yMax), WrOut(varCommon%xMax, varCommon%yMax))
    allocate (seaIceFlagIn(var%xDim, var%yDim), seaIceFlagOut(varCommon%xMax, varCommon%yMax))
    allocate (sstLabel(var%yDim))

    ! Read in SpectraL Coefficient of Topography from ICn
    ! to Ensure that Topography is the Same as Used by Model
    coefTopInFilename = trim(varCommon%dirModelIn) // var%preffix // varCommon%date // var%suffix // trunc
    inFileUnit = openFile(trim(coefTopInFilename), 'unformatted', 'sequential', -1, 'read', 'old') ! lRecIn ?
    if(inFileUnit < 0) return
    read  (unit = inFileUnit) forecastDay, timeOfDay, iCnDate, currentDate
    read  (unit = inFileUnit) coefTopIn
    close (unit = inFileUnit)
    coefTop = real(coefTopIn, p_r8)
    call transp(mnwv2, mEnd1, mEnd2, coefTop)
    call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, coefTop, topog)

    ! Read in Land-Sea Mask Data Set
    landSeaMaskInFilename = trim(varCommon%dirPreOut) // nameLSM // nLats
    inFileUnit = openFile(trim(landSeaMaskInFilename), 'unformatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return
    ! Land Sea Mask : Input
    read  (unit = inFileUnit) landSeaMask
    close (unit = inFileUnit)

    ! Read Mean Weekly 1 deg x 1 deg SST
    inquire (iolength = lRecIn) sstWklIn
    sstWklInFilename = trim(varCommon%dirPreOut) // varName // '.' // varCommon%date(1:8)
    inFileUnit = openFile(trim(sstWklInFilename), 'unformatted', 'direct', lRecIn, 'read', 'old')
    if(inFileUnit < 0) return
    read  (unit = inFileUnit, rec = 1) sstWklIn
    close (unit = inFileUnit)
    if (maxval(sstWklIn) < 100.0_p_r8) sstWklIn = sstWklIn + p_to

    ! Get SSTClim Climatological and Index for High Latitude
    ! Substitution of SST Actual by Climatology
    if (var%climWindow) then
      if (.not. sstClimatological() ) return
      if (maxval(sstClim) < 100.0_p_r8) sstClim = sstClim + p_to
      call sstClimaWindow ()
      jn1 = jn - 1
      js1 = js + 1
      ja = jn
      jb = js
    else
      jn = 0
      js = var%yDim + 1
      jn1 = 0
      js1 = var%yDim + 1
      ja = 1
      jb = var%yDim
    end if

    do j = 1, var%yDim
      if (j >= jn .and. j <= js) then
        sstLabel(j) = 'Weekly'
        do i = 1, var%xDim
          sstIn(i, j) = real(sstWklIn(i, j), p_r8)
        end do
      else
        sstLabel(j) = 'Climat'
        do i = 1, var%xDim
          sstIn(i, j) = sstClim(i, j)
        end do
      end if
    end do
    if (jn1 >= 1) write (unit = p_nfprt, fmt = '(6(I4,1X,A))') (j, sstLabel(j), j = 1, jn1)
    write (unit = p_nfprt, fmt = '(6(I4,1X,A))') (j, sstLabel(j), j = ja, jb)
    if (js1 <= var%yDim) write (unit = p_nfprt, fmt = '(6(I4,1X,A))') (j, sstLabel(j), j = js1, var%yDim)

    ! Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
    ! Over Open Water Set Input SST = MIN of SSTOpenWater
    ! Over Non Ice Points Before Interpolation
    if (var%sstSeaIce < 100.0_p_r8) var%sstSeaIce = var%sstSeaIce + p_to
    do j = 1, var%yDim
      do i = 1, var%xDim
        seaIceFlagIn(i, j) = 0.0_p_r8
        if (sstIn(i, j) < var%sstSeaIce) then
          seaIceFlagIn(i, j) = 1.0_p_r8
        else
          sstIn(i, j) = max(sstIn(i, j), sstOpenWater)
        end if
      end do
    end do
    ! Min And Max Values of Input SST
    rgSSTMax = maxval(sstIn)
    rgSSTMin = minval(sstIn)

    ! Interpolate Flag from 1x1 Grid to Gaussian Grid, Fill SeaIceMask=1
    ! Over Interpolated Points With 50% or More Sea Ice, =0 Otherwise
    if (var%linear) then
      call doLinearInterpolation (seaIceFlagIn, seaIceFlagOut)
    else
      call doAreaInterpolation (seaIceFlagIn, seaIceFlagOut)
    end if
    seaIceMask = int(seaIceFlagOut + 0.5_p_r8)
    where (landSeaMask == 1) seaIceMask = 0

    ! Interpolate SST from 1x1 Grid to Gaussian Grid
    if (var%linear) then
      call doLinearInterpolation (sstIn, sstGaus)
    else
      call doAreaInterpolation (sstIn, sstGaus)
    end if
    ! Min and Max Values of Gaussian Grid
    ggSSTMax = maxval(sstGaus)
    ggSSTMin = minval(sstGaus)

    do j = 1, varCommon%yMax
      do i = 1, varCommon%xMax
        if (landSeaMask(i, j) == 1) then
          ! Set SST = Undef Over Land
          sstGaus(i, j) = p_undef
        else if (seaIceMask(i, j) == 1) then
          ! Set SST Sea Ice Threshold Minus 1 Over Sea Ice
          sstGaus(i, j) = sstSeaIceThreshold - 1.0_p_r8
        else
          ! Correct SST for Topography, Do Not Create or
          ! Destroy Sea Ice Via Topography Correction
          sstGaus(i, j) = sstGaus(i, j) - topog(i, j) * lapseRate
          if (sstGaus(i, j) < sstSeaIceThreshold) &
            sstGaus(i, j) = sstSeaIceThreshold + 0.2_p_r8
        end if
      end do
    end do
    ! Min and Max Values of Corrected Gaussian Grid SST Excluding Land Points
    mgSSTMax = maxval(sstGaus, MASK = sstGaus/=p_undef)
    mgSSTMin = minval(sstGaus, MASK = sstGaus/=p_undef)
    ! Write out Land-Sea Mask and SST Data to Global Model Input
    inquire (iolength = outRecSize) WrOut
    wrOutFilename = trim(varCommon%dirModelIn) // trim(varName) // varCommon%date(1:8) // nLats
    outFileUnit = openFile(trim(wrOutFilename), 'unformatted', 'direct', outRecSize, 'write', 'replace')
    if(outFileUnit < 0) return
    ! Write out Land-Sea Mask to SST Data Set
    ! The LSMask will be Transfered by Model to Post-Processing
    WrOut = real(1 - 2 * landSeaMask, p_r4)
    write (unit = outFileUnit, rec = 1) WrOut
    ! Write out Gaussian Grid Weekly SST
    WrOut = real(sstGaus, p_r4)
    write (unit = outFileUnit, rec = 2) WrOut
    close (unit = outFileUnit)

    write (unit = p_nfprt, fmt = '(/,3(A,I2.2),A,I4)') &
      ' Hour = ', hour, ' Day = ', day, &
      ' Month = ', month, ' Year = ', year

    write (unit = p_nfprt, fmt = '(/,A,3(A,2F8.2,/))') &
      ' Mean Weekly SST Interpolation :', &
      ' Regular  Grid SST: min, max = ', rgSSTMin, rgSSTMax, &
      ' Gaussian Grid SST: min, max = ', ggSSTMin, ggSSTMax, &
      ' Masked G Grid SST: min, max = ', mgSSTMin, mgSSTMax

    if (varCommon%grads) then
      wrOutFilename = trim(varCommon%dirPreOut) // trim(varName) // varCommon%date(1:8) // nLats
      outFileUnit = openFile(trim(wrOutFilename), 'unformatted', 'direct', outRecSize, 'write', 'replace')
      if(outFileUnit < 0) return
      WrOut = real(topog, p_r4)
      write (unit = outFileUnit, rec = 1) WrOut
      WrOut = real(1 - 2 * landSeaMask, p_r4)
      write (unit = outFileUnit, rec = 2) WrOut
      WrOut = real(seaIceMask, p_r4)
      write (unit = outFileUnit, rec = 3) WrOut
      WrOut = real(sstGaus, p_r4)
      write (unit = outFileUnit, rec = 4) WrOut
      close (unit = outFileUnit)

      call generateGrads()
    end if

    call destroyFFT()
    call destroyLegendreObjects()
    
    isExecOk = .true.
  end function generateSSTWeekly


  subroutine generateGrads()
    !# Generates Grads
    !# ---
    !# @info
    !# **Brief:** Generates Grads. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxPathLength) :: ctlPathFileName
    character (len = maxPathLength) :: gradsBaseName
    character (len = maxPathLength) :: ctlFileName
    integer :: ctlFileUnit      ! Temp File Unit

    ! Write GrADS Control File
    gradsBaseName = trim(varName) // varCommon%date(1:8) // nLats
    ctlFileName = trim(gradsBaseName) // '.ctl'
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(ctlFileName)
    ctlFileUnit = openFile(ctlPathFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if(ctlFileUnit < 0) return
    write (unit = ctlFileUnit, fmt = '(A)') 'DSET ' // &
      trim(varCommon%dirPreOut) // varName // varCommon%date(1:8) // nLats
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'TITLE Weekly SST on a Gaussian Grid'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
      'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write (unit = ctlFileUnit, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
    if (var%linear) then
      write (unit = ctlFileUnit, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write (unit = ctlFileUnit, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
    end if
    write (unit = ctlFileUnit, fmt = '(A)') 'ZDEF  1 LEVELS 1000'
    write (unit = ctlFileUnit, fmt = '(A)') 'TDEF  1 LINEAR ' // gradsTime // ' 6HR'
    write (unit = ctlFileUnit, fmt = '(A)') 'VARS  4'
    write (unit = ctlFileUnit, fmt = '(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
    write (unit = ctlFileUnit, fmt = '(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
    write (unit = ctlFileUnit, fmt = '(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
    write (unit = ctlFileUnit, fmt = '(A)') 'SSTW  0 99 Weekly SST Topography Corrected [K]'
    write (unit = ctlFileUnit, fmt = '(A)') 'ENDVARS'

    close (unit = ctlFileUnit)
  end subroutine generateGrads


  function sstClimatological() result(isExecOk)
    !# sstClimatological
    !# ---
    !# @info
    !# **Brief:** sstClimatological. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    ! 1950-1979 1 Degree x 1 Degree SST
    ! Global NCEP OI Monthly Climatology
    ! Grid Orientation (SSTR):
    ! (1,1) = (0.5_p_r8W,89.5_p_r8N)
    ! (xDim,yDim) = (0.5_p_r8E,89.5_p_r8S)
    integer :: m, monthBefore, monthAfter
    real (kind = p_r8) :: dayHour, dayCorrection, FactorBefore, FactorAfter
    integer :: Header(8)
    real (kind = p_r8) :: sstBefore(var%xDim, var%yDim), sstAfter(var%xDim, var%yDim)
    integer :: inFileUnit
    character (len = maxPathLength) :: sstClimInFilename
    logical :: isExecOk

    isExecOk = .false.
    dayHour = real(day, p_r8) + real(hour, p_r8) / 24.0_p_r8
    monthBefore = month - 1
    if (dayHour > (1.0_p_r8 + real(monthLength(month), p_r8) / 2.0_p_r8)) &
      monthBefore = month
    monthAfter = monthBefore + 1
    if (monthBefore < 1) monthBefore = 12
    if (monthAfter > 12) monthAfter = 1
    dayCorrection = real(monthLength(monthBefore), p_r8) / 2.0_p_r8 - 1.0_p_r8
    if (monthBefore == month) dayCorrection = -dayCorrection - 2.0_p_r8
    FactorAfter = 2.0_p_r8 * (dayHour + dayCorrection) / &
      real(monthLength(monthBefore) + monthLength(monthAfter), p_r8)
    FactorBefore = 1.0_p_r8 - FactorAfter

    write (unit = p_nfprt, fmt = '(/,A)') ' From SSTClimatological:'
    write (unit = p_nfprt, fmt = '(/,A,I4,3(A,I2.2))') &
      ' Year = ', year, ' Month = ', month, &
      ' Day = ', day, ' Hour = ', hour
    write (unit = p_nfprt, fmt = '(/,2(A,I2))') &
      ' MonthBefore = ', monthBefore, ' MonthAfter = ', monthAfter
    write (unit = p_nfprt, fmt = '(/,2(A,F9.6),/)') &
      ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter

    sstClimInFilename =  trim(varCommon%dirClmSST) // fileClmSST
    inFileUnit = openFile(trim(sstClimInFilename), 'formatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return
    do m = 1, 12
      read (unit = inFileUnit, fmt = '(8I5)') Header
      write (unit = p_nfprt, fmt = '(/,1X,9I5,/)') m, Header
      read (unit = inFileUnit, fmt = '(16F5.2)') sstClim
      if (m == monthBefore) then
        sstBefore = sstClim
      end if
      if (m == monthAfter) then
        sstAfter = sstClim
      end if
    end do
    close (unit = inFileUnit)

    ! Linear Interpolation in Time for Year, Month, Day and Hour
    ! of the Initial Condition
    sstClim = FactorBefore * sstBefore + FactorAfter * sstAfter
    isExecOk = .true.
  end function sstClimatological


  subroutine sstClimaWindow ()
    !# sstClimaWindow
    !# ---
    !# @info
    !# **Brief:** sstClimaWindow. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    integer :: j
    real (kind = p_r8) :: lat, dLat

    ! Get Indices to Use CLimatological SST Out of LatClimSouth to LatClimNorth
    js = 0
    jn = 0
    dLat = 2.0_p_r8 * p_lat0 / real(var%yDim - 1, p_r8)
    do j = 1, var%yDim
      lat = p_lat0 - real(j - 1, p_r8) * dLat
      if (lat > var%latClimSouth) js = j
      if (lat > var%latClimNorth) jn = j
    end do
    js = js + 1

    write (unit = p_nfprt, fmt = '(/,A,/)')' From SSTClimaWindow:'
    write (unit = p_nfprt, fmt = '(A,I3,A,F7.3)') &
      ' js = ', js, ' LatClimSouth=', p_lat0 - real(js - 1, p_r8) * dLat
    write (unit = p_nfprt, fmt = '(A,I3,A,F7.3,/)') &
      ' jn = ', jn, ' LatClimNorth=', p_lat0 - real(jn - 1, p_r8) * dLat
  end subroutine sstClimaWindow


  subroutine printNameList()
    !# Prints NameList
    !# ---
    !# @info
    !# **Brief:** Prints NameList. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    write (unit = p_nfprt, fmt = '(/,A)')    ' &SSTWeeklyNameList'
    write (unit = p_nfprt, fmt = '(A,I6)')   '          Mend = ', varCommon%mEnd
    write (unit = p_nfprt, fmt = '(A,I6)')   '          Kmax = ', var%zMax
    write (unit = p_nfprt, fmt = '(A,I6)')   '          Idim = ', var%xDim
    write (unit = p_nfprt, fmt = '(A,I6)')   '          Jdim = ', var%yDim
    write (unit = p_nfprt, fmt = '(A,F6.3)') '     SSTSeaIce = ', var%sstSeaIce
    write (unit = p_nfprt, fmt = '(A,F6.1)') '  LatClimSouth = ', var%latClimSouth
    write (unit = p_nfprt, fmt = '(A,F6.1)') '  LatClimNorth = ', var%latClimNorth
    write (unit = p_nfprt, fmt = '(A,L6)')   '    ClimWindow = ', var%climWindow
    write (unit = p_nfprt, fmt = '(A,L6)')   '        Linear = ', var%linear
    write (unit = p_nfprt, fmt = '(A,L6)')   '    LinearGrid = ', var%linearGrid
    write (unit = p_nfprt, fmt = '(A,L6)')   '         GrADS = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)')      '       DateICn = ' // varCommon%date
    write (unit = p_nfprt, fmt = '(A)')      '       Preffix = ' // var%preffix
    write (unit = p_nfprt, fmt = '(A,/)')    ' /'
  end subroutine printNameList

end module Mod_SSTWeekly

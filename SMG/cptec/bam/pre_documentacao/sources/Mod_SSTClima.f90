!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SSTClima </br></br>
!#
!# **Brief**: Module responsible for correcting SST by topography and generate 
!# the sea ice mask (1: Sea Ice 0: No Ice) </br></br>
!# 
!# The sea surface temperature (SST) data sets used by the model are 
!# interpolated to the model's Gaussian grid. Currently, an interpolation in the
!# area is used separately for SST and sea ice. When a region has an average 
!# area of sea ice greater than or equal 50%, this region is defined as a sea 
!# ice point with a value of 270.2 K (-3ºC). Locations without sea ice have an 
!# interpolated SST value or a value of 271.4 K, whichever is higher. Due to 
!# aliasing Gibbs, SST in locations near the continent will often not be 
!# interpreted as originating at sea level. Instead, the SST will be seen as a 
!# surface for the surface spectral topography, transformed to the Gaussian grid 
!# of the model. This can “raise” the SST to elevations of 3 km or more. To 
!# avoid spurious warming, the SST (excluding sea ice points) is adjusted by the
!# topography of the surface re-transformed by the model using the standard 
!# atmospheric lapse rate of -6.5 K km-1. </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ (= pre/dataout/LandSeaMask.GZZZZZ) ???? </br>
!# &bull; pre/databcs/ersst.form </br>
!# &bull; model/datain/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; model/datain/SSTClimaYYYYMMDD.GZZZZZ (Ex.: model/datain/SSTClima20150430.G00450) </br>
!# &bull; pre/dataout/SSTClimaYYYYMMDD.GZZZZZ </br>
!# &bull; pre/dataout/SSTClimaYYYYMMDD.GZZZZZ.ctl
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
!#  <li>01-08-2007 - Tomita          - version: 1.1.1 </li>
!#  <li>02-07-2019 - Eduardo Khamis  - version: 2.0.1 </li>
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

module Mod_SSTClima

  use Mod_FastFourierTransform, only : createFFT, destroyFFT
  use Mod_LegendreTransform, only : createGaussRep, createSpectralRep, createLegTrans, destroyLegendreObjects
  use Mod_SpectralGrid, only : transp, specCoef2Grid
  use Mod_LinearInterpolation, only : gLatsL => latOut, &
    initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, &
    initAreaInterpolation, doAreaInterpolation
  use Mod_Get_xMax_yMax, only : getxMaxyMax
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  ! Reads the 1x1 SST Global Monthly OI Climatology From NCEP,
  ! Interpolates it Using Area Weigth or Bi-Linear Into a Gaussian Grid

  public :: getNameSSTClima, initSSTClima, generateSSTClima, shouldRunSSTClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !  input parameters

  type SSTClimaNameListData
    integer :: zMax
    integer :: xDim
    integer :: yDim
    real (kind = p_r8) :: sstSeaIce = -1.749
    logical :: linear = .false.
    logical :: linearGrid = .false.
    character (len = 128) :: nameLSM = 'ModelLandSeaMask'
    character (len = 7) :: preffix = 'GANLSMT'
    character (len = 6) :: suffix = 'S.unf.'
  end type SSTClimaNameListData

  type(varCommonNameListData) :: varCommon
  type(SSTClimaNameListData)  :: var
  namelist /SSTClimaNameList/    var

  ! internal variables

  integer :: xmx, yMaxHf, &
    mEnd1, mEnd2, mnwv0, mnwv1, mnwv2, mnwv3
  integer :: year, month, day, hour

  real (kind = p_r8) :: lon0, lat0, p_to, sstOpenWater, &
    sstSeaIceThreshold, lapseRate

  logical :: flagInput(5), flagOutput(5)

  character (len = 10) :: trunc = 'T     L   '
  character (len = 7) :: nLats = '.G     '
  character (len = 10) :: mskfmt = '(      I1)'
  character (len = 8) :: varName = 'SSTClima'
  character (len = 11) :: fileClmSST = 'ersst.form'

  integer :: nficn  ! To read Topography from Initial Condition
  integer :: nflsm  ! To read Formatted Land Sea Mask
  integer :: nfclm  ! To read Formatted Climatological SST
  integer :: nfsto  ! To write Unformatted Gaussian Grid SST
  integer :: nfout  ! To write GrADS Topography, Land Sea, Se Ice and Gauss SST
  integer :: nfctl  ! To write GrADS Control file

  integer :: j, i, m, nr, lRec
  integer :: forecastDay
  real (kind = p_r4) :: timeOfDay
  real (kind = p_r8) :: ggSSTMax, ggSSTMin, mgSSTMax, mgSSTMin
  real (kind = p_r8) :: rgSSTMax, rgSSTMin
  integer :: iCnDate(4), currentDate(4), header(5)
  integer, dimension (:, :), allocatable :: landSeaMask, seaIceMask
  real (kind = p_r4), dimension (:), allocatable :: coefTopIn
  real (kind = p_r4), dimension (:, :), allocatable :: wrOut
  real (kind = p_r8), dimension (:), allocatable :: coefTop
  real (kind = p_r8), dimension (:, :), allocatable :: topog, sstIn, sstGaus, &
    seaIceFlagIn, seaIceFlagOut

  character (len = 6), dimension (:), allocatable :: sstLabel


contains


  function getNameSSTClima() result(returnModuleName)
    !# Returns SSTClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SSTClima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    returnModuleName = "SSTClima"
  end function getNameSSTClima


  function shouldRunSSTClima() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it 
    !# does not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunSSTClima


  function getOutFileName() result(sstClimaOutFilename)
    !# Gets SSTClima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SSTClima Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    character(len = maxPathLength) :: sstClimaOutFilename

    sstClimaOutFilename = trim(varCommon%dirModelIn) // varName // varCommon%date(1:8) // nLats
  end function getOutFileName


  subroutine initSSTClima (nameListFileUnit, varCommon_)
    !# Initializes SSTClima
    !# ---
    !# @info
    !# **Brief:** Initializes SSTClima. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo

    !      integer :: ios
    !
    !
    !      mEnd=62               ! Spectral Resolution Horizontal Truncation
    !      zMax=28               ! Number of Layers of the Initial Condition for the Global Model
    !      xDim=360              ! Number of Longitudes For Climatological CO2 data
    !      yDim=180              ! Number of Latitudes For Climatological CO2 data
    !      sstSeaIce=-1.749_p_r8 ! SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
    !      linear=.true.         ! Flag to Bi-linear (T) or Area (F) Interpolation
    !      linearGrid=.false.    ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
    !      grads=.true.          ! Flag for GrADS Outputs
    !      dateICn='yyyymmddhh'  ! Date of the Initial Condition for the Global Model
    !      preffix='GANLCPT'     ! Preffix of the Initial Condition for the Global Model
    !      suffix='S.unf.'       ! Suffix of the Initial Condition for the Global Model
    !      dirMain='./ '         ! Main Datain/Dataout Directory
    !
    !      open (unit=nfinp, file='./'//nameNML, &
    !            FORM='FORMATTED', ACCESS='SEQUENTIAL', &
    !            ACTION='READ', STATUS='OLD', IOSTAT=ios)
    !      if (ios /= 0) then
    !         write (unit=nferr, fmt='(3A,I4)') &
    !               ' ** (Error) ** open file ', &
    !                 './'//nameNML, &
    !               ' returned iostat = ', ios
    !         stop  ' ** (Error) **'
    !      end if
    !      read  (unit=nfinp, NML=SSTClimaNameList)
    !      close (unit=nfinp)


    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
    
    read(unit = nameListFileUnit, nml = SSTClimaNameList)
    varCommon = varCommon_

    call getxMaxyMax (varCommon%mEnd, varCommon%xMax, varCommon%yMax)
    write (nLats(3:7), '(I5.5)') varCommon%yMax

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

    ! For Linear Interpolation
    lon0 = 0.5_p_r8  ! Start Near Greenwhich
    lat0 = 89.5_p_r8 ! Start Near North Pole

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

    if (var%linearGrid) then
      trunc(2:2) = 'L'
    else
      trunc(2:2) = 'Q'
    end if
    write (trunc(3:6), fmt = '(I4.4)') varCommon%mEnd
    write (trunc(8:10), fmt = '(I3.3)') var%zMax

    write (mskfmt(2:7), '(I6)') varCommon%xMax
  end subroutine initSSTClima


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
    
    write (unit = p_nfprt, fmt = '(/,A)')    ' &SSTClimaNameList'
    write (unit = p_nfprt, fmt = '(A,I6)')   '       mEnd = ', varCommon%mEnd
    write (unit = p_nfprt, fmt = '(A,I6)')   '       zMmax = ', var%zMax
    write (unit = p_nfprt, fmt = '(A,I6)')   '       xDim = ', var%xDim
    write (unit = p_nfprt, fmt = '(A,I6)')   '       yDim = ', var%yDim
    write (unit = p_nfprt, fmt = '(A,F6.3)') '  sstSeaIce = ', var%sstSeaIce
    write (unit = p_nfprt, fmt = '(A,L6)')   '     linear = ', var%linear
    write (unit = p_nfprt, fmt = '(A,L6)')   ' linearGrid = ', var%linearGrid
    write (unit = p_nfprt, fmt = '(A,L6)')   '      grads = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)')      '    dateICn = ' // varCommon%date
    write (unit = p_nfprt, fmt = '(A)')      '    preffix = ' // var%preffix
    write (unit = p_nfprt, fmt = '(A)')      '    suffix = ' // var%suffix
    write (unit = p_nfprt, fmt = '(A)')      '    dirPreOut = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, fmt = '(A)')      '    dirModelIn = ' // trim(varCommon%dirModelIn)
    write (unit = p_nfprt, fmt = '(A)')      '    dirClmSST = ' // trim(varCommon%dirClmSST)
    write (unit = p_nfprt, fmt = '(A,/)')    ' /'
  end subroutine printNameList


  function generateSSTClima() result(isExecOk)
    !# Generates SSTClima
    !# ---
    !# @info
    !# **Brief:** Generates SSTClima. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk

    isExecOk = .false.

    write (unit = p_nfprt, fmt = '(/,A,I6)') '   Imax = ', varCommon%xMax
    write (unit = p_nfprt, fmt = '(A,I6,/)') '   Jmax = ', varCommon%yMax

    call createSpectralRep (mEnd1, mEnd2, mnwv1)
    call createGaussRep (varCommon%yMax, yMaxHf)
    call createFFT (varCommon%xMax)
    call createLegTrans (mnwv0, mnwv1, mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf)
    if (var%linear) then
      call initLinearInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, lat0, lon0)
    else
      call initAreaInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if

    allocate (landSeaMask(varCommon%xMax, varCommon%yMax), seaIceMask(varCommon%xMax, varCommon%yMax))
    allocate (coefTopIn(mnwv2))
    allocate (coefTop(mnwv2), topog(varCommon%xMax, varCommon%yMax))
    allocate (sstIn(var%xDim, var%yDim), sstGaus(varCommon%xMax, varCommon%yMax), wrOut(varCommon%xMax, varCommon%yMax))
    allocate (seaIceFlagIn(var%xDim, var%yDim), seaIceFlagOut(varCommon%xMax, varCommon%yMax))
    allocate (sstLabel(var%yDim))

    ! Read in SpectraL Coefficient of Topography from ICn
    ! to Ensure that Topography is the Same as Used by Model
    nficn = openFile(trim(varCommon%dirModelIn) // var%preffix // varCommon%date // var%suffix // trunc, &
      'unformatted', 'sequential', -1, 'read', 'old') 
    if(nficn < 0) return
    read  (unit = nficn) forecastDay, timeOfDay, iCnDate, currentDate
    read  (unit = nficn) coefTopIn
    close (unit = nficn)
    coefTop = real(coefTopIn, p_r8)
    call transp(mnwv2, mEnd1, mEnd2, coefTop)
    call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, coefTop, topog)

    ! Read in Land-Sea Mask Data Set
    nflsm = openFile(trim(varCommon%dirPreOut) // trim(var%nameLSM) // nLats, 'unformatted', 'sequential', -1, 'read', 'old') 
    if(nflsm < 0) return
    read  (unit = nflsm) landSeaMask
    close (unit = nflsm)

    ! Open File for Land-Sea Mask and SST Data to Global Model Input
    inquire (iolength = lRec) wrOut
    nfsto = openFile(getOutFileName(), 'unformatted', 'direct', lRec, 'write', 'replace')
    if(nfsto < 0) return
    ! Write out Land-Sea Mask to SST Data Set
    ! The LSMask will be Transfered by Model to Post-Processing
    wrOut = real(1 - 2 * landSeaMask, p_r4)
    write (unit = nfsto, rec = 1) wrOut

    if (varCommon%grads) then
      nfout = openFile(trim(varCommon%dirPreOut) // varName // varCommon%date(1:8) // nLats ,&
        'unformatted', 'direct', lRec, 'write', 'replace')
        if(nfout < 0) return
    end if

    nfclm = openFile(trim(varCommon%dirClmSST) // fileClmSST, 'formatted', 'sequential', -1, 'read', 'old') 
    if(nfclm < 0) return

    ! Loop Through Months
    do m = 1, 12
      read (unit = nfclm, fmt = '(8I5)') header
      write (unit = p_nfprt, fmt = '(/,1X,9I5,/)') m, header
      read (unit = nfclm, fmt = '(16F5.2)') sstIn

      if (maxval(sstIn) < 100.0_p_r8) sstIn = sstIn + p_to

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
      mgSSTMax = maxval(sstGaus, mask = sstGaus/=p_undef)
      mgSSTMin = minval(sstGaus, mask = sstGaus/=p_undef)

      ! Write out Gaussian Grid Weekly SST
      wrOut = real(sstGaus, p_r4)
      write (unit = nfsto, rec = m + 1) wrOut

      write (unit = p_nfprt, fmt = '(/,3(A,I2.2),A,I4)') &
        ' Hour = ', hour, ' Day = ', day, &
        ' Month = ', month, ' Year = ', year

      write (unit = p_nfprt, fmt = '(/,A,3(A,2F8.2,/))') &
        ' Mean Weekly SST Interpolation :', &
        ' Regular  Grid SST: min, max = ', rgSSTMin, rgSSTMax, &
        ' Gaussian Grid SST: min, max = ', ggSSTMin, ggSSTMax, &
        ' Masked G Grid SST: min, max = ', mgSSTMin, mgSSTMax

      if (varCommon%grads) then
        nr = 1 + 4 * (m - 1)
        wrOut = real(topog, p_r4)
        write (unit = nfout, rec = nr) wrOut
        wrOut = real(1 - 2 * landSeaMask, p_r4)
        write (unit = nfout, rec = nr + 1) wrOut
        wrOut = real(seaIceMask, p_r4)
        write (unit = nfout, rec = nr + 2) wrOut
        wrOut = real(sstGaus, p_r4)
        write (unit = nfout, rec = nr + 3) wrOut
      end if

      ! End Loop Through Months
    end do
    close (unit = nfclm)
    close (unit = nfsto)
    close (unit = nfout)

    if (varCommon%grads) then
      ! Write GrADS Control File
      nfctl = openFile(trim(varCommon%dirPreOut) // varName // varCommon%date(1:8) // nLats // '.ctl', &
        'formatted', 'sequential', -1, 'write', 'replace')
      if(nfctl < 0) return
      write (unit = nfctl, fmt = '(A)') 'DSET ^' // &
        varName // varCommon%date(1:8) // nLats
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'TITLE Monthly Climatological OI SST on a Gaussian Grid'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
        'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
      write (unit = nfctl, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
      if (var%linear) then
        write (unit = nfctl, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
      else
        write (unit = nfctl, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
      end if
      write (unit = nfctl, fmt = '(A)') 'ZDEF  1 LEVELS 1000'
      write (unit = nfctl, fmt = '(A)') 'TDEF 12 LINEAR JAN2007 1MO'
      write (unit = nfctl, fmt = '(A)') 'VARS  4'
      write (unit = nfctl, fmt = '(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      write (unit = nfctl, fmt = '(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      write (unit = nfctl, fmt = '(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      write (unit = nfctl, fmt = '(A)') 'SSTC  0 99 Climatological SST Topography Corrected [K]'
      write (unit = nfctl, fmt = '(A)') 'ENDVARS'

      close (unit = nfctl)
    end if

    call destroyFFT()
    call destroyLegendreObjects()

    isExecOk = .true.
  end function generateSSTClima


end module Mod_SSTClima

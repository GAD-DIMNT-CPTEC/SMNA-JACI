!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_FLUXCO2Clima </br></br>
!#
!# **Brief**: Module responsible for reads the 1x1 CO2 Global Monthly OI 
!# Climatology From NCEP </br>
!#
!! The FLUXCO2Clima.f90 subroutine reads the NCEP global monthly climatological 
!# CO2 set and interpolates them using either the area weight or bi-linear 
!# interpolation in a Gaussian grid of the model. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ </br>
!# &bull; pre/databcs/FluxCO2.bin </br>
!# &bull; pre/databcs/ersst.form </br>
!# &bull; model/datain/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ (Ex.: model/datain/GANLNMC2015043000S.unf.TQ0299L064)
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/FLUXCO2ClimaYYYYMMDD.GZZZZZ  (Ex.: model/datain/FLUXCO2ClimaYYYYMMDD.G00450) </br>
!# &bull; pre/dataout/FLUXCO2ClimaYYYYMMDD.GZZZZZ.ctl </br>
!# &bull; pre/dataout/FLUXCO2ClimaYYYYMMDD.GZZZZZ 
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
!#  <li>04-06-2019 - Eduardo Khamis  - version: 2.0.1 </li>
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

module Mod_FLUXCO2Clima  ! SSTClima + CO2

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

  public :: getNameFLUXCO2Clima, initFLUXCO2Clima, generateFLUXCO2Clima, shouldRunFLUXCO2Clima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !input variables

  type FluxoCO2ClimaNameListData
    integer :: zMax
    integer :: xDim
    integer :: yDim
    real (kind = p_r8) :: sstSeaIce = -1.749_p_r8
    logical :: linear
    logical :: linearGrid
    character (len = 7) :: preffix
    character (len = 6) :: suffix
    character (len = 128) :: nameLSM = 'LandSeaMask'
  end type FluxoCO2ClimaNameListData

  type(varCommonNameListData)     :: varCommon
  type(FluxoCO2ClimaNameListData) :: var
  namelist /FluxoCO2ClimaNameList/   var

  ! internal variables

  integer :: xmx, yMaxHf, mEnd1, mEnd2, mnwv0, mnwv1, mnwv2, mnwv3
  integer :: year, month, day, hour
  real (kind = p_r8) :: lon0, lat0, p_to, co2SeaIce, co2OpenWater, sstOpenWater, &
    co2SeaIceThreshold, sstSeaIceThreshold, lapseRate
  logical :: flagInput(5), flagOutput(5)
  character (len = 10) :: trunc = 'T     L   '
  character (len = 7) :: nLats = '.G     '
  character (len = 10) :: mskfmt = '(      I1)'
  character (len = 12) :: varName = 'FLUXCO2Clima'
  character (len = 11) :: fileClmFluxCO2 = 'FluxCO2.bin'
  character (len = 11) :: fileClmSST = 'ersst.form'

  integer :: j, i, m, nr, lRec
  integer :: forecastDay
  real (kind = p_r4) :: timeOfDay
  real (kind = p_r8) :: ggCO2Max, ggCO2Min, mgCO2Max, mgCO2Min
  real (kind = p_r8) :: rgSSTMax, rgSSTMin
  integer :: iCnDate(4), currentDate(4), header(5)
  integer, dimension (:, :), allocatable :: landSeaMask, seaIceMask
  real (kind = p_r4), dimension (:), allocatable :: coefTopIn
  real (kind = p_r4), dimension (:, :), allocatable :: wrOut, wrIn
  real (kind = p_r8), dimension (:), allocatable :: coefTop
  real (kind = p_r8), dimension (:, :), allocatable :: topog, co2In, sstIn, co2Gaus, sstGaus, &
    seaIceFlagIn, seaIceFlagOut
  character (len = 6), dimension (:), allocatable :: co2Label
  integer, parameter :: lmon(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)


contains


  function getNameFLUXCO2Clima() result(returnModuleName)
    !# Returns FLUXCO2Clima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns FLUXCO2Clima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "FLUXCO2Clima"
  end function getNameFLUXCO2Clima


  function shouldRunFLUXCO2Clima() result(shouldRun)
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
  end function shouldRunFLUXCO2Clima


  function getOutFileName() result(FLUXCO2ClimaOutFilename)
    !# Gets FLUXCO2Clima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets FLUXCO2Clima Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: FLUXCO2ClimaOutFilename

    FLUXCO2ClimaOutFilename = trim(varCommon%dirModelIn) // varName // varCommon%date(1:8) // nLats
  end function getOutFileName


  subroutine initFLUXCO2Clima(nameListFileUnit, varCommon_)
    !# Initializes FLUXCO2Clima
    !# ---
    !# @info
    !# **Brief:** Initializes FLUXCO2Clima. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
      
    !  integer :: ios </br>
    !  mEnd=62               - Spectral Resolution Horizontal Truncation </br>
    !  zMax=28               - Number of Layers of the Initial Condition for the Global Model </br>
    !  xDim=360              - Number of Longitudes For Climatological CO2 data </br>
    !  yDim=180              - Number of Latitudes For Climatological CO2 data </br>
    !  sstSeaIce=-1.749_p_r8 - SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC) </br>
    !  linear=.true.         - Flag to Bi-linear (T) or Area (F) Interpolation </br>
    !  linearGrid=.false.    - Flag to Set Linear (T) or Quadratic Gaussian Grid (T) </br>
    !  grads=.true.          - Flag for GrADS Outputs </br>
    !  dateICn='yyyymmddhh'  - Date of the Initial Condition for the Global Model </br>
    !  preffix='GANLCPT'     - Preffix of the Initial Condition for the Global Model </br>
    !  suffix='S.unf.'       - Suffix of the Initial Condition for the Global Model </br>
    !  dirMain='./ '         - Main Datain/Dataout Directory </br>
    
    !      open (unit=p_nfinp, file='./'//nameNML, &
    !            FORM='FORMATTED', ACCESS='SEQUENTIAL', &
    !            ACTION='READ', STATUS='OLD', IOSTAT=ios)
    !      if (ios /= 0) then
    !         write (unit=p_nferr, fmt='(3A,I4)') &
    !               ' ** (Error) ** open file ', &
    !                 './'//nameNML, &
    !               ' returned iostat = ', ios
    !         stop  ' ** (Error) **'
    !      end if
    !      read  (unit=p_nfinp, NML=FluxoCO2ClimaNameList)
    !      close (unit=p_nfinp)
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = FluxoCO2ClimaNameList)
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
    co2OpenWater = 0.0_p_r8
    co2SeaIceThreshold = 0.0_p_r8
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

  end subroutine initFLUXCO2Clima


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

    write (unit = p_nfprt, fmt = '(/,A)')    ' &FluxoCO2Clima'
    write (unit = p_nfprt, fmt = '(A,I6)')   '       Mend = ', varCommon%mEnd
    write (unit = p_nfprt, fmt = '(A,I6)')   '       Kmax = ', var%zMax
    write (unit = p_nfprt, fmt = '(A,I6)')   '       Idim = ', var%xDim
    write (unit = p_nfprt, fmt = '(A,I6)')   '       Jdim = ', var%yDim
    write (unit = p_nfprt, fmt = '(A,F6.3)') '  SSTSeaIce = ', var%sstSeaIce
    write (unit = p_nfprt, fmt = '(A,L6)')   '     Linear = ', var%linear
    write (unit = p_nfprt, fmt = '(A,L6)')   ' LinearGrid = ', var%linearGrid
    write (unit = p_nfprt, fmt = '(A,L6)')   '      GrADS = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)')      '    DateICn = ' // varCommon%date
    write (unit = p_nfprt, fmt = '(A)')      '    Preffix = ' // var%preffix
    write (unit = p_nfprt, fmt = '(A)')      '    DirPreOut = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, fmt = '(A)')      '    DirModelIn = ' // trim(varCommon%dirModelIn)
    write (unit = p_nfprt, fmt = '(A)')      '    DirClmFluxCO2 = ' // trim(varCommon%dirClmFluxCO2)
    write (unit = p_nfprt, fmt = '(A)')      '    DirClmSST = ' // trim(varCommon%dirClmSST)
    write (unit = p_nfprt, fmt = '(A)')      '    nameLSM = ' // trim(var%nameLSM)
    write (unit = p_nfprt, fmt = '(A,/)')    ' /'

  end subroutine printNameList


  function generateFLUXCO2Clima() result(isExecOk)
    !# Generates FLUXCO2Clima
    !# ---
    !# @info
    !# **Brief:** Generates FLUXCO2Clima. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: icnSpectralCoefOfTopographyFileName
    character(len = maxPathLength) :: landSeaMaskDataSetFileName
    character(len = maxPathLength) :: fileClmSSTFileName
    character(len = maxPathLength) :: fileClmFluxCO2FileNamme
    character(len = maxPathLength) :: landSeaMaskCO2GradsCtlFileName
    character(len = maxPathLength) :: landSeaMaskCO2GradsFileName
    
    integer :: nficn  ! To read Topography from Initial Condition
    integer :: nflsm  ! To read Formatted Land Sea Mask
    integer :: nfclm  ! To read Formatted Climatological CO2
    integer :: nfclm2 ! To read Formatted Climatological SST
    integer :: nfsto  ! To write Unformatted Gaussian Grid CO2
    integer :: nfout  ! To write GrADS Topography, Land Sea, Se Ice and Gauss CO2
    integer :: nfctl  ! To write GrADS Control file
    logical :: isExecOk

    isExecOk = .false.
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
    allocate (co2In(var%xDim, var%yDim), sstIn(var%xDim, var%yDim), wrIn(var%xDim, var%yDim), co2Gaus(varCommon%xMax, varCommon%yMax), sstGaus(varCommon%xMax, varCommon%yMax), wrOut(varCommon%xMax, varCommon%yMax))
    allocate (seaIceFlagIn(var%xDim, var%yDim), seaIceFlagOut(varCommon%xMax, varCommon%yMax))
    allocate (co2Label(var%yDim))

    ! Read in SpectraL Coefficient of Topography from ICn
    ! to Ensure that Topography is the Same as Used by Model
    icnSpectralCoefOfTopographyFileName = trim(varCommon%dirModelIn) // var%preffix // varCommon%date // var%suffix // trunc
    nficn = openFile(trim(icnSpectralCoefOfTopographyFileName), 'unformatted', 'sequential', -1, 'read', 'old')
    if(nficn < 0) return

    read  (unit = nficn) forecastDay, timeOfDay, iCnDate, currentDate
    read  (unit = nficn) coefTopIn
    close (unit = nficn)
    coefTop = real(coefTopIn, p_r8)
    call transp(mnwv2, mEnd1, mEnd2, coefTop)
    call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, coefTop, topog)

    ! Read in Land-Sea Mask Data Set
    landSeaMaskDataSetFileName = trim(varCommon%dirPreOut) // trim(var%nameLSM) // nLats
    nflsm = openFile(trim(landSeaMaskDataSetFileName), 'unformatted', 'sequential', -1, 'read', 'old')
    if(nflsm < 0) return
    read  (unit = nflsm) landSeaMask
    close (unit = nflsm)

    ! Open File for Land-Sea Mask and CO2 Data to Global Model Input
    inquire (iolength = lRec) wrOut
    nfsto = openFile(trim(getOutFileName()), 'unformatted', 'direct', lRec, 'write', 'replace')
    if(nfsto < 0) return

    ! Write out Land-Sea Mask to CO2 Data Set
    ! The LSMask will be Transfered by Model to Post-Processing
    wrOut = real(1 - 2 * landSeaMask, p_r4)
    write (unit = nfsto, rec = 1) wrOut
    if (varCommon%grads) then
      landSeaMaskCO2GradsFileName = trim(varCommon%dirPreOut) // varName // varCommon%date(1:8) // nLats
      nfout = openFile(trim(landSeaMaskCO2GradsFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
      if(nfout < 0) return
    endif

    fileClmSSTFileName = trim(varCommon%dirClmSST) // fileClmSST
    nfclm2 = openFile(trim(fileClmSSTFileName), 'formatted', 'sequential', -1, 'read', 'old')
    if(nfclm2 < 0) return

    inquire(iolength = lrec) wrIn
    fileClmFluxCO2FileNamme = trim(varCommon%dirClmFluxCO2) // fileClmFluxCO2
    nfclm = openFile(trim(fileClmFluxCO2FileNamme), 'unformatted', 'direct', lrec, 'read', 'old')
    if(nfclm < 0) return

    ! Loop Through Months
    do m = 1, 12
      read (unit = nfclm2, fmt = '(8I5)') header
      write (unit = p_nfprt, fmt = '(/,1X,9I5,/)') m, header
      read (unit = nfclm2, fmt = '(16F5.2)') sstIn

      read (nfclm, rec = m) wrIn
      co2In = wrIn
      if (maxval(sstIn) < 100.0_p_r8) sstIn = sstIn + p_to

      ! Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
      ! Over Open Water Set Input CO2 = MIN of CO2OpenWater
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
      ! Min And Max Values of Input CO2
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

      ! Interpolate CO2 from 1x1 Grid to Gaussian Grid
      if (var%linear) then
        call doLinearInterpolation (co2In, co2Gaus)
      else
        call doAreaInterpolation (co2In, co2Gaus)
      end if
      ! Min and Max Values of Gaussian Grid
      ggCO2Max = maxval(co2Gaus)
      ggCO2Min = minval(co2Gaus)

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

      do j = 1, varCommon%yMax
        do i = 1, varCommon%xMax
          if (landSeaMask(i, j) == 1) then
            ! Set CO2 = Undef Over Land
            co2Gaus(i, j) = p_undef
          else if (seaIceMask(i, j) == 1) then
            ! Set CO2 Sea Ice Threshold Minus 1 Over Sea Ice
            co2Gaus(i, j) = co2SeaIceThreshold
          else
            ! Correct CO2 for Topography, Do Not Create or
            ! Destroy Sea Ice Via Topography Correction
            co2Gaus(i, j) = (co2Gaus(i, j) + 0.03333334_p_r8 + 0.00136_p_r8) / (86400.0_p_r8 * real(lmon(m), kind = p_r8))!-topog(i,j)*lapseRate
            !IF (co2Gaus(i,j) < co2SeaIceThreshold) &
            !   co2Gaus(i,j)=co2SeaIceThreshold+0.2_p_r8
          end if
        end do
      end do
      ! Min and Max Values of Corrected Gaussian Grid CO2 Excluding Land Points
      mgCO2Max = maxval(co2Gaus, mask = co2Gaus/=p_undef)
      mgCO2Min = minval(co2Gaus, mask = co2Gaus/=p_undef)

      ! Write out Gaussian Grid Weekly CO2
      wrOut = real(co2Gaus, p_r4)
      write (unit = nfsto, rec = m + 1) wrOut

      write (unit = p_nfprt, fmt = '(/,3(A,I2.2),A,I4)') &
        ' Hour = ', hour, ' Day = ', day, &
        ' Month = ', month, ' Year = ', year

      write (unit = p_nfprt, fmt = '(/,A,3(A,2F8.2,/))') &
        ' Mean Weekly CO2 Interpolation :', &
        ' Regular  Grid SST: min, max = ', rgSSTMin, rgSSTMax, &
        ' Gaussian Grid CO2: min, max = ', ggCO2Min, ggCO2Max, &
        ' Masked G Grid CO2: min, max = ', mgCO2Min, mgCO2Max

      if (varCommon%grads) then
        nr = 1 + 5 * (m - 1)
        wrOut = real(topog, p_r4)
        write (unit = nfout, rec = nr) wrOut
        wrOut = real(1 - 2 * landSeaMask, p_r4)
        write (unit = nfout, rec = nr + 1) wrOut
        wrOut = real(seaIceMask, p_r4)
        write (unit = nfout, rec = nr + 2) wrOut
        wrOut = real(sstGaus, p_r4)
        write (unit = nfout, rec = nr + 3) wrOut
        wrOut = real(co2Gaus, p_r4)
        write (unit = nfout, rec = nr + 4) wrOut
      end if

      ! End Loop Through Months
    end do
    close (unit = nfclm)
    close (unit = nfsto)
    if (varCommon%grads) then
      close (unit = nfout)

      ! Write GrADS Control File
      landSeaMaskCO2GradsCtlFileName = trim(varCommon%dirPreOut) // varName // varCommon%date(1:8) // nLats // '.ctl'
      nfctl = openFile(trim(landSeaMaskCO2GradsCtlFileName), 'formatted', 'sequential', -1, 'write', 'replace')
      if(nfctl < 0) return
      write (unit = nfctl, fmt = '(A)') 'DSET ^' // &
        varName // varCommon%date(1:8) // nLats
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'TITLE Monthly Climatological OI CO2 on a Gaussian Grid'
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
      write (unit = nfctl, fmt = '(A)') 'VARS  5'
      write (unit = nfctl, fmt = '(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      write (unit = nfctl, fmt = '(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      write (unit = nfctl, fmt = '(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      write (unit = nfctl, fmt = '(A)') 'SSTC  0 99 Climatological SST Topography Corrected [K]'
      write (unit = nfctl, fmt = '(A)') 'CO2C  0 99 Climatological CO2 Topography Corrected [kg/m2/s]'
      write (unit = nfctl, fmt = '(A)') 'ENDVARS'

      close (unit = nfctl)
    end if

    deallocate (landSeaMask, seaIceMask)
    deallocate (coefTopIn)
    deallocate (coefTop, topog)
    deallocate (co2In, sstIn, wrIn, co2Gaus, sstGaus, wrOut)
    deallocate (seaIceFlagIn, seaIceFlagOut)
    deallocate (co2Label)

    call destroyFFT()
    call destroyLegendreObjects()

    isExecOk = .true.
  end function generateFLUXCO2Clima


end module Mod_FLUXCO2Clima

!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_TopographyGradient </br></br>
!#
!# **Brief**: Module responsible for reading the topography that is in the 
!# analysis and calculates the derivatives at x and y. </br></br>
!#
!# Note that the analysis is in the file GANLNMCYYYYMMDDHHS.unf.TQYYYYLZZZ, 
!# which contains the model prognostic variables. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/HPRIME.dat
!# &bull; model/datain/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ (Ex.: GANLNMC2015043000S.unf.TQ0299L064)
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/HPRIME.G00450 </br>
!# &bull; model/datain/TopographyGradientYYYYMMDDHH.GZZZZZ (Ex.: TopographyGradient2015043000.G00450) </br>
!# &bull; pre/dataout/HPRIME.G00450.ctl </br>
!# &bull; pre/dataout/TopographyGradientYYYYMMDDHH.GZZZZZ.ctl (Ex.: TopographyGradient2015043000.G00450.ctl)
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita          - version: 1.1.1 </li>
!#  <li>12-02-2020 - Denis Eiras     - version: 2.0.0 </li>
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

module Mod_TopographyGradient

  use Mod_LinearInterpolation, only : gLatsL => LatOut, &
    initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, &
    initAreaInterpolation, doAreaInterpolation
  use Mod_FastFourierTransform, only : createFFT, destroyFFT
  use Mod_LegendreTransform, only : createGaussRep, createSpectralRep, &
    createLegTrans, gLats, destroyLegendreObjects, emRad1
  use Mod_SpectralGrid, only : transp, specCoef2Grid, specCoef2GridD
  use Mod_Get_xMax_yMax, only : getxMaxyMax
  use Mod_FileManager, only: &
    getFileUnit,  &
    openFile, &
    fileExists
  use Mod_Messages, only: &
    msgOut
  
  use Mod_Namelist, only : varCommonNameListData

  implicit none
  private
  include 'pre.h'
  include 'files.h'
  include 'precision.h'

  public :: getNameTopographyGradient, initTopographyGradient
  public :: generateTopographyGradient, shouldRunTopographyGradient

  ! inputparameters

  type TopographyGradientNameListData
    integer :: zMax
    integer :: xDim
    integer :: yDim
    logical :: linear = .true.
    character (len = maxPathLength) :: preffix = 'SMT' 
    !# Preffix for Initial Condition File Name
    character(len = 1) :: grid = 'Q'                   
    !# Grid Type: Q for Quadratic and L for linear Gaussian
  end type TopographyGradientNameListData

  type(varCommonNameListData)          :: varCommon
  type(TopographyGradientNameListData) :: var
  namelist /TopographyGradientNameList/   var

  ! internal variables

  real(kind = p_r8) :: lon0, lat0
  logical :: flagInput(5), flagOutput(5)
  integer :: mEnd1, mEnd2, mEnd3
  integer :: mnwv2, mnwv3, mnwv0, mnwv1
  integer :: xmx, yMaxHf
  integer :: mFactorFourier, mTrigsFourier
  integer(kind = p_r4) :: forecastDay
  real(kind = p_r4) :: timeOfDay
  real(kind = p_r8) :: p_undef_topo, rad, dLon
  character(len = 12) :: tGrads
  character(len = 33) :: fileInp
  character(len = 35) :: fileOut
  character(len = 3), dimension(12) :: months
  integer, dimension(4) :: inputDate
  character(len = 7) :: nLats = '.G     '
  character(len = 10) :: varNameT = 'HPRIME'
  character(len = 10) :: varName = 'HPRIME'
  character(len = *), parameter :: headerMsg = 'Topography Gradient   | '

  ! inputarrays

  integer(kind = p_r4), dimension(:), allocatable :: initialDate, currentDate
  real(kind = p_r4), dimension(:), allocatable :: sigmaInteface, sigmaLayer, qTopoInp
  real(kind = p_r4), dimension(:, :), allocatable :: gradsOut
  real(kind = p_r8), dimension(:), allocatable :: qTopo, qTopoS, cosLatInv
  real(kind = p_r8), dimension(:, :), allocatable :: topo, dTopoDx, dTopoDy

  ! local

  integer :: i, j, ll, mm, nn, lRec, lRecIn, lRecOut
  integer, parameter :: nVar = 14

contains


  function getNameTopographyGradient() result(returnModuleName)
    !# Returns TopographyGradient Module Name
    !# ---
    !# @info
    !# **Brief:** Returns TopographyGradient Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    returnModuleName = "TopographyGradient"
  end function getNameTopographyGradient


  subroutine initTopographyGradient(nameListFileUnit, varCommon_)
    !# Initialization of TopographyGradient module
    !# ---
    !# @info
    !# **Brief:** Initialization of TopographyGradient module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = TopographyGradientNameList)
    varCommon = varCommon_

    months = (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
      'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

    ! Files Units

    !    p_nferr = 0    ! Standard Error Print Out
    !    p_nfinpnml = 5    ! Standard Read In

    !    mEnd = 62  ! Spectral Horizontal Resolution of Output Data
    !    var%zMax = 28  ! Number of Levels of Data
    !    xDim = 2161 ! Number of Longitudes in Navy Topography Data
    !    yDim = 1081 ! Number of Latitudes in Navy Topography Data
    !    grads = .true. ! Flag for GrADS Outputs
    !    linear = .true. ! Flag for linear (T) or Area Weighted (F) Interpolation
    !    iDate = '2004032600'! Date of Initial Condition
    !    preffix = 'NMC' ! Preffix for Initial Condition File Name
    !    grid = 'Q' ! Grid Type: Q for Quadratic and L for linear Gaussian

    ! xMax - Number of Longitudes of Output Gaussian Grid Data
    ! yMax - Number of Latitudes of Output Gaussian Grid Data

    call getxMaxyMax (varCommon%mEnd, varCommon%xMax, varCommon%yMax)

    write (nLats(3:7), fmt = '(I5.5)') varCommon%yMax
    call fillFileInp()
    call fillFileOut()

    mEnd1 = varCommon%mEnd + 1
    mEnd2 = varCommon%mEnd + 2
    mEnd3 = varCommon%mEnd + 3
    mnwv2 = mEnd1 * mEnd2
    mnwv0 = mnwv2 / 2
    mnwv3 = mnwv2 + 2 * mEnd1
    mnwv1 = mnwv3 / 2

    xmx = varCommon%xMax + 2
    yMaxHf = varCommon%yMax / 2

    mFactorFourier = 64
    mTrigsFourier = 3 * varCommon%xMax / 2

    rad = atan(1.0_p_r8) / 45.0_p_r8
    dLon = 360.0_p_r8 / real(varCommon%xMax, p_r8)

    p_undef_topo = -99999.0_p_r8

    ! For linear Interpolation
    lon0 = 0.0_p_r8  ! Start at Prime Meridian
    lat0 = 90.0_p_r8 - 0.5_p_r8 * (360.0_p_r8 / real(var%xDim, p_r8)) ! Start at North Pole

  end subroutine initTopographyGradient


  function shouldRunTopographyGradient() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunTopographyGradient


  subroutine fillFileInp()
    !# Fills Input Data File Name
    !# ---
    !# @info
    !# **Brief:** Fills Input Data File Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none

    fileInp = 'FileInp.dat '     ! Name of File with Input Data
    read(varCommon%date(1:4), '(I4)')inputDate(4)
    read(varCommon%date(5:6), '(I2)')inputDate(3)
    read(varCommon%date(7:8), '(I2)')inputDate(2)
    read(varCommon%date(9:10), '(I2)')inputDate(1)

    fileInp = 'GANL             S.unf.T     L   '
    fileInp(5:7) = trim(var%preffix)
    write(fileInp(8:11), fmt = '(I4.4)') inputDate(4)
    write(fileInp(12:13), fmt = '(I2.2)') inputDate(3)
    write(fileInp(14:15), fmt = '(I2.2)') inputDate(2)
    write(fileInp(16:17), fmt = '(I2.2)') inputDate(1)
    fileInp(25:25) = var%grid
    write(fileInp(26:29), fmt = '(I4.4)') varCommon%mEnd
    write(fileInp(31:33), fmt = '(I3.3)') var%zMax
  end subroutine fillFileInp


  subroutine fillFileOut()
    !# Fills Output Data File Name
    !# ---
    !# @info
    !# **Brief:** Fills Output Data File Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none

    call fillFileInp()
    fileOut = 'TopographyGradient          .G     '
    fileOut(19:28) = fileInp(8:17)
    call getxMaxyMax (varCommon%mEnd, varCommon%xMax, varCommon%yMax)
    write(fileOut(31:35), fmt = '(I5.5)') varCommon%yMax
  end subroutine fillFileOut


  function getOutFileName() result(topographyGradientOutFilename)
    !# Gets TopographyGradient Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets TopographyGradient Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: topographyGradientOutFilename

    call fillFileInp()
    call fillFileOut()

    topographyGradientOutFilename = trim(varCommon%dirModelIn) // trim(fileOut)
    !model/datain/HPRIME.G00450
    !pre/dataout/HPRIME.G00450.ctl
    !model/datain/TopographyGradient2015043000.G00450
    !pre/dataout/TopographyGradient2015043000.G00450.ctl
  end function getOutFileName


  subroutine getArrays
    !# Gets Arrays
    !# ---
    !# @info
    !# **Brief:** Gets Arrays. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none

    allocate(initialDate(4), currentDate(4))
    allocate(gradsOut(varCommon%xMax, varCommon%yMax))
    allocate(sigmaInteface(var%zMax + 1), sigmaLayer(var%zMax), qTopoInp(mnwv2))
    allocate(qTopo(mnwv2), qTopoS(mnwv2), cosLatInv(varCommon%yMax))
    allocate(topo(varCommon%xMax, varCommon%yMax), dTopoDx(varCommon%xMax, varCommon%yMax), dTopoDy(varCommon%xMax, varCommon%yMax))
  end subroutine getArrays


  subroutine clsArrays
    !# Cleans Arrays
    !# ---
    !# @info
    !# **Brief:** Cleans Arrays. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none

    deallocate(initialDate, currentDate)
    deallocate(gradsOut)
    deallocate(sigmaInteface, sigmaLayer, qTopoInp)
    deallocate(qTopo, qTopoS, cosLatInv)
    deallocate(topo, dTopoDx, dTopoDy)
  end subroutine clsArrays


  function readwriteHPRIME() result(isExecOk)
    !# Reads and writes HPRIME
    !# ---
    !# @info
    !# **Brief:** Reads and writes HPRIMECleans Arrays. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    integer :: inFileUnit, outFileUnit, inRecSize, outRecSize
    character (len = maxPathLength) :: hPrimeInFilename, fieldOutFilename
    logical :: isExecOk

    !    integer :: p_nfclm = 10   ! To Read Topography Data (GTOP30 or Navy)
    !    integer :: p_nfoub = 30   ! To Write Intepolated Topography Data

    real(kind = p_r4), dimension(:, :), allocatable :: hPrimeIn
    real(kind = p_r8), dimension(:, :), allocatable :: hPrimeInput
    real(kind = p_r8), dimension(:, :), allocatable :: hPrimeOutput
    real(kind = p_r4), dimension(:, :), allocatable :: fieldOut

    isExecOk = .false.
    if (var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, lat0, lon0)
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

    allocate(hPrimeIn(var%xDim, var%yDim))
    allocate(hPrimeInput(var%xDim, var%yDim))
    allocate(hPrimeOutput(varCommon%xMax, varCommon%yMax))
    allocate(fieldOut(varCommon%xMax, varCommon%yMax))

    ! Read In Input Topo DATA
    inquire(ioLength = inRecSize) hPrimeIn(:, :)
    hPrimeInFilename = trim(varCommon%dirPreOut) // trim(varName) // '.dat'
    inFileUnit = openFile(trim(hPrimeInFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return
    ! write In Input Topo DATA
    inquire(ioLength = outRecSize) fieldOut(:, :)
    fieldOutFilename = trim(varCommon%dirModelIn) // trim(varNameT) // nLats
    outFileUnit = openFile(trim(fieldOutFilename), 'unformatted', 'direct', outRecSize, 'write', 'replace')
    if(outFileUnit < 0) return

    do j = 1, nVar
      read(unit = inFileUnit, REC = j) hPrimeIn(:, :)
      hPrimeInput = real(hPrimeIn, p_r8)
      ! Interpolate Input Regular Grid Albedo To Gaussian Grid on Output
      if (var%linear) then
        call doLinearInterpolation (hPrimeInput(:, :), hPrimeOutput(:, :))
      else
        call doAreaInterpolation   (hPrimeInput(:, :), hPrimeOutput(:, :))
      end if
      fieldOut = real(hPrimeOutput, p_r4)
      write(unit = outFileUnit, REC = j) fieldOut
    end do

    deallocate(hPrimeIn)
    deallocate(hPrimeInput)
    deallocate(hPrimeOutput)
    deallocate(fieldOut)

    close(unit = inFileUnit)
    close(unit = outFileUnit)

    if (varCommon%grads) then
      call generateGradsHPRIME
    end if
    isExecOk = .true.
  end function readwriteHPRIME


  subroutine generateGradsHPRIME()
    !# Generates Grads Control File to the HPRIME
    !# ---
    !# @info
    !# **Brief:** Generates Grads HPRIME. </br>
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
    gradsBaseName = trim(varNameT) // nLats
    ctlFileName = trim(gradsBaseName) // '.ctl'
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(ctlFileName)
    ctlFileUnit = openFile(ctlPathFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if(ctlFileUnit < 0) return
    write (unit = ctlFileUnit, fmt = '(A)') 'DSET ^' // trim(varNameT) // nLats
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef_topo
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'TITLE Topography and Variance on a Gaussian Grid'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)')'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write (unit = ctlFileUnit, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
    if (var%linear) then
      write (unit = ctlFileUnit, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write (unit = ctlFileUnit, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
    end if
    write (unit = ctlFileUnit, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
    write (unit = ctlFileUnit, fmt = '(A)') 'TDEF 1 LINEAR JAN2005 1MO'
    write (unit = ctlFileUnit, fmt = '(A,i5)')'VARS ', nVar
    write (unit = ctlFileUnit, fmt = '(A)')'HSTDV  0 99 standard deviation of orography'
    write (unit = ctlFileUnit, fmt = '(A)')'HCNVX  0 99 Normalized convexity'
    write (unit = ctlFileUnit, fmt = '(A)')'HASYW  0 99 orographic asymmetry in W-E plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HASYS  0 99 orographic asymmetry in S-N plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HASYSW 0 99 orographic asymmetry in SW-NE plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HASYNW 0 99 orographic asymmetry in NW-SE plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HLENW  0 99 orographic length scale in W-E plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HLENS  0 99 orographic length scale in S-N plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HLENSW 0 99 orographic length scale in SW-NE plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HLENNW 0 99 orographic length scale in NW-SE plane'
    write (unit = ctlFileUnit, fmt = '(A)')'HANGL  0 99 angle of the mountain range w/r/t east'
    write (unit = ctlFileUnit, fmt = '(A)')'HSLOP  0 99 slope of orography'
    write (unit = ctlFileUnit, fmt = '(A)')'HANIS  0 99 anisotropy/aspect ratio'
    write (unit = ctlFileUnit, fmt = '(A)')'HZMAX  0 99 max height above mean orography'
    write (unit = ctlFileUnit, fmt = '(A)')'ENDVARS'
    close (unit = ctlFileUnit)
  end subroutine generateGradsHPRIME


  function generateTopographyGradient() result(isExecOk)
    !# Generates Topography Gradient
    !# ---
    !# @info
    !# **Brief:** Generates Topography Gradient. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
  
    implicit none
    !    integer                         :: inFileUnit, outFileUnit, inRecSize, outRecSize
    character (len = maxPathLength) :: fileInpInFilename !, fileOutFilename
    integer :: p_nfinp, p_nfout, p_nfGradsOut
    logical :: isExecOk

    character (len = maxPathLength) :: gradsBaseName, gradsDatPathName

    isExecOk = .false.

    if(.not. readwriteHPRIME()) return
    call getArrays()
    call createSpectralRep(mEnd1, mEnd2, mnwv1)
    call createGaussRep(varCommon%yMax, yMaxHf)
    call createFFT(varCommon%xMax)
    call createLegTrans(mnwv0, mnwv1, mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf)

    ! Write Out Adjusted Interpolated Topography data

    cosLatInV = 1.0_p_r8 / cos(rad * gLats)
    
    fileInpInFilename = trim(varCommon%dirModelIn) // trim(fileInp)
    call msgOut(headerMsg, "Opening as input: " // trim(fileInpInFilename))
    p_nfinp = openFile(fileInpInFilename, 'unformatted', 'sequential', -1, 'read', 'old')
    if(p_nfinp < 0) return

    inquire (ioLength = lRec) gradsOut
    p_nfout = openFile(getOutFileName(), 'unformatted', 'direct', lRec, 'write', 'replace')
    if(p_nfout < 0) return
    call msgOut(headerMsg, "Opening as output: " // trim(getOutFileName()))

    if(varCommon%grads) then
      gradsBaseName = trim(fileOut)
      gradsDatPathName = trim(varCommon%dirPreOut) // trim(gradsBaseName)
      p_nfGradsOut = openFile(trim(gradsDatPathName) , 'unformatted', 'direct', lRec, 'write', 'replace')
      if(p_nfGradsOut < 0) return
      call msgOut(headerMsg, "Opening as output: " // trim(gradsDatPathName))
    endif
    ! Read Topography Coeficients From Model Spectral Initial Condition

    read (unit = p_nfinp) forecastDay, timeOfDay, initialDate, currentDate, &
      sigmaInteface, sigmaLayer
    write (*, *) forecastDay, timeOfDay
    write (*, *) initialDate
    write (*, *) currentDate
    read (unit = p_nfinp) qTopoInp
    qTopo = real(qTopoInp, p_r8)
    call transp(mnwv2, mEnd1, mEnd2, qTopo)
    call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, qTopo, topo)
    do i = 1, varCommon%xMax
      gradsOut(i, varCommon%yMax:1:-1) = real(topo(i, 1:varCommon%yMax), p_r4)
    end do
    write (*, *) minval(gradsout)
    write (*, *) maxval(gradsout)
    write (unit = p_nfout, REC = 1) gradsOut
    if(varCommon%grads) then
      write (unit = p_nfGradsOut, REC = 1) gradsOut
    endif


    ! Zonal Gradient off Topography
    ll = 0
    do mm = 1, mEnd1
      do nn = 1, mEnd2 - mm
        ll = ll + 1
        qTopoS(2 * ll - 1) = -real(nn - 1, p_r8) * qTopo(2 * ll)
        qTopoS(2 * ll) = +real(nn - 1, p_r8) * qTopo(2 * ll - 1)
      end do
    end do
    call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, qTopoS, dTopoDx)
    do i = 1, varCommon%xMax
      gradsOut(i, varCommon%yMax:1:-1) = real(dTopoDx(i, 1:varCommon%yMax) * cosLatInV(1:varCommon%yMax) * emRad1, p_r4)
    end do
    write (unit = p_nfout, REC = 2) gradsOut
    if(varCommon%grads) then
      write (unit = p_nfGradsOut, REC = 2) gradsOut
    endif


    ! Meridional Gradient off Topography
    call specCoef2GridD (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, qTopo, dTopoDy)
    do i = 1, varCommon%xMax
      gradsOut(i, varCommon%yMax:1:-1) = real(dTopoDy(i, 1:varCommon%yMax) * cosLatInV(1:varCommon%yMax) * emRad1, p_r4)
    end do
    write (unit = p_nfout, REC = 3) gradsOut
    if(varCommon%grads) then
      write (unit = p_nfGradsOut, REC = 3) gradsOut
    endif

    close(unit = p_nfinp)
    close(unit = p_nfout)
    if(varCommon%grads) then
      close(unit = p_nfGradsOut)
      call generateGradsCtlTopographyGradient()
    endif

    

    call clsArrays()
    call destroyFFT()
    call destroyLegendreObjects()

    isExecOk = .true.
  end function generateTopographyGradient


  subroutine generateGradsCtlTopographyGradient()
    !# Generates GrADS Control file to the Topography Gradient
    !# ---
    !# @info
    !# **Brief:** Generates GrADS Control file to the Topography Gradient. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxPathLength) :: ctlPathFileName
    character (len = maxPathLength) :: gradsBaseName
    character (len = maxPathLength) :: ctlFileName
    integer :: ctlFileUnit      ! Temp File Unit

    tGrads = '  Z         '
    write (tGrads(1:2), fmt = '(I2.2)') initialDate(1)
    write (tGrads(4:5), fmt = '(I2.2)') initialDate(3)
    tGrads(6:8) = months(initialDate(2))
    write (tGrads(9:12), fmt = '(I4.4)') initialDate(4)

    ! Write GrADS Control File
    gradsBaseName = trim(fileOut)

    ctlFileName = trim(gradsBaseName) // '.ctl'
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(ctlFileName)
    ctlFileUnit = openFile(ctlPathFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if(ctlFileUnit < 0) return

    write (unit = ctlFileUnit, fmt = '(A)') 'DSET ^' // trim(fileOut)
    write (unit = ctlFileUnit, fmt = '(A)') 'options big_endian'
    write (unit = ctlFileUnit, fmt = '(A,1P,G15.7)') 'undef ', p_undef_topo
    write (unit = ctlFileUnit, fmt = '(A)') 'title Topography and its Gradient'
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F10.5)') 'xdef ', varCommon%xMax, ' linear    0.0 ', dLon
    write (unit = ctlFileUnit, fmt = '(A,I5,A,5F10.5)') 'ydef ', varCommon%yMax, ' levels ', gLats(varCommon%yMax:varCommon%yMax - 4:-1)
    write (unit = ctlFileUnit, fmt = '(18X,5F10.5)') gLats(varCommon%yMax - 5:1:-1)
    write (unit = ctlFileUnit, fmt = '(A)') 'zdef 1 levels 1000'
    write (unit = ctlFileUnit, fmt = '(3A)') 'tdef 1 linear ', tGrads, ' 06hr'
    write (unit = ctlFileUnit, fmt = '(A)') 'vars 3'
    write (unit = ctlFileUnit, fmt = '(A)') 'topo 1 99 Topography (m)'
    write (unit = ctlFileUnit, fmt = '(A)') 'dtpx 1 99 Zonal Gradient of Topography (m/m)'
    write (unit = ctlFileUnit, fmt = '(A)') 'dtpy 1 99 Meridional Gradient of Topography (m/m)'
    write (unit = ctlFileUnit, fmt = '(A)') 'endvars'
    close (unit = ctlFileUnit)
  end subroutine generateGradsCtlTopographyGradient


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
    integer :: p_nfprt = 6    ! Standard print Out

    write (unit = p_nfprt, fmt = '(A)') ' '
    write (unit = p_nfprt, fmt = '(A,I5)') '  Imax = ', varCommon%xMax
    write (unit = p_nfprt, fmt = '(A,I5)') '  Jmax = ', varCommon%yMax
    write (unit = p_nfprt, fmt = '(A,I6)') '  Idim = ', var%xDim
    write (unit = p_nfprt, fmt = '(A,I6)') '  Jdim = ', var%yDim
    write (unit = p_nfprt, fmt = '(A,L6)') '  GrADS = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A,L6)') '  Linear = ', var%linear
    write (unit = p_nfprt, fmt = '(/,A)')  '  Numerical Precision (KIND): '
    write (unit = p_nfprt, fmt = '(A,I3)') '  r4  = ', p_r4
    write (unit = p_nfprt, fmt = '(A,I3)') '  r8  = ', p_r8
    write (unit = p_nfprt, fmt = '(A)') ' '

    write (unit = p_nfprt, fmt = '(/,2A)')  '  Input : ', fileInp
    write (unit = p_nfprt, fmt = '(2A)')    '  Output: ', fileOut
    write (unit = p_nfprt, fmt = '(A)') ' '
  end subroutine printNameList

end module Mod_TopographyGradient

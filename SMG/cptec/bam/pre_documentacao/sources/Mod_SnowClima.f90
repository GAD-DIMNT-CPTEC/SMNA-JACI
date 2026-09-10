!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SnowClima </br></br>
!#
!# **Brief**: Module responsible for generating the informations about snow
!# accumulated</br></br>
!# 
!# The Mod_SnowClima.f90 subroutine contains the algorithim that will process 
!# the informations about snow accumulated on the surface from Albedo global 
!# distribuition informations. The algorithim reads the albedo informations and
!# land/sea mask informations in the model resolution and does the fields 
!# composition, creating the thickness snow accumulated.</br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/dataout/Albedo.GZZZZZ  (Ex.: pre/dataout/Albedo.G00450) </br>
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/SnowYYYYMMDDHHS.unf.GZZZZZ (Ex.: model/datain/Snow2015043000S.unf.G00450) </br>
!# &bull; pre/dataout/SnowYYYYMMDDHHS.unf.GZZZZZ.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti   - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita            - version: 1.1.1 </li>
!#  <li>07-05-2019 - Eduardo Khamis    - version: 2.0.0 </li>
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

module Mod_SnowClima

  use Mod_FileManager
  use Mod_Namelist, only : &
    varCommonNameListData

  use Mod_Messages, only : &
    msgOut &
    , msgInLineFormatOut &
    , msgNewLine

  implicit none

  public :: getNameSnowClima, initSnowClima, generateSnowClima, shouldRunSnowClima

  private
  include 'pre.h'
  include 'files.h'
  include 'precision.h'
  include 'messages.h'

  ! input parameters

  type SnowClimaNameListData
    character (len = 128) :: nameLSM = 'LandSeaMask' 
    !# Land Sea mask input file name
  end type SnowClimaNameListData

  type(varCommonNameListData) :: varCommon
  type(SnowClimaNameListData) :: var
  namelist /SnowClimaNameList/   var

  ! internal variables

  integer :: monthBefore, monthAfter

  real (kind = p_r8) :: factorA, factorB

  logical :: icePoints = .false.

  character (len = 5) :: exts = 'S.unf'
  character (len = 9) :: varName = 'SnowClima'
  character (len = 6) :: nameAlb = 'Albedo'
  character (len = 4) :: varNameS = 'Snow'
  character (len = 7) :: nLats = '.G     '
  character (len = 10) :: mskfmt = '(      I1)'
  character (len = 12) :: timeGrADS = '  Z         '

  integer :: i, j
  integer, dimension (:, :), allocatable :: landSeaMask

  real (kind = p_r8), dimension (:), allocatable :: gLats
  real (kind = p_r4), dimension (:, :), allocatable :: snowModel
  real (kind = p_r4), dimension (:, :), allocatable :: albedoBefore
  real (kind = p_r4), dimension (:, :), allocatable :: albedoAfter
  real (kind = p_r8), dimension (:, :), allocatable :: albedo
  real (kind = p_r8), dimension (:, :), allocatable :: snowOut
  integer :: idate(4)

  character(len=*), parameter :: header = 'Snow Clima            | '

contains


  function getNameSnowClima() result(returnModuleName)
    !# Returns SnowClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SnowClima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "SnowClima"
  end function getNameSnowClima


  function shouldRunSnowClima() result(shouldRun)
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
  end function shouldRunSnowClima


  function getOutFileName() result(snowClimaOutFilename)
    !# Gets SnowClima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SnowClima Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo  
    implicit none
    character(len = maxPathLength) :: snowClimaOutFilename

    snowClimaOutFilename = trim(varCommon%dirModelIn) // varNameS // varCommon%date // exts // nLats
  end function getOutFileName


  subroutine initSnowClima (nameListFileUnit, varCommon_)
    !# Initialization of SnowClima module
    !# ---
    !# @info
    !# **Brief:** Initialization of SnowClima module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo    
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    integer :: mon, iadd
    integer, dimension (12) :: monthLength = &
      (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    character (len = 3), dimension (12) :: monthChar = (/ &
      'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
        'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

    !    xMax = 192
    !    yMax = 96
    !    date = '          '
    !    grads = .true.
    !    dirMain = './ '

    read(unit = nameListFileUnit, nml = SnowClimaNameList)
    varCommon = varCommon_

    write (nLats(3:7), '(I5.5)') varCommon%yMax
    write (mskfmt(2:7), '(I6)')  varCommon%xMax

    ! Getting Date
    read (varCommon%date(1:4), '(I4)') idate(4)
    read (varCommon%date(5:6), '(I2)') idate(2)
    read (varCommon%date(7:8), '(I2)') idate(3)
    read (varCommon%date(9:10), '(I2)') idate(1)

    write (timeGrADS(1:2), '(I2.2)') idate(1)
    write (timeGrADS(4:5), '(I2.2)') idate(3)
    write (timeGrADS(6:8), '(A3)') monthChar(idate(2))
    write (timeGrADS(9:12), '(I4.4)') idate(4)

    ! Linear Time Interpolation Factors A and B
    mon = idate(2)
    if (mod(idate(4), 4) == 0) monthLength(2) = 29
    monthBefore = mon - 1
    if (idate(3) > monthLength(mon) / 2) monthBefore = mon
    monthAfter = monthBefore + 1
    if (monthBefore < 1) monthBefore = 12
    if (monthAfter > 12) monthAfter = 1
    iadd = monthLength(monthBefore) / 2
    if (monthBefore == mon) iadd = -iadd
    factorB = 2.0_p_r8 * real(idate(3) + iadd, p_r8) / &
      real(monthLength(monthBefore) + monthLength(monthAfter), p_r8)
    factorA = 1.0_p_r8 - factorB

  end subroutine initSnowClima


  function generateSnowClima() result(isExecOk)
    !# Generates SnowClima
    !# ---
    !# @info
    !# **Brief:** Generates SnowClima. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo   
    implicit none
    integer :: inFileUnit, outFileUnit, lRecIn, lRecOut
    character (len = maxPathLength) :: albedoBeforeInFilename, landSeaMaskInFilename
    logical :: isExecOk

    isExecOk = .false.

    call msgInLineFormatOut(' ' // header, '(A)')
    call msgInLineFormatOut(' Climatological Snow for idate : ', '(A)')
    call msgInLineFormatOut(idate, '(3I3,I5)')
    call msgInLineFormatOut(' MonthBefore = ', '(A)')
    call msgInLineFormatOut(monthBefore, '(I3)')
    call msgInLineFormatOut( ' MonthAfter = ', '(A)')
    call msgInLineFormatOut(monthAfter, '(I3)')
    call msgInLineFormatOut(' FactorA = ', '(A)')
    call msgInLineFormatOut(factorA, '(F8.5)')
    call msgInLineFormatOut(' FactorB = ', '(A)')
    call msgInLineFormatOut(factorB, '(F8.5)')
    call msgNewLine()
   
    call allocateData()

    landSeaMaskInFilename = trim(varCommon%dirPreOut) // trim(var%nameLSM) // nLats
    inFileUnit = openFile(trim(landSeaMaskInFilename), 'unformatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return
    ! Land Sea Mask : Input
    read  (unit = inFileUnit) landSeaMask
    close (unit = inFileUnit)

    ! Climatological Albedo at Model Grid : Input
    inquire (iolength = lRecIn) albedoBefore
    albedoBeforeInFilename = trim(varCommon%dirPreOut) // nameAlb // nLats
    inFileUnit = openFile(trim(albedoBeforeInFilename), 'unformatted', 'direct', lRecIn, 'read', 'old')
    if(inFileUnit < 0) return

    ! Land Sea Mask : Input
    read  (unit = inFileUnit, rec = monthBefore) albedoBefore
    read  (unit = inFileUnit, rec = monthAfter) albedoAfter
    close (unit = inFileUnit)

    ! Computes Climatological Snow Based on Climatological Albedo
    do j = 1, varCommon%yMax
      do i = 1, varCommon%xMax
        snowOut(i, j) = 0.0_p_r8
        ! Linear Interpolation of Albedo in Time
        albedo(i, j) = factorA * real(albedoBefore(i, j), p_r8) + &
          factorB * real(albedoAfter(i, j), p_r8)
        if (landSeaMask(i, j) /= 0) then
          if (.not.icePoints) then
            if (albedo(i, j) >= 0.40_p_r8 .and. albedo(i, j) < 0.49_p_r8) then
              snowOut(i, j) = 5.0_p_r8
            end if
            if (albedo(i, j) >= 0.49_p_r8 .and. albedo(i, j) < 0.69_p_r8) then
              snowOut(i, j) = 10.0_p_r8
            end if
            if (albedo(i, j) >= 0.69_p_r8 .and. albedo(i, j) < 0.75_p_r8) then
              snowOut(i, j) = 20.0_p_r8
            end if
          end if
          if (albedo(i, j) >= 0.75_p_r8) then
            snowOut(i, j) = 3000.0_p_r8
          end if
        end if
      end do
    end do

    ! IEEE-32 Bits Output
    snowModel = real(snowOut, p_r4)

    ! Climatological Snow at Model Grid and
    ! Time of Initial Condititon : Model Input
    inquire (iolength = lRecOut) snowModel
    outFileUnit = openFile(getOutFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(outFileUnit < 0) return
    write (unit = outFileUnit, REC = 1) snowModel
    close (unit = outFileUnit)

    if (varCommon%grads) then
      call generateGrads()
    end if

    call deallocateData()

    isExecOk = .true.
  end function generateSnowClima


  subroutine gaussianLatitudes ()
    !# Gaussian Latitudes
    !# ---
    !# @info
    !# **Brief:** Gaussian Latitudes. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo  
    implicit none
    integer :: yMaxH, j
    real (kind = p_r8) :: eps, rad, dGcolIn, gcol, dGcol, p1, p2

    eps = epsilon(1.0_p_r8) * 100.0_p_r8
    yMaxH = varCommon%yMax / 2
    dGcolIn = atan(1.0_p_r8) / real(varCommon%yMax, p_r8)
    rad = 45.0_p_r8 / atan(1.0_p_r8)
    gcol = 0.0_p_r8
    do j = 1, yMaxH
      dGcol = dGcolIn
      do
        call legendrePolynomial (varCommon%yMax, gcol, p2)
        do
          p1 = p2
          gcol = gcol + dGcol
          call legendrePolynomial (varCommon%yMax, gcol, p2)
          if (sign(1.0_p_r8, p1) /= sign(1.0_p_r8, p2)) exit
        end do
        if (dGcol <= eps) exit
        gcol = gcol - dGcol
        dGcol = dGcol * 0.25_p_r8
      end do
      gLats(j) = 90.0_p_r8 - rad * gcol
      gLats(varCommon%yMax - j + 1) = -gLats(j)
    end do
  end subroutine gaussianLatitudes


  subroutine legendrePolynomial (n, colatitude, pln)
    !# Legendre Polynomial
    !# ---
    !# @info
    !# **Brief:** Legendre Polynomial. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo 
    implicit none
    integer, intent(in) :: n
    real (kind = p_r8), intent(in) :: colatitude
    real (kind = p_r8), intent(out) :: pln

    integer :: i
    real (kind = p_r8) :: x, y1, y2, y3, g

    x = cos(colatitude)
    y1 = 1.0_p_r8
    y2 = x
    do i = 2, n
      g = x * y2
      y3 = g - y1 + g - (g - y1) / real(i, p_r8)
      y1 = y2
      y2 = y3
    end do
    pln = y3
  end subroutine legendrePolynomial


  subroutine allocateData()
    !# Allocates Data
    !# ---
    !# @info
    !# **Brief:** Allocates Data. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo 
    implicit none
    allocate (landSeaMask(varCommon%xMax, varCommon%yMax))
    allocate (snowModel(varCommon%xMax, varCommon%yMax))
    allocate (albedoBefore(varCommon%xMax, varCommon%yMax), albedoAfter(varCommon%xMax, varCommon%yMax))
    allocate (snowOut(varCommon%xMax, varCommon%yMax), albedo(varCommon%xMax, varCommon%yMax))
    if (varCommon%grads) then
      allocate (gLats(varCommon%yMax))
    endif
  end subroutine allocateData


  subroutine deallocateData()
    !# Deallocates Data
    !# ---
    !# @info
    !# **Brief:** Deallocates Data. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo 
    implicit none
    deallocate (landSeaMask)
    deallocate (snowModel)
    deallocate (albedoBefore, albedoAfter)
    deallocate (snowOut, albedo)
    if (varCommon%grads) then
      deallocate (gLats)
    endif
  end subroutine deallocateData


  subroutine generateGrads()
    !# Generates Grads
    !# ---
    !# @info
    !# **Brief:** Generates Grads. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo 
    
    ! Climatological Snow at Model Grid and
    ! Time of Initial Condititon : GrADS
    ! Getting Gaussian Latitudes
    implicit none
    character (len = maxPathLength) :: ctlPathFileName
    character (len = maxPathLength) :: gradsBaseName
    character (len = maxPathLength) :: ctlFileName
    integer :: ctlFileUnit      ! Temp File Unit

    call gaussianLatitudes ()

    ! Write GrADS Control File
    gradsBaseName = trim(varName) // varCommon%date // exts // nLats
    ctlFileName = trim(gradsBaseName) // '.ctl'
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(ctlFileName)
    ctlFileUnit = openFile(ctlPathFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if(ctlFileUnit < 0) return
    write (unit = ctlFileUnit, fmt = '(A)') 'DSET ' // &
      trim(varCommon%dirModelIn) // varNameS // varCommon%date // exts // nLats
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'TITLE Climatological Snow on a Gaussian Grid'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
      'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write (unit = ctlFileUnit, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
    write (unit = ctlFileUnit, fmt = '(8F10.5)') gLats(varCommon%yMax:1:-1)
    write (unit = ctlFileUnit, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
    write (unit = ctlFileUnit, fmt = '(A)') 'TDEF 1 LINEAR ' // timeGrADS // ' 6HR'
    write (unit = ctlFileUnit, fmt = '(A)') 'VARS 1'
    write (unit = ctlFileUnit, fmt = '(A)') 'SNOW 0 99 Climatological Snow Depth [kg/m2]'
    write (unit = ctlFileUnit, fmt = '(A)') 'ENDVARS'
    close (unit = ctlFileUnit)
  end subroutine generateGrads


  subroutine printNameList
    !# Prints NameList
    !# ---
    !# @info
    !# **Brief:** Prints NameList. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: may/2019 </br>
    !# @endinfo 
    implicit none

    write (unit = p_nfprt, fmt = '(/,A)')  ' &SnowClimaNamelist'
    write (unit = p_nfprt, fmt = '(A,I6)') '     xMax = ', varCommon%xMax
    write (unit = p_nfprt, fmt = '(A,I6)') '     yMax = ', varCommon%yMax
    write (unit = p_nfprt, fmt = '(A)')    '     date = ' // varCommon%date
    write (unit = p_nfprt, fmt = '(A,L6)') '    grads = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)')    '  dirPreOut = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, fmt = '(A)')    '  dirModelIn = ' // trim(varCommon%dirModelIn)
    write (unit = p_nfprt, fmt = '(A)')    '  nameLSM = ' // trim(var%nameLSM)
  end subroutine printNameList

end module Mod_SnowClima

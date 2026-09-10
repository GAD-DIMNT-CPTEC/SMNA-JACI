!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_AlbedoClima </br></br>
!#
!# **Brief**: Module responsible for generating AlbedoClima files </br></br>
!#
!# To create the accumulated snow field you must have a climatological albedo
!# field. This albedo data set will be created by the AlbedoClima.f90 subroutine.
!# In the first step, the program reads the albedo.form file with the albedo
!# information and creates 2 files, one with the global map with the albedo
!# distribution, this data is placed in the AlbedoClima.dat file and a AlbedoClima.ctl
!# descriptor file for the GrADS. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/databcs/albedo.form </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/AlbedoClima.dat </br>
!# &bull; pre/dataout/AlbedoClima.ctl
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
!#  <li>01-08-2007 - Tomita - version: 1.1.1.1 </li>
!#  <li>10-12-2018 - Daniel Lamosa - version: 2.0.0 - Module creation </li>
!#  <li>25-03-2019 - Denis Eiras - version: 2.1.0 - Parallel Pre Program version </li>
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

module Mod_AlbedoClima

  use Mod_LinearInterpolation, only : gLatsL => latOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, initAreaInterpolation, doAreaInterpolation
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgOut, msgInLineFormatOut, msgNewLine

  implicit none

  public :: generateAlbedoClima
  public :: getNameAlbedoClima
  public :: initAlbedoClima
  public :: shouldRunAlbedoClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !parameters
  character (len = 11), parameter :: fileBCs = 'albedo.form'
  !# input file name

  !input variables
  type AlbedoClimaNameListData
  !# Type for namelist
    integer :: xDim
    !# Number of Longitudes
    integer :: yDim
    !# Number of Latitudes
    ! set values are default
    character(len = maxPathLength) :: varName = 'AlbedoClima'
    !# output prefix file name
  end type AlbedoClimaNameListData

  type(varCommonNameListData)   :: varCommon
  !# common namelist variable
  type(AlbedoClimaNameListData) :: var
  !# variable for use in namelist
  namelist /AlbedoClimaNameList/   var

  ! internal variables
  real (kind = p_r4), dimension (:, :), allocatable :: Albedo
  character(len=*), parameter :: header = 'Albedo Clima          | '


contains


  function getNameAlbedoClima() result(returnModuleName)
    !# Returns Albedo Clima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns Albedo Clima Module Name </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    !# name of the module

    returnModuleName = "AlbedoClima"
  end function getNameAlbedoClima


  subroutine initAlbedoClima(nameListFileUnit, varCommon_)
    !# Initializes Albedo Clima module
    !# ---
    !# @info
    !# **Brief:** Initializes AlbedoClima module, defined in PRE_run.nml </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    !# namelist file Unit comming from Pre_Program
    type(varCommonNameListData), intent(in) :: varCommon_
    !# common namelist variable
    read(unit = nameListFileUnit, nml = AlbedoClimaNameList)
    varCommon = varCommon_
  end subroutine initAlbedoClima


  function shouldRunAlbedoClima() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun
    !# true if module should run (there's output files)

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunAlbedoClima


  function getOutFileName() result(albedoOutFilename)
    !# Gets Albedo Clima Out file name
    !# ---
    !# @info
    !# **Brief:** Gets Albedo Clima Out file name </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: albedoOutFilename
    !# output file name

    albedoOutFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
  end function getOutFileName


  function generateAlbedoClima() result(isExecOk)
    !# Generates Albedo Clima
    !# ---
    !# @info
    !# **Brief:** Generates Albedo Clima output </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk
    !# true if execution was sucessfull
    character(len = maxPathLength) :: albedoInFilename, albedoOutFilename
    integer :: outFileUnit, inFileUnit, outRecSize, month

    isExecOk = .false.
    allocate (Albedo(var%xDim, var%yDim))

    albedoInFilename = trim(varCommon%dirBCs) // trim(fileBCs)
    inFileUnit = openFile(trim(albedoInFilename), 'formatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return

    inquire (iolength = outRecSize) Albedo
    albedoOutFilename = getOutFileName()
    outFileUnit = openFile(trim(albedoOutFilename), 'unformatted', 'direct', outRecSize, 'write', 'replace')
    if(outFileUnit < 0) return

    do month = 1, 12
      read (unit = inFileUnit, fmt = '(5E15.8)') Albedo
      call FlipMatrix (var%xDim, var%yDim, Albedo)
      write (unit = outFileUnit, rec = month) Albedo
    enddo
    close (unit = inFileUnit)
    close (unit = outFileUnit)

    call msgInLineFormatOut(' ' // header // 'Outputfile generated: ', '(A)')
    call msgInLineFormatOut(albedoOutFilename,'(A)')
    call msgNewLine();

    if(varCommon%grads) then
      if (.not. generateGrads()) return
    endif

    deallocate(Albedo)

    isExecOk = .true.

  end function generateAlbedoClima


  function generateGrads() result(isExecOk)
    !# Generates Grads files
    !# ---
    !# @info
    !# **Brief:** Generates .ctl and .dat files to check output Grads files </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    implicit none
    character (len = maxPathLength) :: ctlPathFileName
    character (len = maxPathLength) :: gradsBaseName
    character (len = maxPathLength) :: ctlFileName
    integer :: ctlFileUnit
    logical :: isExecOk
    !# Temp File Unit

    isExecOk = .false.
    gradsBaseName = trim(var%varName)
    ctlFileName = trim(gradsBaseName) // '.ctl'
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(ctlFileName)
    
    ctlFileUnit = openFile(ctlPathFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if(ctlFileUnit < 0) return
    
    write (unit = ctlFileUnit, fmt = '(A)') 'DSET ^' // trim(gradsBaseName) // '.dat'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,1PG12.5)') 'UNDEF', -999.0_p_r8
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'TITLE Climatological Albedo'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
    'XDEF ', var%xDim, ' LINEAR ', 0.0_p_r4, 360.0_p_r4 / real(var%xDim, p_r4)
    write (unit = ctlFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
    'YDEF ', var%yDim, ' LINEAR ', -90.0_p_r4, 180.0_p_r4 / real(var%yDim - 1, p_r4)
    write (unit = ctlFileUnit, fmt = '(A)') 'ZDEF  1 LEVELS 1000'
    write (unit = ctlFileUnit, fmt = '(A)') 'TDEF 12 LINEAR JAN2005 1MO'
    write (unit = ctlFileUnit, fmt = '(A)') '*'
    write (unit = ctlFileUnit, fmt = '(A)') 'VARS  1'
    write (unit = ctlFileUnit, fmt = '(A)') 'ALBE  0 99 Albedo [No Dim]'
    write (unit = ctlFileUnit, fmt = '(A)') 'ENDVARS'
    close(unit = ctlFileUnit)

    call msgInLineFormatOut(' ' // header // 'Grads Ctl File generated: ', '(A)')
    call msgInLineFormatOut(ctlPathFileName,'(A)')
    call msgNewLine()

    isExecOk = .true.

  end function generateGrads


  subroutine printNameList()
    !# Prints Namelist of AlbedoClima
    !# ---
    !# @info
    !# **Brief:** Prints Namelist of Albedo Clima module, defined in PRE_run.nml </br>
    !# **Author**: Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    write (unit = p_nfprt, fmt = '(/,A)')  ' &AlbedoClimaNameList'
    write (unit = p_nfprt, fmt = '(A,I6)') '      xDim = ', var%xDim
    write (unit = p_nfprt, fmt = '(A,I6)') '      yDim = ', var%yDim
    write (unit = p_nfprt, fmt = '(A,L6)') '     grads = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)') '    dirBCs = ' // trim(varCommon%dirBCs)
    write (unit = p_nfprt, fmt = '(A)') ' dirPreOut = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, fmt = '(A)') '   varName = ' // trim(var%varName)
  end subroutine printNameList


  subroutine FlipMatrix (xDim, yDim, h)
    !# Flips Matrix
    !# ---
    !# @info
    !# **Brief:** Flips a Matrix Over I.D.L. and Greenwitch and Southern and Northern Hemispheres </br>
    !# **Authors**:</br>
    !# &bull; Tomita S. </br>
    !# **Date**: aug/2007 </br>
    !# @endinfo
    integer, intent(in) :: xDim
    !# X Column Dimension of h(xDim,yDim)
    integer, intent(in) :: yDim
    !# Y Column Dimension of h(xDim,yDim)
    real (kind = p_r4), intent(inout) :: h(xDim, yDim)
    !# Matrix to be Flipped

    integer :: xDimd, xDimd1
    real (kind = p_r4) :: wk(xDim, yDim)

    xDimd = xDim / 2
    xDimd1 = xDimd + 1

    wk = h
    h(1:xDimd, :) = wk(xDimd1:xDim, :)
    h(xDimd1:xDim, :) = wk(1:xDimd, :)

    wk = h
    h(:, 1:yDim) = wk(:, yDim:1:-1)
  end subroutine FlipMatrix

end module Mod_AlbedoClima

!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SoilMoistureWeekly </br></br>
!#
!# **Brief**: Module responsible for using daily soil moisture generated in 
!# CPTEC as another option for model initialization.</br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/datain/GL_SM.GPNR.YYYYMMDDHH.vfm (Ex.: pre/datain/GL_SM.GPNR.2015043000.vfm)
!# </br></br>
!# 
!# **intermediate files:**
!#
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ  (Ex.: pre/dataout/ModelLandSeaMask.G00450) </br>
!# &bull; temporary/SoilMoistureWeekly.YYYYMMDD </br>
!# &bull; temporary/SoilMoistureWeekly.YYYYMMDD.ctl
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/SoilMoistureWeeklyYYYYMMDD.GZZZZZ (Ex.: model/datain/SoilMoistureWeekly20150430.G00450) </br>
!# &bull; pre/dataout/SoilMoistureWeeklyYYYYMMDD.GZZZZZ </br>
!# &bull; pre/dataout/SoilMoistureWeeklyYYYYMMDD.GZZZZZ.ctl
!# </br></br>
!# 
!# The SoilMoistureWeeklyYYYYMMDD.GZZZZZ file that is in the model/datain
!# directory is the file that the model will use in its execution, has the 
!# specific format (unformatted, direct access, 64-byte integer) and the 
!# SoilMoistureWeeklyYYYYMMDD.GZZZZZ file that is in the pre/dataout directory
!# is specific to viewing in GrADS software.</br></br>
!# 
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti         - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita                  - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel M. Lamosa/Barbara A. G. P. Yamada - version: 2.0.0 </li>
!#  <li>18-12-2018 - Denis Eiras             - version: 2.0.1 </li>
!#  <li>20-05-2019 - Eduardo Khamis          - version: 2.1.0 </li>
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

module Mod_SoilMoistureWeekly

  use Mod_LinearInterpolation, only : gLatsL => latOut, initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only : gLatsA => gLats, initAreaInterpolation, doAreaInterpolation
  use Mod_Vfm, only : vfirec
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  public :: generateSoilMoistureWeekly
  public :: getNameSoilMoistureWeekly
  public :: initSoilMoistureWeekly
  public :: shouldRunSoilMoistureWeekly

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !input variables

  type SoilMoistureWeeklyNameListData
    integer :: xDim  
    !# Number of Longitudes at the NCEP SoilMoisture Weekly Grid
    integer :: yDim  
    !# Number of Latitudes at the NCEP SoilMoisture Weekly Grid
    integer :: zDim  
    !# Number of Levels at the NCEP SoilMoisture Weekly Grid
    logical :: linear = .true. 
    !# Flag for Linear (.TRUE.) or Area Weighted (.FALSE.) Interpolation
    character(len = maxPathLength) :: varName = 'SoilMoistureWeekly' 
    !# output prefix file name
  end type SoilMoistureWeeklyNameListData

  type(varCommonNameListData)          :: varCommon
  type(SoilMoistureWeeklyNameListData) :: var
  namelist /SoilMoistureWeeklyNameList/   var

  !internal variables

  !real(kind = p_r4), dimension(:, :   ), allocatable :: soilMoistureWeeklyIn
  real(kind = p_r8), dimension(:, :, :), allocatable :: soilMoistureWeeklyInput
  real(kind = p_r4), dimension(:, :, :), allocatable :: soilMoistureWeeklyOut
  real(kind = p_r8), dimension(:, :, :), allocatable :: soilMoistureWeeklyOutput
  logical :: flagInput(5)  
  !# Input  grid flags
  logical :: flagOutput(5) 
  !# Output grid flags

  ! other variables
  character (len = 7) :: nLats = '.G     '


contains


  function getNameSoilMoistureWeekly() result(returnModuleName)
    !# Returns Soil Moisture Weekly Module Name
    !# ---
    !# @info
    !# **Brief:** Returns Soil Moisture Weekly Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo 
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "SoilMoistureWeekly"
  end function getNameSoilMoistureWeekly


  function shouldRunSoilMoistureWeekly() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo 
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunSoilMoistureWeekly


  function getOutFileName() result(soilMoistureWeeklyFileName)
    !# Gets SoilMoistureWeekly Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SoilMoistureWeekly Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo 
    implicit none
    character(len = maxPathLength) :: soilMoistureWeeklyFileName

    soilMoistureWeeklyFileName = trim(varCommon%dirModelIn) // trim(var%varName) // '.' // trim(varCommon%date(1:8)) // nLats
  end function getOutFileName


  subroutine initSoilMoistureWeekly(nameListFileUnit, varCommon_)
    !# Initialization of SoilMoistureWeekly module
    !# ---
    !# @info
    !# **Brief:** Initialization of SoilMoistureWeekly module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo 
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = SoilMoistureWeeklyNameList)
    varCommon = varCommon_

    write (nLats(3:7), '(I5.5)') varCommon%yMax
  end subroutine initSoilMoistureWeekly


  function generateSoilMoistureWeekly() result(isExecOk)
    !# Generates Soil Moisture Weekly
    !# ---
    !# @info
    !# **Brief:** Generates Soil Moisture Weekly output. This subroutine is the
    !# main method for use this module. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    
    implicit none
    ! Linear (.TRUE.) or Area Weighted (.FALSE.) Interpolation
    ! Interpolatp_Lat0e Regular To Gaussian
    ! Regular Input Data is Assumed to be Oriented with
    ! the North Pole and Greenwich as the First Point
    ! Gaussian Output Data is Interpolated to be Oriented with
    ! the North Pole and Greenwich as the First Point
    ! Input for the AGCM is Assumed IEEE 32 Bits Big Endian

    real(kind = p_r8), parameter :: p_Lat0 = 90.0_p_r8    
    !# Start at North Pole
    real(kind = p_r8), parameter :: p_Lon0 = 0.0_p_r8     
    !# Start at Prime Meridian
    logical :: isExecOk
    INTEGER :: k,i,j
    isExecOk = .false.

    if (var%linear) then
      call initLinearInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_Lat0, p_Lon0)
    else
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
      call initAreaInterpolation(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if
    call allocateData()
    call readSoilMoistureWeekly()

    DO k=1,var%zDim

       ! Interpolate Input Regular Rough Length To Gaussian Grid on Output
       if(var%linear) then
          call doLinearInterpolation(soilMoistureWeeklyInput( :, :, k), soilMoistureWeeklyOutput( :, :, k))
       else
          call doAreaInterpolation(soilMoistureWeeklyInput  ( :, :, k), soilMoistureWeeklyOutput( :, :, k))
       end if
       soilMoistureWeeklyOut( :, :, k) = real(soilMoistureWeeklyOutput( :, :, k), p_r4)
       DO j=1,varCommon%yMax
          DO i=1,varCommon%xMax
             IF (SoilMoistureWeeklyOut(i,j,k) <=0.0_p_r4)THEN
                SoilMoistureWeeklyOut(i,j,k) = 0.9_p_r4
             END IF
          END DO
       END DO
    END DO
    call writeSoilMoistureWeekly(getOutFileName())

    if(varCommon%grads) then
      call generateGrads()
    end if
    call deallocateData()

    isExecOk = .true.
  end function generateSoilMoistureWeekly


  subroutine generateGrads()
    !# Generates Grads files
    !# ---
    !# @info
    !# **Brief:** Generates .ctl and .dat files to check output. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    implicit none
    character (len = maxPathLength) :: gradsPathFileName
    character (len = maxPathLength) :: gradsBaseName
    integer :: fileUnit ! Temp File Unit

    character (len = 12) :: TimeGrADS = '  Z         '
    character (len = 3), dimension(12) :: MonthChar = &
      (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
        'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC' /)
    integer :: month

    gradsBaseName = trim(var%varName) // '.' // trim(varCommon%date(1:8)) // nLats
    gradsPathFileName = trim(varCommon%dirPreOut) // trim(gradsBaseName)

    ! Save the same file in temporary directory to view in grads
    call writeSoilMoistureWeekly(trim(gradsPathFileName))

    fileUnit = openFile(trim(gradsPathFileName) // '.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if(fileUnit < 0) return
    write(unit = fileUnit, fmt = '(a)') 'DSET ^' // trim(gradsBaseName)
    write(unit = fileUnit, fmt = '(a)') '*'
    write(unit = fileUnit, fmt = '(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit = fileUnit, fmt = '(a)') '*'
    write(unit = fileUnit, fmt = '(a)') 'UNDEF -999.0'
    write(unit = fileUnit, fmt = '(a)') '*'
    write(unit = fileUnit, fmt = '(a)') 'TITLE SoilMoistureWeekly on a Gaussian Grid'
    write(unit = fileUnit, fmt = '(a)') '*'
    write(unit = fileUnit, fmt = '(a,i5,a,f8.3,f15.10)') 'XDEF ', varCommon%xMax, ' LINEAR ', &
      0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write(unit = fileUnit, fmt = '(a,i5,a)') 'YDEF ', varCommon%yMax, ' LEVELS '
    if(var%linear) then
      write(unit = fileUnit, fmt = '(8f10.5)') gLatsL(varCommon%yMax:1:-1)
    else
      write(unit = fileUnit, fmt = '(8f10.5)') gLatsA(varCommon%yMax:1:-1)
    end if

    TimeGrADS(1:2) = varCommon%date(9:10)
    TimeGrADS(4:5) = varCommon%date(7:8)
    TimeGrADS(9:12) = varCommon%date(1:4)
    read (varCommon%date(5:6), FMT = '(I2)') Month
    TimeGrADS(6:8) = MonthChar(Month)

    write(unit = fileUnit, fmt = '(a)') 'ZDEF 8 LEVELS 1 2 3 4 5 6 7 8'
    write(unit = fileUnit, fmt = '(a)') 'TDEF 1 LINEAR ' // TimeGrADS // ' 1MO'
    write(unit = fileUnit, fmt = '(a)') 'VARS 1'
    write(unit = fileUnit, fmt = '(a)') 'SoilMoisture 8 99 SoilMoistureWeekly [kg/m3]'
    write(unit = fileUnit, fmt = '(a)') 'ENDVARS'
    close(unit = fileUnit)
  end subroutine generateGrads


  subroutine writeSoilMoistureWeekly(outputFile)
    !# Writes Soil Moisture Weekly
    !# ---
    !# @info
    !# **Brief:** Writes Soil Moisture Weekly output file. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: dec/2018 </br>
    !# @endinfo
    
    implicit none
    character(len = *), intent(in) :: outputFile 
    !# complete path to write file
    integer :: lRecOut
    integer :: fileUnit ! Temp File Unit
    integer :: k

    inquire(iolength = lRecOut) SoilMoistureWeeklyOut(:,:,1)
    fileUnit = openFile(outputFile, 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(fileUnit < 0) return
    DO k=1,var%zDim
       write(unit = fileUnit, rec = k) SoilMoistureWeeklyOut(:,:,k)
    END DO
    close(unit = fileUnit)
  end subroutine writeSoilMoistureWeekly


  subroutine readSoilMoistureWeekly()
    implicit none

    integer :: lRecIn, i, j, k
    character (len = maxPathLength) :: soilMoistureWeeklyCPTECFilenameIn
    character (len = 25) :: fileSoilMoisture = 'GL_SM.GPNR.          .vfm'
    real (kind = p_r8), dimension (:, :), allocatable :: SoilMoisture
    real (kind = p_r8), dimension (:, :, :), allocatable :: SoilMoisture_aux
    integer :: fileUnit

    fileSoilMoisture(12:21) = varCommon%date(1:10)! hour
    soilMoistureWeeklyCPTECFilenameIn = trim(varCommon%dirPreIn) // trim(fileSoilMoisture)
    fileUnit = openFile(soilMoistureWeeklyCPTECFilenameIn, 'formatted', 'sequential', -1, 'read', 'old')
    if(fileUnit < 0) return

    allocate (SoilMoisture(var%xDim    , var%yDim))
    allocate (SoilMoisture_aux(var%xDim, var%yDim, var%zDim))
    call vfirec(fileUnit, SoilMoisture_aux(1, 1, 1), var%xDim * var%yDim * var%zDim, 'LIN')
    !SoilMoisture = SoilMoisture_aux
    close (unit = fileUnit)

    !inquire (iolength = lRecIn) SoilMoisture (1:var%xDim, 1:var%yDim, 1)
    do k = 1, var%zDim
      do j = 1, var%yDim
        do i = 1, var%xDim
            SoilMoisture(i, j) = max(SoilMoisture_aux(i, var%yDim + 1 - j, k), 0.0_p_r8)
        end do
      end do
      soilMoistureWeeklyInput(:, :, k) = cshift (SoilMoisture, shift = size(SoilMoisture, 1) / 2, dim = 1)
    end do

  end subroutine readSoilMoistureWeekly


  subroutine printNameList()
    !# Prints namelist of SoilMoiustureWeekly
    !# ---
    !# @info
    !# **Brief:** Print namelist SoilMoiustureWeekly, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    write (unit = p_nfprt, FMT = '(/,A)')  ' &SoilMoistureWeeklyNameList'
    write (unit = p_nfprt, FMT = '(A,I6)') '      xDim = ', var%xDim
    write (unit = p_nfprt, FMT = '(A,I6)') '      yDim = ', var%yDim
    write (unit = p_nfprt, FMT = '(A,I6)') '      zDim = ', var%zDim
    write (unit = p_nfprt, FMT = '(A,I6)') '      xMax = ', varCommon%xMax
    write (unit = p_nfprt, FMT = '(A,I6)') '      yMax = ', varCommon%yMax
    write (unit = p_nfprt, FMT = '(A,L6)') '     grads = ', varCommon%grads
    write (unit = p_nfprt, FMT = '(A,L6)') '    linear = ', var%linear

    write (unit = p_nfprt, FMT = '(A)') ' dirPreIn   = ' // trim(varCommon%dirPreIn)
    write (unit = p_nfprt, FMT = '(A)') ' dirPreOut  = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, FMT = '(A)') ' dirModelIn = ' // trim(varCommon%dirModelIn)
    write (unit = p_nfprt, FMT = '(A)') '    varName = ' // trim(var%varName)
    write (unit = p_nfprt, FMT = '(A)') '       date = ' // trim(varCommon%date)
  end subroutine printNameList


  subroutine allocateData()
    implicit none
    !allocate(soilMoistureWeeklyIn(var%xDim, var%yDim))
    allocate(soilMoistureWeeklyInput(var%xDim, var%yDim, var%zDim))
    allocate(soilMoistureWeeklyOut(varCommon%xMax, varCommon%yMax, var%zDim))
    allocate(soilMoistureWeeklyOutput(varCommon%xMax, varCommon%yMax, var%zDim))
  end subroutine allocateData


  subroutine deallocateData()
    implicit none
    !deallocate(soilMoistureWeeklyIn)
    deallocate(soilMoistureWeeklyInput)
    deallocate(soilMoistureWeeklyOut)
    deallocate(soilMoistureWeeklyOutput)
  end subroutine deallocateData


end module Mod_SoilMoistureWeekly


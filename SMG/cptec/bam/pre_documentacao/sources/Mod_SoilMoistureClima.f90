!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SoilMoistureClima </br></br>
!#
!# **Brief**: Module responsible for reading the file soilms.form which contains
!# the climatological data of soil moisture and then generates the global map of
!# moisture distribution in the grid of 1 X 1 degree. </br></br>
!# 
!# First Point of Initial Data is at North Pole and I. D. Line. First
!# Point of Output  Data is at North Pole and Greenwhich. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/databcs/soilms.form
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/SoilMoistureClima.dat
!# &bull; pre/dataout/SoilMoistureClima.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <table>
!# <ul type="disc">
!# <tr><td><li>18-09-2007 - </li></td> <td> Jose P. Bonatti  </td> <td> - version: 1.2.0 <td></tr>
!# <tr><td><li>01-04-2018 - </li></td> <td> Daniel M. Lamosa </td> <td> - version: 2.0.0 <td></tr>
!# <tr><td><li>12-03-2020 - </li></td> <td> Eduardo Khamis   </td> <td> - version: 2.1.0 <td></tr>
!# </ul>
!# </table>
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

module Mod_SoilMoistureClima

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut

  implicit none
  
  public :: initSoilMoistureClima
  public :: generateSoilMoistureClima
  public :: getnameSoilMoistureClima
  public :: shouldRunSoilMoistureClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !input variables
  type SoilMoistureClimaNameListData
    integer :: xDim  
    integer :: yDim  
    character(len=maxPathLength) :: varName='SoilMoistureClima'    
    !# output prefix file name
    character(len=maxPathLength) :: fileBCs='soilms.form'          
  end type SoilMoistureClimaNameListData

  type(varCommonNameListData)         :: varCommon
  type(SoilMoistureClimaNameListData) :: var
  namelist /SoilMoistureClimaNameList/   var  

  !internal variables
  real(kind=p_r4), dimension(:,:), allocatable :: soilMoisture
  
  character(len=*), parameter :: header = 'Soil Moisture Clima          | '


  contains

  
  function getNameSoilMoistureClima() result(returnModuleName)
    !# Returns SoilMoistureClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SoilMoistureClima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo 
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "SoilMoistureClima"
  end function getNameSoilMoistureClima


  subroutine initSoilMoistureClima(nameListFileUnit, varCommon_)
    !# Initialization of SoilMoistureClima module
    !# ---
    !# @info
    !# **Brief:** Initialization of SoilMoistureClima module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = SoilMoistureClimaNameList)
    varCommon = varCommon_
  end subroutine initSoilMoistureClima


  function shouldRunSoilMoistureClima() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunSoilMoistureClima


  function getOutFileName() result(soilMoistureClimaOutFilename)
    !# Gets SoilMoistureClima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SoilMoistureClima Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: mar/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: soilMoistureClimaOutFilename

    soilMoistureClimaOutFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
  end function getOutFileName

  
  subroutine flipMatrix(xDim, yDim, matrix)
    !# Flips Matrix
    !# ---
    !# @info
    !# **Brief:** Flips a Matrix Over I.D.L. and Greenwitch. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    integer, intent(in)            :: xDim               
    !# Column Dimension
    integer, intent(in)            :: yDim               
    !# Row Dimension
    real(kind=p_r4), intent(inout) :: matrix(xDim, yDim) 
    !# Matrix to be Flipped
 
    real (kind=p_r4) :: temp(xDim, yDim) 
    !# temporary matrix

    integer :: xDimd  
    integer :: xDimd1 
    !# xdimd position plus 1
    
    xDimd  = xDim / 2
    xDimd1 = xDimd + 1

    temp = matrix
    matrix(1:xDimd,:)     = temp(xDimd1:xDim,:)
    matrix(xDimd1:xDim,:) = temp(1:xDimd,:)
  end subroutine flipMatrix

  
  function generateGrads() result(isGradsOk)
    !# Generates Grads files 
    !# ---
    !# @info
    !# **Brief:** Generates .ctl file to check output. This routine differ to
    !# another, because the generation output file differ too. Ctl file put in
    !# output directory. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isGradsOk
    integer :: gradsFileUnit
    
    isGradsOk = .false.

    ! Write GrADS Control File
    gradsFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varName)//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if (gradsFileUnit < 0) return

    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'//trim(var%varName)//'.dat'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'UNDEF -999.0'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE CLimatological Soil Moisture'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',var%xDim, &
               ' LINEAR ', 0.0_p_r4, 360.0_p_r4 / real(var%xDim, p_r4)
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'YDEF ',var%yDim, &
               ' LINEAR ', -90.0_p_r4, 180.0_p_r4 / real(var%yDim-1, p_r4)
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF  1 LEVELS 1000'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 12 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS  1'
    write(unit=gradsFileUnit, fmt='(a)') 'SOMO  0 99 SoilMoisture [cm]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)

    isGradsOk = .true.
  end function generateGrads
  
  
  function generateSoilMoistureClima() result(isExecOk)
    !# Generates Soil Moisture Clima 
    !# ---
    !# @info
    !# **Brief:** Generates Soil Moisture Clima output. This subroutine is the 
    !# main method for use this module. Only file name of namelist is needed to
    !# use it. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk
    integer :: inFileUnit, outFileUnit
    
    !internal variables
    integer :: lRec  
    integer :: month 
    !# counter of loop

    isExecOk = .false.
    
    allocate(soilMoisture(var%xDim, var%yDim))   
    
    ! open soilms file to read -----------------------------------------------------
    inFileUnit = openFile(trim(varCommon%dirBCs)//trim(var%fileBCs), 'formatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return
    ! ------------------------------------------------------------------------------
    
    ! Open File to write -----------------------------------------------------------
    inquire(iolength=lRec) soilMoisture
    outFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varName)//'.dat', 'unformatted', 'direct', lRec, 'write', 'replace')
    if(outFileUnit < 0) return
    ! -------------------------------------------------------------------------------
    
    do month=1, 12
      read(unit=inFileUnit, fmt='(5e15.8)') soilMoisture
      call flipMatrix(var%xDim, var%yDim, soilMoisture)
      write(unit=outFileUnit, rec=month) soilMoisture
    end do

    ! Close files
    close (unit=inFileUnit)
    close (unit=outFileUnit)
    
    if(varCommon%grads) then
      if (.not. generateGrads()) return
    end if
    
    deallocate(soilMoisture)
    isExecOk = .true.
  end function generateSoilMoistureClima

end module Mod_SoilMoistureClima

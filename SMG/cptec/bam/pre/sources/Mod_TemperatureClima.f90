!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_TemperatureClima </br></br>
!#
!# **Brief**: Module responsible for reading the temperature data from the
!# clmt.form file and creating the global temperature map in a pre-defined
!# resolution of 0.5 degrees </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/databcs/clmt.form
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; pre/dataout/TemperatureClima.dat </br>
!# &bull; pre/dataout/TemperatureClima.ctl
!# </br></br>
!# 
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti   - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita            - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel M. Lamosa  - version: 2.0.0 </li>
!#  <li>04-02-2020 - Denis Eiras       - version: 2.1.0 </li>
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

module Mod_TemperatureClima

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  public generateTemperatureClima
  public initTemperatureClima
  public getNameTemperatureClima
  public shouldRunTemperatureClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'

  !input variables
  type TemperatureClimaNameListData
    integer :: xDim  
    !# Number of Longitudes
    integer :: yDim  
    !# Number of Latitudes
    character(len=maxPathLength) :: varName='TemperatureClima'  
    !# output prefix file name
    character(len=maxPathLength) :: fileBCs='clmt.form'         
    !# climate temperature filename
  end type TemperatureClimaNameListData

  type(varCommonNameListData) :: varCommon
  type(TemperatureClimaNameListData)    :: var
  namelist /TemperatureClimaNameList/      var
  !------------------------------------------------------------------------------------------
  
  real(kind=p_r4), dimension (:,:,:), allocatable :: climateTemp
  

  contains
  

  function getNameTemperatureClima() result(returnModuleName)
    !# Returns TemperatureClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns TemperatureClima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName 
    ! variable for store module name

    returnModuleName = "TemperatureClima"
  end function getNameTemperatureClima


  subroutine initTemperatureClima(nameListFileUnit, varCommon_)
    !# Returns TemperatureClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns TemperatureClima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    
    implicit none
    integer, intent(in) :: nameListFileUnit 
    !# file unit of namelist PRE_run.nml
    type(varCommonNameListData), intent(in) :: varCommon_ 
    !# variable of type varCommonNameListData for managing common variables at PRE_run.nml

    read(unit = nameListFileUnit, nml = TemperatureClimaNameList)
    varCommon = varCommon_
  end subroutine initTemperatureClima


  function shouldRunTemperatureClima() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    
    implicit none
    logical :: shouldRun 
    ! result variable for store if method should run

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunTemperatureClima


  function getOutFileName() result(outFileName)
    !# Gets TemperatureClima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets TemperatureClima Unformatted Climatological PorcClay Mask file name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    
    implicit none
    character(len = maxPathLength) :: outFileName  
    ! Unformatted Climatological PorcClay Mask

    outFileName = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
  end function getOutFileName


  subroutine allocateData()
    !# Allocates data
    !# ---
    !# @info
    !# **Brief:** Allocates matrixes and variables based in namelist. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    allocate(climateTemp(var%xDim, var%yDim, 12))
    climateTemp = 0.0_p_r4
    !var%xDim = var%xDim!-1   ! Commented minus 1, I don't know
  end subroutine allocateData
  

  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Dellocates matrixes. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    deallocate(climateTemp)
  end subroutine deallocateData


  function readTemperatureClima() result(isReadOk)
    !# Reads Delta Temp Coldest C
    !# ---
    !# @info
    !# **Brief:** Reads Delta Temp Coldest C. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    logical :: isReadOk 
    ! return value (if function was executed sucessfully)
    integer :: inRecSize 
    !# size of file record
    character (len = maxPathLength) :: inFilename 
    !# input filename
    integer :: inFileUnit 
    !# To Read Formatted Climatological PorcClay Mask

    isReadOk = .false.
    inquire (iolength=inRecSize) climateTemp(:,:,:)

    inFilename = trim(varCommon%dirBCs) // trim(var%fileBCs)
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return

    read(unit=inFileUnit, rec=1) climateTemp
    close(unit=inFileUnit)
    isReadOk = .true.
  end function readTemperatureClima


  function writeTemperatureClima() result(isWriteOk)
    !# Writes Delta Temp Coldest C
    !# ---
    !# @info
    !# **Brief:** Writes Delta Temp Coldest C. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    logical :: isWriteOk  
    ! return value (if function was executed sucessfully)
    integer :: lRecOut  
    !# size of file record
    integer :: datFileUnit  
    !# To Write Unformatted Climatological PorcClay Mask

    isWriteOk = .false.
    inquire(iolength=lRecOut) climateTemp(:,:,:)
    datFileUnit = openFile(getOutFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(datFileUnit < 0) return

    write(unit=datFileUnit, rec=1) climateTemp(1:var%xDim,var%yDim:1:-1,1:12)
    close(unit=datFileUnit)
    isWriteOk = .true.

  end function writeTemperatureClima

 
  function generateGrads() result(isGradsOk)
    !# Generates Grads files
    !# ---
    !# @info
    !# **Brief:** Generates .ctl and .dat files to check output. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    logical :: isGradsOk 
    ! return value (if function was executed sucessfully)
    integer :: lRecGad 
    !# size of file record
    real(kind=p_r4), dimension (:,:,:), allocatable :: clayMaskGad   
    !# PorcClay Mask in real values
    character (len = maxPathLength) :: ctlPathFileName 
    !# ctl Grads File name with path
    character (len = maxPathLength) :: gradsBaseName 
    !# dat GrADS base name without extension
    integer :: gradsFileUnit ! File Unit for all grads files

    ! dat file is the same of getOutputFile; not generated here
    isGradsOk = .false.
    gradsBaseName = trim(var%varName)
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.ctl'
    
    ! Write .ctl file -------------------------------------------------------------
    gradsFileUnit = openFile(trim(ctlPathFileName), 'formatted', 'sequential', -1, 'write', 'replace')
    if(gradsFileUnit < 0) return

    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'// trim(gradsBaseName) // '.dat'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'UNDEF -999.0'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE CLimatological Temp'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',var%xDim,' LINEAR ', &
               0.0_p_r4, 360.0_p_r4 / real(var%xDim, p_r4)
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'YDEF ',var%yDim,' LINEAR ',  &
               -90.0_p_r4, 180.0_p_r4 / real(var%yDim, p_r4)
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF 1 LEVELS 1000'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 12 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS 1'
    write(unit=gradsFileUnit, fmt='(a)') 'clmt 0 99 TemperatureClima [%]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)
    isGradsOk = .true.
  end function generateGrads
  

  function generateTemperatureClima() result(isExecOk)
    !# Generates Delta Temp Coldest Clima
    !# ---
    !# @info
    !# **Brief:** Generates Delta Temp Coldest Clima output. This subroutine is
    !# the main method for use this module. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    logical :: isExecOk 
    ! return value (if function was executed sucessfully)
    integer :: i

    isExecOk = .false.
    call allocateData()
    if (.not. readTemperatureClima()) return
    if(.not. writeTemperatureClima()) return
    ! Generate grads output for debug or test
    if(varCommon%grads) then
      if (.not. generateGrads()) return
    end if
    call deallocateData()

    isExecOk = .true.
  end function generateTemperatureClima


end module Mod_TemperatureClima

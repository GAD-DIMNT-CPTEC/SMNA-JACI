!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_PorceClayMaskIBISClima </br></br>
!#
!# **Brief**: Module responsible for generating a global clay percentage field
!# in the soil </br></br>
!# 
!# This subroutine reads the file claymsk.form and generates the global 
!# distribution of the clay percentage in the soil and places the file 
!# PorceClayMaskIBISClima.dat in the pre/dataout directory with the information.
!# The same data is also generated in the PorceClayMaskIBISClimaG.dat and
!# PorceClayMaskIBISClimaG.ctl files for viewing in the GrADS software. </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/databcs/claymsk.form </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/PorceClayMaskIBISClima.dat </br>
!# &bull; pre/dataout/PorceClayMaskIBISClimaG.dat </br>
!# &bull; pre/dataout/PorceClayMaskIBISClimaG.ctl
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
!#  <li>20-01-2020 - Denis Eiras       - version: 2.1.0 </li>
!#  <li>15-10-2020 - Eduardo Khamis    - version: 2.1.1 </li>
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

module Mod_PorceClayMaskIBISClima

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  public initPorceClayMaskIBISClima
  public generatePorceClayMaskIBISClima
  public getNamePorceClayMaskIBISClima
  public shouldRunPorceClayMaskIBISClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'

  !input variables
  type PorceClayMaskIBISClimaNameListData
    integer :: xDim  
    !# Number of Longitudes
    integer :: yDim  
    !# Number of Latitudes
    integer :: layer
    !# set values are default
    character(len = maxPathLength) :: varName='PorceClayMaskIBISClima'     
    !# output prefix file name
    character(len = maxPathLength) :: varNameG='PorceClayMaskIBISClimaG'   
    !# output grads file name
    character(len = maxPathLength) :: fileBCs='claymsk.form' 
  end type PorceClayMaskIBISClimaNameListData

  type(varCommonNameListData) :: varCommon
  type(PorceClayMaskIBISClimaNameListData)    :: var
  namelist /PorceClayMaskIBISClimaNameList/      var

  !Internal variables
  integer,         dimension (:,:,:), allocatable :: porcClayMask  
  real(kind=p_r4), dimension (:,:,:), allocatable :: porcClayMask2 


  contains


  function getNamePorceClayMaskIBISClima() result(returnModuleName)
    !# Returns PorceClayMaskIBISClima Module Name
    !# ---
    !# @info
    !# **Brief:** Returns PorceClayMaskIBISClima Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName ! variable for store module name

    returnModuleName = "PorceClayMaskIBISClima"
  end function getNamePorceClayMaskIBISClima


  subroutine initPorceClayMaskIBISClima(nameListFileUnit, varCommon_)
    !# Initializes PorceClayMaskIBISClima module
    !# ---
    !# @info
    !# **Brief:** Initializes PorceClayMaskIBISClima module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit 
    !# file unit of namelist PRE_run.nml
    type(varCommonNameListData), intent(in) :: varCommon_ 
    !# variable of type varCommonNameListData for managing common variables at PRE_run.nml

    read(unit = nameListFileUnit, nml = PorceClayMaskIBISClimaNameList)
    varCommon = varCommon_
  end subroutine initPorceClayMaskIBISClima


  function shouldRunPorceClayMaskIBISClima() result(shouldRun)
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
    logical :: shouldRun ! result variable for store if method should run

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunPorceClayMaskIBISClima


  function getOutFileName() result(outFileName)
    !# Gets PorceClayMaskIBISClima Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets PorceClayMaskIBISClima Unformatted Climatological PorcClay
    !# Mask file name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: outFileName  ! Unformatted Climatological PorcClay Mask

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
    implicit none
    allocate(porcClayMask(var%xDim,var%yDim,var%layer))
    allocate(porcClayMask2(var%xDim,var%yDim,var%layer))
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
    implicit none
    deallocate(porcClayMask)
    deallocate(porcClayMask2)
  end subroutine deallocateData

  
  function readPorcClayMask() result(isReadOk)
    !# Reads PorcClay Mask
    !# ---
    !# @info
    !# **Brief:** Reads PorcClay Mask file and load in porcClayMask2. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReadOk ! return value (if function was executed sucessfully)
    integer :: inRecSize 
    !# size of file record
    character (len = maxPathLength) :: inFilename 
    !# input filename
    integer :: inFileUnit 
    !# To Read Formatted Climatological PorcClay Mask

    isReadOk = .false.
    inquire (iolength=inRecSize) porcClayMask2

    inFilename = trim(varCommon%dirBCs) // trim(var%fileBCs)
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return

    read(unit=inFileUnit, rec=1) porcClayMask2
    where(porcClayMask2 > 100)
      porcClayMask2=0.0
    end where
    porcClayMask = int(porcClayMask2)
    
    close(unit=inFileUnit)
    isReadOk = .true.
  end function readPorcClayMask
  

  function writePorcClayMask() result(isWriteOk)
    !# Writes PorcClay Mask
    !# ---
    !# @info
    !# **Brief:** Writes a file with PorcClay Mask. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isWriteOk  !return value (if function was executed sucessfully)
    integer :: lRecOut  
    !# size of file record
    integer :: datFileUnit  
    !# To Write Unformatted Climatological PorcClay Mask

    isWriteOk = .false.
    inquire(iolength=lRecOut) porcClayMask
    datFileUnit = openFile(getOutFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(datFileUnit < 0) return

    write(unit=datFileUnit, rec=1) porcClayMask
    close(unit=datFileUnit)
    isWriteOk = .true.

  end function writePorcClayMask
  

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
    logical :: isGradsOk !return value (if function was executed sucessfully)
    integer :: lRecGad 
    !# size of file record
    real(kind=p_r4), dimension (:,:,:), allocatable :: clayMaskGad   
    !# PorcClay Mask in real values
    character (len = maxPathLength) :: ctlPathFileName 
    !# ctl Grads File name with path
    character (len = maxPathLength) :: gradsBaseName 
    !# dat GrADS base name without extension
    integer :: gradsFileUnit 
    !# File Unit for all grads files

    isGradsOk = .false.
    gradsBaseName = trim(var%varNameG)
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.ctl'

    allocate(clayMaskGad(var%xDim,var%yDim,var%layer))

    ! write .dat file --------------------------------------------------------------
    inquire(iolength=lRecGad) clayMaskGad
    gradsFileUnit = openFile(trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.dat', 'unformatted', 'direct', lRecGad, 'write', 'replace')
    if(gradsFileUnit < 0) return

    clayMaskGad = real(porcClayMask, p_r4)
    write(unit=gradsFileUnit, rec=1) clayMaskGad
    close(unit=gradsFileUnit)
    ! -----------------------------------------------------------------------------

    ! Write .ctl file -------------------------------------------------------------
    gradsFileUnit = openFile(trim(ctlPathFileName), 'formatted', 'sequential', -1, 'write', 'replace')
    if(gradsFileUnit < 0) return

    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'// trim(gradsBaseName) // '.dat'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'UNDEF -999.0'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE IBIS PorcClay Mask'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',var%xDim,' LINEAR ', &
               0.0_p_r4, 360.0_p_r4 / real(var%xDim, p_r4)
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'YDEF ',var%yDim,' LINEAR ', &
               -89.5_p_r4, 179.0_p_r4 / real(var%yDim-1, p_r4)
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF 6 LEVELS 0 1 2 3 4 5'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 1 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS 1'
    write(unit=gradsFileUnit, fmt='(a)') 'VEGM 6 99 PorcClay Mask [No Dim]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)
    ! -----------------------------------------------------------------------------
    isGradsOk = .true.
  end function generateGrads
  

  subroutine flipMatrix(xDim, yDim, layer, h)
    !# Flips matrix
    !# ---
    !# @info
    !# **Brief:** Flips over the rows of a matrix, after flips over I.D.L. and
    !# Greenwitch. Matrix to flip is porcClayMask. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    ! Eduardo Khamis - bug fix - 15/10/2020
    integer, intent(in) :: xDim
    integer, intent(in) :: yDim
    integer, intent(in) :: layer
    integer, intent(inout) :: h(xDim, yDim, layer)

    integer :: xDimd  !start position of new dimension
    integer :: xDimd1 
    !# xdimd position plus 1
    integer :: temp(xDim, yDim, layer) 
    !# temporary matrix

    xDimd  = xDim / 2
    xDimd1 = xDimd + 1

    temp = h
    h(1:xDimd,:,:)     = temp(xDimd1:xDim,:,:)
    h(xDimd1:xDim,:,:) = temp(1:xDimd,:,:)

    temp = h
    ! *************************************************
    ! ******************** WARNNIG ********************
    ! *************************************************
    ! Comment was preserved below.
    ! *************************************************
    !h(:,1:yDim,:)=temp(:,yDim:1:-1,:)

  end subroutine flipMatrix
  

  function generatePorceClayMaskIBISClima() result(isExecOk)
    !# Generates PorceClay Mask IBIS Clima 
    !# ---
    !# @info
    !# **Brief:** Generates PorceClay Mask IBIS Clima output. This subroutine is
    !# the main method for use this module. Only file name of namelist is needed
    !# to use it.. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk ! return value (if function was executed sucessfully)

    isExecOk = .false.
    call allocateData()
    if (.not. readPorcClayMask()) return
    call flipMatrix(var%xDim, var%yDim, var%layer, porcClayMask)
    if (.not. writePorcClayMask()) return
    ! Generate grads output for debug or test
    if(varCommon%grads) then
      if (.not. generateGrads()) return
    end if
    call deallocateData()

    isExecOk = .true.
  end function generatePorceClayMaskIBISClima
  

end module Mod_PorceClayMaskIBISClima

!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_PorceSandMaskIBISClima </br></br>
!#
!# **Brief**: Module responsible for reading the file sandmsk.form that contains
!# the data of percentage of sand in the soil and generating the global map of
!# the distribution of this field </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/databcs/sandmsk.form </br></br>
!# 
!# **Files out:**
!#
!# &bull; pre/dataout/PorceSandMaskIBISClima.dat
!# &bull; pre/dataout/PorceSandMaskIBISClimaG.dat
!# &bull; pre/dataout/PorceSandMaskIBISClimaG.ctl
!# </br></br>
!# 
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>           - Jose P. Bonatti   - version: 1.0.0</li>
!#  <li>01-08-2007 - Tomita            - version: 1.1.1.1</li>
!#  <li>01-04-2018 - Daniel M. Lamosa  - version: 2.0.0</li>
!#  <li>28-01-2020 - Denis Eiras       - version: 2.1.0</li>
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
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2>
!# **&#9993;**<mailto:atende.cptec@inpe.br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---

module Mod_PorceSandMaskIBISClima

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  public initPorceSandMaskIBISClima
  public generatePorceSandMaskIBISClima
  public getNamePorceSandMaskIBISClima
  public shouldRunPorceSandMaskIBISClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'

  !input variables
  type PorceSandMaskIBISClimaNameListData
    integer :: xDim  !< Number of Longitudes
    integer :: yDim  !< Number of Latitudes
    integer :: layer  !
    ! set values are default
  !  character(len=128) :: dirPreIn='./'                        !< input data directory
  !  character(len=128) :: dirPreTemp='./'                      !< temporary data directory
  !  character(len=128) :: dirPreOut='./'                       !< output data directory
    character(len = maxPathLength) :: varName='PorceSandMaskIBISClima'     !< output prefix file name
    character(len = maxPathLength) :: varNameG='PorceSandMaskIBISClimaG'   !< output grads file name
    character(len = maxPathLength) :: fileBCs='sandmsk.form'               !< ???
  end type PorceSandMaskIBISClimaNameListData

  type(varCommonNameListData) :: varCommon
  type(PorceSandMaskIBISClimaNameListData)    :: var
  namelist /PorceSandMaskIBISClimaNameList/      var

  !Internal variables
  integer,         dimension (:,:,:), allocatable :: porcSandMask  !< ???
  real(kind=p_r4), dimension (:,:,:), allocatable :: porcSandMask2 !< ???


  contains

  ! ---------------------------------------------------------------------------
  !> @brief Returns PorceSandMaskIBISClima Module Name
  !!
  !! @details Returns PorceSandMaskIBISClima Module Name.\n
  !!
  !! @author Denis Eiras
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  function getNamePorceSandMaskIBISClima() result(returnModuleName)
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName !< variable for store module name

    returnModuleName = "PorceSandMaskIBISClima"
  end function getNamePorceSandMaskIBISClima

  ! ---------------------------------------------------------------------------
  !> @brief Initialization of PorceSandMaskIBISClima module
  !!
  !! @details Initialization of PorceSandMaskIBISClima module, defined in PRE_run.nml
  !!
  !! @author Denis Eiras
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  subroutine initPorceSandMaskIBISClima(nameListFileUnit, varCommon_)
    implicit none
    integer, intent(in) :: nameListFileUnit !< file unit of namelist PRE_run.nml
    type(varCommonNameListData), intent(in) :: varCommon_ !< variable of type varCommonNameListData for managing common variables at PRE_run.nml

    read(unit = nameListFileUnit, nml = PorceSandMaskIBISClimaNameList)
    varCommon = varCommon_
  end subroutine initPorceSandMaskIBISClima

  ! ---------------------------------------------------------------------------
  !> @brief Returns true if Module Should Run as a dependency
  !!
  !! @details Returns true if Module Should Run as a dependency, when it does not generated its out files and was not\n
  !! marked to run
  !!
  !! @author Denis Eiras
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  function shouldRunPorceSandMaskIBISClima() result(shouldRun)
    implicit none
    logical :: shouldRun !< result variable for store if method should run

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunPorceSandMaskIBISClima

  ! ---------------------------------------------------------------------------
  !> @brief Get PorceSandMaskIBISClima Out Filename
  !!
  !! @details Get PorceSandMaskIBISClima Unformatted Climatological PorcSand Mask file name
  !!
  !! @author Denis Eiras
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  function getOutFileName() result(outFileName)
    implicit none
    character(len = maxPathLength) :: outFileName  !< Unformatted Climatological PorcSand Mask

    outFileName = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
  end function getOutFileName

  ! ---------------------------------------------------------------------------
  !> @brief Allocate data
  !!
  !! @details Allocate matrixes and variables based in namelist
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine allocateData()
    implicit none
    allocate(porcSandMask(var%xDim,var%yDim,var%layer))
    allocate(porcSandMask2(var%xDim,var%yDim,var%layer))
  end subroutine allocateData                             
                          
  ! ---------------------------------------------------------------------------
  !> @brief Deallocate data
  !!
  !! @details dellocate matrixes
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine deallocateData()
    implicit none
    deallocate(porcSandMask)
    deallocate(porcSandMask2)
  end subroutine deallocateData

  
  ! ---------------------------------------------------------------------------
  !> @brief Read PorcSand Mask
  !!
  !! @details read PorcSand Mask file and load in porcSandMask2
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function readPorcSandMask() result(isReadOk)
    implicit none

    logical :: isReadOk !< return value (if function was executed sucessfully)
    integer :: inRecSize !< size of file record
    character (len = maxPathLength) :: inFilename !< input filename
    integer :: inFileUnit !< To Read Formatted Climatological PorcSand Mask

    isReadOk = .false.
    inquire (iolength=inRecSize) porcSandMask2

    inFilename = trim(varCommon%dirBCs) // trim(var%fileBCs)
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return

    read(unit=inFileUnit, rec=1) porcSandMask2
    where(porcSandMask2 > 100)
      porcSandMask2=0.0
    end where
    porcSandMask = int(porcSandMask2)
    
    close(unit=inFileUnit)
    isReadOk = .true.
  end function readPorcSandMask
  
  ! ---------------------------------------------------------------------------
  !> @brief Write porcSand Mask
  !!
  !! @details write a file with porcSand Mask
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function writePorcSandMask() result(isWriteOk)
    implicit none
    logical :: isWriteOk  !< return value (if function was executed sucessfully)
    integer :: lRecOut  !< size of file record
    integer :: datFileUnit  !< To Write Unformatted Climatological PorcSand Mask

    isWriteOk = .false.
    inquire(iolength=lRecOut) porcSandMask
    datFileUnit = openFile(getOutFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(datFileUnit < 0) return

    write(unit=datFileUnit, rec=1) porcSandMask
    close(unit=datFileUnit)
    isWriteOk = .true.

  end function writePorcSandMask
  
  ! ---------------------------------------------------------------------------
  !> @brief Generate Grads files 
  !!
  !! @details Generate .ctl and .dat files to check output
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function generateGrads() result(isGradsOk)
    implicit none
    logical :: isGradsOk !< return value (if function was executed sucessfully)
    integer :: lRecGad !< size of file record
    real(kind=p_r4), dimension (:,:,:), allocatable :: sandMaskGad   !< PorcSand Mask in real values
    character (len = maxPathLength) :: ctlPathFileName !< ctl Grads File name with path
    character (len = maxPathLength) :: gradsBaseName !< dat GrADS base name without extension
    integer :: gradsFileUnit ! File Unit for all grads files

    isGradsOk = .false.
    gradsBaseName = trim(var%varNameG)
    ctlPathFileName = trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.ctl'

    allocate(sandMaskGad(var%xDim,var%yDim,var%layer))

    ! write .dat file --------------------------------------------------------------
    inquire(iolength=lRecGad) sandMaskGad
    gradsFileUnit = openFile(trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.dat', 'unformatted', 'direct', lRecGad, 'write', 'replace')
    if(gradsFileUnit < 0) return

    sandMaskGad = real(porcSandMask, p_r4)
    write(unit=gradsFileUnit, rec=1) sandMaskGad
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
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE IBIS PorcSand Mask'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',var%xDim,' LINEAR ', &
               0.0_p_r4, 360.0_p_r4 / real(var%xDim, p_r4)
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'YDEF ',var%yDim,' LINEAR ', &
               -89.5_p_r4, 179.0_p_r4 / real(var%yDim-1, p_r4)
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF 6 LEVELS 0 1 2 3 4 5'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF 1 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS 1'
    write(unit=gradsFileUnit, fmt='(a)') 'VEGM 6 99 PorcSand Mask [No Dim]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)
    ! -----------------------------------------------------------------------------
    isGradsOk = .true.
  end function generateGrads
  
  ! ---------------------------------------------------------------------------
  !> @brief Flip matrix
  !!
  !! @details Flips over the rows of a matrix, after flips over I.D.L. and
  !! Greenwitch \n
  !! Matrix to flip is porcSandMask
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine flipMatrix()
    implicit none
    integer :: xDimd  !< start position of new dimension
    integer :: xDimd1 !< xdimd position plus 1
    integer :: temp(var%xDim,var%yDim,var%layer) !< temporary matrix

    xDimd  = var%xDim / 2
    xDimd1 = xDimd + 1

    temp = porcSandMask
    porcSandMask(1:xDimd,:,:)     = temp(xDimd1:var%xDim,:,:)
    porcSandMask(xDimd1:var%xDim,:,:) = temp(1:xDimd,:,:)

    temp = porcSandMask
    ! *************************************************
    ! ******************** WARNNIG ********************
    ! *************************************************
    ! Comment was preserved below.
    ! *************************************************
    !porcSandMask(:,1:yDim,:)=temp(:,yDim:1:-1,:)

  end subroutine flipMatrix
  
  ! ---------------------------------------------------------------------------
  !> @brief Generate PorceSand Mask IBIS Clima
  !!
  !! @details Generate PorceSand Mask IBIS Clima output. \n
  !! This subroutine is the main method for use this module. 
  !! Only file name of namelist is needed to use it.
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function generatePorceSandMaskIBISClima() result(isExecOk)
    implicit none
    logical :: isExecOk !< return value (if function was executed sucessfully)

    isExecOk = .false.
    call allocateData()
    if (.not. readPorcSandMask()) return
    call flipMatrix()
    if (.not. writePorcSandMask()) return
    ! Generate grads output for debug or test
    if(varCommon%grads) then
      if (.not. generateGrads()) return
    end if
    call deallocateData()

    isExecOk = .true.
  end function generatePorceSandMaskIBISClima
  

end module Mod_PorceSandMaskIBISClima

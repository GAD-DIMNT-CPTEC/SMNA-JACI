!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_VegetationMaskSSiB </br></br>
!#
!# **Brief**: Module responsible for generating the global mask of vegetation cover</br></br>
!# 
!# This subroutine generates the global mask of vegetation cover according with 
!# climatological data of the SsiB vegetation following Kuchler ranking. The 
!# input data for this subroutine is in the sibmsk.form file, which is in a 
!# special resolution of 1°.</br></br>
!#
!# **Files in:**
!#
!# &bull; pre/databcs/sibmsk.form</br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/VegetationMaskClima.dat
!# &bull; pre/dataout/VegetationMaskClimaG.dat
!# &bull; pre/dataout/VegetationMaskClimaG.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>??-??-???? - Jose P. Bonatti         - version: 1.0.0</li>
!#  <li>01-08-2007 - Tomita                  - version: 1.1.1.1</li>
!#  <li>01-04-2018 - Daniel M. Lamosa        - version: 2.0.0</li>
!#  <li>           - Barbara A. G. P. Yamada</li>
!#  <li>04-02-2020 - Eduardo Khamis          - version: 2.1.0</li>
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

module Mod_VegetationMaskSSiB

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut

  implicit none
  
  public :: generateVegetationMaskSSiB
  public :: getNameVegetationMaskSSiB
  public :: initVegetationMaskSSiB
  public :: shouldRunVegetationMaskSSiB

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'
  
  !input variables ---------------------------------------------------------------------------
  type VegetationMaskSSiBNameListData
    integer :: xDim  !< Number of Longitudes in SSiB Vegetation Mask Data
    integer :: yDim  !< Number of Latitudes in SSiB Vegetation Mask Data
    character(len=maxPathLength) :: varName='VegetationMaskClima'   !< VegetationMask SIB output file name prefix
    character(len=maxPathLength) :: varNameG='VegetationMaskClimaG' !< output grads file name
    character(len=maxPathLength) :: fileBCs='sibmsk.form'           !< SIB Vegetation mask input file name
  end type VegetationMaskSSiBNameListData

  type(varCommonNameListData)          :: varCommon
  type(VegetationMaskSSiBNameListData) :: var
  namelist /VegetationMaskSSiBNameList/   var

  !Internal variables
  integer,         dimension (:,:), allocatable :: vegetationMask !< ???
  real(kind=p_r4), dimension (:,:), allocatable :: vegMaskGad     !< ???

  character(len=*), parameter :: header = 'Vegetation Mask SSiB          | '


  contains  


  ! ---------------------------------------------------------------------------
  !> @brief Returns VegetationMaskSSiB Module Name
  !!
  !! @details Returns VegetationMaskSSiB Module Name.\n
  !!
  !! @author Eduardo Khamis
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  function getNameVegetationMaskSSiB() result(returnModuleName)
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "VegetationMaskSSiB"
  end function getNameVegetationMaskSSiB


  ! ---------------------------------------------------------------------------
  !> @brief Initialization of VegetationMaskSSiB module
  !!
  !! @details Initialization of VegetationMaskSSiB module, defined in PRE_run.nml
  !!
  !! @author Eduardo Khamis
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  subroutine initVegetationMaskSSiB(nameListFileUnit, varCommon_)
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = VegetationMaskSSiBNameList)
    varCommon = varCommon_
  end subroutine initVegetationMaskSSiB


  ! ---------------------------------------------------------------------------
  !> @brief Returns true if Module Should Run as a dependency
  !!
  !! @details Returns true if Module Should Run as a dependency, when it does not generated its out files and was not\n
  !! marked to run
  !!
  !! @author Eduardo Khamis
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  function shouldRunVegetationMaskSSiB() result(shouldRun)
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunVegetationMaskSSiB


  ! ---------------------------------------------------------------------------
  !> @brief Get VegetationMaskSSiB Out Filename
  !!
  !! @details Get VegetationMaskSSiB Out Filename
  !!
  !! @author Eduardo Khamis
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  function getOutFileName() result(vegetationMaskSSiBOutFilename)
    implicit none
    character(len = maxPathLength) :: vegetationMaskSSiBOutFilename

    vegetationMaskSSiBOutFilename = trim(varCommon%dirPreOut) // trim(var%varName) // '.dat'
  end function getOutFileName


  ! ---------------------------------------------------------------------------
  !> @brief Flip matrix
  !!
  !! @details Flips over the rows of a matrix, after flips over I.D.L. and
  !! Greenwitch \n
  !! Matrix to flip is VegetationMask
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine flipMatrix()
    implicit none
    integer :: xDimd  !< start position of new dimension
    integer :: xDimd1 !< xdimd position plus 1
    integer :: temp(var%xDim,var%yDim) !< temporary matrix

   xDimd  = var%xDim / 2
   xDimd1 = xDimd + 1

   temp = vegetationMask
   vegetationMask(1:xDimd,:) = temp(xDimd1:var%xDim,:)
   vegetationMask(xDimd1:var%xDim,:) = temp(1:xDimd,:)

   temp = vegetationMask
   vegetationMask(:,1:var%yDim) = temp(:,var%yDim:1:-1)

  end subroutine flipMatrix

  
  ! ---------------------------------------------------------------------------
  !> @brief Read sib mask file
  !!
  !! @details read sib mask file and load in vegetationMask
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function readSibMask() result(isReadOk)
    implicit none
    logical :: isReadOk
    integer :: lRecOut 
    integer :: inFileUnit 
    character(len=maxPathLength) :: sibMaskFile !< complete path do sib mask file
    
    isReadOk = .false.
    sibMaskFile=trim(varCommon%dirBCs)//trim(var%fileBCs)
    inFileUnit = openFile(trim(sibMaskFile), 'formatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return
    read(unit=inFileUnit, fmt='(40i2)') vegetationMask
    close(unit=inFileUnit)
    isReadOk = .true.
 
  end function readSibMask

  
  ! ---------------------------------------------------------------------------
  !> @brief Allocate matrixes 
  !!
  !! @details Allocate matrixes vetegetationMask and vegMaskGad
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine allocateMatrixes()
    implicit none
    allocate(vegetationMask(var%xDim, var%yDim))
    allocate(vegMaskGad(var%xDim, var%yDim))
  end subroutine allocateMatrixes
  

  ! ---------------------------------------------------------------------------
  !> @brief Deallocate matrixes 
  !!
  !! @details Deallocate matrixes vetegetationMask and vegMaskGad
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine deallocateMatrixes()
    implicit none
    deallocate(vegetationMask)
    deallocate(vegMaskGad)
  end subroutine deallocateMatrixes
  
  
  ! ---------------------------------------------------------------------------
  !> @brief Write Vegetation Mask 
  !!
  !! @details write a file with vetegetation Mask
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function writeVegMask() result(isWriteOk)
    implicit none
    logical :: isWriteOk
    integer :: lRecOut
    integer :: outFileUnit
    character(len=maxPathLength) :: datFileName !< output filename

    isWriteOk = .false.
    datFileName=trim(getOutFileName())
    inquire(iolength=lRecOut) vegetationMask
    outFileUnit = openFile(trim(datFilename), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(outFileUnit < 0) return
    write(unit=outFileUnit, rec=1) vegetationMask
    close(unit=outFileUnit)
    isWriteOk = .true.

  end function writeVegMask

  
  ! ---------------------------------------------------------------------------
  !> @brief Generate Grads files 
  !!
  !! @details Generate .ctl and .dat files to check output
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function generateGrads() result(isGradsOk)
    implicit none
    logical :: isGradsOk
    integer :: lRecGad
    integer :: gradsFileUnit
    integer :: ctlFileUnit
    character(len=maxPathLength) :: datFileName !< output filename
    character(len=maxPathLength) :: graFileName !< gra filename

    isGradsOk = .false.
 
    datFileName = trim(varCommon%dirPreOut)//trim(var%varNameG)//'.dat'
    graFileName = trim(varCommon%dirPreOut)//trim(var%varNameG)//'.ctl'
    ! create .dat file --------------------------------------------------------
    inquire(iolength=lRecGad) vegMaskGad
    gradsFileUnit = openFile(datFileName, 'unformatted', 'direct', lRecGad, 'write', 'replace')
    if (gradsFileUnit < 0) return
    vegMaskGad = real(vegetationMask,p_r4)
    write(unit=gradsFileUnit, rec=1) vegMaskGad
    close(unit=gradsFileUnit)
    ! -------------------------------------------------------------------------
    
    ! create .ctl file --------------------------------------------------------
    ctlFileUnit = openFile(graFileName, 'formatted', 'sequential', -1, 'write', 'replace')
    if (ctlFileUnit < 0) return
    
    write(unit=ctlFileUnit, fmt='(a)') 'DSET ^'//trim(var%varNameG)//'.dat'
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a)') 'UNDEF -999.0'
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a)') 'TITLE SSiB Vegetation Mask'
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a,i5,a,f8.3,f15.10)') &
                                 'XDEF ',var%xDim,' LINEAR ', 0.0_p_r4,360.0_p_r4/real(var%xDim,p_r4) 
    write(unit=ctlFileUnit, fmt='(a,i5,a,f8.3,f15.10)') &
                                 'YDEF ',var%yDim,' LINEAR ', -89.5_p_r4,179.0_p_r4/real(var%yDim-1,p_r4)
    write(unit=ctlFileUnit, fmt='(a)') 'ZDEF 1 LEVELS 1000'
    write(unit=ctlFileUnit, fmt='(a)') 'TDEF 1 LINEAR JAN2005 1MO'
    write(unit=ctlFileUnit, fmt='(a)') '*'
    write(unit=ctlFileUnit, fmt='(a)') 'VARS 1'
    write(unit=ctlFileUnit, fmt='(a)') 'VEGM 0 99 Vegetation Mask [No Dim]'
    write(unit=ctlFileUnit, fmt='(a)') 'ENDVARS'

    close(unit=ctlFileUnit) 
    ! -------------------------------------------------------------------------

    isGradsOk = .true.
  end function generateGrads

  
  ! ---------------------------------------------------------------------------
  !> @brief Generate Vegetation Mask SSib 
  !!
  !! @details Generate Vegetation Mask output. \n
  !! This subroutine is the main method for use this module. 
  !! Only file name of namelist is needed to use it.
  !!
  !! @authors Daniel M. Lamosa \n
  !> Barbara A. G. P. Yamada
  !!
  !! @date mar/2018
  !!       jan/2020 - Eduardo Khamis   - changing subroutine to function
  ! ---------------------------------------------------------------------------
  function generateVegetationMaskSSib() result(isExecOk)
    implicit none
    logical :: isExecOk

    isExecOk = .false.    

    call allocateMatrixes()

    if (.not. readSibMask()) then
      call msgWarningOut(header, "Error reading sibmsk.form file")
      return
    end if
    call flipMatrix()
    if (.not. writeVegMask()) then
      call msgWarningOut(header, "Error writing VegetationMaskClima file")
      return
    end if
    
    ! Generate grads output for debug or test
    if(varCommon%grads) then
      if (.not. generateGrads()) then
        call msgWarningOut(header, "Error while generating grads files")
        return
      end if
    end if

    call deallocateMatrixes()

    isExecOk = .true.
  end function generateVegetationMaskSSib


end module Mod_VegetationMaskSSiB

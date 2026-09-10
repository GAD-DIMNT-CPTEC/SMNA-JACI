!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_VegetationMask </br></br>
!#
!# **Brief**: Module responsible for interpolating the vegetation mask in the 
!# model resolution and correcting the land/sea mask of the model </br></br>
!# 
!# The VegetationMask.GZZZZ file that is in the model/datain directory is the 
!# file that the model will use in its execution, has the specific format 
!# (unformatted, direct access, 64-byte integer) and the VegetationMask.GZZZZ 
!# file that is in the pre/dataout directory is specific to viewing in GrADS 
!# software. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/VegetationMaskClima.dat </br>
!# &bull; pre/dataout/LandSeaMask.GZZZZZ (Ex.: pre/dataout/LandSeaMask.G00450)
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ </br>
!# &bull; model/datain/VegetationMask.GZZZZ </br>
!# &bull; pre/dataout/VegetationMask.GZZZZZ </br>
!# &bull; pre/dataout/VegetationMask.GZZZZZ.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti  - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita           - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel M. Lamosa - version: 2.0.0 </li>
!#  <li>04-02-2020 - Eduardo Khamis   - version: 2.1.0 </li>
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

module Mod_VegetationMask

  use Mod_AreaIntegerInterp, only: gLats, initAreaIntegerInterp, &
                                   doAreaIntegerInterp
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgWarningOut

  implicit none
  
  public :: generateVegetationMask
  public :: getNameVegetationMask
  public :: initVegetationMask
  public :: shouldRunVegetationMask

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  integer, parameter :: p_numVegClasses = 13
  integer, dimension (p_numVegClasses) :: vegClass = &
            (/ 1, 1, 1, 1, 1, 1, 2, 2, 3, 2, 3, 4, 5 /)
  integer, parameter :: p_undefV          = 0 
  !# Undefined value which if found in input array causes that location to be ignored in interpolation.  Used as the output value for output points with no defined and/or unmasked data

  !input variables
  type VegetationMaskNameListData
    integer :: xDim          
    !# Number of Longitudes in SSiB Vegetation Mask Data
    integer :: yDim          
    !# Number of Latitudes in SSiB Vegetation Mask Data
    character(len=maxPathLength) :: nameLSM='LandSeaMask'           
    !# Land Sea Mask input file name prefix
    character(len=maxPathLength) :: varName='VegetationMaskClima'   
    !# Vegetation mask clima input file name prefix
    character(len=maxPathLength) :: nameLSMSSiB='ModelLandSeaMask'  
    !# Land Sea Mask Sib output file name prefix
    character(len=maxPathLength) :: varNameVeg='VegetationMask'     
    !# Interpolated Vegetation Mask output file name prefix
  end type VegetationMaskNameListData

  type(varCommonNameListData)      :: varCommon
  type(VegetationMaskNameListData) :: var
  namelist /VegetationMaskNameList/   var
  
  !internal variables
  character(len=10)                    :: mskfmt = '(      i1)'
  !# format of landSeaMask file
  character(len=7)   :: nLats='.G     '                 
  !# posfix land sea mask file
  integer, dimension(:,:), allocatable :: vegetationMaskIn    
  integer, dimension(:,:), allocatable :: vegetationMaskOut   
  integer, dimension(:,:), allocatable :: vegMaskSave         
  integer, dimension(:,:), allocatable :: landSeaMask         
  integer, dimension(:,:), allocatable :: lSMaskSave           
  integer, dimension(:,:), allocatable :: vegetationMaskInput  
  integer, dimension(:,:), allocatable :: vegetationMaskOutput 
  logical :: flagInput(5), flagOutput(5)
  
  character(len=*), parameter :: header = 'Vegetation Mask          | '


  contains


  function getNameVegetationMask() result(returnModuleName)
    !# Returns VegetationMask Module Name
    !# ---
    !# @info
    !# **Brief:** Returns VegetationMask Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "VegetationMask"
  end function getNameVegetationMask


  subroutine initVegetationMask(nameListFileUnit, varCommon_)
    !# Initialization of VegetationMask module
    !# ---
    !# @info
    !# **Brief:** Initialization of VegetationMask module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = VegetationMaskNameList)
    varCommon = varCommon_

    write (nLats(3:7), '(I5.5)') varCommon%yMax
    write (mskfmt(2:7),  '(I6)') varCommon%xMax

    flagInput(1)=.true.   ! Start at North Pole
    flagInput(2)=.true.   ! Start at Prime Meridian
    flagInput(3)=.false.  ! Latitudes Are at North Edge
    flagInput(4)=.false.  ! Longitudes Are at Western Edge
    flagInput(5)=.false.  ! Regular Grid
    flagOutput(1)=.true.  ! Start at North Pole
    flagOutput(2)=.true.  ! Start at Prime Meridian
    flagOutput(3)=.false. ! Latitudes Are at North Edge of Box
    flagOutput(4)=.true.  ! Longitudes Are at Center of Box
    flagOutput(5)=.true.  ! Gaussian Grid

  end subroutine initVegetationMask


  function shouldRunVegetationMask() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunVegetationMask


  function getOutFileName() result(vegetationMaskOutFilename)
    !# Gets VegetationMask Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets VegetationMask Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: vegetationMaskOutFilename

    vegetationMaskOutFilename = trim(varCommon%dirModelIn) // trim(var%varNameVeg) // nLats
  end function getOutFileName

  
  function readlandSeaMask() result(isReadOk)
    !# Reads Land Sea Mask
    !# ---
    !# @info
    !# **Brief:** Reads Land Sea Mask, input file. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReadOk 
    character(len = maxPathLength) :: inFilename 
    integer :: lRecIn 
    integer :: inFileUnit 

    isReadOk = .false.

    inFileName = trim(varCommon%dirPreOut)//trim(var%nameLSM)//nLats
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'sequential', -1, 'read', 'old')
    if (inFileUnit < 0) return
    read(unit=inFileUnit) landSeaMask
    lSMaskSave = landSeaMask
    close(unit=inFileUnit)

    isReadOk = .true.
  end function readlandSeaMask
 
 
  function readVegetationMask() result(isReadOk)
    !# Reads Vegetation Mask
    !# ---
    !# @info
    !# **Brief:** Reads Vegetation Mask, input file. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReadOk 
    character(len = maxPathLength) :: inFilename 
    integer :: lRecIn 
    integer :: inFileUnit 

    isReadOk = .false.
    inquire(iolength=lRecIn) vegetationMaskIn
    inFileName = trim(varCommon%dirPreOut)//trim(var%varName)//'.dat'
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', lRecIn, 'read', 'old')
    if (inFileUnit < 0) return
    read(unit=inFileUnit, rec=1) vegetationMaskIn
    close(unit=inFileUnit)

    isReadOk = .true.
  end function readVegetationMask
  

  function writeInterpolatedVegetationMask() result(isWriteOk)
    !# Writes Interpolated Vegetation Mask
    !# ---
    !# @info
    !# **Brief:** Writes Out Adjusted Interpolated VegetationMask. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    logical :: isWriteOk 
    character(len=maxPathLength) :: outputFile 
    !# complete path to write file
    integer :: lRecOut
    integer :: outFileUnit
 
    isWriteOk = .false.   
    
    inquire(iolength=lRecOut) vegetationMaskOut
    outputFile = trim(varCommon%dirModelIn)//trim(var%varNameVeg)//nLats
    outFileUnit = openFile(outputFile, 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if (outFileUnit < 0) return
    write(unit=outFileUnit, rec=1) vegetationMaskOut
    close(unit=outFileUnit)

    isWriteOk = .true.
  end function writeInterpolatedVegetationMask
  

  function writeLandSeaMask() result(isWriteOk)
    !# Writes Land Sea Mask
    !# ---
    !# @info
    !# **Brief:** Writes Land Sea Mask. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    logical :: isWriteOk 
    character(len=maxPathLength) :: outputFile 
    !# complete path to write file
    integer :: outFileUnit
 
    isWriteOk = .false.   
    
    outputFile = trim(varCommon%dirPreOut)//trim(var%nameLSMSSiB)//nLats
    outFileUnit = openFile(outputFile, 'unformatted', 'sequential', -1, 'write', 'replace')
    if (outFileUnit < 0) return
    write(unit=outFileUnit) landSeaMask
    close(unit=outFileUnit)

    isWriteOk = .true.
  end function writeLandSeaMask
  

  subroutine allocateData()
    !# Allocates data
    !# ---
    !# @info
    !# **Brief:** Allocates matrixes based in namelist. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    allocate(vegetationMaskInput(var%xDim,var%yDim))
    allocate(vegetationMaskOutput(varCommon%xMax,varCommon%yMax))
    allocate(vegetationMaskIn(var%xDim,var%yDim))
    allocate(vegetationMaskOut(varCommon%xMax,varCommon%yMax))
    allocate(vegMaskSave(varCommon%xMax,varCommon%yMax))
    allocate(landSeaMask(varCommon%xMax,varCommon%yMax))
    allocate(lSMaskSave(varCommon%xMax,varCommon%yMax))

  end subroutine allocateData
  

  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Deallocates matrixes. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    deallocate(vegetationMaskInput)
    deallocate(vegetationMaskOutput)
    deallocate(vegetationMaskIn)
    deallocate(vegetationMaskOut)
    deallocate(vegMaskSave)
    deallocate(landSeaMask)
    deallocate(lSMaskSave)

  end subroutine deallocateData
   

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
    integer :: lRecGad
    integer :: gradsFileUnit
    real(kind=p_r4), dimension(:,:), allocatable :: gad
    real(kind=p_r4), parameter :: p_undefG = -99.0     
    
    isGradsOk = .false.

    allocate(gad(varCommon%xMax,varCommon%yMax))
    
    ! Write .dat file --------------------------------------------------------
    inquire(iolength=lRecGad) gad
    gradsFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varNameVeg)//nLats, 'unformatted', 'direct', lRecGad, 'write', 'replace')
    if (gradsFileUnit < 0) return

    gad=real(lSMaskSave, p_r4)
    write(unit=gradsFileUnit, rec=1) gad
    
    gad=real(landSeaMask, p_r4)
    write(unit=gradsFileUnit, rec=2) gad
    
    gad=real(vegMaskSave, p_r4)
    write(unit=gradsFileUnit, rec=3) gad
    
    gad=real(vegetationMaskOut, p_r4)
    write(unit=gradsFileUnit, rec=4) gad
    close(unit=gradsFileUnit)
    ! -------------------------------------------------------------------------
    
    ! Write .ctl file ---------------------------------------------------------
    gradsFileUnit = openFile(trim(varCommon%dirPreOut)//trim(var%varNameVeg)//nLats//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if (gradsFileUnit < 0) return 

    write(unit=gradsFileUnit, fmt='(a)') 'DSET ^'//trim(var%varNameVeg)//nLats//'.dat'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,1pg12.5)') 'UNDEF ', p_undefG
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a)') 'TITLE Vegetation Mask on a Gaussian Grid'
    write(unit=gradsFileUnit, fmt='(a)') '*'
    write(unit=gradsFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',varCommon%xMax,' LINEAR ', &
          0.0_p_r8,360.0_p_r8/real(varCommon%xMax,p_r8)
    write(unit=gradsFileUnit, fmt='(a,i5,a)') 'YDEF ',varCommon%yMax,' LEVELS '
    write(unit=gradsFileUnit, fmt='(8f10.5)') gLats(varCommon%yMax:1:-1)
    write(unit=gradsFileUnit, fmt='(a)') 'ZDEF  1 LEVELS 1000'
    write(unit=gradsFileUnit, fmt='(a)') 'TDEF  1 LINEAR JAN2005 1MO'
    write(unit=gradsFileUnit, fmt='(a)') 'VARS  4'
    write(unit=gradsFileUnit, fmt='(a)') 'LSMO  0 99 Land Sea Mask Before Fix  [0-Sea 1-Land]'
    write(unit=gradsFileUnit, fmt='(a)') 'LSMK  0 99 Land Sea Mask For Model   [0-Sea 1-Land]'
    write(unit=gradsFileUnit, fmt='(a)') 'VGMO  0 99 Vegetation Mask Before Fix [0 to 13]'
    write(unit=gradsFileUnit, fmt='(a)') 'VEGM  0 99 Vegetation Mask For Model  [0 to 13]'
    write(unit=gradsFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsFileUnit)
      
    deallocate(gad)

    isGradsOk = .true.

  end function generateGrads
  

  function generateVegetationMask() result(isExecOk)
    !# Generates Vegetation Mask
    !# ---
    !# @info
    !# **Brief:** Generates Vegetation Mask output. This subroutine is the main
    !# method for use this module. Only file name of namelist is needed to use it. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Eduardo Khamis - changing subroutine to function </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk
    
    logical, parameter :: dumpLocal=.false.

    isExecOk = .false.
    call allocateData()
    if (.not. readlandSeaMask()) then
      call msgWarningOut(header, "Error reading LandSeaMask file")
      return
    end if
    if (.not. readVegetationMask()) then
      call msgWarningOut(header, "Error reading VegetationMask file")
      return
    end if
    call initAreaIntegerInterp(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_numVegClasses, vegClass, flagInput, flagOutput)
    
    ! Interpolate Input Regular VegetationMask To Gaussian Grid Output
    vegetationMaskInput=int(vegetationMaskIn)
    if(.not. doAreaIntegerInterp(vegetationMaskInput, vegetationMaskOutput)) return

    vegMaskSave=int(vegetationMaskOutput)
    
    ! Fix Problems on Vegetation Mask and Land Sea Mask
    call vegetationMaskCheck(dumpLocal)
    vegetationMaskOut = int(vegetationMaskOutput)
      
    if (.not. writeInterpolatedVegetationMask()) then
      call msgWarningOut(header, "Error writing VegetationMask file")
      return
    end if
    if (.not. writeLandSeaMask()) then
      call msgWarningOut(header, "Error writing ModelLandSeaMask file")
      return
    end if
    
    ! Generate grads output for debug or test
    if(varCommon%grads) then
      if (.not. generateGrads()) then
        call msgWarningOut(header, "Error while generating grads files")
        return
      end if
    end if
    
    call deallocateData()

    isExecOk = .true.
  end function generateVegetationMask
  
  
  subroutine vegetationMaskCheck(dumpLocal)
    !# Vegetation Mask Check
    !# ---
    !# @info
    !# **Brief:** Vegetation Mask Check. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none

    logical, intent(in) :: dumpLocal
    
    logical :: flag=.false. 
    !# exit loop control
    integer :: i
    integer :: im
    integer :: ip
    integer :: j
    integer :: jm
    integer :: jp
    integer :: k1
    integer :: k2
    integer :: k3
    integer :: k4
    integer :: k5
    integer :: k6
    integer :: k7
    integer :: k8
    integer :: kk
    
    
    ! infinite loop
    do 
      mainLoop: do j=1, varCommon%yMax
        do i=1, varCommon%xMax
          if(landSeaMask(i, j) /= 0 .and. vegetationMaskOutput(i, j) == p_undefV) then

            if(im == 0) then
              im = varCommon%xMax
            else
              im = i - 1
            end if
            if(ip > varCommon%xMax) then
              ip = 1
            else
              ip = i + 1
            end if
            if(jm == 0) then
              jm = 1
            else
              jm = j - 1
            end if
            if(jp > varCommon%yMax) then
              jp = varCommon%yMax
            else
              jp = j + 1
            end if
            
            if(landSeaMask(im, jm) /= 0) then
              k1 = vegetationMaskOutput(im, jm)
            else
              k1 = 0
            end if
            
            if(landSeaMask(i, jm) /= 0) then
              k2 = vegetationMaskOutput(i, jm)
            else
              k2 = 0
            end if
            
            if(landSeaMask(ip, jm) /= 0) then
              k3 = vegetationMaskOutput(ip, jm)
            else
              k3 = 0
            end if
            
            if(landSeaMask(im, j) /= 0) then
              k4 = vegetationMaskOutput(im, j)
            else
              k4 = 0
            end if
            
            if(landSeaMask(ip, j) /= 0) then
              k5 = vegetationMaskOutput(ip, j)
            else
              k5 = 0
            end if
            
            if(landSeaMask(im, jp) /= 0) then
              k6 = vegetationMaskOutput(im, jp)
            else
              k6 = 0
            end if
            
            if(landSeaMask(i, jp) /= 0) then
              k7 = vegetationMaskOutput(i, jp)
            else
              k7 = 0
            end if
          
            if(landSeaMask(ip, jp) /= 0) then
              k8 = vegetationMaskOutput(ip, jp)
            else
              k8 = 0
            end if

            if(k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8 == 0) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if
          
            if(k1 /= 0) then
              kk = k1
            else
              kk = -1
            end if
          
            if(k2 /= 0 .and. kk == -1) then
              kk = k2
            else if(k2 /= 0 .and. kk /= k2) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if
          
            if(k3 /= 0 .and. kk == -1) then
              kk = k3
            else if(k3 /= 0 .and. kk /= k3) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if

            if(k4 /= 0 .and. kk == -1) then
              kk = k4
            else if(k4 /= 0 .and. kk /= k4) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if

            if(k5 /= 0 .and. kk == -1) then
              kk = k5
            else if(k5 /= 0 .and. kk /= k5) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if

            if(k6 /= 0 .and. kk == -1) then
              kk = k6
            else if(k6 /= 0 .and. kk /= k6) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if

            if(k7 /= 0 .and. kk == -1) then
              kk = k7
            else if(k7 /= 0 .and. kk /= k7) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if

            if(k8 /= 0 .and. kk == -1) then
              kk = k8
            else if(k8 /= 0 .and. kk /= k8) then
              call vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
              exit mainLoop
            end if

            vegetationMaskOutput(i, j) = kk
            if(dumpLocal) write(unit=*, fmt='(a,2i4,a,i2)') &
                              ' Undefined Location i, j = ', i, j, &
                              ' Filled with Nearby Value = ', kk
          end if

          if(landSeaMask(i, j) == 0) vegetationMaskOutput(i, j) = p_undefV
        end do
        ! Normally ended nested loops flag exit infinite loop
        if( i == varCommon%xMax + 1 .and. j == varCommon%yMax) flag = .true. 
      end do mainLoop
      if(flag) exit
    end do

  end subroutine vegetationMaskCheck


  subroutine vegetationMaskFix(i, j, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
    !# Vegetation Mask Fix
    !# ---
    !# @info
    !# **Brief:** Vegetation Mask Fix. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    implicit none
    integer, intent(in) :: i 
    !# i index to check
    integer, intent(in) :: j 
    !# j index to check
    integer, intent(in) :: k1
    integer, intent(in) :: k2
    integer, intent(in) :: k3
    integer, intent(in) :: k4
    integer, intent(in) :: k5
    integer, intent(in) :: k6
    integer, intent(in) :: k7
    integer, intent(in) :: k8
    logical, intent(in) :: dumpLocal
    
    if(dumpLocal) write(unit=*, fmt='(a,2i4,/,a,/,(11x,3i3))') &
                       ' At Land Point Value Undefined at i, j = ', i, j, &
                       ' With no Unambigous Nearby Values.  Local Area:', &
                       k1, k2, k3, k4, VegetationMaskOutput(i,j), k5, k6, k7, k8

    if(landSeaMask(i, j) == 1) then
      landSeaMask(i, j) = 0
      if(dumpLocal) write(unit=*, fmt='(a)') ' Land Sea Mask Changed from 1 to 0'
    else
      landSeaMask(i,j)=1
      if(dumpLocal) write(unit=*, fmt='(a)') ' Land Sea Mask Changed from 0 to 1'
    end if

  end subroutine vegetationMaskFix


END module Mod_VegetationMask

!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_PorceSandMaskIBIS </br></br>
!#
!# **Brief**: Module responsible for reading the PorceSandMaskIBISCLima.dat file
!# and the land/sea mask in the model resolution and making the intersection and
!# interpolation of the fields in the horizontal and vertical resolutions in 
!# which the model will be executed</br></br>
!# 
!# In the interpolation process, the processing algorithm defines a value of 60%
!# sand where there is no data in the regions of Antarctica and Greenland.</br></br>
!# 
!# 
!# **Files in:**
!#
!# &bull; pre/dataout/LandSeaMask.GZZZZZ  (Ex.: pre/dataout/LandSeaMask.G00450)
!# &bull; pre/dataout/PorceSandMaskIBISClima.dat
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; model/datain/PorceSandMaskIBIS.GZZZZZ
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ
!# &bull; pre/dataout/PorceSandMaskIBIS.GZZZZZ
!# &bull; pre/dataout/orceSandMaskIBIS.GZZZZZ.ctl
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
!#  <li>22-01-2020 - Denis Eiras       - version: 2.1.0</li>
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

module Mod_PorceSandMaskIBIS

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_AreaIntegerInterp, only: gLats, initAreaIntegerInterp, &
                                   doAreaIntegerInterp
  use Mod_Messages, only : msgInLineFormatOut, msgNewLine, msgOut, msgWarningOut

  implicit none

  public initPorceSandMaskIBIS
  public generatePorceSandMaskIBIS
  public getNamePorceSandMaskIBIS
  public shouldRunPorceSandMaskIBIS

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'
  include 'porceSand.h'

  !input variables
  type PorceSandMaskIBISNameListData
    integer :: xDim  !< Number of Longitudes
    integer :: yDim  !< Number of Latitudes
    integer :: layer  !
    character(len = maxPathLength) :: nameLSM='LandSeaMask'            !< Read  - Formatted Land Sea Mask
    character(len = maxPathLength) :: varName='PorceSandMaskIBISClima' !< Read  - Input
    character(len = maxPathLength) :: nameLSMSSiB='ModelLandSeaMask'   !< Write - Formatted Land Sea Mask (same as nameLSM)
    character(len = maxPathLength) :: varNameVeg='PorceSandMaskIBIS'   !< Write - Adjusted Interpolated PorceSandMaskIBIS
  end type PorceSandMaskIBISNameListData

  type(varCommonNameListData) :: varCommon
  type(PorceSandMaskIBISNameListData)    :: var
  namelist /PorceSandMaskIBISNameList/ var

  !internal variables
  character(len=7)   :: nLats='.G     '                  !< posfix land sea mask file
  character(len=10)  :: mskfmt = '(      i1)' !< format of landSeaMask file
  integer, dimension(:,:,:), allocatable :: porceSandMaskIBISIn     !< ???
  integer, dimension(:,:,:), allocatable :: porceSandMaskIBISOut    !< ???
  integer, dimension(:,:,:), allocatable :: vegMaskSave             !< ???
  integer, dimension(:,:),   allocatable :: landSeaMask             !< ???
  integer, dimension(:,:),   allocatable :: lSMaskSave              !< ???
  integer, dimension(:,:,:), allocatable :: porceSandMaskIBISInput  !< ???
  integer, dimension(:,:,:), allocatable :: porceSandMaskIBISOutput !< ???
  logical :: flagInput(5), flagOutput(5)

  character(len=*), parameter :: header = 'PorceSandMaskIBIS     | '

  contains

  ! ---------------------------------------------------------------------------
  !> @brief Returns PorceSandMaskIBIS Module Name
  !!
  !! @details Returns PorceSandMaskIBIS Module Name.\n
  !!
  !! @author Denis Eiras
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  function getNamePorceSandMaskIBIS() result(returnModuleName)
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName !< variable for store module name

    returnModuleName = "PorceSandMaskIBIS"
  end function getNamePorceSandMaskIBIS

  ! ---------------------------------------------------------------------------
  !> @brief Initialization of PorceSandMaskIBIS module
  !!
  !! @details Initialization of PorceSandMaskIBIS module, defined in PRE_run.nml
  !!
  !! @author Denis Eiras
  !!
  !! @date jan/2020
  ! ---------------------------------------------------------------------------
  subroutine initPorceSandMaskIBIS(nameListFileUnit, varCommon_)
    implicit none
    integer, intent(in) :: nameListFileUnit !< file unit of namelist PRE_run.nml
    type(varCommonNameListData), intent(in) :: varCommon_ !< variable of type varCommonNameListData for managing common variables at PRE_run.nml

    read(unit = nameListFileUnit, nml = PorceSandMaskIBISNameList)
    varCommon = varCommon_

    write (nLats(3:7), '(I5.5)') varCommon%yMax
    write (mskfmt(2:7), '(I6)') varCommon%xMax

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

  end subroutine initPorceSandMaskIBIS

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
  function shouldRunPorceSandMaskIBIS() result(shouldRun)
    implicit none
    logical :: shouldRun !< result variable for store if method should run

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunPorceSandMaskIBIS

  ! ---------------------------------------------------------------------------
  !> @brief Get PorceSandMaskIBIS Out Filename
  !!
  !! @details Get PorceSandMaskIBIS Unformatted Climatological PorcSand Mask file name
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
  !! @details Allocate matrixes based in namelist
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine allocateData()
    implicit none
    allocate(porceSandMaskIBISInput(var%xDim,var%yDim,var%layer))
    allocate(porceSandMaskIBISOutput(varCommon%xMax,varCommon%yMax,var%layer))
    allocate(porceSandMaskIBISIn(var%xDim,var%yDim,var%layer))
    allocate(porceSandMaskIBISOut(varCommon%xMax,varCommon%yMax,var%layer))
    allocate(vegMaskSave(varCommon%xMax,varCommon%yMax,var%layer))
    allocate(landSeaMask(varCommon%xMax,varCommon%yMax))
    allocate(lSMaskSave(varCommon%xMax,varCommon%yMax))
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
    deallocate(porceSandMaskIBISInput)
    deallocate(porceSandMaskIBISOutput)
    deallocate(porceSandMaskIBISIn)
    deallocate(porceSandMaskIBISOut)
    deallocate(vegMaskSave)
    deallocate(landSeaMask)
    deallocate(lSMaskSave)
  end subroutine deallocateData
  
  ! ---------------------------------------------------------------------------
  !> @brief Read Land Sea Mask
  !!
  !! @details read Land Sea Mask, input file
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function readlandSeaMask() result(isReturnOk)
    implicit none
    logical :: isReturnOk !< return value (if function was executed sucessfully)
    character (len = maxPathLength) :: inFilename !< input filename
    integer :: inFileUnit !< To Read Formatted Land Sea Mask

    isReturnOk = .false.
    inFilename = trim(varCommon%dirPreOut)//trim(var%nameLSM)//nLats
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return

    read(unit=inFileUnit) landSeaMask
    lSMaskSave = landSeaMask
    close(unit=inFileUnit)
    isReturnOk = .true.
  end function readlandSeaMask
  
  ! ---------------------------------------------------------------------------
  !> @brief Read PorceSandMaskIBIS
  !!
  !! @details read PorceSandMaskIBIS, input file
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function readPorceSandMaskIBIS() result(isReturnOk)
    implicit none
    logical :: isReturnOk !< return value (if function was executed sucessfully)
    character (len = maxPathLength) :: inFilename !< input filename
    integer :: inFileUnit !< To Read Unformatted Climatological Vegetation Mask
    integer :: inRecSize !< size of file record

    isReturnOk = .false.
    inquire(iolength=inRecSize) porceSandMaskIBISIn
    inFilename = trim(varCommon%dirPreOut)//trim(var%varName)//'.dat'
    inFileUnit = openFile(trim(inFilename), 'unformatted', 'direct', inRecSize, 'read', 'old')
    if(inFileUnit < 0) return
    read(unit=inFileUnit, rec=1) porceSandMaskIBISIn
    close(unit=inFileUnit)
    isReturnOk = .true.
  end function readPorceSandMaskIBIS

  ! ---------------------------------------------------------------------------
  !> @brief write Interpolated PorceSandMaskIBIS
  !!
  !! @details Write Out Adjusted Interpolated PorceSandMaskIBIS
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function writeInterpolatedPorceSandMaskIBIS() result(isReturnOk)
    implicit none
    logical :: isReturnOk !< return value (if function was executed sucessfully)
    character (len = maxPathLength) :: outFilename !< output Adjusted Interpolated PorceSandMaskIBIS
    integer :: outFileUnit !< out file unit
    integer :: lRecOut !< size of file record

    isReturnOk = .false.
    inquire(iolength=lRecOut) porceSandMaskIBISOut
    outFilename = trim(varCommon%dirModelIn)//trim(var%varNameVeg)//nLats
    outFileUnit = openFile(trim(outFilename), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(outFileUnit < 0) return

    write(unit=outFileUnit, rec=1) porceSandMaskIBISOut
    close(unit=outFileUnit)
    isReturnOk = .true.
  end function writeInterpolatedPorceSandMaskIBIS
  
  ! ---------------------------------------------------------------------------
  !> @brief Writes out Formatted Land Sea Mask Modified by Vegetation
  !!
  !! @details Writes out Formatted Land Sea Mask Modified by Vegetation
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function writeLandSeaMask() result(isReturnOk)
    implicit none
    logical :: isReturnOk !< return value (if function was executed sucessfully)
    character (len = maxPathLength) :: outFilename !< output file name of Formatted Land Sea Mask Modified by Vegetation
    integer :: outFileUnit !< out file unit

    isReturnOk = .false.
    outFilename = trim(varCommon%dirPreOut)//trim(var%nameLSMSSiB)//nLats
    outFileUnit = openFile(trim(outFilename), 'unformatted', 'sequential', -1, 'write', 'replace')
    if(outFileUnit < 0) return
    write(unit=outFileUnit) landSeaMask
    close(unit=outFileUnit)
    isReturnOk = .true.
  end function writeLandSeaMask
  
  ! ---------------------------------------------------------------------------
  !> @brief Generate Grads files 
  !!
  !! @details Generate .ctl and .dat files to check output
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function generateGrads() result(isReturnOk)
    implicit none
    logical :: isReturnOk !< return value (if function was executed sucessfully)
    integer :: lRecGad !< size of file record
    integer :: k       !< temporary variable
    character (len = maxPathLength) :: gradsBaseName !< dat GrADS base name without extension
    character (len = maxPathLength) :: gradsPathFileName !< dat Grads File name with path
    character (len = maxPathLength) :: ctlPathFileName !< ctl Grads File name with path
    integer :: gradsDatFileUnit, gradsCtlFileUnit !< File Units for grads files
 
    real(kind=p_r4), dimension(:,:,:), allocatable :: gad !< ???
    real(kind=p_r4), parameter :: p_undefG = -99.0      !< ???

    isReturnOk = .false.
    gradsBaseName = trim(var%varNameVeg) // nlats
    gradsPathFileName = trim(varCommon%dirPreOut) // trim(gradsBaseName)
    ctlPathFileName   = trim(varCommon%dirPreOut) // trim(gradsBaseName) // '.ctl'
    allocate(gad(varCommon%xMax,varCommon%yMax,var%layer))

    ! Write .dat file --------------------------------------------------------
    inquire(iolength=lRecGad) gad(:,:,1)
    gradsDatFileUnit = openFile(trim(gradsPathFileName), 'unformatted', 'direct', lRecGad, 'write', 'replace')
    if(gradsDatFileUnit < 0) return

    gad(:,:,1) = real(lSMaskSave,p_r4)
    write(unit=gradsDatFileUnit, rec=1) gad(:,:,1)
    gad(:,:,1) = real(landSeaMask,p_r4)
    write(unit=gradsDatFileUnit, rec=2) gad(:,:,1)
    do k=1,var%layer
      gad(:,:,k) = real(vegMaskSave(:,:,k), p_r4)
      write(unit=gradsDatFileUnit, rec=2+k) gad(:,:,k)
    end do 
    do k=1,var%layer
      gad(:,:,k) = real(porceSandMaskIBISOut(:,:,k),p_r4)
      write(unit=gradsDatFileUnit, rec=2+var%layer+k) gad(:,:,k)
    end do
    close (unit=gradsDatFileUnit)
    ! -------------------------------------------------------------------------
    
    ! Write GrADS Control File ---------------------------------------------------------
    gradsCtlFileUnit = openFile(trim(ctlPathFileName), 'formatted', 'sequential', -1, 'write', 'replace')
    if(gradsCtlFileUnit < 0) return

    write(unit=gradsCtlFileUnit, fmt='(a)') 'DSET ^' // trim(gradsBaseName)
    write(unit=gradsCtlFileUnit, fmt='(a)') '*'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'OPTIONS YREV BIG_ENDIAN'
    write(unit=gradsCtlFileUnit, fmt='(a)') '*'
    write(unit=gradsCtlFileUnit, fmt='(a,1pg12.5)') 'UNDEF ', p_undefG
    write(unit=gradsCtlFileUnit, fmt='(a)') '*'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'TITLE Vegetation Mask on a Gaussian Grid'
    write(unit=gradsCtlFileUnit, fmt='(a)') '*'
    write(unit=gradsCtlFileUnit, fmt='(a,i5,a,f8.3,f15.10)') 'XDEF ',varCommon%xMax,' LINEAR ', &
               0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
    write(unit=gradsCtlFileUnit, fmt='(a,i5,a)') 'YDEF ',varCommon%yMax,' LEVELS '
    write(unit=gradsCtlFileUnit, fmt='(8f10.5)') gLats(varCommon%yMax:1:-1)
    write(unit=gradsCtlFileUnit, fmt='(a)') 'ZDEF  6 LEVELS 0 1 2 3 4 5'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'TDEF  1 LINEAR JAN2005 1MO'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'VARS  4'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'LSMO  0 99 Land Sea Mask Before Fix  [0-Sea 1-Land]'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'LSMK  0 99 Land Sea Mask For Model   [0-Sea 1-Land]'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'VGMO  6 99 Vegetation Mask Before Fix [0 to 13]'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'VEGM  6 99 Vegetation Mask For Model  [0 to 13]'
    write(unit=gradsCtlFileUnit, fmt='(a)') 'ENDVARS'
    close(unit=gradsCtlFileUnit)
    ! -------------------------------------------------------------------------------------
    deallocate(gad)
    isReturnOk = .true.
  end function generateGrads
  
  ! ---------------------------------------------------------------------------
  !> @brief Generate PorceSandMaskIBIS
  !!
  !! @details Generate PorceSandMaskIBIS output. \n
  !! This subroutine is the main method for use this module. 
  !! Only file name of namelist is needed to use it.
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  function generatePorceSandMaskIBIS() result(isReturnOk)
    implicit none
    logical :: isReturnOk !< return value (if function was executed sucessfully)
    logical, parameter :: dumpLocal=.false.
    integer :: i
    integer :: j
    integer :: k

    isReturnOk = .false.
    call allocateData()
    if (.not. readlandSeaMask() ) then
      call msgWarningOut(header, "Error reading LandSeaMask file")
      return
    end if
    if (.not. readPorceSandMaskIBIS() ) then
      call msgWarningOut(header, "Error reading PorceSandMaskIBISClima file")
      return
    end if

    call initAreaIntegerInterp(var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_numVegClasses, vegClasses, flagInput, flagOutput)

    ! Interpolate Input Regular PorceSandMaskIBIS To Gaussian Grid Output
    porceSandMaskIBISInput = int(porceSandMaskIBISIn)
    do k=1, var%layer
      if (.not. doAreaIntegerInterp(porceSandMaskIBISInput(:,:,k), porceSandMaskIBISOutput(:,:,k))) return
    end do

    do k=1, var%layer
      do j=1, varCommon%yMax
        do i=1, varCommon%xMax
          if(lSMaskSave(i,j) /= 0 .and.porceSandMaskIBISOutput(i,j,k) == 0.0) then
            porceSandMaskIBISOutput(i,j,k) = 60
          end if
        end do
      end do
    end do

    vegMaskSave = int(porceSandMaskIBISOutput)
   
    ! Fix Problems on Vegetation Mask and Land Sea Mask
    call vegetationMaskCheck(dumpLocal)
    porceSandMaskIBISOut = int(porceSandMaskIBISOutput)

    if (.not. writeInterpolatedPorceSandMaskIBIS() ) then
      call msgWarningOut(header, "Error writing PorceSandMaskIBIS file")
      return
    end if
    if (.not. writeLandSeaMask() ) then
      call msgWarningOut(header, "Error writing ModelLandSeaMask file")
      return
    end if

    ! Generate grads output for debug or test
    if(varCommon%grads .and. .not. generateGrads()) then
      call msgWarningOut(header, "Error while generating grads files")
      return
    end if

    call deallocateData()
    isReturnOk = .true.
  end function generatePorceSandMaskIBIS

  ! ---------------------------------------------------------------------------
  !> @brief Vegetation Mask Check 
  !!
  !! @details ???
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine vegetationMaskCheck(dumpLocal)
    logical, intent(in) :: dumpLocal !< print debug information
    
    ! local variables
    logical :: flag=.false. !< exit loop control
    integer :: i            !<
    integer :: im           !<
    integer :: ip           !<
    integer :: j            !<
    integer :: jm           !<
    integer :: jp           !<
    integer :: k            !< 
    integer :: k1           !<
    integer :: k2           !<
    integer :: k3           !<
    integer :: k4           !<
    integer :: k5           !<
    integer :: k6           !<
    integer :: k7           !<
    integer :: k8           !<
    integer :: kk           !<

    integer :: p_undef_int_zero = 0
  
    do k=1, var%layer
      im=0
      ip=0
      jm=0
      jp=0
      k1=0
      k2=0
      k3=0
      k4=0
      k5=0
      k6=0
      k7=0
      k8=0
      kk=0
      ! infinite loop
      do 
        mainLoop: do j=1, varCommon%yMax
          do i=1, varCommon%xMax
            if(landSeaMask(i,j) /= 0 .and. porceSandMaskIBISOutput(i,j,k) == p_undef_int_zero) then
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
              
              im = max(min(varCommon%xMax,im),1)
              jm = max(min(varCommon%yMax,jm),1)
              if(landSeaMask(im,jm) /= 0) then
                k1 = porceSandMaskIBISOutput(im,jm,k)
              else
                k1 = 0
              end if
              if(landSeaMask(i,jm) /= 0) then
                k2 = porceSandMaskIBISOutput(i,jm,k)
              else
                k2=0
              end if
              
              ip = max(min(varCommon%xMax,ip),1)
              if(landSeaMask(ip,jm) /= 0) then
                k3 = porceSandMaskIBISOutput(ip,jm,k)
              else
                k3 = 0
              end if
              if(landSeaMask(im,j) /= 0) then
                k4 = porceSandMaskIBISOutput(im,j,k)
              else
                k4 = 0
              end if
              if(landSeaMask(ip,j) /= 0) then
                k5 = porceSandMaskIBISOutput(ip,j,k)
              else
                k5 = 0
              end if
               
              jp = max(min(varCommon%yMax,jp),1)
              if(landSeaMask(im,jp) /= 0) then
                k6 = porceSandMaskIBISOutput(im,jp,k)
              else
                k6 = 0
              end if
              if(landSeaMask(i,jp) /= 0) then
                k7 = porceSandMaskIBISOutput(i,jp,k)
              else
                k7 = 0
              end if
              if(landSeaMask(ip,jp) /= 0) then
                k8 = porceSandMaskIBISOutput(ip,jp,k)
              else
                k8 = 0
              end if

              if(k1+k2+k3+k4+k5+k6+k7+k8 == 0) then
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
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
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
                exit mainLoop
              end if

              if(k3 /= 0 .and. kk == -1) then
                kk = k3
              else if (k3 /= 0 .and. kk /= k3) then
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
                exit mainLoop
              end if

              if(k4 /= 0 .and. kk == -1) then
                kk = k4
              else if(k4 /= 0 .and. kk /= k4) then
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
                exit mainLoop
              end if

              if(k5 /= 0 .and. kk == -1) then
                kk = k5
              else if(k5 /= 0 .and. kk /= k5) then
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
                exit mainLoop
              end if

              if(k6 /= 0 .and. kk == -1) then
                kk = k6
              else if(k6 /= 0 .and. kk /= k6) then
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
                exit mainLoop
              end if

              if(k7 /= 0 .and. kk == -1) then
                kk = k7
              else if (k7 /= 0 .and. kk /= k7) then
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
                exit mainLoop
              end if

              if(k8 /= 0 .and. kk == -1) then
                kk = k8
              else if(k8 /= 0 .and. kk /= k8) then
                call vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
                exit mainLoop
              end if

              porceSandMaskIBISOutput(i,j,k) = kk
              if(dumpLocal) then
                call msgInLineFormatOut(' ' // header // ' Undefined Location i, j = (', '(A)')
                call msgInLineFormatOut(i, '(i4)')
                call msgInLineFormatOut(', ', '(A)')
                call msgInLineFormatOut(j, '(i4)')
                call msgInLineFormatOut(') Filled with Nearby Value = ', '(A)')
                call msgInLineFormatOut(kk, '(i4)')
                call msgNewLine();
              end if
            end if

            if (landSeaMask(i,j) == 0) porceSandMaskIBISOutput(i,j,k) = p_undef_int_zero
          end do ! end i

          ! Normally ended nested loops flag exit infinite loop
          if( i == varCommon%xMax + 1 .and. j == varCommon%yMax) flag = .true.
          
        end do mainLoop
        if(flag) exit
      end do ! infinite loop
    end do ! end k
    
end subroutine vegetationMaskCheck


  ! ---------------------------------------------------------------------------
  !> @brief Vegetation Mask Fix 
  !!
  !! @details ???
  !!
  !! @author Daniel M. Lamosa
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine vegetationMaskFix(i, j, k, k1, k2, k3, k4 ,k5, k6, k7, k8, dumpLocal)
    integer, intent(in) :: i !< i index to check
    integer, intent(in) :: j !< j index to check
    integer, intent(in) :: k !< k index to check
    integer, intent(in) :: k1
    integer, intent(in) :: k2
    integer, intent(in) :: k3
    integer, intent(in) :: k4
    integer, intent(in) :: k5
    integer, intent(in) :: k6
    integer, intent(in) :: k7
    integer, intent(in) :: k8
    logical, intent(in) :: dumpLocal
    
    if(dumpLocal) then
      call msgInLineFormatOut(' ' // header // ' At Land Point Value Undefined at i, j, k = ', '(A)')
      call msgInLineFormatOut(i, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(j, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k, '(i8)')
      call msgNewLine();
      call msgInLineFormatOut(' ' // header // ' With no Unambigous Nearby Values.  Local Area: ', '(A)')

      call msgInLineFormatOut(k1, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k2, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k3, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k4, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(porceSandMaskIBISOutput(i,j,k), '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k5, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k6, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k7, '(i8)')
      call msgInLineFormatOut(', ', '(A)')
      call msgInLineFormatOut(k8, '(i8)')
      call msgNewLine();
    end if
    
    if(landSeaMask(i, j) == 1) then
      landSeaMask(i, j) = 0
      if(dumpLocal) call msgOut(header, ' Land Sea Mask Changed from 1 to 0')
    else
      landSeaMask(i,j) = 1
      if(dumpLocal) call msgOut(header, ' Land Sea Mask Changed from 0 to 1')
    end if

  end subroutine vegetationMaskFix

end module Mod_PorceSandMaskIBIS

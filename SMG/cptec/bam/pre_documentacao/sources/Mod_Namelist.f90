!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Namelist </br></br>
!#
!# **Brief**: Module used for file manipulating namelists defined in 
!# PRE_run.nml and PRE_cfg.nml </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!#
!# **Author**: Denis Eiras </br>
!#
!# **Version**: 1.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>26-04-2019 - Denis Eiras  -  version: 1.0.0 </li>
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

module Mod_Namelist

  use Mod_FileManager, only: &
    openFileErrorMsg &
    , openFile

  private
  include "files.h"

  public :: oneConfigFile
  public :: createConfigFile
  public :: destroyConfigFile
  public :: getConfigFileName
  public :: openConfigFile
  public :: closeConfigFile

  public :: oneNamelistFile
  public :: createNamelistFile
  public :: destroyNamelistFile
  public :: getNamelistFileName
  public :: openNamelistFile
  public :: closeNamelistFile

  public :: varCommonNameListData
  type namelistFile
    character(len = maxPathLength) :: fileName
  end type namelistFile

  type(namelistFile), pointer :: oneConfigFile
  type(namelistFile), pointer :: oneNamelistFile
  ! variables common to all modules
  type varCommonNameListData
    integer                        :: xMax                     
    !# max value of longitude
    integer                        :: yMax                     
    !# max value of latitude
    integer                        :: mEnd
    character(len = maxPathLength) :: dirPreIn      = './'     
    !# input data directory
    character(len = maxPathLength) :: dirPreOut     = './'     
    !# output data directory
    character(len = maxPathLength) :: dirBCs        = './'     
    !# Boundary Conditions directory
    character(len = maxPathLength) :: dirClmCO2     = './'     
    character(len = maxPathLength) :: dirClmSST     = './'     
    character(len = maxPathLength) :: dirObsSST     = './'     
    !# Climatological SST Datain Directory
    character(len = maxPathLength) :: dirClmOCM     = './'     
    !# Climatological OCM Datain Directory
    character(len = maxPathLength) :: dirClmFluxCO2 = './'     
    !# Climatological CO2 Datain Directory
    character(len = maxPathLength) :: dirModelIn    = './'
    character(len = 10)            :: date
    logical                        :: grads                    
    !# flag for grads output (true = yes, false = no)
  end type varCommonNameListData


contains


  subroutine createConfigFile()
    !# Creates pointer for type namelistfile configuration file PRE_cfg.nml
    !# ---
    !# @info
    !# **Brief:** Creates pointer for type namelistfile configuration file PRE_cfg.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    if (associated(oneConfigFile)) then
      deallocate(oneConfigFile)
    end if
    allocate(oneConfigFile)
  end subroutine createConfigFile


  subroutine destroyConfigFile()
    !# Deallocates pointer for type namelistfile configuration file PRE_cfg.nml
    !# ---
    !# @info
    !# **Brief:** Deallocates pointer for type namelistfile configuration file PRE_cfg.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    if (associated(oneConfigFile)) then
      deallocate(oneConfigFile)
    end if
    nullify(oneConfigFile)
  end subroutine destroyConfigFile


  function getConfigFileName() result(fileName)
    !# Gets filename of namelistfile configuration file PRE_cfg.nml
    !# ---
    !# @info
    !# **Brief:** Gets filename of namelistfile configuration file PRE_cfg.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none

    integer :: nargs
    integer :: iarg
    integer :: lenArg
    integer :: status
    logical :: flagName ! true iff arg="-f"; next arg is file name
    character(len = maxPathLength) :: arg
    character(len = maxPathLength) :: fileName

    oneConfigFile%fileName = "PRE_cfg.nml" ! default namelist
    ! search command line for "-f " <namelist file name>
    ! return default if not found
    nargs = command_argument_count()
    if (nargs >= 0) then
      flagName = .false.
      do iarg = 0, nargs
        call get_command_argument(iarg, arg, lenArg, status)
        if (status == 0) then
          if (flagName) then
            oneConfigFile%fileName = arg(1:lenArg)
            exit
          else
            flagName = arg(1:lenArg) == "-c"
          end if
        end if
      end do
    end if
    fileName = oneConfigFile%fileName
  end function getConfigFileName


  function openConfigFile() result(fu_namelist)
    !# Opens configuration file PRE_cfg.nml
    !# ---
    !# @info
    !# **Brief:** Opens configuration file PRE_cfg.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    integer :: fu_namelist

    fu_namelist = openFile(trim(oneConfigFile%fileName), 'formatted', 'sequential', -1, 'read', 'old')

  end function openConfigFile


  subroutine closeConfigFile(fu_namelist)
    !# Closes configuration file PRE_cfg.nml
    !# ---
    !# @info
    !# **Brief:** Closes configuration file PRE_cfg.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    integer, intent(in) :: fu_namelist
    integer :: err

    close(fu_namelist, iostat = err)
    if (err /= 0) then
      call openFileErrorMsg(err, oneConfigFile%fileName)
    end if

  end subroutine closeConfigFile


  subroutine createNamelistFile()
    !# Creates pointer for type namelistfile configuration file PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Creates pointer for type namelistfile configuration file PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    if (associated(oneNamelistFile)) then
      deallocate(oneNamelistFile)
    end if
    allocate(oneNamelistFile)
  end subroutine createNamelistFile


  subroutine destroyNamelistFile()
    !# Deallocates pointer for type namelistfile configuration file PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Deallocates pointer for type namelistfile configuration file 
    !# PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    if (associated(oneNamelistFile)) then
      deallocate(oneNamelistFile)
    end if
    nullify(oneNamelistFile)
  end subroutine destroyNamelistFile


  subroutine setNamelistFileName(fileName)
    !# Sets filename of namelistfile configuration file PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Sets filename of namelistfile configuration file PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength), intent(in) :: fileName

    oneNamelistFile%fileName = fileName
  end subroutine setNamelistFileName


  function getNamelistFileName() result(fileName)
    !# Gets filename of namelistfile configuration file PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Gets filename of namelistfile configuration file PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none

    integer :: nargs
    integer :: iarg
    integer :: lenArg
    integer :: status
    logical :: flagName ! true iff arg="-f"; next arg is file name
    character(len = maxPathLength) :: arg
    character(len = maxPathLength) :: fileName

    oneNamelistFile%fileName = "PRE_run.nml" ! default namelist
    ! search command line for "-f " <namelist file name>
    ! return default if not found
    nargs = command_argument_count()
    if (nargs >= 0) then
      flagName = .false.
      do iarg = 0, nargs
        call get_command_argument(iarg, arg, lenArg, status)
        if (status == 0) then
          if (flagName) then
            oneNamelistFile%fileName = arg(1:lenArg)
            exit
          else
            flagName = arg(1:lenArg) == "-f"
          end if
        end if
      end do
    end if
    fileName = oneNamelistFile%fileName
  end function getNamelistFileName


  function openNamelistFile() result(fu_namelist)    
    !# Opens run configuration file PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Opens run configuration file PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    integer :: fu_namelist

    fu_namelist = openFile(trim(oneNamelistFile%fileName), 'formatted', 'sequential', -1, 'read', 'old')

  end function openNamelistFile


  subroutine closeNamelistFile(fu_namelist)
    !# Closes run configuration file PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Closes run configuration file PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo

    integer, intent(in) :: fu_namelist
    integer :: err

    close(fu_namelist, iostat = err)
    if (err /= 0) then
      call openFileErrorMsg(err, oneNamelistFile%fileName)
    end if

  end subroutine closeNamelistFile


end module Mod_Namelist

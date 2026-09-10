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

  use Mod_String_Functions, only : &
    intToStr

  use Mod_Messages_Parallel, only : &
    msgOutMaster

  use Mod_Parallelism, only : &
    fatalError

  private
  include "files.h"
  include "pre.h"

  character(len = *), parameter :: headerMsg = 'Namelist Module  | '

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

  ! ModulesCommonNameList - values common to all modules
  type(varCommonNameListData), pointer  ::  varCommon
  namelist /ModulesCommonNameList/ varCommon


  integer :: fu_namelist, fu_config
  character(len = maxPathLength) :: nameListFileName, configFileName

  ! RunNameList - namelist which contains the list of modules that will run
  ! max modules to run, each with max characters in namelist name, defined in pre.h
  integer :: numberOfModulesToRun ! size of modulesToRun
  character (len = maxModuleNameLength), dimension(maxNumberOfModules) :: modulesToRun = ""
  namelist /RunNameList/ modulesToRun

  ! ModulesNameList
  ! List of all modules available, containing:
  ! - moduleName: name fo the module;
  ! - dependencies: list of modules which module depends. Each module aims to configure the modules  which it depends.
  !     Then, the dependencies will run automatically before the module declared as element, in cases:
  !     - (1) The dependency is configured to run in namelist "RunNameList" (PRE_run.nml). Always runs;
  !     - (2) The output files of the dependency does not exists. Will run even wheter is not configured to run in
  !       RunNameList.
  ! - processors: number of processors that module runs. Normally equals 1. Only Chopping accepts more than 1, until now.
  integer :: numberOfModules ! size of all modules
  type ModuleConfig
    character(len = maxModuleNameLength) :: moduleName = ""
    character(len = maxModuleNameLength) :: dependencies(maxNumberOfModules - 1) = ""
    integer :: processors = 1
    logical :: runsAlone = .false.
  end type ModuleConfig
  type(ModuleConfig), pointer :: modCfg(:)
  namelist /ModulesNameList/ modCfg

  ! PUBLIC section

  public :: createConfigFile
  public :: createNamelistFile
  public :: setConfigFileName
  public :: setNamelistFileName
  public :: initializeNamelists
  public :: initializeNamelistsDefaultFileNames
  public :: getNumberOfModules
  public :: numberOfModulesToRun
  !# getter and setter for Module Namelist (module manipulator)
  public :: modulesToRun
  !# getter for Module Pre (cannot manipulate modulesToRun)
  public :: getModuleToRun
  public :: getVarCommon
  public :: ModuleConfig
  public :: getModCfg
  public :: getModCfgAll
  public :: varCommonNameListData
  public :: getModuleIndex
  public :: get_fu_namelist
  public :: finalizeNamelists

contains


  subroutine initializeNamelistsDefaultFileNames()
    !# Initialize namelists from PRE_cfg.nml and PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Initialize namelists from PRE_cfg.nml and PRE_run.nml using default names or command line options names </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: aug/2021 </br>
    !# @endinfo
    implicit none

    call createConfigFile()
    configFileName = getConfigFileName()
    call createNamelistFile()
    nameListFileName = getNamelistFileName()

    call initializeNamelists()

  end subroutine initializeNamelistsDefaultFileNames


  subroutine initializeNamelists()
    !# Initialize namelists from PRE_cfg.nml and PRE_run.nml
    !# ---
    !# @info
    !# **Brief:** Initialize namelists from PRE_cfg.nml and PRE_run.nml </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: aug/2021 </br>
    !# @endinfo

    implicit none

    integer :: ioError, moduleIndex, runIdx
    character (len = maxModuleNameLength) :: moduleName_aux, moduleNameWithProcs

    allocate(modCfg(maxNumberOfModules))
    allocate(varCommon)

    call msgOutMaster("$" // headerMsg, "Reading Modules Configuration File ...")
    fu_config = openConfigFile()
    if(fu_config < 0) call fatalError(headerMsg, "Configuration Namelist file does not exists: " // configFileName)

    ! calculates numberOfModules through ModulesNameList modulesNames array size
    read(unit = fu_config, nml = ModulesNameList, iostat=ioError)
    if (ioError < 0) then
      call fatalError(headerMsg, 'Error reading ModulesNameList in PRE_cfg.nml')
    endif
    call msgOutMaster("$" // headerMsg, "Modules implemented (moduleName, processors):")
    do moduleIndex = 1, maxNumberOfModules
      moduleName_aux = modCfg(moduleIndex)%moduleName
      if(len_trim(moduleName_aux) .eq. 0) then
        numberOfModules = moduleIndex - 1
        exit
      endif
      if(modCfg(moduleIndex)%processors .eq. 0) then
        modCfg(moduleIndex)%processors = 1
      endif
      call msgOutMaster("", " => " // trim(moduleName_aux) // ", " // intToStr(modCfg(moduleIndex)%processors))
    enddo

    call msgOutMaster("$" // headerMsg, "Reading Namelist File ...")
    modulesToRun(:) = ""
    fu_namelist = openNamelistFile()
    if(fu_namelist < 0) call fatalError(headerMsg, "Run Namelist file does not exists: " // nameListFileName)

    ! calculates numberOfModulesToRun
    read(unit = fu_namelist, nml = RunNameList, iostat=ioError)
    if (ioError < 0) then
      call fatalError(headerMsg, 'Error reading RunNameList in PRE_run.nml')
    endif
    do runIdx = 1, maxNumberOfModules
      moduleName_aux = trim(modulesToRun(runIdx))
      if(len(trim(moduleName_aux)) .eq. 0) then
        numberOfModulesToRun = runIdx - 1
        exit
      endif
    enddo

    if(numberOfModulesToRun .eq. 0) then
      call fatalError("$" // headerMsg, "None modules configured to run in RunNameList. Exiting ...")
    endif


    ! reads Modules common namelist
    read(unit = fu_namelist, nml = ModulesCommonNameList, iostat=ioError)
    if (ioError < 0) then
      call fatalError(headerMsg, 'Error reading ModulesCommonNameList in PRE_run.nml')
    endif

  end subroutine initializeNamelists


  function get_fu_namelist()
    implicit none
    integer :: get_fu_namelist
    get_fu_namelist = fu_namelist
  end function get_fu_namelist


  function getVarCommon()
    implicit none
    type(varCommonNameListData), pointer  ::  getVarCommon

    if (associated(varCommon)) then
      allocate(getVarCommon)
      getVarCommon = varCommon
    else
      getVarCommon => NULL()
    end if
  end function getVarCommon


  function getModCfg(modIndex)
    implicit none
    integer, intent(in) :: modIndex
    type(ModuleConfig) :: getModCfg

    getModCfg = modCfg(modIndex)
  end function getModCfg


  subroutine getModCfgAll(modCfgAll)
    implicit none
    type(ModuleConfig), pointer, intent(out) :: modCfgAll(:)

    if (associated(modCfg)) then
      allocate(modCfgAll(maxNumberOfModules))
      modCfgAll = modCfg
    else
      modCfgAll => NULL()
    end if
  end subroutine getModCfgAll


  function getModuleToRun(idx)
    implicit none
    integer, intent(in) :: idx
    character (len = maxModuleNameLength) :: getModuleToRun
    getModuleToRun = modulesToRun(idx)
  end function getModuleToRun


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


  function getNumberOfModules() result(numOfModules)
    implicit none
    integer :: numOfModules
    numOfModules = numberOfModules

  end function getNumberOfModules


  function getModuleIndex(moduleName) result(moduleIndex)
    !# Gets module index, given module name
    !# ---
    !# @info
    !# **Brief:** Gets module index, given module name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 </br>
    !# @endinfo
    implicit none
    character(len = *), intent(in) :: moduleName
    integer :: moduleIndex

    do moduleIndex = 1, numberOfModules
      if(trim(moduleName) .eq. trim(modCfg(moduleIndex)%moduleName)) return
    enddo
    moduleIndex = - 1
  end function getModuleIndex


  subroutine setConfigFileName(fileName)
    !# Sets filename of namelistfile configuration file PRE_cfg.nml
    !# ---
    !# @info
    !# **Brief:** Sets filename of namelistfile configuration file PRE_cfg.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: aug/2021 </br>
    !# @endinfo
    implicit none
    character(len = *), intent(in) :: fileName

    oneConfigFile%fileName = trim(fileName)

  end subroutine setConfigFileName


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
    character(len = *), intent(in) :: fileName

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


  subroutine finalizeNamelists()
    call msgOutMaster(headerMsg, "Closing Namelist File ...")
    call closeNamelistFile(fu_namelist)
    call destroyNamelistFile()
    call msgOutMaster(headerMsg, "Closing Config File ...")
    call closeNamelistFile(fu_config)
    call destroyConfigFile()
    numberOfModules = 0
    numberOfModulesToRun = 0
    modulesToRun(:) = ""
    deallocate(varCommon)
    deallocate(modCfg)
  end subroutine finalizeNamelists


end module Mod_Namelist

!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Pre </br></br>
!#
!# **Brief**: This is the main module, which invokes all the pre modules</br></br>
!# 
!# **Author**: Denis Eiras </br>
!#
!# **Version**: 2.0.4 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>26-04-2019 - Denis Eiras    - version: 1.0.0</li></br>
!#  <li>20-05-2019 - Eduardo Khamis - version: 2.0.0</li></br>
!#  <li>17-06-2019 - Eduardo Khamis - version: 2.0.1 - adding modules</li></br>
!#  <li>02-07-2019 - Eduardo Khamis - version: 2.0.2 - adding SSTClima</li></br>
!#  <li>23-01-2020 - Eduardo Khamis - version: 2.0.3 - adding Vegetation modules</li></br>
!#  <li>27-01-2020 - Eduardo Khamis - version: 2.0.4 - adding DeepSoilTemperature* and RoughnessLength* modules</li></br></br>
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
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2></br>
!# **&#9993;**<mailto:atende.cptec@inpe.br></br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---

module Mod_Pre

  use Mod_Albedo, only : &
    generateAlbedo &
    , getNameAlbedo &
    , initAlbedo &
    , shouldRunAlbedo

  use Mod_AlbedoClima, only : &
    generateAlbedoClima &
    , getNameAlbedoClima &
    , initAlbedoClima &
    , shouldRunAlbedoClima

  use Mod_Chopping, only : &
    generateChopping &
    , getNameChopping &
    , initChopping &
    , shouldRunChopping

  use Mod_CO2MonthlyDirec, only : &
    generateCO2MonthlyDirec &
    , getNameCO2MonthlyDirec &
    , initCO2MonthlyDirec &
    , shouldRunCO2MonthlyDirec

  use Mod_DeepSoilTemperature, only : &
    generateDeepSoilTemperature &
    , getNameDeepSoilTemperature &
    , initDeepSoilTemperature &
    , shouldRunDeepSoilTemperature

  use Mod_DeepSoilTemperatureClima, only : &
    generateDeepSoilTemperatureClima &
    , getNameDeepSoilTemperatureClima &
    , initDeepSoilTemperatureClima &
    , shouldRunDeepSoilTemperatureClima

  use Mod_DeltaTempColdest, only : &
    generateDeltaTempColdest &
    , getNameDeltaTempColdest &
    , initDeltaTempColdest &
    , shouldRunDeltaTempColdest
  
  use Mod_DeltaTempColdestClima, only : &
    generateDeltaTempColdestClima &
    , getNameDeltaTempColdestClima &
    , initDeltaTempColdestClima &
    , shouldRunDeltaTempColdestClima

  use Mod_FLUXCO2Clima, only : &
    generateFLUXCO2Clima &
    , getNameFLUXCO2Clima &
    , initFLUXCO2Clima &
    , shouldRunFLUXCO2Clima

  use Mod_LandSeaMask, only : &
    generateLandSeaMask &
    , getNameLandSeaMask &
    , initLandSeaMask &
    , shouldRunLandSeaMask

  use Mod_OCMClima, only : &
    generateOCMClima &
    , getNameOCMClima &
    , initOCMClima &
    , shouldRunOCMClima

  use Mod_PorceClayMaskIBIS, only : &
    generatePorceClayMaskIBIS &
    , getNamePorceClayMaskIBIS &
    , initPorceClayMaskIBIS &
    , shouldRunPorceClayMaskIBIS

  use Mod_PorceClayMaskIBISClima, only : &
    generatePorceClayMaskIBISClima &
    , getNamePorceClayMaskIBISClima &
    , initPorceClayMaskIBISClima &
    , shouldRunPorceClayMaskIBISClima

  use Mod_PorceSandMaskIBIS, only : &
    generatePorceSandMaskIBIS &
    , getNamePorceSandMaskIBIS &
    , initPorceSandMaskIBIS &
    , shouldRunPorceSandMaskIBIS

  use Mod_PorceSandMaskIBISClima, only : &
    generatePorceSandMaskIBISClima &
    , getNamePorceSandMaskIBISClima &
    , initPorceSandMaskIBISClima &
    , shouldRunPorceSandMaskIBISClima

  use Mod_RoughnessLength, only : &
    generateRoughnessLength &
    , getNameRoughnessLength &
    , initRoughnessLength &
    , shouldRunRoughnessLength

  use Mod_RoughnessLengthClima, only : &
    generateRoughnessLengthClima &
    , getNameRoughnessLengthClima &
    , initRoughnessLengthClima &
    , shouldRunRoughnessLengthClima

  use Mod_SnowClima, only : &
    generateSnowClima &
    , getNameSnowClima &
    , initSnowClima &
    , shouldRunSnowClima

  use Mod_SnowWeeklyNCEP, only : &
    generateSnowWeeklyNCEP &
    , getNameSnowWeeklyNCEP &
    , initSnowWeeklyNCEP &
    , shouldRunSnowWeeklyNCEP

  use Mod_SoilMoisture, only : &
    generateSoilMoisture &
    , getNameSoilMoisture &
    , initSoilMoisture &
    , shouldRunSoilMoisture

  use Mod_SoilMoistureClima, only : &
    generateSoilMoistureClima &
    , getNameSoilMoistureClima &
    , initSoilMoistureClima &
    , shouldRunSoilMoistureClima

  use Mod_SoilMoistureWeekly, only : &
    generateSoilMoistureWeekly &
    , getNameSoilMoistureWeekly &
    , initSoilMoistureWeekly &
    , shouldRunSoilMoistureWeekly

  use Mod_SSTClima, only : &
    generateSSTClima &
    , getNameSSTClima &
    , initSSTClima &
    , shouldRunSSTClima

  use Mod_SSTDailyDirec, only : &
    generateSSTDailyDirec &
    , getNameSSTDailyDirec &
    , initSSTDailyDirec &
    , shouldRunSSTDailyDirec

  use Mod_SSTMonthlyDirec, only : &
    generateSSTMonthlyDirec &
    , getNameSSTMonthlyDirec &
    , initSSTMonthlyDirec &
    , shouldRunSSTMonthlyDirec

  use Mod_SSTSeasonDirec, only : &
    generateSSTSeasonDirec &
    , getNameSSTSeasonDirec &
    , initSSTSeasonDirec &
    , shouldRunSSTSeasonDirec

  use Mod_SSTWeekly, only : &
    generateSSTWeekly &
    , getNameSSTWeekly &
    , initSSTWeekly &
    , shouldRunSSTWeekly

  use Mod_SSTWeeklyNCEP, only : &
    generateSSTWeeklyNCEP &
    , getNameSSTWeeklyNCEP &
    , initSSTWeeklyNCEP &
    , shouldRunSSTWeeklyNCEP

  use Mod_Temperature, only : &
    generateTemperature &
    , getNameTemperature &
    , initTemperature &
    , shouldRunTemperature
  
  use Mod_TemperatureClima, only : &
    generateTemperatureClima &
    , getNameTemperatureClima &
    , initTemperatureClima &
    , shouldRunTemperatureClima
  
  use Mod_TopographyGradient, only : &
    generateTopographyGradient &
    , getNameTopographyGradient &
    , initTopographyGradient &
    , shouldRunTopographyGradient

  use Mod_TopoSpectral, only : &
    generateTopoSpectral &
    , getNameTopoSpectral &
    , initTopoSpectral &
    , shouldRunTopoSpectral

  use Mod_TopoWaterPercGT30, only : &
    generateTopoWaterPercGT30 &
    , getNameTopoWaterPercGT30 &
    , initTopoWaterPercGT30 &
    , shouldRunTopoWaterPercGT30

  use Mod_TopoWaterPercNavy, only : &
    generateTopoWaterPercNavy &
    , getNameTopoWaterPercNavy &
    , initTopoWaterPercNavy &
    , shouldRunTopoWaterPercNavy

  use Mod_VarTopo, only : &
    generateVarTopo &
    , getNameVarTopo &
    , initVarTopo &
    , shouldRunVarTopo

  use Mod_VegetationAlbedoSSiB, only : &
    generateVegetationAlbedoSSiB &
    , getNameVegetationAlbedoSSiB &
    , initVegetationAlbedoSSiB &
    , shouldRunVegetationAlbedoSSiB

  use Mod_VegetationMask, only : &
    generateVegetationMask &
    , getNameVegetationMask &
    , initVegetationMask &
    , shouldRunVegetationMask

  use Mod_VegetationMaskIBIS, only : &
    generateVegetationMaskIBIS &
    , getNameVegetationMaskIBIS &
    , initVegetationMaskIBIS &
    , shouldRunVegetationMaskIBIS

  use Mod_VegetationMaskIBISClima, only : &
    generateVegetationMaskIBISClima &
    , getNameVegetationMaskIBISClima &
    , initVegetationMaskIBISClima &
    , shouldRunVegetationMaskIBISClima

  use Mod_VegetationMaskSSiB, only : &
    generateVegetationMaskSSiB &
    , getNameVegetationMaskSSiB &
    , initVegetationMaskSSiB &
    , shouldRunVegetationMaskSSiB

  use Mod_Namelist

  use Mod_String_Functions

  use Mod_Parallelism, only : &
    myId &
    , parf_bcast_int_scalar &
    , parf_bcast_int_1d &
    , parf_bcast_char &
    , parf_barrier &
    , isMasterProc &
    , mpiMasterProc &
    , getMyIdString &
    , parf_send_char &
    , parf_get_char &
    , fatalError &
    , parf_get_noblock_char &
    , parf_send_noblock_char &
    , parf_wait_any_nostatus &
    , maxNodes

  use Mod_Parallelism_Group, only : &
    createMpiCommGroup

  use Mod_Messages, only : &
    msgOut &
    , setDebugMode &
    , msgInLineFormatOut &
    , msgNewLine &
    , msgWarningOut

  use Mod_Messages_Parallel, only : &
    msgOutMaster &
    , msgWarningOutMaster

  implicit none

  include 'pre.h'
  include 'messages.h'
  include 'files.h'
  include 'precision.h'

  public :: runPre

  private
  character(len = *), parameter :: headerMsg = 'Master Pre Program    | '

  ! RunNameList - namelist which contains the list of modules that will run
  ! max modules to run, each with max characters in namelist name, defined in pre.h
  integer :: numberOfModulesToRun ! size of modulesToRun
  character (len = maxModuleNameLength), dimension(maxNumberOfModules) :: modulesToRun = ""
  namelist /RunNameList/ modulesToRun

  ! ModulesCommonNameList - values common to all modules
  type(varCommonNameListData)  ::  varCommon
  namelist /ModulesCommonNameList/ varCommon

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
  type(ModuleConfig) :: modCfg(maxNumberOfModules)
  namelist /ModulesNameList/ modCfg

  character (len = *), parameter :: fixedStringProc = "__PROC_#__"
  character (len = *), parameter :: fixedStringProcGroup = "__STAND_ALONE_TOTAL_"
  character (len = *), parameter :: fixedStringMasterProc = "__TOTAL__"


contains


  subroutine runPre()
    !# Main subroutine to run Pre
    !# ---
    !# @info
    !# **Brief:** This subroutine is the main method for use in this module.
    !# The namelist PRE_run.nml is read to configure:
    !# - $RunNameList - array of strings that contais the the modules names
    !#   which will run
    !# - $ModulesCommonNameList - common parameters to all modules
    !# - The next Namelists are specific for each module
    !# Then the modules configured to run, runs in parrallel, each one in one
    !# processor (until the module be converted to its parallel version), and
    !# the modules which runs in parallel uses more than one processors (like
    !# Chopping). </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 </br>
    !# @endinfo

    implicit none
    integer :: fu_namelist, fu_config
    character(len = maxPathLength) :: nameListFileName, configFileName
    integer :: moduleIndex  ! temporary variables for iteration of modules
    integer :: runIdx ! temporary variable for iteration of run modules
    integer :: procRank        ! temporary variable for mpi rank
    integer :: dependenciesShouldRunCount
    integer :: processorsNum_aux
    character (len = maxModuleNameLength) :: moduleName_aux, moduleNameWithProcs

    ! variables used in execution of master processor
    character(len = maxModuleNameLength), allocatable :: modulesToExecuteRemaining(:)
    character(len = maxModuleNameLength), allocatable :: modulesExecutedOnTurn(:)
    character(len = maxModuleNameLength), allocatable :: modulesProcessesToExecuteRemaining(:)
    character(len = maxModuleNameLength), allocatable :: modulesExecuted(:)

    integer, allocatable :: requestRecv(:), requestSend(:)
    integer :: tagRecv, tagSend, numberOfModulesToRunAtOnce, &
      dummyRecNum, runIndex, totalProcessors, numRemainingMods
    integer, parameter :: tagSendInitialValue = 10000
    logical :: isStandAloneRun
    logical, parameter :: isDebugMessageType = .true.
    integer :: ioError   
    !# Variable to handle error on manipulating files

    call setDebugMode(.false.)

    call msgOutMaster(headerMsg, "Initializing PRE ...")
    mpiMasterProc = maxNodes - 1

    call msgOutMaster("$" // headerMsg, "Reading Modules Configuration File ...")
    call createConfigFile()
    configFileName = getConfigFileName()
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
    call createNamelistFile()
    nameListFileName = getNamelistFileName()
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

    ! TODO - load namelist only of modules will run and dependencies (problem 
    ! with calling dependencyShouldRun without reading namelist initialize
    ! variables and check processors avaiable
    do moduleIndex = 1, numberOfModules
      moduleName_aux = modCfg(moduleIndex)%moduleName
      processorsNum_aux = modCfg(moduleIndex)%processors
      ! Checks
      if(modCfg(moduleIndex)%runsAlone) then
        if(processorsNum_aux > maxNodes) then
          call fatalError("$" // headerMsg, "Module " // trim(moduleName_aux) // " needs at least " // intToStr(processorsNum_aux) // " processors to run at once !")
        endif
      else
        if(processorsNum_aux > maxNodes - 1) then
          call fatalError("$" // headerMsg, "Module " // trim(moduleName_aux) // " needs at least " // intToStr(processorsNum_aux) // " processors (+1 master) to run at once !")
        endif
      endif
      call initModule(trim(moduleName_aux), fu_namelist)
    enddo

    ! Do checkings ...
    if(isMasterProc()) then
      call msgOutMaster("$" // headerMsg, "The following modules in RunNameList will run:")
      do runIdx = 1, numberOfModulesToRun
        moduleName_aux = modulesToRun(runIdx)
        moduleIndex = getModuleIndex(trim(moduleName_aux))
        call msgOutMaster("", " => " // trim(moduleName_aux) // ", " // trim(intToStr(modCfg(moduleIndex)%processors)) // " processors")
      enddo

      dependenciesShouldRunCount = 0
      ! Iterates over modules to run, check dependecies should run and adds to the run list
      call msgOutMaster("$" // headerMsg, "The following modules will run even if not in RunNameList, because modules in RunNameList depends on its output files:")
      do runIdx = 1, numberOfModulesToRun
        call addDependenciesInModulesToRun(trim(modulesToRun(runIdx)), dependenciesShouldRunCount)
      enddo

      if(dependenciesShouldRunCount == 0) then
        call msgOutMaster("", " => None modules will run");
      endif

      if(maxNodes < 2) then
        call fatalError("$" // headerMsg, "Pre program needs at least two processors to run, one for master and one for processing modules ... ")
      endif

      ! check if there are avaiable processors to run at once each module
      totalProcessors = 0
      do runIdx = 1, numberOfModulesToRun
        moduleIndex = getModuleIndex(modulesToRun(runIdx))
        totalProcessors = totalProcessors + modCfg(moduleIndex)%processors
      enddo

      if(numberOfModulesToRun .eq. 0) then
        call fatalError("$" // headerMsg, "No modules defined in RunNameList ! Exiting ...")
      endif
    endif

    call parf_barrier(1)
    ! Begin execution ...
    call createParallelismGroupsForModulesToRun()
    if(isMasterProc()) then
      allocate(modulesExecuted(0))

      ! while there are modules to run
      if(maxNodes > totalProcessors + 1) then
        call msgOutMaster("$" // headerMsg, "Warning:The number of processors in run is greather than the number of modules processors sum.")
        call msgOutMaster(headerMsg, "You are wasting some processors ...")
        call msgOutMaster(headerMsg, "Total processors required (Maximum, if all modules runs at same time): " // trim(intToStr(totalProcessors + 1)))
        call msgOutMaster(headerMsg, "Total processors used: " // trim(intToStr(maxNodes)))
      endif

      allocate(modulesToExecuteRemaining(0))
      numRemainingMods = getNextModulesToRun(modulesExecuted, modulesToExecuteRemaining)
      if(numRemainingMods == 0) then
        call fatalError("$" // headerMsg, "No process can execute due to dependencies not executed.&Check dependencies needed by modules. Maybe there are cyclic dependencies ")
      endif

      call msgOutMaster("$" // headerMsg, "================================ Execution of modules Start  =============================")
      call msgOutMaster(headerMsg, "==========================================================================================")
      call msgOutMaster(headerMsg, "Total Processors in run: " // trim(intToStr(maxNodes)))


      ! allocates modules to run in processors avaiable. Run first modules without dependencies
      do while(numRemainingMods > 0)

        call getNextModulesProcessToRun(modulesToExecuteRemaining, modulesProcessesToExecuteRemaining)
        isStandAloneRun = size(modulesToExecuteRemaining) == 1 .and. modCfg(getModuleIndex(modulesToExecuteRemaining(1)))%runsAlone
        numberOfModulesToRunAtOnce = min(size(modulesProcessesToExecuteRemaining), maxNodes - 1) ! -1 to remove master

        call msgOutMaster("$" // headerMsg, "============================ Processes of modules will run now ============================", isDebugMessageType)
        do runIndex = 1, numberOfModulesToRunAtOnce
          call msgOutMaster("", " => " // trim(modulesProcessesToExecuteRemaining(runIndex)), isDebugMessageType)
        enddo

        allocate(modulesExecutedOnTurn(numberOfModulesToRunAtOnce))
        allocate(requestRecv(numberOfModulesToRunAtOnce))
        allocate(requestSend(numberOfModulesToRunAtOnce))

        ! Receives finished messages from slaves
        do runIndex = 1, numberOfModulesToRunAtOnce
          procRank = runIndex - 1
          tagRecv = procRank
          call parf_get_noblock_char(modulesExecutedOnTurn(runIndex), maxModuleNameLength, procRank, tagRecv, requestRecv(runIndex))
        enddo

        ! Send signal to slaves start running
        do runIndex = 1, numberOfModulesToRunAtOnce
          moduleName_aux = modulesProcessesToExecuteRemaining(runIndex)
          procRank = runIndex - 1
          tagSend = procRank + tagSendInitialValue
          call parf_send_noblock_char(moduleName_aux, maxModuleNameLength, procRank, tagSend, requestSend(runIndex))
          call msgOutMaster(headerMsg, "Master Processor send to processor #" // trim(intToStr(procRank)) // " module: " // trim(moduleName_aux) // ". Total sends = " // trim(intToStr(runIndex)), isDebugMessageType)
        enddo

        ! Runs the module even in the master
        if(isStandAloneRun .and. size(modulesProcessesToExecuteRemaining) == maxNodes) then
          call msgOutMaster("", " => " // trim(modulesToExecuteRemaining(1)) // fixedStringProc // trim(intToStr(maxNodes-1)) // fixedStringProcGroup // trim(intToStr(maxNodes)))
          call runModule(modulesToExecuteRemaining(1))
        endif

        ! Waits receives of finished messages from slaves
        do runIndex = 1, numberOfModulesToRunAtOnce
          call parf_wait_any_nostatus(numberOfModulesToRunAtOnce, requestRecv, dummyRecNum)
          call msgOutMaster(headerMsg, "Master Processor Received finished signal from " // trim(intToStr(runIndex)) // " processor of " // trim(intToStr(numberOfModulesToRunAtOnce)) // " processors", isDebugMessageType)
        enddo

        call markModulesAsRan(modulesExecutedOnTurn, modulesExecuted)
        ! TODO check modulesExecutedOnTurn == modulesToExecuteRemaining before mark modules (maybe stop/warning)

        ! Waits send signal to slaves start running
        do runIndex = 1, numberOfModulesToRunAtOnce
          call msgOutMaster(headerMsg, "Waiting for " // trim(intToStr(numberOfModulesToRunAtOnce - runIndex + 1)) // " processors of " // trim(intToStr(numberOfModulesToRunAtOnce)) // " to finish", isDebugMessageType)
          call parf_wait_any_nostatus(numberOfModulesToRunAtOnce, requestSend, dummyRecNum)
          call msgOutMaster(headerMsg, "Finished " // trim(intToStr(numberOfModulesToRunAtOnce - runIndex + 1)) // " processors of " // trim(intToStr(numberOfModulesToRunAtOnce)), isDebugMessageType)
        enddo

        deallocate(modulesToExecuteRemaining)
        deallocate(modulesExecutedOnTurn)
        deallocate(modulesProcessesToExecuteRemaining)
        deallocate(requestRecv)
        deallocate(requestSend)

        allocate(modulesToExecuteRemaining(0))
        numRemainingMods = getNextModulesToRun(modulesExecuted, modulesToExecuteRemaining)

      enddo

      call msgOutMaster("$" // headerMsg, "Sending end of processing signal to all processors ...", isDebugMessageType)
      ! Send signal "*" to kill slaves
      allocate(requestSend(maxNodes - 1))
      do runIndex = 1, maxNodes - 1
        procRank = runIndex - 1
        tagSend = procRank + tagSendInitialValue
        call flush(p_nfprt)
        call parf_send_noblock_char("*", maxModuleNameLength, procRank, tagSend, requestSend(runIndex))
      enddo

      ! ... processors receive modules and execute them
    else
      ! while there are modules to run
      do while (.true.)
        moduleNameWithProcs = ""
        call parf_get_char(moduleNameWithProcs, maxModuleNameLength, mpiMasterProc, myId + tagSendInitialValue)
        if(moduleNameWithProcs(1:1) .eq. "*") exit
        moduleName_aux = moduleNameWithProcs(1:index(trim(moduleNameWithProcs), fixedStringProc) - 1)
        call runModule(moduleName_aux)
        call parf_send_char(moduleNameWithProcs, maxModuleNameLength, mpiMasterProc, myId)
      enddo
    endif

    call flush(p_nfprt)
    call msgOutMaster(headerMsg, "All processors are free")

    call parf_barrier(1)
    call msgOutMaster(headerMsg, "Closing Namelist File ...")
    call closeNamelistFile(fu_namelist)
    call destroyNamelistFile()
    call msgOutMaster(headerMsg, "Closing Config File ...")
    call closeNamelistFile(fu_config)
    call destroyConfigFile()
    call msgOutMaster(headerMsg, "PRE execution ends")
  end subroutine runPre


  function isValidModuleName(moduleName) result(isValid)
    !# Checks if exists module, given name
    !# ---
    !# @info
    !# **Brief:** Checks if exists module, given name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 </br>
    !# @endinfo
    implicit none
    character(len = *), intent(in) :: moduleName
    logical :: isValid

    isValid = getModuleIndex(moduleName) > 0
  end function isValidModuleName


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


  function dependencyShouldRun(moduleName) result(shouldRun)
    !# Checks if dependency modules should run
    !# ---
    !# @info
    !# **Brief:** Checks if dependency modules should run, even if dependency was
    !# not in run list. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 </br>
    !# @endinfo
    character (len = *), intent(in) :: moduleName
    logical :: shouldRun

    shouldRun = .false.
    ! TODO Implement should run modules
    if (trim(moduleName) .eq. getNameAlbedo()) then
      shouldRun = shouldRunAlbedo()
    elseif(trim(moduleName) .eq. getNameAlbedoClima()) then
      shouldRun = shouldRunAlbedoClima()
    elseif(trim(moduleName) .eq. getNameChopping()) then
      shouldRun = shouldRunChopping()
    elseif(trim(moduleName) .eq. getNameCO2MonthlyDirec()) then
      shouldRun = shouldRunCO2MonthlyDirec()
    elseif (trim(moduleName) .eq. getNameDeepSoilTemperature()) then
      shouldRun = shouldRunDeepSoilTemperature()
    elseif (trim(moduleName) .eq. getNameDeepSoilTemperatureClima()) then
      shouldRun = shouldRunDeepSoilTemperatureClima()
    elseif (trim(moduleName) .eq. getNameDeltaTempColdestClima()) then
      shouldRun = shouldRunDeltaTempColdestClima()
    elseif (trim(moduleName) .eq. getNameDeltaTempColdest()) then
      shouldRun = shouldRunDeltaTempColdest()
    elseif(trim(moduleName) .eq. getNameFLUXCO2Clima()) then
      shouldRun = shouldRunFLUXCO2Clima()
    elseif(trim(moduleName) .eq. getNameLandSeaMask()) then
      shouldRun = shouldRunLandSeaMask()
    elseif(trim(moduleName) .eq. getNameOCMClima()) then
      shouldRun = shouldRunOCMClima()
    elseif(trim(moduleName) .eq. getNamePorceClayMaskIBIS()) then
      shouldRun = shouldRunPorceClayMaskIBIS()
    elseif(trim(moduleName) .eq. getNamePorceClayMaskIBISClima()) then
      shouldRun = shouldRunPorceClayMaskIBISClima()
    elseif(trim(moduleName) .eq. getNamePorceSandMaskIBIS()) then
      shouldRun = shouldRunPorceSandMaskIBIS()
    elseif(trim(moduleName) .eq. getNamePorceSandMaskIBISClima()) then
      shouldRun = shouldRunPorceSandMaskIBISClima()
    elseif (trim(moduleName) .eq. getNameRoughnessLength()) then
      shouldRun = shouldRunRoughnessLength()
    elseif (trim(moduleName) .eq. getNameRoughnessLengthClima()) then
      shouldRun = shouldRunRoughnessLengthClima()
    elseif(trim(moduleName) .eq. getNameSnowClima()) then
      shouldRun = shouldRunSnowClima()
    elseif(trim(moduleName) .eq. getNameSnowWeeklyNCEP()) then
      shouldRun = shouldRunSnowWeeklyNCEP()
    elseif (trim(moduleName) .eq. getNameSoilMoisture()) then
      shouldRun = shouldRunSoilMoisture()
    elseif (trim(moduleName) .eq. getNameSoilMoistureClima()) then
      shouldRun = shouldRunSoilMoistureClima()
    elseif (trim(moduleName) .eq. getNameSoilMoistureWeekly()) then
      shouldRun = shouldRunSoilMoistureWeekly()
    elseif(trim(moduleName) .eq. getNameSSTClima()) then
      shouldRun = shouldRunSSTClima()
    elseif(trim(moduleName) .eq. getNameSSTDailyDirec()) then
      shouldRun = shouldRunSSTDailyDirec()
    elseif(trim(moduleName) .eq. getNameSSTMonthlyDirec()) then
      shouldRun = shouldRunSSTMonthlyDirec()
    elseif(trim(moduleName) .eq. getNameSSTSeasonDirec()) then
      shouldRun = shouldRunSSTSeasonDirec()
    elseif(trim(moduleName) .eq. getNameSSTWeekly()) then
      shouldRun = shouldRunSSTWeekly()
    elseif(trim(moduleName) .eq. getNameSSTWeeklyNCEP()) then
      shouldRun = shouldRunSSTWeeklyNCEP()
    elseif (trim(moduleName) .eq. getNameTemperature()) then
      shouldRun = shouldRunTemperature()
    elseif (trim(moduleName) .eq. getNameTemperatureClima()) then
      shouldRun = shouldRunTemperatureClima()
    elseif(trim(moduleName) .eq. getNameTopographyGradient()) then
      shouldRun = shouldRunTopographyGradient()
    elseif(trim(moduleName) .eq. getNameTopoSpectral()) then
      shouldRun = shouldRunTopoSpectral()
    elseif(trim(moduleName) .eq. getNameTopoWaterPercGT30()) then
      shouldRun = shouldRunTopoWaterPercGT30()
    elseif(trim(moduleName) .eq. getNameTopoWaterPercNavy()) then
      shouldRun = shouldRunTopoWaterPercNavy()
    elseif(trim(moduleName) .eq. getNameVarTopo()) then
      shouldRun = shouldRunVarTopo()
    elseif(trim(moduleName) .eq. getNameVegetationAlbedoSSiB()) then
      shouldRun = shouldRunVegetationAlbedoSSiB()
    elseif(trim(moduleName) .eq. getNameVegetationMask()) then
      shouldRun = shouldRunVegetationMask()
    elseif(trim(moduleName) .eq. getNameVegetationMaskIBIS()) then
      shouldRun = shouldRunVegetationMaskIBIS()
    elseif(trim(moduleName) .eq. getNameVegetationMaskIBISClima()) then
      shouldRun = shouldRunVegetationMaskIBISClima()
    elseif(trim(moduleName) .eq. getNameVegetationMaskSSiB()) then
      shouldRun = shouldRunVegetationMaskSSiB()
    endif

  end function dependencyShouldRun


  subroutine initModule(moduleName, fu_namelist)
    !# Reads namelist, given module name
    !# ---
    !# @info
    !# **Brief:** Reads namelist, given module name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 </br>
    !# @endinfo
    character (len = *), intent(in) :: moduleName
    integer, intent(in) :: fu_namelist


    if (trim(moduleName) .eq. getNameAlbedo()) then
      call initAlbedo(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameAlbedoClima()) then
      call initAlbedoClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameChopping()) then
      call initChopping(fu_namelist)
    elseif(trim(moduleName) .eq. getNameCO2MonthlyDirec()) then
      call initCO2MonthlyDirec(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameDeepSoilTemperature()) then
      call initDeepSoilTemperature(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameDeepSoilTemperatureClima()) then
      call initDeepSoilTemperatureClima(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameDeltaTempColdestClima()) then
      call initDeltaTempColdestClima(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameDeltaTempColdest()) then
      call initDeltaTempColdest(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameFLUXCO2Clima()) then
      call initFLUXCO2Clima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameLandSeaMask()) then
      call initLandSeaMask(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameOCMClima()) then
      call initOCMClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNamePorceClayMaskIBIS()) then
      call initPorceClayMaskIBIS(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNamePorceClayMaskIBISClima()) then
      call initPorceClayMaskIBISClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNamePorceSandMaskIBIS()) then
      call initPorceSandMaskIBIS(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNamePorceSandMaskIBISClima()) then
      call initPorceSandMaskIBISClima(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameRoughnessLength()) then
      call initRoughnessLength(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameRoughnessLengthClima()) then
      call initRoughnessLengthClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSnowClima()) then
      call initSnowClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSnowWeeklyNCEP()) then
      call initSnowWeeklyNCEP(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameSoilMoisture()) then
      call initSoilMoisture(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameSoilMoistureClima()) then
      call initSoilMoistureClima(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameSoilMoistureWeekly()) then
      call initSoilMoistureWeekly(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSSTClima()) then
      call initSSTClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSSTDailyDirec()) then
      call initSSTDailyDirec(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSSTMonthlyDirec()) then
      call initSSTMonthlyDirec(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSSTSeasonDirec()) then
      call initSSTSeasonDirec(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSSTWeekly()) then
      call initSSTWeekly(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameSSTWeeklyNCEP()) then
      call initSSTWeeklyNCEP(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameTemperature()) then
      call initTemperature(fu_namelist, varCommon)
    elseif (trim(moduleName) .eq. getNameTemperatureClima()) then
      call initTemperatureClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameTopographyGradient()) then
      call initTopographyGradient(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameTopoSpectral()) then
      call initTopoSpectral(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameTopoWaterPercGT30()) then
      call initTopoWaterPercGT30(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameTopoWaterPercNavy()) then
      call initTopoWaterPercNavy(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameVarTopo()) then
      call initVarTopo(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameVegetationAlbedoSSiB()) then
      call initVegetationAlbedoSSiB(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameVegetationMask()) then
      call initVegetationMask(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameVegetationMaskIBIS()) then
      call initVegetationMaskIBIS(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameVegetationMaskIBISClima()) then
      call initVegetationMaskIBISClima(fu_namelist, varCommon)
    elseif(trim(moduleName) .eq. getNameVegetationMaskSSiB()) then
      call initVegetationMaskSSiB(fu_namelist, varCommon)
    else
      call FatalError("$" // headerMsg, "Wrong Module Name for reading namelist: " // trim(moduleName))
    endif
  end subroutine initModule


  subroutine runModule(moduleName)
    !# Runs module, given module name
    !# ---
    !# @info
    !# **Brief:** Runs module, given module name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 </br>
    !# @endinfo
    implicit none

    character (len = *), intent(in) :: moduleName
    logical :: isParallelModule
    real (kind = p_r8) :: t1, t2, tdiff
    logical :: isModuleExecOk
    integer :: ierror

    isModuleExecOk = .false.
    call msgOut(headerMsg, " ~ ~ ~ ~ ~ Starting execution of module " // trim(moduleName) // " at processor " // trim(getMyIdString()))

    isParallelModule = modCfg(getModuleIndex(moduleName))%processors > 1
    if(isParallelModule) then
      if (isMasterProc()) call cpu_time(t1)
    else
      call cpu_time(t1)
    endif
    if (trim(moduleName) .eq. getNameAlbedo()) then
      isModuleExecOk = generateAlbedo()
    elseif (trim(moduleName) .eq. getNameAlbedoClima()) then
      isModuleExecOk = generateAlbedoClima()
    elseif(trim(moduleName) .eq. getNameChopping()) then
      isModuleExecOk = generateChopping()
    elseif(trim(moduleName) .eq. getNameCO2MonthlyDirec()) then
      isModuleExecOk = generateCO2MonthlyDirec()
    elseif (trim(moduleName) .eq. getNameDeepSoilTemperature()) then
      isModuleExecOk = generateDeepSoilTemperature()
    elseif (trim(moduleName) .eq. getNameDeepSoilTemperatureClima()) then
      isModuleExecOk = generateDeepSoilTemperatureClima()
    elseif (trim(moduleName) .eq. getNameDeltaTempColdestClima()) then
      isModuleExecOk = generateDeltaTempColdestClima()
    elseif (trim(moduleName) .eq. getNameDeltaTempColdest()) then
      isModuleExecOk = generateDeltaTempColdest()
    elseif(trim(moduleName) .eq. getNameFLUXCO2Clima()) then
      isModuleExecOk = generateFLUXCO2Clima()
    elseif(trim(moduleName) .eq. getNameLandSeaMask()) then
      isModuleExecOk = generateLandSeaMask()
    elseif(trim(moduleName) .eq. getNameOCMClima()) then
      isModuleExecOk = generateOCMClima()
    elseif(trim(moduleName) .eq. getNamePorceClayMaskIBIS()) then
      isModuleExecOk = generatePorceClayMaskIBIS()
    elseif(trim(moduleName) .eq. getNamePorceClayMaskIBISClima()) then
      isModuleExecOk = generatePorceClayMaskIBISClima()
    elseif(trim(moduleName) .eq. getNamePorceSandMaskIBIS()) then
      isModuleExecOk = generatePorceSandMaskIBIS()      
    elseif(trim(moduleName) .eq. getNamePorceSandMaskIBISClima()) then
      isModuleExecOk = generatePorceSandMaskIBISClima()      
    elseif (trim(moduleName) .eq. getNameRoughnessLength()) then
      isModuleExecOk = generateRoughnessLength()
    elseif (trim(moduleName) .eq. getNameRoughnessLengthClima()) then
      isModuleExecOk = generateRoughnessLengthClima()
    elseif(trim(moduleName) .eq. getNameSnowClima()) then
      isModuleExecOk = generateSnowClima()
    elseif(trim(moduleName) .eq. getNameSnowWeeklyNCEP()) then
      isModuleExecOk = generateSnowWeeklyNCEP()
    elseif (trim(moduleName) .eq. getNameSoilMoisture()) then
      isModuleExecOk = generateSoilMoisture()
    elseif (trim(moduleName) .eq. getNameSoilMoistureClima()) then
      isModuleExecOk = generateSoilMoistureClima()
    elseif (trim(moduleName) .eq. getNameSoilMoistureWeekly()) then
      isModuleExecOk = generateSoilMoistureWeekly()
    elseif(trim(moduleName) .eq. getNameSSTClima()) then
      isModuleExecOk = generateSSTClima()
    elseif(trim(moduleName) .eq. getNameSSTDailyDirec()) then
      isModuleExecOk = generateSSTDailyDirec()
    elseif(trim(moduleName) .eq. getNameSSTMonthlyDirec()) then
      isModuleExecOk = generateSSTMonthlyDirec()
    elseif(trim(moduleName) .eq. getNameSSTSeasonDirec()) then
      isModuleExecOk = generateSSTSeasonDirec()
    elseif(trim(moduleName) .eq. getNameSSTWeekly()) then
      isModuleExecOk = generateSSTWeekly()
    elseif(trim(moduleName) .eq. getNameSSTWeeklyNCEP()) then
      isModuleExecOk = generateSSTWeeklyNCEP()
    elseif (trim(moduleName) .eq. getNameTemperature()) then
      isModuleExecOk = generateTemperature()
    elseif (trim(moduleName) .eq. getNameTemperatureClima()) then
      isModuleExecOk = generateTemperatureClima()
    elseif(trim(moduleName) .eq. getNameTopographyGradient()) then
      isModuleExecOk = generateTopographyGradient()
    elseif(trim(moduleName) .eq. getNameTopoSpectral()) then
      isModuleExecOk = generateTopoSpectral()
    elseif(trim(moduleName) .eq. getNameTopoWaterPercGT30()) then
      isModuleExecOk = generateTopoWaterPercGT30()
    elseif(trim(moduleName) .eq. getNameTopoWaterPercNavy()) then
      isModuleExecOk = generateTopoWaterPercNavy()
    elseif(trim(moduleName) .eq. getNameVarTopo()) then
      isModuleExecOk = generateVarTopo()
    elseif(trim(moduleName) .eq. getNameVegetationAlbedoSSiB()) then
      isModuleExecOk = generateVegetationAlbedoSSiB()
    elseif(trim(moduleName) .eq. getNameVegetationMask()) then
      isModuleExecOk = generateVegetationMask()
    elseif(trim(moduleName) .eq. getNameVegetationMaskIBIS()) then
      isModuleExecOk = generateVegetationMaskIBIS()
    elseif(trim(moduleName) .eq. getNameVegetationMaskIBISClima()) then
      isModuleExecOk = generateVegetationMaskIBISClima()
    elseif(trim(moduleName) .eq. getNameVegetationMaskSSiB()) then
      isModuleExecOk = generateVegetationMaskSSiB()
    else
      call msgWarningOutMaster(headerMsg, "Wrong Module Name for running: " // trim(moduleName) // ". Program wil continue ")
    endif
    if(isParallelModule .and. .not. isMasterProc()) return
    call cpu_time(t2)
    
    if(isModuleExecOk) then
      tdiff = t2-t1
      call msgInLineFormatOut(' ' // headerMsg // "CPU time (seconds) " // trim(moduleName), '(A)')
      call msgInLineFormatOut(tdiff, '(F10.2)')
      call msgNewLine()
      call msgOut(headerMsg, " ~ ~ ~ ~ ~ Execution of module " // trim(moduleName) // " was finished successfully !!! ")
    else
      call msgWarningOut(headerMsg, " Execution of module " // trim(moduleName) // " finished anormally! Program will continue. Check results of dependent modules")
    endif

  end subroutine runModule


  subroutine createParallelismGroupsForModulesToRun()
    !# Creates Parallelism Group if necessary
    !# ---
    !# @info
    !# **Brief:** Creates Parallelism Group if necessary. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none

    integer :: processors, modIdx
    character (len = maxModuleNameLength) :: moduleName_aux

    do modIdx = 1, maxNumberOfModules
      moduleName_aux = trim(modCfg(modIdx)%moduleName)
      if(modIdx < 0 ) return
      if(modCfg(modIdx)%runsAlone) then
        processors = modCfg(modIdx)%processors
        ! defines mpi communicator for using more than one processor
        call createMpiCommGroup(trim(moduleName_aux), 0, processors)
      endif
    enddo

  end subroutine createParallelismGroupsForModulesToRun


  subroutine getModuleDependencies(moduleName, dependencies)
    !# Gets a list of module dependencies, given module name
    !# ---
    !# @info
    !# **Brief:** Gets a list of module dependencies, given module name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 </br>
    !# @endinfo
    implicit none
    character (len = *), intent(in) :: moduleName
    character(len = maxModuleNameLength), allocatable, intent(out) :: dependencies(:)
    character(len = maxModuleNameLength) :: allDependencies(maxNumberOfModules - 1)

    integer :: moduleIndex, depIdx, dependenciesSize

    moduleIndex = getModuleIndex(moduleName)
    allDependencies = modCfg(moduleIndex)%dependencies

    dependenciesSize = 0
    do depIdx = 1, maxNumberOfModules - 1
      if(isValidModuleName(allDependencies(depIdx))) then
        dependenciesSize = dependenciesSize + 1
      else
        exit
      endif
    end do

    allocate(dependencies(dependenciesSize))
    do depIdx = 1, dependenciesSize
      dependencies(depIdx) = allDependencies(depIdx)
    end do

  end subroutine getModuleDependencies


  function getNextModulesToRun(modulesExecuted, modulesToExecuteRemaining) result(modulesNameSize)
    !# Gets remaining modules to run, which had dependencies ran before
    !# ---
    !# @info
    !# **Brief:** Gets remaining modules to run, which had dependencies ran before. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxModuleNameLength), allocatable, intent(in) :: modulesExecuted(:)
    character(len = maxModuleNameLength), allocatable, intent(inout) :: modulesToExecuteRemaining(:)
    integer :: modulesNameSize

    character(len = maxModuleNameLength), allocatable :: modulesNameLast(:)
    logical :: isCandidate
    integer :: modIndex, execIndex, depIndex, execDepIndex, numDependeciesExecuted
    character(len = maxModuleNameLength), allocatable :: dependencies(:)

    modulesNameSize = 0
    allocate(modulesNameLast(0))

    do modIndex = 1, numberOfModulesToRun
      isCandidate = .true.

      do execIndex = 1, size(modulesExecuted)
        ! check if module was executed
        if(trim(modulesToRun(modIndex)) .eq. trim(modulesExecuted(execIndex))) then
          isCandidate = .false.
          exit
        endif
      enddo
      if(isCandidate) then
        ! check if module will run
        do execIndex = 1, size(modulesToExecuteRemaining)
          if(trim(modulesToRun(modIndex)) .eq. trim(modulesToExecuteRemaining(execIndex))) then
            isCandidate = .false.
            exit
          endif
        enddo

        if(isCandidate) then
          numDependeciesExecuted = 0
          call getModuleDependencies(modulesToRun(modIndex), dependencies)

          do depIndex = 1, size(dependencies)
            if(.not. dependencyShouldRun(trim(dependencies(depIndex)))) then
              numDependeciesExecuted = numDependeciesExecuted + 1
            else
              do execDepIndex = 1, size(modulesExecuted)
                if(trim(dependencies(depIndex)) .eq. trim(modulesExecuted(execDepIndex))) then
                  numDependeciesExecuted = numDependeciesExecuted + 1
                endif
              enddo
            endif
          enddo
          if(numDependeciesExecuted .eq. size(dependencies)) then
            ! all dependencies executed !
            if(modCfg(getModuleIndex(modulesToRun(modIndex)))%runsAlone) then
              deallocate(modulesToExecuteRemaining)
              allocate(modulesToExecuteRemaining(1))
              modulesToExecuteRemaining(1) = trim(modulesToRun(modIndex))
              modulesNameSize = modulesNameSize + 1
              return
            endif
            deallocate(modulesNameLast)
            allocate(modulesNameLast(modulesNameSize))
            modulesNameLast = modulesToExecuteRemaining
            deallocate(modulesToExecuteRemaining)
            allocate(modulesToExecuteRemaining(modulesNameSize + 1))
            modulesToExecuteRemaining(1:modulesNameSize) = modulesNameLast
            modulesToExecuteRemaining(modulesNameSize + 1) = trim(modulesToRun(modIndex))
            modulesNameSize = modulesNameSize + 1
          endif
        endif
      endif
    enddo

  end function getNextModulesToRun


  recursive subroutine addDependenciesInModulesToRun(moduleName_aux, dependenciesShouldRunCount)
    !# Adds dependencies in modules to run
    !# ---
    !# @info
    !# **Brief:** Adds dependencies in modules to run. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(inout) :: dependenciesShouldRunCount
    character (len = *), intent(in) :: moduleName_aux

    character(len = maxModuleNameLength), allocatable :: dependsOnModules_aux(:)
    logical :: moduleAlreadyAddedToTun
    integer :: dependentIndex, existsIdx

    ! Iterates over modules which moduleName_aux depends and add modules should run
    call getModuleDependencies(moduleName_aux, dependsOnModules_aux)
    do dependentIndex = 1, size(dependsOnModules_aux)
      moduleAlreadyAddedToTun = .false.
      if(isValidModuleName(dependsOnModules_aux(dependentIndex)) .and. dependencyShouldRun(dependsOnModules_aux(dependentIndex))) then
        do existsIdx = 1, numberOfModulesToRun
          if(trim(modulesToRun(existsIdx)) .eq. dependsOnModules_aux(dependentIndex)) then
            moduleAlreadyAddedToTun = .true.
            exit
          endif
        enddo
        if(.not. moduleAlreadyAddedToTun) then
          call addDependenciesInModulesToRun(dependsOnModules_aux(dependentIndex), dependenciesShouldRunCount)
          numberOfModulesToRun = numberOfModulesToRun + 1
          modulesToRun(numberOfModulesToRun) = trim(dependsOnModules_aux(dependentIndex))
          dependenciesShouldRunCount = dependenciesShouldRunCount + 1
          call msgOutMaster("", " => " // trim(modulesToRun(numberOfModulesToRun)))
        endif
      endif
    enddo

  end subroutine addDependenciesInModulesToRun


  subroutine getNextModulesProcessToRun(modsToExecOnTurn, modulesProcessToExecuteOnTurn)
    !# Gets next modules process to run
    !# ---
    !# @info
    !# **Brief:** Gets next modules process to run. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxModuleNameLength), allocatable, intent(in) :: modsToExecOnTurn(:)
    character(len = maxModuleNameLength), allocatable, intent(inout) :: modulesProcessToExecuteOnTurn(:)
    character(len = maxModuleNameLength), allocatable :: modulesToExecuteOnTurnLast(:)
    character(len = maxModuleNameLength), allocatable :: modulesToExecuteOnTurn(:)

    character(len = maxModuleNameLength) :: modName
    integer :: procRank, modIdx, procIdx, slotsRemaining, modProcsNum, masterProcNum, modIdx_aux

    if(size(modsToExecOnTurn) == 1 .and. modCfg(getModuleIndex(modsToExecOnTurn(1)))%runsAlone) then
      modIdx_aux = getModuleIndex(modsToExecOnTurn(1))
      modProcsNum = modCfg(modIdx_aux)%processors
      modName = modCfg(modIdx_aux)%moduleName
      allocate(modulesProcessToExecuteOnTurn(modProcsNum))
      procRank = 0
      do procIdx = 1, modProcsNum
        modulesProcessToExecuteOnTurn(procRank + 1) = trim(modName) // fixedStringProc // trim(intToStr(procRank)) // fixedStringProcGroup // trim(intToStr(modProcsNum))
        procRank = procRank + 1
      enddo
    else
      slotsRemaining = maxNodes - 1
      allocate(modulesToExecuteOnTurn(0))
      allocate(modulesToExecuteOnTurnLast(0))
      ! allocate modules to execute on turn anc update slots remaining
      do modIdx = 1, size(modsToExecOnTurn)
        modName = modsToExecOnTurn(modIdx)
        modIdx_aux = getModuleIndex(modName)
        if(modCfg(modIdx_aux)%processors <= slotsRemaining) then
          deallocate(modulesToExecuteOnTurnLast)
          allocate(modulesToExecuteOnTurnLast(size(modulesToExecuteOnTurn)))
          modulesToExecuteOnTurnLast = modulesToExecuteOnTurn
          deallocate(modulesToExecuteOnTurn)
          allocate(modulesToExecuteOnTurn(size(modulesToExecuteOnTurnLast) + 1))
          modulesToExecuteOnTurn(1:size(modulesToExecuteOnTurnLast)) = modulesToExecuteOnTurnLast
          modulesToExecuteOnTurn(size(modulesToExecuteOnTurnLast) + 1) = modName
          slotsRemaining = slotsRemaining - modCfg(modIdx_aux)%processors
        endif
      enddo
      allocate(modulesProcessToExecuteOnTurn(maxNodes - 1 - slotsRemaining))

      procRank = 0
      do modIdx = 1, size(modulesToExecuteOnTurn)
        modName = modulesToExecuteOnTurn(modIdx)
        modProcsNum = modCfg(getModuleIndex(modName))%processors
        masterProcNum = procRank
        do procIdx = 1, modProcsNum
          modulesProcessToExecuteOnTurn(procRank + 1) = trim(modName) // fixedStringProc // trim(intToStr(procRank)) // fixedStringMasterProc // trim(intToStr(masterProcNum))
          procRank = procRank + 1
        end do
      end do
    endif

  end subroutine getNextModulesProcessToRun


  subroutine markModulesAsRan(modulesExecutedOnTurn, modulesExecuted)
    !# Marks modules as Ran
    !# ---
    !# @info
    !# **Brief:** Marks modules as Ran </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxModuleNameLength), allocatable, intent(in) :: modulesExecutedOnTurn(:)
    character(len = maxModuleNameLength), allocatable, intent(inout) :: modulesExecuted(:)
    character(len = maxModuleNameLength), allocatable :: modulesExecutedLast(:)
    character(len = maxModuleNameLength) :: moduleName_aux
    integer :: modIdx, modExecIdx, newIdx
    logical :: isExecuted

    newIdx = 0
    do modIdx = 1, size(modulesExecutedOnTurn)
      moduleName_aux = modulesExecutedOnTurn(modIdx)
      moduleName_aux = moduleName_aux(1:index(trim(moduleName_aux), fixedStringProc) - 1)

      ! check already executed
      isExecuted = .false.
      do modExecIdx = 1, size(modulesExecuted)
        if(trim(moduleName_aux) .eq. trim(modulesExecuted(modExecIdx))) then
          isExecuted = .true.
          exit
        endif
      enddo

      if(.not. isExecuted) then
        allocate(modulesExecutedLast(size(modulesExecuted)))
        modulesExecutedLast = modulesExecuted
        deallocate(modulesExecuted)
        allocate(modulesExecuted(size(modulesExecutedLast) + 1))
        modulesExecuted(1:size(modulesExecutedLast)) = modulesExecutedLast
        modulesExecuted(size(modulesExecutedLast) + 1) = trim(moduleName_aux)
        deallocate(modulesExecutedLast)
      endif
    enddo

  end subroutine markModulesAsRan


end module Mod_Pre

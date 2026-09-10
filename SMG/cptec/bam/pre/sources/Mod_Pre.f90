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

   use Mod_Namelist, only : &
         initializeNamelistsDefaultFileNames &
         , finalizeNamelists &
         , getVarCommon &
         , getmodcfg &
         , getNumberOfModules &
         , numberOfModulesToRun &
         , getModuleToRun &
         , getModuleIndex &
         , ModuleConfig

   use Mod_String_Functions, only : &
         intToStr

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

   use Mod_ModulesToRun, only : &
         markModulesAsRan &
         , getNextModulesToRun &
         , addAllDependenciesInModulesToRun &
         , getNextModulesProcessToRun &
         , initializeModules &
         , getTotalProcessors

   use Mod_AbstractModulesFacade, only : &
         AbstractModulesFacade

   implicit none

   include 'pre.h'
   include 'messages.h'
   include 'files.h'
   include 'precision.h'

   public :: runPre, runPreInitializedModulesAndNamelists

   private
   character(len = *), parameter :: headerMsg = 'Master Pre Program    | '


contains


   subroutine runPre(moduleFacade)
      implicit none
      class(AbstractModulesFacade), intent(inout) :: moduleFacade

            !# Variable to handle error on manipulating files
      call setDebugMode(.true.)
      call msgOutMaster(headerMsg, "Initializing PRE ...")
      call initializeNamelistsDefaultFileNames()
      call initializeModules(moduleFacade)
      call runPreInitializedModulesAndNamelists(moduleFacade, 2)
   end subroutine runPre


   subroutine runPreInitializedModulesAndNamelists(moduleFacade, execRoundSleepTime)
   
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
      class(AbstractModulesFacade), intent(inout) :: moduleFacade
      integer, intent(in) :: execRoundSleepTime

      integer :: procRank        ! temporary variable for mpi rank
      character (len = maxModuleNameLength) :: moduleName_aux, moduleNameWithProcs

      ! variables used in execution of master processor
      character(len = maxModuleNameLength), allocatable :: modulesToExecuteRemaining(:)
      character(len = maxModuleNameLength), allocatable :: modulesExecutedOnTurn(:)
      character(len = maxModuleNameLength), allocatable :: modulesProcessesToExecuteRemaining(:)
      character(len = maxModuleNameLength), allocatable :: modulesExecuted(:)

      integer, allocatable :: requestRecv(:)
      integer :: tagRecv, tagSend, numberOfModulesToRunAtOnce, &
            dummyRecNum, runIndex, numRemainingMods, execRound
      integer, parameter :: tagSendInitialValue = 1000
      logical :: isStandAloneRun
      logical, parameter :: isDebugMessageType = .true.

      type(ModuleConfig) :: modCfg


      call parf_barrier(1)
      call msgOutMaster(headerMsg, "Running PRE ...")
      ! Begin execution ...
      call createParallelismGroupsForModulesToRun()
      if(isMasterProc()) then
         allocate(modulesExecuted(0))

         ! while there are modules to run
         if(maxNodes > getTotalProcessors() + 1) then
            call msgOutMaster("$" // headerMsg, "Warning:The number of processors in run is greather than the number of modules processors sum.")
            call msgOutMaster(headerMsg, "You are wasting some processors ...")
            call msgOutMaster(headerMsg, "Total processors required (Maximum, if all modules runs at same time): " // trim(intToStr(getTotalProcessors() + 1)))
            call msgOutMaster(headerMsg, "Total processors used: " // trim(intToStr(maxNodes)))
         endif

         allocate(modulesToExecuteRemaining(0))
         numRemainingMods = getNextModulesToRun(moduleFacade, modulesExecuted, modulesToExecuteRemaining)
         if(numRemainingMods == 0) then
            call fatalError("$" // headerMsg, "No process can execute due to dependencies not executed.&Check dependencies needed by modules. Maybe there are cyclic dependencies ")
         endif

         call msgOutMaster("$" // headerMsg, "================================ Execution of modules Start  =============================")
         call msgOutMaster(headerMsg, "==========================================================================================")
         call msgOutMaster(headerMsg, "Total Processors in run: " // trim(intToStr(maxNodes)))


         ! TODO refactor this code to use the same code for Mod_Pre and tests
         ! allocates modules to run in processors avaiable. Run first modules without dependencies
         execRound = 0
         do while(numRemainingMods > 0)

            call getNextModulesProcessToRun(modulesToExecuteRemaining, modulesProcessesToExecuteRemaining)
            modCfg = getModCfg(getModuleIndex(modulesToExecuteRemaining(1)))
            isStandAloneRun = size(modulesToExecuteRemaining) == 1 .and. modCfg%runsAlone
            numberOfModulesToRunAtOnce = min(size(modulesProcessesToExecuteRemaining), maxNodes - 1) ! -1 to remove master
            
            execRound = execRound + 1
            call msgOutMaster("$" // headerMsg, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Execution Round " // trim(intToStr(execRound)) // " will execute: ")
            do runIndex = 1, numberOfModulesToRunAtOnce
               call msgOutMaster("", " => " // trim(modulesProcessesToExecuteRemaining(runIndex)))
            enddo

            allocate(modulesExecutedOnTurn(numberOfModulesToRunAtOnce))
            allocate(requestRecv(numberOfModulesToRunAtOnce))

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
               call parf_send_char(moduleName_aux, maxModuleNameLength, procRank, tagSend)
               call msgOutMaster(headerMsg, "Master Processor send to processor #" // trim(intToStr(procRank)) // " tag: " //trim(intToStr(tagSend)) // " module: " & 
                              // trim(moduleName_aux) // ". Total sends = " // trim(intToStr(runIndex)), isDebugMessageType)
            enddo

            ! if runsAlone (chopping does ) then runs the module even in the master
            if(isStandAloneRun .and. size(modulesProcessesToExecuteRemaining) == maxNodes) then
               call msgOutMaster("", " => " // trim(modulesToExecuteRemaining(1)) // fixedStringProc // trim(intToStr(maxNodes - 1)) // fixedStringProcGroup // trim(intToStr(maxNodes)), isDebugMessageType)
               call moduleFacade%runModule(modulesToExecuteRemaining(1))
            endif

            ! Waits receives of finished messages from slaves
            do runIndex = 1, numberOfModulesToRunAtOnce
               call parf_wait_any_nostatus(numberOfModulesToRunAtOnce, requestRecv, dummyRecNum)
               call msgOutMaster(headerMsg, "Master Processor Received finished signal from " // trim(intToStr(runIndex)) // " processors of " // trim(intToStr(numberOfModulesToRunAtOnce)) & 
                                 // " processors", isDebugMessageType)
            enddo

            call markModulesAsRan(modulesExecutedOnTurn, modulesExecuted)
            ! TODO check modulesExecutedOnTurn == modulesToExecuteRemaining before mark modules (maybe stop/warning)

            deallocate(modulesToExecuteRemaining)
            deallocate(modulesExecutedOnTurn)
            deallocate(modulesProcessesToExecuteRemaining)
            deallocate(requestRecv)

            call msgOutMaster(headerMsg, "Waiting " // trim(intToStr(execRoundSleepTime)) // " seconds for round ends, due to IO time. " )
            
            call flush(p_nfprt)
            call sleep(execRoundSleepTime)
            allocate(modulesToExecuteRemaining(0))
            numRemainingMods = getNextModulesToRun(moduleFacade, modulesExecuted, modulesToExecuteRemaining)

         enddo

         call msgOutMaster("$" // headerMsg, "Sending end of processing signal to all processors ...", isDebugMessageType)
         ! Send signal "*" to kill slaves
         do runIndex = 1, maxNodes - 1
            procRank = runIndex - 1
            tagSend = procRank + tagSendInitialValue
            call flush(p_nfprt)
            call parf_send_char("*", maxModuleNameLength, procRank, tagSend)
         enddo

         ! ... processors receive modules and execute them
      else
         ! while there are modules to run
         do while (.true.)
            moduleNameWithProcs = ""
            call parf_get_char(moduleNameWithProcs, maxModuleNameLength, mpiMasterProc, myId + tagSendInitialValue)
            if(moduleNameWithProcs(1:1) .eq. "*") exit
            call msgOut(headerMsg, "Process " // trim(moduleNameWithProcs) // " of rank " // trim(intToStr(myId)) // " tag " // trim(intToStr(myId + tagSendInitialValue)) &
                        // " will execute " , isDebugMessageType)
            moduleName_aux = moduleNameWithProcs(1:index(trim(moduleNameWithProcs), fixedStringProc) - 1)
            call moduleFacade%runModule(moduleName_aux)
            call parf_send_char(moduleNameWithProcs, maxModuleNameLength, mpiMasterProc, myId)
         enddo
      endif

      call flush(p_nfprt)
      call msgOutMaster(headerMsg, "All processors are free")

      call parf_barrier(1)
      call finalizeNamelists()
      call msgOutMaster(headerMsg, "PRE execution ends")
   end subroutine runPreInitializedModulesAndNamelists


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

      integer :: modIdx
      character (len = maxModuleNameLength) :: moduleName_aux
      type(ModuleConfig) :: modCfg

      do modIdx = 1, maxNumberOfModules
         modCfg = getModCfg(modIdx)
         moduleName_aux = trim(modCfg%moduleName)
         if(modIdx < 0) return
         if(modCfg%runsAlone) then
            ! defines mpi communicator for using more than one processor
            call createMpiCommGroup(trim(moduleName_aux), 0, modCfg%processors)
         endif
      enddo

   end subroutine createParallelismGroupsForModulesToRun


end module Mod_Pre


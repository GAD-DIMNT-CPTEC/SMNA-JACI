!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_ModulesToRun </br></br>
!#
!# **Brief**: Selects modules to run</br></br>
!# 
!# **Author**: Denis Eiras </br>
!#
!# **Version**: 1.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>06-08-2021 - Denis Eiras    - version: 1.0.0</li></br>
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

module Mod_ModulesToRun

  use Mod_String_Functions, only : intToStr
  use Mod_Parallelism, only : maxNodes, isMasterProc, fatalError
  use Mod_Messages_Parallel, only : msgOutMaster
  use Mod_AbstractModulesFacade, only: AbstractModulesFacade
  use Mod_Namelist, only: getModuleIndex, getModCfg, ModuleConfig, modulesToRun &
    , numberOfModulesToRun, getModuleToRun, getNumberOfModules

  implicit none
  include 'pre.h'
!  include 'messages.h'
!  include 'files.h'
!  include 'precision.h'

  public :: markModulesAsRan, getNextModulesToRun, addAllDependenciesInModulesToRun &
    , getNextModulesProcessToRun, getmoduledependencies, initializeModules, getTotalProcessors

  private
  character(len = *), parameter :: headerMsg = 'ModulesToRun Module  | '
  integer :: totalProcessors = 0


contains


   function getTotalProcessors()
     implicit none
     integer :: getTotalProcessors

     getTotalProcessors = totalProcessors
   end function getTotalProcessors


   subroutine initializeModules(moduleFacade)
      implicit none

      class(AbstractModulesFacade), intent(inout) :: moduleFacade  ! Facade for Concrete or Mock implementation
      type(ModuleConfig) :: modCfg  ! auxiliar variable of ModuleConfig type
      integer :: processorsNum_aux  ! temporary variable
      character (len = maxModuleNameLength) :: moduleName_aux ! moduleName_aux
      integer :: moduleIndex  ! temporary variables for iteration of modules
      integer :: runIdx ! temporary variable for iteration of run modules
      integer :: dependenciesShouldRunCount  ! count number of dependencies


      ! TODO - load namelist only of modules will run and dependencies (problem
      ! with calling dependencyShouldRun without reading namelist initialize
      ! variables and check processors avaiable
      do moduleIndex = 1, getNumberOfModules()
         modCfg = getModCfg(moduleIndex)
         moduleName_aux = modCfg%moduleName
         processorsNum_aux = modCfg%processors
         ! Checks
         if(modCfg%runsAlone) then
            if(processorsNum_aux > maxNodes) then
               call fatalError("$" // headerMsg, "Module " // trim(moduleName_aux) // " needs at least " // intToStr(processorsNum_aux) // " processors to run at once !")
            endif
         else
            if(processorsNum_aux > maxNodes - 1) then
               call fatalError("$" // headerMsg, "Module " // trim(moduleName_aux) // " needs at least " // intToStr(processorsNum_aux) // " processors (+1 master) to run at once !")
            endif
         endif
         call moduleFacade%initModule(trim(moduleName_aux))
      enddo

      if(isMasterProc()) then
         call msgOutMaster("$" // headerMsg, "The following modules in RunNameList are set to run:")
         do runIdx = 1, numberOfModulesToRun
            moduleName_aux = getModuleToRun(runIdx)
            moduleIndex = getModuleIndex(trim(moduleName_aux))
            modCfg = getModCfg(moduleIndex)
            call msgOutMaster("", " => " // trim(moduleName_aux) // ", " // trim(intToStr(modCfg%processors)) // " processors")
         enddo

         call addAllDependenciesInModulesToRun(moduleFacade, dependenciesShouldRunCount)

         if(maxNodes < 2) then
            call fatalError("$" // headerMsg, "Pre program needs at least two processors to run, one for master and one for processing modules ... ")
         endif

         ! check if there are avaiable processors to run at once each module
         totalProcessors = 0
         do runIdx = 1, numberOfModulesToRun
            moduleIndex = getModuleIndex(getModuleToRun(runIdx))
            modCfg = getModCfg(moduleIndex)
            totalProcessors = totalProcessors + modCfg%processors
         enddo

         if(numberOfModulesToRun .eq. 0) then
            call fatalError("$" // headerMsg, "No modules defined in RunNameList ! Exiting ...")
         endif
      endif
   end subroutine initializeModules


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

    integer :: depIdx, dependenciesSize
    type(ModuleConfig) :: modCfg

    modCfg=getModCfg(getModuleIndex(moduleName))
    allDependencies = modCfg%dependencies

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


  function getNextModulesToRun(moduleFacade, modulesExecuted, modulesToExecuteRemaining) result(modulesNameSize)
    !# Gets remaining modules to run, which had dependencies ran before
    !# ---
    !# @info
    !# **Brief:** Gets remaining modules to run, which had dependencies ran before. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    class(AbstractModulesFacade), intent(inout) :: moduleFacade

    character(len = maxModuleNameLength), allocatable, intent(in) :: modulesExecuted(:)
    character(len = maxModuleNameLength), allocatable, intent(inout) :: modulesToExecuteRemaining(:)
    integer :: modulesNameSize

    character(len = maxModuleNameLength), allocatable :: modulesNameLast(:)
    logical :: isCandidate, depWillRunInRound, depWillRun
    integer :: modIndex, execIndex, depIndex, execDepIndex, numDependeciesExecuted
    character(len = maxModuleNameLength), allocatable :: dependencies(:)
    type(ModuleConfig) :: modCfg

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
            ! check if dependency will run ( to not run with dependent)
            depWillRunInRound = .false.
            do execIndex = 1, size(modulesToExecuteRemaining)
              if(trim(dependencies(depIndex)) .eq. trim(modulesToExecuteRemaining(execIndex))) then
                depWillRunInRound = .true.
                exit
              endif
            enddo
            if(.not. depWillRunInRound) then
              ! check if dependency is in list to run ... (A)
              depWillRun = .false.
              do execIndex = 1, numberOfModulesToRun
                if(trim(dependencies(depIndex)) .eq. trim(modulesToRun(execIndex))) then
                  depWillRun = .true.
                endif
              enddo
              if (depWillRun) then
                ! (A)... and if dependency already ran, sum as executed 
                do execDepIndex = 1, size(modulesExecuted)
                  if(trim(dependencies(depIndex)) .eq. trim(modulesExecuted(execDepIndex))) then
                    numDependeciesExecuted = numDependeciesExecuted + 1
                  endif
                enddo
              else
                ! if dependency should not run, sum as executed 
                if(.not. moduleFacade%dependencyShouldRun(trim(dependencies(depIndex)))) then
                  numDependeciesExecuted = numDependeciesExecuted + 1
                endif

              endif
            endif

          enddo
          if(numDependeciesExecuted .eq. size(dependencies)) then
            ! all dependencies executed !
            modCfg=getModCfg(getModuleIndex(modulesToRun(modIndex)))
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


  subroutine addAllDependenciesInModulesToRun(moduleFacade, dependenciesShouldRunCount)
    implicit none
    class(AbstractModulesFacade), intent(inout) :: moduleFacade
    integer :: runIdx, dependenciesShouldRunCount
    character(len = maxModuleNameLength), allocatable :: dependencies(:)

    dependenciesShouldRunCount = 0
    ! Iterates over modules to run, check dependecies should run and adds to the run list
    call msgOutMaster("$" // headerMsg, "The following modules will run even if not in RunNameList, because modules in RunNameList depends on its output files:")



    do runIdx = 1, numberOfModulesToRun
      call addDependenciesInModulesToRun(moduleFacade, trim(getModuleToRun(runIdx)), dependenciesShouldRunCount)
    enddo

    ! print*, "numberOfModulesToRunFinal=", numberOfModulesToRun
    ! print*, "modulesToRunFinal", modulesToRun(1:numberOfModulesToRun)

    if(dependenciesShouldRunCount == 0) then
      call msgOutMaster("", " => None modules will run");
    endif

  end subroutine addAllDependenciesInModulesToRun


  recursive subroutine addDependenciesInModulesToRun(moduleFacade, moduleName_aux, dependenciesShouldRunCount)
    !# Adds dependencies in modules to run
    !# ---
    !# @info
    !# **Brief:** Adds dependencies in modules to run. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    class(AbstractModulesFacade), intent(inout) :: moduleFacade

    character (len = *), intent(in) :: moduleName_aux
    integer, intent(inout) :: dependenciesShouldRunCount

    character(len = maxModuleNameLength), allocatable :: dependsOnModules_aux(:)
    logical :: moduleAlreadyAddedToRun
    integer :: dependentIndex, existsIdx

    ! Iterates over modules which moduleName_aux depends and add modules should run
    call getModuleDependencies(moduleName_aux, dependsOnModules_aux)
    do dependentIndex = 1, size(dependsOnModules_aux)
      moduleAlreadyAddedToRun = .false.
      if(isValidModuleName(dependsOnModules_aux(dependentIndex))) then
        do existsIdx = 1, numberOfModulesToRun
          if(trim(modulesToRun(existsIdx)) .eq. dependsOnModules_aux(dependentIndex)) then
            moduleAlreadyAddedToRun = .true.
            exit
          endif
        enddo
        if(.not. moduleAlreadyAddedToRun .and. moduleFacade%dependencyShouldRun(trim(dependsOnModules_aux(dependentIndex))) ) then
          numberOfModulesToRun = numberOfModulesToRun + 1
          modulesToRun(numberOfModulesToRun) = trim(dependsOnModules_aux(dependentIndex))
          dependenciesShouldRunCount = dependenciesShouldRunCount + 1
          call msgOutMaster("", " => " // trim(modulesToRun(numberOfModulesToRun)))
        endif
        call addDependenciesInModulesToRun(moduleFacade, dependsOnModules_aux(dependentIndex), dependenciesShouldRunCount)
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
    type(ModuleConfig) :: modCfg

    modCfg=getModCfg(getModuleIndex(modsToExecOnTurn(1)))
    if(size(modsToExecOnTurn) == 1 .and. modCfg%runsAlone) then
      modIdx_aux = getModuleIndex(modsToExecOnTurn(1))
      modCfg=getModCfg(modIdx_aux)
      modProcsNum = modCfg%processors
      modName = modCfg%moduleName
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
        modCfg=getModCfg(modIdx_aux)
        if(modCfg%processors <= slotsRemaining) then
          deallocate(modulesToExecuteOnTurnLast)
          allocate(modulesToExecuteOnTurnLast(size(modulesToExecuteOnTurn)))
          modulesToExecuteOnTurnLast = modulesToExecuteOnTurn
          deallocate(modulesToExecuteOnTurn)
          allocate(modulesToExecuteOnTurn(size(modulesToExecuteOnTurnLast) + 1))
          modulesToExecuteOnTurn(1:size(modulesToExecuteOnTurnLast)) = modulesToExecuteOnTurnLast
          modulesToExecuteOnTurn(size(modulesToExecuteOnTurnLast) + 1) = modName
          slotsRemaining = slotsRemaining - modCfg%processors
        endif
      enddo
      allocate(modulesProcessToExecuteOnTurn(maxNodes - 1 - slotsRemaining))

      procRank = 0
      do modIdx = 1, size(modulesToExecuteOnTurn)
        modName = modulesToExecuteOnTurn(modIdx)
        modCfg=getModCfg(getModuleIndex(modName))
        modProcsNum = modCfg%processors
        masterProcNum = procRank
        do procIdx = 1, modProcsNum
          modulesProcessToExecuteOnTurn(procRank + 1) = trim(modName) // fixedStringProc // trim(intToStr(procRank)) // fixedStringMasterProc // trim(intToStr(masterProcNum))
          procRank = procRank + 1
        end do
      end do
    endif

  end subroutine getNextModulesProcessToRun


  subroutine markModulesAsRan(modulesProcsExecutedOnTurn, modulesExecuted)
    !# Marks modules as Ran
    !# ---
    !# @info
    !# **Brief:** Marks modules as Ran </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxModuleNameLength), allocatable, intent(in) :: modulesProcsExecutedOnTurn(:)
    character(len = maxModuleNameLength), allocatable, intent(inout) :: modulesExecuted(:)
    character(len = maxModuleNameLength), allocatable :: modulesExecutedLast(:)
    character(len = maxModuleNameLength) :: moduleName_aux
    integer :: modIdx, modExecIdx, newIdx
    logical :: isExecuted

    newIdx = 0
    do modIdx = 1, size(modulesProcsExecutedOnTurn)
      moduleName_aux = modulesProcsExecutedOnTurn(modIdx)
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


end module Mod_ModulesToRun

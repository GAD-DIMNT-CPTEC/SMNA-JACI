!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_ModulesFacade </br></br>
!#
!# **Brief**: Facade for acess Pre Modules</br></br>
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

module Mod_ModulesFacade

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


   ! Another modules
   use Mod_Namelist, only : &
         get_fu_namelist &
         , getVarCommon &
         , getModuleIndex &
         , getModCfg &
         , ModuleConfig &
         , getNumberOfModules

   use Mod_Parallelism, only : &
         getmyidstring &
         , isMasterProc &
         , fatalError &
         , maxNodes

   use Mod_Messages_Parallel, only : &
         msgOutMaster &
         , msgwarningoutmaster

   use Mod_Messages, only : &
         msgwarningout &
         , msgnewline &
         , msginlineformatout &
         , msgout

   use Mod_AbstractModulesFacade, only : &
         AbstractModulesFacade

   use Mod_String_Functions, only : &
         intToStr

   implicit none
   private
   include 'pre.h'
   include 'precision.h'

   type, extends(AbstractModulesFacade) :: ConcreteModulesFacade
   contains
      procedure :: runModule => concreteRunModule
      procedure :: initModule => concreteInitModule
      procedure :: dependencyShouldRun => concreteDependencyShouldRun
   end type ConcreteModulesFacade

   character(len = *), parameter :: headerMsg = 'ModulesFacade Module  | '

   public :: ConcreteModulesFacade


contains

   ! ContreteModulesFacade methods ==============================

   function concreteDependencyShouldRun(this, moduleName) result(shouldRun)
      !# Checks if dependency modules should run
      !# ---
      !# @info
      !# **Brief:** Checks if dependency modules should run, even if dependency was
      !# not in run list. </br>
      !# **Authors**: </br>
      !# &bull; Denis Eiras </br>
      !# **Date**: apr/2019 </br>
      !# @endinfo
      class(ConcreteModulesFacade) :: this
      character (len = *), intent(in) :: moduleName
      logical :: shouldRun

      shouldRun = .false.
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

   end function concreteDependencyShouldRun


   subroutine concreteInitModule(this, moduleName)
      !# Reads namelist, given module name
      !# ---
      !# @info
      !# **Brief:** Reads namelist, given module name. </br>
      !# **Authors**: </br>
      !# &bull; Denis Eiras </br>
      !# **Date**: apr/2019 </br>
      !# @endinfo
      class(ConcreteModulesFacade) :: this
      character (len = *), intent(in) :: moduleName

      if (trim(moduleName) .eq. getNameAlbedo()) then
         call initAlbedo(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameAlbedoClima()) then
         call initAlbedoClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameChopping()) then
         call initChopping(get_fu_namelist())
      elseif(trim(moduleName) .eq. getNameCO2MonthlyDirec()) then
         call initCO2MonthlyDirec(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameDeepSoilTemperature()) then
         call initDeepSoilTemperature(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameDeepSoilTemperatureClima()) then
         call initDeepSoilTemperatureClima(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameDeltaTempColdestClima()) then
         call initDeltaTempColdestClima(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameDeltaTempColdest()) then
         call initDeltaTempColdest(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameFLUXCO2Clima()) then
         call initFLUXCO2Clima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameLandSeaMask()) then
         call initLandSeaMask(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameOCMClima()) then
         call initOCMClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNamePorceClayMaskIBIS()) then
         call initPorceClayMaskIBIS(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNamePorceClayMaskIBISClima()) then
         call initPorceClayMaskIBISClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNamePorceSandMaskIBIS()) then
         call initPorceSandMaskIBIS(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNamePorceSandMaskIBISClima()) then
         call initPorceSandMaskIBISClima(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameRoughnessLength()) then
         call initRoughnessLength(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameRoughnessLengthClima()) then
         call initRoughnessLengthClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSnowClima()) then
         call initSnowClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSnowWeeklyNCEP()) then
         call initSnowWeeklyNCEP(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameSoilMoisture()) then
         call initSoilMoisture(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameSoilMoistureClima()) then
         call initSoilMoistureClima(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameSoilMoistureWeekly()) then
         call initSoilMoistureWeekly(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSSTClima()) then
         call initSSTClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSSTDailyDirec()) then
         call initSSTDailyDirec(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSSTMonthlyDirec()) then
         call initSSTMonthlyDirec(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSSTSeasonDirec()) then
         call initSSTSeasonDirec(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSSTWeekly()) then
         call initSSTWeekly(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameSSTWeeklyNCEP()) then
         call initSSTWeeklyNCEP(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameTemperature()) then
         call initTemperature(get_fu_namelist(), getVarCommon())
      elseif (trim(moduleName) .eq. getNameTemperatureClima()) then
         call initTemperatureClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameTopographyGradient()) then
         call initTopographyGradient(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameTopoSpectral()) then
         call initTopoSpectral(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameTopoWaterPercGT30()) then
         call initTopoWaterPercGT30(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameTopoWaterPercNavy()) then
         call initTopoWaterPercNavy(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameVarTopo()) then
         call initVarTopo(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameVegetationAlbedoSSiB()) then
         call initVegetationAlbedoSSiB(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameVegetationMask()) then
         call initVegetationMask(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameVegetationMaskIBIS()) then
         call initVegetationMaskIBIS(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameVegetationMaskIBISClima()) then
         call initVegetationMaskIBISClima(get_fu_namelist(), getVarCommon())
      elseif(trim(moduleName) .eq. getNameVegetationMaskSSiB()) then
         call initVegetationMaskSSiB(get_fu_namelist(), getVarCommon())
      else
         call FatalError("$" // headerMsg, "Wrong Module Name for reading namelist: " // trim(moduleName))
      endif
   end subroutine concreteInitModule


   subroutine concreteRunModule(this, moduleName)
      !# Runs module, given module name
      !# ---
      !# @info
      !# **Brief:** Runs module, given module name. </br>
      !# **Authors**: </br>
      !# &bull; Denis Eiras </br>
      !# **Date**: apr/2019 </br>
      !# @endinfo
      implicit none
      class(ConcreteModulesFacade) :: this
      character (len = *), intent(in) :: moduleName
      logical :: isParallelModule
      real (kind = p_r8) :: t1, t2, tdiff
      logical :: isModuleExecOk
      integer :: ierror
      type(ModuleConfig) :: modCfg

      isModuleExecOk = .false.
      call msgOut(headerMsg, "~ ~ ~ " // trim(moduleName) // " started at processor " // trim(getMyIdString()))

      modCfg = getModCfg(getModuleIndex(moduleName))
      isParallelModule = modCfg%processors > 1

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
         tdiff = t2 - t1
         call msgInLineFormatOut(' ' // headerMsg // "o o o " // trim(moduleName) // " sucessfullt executed in ", '(A)')
         call msgInLineFormatOut(tdiff, '(F10.2)')
         call msgInLineFormatOut(" seconds ", '(A)')
         call msgnewline()
      else
         call msgWarningOut(headerMsg, " Execution of module " // trim(moduleName) // " finished anormally! Program will continue. Check results of dependent modules")
      endif

   end subroutine concreteRunModule


end module Mod_ModulesFacade

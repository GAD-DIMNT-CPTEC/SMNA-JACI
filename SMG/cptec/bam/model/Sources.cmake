#-----------------------------------------------------------------------------!
# Project: BAM_Model
#-----------------------------------------------------------------------------!
#BOI
# !TITLE: Sources.cmake - Source Files for BAM_Model
# !AUTHORS: João Gerd Zell de Mattos
# !AFFILIATION: CPTEC/INPE - Data Assimilation Development Group (GDAD)
# !DATE: 2025-02-26
# !INTRODUCTION:
# This script defines a list variable with all the Fortran source files
# for the BAM_Model project. Use absolute or relative paths consistently.
# It is assumed that this file is located in the same directory as your
# CMakeLists.txt, and that the source code is under the "source/<etc>" directory.
#
# If this file is located at the root of the project (same directory as CMakeLists.txt),
# CMAKE_CURRENT_LIST_DIR is used to construct absolute paths in a robust way.
#EOI


set(BAM_MODEL_SOURCES
    # Files from source/Assimilation
    ${CMAKE_CURRENT_LIST_DIR}/source/Assimilation/GridDump.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Assimilation/SpecDump.f90
    
    # Files from source/Diagnostics
    ${CMAKE_CURRENT_LIST_DIR}/source/Diagnostics/Diagnostics.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Diagnostics/GridHistory.f90
    
    # Files from source/Dynamics
    ${CMAKE_CURRENT_LIST_DIR}/source/Dynamics/GridDynamics.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Dynamics/SemiLagrangian.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Dynamics/SpecDynamics.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Dynamics/TimeStep.f90
    
    # Files from source/Fields
    ${CMAKE_CURRENT_LIST_DIR}/source/Fields/FieldsDynamics.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Fields/FieldsPhysics.f90
    
    # Files from source/Initialization
    ${CMAKE_CURRENT_LIST_DIR}/source/Initialization/Init.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Initialization/NonLinearNMI.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Initialization/Options.f90
    
    # Files from source/InputOutput
    ${CMAKE_CURRENT_LIST_DIR}/source/InputOutput/Dumpgraph.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/InputOutput/InputOutput.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/InputOutput/IOLowLevel.f90
    
    # Files from source/Main
    ${CMAKE_CURRENT_LIST_DIR}/source/Main/Atmos_Model.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Main/Model.f90
    
    # Files from source/Parallelism
    ${CMAKE_CURRENT_LIST_DIR}/source/Parallelism/Communications.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Parallelism/Parallelism.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Parallelism/Sizes.f90

    # Files from source/Physics/PhysicsDriver
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/PhysicsDriver.f90

    # Files from source/Physics/BoundaryLayer/PblDriver
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/BoundaryLayer/PblDriver.f90

    # Files from source/Physics/BoundaryLayer
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/BoundaryLayer/HostlagBoville/Pbl_HostlagBoville.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/BoundaryLayer/MellorYamada0/Pbl_MellorYamada0.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/BoundaryLayer/MellorYamada1/Pbl_MellorYamada1.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/BoundaryLayer/ParkBretherton/Pbl_UniversityWashington.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/BoundaryLayer/PBL_Entrain.f90
    
    # Files from source/Physics/Convection
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/Convection.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/DeepConvection.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/GDM_BAM/Cu_GDM_BAM.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/GrellEnsCPTEC/Cu_Grellens_CPTEC.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/GrellEns/Cu_Grellens.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/Kuo/Cu_Kuolcl.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/RAS3PHASE/Cu_RAS3PHASE.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/Ras/Cu_RAS.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/DeepConvection/Zhang/Cu_ZhangMcFarlane.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/Ferrier/Micro_Ferrier.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/Hack/Micro_Hack.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/HWRF/Micro_HWRF.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/LrgScl/Micro_LrgScl.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/MicroPhysics.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/MORRISON_AERO/Micro_HugMorr.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/MORRISON/Micro_MORR.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/UKME/Micro_UKME.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/MicroPhysics/UKME/StratCloudFraction.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/ShallowConvection/JHack/Shall_JHack.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/ShallowConvection/MasFlux/Shall_MasFlux.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/ShallowConvection/ShallowConvection.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/ShallowConvection/Souza/Shall_Souza.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/ShallowConvection/Tied/Shall_Tied.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Convection/ShallowConvection/UWShCu/Shall_UWShCu.f90

    # Files from source/Physics/GravityWaveDra/GwddDriver
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/GravityWaveDrag/GwddDriver.f90

    # Files from source/Physics/GravityWaveDra
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/GravityWaveDrag/Alpert/GwddSchemeAlpert.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/GravityWaveDrag/CAM/GwddSchemeCAM.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/GravityWaveDrag/ECMWF/Gwdd_ECMWF.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/GravityWaveDrag/UKMET/GwddSchemeCPTEC.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/GravityWaveDrag/UKMET/GwddSchemeUSSP.f90
    
    # Files from source/Physics/Radiation/RadiationDriver
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/RadiationDriver.f90

    # Files from source/Physics/Radiation
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/Clirad/Rad_Clirad.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/Clirad/Rad_CliRadLW.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/CliradTarasova/Rad_CliRadLWTarasova.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/CliradTarasova/Rad_CliradTarasova.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/CloudOpticalProperty.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/COLA/Rad_COLA.f90
#    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/RRTMG/Rad_RRTMG.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Radiation/UKMET/Rad_UKMO.f90

    # Files from source/Physics/Surface/SfcPBLDriver 
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/SfcPBLDriver.f90

    # Files from source/Physics/Surface
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Land/IBIS2.6/Sfc_Ibis_BioGeoChemistry.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Land/IBIS2.6/Sfc_Ibis_BioGeoPhysics.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Land/IBIS2.6/Sfc_Ibis_Fiels.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Land/IBIS2.6/Sfc_Ibis_Interface.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Land/SiB2.5/Sfc_SiB2.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Land/SSiB/Sfc_SSiB.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Ocean/SeaFlux_COLA/Sfc_SeaFlux_COLA_Model.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Ocean/SeaFlux_UKME/Sfc_SeaFlux_UKME_Model.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Ocean/SeaFlux_WGFS/Sfc_SeaFlux_WGFS_Model.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Ocean/SeaIceFlux_WRF/Sfc_SeaIceFlux_WRF_Model.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Ocean/Sfc_SeaFlux_Interface.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Ocean/SLAB/SlabOceanModel.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/sfcpbl/Sfc_MellorYamada0.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/sfcpbl/Sfc_MellorYamada1.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/Surface/Surface.f90
    
    # Files from source/Physics/ThermCellV0
    ${CMAKE_CURRENT_LIST_DIR}/source/Physics/ThermCellV0/ModThermalCell.f90
    
    # Files from source/Transform
    ${CMAKE_CURRENT_LIST_DIR}/source/Transform/Transform.f90
    
    # Files from source/Utils
    ${CMAKE_CURRENT_LIST_DIR}/source/Utils/Constants.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Utils/PhysicalFunctions.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Utils/Utils.f90
    ${CMAKE_CURRENT_LIST_DIR}/source/Utils/Watches.f90
)

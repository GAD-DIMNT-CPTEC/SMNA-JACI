# ----------------------------------------------------------------------
# RRTMG_Ext: CMakeLists
# This file compiles the entire RRTMG package (LW, SW, MCICA) as one library.
# ----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.16)
project(RRTMG_Ext LANGUAGES Fortran)

# --- 1) csm_share library ---
add_library(rrtmg_share
    ${CMAKE_CURRENT_SOURCE_DIR}/../csm_share/shr/abortutils.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../csm_share/shr/shr_const_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../csm_share/shr/shr_kind_mod.f90
)

# --- 2) Merge LW, SW, and MCICA into One Library ---
file(GLOB RRTMG_ALL_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/rrtmg_lw/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/rrtmg_sw/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/rrtmg_mcica/*.f90"
)

# Add the driver file `Rad_RRTMG.f90`
list(APPEND RRTMG_ALL_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/../Rad_RRTMG.f90"
)

list(REMOVE_ITEM RRTMG_ALL_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/rrtmg_lw/lixo.f90"
)

# Define the RRTMG all-in-one library
add_library(rrtmg_all_lib ${RRTMG_ALL_SRCS})

# Ensure module files are placed in a single directory
set_target_properties(rrtmg_all_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)

# Link with rrtmg_share to resolve dependencies
target_link_libraries(rrtmg_all_lib PRIVATE rrtmg_share)


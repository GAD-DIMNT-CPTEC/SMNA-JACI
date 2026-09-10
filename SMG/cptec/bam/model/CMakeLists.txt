cmake_minimum_required(VERSION 3.16)
project(BAM_Model LANGUAGES Fortran)
#-----------------------------------------------------------------------------!
# Project: BAM_Model
#-----------------------------------------------------------------------------!
#BOI
# !TITLE: CMake Configuration for BAM_Model
# !AUTHORS: João Gerd Zell de Mattos
# !AFFILIATION: CPTEC/INPE - Data Assimilation Development Group (GDAD)
# !DATE: 2025-02-26
# !INTRODUCTION:
#    This CMake configuration file is designed to build the BAM_Model
#    application using CMake. It automatically detects Fortran source files,
#    sets up compiler flags based on the detected compiler, and defines an 
#    executable target. The script is fully compatible with GNU, Intel, Cray, 
#    and PGI compilers.
#
#    Supported build modes:
#      - Debug:         Includes additional checks and debugging symbols.
#      - Release:       Optimized for maximum performance.
#      - RelWithDebInfo: Combines optimizations with debugging information.
#      - MinSizeRel:    Optimized for the smallest executable size.
#
#    How to manually select the compiler:
#      To specify which compiler to use, run CMake with the following option:
#        $ cmake -DCMAKE_Fortran_COMPILER=<compiler> ..
#      
#      Examples:
#        - GNU Fortran:  $ cmake -DCMAKE_Fortran_COMPILER=gfortran ..
#        - Intel Fortran: $ cmake -DCMAKE_Fortran_COMPILER=ifort ..
#        - Cray Fortran:  $ cmake -DCMAKE_Fortran_COMPILER=crayftn ..
#        - NVIDIA (PGI):  $ cmake -DCMAKE_Fortran_COMPILER=nvfortran ..
#      
#      Alternatively, you can set the FC environment variable before running CMake:
#        $ export FC=gfortran && cmake ..
#
#    How to choose the build mode:
#      CMake allows different compilation modes using the -DCMAKE_BUILD_TYPE option:
#
#        - Debug mode (for debugging):
#          $ cmake -DCMAKE_BUILD_TYPE=Debug ..
#
#        - Release mode (optimized for performance):
#          $ cmake -DCMAKE_BUILD_TYPE=Release ..
#
#        - RelWithDebInfo mode (optimized + debug symbols):
#          $ cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
#
#        - MinSizeRel mode (optimized for minimal executable size):
#          $ cmake -DCMAKE_BUILD_TYPE=MinSizeRel ..
#
#    Example usage:
#      $ rm -fr build && mkdir build && cd build
#      $ cmake -DCMAKE_BUILD_TYPE=Debug ..
#      $ make
#EOI

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

#BOP
# !SECTION: Source File Detection
# !DESCRIPTION:
#    This section manually includes the Fortran source files from a dedicated
#    file (Sources.cmake). Maintaining a manual list of sources is a
#    recommended practice in CMake, as it ensures that any addition, removal,
#    or renaming of files will explicitly trigger a reconfiguration.
#    This avoids the potential pitfalls of using file(GLOB_RECURSE ...),
#    which does not automatically reconfigure on file system changes.
#EOP

# Include the list of source files defined in Sources.cmake
include(${CMAKE_SOURCE_DIR}/Sources.cmake)

# Check if the list is empty (optional safety check)
if(NOT BAM_MODEL_SOURCES)
    message(FATAL_ERROR
        "No Fortran source files found (BAM_MODEL_SOURCES is empty). "
        "Check Sources.cmake.")
endif()

# Criando o executável
add_executable(ParModel_MPI ${BAM_MODEL_SOURCES})

# Add the RRTMG Ext subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/source/Physics/Radiation/RRTMG/ext)

# Now we have the targets:
#   rrtmg_share
#   rrtmg_lw_lib
#   rrtmg_sw_lib
#   rrtmg_mcica_lib
#   (possibly rrtmg_driver, if you made that library)

# Ensure bam_model finds module files from rrtmg_all_lib
target_include_directories(ParModel_MPI PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/source/Physics/Radiation/RRTMG/ext/mod
)

# Ensure bam_model links to rrtmg_all_lib
target_link_libraries(ParModel_MPI PRIVATE rrtmg_all_lib rrtmg_share)


#BOP
# !SECTION: Compiler Flags
# !DESCRIPTION:
#    This section defines compiler-specific flags for different compilers
#    (GNU, Intel, Cray, PGI). Flags are set for free-form Fortran 2008 and
#    include options for bounds checking and optimization.
#EOP
# Detect Fortran compiler and set appropriate flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    #! GNU Fortran Compiler
    #!  - Aggressive optimizations: -O3, -ffast-math, -funroll-loops
    #!  - Big-endian support
    #!  - No line length limit
    set(dialect -ffree-line-length-none)
    set(optimization -ffast-math -funroll-loops -ftree-vectorizer-verbose=2)
    set(bounds -fbounds-check -fcheck=all)
    set(debug -g -Wall -fbacktrace)
    set(byte_order -fconvert=big-endian)

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    #! Intel Fortran Compiler
    #!  - Advanced optimizations: -O3, -xHost
    #!  - Additional checks in Debug mode: -check all
    #!  - Big-endian support: -convert big_endian
    set(dialect -free)
    set(optimization -O3 -xHost)
    set(bounds -check all)
    set(debug -g -traceback)
    set(byte_order -convert big_endian)

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    #! Intel Fortran Compiler
    #!  - Advanced optimizations: -O3, -xHost
    #!  - Additional checks in Debug mode: -check all
    #!  - Big-endian support: -convert big_endian
    set(dialect -free)
    set(optimization -O3 -xHost)
    set(bounds -check all)
    set(debug -g -traceback)
    set(byte_order -convert big_endian)

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    #! Cray Fortran Compiler
    #!  - Optimizations for Cray architecture: -h vector0,cache0,fp2
    #!  - Bounds checking in Debug mode: -h bounds
    set(dialect -free)
    set(optimization -O3 -h vector0,cache0,fp2)
    set(bounds -h bounds)
    set(debug -g -O0 -h traceback)
    set(byte_order)

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    #! PGI/NVIDIA Fortran Compiler
    #!  - Aggressive optimizations: -O3, -Mbyteswapio, -fast
    #!  - Additional checks in Debug mode
    set(dialect -free)
    set(optimization -O3 -Mbyteswapio -fast)
    set(bounds -Mbounds -Mchkstk -Mchkptr -Mchkfpstk)
    set(debug -g)
    set(byte_order -Mbyteswapio)
endif()



#BOP
# !SECTION: Apply Compiler Flags
# !DESCRIPTION:
#    The compiler flags are applied to the target using `target_compile_options`.
#    This ensures that only the relevant target receives the flags instead of
#    modifying global CMake variables.
#    Create variables for each configuration so we can:
#    1) Print them easily.
#    2) Keep code clear and organized.
#EOP

set(DEBUG_OPTIONS          ${dialect} ${debug} ${byte_order} ${bounds})
set(RELEASE_OPTIONS        ${dialect} ${optimization} ${byte_order})
set(RELWITHDEBINFO_OPTIONS ${dialect} ${optimization} ${byte_order} ${debug})
set(MINSIZEREL_OPTIONS     ${dialect} ${optimization} ${byte_order})

# Here is your default fallback:
set(DEFAULT_OPTIONS        ${dialect} ${byte_order})

# Print for confirmation
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "DEBUG flags:          ${DEBUG_OPTIONS}")
message(STATUS "RELEASE flags:        ${RELEASE_OPTIONS}")
message(STATUS "RELWITHDEBINFO flags: ${RELWITHDEBINFO_OPTIONS}")
message(STATUS "MINSIZEREL flags:     ${MINSIZEREL_OPTIONS}")
message(STATUS "DEFAULT (no-config) : ${DEFAULT_OPTIONS}")

# Apply compiler flags conditionally
target_compile_options(ParModel_MPI PRIVATE
    # Only used if CMAKE_BUILD_TYPE=Debug
    $<$<STREQUAL:$<CONFIG>,Debug>:${DEBUG_OPTIONS}>

    # Only used if CMAKE_BUILD_TYPE=Release
    $<$<STREQUAL:$<CONFIG>,Release>:${RELEASE_OPTIONS}>

    # Only used if CMAKE_BUILD_TYPE=RelWithDebInfo
    $<$<STREQUAL:$<CONFIG>,RelWithDebInfo>:${RELWITHDEBINFO_OPTIONS}>

    # Only used if CMAKE_BUILD_TYPE=MinSizeRel
    $<$<STREQUAL:$<CONFIG>,MinSizeRel>:${MINSIZEREL_OPTIONS}>

    # Fallback if no build type is set (CMAKE_BUILD_TYPE is empty)
    $<$<STREQUAL:$<CONFIG>,>:${DEFAULT_OPTIONS}>
)
#BOP
# !SECTION: Installation
# !DESCRIPTION:
#    This section installs the compiled executable into the `exec` directory
#    inside `CMAKE_CURRENT_LIST_DIR`.
#    This ensures a structured deployment, keeping binaries separate from other files.
#EOP

# Set the installation directory to be inside the project folder.
# All installed files (binaries, libraries, headers) will be placed inside
# the 'exec' directory next to this CMakeLists.txt file.
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}" CACHE PATH "Installation directory" FORCE)

# Install the executable in the 'bin' directory within 'exec'
install(TARGETS ParModel_MPI DESTINATION exec)


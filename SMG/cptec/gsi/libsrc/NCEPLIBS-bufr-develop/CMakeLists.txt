# ==============================================================================
# File: CMakeLists.txt
# Project: NCEPLIBS-bufr
#
# Description:
#   Main configuration file for NCEPLIBS-bufr.
#   This file defines the project, user options, build type, compiler flags,
#   and includes all necessary subdirectories (src, utils, python, tables, test, docs).
#
# Authors:
#   Mark Potts, Kyle Gerheiser, Rahul Mahajan, Ed Hartnett, Jeff Ator
#
# Revision History:
#   2025-02-24: Added bufr_types.F90 to handle Fortran kinds. Inserted the "-cpp" 
#               flag to ensure proper processing of preprocessor directives and 
#               added the "-DBUFR_INT8" directive for enabling 8-byte integer precision.
#   2025-02-24: Introduced the "KINDS" cache variable to allow users to specify,
#               via the command line, the list of kinds to compile (default remains
#               "4" if not specified).

# ==============================================================================

cmake_minimum_required(VERSION 3.15)

#==============================================================================
# Read Version Information
#------------------------------------------------------------------------------
file(STRINGS "VERSION" pVersion)

#==============================================================================
# Define the Project
#------------------------------------------------------------------------------
project(
    bufr
    VERSION ${pVersion}
    LANGUAGES C Fortran)

#==============================================================================
# User Options
#------------------------------------------------------------------------------
# ENABLE_DOCS:
#   Enable generation of doxygen-based documentation.
#   Default: OFF
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)

# ENABLE_PYTHON:
#   Enable building the Python module 'ncepbufr'.
#   Default: OFF
option(ENABLE_PYTHON "Enable building python module 'ncepbufr'" OFF)

# BUILD_SHARED_LIBS:
#   Enable building shared libraries.
#   Default: OFF
option(BUILD_SHARED_LIBS "Enable building shared libraries" OFF)

# BUILD_TESTS:
#   Enable building tests for the project.
#   Default: ON
option(BUILD_TESTS "Build tests" ON)

# TEST_FILE_DIR:
#   Specifies a local directory containing test files.
#   If present, files will be copied from here instead of fetched via FTP.
#   Default: "."
set(TEST_FILE_DIR "." CACHE STRING "Check this directory for test files before using FTP.")
message(STATUS "Finding test data files in directory ${TEST_FILE_DIR}.")

#==============================================================================
# BUFR Master Table Directory Configuration
#------------------------------------------------------------------------------
# MASTER_TABLE_DIR:
#   Defines the installation location for the Master BUFR Tables.
set(MASTER_TABLE_DIR "${CMAKE_INSTALL_PREFIX}/tables" CACHE STRING "Installation location of Master BUFR Tables")

# Format the directory string for Fortran90 and C source files, ensuring that
# each line does not exceed the maximum number of characters (_lenslmax).
string(LENGTH "${MASTER_TABLE_DIR}" _lenmtd)
set(_lenslmax 60)  # Maximum number of characters per line
if(${_lenmtd} LESS_EQUAL ${_lenslmax})
  set(MASTER_TABLE_DIR_F "${MASTER_TABLE_DIR}")
  set(MASTER_TABLE_DIR_C "${MASTER_TABLE_DIR}")
else()
  # Split the directory string into chunks of _lenslmax characters for inclusion in source code
  set(_mtd_f "")
  set(_mtd_c "")
  set(_len2go ${_lenmtd})
  set(_start 0)
  while(${_len2go} GREATER ${_lenslmax})
    string(SUBSTRING "${MASTER_TABLE_DIR}" ${_start} ${_lenslmax} _chunk)
    string(APPEND _mtd_f "${_chunk}" "' // &\n'")
    string(APPEND _mtd_c "${_chunk}" "\" \n\t\t\"")
    math(EXPR _len2go "${_len2go} - ${_lenslmax}")
    math(EXPR _start "${_start} + ${_lenslmax}")
  endwhile()
  string(SUBSTRING "${MASTER_TABLE_DIR}" ${_start} -1 _chunk)
  string(APPEND _mtd_f "${_chunk}")
  string(APPEND _mtd_c "${_chunk}")
  set(MASTER_TABLE_DIR_F "${_mtd_f}")
  set(MASTER_TABLE_DIR_C "${_mtd_c}")
endif()

#==============================================================================
# Build Type Configuration
#------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#==============================================================================
# Common Inclusions and Global Settings
#------------------------------------------------------------------------------
include(GNUInstallDirs)
enable_testing()

# Ensure position independent code is generated (required for shared libraries)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#==============================================================================
# Python Support Configuration
#------------------------------------------------------------------------------
if(ENABLE_PYTHON)
  if(APPLE)
    set(Python3_FIND_FRAMEWORK "LAST")
  endif()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

#==============================================================================
# macOS Specific Settings
#------------------------------------------------------------------------------
if(APPLE AND NOT "${CMAKE_RANLIB}" MATCHES "^.*(llvm-ranlib)$")
  # macOS linker does not include 'common symbols' by default.
  # Passing the -c flag forces their inclusion, resolving undefined symbol errors.
  set(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
  set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
endif()


#==============================================================================
# Define user-specified kinds (for example, "4:8:d")
# If not specified via command line, it defaults to "4".
set(KINDS "4" CACHE STRING "Semicolon-separated list of kinds to compile(e.g.,'4:8:d')")
string(REPLACE ":" " " KINDS_LIST "${KINDS}")
separate_arguments(KINDS_LIST)
set(kinds ${KINDS_LIST})
#--------------------------------------------------------------------------

#==============================================================================
# Library Version Configuration
#------------------------------------------------------------------------------
# Use the 'kinds' variable defined above.
# list(APPEND kinds "4") is no longer needed.

#==============================================================================
# Compiler Flags Configuration
#------------------------------------------------------------------------------
# C Compiler Flags
if(CMAKE_C_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -traceback")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -DUNDERSCORE -O3")
  set(CMAKE_C_FLAGS_DEBUG "-ggdb -Wall -O0")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS}")
endif()

# Fortran Compiler Flags
# Note: Even when building only the _4 version, _8 and _d flags are required
#       for compiling test codes.
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(fortran_8_flags "-r8 -i8 -DBUFR_INT8")
  set(fortran_d_flags "-r8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(fortran_8_flags "-fdefault-real-8 -fdefault-integer-8 -DBUFR_INT8")
  set(fortran_d_flags "-fdefault-real-8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -O3 -g -fbacktrace -fno-second-underscore -fno-range-check -fconvert=big-endian")
  set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -Wall -O0")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}")
endif()

#==============================================================================
# Include Subdirectories
#------------------------------------------------------------------------------
# Main library source code
add_subdirectory(src)

# Utility tools and additional code
add_subdirectory(utils)

# Python module (optional)
if(ENABLE_PYTHON)
  add_subdirectory(python)
endif()

# BUFR Tables required for operation
add_subdirectory(tables)

# Test suite (optional)
if(BUILD_TESTS)
  add_subdirectory(test)
endif()

#==============================================================================
# Documentation Generation
#------------------------------------------------------------------------------
# If ENABLE_DOCS is turned on, configure documentation generation using Doxygen.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  add_subdirectory(docs)
endif()


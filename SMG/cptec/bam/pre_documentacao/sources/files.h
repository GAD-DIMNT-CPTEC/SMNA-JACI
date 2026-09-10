! **Include File:** files.h
!
! **Brief:** Defines constants for Module Mod_FileManager
!
! **Author:** Denis Eiras; denis.eiras@inpe.br
!
! **Version:** 1.0.0
!
! **Changes:**
! - 01-05-2019 - Denis Eiras             - version: 1.0
!
! **Bugs:**
! - No items at this time
!
! **ToDo:**
! - No items at this time
!
! **Documentation:**
! For theoretical information, please visit the following link:
! http://urlib.net/8JMKD3MGP3W34R/3SME6J2
! Copyright Under GLP-3.0
! &copy; https://opensource.org/licenses/GPL-3.0


integer, parameter :: maxPathLength = 4096
!# Max size of path in Linux
integer, parameter :: firstUnit = 20
!# lowest io unit number available
integer, parameter :: maxUnitsPerRank = 100
!# Number of units per rank mpi

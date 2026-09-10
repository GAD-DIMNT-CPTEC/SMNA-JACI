! **Include File:** pre.h
!
! **Brief:** Defines constants for Module Mod_Pre
!
! **Author:** Denis Eiras; denis.eiras@inpe.br
!
! **Version:** 1.0.0
!
! **Changes:**
! - No items at this time
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
! https://opensource.org/licenses/GPL-3.0


integer, parameter :: maxModuleNameLength = 80
!# max length of a module name
integer, parameter :: maxNumberOfModules = 100
!# max number of modules may be implemented (Just raise this number if necessary)

! MPI Process compound string names
character (len = *), parameter :: fixedStringProc = "__PROC_#__"
character (len = *), parameter :: fixedStringProcGroup = "__STAND_ALONE_TOTAL_"
character (len = *), parameter :: fixedStringMasterProc = "__TOTAL__"

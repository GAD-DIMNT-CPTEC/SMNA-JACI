! **Include File:** porceClay.h
!
! **Brief:** Defines constants for Module Mod_PorceClayMaskIBIS
!
! **Author:** Denis Eiras; denis.eiras@inpe.br
!
! **Version:** 1.0.0
!
! **Changes:**
!  - 01-05-2019 - Denis Eiras             - version: 1.0
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

integer, parameter :: p_numVegClasses = 100
!# Number of vegetation classes
integer, dimension (p_numVegClasses) :: vegClasses = &
		(/1,1,1,1,1,1,1,1,1,1,&
		  2,2,2,2,2,2,2,2,2,2,&
		  3,3,3,3,3,3,3,3,3,3,&
		  4,4,4,4,4,4,4,4,4,4,&
		  5,5,5,5,5,5,5,5,5,5,&
		  6,6,6,6,6,6,6,6,6,6,&
		  7,7,7,7,7,7,7,7,7,7,&
		  8,8,8,8,8,8,8,8,8,8,&
		  9,9,9,9,9,9,9,9,9,9,&
		  2,2,2,2,2,2,2,2,2,2/)
!# Vegetation classes
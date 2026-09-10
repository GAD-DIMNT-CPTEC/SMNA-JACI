!------------------------------------------------------------------------------------------------!
!                           BUFR Types Module - Version 1.0                                      !
!------------------------------------------------------------------------------------------------!
!BOP                                                                                             !
! Description:                                                                                   !
!     This module defines a standard integer type (bufr_int) for BUFR libraries.                 !
!     It ensures compatibility with both 4-byte (INTEGER(4)) and 8-byte (INTEGER(8))             !
!     integer precisions, depending on the compilation flag BUFR_INT8.                           !
!                                                                                                !  
! Usage:                                                                                         !
!     - Include this module in Fortran source files that require precision control for integers. !
!     - The integer type bufr_int is set dynamically based on the compilation configuration.     !
!                                                                                                !
! Example:                                                                                       !
!     use bufr_types                                                                             !
!     integer(bufr_int) :: var                                                                   !
!                                                                                                !
! Compilation Options:                                                                           !
!     - Default: bufr_int is set to c_int (4-byte integers).                                     !
!     - If compiled with -DBUFR_INT8, bufr_int is set to c_long (8-byte integers).               !
!                                                                                                !
! Compilation Example:                                                                           !
!     gfortran -c myfile.F90                ! Uses 4-byte integers                               !
!     gfortran -DBUFR_INT8 -c myfile.F90     ! Uses 8-byte integers                              !
!                                                                                                !
! Notes:                                                                                         !
!     - This module allows BUFR libraries to be compiled for different precision levels.         !
!     - Ensures compatibility across multiple architectures.                                     !
!                                                                                                !
! Revisions:                                                                                     !
!     * 22-02-2025: de Mattos, J. G. - Initial module creation                                   !
!EOP                                                                                             !
!------------------------------------------------------------------------------------------------!

module bufr_types
  use iso_c_binding
#ifdef BUFR_INT8
  integer, parameter :: bufr_int = c_long  ! 8-byte integer when BUFR_INT8 is defined
#else
  integer, parameter :: bufr_int = c_int   ! Default: 4-byte integer
#endif
end module bufr_types

! **Include File:** precision.h
!
! **Brief:** Defines constants of number's precision
!
! **Author:** Denis Eiras; denis.eiras@inpe.br
!
! **Version:** 1.0.0
!
! **Changes:**
!  - 01-05-2019 - Denis Eiras - version: 1.0
!  - 08-09-2019 - Denis Eiras - version: 1.1
!  - 12-11-2019 - Denis Eiras - Task 6061 - precision adjustments: fixes TopoWaterGT30 reprodutibility - version: 1.2
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


! integer range 2^(4bytes*8bits) = 2^32 => smallest = -2^32/2, biggest = 2^32/2 -1 (-1 due to zero)
! integer range 2^(8bytes*8bits) = 2^64 => smallest = -2^64/2, biggest = 2^64/2 -1 (-1 due to zero)
! real range : first parameter represents minimum guaranteed size of precision, the second represents the minimum guaranteed size.
! PS: The origin of the base numbers below was not studied

! Values below are calculated on a intel64 i7
integer, parameter, public :: p_i4 = selected_int_kind(4)
!# smallest value = -2147483648              biggest value = +2147483647
integer, parameter, public :: p_i8 = selected_int_kind(8)
!# smallest value = -9223372036854775808     biggest value = +9223372036854775807
integer, parameter, public :: p_r4 = selected_real_kind(6,37)
!# smallest positive value = 1.17549435E-38           biggest value= 3.40282347E+38
integer, parameter, public :: p_r8 = selected_real_kind(15,307)
!# smallest positive value = 2.2250738585072014E-308  biggest value= 1.7976931348623157E+308

! Tupa + Cray doesnt accepts p_r16 = selected_real_kind(33, 4931). Legacy code always used p_r16 = p_r8
! integer, parameter, public :: p_r16 = selected_real_kind(33, 4931) ! smallest positive value = 3.36210314311209350626267781732175260E-4932  biggest value=   1.18973149535723176508575932662800702E+4932
integer, parameter, public :: p_r16 = selected_real_kind(15,307)
!# smallest positive value = 2.2250738585072014E-308  biggest value= 1.7976931348623157E+308

real (kind = p_r8), parameter , public :: p_undef = -999.0_p_r8
!# undef value for grADS ctl

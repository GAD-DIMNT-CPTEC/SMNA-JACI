!===============================================================================
! SVN $Id: shr_kind_mod.F90 41285 2012-10-26 01:46:45Z sacks $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_130528/shr/shr_kind_mod.F90 $
!===============================================================================

MODULE shr_kind_mod

   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   PUBLIC
   INTEGER,PARAMETER :: SHR_KIND_R8 = selected_real_kind(12) ! 8 byte real
   INTEGER,PARAMETER :: SHR_KIND_R4 = selected_real_kind( 6) ! 4 byte real
   INTEGER,PARAMETER :: SHR_KIND_RN = kind(1.0)              ! native real
   INTEGER,PARAMETER :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
   INTEGER,PARAMETER :: SHR_KIND_I4 = selected_int_kind ( 6) ! 4 byte integer
   INTEGER,PARAMETER :: SHR_KIND_IN = kind(1)                ! native integer
   INTEGER,PARAMETER :: SHR_KIND_CS = 80                     ! short char
   INTEGER,PARAMETER :: SHR_KIND_CL = 256                    ! long char
   INTEGER,PARAMETER :: SHR_KIND_CX = 512                    ! extra-long char
   INTEGER,PARAMETER :: SHR_KIND_CXX= 4096                   ! extra-extra-long char

END MODULE shr_kind_mod

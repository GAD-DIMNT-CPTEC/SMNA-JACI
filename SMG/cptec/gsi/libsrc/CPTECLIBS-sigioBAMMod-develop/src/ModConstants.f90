module ModConstants
   implicit none
   private

   !
   !precisao dos dados do modelo BAM
   integer, public, parameter :: I4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
   integer, public, parameter :: I8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
   integer, public, parameter :: R4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
   integer, public, parameter :: R8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers

   !Logical Units 
   integer,  public, parameter :: stderr = 0 ! Error Unit
   integer,  public, parameter :: stdinp = 5 ! Input Unit
   integer,  public, parameter :: stdout = 6 ! Output Unit


   ! Constants
   real(kind=r8), public, parameter :: rd      = 45.0/ATAN(1.0) ! convert to radian
   real(kind=r8), public, parameter :: emRad   = 6.37E6_r8
   real(kind=r8), public, parameter :: emRad2  = emRad*emRad
   real(kind=r8), public, parameter :: emRad1  = 1.0_r8/emRad
   real(kind=r8), public, parameter :: emRad12 = emRad1*emRad1
  
   real(kind=r8), public, parameter :: rad  = ATAN(1.0_r8)/45.0_r8

   integer, public, parameter :: strlen = 512!1024

   real(kind=r4), public, parameter :: undef = 1.0E-20
   
end module ModConstants

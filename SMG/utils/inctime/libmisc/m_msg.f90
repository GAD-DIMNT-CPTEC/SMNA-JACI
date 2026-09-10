!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!BOP
!  !MODULE: m_msg - a module to print some messages
!
!  !DESCRIPTION: 
!     TO DO
!
!  !INTERFACE:
!
! !REVISION HISTORY:
!
! 	19 may 2022 - J. G. Z. de Mattos - Initial Version
!
!----------------------------------------------------------------------!

MODULE m_msg
   use m_stdio, only: stdout, stderr

   implicit none
   private

   public :: perr
   public :: pwrn

   interface perr
      module procedure perr1_, perr2_
   end interface

   interface pwrn
      module procedure pwrn1_, pwrn2_
   end interface

   contains

   subroutine perr1_(where, message)
      implicit none
      character(len=*), intent(in) :: where
      character(len=*), intent(in) :: message
 
      write(stderr,'(4A)')'Error at ',trim(where),' : ', trim(message)
 
   end subroutine
 
   subroutine perr2_(where, message, cod)
      implicit none
      character(len=*), intent(in) :: where
      character(len=*), intent(in) :: message
      integer,          intent(in) :: cod
 
      write(stderr,'(4A,1x,I4)')'Error at ',trim(where),' : ', trim(message),cod
 
   end subroutine
 
   subroutine pwrn1_(where, message)
      implicit none
      character(len=*), intent(in) :: where
      character(len=*), intent(in) :: message
 
      write(stdout,'(4A)')'WARNING ',trim(where),' : ', trim(message)
 
   end subroutine
 
   subroutine pwrn2_(where, message1, message2)
      implicit none
      character(len=*), intent(in) :: where
      character(len=*), intent(in) :: message1, message2
 
      integer :: isize
      character(len=60) :: fmt
 
      isize=len_trim(where)
      write(fmt,'(A,I3.1,A)') '(A7,1x,A',isize,',1x,A,/,11x,A,/,11x,A)'
      write(stdout,fmt)'WARNING ',trim(where),' : ', trim(message1),trim(message2)
      
   end subroutine

END MODULE

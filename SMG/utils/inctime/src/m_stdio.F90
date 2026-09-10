module m_stdio
!
! A F90 module defines std. I/O parameters
!
! Description
!   Define system dependent I/O parameters.
!
! History
!   * 23 Mar 2011 - J.G. de Mattos - Initial Code
!
  implicit none
  private

  public   :: stdin    ! a unit linked to UNIX stdin
  public   :: stdout   ! a unit linked to UNIX stdout
  public   :: stderr   ! a unit linked to UNIX stderr

  ! Defines standar i/o units.

  integer, parameter :: stdin  = 5
  integer, parameter :: stdout = 6

  ! Generic setting for UNIX other than HP-UX

  integer, parameter :: stderr = 0

end module m_stdio

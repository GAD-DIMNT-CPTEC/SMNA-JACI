!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Get_xMax_yMax </br></br>
!#
!# **Brief**: Created to avoid unnecessary repetition of the GetImaxJmax subroutine in several InputParameters.f90 files that were removed (kept only in the Chopping folder). </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita          - version: 1.1.1 </li>
!#  <li>08-08-2019 - Eduardo Khamis  - version: 2.0.0 </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endtodo
!#
!# @documentation
!#
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2> </br>
!# **&#9993;**<mailto:atende.cptec@inpe.br> </br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---


module Mod_Get_xMax_yMax

  implicit none

  private

  public :: getxMaxyMax

  integer, parameter :: p_r8 = selected_real_kind(15)


contains


  subroutine getxMaxyMax (mEnd, xMax, yMax, linG)
    !# Gets xMax yMax
    !# ---
    !# @info
    !# **Brief:** Used to compute the value of xMax and yMax (old Imax e Jmax,
    !# renamed), when the process needs to make transforms.  </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: ago/2019 </br>
    !# @endinfo   
    implicit none

    integer, intent(in) :: mEnd
    integer, intent(out) :: xMax, yMax
    logical, intent(in), optional :: linG       
    !# Flag for linear (T) or quadratic (F) triangular truncation

    logical :: linearGrid
    integer :: nx, nm, n2m, n3m, n5m, n2, n3, n5, j, n, check, yfft
    integer, save :: lfft = 40000

    integer, dimension(:), allocatable, save :: xfft

    !    real(kind=p_r8) :: dl, dx, dKm = 112.0_p_r8

    if (present(linG)) then
      linearGrid = linG
    else
      linearGrid = .false.
    end if

    n2m = ceiling(log(real(lfft, p_r8)) / log(2.0_p_r8))
    n3m = ceiling(log(real(lfft, p_r8)) / log(3.0_p_r8))
    n5m = ceiling(log(real(lfft, p_r8)) / log(5.0_p_r8))
    nx = n2m * (n3m + 1) * (n5m + 1)

    allocate(xfft (nx))
    xfft = 0

    n = 0
    do n2 = 1, n2m
      yfft = (2**n2)
      if (yfft > lfft) exit
      do n3 = 0, n3m
        yfft = (2**n2) * (3**n3)
        if (yfft > lfft) exit
        do n5 = 0, n5m
          yfft = (2**n2) * (3**n3) * (5**n5)
          if (yfft > lfft) exit
          n = n + 1
          xfft(n) = yfft
        end do
      end do
    end do
    nm = n

    n = 0
    do
      check = 0
      n = n + 1
      do j = 1, nm - 1
        if (xfft(j) > xfft(j + 1)) then
          yfft = xfft(j)
          xfft(j) = xfft(j + 1)
          xfft(j + 1) = yfft
          check = 1
        end if
      end do
      if (check == 0) exit
    end do

    if (linearGrid) then
      yfft = 2
    else
      yfft = 3
    end if
    xMax = yfft * mEnd + 1
    do n = 1, nm
      if (xfft(n) >= xMax) then
        xMax = xfft(n)
        exit
      end if
    end do
    yMax = xMax / 2
    if (mod(yMax, 2) /= 0) yMax = yMax + 1

    deallocate(xfft)


    !!   For debuging :
    ! 
    !    if (linearGrid) then
    !       write(unit=*, fmt='(/,A)') &
    !             ' Linear Triangular Truncation : '
    !    else
    !       write(unit=*, fmt='(/,A)') &
    !             ' Quadratic Triangular Truncation : '
    !    end if
    !
    !    dl = 360.0_p_r8 / real(xMax, p_r8)
    !    dx = dl * dKm
    !    write(unit=*, fmt='(/,3(A,I5,/))') &
    !          ' mEnd : ', mEnd, ' xMax : ', xMax, ' yMax : ', yMax
    !    write(unit=*, fmt='(A,F13.9,A)')   ' dl: ', dl, ' degrees'
    !    write(unit=*, fmt='(A,F13.2,A,/)') ' dx: ', dx, ' km'

  end subroutine getxMaxyMax


end module Mod_Get_xMax_yMax

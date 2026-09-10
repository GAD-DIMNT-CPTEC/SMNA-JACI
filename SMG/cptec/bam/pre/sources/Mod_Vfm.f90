!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Vfm </br></br>
!#
!# **Brief**: Module responsible for reading vformat file</br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 1.0.1
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>??-??-???? - Jose P. Bonatti - version: 1.0.0</li>
!#  <li>26-04-2019 - Denis Eiras     - version: 1.0.1</li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time</li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time</li>
!# </ul>
!# @endtodo
!#
!# @documentation
!#
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2>
!# **&#9993;**<mailto:atende.cptec@inpe.br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---

module Mod_Vfm

  implicit none

  public :: vfirec
  private
  include 'precision.h'


contains

  ! ---------------------------------------------------------------------------
  !> @brief read vformat file
  !!
  !! @details read vformat file
  !!
  !! @param iunit vformat file unit for reading
  !! @param aa output array containing the file read
  !! @param n size of record
  !! @param type "LIN" - linear; "LOG" - logarithmic
  !!
  !! @author Denis M. A. Eiras
  !!
  !! @date mar/2019
  ! ---------------------------------------------------------------------------
  subroutine vfirec(iunit, a, n, type)

    integer, intent(in) :: iunit  !#TO deve ser kind default
    integer, intent(in) :: n
    real(kind = p_r8), intent(OUT) :: a(n)
    character(len = *), intent(in) :: type
    !
    ! local
    !
    character(len = 1) :: vc(0:63)
    character(len = 80) :: line
    character(len = 1) :: cs
    integer :: ich0
    integer :: ich9
    integer :: ichcz
    integer :: ichca
    integer :: ichla
    integer :: ichlz
    integer :: i
    integer :: nvalline
    integer :: nchs
    integer :: ic
    integer :: ii
    integer :: isval
    integer :: iii
    integer :: ics
    integer :: nn
    integer :: nbits
    integer :: nc
    real(kind = p_r8) :: bias
    real(kind = p_r8) :: fact
    real(kind = p_r8) :: facti
    real(kind = p_r8) :: scfct

    vc = '0'
    if (vc(0).ne.'0') call vfinit(vc)

    ich0 = ichar('0')
    ich9 = ichar('9')
    ichcz = ichar('Z')
    ichlz = ichar('z')
    ichca = ichar('A')
    ichla = ichar('a')

    read (iunit, '(2i8,2e20.10)')nn, nbits, bias, fact

    if (nn.ne.n) then
      print*, ' Word count mismatch on vfirec record '
      print*, ' Words on record - ', nn
      print*, ' Words expected  - ', n
      stop 'vfirec'
    end if

    nvalline = (78 * 6) / nbits
    nchs = nbits / 6

    do i = 1, n, nvalline
      read(iunit, '(a78)') line
      ic = 0
      do ii = i, i + nvalline - 1
        isval = 0
        if(ii.gt.n) EXIT
        do iii = 1, nchs
          ic = ic + 1
          cs = line(ic:ic)
          ics = ichar(cs)
          if (ics.le.ich9) then
            nc = ics - ich0
          else if (ics.le.ichcz) then
            nc = ics - ichca + 10
          else
            nc = ics - ichla + 36
          end if
          isval = ior(ishft(nc, 6 * (nchs - iii)), isval)
        end do ! loop iii
        a(ii) = isval
      end do ! loop ii

    end do ! loop i

    facti = 1.0_p_r8 / fact

    if (type.eq.'LIN') then
      do i = 1, n

        a(i) = a(i) * facti - bias

        !print*,'VFM=',i,a(i)
      end do
    else if (type.eq.'LOG') then
      scfct = 2.0_p_r8**(nbits - 1)
      do i = 1, n
        a(i) = sign(1.0_p_r8, a(i) - scfct)  &
          * (10.0_p_r8**(abs(20.0_p_r8 * (a(i) / scfct - 1.0_p_r8)) - 10.0_p_r8))
      end do
    end if
  end subroutine vfirec


  !--------------------------------------------------------
  subroutine vfinit(vc)
    character(len = 1), intent(OUT) :: vc   (*)
    character(len = 1) :: vcscr(0:63)
    integer :: n

    data vcscr/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'   &
      , 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'  &
      , 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T'  &
      , 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd'  &
      , 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n'  &
      , 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x'  &
      , 'y', 'z', '{', '|'/

    do n = 0, 63
      vc(n) = vcscr(n)
    end do
  end subroutine vfinit

end module Mod_Vfm

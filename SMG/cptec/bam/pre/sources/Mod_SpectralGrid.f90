!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SpectralGrid </br></br>
!#
!# **Brief**: Module responsible for generating the spectral representation of a
!# global field diagonalwise storage</br></br>
!# 
!# **Files in:**
!#
!# &bull; ?
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; ?
!# </br></br>
!# 
!# **Author**: Eduardo Khamis </br>
!#
!# **Version**: 1.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>25-09-2019 - Eduardo Khamis    - version: 1.0.0 </li>
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

module Mod_SpectralGrid

  use Mod_FastFourierTransform, only : invFFT, dirFFT

  use Mod_LegendreTransform, only : la0, spec2Four, four2Spec

  implicit none

  private

  public :: transp, specCoef2Grid, specCoef2GridD, grid2SpecCoef

  !parameters
  integer, parameter :: p_r8 = selected_real_kind(15)  
  !# kind for 64bits


contains


  subroutine transp(mnwv2In, mEnd1In, mEnd2In, qf)
    !# Transp
    !# ---
    !# @info
    !# **Brief:** Transp. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo 
    implicit none
    ! qf(Mnwv2) input: spectral representation of a
    !                  global field coluMnwise storage
    !          output: spectral representation of a
    !                  global field diagonalwise storage
    integer, intent(in) :: mnwv2In
    integer, intent(in) :: mEnd1In
    integer, intent(in) :: mEnd2In
    real(kind = p_r8), intent(inout) :: qf(mnwv2In)

    real(kind = p_r8) :: qw(mnwv2In)
    integer :: l
    integer :: lx
    integer :: mn
    integer :: mm
    integer :: nlast
    integer :: nn

    l = 0
    do mm = 1, mEnd1In
      nlast = mEnd2In - mm
      do nn = 1, nlast
        l = l + 1
        lx = la0(mm, nn)
        qw(2 * lx - 1) = qf(2 * l - 1)
        qw(2 * lx) = qf(2 * l)
      end do
    end do
    do mn = 1, mnwv2In
      qf(mn) = qw(mn)
    end do
  end subroutine transp


  subroutine specCoef2Grid(mnwv2In, mnwv3In, mEnd1In, mEnd2In, xMaxIn, yMaxIn, xmxIn, yMaxHfIn, qf, gf)
    !# Spectral Coeficient to Grid
    !# ---
    !# @info
    !# **Brief:** Spectral Coeficient to Grid. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo 
    implicit none
    integer, intent(in) :: mnwv2In
    integer, intent(in) :: mnwv3In
    integer, intent(in) :: mEnd1In
    integer, intent(in) :: mEnd2In
    integer, intent(in) :: xMaxIn
    integer, intent(in) :: yMaxIn
    integer, intent(in) :: xmxIn
    integer, intent(in) :: yMaxHfIn
    real(kind = p_r8), intent(in) :: qf(mnwv2In)
    real(kind = p_r8), intent(out) :: gf(xMaxIn, yMaxIn)
    real(kind = p_r8) :: four(xmxIn, yMaxIn)

    call spec2Four(mnwv2In, mnwv3In, mEnd1In, mEnd2In, yMaxHfIn, qf, four, .false.)
    call invFFT(four, gf)
  end subroutine SpecCoef2Grid


  subroutine specCoef2GridD(mnwv2In, mnwv3In, mEnd1In, mEnd2In, xMaxIn, yMaxIn, xmxIn, yMaxHfIn, qf, gf)
    !# Spectral Coeficient to GridD
    !# ---
    !# @info
    !# **Brief:** Spectral Coeficient to GridD. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo   
    implicit none
    integer, intent(in) :: mnwv2In
    integer, intent(in) :: mnwv3In
    integer, intent(in) :: mEnd1In
    integer, intent(in) :: mEnd2In
    integer, intent(in) :: xMaxIn
    integer, intent(in) :: yMaxIn
    integer, intent(in) :: xmxIn
    integer, intent(in) :: yMaxHfIn
    real(kind = p_r8), intent(in) :: qf(mnwv2In)
    real(kind = p_r8), intent(out) :: gf(xMaxIn, yMaxIn)
    real(kind = p_r8) :: four(xmxIn, yMaxIn)

    call spec2Four(mnwv2In, mnwv3In, mEnd1In, mEnd2In, yMaxHfIn, qf, four, .true.)
    call invFFT(four, gf)
    gf(:, yMaxHfIn + 1:yMaxIn) = -gf(:, yMaxHfIn + 1:yMaxIn)
  end subroutine specCoef2GridD


  subroutine grid2SpecCoef (mnwv2In, mEnd1In, xMaxIn, yMaxIn, xmxIn, yMaxHfIn, gf, qf)
    !# Grid to Spectral Coeficient
    !# ---
    !# @info
    !# **Brief:** Grid to Spectral Coeficient. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo  
    implicit none
    integer, intent(in) :: mnwv2In
    integer, intent(in) :: mEnd1In
    integer, intent(in) :: xMaxIn
    integer, intent(in) :: yMaxIn
    integer, intent(in) :: xmxIn
    integer, intent(in) :: yMaxHfIn
    real (kind = p_r8), intent(in) :: gf(xMaxIn, yMaxIn)
    real (kind = p_r8), intent(out) :: qf(mnwv2In)
    real (kind = p_r8) :: four(xmxIn, yMaxIn)

    call dirFFT (gf, Four)
    call four2Spec(mnwv2In, mEnd1In, yMaxHfIn, four, qf)
  end subroutine grid2SpecCoef


end module Mod_SpectralGrid

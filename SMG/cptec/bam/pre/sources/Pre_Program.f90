!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Program**: Global Pre Program - Main Program </br></br>
!#
!# **Brief**: Program used for generating input data for CPTEC Global Model </br></br>
!#
!#
!# **Author**: Eiras, Denis .M.A. </br>
!#
!# **Version**: 1.0.0 </br></br> 
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-02-2020 - Denis Eiras - version: 1.0.0 - Parallel Pre Program creation </li>
!# </ul>
!# 
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>No active bugs reported now </li>
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

program Pre_Program

  use Mod_Pre
  use Mod_Parallelism, only : &
    CreateParallelism, &
    DestroyParallelism
  use Mod_ModulesFacade, only : &
    ConcreteModulesFacade

  implicit none
  type(ConcreteModulesFacade) :: moduleFacade

  call CreateParallelism(.true.)
  call runPre(moduleFacade)
  call DestroyParallelism("Pre Program Parallelism Destroyed")
  call print_size()

contains

subroutine print_size
!Cray-JNT!
! Esta subrotina reporta o Virtual Memory High Water Mark (VmHWM)
! presente no arquivo /proc/self/status do sistema Linux.
!
integer :: iunit
character(len=80)            :: linha
  iunit = 7
  open(iunit,file="/proc/self/status")
  do
    read(iunit,'(a)',END=999) linha
    if (index(linha,'VmHWM').eq.1) then
       print *,linha
    endif
  enddo
999 close(iunit)
end subroutine print_size


end program Pre_Program

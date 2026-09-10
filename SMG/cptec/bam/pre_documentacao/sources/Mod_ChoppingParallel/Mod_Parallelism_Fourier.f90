!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Parallelism_Fourier </br></br>
!#
!# **Brief**: Deals with parallelism of Fourier in Chopping. </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!# 
!# **Author**: Denis Eiras </br>
!#
!# **Version**: 1.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>03-06-2019 - Denis Eiras - version 1.0.0 - Moving Funcionality from Mod_Parallelism to this module </li>
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

module Mod_Parallelism_Fourier

  use Mod_Parallelism_Group_Chopping, only : &
    mpiCommGroup

  implicit none

  private

  ! public data
  integer, public :: maxNodes_four 
  !# MPI processes in fourier group
  integer, public :: myId_four     
  !# MPI process rank in fourier group
  integer, public :: mygroup_four  
  !# fourier group
  integer, public :: COMM_FOUR     
  !# Communicator of Fourier Group

  ! public procedures
  public :: CreateFourierGroup

  include 'mpif.h'

contains


  subroutine CreateFourierGroup(mygroup_four, myid_four)
    !# Creates MPI Fourier Group
    !# ---
    !# @info
    !# **Brief:** Creates MPI Fourier Group. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 <br>
    !# @endin
    integer, intent(in) :: mygroup_four
    integer, intent(in) :: myid_four

    integer :: ierror

    call MPI_COMM_SPLIT(mpiCommGroup, mygroup_four, myid_four, COMM_FOUR, ierror)

  end subroutine CreateFourierGroup

end module Mod_Parallelism_Fourier


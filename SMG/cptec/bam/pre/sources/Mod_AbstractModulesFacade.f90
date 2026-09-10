!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_AbstractModulesFacade </br></br>
!#
!# **Brief**: Abstract Facade for overriding facades Pre Modules</br></br>
!# 
!# **Author**: Denis Eiras </br>
!#
!# **Version**: 1.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>27-08-2021 - Denis Eiras    - version: 1.0.0</li></br>
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
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2></br>
!# **&#9993;**<mailto:atende.cptec@inpe.br></br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---

module Mod_AbstractModulesFacade

   implicit none
   private

   !# Abstract class. Must be implemented for real (ConcreteModulesFacade)  and for tests (MockModulesFacade)
   type AbstractModulesFacade
   contains
      procedure :: runModule => abstractRunModule
      procedure :: initModule => abstractInitModule
      procedure :: dependencyShouldRun => abstractDependencyShouldRun
   end type AbstractModulesFacade

   character(len = *), parameter :: headerMsg = 'AbstractModulesFacade Module  | '


   public :: AbstractModulesFacade


contains


   subroutine abstractInitModule(this, moduleName)
      implicit none
      class(AbstractModulesFacade) :: this
      character (len = *), intent(in) :: moduleName

      print*, "$" // headerMsg // "Abstract InitModule must be implemented"
      stop
   end subroutine abstractInitModule


   subroutine abstractRunModule(this, moduleName)
      implicit none
      class(AbstractModulesFacade) :: this
      character (len = *), intent(in) :: moduleName

      print*, "$" // headerMsg // "Abstract runModule must be implemented"

   end subroutine abstractRunModule


   function abstractDependencyShouldRun(this, moduleName) result(shouldRun)
      implicit none
      !# Super method to be implemented by subclasses
      !# ---
      !# @info
      !# **Brief:**  Super method to be implemented by subclasses
      !# Checks if dependency modules should run, even if dependency was
      !# not in run list. </br>
      !# **Authors**: </br>
      !# &bull; Denis Eiras </br>
      !# **Date**: apr/2019 </br>
      !# @endinfo
      class(AbstractModulesFacade) :: this
      character (len = *), intent(in) :: moduleName
      logical :: shouldRun

      print*, "$" // headerMsg // "Abstract dependencyShouldRun must be implemented"
      stop
   end function abstractDependencyShouldRun


end module Mod_AbstractModulesFacade

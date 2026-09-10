!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_InputArrays </br></br>
!#
!# **Brief**: Module used for Get Input Arrays. </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!# 
!# **Author**: Paulo Bonatti </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Paulo Bonatti  - version: 1.0.0 </li>
!#  <li>01-08-2007 - Simone Tomita  - version: 1.1.1 </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 2.0.0 - some adaptations for modularizing Chopping </li>
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


module Mod_InputArrays

  use Mod_InputParameters, only : KmaxInpp, NTracers, &
    ImaxOut, JmaxOut, KmaxOutp, ImaxInp, jmaxInp
  use Mod_Sizes, only : ibMax, jbMax, kMax, kMaxloc_out, kmaxloc_in, &
    iMax, jMax, jMaxHalf, jMinPerM, iMaxPerJ, &
    mMax, nExtMax, mnMax, mnExtMax, mnExtMap, &
    mymnMax, mymnExtMax, mnMax_out, mnMax_out

  implicit none

  private
  include 'precision.h'

  public :: GetArrays, ClsArrays, GetSpHuTracers, ClsSpHuTracers

  integer (kind = p_i8), dimension (4), public :: DateInitial, DateCurrent

  real (kind = p_r4), dimension (:), allocatable, public :: &
    DelSInp, SigIInp, SigLInp, SigIOut, SigLOut

  real (kind = p_r4), dimension (:), allocatable, public :: qWorkInp, qWorkprOut, qtorto
  real (kind = p_r4), dimension (:,:,:), allocatable, public :: qWorkInp3D
  real (kind = p_r8), dimension (:, :), allocatable, public :: qWorkOut, qWorkOut1, qWorkInOut, qWorkInOut1

  real (kind = p_r8), dimension (:), allocatable, public :: &
    DelSigmaInp, SigInterInp, SigLayerInp, &
    DelSigmaOut, SigInterOut, SigLayerOut

  real (kind = p_r8), dimension (:), allocatable, public :: &
    qTopoInp, qLnPsInp, qTopoOut, qTopoOutSpec, qLnPsOut

  real (kind = p_r8), dimension (:, :), allocatable, public :: &
    qDivgInp, qVortInp, qTvirInp, gWorkprInp, &
    qDivgOut, qVortOut, qTvirOut, &
    qUvelInp, qVvelInp, qUvelOut, qVvelOut

  real (kind = p_r8), dimension (:, :, :), allocatable, public :: &
    qSpHuInp, qSpHuOut

  real (kind = p_r8), dimension (:, :, :), allocatable, public :: gWorkOut
  real (kind = p_r4), dimension (:, :), allocatable, public :: gWorkprout

  real (kind = p_r8), dimension (:, :), allocatable, public :: &
    gTopoInp, gTopoOut, gTopoOutGaus, gTopoOutGaus8, gTopoDel, &
    gLnPsInp, gPsfcInp, gLnPsOut, gPsfcOut, gpresaux

  real (kind = p_r8), dimension (:, :, :), allocatable, public :: &
    gUvelInp, gVvelInp, gTvirInp, &
    gDivgInp, gVortInp, gPresInp, gPresInpp, &
    gUvelOut, gVvelOut, gTvirOut, &
    gPresOut

  real (kind = p_r8), dimension (:, :, :, :), allocatable, public :: &
    gSpHuInp, gSpHuOut

  integer :: Mnwv2Rec, Mnwv2Out, Mnwv2Inp, Mnwv3Inp, Mnwv3Out

contains


  subroutine GetArrays(kMaxInp, kMaxOut)
    !# Gets Input Arrays
    !# ---
    !# @info
    !# **Brief:** Gets Input Arrays. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    implicit none
    integer, intent(in) :: kMaxInp
    integer, intent(in) :: kMaxOut

    Mnwv2Inp = 2 * mymnMax
    Mnwv2Out = Mnwv2Inp
    Mnwv3Inp = 2 * mymnExtMax
    Mnwv3Out = Mnwv3Inp
    Mnwv2Rec = Mnwv2Inp
    !  ALLOCATE (qTopoRec    (Mnwv2Rec))
    allocate (qWorkInp    (2 * mnmax))
    allocate (qWorkInp3D  (2 * mnmax,kMaxInp,2))
    !  ALLOCATE (qWorkInpAux (2*mnMax_out))
    allocate (qWorkOut    (2 * mnmax_Out, kMaxOut))
    allocate (qWorkOut1   (2 * mnmax_Out, kMaxOut))
    allocate (qWorkInOut  (2 * mnmax_Out, kMaxInp))
    allocate (qWorkInOut1 (2 * mnmax_Out, kMaxInp))
    allocate (qWorkprOut  (2 * mnmax_Out))
    allocate (qtorto      (2 * mnmax_Out))
    allocate (SigIInp     (KmaxInpp), SigLInp (KmaxInpp), DelSInp (kMaxInp))
    allocate (SigIOut     (KmaxOutp), SigLOut (KmaxOutp))
    allocate (DelSigmaInp (KmaxInpp), SigInterInp (KmaxInpp), SigLayerInp (KmaxInpp))
    allocate (DelSigmaOut (kMaxOut), SigInterOut (KmaxOutp), SigLayerOut (KmaxOutp))
    allocate (qTopoInp    (Mnwv2Inp), qLnPsInp (Mnwv2Inp))
    allocate (qTopoOut    (Mnwv2Out), qTopoOutSpec(Mnwv2Out), qLnPsOut (Mnwv2Out))
    allocate (qDivgInp (Mnwv2Inp, Kmaxloc_In), &
      qVortInp (Mnwv2Inp, Kmaxloc_In))
    allocate (qDivgOut (Mnwv2Out, Kmaxloc_Out), &
      qVortOut (Mnwv2Out, Kmaxloc_Out))
    allocate (qUvelInp (Mnwv3Out, Kmaxloc_In), &
      qVvelInp (Mnwv3Out, Kmaxloc_In))
    allocate (qUvelOut (Mnwv3Out, Kmaxloc_Out), &
      qVvelOut (Mnwv3Out, Kmaxloc_Out))
    allocate (qTvirInp (Mnwv2Inp, Kmaxloc_In))
    allocate (qTvirOut (Mnwv2Out, Kmaxloc_Out))
    allocate (gWorkOut (ImaxOut, JmaxOut, max(kMaxOut, kMaxInp)))
    allocate (gWorkprOut (ImaxOut, JmaxOut), &
      gWorkprInp (ImaxInp, jmaxInp))
    allocate (gTopoInp (Ibmax, jbmax), &
      gTopoOut (Ibmax, jbmax), &
      gTopoOutGaus(Ibmax, jbmax), &
      gTopoOutGaus8(ImaxOut, JmaxOut), &
      gTopoDel (Ibmax, jbmax), &
      gLnPsInp (Ibmax, jbmax), &
      gLnPsOut (Ibmax, jbmax), &
      gPsfcInp (Ibmax, jbmax), &
      gPsfcOut (Ibmax, jbmax))
    allocate (gpresaux (Ibmax, KmaxOutp))
    allocate (gUvelInp (Ibmax, kMaxInp, Jbmax), &
      gVvelInp (Ibmax, kMaxInp, Jbmax))
    allocate (gDivgInp (Ibmax, kMaxInp, Jbmax), &
      gVortInp (Ibmax, kMaxInp, Jbmax))
    allocate (gTvirInp (Ibmax, kMaxInp, Jbmax), &
      gPresInp (Ibmax, kMaxInp, Jbmax), &
      gPresInpp(Ibmax, KmaxInpp, Jbmax))
    allocate (gUvelOut (Ibmax, kMaxOut, Jbmax), &
      gVvelOut (Ibmax, kMaxOut, Jbmax))
    allocate (gTvirOut (Ibmax, kMaxOut, Jbmax), &
      gPresOut (Ibmax, kMaxOut, Jbmax))

    gPsfcOut = 0.0_p_r8

    !  DateInitial=0;DateCurrent=0

    !  DelSInp=0.0_r4; SigIInp=0.0_r4;SigLInp=0.0_r4;SigIOut=0.0_r4; SigLOut=0.0_r4

    !  qWorkInp=0.0_r4; qWorkprOut=0.0_r4; qtorto=0.0_r4;
    !  qWorkOut=0.0_r8; qWorkOut1=0.0_r8;qWorkInOut=0.0_r8; qWorkInOut1=0.0_r8

    !  DelSigmaInp=0.0_r8; SigInterInp=0.0_r8; SigLayerInp=0.0_r8
    !  DelSigmaOut=0.0_r8; SigInterOut=0.0_r8; SigLayerOut=0.0_r8

    !  qTopoInp=0.0_r8; qLnPsInp=0.0_r8; qTopoOut=0.0_r8;qTopoOutSpec=0.0_r8
    !  qLnPsOut=0.0_r8

    !  qDivgInp=0.0_r8; qVortInp=0.0_r8; qTvirInp=0.0_r8; gWorkprInp=0.0_r8
    !  qDivgOut=0.0_r8; qVortOut=0.0_r8
    !  qTvirOut=0.0_r8
    !  qUvelInp=0.0_r8; qVvelInp=0.0_r8; qUvelOut=0.0_r8
    !  qVvelOut=0.0_r8


    !  gWorkOut=0.0_r8
    !  gWorkprout=0.0_r4

    !  gTopoInp=0.0_r8; gTopoOut=0.0_r8;gTopoOutGaus=0.0_r8
    !  gTopoOutGaus8=0.0_r8;gTopoDel=0.0_r8
    !  gLnPsInp=0.0_r8; gPsfcInp=0.0_r8; gLnPsOut=0.0_r8; gPsfcOut=0.0_r8; gpresaux=0.0_r8

    !  gUvelInp=0.0_r8; gVvelInp=0.0_r8; gTvirInp=0.0_r8
    !  gDivgInp=0.0_r8; gVortInp=0.0_r8; gPresInp=0.0_r8; gPresInpp=0.0_r8
    !  gUvelOut=0.0_r8; gVvelOut=0.0_r8; gTvirOut=0.0_r8
    !  gPresOut=0.0_r8;

  end subroutine GetArrays


  subroutine ClsArrays
    !# Cleans Arrays
    !# ---
    !# @info
    !# **Brief:** Cleans Arrays. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    implicit none

    !  DEALLOCATE (qTopoRec )
    deallocate (qWorkInp)
    deallocate (qWorkInp3D) 
    !  DEALLOCATE (qWorkInpAux )
    deallocate (qWorkOut)
    deallocate (qWorkOut1)
    deallocate (qWorkInOut)
    deallocate (qWorkInOut1)
    deallocate (qWorkprOut)
    deallocate (qtorto)
    deallocate (SigIInp, SigLInp, DelSInp)
    deallocate (SigIOut, SigLOut)
    deallocate (DelSigmaInp, SigInterInp, SigLayerInp)
    deallocate (DelSigmaOut, SigInterOut, SigLayerOut)
    deallocate (qTopoInp, qLnPsInp)
    deallocate (qTopoOut, qTopoOutSpec, qLnPsOut)
    deallocate (qDivgInp, &
      qVortInp)
    deallocate (qDivgOut, &
      qVortOut)
    deallocate (qUvelInp, &
      qVvelInp)
    deallocate (qUvelOut, &
      qVvelOut)
    deallocate (qTvirInp)
    deallocate (qTvirOut)
    deallocate (gWorkOut)
    deallocate (gWorkprOut, &
      gWorkprInp)
    deallocate (gTopoInp, &
      gTopoOut, &
      gTopoOutGaus, &
      gTopoOutGaus8, &
      gTopoDel, &
      gLnPsInp, &
      gLnPsOut, &
      gPsfcInp, &
      gPsfcOut)
    deallocate (gpresaux)
    deallocate (gUvelInp, &
      gVvelInp)
    deallocate (gDivgInp, &
      gVortInp)
    deallocate (gTvirInp, &
      gPresInp, &
      gPresInpp)
    deallocate (gUvelOut, &
      gVvelOut)
    deallocate (gTvirOut, &
      gPresOut)

  end subroutine ClsArrays


  subroutine GetSpHuTracers(kMaxInp, kMaxOut)
    !# Gets SpHu Tracers
    !# ---
    !# @info
    !# **Brief:** Gets SpHu Tracers. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    implicit none
    integer, intent(in) :: kMaxInp
    integer, intent(in) :: kMaxOut

    if (.not.allocated(qSpHuInp))then
      allocate (qSpHuInp (Mnwv2Inp, kMaxInp, NTracers + 2))
    else
      deallocate (qSpHuInp)
      allocate (qSpHuInp (Mnwv2Inp, kMaxInp, NTracers + 2))
    end if

    if (.not.allocated(qSpHuOut))then
      allocate (qSpHuOut (Mnwv2Out, kMaxOut, NTracers + 2))
    else
      deallocate (qSpHuOut)
      allocate (qSpHuOut (Mnwv2Out, kMaxOut, NTracers + 2))
    end if

    if (.not.allocated(gSpHuInp))then
      allocate (gSpHuInp (Ibmax, kMaxInp, Jbmax, NTracers + 2))
    else
      deallocate (gSpHuInp)
      allocate (gSpHuInp (Ibmax, kMaxInp, Jbmax, NTracers + 2))
    end if

    if (.not.allocated(gSpHuOut))then
      allocate (gSpHuOut (Ibmax, kMaxOut, Jbmax, NTracers + 2))
    else
      deallocate (gSpHuOut)
      allocate (gSpHuOut (Ibmax, kMaxOut, Jbmax, NTracers + 2))
    end if

  end subroutine GetSpHuTracers


  subroutine ClsSpHuTracers
    !# Cleans SpHu Tracers
    !# ---
    !# @info
    !# **Brief:** Cleans SpHu Tracers. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    implicit none
    if (allocated(qSpHuInp))then
      deallocate (qSpHuInp)
    end if

    if (allocated(qSpHuOut))then
      deallocate (qSpHuOut)
    end if

    if (allocated(gSpHuInp))then
      deallocate (gSpHuInp)
    end if

    if (allocated(gSpHuOut))then
      deallocate (gSpHuOut)
    end if

  end subroutine ClsSpHuTracers


end module Mod_InputArrays

! ------------------------- ------------------------------------------------------
! INPE/CPTEC, DMD, Modelling and Development Division
! -------------------------------------------------------------------------------
!
! Include File: pre.h
!
! REVISION HISTORY:
! 01-05-2019 - Denis Eiras             - version: 1.0
!
!> @author
!> Denis M. A. Eiras (last revision) \n
!!
!> @brief Type used in ChoppingNameList, configured in PRE.nml, which contains all Chopping parameters
!!
!! @version 1.0
!! @date 01-04-2019
!!
!! @copyright Under GLP-3.0
!! @link: https://opensource.org/licenses/GPL-3.0
! -------------------------------------------------------------------------------
! ---------------------------------------------------------------------------

include 'pre.h'
type ChoppingNameListData
  integer :: mEndInp, kMaxInp, mEndOut, kMaxOut, mEndMin, mEndCut, iter, nProc_vert, ibdim_size, tamBlock
  real (kind = p_r8) :: smthPerCut
  logical :: getOzone, getTracers, grADS, grADSOnly, gdasOnly, smoothTopo, rmGANL, linearGrid, givenfouriergroups
  character (len = 10) :: dateLabel
  character (len = 2) :: utc
  character (len = 16) :: nCepName
  character (len = 500) :: dataGDAS, gdasInp
  character(len = maxPathLength) :: dirInp, dirOut, dirTop, dirSig, dirGrd, dGDInp
  character(len = maxPathLength) :: dirPreOut = './'                       !< output data directory
end type ChoppingNameListData
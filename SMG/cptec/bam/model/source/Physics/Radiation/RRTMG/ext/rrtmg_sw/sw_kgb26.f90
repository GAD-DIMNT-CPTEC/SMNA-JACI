!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:13 $

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! **************************************************************************
!      subroutine sw_kgbnn
! **************************************************************************
!  RRTM Shortwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original by J.Delamere, Atmospheric & Environmental Research.
!  Reformatted for F90: JJMorcrette, ECMWF
!  Further F90 and GCM revisions:  MJIacono, AER, July 2002
!
!  This file contains 14 subroutines that include the 
!  absorption coefficients and other data for each of the 14 shortwave
!  spectral bands used in RRTM_SW.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in routine RRTMG_SW_INIT to reduce
!  the total number of g-points from 224 to 112 for use in the GCM.
! **************************************************************************
! **************************************************************************
      subroutine sw_kgb26
! **************************************************************************

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb 
      use rrsw_kg26, only : sfluxrefo, raylo

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:) = (/ &
!         &     129.462_r8, 15*0._r8 /)
        &   29.0079_r8,  28.4088_r8,     20.3099_r8,  13.0283_r8 &
        &,  11.8619_r8,  9.95840_r8,     6.68696_r8,  5.38987_r8 &
        &,  3.49829_r8, 0.407693_r8,    0.299027_r8, 0.236827_r8 &
        &, 0.188502_r8, 0.163489_r8, 4.64335e-02_r8, 2.72662e-03_r8 /)

! Rayleigh extinction coefficient at all v 
      raylo(:) = (/ &
        &  1.21263e-06_r8,1.43428e-06_r8,1.67677e-06_r8,1.93255e-06_r8 &
        &, 2.19177e-06_r8,2.44195e-06_r8,2.66926e-06_r8,2.85990e-06_r8 &
        &, 3.00380e-06_r8,3.06996e-06_r8,3.08184e-06_r8,3.09172e-06_r8 &
        &, 3.09938e-06_r8,3.10456e-06_r8,3.10727e-06_r8,3.10818e-06_r8 /)

      end subroutine sw_kgb26

!**************************************************************************

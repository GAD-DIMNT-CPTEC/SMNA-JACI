!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/04/24 16:17:27 $
!

       module rrtmg_lw_rad

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                              RRTMG_LW                                    *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                   a rapid radiative transfer model                       *
! *                       for the longwave region                            * 
! *             for application to general circulation models                *
! *                                                                          *
! *                                                                          *
! *            Atmospheric and Environmental Research, Inc.                  *
! *                        131 Hartwell Avenue                               *
! *                        Lexington, MA 02421                               *
! *                                                                          *
! *                                                                          *
! *                           Eli J. Mlawer                                  *
! *                        Jennifer S. Delamere                              *
! *                         Michael J. Iacono                                *
! *                         Shepard A. Clough                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                       email:  miacono@aer.com                            *
! *                       email:  emlawer@aer.com                            *
! *                       email:  jdelamer@aer.com                           *
! *                                                                          *
! *        The authors wish to acknowledge the contributions of the          *
! *        following people:  Steven J. Taubman, Karen Cady-Pereira,         *
! *        Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.  *
! *                                                                          *
! ****************************************************************************

! -------- Modules --------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb 
      use rrlw_vsn
      use mcica_subcol_gen_lw, only: mcica_subcol_lw
      use rrtmg_lw_cldprop, only: cldprop
      use rrtmg_lw_cldprmc, only: cldprmc
! Move call to rrtmg_lw_ini and following use association to 
! GCM initialization area
!      use rrtmg_lw_init, only: rrtmg_lw_ini
      use rrtmg_lw_rtrn, only: rtrn
      use rrtmg_lw_rtrnmr, only: rtrnmr
      use rrtmg_lw_rtrnmc, only: rtrnmc
      use rrtmg_lw_setcoef, only: setcoef
      use rrtmg_lw_taumol, only: taumol

      implicit none

! public interfaces/functions/subroutines
      public :: rrtmg_lw, inatm

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------
!  call rrtmg_lw(&
!       imca                                              , &
!       pcols                                             , &
!       ncol                                              , &
!       rrtmg_levs                                        , &
!       icld                                              , &
!       r_state_pmidmb  (1:pcols,1:pver)                  , &
!       r_state_pintmb  (1:pcols,1:pver+1)                , &
!       r_state_tlay    (1:pcols,1:pver)                  , &
!       r_state_tlev    (1:pcols,1:pver+1)                , &
!       tsfc            (1:pcols)                         , &
!       r_state_h2ovmr  (1:pcols,1:pver)                  , &
!       r_state_o3vmr   (1:pcols,1:pver)                  , &
!       r_state_co2vmr  (1:pcols,1:pver)                  , &
!       r_state_ch4vmr  (1:pcols,1:pver)                  , &
!       r_state_o2vmr   (1:pcols,1:pver)                  , &
!       r_state_n2ovmr  (1:pcols,1:pver)                  , &
!       r_state_cfc11vmr(1:pcols,1:pver)                  , &
!       r_state_cfc12vmr(1:pcols,1:pver)                  , &
!       r_state_cfc22vmr(1:pcols,1:pver)                  , &
!       r_state_ccl4vmr (1:pcols,1:pver)                  , &
!       emis            (1:pcols,1:nbndlw)                , &
!       inflglw                                           , &
!       iceflglw                                          , &
!       liqflglw                                          , &
!       cld_stolw       (1:nsubclw,1:pcols,1:rrtmg_levs-1), &
!       tauc_stolw      (1:nsubclw,1:pcols,1:rrtmg_levs-1), &
!       cicewp_stolw    (1:nsubclw,1:pcols,1:rrtmg_levs-1), &
!       cliqwp_stolw    (1:nsubclw,1:pcols,1:rrtmg_levs-1), &
!       cld             (1:pcols,1:pver)                  , &
!       cicewp          (1:pcols,1:rrtmg_levs-1)          , &
!       cliqwp          (1:pcols,1:rrtmg_levs-1)          , &
!       tauc_lw         (1:nbndlw,1:pcols,1:pver)         , &
!       rei             (1:pcols,1:rrtmg_levs-1)          , &
!       rel             (1:pcols,1:rrtmg_levs-1)          , &
!       taua_lw         (1:pcols,1:rrtmg_levs-1,1:nbndlw) , &
!       uflx            (1:pcols,1:rrtmg_levs+1)          , &
!       dflx            (1:pcols,1:rrtmg_levs+1)          , &
!       hr              (1:pcols,1:rrtmg_levs)            , &
!       uflxc           (1:pcols,1:rrtmg_levs+1)          , &
!       dflxc           (1:pcols,1:rrtmg_levs+1)          , &
!       hrc             (1:pcols,1:rrtmg_levs)            , &
!       lwuflxs         (1:nbndlw,1:pcols,1:pverp+1)      , &
!       lwdflxs         (1:nbndlw,1:pcols,1:pverp+1)        )

      subroutine rrtmg_lw(&
       imca            , &!integer      , intent(in   ) :: imca                          ! flag for mcica [0=off, 1=on]
       pcols           , &!integer      , intent(in   ) :: pcols
       ncol            , &!integer      , intent(in   ) :: ncol              ! Number of horizontal columns
       nlay            , &!integer      , intent(in   ) :: nlay              ! Number of model layers
       icld            , &!integer      , intent(inout) :: icld              ! Cloud overlap method
       play            , &!real(kind=r8), intent(in   ) :: play(ncol,nlay)   ! Layer pressures (hPa, mb)
       plev            , &!real(kind=r8), intent(in   ) :: plev(ncol,nlay+1) ! Interface pressures (hPa, mb)
       tlay            , &!real(kind=r8), intent(in   ) :: tlay(ncol,nlay)   ! Layer temperatures (K)
       tlev            , &!real(kind=r8), intent(in   ) :: tlev(ncol,nlay+1)    ! Interface temperatures (K)
       tsfc            , &!real(kind=r8), intent(in   ) :: tsfc(ncol)        ! Surface temperature (K)
       h2ovmr          , &!real(kind=r8), intent(in) :: h2ovmr (ncol,nlay)          ! H2O volume mixing ratio
       o3vmr           , &! real(kind=r8), intent(in) :: o3vmr (ncol,nlay)           ! O3 volume mixing ratio
       co2vmr          , &!real(kind=r8), intent(in) :: co2vmr (ncol,nlay)          ! CO2 volume mixing ratio
       ch4vmr          , &!real(kind=r8), intent(in) :: ch4vmr (ncol,nlay)         ! Methane volume mixing ratio
       o2vmr           , &!real(kind=r8), intent(in) :: o2vmr (ncol,nlay)          ! O2 volume mixing ratio
       n2ovmr          , &!real(kind=r8), intent(in) :: n2ovmr (ncol,nlay)          ! Nitrous oxide volume mixing ratio
       cfc11vmr        , &!real(kind=r8), intent(in) :: cfc11vmr (ncol,nlay)        ! CFC11 volume mixing ratio
       cfc12vmr        , &!real(kind=r8), intent(in) :: cfc12vmr (ncol,nlay)        ! CFC12 volume mixing ratio
       cfc22vmr        , &!real(kind=r8), intent(in) :: cfc22vmr (ncol,nlay)        ! CFC22 volume mixing ratio
       ccl4vmr         , &!real(kind=r8), intent(in) :: ccl4vmr (ncol,nlay)         ! CCL4 volume mixing ratio
       emis            , &!real(kind=r8), intent(in) :: emis (ncol,nbndlw)          ! Surface emissivity
       inflglw         , &!integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
       iceflglw        , &!integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
       liqflglw        , &!integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification
       cldfmcl         , &!real(kind=r8), intent(in) :: cldfmcl(ngptlw,ncol,nlay-1)       ! Cloud fraction
       taucmcl         , &!real(kind=r8), intent(in) :: taucmcl (ngptlw,ncol,nlay-1)       ! Cloud optical depth
       ciwpmcl         , &!real(kind=r8), intent(in) :: ciwpmcl(ngptlw,ncol,nlay-1)       ! Cloud ice water path (g/m2)
       clwpmcl         , &!real(kind=r8), intent(in) :: clwpmcl (ngptlw,ncol,nlay-1)      ! Cloud liquid water path (g/m2)
       cld             , &!real(kind=r8), intent(in) :: cld(ncol,nlay)       ! layer cloud fraction
       cicewp          , &!real(kind=r8), intent(in) :: cicewp(ncol,nlay-1) 
       cliqwp          , &!real(kind=r8), intent(in) :: cliqwp(ncol,nlay-1)  
       tauc_lw         , &!real(kind=r8), intent(in) :: tauc_lw(nbndlw,ncol,nlay)        ! cloud optical depth
       reicmcl         , &!real(kind=r8), intent(in) :: reicmcl (ncol,nlay-1)         ! Cloud ice effective radius (microns)
       relqmcl         , &!real(kind=r8), intent(in) :: relqmcl (ncol,nlay-1)         ! Cloud water drop effective radius (microns)
       tauaer          , &!real(kind=r8), intent(in) :: tauaer (ncol,nlay-1,nbndlw)        ! aerosol optical depth
       uflx            , &!real(kind=r8), intent(inout) :: uflx(ncol,nlay+1)           ! Total sky longwave upward flux (W/m2)
       dflx            , &!real(kind=r8), intent(inout) :: dflx(ncol,nlay+1)           ! Total sky longwave downward flux (W/m2)
       hr              , &!real(kind=r8), intent(inout) :: hr (ncol,nlay)             ! Total sky longwave radiative heating rate (K/d)
       uflxc           , &! real(kind=r8), intent(inout) :: uflxc (ncol,nlay+1)          ! Clear sky longwave upward flux (W/m2)
       dflxc           , &!real(kind=r8), intent(inout) :: dflxc (ncol,nlay+1)          ! Clear sky longwave downward flux (W/m2)
       hrc             , &!real(kind=r8), intent(inout) :: hrc(ncol,nlay)           ! Clear sky longwave radiative heating rate (K/d)
       uflxs           , &!real(kind=r8), intent(inout) :: uflxs (nbndlw,ncol,nlay+1)        ! Total sky longwave upward flux spectral (W/m2)
       dflxs             )!real(kind=r8), intent(inout) :: dflxs (nbndlw,ncol,nlay+1)        ! Total sky longwave downward flux spectral (W/m2)

!      subroutine rrtmg_lw &
!            (ncol    ,nlay    ,icld    ,                   &
!             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
!             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,&
!             cfc11vmr,cfc12vmr, &
!             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
!             cldfmcl ,taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
!             tauaer  , &
!             uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc     , &
!             uflxs   ,dflxs   , &
!             pcols)

! -------- Description --------

! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
! model for application to GCMs, that has been adapted from RRTM_LW for
! improved efficiency.
!
! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
!  area, since this has to be called only once. 
!
! This routine:
!    a) calls INATM to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    b) calls CLDPRMC to set cloud optical depth for McICA based 
!       on input cloud properties 
!    c) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    d) calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands
!    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation using McICA, the Monte-Carlo 
!       Independent Column Approximation, to represent sub-grid scale 
!       cloud variability
!    f) passes the necessary fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
!     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!
! This call to RRTMG_LW must be preceeded by a call to the module
!     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngpt) dimension.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
!     and subroutine rrtmg_lw_cldprop.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
!    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflglw and liqflglw
!
! One method of aerosol property input is possible:
!     Aerosol properties can be input in only one way (controlled by input 
!     flag iaer, see text file rrtmg_lw_instructions for further details):
!
!    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
!       band average optical depth at the mid-point of each spectral band.
!       RRTMG_LW currently treats only aerosol absorption;
!       scattering capability is not presently available. 
!
!
! ------- Modifications -------
!
! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
! set of g-points for application to GCMs.  
!
!-- Original version (derived from RRTM_LW), reduction of g-points, other
!   revisions for use with GCMs.  
!     1999: M. J. Iacono, AER, Inc.
!-- Adapted for use with NCAR/CAM.
!     May 2004: M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability. 
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Conversion to F90 formatting for consistency with rrtmg_sw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays.
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to add longwave aerosol absorption.
!     Apr 2008: M. J. Iacono, AER, Inc.

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi
      use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
      integer      , intent(in   ) :: imca                          ! flag for mcica [0=off, 1=on]
      integer      , intent(in   ) :: pcols
      integer      , intent(in   ) :: ncol              ! Number of horizontal columns
      integer      , intent(in   ) :: nlay              ! Number of model layers
      integer      , intent(inout) :: icld              ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum
      real(kind=r8), intent(in) :: play(ncol,nlay)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(ncol,nlay+1)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(ncol,nlay)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(ncol,nlay+1)    ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(ncol)             ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr (ncol,nlay)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr (ncol,nlay)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr (ncol,nlay)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr (ncol,nlay)         ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr (ncol,nlay)          ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr (ncol,nlay)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc11vmr (ncol,nlay)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc12vmr (ncol,nlay)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc22vmr (ncol,nlay)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ccl4vmr (ncol,nlay)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: emis (ncol,nbndlw)          ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(ngptlw,ncol,nlay-1)       ! Cloud fraction
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(ngptlw,ncol,nlay-1)       ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl (ngptlw,ncol,nlay-1)       ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: cld(ncol,nlay)  

      real(kind=r8), intent(in) :: cicewp(ncol,nlay-1) 
      real(kind=r8), intent(in) :: cliqwp(ncol,nlay-1)  
      
      real(kind=r8), intent(in) :: tauc_lw(nbndlw,ncol,nlay)        ! cloud optical depth
      
      real(kind=r8), intent(in) :: reicmcl (ncol,nlay-1)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl (ncol,nlay-1)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl (ngptlw,ncol,nlay-1)       ! Cloud optical depth
                                                        !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=r8), intent(in) :: ssacmcl(:,:,:)      ! Cloud single scattering albedo
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   for future expansion
                                                        !   lw scattering not yet available
!      real(kind=r8), intent(in) :: asmcmcl(:,:,:)      ! Cloud asymmetry parameter
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   for future expansion
                                                        !   lw scattering not yet available
      real(kind=r8), intent(in) :: tauaer (ncol,nlay-1,nbndlw)        ! aerosol optical depth
                                                        !   at mid-point of LW spectral bands
                                                        !    Dimensions: (ncol,nlay,nbndlw)
!      real(kind=r8), intent(in) :: ssaaer(:,:,:)       ! aerosol single scattering albedo
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)
!      real(kind=r8), intent(in) :: asmaer(:,:,:)       ! aerosol asymmetry parameter
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)

! ----- Output -----

      real(kind=r8), intent(inout) :: uflx(ncol,nlay+1)           ! Total sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: dflx(ncol,nlay+1)           ! Total sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: hr (ncol,nlay)             ! Total sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(inout) :: uflxc (ncol,nlay+1)          ! Clear sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: dflxc (ncol,nlay+1)          ! Clear sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: hrc(ncol,nlay)           ! Clear sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(inout) :: uflxs (nbndlw,ncol,nlay+1)        ! Total sky longwave upward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)
      real(kind=r8), intent(inout) :: dflxs (nbndlw,ncol,nlay+1)        ! Total sky longwave downward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)

! ----- Local -----
! ----- Output -----

      real(kind=r8) :: uflx_local(ncol,nlay+1)          ! Total sky longwave upward flux (W/m2)
                                         !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: dflx_local(ncol,nlay+1)          ! Total sky longwave downward flux (W/m2)
                                         !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: hr_local (ncol,nlay)                ! Total sky longwave radiative heating rate (K/d)
                                         !    Dimensions: (ncol,nlay)
      real(kind=r8) :: uflxc_local (ncol,nlay+1)           ! Clear sky longwave upward flux (W/m2)
                                         !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: dflxc_local (ncol,nlay+1)           ! Clear sky longwave downward flux (W/m2)
                                         !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: hrc_local(ncol,nlay)              ! Clear sky longwave radiative heating rate (K/d)
                                         !    Dimensions: (ncol,nlay)
      real(kind=r8) :: uflxs_local (nbndlw,ncol,nlay+1)        ! Total sky longwave upward flux spectral (W/m2)
                                         !    Dimensions: (nbndlw,ncol,nlay+1)
      real(kind=r8) :: dflxs_local (nbndlw,ncol,nlay+1)        ! Total sky longwave downward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)


! Control
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: iout                           ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
      integer :: iplon                          ! column loop index
      integer :: ims                            ! value for changing mcica permute seed
      integer :: k                              ! layer loop index
      integer :: ig                             ! g-point loop index
      integer :: i                              ! layer loop index
      integer :: j                              ! layer loop index

! Atmosphere
      real(kind=r8) :: pavel(nlay)              ! layer pressures (mb) 
      real(kind=r8) :: tavel(nlay)              ! layer temperatures (K)
      real(kind=r8) :: pz(0:nlay)               ! level (interface) pressures (hPa, mb)
      real(kind=r8) :: tz(0:nlay)               ! level (interface) temperatures (K)
      real(kind=r8) :: tbound                   ! surface temperature (K)
      real(kind=r8) :: coldry(nlay)             ! dry air column density (mol/cm2)
      real(kind=r8) :: wbrodl(nlay)             ! broadening gas column density (mol/cm2)
      real(kind=r8) :: wkl(mxmol,nlay)          ! molecular amounts (mol/cm-2)
      real(kind=r8) :: wx(maxxsec,nlay)         ! cross-section amounts (mol/cm-2)
      real(kind=r8) :: pwvcm                    ! precipitable water vapor (cm)
      real(kind=r8) :: semiss(nbndlw)           ! lw surface emissivity
      real(kind=r8) :: fracs(nlay,ngptlw)       ! 
      real(kind=r8) :: taug(nlay,ngptlw)        ! gaseous optical depths
      real(kind=r8) :: taut(nlay,ngptlw)        ! gaseous + aerosol optical depths

      real(kind=r8) :: taua(nlay,nbndlw)        ! aerosol optical depth
!      real(kind=r8) :: ssaa(nlay,nbndlw)        ! aerosol single scattering albedo
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)
!      real(kind=r8) :: asma(nlay+1,nbndlw)      ! aerosol asymmetry parameter
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)
      real(kind=r8) :: taucloud(nlay,nbndlw) ! in-cloud optical depth; delta scaled
      real(kind=r8) :: cldfrac(nlay)            ! layer cloud fraction
      real(kind=r8) :: ciwp(nlay)               ! cloud ice water path
      real(kind=r8) :: clwp(nlay)               ! cloud liquid water path
      real(kind=r8) :: rei(nlay)                ! cloud ice particle size
      real(kind=r8) :: rel(nlay)                ! cloud liquid particle size

! Atmosphere - setcoef
      integer :: laytrop                          ! tropopause layer index
      integer :: jp(nlay)                         ! lookup table index 
      integer :: jt(nlay)                         ! lookup table index 
      integer :: jt1(nlay)                        ! lookup table index 
      real(kind=r8) :: planklay(nlay,nbndlw)      ! 
      real(kind=r8) :: planklev(0:nlay,nbndlw)    ! 
      real(kind=r8) :: plankbnd(nbndlw)           ! 

      real(kind=r8) :: colh2o(nlay)               ! column amount (h2o)
      real(kind=r8) :: colco2(nlay)               ! column amount (co2)
      real(kind=r8) :: colo3(nlay)                ! column amount (o3)
      real(kind=r8) :: coln2o(nlay)               ! column amount (n2o)
      real(kind=r8) :: colco(nlay)                ! column amount (co)
      real(kind=r8) :: colch4(nlay)               ! column amount (ch4)
      real(kind=r8) :: colo2(nlay)                ! column amount (o2)
      real(kind=r8) :: colbrd(nlay)               ! column amount (broadening gases)

      integer :: indself(nlay)
      integer :: indfor(nlay)
      real(kind=r8) :: selffac(nlay)
      real(kind=r8) :: selffrac(nlay)
      real(kind=r8) :: forfac(nlay)
      real(kind=r8) :: forfrac(nlay)

      integer :: indminor(nlay)
      real(kind=r8) :: minorfrac(nlay)
      real(kind=r8) :: scaleminor(nlay)
      real(kind=r8) :: scaleminorn2(nlay)

      real(kind=r8) :: &                          !
                         fac00(nlay), fac01(nlay), &
                         fac10(nlay), fac11(nlay) 
      real(kind=r8) :: &                          !
                         rat_h2oco2(nlay),rat_h2oco2_1(nlay), &
                         rat_h2oo3(nlay),rat_h2oo3_1(nlay), &
                         rat_h2on2o(nlay),rat_h2on2o_1(nlay), &
                         rat_h2och4(nlay),rat_h2och4_1(nlay), &
                         rat_n2oco2(nlay),rat_n2oco2_1(nlay), &
                         rat_o3co2(nlay),rat_o3co2_1(nlay)

! Atmosphere/clouds - cldprop
      integer :: ncbands                          ! number of cloud spectral bands
      integer :: inflag                           ! flag for cloud property method
      integer :: iceflag                          ! flag for ice cloud properties
      integer :: liqflag                          ! flag for liquid cloud properties

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=r8) :: cldfmc(ngptlw,nlay)      ! cloud fraction [mcica]
      real(kind=r8) :: ciwpmc(ngptlw,nlay)      ! cloud ice water path [mcica]
      real(kind=r8) :: clwpmc(ngptlw,nlay)      ! cloud liquid water path [mcica]
      real(kind=r8) :: relqmc(nlay)             ! liquid particle size (microns)
      real(kind=r8) :: reicmc(nlay)             ! ice particle effective radius (microns)
      real(kind=r8) :: dgesmc(nlay)             ! ice particle generalized effective size (microns)
      real(kind=r8) :: taucmc(ngptlw,nlay)      ! cloud optical depth [mcica]
!      real(kind=r8) :: ssacmc(ngptlw,nlay)     ! cloud single scattering albedo [mcica]
                                                !   for future expansion 
                                                !   (lw scattering not yet available)
!      real(kind=r8) :: asmcmc(ngptlw,nlay)     ! cloud asymmetry parameter [mcica]
                                                !   for future expansion 
                                                !   (lw scattering not yet available)

! Output
      real(kind=r8) :: totuflux(0:nlay)         ! upward longwave flux (w/m2)
      real(kind=r8) :: totdflux(0:nlay)         ! downward longwave flux (w/m2)
      real(kind=r8) :: totufluxs(nbndlw,0:nlay) ! upward longwave flux spectral (w/m2)
      real(kind=r8) :: totdfluxs(nbndlw,0:nlay) ! downward longwave flux spectral (w/m2)
      real(kind=r8) :: fnet(0:nlay)             ! net longwave flux (w/m2)
      real(kind=r8) :: htr(0:nlay)              ! longwave heating rate (k/day)
      real(kind=r8) :: totuclfl(0:nlay)         ! clear sky upward longwave flux (w/m2)
      real(kind=r8) :: totdclfl(0:nlay)         ! clear sky downward longwave flux (w/m2)
      real(kind=r8) :: fnetc(0:nlay)            ! clear sky net longwave flux (w/m2)
      real(kind=r8) :: htrc(0:nlay)             ! clear sky longwave heating rate (k/day)

! Initializations


      uflx_local=0.0_r8 ;   dflx_local=0.0_r8;    hr_local=0.0_r8;   uflxc_local=0.0_r8 ;  
      dflxc_local=0.0_r8;   hrc_local=0.0_r8;     uflxs_local=0.0_r8;dflxs_local=0.0_r8;
      
! Atmosphere
      pavel = 0.0_r8;     tavel = 0.0_r8;     pz = 0.0_r8;        tz = 0.0_r8;        
      tbound= 0.0_r8;     coldry = 0.0_r8;    wbrodl = 0.0_r8;    wkl = 0.0_r8;       
      wx = 0.0_r8;        pwvcm= 0.0_r8;      semiss = 0.0_r8;    
      fracs = 0.0_r8 ;    taug = 0.0_r8  ;    taut = 0.0_r8  ;    taua= 0.0_r8;       

! Atmosphere - setcoef
      laytrop=0;  ! tropopause layer index
      jp = 0;jt = 0;jt1 = 0;  planklay = 0.0_r8;  planklev = 0.0_r8;  plankbnd = 0.0_r8;  
       colh2o = 0.0_r8;    colco2 = 0.0_r8;    colo3 = 0.0_r8;     coln2o = 0.0_r8;    colco = 0.0_r8;     
       colch4 = 0.0_r8;    colo2 = 0.0_r8;     colbrd = 0.0_r8;    
       indself = 0;indfor = 0;   selffac = 0.0_r8;    selffrac = 0.0_r8;   forfac = 0.0_r8;     forfrac = 0.0_r8;    
      indminor = 0;        minorfrac = 0.0_r8;           
      scaleminor = 0.0_r8;          scaleminorn2 = 0.0_r8;        fac00 = 0.0_r8; fac01 = 0.0_r8;
      fac10 = 0.0_r8; fac11 = 0.0_r8;rat_h2oco2 = 0.0_r8;  rat_h2oco2_1 = 0.0_r8;
      rat_h2oo3 = 0.0_r8;   rat_h2oo3_1 = 0.0_r8;rat_h2on2o = 0.0_r8;  rat_h2on2o_1 = 0.0_r8;
      rat_h2och4 = 0.0_r8;  rat_h2och4_1 = 0.0_r8;rat_n2oco2 = 0.0_r8;  rat_n2oco2_1 = 0.0_r8;
      rat_o3co2 = 0.0_r8;   rat_o3co2_1 = 0.0_r8;
! Atmosphere/clouds - cldprop
      ncbands=0;inflag=0;iceflag=0;liqflag=0;
! Atmosphere/clouds - cldprmc [mcica]
      cldfmc = 0.0_r8; ciwpmc = 0.0_r8; clwpmc = 0.0_r8; relqmc = 0.0_r8; 
      reicmc = 0.0_r8; dgesmc = 0.0_r8; taucmc = 0.0_r8; 
      taucloud=0.0_r8; cldfrac= 0.0_r8; ciwp = 0.0_r8; ;clwp= 0.0_r8;
      rei= 0.0_r8; rel= 0.0_r8;
! Output
      totuflux = 0.0_r8;   totdflux = 0.0_r8;   totufluxs = 0.0_r8;  
      totdfluxs = 0.0_r8;  fnet = 0.0_r8;       htr = 0.0_r8;        
      totuclfl = 0.0_r8;   totdclfl = 0.0_r8;   fnetc = 0.0_r8;      htrc = 0.0_r8;       

! Atmosphere/clouds - cldprop
      
      
      oneminus = 1._r8 - 1.e-6_r8
      pi = 2._r8 * asin(1._r8)
      fluxfac = pi * 2.e4_r8                    ! orig:   fluxfac = pi * 2.d4  
      istart = 1
      iend = 16
      iout = 0
      ims = 1

! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud overlap method  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap
! icld = 2, with clouds using maximum/random cloud overlap
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 10, input total aerosol optical depth (tauaer) directly 
      iaer = 10

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 256 to 140 for input absorption coefficient 
! data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_lw_ini

!  This is the main longitude/column loop within RRTMG.
      do iplon = 1, ncol

!  Prepare atmospheric profile from GCM for use in RRTMG, and define
!  other input parameters.  

         call inatm ( &
              iplon                                   , &
              ncol                                    , &
              nlay                                    , &
              icld                                    , &
              iaer                                    , &
              play      (1:ncol,1:nlay)               , &
              plev      (1:ncol,1:nlay+1)             , &
              tlay      (1:ncol,1:nlay)               , &
              tlev      (1:ncol,1:nlay+1)             , &
              tsfc      (1:ncol)                      , &
              h2ovmr    (1:ncol,1:nlay)               , &
              o3vmr     (1:ncol,1:nlay)               , &
              co2vmr    (1:ncol,1:nlay)               , &
              ch4vmr    (1:ncol,1:nlay)               , &
              o2vmr     (1:ncol,1:nlay)               , &
              n2ovmr    (1:ncol,1:nlay)               , &
              cfc11vmr  (1:ncol,1:nlay)               , &
              cfc12vmr  (1:ncol,1:nlay)               , &
              cfc22vmr  (1:ncol,1:nlay)               , &
              ccl4vmr   (1:ncol,1:nlay)               , &
              emis      (1:ncol,1:nbndlw)             , &
              inflglw                                 , &
              iceflglw                                , &
              liqflglw                                , &
              cldfmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              taucmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              ciwpmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              clwpmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              reicmcl   (1:ncol  ,1:nlay-1)           , &
              relqmcl   (1:ncol  ,1:nlay-1)           , &
              tauaer    (1:ncol  ,1:nlay-1 ,1:nbndlw) , &
              pavel     (1:nlay)                      , &
              pz        (0:nlay)                      , &
              tavel     (1:nlay)                      , &
              tz        (0:nlay)                      , &
              tbound                                  , &
              semiss    (1:nbndlw)                    , &
              coldry    (1:nlay)                      , &
              wkl       (1:mxmol,1:nlay)              , &
              wbrodl    (1:nlay)                      , &
              wx        (1:maxxsec,1:nlay)            , &
              pwvcm                                   , &
              inflag                                  , &
              iceflag                                 , &
              liqflag                                 , &
              cldfmc    (1:ngptlw,1:nlay)             , &
              taucmc    (1:ngptlw,1:nlay)             , &
              ciwpmc    (1:ngptlw,1:nlay)             , &
              clwpmc    (1:ngptlw,1:nlay)             , &
              reicmc    (1:nlay)                      , &
              dgesmc    (1:nlay)                      , &
              relqmc    (1:nlay)                      , &
              taua      (1:nlay,1:nbndlw)                )

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprop.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprop.  

!  If McICA is requested use cloud fraction and cloud physical properties 
!  generated by sub-column cloud generator above. 
         cldfrac(1:nlay)     =  cld(iplon,1:nlay) 

         ciwp   (1:nlay-1) = cicewp   (iplon,1:nlay-1)
         ciwp   (nlay)     = cicewp   (iplon,nlay-1)

         clwp   (1:nlay-1) = cliqwp   (iplon,1:nlay-1)
         clwp   (nlay)     = cliqwp   (iplon,nlay-1)

         rei    (1:nlay-1) = reicmcl  (iplon,1:nlay-1) 
         rei(nlay)         = reicmcl  (iplon,nlay-1)
 
         rel(1:nlay-1)     = relqmcl  (iplon,1:nlay-1)
         rel(nlay)          =relqmcl  (iplon,nlay-1)

         if (imca.eq.0) then
             do i = 1, nlay
               cldfrac(i)=MIN(MAX(cldfrac(i),0.00001_r8),0.9999_r8)
             enddo
            call cldprop( &
                         nlay                              , &
                         inflag                            , &
                         iceflag                           , &
                         liqflag                           , &
                         cldfrac  (1:nlay)                 , &
                         tauc_lw  (1:nbndlw,iplon,1:nlay)  , &
                         ciwp     (1:nlay)                 , &
                         clwp     (1:nlay)                 , &
                         rei      (1:nlay)                 , &
                         2*rei    (1:nlay)                 , &
                         rel      (1:nlay)                 , &
                         ncbands                           , &
                         taucloud (1:nlay,1:nbndlw)          )
         else
            call cldprmc( &
                         nlay                      , &
                         inflag                    , &
                         iceflag                   , &
                         liqflag                   , &
                         cldfmc   (1:ngptlw,1:nlay), &
                         ciwpmc   (1:ngptlw,1:nlay), &
                         clwpmc   (1:ngptlw,1:nlay), &
                         reicmc   (1:nlay)         , &
                         dgesmc   (1:nlay)         , &
                         relqmc   (1:nlay)         , &
                         ncbands                   , &
                         taucmc   (1:ngptlw,1:nlay)  )
         endif

! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the 
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres. 

         call setcoef( &
                      nlay                         , &
                      istart                       , &
                      pavel       (1:nlay)         , &
                      tavel       (1:nlay)         , &
                      tz          (0:nlay)         , &
                      tbound                       , &
                      semiss      (1:nbndlw)       , &
                      coldry      (1:nlay)         , &
                      wkl         (1:mxmol,1:nlay) , &
                      wbrodl      (1:nlay)         , &
                      laytrop                      , &
                      jp          (1:nlay)         , &
                      jt          (1:nlay)         , &
                      jt1         (1:nlay)         , &
                      planklay    (1:nlay,1:nbndlw), &
                      planklev    (0:nlay,1:nbndlw), &
                      plankbnd    (1:nbndlw)       , &
                      colh2o      (1:nlay)         , &
                      colco2      (1:nlay)         , &
                      colo3       (1:nlay)         , &
                      coln2o      (1:nlay)         , &
                      colco       (1:nlay)         , &
                      colch4      (1:nlay)         , &
                      colo2       (1:nlay)         , &
                      colbrd      (1:nlay)         , &
                      fac00       (1:nlay)         , &
                      fac01       (1:nlay)         , &
                      fac10       (1:nlay)         , &
                      fac11       (1:nlay)         , &
                      rat_h2oco2  (1:nlay)         , &
                      rat_h2oco2_1(1:nlay)         , &
                      rat_h2oo3   (1:nlay)         , &
                      rat_h2oo3_1 (1:nlay)         , &
                      rat_h2on2o  (1:nlay)         , &
                      rat_h2on2o_1(1:nlay)         , &
                      rat_h2och4  (1:nlay)         , &
                      rat_h2och4_1(1:nlay)         , &
                      rat_n2oco2  (1:nlay)         , &
                      rat_n2oco2_1(1:nlay)         , &
                      rat_o3co2   (1:nlay)         , &
                      rat_o3co2_1 (1:nlay)         , &
                      selffac     (1:nlay)         , &
                      selffrac    (1:nlay)         , &
                      indself     (1:nlay)         , &
                      forfac      (1:nlay)         , &
                      forfrac     (1:nlay)         , &
                      indfor      (1:nlay)         , &
                      minorfrac   (1:nlay)         , &
                      scaleminor  (1:nlay)         , &
                      scaleminorn2(1:nlay)         , &
                      indminor    (1:nlay)           )

!  Calculate the gaseous optical depths and Planck fractions for 
!  each longwave spectral band.

         call taumol( &
                     nlay                           , &
                     pavel        (1:nlay)          , &
                     wx           (1:maxxsec,1:nlay), &
                     coldry       (1:nlay)          , &
                     laytrop                        , &
                     jp           (1:nlay)          , &
                     jt           (1:nlay)          , &
                     jt1          (1:nlay)          , &
                     planklay     (1:nlay,1:nbndlw) , &
                     planklev     (0:nlay,1:nbndlw) , &
                     plankbnd     (1:nbndlw)        , &
                     colh2o       (1:nlay)          , &
                     colco2       (1:nlay)          , &
                     colo3        (1:nlay)          , &
                     coln2o       (1:nlay)          , &
                     colco        (1:nlay)          , &
                     colch4       (1:nlay)          , &
                     colo2        (1:nlay)          , &
                     colbrd       (1:nlay)          , &
                     fac00        (1:nlay)          , &
                     fac01        (1:nlay)          , &
                     fac10        (1:nlay)          , &
                     fac11        (1:nlay)          , &
                     rat_h2oco2   (1:nlay)          , &
                     rat_h2oco2_1 (1:nlay)          , &
                     rat_h2oo3    (1:nlay)          , &
                     rat_h2oo3_1  (1:nlay)          , &
                     rat_h2on2o   (1:nlay)          , &
                     rat_h2on2o_1 (1:nlay)          , &
                     rat_h2och4   (1:nlay)          , &
                     rat_h2och4_1 (1:nlay)          , &
                     rat_n2oco2   (1:nlay)          , &
                     rat_n2oco2_1 (1:nlay)          , &
                     rat_o3co2    (1:nlay)          , &
                     rat_o3co2_1  (1:nlay)          , &
                     selffac      (1:nlay)          , &
                     selffrac     (1:nlay)          , &
                     indself      (1:nlay)          , &
                     forfac       (1:nlay)          , &
                     forfrac      (1:nlay)          , &
                     indfor       (1:nlay)          , &
                     minorfrac    (1:nlay)          , &
                     scaleminor   (1:nlay)          , &
                     scaleminorn2 (1:nlay)          , &
                     indminor     (1:nlay)          , &
                     fracs        (1:nlay,1:ngptlw) , &
                     taug         (1:nlay,1:ngptlw)  )



! Combine gaseous and aerosol optical depths, if aerosol active
         if (iaer .eq. 0) then
            do k = 1, nlay
               do ig = 1, ngptlw 
                  taut(k,ig) = taug(k,ig)
               enddo
            enddo
         elseif (iaer .eq. 10) then
            do k = 1, nlay
               do ig = 1, ngptlw 
                  taut(k,ig) = taug(k,ig) + taua(k,ngb(ig))
               enddo
            enddo
         endif

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.
! For McICA, RTRNMC is called for clear and cloudy calculations.
        if (imca .eq. 0) then
           if (icld .eq. 1) then
              call rtrn( &
                     nlay                          , &
                     istart                        , &
                     iend                          , &
                     iout                          , &
                     pz        (0:nlay)            , &
                     semiss    (1:nbndlw)          , &
                     ncbands                       , &
                     cldfrac   (1:nlay)            , &
                     taucloud  (1:nlay ,1:nbndlw)  , &
                     planklay  (1:nlay ,1:nbndlw)  , &
                     planklev  (0:nlay ,1:nbndlw)  , &
                     plankbnd  (1:nbndlw)          , &
                     pwvcm                         , &
                     fracs     (1:nlay,1:ngptlw)   , &
                     taut      (1:nlay,1:ngptlw)   , &
                     totuflux  (0:nlay)            , &
                     totdflux  (0:nlay)            , &
                     fnet      (0:nlay)            , &
                     htr       (0:nlay)            , &
                     totuclfl  (0:nlay)            , &
                     totdclfl  (0:nlay)            , &
                     fnetc     (0:nlay)            , &
                     htrc      (0:nlay)            , &
                     totufluxs (1:nbndlw,0:nlay)   , &
                     totdfluxs (1:nbndlw,0:nlay)     )
           else
              call rtrnmr( &
                     nlay                          , &
                     istart                        , &
                     iend                          , &
                     iout                          , &
                     pz         (0:nlay)           , &
                     semiss     (1:nbndlw)         , &
                     ncbands                       , &
                     cldfrac    (1:nlay)           , &
                     taucloud   (1:nlay ,1:nbndlw) , &
                     planklay   (1:nlay ,1:nbndlw) , &
                     planklev   (0:nlay ,1:nbndlw) , &
                     plankbnd   (1:nbndlw)         , &
                     pwvcm                         , &
                     fracs      (1:nlay,1:ngptlw)  , &
                     taut       (1:nlay,1:ngptlw)  , &
                     totuflux   (0:nlay)           , &
                     totdflux   (0:nlay)           , &
                     fnet       (0:nlay)           , &
                     htr        (0:nlay)           , &
                     totuclfl   (0:nlay)           , &
                     totdclfl   (0:nlay)           , &
                     fnetc      (0:nlay)           , &
                     htrc       (0:nlay)           , &
                     totufluxs  (1:nbndlw,0:nlay)  , &
                     totdfluxs  (1:nbndlw,0:nlay)     )

 
           endif
        elseif (imca .eq. 1) then

              call rtrnmc( &
                     nlay                       , &
                     istart                     , &
                     iend                       , &
                     iout                       , &
                     pz        (0:nlay)         , &
                     semiss    (1:nbndlw)       , &
                     ncbands                    , &
                     cldfmc    (1:ngptlw,1:nlay), &
                     taucmc    (1:ngptlw,1:nlay), &
                     planklay  (1:nlay,1:nbndlw), &
                     planklev  (0:nlay,1:nbndlw), &
                     plankbnd  (1:nbndlw)       , &
                     pwvcm                      , &
                     fracs     (1:nlay,1:ngptlw), &
                     taut      (1:nlay,1:ngptlw), &
                     totuflux  (0:nlay)         , &
                     totdflux  (0:nlay)         , &
                     fnet      (0:nlay)         , &
                     htr       (0:nlay)         , &
                     totuclfl  (0:nlay)         , &
                     totdclfl  (0:nlay)         , &
                     fnetc     (0:nlay)         , &
                     htrc      (0:nlay)         , &
                     totufluxs (1:nbndlw,0:nlay), &
                     totdfluxs (1:nbndlw,0:nlay)  )
        endif

!  Transfer up and down fluxes and heating rate to output arrays.
!  Vertical indexing goes from bottom to top

         do k = 0, nlay
            uflx_local(iplon,k+1) = totuflux(k)
            dflx_local(iplon,k+1) = totdflux(k)
            uflxc_local(iplon,k+1) = totuclfl(k)
            dflxc_local(iplon,k+1) = totdclfl(k)
            uflxs_local(:,iplon,k+1) = totufluxs(:,k)
            dflxs_local(:,iplon,k+1) = totdfluxs(:,k)
         enddo
         do k = 0, nlay-1
            hr_local(iplon,k+1) = htr(k)
            hrc_local(iplon,k+1) = htrc(k)
         enddo
!***********************************************************************************************************************

! ----- Output -----
         do k = 0, nlay
            uflx (iplon,k+1) = uflx (iplon,k+1) + uflx_local (iplon,k+1)! Total sky longwave upward flux (W/m2)  !    Dimensions: (ncol,nlay+1)
            dflx (iplon,k+1) = dflx (iplon,k+1) + dflx_local (iplon,k+1)! Total sky longwave downward flux (W/m2)  !    Dimensions: (ncol,nlay+1)
            uflxc(iplon,k+1) = uflxc(iplon,k+1) + uflxc_local(iplon,k+1)! Clear sky longwave upward flux (W/m2)  !    Dimensions: (ncol,nlay+1)
            dflxc(iplon,k+1) = dflxc(iplon,k+1) + dflxc_local(iplon,k+1)! Clear sky longwave downward flux (W/m2)  !    Dimensions: (ncol,nlay+1)
         enddo
         do k = 0, nlay
            DO j=1,nbndlw
               uflxs(j,iplon,k+1) = uflxs(j,iplon,k+1) + uflxs_local(j,iplon,k+1) ! Total sky longwave upward flux spectral (W/m2)  !    Dimensions: (nbndlw,ncol,nlay+1)
               dflxs(j,iplon,k+1) = dflxs(j,iplon,k+1) + dflxs_local(j,iplon,k+1) ! Total sky longwave downward flux spectral (W/m2)!    Dimensions: (nbndlw,ncol,nlay+1)
            END DO
         END DO    
         do k = 0, nlay-1
            hr   (iplon,k+1) = hr   (iplon,k+1) + hr_local   (iplon,k+1)  ! Total sky longwave radiative heating rate (K/d) ! Dimensions: (ncol,nlay)
            hrc  (iplon,k+1) = hrc  (iplon,k+1) + hrc_local  (iplon,k+1)  ! Clear sky longwave radiative heating rate (K/d) ! Dimensions: (ncol,nlay)
         END DO    

! End longitude loop

      enddo

      end subroutine rrtmg_lw

!***************************************************************************
      subroutine inatm (  &!                 call inatm ( &
              iplon     , &!                      iplon                                   , &
              ncol      , &!                      ncol                                    , &
              nlay      , &!                      nlay                                    , &
              icld      , &!                      icld                                    , &
              iaer      , &!                      iaer                                    , &
              play      , &!                      play      (1:ncol,1:nlay)               , &
              plev      , &!                      plev      (1:ncol,1:nlay+1)             , &
              tlay      , &!                      tlay      (1:ncol,1:nlay)               , &
              tlev      , &!                      tlev      (1:ncol,1:nlay+1)             , &
              tsfc      , &!                      tsfc      (1:ncol)                      , &
              h2ovmr    , &!                      h2ovmr    (1:ncol,1:nlay)               , &
              o3vmr     , &!                      o3vmr     (1:ncol,1:nlay)               , &
              co2vmr    , &!                      co2vmr    (1:ncol,1:nlay)               , &
              ch4vmr    , &!                      ch4vmr    (1:ncol,1:nlay)               , &
              o2vmr     , &!                      o2vmr     (1:ncol,1:nlay)               , &
              n2ovmr    , &!                      n2ovmr    (1:ncol,1:nlay)               , &
              cfc11vmr  , &!                      cfc11vmr  (1:ncol,1:nlay)               , &
              cfc12vmr  , &!                      cfc12vmr  (1:ncol,1:nlay)               , &
              cfc22vmr  , &!                      cfc22vmr  (1:ncol,1:nlay)               , &
              ccl4vmr   , &!                      ccl4vmr   (1:ncol,1:nlay)               , &
              emis      , &!                      emis      (1:ncol,1:nbndlw)             , &
              inflglw   , &!                      inflglw                                 , &
              iceflglw  , &!                      iceflglw                                , &
              liqflglw  , &!                      liqflglw                                , &
              cldfmcl   , &!                      cldfmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              taucmcl   , &!                      taucmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              ciwpmcl   , &!                      ciwpmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              clwpmcl   , &!                      clwpmcl   (1:ngptlw,1:ncol,1:nlay-1)    , &
              reicmcl   , &!                      reicmcl   (1:ncol  ,1:nlay-1)           , &
              relqmcl   , &!                      relqmcl   (1:ncol  ,1:nlay-1)           , &
              tauaer    , &!                      tauaer    (1:ncol  ,1:nlay-1 ,1:nbndlw) , &
              pavel     , &!                      pavel     (1:nlay)                      , &
              pz        , &!                      pz        (0:nlay)                      , &
              tavel     , &!                      tavel     (1:nlay)                      , &
              tz        , &!                      tz        (0:nlay)                      , &
              tbound    , &!                      tbound                                  , &
              semiss    , &!                      semiss    (1:nbndlw)                    , &
              coldry    , &!                      coldry    (1:nlay)                      , &
              wkl       , &!                      wkl       (1:mxmol,1:nlay)              , &
              wbrodl    , &!                      wbrodl    (1:nlay)                      , &
              wx        , &!                      wx        (1:maxxsec,1:nlay)            , &
              pwvcm     , &!                      pwvcm                                   , &
              inflag    , &!                      inflag                                  , &
              iceflag   , &!                      iceflag                                 , &
              liqflag   , &!                      liqflag                                 , &
              cldfmc    , &!                      cldfmc    (1:ngptlw,1:nlay)             , &
              taucmc    , &!                      taucmc    (1:ngptlw,1:nlay)             , &
              ciwpmc    , &!                      ciwpmc    (1:ngptlw,1:nlay)             , &
              clwpmc    , &!                      clwpmc    (1:ngptlw,1:nlay)             , &
              reicmc    , &!                      reicmc    (1:nlay)                      , &
              dgesmc    , &!                      dgesmc    (1:nlay)                      , &
              relqmc    , &!                      relqmc    (1:nlay)                      , &
              taua        )!                      taua      (1:nlay,1:nbndlw)                )
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_LW.
!  Set other RRTMG_LW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, nmol, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrlw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave, ixindx

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: ncol
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(in) :: icld                       ! clear/cloud and cloud overlap flag
      integer, intent(in) :: iaer                       ! aerosol option flag

      real(kind=r8), intent(in) :: play(1:ncol,1:nlay)             ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(1:ncol,1:nlay+1)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(1:ncol,1:nlay)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(1:ncol,1:nlay+1)            ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(1:ncol)              ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(1:ncol,1:nlay)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(1:ncol,1:nlay)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(1:ncol,1:nlay)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(1:ncol,1:nlay)          ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(1:ncol,1:nlay)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(1:ncol,1:nlay)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc11vmr(1:ncol,1:nlay)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc12vmr(1:ncol,1:nlay)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc22vmr(1:ncol,1:nlay)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ccl4vmr(1:ncol,1:nlay)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: emis(1:ncol,1:nbndlw)             ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(1:ngptlw,1:ncol,1:nlay-1)       ! Cloud fraction
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(1:ngptlw,1:ncol,1:nlay-1)       ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl(1:ngptlw,1:ncol,1:nlay-1)       ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: reicmcl(1:ncol  ,1:nlay-1)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl(1:ncol  ,1:nlay-1)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl(1:ngptlw,1:ncol,1:nlay-1)       ! Cloud optical depth
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: tauaer (1:ncol  ,1:nlay-1 ,1:nbndlw)        ! Aerosol optical depth
                                                        !    Dimensions: (ncol,nlay,nbndlw)

! ----- Output -----
! Atmosphere
      real(kind=r8), intent(out) :: pavel(1:nlay)            ! layer pressures (mb) 
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: tavel(1:nlay)             ! layer temperatures (K)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: pz(0:nlay)              ! level (interface) pressures (hPa, mb)
                                                        !    Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tz(0:nlay)              ! level (interface) temperatures (K)
                                                        !    Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tbound              ! surface temperature (K)
      real(kind=r8), intent(out) :: coldry(1:nlay)           ! dry air column density (mol/cm2)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: wbrodl(1:nlay)            ! broadening gas column density (mol/cm2)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: wkl(1:mxmol,1:nlay)             ! molecular amounts (mol/cm-2)
                                                        !    Dimensions: (mxmol,nlay)
      real(kind=r8), intent(out) :: wx(1:maxxsec,1:nlay)             ! cross-section amounts (mol/cm-2)
                                                        !    Dimensions: (maxxsec,nlay)
      real(kind=r8), intent(out) :: pwvcm               ! precipitable water vapor (cm)
      real(kind=r8), intent(out) :: semiss(1:nbndlw)           ! lw surface emissivity
                                                        !    Dimensions: (nbndlw)

! Atmosphere/clouds - cldprop
      integer, intent(out) :: inflag                    ! flag for cloud property method
      integer, intent(out) :: iceflag                   ! flag for ice cloud properties
      integer, intent(out) :: liqflag                   ! flag for liquid cloud properties

      real(kind=r8), intent(out) :: cldfmc(1:ngptlw,1:nlay)          ! cloud fraction [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: ciwpmc(1:ngptlw,1:nlay)          ! cloud ice water path [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: clwpmc(1:ngptlw,1:nlay)          ! cloud liquid water path [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: relqmc(1:nlay)           ! liquid particle effective radius (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: reicmc(1:nlay)           ! ice particle effective radius (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: dgesmc(1:nlay)           ! ice particle generalized effective size (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: taucmc(1:ngptlw,1:nlay)          ! cloud optical depth [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: taua(1:nlay,1:nbndlw)           ! Aerosol optical depth
                                                        ! Dimensions: (nlay,nbndlw)


! ----- Local -----
      real(kind=r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
      real(kind=r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)
!      real(kind=r8), parameter :: amc = 44.0098_r8      ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=r8), parameter :: amo = 47.9998_r8      ! Molecular weight of ozone (g/mol)
!      real(kind=r8), parameter :: amo2 = 31.9999_r8     ! Molecular weight of oxygen (g/mol)
!      real(kind=r8), parameter :: amch4 = 16.0430_r8    ! Molecular weight of methane (g/mol)
!      real(kind=r8), parameter :: amn2o = 44.0128_r8    ! Molecular weight of nitrous oxide (g/mol)
!      real(kind=r8), parameter :: amc11 = 137.3684_r8   ! Molecular weight of CFC11 (g/mol) - CCL3F
!      real(kind=r8), parameter :: amc12 = 120.9138_r8   ! Molecular weight of CFC12 (g/mol) - CCL2F2
!      real(kind=r8), parameter :: amc22 = 86.4688_r8    ! Molecular weight of CFC22 (g/mol) - CHCLF2
!      real(kind=r8), parameter :: amcl4 = 153.823_r8    ! Molecular weight of CCL4 (g/mol) - CCL4

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
      real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
      real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
      real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
      real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
      real(kind=r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
      real(kind=r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

      real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)

      integer :: isp, l, ix, n, imol, ib, ig            ! Loop indices
      real(kind=r8) :: amm, amttl, wvttl, wvsh, summol  

!  Initialize all molecular amounts and cloud properties to zero here, then pass input amounts
!  into RRTM arrays below.

      wkl(:,:) = 0.0_r8
      wx(:,:) = 0.0_r8
      cldfmc(:,:) = 0.0_r8
      taucmc(:,:) = 0.0_r8
      ciwpmc(:,:) = 0.0_r8
      clwpmc(:,:) = 0.0_r8
      reicmc(:) = 0.0_r8
      dgesmc(:) = 0.0_r8
      relqmc(:) = 0.0_r8
      taua(:,:) = 0.0_r8
      amttl = 0.0_r8
      wvttl = 0.0_r8
 
!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_LW arrays for pressure, temperature,
!  and molecular amounts.  
!  Pressures are input in mb, or are converted to mb here.
!  Molecular amounts are input in volume mixing ratio, or are converted from 
!  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
!  here. These are then converted to molecular amount (molec/cm2) below.  
!  The dry air column COLDRY (in molec/cm2) is calculated from the level 
!  pressures, pz (in mb), based on the hydrostatic equation and includes a 
!  correction to account for h2o in the layer.  The molecular weight of moist 
!  air (amm) is calculated for each layer.  
!  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
!  assumes GCM input fields are also bottom to top. Input layer indexing
!  from GCM fields should be reversed here if necessary.

      pz(0) = plev(iplon,nlay+1)
      tz(0) = tlev(iplon,nlay+1)
      do l = 1, nlay
         pavel(l) = play(iplon,nlay-l+1)
         tavel(l) = tlay(iplon,nlay-l+1)
         pz(l) = plev(iplon,nlay-l+1)
         tz(l) = tlev(iplon,nlay-l+1)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,nlay-l+1)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,nlay-l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,nlay-l)/(1._r8 - h2o(iplon,nlay-l)))*amdw
         wkl(2,l) = co2vmr(iplon,nlay-l+1)
         wkl(3,l) = o3vmr(iplon,nlay-l+1)
         wkl(4,l) = n2ovmr(iplon,nlay-l+1)
         wkl(6,l) = ch4vmr(iplon,nlay-l+1)
         wkl(7,l) = o2vmr(iplon,nlay-l+1)

         amm = (1._r8 - wkl(1,l)) * amd + wkl(1,l) * amw            

         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_r8 * avogad / &
                     (1.e2_r8 * grav * amm * (1._r8 + wkl(1,l)))

! Set cross section molecule amounts from input; convert to vmr if necessary
         wx(1,l) = ccl4vmr(iplon,nlay-l+1)
         wx(2,l) = cfc11vmr(iplon,nlay-l+1)
         wx(3,l) = cfc12vmr(iplon,nlay-l+1)
         wx(4,l) = cfc22vmr(iplon,nlay-l+1)

      enddo

      coldry(nlay) = (pz(nlay-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(1,nlay-1)))

! At this point all molecular amounts in wkl and wx are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  also, compute precipitable
! water vapor for diffusivity angle adjustments in rtrn and rtrnmr.

      do l = 1, nlay
         summol = 0.0_r8
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         wbrodl(l) = coldry(l) * (1._r8 - summol)
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
         amttl = amttl + coldry(l)+wkl(1,l)
         wvttl = wvttl + wkl(1,l)
         do ix = 1,maxxsec
            if (ixindx(ix) .ne. 0) then
               wx(ixindx(ix),l) = coldry(l) * wx(ix,l) * 1.e-20_r8
            endif
         enddo
      enddo

      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm = wvsh * (1.e3_r8 * pz(0)) / (1.e2_r8 * grav)

! Set spectral surface emissivity for each longwave band.  

      do n=1,nbndlw
         semiss(n) = emis(iplon,n)
!          semiss(n) = 1.0_r8
      enddo

! Transfer aerosol optical properties to RRTM variable;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlay-1
            do ib = 1, nbndlw
               taua(l,ib) = tauaer(iplon,nlay-l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables,
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflglw
         iceflag = iceflglw
         liqflag = liqflglw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

         do l = 1, nlay-1
            do ig = 1, ngptlw
               cldfmc(ig,l) = cldfmcl(ig,iplon,nlay-l)
               taucmc(ig,l) = taucmcl(ig,iplon,nlay-l)
               ciwpmc(ig,l) = ciwpmcl(ig,iplon,nlay-l)
               clwpmc(ig,l) = clwpmcl(ig,iplon,nlay-l)
            enddo
            reicmc(l) = reicmcl(iplon,nlay-l)
            if (iceflag .eq. 3) then
               dgesmc(l) = 1.5396_r8 * reicmcl(iplon,nlay-l)
            endif
            relqmc(l) = relqmcl(iplon,nlay-l)
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

         cldfmc(:,nlay) = 0.0_r8
         taucmc(:,nlay) = 0.0_r8
         ciwpmc(:,nlay) = 0.0_r8
         clwpmc(:,nlay) = 0.0_r8
         reicmc(nlay) = 0.0_r8
         dgesmc(nlay) = 0.0_r8
         relqmc(nlay) = 0.0_r8
         taua(nlay,:) = 0.0_r8

      endif
      
      end subroutine inatm

      end module rrtmg_lw_rad


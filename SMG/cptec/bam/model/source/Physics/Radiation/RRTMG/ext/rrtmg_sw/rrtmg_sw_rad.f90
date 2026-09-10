!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/01/03 21:35:35 $
!
!  Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) (coarse mode) included
!  Fine aerosol mode is included
!  Modifications are marked by
!  !tar begin.....!tar end


       module rrtmg_sw_rad

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
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! --------- Modules ---------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb
      use rrsw_vsn
      use mcica_subcol_gen_sw, only: mcica_subcol_sw
      use rrtmg_sw_cldprop, only: cldprop_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
! Move call to rrtmg_sw_ini and following use association to 
! GCM initialization area
!      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvmc, only: spcvmc_sw
      use rrtmg_sw_spcvrt, only: spcvrt_sw

!      use perf_mod

      implicit none

! public interfaces/functions/subroutines
!      public :: rrtmg_sw, inatm_sw, earth_sun
      public :: rrtmg_sw

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------
!      subroutine rrtmg_sw &
!            ( ncol    ,nlay    ,icld    ,          &
!             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
!             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  , &
!             asdir   ,asdif   ,aldir   ,aldif   , &
!             coszen  ,adjes   ,dyofyr  ,solvar, &
!             inflgsw ,iceflgsw,liqflgsw, &
!             cldfmcl ,taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl, &
!             ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
!             tauaer  ,ssaaer  ,asmaer  , &
!             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
!             dirdnuv, dirdnir, difdnuv, difdnir, &
!             dirdnuvc,difdnuvc,dirdnirc,difdnirc,ninflx, ninflxc, &
!             swuflxs, swdflxs)

      subroutine rrtmg_sw( &
             imca    , & !                ! flag for mcica [0=off, 1=on]
             pcols   , & !integer, intent(in) :: pcols                      ! Total Number of horizontal columns 
             ncol    , & !integer, intent(in) :: ncol                       ! Number of horizontal columns     
             nlay    , & !integer, intent(in) :: nlay                       ! Number of model layers
             icld    , & !integer, intent(inout) :: icld                    ! Cloud overlap method
             play    , & !real(kind=r8), intent(in) :: play(pcols,nlay)            ! Layer pressures (hPa, mb)
             plev    , & !real(kind=r8), intent(in) :: plev (pcols,nlay+1)            ! Interface pressures (hPa, mb)
             tlay    , & !real(kind=r8), intent(in) :: tlay(pcols,nlay)            ! Layer temperatures (K)
             tlev    , & !real(kind=r8), intent(in) :: tlev (pcols,nlay+1)            ! Interface temperatures (K)
             tsfc    , & !real(kind=r8), intent(in) :: tsfc(pcols)              ! Surface temperature (K)
             h2ovmr  , & ! real(kind=r8), intent(in) :: h2ovmr(pcols,nlay)          ! H2O volume mixing ratio
             o3vmr   , & !real(kind=r8), intent(in) :: o3vmr(pcols,nlay)           ! O3 volume mixing ratio
             co2vmr  , & !real(kind=r8), intent(in) :: co2vmr(pcols,nlay)          ! CO2 volume mixing ratio
             ch4vmr  , & ! real(kind=r8), intent(in) :: ch4vmr(pcols,nlay)          ! Methane volume mixing ratio
             o2vmr   , & !real(kind=r8), intent(in) :: o2vmr(pcols,nlay)           ! O2 volume mixing ratio
             n2ovmr  , & !real(kind=r8), intent(in) :: n2ovmr(pcols,nlay)          ! Nitrous oxide volume mixing ratio
             asdir   , & !real(kind=r8), intent(in) :: asdir(pcols)             ! UV/vis surface albedo direct rad
             asdif   , & !real(kind=r8), intent(in) :: asdif(pcols)             ! UV/vis surface albedo: diffuse rad
             aldir   , & !real(kind=r8), intent(in) :: aldir(pcols)             ! Near-IR surface albedo direct rad
             aldif   , & !real(kind=r8), intent(in) :: aldif(pcols)             ! Near-IR surface albedo: diffuse rad
             coszen  , & !real(kind=r8), intent(in) :: coszen(pcols)            ! Cosine of solar zenith angle
             adjes   , & !real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
             dyofyr  , & ! integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
             solvar  , & !real(kind=r8), intent(in) :: solvar(1:nbndsw)     ! Solar constant (Wm-2) scaling per band
             inflgsw , & !integer, intent(in) :: inflgsw                    ! Flag for cloud optical properties
             iceflgsw, & !integer, intent(in) :: iceflgsw                   ! Flag for ice particle specification
             liqflgsw, & !integer, intent(in) :: liqflgsw                   ! Flag for liquid droplet specification
             cldfmcl , & !real(kind=r8), intent(in) :: cldfmcl(ngptsw,pcols,nlay-1)       ! Cloud fraction
             taucmcl , & !real(kind=r8), intent(in) :: taucmcl(ngptsw,pcols,nlay-1)      ! Cloud optical depth
             ssacmcl , & !real(kind=r8), intent(in) :: ssacmcl(ngptsw,pcols,nlay-1)       ! Cloud single scattering albedo
             asmcmcl , & !real(kind=r8), intent(in) :: asmcmcl(ngptsw,pcols,nlay-1)       ! Cloud asymmetry parameter
             fsfcmcl , & !real(kind=r8), intent(in) :: fsfcmcl(ngptsw,pcols,nlay-1)       ! Cloud forward scattering parameter
             ciwpmcl , & !real(kind=r8), intent(in) :: ciwpmcl(ngptsw,pcols,nlay-1)      ! Cloud ice water path (g/m2)
             clwpmcl , & !real(kind=r8), intent(in) :: clwpmcl(ngptsw,pcols,nlay-1)      ! Cloud liquid water path (g/m2)
             cld                , &! Fractional cloud cover
             cicewp             , &! in-cloud cloud ice water path
             cliqwp             , &! in-cloud cloud liquid water path
             tauc_sw            , &! cloud optical depth  
             ssac_sw            , &! cloud single scat. albedo   
             asmc_sw           , &! cloud asymmetry parameter  
             fsfc_sw            , &! cloud forward scattering fraction

             reicmcl , & !real(kind=r8), intent(in) :: reicmcl(pcols,nlay-1)         ! Cloud ice effective radius (microns)
             relqmcl , & !real(kind=r8), intent(in) :: relqmcl (pcols,nlay-1)         ! Cloud water drop effective radius (microns)
             tauaer  , & !real(kind=r8), intent(inout) :: tauaer(pcols,nlay-1,nbndsw)        ! Aerosol optical depth (iaer=10 only)
             ssaaer  , & !real(kind=r8), intent(inout) :: ssaaer(pcols,nlay-1,nbndsw)       ! Aerosol single scattering albedo (iaer=10 only)
             asmaer  , & !real(kind=r8), intent(inout) :: asmaer(pcols,nlay-1,nbndsw)        ! Aerosol asymmetry parameter (iaer=10 only)
             swuflx  , & !real(kind=r8), intent(out) :: swuflx(pcols,nlay+1)        ! Total sky shortwave upward flux (W/m2)
             swdflx  , & !real(kind=r8), intent(out) :: swdflx(pcols,nlay+1)        ! Total sky shortwave downward flux (W/m2)
             swhr    , & !real(kind=r8), intent(out) :: swhr (pcols,nlay)           ! Total sky shortwave radiative heating rate (K/d)
             swuflxc , & !real(kind=r8), intent(out) :: swuflxc(pcols,nlay+1)   ! Clear sky shortwave upward flux (W/m2)
             swdflxc , & !real(kind=r8), intent(out) :: swdflxc(pcols,nlay+1)        ! Clear sky shortwave downward flux (W/m2)
             swhrc   , & !real(kind=r8), intent(out) :: swhrc(pcols,nlay)          ! Clear sky shortwave radiative heating rate (K/d)
             dirdnuv , & !real(kind=r8), intent(out) :: dirdnuv(pcols,nlay+1)        ! Direct downward shortwave flux, UV/vis
             dirdnir , & !real(kind=r8), intent(out) :: dirdnir(pcols,nlay+1)        ! Direct downward shortwave flux, near-IR
             difdnuv , & !real(kind=r8), intent(out) :: difdnuv(pcols,nlay+1)        ! Diffuse downward shortwave flux, UV/vis
             difdnir , & !real(kind=r8), intent(out) :: difdnir(pcols,nlay+1)        ! Diffuse downward shortwave flux, near-IR
             dirdnuvc, & !real(kind=r8), intent(out) :: dirdnuvc(pcols,nlay+1)       ! Direct downward shortwave flux, UV/vis
             difdnuvc, & !real(kind=r8), intent(out) :: difdnuvc(pcols,nlay+1)       ! Diffuse downward shortwave flux, UV/vis
             dirdnirc, & !real(kind=r8), intent(out) :: dirdnirc(pcols,nlay+1)       ! Direct downward shortwave flux, near-IR
             difdnirc, & !real(kind=r8), intent(out) :: difdnirc(pcols,nlay+1)       ! Diffuse downward shortwave flux, near-IR
             ninflx  , & !real(kind=r8), intent(out) :: ninflx(pcols,nlay+1)         ! Net shortwave flux, near-IR
             ninflxc , & !real(kind=r8), intent(out) :: ninflxc(pcols,nlay+1)        ! Net clear sky shortwave flux, near-IR
             swuflxs , & !real(kind=r8), intent(out)  :: swuflxs(nbndsw,pcols,nlay+1)   ! shortwave spectral flux up
             swdflxs , & !real(kind=r8), intent(out)  :: swdflxs(nbndsw,pcols,nlay+1)   ! shortwave spectral flux down
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013)
               ifaeros,aodsol,asysol,ssasol,z_aersol,topogsol, &
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013)
                aodFsol,asyFsol,ssaFsol,z_aerFsol )
!tar end



! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! Note: The call to RRTMG_SW_INI should be moved to the GCM initialization 
!  area, since this has to be called only once. 
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPRMC_SW to set cloud optical depth for McICA based
!       on input cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVMC to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer
!       using McICA, the Monte-Carlo Independent Column Approximation,
!       to represent sub-grid scale cloud variability
!    f) passes the calculated fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_sw.nomcica.f90 (to not use
!     McICA) or rrtmg_sw.f90 (to use McICA) to interface with a GCM.
!
!    1) Standard, single forward model calculation (imca = 0); this is 
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       This method is valid for clear sky or partial cloud conditions.
!
! This call to RRTMG_SW must be preceeded by a call to the module
!     mcica_subcol_gen_sw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngptsw) dimension.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprop.f90 for further details):
!
!    1) Input cloud fraction, cloud optical depth, single scattering albedo 
!       and asymmetry parameter directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties: ice fracion,
!       ice and liquid particle sizes (inflgsw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflgsw and liqflgsw
!
! Two methods of aerosol property input are possible:
!     Aerosol properties can be input in one of two ways (controlled by input 
!     flag iaer, see text file rrtmg_sw_instructions for further details):
!
!    1) Input aerosol optical depth, single scattering albedo and asymmetry
!       parameter directly by layer and spectral band (iaer=10)
!    2) Input aerosol optical depth and 0.55 micron directly by layer and use
!       one or more of six ECMWF aerosol types (iaer=6)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs. 
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90 
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse 
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays. 
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to output direct and diffuse fluxes either with or without
!   delta scaling based on setting of idelm flag
!     Dec 2008: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2

! ------- Declarations

! ----- Input -----
      integer, intent(in) :: imca                       !                ! flag for mcica [0=off, 1=on]
      integer, intent(in) :: pcols                      ! Total Number of horizontal columns     
      integer, intent(in) :: ncol                       ! Number of horizontal columns     
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(inout) :: icld                    ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum
      real(kind=r8), intent(in) :: play(pcols,nlay)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev (pcols,nlay+1)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(pcols,nlay)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev (pcols,nlay+1)            ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(pcols)              ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(pcols,nlay)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(pcols,nlay)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(pcols,nlay)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(pcols,nlay)          ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(pcols,nlay)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(pcols,nlay)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: asdir(pcols)             ! UV/vis surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: aldir(pcols)             ! Near-IR surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: asdif(pcols)             ! UV/vis surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: aldif(pcols)             ! Near-IR surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)

      integer      , intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
      real(kind=r8), intent(in) :: coszen(pcols)            ! Cosine of solar zenith angle
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: solvar(1:nbndsw)     ! Solar constant (Wm-2) scaling per band

      integer, intent(in) :: inflgsw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflgsw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(ngptsw,pcols,nlay-1)       ! Cloud fraction
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl(ngptsw,pcols,nlay-1)      ! Cloud optical depth
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ssacmcl(ngptsw,pcols,nlay-1)       ! Cloud single scattering albedo
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: asmcmcl(ngptsw,pcols,nlay-1)       ! Cloud asymmetry parameter
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: fsfcmcl(ngptsw,pcols,nlay-1)       ! Cloud forward scattering parameter
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(ngptsw,pcols,nlay-1)      ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl(ngptsw,pcols,nlay-1)      ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptsw,ncol,nlay)

      real(kind=r8), intent(in) :: cld     (1:pcols,1:nlay-1)       ! Fractional cloud cover
      real(kind=r8), intent(in) :: cicewp  (1:pcols,1:nlay-1)       ! in-cloud cloud ice water path
      real(kind=r8), intent(in) :: cliqwp  (1:pcols,1:nlay-1)       ! in-cloud cloud liquid water path
      real(kind=r8), intent(in) :: tauc_sw (1:nbndsw,1:pcols,1:nlay-1) ! cloud optical depth  
      real(kind=r8), intent(in) :: ssac_sw (1:nbndsw,1:pcols,1:nlay-1) ! cloud single scat. albedo   
      real(kind=r8), intent(in) :: asmc_sw (1:nbndsw,1:pcols,1:nlay-1) ! cloud asymmetry parameter  
      real(kind=r8), intent(in) :: fsfc_sw (1:nbndsw,1:pcols,1:nlay-1) ! cloud forward scattering fraction

      real(kind=r8), intent(in) :: reicmcl(pcols,nlay-1)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl (pcols,nlay-1)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(inout) :: tauaer(pcols,nlay-1,nbndsw)        ! Aerosol optical depth (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=r8), intent(inout) :: ssaaer(pcols,nlay-1,nbndsw)       ! Aerosol single scattering albedo (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=r8), intent(inout) :: asmaer(pcols,nlay-1,nbndsw)        ! Aerosol asymmetry parameter (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=r8) :: ecaer(ncol,nlay,naerec)        ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                        !    Dimensions: (ncol,nlay,naerec)
                                                        ! (non-delta scaled)      
!tar begin
! input aerosol parameters
!
    INTEGER,    INTENT(in   ) :: ifaeros
    REAL(KIND=r8),    INTENT(in   )  :: aodsol(pcols,14) 
    REAL(KIND=r8),    INTENT(in   )  :: asysol(pcols,14)
    REAL(KIND=r8),    INTENT(in   )  :: ssasol(pcols,14)
    REAL(KIND=r8),    INTENT(in   )  :: z_aersol(pcols,40)       
    REAL(KIND=r8),    INTENT(in   )  :: topogsol (pcols)   
    REAL(KIND=r8),    INTENT(in   )  :: aodFsol(pcols,14) 
    REAL(KIND=r8),    INTENT(in   )  :: asyFsol(pcols,14)
    REAL(KIND=r8),    INTENT(in   )  :: ssaFsol(pcols,14)
    REAL(KIND=r8),    INTENT(in   )  :: z_aerFsol(pcols,40) 
!    
!tar end  
! ----- Output -----

      real(kind=r8), intent(inout) :: swuflx(pcols,nlay+1)        ! Total sky shortwave upward flux (W/m2)
                                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: swdflx(pcols,nlay+1)        ! Total sky shortwave downward flux (W/m2)
                                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: swhr (pcols,nlay)           ! Total sky shortwave radiative heating rate (K/d)
                                                                      !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(inout) :: swuflxc(pcols,nlay+1)       ! Clear sky shortwave upward flux (W/m2)
                                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: swdflxc(pcols,nlay+1)       ! Clear sky shortwave downward flux (W/m2)
                                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(inout) :: swhrc(pcols,nlay)          ! Clear sky shortwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)

      real(kind=r8), intent(inout) :: dirdnuv(pcols,nlay+1)        ! Direct downward shortwave flux, UV/vis
      real(kind=r8), intent(inout) :: difdnuv(pcols,nlay+1)        ! Diffuse downward shortwave flux, UV/vis
      real(kind=r8), intent(inout) :: dirdnir(pcols,nlay+1)        ! Direct downward shortwave flux, near-IR
      real(kind=r8), intent(inout) :: difdnir(pcols,nlay+1)        ! Diffuse downward shortwave flux, near-IR

      real(kind=r8), intent(inout) :: dirdnuvc(pcols,nlay+1)        ! Direct downward shortwave flux, UV/vis
      real(kind=r8), intent(inout) :: difdnuvc(pcols,nlay+1)        ! Diffuse downward shortwave flux, UV/vis
      real(kind=r8), intent(inout) :: dirdnirc(pcols,nlay+1)        ! Direct downward shortwave flux, near-IR
      real(kind=r8), intent(inout) :: difdnirc(pcols,nlay+1)        ! Diffuse downward shortwave flux, near-IR

      real(kind=r8), intent(inout) :: ninflx(pcols,nlay+1)         ! Net shortwave flux, near-IR
      real(kind=r8), intent(inout) :: ninflxc(pcols,nlay+1)        ! Net clear sky shortwave flux, near-IR

      real(kind=r8), intent(inout) :: swuflxs(nbndsw,pcols,nlay+1)   ! shortwave spectral flux up
      real(kind=r8), intent(inout) :: swdflxs(nbndsw,pcols,nlay+1)   ! shortwave spectral flux down

! ----- Output -----

      real(kind=r8) :: swuflx_local(pcols,nlay+1)         ! Total sky shortwave upward flux (W/m2)
                                                            !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: swdflx_local(pcols,nlay+1)         ! Total sky shortwave downward flux (W/m2)
                                                            !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: swhr_local (pcols,nlay)           ! Total sky shortwave radiative heating rate (K/d)
                                                            !    Dimensions: (ncol,nlay)
      real(kind=r8) :: swuflxc_local(pcols,nlay+1)         ! Clear sky shortwave upward flux (W/m2)
                                                            !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: swdflxc_local(pcols,nlay+1)         ! Clear sky shortwave downward flux (W/m2)
                                                            !    Dimensions: (ncol,nlay+1)
      real(kind=r8) :: swhrc_local(pcols,nlay)          ! Clear sky shortwave radiative heating rate (K/d)
                                              !        Dimensions: (ncol,nlay)

      real(kind=r8) :: dirdnuv_local(pcols,nlay+1)          ! Direct downward shortwave flux, UV/vis
      real(kind=r8) :: difdnuv_local(pcols,nlay+1)          ! Diffuse downward shortwave flux, UV/vis
      real(kind=r8) :: dirdnir_local(pcols,nlay+1)          ! Direct downward shortwave flux, near-IR
      real(kind=r8) :: difdnir_local(pcols,nlay+1)          ! Diffuse downward shortwave flux, near-IR

      real(kind=r8) :: dirdnuvc_local(pcols,nlay+1)           ! Direct downward shortwave flux, UV/vis
      real(kind=r8) :: difdnuvc_local(pcols,nlay+1)           ! Diffuse downward shortwave flux, UV/vis
      real(kind=r8) :: dirdnirc_local(pcols,nlay+1)           ! Direct downward shortwave flux, near-IR
      real(kind=r8) :: difdnirc_local(pcols,nlay+1)           ! Diffuse downward shortwave flux, near-IR

      real(kind=r8) :: ninflx_local(pcols,nlay+1)          ! Net shortwave flux, near-IR
      real(kind=r8) :: ninflxc_local(pcols,nlay+1)          ! Net clear sky shortwave flux, near-IR

      real(kind=r8)  :: swuflxs_local(nbndsw,pcols,nlay+1)   ! shortwave spectral flux up
      real(kind=r8)  :: swdflxs_local(nbndsw,pcols,nlay+1)   ! shortwave spectral flux down

! ----- Local -----

! Control
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: icpr                           ! cldprop/cldprmc use flag
      integer :: iout = 0                       ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
      integer :: idelm                          ! delta-m scaling flag
                                                ! [0 = direct and diffuse fluxes are unscaled]
                                                ! [1 = direct and diffuse fluxes are scaled]
                                                ! (total downward fluxes are always delta scaled)
      integer :: isccos                         ! instrumental cosine response flag (inactive)
      integer :: iplon                          ! column loop index
      integer :: i                              ! layer loop index                       ! jk
      integer :: ib                             ! band loop index                        ! jsw
      integer :: ia, ig                         ! indices
      integer :: k                              ! layer loop index
      integer :: ims                            ! value for changing mcica permute seed

      real(kind=r8) :: zepsec, zepzen           ! epsilon
      real(kind=r8) :: zdpgcp                   ! flux to heating conversion ratio

! Atmosphere
      real(kind=r8) :: pavel(nlay)            ! layer pressures (mb) 
      real(kind=r8) :: tavel(nlay)            ! layer temperatures (K)
      real(kind=r8) :: pz(0:nlay)             ! level (interface) pressures (hPa, mb)
      real(kind=r8) :: tz(0:nlay)             ! level (interface) temperatures (K)
      real(kind=r8) :: tbound                   ! surface temperature (K)
      real(kind=r8) :: pdp(nlay)              ! layer pressure thickness (hPa, mb)
      real(kind=r8) :: coldry (nlay)           ! dry air column amount
      real(kind=r8) :: wkl    (mxmol,nlay)        ! molecular amounts (mol/cm-2)

!      real(kind=r8) :: earth_sun               ! function for Earth/Sun distance factor
      real(kind=r8) :: cossza                   ! Cosine of solar zenith angle
      real(kind=r8) :: adjflux(jpband)          ! adjustment for current Earth/Sun distance
!      real(kind=r8) :: solvar(jpband)           ! solar constant scaling factor from rrtmg_sw
                                                !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=r8) :: albdir(nbndsw)           ! surface albedo, direct          ! zalbp
      real(kind=r8) :: albdif(nbndsw)           ! surface albedo, diffuse         ! zalbd

      real(kind=r8) :: taua(nlay,nbndsw)      ! Aerosol optical depth
      real(kind=r8) :: ssaa(nlay,nbndsw)      ! Aerosol single scattering albedo
      real(kind=r8) :: asma(nlay,nbndsw)      ! Aerosol asymmetry parameter

! Atmosphere - setcoef
      integer :: laytrop                        ! tropopause layer index
      integer :: layswtch                       ! 
      integer :: laylow                         ! 
      integer :: jp(nlay)                     ! 
      integer :: jt(nlay)                     !
      integer :: jt1(nlay)                    !

      real(kind=r8) :: colh2o(nlay)           ! column amount (h2o)
      real(kind=r8) :: colco2(nlay)           ! column amount (co2)
      real(kind=r8) :: colo3(nlay)            ! column amount (o3)
      real(kind=r8) :: coln2o(nlay)           ! column amount (n2o)
      real(kind=r8) :: colch4(nlay)           ! column amount (ch4)
      real(kind=r8) :: colo2(nlay)            ! column amount (o2)
      real(kind=r8) :: colmol(nlay)           ! column amount
      real(kind=r8) :: co2mult(nlay)          ! column amount 

      integer :: indself(nlay)
      integer :: indfor(nlay)
      real(kind=r8) :: selffac(nlay)
      real(kind=r8) :: selffrac(nlay)
      real(kind=r8) :: forfac(nlay)
      real(kind=r8) :: forfrac(nlay)
      real(kind=r8) :: cldfrac(nlay)         ! layer cloud fraction

      real(kind=r8) :: &                        !
                         fac00(nlay), fac01(nlay), &
                         fac10(nlay), fac11(nlay) 

! Atmosphere/clouds - cldprop
      integer :: ncbands                        ! number of cloud spectral bands
      integer :: inflag                         ! flag for cloud property method
      integer :: iceflag                        ! flag for ice cloud properties
      integer :: liqflag                        ! flag for liquid cloud properties

!      real(kind=r8) :: cldfrac(nlay)            ! layer cloud fraction
!      real(kind=r8) :: tauc(nlay)               ! cloud optical depth (non-delta scaled)
!      real(kind=r8) :: ssac(nlay)               ! cloud single scattering albedo (non-delta scaled)
!      real(kind=r8) :: asmc(nlay)               ! cloud asymmetry parameter (non-delta scaled)
!      real(kind=r8) :: ciwp(nlay)               ! cloud ice water path
!      real(kind=r8) :: clwp(nlay)               ! cloud liquid water path
!      real(kind=r8) :: rei(nlay)                ! cloud ice particle size
!      real(kind=r8) :: rel(nlay)                ! cloud liquid particle size

      real(kind=r8) :: taucloud  (nlay,jpband)    ! cloud optical depth
      real(kind=r8) :: taucldorig(nlay,jpband)  ! cloud optical depth (non-delta scaled)
      real(kind=r8) :: ssacloud  (nlay,jpband)    ! cloud single scattering albedo
      real(kind=r8) :: asmcloud  (nlay,jpband)    ! cloud asymmetry parameter

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=r8) :: cldfmc(ngptsw,nlay)    ! cloud fraction [mcica]
      real(kind=r8) :: ciwpmc(ngptsw,nlay)    ! cloud ice water path [mcica]
      real(kind=r8) :: clwpmc(ngptsw,nlay)    ! cloud liquid water path [mcica]
      real(kind=r8) :: relqmc(nlay)           ! liquid particle size (microns)
      real(kind=r8) :: reicmc(nlay)           ! ice particle effective radius (microns)
      real(kind=r8) :: dgesmc(nlay)           ! ice particle generalized effective size (microns)
      real(kind=r8) :: taucmc(ngptsw,nlay)    ! cloud optical depth [mcica]
      real(kind=r8) :: taormc(ngptsw,nlay)    ! unscaled cloud optical depth [mcica]
      real(kind=r8) :: ssacmc(ngptsw,nlay)    ! cloud single scattering albedo [mcica]
      real(kind=r8) :: asmcmc(ngptsw,nlay)    ! cloud asymmetry parameter [mcica]
      real(kind=r8) :: fsfcmc(ngptsw,nlay)    ! cloud forward scattering fraction [mcica]

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=r8) :: ztauc(nlay,nbndsw)     ! cloud optical depth
      real(kind=r8) :: ztaucorig(nlay,nbndsw) ! unscaled cloud optical depth
      real(kind=r8) :: zasyc(nlay,nbndsw)     ! cloud asymmetry parameter 
                                                !  (first moment of phase function)
      real(kind=r8) :: zomgc(nlay,nbndsw)     ! cloud single scattering albedo
      real(kind=r8) :: ztaua(nlay,nbndsw)     ! total aerosol optical depth
      real(kind=r8) :: zasya(nlay,nbndsw)     ! total aerosol asymmetry parameter 
      real(kind=r8) :: zomga(nlay,nbndsw)     ! total aerosol single scattering albedo

      real(kind=r8) :: zcldfmc(nlay,ngptsw)   ! cloud fraction [mcica]
      real(kind=r8) :: ztaucmc(nlay,ngptsw)   ! cloud optical depth [mcica]
      real(kind=r8) :: ztaormc(nlay,ngptsw)   ! unscaled cloud optical depth [mcica]
      real(kind=r8) :: zasycmc(nlay,ngptsw)   ! cloud asymmetry parameter [mcica] 
      real(kind=r8) :: zomgcmc(nlay,ngptsw)   ! cloud single scattering albedo [mcica]

      real(kind=r8) :: zbbfu(nlay+2)          ! temporary upward shortwave flux (w/m2)
      real(kind=r8) :: zbbfd(nlay+2)          ! temporary downward shortwave flux (w/m2)
      real(kind=r8) :: zbbcu(nlay+2)          ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=r8) :: zbbcd(nlay+2)          ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=r8) :: zbbfddir(nlay+2)       ! temporary downward direct shortwave flux (w/m2)
      real(kind=r8) :: zbbcddir(nlay+2)       ! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=r8) :: zuvfd(nlay+2)          ! temporary UV downward shortwave flux (w/m2)
      real(kind=r8) :: zuvcd(nlay+2)          ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=r8) :: zuvfddir(nlay+2)       ! temporary UV downward direct shortwave flux (w/m2)
      real(kind=r8) :: zuvcddir(nlay+2)       ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=r8) :: znifd(nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=r8) :: znicd(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=r8) :: znifddir(nlay+2)       ! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=r8) :: znicddir(nlay+2)       ! temporary clear sky near-IR downward direct shortwave flux (w/m2)
! Added for near-IR flux diagnostic
      real(kind=r8) :: znifu(nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=r8) :: znicu(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)

! Optional output fields 
      real(kind=r8) :: swnflx(nlay+2)         ! Total sky shortwave net flux (W/m2)
      real(kind=r8) :: swnflxc(nlay+2)        ! Clear sky shortwave net flux (W/m2)
      real(kind=r8) :: dirdflux(nlay+2)       ! Direct downward shortwave surface flux
      real(kind=r8) :: difdflux(nlay+2)       ! Diffuse downward shortwave surface flux
      real(kind=r8) :: uvdflx(nlay+2)         ! Total sky downward shortwave flux, UV/vis   
      real(kind=r8) :: nidflx(nlay+2)         ! Total sky downward shortwave flux, near-IR  
      real(kind=r8) :: zbbfsu(nbndsw,nlay+2)  ! temporary upward shortwave flux spectral (w/m2)
      real(kind=r8) :: zbbfsd(nbndsw,nlay+2)  ! temporary downward shortwave flux spectral (w/m2)

! Output - inactive
!      real(kind=r8) :: zuvfu(nlay+2)         ! temporary upward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvfd(nlay+2)         ! temporary downward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvcu(nlay+2)         ! temporary clear sky upward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvcd(nlay+2)         ! temporary clear sky downward UV shortwave flux (w/m2)
!      real(kind=r8) :: zvsfu(nlay+2)         ! temporary upward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvsfd(nlay+2)         ! temporary downward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvscu(nlay+2)         ! temporary clear sky upward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvscd(nlay+2)         ! temporary clear sky downward visible shortwave flux (w/m2)
!      real(kind=r8) :: znifu(nlay+2)         ! temporary upward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znifd(nlay+2)         ! temporary downward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znicu(nlay+2)         ! temporary clear sky upward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znicd(nlay+2)         ! temporary clear sky downward near-IR shortwave flux (w/m2)


!tar begin
! Climate aerosol parameters
    REAL(KIND=r8), DIMENSION(ncol,nlay,14)  :: taual,ssaal,asyal
!
!Fine mode
    REAL(KIND=r8), DIMENSION(ncol,nlay,14)  :: taualF,ssaalF,asyalF
!Coarse mode    
    REAL(KIND=r8), DIMENSION(ncol,nlay,14)  :: taualC,ssaalC,asyalC    
!
! Aerosol local variables 
    REAL(KIND=r8) :: pr(nlay+1)     !levels pressure from surface in Pa
    REAL(KIND=r8) :: tmp(nlay)       ! layer temperature from surface in K
    REAL(KIND=r8) :: tpg          !grid topography
    REAL(KIND=r8) :: zlv(nlay+1)   !level height from surface   (m) 
    REAL(KIND=r8) :: dz(nlay)       ! layer depth from surface  (m)
!    REAL(KIND=r8) :: dzint(ncol,nlay)       ! layer depth from surface  (m)
    INTEGER :: j
    REAL(KIND=r8) :: ep(40),f(42),ha(41)  ! extinction in CA layers of coarse mode
    REAL(KIND=r8) :: fp(nlay)      !extinction in model layers from surface of coarse mode
    REAL(KIND=r8) :: od_n(nlay)    ! AOD in model layers from surface  of coarse mode           
    REAL(KIND=r8) :: aod_norm(ncol,nlay) ! extinction in model layers from top  of coarse mode     
    REAL(KIND=r8) :: epF(40),fF(42),haF(41)  ! extinction in CA layers of fine mode
    REAL(KIND=r8) :: fpF(nlay)      !extinction in model layers from surface of fine mode
    REAL(KIND=r8) :: od_nF(nlay)    ! AOD in model layers from surface of fine mode           
    REAL(KIND=r8) :: aod_normF(ncol,nlay) ! extinction in model layers from top of fine mode 
      
!tar end 

! Initializations
       swuflx_local=0.0_r8;     swdflx_local=0.0_r8;     swhr_local =0.0_r8;      swuflxc_local=0.0_r8;    
       swdflxc_local=0.0_r8;    swhrc_local=0.0_r8;      dirdnuv_local=0.0_r8;    difdnuv_local=0.0_r8;   
       dirdnir_local=0.0_r8;    difdnir_local=0.0_r8;    
       dirdnuvc_local=0.0_r8;   difdnuvc_local=0.0_r8;   dirdnirc_local=0.0_r8;   difdnirc_local=0.0_r8;   
       ninflx_local=0.0_r8;     ninflxc_local=0.0_r8;    swuflxs_local=0.0_r8;    swdflxs_local=0.0_r8;    

       indself=0
       indfor=0
       pavel=0.0_r8;tavel=0.0_r8;pz=0.0_r8;tz=0.0_r8
       pdp=0.0_r8;coldry=0.0_r8;wkl=0.0_r8;adjflux=0.0_r8;albdir=0.0_r8;albdif=0.0_r8;
       taua=0.0_r8;ssaa=0.0_r8;asma=0.0_r8;colh2o=0.0_r8;colco2=0.0_r8;colo3=0.0_r8;coln2o=0.0_r8;
       colch4=0.0_r8;colo2=0.0_r8;colmol=0.0_r8;co2mult=0.0_r8;

       selffac=0.0_r8;selffrac=0.0_r8;forfac=0.0_r8;forfrac=0.0_r8;


       fac00=0.0_r8;fac01=0.0_r8;fac10=0.0_r8;fac11=0.0_r8;


       cldfmc=0.0_r8;ciwpmc=0.0_r8;clwpmc=0.0_r8;relqmc=0.0_r8;reicmc=0.0_r8;dgesmc=0.0_r8;taucmc=0.0_r8;
       taormc=0.0_r8;ssacmc=0.0_r8;asmcmc=0.0_r8;fsfcmc=0.0_r8;
       taucldorig=0.0_r8;taucloud  =0.0_r8;
       ssacloud  =0.0_r8;asmcloud  =0.0_r8;cldfrac=0.0_r8
! Atmosphere/clouds/aerosol - spcvrt,spcvmc
       ztauc=0.0_r8;ztaucorig=0.0_r8;zasyc=0.0_r8;
                            
       zomgc=0.0_r8;ztaua=0.0_r8;zasya=0.0_r8;zomga=0.0_r8;

       zcldfmc=0.0_r8;ztaucmc=0.0_r8;ztaormc=0.0_r8;zasycmc=0.0_r8;zomgcmc=0.0_r8;zbbfu=0.0_r8;zbbfd=0.0_r8;
       zbbcu=0.0_r8;zbbcd=0.0_r8;zbbfddir=0.0_r8;zbbcddir=0.0_r8;zuvfd=0.0_r8;zuvcd=0.0_r8;
       zuvfddir=0.0_r8;zuvcddir=0.0_r8;znifd=0.0_r8;znicd=0.0_r8;znifddir=0.0_r8;znicddir=0.0_r8;
! Added for near-IR flux diagnostic
       znifu=0.0_r8;
       znicu=0.0_r8;

! Optional output fields 
       swnflx=0.0_r8;    
       swnflxc=0.0_r8;   
       dirdflux=0.0_r8;  
       difdflux=0.0_r8;  
       uvdflx=0.0_r8;    
       nidflx=0.0_r8;    
       zbbfsu=0.0_r8;    
       zbbfsd=0.0_r8;    












      zepsec = 1.e-06_r8
      zepzen = 1.e-10_r8
      oneminus = 1.0_r8 - zepsec
      pi = 2._r8 * asin(1._r8)

      istart = jpb1
      iend = jpb2
      icpr = 0
      ims = 2
!      
!tar begin
!
! Climate aerosol (coarse + fine modes)  (Kinne, 2013)  
      IF(ifaeros==2) THEN 

! Initialize local vectors and output variables
!    
         taual=0.0_r8; ssaal=0.0_r8; asyal=0.0_r8
!
! Fine mode
         taualF=0.0_r8; ssaalF=0.0_r8; asyalF=0.0_r8
! Coarse mode    
         taualC=0.0_r8; ssaalC=0.0_r8; asyalC=0.0_r8 
!Initialize new aerosol local vectors
         pr=0.0_r8; tmp=0.0_r8; tpg=0.0_r8; zlv=0.0_r8; dz=0.0_r8; ep=0.0_r8
         f=0.0_r8; ha=0.0_r8; fp=0.0_r8; od_n=0.0_r8; aod_norm=0.0_r8
!Initialize new aerosol local vectors of fine mode 
         epF=0.0_r8;fF=0.0_r8; haF=0.0_r8; fpF=0.0_r8; od_nF=0.0_r8; aod_normF=0.0_r8
!       
!   
        DO i=1,ncol
!
           DO k=1,nlay
              tmp(k)=tlay(i,nlay-k+1)  !Flipping from surface
           END DO 
!
           DO k=1,nlay+1     
              pr(k)=plev(i,nlay-k+2)*100.0_r8  !Flipping from surface in Pa   
           END DO 
!
           tpg=topogsol(i)  !topography
!
           DO j=1,40
              ep(j)=z_aersol(i,j)  ! Extinction profile from surface, coarse mode
!
              epF(j)=z_aerFsol(i,j)  ! Extinction profile from surface, fine mode
!
           END DO
!
           CALL zplev(pr,tmp,tpg,nlay, zlv,dz)
! coarse aerosol mode
           CALL zaero(ep,40, f,ha)
!
           CALL aeros_interp(zlv,ha,f,40,nlay, fp)
!
           CALL aod_n(fp,dz,nlay, od_n)
!
! fine aerosol mode
           CALL zaero(epF,40, fF,haF)
!
           CALL   aeros_interp(zlv,haF,fF,40,nlay, fpF)
!
           CALL   aod_n(fpF,dz,nlay, od_nF)
!
!
!
           DO k=1,nlay
              aod_norm(i,k)=od_n(nlay-k+1)  !Flipping from top, coarse mode
!
              aod_normF(i,k)=od_nF(nlay-k+1)  !Flipping from top, fine mode
!
           END DO
!Layer depth in km
!           DO k=1,nlay        
!              dzint(i,k)=dz(nlay-k+1)*0.001_r8 !Flipping from top 
!           END DO
!
        END DO   !do i=1,ncol


        DO k=1,nlay
           DO ib=1,14
              DO i=1,ncol
                 taualC(i,k,ib) = aodsol(i,ib)*aod_norm(i,k)
! 
! Fine aerosol mode included

                 taualF(i,k,ib) = aodFsol(i,ib)*aod_normF(i,k)
!
                 IF(ssasol(i,ib).GE.0.99_r8) THEN
                    ssaalC(i,k,ib)=0.99_r8
                 ELSE
                    ssaalC(i,k,ib) =ssasol(i,ib)
                 END IF
!
                 IF(ssaFsol(i,ib).GE.0.99_r8) THEN
                    ssaalF(i,k,ib)=0.99_r8
                 ELSE
                    ssaalF(i,k,ib) =ssaFsol(i,ib)
                 END IF
!
                 asyalC(i,k,ib) =asysol(i,ib) 
!
                 asyalF(i,k,ib) =asyFsol(i,ib)
!
!                taual(i,k,ib)=taualF(i,k,ib)
!                ssaal(i,k,ib)=ssaalF(i,k,ib)
!                asyal(i,k,ib)=asyalF(i,k,ib)
!
                 taual(i,k,ib)=taualF(i,k,ib)+taualC(i,k,ib)  
!  
                 IF (taual(i,k,ib)==0.0_r8) THEN
                    ssaal(i,k,ib)=0.99_r8
                    asyal(i,k,ib)=0.7_r8
!    
                 ELSE
                    ssaal(i,k,ib)=(ssaalF(i,k,ib)*taualF(i,k,ib)+ &
                    ssaalC(i,k,ib)*taualC(i,k,ib))/ taual(i,k,ib)
!  
                    asyal(i,k,ib)=(asyalF(i,k,ib)*ssaalF(i,k,ib)*taualF(i,k,ib)+   &
                    asyalC(i,k,ib)*ssaalC(i,k,ib)*taualC(i,k,ib))/  &  
                    (ssaalF(i,k,ib)*taualF(i,k,ib)+ssaalC(i,k,ib)*taualC(i,k,ib))
! 
                 END IF
!
!
              END DO
           END DO
        END DO 
!
        DO k=1,nlay-1
           DO ib=1,14
              DO i=1,ncol 
                 tauaer(i,k,ib)=taual(i,nlay-k+1,ib)  !Flipping from surface
                 ssaaer(i,k,ib)=ssaal(i,nlay-k+1,ib)  !Flipping from surface
!                 asyaer(i,k,ib)=asyal(i,nlay-k+1,ib)  !Flipping from surface
              END DO
           END DO
        END DO
   END IF    !IF(ifaeros==2) THEN 

! tar end      
      
      

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation (clear and overcast only)
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!            (clear, overcast or partial cloud conditions)

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud 
! overlap method (read by subroutine readprof from input file INPUT_RRTM):  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 6, use six ECMWF aerosol types
!           input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
! iaer = 10, input total aerosol optical depth, single scattering albedo 
!            and asymmetry parameter (tauaer, ssaaer, asmaer) directly
!PK      iaer = 10
       iaer = 10
! Set idelm to select between delta-M scaled or unscaled output direct and diffuse fluxes
! NOTE: total downward fluxes are always delta scaled
! idelm = 0, output direct and diffuse flux components are not delta scaled
!            (direct flux does not include forward scattering peak)
! idelm = 1, output direct and diffuse flux components are delta scaled (default)
!            (direct flux includes part or most of forward scattering peak)
      idelm = 1

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_sw_ini

! This is the main longitude/column loop in RRTMG.
! Modify to loop over all columns (nlon) or over daylight columns

      do iplon = 1, ncol

! Prepare atmosphere profile from GCM for use in RRTMG, and define
! other input parameters

         call inatm_sw (    &
             iplon                                   , &
             pcols                                   , &
             nlay                                    , &
             icld                                    , &
             iaer                                    , &
             play         (1:pcols,1:nlay)           , &
             plev         (1:pcols,1:nlay+1)         , &
             tlay         (1:pcols,1:nlay)           , &
             tlev         (1:pcols,1:nlay+1)         , &
             tsfc         (1:pcols)                  , &
             h2ovmr       (1:pcols,1:nlay)           , &
             o3vmr        (1:pcols,1:nlay)           , &
             co2vmr       (1:pcols,1:nlay)           , &
             ch4vmr       (1:pcols,1:nlay)           , &
             o2vmr        (1:pcols,1:nlay)           , &
             n2ovmr       (1:pcols,1:nlay)           , &
             adjes                                   , &
             dyofyr                                  , &
             solvar       (1:nbndsw)                 , &
             inflgsw                                 , &
             iceflgsw                                , &
             liqflgsw                                , &
             cldfmcl      (1:ngptsw,1:pcols,1:nlay-1), &
             taucmcl      (1:ngptsw,1:pcols,1:nlay-1), &
             ssacmcl      (1:ngptsw,1:pcols,1:nlay-1), &
             asmcmcl      (1:ngptsw,1:pcols,1:nlay-1), &
             fsfcmcl      (1:ngptsw,1:pcols,1:nlay-1), &
             ciwpmcl      (1:ngptsw,1:pcols,1:nlay-1), &
             clwpmcl      (1:ngptsw,1:pcols,1:nlay-1), &
             reicmcl      (1:pcols,1:nlay-1)         , &
             relqmcl      (1:pcols,1:nlay-1)         , &
             tauaer       (1:pcols,1:nlay-1,1:nbndsw), &
             ssaaer       (1:pcols,1:nlay-1,1:nbndsw), &
             asmaer       (1:pcols,1:nlay-1,1:nbndsw), &
             pavel        (1:nlay)                   , &
             pz           (0:nlay)                   , &
             pdp          (1:nlay)                   , &
             tavel        (1:nlay)                   , &
             tz           (0:nlay)                   , &
             tbound                                  , &
             coldry       (1:nlay)                   , &
             wkl          (1:mxmol,1:nlay)           , &
             adjflux      (1:jpband)                 , &
             inflag                                  , &
             iceflag                                 , &
             liqflag                                 , &
             cldfmc       (1:ngptsw,1:nlay)          , &
             taucmc       (1:ngptsw,1:nlay)          , &
             ssacmc       (1:ngptsw,1:nlay)          , &
             asmcmc       (1:ngptsw,1:nlay)          , &
             fsfcmc       (1:ngptsw,1:nlay)          , &
             ciwpmc       (1:ngptsw,1:nlay)          , &
             clwpmc       (1:ngptsw,1:nlay)          , &
             reicmc       (1:nlay)                   , &
             dgesmc       (1:nlay)                   , &
             relqmc       (1:nlay)                   , &
             taua         (1:nlay,1:nbndsw)          , &
             ssaa         (1:nlay,1:nbndsw)          , &
             asma         (1:nlay,1:nbndsw)            )

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  

!  If McICA is requested use cloud fraction and cloud physical properties 
!  generated by sub-column cloud generator above. 
      cldfrac(1:nlay-1)=cld(iplon,1:nlay-1)

         if (imca.eq.0) then
            do i = 1, nlay
               if (cldfrac(i).ne.0.0_r8 .and. cldfrac(i).ne.1.0_r8) then
                     cldfrac(i)=MIN(MAX(cldfrac(i),0.00001_r8),0.9999_r8)
                     !stop 'PARTIAL CLOUD NOT ALLOWED FOR IMCA=0'
               endif
            enddo
               call cldprop_sw(nlay-1   , &
                               inflag   , &
                               iceflag  , &
                               liqflag  , &
                               cldfrac   (1:nlay-1)                 , &
                               tauc_sw   (1:nbndsw,iplon,1:nlay-1)  , &
                               ssac_sw   (1:nbndsw,iplon,1:nlay-1)  , &
                               asmc_sw   (1:nbndsw,iplon,1:nlay-1)  , &
                               cicewp    (iplon,1:nlay-1)           , &
                               cliqwp    (iplon,1:nlay-1)           , &
                               reicmcl   (iplon,1:nlay-1)           , &
                               2*reicmcl (iplon,1:nlay-1)           , &
                               relqmcl   (iplon,1:nlay-1)           , &
                               taucldorig(1:nlay-1,1:jpband)          , &
                               taucloud  (1:nlay-1,1:jpband)          , &
                               ssacloud  (1:nlay-1,1:jpband)          , &
                               asmcloud  (1:nlay-1,1:jpband)            )
               icpr = 1

         else
            call cldprmc_sw( &
                         nlay                                , &
                         inflag                              , &
                         iceflag                             , &
                         liqflag                             , &
                         cldfmc          (1:ngptsw,1:nlay)   , &
                         ciwpmc          (1:ngptsw,1:nlay)   , &
                         clwpmc          (1:ngptsw,1:nlay)   , &
                         reicmc          (1:nlay)            , &
                         dgesmc          (1:nlay)            , &
                         relqmc          (1:nlay)            , &
                         taormc          (1:ngptsw,1:nlay)   , & 
                         taucmc          (1:ngptsw,1:nlay)   , & 
                         ssacmc          (1:ngptsw,1:nlay)   , & 
                         asmcmc          (1:ngptsw,1:nlay)   , & 
                         fsfcmc          (1:ngptsw,1:nlay)   ) 
            icpr = 1

         endif
! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

         call setcoef_sw( &
                         nlay                  , &
                         pavel         (1:nlay), &
                         tavel         (1:nlay), &
                         pz            (0:nlay), &
                         tz            (0:nlay)        , &
                         tbound                        , &
                         coldry        (1:nlay)        , &
                         wkl           (1:mxmol,1:nlay), &
                         laytrop               , &
                         layswtch              , &
                         laylow                , &
                         jp            (1:nlay), &
                         jt            (1:nlay), &
                         jt1           (1:nlay), &
                         co2mult       (1:nlay), &
                         colch4        (1:nlay), &
                         colco2        (1:nlay), &
                         colh2o        (1:nlay), &
                         colmol        (1:nlay), &
                         coln2o        (1:nlay), &
                         colo2         (1:nlay), &
                         colo3         (1:nlay), &
                         fac00         (1:nlay), &
                         fac01         (1:nlay), &
                         fac10         (1:nlay), &
                         fac11         (1:nlay), &
                         selffac       (1:nlay), &
                         selffrac      (1:nlay), &
                         indself       (1:nlay), &
                         forfac        (1:nlay), &
                         forfrac       (1:nlay), &
                         indfor        (1:nlay)  )

! Cosine of the solar zenith angle 
!  Prevent using value of zero; ideally, SW model is not called from host model when sun 
!  is below horizon

         cossza = coszen(iplon)
         if (cossza .lt. zepzen) cossza = zepzen

! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 

! Surface albedo
!  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
!         do ib=1,9
         do ib=1,8
            albdir(ib) = aldir(iplon)
            albdif(ib) = aldif(iplon)
         enddo
         albdir(nbndsw) = aldir(iplon)
         albdif(nbndsw) = aldif(iplon)
!  Set band 24 (or, band 9 counting from 1) to use linear average of UV/visible
!  and near-IR values, since this band straddles 0.7 microns: 
         albdir(9) = 0.5*(aldir(iplon) + asdir(iplon))
         albdif(9) = 0.5*(aldif(iplon) + asdif(iplon))
!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
         do ib=10,13
            albdir(ib) = asdir(iplon)
            albdif(ib) = asdif(iplon)
         enddo


! Clouds
         if (icld.eq.0) then

            ztauc(:,:) = 0.0_r8
            ztaucorig(:,:) = 0.0_r8
            zasyc(:,:) = 0.0_r8
            zomgc(:,:) = 1.0_r8
            zcldfmc(:,:) = 0._r8
            ztaucmc(:,:) = 0._r8
            ztaormc(:,:) = 0._r8
            zasycmc(:,:) = 0._r8
            zomgcmc(:,:) = 1._r8

         elseif (icld.ge.1) then
            if (imca.eq.0) then
               do i=1,nlay
                  do ib=1,nbndsw
                     if (cldfrac(i) .ge. zepsec) then
                        ztauc(i,ib) = taucloud(i,jpb1-1+ib)
                        ztaucorig(i,ib) = taucldorig(i,jpb1-1+ib)
                        zasyc(i,ib) = asmcloud(i,jpb1-1+ib)
                        zomgc(i,ib) = ssacloud(i,jpb1-1+ib)
                     endif
                  enddo
               enddo
            else
               do i=1,nlay
                  do ig=1,ngptsw
                     zcldfmc(i,ig) = cldfmc(ig,i)
                     ztaucmc(i,ig) = taucmc(ig,i)
                     ztaormc(i,ig) = taormc(ig,i)
                     zasycmc(i,ig) = asmcmc(ig,i)
                     zomgcmc(i,ig) = ssacmc(ig,i)
                  enddo
               enddo
            endif   
         endif   

! Aerosol
! IAER = 0: no aerosols
         if (iaer.eq.0) then

            ztaua(:,:) = 0._r8
            zasya(:,:) = 0._r8
            zomga(:,:) = 1._r8

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer), 
! or set manually here for each aerosol and layer.
         elseif (iaer.eq.6) then

            do i = 1, nlay
               do ia = 1, naerec
                  ecaer(iplon,i,ia) =0.0e-15_r8
               enddo
            enddo

            do i = 1, nlay
               do ib = 1, nbndsw
                  ztaua(i,ib) = 0._r8
                  zasya(i,ib) = 0._r8
                  zomga(i,ib) = 1._r8
                  do ia = 1, naerec
                     ztaua(i,ib) = ztaua(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia)
                    zomga(i,ib) = zomga(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
                                   rsrpiza(ib,ia)
                     zasya(i,ib) = zasya(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
                                   rsrpiza(ib,ia) * rsrasya(ib,ia)
                  enddo
                  if (zomga(i,ib) /= 0._r8) then
                     zasya(i,ib) = zasya(i,ib) / zomga(i,ib)
                  endif
                  if (ztaua(i,ib) /= 0._r8) then
                     zomga(i,ib) = zomga(i,ib) / ztaua(i,ib)
                  endif
               enddo
            enddo

! IAER=10: Direct specification of aerosol optical properties from GCM
         elseif (iaer.eq.10) then

            do i = 1 ,nlay
               do ib = 1 ,nbndsw
                  ztaua(i,ib) = taua(i,ib)
                  zasya(i,ib) = asma(i,ib)
                  zomga(i,ib) = ssaa(i,ib)
               enddo
            enddo

         endif


! Call the 2-stream radiation transfer model

         do i=1,nlay+1
            zbbcu(i) = 0._r8
            zbbcd(i) = 0._r8
            zbbfu(i) = 0._r8
            zbbfd(i) = 0._r8
            zbbcddir(i) = 0._r8
            zbbfddir(i) = 0._r8
            zuvcd(i) = 0._r8
            zuvfd(i) = 0._r8
            zuvcddir(i) = 0._r8
            zuvfddir(i) = 0._r8  
            znicd(i) = 0._r8
            znifd(i) = 0._r8
            znicddir(i) = 0._r8
            znifddir(i) = 0._r8
    
            znicu(i) = 0._r8
            znifu(i) = 0._r8
            zbbfsu(:,i) = 0._r8
            zbbfsd(:,i) = 0._r8
         enddo

         if (imca.eq.0) then
            call spcvrt_sw (&
              nlay                                , &! nlay                             ,&!
              istart                              , &! istart                           ,&!
              iend                                , &! iend                             ,&!
              icpr                                , &! icpr                             ,&!
              idelm                               , &! idelm                            ,&!
              iout                                , &! iout                             ,&!
              pavel         (1:nlay)              , &! pavel      (1:nlay)              ,&!
              tavel         (1:nlay)              , &! tavel      (1:nlay)              ,&!
              pz            (0:nlay)              , &! pz         (0:nlay)              ,&!
              tz            (0:nlay)              , &! tz         (0:nlay)              ,&!
              tbound                              , &! tbound                           ,&!
              albdif        (1:nbndsw)            , &! albdif     (1:nbndsw)            ,&!
              albdir        (1:nbndsw)            , &! albdir     (1:nbndsw)            ,&!
              cldfrac       (1:nlay)              , &! zcldfmc    (1:nlay,1:ngptsw)     ,&!
              ztauc         (1:nlay,1:nbndsw)     , &! ztaucmc    (1:nlay,1:ngptsw)     ,&!
              zasyc         (1:nlay,1:nbndsw)     , &! zasycmc    (1:nlay,1:ngptsw)     ,&!
              zomgc         (1:nlay,1:nbndsw)     , &! zomgcmc    (1:nlay,1:ngptsw)     ,&!
              ztaucorig     (1:nlay,1:nbndsw)     , &! ztaormc    (1:nlay,1:ngptsw)     ,&!
              ztaua         (1:nlay,1:nbndsw)     , &! ztaua      (1:nlay,1:nbndsw)     ,&!
              zasya         (1:nlay,1:nbndsw)     , &! zasya      (1:nlay,1:nbndsw)     ,&!
              zomga         (1:nlay,1:nbndsw)     , &! zomga      (1:nlay,1:nbndsw)     ,&!
              cossza                              , &! cossza                           ,&!
              coldry        (1:nlay)              , &! coldry     (1:nlay)              ,&!
              wkl           (1:mxmol,1:nlay)      , &! wkl        (1:mxmol,1:nlay)      ,&!
              adjflux       (1:jpband)            , &! adjflux    (1:jpband)            ,&!
              laytrop                             , &! laytrop                          ,&!
              layswtch                            , &! layswtch                         ,&!
              laylow                              , &! laylow                           ,&!
              jp            (1:nlay)              , &! jp          (1:nlay)             ,&!
              jt            (1:nlay)              , &! jt          (1:nlay)             ,&!
              jt1           (1:nlay)              , &! jt1         (1:nlay)             ,&!
              co2mult       (1:nlay)              , &! co2mult     (1:nlay)             ,&!
              colch4        (1:nlay)              , &! colch4      (1:nlay)             ,&!
              colco2        (1:nlay)              , &! colco2      (1:nlay)             ,&!
              colh2o        (1:nlay)              , &! colh2o      (1:nlay)             ,&!
              colmol        (1:nlay)              , &! colmol      (1:nlay)             ,&!
              coln2o        (1:nlay)              , &! coln2o      (1:nlay)             ,&!
              colo2         (1:nlay)              , &! colo2       (1:nlay)             ,&!
              colo3         (1:nlay)              , &! colo3       (1:nlay)             ,&!
              fac00         (1:nlay)              , &! fac00       (1:nlay)             ,&!
              fac01         (1:nlay)              , &! fac01       (1:nlay)             ,&!
              fac10         (1:nlay)              , &! fac10       (1:nlay)             ,&!
              fac11         (1:nlay)              , &! fac11       (1:nlay)             ,&!
              selffac       (1:nlay)              , &! selffac     (1:nlay)             ,&!
              selffrac      (1:nlay)              , &! selffrac    (1:nlay)             ,&!
              indself       (1:nlay)              , &! indself     (1:nlay)             ,&!
              forfac        (1:nlay)              , &! forfac      (1:nlay)             ,&!
              forfrac       (1:nlay)              , &! forfrac     (1:nlay)             ,&!
              indfor        (1:nlay)              , &! indfor      (1:nlay)             ,&!
              zbbfd         (1:nlay+2)            , &! zbbfd       (1:nlay+2)           ,&!
              zbbfu         (1:nlay+2)            , &! zbbfu       (1:nlay+2)           ,&!
              zbbcd         (1:nlay+2)            , &! zbbcd       (1:nlay+2)           ,&!
              zbbcu         (1:nlay+2)            , &! zbbcu       (1:nlay+2)           ,&!
              zuvfd         (1:nlay+2)            , &! zuvfd       (1:nlay+2)           ,&!
              zuvcd         (1:nlay+2)            , &! zuvcd       (1:nlay+2)           ,&!
              znifd         (1:nlay+2)            , &! znifd       (1:nlay+2)           ,&!
              znicd         (1:nlay+2)            , &! znicd       (1:nlay+2)           ,&!
              znifu         (1:nlay+2)            , &! znifu       (1:nlay+2)           ,&!
              znicu         (1:nlay+2)            , &! znicu       (1:nlay+2)           ,&!
              zbbfddir      (1:nlay+2)            , &! zbbfddir    (1:nlay+2)           ,&!
              zbbcddir      (1:nlay+2)            , &! zbbcddir    (1:nlay+2)           ,&!
              zuvfddir      (1:nlay+2)            , &! zuvfddir    (1:nlay+2)           ,&!
              zuvcddir      (1:nlay+2)            , &! zuvcddir    (1:nlay+2)           ,&!
              znifddir      (1:nlay+2)            , &! znifddir    (1:nlay+2)           ,&!
              znicddir      (1:nlay+2)              )! znicddir    (1:nlay+2)           ,&!
         else
            call spcvmc_sw( &
              iplon                          ,&!
              nlay                           ,&!
              istart                         ,&!
              iend                           ,&!
              icpr                           ,&!
              idelm                          ,&!
              iout                           ,&!
              pavel      (1:nlay)            ,&!
              tavel      (1:nlay)            ,&!
              pz         (0:nlay)            ,&!
              tz         (0:nlay)            ,&!
              tbound                         ,&!
              albdif     (1:nbndsw)          ,&!
              albdir     (1:nbndsw)          ,&!
              zcldfmc    (1:nlay,1:ngptsw)   ,&!
              ztaucmc    (1:nlay,1:ngptsw)   ,&!
              zasycmc    (1:nlay,1:ngptsw)   ,&!
              zomgcmc    (1:nlay,1:ngptsw)   ,&!
              ztaormc    (1:nlay,1:ngptsw)   ,&!
              ztaua      (1:nlay,1:nbndsw)   ,&!
              zasya      (1:nlay,1:nbndsw)   ,&!
              zomga      (1:nlay,1:nbndsw)   ,&!
              cossza                         ,&!
              coldry     (1:nlay)            ,&!  
              wkl        (1:mxmol,1:nlay)    ,&! 
              adjflux    (1:jpband)          ,&!
              laytrop                        ,&!
              layswtch                       ,&!
              laylow                         ,&!
              jp         (1:nlay)            ,&! 
              jt         (1:nlay)            ,&! 
              jt1        (1:nlay)            ,&!
              co2mult    (1:nlay)            ,&!
              colch4     (1:nlay)            ,&!
              colco2     (1:nlay)            ,&!
              colh2o     (1:nlay)            ,&!
              colmol     (1:nlay)            ,&!
              coln2o     (1:nlay)            ,&!
              colo2      (1:nlay)            ,&!
              colo3      (1:nlay)            ,&!
              fac00      (1:nlay)            ,&!
              fac01      (1:nlay)            ,&!
              fac10      (1:nlay)            ,&!
              fac11      (1:nlay)            ,&!
              selffac    (1:nlay)            ,&!
              selffrac   (1:nlay)            ,&!
              indself    (1:nlay)            ,&!
              forfac     (1:nlay)            ,&!
              forfrac    (1:nlay)            ,&!
              indfor     (1:nlay)            ,&!
              zbbfd      (1:nlay+2)          ,&! zbbfd(nlay+2)          ! temporary downward shortwave flux (w/m2)
              zbbfu      (1:nlay+2)          ,&! zbbfu(nlay+2)          ! temporary upward shortwave flux (w/m2)
              zbbcd      (1:nlay+2)          ,&! zbbcd(nlay+2)          ! temporary clear sky downward shortwave flux (w/m2)
              zbbcu      (1:nlay+2)          ,&! zbbcu(nlay+2)          ! temporary clear sky upward shortwave flux (w/m2)
              zuvfd      (1:nlay+2)          ,&! zuvfd(nlay+2)          ! temporary UV downward shortwave flux (w/m2)
              zuvcd      (1:nlay+2)          ,&! zuvcd(nlay+2)          ! temporary clear sky UV downward shortwave flux (w/m2)
              znifd      (1:nlay+2)          ,&! znifd(nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
              znicd      (1:nlay+2)          ,&! znicd(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)
              znifu      (1:nlay+2)          ,&! znifu(nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
              znicu      (1:nlay+2)          ,&! znicu(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)
              zbbfddir   (1:nlay+2)          ,&! zbbfddir(nlay+2)       ! temporary downward direct shortwave flux (w/m2)
              zbbcddir   (1:nlay+2)          ,&! zbbcddir(nlay+2)       ! temporary clear sky downward direct shortwave flux (w/m2)
              zuvfddir   (1:nlay+2)          ,&! zuvfddir(nlay+2)       ! temporary UV downward direct shortwave flux (w/m2)
              zuvcddir   (1:nlay+2)          ,&! zuvcddir(nlay+2)       ! temporary clear sky UV downward direct shortwave flux (w/m2)
              znifddir   (1:nlay+2)          ,&! znifddir(nlay+2)       ! temporary near-IR downward direct shortwave flux (w/m2)
              znicddir   (1:nlay+2)          ,&! znicddir(nlay+2)       ! temporary clear sky near-IR downward direct shortwave flux (w/m2)
              zbbfsu     (1:nbndsw,1:nlay+2) ,&! zbbfsu(nbndsw,nlay+2)  ! temporary upward shortwave flux spectral (w/m2)
              zbbfsd     (1:nbndsw,1:nlay+2)  )! zbbfsd(nbndsw,nlay+2)  ! temporary downward shortwave flux spectral (w/m2)
         endif

! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top
!      real(kind=r8), intent(out) :: dirdnuv_localc(:,:)        ! Direct downward shortwave flux, UV/vis
!      real(kind=r8), intent(out) :: difdnuvc_local(:,:)        ! Diffuse downward shortwave flux, UV/vis
!      real(kind=r8), intent(out) :: dirdnirc_local(:,:)        ! Direct downward shortwave flux, near-IR
!      real(kind=r8), intent(out) :: difdnirc_local(:,:)        ! Diffuse downward shortwave flux, near-IR
 
!      real(kind=r8) :: zuvcd(nlay+2)          ! temporary clear sky UV downward shortwave flux (w/m2)
!      real(kind=r8) :: znifd(nlay+2)          ! temporary           near-IR downward shortwave flux (w/m2)
!      real(kind=r8) :: znicd(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)

         do i = 1, nlay+1
            swuflxc_local(iplon,i) = zbbcu(i)
            swdflxc_local(iplon,i) = zbbcd(i)
            swuflx_local(iplon,i) = zbbfu(i)
            swdflx_local(iplon,i) = zbbfd(i)
            swuflxs_local(:,iplon,i) = zbbfsu(:,i)
            swdflxs_local(:,iplon,i) = zbbfsd(:,i)
            uvdflx(i) = zuvfd(i)
            nidflx(i) = znifd(i)
!  Direct/diffuse fluxes
            dirdflux(i) = zbbfddir(i)
            difdflux(i) = swdflx_local(iplon,i) - dirdflux(i)
!  UV/visible direct/diffuse fluxes
            dirdnuvc_local(iplon,i) =zuvcddir(i) ! temporary clear sky UV downward direct shortwave flux (w/m2)
            dirdnuv_local(iplon,i) = zuvfddir(i) ! temporary           UV downward direct shortwave flux (w/m2)
            difdnuv_local(iplon,i) = zuvfd(i) - dirdnuv_local(iplon,i)
            difdnuvc_local(iplon,i)= zuvcd(i) - dirdnuvc_local(iplon,i)
!  Near-IR direct/diffuse fluxes
            dirdnirc_local(iplon,i) = znicddir(i)       ! temporary clear sky near-IR downward direct shortwave flux (w/m2)
            dirdnir_local(iplon,i)  = znifddir(i)       ! temporary           near-IR downward direct shortwave flux (w/m2)
            difdnir_local(iplon,i)  = znifd(i) - dirdnir_local(iplon,i)
            difdnirc_local(iplon,i) = znicd(i) - dirdnirc_local(iplon,i)
!  Added for net near-IR diagnostic
            ninflx_local(iplon,i) = znifd(i) - znifu(i)
            ninflxc_local(iplon,i) = znicd(i) - znicu(i)
         enddo

!  Total and clear sky net fluxes
         do i = 1, nlay+1
            swnflxc(i) = swdflxc_local(iplon,i) - swuflxc_local(iplon,i)
            swnflx(i) = swdflx_local(iplon,i) - swuflx_local(iplon,i)
         enddo

!  Total and clear sky heating rates
!  Heating units are in K/d. Flux units are in W/m2.
         do i = 1, nlay
            zdpgcp = heatfac / pdp(i)
            swhrc_local(iplon,i) = (swnflxc(i+1) - swnflxc(i)) * zdpgcp
            swhr_local(iplon,i) = (swnflx(i+1) - swnflx(i)) * zdpgcp
         enddo
         swhrc_local(iplon,nlay) = 0._r8
         swhr_local(iplon,nlay) = 0._r8
!***********************************************************************************************************************
         do i = 1, nlay+1
            swuflx  (iplon,i) = swuflx  (iplon,i) +  swuflx_local  (iplon,i) ! Total sky shortwave upward flux (W/m2)!    Dimensions: (ncol,nlay+1)
            swdflx  (iplon,i) = swdflx  (iplon,i) +  swdflx_local  (iplon,i) ! Total sky shortwave downward flux (W/m2)!    Dimensions: (ncol,nlay+1)
            swuflxc (iplon,i) = swuflxc (iplon,i) +  swuflxc_local (iplon,i) ! Clear sky shortwave upward flux (W/m2)!    Dimensions: (ncol,nlay+1)
            swdflxc (iplon,i) = swdflxc (iplon,i) +  swdflxc_local (iplon,i) ! Clear sky shortwave downward flux (W/m2)!    Dimensions: (ncol,nlay+1)

            dirdnuv (iplon,i) = dirdnuv (iplon,i) +  dirdnuv_local (iplon,i) ! Direct downward shortwave flux, UV/vis
            difdnuv (iplon,i) = difdnuv (iplon,i) +  difdnuv_local (iplon,i) ! Diffuse downward shortwave flux, UV/vis
            dirdnir (iplon,i) = dirdnir (iplon,i) +  dirdnir_local (iplon,i) ! Direct downward shortwave flux, near-IR
            difdnir (iplon,i) = difdnir (iplon,i) +  difdnir_local (iplon,i) ! Diffuse downward shortwave flux, near-IR

            dirdnuvc(iplon,i) = dirdnuvc(iplon,i) +  dirdnuvc_local(iplon,i) ! Direct downward shortwave flux, UV/vis
            difdnuvc(iplon,i) = difdnuvc(iplon,i) +  difdnuvc_local(iplon,i) ! Diffuse downward shortwave flux, UV/vis
            dirdnirc(iplon,i) = dirdnirc(iplon,i) +  dirdnirc_local(iplon,i) ! Direct downward shortwave flux, near-IR
            difdnirc(iplon,i) = difdnirc(iplon,i) +  difdnirc_local(iplon,i) ! Diffuse downward shortwave flux, near-IR

            ninflx  (iplon,i) = ninflx  (iplon,i) +  ninflx_local  (iplon,i) ! Net shortwave flux, near-IR
            ninflxc (iplon,i) = ninflxc (iplon,i) +  ninflxc_local (iplon,i) ! Net clear sky shortwave flux, near-IR
         END DO
         do i = 1, nlay+1
            DO k=1,nbndsw
               swuflxs (k,iplon,i) = swuflxs (k,iplon,i) + swuflxs_local (k,iplon,i)! shortwave spectral flux up
               swdflxs (k,iplon,i) = swdflxs (k,iplon,i) + swdflxs_local (k,iplon,i)! shortwave spectral flux down
            END DO
         END DO
         do i = 1, nlay
            swhrc   (iplon,i) = swhrc   (iplon,i) + swhrc_local   (iplon,i) ! Clear sky shortwave radiative heating rate (K/d)!    Dimensions: (ncol,nlay) 
            swhr    (iplon,i) = swhr    (iplon,i) + swhr_local    (iplon,i) ! Total sky shortwave radiative heating rate (K/d)!    Dimensions: (ncol,nlay)
         END DO
! End longitude loop
      enddo

      end subroutine rrtmg_sw

!*************************************************************************
      real(kind=r8) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use rrsw_con, only : pi

      integer, intent(in) :: idn

      real(kind=r8) :: gamma

      gamma = 2._r8*pi*(idn-1)/365._r8

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_r8 + .034221_r8 * cos(gamma) + .001289_r8 * sin(gamma) + &
                   .000719_r8 * cos(2._r8*gamma) + .000077_r8 * sin(2._r8*gamma)

      end function earth_sun

!***************************************************************************
      subroutine inatm_sw ( &        !  call inatm_sw (    &
            iplon      , &!              iplon                                     , &
            pcols      , &!              pcols
            nlay       , &!              nlay                                      , &
            icld       , &!              icld                                      , &
            iaer       , &!              iaer                                      , &
            play       , &!              play           (1:pcols,1:nlay)           , &
            plev       , &!              plev           (1:pcols,1:nlay+1)         , &
            tlay       , &!              tlay           (1:pcols,1:nlay)           , &
            tlev       , &!              tlev           (1:pcols,1:nlay+1)         , &
            tsfc       , &!              tsfc           (1:pcols)                  , &
            h2ovmr     , &!              h2ovmr         (1:pcols,1:nlay)           , &
            o3vmr      , &!              o3vmr          (1:pcols,1:nlay)           , &
            co2vmr     , &!              co2vmr         (1:pcols,1:nlay)           , &
            ch4vmr     , &!              ch4vmr         (1:pcols,1:nlay)           , &
            o2vmr      , &!              o2vmr          (1:pcols,1:nlay)           , &
            n2ovmr     , &!              n2ovmr         (1:pcols,1:nlay)           , &
            adjes      , &!              adjes                                     , &
            dyofyr     , &!              dyofyr                                    , &
            solvar     , &!              solvar         (1:nbndsw)                 , &
            inflgsw    , &!              inflgsw                                   , &
            iceflgsw   , &!              iceflgsw                                  , &
            liqflgsw   , &!              liqflgsw                                  , &
            cldfmcl    , &!              cldfmcl        (1:ngptsw,1:pcols,1:nlay-1), &
            taucmcl    , &!              taucmcl        (1:ngptsw,1:pcols,1:nlay-1), &
            ssacmcl    , &!              ssacmcl        (1:ngptsw,1:pcols,1:nlay-1), &
            asmcmcl    , &!              asmcmcl        (1:ngptsw,1:pcols,1:nlay-1), &
            fsfcmcl    , &!              fsfcmcl        (1:ngptsw,1:pcols,1:nlay-1), &
            ciwpmcl    , &!              ciwpmcl        (1:ngptsw,1:pcols,1:nlay-1), &
            clwpmcl    , &!              clwpmcl        (1:ngptsw,1:pcols,1:nlay-1), &
            reicmcl    , &!              reicmcl        (1:pcols,1:nlay-1)         , &
            relqmcl    , &!              relqmcl        (1:pcols,1:nlay-1)         , &
            tauaer     , &!              tauaer         (1:pcols,1:nlay-1,1:nbndsw), &
            ssaaer     , &!              ssaaer         (1:pcols,1:nlay-1,1:nbndsw), &
            asmaer     , &!              asmaer         (1:pcols,1:nlay-1,1:nbndsw), &
            pavel      , &!              pavel          (1:nlay)                   , &
            pz         , &!              pz             (0:nlay)                   , &
            pdp        , &!              pdp            (1:nlay)                   , &
            tavel      , &!              tavel          (1:nlay)                   , &
            tz         , &!              tz             (0:nlay)                   , &
            tbound     , &!              tbound                                    , &
            coldry     , &!              coldry         (1:nlay)                   , &
            wkl        , &!              wkl            (1:mxmol,1:nlay)           , &
            adjflux    , &!              adjflux           (1:jpband)              , &
            inflag     , &!              inflag                                    , &
            iceflag    , &!              iceflag                                   , &
            liqflag    , &!              liqflag                                   , &
            cldfmc     , &!              cldfmc           (1:ngptsw,1:nlay)        , &
            taucmc     , &!              taucmc           (1:ngptsw,1:nlay)        , &
            ssacmc     , &!              ssacmc           (1:ngptsw,1:nlay)        , &
            asmcmc     , &!              asmcmc           (1:ngptsw,1:nlay)        , &
            fsfcmc     , &!              fsfcmc           (1:ngptsw,1:nlay)        , &
            ciwpmc     , &!              ciwpmc           (1:ngptsw,1:nlay)        , &
            clwpmc     , &!              clwpmc           (1:ngptsw,1:nlay)        , &
            reicmc     , &!              reicmc           (1:nlay)                 , &
            dgesmc     , &!              dgesmc           (1:nlay)                 , &
            relqmc     , &!              relqmc           (1:nlay)                 , &
            taua       , &!              taua           (1:nlay,1:nbndsw)          , &
            ssaa       , &!              ssaa           (1:nlay,1:nbndsw)          , &
            asma         )!              asma           (1:nlay,1:nbndsw)            )
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_SW.
!  Set other RRTMG_SW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrsw, only : nbndsw, ngptsw, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_con, only : heatfac, oneminus, pi, grav, avogad
      use rrsw_wvn, only : ng, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: pcols
      integer, intent(in) :: nlay                       ! number of model layers
      integer, intent(in) :: icld                       ! clear/cloud and cloud overlap flag
      integer, intent(in) :: iaer                       ! aerosol option flag

      real(kind=r8), intent(in) :: play(1:pcols,1:nlay)            ! Layer pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(1:pcols,1:nlay+1)            ! Interface pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(1:pcols,1:nlay)             ! Layer temperatures (K)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(1:pcols,1:nlay+1)            ! Interface temperatures (K)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(1:pcols)              ! Surface temperature (K)
                                                        ! Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(1:pcols,1:nlay)           ! H2O volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(1:pcols,1:nlay)            ! O3 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(1:pcols,1:nlay)           ! CO2 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(1:pcols,1:nlay)           ! Methane volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(1:pcols,1:nlay)            ! O2 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(1:pcols,1:nlay)          ! Nitrous oxide volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)

      integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
      real(kind=r8), intent(in) :: solvar(jpb1:jpb2)    ! Solar constant (Wm-2) scaling per band

      integer, intent(in) :: inflgsw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflgsw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(1:ngptsw,1:pcols,1:nlay-1)       ! Cloud fraction
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl(1:ngptsw,1:pcols,1:nlay-1)       ! Cloud optical depth (optional)
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ssacmcl(1:ngptsw,1:pcols,1:nlay-1)       ! Cloud single scattering albedo
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: asmcmcl(1:ngptsw,1:pcols,1:nlay-1)       ! Cloud asymmetry parameter
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: fsfcmcl(1:ngptsw,1:pcols,1:nlay-1)       ! Cloud forward scattering fraction
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(1:ngptsw,1:pcols,1:nlay-1)       ! Cloud ice water path (g/m2)
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl(1:ngptsw,1:pcols,1:nlay-1)       ! Cloud liquid water path (g/m2)
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: reicmcl(1:pcols,1:nlay-1)         ! Cloud ice effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl(1:pcols,1:nlay-1)         ! Cloud water drop effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)

      real(kind=r8), intent(in) :: tauaer(1:pcols,1:nlay-1,1:nbndsw)        ! Aerosol optical depth
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=r8), intent(in) :: ssaaer(1:pcols,1:nlay-1,1:nbndsw)        ! Aerosol single scattering albedo
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=r8), intent(in) :: asmaer(1:pcols,1:nlay-1,1:nbndsw)        ! Aerosol asymmetry parameter
                                                        ! Dimensions: (ncol,nlay,nbndsw)

! Atmosphere

      real(kind=r8), intent(out) :: pavel(1:nlay)            ! layer pressures (mb) 
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: tavel (1:nlay)            ! layer temperatures (K)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: pz(0:nlay)              ! level (interface) pressures (hPa, mb)
                                                        ! Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tz(0:nlay)              ! level (interface) temperatures (K)
                                                        ! Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tbound              ! surface temperature (K)
      real(kind=r8), intent(out) :: pdp (1:nlay)              ! layer pressure thickness (hPa, mb)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: coldry(1:nlay)           ! dry air column density (mol/cm2)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: wkl(1:mxmol,1:nlay)            ! molecular amounts (mol/cm-2)
                                                        ! Dimensions: (mxmol,nlay)

      real(kind=r8), intent(out) :: adjflux(1:jpband)           ! adjustment for current Earth/Sun distance
                                                        ! Dimensions: (jpband)
!      real(kind=r8), intent(out) :: solvar(:)           ! solar constant scaling factor from rrtmg_sw
                                                        ! Dimensions: (jpband)
                                                        !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=r8), intent(out) :: taua(1:nlay,1:nbndsw)            ! Aerosol optical depth
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=r8), intent(out) :: ssaa(1:nlay,1:nbndsw)            ! Aerosol single scattering albedo
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=r8), intent(out) :: asma(1:nlay,1:nbndsw)            ! Aerosol asymmetry parameter
                                                        ! Dimensions: (nlay,nbndsw)

! Atmosphere/clouds - cldprop
      integer, intent(out) :: inflag                    ! flag for cloud property method
      integer, intent(out) :: iceflag                   ! flag for ice cloud properties
      integer, intent(out) :: liqflag                   ! flag for liquid cloud properties

      real(kind=r8), intent(out) :: cldfmc(1:ngptsw,1:nlay)         ! layer cloud fraction
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: taucmc(1:ngptsw,1:nlay)         ! cloud optical depth (non-delta scaled)
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: ssacmc(1:ngptsw,1:nlay)         ! cloud single scattering albedo (non-delta-scaled)
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: asmcmc(1:ngptsw,1:nlay)         ! cloud asymmetry parameter (non-delta scaled)
      real(kind=r8), intent(out) :: fsfcmc(1:ngptsw,1:nlay)         ! cloud forward scattering fraction (non-delta scaled)
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: ciwpmc(1:ngptsw,1:nlay)         ! cloud ice water path
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: clwpmc(1:ngptsw,1:nlay)         ! cloud liquid water path
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: reicmc(1:nlay)           ! cloud ice particle effective radius
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: dgesmc(1:nlay)           ! cloud ice particle effective radius
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: relqmc(1:nlay)           ! cloud liquid particle size
                                                        ! Dimensions: (nlay)

! ----- Local -----
      real(kind=r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
      real(kind=r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)
!      real(kind=r8), parameter :: amc = 44.0098_r8      ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=r8), parameter :: amo = 47.9998_r8      ! Molecular weight of ozone (g/mol)
!      real(kind=r8), parameter :: amo2 = 31.9999_r8     ! Molecular weight of oxygen (g/mol)
!      real(kind=r8), parameter :: amch4 = 16.0430_r8    ! Molecular weight of methane (g/mol)
!      real(kind=r8), parameter :: amn2o = 44.0128_r8    ! Molecular weight of nitrous oxide (g/mol)

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
      real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
      real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
      real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
      real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide

      real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)

      integer :: isp, l, ix, n, imol, ib, ig   ! Loop indices
      real(kind=r8) :: amm, summol                      ! 
      real(kind=r8) :: adjflx                           ! flux adjustment for Earth/Sun distance
!      real(kind=r8) :: earth_sun                        ! function for Earth/Sun distance adjustment
!      real(kind=r8) :: solar_band_irrad(jpb1:jpb2) ! rrtmg assumed-solar irradiance in each sw band

!  Initialize all molecular amounts to zero here, then pass input amounts
!  into RRTM array WKL below.

       wkl(:,:) = 0.0_r8
       cldfmc(:,:) = 0.0_r8
       taucmc(:,:) = 0.0_r8
       ssacmc(:,:) = 1.0_r8
       asmcmc(:,:) = 0.0_r8
       fsfcmc(:,:) = 0.0_r8
       ciwpmc(:,:) = 0.0_r8
       clwpmc(:,:) = 0.0_r8
       reicmc(:) = 0.0_r8
       dgesmc(:) = 0.0_r8
       relqmc(:) = 0.0_r8
       taua(:,:) = 0.0_r8
       ssaa(:,:) = 1.0_r8
       asma(:,:) = 0.0_r8
 
! Set flux adjustment for current Earth/Sun distance (two options).
! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
      adjflx = adjes
!
! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
!    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 
      if (dyofyr .gt. 0) then
         adjflx = earth_sun(dyofyr)
      endif

! Set incoming solar flux adjustment to include adjustment for
! current Earth/Sun distance (ADJFLX) and scaling of default internal
! solar constant (rrsw_scon = 1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
! to a single scaling factor as needed, or to a different value in each 
! band, which may be necessary for paleoclimate simulations. 
! 

      adjflux(:) = 0._r8
      do ib = jpb1,jpb2
         adjflux(ib) = adjflx * solvar(ib)
      enddo

!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_SW arrays for pressure, temperature,
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
         pdp(l) = pz(l-1) - pz(l)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,nlay-l+1)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,nlayers-l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,nlayers-l)/(1._r8 - h2o(iplon,nlayers-l)))*amdw
         wkl(2,l) = co2vmr(iplon,nlay-l+1)
         wkl(3,l) = o3vmr(iplon,nlay-l+1)
         wkl(4,l) = n2ovmr(iplon,nlay-l+1)
         wkl(6,l) = ch4vmr(iplon,nlay-l+1)
         wkl(7,l) = o2vmr(iplon,nlay-l+1) 
         amm = (1._r8 - wkl(1,l)) * amd + wkl(1,l) * amw            
         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_r8 * avogad / &
                     (1.e2_r8 * grav * amm * (1._r8 + wkl(1,l)))
      enddo

      coldry(nlay) = (pz(nlay-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(1,nlay-1)))

! At this point all molecular amounts in wkl are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  

      do l = 1, nlay
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
      enddo

! Transfer aerosol optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlay-1
            do ib = 1, nbndsw
               taua(l,ib) = tauaer(iplon,nlay-l,ib)
               ssaa(l,ib) = ssaaer(iplon,nlay-l,ib)
               asma(l,ib) = asmaer(iplon,nlay-l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflgsw
         iceflag = iceflgsw
         liqflag = liqflgsw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

         do l = 1, nlay-1
            do ig = 1, ngptsw
               cldfmc(ig,l) = cldfmcl(ig,iplon,nlay-l)
               taucmc(ig,l) = taucmcl(ig,iplon,nlay-l)
               ssacmc(ig,l) = ssacmcl(ig,iplon,nlay-l)
               asmcmc(ig,l) = asmcmcl(ig,iplon,nlay-l)
               fsfcmc(ig,l) = fsfcmcl(ig,iplon,nlay-l)
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
         ssacmc(:,nlay) = 1.0_r8
         asmcmc(:,nlay) = 0.0_r8
         fsfcmc(:,nlay) = 0.0_r8
         ciwpmc(:,nlay) = 0.0_r8
         clwpmc(:,nlay) = 0.0_r8
         reicmc(nlay) = 0.0_r8
         dgesmc(nlay) = 0.0_r8
         relqmc(nlay) = 0.0_r8
         taua(nlay,:) = 0.0_r8
         ssaa(nlay,:) = 1.0_r8
         asma(nlay,:) = 0.0_r8
     
      endif

      end subroutine inatm_sw
!tar begin
!  Vertical interpolation of Climate Aerosol (Kinne, 2013) (coarse mode)
!   4 subroutines:
!  
          subroutine zplev(pl,tl,top,n, zlev,dz)
!
!  written by Tarasova, Pisnitchenko, July 2015
!  Calculates heights and depth of the model layers 
!  using model level pressure and layer temperature  
!  INPUT:
!  pl...model level pressure in Pa from surface
!  tl...model layer temperature in K from surface
!  top...model topography in m 
!  n....model number of layers
!  OUTPUT:
!  zlev...model level heights in m from surface
!  dz...model layer depth (m) from surface
!
      Implicit none
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15)      
  INTEGER         , INTENT(IN   ) :: n
  REAL(KIND=r8)   , INTENT(IN   ) :: top
  REAL(KIND=r8)   , INTENT(IN   ) :: pl(n+1)  
  REAL(KIND=r8)   , INTENT(IN   ) :: tl(n)    
  REAL(KIND=r8)   , INTENT(OUT  ) :: zlev(n+1) 
  REAL(KIND=r8)   , INTENT(OUT  ) :: dz(n)     
! working  var  
  REAL(KIND=r8)   , PARAMETER ::  RG=287.04
  REAL(KIND=r8)   , PARAMETER ::  gt=9.80665
  INTEGER :: k
!
!  
  zlev(1)=top
  do k=2,n+1
   dz(k-1)= (RG*tl(k-1)/gt)*(log(pl(k-1))-log(pl(k)))
      zlev(k)=zlev(k-1)+ dz(k-1)
  enddo 
    end subroutine zplev      
!
!
!
    subroutine zaero(ep,m, f,h)  
!
! written by Tarasova, Pisnitchenko, July 2015
! Calculates Climate Aerosol extinction at new levels
! from Climate Aerosol extinction of (Kinne, 2013)
! INPUT:
! ep.....Climate Aerosol extinction in 40 layers from bottom
! m=40... Number of layers in Climate aerosol of Kinne, 2013
!
!  OUTPUT:
!  f......Climate Aerosol extinction at 42 levels from bottom
!  h......Climate Aerosol Height of 41 levels from bottom
! 
    Implicit none
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15)  
  REAL(KIND=r8)   , INTENT(IN   ) :: ep(m) 
  INTEGER         , INTENT(IN   ) :: m
  REAL(KIND=r8)   , INTENT(OUT   ) :: f(m+2) 
  REAL(KIND=r8)   , INTENT(OUT   ) :: h(m+1)     
! working var  
  INTEGER :: i 
  REAL(KIND=r8)   , PARAMETER ::  dd=500.0
!
   f(1)=ep(1)
  do i=1,m
   f(i+1)=ep(i)
  enddo
!   f(m+2)=ep(m)
    f(m+2)=0.0
!   
    h(1)=0.0
   do i=1,m
    h(i+1)=h(i)+dd
   enddo 
   end subroutine zaero
!
!
!     
         subroutine aeros_interp(zlev,h,f,m,n, fp)
!
!    written by Tarasova, Pisnitchenko, July 2015
! Calculates interpolation of aerosol extinction for the model layers
!   INPUT:
!   zlev....height (m) of model levels from surface
!   h......height (m) of Climate Aerosol at 41 levels from surface
!   f......Climate Aerosol extinction at 42 levels from surface
!   m......=40 number of Climate Aerosol layers (Kinne, 2013
!   n......number of model layers
!   OUTPUT:
!   fp.....aerosol extinction (1/m) in model layers from surface
!
       Implicit none      
!
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15)
  REAL(KIND=r8)   , INTENT(IN   ) :: zlev(n+1) 
  REAL(KIND=r8)   , INTENT(IN   ) :: h(m+1)
  REAL(KIND=r8)   , INTENT(IN   ) :: f(m+2) 
  INTEGER         , INTENT(IN   ) :: m,n  
  REAL(KIND=r8)   , INTENT(OUT  ) :: fp(n) 
! working var                   
  REAL(KIND=r8)   , PARAMETER ::  dd=500.0
  REAL(KIND=r8) ::dz1,dz2,di,fi1,fi2,fi
  INTEGER :: j,i,k,i1,i2,ik
!
  i=1
 do j=1,n
  If(zlev(j) > h(m+1)) then
  fp(j)=f(m+2)
   else
  do while(zlev(j) > h(i))
  i=i+1
   enddo
    i1=i
!    
   If(zlev(j+1) > h(m+1)) then
   dz1=(h(i2)-zlev(j))/dd
   dz2=(zlev(j+1)-h(m+1))/dd
   fi1=f(i2)*dz1
   fi=0.0
   do k=1,(m+1)-i2
    fi=fi+f(i2+k)
   enddo
   fi2=f(m+2)*dz2
   fi=fi+fi1+fi2
   di=dz1+dz2+(m+1-i2)*1.0
   fp(j)=fi/di
   else   
!       
  do while(zlev(j+1) > h(i))
    i=i+1
   enddo
   i2=i
!   
  dz1=h(i1)-zlev(j)
  dz2=h(i2)-zlev(j+1)
  dz2=dd-dz2
  dz1=dz1/dd
  dz2=dz2/dd
  fi1=f(i1)*dz1
    fi=0.0
    ik=0
   do k=1,i2-i1-1
    fi=fi+f(i1+k)
    ik=ik+1
   enddo
   fi2=f(i2)*dz2
   fi=fi+fi1+fi2
   di=dz1+dz2+ik*1.0
   fp(j)=fi/di
!
    endif
   endif 
  enddo
  end subroutine aeros_interp
!
!  
          subroutine  aod_n(fp,dz,n,od_n)
!
!    written by Tarasova, Pisnitchenko, July 2015
! Calculates total aerosol optical depth in all layers
! and normalized optical depth in each model layer from surface
!
!   INPUT:
!  fp...extinction (1/m) in each model layer from surface
!  dz...depth (m) of each model layer from surface
!  n...number of model layers 
!
!   OUTPUT:
!  od_n...normalized aerosol optical depth in each layer from surface 
!
         Implicit none
!         
!
 INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15)       
 REAL(KIND=r8)   , INTENT(IN   ) :: fp(n)  
 REAL(KIND=r8)   , INTENT(IN   ) :: dz(n) 
 REAL(KIND=r8)   , INTENT(OUT   ) :: od_n(n) 
 INTEGER   , INTENT(IN   )  :: n
! working var
  REAL(KIND=r8) :: a
  Integer :: i
!
!
    a=0.0
  do i=1,n
  a=a+fp(i)*dz(i)
   enddo
!      
     do i=1,n
   od_n(i)=(fp(i)*dz(i))/a
   enddo
   end subroutine aod_n

!tar end
      end module rrtmg_sw_rad



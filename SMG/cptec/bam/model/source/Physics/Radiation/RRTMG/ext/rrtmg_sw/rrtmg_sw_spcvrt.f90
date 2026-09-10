!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_spcvrt.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:14 $

      module rrtmg_sw_spcvrt

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use ppgrid,       only: pcols, begchunk, endchunk

!      use parkind, only : jpim, jprb
      use parrrsw, only : nbndsw, ngptsw, mxmol, jpband
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrspv, hnamspv
      use rrsw_wvn, only : ngc, ngs
      use rrtmg_sw_reftra, only: reftra_sw
      use rrtmg_sw_taumol, only: taumol_sw
      use rrtmg_sw_vrtqdr, only: vrtqdr_sw

      implicit none
!      INTEGER, PARAMETER :: pcols=1
!      INTEGER, PARAMETER :: begchunk=1
!      INTEGER, PARAMETER :: endchunk=1
!      real(kind=r8), allocatable, save :: ztaugst(:,:,:,:)      ! Gaseous optical depth storage array
!      real(kind=r8), allocatable, save :: ztaurst(:,:,:,:)      ! Rayleigh optical depth storage array
!      real(kind=r8), allocatable, save :: zsflxzenst(:,:,:)     ! Solar source function storage array

      contains

! ---------------------------------------------------------------------------
      subroutine spcvrt_sw(   & !           call spcvrt_sw (&
             nlayers        , & !             nlayers                                  , &! nlay                                ,&!
             istart         , & !             istart                                   , &! istart                              ,&!
             iend           , & !             iend                                     , &! iend                                ,&!
             icpr           , & !             icpr                                     , &! icpr                                ,&!
             idelm          , & !             idelm                                    , &! idelm                               ,&!
             iout           , & !             iout                                     , &! iout                                ,&!
             pavel          , & !             pavel         (1:nlayers)                , &! pavel         (1:nlay)              ,&!
             tavel          , & !             tavel         (1:nlayers)                , &! tavel         (1:nlay)              ,&!
             pz             , & !             pz            (0:nlayers)                , &! pz            (0:nlay)              ,&!
             tz             , & !             tz            (0:nlayers)                , &! tz            (0:nlay)              ,&!
             tbound         , & !             tbound                                   , &! tbound                              ,&!
             palbd          , & !             albdif        (1:nbndsw)                 , &! albdif        (1:nbndsw)            ,&!
             palbp          , & !             albdir        (1:nbndsw)                 , &! albdir        (1:nbndsw)            ,&!
             pclfr          , & !             cldfrac       (1:nlayers)                , &! zcldfmc       (1:nlay,1:ngptsw)     ,&!
             ptauc          , & !             ztauc         (1:nlayers,1:nbndsw)       , &! ztaucmc       (1:nlay,1:ngptsw)     ,&!
             pasyc          , & !             zasyc         (1:nlayers,1:nbndsw)       , &! zasycmc       (1:nlay,1:ngptsw)     ,&!
             pomgc          , & !             zomgc         (1:nlayers,1:nbndsw)       , &! zomgcmc       (1:nlay,1:ngptsw)     ,&!
             ptaucorig      , & !             ztaucorig     (1:nlayers,1:nbndsw)       , &! ztaormc       (1:nlay,1:ngptsw)     ,&!
             ptaua          , & !             ztaua         (1:nlayers,1:nbndsw)       , &! ztaua         (1:nlay,1:nbndsw)     ,&!
             pasya          , & !             zasya         (1:nlayers,1:nbndsw)       , &! zasya         (1:nlay,1:nbndsw)     ,&!
             pomga          , & !             zomga         (1:nlayers,1:nbndsw)       , &! zomga         (1:nlay,1:nbndsw)     ,&!
             prmu0          , & !             cossza                                   , &! cossza                              ,&!
             coldry         , & !             coldry        (1:nlayers)                , &! coldry        (1:nlay)              ,&!
             wkl            , & !             wkl           (1:mxmol,1:nlayers)        , &! wkl           (1:mxmol,1:nlay)      ,&!
             adjflux        , & !             adjflux       (1:jpband)                 , &! adjflux       (1:jpband)            ,&!
             laytrop        , & !             laytrop                                  , &! laytrop                             ,&!
             layswtch       , & !             layswtch                                 , &! layswtch                            ,&!
             laylow         , & !             laylow                                   , &! laylow                              ,&!
             jp             , & !             jp            (1:nlayers)                , &! jp            (1:nlay)              ,&!
             jt             , & !             jt            (1:nlayers)                , &! jt            (1:nlay)              ,&!
             jt1            , & !             jt1           (1:nlayers)                , &! jt1           (1:nlay)              ,&!
             co2mult        , & !             co2mult       (1:nlayers)                , &! co2mult       (1:nlay)              ,&!
             colch4         , & !             colch4        (1:nlayers)                , &! colch4        (1:nlay)               ,&!
             colco2         , & !             colco2        (1:nlayers)                , &! colco2        (1:nlay)               ,&!
             colh2o         , & !             colh2o        (1:nlayers)                , &! colh2o        (1:nlay)               ,&!
             colmol         , & !             colmol        (1:nlayers)                , &! colmol        (1:nlay)               ,&!
             coln2o         , & !             coln2o        (1:nlayers)                , &! coln2o        (1:nlay)               ,&!
             colo2          , & !             colo2         (1:nlayers)                , &! colo2         (1:nlay)                ,&!
             colo3          , & !             colo3         (1:nlayers)                , &! colo3         (1:nlay)                ,&!
             fac00          , & !             fac00         (1:nlayers)                , &! fac00         (1:nlay)                ,&!
             fac01          , & !             fac01         (1:nlayers)                , &! fac01         (1:nlay)                ,&!
             fac10          , & !             fac10         (1:nlayers)                , &! fac10         (1:nlay)                ,&!
             fac11          , & !             fac11         (1:nlayers)                , &! fac11         (1:nlay)                ,&!
             selffac        , & !             selffac       (1:nlayers)                , &! selffac       (1:nlay)                ,&!
             selffrac       , & !             selffrac      (1:nlayers)                , &! selffrac      (1:nlay)                ,&!
             indself        , & !             indself       (1:nlayers)                , &! indself       (1:nlay)                ,&!
             forfac         , & !             forfac        (1:nlayers)                , &! forfac        (1:nlay)                ,&!
             forfrac        , & !             forfrac       (1:nlayers)                , &! forfrac       (1:nlay)                ,&!
             indfor         , & !             indfor        (1:nlayers)                , &! indfor        (1:nlay)                ,&!
             pbbfd          , & !             zbbfd         (1:nlayers+2)              , &! zbbfd         (1:nlay+2)              ,&!
             pbbfu          , & !             zbbfu         (1:nlayers+2)              , &! zbbfu         (1:nlay+2)              ,&!
             pbbcd          , & !             zbbcd         (1:nlayers+2)              , &! zbbcd         (1:nlay+2)              ,&!
             pbbcu          , & !             zbbcu         (1:nlayers+2)              , &! zbbcu         (1:nlay+2)              ,&!
             puvfd          , & !             zuvfd         (1:nlayers+2)              , &! zuvfd         (1:nlay+2)              ,&!
             puvcd          , & !             zuvcd         (1:nlayers+2)              , &! zuvcd         (1:nlay+2)              ,&!
             pnifd          , & !             znifd         (1:nlayers+2)              , &! znifd         (1:nlay+2)              ,&!
             pnicd          , & !             znicd         (1:nlayers+2)              , &! znicd         (1:nlay+2)              ,&!
             pnifu          , & !             znifu         (1:nlayers+2)              , &! znifu         (1:nlay+2)              ,&!
             pnicu          , & !             znicu         (1:nlayers+2)              , &! znicu         (1:nlay+2)              ,&!
             pbbfddir       , & !             zbbfddir      (1:nlayers+2)              , &! zbbfddir      (1:nlay+2)              ,&!
             pbbcddir       , & !             zbbcddir      (1:nlayers+2)              , &! zbbcddir      (1:nlay+2)              ,&!
             puvfddir       , & !             zuvfddir      (1:nlayers+2)              , &! zuvfddir      (1:nlay+2)              ,&!
             puvcddir       , & !             zuvcddir      (1:nlayers+2)              , &! zuvcddir      (1:nlay+2)              ,&!
             pnifddir       , & !             znifddir      (1:nlayers+2)              , &! znifddir      (1:nlay+2)              ,&!
             pnicddir         ) !            znicddir       (1:nlayers+2)                )! znicddir      (1:nlay+2)              ,&!
! ---------------------------------------------------------------------------
!
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker. 
!
! Interface:  *spcvrt_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:
!
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 
!
! ------------------------------------------------------------------

! ------- Declarations ------

! -------- Input -------

!      integer, intent(in) :: lchnk
!      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: nlayers
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      integer, intent(in) :: icpr
      integer, intent(in) :: idelm     ! delta-m scaling flag
                                       ! [0 = direct and diffuse fluxes are unscaled]
                                       ! [1 = direct and diffuse fluxes are scaled]
      integer, intent(in) :: iout

      integer, intent(in) :: laytrop
      integer, intent(in) :: layswtch
      integer, intent(in) :: laylow

      integer, intent(in) :: indfor (1:nlayers)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: indself (1:nlayers)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: jp (1:nlayers)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: jt (1:nlayers)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: jt1 (1:nlayers)
                                                                 !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: pavel(1:nlayers)                    ! layer pressure (hPa, mb) 
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: tavel(1:nlayers)                   ! layer temperature (K)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: pz(0:nlayers)                       ! level (interface) pressure (hPa, mb)
                                                                 !   Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: tz(0:nlayers)                       ! level temperatures (hPa, mb)
                                                                 !   Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: tbound                      ! surface temperature (K)
      real(kind=r8), intent(in) :: wkl(1:mxmol,1:nlayers)                    ! molecular amounts (mol/cm2) 
                                                                 !   Dimensions: (mxmol,nlayers)
      real(kind=r8), intent(in) :: coldry(1:nlayers)                    ! dry air column density (mol/cm2)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colmol(1:nlayers) 
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: adjflux(1:jpband)                   ! Earth/Sun distance adjustment
                                                                 !   Dimensions: (jpband)

      real(kind=r8), intent(in) :: palbd(1:nbndsw)                    ! surface albedo (diffuse)
                                                                 !   Dimensions: (nbndsw)
      real(kind=r8), intent(in) :: palbp(1:nbndsw)                    ! surface albedo (direct)
                                                                 !   Dimensions: (nbndsw)
      real(kind=r8), intent(in) :: prmu0                       ! cosine of solar zenith angle
      real(kind=r8), intent(in) :: pclfr(1:nlayers)                    ! cloud fraction
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: ptauc(1:nlayers,1:nbndsw)                   ! cloud optical depth
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=r8), intent(in) :: pasyc(1:nlayers,1:nbndsw)                  ! cloud asymmetry parameter
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=r8), intent(in) :: pomgc(1:nlayers,1:nbndsw)                   ! cloud single scattering albedo
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=r8), intent(in) :: ptaucorig(1:nlayers,1:nbndsw)               ! cloud optical depth, non-delta scaled
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=r8), intent(in) :: ptaua(1:nlayers,1:nbndsw)                   ! aerosol optical depth
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=r8), intent(in) :: pasya(1:nlayers,1:nbndsw)                  ! aerosol asymmetry parameter
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=r8), intent(in) :: pomga(1:nlayers,1:nbndsw)                   ! aerosol single scattering albedo
                                                                 !   Dimensions: (nlayers,nbndsw)

      real(kind=r8), intent(in) :: colh2o(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colco2(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colch4(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: co2mult(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colo3(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colo2(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: coln2o(1:nlayers)
                                                                 !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: forfac(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: forfrac(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: selffac(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: selffrac(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac00(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac01(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac10(1:nlayers)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac11(1:nlayers)
                                                                 !   Dimensions: (nlayers)

! ------- Output -------
                                                                 !   All Dimensions: (nlayers+1)
      real(kind=r8), intent(out) :: pbbcd(1:nlayers+2)
      real(kind=r8), intent(out) :: pbbcu(1:nlayers+2)
      real(kind=r8), intent(out) :: pbbfd(1:nlayers+2)
      real(kind=r8), intent(out) :: pbbfu(1:nlayers+2)
      real(kind=r8), intent(out) :: pbbfddir(1:nlayers+2)
      real(kind=r8), intent(out) :: pbbcddir(1:nlayers+2)

      real(kind=r8), intent(out) :: puvcd(1:nlayers+2)
      real(kind=r8), intent(out) :: puvfd(1:nlayers+2)
      real(kind=r8), intent(out) :: puvcddir(1:nlayers+2)
      real(kind=r8), intent(out) :: puvfddir(1:nlayers+2)

      real(kind=r8), intent(out) :: pnicd(1:nlayers+2)
      real(kind=r8), intent(out) :: pnifd(1:nlayers+2)
      real(kind=r8), intent(out) :: pnicddir(1:nlayers+2)
      real(kind=r8), intent(out) :: pnifddir(1:nlayers+2)

! Added for net near-IR flux diagnostic
      real(kind=r8), intent(out) :: pnicu(1:nlayers+2)
      real(kind=r8), intent(out) :: pnifu(1:nlayers+2)

! Output - inactive                                              !   All Dimensions: (nlayers+1)
!      real(kind=r8), intent(out) :: puvcu(:)
!      real(kind=r8), intent(out) :: puvfu(:)
!      real(kind=r8), intent(out) :: pnicu(:)
!      real(kind=r8), intent(out) :: pnifu(:)
!      real(kind=r8), intent(out) :: pvscd(:)
!      real(kind=r8), intent(out) :: pvscu(:)
!      real(kind=r8), intent(out) :: pvsfd(:)
!      real(kind=r8), intent(out) :: pvsfu(:)


! ------- Local -------

      logical :: lrtchkclr(nlayers),lrtchkcld(nlayers)

      integer  :: klev
      integer :: ib1, ib2, ibm, igt, ikl, ikp, ikx
      integer :: iw, jb, jg, jl, jk
!      integer, parameter :: nuv = ?? 
!      integer, parameter :: nvs = ?? 
      integer :: itind

      real(kind=r8) :: tblind, ze1
      real(kind=r8) :: zclear, zcloud
      real(kind=r8) :: zdbt(nlayers+1), zdbt_nodel(nlayers+1)
      real(kind=r8) :: zgc(nlayers), zgcc(nlayers), zgco(nlayers)
      real(kind=r8) :: zomc(nlayers), zomcc(nlayers), zomco(nlayers)
      real(kind=r8) :: zrdnd(nlayers+1), zrdndc(nlayers+1)
      real(kind=r8) :: zref(nlayers+1), zrefc(nlayers+1), zrefo(nlayers+1)
      real(kind=r8) :: zrefd(nlayers+1), zrefdc(nlayers+1), zrefdo(nlayers+1)
      real(kind=r8) :: zrup(nlayers+1), zrupd(nlayers+1)
      real(kind=r8) :: zrupc(nlayers+1), zrupdc(nlayers+1)
      real(kind=r8) :: zs1(nlayers+1)
      real(kind=r8) :: ztauc(nlayers), ztauo(nlayers)
      real(kind=r8) :: ztdn(nlayers+1), ztdnd(nlayers+1), ztdbt(nlayers+1)
      real(kind=r8) :: ztoc(nlayers), ztor(nlayers)
      real(kind=r8) :: ztra(nlayers+1), ztrac(nlayers+1), ztrao(nlayers+1)
      real(kind=r8) :: ztrad(nlayers+1), ztradc(nlayers+1), ztrado(nlayers+1)
      real(kind=r8) :: zdbtc(nlayers+1), ztdbtc(nlayers+1)
      real(kind=r8) :: zincflx(ngptsw), zdbtc_nodel(nlayers+1) 
      real(kind=r8) :: ztdbt_nodel(nlayers+1), ztdbtc_nodel(nlayers+1)

      real(kind=r8) :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real(kind=r8) :: zwf, tauorig, repclc
!     real(kind=r8) :: zincflux                                   ! inactive

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=r8) :: ztaug(nlayers,16), ztaur(nlayers,16)
!      real(kind=r8) :: zsflxzen(16)
      real(kind=r8) :: ztaug(nlayers,ngptsw), ztaur(nlayers,ngptsw)
      real(kind=r8) :: zsflxzen(ngptsw)

! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=r8) :: zcd(nlayers+1,ngptsw), zcu(nlayers+1,ngptsw)
      real(kind=r8) :: zfd(nlayers+1,ngptsw), zfu(nlayers+1,ngptsw)

! Inactive arrays
!     real(kind=r8) :: zbbcd(nlayers+1), zbbcu(nlayers+1)
!     real(kind=r8) :: zbbfd(nlayers+1), zbbfu(nlayers+1)
!     real(kind=r8) :: zbbfddir(nlayers+1), zbbcddir(nlayers+1)

! Local storage arrays for zsflxzen, ztaug, and ztaur 

!PK      real(kind=r8), allocatable, save :: ztaugst(:,:,:,:)      ! Gaseous optical depth storage array
!PK      real(kind=r8), allocatable, save :: ztaurst(:,:,:,:)      ! Rayleigh optical depth storage array
!PK      real(kind=r8), allocatable, save :: zsflxzenst(:,:,:)     ! Solar source function storage array

! Allocate storage arrays
!      if (.not.allocated(ztaugst)) allocate (ztaugst(pcols,nlayers,ngptsw,begchunk:endchunk))
!      if (.not.allocated(ztaurst)) allocate (ztaurst(pcols,nlayers,ngptsw,begchunk:endchunk))
!      if (.not.allocated(zsflxzenst)) allocate (zsflxzenst(pcols,ngptsw,begchunk:endchunk))
       
! ------------------------------------------------------------------

! Initializations

      ib1 = istart
      ib2 = iend
      klev = nlayers
      iw = 0
      repclc = 1.e-12_r8
!      zincflux = 0.0_r8

      do jk=1,klev+1
         pbbcd(jk)=0._r8
         pbbcu(jk)=0._r8
         pbbfd(jk)=0._r8
         pbbfu(jk)=0._r8
         pbbcddir(jk)=0._r8
         pbbfddir(jk)=0._r8
         puvcd(jk)=0._r8
         puvfd(jk)=0._r8
         puvcddir(jk)=0._r8
         puvfddir(jk)=0._r8
         pnicd(jk)=0._r8
         pnifd(jk)=0._r8
         pnicddir(jk)=0._r8
         pnifddir(jk)=0._r8
         pnicu(jk)=0._r8
         pnifu(jk)=0._r8
      enddo


! Perform calculations of optical depth and solar source function at interval
! specified by dotau and store output arrays zsflxzen, ztaug, and ztaur in memory.
! In intervening time steps, obtain zsflxzen, ztaug, and ztaur from stored arrays. 
! Memory storage should be replaced with I/O to netCDF file.

!PK         if (dotau) then
 
! Calculate the optical depths for gaseous absorption and Rayleigh scattering

            call taumol_sw(&
                     nlayers                       , &
                     ngptsw                        , &
                     colh2o    (1:nlayers)         , &
                     colco2    (1:nlayers)         , &
                     colch4    (1:nlayers)         , &
                     colo2     (1:nlayers)         , &
                     colo3     (1:nlayers)         , &
                     colmol    (1:nlayers)         , &
                     laytrop                       , &
                     jp        (1:nlayers)         , &
                     jt        (1:nlayers)         , &
                     jt1       (1:nlayers)         , &
                     fac00     (1:nlayers)         , &
                     fac01     (1:nlayers)         , &
                     fac10     (1:nlayers)         , &
                     fac11     (1:nlayers)         , &
                     selffac   (1:nlayers)         , &
                     selffrac  (1:nlayers)         , &
                     indself   (1:nlayers)         , &
                     forfac    (1:nlayers)         , &
                     forfrac   (1:nlayers)         , &
                     indfor    (1:nlayers)         , &
                     zsflxzen  (1:ngptsw)          , &
                     ztaug     (1:nlayers,1:ngptsw), &
                     ztaur     (1:nlayers,1:ngptsw))

! Store ztaug, ztaur, and zsflxzen for use when dotau is false
!PK            ztaugst(iplon,:nlayers,:ngptsw,lchnk) = ztaug(:nlayers,:ngptsw)
!PK            ztaurst(iplon,:nlayers,:ngptsw,lchnk) = ztaur(:nlayers,:ngptsw)
!PK            zsflxzenst(iplon,:ngptsw,lchnk) = zsflxzen(:ngptsw)

!PK         else

! Restore ztaug, ztaur, and zsflxzen from storage when dotau is false
!PK            ztaug(:nlayers,:ngptsw) = ztaugst(iplon,:nlayers,:ngptsw,lchnk)
!PK            ztaur(:nlayers,:ngptsw) = ztaurst(iplon,:nlayers,:ngptsw,lchnk)
!PK            zsflxzen(:ngptsw) = zsflxzenst(iplon,:ngptsw,lchnk)

!PK         endif

! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

      jb = ib1-1                  ! ???
      do jb = ib1, ib2
         ibm = jb-15
         igt = ngc(ibm)

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.ibm.ge.2) iw = ngs(ibm-1)

!        do jk=1,klev+1
!           zbbcd(jk)=0.0_r8
!           zbbcu(jk)=0.0_r8
!           zbbfd(jk)=0.0_r8
!           zbbfu(jk)=0.0_r8
!        enddo

! Top of g-point interval loop within each band (iw is cumulative counter) 
         do jg = 1,igt
            iw = iw+1

! Apply adjustments for correct Earth/Sun distance and zenith angle to incoming solar flux
            zincflx(iw) = adjflux(jb) * zsflxzen(iw) * prmu0
!             zincflux = zincflux + adjflux(jb) * zsflxzen(iw) * prmu0           ! inactive

! Compute layer reflectances and transmittances for direct and diffuse sources, 
! first clear then cloudy

! zrefc(jk)  direct albedo for clear
! zrefo(jk)  direct albedo for cloud
! zrefdc(jk) diffuse albedo for clear
! zrefdo(jk) diffuse albedo for cloud
! ztrac(jk)  direct transmittance for clear
! ztrao(jk)  direct transmittance for cloudy
! ztradc(jk) diffuse transmittance for clear
! ztrado(jk) diffuse transmittance for cloudy
!  
! zref(jk)   direct reflectance
! zrefd(jk)  diffuse reflectance
! ztra(jk)   direct transmittance
! ztrad(jk)  diffuse transmittance
!
! zdbtc(jk)  clear direct beam transmittance
! zdbto(jk)  cloudy direct beam transmittance
! zdbt(jk)   layer mean direct beam transmittance
! ztdbt(jk)  total direct beam transmittance at levels

! Clear-sky    
!   TOA direct beam    
            ztdbtc(1)=1.0_r8
            ztdbtc_nodel(1)=1.0_r8
!   Surface values
            zdbtc(klev+1) =0.0_r8
            ztrac(klev+1) =0.0_r8
            ztradc(klev+1)=0.0_r8
            zrefc(klev+1) =palbp(ibm)
            zrefdc(klev+1)=palbd(ibm)
            zrupc(klev+1) =palbp(ibm)
            zrupdc(klev+1)=palbd(ibm)
           
! Total sky    
!   TOA direct beam    
            ztdbt(1)=1.0_r8
            ztdbt_nodel(1)=1.0_r8
!   Surface values
            zdbt(klev+1) =0.0_r8
            ztra(klev+1) =0.0_r8
            ztrad(klev+1)=0.0_r8
            zref(klev+1) =palbp(ibm)
            zrefd(klev+1)=palbd(ibm)
            zrup(klev+1) =palbp(ibm)
            zrupd(klev+1)=palbd(ibm)
    
    
! Top of layer loop
            do jk=1,klev

! Note: two-stream calculations proceed from top to bottom; 
!   RRTMG_SW quantities are given bottom to top and are reversed here

               ikl=klev+1-jk

! Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
               lrtchkclr(jk)=.true.

!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
               lrtchkcld(jk)=.false.
               lrtchkcld(jk)=(pclfr(ikl) > repclc)

! Clear-sky optical parameters - this section inactive     
!   Original
!               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw)
!               zomcc(jk) = ztaur(ikl,iw) / ztauc(jk)
!               zgcc(jk) = 0.0001_r8
!   Total sky optical parameters        
!               ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptauc(ikl,ibm)
!               zomco(jk) = ptauc(ikl,ibm) * pomgc(ikl,ibm) + ztaur(ikl,iw)
!               zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
!                           ztaur(ikl,iw) * 0.0001_r8) / zomco(jk)
!               zomco(jk) = zomco(jk) / ztauo(jk)

! Clear-sky optical parameters including aerosols
               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm)
               zomcc(jk) = ztaur(ikl,iw) * 1.0_r8 + ptaua(ikl,ibm) * pomga(ikl,ibm)
               zgcc(jk) = pasya(ikl,ibm) * pomga(ikl,ibm) * ptaua(ikl,ibm) / zomcc(jk)
               zomcc(jk) = zomcc(jk) / ztauc(jk)

! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
!   \/\/\/ This block of code is only needed for unscaled direct beam calculation
               if (idelm .eq. 0) then
!     
                  zclear = 1.0_r8 - pclfr(ikl)
                  zcloud = pclfr(ikl)

! Clear
!                   zdbtmc = exp(-ztauc(jk) / prmu0)
 
! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ze1 = ztauc(jk) / prmu0
                  if (ze1 .le. od_lo) then
                     zdbtmc = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
                  else 
                     tblind = ze1 / (bpade + ze1)
                     itind = tblint * tblind + 0.5_r8
                     zdbtmc = exp_tbl(itind)
                  endif

                  zdbtc_nodel(jk) = zdbtmc
                  ztdbtc_nodel(jk+1) = zdbtc_nodel(jk) * ztdbtc_nodel(jk)

! Clear + Cloud
                  tauorig = ztauc(jk) + ptaucorig(ikl,ibm)
!                   zdbtmo = exp(-tauorig / prmu0)

! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ze1 = tauorig / prmu0
                  if (ze1 .le. od_lo) then
                     zdbtmo = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
                  else
                     tblind = ze1 / (bpade + ze1)
                     itind = tblint * tblind + 0.5_r8
                     zdbtmo = exp_tbl(itind)
                  endif

                  zdbt_nodel(jk) = zclear * zdbtmc + zcloud * zdbtmo
                  ztdbt_nodel(jk+1) = zdbt_nodel(jk) * ztdbt_nodel(jk)

               endif
!   /\/\/\ Above code only needed for unscaled direct beam calculation


! Delta scaling - clear   
               zf = zgcc(jk) * zgcc(jk)
               zwf = zomcc(jk) * zf
               ztauc(jk) = (1.0_r8 - zwf) * ztauc(jk)
               zomcc(jk) = (zomcc(jk) - zwf) / (1.0_r8 - zwf)
               zgcc (jk) = (zgcc(jk) - zf) / (1.0_r8 - zf)

! Total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
               if (icpr .ge. 1) then
                  ztauo(jk) = ztauc(jk) + ptauc(ikl,ibm)
                  zomco(jk) = ztauc(jk) * zomcc(jk) + ptauc(ikl,ibm) * pomgc(ikl,ibm) 
                  zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
                              ztauc(jk) * zomcc(jk) * zgcc(jk)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
               elseif (icpr .eq. 0) then
                  ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm) + ptauc(ikl,ibm)
                  zomco(jk) = ptaua(ikl,ibm) * pomga(ikl,ibm) + ptauc(ikl,ibm) * pomgc(ikl,ibm) + &
                              ztaur(ikl,iw) * 1.0_r8
                  zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
                              ptaua(ikl,ibm)*pomga(ikl,ibm)*pasya(ikl,ibm)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Delta scaling - clouds 
!   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
                  zf = zgco(jk) * zgco(jk)
                  zwf = zomco(jk) * zf
                  ztauo(jk) = (1._r8 - zwf) * ztauo(jk)
                  zomco(jk) = (zomco(jk) - zwf) / (1.0_r8 - zwf)
                  zgco (jk) = (zgco(jk) - zf) / (1.0_r8 - zf)
               endif 

! End of layer loop
            enddo    


! Clear sky reflectivities
            call reftra_sw ( &
                            nlayers               , &
                            lrtchkclr(1:nlayers)  , &
                            zgcc     (1:nlayers)  , &
                            prmu0                 , &
                            ztauc    (1:nlayers)  , &
                            zomcc    (1:nlayers)  , &
                            zrefc    (1:nlayers+1), &
                            zrefdc   (1:nlayers+1), &
                            ztrac    (1:nlayers+1), &
                            ztradc   (1:nlayers+1)  )

! Total sky reflectivities      
            call reftra_sw ( &
                            nlayers               , &
                            lrtchkcld(1:nlayers)  , &
                            zgco     (1:nlayers)  , &
                            prmu0                 , &
                            ztauo    (1:nlayers)  , &
                            zomco    (1:nlayers)  , &
                            zrefo    (1:nlayers+1), &
                            zrefdo   (1:nlayers+1), &
                            ztrao    (1:nlayers+1), &
                            ztrado   (1:nlayers+1)  )


            do jk=1,klev

! Combine clear and cloudy contributions for total sky
               ikl = klev+1-jk 
               zclear = 1.0_r8 - pclfr(ikl)
               zcloud = pclfr(ikl)

               zref(jk) = zclear*zrefc(jk) + zcloud*zrefo(jk)
               zrefd(jk)= zclear*zrefdc(jk) + zcloud*zrefdo(jk)
               ztra(jk) = zclear*ztrac(jk) + zcloud*ztrao(jk)
               ztrad(jk)= zclear*ztradc(jk) + zcloud*ztrado(jk)

! Direct beam transmittance        

! Clear
!                zdbtmc = exp(-ztauc(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauc(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmc = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  zdbtmc = exp_tbl(itind)
               endif

               zdbtc(jk) = zdbtmc
               ztdbtc(jk+1) = zdbtc(jk)*ztdbtc(jk)

! Clear + Cloud
!                zdbtmo = exp(-ztauo(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauo(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmo = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  zdbtmo = exp_tbl(itind)
               endif

               zdbt(jk) = zclear*zdbtmc + zcloud*zdbtmo
               ztdbt(jk+1) = zdbt(jk)*ztdbt(jk)
        
            enddo           
                 
! Vertical quadrature for clear-sky fluxes

            call vrtqdr_sw ( &
                            ngptsw                        , &
                            nlayers                       , &
                            iw                            , &
                            zrefc   (1:nlayers+1)         , &
                            zrefdc  (1:nlayers+1)         , &
                            ztrac   (1:nlayers+1)         , &
                            ztradc  (1:nlayers+1)         , &
                            zdbtc   (1:nlayers+1)         , &
                            zrdndc  (1:nlayers+1)         , &
                            zrupc   (1:nlayers+1)         , &
                            zrupdc  (1:nlayers+1)         , &
                            ztdbtc  (1:nlayers+1)         , &
                            zcd     (1:nlayers+1,1:ngptsw), &
                            zcu     (1:nlayers+1,1:ngptsw)  )
      
! Vertical quadrature for cloudy fluxes

            call vrtqdr_sw ( &
                            ngptsw                        , &
                            nlayers                       , &
                            iw                            , &
                            zref    (1:nlayers+1)         , &
                            zrefd   (1:nlayers+1)         , &
                            ztra    (1:nlayers+1)         , &
                            ztrad   (1:nlayers+1)         , &
                            zdbt    (1:nlayers+1)         , &
                            zrdnd   (1:nlayers+1)         , &
                            zrup    (1:nlayers+1)         , &
                            zrupd   (1:nlayers+1)         , &
                            ztdbt   (1:nlayers+1)         , &
                            zfd     (1:nlayers+1,1:ngptsw), &
                            zfu     (1:nlayers+1,1:ngptsw)  )

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays

            do jk=1,klev+1
               ikl=klev+2-jk

! Accumulate spectral fluxes over bands - inactive
!               zbbfu(ikl) = zbbfu(ikl) + zincflx(iw)*zfu(jk,iw)  
!               zbbfd(ikl) = zbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
!               zbbcu(ikl) = zbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
!               zbbcd(ikl) = zbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
!               zbbfddir(ikl) = zbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!               zbbcddir(ikl) = zbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)

! Accumulate spectral fluxes over whole spectrum  
               pbbfu(ikl) = pbbfu(ikl) + zincflx(iw)*zfu(jk,iw)
               pbbfd(ikl) = pbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
               pbbcu(ikl) = pbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
               pbbcd(ikl) = pbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
               if (idelm .eq. 0) then
                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
               elseif (idelm .eq. 1) then
                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt(jk)
                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc(jk)
               endif

! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  puvcd(ikl) = puvcd(ikl) + zincflx(iw)*zcd(jk,iw)
                  puvfd(ikl) = puvfd(ikl) + zincflx(iw)*zfd(jk,iw)
                  if (idelm .eq. 0) then
                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  elseif (idelm .eq. 1) then
                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt(jk)
                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc(jk)
                  endif
! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                  pnicd(ikl) = pnicd(ikl) + zincflx(iw)*zcd(jk,iw)
                  pnifd(ikl) = pnifd(ikl) + zincflx(iw)*zfd(jk,iw)
                  if (idelm .eq. 0) then
                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  elseif (idelm .eq. 1) then
                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt(jk)
                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc(jk)
                  endif
! Added for net near-IR flux diagnostic 
                  pnicu(ikl) = pnicu(ikl) + zincflx(iw)*zcu(jk,iw)
                  pnifu(ikl) = pnifu(ikl) + zincflx(iw)*zfu(jk,iw)
               endif

            enddo

! End loop on jg, g-point interval
         enddo             

! End loop on jb, spectral band
      enddo                    

      end subroutine spcvrt_sw

      end module rrtmg_sw_spcvrt



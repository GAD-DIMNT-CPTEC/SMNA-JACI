!
!
!
!
!  radlw_init___rrtmg_lw_ini
!                          |
!                          |__ lwdatinit
!                          |
!                          |__ lwcmbdat       ! g-point interval reduction data
!                          |
!                          |__ lwcldpr        ! cloud optical properties
!                          |
!                          |__ lwatmref       ! reference MLS profile
!                          |
!                          |__ lwavplank      ! Planck function 
!                          |
!                          |__ lw_kgb01       ! molecular absorption coefficients
!                          |
!                          |__ lw_kgb02
!                          |
!                          |__ lw_kgb03
!                          |
!                          |__ lw_kgb04
!                          |
!                          |__ lw_kgb05
!                          |
!                          |__ lw_kgb06
!                          |
!                          |__ lw_kgb07
!                          |
!                          |__ lw_kgb08
!                          |
!                          |__ lw_kgb09
!                          |
!                          |__ lw_kgb10
!                          |
!                          |__ lw_kgb11
!                          |
!                          |__ lw_kgb12
!                          |
!                          |__ lw_kgb13
!                          |
!                          |__ lw_kgb14
!                          |
!                          |__ lw_kgb15
!                          |
!                          |__ lw_kgb16
!                          |
!                          |__ cmbgb1
!                          |
!                          |__ cmbgb2
!                          |
!                          |__ cmbgb3
!                          |
!                          |__ cmbgb4
!                          |
!                          |__ cmbgb5
!                          |
!                          |__ cmbgb6
!                          |
!                          |__ cmbgb7
!                          |
!                          |__ cmbgb8
!                          |
!                          |__ cmbgb9
!                          |
!                          |__ cmbgb10
!                          |
!                          |__ cmbgb11
!                          |
!                          |__ cmbgb12
!                          |
!                          |__ cmbgb13
!                          |
!                          |__ cmbgb14
!                          |
!                          |__ cmbgb15
!                          |
!                          |__ cmbgb16

!
!  rad_rrtmg_lw____
!              |
!              |__ mcica_subcol_lw __ generate_stochastic_clouds __ kissvec
!              |
!              |__ rrtmg_lw __ inatm
!                          |
!                          |__ cldprmc
!                          |
!                          |__ setcoef
!                          |
!                          |__ taumol
!                          |          |
!                          |          |__ taugb1
!                          |          |
!                          |          |__ taugb2
!                          |          |
!                          |          |__ taugb3
!                          |          |
!                          |          |__ taugb4
!                          |          |
!                          |          |__ taugb5
!                          |          |
!                          |          |__ taugb6
!                          |          |
!                          |          |__ taugb7
!                          |          |
!                          |          |__ taugb8
!                          |          |
!                          |          |__ taugb9
!                          |          |
!                          |          |__ taugb10
!                          |          |
!                          |          |__ taugb11
!                          |          |
!                          |          |__ taugb12
!                          |          |
!                          |          |__ taugb13
!                          |          |
!                          |          |__ taugb14
!                          |          |
!                          |          |__ taugb15
!                          |          |
!                          |          |__ taugb16
!                          |
!                          |__ rtrnmc
!
!
!
!
module radlw
!-----------------------------------------------------------------------
! 
! Purpose: Longwave radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,      only: r8 => shr_kind_r8
!use ppgrid,            only: pcols, pver, pverp
!use scamMod,           only: single_column, scm_crm_mode
use parrrtm,           only: nbndlw, ngptlw
use rrtmg_lw_init,     only: rrtmg_lw_ini
use rrtmg_lw_rad,      only: rrtmg_lw
!use spmd_utils,        only: masterproc
!use perf_mod,          only: t_startf, t_stopf
!use cam_logfile,       only: iulog
!use abortutils,        only: endrun
use radconstants,      only: nlwbands

implicit none

private
save

! Public methods

public ::&
   radlw_init,   &! initialize constants
   rad_rrtmg_lw   ! driver for longwave radiation code
   
! Private data
integer :: ntoplw    ! top level to solve for longwave cooling

!===============================================================================
CONTAINS
!===============================================================================

subroutine rad_rrtmg_lw(&
           imca,ncol            ,rrtmg_levs      , &
           r_state_pmidmb  ,r_state_pintmb  ,r_state_tlay    , &
           r_state_tlev    ,r_state_h2ovmr  ,r_state_o3vmr   ,r_state_co2vmr  , &
           r_state_ch4vmr  ,r_state_o2vmr   ,r_state_n2ovmr  ,r_state_cfc11vmr, &
           r_state_cfc12vmr,r_state_cfc22vmr,r_state_ccl4vmr ,       &
           pmid            ,aer_lw_abs      ,cld             ,cicewp_in          , &
           cliqwp_in          ,rei_in             ,rel_in             ,tauc_lw         , &
           qrl             ,qrlc            ,flns            ,flnt            , &
           flnsc           ,flntc           ,flwds           ,flut            , &
           flutc           ,fnl             ,fcnl            ,fldsc           , &
           lu              ,ld              , FeedBackOptics_cld,&
           pcols           ,pver            ,pverp)

!-----------------------------------------------------------------------
!   use cam_history,         only: outfld
   use mcica_subcol_gen_lw, only: mcica_subcol_lw
   use shr_const_mod, only: cpair       => shr_const_cpdair 

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: imca
   integer, intent(in) :: pcols, pver, pverp
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: rrtmg_levs            ! number of levels rad is applied
   CHARACTER(LEN=*), INTENT(IN   ) :: FeedBackOptics_cld

!
! Input arguments which are only passed to other routines
!
    real(r8), intent(in) :: r_state_tlev    (pcols,rrtmg_levs+1)
    real(r8), intent(in) :: r_state_pmidmb  (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_pintmb  (pcols,rrtmg_levs+1)
    real(r8), intent(in) :: r_state_tlay    (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_h2ovmr  (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_o3vmr   (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_co2vmr  (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_ch4vmr  (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_o2vmr   (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_n2ovmr  (pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_cfc11vmr(pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_cfc12vmr(pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_cfc22vmr(pcols,rrtmg_levs)
    real(r8), intent(in) :: r_state_ccl4vmr (pcols,rrtmg_levs)

   real(r8), intent(in) :: pmid(pcols,rrtmg_levs)     ! Level pressure (Pascals)

   real(r8), intent(in) :: aer_lw_abs (pcols,rrtmg_levs,nbndlw) ! aerosol absorption optics depth (LW)

   real(r8), intent(in) :: cld   (pcols,rrtmg_levs)      ! Cloud cover
   real(r8), intent(in) :: cicewp_in(pcols,rrtmg_levs)   ! in-cloud cloud ice water path
   real(r8), intent(in) :: cliqwp_in(pcols,rrtmg_levs)   ! in-cloud cloud liquid water path
   real(r8), intent(in) :: rei_in   (pcols,rrtmg_levs)      ! ice particle effective radius (microns)
   real(r8), intent(in) :: rel_in   (pcols,rrtmg_levs)      ! liquid particle radius (micron)

   real(r8), intent(in) :: tauc_lw(nbndlw,pcols,rrtmg_levs)   ! Cloud longwave optical depth by band
   real(r8) :: cicewp(pcols,rrtmg_levs)   ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,rrtmg_levs)   ! in-cloud cloud liquid water path
   real(r8) :: rei   (pcols,rrtmg_levs)      ! ice particle effective radius (microns)
   real(r8) :: rel   (pcols,rrtmg_levs)      ! liquid particle radius (micron)

!
! Output arguments
!
   real(r8), intent(out) :: qrl (pcols,rrtmg_levs)     ! Longwave heating rate
   real(r8), intent(out) :: qrlc(pcols,rrtmg_levs)     ! Clearsky longwave heating rate
   real(r8), intent(out) :: flns(pcols)          ! Surface cooling flux
   real(r8), intent(out) :: flnt(pcols)          ! Net outgoing flux
   real(r8), intent(out) :: flut(pcols)          ! Upward flux at top of model
   real(r8), intent(out) :: flnsc(pcols)         ! Clear sky surface cooing
   real(r8), intent(out) :: flntc(pcols)         ! Net clear sky outgoing flux
   real(r8), intent(out) :: flutc(pcols)         ! Upward clear-sky flux at top of model
   real(r8), intent(out) :: flwds(pcols)         ! Down longwave flux at surface
   real(r8), intent(out) :: fldsc(pcols)         ! Down longwave clear flux at surface
   real(r8), intent(out) :: fcnl(pcols,rrtmg_levs+1)    ! clear sky net flux at interfaces
   real(r8), intent(out) :: fnl(pcols,rrtmg_levs+1)     ! net flux at interfaces

   real(r8), intent(out)  :: lu( pCols,rrtmg_levs+1,nbndlw) ! longwave spectral flux up
   real(r8), intent(out)  :: ld( pCols,rrtmg_levs+1,nbndlw) ! longwave spectral flux down
   
!
!---------------------------Local variables-----------------------------
!
   integer :: i,j, k, kk, nbnd         ! indices

   real(r8) :: ful(pcols,rrtmg_levs+1)     ! Total upwards longwave flux
   real(r8) :: fsul(pcols,rrtmg_levs+1)    ! Clear sky upwards longwave flux
   real(r8) :: fdl(pcols,rrtmg_levs+1)     ! Total downwards longwave flux
   real(r8) :: fsdl(pcols,rrtmg_levs+1)    ! Clear sky downwards longwv flux

   integer :: inflglw               ! Flag for cloud parameterization method
   integer :: iceflglw              ! Flag for ice cloud param method
   integer :: liqflglw              ! Flag for liquid cloud param method
   integer :: icld                  ! Flag for cloud overlap method
                                 ! 0=clear, 1=random, 2=maximum/random, 3=maximum

   real(r8) :: tsfc(pcols)          ! surface temperature
   real(r8) :: emis(pcols,nbndlw)   ! surface emissivity

   real(r8) :: taua_lw(pcols,rrtmg_levs-1,nbndlw)     ! aerosol optical depth by band

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day

   ! Cloud arrays for McICA 
   integer, parameter :: nsubclw = ngptlw       ! rrtmg_lw g-point (quadrature point) dimension
   integer :: permuteseed                       ! permute seed for sub-column generator


   real(r8) :: cld_stolw    (nsubclw, pcols, rrtmg_levs-1)     ! cloud fraction (mcica)
   real(r8) :: cicewp_stolw (nsubclw, pcols, rrtmg_levs-1)  ! cloud ice water path (mcica)
   real(r8) :: cliqwp_stolw (nsubclw, pcols, rrtmg_levs-1)  ! cloud liquid water path (mcica)
   real(r8) :: rei_stolw    (pcols,rrtmg_levs-1)               ! ice particle size (mcica)
   real(r8) :: rel_stolw    (pcols,rrtmg_levs-1)               ! liquid particle size (mcica)
   real(r8) :: tauc_stolw   (nsubclw, pcols, rrtmg_levs-1)    ! cloud optical depth (mcica - optional)

   ! Includes extra layer above model top
   real(r8) :: uflx   (pcols,rrtmg_levs+1)  ! Total upwards longwave flux
   real(r8) :: uflxc  (pcols,rrtmg_levs+1) ! Clear sky upwards longwave flux
   real(r8) :: dflx   (pcols,rrtmg_levs+1)  ! Total downwards longwave flux
   real(r8) :: dflxc  (pcols,rrtmg_levs+1) ! Clear sky downwards longwv flux
   real(r8) :: hr     (pcols,rrtmg_levs)      ! Longwave heating rate (K/d)
   real(r8) :: hrc    (pcols,rrtmg_levs)     ! Clear sky longwave heating rate (K/d)
   real(r8) :: lwuflxs(nbndlw,pcols,rrtmg_levs+1)  ! Longwave spectral flux up
   real(r8) :: lwdflxs(nbndlw,pcols,rrtmg_levs+1)  ! Longwave spectral flux down
   !-----------------------------------------------------------------------
   INTEGER :: nmca
   INTEGER :: ims

   ! mji/rrtmg
   qrl  = 0.0_r8  ;   qrlc = 0.0_r8  ;   flns = 0.0_r8  ;   flnt = 0.0_r8  ;
   flut = 0.0_r8  ;   flnsc = 0.0_r8 ;   flntc = 0.0_r8 ;   flutc = 0.0_r8 ;
   flwds = 0.0_r8 ;   fldsc = 0.0_r8 ;   fcnl = 0.0_r8  ;   fnl = 0.0_r8   ;
   lu = 0.0_r8    ;   ld = 0.0_r8    ;
   
!
!---------------------------Local variables-----------------------------
!

   ful = 0.0_r8  ;            fsul = 0.0_r8  ;           fdl = 0.0_r8        ;      fsdl = 0.0_r8   ;       
   tsfc = 0.0_r8 ;            emis = 0.0_r8  ;           taua_lw = 0.0_r8       ;


   ! Cloud arrays for McICA 

    cicewp = cicewp_in!0.0_r8             
    cliqwp = cliqwp_in! 0.0_r8             
    rei = rei_in    !0.0_r8                
    rel = rel_in    !0.0_r8                

    cld_stolw  = 0.0_r8      
    cicewp_stolw  = 0.0_r8      
    cliqwp_stolw  = 0.0_r8      
    rei_stolw  = 0.0_r8      
    rel_stolw  = 0.0_r8      
    tauc_stolw    = 0.0_r8      

   ! Includes extra layer above model top
   uflx    = 0.0_r8             
   uflxc   = 0.0_r8             
   dflx    = 0.0_r8             
   dflxc   = 0.0_r8             
   hr   = 0.0_r8             
   hrc     = 0.0_r8             
   lwuflxs = 0.0_r8             
   lwdflxs = 0.0_r8             
   !-----------------------------------------------------------------------

   ! Calculate cloud optical properties here if using CAM method, or if using one of the
   ! methods in RRTMG_LW, then pass in cloud physical properties and zero out cloud optical 
   ! properties here
   
   ! Zero optional cloud optical depth input array tauc_lw, 
   ! if inputting cloud physical properties into RRTMG_LW
   !          tauc_lw(:,:,:) = 0.
   ! Or, pass in CAM cloud longwave optical depth to RRTMG_LW
   ! do nbnd = 1, nbndlw
   !    tauc_lw(nbnd,:ncol,:pver) = cldtau(:ncol,:pver)
   ! end do

   ! Call mcica sub-column generator for RRTMG_LW

   ! Call sub-column generator for McICA in radiation
   !call t_startf('mcica_subcol_lw')

   ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
   icld = 2
   ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
   ! effective randomization)
   permuteseed = 150

   ! These fields are no longer supplied by CAM.
 !  cicewp = 0.0_r8
 !  cliqwp = 0.0_r8
 !  rei = 0.0_r8
 !  rel = 0.0_r8

! Set nmca to sample size for Monte Carlo calculation
   IF(TRIM(FeedBackOptics_cld) == 'WRF')THEN
      if (imca.eq.0) nmca = 1
      if (imca.eq.1) nmca = 1!200
   ELSE
      if (imca.eq.0) nmca = 1
      if (imca.eq.1) nmca = 2!200
   END IF
! This is the statistical sampling loop for McICA

   do ims = 1, nmca


      ! Call sub-colum cloud generator for McICA calculations
      ! Output will be written for all nmca samples.  This will be excessive if
      ! band output (iout=99) option is selected. 

      if (imca.eq.1) then
         call mcica_subcol_lw( ncol                                                 , & !integer, intent(in) :: ncol        ! number of columns
                               rrtmg_levs-1                                         , & !integer, intent(in) :: nlay        ! number of model layers
                               icld                                                 , & !integer, intent(in) :: icld        ! clear/cloud, cloud overlap flag
                               permuteseed                                          , & !integer, intent(in) :: permuteseed    ! if the cloud generator is called multiple times, 
                               pmid        (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(in) :: play(ncol,nlay)           ! layer pressures (mb) 
                               cld         (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(in) :: cldfrac(ncol,nlay)       ! layer cloud fraction
                               cicewp      (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(in) :: ciwp (ncol,nlay)      ! cloud ice water path
                               cliqwp      (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(in) :: clwp (ncol,nlay)      ! cloud liquid water path
                               rei         (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(in) :: rei (ncol,nlay)        ! cloud ice particle size
                               rel         (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(in) :: rel (ncol,nlay)           ! cloud liquid particle size
                               tauc_lw     (1:nbndlw ,1: pcols,1: rrtmg_levs-1)     , & !real(kind=r8), intent(in) :: tauc (nbndlw,ncol,nlay)         ! cloud optical depth
                               cld_stolw   (1:nsubclw,1: pcols,1: rrtmg_levs-1)     , & !real(kind=r8), intent(out) :: cldfmcl(ngptlw,ncol,nlay)    ! cloud fraction [mcica]
                               cicewp_stolw(1:nsubclw,1: pcols,1: rrtmg_levs-1)     , & !real(kind=r8), intent(out) :: ciwpmcl(ngptlw,ncol,nlay)    ! cloud ice water path [mcica]
                               cliqwp_stolw(1:nsubclw,1: pcols,1: rrtmg_levs-1)     , & !real(kind=r8), intent(out) :: clwpmcl(ngptlw,ncol,nlay)    ! cloud liquid water path [mcica]
                               rei_stolw   (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(out) :: reicmcl (ncol,nlay)    ! ice partcle size (microns)
                               rel_stolw   (1:pcols,1:rrtmg_levs-1)                 , & !real(kind=r8), intent(out) :: relqmcl (ncol,nlay)    ! liquid particle size (microns)
                               tauc_stolw  (1:nsubclw,1: pcols,1: rrtmg_levs-1)       ) !real(kind=r8), intent(out) :: taucmcl(ngptlw,ncol,nlay)    ! cloud optical depth [mcica]
      endif

   !call t_stopf('mcica_subcol_lw')

   
   !call t_startf('rrtmg_lw')

   !
   ! Call RRTMG_LW model
   !
   ! Set input flags for cloud parameterizations
   ! Use separate specification of ice and liquid cloud optical depth.
   ! Use either Ebert and Curry ice parameterization (iceflglw = 0 or 1), 
   ! or use Key (Streamer) approach (iceflglw = 2), or use Fu method
   ! (iceflglw = 3), and Hu/Stamnes for liquid (liqflglw = 1).
   ! For use in Fu method (iceflglw = 3), rei is converted in RRTMG_LW
   ! from effective radius to generalized effective size using the
   ! conversion of D. Mitchell, JAS, 2002.  For ice particles outside
   ! the effective range of either the Key or Fu approaches, the 
   ! Ebert and Curry method is applied. 
   IF(TRIM(FeedBackOptics_cld) == 'WRF')THEN
   !   !WRFMODEL ! For passing in cloud physical properties; cloud optics parameterized in RRTMG:
   !   !WRFMODEL       icld = 2
   !   !WRFMODEL       inflglw = 2
   !   !WRFMODEL       iceflglw = 3
   !   !WRFMODEL       liqflglw = 1
       inflglw = 0
       iceflglw = 0
       liqflglw = 0
   ELSE IF(TRIM(FeedBackOptics_cld) == 'CAM5')THEN
      inflglw = 0
      iceflglw = 0
      liqflglw = 0
   ELSE
      ! IF(ims == 10)THEN
      !   ! Input CAM cloud optical depth directly
      !   !PK inflglw = 0
      !   !PK iceflglw = 0
      !   !PK liqflglw = 0
      inflglw = 0
      iceflglw = 0
      liqflglw = 0
       !ELSE
       !   ! Use E&C approach for ice to mimic CAM3
       !   !   inflglw = 2
       !   !   iceflglw = 1
       !   !   liqflglw = 1
       !   ! Use merged Fu and E&C params for ice
       !   !   inflglw = 2
       !   !   iceflglw = 3
       !   !   liqflglw = 1
       !   !WRFMODEL ! For passing in cloud physical properties; cloud optics parameterized in RRTMG:
       !   !WRFMODEL       icld = 2
       !   !WRFMODEL       inflglw = 2
       !   !WRFMODEL       iceflglw = 3
       !   !WRFMODEL       liqflglw = 1
       !    inflglw = 2
       !    iceflglw = 3
       !    liqflglw = 1
       !END IF 
   END IF
   ! Convert incoming water amounts from specific humidity to vmr as needed;
   ! Convert other incoming molecular amounts from mmr to vmr as needed;
   ! Convert pressures from Pa to hPa;
   ! Set surface emissivity to 1.0 here, this is treated in land surface model;
   ! Set surface temperature
   ! Set aerosol optical depth to zero for now

   emis(:ncol,:nbndlw) = 1._r8
   tsfc(:ncol) = r_state_tlev(:ncol,rrtmg_levs+1)
   taua_lw(:ncol, 1:rrtmg_levs-1, :nbndlw) = aer_lw_abs(:ncol,pverp-rrtmg_levs+1:pverp-1,:nbndlw)

    lu(1:ncol,:,:) = 0.0_r8
    ld(1:ncol,:,:) = 0.0_r8

   call rrtmg_lw(&
        imca                                              , & !integer      , intent(in   ) :: imca
        pcols                                             , & !integer      , intent(in   ) :: pcols
        ncol                                              , & !integer      , intent(in   ) :: ncol                   ! Number of horizontal columns
        rrtmg_levs                                        , & !integer      , intent(in   ) :: nlay                   ! Number of model layers
        icld                                              , & !integer      , intent(inout) :: icld                   ! Cloud overlap method
        r_state_pmidmb  (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: play     (ncol,nlay)   ! Layer pressures (hPa, mb)
        r_state_pintmb  (1:pcols,1:rrtmg_levs+1)          , & !real(kind=r8), intent(in   ) :: plev     (ncol,nlay+1) ! Interface pressures (hPa, mb)
        r_state_tlay    (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: tlay     (ncol,nlay)   ! Layer temperatures (K)
        r_state_tlev    (1:pcols,1:rrtmg_levs+1)          , & !real(kind=r8), intent(in   ) :: tlev     (ncol,nlay+1)    ! Interface temperatures (K)
        tsfc            (1:pcols)                         , & !real(kind=r8), intent(in   ) :: tsfc     (ncol)        ! Surface temperature (K)
        r_state_h2ovmr  (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: h2ovmr   (ncol,nlay)     ! H2O volume mixing ratio
        r_state_o3vmr   (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: o3vmr    (ncol,nlay)      ! O3 volume mixing ratio
        r_state_co2vmr  (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: co2vmr   (ncol,nlay)     ! CO2 volume mixing ratio
        r_state_ch4vmr  (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: ch4vmr   (ncol,nlay)         ! Methane volume mixing ratio
        r_state_o2vmr   (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: o2vmr    (ncol,nlay)          ! O2 volume mixing ratio
        r_state_n2ovmr  (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: n2ovmr   (ncol,nlay)     ! Nitrous oxide volume mixing ratio
        r_state_cfc11vmr(1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: cfc11vmr (ncol,nlay)   ! CFC11 volume mixing ratio
        r_state_cfc12vmr(1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: cfc12vmr (ncol,nlay)   ! CFC12 volume mixing ratio
        r_state_cfc22vmr(1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: cfc22vmr (ncol,nlay)   ! CFC22 volume mixing ratio
        r_state_ccl4vmr (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: ccl4vmr  (ncol,nlay)    ! CCL4 volume mixing ratio
        emis            (1:pcols,1:nbndlw)                , & !real(kind=r8), intent(in   ) :: emis     (ncol,nbndlw)     ! Surface emissivity
        inflglw                                           , & !integer      , intent(in   ) :: inflglw                          ! Flag for cloud optical properties
        iceflglw                                          , & !integer      , intent(in   ) :: iceflglw                         ! Flag for ice particle specification
        liqflglw                                          , & !integer      , intent(in   ) :: liqflglw                         ! Flag for liquid droplet specification
        cld_stolw       (1:nsubclw,1:pcols,1:rrtmg_levs-1), & !real(kind=r8), intent(in   ) :: cldfmcl (ngptlw,ncol,nlay-1)       ! Cloud fraction
        tauc_stolw      (1:nsubclw,1:pcols,1:rrtmg_levs-1), & !real(kind=r8), intent(in   ) :: taucmcl (ngptlw,ncol,nlay-1)       ! Cloud optical depth
        cicewp_stolw    (1:nsubclw,1:pcols,1:rrtmg_levs-1), & !real(kind=r8), intent(in   ) :: ciwpmcl (ngptlw,ncol,nlay-1)       ! Cloud ice water path (g/m2)
        cliqwp_stolw    (1:nsubclw,1:pcols,1:rrtmg_levs-1), & !real(kind=r8), intent(in   ) :: clwpmcl (ngptlw,ncol,nlay-1)      ! Cloud liquid water path (g/m2)
        cld             (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(in   ) :: cld     (ncol,nlay)    ! layer cloud fraction
        cicewp          (1:pcols,1:rrtmg_levs-1)          , & !real(kind=r8), intent(in   ) :: cicewp  (ncol,nlay-1) 
        cliqwp          (1:pcols,1:rrtmg_levs-1)          , & !real(kind=r8), intent(in   ) :: cliqwp  (ncol,nlay-1)  
        tauc_lw         (1:nbndlw,1:pcols,1:rrtmg_levs)   , & !real(kind=r8), intent(in   ) :: tauc_lw (nbndlw,ncol,nlay)        ! cloud optical depth
        rei             (1:pcols,1:rrtmg_levs-1)          , & !real(kind=r8), intent(in   ) :: reicmcl (ncol,nlay-1)     ! Cloud ice effective radius (microns)
        rel             (1:pcols,1:rrtmg_levs-1)          , & !real(kind=r8), intent(in   ) :: relqmcl (ncol,nlay-1)     ! Cloud water drop effective radius (microns)
        taua_lw         (1:pcols,1:rrtmg_levs-1,1:nbndlw) , & !real(kind=r8), intent(in   ) :: tauaer  (ncol,nlay-1,nbndlw)        ! aerosol optical depth
        uflx            (1:pcols,1:rrtmg_levs+1)          , & !real(kind=r8), intent(inout) :: uflx    (ncol,nlay+1)      ! Total sky longwave upward flux (W/m2)
        dflx            (1:pcols,1:rrtmg_levs+1)          , & !real(kind=r8), intent(inout) :: dflx    (ncol,nlay+1)      ! Total sky longwave downward flux (W/m2)
        hr              (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(inout) :: hr      (ncol,nlay)        ! Total sky longwave radiative heating rate (K/d)
        uflxc           (1:pcols,1:rrtmg_levs+1)          , & !real(kind=r8), intent(inout) :: uflxc   (ncol,nlay+1)        ! Clear sky longwave upward flux (W/m2)
        dflxc           (1:pcols,1:rrtmg_levs+1)          , & !real(kind=r8), intent(inout) :: dflxc   (ncol,nlay+1)       ! Clear sky longwave downward flux (W/m2)
        hrc             (1:pcols,1:rrtmg_levs)            , & !real(kind=r8), intent(inout) :: hrc     (ncol,nlay)      ! Clear sky longwave radiative heating rate (K/d)
        lwuflxs         (1:nbndlw,1:pcols,1:rrtmg_levs+1) , & !real(kind=r8), intent(inout) :: uflxs   (nbndlw,ncol,nlay+1)    ! Total sky longwave upward flux spectral (W/m2)
        lwdflxs         (1:nbndlw,1:pcols,1:rrtmg_levs+1)   ) !real(kind=r8), intent(inout) :: dflxs   (nbndlw,ncol,nlay+1)    ! Total sky longwave downward flux spectral (W/m2)
!      subroutine rrtmg_lw(&
!       pcols           ,ncol            ,nlay            ,icld            , &
!       play            ,plev            ,tlay            ,tlev            , &
!       tsfc            ,h2ovmr          ,o3vmr           ,co2vmr          , &
!       ch4vmr          ,o2vmr           ,n2ovmr          ,cfc11vmr        , &
!       cfc12vmr        ,cfc22vmr        ,ccl4vmr         ,emis            , &
!       inflglw         ,iceflglw        ,liqflglw        ,cldfmcl         , &
!       taucmcl         ,ciwpmcl         ,clwpmcl         ,reicmcl         , &
!       relqmcl         ,tauaer          ,uflx            ,dflx            , &
!       hr              ,uflxc           ,dflxc           ,hrc             , &
!       uflxs           ,dflxs            )

            if (imca .eq. 0) then
               if (ims .eq. nmca) then
                  DO k=1,rrtmg_levs+1
                     DO i=1,ncol
                        uflx     (i,k) = uflx     (i,k)  /nmca
                        dflx     (i,k) = dflx     (i,k)  /nmca
                        uflxc    (i,k) = uflxc    (i,k)  /nmca
                        dflxc    (i,k) = dflxc    (i,k)  /nmca
                     END DO
                  END DO 
                  DO k=1,rrtmg_levs
                     DO i=1,ncol
                         hr          (i,k) = hr     (i,k)  /nmca
                         hrc         (i,k) = hrc    (i,k)  /nmca
                     END DO
                  END DO 
                  DO k=1,rrtmg_levs+1
                     DO i=1,ncol
                        DO j=1,nbndlw
                           lwuflxs       (j,i,k)    =lwuflxs    (j,i,k)/nmca
                           lwdflxs       (j,i,k)    =lwdflxs    (j,i,k)/nmca
                        END DO
                      END DO
                  END DO 
               endif
            elseif (imca .eq. 1) then
               if (ims .eq. nmca) then
                  DO k=1,rrtmg_levs+1
                     DO i=1,ncol
                        uflx     (i,k) = uflx     (i,k)  /nmca
                        dflx     (i,k) = dflx     (i,k)  /nmca
                        uflxc    (i,k) = uflxc    (i,k)  /nmca
                        dflxc    (i,k) = dflxc    (i,k)  /nmca
                     END DO
                  END DO 
                  DO k=1,rrtmg_levs
                     DO i=1,ncol
                         hr          (i,k) = hr     (i,k)  /nmca
                         hrc         (i,k) = hrc    (i,k)  /nmca
                     END DO
                  END DO 
                  DO k=1,rrtmg_levs+1
                     DO i=1,ncol
                        DO j=1,nbndlw
                           lwuflxs       (j,i,k)    =lwuflxs    (j,i,k)/nmca
                           lwdflxs       (j,i,k)    =lwdflxs    (j,i,k)/nmca
                        END DO
                      END DO
                  END DO 
               endif
            endif

! End statistical loop for McICA

  END DO 
   !
   !----------------------------------------------------------------------
   ! All longitudes: store history tape quantities
   ! Flux units are in W/m2 on output from rrtmg_lw and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are in K/d on output from RRTMG and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are converted to J/kg/s below for use in CAM. 

   flwds(:ncol) = dflx (:ncol,1)
   fldsc(:ncol) = dflxc(:ncol,1)
   flns(:ncol)  = uflx (:ncol,1) - dflx (:ncol,1)
   flnsc(:ncol) = uflxc(:ncol,1) - dflxc(:ncol,1)
   flnt(:ncol)  = uflx (:ncol,rrtmg_levs) - dflx (:ncol,rrtmg_levs)
   flntc(:ncol) = uflxc(:ncol,rrtmg_levs) - dflxc(:ncol,rrtmg_levs)
   flut(:ncol)  = uflx (:ncol,rrtmg_levs)
   flutc(:ncol) = uflxc(:ncol,rrtmg_levs)

   !
   ! Reverse vertical indexing here for CAM arrays to go from top to bottom.
   !
   ful = 0._r8
   fdl = 0._r8
   fsul = 0._r8
   fsdl = 0._r8
   ful (:ncol,pverp-rrtmg_levs+1:pverp)= uflx(:ncol,rrtmg_levs:1:-1)
   fdl (:ncol,pverp-rrtmg_levs+1:pverp)= dflx(:ncol,rrtmg_levs:1:-1)
   fsul(:ncol,pverp-rrtmg_levs+1:pverp)=uflxc(:ncol,rrtmg_levs:1:-1)
   fsdl(:ncol,pverp-rrtmg_levs+1:pverp)=dflxc(:ncol,rrtmg_levs:1:-1)

!   if (single_column.and.scm_crm_mode) then
!      call outfld('FUL     ',ful,pcols,lchnk)
!      call outfld('FDL     ',fdl,pcols,lchnk)
!      call outfld('FULC    ',fsul,pcols,lchnk)
!      call outfld('FDLC    ',fsdl,pcols,lchnk)
!   endif
   
   fnl(:ncol,:) = ful(:ncol,:) - fdl(:ncol,:)
   ! mji/ cam excluded this?
   fcnl(:ncol,:) = fsul(:ncol,:) - fsdl(:ncol,:)

   ! Pass longwave heating to CAM arrays and convert from K/d to J/kg/s
  ! qrl = 0._r8
  ! qrlc = 0._r8
  ! qrl (:ncol,pverp-rrtmg_levs+1:pver)=hr (:ncol,rrtmg_levs-1:1:-1)*cpair*dps
  ! qrlc(:ncol,pverp-rrtmg_levs+1:pver)=hrc(:ncol,rrtmg_levs-1:1:-1)*cpair*dps
  
   ! Pass longwave heating to CAM arrays and convert from K/d to K/s
   qrl = 0._r8
   qrlc = 0._r8
   qrl (:ncol,pverp-rrtmg_levs+1:pver)=hr (:ncol,rrtmg_levs-1:1:-1)/86400.0_r8 !*cpair*dps
   qrlc(:ncol,pverp-rrtmg_levs+1:pver)=hrc(:ncol,rrtmg_levs-1:1:-1)/86400.0_r8 !*cpair*dps

   ! Return 0 above solution domain
   if ( ntoplw > 1 )then
      qrl(:ncol,:ntoplw-1) = 0._r8
      qrlc(:ncol,:ntoplw-1) = 0._r8
   end if

   ! Pass spectral fluxes, reverse layering
   ! order=(/3,1,2/) maps the first index of lwuflxs to the third index of lu.
   !if (associated(lu)) then
      lu(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwuflxs(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   !end if
   
   !if (associated(ld)) then
      ld(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwdflxs(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   !end if
   
   !call t_stopf('rrtmg_lw')

end subroutine rad_rrtmg_lw

!-------------------------------------------------------------------------------

subroutine radlw_init(pref_mid,pver)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme.
!
!-----------------------------------------------------------------------

!   use ref_pres, only : pref_mid
   integer, intent(in) :: pver
   real(r8), intent(in) :: pref_mid(:)  ! reference pressure at layer midpoints (Pa)

   integer :: k
!   do k = 1, plev
!      pref_edge(k) = hypi(k)
!      pref_mid(k)  = hypm(k)
!   end do

   ! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
   ! longwave cooling to about 80 km (1 Pa)
   if (pref_mid(1) .lt. 0.1_r8) then
      do k = 1, pver
         if (pref_mid(k) .lt. 1._r8) ntoplw  = k
      end do
   else
      ntoplw  = 1
   end if
!   if (masterproc) then
!      write(iulog,*) 'radlw_init: ntoplw =',ntoplw
!   endif

   call rrtmg_lw_ini

end subroutine radlw_init

!-------------------------------------------------------------------------------

end module radlw
